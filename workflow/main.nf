#!/usr/bin/env nextflow
/*
 * AAFTF internal benchmarking suite (see aaftf_benchmark_specification.md)
 *
 *   1. Normalise reference headers (FUNGUS_*, CONTAM_<taxon>_*) and build one combined FASTA per sample
 *      (sample + all contaminants, data/sim/references/combined/<sample>_with_contaminants.fa)
 *   2. Simulate reads at 0/1/5/15% contaminant spike-in: InSilicoSeq Illumina and/or Badread ONT (as --tracks needs)
 *   3. Ablation checkpoints per dataset (track 'short'):
 *        cp00_raw_baseline        spades.py --isolate on raw reads
 *        cp01_spades_asm          AAFTF trim → mito → filter (CLEAN_READS), then assemble with SPAdes
 *        cp02_vecscreen           AAFTF vecscreen      (comparison branch, stops here: cp02_fcs_screen, AAFTF fcs_screen)
 *        cp03_taxpurge            AAFTF taxpurge       (comparison branch, stops here: cp03_fcsgx, AAFTF fcs_gx_purge)
 *        cp04_covpurge            AAFTF covpurge (filtered Illumina reads)
 *        cp05_rmdup               AAFTF rmdup
 *      --tracks (short,long,hybrid) picks which tracks run per sample × spike level (ONT assembled with AAFTF assemble --method flye):
 *        long:   cp01_flye_asm → cp02 → cp03 → cp05 → cp06_polished (racon); no cp04 (no Illumina reads)
 *        hybrid: cp01-cp03 shared with 'long'; then cp04 covpurge and cp05 rmdup with its own Illumina reads,
 *                cp06_polished = racon, then polypolish with the filtered Illumina reads of the same dataset
 * Each checkpoint FASTA is named <dataset>.<checkpoint>.fasta; one end-to-end path per dataset (vecscreen → taxpurge).
 *   4. Evaluate every checkpoint: QUAST, BUSCO (miniprot), minimap2 vs the combined FASTA,
 *      minimap2 PAF vs the reference for D-GENIES dot plots
 *   5. Aggregate into results/final_report.tsv
 * Every AAFTF / assembler call is wrapped in GNU time (scripts/timeit.sh).
 */
nextflow.enable.dsl = 2

params.sim_samples = 'config/sim_samples.tsv'
params.contaminants = 'config/contaminants.tsv'
params.outdir = 'results'
// simulated datasets (references, reads) and real datasets (same-strain references, sequenced reads)
params.sim_dir = 'data/sim'
params.real_dir = 'data/real'
// real datasets (same-strain reference + sequenced reads); skipped when the file does not exist
params.real_samples = 'config/real_samples.tsv'
params.short_coverage = 50
params.read_length = 150
params.iss_model = 'novaseq'
params.spike_levels = '0,1,5,15'
params.group_weights = 'bacteria:0.65,host:0.15,vector:0.1,fungal:0.1'
// any of: short,long,hybrid
params.tracks = 'short'
// use reads already in <sim_dir>/reads when their .sim_hash matches
params.reuse_reads = true
params.badread_args = '--length 15000,8000 --error_model nanopore2023 --qscore_model nanopore2023 --identity 98.5,99.8,1 --chimeras 1.0 --glitches 10000,25,25'
params.long_coverage = 35
params.longread_type = 'nano-hq'
params.trim_method = 'bbduk'
params.filter_aligner = 'bbduk'
params.min_contig = 500
// covpurge: remove contigs below this % of the N50 contigs' mean coverage (AAFTF default 5)
params.covpurge_mincovpct = 5
params.fcs_gx_db = null
// NCBI FCS-adaptor (AAFTF fcs_screen) as a comparison branch to vecscreen at cp02 (needs Singularity)
params.fcs_screen = true
// BUSCO download folder (<dir>/lineages/<lineage>; read-only is fine, used with --offline)
params.busco_downloads = null

// (sample, spike level) links the short, long and hybrid datasets built from the same simulation
// file/header-safe key for an NCBI organism name (contaminant taxa and sample names), e.g. 'Cutibacterium acnes subsp. acnes' -> 'Cutibacterium_acnes_subsp_acnes'
def safeKey(name) {
    name.trim().replaceAll(/[^A-Za-z0-9]+/, '_').replaceAll(/^_|_$/, '')
}

// one meta map + reference per row of sim_samples.tsv (origin 'sim') or real_samples.tsv (origin 'real')
// AAFTF commit pinned in pixi.lock. Referenced as params.aaftf_commit in every AAFTF job script, so changing
// the pin (pixi.toml + pixi lock) reruns exactly the AAFTF jobs and what depends on them on the next -resume.
def aaftfCommit() {
    def lock = file("${projectDir}/../pixi.lock")
    def m = lock.exists() ? (lock.text =~ /AAFTF\.git[^#\s]*#([0-9a-f]{40})/) : null
    m ? m[0][1] : 'unknown'
}

// TSV rows whose first column (organism) starts with '#' are commented out and skipped
def activeRow(r) {
    !(r.organism ?: '').trim().startsWith('#')
}

def optFile(x) {
    x ? file(x) : null
}

def sampleRow(r, origin) {
    def meta = [
        sample: safeKey(r.organism), name: r.organism, origin: origin, genome_size: r.genome_size, phylum: r.phylum,
        taxid: r.taxid, lineage: r.lineage, mito_seed: optFile(r.mito_seed), source_reference: r.reference,
        accessions: r.accessions ?: '', read1: optFile(r.read1), read2: optFile(r.read2), longreads: optFile(r.longreads), longread_type: r.longread_type ?: null,
    ]
    tuple(meta, file(r.reference))
}

// order- and name-independent fingerprint of an assembly: md5 of the sorted per-contig sequence md5s
def seqFingerprint(f) {
    def records = file(f).text.split(/(?m)^>/).findAll { r -> r.trim() }
    def digests = records.collect { r -> r.substring(r.indexOf('\n') + 1).replaceAll(/\s+/, '').toUpperCase().md5() }
    digests.sort().join(',').md5()
}

def trackLinkKey(m) {
    [m.sample, m.spike_pct]
}
def toHybridMeta(m) {
    m + [id: m.id.replaceAll('_long$', '_hybrid'), track: 'hybrid']
}

// Fingerprint of everything that determines a simulated read set: input FASTAs (path, size,
// mtime, as Nextflow's own cache does), the contaminant table, the simulation params and the
// helper scripts. Stored as <dataset>_<kind>.sim_hash next to the reads.
def fileSig(f) {
    def p = file(f)
    p.exists() ? "${p.toAbsolutePath()}:${p.size()}:${p.lastModified()}" : "${f}:missing"
}

def simHash(m, pct, kind) {
    // contaminant FASTAs: the `reference` column, looked up by name (column order may change)
    def lines = file(params.contaminants).readLines().findAll { l -> l.trim() && !l.trim().startsWith('#') }
    def ref_col = lines[0].split('\t', -1).toList().indexOf('reference')
    def contam_refs = lines.drop(1).collect { l -> l.split('\t', -1)[ref_col] }
    def parts = [kind, pct, fileSig(m.source_reference), fileSig(params.contaminants)] + contam_refs.collect { f -> fileSig(f) } + [params.group_weights, file("${projectDir}/scripts/format_reference_headers.py").text]
    parts += kind == 'short'
        ? [params.short_coverage, params.read_length, params.iss_model, file("${projectDir}/scripts/make_readcounts.py").text]
        : [params.long_coverage, params.badread_args, file("${projectDir}/scripts/make_long_quantities.py").text]
    parts.join('|').md5()
}

// previously simulated files for a sample × spike level ('short' or 'long'); [] unless all exist
// and were made with the same settings (matching .sim_hash)
def reusableReads(m, pct, kind) {
    def d = "${params.sim_dir}/reads/${m.sample}_spike${pct}"
    def fs = kind == 'short'
        ? [file("${d}_R1.fastq.gz"), file("${d}_R2.fastq.gz"), file("${d}_short.composition.tsv")]
        : [file("${d}_ont.fastq.gz"), file("${d}_long.composition.tsv")]
    def h = file("${d}_${kind}.sim_hash")
    def ok = params.reuse_reads.toString() == 'true' && fs.every { f -> f.exists() && f.size() > 0 } && h.exists() && h.text.trim() == simHash(m, pct, kind)
    ok ? fs : []
}

// ================================================================ data generation

params.aaftf_commit = aaftfCommit()
// pixi manifest, used by TOOL_VERSIONS to query the QUAST / BUSCO environments
params.pixi_manifest = "${projectDir}/../pixi.toml"

process TAG_FUNGAL_REF {
    tag "${meta.sample}"
    publishDir path: { meta.origin == 'real' ? "${params.real_dir}/references" : "${params.sim_dir}/references/fungus" }, mode: 'copy'

    input:
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("${meta.sample}.fa")

    script:
    """
    format_reference_headers.py ${ref} --prefix FUNGUS --accessions '${meta.accessions}' > ${meta.sample}.fa
    """
}

process TAG_CONTAM_REF {
    tag "${taxon}"
    publishDir "${params.sim_dir}/references/contaminants", mode: 'copy'

    input:
    tuple val(taxon), val(group), val(accessions), path(ref)

    output:
    tuple val(taxon), val(group), path("CONTAM_${taxon}.fa")

    script:
    """
    format_reference_headers.py ${ref} --prefix CONTAM_${taxon} \\
        --accessions '${accessions}' > CONTAM_${taxon}.fa
    """
}

// sample genome + all contaminants: minimap2 labels each contig FUNGUS_ / CONTAM_ by its best match
process BUILD_COMBINED_REF {
    tag "${meta.sample}"
    publishDir "${params.sim_dir}/references/combined", mode: 'copy'

    input:
    tuple val(meta), path(fungus), path(contams)

    output:
    tuple val(meta), path("${meta.sample}_with_contaminants.fa")

    script:
    """
    cat ${fungus} ${contams} > ${meta.sample}_with_contaminants.fa
    """
}

process SIMULATE_ILLUMINA {
    tag "${meta.sample}:${pct}%"
    publishDir "${params.sim_dir}/reads", mode: 'copy'

    input:
    tuple val(meta), path(fungus), val(pct), val(groups), path(contams), val(hash)

    output:
    tuple val(meta), val(pct), path("${meta.sample}_spike${pct}_R{1,2}.fastq.gz"), emit: reads
    path ("${meta.sample}_spike${pct}_short.composition.tsv"), emit: composition
    path "${meta.sample}_spike${pct}_short.sim_hash"

    script:
    def ds = "${meta.sample}_spike${pct}"
    def ds_track = "${ds}_short"
    """
    make_readcounts.py --fungus ${fungus} --contam ${contams} --groups ${groups.join(' ')} \\
        --weights ${params.group_weights} --coverage ${params.short_coverage} \\
        --read_length ${params.read_length} --percent ${pct} \\
        --out_fungus fungus.readcount.tsv --out_contam contam.readcount.tsv

    iss generate --genomes ${fungus} --readcount_file fungus.readcount.tsv --model ${params.iss_model} \\
        --output fungus --cpus ${task.cpus} --seed 42 --compress
    cp fungus_R1.fastq.gz ${ds}_R1.fastq.gz; cp fungus_R2.fastq.gz ${ds}_R2.fastq.gz
    if [ -s contam.readcount.tsv ]; then
        cat ${contams} > contam_all.fa
        iss generate --genomes contam_all.fa --readcount_file contam.readcount.tsv \\
            --model ${params.iss_model} --output contam --cpus ${task.cpus} --seed 43 --compress
        cat contam_R1.fastq.gz >> ${ds}_R1.fastq.gz
        cat contam_R2.fastq.gz >> ${ds}_R2.fastq.gz
    fi

    # realised ground truth (read pairs per source)
    nf=\$(( \$(zcat fungus_R1.fastq.gz | wc -l) / 4 ))
    nt=\$(( \$(zcat ${ds}_R1.fastq.gz | wc -l) / 4 ))
    # never hand an empty read set to the assembly steps
    [ "\$nf" -gt 0 ] || { echo "ERROR: InSilicoSeq produced no fungal reads" >&2; exit 1; }
    [ "${pct}" = "0" ] || [ "\$nt" -gt "\$nf" ] || { echo "ERROR: no contaminant reads at ${pct}%" >&2; exit 1; }
    awk -v ds=${ds_track} -v nf=\$nf -v nt=\$nt 'BEGIN { printf "dataset\\tfungal_pairs\\tcontam_pairs\\tcontam_pct\\n%s\\t%d\\t%d\\t%.3f\\n", ds, nf, nt-nf, (nt > 0 ? (nt-nf)*100/nt : 0) }' \\
        > ${ds_track}.composition.tsv
    echo ${hash} > ${ds_track}.sim_hash
    """
}

process SIMULATE_ONT_FUNGUS {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(fungus)

    output:
    tuple val(meta), path("${meta.sample}_ont_fungus.fastq.gz")

    script:
    """
    badread simulate --reference ${fungus} --quantity ${params.long_coverage}x \\
        ${params.badread_args} --seed 42 \\
        | gzip > ${meta.sample}_ont_fungus.fastq.gz
    # never hand an empty read set to the next steps
    [ "\$(zcat ${meta.sample}_ont_fungus.fastq.gz | head -c 1 | wc -c)" -gt 0 ] || { echo "ERROR: Badread produced no reads" >&2; exit 1; }
    """
}

process SIMULATE_ONT_SPIKE {
    tag "${meta.sample}:${pct}%"
    publishDir "${params.sim_dir}/reads", mode: 'copy'

    input:
    tuple val(meta), path(fungus), path(fungus_reads), val(pct), val(groups), path(contams), val(hash)

    output:
    tuple val(meta), val(pct), path("${meta.sample}_spike${pct}_ont.fastq.gz"), emit: reads
    path ("${meta.sample}_spike${pct}_long.composition.tsv"), emit: composition
    path "${meta.sample}_spike${pct}_long.sim_hash"

    script:
    def ds = "${meta.sample}_spike${pct}_long"
    """
    cp ${fungus_reads} ${meta.sample}_spike${pct}_ont.fastq.gz
    make_long_quantities.py --fungus ${fungus} --coverage ${params.long_coverage} --percent ${pct} \\
        --contam ${contams} --groups ${groups.join(' ')} --weights ${params.group_weights} > quantities.tsv
    while read fa bases; do
        [ "\$bases" -gt 0 ] || continue
        badread simulate --reference \$fa --quantity \$bases \\
            ${params.badread_args} --seed 43 \\
            | gzip >> contam_ont.fastq.gz
    done < quantities.tsv
    [ -s contam_ont.fastq.gz ] && cat contam_ont.fastq.gz >> ${meta.sample}_spike${pct}_ont.fastq.gz

    # realised ground truth (bases per source)
    fb=\$(zcat ${fungus_reads} | awk 'NR%4==2 {n+=length(\$0)} END {print n+0}')
    tb=\$(zcat ${meta.sample}_spike${pct}_ont.fastq.gz | awk 'NR%4==2 {n+=length(\$0)} END {print n+0}')
    [ "\$fb" -gt 0 ] || { echo "ERROR: no fungal long reads" >&2; exit 1; }
    [ "${pct}" = "0" ] || [ "\$tb" -gt "\$fb" ] || { echo "ERROR: no contaminant long reads at ${pct}%" >&2; exit 1; }
    awk -v ds=${ds} -v f=\$fb -v t=\$tb 'BEGIN { printf "dataset\\tfungal_bases\\tcontam_bases\\tcontam_pct\\n%s\\t%d\\t%d\\t%.3f\\n", ds, f, t-f, (t > 0 ? (t-f)*100/t : 0) }' \\
        > ${ds}.composition.tsv
    echo ${hash} > ${ds}.sim_hash
    """
}

// ================================================================ Track 1: ablation

process CP00_RAW_BASELINE {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), val('cp00_raw_baseline'), path("${meta.id}.cp00_raw_baseline.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    timeit.sh ${meta.id} cp00_raw_baseline spades \\
        spades.py -1 ${reads[0]} -2 ${reads[1]} --isolate -o spades -t ${task.cpus} -m ${task.memory.toGiga()}
    cp spades/scaffolds.fasta ${meta.id}.cp00_raw_baseline.fasta
    """
}

// cp01 is split in two so the hybrid track can reuse the cleaned Illumina reads without SPAdes
process CLEAN_READS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("reads_filtered_filtered_{1,2}.fastq.gz"), emit: reads
    path ("*.time.tsv"), emit: times

    script:
    def cp = 'cp01_spades_asm'
    def seed = meta.mito_seed ? "-s ${meta.mito_seed}" : ''
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} ${cp} trim \\
        AAFTF trim --method ${params.trim_method} -c ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o reads_trimmed
    # NOVOPlasty fails when the reference has no mitochondrial contig; then filter runs without it
    timeit.sh ${meta.id} ${cp} mito \\
        AAFTF mito -1 reads_trimmed_1P.fastq.gz -2 reads_trimmed_2P.fastq.gz ${seed} -o mito_assembly.fasta \\
        || echo "AAFTF mito failed; continuing without mitochondrial screen" >&2
    screen=\$( [ -s mito_assembly.fasta ] && echo "-s mito_assembly.fasta" || true )
    timeit.sh ${meta.id} ${cp} filter \\
        AAFTF filter -1 reads_trimmed_1P.fastq.gz -2 reads_trimmed_2P.fastq.gz \\
        --aligner ${params.filter_aligner} -c ${task.cpus} -o reads_filtered \$screen
    """
}

process CP01_SPADES_ASM {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), val('cp01_spades_asm'), path("${meta.id}.cp01_spades_asm.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} cp01_spades_asm assemble \\
        AAFTF assemble -1 ${reads[0]} -2 ${reads[1]} \\
        -c ${task.cpus} -m ${task.memory.toGiga()} --method spades -o ${meta.id}.cp01_spades_asm.fasta
    """
}

process CP01_FLYE_ASM {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), val('cp01_flye_asm'), path("${meta.id}.cp01_flye_asm.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} cp01_flye_asm assemble \\
        AAFTF assemble --method flye -lr ${reads} --longread_type ${meta.longread_type ?: params.longread_type} \\
        --genome_size ${meta.genome_size} -c ${task.cpus} -m ${task.memory.toGiga()} -o ${meta.id}.cp01_flye_asm.fasta
    """
}

process CP02_VECSCREEN {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm)

    output:
    tuple val(meta), val('cp02_vecscreen'), path("${meta.id}.cp02_vecscreen.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} cp02_vecscreen vecscreen \\
        AAFTF vecscreen -i ${asm} -c ${task.cpus} -o ${meta.id}.cp02_vecscreen.fasta
    """
}

// NCBI FCS-adaptor screen, compared with vecscreen at cp02; this branch is not continued to cp03-cp06.
// The run_fcsadaptor.sh script and container image come from the AAFTF databases (AAFTF_DB).
process CP02_FCS_SCREEN {
    tag "${meta.id}"
    module 'singularity'
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm)

    output:
    tuple val(meta), val('cp02_fcs_screen'), path("${meta.id}.cp02_fcs_screen.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} cp02_fcs_screen fcs_screen \\
        AAFTF fcs_screen -i ${asm} -o ${meta.id}.cp02_fcs_screen.fasta
    """
}

process CP03_TAXPURGE {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm)

    output:
    tuple val(meta), val('cp03_taxpurge'), path("${meta.id}.cp03_taxpurge.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} cp03_taxpurge taxpurge \\
        AAFTF taxpurge -i ${asm} -c ${task.cpus} -p ${meta.phylum} -o ${meta.id}.cp03_taxpurge.fasta
    """
}

// FCS-GX for all cp02 assemblies of the run in ONE job: the ~470 GB database is copied once to node-local
// $TMPDIR (unless it is already on tmpfs, e.g. /dev/shm) and loaded once, then `AAFTF fcs_gx_purge` runs once
// per taxid (it takes one -t per run) on all assemblies of that taxid. fcs_gx_db is the database prefix,
// e.g. /srv/projects/db/ncbi-fcs/0.5.4/gxdb/all.
process CP03_FCSGX {
    tag "${jobs.size()} assemblies"
    publishDir params.outdir, mode: 'copy', pattern: '*.fasta',
        saveAs: { fn -> "01_checkpoints/${fn.replaceAll(/\.cp03_fcsgx\.fasta$/, '')}/${fn}" }

    input:
    val jobs                               // [[dataset id, taxid], ...]
    path asms, stageAs: 'in/*'             // <dataset id>.cp02_vecscreen.fasta, same order as jobs

    output:
    path "*.cp03_fcsgx.fasta", emit: fasta
    path "*.time.tsv", emit: times

    script:
    def db_name = file(params.fcs_gx_db).name
    def runs = jobs.groupBy { j -> j[1] }.collect { taxid, js ->
        def inputs = js.collect { j -> "in/${j[0]}.cp02_vecscreen.fasta" }.join(' ')
        "timeit.sh fcsgx_batch cp03_fcsgx fcs_gx_purge_taxid${taxid} \\\n" +
        "        AAFTF fcs_gx_purge -i ${inputs} -c ${task.cpus} -d \"\$db\" -t ${taxid} -o fcsgx"
    }.join('\n    ')
    """
    # AAFTF ${params.aaftf_commit}
    if [ "\$(stat -f -c %T "\$(dirname ${params.fcs_gx_db})")" = tmpfs ]; then
        db=${params.fcs_gx_db}
    else
        db_dir=\$(mktemp -d "\${TMPDIR:-\$PWD}/gxdb.XXXXXX")
        trap 'rm -rf "\$db_dir"' EXIT
        timeit.sh fcsgx_batch cp03_fcsgx gxdb_copy cp ${params.fcs_gx_db}.* "\$db_dir"/
        db="\$db_dir/${db_name}"
    fi
    ${runs}
    for f in fcsgx/*.cp02_vecscreen.fcsgx_cleaned.fasta; do
        mv "\$f" "\$(basename "\$f" .cp02_vecscreen.fcsgx_cleaned.fasta).cp03_fcsgx.fasta"
    done
    """
}

// coverage purge with the filtered Illumina reads (short and hybrid; as AAFTF pipeline_short / pipeline_hybrid)
process CP04_COVPURGE {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm), path(reads)

    output:
    tuple val(meta), val('cp04_covpurge'), path("${meta.id}.cp04_covpurge.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} cp04_covpurge covpurge \\
        AAFTF covpurge -i ${asm} -1 ${reads[0]} -2 ${reads[1]} -mc ${params.covpurge_mincovpct} -c ${task.cpus} \\
        -o ${meta.id}.cp04_covpurge.fasta
    """
}

process CP05_RMDUP {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm)

    output:
    tuple val(meta), val('cp05_rmdup'), path("${meta.id}.cp05_rmdup.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} cp05_rmdup rmdup \\
        AAFTF rmdup -i ${asm} -c ${task.cpus} -ml ${params.min_contig} -o ${meta.id}.cp05_rmdup.fasta
    """
}

// long: racon with the ONT reads; hybrid: racon, then polypolish with filtered Illumina reads
process CP06_POLISH {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm), path(longreads), path(short_reads)

    output:
    tuple val(meta), val('cp06_polished'), path("${meta.id}.cp06_polished.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    def c = "cp06_polished"
    def hybrid = meta.track == 'hybrid'
    """
    # AAFTF ${params.aaftf_commit}
    timeit.sh ${meta.id} ${c} polish_racon \\
        AAFTF polish -i ${asm} -lr ${longreads} --method racon -c ${task.cpus} -m ${task.memory.toGiga()} -o racon.fasta
    ${hybrid
        ? """timeit.sh ${meta.id} ${c} polish_polypolish \\
        AAFTF polish -i racon.fasta -1 ${short_reads[0]} -2 ${short_reads[1]} --method polypolish \\
        -c ${task.cpus} -m ${task.memory.toGiga()} -o polished.fasta"""
        : 'mv racon.fasta polished.fasta'}
    mv polished.fasta ${meta.id}.cp06_polished.fasta
    """
}

// ================================================================ evaluation

process QUAST {
    tag "${meta.id}:${cp}"
    publishDir "${params.outdir}/02_quast", mode: 'copy'

    input:
    tuple val(meta), val(cp), path(asm), path(ref)

    output:
    path "${meta.id}__${cp}"

    script:
    """
    quast.py -r ${ref} --fungus --threads ${task.cpus} --output-dir ${meta.id}__${cp} ${asm}
    """
}

process BUSCO {
    tag "${meta.id}:${cp}"
    publishDir "${params.outdir}/03_completeness/busco", mode: 'copy'

    input:
    tuple val(meta), val(cp), path(asm)

    output:
    path "${meta.id}__${cp}"

    script:
    // shared, read-only BUSCO download folder (<dir>/lineages/<lineage>); --offline: never download or write there
    def dl = params.busco_downloads ? "--download_path ${params.busco_downloads} --offline" : ''
    // --opt-out-run-stats: skip BUSCO 6's usage-statistics upload, which fails on the cluster (HTTP 502) after a
    // ~16 min wait per job
    """
    busco -i ${asm} -m genome --miniprot -l ${meta.lineage} -o ${meta.id}__${cp} --out_path . \\
        --cpu ${task.cpus} ${dl} --opt-out-run-stats
    # drop ~3.5 GB of intermediates per run (miniprot alignments, a log duplicating them, hmmsearch output,
    # the unpacked protein database); keep short summaries, full_table.tsv, missing list, busco_sequences
    rm -rf ${meta.id}__${cp}/tmp ${meta.id}__${cp}/run_*/miniprot_output ${meta.id}__${cp}/run_*/hmmer_output \\
        ${meta.id}__${cp}/logs/miniprot_align_*_out.log
    """
}

// alignment of a checkpoint (query) to the sample's reference (target) for interactive dot plots in
// D-GENIES "Plot alignment": upload <name>.paf, <name>.target.idx and <name>.query.idx (D-GENIES index
// format: sample name, then one "sequence<TAB>length" line per sequence; small, unlike the FASTA files)
process DOTPLOT_PAF {
    tag "${meta.id}:${cp}"
    publishDir "${params.outdir}/05_dotplots", mode: 'copy'

    input:
    tuple val(meta), val(cp), path(asm), path(ref)

    output:
    path "${meta.id}__${cp}.{paf,target.idx,query.idx}"

    script:
    def n = "${meta.id}__${cp}"
    """
    minimap2 -x asm5 -t ${task.cpus} ${ref} ${asm} > ${n}.paf
    fasta_idx() { awk -v name="\$2" 'BEGIN { print name } /^>/ { if (id) print id "\t" len; id = substr(\$1, 2); len = 0; next } { len += length(\$0) } END { if (id) print id "\t" len }' "\$1"; }
    fasta_idx ${ref} ${meta.sample} > ${n}.target.idx
    fasta_idx ${asm} ${n} > ${n}.query.idx
    """
}

process MAP_TO_COMBINED_REF {
    tag "${meta.id}:${cp}"
    publishDir "${params.outdir}/04_decontamination/paf", mode: 'copy'

    input:
    tuple val(meta), val(cp), path(asm), path(combined)

    output:
    tuple val(meta), val(cp), path(asm), path("${meta.id}__${cp}.paf")

    script:
    """
    minimap2 -x asm5 --secondary=no -t ${task.cpus} ${combined} ${asm} > ${meta.id}__${cp}.paf
    """
}

process DECON_METRICS {
    tag "${meta.id}:${cp}"
    publishDir "${params.outdir}/04_decontamination/metrics", mode: 'copy'

    input:
    tuple val(meta), val(cp), path(asm, stageAs: 'after/*'), path(paf, stageAs: 'after/*'), path(before_asm, stageAs: 'before/*'), path(before_paf, stageAs: 'before/*')

    output:
    path "${meta.id}__${cp}.precision_recall_f1.tsv"

    script:
    // confusion matrix only for the purge step (cp02 → cp03); composition for every checkpoint
    // confusion matrix for the screens (cp02 vecscreen / fcs_screen, vs cp01) and the purge steps (cp03 taxpurge /
    // fcs_gx_purge, cp04 covpurge; cumulative from cp02); the 'before' input is chosen in the workflow
    def before = cp ==~ /cp0[234]_.*/ ? "--before_fasta ${before_asm} --before_paf ${before_paf}" : ''
    """
    calculate_decon_metrics.py --dataset ${meta.id} --checkpoint ${cp} --fasta ${asm} --paf ${paf} ${before} \\
        -o ${meta.id}__${cp}.precision_recall_f1.tsv
    """
}

// versions of the tools that produced the results (for the methods section), written once per run
process TOOL_VERSIONS {
    publishDir params.outdir, mode: 'copy'

    output:
    path 'versions.tsv'

    script:
    """
    # first version-like token of a command's output, NA if the tool is missing
    ver() { local out; out=\$("\$@" 2>&1 || true); awk 'match(\$0, /[0-9]+\\.[0-9]+[0-9A-Za-z.+_-]*/) { print substr(\$0, RSTART, RLENGTH); found = 1; exit } END { if (!found) print "NA" }' <<< "\$out"; }
    # tools of the eval / busco pixi environments (QUAST, BUSCO): only the line naming the tool, since pixi
    # prints its own messages first; NA when pixi is not used
    ver_env() { local env=\$1 name=\$2; shift 2; ver grep -i -m1 "\$name" <<< "\$(pixi run --manifest-path ${params.pixi_manifest} -e \$env "\$@" 2>&1 || true)"; }
    {
        printf "tool\\tversion\\n"
        printf "AAFTF commit (pixi.lock)\\t%s\\n" "${params.aaftf_commit}"
        printf "AAFTF\\t%s\\n" "\$(ver AAFTF --version)"
        printf "Nextflow\\t%s\\n" "${nextflow.version}"
        printf "SPAdes\\t%s\\n" "\$(ver spades.py --version)"
        printf "Flye\\t%s\\n" "\$(ver flye --version)"
        printf "minimap2\\t%s\\n" "\$(ver minimap2 --version)"
        printf "Racon\\t%s\\n" "\$(ver racon --version)"
        printf "Polypolish\\t%s\\n" "\$(ver polypolish --version)"
        printf "BBTools (bbduk)\\t%s\\n" "\$(ver bbduk.sh --version)"
        printf "BLAST+\\t%s\\n" "\$(ver blastn -version)"
        printf "sourmash\\t%s\\n" "\$(ver sourmash --version)"
        printf "samtools\\t%s\\n" "\$(ver samtools --version)"
        printf "InSilicoSeq\\t%s\\n" "\$(ver iss --version)"
        printf "Badread\\t%s\\n" "\$(ver badread --version)"
        printf "QUAST\\t%s\\n" "\$(ver_env eval quast quast.py --version)"
        printf "BUSCO\\t%s\\n" "\$(ver_env busco busco busco --version)"
    } > versions.tsv
    """
}

process AGGREGATE_RESULTS {
    publishDir params.outdir, mode: 'copy'

    input:
    path aliases, stageAs: 'eval_aliases.tsv'
    path quast, stageAs: 'quast/*'
    path busco, stageAs: 'busco/*'
    path decon, stageAs: 'decon/*'
    path times, stageAs: 'times.tsv'
    path composition, stageAs: 'comp/*'

    output:
    path 'final_report.tsv'
    path 'precision_recall_f1.tsv'
    path 'time_summary.tsv'
    path 'read_composition.tsv'

    script:
    """
    aggregate_benchmark_results.py --quast quast/* --busco busco/* \\
        --decon decon/* --times times.tsv --aliases eval_aliases.tsv -o final_report.tsv
    # cp03 rows have extra columns (TP/FN/FP/TN, precision, recall, F1): merge by column name
    python -c "import glob, pandas as pd; pd.concat([pd.read_csv(f, sep='\\t') for f in sorted(glob.glob('decon/*'))]).to_csv('precision_recall_f1.tsv', sep='\\t', index=False)"
    cp times.tsv time_summary.tsv
    awk 'FNR>1 || NR==1' comp/*   > read_composition.tsv
    """
}

// ================================================================ workflow

// Channel notation used in the comments below:
//   meta  sample metadata map, one per sample (sampleRow): sample (file-safe key), name (NCBI organism
//         name), origin ('sim' | 'real'), genome_size, phylum, taxid, lineage, mito_seed, source_reference
//         (downloaded FASTA), accessions, read1/read2/longreads/longread_type (real datasets only). Per
//         dataset it also gets id (<sample>_spike<pct>_<track> or <sample>_real_<track>), spike_pct ('real'
//         for real datasets), track, tagged_ref (tagged fungal FASTA) and combined_ref (sample +
//         contaminants FASTA; null for real datasets).
//   cp    checkpoint name, e.g. 'cp03_taxpurge'
//   (a, b, ...) one channel item (tuple)
workflow {
    // placeholder for optional path inputs (e.g. no Illumina reads for long-only polishing)
    NO_FILE = file("${projectDir}/assets/NO_FILE")

    // tracks to run (--tracks): which read types are needed follows from them
    def tracks = params.tracks.toString().tokenize(',')*.trim()
    def unknown_tracks = tracks - ['short', 'long', 'hybrid']
    if (unknown_tracks) {
        error("Unknown track(s) in --tracks: ${unknown_tracks.join(', ')} (use short, long, hybrid)")
    }
    def want_short = 'short' in tracks
    def want_long = 'long' in tracks
    def want_hybrid = 'hybrid' in tracks
    def need_illumina = want_short || want_hybrid
    def need_ont = want_long || want_hybrid

    // ---- references: simulated samples (sim_samples.tsv) and real datasets (real_samples.tsv)
    // sim_rows / real_rows: (meta, downloaded reference FASTA), one per row of sim_samples.tsv / real_samples.tsv
    sim_rows = channel.fromPath(params.sim_samples).splitCsv(header: true, sep: '\t').filter { r -> activeRow(r) }.map { r -> sampleRow(r, 'sim') }
    def real_tsv = params.real_samples ? file(params.real_samples) : null
    real_rows = real_tsv && real_tsv.exists()
        ? channel.fromPath(real_tsv).splitCsv(header: true, sep: '\t').filter { r -> activeRow(r) }.map { r -> sampleRow(r, 'real') }
        : channel.empty()
    // tagged_refs: (meta, tagged reference FASTA, headers FUNGUS_*) for every sample; a process can only be
    // called once, so simulated and real samples go through TAG_FUNGAL_REF together and are split again:
    // sim_refs (simulated samples) and real_refs (real datasets)
    tagged_refs = TAG_FUNGAL_REF(sim_rows.mix(real_rows))
    sim_refs = tagged_refs.filter { m, _f -> m.origin == 'sim' }
    real_refs = tagged_refs.filter { m, _f -> m.origin == 'real' }

    // contam_refs: (taxon key, group, tagged contaminant FASTA, headers CONTAM_<taxon>_*), one per contaminants.tsv row
    contam_refs = TAG_CONTAM_REF(
        channel.fromPath(params.contaminants).splitCsv(header: true, sep: '\t').filter { r -> activeRow(r) }.map { r -> tuple(safeKey(r.organism), r.group, r.accessions ?: '', file(r.reference)) }
    )
    // contam_panel: a single item (['taxon:group', ...], [contaminant FASTAs]) with the whole panel,
    // combined with every simulation job
    contam_panel = contam_refs.toList().map { l -> tuple(l.collect { t -> "${t[0]}:${t[1]}" }, l.collect { t -> t[2] }) }

    // combined: (meta, <sample>_with_contaminants.fa) per simulated sample; minimap2 target for contig labelling
    combined_refs = BUILD_COMBINED_REF(sim_refs.combine(contam_panel).map { m, f, _g, c -> tuple(m, f, c) })

    // levels: one item per spike level ('0', '1', '5', '15')
    spike_levels = channel.fromList(params.spike_levels.toString().tokenize(','))
    // sim_refs_combined: (meta, tagged fungal FASTA, sample_with_contaminants.fa) per simulated sample
    sim_refs_combined = sim_refs.join(combined_refs)
    // Channels filled step by step (mix) as the tracks below run:
    //   timing_tsvs:      timing TSVs from timeit.sh, one per timed command
    //   composition_tsvs: <dataset>.composition.tsv, realised fungal / contaminant amounts per simulated dataset
    //   cp01_assemblies:  (meta, cp, FASTA) of the first assembly of each dataset (cp01), input to the cleanup
    //   checkpoints:      (meta, cp, FASTA) of every checkpoint, input to the evaluation
    timing_tsvs = channel.empty()
    composition_tsvs = channel.empty()
    cp01_assemblies = channel.empty()
    checkpoints = channel.empty()
    // clean_illumina: (meta of short dataset, [filtered R1, R2]) from CLEAN_READS; used by cp01 SPAdes, covpurge and polypolish
    clean_illumina = channel.empty()
    // ont_reads: (meta of long dataset, ONT reads); used by flye (cp01) and racon (cp06)
    ont_reads = channel.empty()

    // ---- Illumina datasets: needed by the short track and for hybrid polishing
    if (need_illumina) {
        // illumina_jobs: (meta, fungal FASTA, pct, groups, contaminant FASTAs) per sample × spike level, branched
        // into 'reuse' (reads with a matching .sim_hash already in sim_dir/reads) and 'simulate'
        illumina_jobs = sim_refs
            .combine(spike_levels)
            .combine(contam_panel)
            .branch { m, _f, pct, _g, _c ->
                reuse: reusableReads(m, pct, 'short')
                simulate: true
            }
        // sim_illumina.reads: (meta, pct, [R1, R2]) newly simulated; sim_illumina.composition: composition TSVs
        sim_illumina = SIMULATE_ILLUMINA(illumina_jobs.simulate.map { m, f, pct, g, c -> tuple(m, f, pct, g, c, simHash(m, pct, 'short')) })
        // reused_illumina: (meta, pct, [R1, R2], composition TSV) from an earlier run
        reused_illumina = illumina_jobs.reuse.map { m, _f, pct, _g, _c ->
            def fs = reusableReads(m, pct, 'short')
            tuple(m, pct, fs[0..1], fs[2])
        }
        composition_tsvs = composition_tsvs.mix(sim_illumina.composition, reused_illumina.map { _m, _p, _r, c -> c })
        // sim_short_ds: (dataset meta, [R1, R2]) per simulated short-read dataset; dataset fields added to meta
        sim_short_ds = sim_illumina.reads
            .mix(reused_illumina.map { m, pct, r, _c -> tuple(m, pct, r) })
            .combine(sim_refs_combined, by: 0)
            .map { m, pct, reads, ref, comb ->
                tuple(m + [id: "${m.sample}_spike${pct}_short".toString(), spike_pct: pct, track: 'short', tagged_ref: ref, combined_ref: comb], reads)
            }
        // real: no contaminant answer key, so no combined FASTA (decontamination metrics are skipped)
        real_short_ds = real_refs
            .filter { m, _f -> m.read1 }
            .map { m, ref ->
                tuple(m + [id: "${m.sample}_real_short".toString(), spike_pct: 'real', track: 'short', tagged_ref: ref, combined_ref: null], [m.read1, m.read2].findAll { r -> r })
            }
        // short_ds: (dataset meta, [R1, R2]) for all short-read datasets, simulated and real
        short_ds = sim_short_ds.mix(real_short_ds)
        // cleaned.reads: (dataset meta, [filtered R1, R2]) after AAFTF trim → mito → filter
        cleaned = CLEAN_READS(short_ds)
        clean_illumina = cleaned.reads
        timing_tsvs = timing_tsvs.mix(cleaned.times)

        // ---- short track: cp00 baseline + cp01 SPAdes on cleaned reads
        if (want_short) {
            // cp00.fasta / cp01_short.fasta: (dataset meta, cp, FASTA)
            cp00 = CP00_RAW_BASELINE(short_ds)
            cp01_short = CP01_SPADES_ASM(clean_illumina)
            cp01_assemblies = cp01_assemblies.mix(cp01_short.fasta)
            checkpoints = checkpoints.mix(cp00.fasta, cp01_short.fasta)
            timing_tsvs = timing_tsvs.mix(cp00.times, cp01_short.times)
        }
    }

    // ---- ONT datasets: one flye assembly per dataset, shared by the long and hybrid
    // tracks until polishing (those shared checkpoints are labelled '..._long')
    if (need_ont) {
        // ont_jobs: as illumina_jobs, for ONT reads
        ont_jobs = sim_refs
            .combine(spike_levels)
            .combine(contam_panel)
            .branch { m, _f, pct, _g, _c ->
                reuse: reusableReads(m, pct, 'long')
                simulate: true
            }
        // fungal_ont_reads: (meta, fungal ONT reads), one per sample (shared by all its spike levels); badread runs
        // only for samples with at least one spike level still to simulate
        fungal_ont_reads = SIMULATE_ONT_FUNGUS(ont_jobs.simulate.map { m, f, _p, _g, _c -> tuple(m, f) }.unique())
        // sim_ont.reads: (meta, pct, ONT reads) with contaminant reads added; sim_ont.composition: composition TSVs
        sim_ont = SIMULATE_ONT_SPIKE(
            ont_jobs.simulate.combine(fungal_ont_reads, by: 0).map { m, f, pct, g, c, fr -> tuple(m, f, fr, pct, g, c, simHash(m, pct, 'long')) }
        )
        // reused_ont: (meta, pct, ONT reads, composition TSV) from an earlier run
        reused_ont = ont_jobs.reuse.map { m, _f, pct, _g, _c ->
            def fs = reusableReads(m, pct, 'long')
            tuple(m, pct, fs[0], fs[1])
        }
        composition_tsvs = composition_tsvs.mix(sim_ont.composition, reused_ont.map { _m, _p, _r, c -> c })
        // sim_long_ds: (dataset meta, ONT reads) per simulated long-read dataset
        sim_long_ds = sim_ont.reads
            .mix(reused_ont.map { m, pct, r, _c -> tuple(m, pct, r) })
            .combine(sim_refs_combined, by: 0)
            .map { m, pct, reads, ref, comb ->
                tuple(m + [id: "${m.sample}_spike${pct}_long".toString(), spike_pct: pct, track: 'long', tagged_ref: ref, combined_ref: comb], reads)
            }
        real_long_ds = real_refs
            .filter { m, _f -> m.longreads }
            .map { m, ref ->
                tuple(m + [id: "${m.sample}_real_long".toString(), spike_pct: 'real', track: 'long', tagged_ref: ref, combined_ref: null], m.longreads)
            }
        // long_ds: (dataset meta, ONT reads) for all long-read datasets, simulated and real
        long_ds = sim_long_ds.mix(real_long_ds)
        // cp01_long.fasta: (dataset meta, 'cp01_flye_asm', FASTA)
        cp01_long = CP01_FLYE_ASM(long_ds)
        cp01_assemblies = cp01_assemblies.mix(cp01_long.fasta)
        checkpoints = checkpoints.mix(cp01_long.fasta)
        timing_tsvs = timing_tsvs.mix(cp01_long.times)
        ont_reads = long_ds
    }

    // ---- cleanup: cp02 vecscreen → cp03 purge (shared by long and hybrid) → cp04 covpurge (short, hybrid)
    // → cp05 rmdup → cp06 polish (long, hybrid). The long track has no Illumina reads, so no cp04.
    // cp02.fasta, cp03_taxpurge.fasta: (dataset meta, cp, FASTA)
    cp02 = CP02_VECSCREEN(cp01_assemblies)
    // cp02_fcs: (dataset meta, 'cp02_fcs_screen', FASTA); FCS-adaptor on the same cp01 assemblies, compared
    // with vecscreen at cp02 only (cp03-cp06 continue from vecscreen)
    // toString(): --fcs_screen false on the command line arrives as the text 'false'
    if (params.fcs_screen.toString() == 'true') {
        cp02_fcs = CP02_FCS_SCREEN(cp01_assemblies)
        checkpoints = checkpoints.mix(cp02_fcs.fasta)
        timing_tsvs = timing_tsvs.mix(cp02_fcs.times)
    }
    cp03_taxpurge = CP03_TAXPURGE(cp02.fasta)
    // cp03: the taxpurge result, which the later steps (cp04-cp06) continue from
    cp03 = cp03_taxpurge.fasta
    timing_tsvs = timing_tsvs.mix(cp02.times, cp03_taxpurge.times)
    checkpoints = checkpoints.mix(cp02.fasta, cp03)

    // FCS-GX branch: only cp03, evaluated and compared with taxpurge; not continued to cp04-cp06.
    // One batched job for all cp02 assemblies (one database copy), so it starts once every cp02 exists.
    if (params.fcs_gx_db) {
        // fcsgx_batch: one item, the list of all (meta, cp, cp02 FASTA)
        fcsgx_batch = cp02.fasta.toList().filter { l -> l }
        cp03_fcsgx_job = CP03_FCSGX(
            fcsgx_batch.map { l -> l.collect { m, _cp, _a -> [m.id, m.taxid] } },
            fcsgx_batch.map { l -> l.collect { _m, _cp, a -> a } },
        )
        // cp03_fcsgx: (dataset meta, 'cp03_fcsgx', FASTA), matched back to each dataset by file name
        cp03_fcsgx = cp03_fcsgx_job.fasta
            .flatten()
            .map { f -> tuple(f.name.replaceAll(/\.cp03_fcsgx\.fasta$/, ''), f) }
            .join(cp02.fasta.map { m, _cp, _a -> tuple(m.id, m) })
            .map { _id, f, m -> tuple(m, 'cp03_fcsgx', f) }
        checkpoints = checkpoints.mix(cp03_fcsgx)
        timing_tsvs = timing_tsvs.mix(cp03_fcsgx_job.times)
    }

    // clean_illumina_by_key: ([sample, spike_pct], [filtered R1, R2]); trackLinkKey [sample, spike_pct] links datasets of the same sample and
    // spike level across tracks (short ↔ long ↔ hybrid)
    clean_illumina_by_key = clean_illumina.map { m, r -> tuple(trackLinkKey(m), r) }
    // covpurge_in_short: (short meta, cp, cp03 FASTA, [filtered R1, R2])
    covpurge_in_short = cp03.filter { m, _cp, _a -> m.track == 'short' }.combine(clean_illumina, by: 0)
    // covpurge_in_hybrid: (hybrid meta, cp, long cp03 FASTA, [filtered R1, R2]); hybrid branches off the long assembly
    // after cp03, where Illumina reads exist for the same sample and spike level
    covpurge_in_hybrid = want_hybrid
        ? cp03.filter { m, _cp, _a -> m.track == 'long' }
            .map { m, cp, a -> tuple(trackLinkKey(m), m, cp, a) }
            .combine(clean_illumina_by_key, by: 0)
            .map { _k, m, cp, a, r -> tuple(toHybridMeta(m), cp, a, r) }
        : channel.empty()
    // cp04.fasta: (short / hybrid meta, 'cp04_covpurge', FASTA)
    cp04 = CP04_COVPURGE(covpurge_in_short.mix(covpurge_in_hybrid))
    timing_tsvs = timing_tsvs.mix(cp04.times)
    checkpoints = checkpoints.mix(cp04.fasta)

    // rmdup_in_long: (long meta, cp, FASTA); the long track goes from cp03 straight to rmdup
    rmdup_in_long = want_long ? cp03.filter { m, _cp, _a -> m.track == 'long' } : channel.empty()
    // cp05.fasta: (meta, 'cp05_rmdup', FASTA) for short, long and hybrid datasets
    cp05 = CP05_RMDUP(cp04.fasta.mix(rmdup_in_long))
    timing_tsvs = timing_tsvs.mix(cp05.times)
    checkpoints = checkpoints.mix(cp05.fasta)

    // ---- cp06: long (racon) and hybrid (racon + polypolish)
    if (need_ont) {
        // ont_reads_by_key: ([sample, spike_pct], ONT reads)
        ont_reads_by_key = ont_reads.map { m, r -> tuple(trackLinkKey(m), r) }
        // polish_in_long: (long meta, cp, cp05 FASTA, ONT reads, NO_FILE): racon only
        polish_in_long = cp05.fasta
            .filter { m, _cp, _a -> m.track == 'long' }
            .map { m, cp, a -> tuple(trackLinkKey(m), m, cp, a) }
            .combine(ont_reads_by_key, by: 0)
            .map { _k, m, cp, a, lr -> tuple(m, cp, a, lr, NO_FILE) }
        // polish_in_hybrid: (hybrid meta, cp, cp05 FASTA, ONT reads, [filtered R1, R2]): racon, then polypolish
        polish_in_hybrid = cp05.fasta
            .filter { m, _cp, _a -> m.track == 'hybrid' }
            .map { m, cp, a -> tuple(trackLinkKey(m), m, cp, a) }
            .combine(ont_reads_by_key, by: 0)
            .combine(clean_illumina_by_key, by: 0)
            .map { _k, m, cp, a, lr, sr -> tuple(m, cp, a, lr, sr) }
        // cp06.fasta: (long / hybrid meta, 'cp06_polished', FASTA)
        cp06 = CP06_POLISH(polish_in_long.mix(polish_in_hybrid))
        timing_tsvs = timing_tsvs.mix(cp06.times)
        checkpoints = checkpoints.mix(cp06.fasta)
    }

    // ---- evaluation
    // QUAST / BUSCO run once per unique assembly (same sample + same sequences, ignoring
    // contig names and order), e.g. when a step removes nothing; the first checkpoint seen is
    // evaluated at once and the others get its results in AGGREGATE_RESULTS (eval_aliases.tsv)
    // checkpoints_fp: ('<sample>:<sequence fingerprint>', meta, cp, FASTA) for every checkpoint
    checkpoints_fp = checkpoints.map { m, cp, a -> tuple("${m.sample}:${seqFingerprint(a)}".toString(), m, cp, a) }
    // unique_checkpoints: (meta, cp, FASTA), the first checkpoint seen for each fingerprint
    unique_checkpoints = checkpoints_fp.unique { t -> t[0] }.map { _k, m, cp, a -> tuple(m, cp, a) }
    // eval_aliases: one TSV (key, <dataset>__<cp>) listing every checkpoint with its fingerprint
    eval_aliases = checkpoints_fp
        .map { k, m, cp, _a -> "${k}\t${m.id}__${cp}\n".toString() }
        .collectFile(name: 'eval_aliases.tsv', seed: "key\tname\n", sort: true)
    // quast_dirs / busco_dirs: one result folder <dataset>__<cp> per unique assembly
    quast_dirs = QUAST(unique_checkpoints.map { m, cp, a -> tuple(m, cp, a, m.tagged_ref) })
    // PAF of every checkpoint against the sample's reference, for D-GENIES (results/05_dotplots/<dataset>__<cp>.paf)
    DOTPLOT_PAF(checkpoints.map { m, cp, a -> tuple(m, cp, a, m.tagged_ref) })
    busco_dirs = BUSCO(unique_checkpoints)
    // combined_paf: (meta, cp, FASTA, PAF vs the combined FASTA) per checkpoint of a simulated dataset (real datasets have no combined FASTA,
    // so no decontamination metrics)
    combined_paf = MAP_TO_COMBINED_REF(checkpoints.filter { m, _cp, _a -> m.combined_ref }.map { m, cp, a -> tuple(m, cp, a, m.combined_ref) })
    // cp02_combined_paf: (dataset id, cp02 FASTA, cp02 PAF), the 'before' state for every checkpoint of that dataset;
    // hybrid datasets share cp01-cp03 with their long dataset, so they are compared with its cp02
    cp02_combined_paf = combined_paf.filter { _m, cp, _a, _p -> cp == 'cp02_vecscreen' }.map { m, _cp, a, p -> tuple(m.id, a, p) }
    // cp01_combined_paf: (dataset id, cp01 FASTA, cp01 PAF), the 'before' state of the two cp02 screens
    // (vecscreen, fcs_screen), so both are scored on the same input
    cp01_combined_paf = combined_paf.filter { _m, cp, _a, _p -> cp ==~ /cp01_.*/ }.map { m, _cp, a, p -> tuple(m.id, a, p) }
    // decon_tsvs: one <dataset>__<cp>.precision_recall_f1.tsv per checkpoint
    decon_tsvs = DECON_METRICS(
        combined_paf.filter { _m, cp, _a, _p -> cp ==~ /cp02_.*/ }
            .map { m, cp, a, p -> tuple(m.id, m, cp, a, p) }
            .combine(cp01_combined_paf, by: 0)
            .mix(
                combined_paf.filter { _m, cp, _a, _p -> !(cp ==~ /cp02_.*/) }
                    .map { m, cp, a, p -> tuple(m.id.replaceAll('_hybrid$', '_long'), m, cp, a, p) }
                    .combine(cp02_combined_paf, by: 0)
            )
            .map { _id, m, cp, a, p, before_a, before_p -> tuple(m, cp, a, p, before_a, before_p) }
    )

    // AGGREGATE_RESULTS waits for everything above and writes final_report.tsv and the summary tables
    TOOL_VERSIONS()

    AGGREGATE_RESULTS(
        eval_aliases,
        quast_dirs.collect(),
        busco_dirs.collect(),
        decon_tsvs.collect().ifEmpty([]),
        // one table: per-step timing files share names (e.g. vecscreen.time.tsv) across datasets
        timing_tsvs.collectFile(name: 'times.tsv', keepHeader: true, skip: 1, sort: true),
        composition_tsvs.collect().ifEmpty([]),
    )
}
