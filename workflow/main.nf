#!/usr/bin/env nextflow
/*
 * AAFTF internal benchmarking suite (see aaftf_benchmark_specification.md)
 *
 *   1. Normalise reference headers (FUNGUS_*, CONTAM_<taxon>_*) and build one combined FASTA per sample
 *      (sample + all contaminants, data/sim/references/combined/<sample>_with_contaminants.fa)
 *   2. Simulate reads at 0/1/5/15% contaminant spike-in: InSilicoSeq Illumina and/or Badread ONT (as --tracks needs)
 *   3. Six-checkpoint ablation per dataset (track 'short'):
 *        cp00_raw_baseline        spades.py --isolate on raw reads
 *        cp01_mito_filtered_asm   AAFTF trim → mito → filter → assemble
 *        cp02_vecscreen           AAFTF vecscreen
 *        cp03_purged_{taxpurge,fcsgx}   branch A: AAFTF taxpurge | branch B: AAFTF fcs_gx_purge
 *        cp04_rmdup_{branch}      AAFTF rmdup
 *        cp05_sorted_{branch}     AAFTF sort (no polishing, as in AAFTF pipeline_short)
 *      --tracks (short,long,hybrid) picks which tracks run per sample × spike level (ONT assembled with AAFTF assemble --method flye):
 *        long:   cp01 flye → cp02 vecscreen → cp03 purge → cp04 rmdup → cp05 polish (racon) → sort
 *        hybrid: same assembly and cleanup as 'long'; cp05 = racon, then polypolish with the
 *                filtered Illumina reads of the same dataset (cp01 of the short track)
 *   4. Evaluate every checkpoint: QUAST, compleasm, BUSCO (miniprot), minimap2 vs the combined FASTA
 *   5. Aggregate into results/final_report.tsv
 * Every AAFTF / assembler call is wrapped in GNU time (scripts/timeit.sh).
 */
nextflow.enable.dsl = 2

params.samples = 'config/samples.tsv'
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
params.badread_args = '--length 15000,8000 --error_model nanopore2023 --qscore_model nanopore2023 --identity 98.5,99.8,99.2 --chimeras 1.0 --glitches 10000,25,25'
params.long_coverage = 35
params.longread_type = 'nano-hq'
params.trim_method = 'bbduk'
params.filter_aligner = 'bbduk'
params.min_contig = 500
params.fcs_gx_db = null
params.busco_downloads = null

// (sample, spike level) links the short, long and hybrid datasets built from the same simulation
// file/header-safe key for an NCBI organism name (contaminant taxa and sample names), e.g. 'Cutibacterium acnes subsp. acnes' -> 'Cutibacterium_acnes_subsp_acnes'
def taxonKey(name) {
    name.trim().replaceAll(/[^A-Za-z0-9]+/, '_').replaceAll(/^_|_$/, '')
}

// one meta map + reference per row of samples.tsv (origin 'sim') or real_samples.tsv (origin 'real')
def optFile(x) {
    x ? file(x) : null
}

def sampleRow(r, origin) {
    def meta = [
        sample: taxonKey(r.organism), name: r.organism, origin: origin, genome_size: r.genome_size, phylum: r.phylum,
        taxid: r.taxid, lineage: r.lineage, mito_seed: optFile(r.mito_seed), reference: r.reference,
        accessions: r.accessions ?: '', read1: optFile(r.read1), read2: optFile(r.read2), longreads: optFile(r.longreads), longread_type: r.longread_type ?: null,
    ]
    tuple(meta, file(r.reference))
}

def datasetKey(m) {
    [m.sample, m.pct]
}
def asHybrid(m) {
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
    def lines = file(params.contaminants).readLines().findAll { l -> l.trim() }
    def ref_col = lines[0].split('\t', -1).toList().indexOf('reference')
    def contam_refs = lines.drop(1).collect { l -> l.split('\t', -1)[ref_col] }
    def parts = [kind, pct, fileSig(m.reference), fileSig(params.contaminants)] + contam_refs.collect { f -> fileSig(f) } + [params.group_weights, file("${projectDir}/scripts/format_reference_headers.py").text]
    parts += kind == 'short'
        ? [params.short_coverage, params.read_length, params.iss_model, file("${projectDir}/scripts/make_readcounts.py").text]
        : [params.long_coverage, params.badread_args, file("${projectDir}/scripts/make_long_quantities.py").text]
    parts.join('|').md5()
}

// previously simulated files for a sample × spike level ('short' or 'long'); [] unless all exist
// and were made with the same settings (matching .sim_hash)
def existingReads(m, pct, kind) {
    def d = "${params.sim_dir}/reads/${m.sample}_spike${pct}"
    def fs = kind == 'short'
        ? [file("${d}_R1.fastq.gz"), file("${d}_R2.fastq.gz"), file("${d}_short.composition.tsv")]
        : [file("${d}_ont.fastq.gz"), file("${d}_long.composition.tsv")]
    def h = file("${d}_${kind}.sim_hash")
    def ok = params.reuse_reads.toString() == 'true' && fs.every { f -> f.exists() && f.size() > 0 } && h.exists() && h.text.trim() == simHash(m, pct, kind)
    ok ? fs : []
}

// ================================================================ data generation

process FORMAT_FUNGUS {
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

process FORMAT_CONTAM {
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
process COMBINED_REFERENCE {
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

process SIMULATE_SHORT {
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
    awk -v ds=${ds_track} -v nf=\$nf -v nt=\$nt 'BEGIN { printf "dataset\\tfungal_pairs\\tcontam_pairs\\tcontam_pct\\n%s\\t%d\\t%d\\t%.3f\\n", ds, nf, nt-nf, (nt-nf)*100/nt }' \\
        > ${ds_track}.composition.tsv
    echo ${hash} > ${ds_track}.sim_hash
    """
}

process SIMULATE_LONG_FUNGUS {
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
    """
}

process SIMULATE_LONG_SPIKE {
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
    awk -v ds=${ds} -v f=\$fb -v t=\$tb 'BEGIN { printf "dataset\\tfungal_bases\\tcontam_bases\\tcontam_pct\\n%s\\t%d\\t%d\\t%.3f\\n", ds, f, t-f, (t-f)*100/t }' \\
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
    tuple val(meta), val('cp00_raw_baseline'), path("${meta.id}.cp00.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    timeit.sh ${meta.id} cp00_raw_baseline spades \\
        spades.py -1 ${reads[0]} -2 ${reads[1]} --isolate -o spades -t ${task.cpus} -m ${task.memory.toGiga()}
    cp spades/scaffolds.fasta ${meta.id}.cp00.fasta
    """
}

// cp01 is split in two so the hybrid track can reuse the cleaned Illumina reads without SPAdes
process READ_PREP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("reads_filtered_filtered_{1,2}.fastq.gz"), emit: reads
    path ("*.time.tsv"), emit: times

    script:
    def cp = 'cp01_mito_filtered_asm'
    def seed = meta.mito_seed ? "-s ${meta.mito_seed}" : ''
    """
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

process CP01_MITO_FILTERED_ASM {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), val('cp01_mito_filtered_asm'), path("${meta.id}.cp01.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    timeit.sh ${meta.id} cp01_mito_filtered_asm assemble \\
        AAFTF assemble -1 ${reads[0]} -2 ${reads[1]} \\
        -c ${task.cpus} -m ${task.memory.toGiga()} --method spades -o ${meta.id}.cp01.fasta
    """
}

process CP01_LONG_ASM {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), val('cp01_flye_asm'), path("${meta.id}.cp01.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    timeit.sh ${meta.id} cp01_flye_asm assemble \\
        AAFTF assemble --method flye -lr ${reads} --longread_type ${meta.longread_type ?: params.longread_type} \\
        --genome_size ${meta.genome_size} -c ${task.cpus} -m ${task.memory.toGiga()} -o ${meta.id}.cp01.fasta
    """
}

process CP02_VECSCREEN {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm)

    output:
    tuple val(meta), val('cp02_vecscreen'), path("${meta.id}.cp02.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    timeit.sh ${meta.id} cp02_vecscreen vecscreen \\
        AAFTF vecscreen -i ${asm} -c ${task.cpus} -o ${meta.id}.cp02.fasta
    """
}

process CP03_TAXPURGE {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm)

    output:
    tuple val(meta), val('cp03_purged_taxpurge'), path("${meta.id}.cp03_taxpurge.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    timeit.sh ${meta.id} cp03_purged_taxpurge taxpurge \\
        AAFTF taxpurge -i ${asm} -c ${task.cpus} -p ${meta.phylum} -o ${meta.id}.cp03_taxpurge.fasta
    """
}

process CP03_FCSGX {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(_cp), path(asm)

    output:
    tuple val(meta), val('cp03_purged_fcsgx'), path("${meta.id}.cp03_fcsgx.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    """
    timeit.sh ${meta.id} cp03_purged_fcsgx fcs_gx_purge \\
        AAFTF fcs_gx_purge -i ${asm} -c ${task.cpus} -d ${params.fcs_gx_db} -t ${meta.taxid} -o fcsgx
    mv fcsgx/*.fcsgx_cleaned.fasta ${meta.id}.cp03_fcsgx.fasta
    """
}

process CP04_RMDUP {
    tag "${meta.id}:${branch}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(cp), path(asm)

    output:
    tuple val(meta), val("cp04_rmdup_${branch}"), path("${meta.id}.cp04_${branch}.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    branch = cp.tokenize('_')[-1]
    """
    timeit.sh ${meta.id} cp04_rmdup_${branch} rmdup \\
        AAFTF rmdup -i ${asm} -c ${task.cpus} -ml ${params.min_contig} -o ${meta.id}.cp04_${branch}.fasta
    """
}

// AAFTF pipeline_short does not polish (re-polishing with the assembly's own Illumina reads rarely helps),
// so the short-track cp05 is sort only
process CP05_SORTED {
    tag "${meta.id}:${branch}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(cp), path(asm)

    output:
    tuple val(meta), val("cp05_sorted_${branch}"), path("${meta.id}.cp05_${branch}.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    branch = cp.tokenize('_')[-1]
    """
    timeit.sh ${meta.id} cp05_sorted_${branch} sort \\
        AAFTF sort -i ${asm} -ml ${params.min_contig} -o ${meta.id}.cp05_${branch}.fasta
    """
}

// long: racon with the ONT reads; hybrid: racon, then polypolish with filtered Illumina reads
process CP05_LONG_POLISHED_SORTED {
    tag "${meta.id}:${branch}"
    publishDir path: { "${params.outdir}/01_checkpoints/${meta.id}" }, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(meta), val(cp), path(asm), path(longreads), path(short_reads)

    output:
    tuple val(meta), val("cp05_polished_${branch}"), path("${meta.id}.cp05_${branch}.fasta"), emit: fasta
    path ("*.time.tsv"), emit: times

    script:
    branch = cp.tokenize('_')[-1]
    def c = "cp05_polished_${branch}"
    def hybrid = meta.track == 'hybrid'
    """
    timeit.sh ${meta.id} ${c} polish_racon \\
        AAFTF polish -i ${asm} -lr ${longreads} --method racon -c ${task.cpus} -m ${task.memory.toGiga()} -o racon.fasta
    ${hybrid
        ? """timeit.sh ${meta.id} ${c} polish_polypolish \\
        AAFTF polish -i racon.fasta -1 ${short_reads[0]} -2 ${short_reads[1]} --method polypolish \\
        -c ${task.cpus} -m ${task.memory.toGiga()} -o polished.fasta"""
        : 'mv racon.fasta polished.fasta'}
    timeit.sh ${meta.id} ${c} sort \\
        AAFTF sort -i polished.fasta -ml ${params.min_contig} -o ${meta.id}.cp05_${branch}.fasta
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

process COMPLEASM {
    tag "${meta.id}:${cp}"
    publishDir "${params.outdir}/03_completeness/compleasm", mode: 'copy'

    input:
    tuple val(meta), val(cp), path(asm)

    output:
    path "${meta.id}__${cp}"

    script:
    def lib = params.busco_downloads ? "-L ${params.busco_downloads}" : ''
    """
    compleasm run -a ${asm} -o ${meta.id}__${cp} -l ${meta.lineage.replaceAll('_odb\\d+$', '')} -t ${task.cpus} ${lib}
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
    def dl = params.busco_downloads ? "--download_path ${params.busco_downloads}" : ''
    """
    busco -i ${asm} -m genome --miniprot -l ${meta.lineage} -o ${meta.id}__${cp} --out_path . \\
        --cpu ${task.cpus} ${dl}
    """
}

process MINIMAP2_COMBINED {
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
    def before = cp.startsWith('cp03') ? "--before_fasta ${before_asm} --before_paf ${before_paf}" : ''
    """
    calculate_decon_metrics.py --dataset ${meta.id} --checkpoint ${cp} --fasta ${asm} --paf ${paf} ${before} \\
        -o ${meta.id}__${cp}.precision_recall_f1.tsv
    """
}

process AGGREGATE {
    publishDir params.outdir, mode: 'copy'

    input:
    path quast, stageAs: 'quast/*'
    path compleasm, stageAs: 'compleasm/*'
    path busco, stageAs: 'busco/*'
    path decon, stageAs: 'decon/*'
    path times, stageAs: 'times/*'
    path composition, stageAs: 'comp/*'

    output:
    path 'final_report.tsv'
    path 'precision_recall_f1.tsv'
    path 'time_summary.tsv'
    path 'read_composition.tsv'

    script:
    """
    aggregate_benchmark_results.py --quast quast/* --compleasm compleasm/* --busco busco/* \\
        --decon decon/* --times times/* -o final_report.tsv
    awk 'FNR>1 || NR==1' decon/*  > precision_recall_f1.tsv
    awk 'FNR>1 || NR==1' times/*  > time_summary.tsv
    awk 'FNR>1 || NR==1' comp/*   > read_composition.tsv
    """
}

// ================================================================ workflow

workflow {
    NO_FILE = file("${projectDir}/assets/NO_FILE")

    def tracks = params.tracks.toString().tokenize(',')*.trim()
    def bad = tracks - ['short', 'long', 'hybrid']
    if (bad) {
        error("Unknown track(s) in --tracks: ${bad.join(', ')} (use short, long, hybrid)")
    }
    def want_short = 'short' in tracks
    def want_long = 'long' in tracks
    def want_hybrid = 'hybrid' in tracks
    def need_illumina = want_short || want_hybrid
    def need_ont = want_long || want_hybrid

    // ---- references: simulated samples (samples.tsv) and real datasets (real_samples.tsv)
    sim_rows = channel.fromPath(params.samples).splitCsv(header: true, sep: '\t').map { r -> sampleRow(r, 'sim') }
    def real_tsv = params.real_samples ? file(params.real_samples) : null
    real_rows = real_tsv && real_tsv.exists()
        ? channel.fromPath(real_tsv).splitCsv(header: true, sep: '\t').map { r -> sampleRow(r, 'real') }
        : channel.empty()
    fungus_all = FORMAT_FUNGUS(sim_rows.mix(real_rows))
    fungus = fungus_all.filter { m, _f -> m.origin == 'sim' }
    real_fungus = fungus_all.filter { m, _f -> m.origin == 'real' }

    contams = FORMAT_CONTAM(
        channel.fromPath(params.contaminants).splitCsv(header: true, sep: '\t').map { r -> tuple(taxonKey(r.organism), r.group, r.accessions ?: '', file(r.reference)) }
    )
    contam_set = contams.toList().map { l -> tuple(l.collect { t -> "${t[0]}:${t[1]}" }, l.collect { t -> t[2] }) }

    combined = COMBINED_REFERENCE(fungus.combine(contam_set).map { m, f, _g, c -> tuple(m, f, c) })

    levels = channel.fromList(params.spike_levels.toString().tokenize(','))
    // (meta, fungus.fa, sample_with_contaminants.fa)
    refs = fungus.join(combined)
    times = channel.empty()
    composition = channel.empty()
    assemblies = channel.empty()
    checkpoints = channel.empty()
    // (meta of short dataset, filtered Illumina reads)
    illumina = channel.empty()
    // (meta of long dataset, ONT reads)
    ont = channel.empty()

    // ---- Illumina datasets: needed by the short track and for hybrid polishing
    if (need_illumina) {
        short_jobs = fungus
            .combine(levels)
            .combine(contam_set)
            .branch { m, _f, pct, _g, _c ->
                reuse: existingReads(m, pct, 'short')
                simulate: true
            }
        sim = SIMULATE_SHORT(short_jobs.simulate.map { m, f, pct, g, c -> tuple(m, f, pct, g, c, simHash(m, pct, 'short')) })
        reused = short_jobs.reuse.map { m, _f, pct, _g, _c ->
            def fs = existingReads(m, pct, 'short')
            tuple(m, pct, fs[0..1], fs[2])
        }
        composition = composition.mix(sim.composition, reused.map { _m, _p, _r, c -> c })
        sim_short = sim.reads
            .mix(reused.map { m, pct, r, _c -> tuple(m, pct, r) })
            .combine(refs, by: 0)
            .map { m, pct, reads, ref, comb ->
                tuple(m + [id: "${m.sample}_spike${pct}_short".toString(), pct: pct, track: 'short', ref: ref, combined: comb], reads)
            }
        // real: no contaminant answer key, so no combined FASTA (decontamination metrics are skipped)
        real_short = real_fungus
            .filter { m, _f -> m.read1 }
            .map { m, ref ->
                tuple(m + [id: "${m.sample}_real_short".toString(), pct: 'real', track: 'short', ref: ref, combined: null], [m.read1, m.read2].findAll { r -> r })
            }
        datasets = sim_short.mix(real_short)
        prep = READ_PREP(datasets)
        illumina = prep.reads
        times = times.mix(prep.times)

        // ---- short track: cp00 baseline + cp01 SPAdes on cleaned reads
        if (want_short) {
            cp00 = CP00_RAW_BASELINE(datasets)
            cp01 = CP01_MITO_FILTERED_ASM(illumina)
            assemblies = assemblies.mix(cp01.fasta)
            checkpoints = checkpoints.mix(cp00.fasta, cp01.fasta)
            times = times.mix(cp00.times, cp01.times)
        }
    }

    // ---- ONT datasets: one flye assembly per dataset, shared by the long and hybrid
    // tracks until polishing (those shared checkpoints are labelled '..._long')
    if (need_ont) {
        long_jobs = fungus
            .combine(levels)
            .combine(contam_set)
            .branch { m, _f, pct, _g, _c ->
                reuse: existingReads(m, pct, 'long')
                simulate: true
            }
        // badread on the fungus only for samples with at least one spike level to simulate
        fungus_ont = SIMULATE_LONG_FUNGUS(long_jobs.simulate.map { m, f, _p, _g, _c -> tuple(m, f) }.unique())
        lsim = SIMULATE_LONG_SPIKE(
            long_jobs.simulate.combine(fungus_ont, by: 0).map { m, f, pct, g, c, fr -> tuple(m, f, fr, pct, g, c, simHash(m, pct, 'long')) }
        )
        lreused = long_jobs.reuse.map { m, _f, pct, _g, _c ->
            def fs = existingReads(m, pct, 'long')
            tuple(m, pct, fs[0], fs[1])
        }
        composition = composition.mix(lsim.composition, lreused.map { _m, _p, _r, c -> c })
        sim_long = lsim.reads
            .mix(lreused.map { m, pct, r, _c -> tuple(m, pct, r) })
            .combine(refs, by: 0)
            .map { m, pct, reads, ref, comb ->
                tuple(m + [id: "${m.sample}_spike${pct}_long".toString(), pct: pct, track: 'long', ref: ref, combined: comb], reads)
            }
        real_long = real_fungus
            .filter { m, _f -> m.longreads }
            .map { m, ref ->
                tuple(m + [id: "${m.sample}_real_long".toString(), pct: 'real', track: 'long', ref: ref, combined: null], m.longreads)
            }
        long_ds = sim_long.mix(real_long)
        cp01l = CP01_LONG_ASM(long_ds)
        assemblies = assemblies.mix(cp01l.fasta)
        checkpoints = checkpoints.mix(cp01l.fasta)
        times = times.mix(cp01l.times)
        ont = long_ds
    }

    // ---- shared cleanup: cp02 vecscreen → cp03 purge → cp04 rmdup
    cp02 = CP02_VECSCREEN(assemblies)
    cp03a = CP03_TAXPURGE(cp02.fasta)
    cp03 = cp03a.fasta
    times = times.mix(cp02.times, cp03a.times)
    if (params.fcs_gx_db) {
        cp03b = CP03_FCSGX(cp02.fasta)
        cp03 = cp03.mix(cp03b.fasta)
        times = times.mix(cp03b.times)
    }
    cp04 = CP04_RMDUP(cp03)
    times = times.mix(cp04.times)
    checkpoints = checkpoints.mix(cp02.fasta, cp03, cp04.fasta)

    // ---- cp05 short: sort only
    if (want_short) {
        cp05 = CP05_SORTED(cp04.fasta.filter { m, _cp, _a -> m.track == 'short' })
        times = times.mix(cp05.times)
        checkpoints = checkpoints.mix(cp05.fasta)
    }

    // ---- cp05 long (racon) and/or hybrid (racon + polypolish)
    if (need_ont) {
        // (sample, pct) links a long dataset to the Illumina reads of the same sample and spike level
        long_cp04 = cp04.fasta
            .filter { m, _cp, _a -> m.track == 'long' }
            .map { m, cp, a -> tuple(datasetKey(m), m, cp, a) }
            .combine(ont.map { m, r -> tuple(datasetKey(m), r) }, by: 0)
        long_in = want_long
            ? long_cp04.map { _k, m, cp, a, lr -> tuple(m, cp, a, lr, NO_FILE) }
            : channel.empty()
        // hybrid only where Illumina reads exist for the same dataset
        hybrid_in = want_hybrid
            ? long_cp04.combine(illumina.map { m, r -> tuple(datasetKey(m), r) }, by: 0).map { _k, m, cp, a, lr, sr -> tuple(asHybrid(m), cp, a, lr, sr) }
            : channel.empty()
        cp05l = CP05_LONG_POLISHED_SORTED(long_in.mix(hybrid_in))
        times = times.mix(cp05l.times)
        checkpoints = checkpoints.mix(cp05l.fasta)
    }

    // ---- evaluation
    quast = QUAST(checkpoints.map { m, cp, a -> tuple(m, cp, a, m.ref) })
    comp = COMPLEASM(checkpoints)
    busco = BUSCO(checkpoints)
    // decontamination metrics need the combined FASTA, i.e. simulated datasets only
    paf = MINIMAP2_COMBINED(checkpoints.filter { m, _cp, _a -> m.combined }.map { m, cp, a -> tuple(m, cp, a, m.combined) })
    cp02_paf = paf.filter { _m, cp, _a, _p -> cp == 'cp02_vecscreen' }.map { m, _cp, a, p -> tuple(m, a, p) }
    decon = DECON_METRICS(paf.combine(cp02_paf, by: 0))

    AGGREGATE(
        quast.collect(),
        comp.collect(),
        busco.collect(),
        decon.collect().ifEmpty([]),
        times.collect(),
        composition.collect().ifEmpty([]),
    )
}
