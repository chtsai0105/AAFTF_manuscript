# Specification: AAFTF Assembly Pipeline Internal Benchmarking Suite

## 1. Project Overview & Scope
This workflow implements an automated, reproducible benchmarking harness for the fungal genome assembly pipeline **AAFTF** (Automatic Assembly For The Fungi).

### Benchmark Design Principles
1. **Self-Contained Evaluation:** Do NOT compare against external wrapper pipelines (e.g., Shovill, fspassemblypipeline). All benchmarks compare AAFTF intermediate stages against a naive raw assembler baseline (short reads) or evaluate interchangeable internal modules.
2. **Deterministic Ground Truth:** All read datasets are simulated from finished, chromosome-level reference genomes using `InSilicoSeq` (short reads) and `Badread` (long reads), with fixed seeds and exact, recorded contaminant fractions.
3. **Decoupled Verification:** Evaluation relies strictly on nucleotide-level validation (`QUAST`, `compleasm`, `BUSCO -m genome --miniprot`, and `minimap2` alignment against known reference targets). Do NOT invoke *ab initio* gene annotation pipelines (e.g., Augustus, BRAKER, funannotate).

### Sequencing Tracks
Each sample × spike level is run on up to three tracks, chosen with `tracks` (`short`, `long`, `hybrid`; default `short`):
* **short:** Illumina reads only (mirrors `AAFTF pipeline_short`).
* **long:** ONT reads only (mirrors `AAFTF pipeline_long`).
* **hybrid:** ONT assembly polished with the Illumina reads (mirrors `AAFTF pipeline_hybrid`).

---

## 2. Directory Structure & Naming Conventions
All scripts, workflows, and outputs follow this layout (repository root):

```text
./
├── aaftf_benchmark_specification.md   # this document
├── README.md
├── pixi.toml / pixi.lock              # tool environments (default, eval, busco)
├── config/
│   ├── samples.tsv                 # simulated samples: organism, reference, accessions, genome_size, phylum, taxid, lineage, mito_seed
│   ├── real_samples.tsv            # real datasets: samples.tsv columns + read1, read2, longreads, longread_type, sra_short, sra_long
│   ├── contaminants.tsv            # organism, reference, accessions, group, assembly (source assembly, informational)
│   └── params.yaml                 # coverage, spike levels, group weights, tracks, partitions, ...
├── data/                           # not tracked in git
│   ├── sim/                        # simulated datasets
│   │   ├── references/
│   │   │   ├── fungus/             # downloaded genomes; tagged copies with headers >FUNGUS_{seqid}
│   │   │   ├── contaminants/       # downloaded genomes; tagged copies with headers >CONTAM_{taxon}_{seqid}
│   │   │   └── combined/           # {sample}_with_contaminants.fa: sample + all contaminants, for minimap2 contig labelling
│   │   └── reads/                  # simulated reads, composition tables and .sim_hash fingerprints
│   └── real/                       # real datasets
│       ├── references/             # same-strain reference genomes (+ tagged copies)
│       └── reads/                  # sequenced reads (Illumina and/or ONT/PacBio)
├── workflow/
│   ├── main.nf                     # Nextflow orchestration workflow
│   ├── nextflow.config             # resources and profiles (slurm, pixi)
│   ├── bin -> scripts              # put on PATH by Nextflow
│   ├── assets/NO_FILE              # placeholder for optional inputs
│   └── scripts/
│       ├── fetch_samples.sh                # download fungal references (samples.tsv)
│       ├── fetch_contaminants.sh           # download contaminant references (contaminants.tsv)
│       ├── fetch_reads.sh                  # download real reads from SRA (real_samples.tsv)
│       ├── format_reference_headers.py     # FUNGUS_/CONTAM_ headers + depth= copy-number tags
│       ├── make_readcounts.py              # exact InSilicoSeq read counts per sequence
│       ├── make_long_quantities.py         # Badread quantity (bases) per contaminant
│       ├── timeit.sh                       # GNU time wrapper for every assembler / AAFTF call
│       ├── calculate_decon_metrics.py
│       └── aggregate_benchmark_results.py
└── results/                        # not tracked in git
    ├── 01_checkpoints/             # intermediate assembly FASTAs per dataset and ablation stage
    ├── 02_quast/                   # QUAST assessment directories
    ├── 03_completeness/            # compleasm and BUSCO outputs
    ├── 04_decontamination/         # PAF alignments and per-checkpoint metrics
    ├── pipeline_info/              # Nextflow trace and report
    ├── final_report.tsv            # master aggregated results matrix (tidy, long format)
    ├── precision_recall_f1.tsv     # decontamination metrics for all checkpoints
    ├── time_summary.tsv            # wall time, peak RSS, CPU %, exit status per command
    └── read_composition.tsv        # realised fungal / contaminant amounts per dataset
```

**Columns.** All three TSVs start with `organism` (NCBI organism name), `reference` (local FASTA) and `accessions` (NCBI accessions, optional `:copies`), followed by file-specific columns.

**Names.** Sample and contaminant names are the NCBI organism names (e.g. `Neurospora crassa OR74A`); file names and FASTA headers use a safe key with non-alphanumeric runs replaced by `_` (e.g. `Neurospora_crassa_OR74A`). Datasets are named `{sample}_spike{pct}_{track}` (simulated) or `{sample}_real_{track}` (real reads).

---

## 3. Dataset Generation Protocol

### 3.1 Reference Genomes and Header Normalization
References are listed in the TSVs with NCBI accessions and downloaded by `fetch_samples.sh` (GCF_/GCA_ assemblies via NCBI Datasets, single sequences via E-utilities) and `fetch_contaminants.sh`. An accession may carry a copy number per nuclear genome as `acc:copies` (default 1), e.g. a mitochondrion `NC_012920.1:200`.

| organism (NCBI organism name) | assembly | phylum | BUSCO lineage |
|---|---|---|---|
| *Neurospora crassa* OR74A | GCF_000182925.2 (incl. mito NC_026614.1) | Ascomycota | sordariomycetes_odb12.2 |
| *Saccharomyces cerevisiae* S288C | GCF_000146045.2 (incl. mito NC_001224.1) | Ascomycota | saccharomycetes_odb12.2 |
| *Cryptococcus neoformans* H99 | GCF_000149245.1 (incl. mito NC_018792.1) | Basidiomycota | tremellomycetes_odb12.2 |

Before simulating reads, every FASTA header is tagged to make downstream mapping and decontamination evaluation deterministic:
* Fungal reference FASTA: rename all sequences to `^FUNGUS_.*`
* Contaminant reference FASTA (bacteria, host, vector, fungal cross-contaminant): rename to `^CONTAM_{TAXON}_.*`
* Copy numbers are written as `depth=N` after the ID (read natively by Badread; used by `make_readcounts.py` for InSilicoSeq).
* Combined FASTA (evaluation only, never used for assembly): concatenate tagged fungal and contaminant references:
  `cat data/sim/references/fungus/{sample}.fa data/sim/references/contaminants/*.fa > data/sim/references/combined/{sample}_with_contaminants.fa`

### 3.2 Real Datasets
Real runs of the exact strains with a reference genome are listed in `config/real_samples.tsv` (reference accessions as in `samples.tsv`, plus local read paths and SRA run accessions) and downloaded with `fetch_samples.sh config/real_samples.tsv` and `fetch_reads.sh`. The available reads decide the tracks (short: read1/read2; long: longreads; hybrid: both). Real datasets run the same checkpoints and the QUAST, completeness and timing evaluation; contaminant precision/recall (§5.3) needs the simulated answer key and is computed for simulated datasets only.

### 3.3 Read Simulation

#### Track A: Short-Read Simulation (`InSilicoSeq`)
Paired-end NovaSeq reads (2 × 150 bp) at `short_coverage` (default 50×) on the fungus. Exact read counts per sequence (length × copy number) are given with `--readcount_file`; `iss --coverage` names a coverage distribution, not a depth.
```bash
make_readcounts.py --fungus {fungus.fa} --contam {contaminants} --coverage 50 --percent {pct} ...
iss generate --genomes {fungus.fa} --readcount_file fungus.readcount.tsv --model novaseq \
  --output fungus --cpus {threads} --seed 42 --compress
```
Output: `data/sim/reads/{sample}_spike{pct}_R1.fastq.gz`, `_R2.fastq.gz`.

#### Track B: Long-Read Simulation (`Badread`)
ONT R10.4.1/Q20+-like reads at `long_coverage` (default 35×), settings in `badread_args`:
```bash
badread simulate --reference {fungus.fa} --quantity 35x \
  --length 15000,8000 --error_model nanopore2023 --qscore_model nanopore2023 \
  --identity 98.5,99.8,99.2 --chimeras 1.0 --glitches 10000,25,25 --seed 42 | gzip > fungus_ont.fastq.gz
```
Output: `data/sim/reads/{sample}_spike{pct}_ont.fastq.gz`. All sequences are simulated as linear.

#### Track C: Synthetic Decontamination Spike-In Panel
Contaminant reads are simulated from the tagged contaminant genomes and added to the fungal reads, so the fungal depth stays fixed:
contaminant amount = fungal amount × p / (100 − p) (reads for short, bases for long).
* **Spike-in levels:** `0%` (clean control), `1%`, `5%`, `15%` (`spike_levels`).
* **Group weights** (`group_weights`): 65% bacteria, 15% host, 10% vector, 10% fungal cross-contaminant; equal split between taxa of a group; length × copy number within an organism.

| group | organism (NCBI organism name) | accessions | rationale |
|---|---|---|---|
| bacteria | *Ralstonia pickettii* (K-288, type material) | NZ_CP066771.1, NZ_CP066772.2, NZ_CP066773.1 | DNA-extraction kit / reagent contaminant |
| bacteria | *Cutibacterium acnes* subsp. *acnes* (NBRC 107605) | NZ_AP019723.1 | skin / handling contaminant |
| host | *Homo sapiens* | NC_000021.9 (chr21), NC_012920.1:200 (mito) | handling contamination |
| vector | Escherichia phage phiX174 | NC_001422.1 | Illumina spike-in |
| vector | Cloning vector pUC19c | L09137.2 | cloning vector |
| fungal | *Penicillium chrysogenum* (P2niaD18) | CM002801.1 (chr IV), CM002802.1:50 (mito) | culture cross-contamination; same phylum as the Ascomycota samples |

No plant host is included: the samples are pure-culture lab strains. Organelle copy numbers are placeholders until measured from public WGS data (organelle:nuclear depth ratio).
* The realised composition is recorded in `results/read_composition.tsv`.
* Simulated reads are reused on reruns when their `.sim_hash` fingerprint (inputs, simulation settings, helper scripts) matches (`reuse_reads`).

---

## 4. Execution Tracks & Checkpoint Definitions

### Track 1: Sequential Ablation Pipeline
For each dataset, the intermediate FASTA after every step is saved to `results/01_checkpoints/{dataset}/` to isolate the value of each AAFTF command. cp02–cp04 are shared by all tracks; the long and hybrid tracks share one assembly and cleanup (labelled `_long`) and split at cp05.

1. **`cp00_raw_baseline`** (short only): standard raw assembler without AAFTF pre-filtering or post-processing:
   ```bash
   spades.py -1 {R1} -2 {R2} --isolate -o {outdir} -t {threads} -m {mem}
   ```
   FASTA output: `{outdir}/scaffolds.fasta`. No long-read baseline: AAFTF does no read cleanup before Flye, so a raw Flye run equals cp01.

2. **`cp01_mito_filtered_asm`** (short): read sanitation followed by core assembly:
   * `AAFTF trim --method bbduk -c {threads} -1 {R1} -2 {R2} -o reads_trimmed` → `reads_trimmed_1P/_2P.fastq.gz`
   * `AAFTF mito -1 reads_trimmed_1P.fastq.gz -2 reads_trimmed_2P.fastq.gz [-s {mito_seed}] -o mito_assembly.fasta` (a failure is tolerated)
   * `AAFTF filter -1 ... -2 ... --aligner bbduk -c {threads} -o reads_filtered [-s mito_assembly.fasta]` → `reads_filtered_filtered_1/_2.fastq.gz`
   * `AAFTF assemble -1 reads_filtered_filtered_1.fastq.gz -2 ..._2.fastq.gz -c {threads} --method spades -o {dataset}.cp01.fasta`

   **`cp01_flye_asm`** (long / hybrid):
   * `AAFTF assemble --method flye -lr {ont.fastq.gz} --longread_type nano-hq --genome_size {genome_size} -c {threads} -o {dataset}.cp01.fasta`

3. **`cp02_vecscreen`**:
   * `AAFTF vecscreen -i {cp01.fasta} -c {threads} -o {dataset}.cp02.fasta`

4. **`cp03_purged_{taxpurge,fcsgx}`**:
   * Branch A: `AAFTF taxpurge -i {cp02.fasta} -c {threads} -p {phylum} -o {dataset}.cp03_taxpurge.fasta`
   * Branch B (only with `fcs_gx_db`): `AAFTF fcs_gx_purge -i {cp02.fasta} -c {threads} -d {fcs_gx_db} -t {taxid} -o fcsgx/` (runs on the `fcsgx_partition`)

5. **`cp04_rmdup_{branch}`**:
   * `AAFTF rmdup -i {cp03.fasta} -c {threads} -ml 500 -o {dataset}.cp04_{branch}.fasta`

6. **cp05** (per branch):
   * short — **`cp05_sorted_{branch}`**: `AAFTF sort -i {cp04.fasta} -ml 500 -o {dataset}.cp05_{branch}.fasta` (no polishing, as in `AAFTF pipeline_short`)
   * long — **`cp05_polished_{branch}`**: `AAFTF polish --method racon -lr {ont}` → `AAFTF sort`
   * hybrid — **`cp05_polished_{branch}`**: `AAFTF polish --method racon -lr {ont}` → `AAFTF polish --method polypolish -1/-2 {filtered Illumina reads}` → `AAFTF sort`

---

## 5. Quantitative Evaluation Suite

Automated processes run four evaluation modules on every checkpoint FASTA:

### 5.1 Contiguity & Misassembly (`QUAST`)
Reference-guided `quast.py` with fungal parameters, against the tagged fungal reference:
```bash
quast.py -r data/{sim,real}/references/.../{sample}.fa --fungus --threads {threads} \
  --output-dir results/02_quast/{dataset}__{checkpoint} {assembly_checkpoint.fasta}
```
**Required Metrics to Extract (`transposed_report.tsv`, `contigs_reports/transposed_report_misassemblies.tsv`):**
* `Genome fraction (%)`
* `N50` and `NGA50`
* `auN` and `auNGA`
* `Total contigs` (total and ≥ 1000 bp)
* `Misassemblies` (Relocations, Translocations, Inversions)
* `Unaligned length`

### 5.2 Biological Completeness (`compleasm` & `BUSCO`)
Both tools in direct nucleotide mode against the sample's OrthoDB v12.2 lineage dataset, using the same lineage files: BUSCO reads the shared download folder (`busco_downloads`, read-only, `--offline`), and compleasm a writable library (`compleasm_library`) whose lineage folders are symlinks into it (`pixi run setup-compleasm`).

```bash
# Tool 1: compleasm (miniprot mode); lineage name without suffix, OrthoDB version from it
compleasm run -a {assembly_checkpoint.fasta} -o results/03_completeness/compleasm/{dataset}__{checkpoint} \
  -l {lineage} --odb odb12.2 -t {threads} -L {compleasm_library}

# Tool 2: BUSCO v6 (genome mode, miniprot engine)
busco -i {assembly_checkpoint.fasta} -m genome --miniprot -l {lineage}_odb12.2 \
  -o {dataset}__{checkpoint} --out_path results/03_completeness/busco --cpu {threads} \
  --download_path {busco_downloads} --offline
```
**Required Metrics to Extract:**
* Complete Single-Copy (S%)
* Complete Duplicated (D%) — *track reduction from checkpoint 00 to 04*
* Fragmented (F%)
* Missing (M%)

### 5.3 Decontamination Efficacy & Confusion Matrix (`minimap2`)
Map the assembly contigs of simulated datasets directly against `data/sim/references/combined/{sample}_with_contaminants.fa` to establish contig taxonomy without BLAST:
```bash
minimap2 -x asm5 --secondary=no -t {threads} \
  data/sim/references/combined/{sample}_with_contaminants.fa \
  {assembly_checkpoint.fasta} > results/04_decontamination/paf/{dataset}__{checkpoint}.paf
```

#### Parsing Logic for `calculate_decon_metrics.py`:
Parse the PAF file (`col 1: query_contig`, `col 3–4: query start/end`, `col 6: target_seq`):
* Assign each contig the class with the most aligned bases: `^FUNGUS_` → **Fungal**, `^CONTAM_` → **Contaminant**; contigs without hits → **Unassigned** (excluded from the confusion matrix).
* Every checkpoint: report fungal, contaminant and unassigned bp and the number of contaminant contigs.
* Compare the contig set before purging (`cp02`) vs. after purging (`cp03`), matched by contig name:
  * **True Positive (TP, bp):** contaminant contig bases successfully removed.
  * **False Negative (FN, bp):** contaminant contig bases remaining in the purged assembly.
  * **False Positive (FP, bp):** fungal contig bases mistakenly removed.
  * **True Negative (TN, bp):** fungal contig bases retained.
* Compute:
  $$\text{Precision} = \frac{TP}{TP + FP}$$
  $$\text{Recall} = \frac{TP}{TP + FN}$$
  $$F_1 = \frac{2 \cdot \text{Precision} \cdot \text{Recall}}{\text{Precision} + \text{Recall}}$$

Expected limitation: `taxpurge -p <phylum>` keeps every contig of the sample's phylum, so it should miss the *Penicillium* cross-contaminant in the Ascomycota samples but not in *C. neoformans*.

### 5.4 Computational Efficiency
Every assembler and AAFTF command is wrapped in GNU `/usr/bin/time -v` (`timeit.sh`); the Nextflow trace records the same per task:
* `Wall-clock runtime (seconds)`
* `Peak Resident Set Size (RSS memory, MB)`
* `CPU usage percentage`
* exit status

---

## 6. Implementation Deliverables

1. **Nextflow workflow (`workflow/main.nf`)** that parameterizes the ablation checkpoints, the three tracks and the evaluation steps across all samples and spike levels, with `slurm` (partition `epyc`; FCS-GX on `highmem`) and `pixi` profiles.
2. **`calculate_decon_metrics.py`** in Python using `pandas` to process the PAF files and emit `precision_recall_f1.tsv`.
3. **`aggregate_benchmark_results.py`** to combine QUAST `transposed_report.tsv`, compleasm/BUSCO summary files, decontamination metrics and execution times into a single tidy `final_report.tsv` (`dataset, sample, spike_pct, track, checkpoint, tool, metric, value`) suitable for plotting in R (`ggplot2`) or Python (`seaborn`).
4. **`pixi.toml` + `pixi.lock`** with pinned versions of `aaftf`, `insilicoseq`, `badread`, `quast`, `compleasm`, `busco` (v6), `minimap2`, `samtools` and the AAFTF tool chain, split into three environments (`default`, `eval`, `busco`) because they cannot be solved together.
