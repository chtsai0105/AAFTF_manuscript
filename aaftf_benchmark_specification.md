# Specification: AAFTF Assembly Pipeline Internal Benchmarking Suite

## 1. Project Overview & Scope
This workflow implements an automated, reproducible benchmarking harness for the fungal genome assembly pipeline **AAFTF** (Automatic Assembly For The Fungi). 

### Benchmark Design Principles
1. **Self-Contained Evaluation:** Do NOT compare against external wrapper pipelines (e.g., Shovill, fspassemblypipeline). All benchmarks compare AAFTF intermediate stages against a naive raw assembler baseline or evaluate interchangeable internal modules.
2. **Deterministic Ground Truth:** All read datasets are simulated from finished, chromosome-level reference genomes using `InSilicoSeq` (short reads) and `Badread` (long reads).
3. **Decoupled Verification:** Evaluation relies strictly on nucleotide-level validation (`QUAST`, `compleasm`, `BUSCO -m genome --miniprot`, and `minimap2` alignment against known reference targets). Do NOT invoke *ab initio* gene annotation pipelines (e.g., Augustus, BRAKER, funannotate).

---

## 2. Directory Structure & Naming Conventions
All scripts, workflows, and outputs must follow this standardized layout:

```text
benchmark_suite/
├── config/
│   ├── samples.tsv                 # Sample IDs, phylum, reference paths, and genome sizes
│   ├── contaminants.tsv            # Paths and taxonomy of spike-in contaminant references
│   └── params.yaml                 # Core configuration parameters (threads, k-mer sizes, etc.)
├── data/
│   ├── references/                 # Source T2T/finished FASTA genomes
│   │   ├── fungus/                 # Formatted headers: >FUNGUS_{seqid}
│   │   └── contaminants/           # Formatted headers: >CONTAM_{taxon}_{seqid}
│   ├── reads/                      # Simulated reads (paired-end Illumina & ONT/HiFi)
│   └── truth_references/          # Combined concatenated truth genomes for minimap2
├── workflow/
│   ├── Snakefile                   # (or main.nf) Orchestration workflow
│   ├── envs/                       # Conda / Pixi environment YAMLs
│   │   └── benchmark_env.yaml
│   └── scripts/
│       ├── format_reference_headers.py
│       ├── calculate_decon_metrics.py
│       └── aggregate_benchmark_results.py
└── results/
    ├── 01_checkpoints/             # Intermediate assembly FASTAs per ablation stage
    ├── 02_quast/                   # QUAST assessment directories
    ├── 03_completeness/            # compleasm and BUSCO outputs
    ├── 04_decontamination/         # PAF alignments and Confusion Matrix metrics
    └── final_report.tsv            # Master aggregated results matrix
```

---

## 3. Dataset Generation Protocol

### 3.1 Reference Genome Header Normalization
Before simulating reads, every FASTA header must be explicitly tagged to make downstream mapping and decontamination evaluation deterministic:
* Fungal reference FASTA: Rename all sequences to match pattern: `^FUNGUS_.*`
* Contaminant reference FASTA (Bacteria, Host Plant, Vector): Rename to match pattern: `^CONTAM_{TAXON}_.*`
* Combined truth FASTA: Concatenate tagged fungal and contaminant references:
  `cat data/references/fungus/{sample}.fa data/references/contaminants/*.fa > data/truth_references/{sample}_combined_truth.fa`

### 3.2 Read Simulation

#### Track A: Short-Read Simulation (`InSilicoSeq`)
Generate paired-end NovaSeq reads ($2 \times 150$ bp) at target coverage depth:
```bash
iss generate \
  --genomes data/references/fungus/{sample}.fa \
  --model novaseq \
  --coverage 50 \
  --output data/reads/{sample}_fungus_50x \
  --cpus {threads}
```

#### Track B: Long-Read Simulation (`Badread`)
Generate long reads (ONT R10.4.1/Q20+ profile) with realistic structural artifacts:
```bash
badread simulate \
  --reference data/references/fungus/{sample}.fa \
  --quantity 35x \
  --length 15000,8000 \
  --error_model nanopore2023 \
  --qscore_model nanopore2023 \
  --identity 98.5,99.8,99.2 \
  --chimeras 1.0 \
  --glitches 10000,25,25 | gzip > data/reads/{sample}_ont_35x.fastq.gz
```

#### Track C: Synthetic Decontamination Spike-In Panel
Simulate contaminant reads from tagged contaminant genomes and mix with the target fungus at exact read ratios:
* **Spike-in levels:** `1%`, `5%`, and `15%` contaminant read abundance.
* **Composition:** 70% Gram-negative/Gram-positive bacteria (*Pseudomonas*, *Burkholderia*, *Staphylococcus*), 20% host plant/human (*Arabidopsis* / *Homo sapiens*), 10% cloning vector / PhiX.
* Use `InSilicoSeq` with `--abundance_file` or simulate separately and merge with fixed read counts to ensure the ground truth percentage is exact.

---

## 4. Execution Tracks & Checkpoint Definitions

### Track 1: 6-Stage Sequential Ablation Pipeline
For each sample, capture the intermediate FASTA outputs to isolate the value of each AAFTF command:

1. **`checkpoint_00_raw_baseline`**:
   Execute standard raw assembler without AAFTF pre-filtering or post-processing:
   ```bash
   spades.py -1 {R1} -2 {R2} --isolate -o {outdir} -t {threads} -m {mem}
   ```
   FASTA output: `{outdir}/scaffolds.fasta`

2. **`checkpoint_01_mito_filtered_asm`**:
   Run pre-assembly read sanitation followed by core assembly:
   * `AAFTF trim --method bbduk -c {threads} -l {R1} -r {R2} -o reads_trimmed`
   * `AAFTF mito -l reads_trimmed_1.fq.gz -r reads_trimmed_2.fq.gz -o mito_assembly -c {threads}`
   * `AAFTF filter -l reads_trimmed_1.fq.gz -r reads_trimmed_2.fq.gz --aligner bbmap -c {threads} -o reads_filtered`
   * `AAFTF assemble -l reads_filtered_1.fq.gz -r reads_filtered_2.fq.gz -c {threads} --method spades -o 01_assembled.fasta`

3. **`checkpoint_02_vecscreen`**:
   * `AAFTF vecscreen -i 01_assembled.fasta -c {threads} -o 02_vecscreen.fasta`

4. **`checkpoint_03_purged`**:
   * Branch A: `AAFTF sourpurge -i 02_vecscreen.fasta -c {threads} -o 03_sourpurge.fasta`
   * Branch B: `AAFTF fcs_gx_purge -i 02_vecscreen.fasta -c {threads} -o 03_fcsgx_purge.fasta`

5. **`checkpoint_04_rmdup`**:
   * `AAFTF rmdup -i 03_purged.fasta -c {threads} -o 04_rmdup.fasta`

6. **`checkpoint_05_polished_sorted`**:
   * `AAFTF polish -i 04_rmdup.fasta -l reads_filtered_1.fq.gz -r reads_filtered_2.fq.gz --method polca -c {threads} -o 05_polished.fasta`
   * `AAFTF sort -i 05_polished.fasta -o 05_final_sorted.fasta`

---

## 5. Quantitative Evaluation Suite

Implement automated rules to run four evaluation modules on every checkpoint FASTA:

### 5.1 Contiguity & Misassembly (`QUAST`)
Use standard reference-guided `quast.py` with fungal parameters:
```bash
quast.py \
  -r data/references/fungus/{sample}.fa \
  --fungus \
  --threads {threads} \
  --output-dir results/02_quast/{sample}_{checkpoint} \
  {assembly_checkpoint.fasta}
```
**Required Metrics to Extract (`report.tsv`):**
* `Genome fraction (%)`
* `N50` and `NGA50`
* `auN` and `auNGA`
* `Total contigs` (total and $\ge 1000$ bp)
* `Misassemblies` (Relocations, Translocations, Inversions)
* `Unaligned length`

### 5.2 Biological Completeness (`compleasm` & `BUSCO`)
Run both tools in direct nucleotide mode against the relevant lineage OrthoDB dataset:

```bash
# Tool 1: compleasm (miniprot mode)
compleasm run \
  -a {assembly_checkpoint.fasta} \
  -o results/03_completeness/compleasm/{sample}_{checkpoint} \
  -l {lineage_odb10} \
  -t {threads}

# Tool 2: BUSCO v5 (genome mode, miniprot engine)
busco \
  -i {assembly_checkpoint.fasta} \
  -m genome \
  -l {lineage_odb10} \
  -o {sample}_{checkpoint} \
  --out_path results/03_completeness/busco \
  --cpu {threads} \
  --download_path {busco_downloads}
```
**Required Metrics to Extract:**
* Complete Single-Copy ($S\%$)
* Complete Duplicated ($D\%$) — *Track reduction from checkpoint 00 to 04*
* Fragmented ($F\%$)
* Missing ($M\%$)

### 5.3 Decontamination Efficacy & Confusion Matrix (`minimap2`)
Map the assembly contigs directly against `combined_truth_ref.fasta` to establish contig taxonomy without BLAST:
```bash
minimap2 -x asm5 --secondary=no -t {threads} \
  data/truth_references/{sample}_combined_truth.fa \
  {assembly_checkpoint.fasta} > results/04_decontamination/{sample}_{checkpoint}.paf
```

#### Parsing Logic for `calculate_decon_metrics.py`:
Parse the resulting PAF file (`col 1: query_contig`, `col 2: query_len`, `col 6: target_seq`):
* Classify contigs matching `^FUNGUS_` as **Fungal**, and `^CONTAM_` as **Contaminant**.
* Compare contig set before purging (`checkpoint_02`) vs. after purging (`checkpoint_03`):
  * **True Positive ($TP$ in bp):** Contaminant contig bases successfully removed.
  * **False Negative ($FN$ in bp):** Contaminant contig bases remaining in purged assembly.
  * **False Positive ($FP$ in bp):** Fungal contig bases mistakenly removed.
  * **True Negative ($TN$ in bp):** Fungal contig bases retained.
* Compute formulas:
  $$\text{Precision} = \frac{TP}{TP + FP}$$
  $$\text{Recall} = \frac{TP}{TP + FN}$$
  $$F_1 = \frac{2 \cdot \text{Precision} \cdot \text{Recall}}{\text{Precision} + \text{Recall}}$$

### 5.4 Computational Efficiency
Wrap all subcommand executions with benchmark resource profilers (e.g., Snakemake `benchmark:` directive or GNU `/usr/bin/time -v`):
* `Wall-clock runtime (seconds / minutes)`
* `Peak Resident Set Size (RSS Memory in MB/GB)`
* `CPU Usage Percentage`

---

## 6. Implementation Deliverables for the Agent

When writing the code for this repository:
1. **Provide a working Snakemake (`Snakefile`) or Nextflow (`main.nf`) script** that parameterizes the 6 ablation checkpoints and evaluation steps across all samples.
2. **Implement `calculate_decon_metrics.py`** in Python using `polars` or `pandas` to process the PAF file and emit `precision_recall_f1.tsv`.
3. **Implement `aggregate_benchmark_results.py`** to combine QUAST `transposed_report.tsv`, compleasm/BUSCO summary files, decontamination metrics, and execution times into a single tidy `final_report.tsv` suitable for plotting in R (`ggplot2`) or Python (`seaborn`).
4. **Include a `conda_env.yaml` / `pixi.toml`** containing exact versions of `aaftf`, `insilicoseq`, `badread`, `quast`, `compleasm`, `busco`, `minimap2`, and `samtools`.