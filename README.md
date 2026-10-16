# AAFTF internal benchmarking suite

Implements [aaftf_benchmark_specification.md](aaftf_benchmark_specification.md) in Nextflow.

```
pixi install --all                       # tools: see pixi.toml (environments default, eval, busco)
pixi run fetch-samples && pixi run fetch-contaminants
pixi run nextflow run workflow/main.nf -params-file config/params.yaml -profile pixi,slurm
```

With `-profile slurm`, jobs go to the `partition` set in `config/params.yaml` (default `epyc`), and FCS-GX
(`CP03_FCSGX`, ~470 GB database in RAM) goes to `fcsgx_partition` (default `highmem`). Override either on the
command line, e.g. `--partition batch`.

1. Fungal references are defined in `config/samples.tsv`; `workflow/scripts/fetch_samples.sh` downloads each row's
   `accessions` into its `reference` path. A `GCF_`/`GCA_` accession fetches the whole assembly (NCBI Datasets);
   other accessions fetch single sequences, e.g. a mitogenome missing from the assembly (skipped when the assembly
   already has it). `acc:copies` sets a copy number as for contaminants: the mitochondria are at a placeholder 50.
   Both fetch scripts skip existing files; `--force` re-downloads. The contaminant panel is defined in `config/contaminants.tsv`
   (`taxon, group, reference, accessions`); `workflow/scripts/fetch_contaminants.sh` downloads each row's NCBI
   nuccore accessions into its `reference` path (~57 Mb: *Ralstonia pickettii*, *Cutibacterium acnes*, human chr21 + mito,
   PhiX, pUC19, *Penicillium chrysogenum* chr IV + mito). No plant host: the samples are pure-culture lab strains, so
   plant DNA would only be realistic for a plant pathogen sequenced from host tissue.
   Edit the TSV to change the panel; rows without accessions must point to a local FASTA.
   Eukaryotic contaminants carry their organelles in the same row (human chr21 + mitochondrion; *Penicillium*
   chr IV + mitochondrion), because within one organism the nuclear:organelle read ratio is set by
   biology, not by the scenario. An accession can carry a copy number per nuclear genome as `acc:copies`
   (default 1), e.g. `NC_000021.9,NC_012920.1:200`. It is written as `depth=` in the FASTA
   header; Badread reads it natively and `make_readcounts.py` scales InSilicoSeq read counts by it.
   **The current copy numbers (human mito 200, fungal mitos 50) are placeholders**;
   replace them with organelle:nuclear depth ratios measured from real WGS data before the final run.
   All sequences are simulated as linear (no reads span a circular genome's origin).
   The fungal cross-contaminant (*Penicillium chrysogenum* P2niaD18 chr IV + mitochondrion) is in the same phylum
   (Ascomycota, different class) as *N. crassa* and *S. cerevisiae*, and a different phylum from *C. neoformans*.
   `taxpurge -p <phylum>` keeps every contig of the sample's phylum, so it is expected to miss this contaminant in
   the Ascomycota samples; this tests the limit of phylum-level purging against FCS-GX.
   The headers are rewritten to `FUNGUS_*` / `CONTAM_<taxon>_*` automatically.
2. Reads are simulated with InSilicoSeq at 0/1/5/15% contaminant reads, split 65/15/10/10 between bacteria, host (human), vector and a fungal cross-contaminant (`group_weights`).
   Exact read counts per sequence go through `--readcount_file`; the realised mix is written to `results/read_composition.tsv`.
3. Ablation checkpoints `cp00`–`cp05` (branch A `taxpurge`, branch B `fcs_gx_purge` when `fcs_gx_db` is set).
   Short-track cp05 is `AAFTF sort` only (`cp05_sorted_*`), matching `AAFTF pipeline_short`, which does not polish.
   `tracks` (params.yaml or `--tracks`) picks what runs, as a comma list:
   - `short` (default): the six checkpoints above
   - `long`: Badread ONT reads → `cp01_flye_asm` (AAFTF assemble --method flye) → cp02–cp04 → cp05 racon → sort
   - `hybrid`: same assembly and cleanup as `long`; cp05 = racon, then polypolish with cleaned Illumina reads
   `long` and `hybrid` share one flye assembly and cleanup (labelled `<sample>_spike<pct>_long`).
   Illumina reads are only simulated when `short` or `hybrid` is chosen; with `hybrid` but not `short`, they are
   cleaned (trim/mito/filter) but not assembled. No cp00 for long reads: a raw flye run would equal cp01.
   Long contaminant reads use the same % (by bases) and group weights.

   ```
   nextflow run workflow/main.nf -params-file config/params.yaml --tracks long
   ```
   Simulated reads are saved in `data/reads/` with a `<sample>_spike<pct>_<short|long>.sim_hash` fingerprint of
   the settings that made them. With `reuse_reads: true` (default) a rerun reuses a set only when all of its files
   exist and the stored hash matches the current one; otherwise that set is simulated again and overwritten.
   The hash covers: the fungal and contaminant FASTAs (path, size, modification time), `contaminants.tsv`,
   `group_weights`, the spike level, `short_coverage`/`read_length`/`iss_model` (short) or
   `long_coverage`/`badread_args` (long), and the helper scripts that compute read counts and rename headers.
   It does not cover the simulator versions or the fixed seeds; set `reuse_reads: false` to force new reads.
4. Every checkpoint is evaluated with QUAST, compleasm, BUSCO (`--miniprot`) and minimap2 against the truth reference.

Outputs (`results/`):
- `final_report.tsv`: tidy long table (`dataset, sample, spike_pct, track, checkpoint, tool, metric, value`); datasets are `<sample>_spike<pct>_<track>`
- `precision_recall_f1.tsv`: contig composition for every checkpoint, plus TP/FN/FP/TN, precision, recall and F1 for cp02 → cp03
- `time_summary.tsv`: wall time, peak RSS, CPU % and exit status for every AAFTF / SPAdes call
- `02_quast/`, `03_completeness/`, `04_decontamination/`: raw tool outputs

## Differences from the spec (current AAFTF CLI)
- Reads use `-1/-2` (not `-l/-r`). trim writes `_1P/_2P.fastq.gz`; filter writes `_filtered_1/2.fastq.gz`.
- `filter --aligner bbmap` is not a valid choice, so `bbduk` is the default (`filter_aligner`).
- `sourpurge` is now `taxpurge`, which needs `-p <phylum>` from samples.tsv.
- `fcs_gx_purge` needs `-d <db> -t <taxid>` and writes into an output directory.
- `mito` output is passed to `filter -s` to screen mitochondrial reads. If NOVOPlasty fails (for example, the reference has no mitogenome), the run continues without it.
- `iss --coverage 50` is a distribution name in InSilicoSeq, not a depth; depth is set with `--readcount_file`.
