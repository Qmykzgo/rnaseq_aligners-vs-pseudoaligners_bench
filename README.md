# RNA-seq Alignment & Pseudoalignment Benchmarking Pipeline

This project evaluates the performance, accuracy, and reproducibility of alignment-based pipelines (**STAR**, **HISAT2** + **featureCounts**) against alignment-free pseudoaligners (**Salmon**, **kallisto**) across three distinct benchmark Aims.

## System Requirements & Prerequisites

> [!IMPORTANT]
> **Linux/WSL Environment Required:**
> Bioinformatics packages hosted on Bioconda (STAR, HISAT2, Salmon, kallisto, featureCounts) compile only for Linux and macOS. If you are on Windows, you **must** run these scripts in **WSL (Windows Subsystem for Linux)** or a Linux VM.

- **Conda package manager** (Miniconda or Anaconda installed and on your PATH)
- **GNU time utility** (installed at `/usr/bin/time` to measure CPU/RAM usage). On Debian/Ubuntu/WSL, run:
  ```bash
  sudo apt-get update && sudo apt-get install time -y
  ```

---

## Benchmark Structure

- **Aim 1: Technical Replicate Reproducibility**
  - Evaluates within-sample consistency (Pearson & Spearman correlations) on technical replicates from SEQC/MAQC-III Sample A and Sample B.
- **Aim 2: Quantification Accuracy**
  - Compares the log2 fold-change calculated by each tool (MAQCA vs. MAQCB) against the Everaert et al. 2017 whole-transcriptome RT-qPCR normalized reference dataset (**GSE83402**).
- **Aim 3: Novel Splice-Junction Recovery**
  - Assesses precision, recall, and F1-score of STAR (1-pass vs. 2-pass) and HISAT2 on chromosome 21 novel junctions using an indexing mask strategy.

---

## Execution Guide

Run the pipeline scripts sequentially:

### Step 0: Create Conda Environments
Creates the tool execution environment (`rnaseq_bench_env`) and metrics reporting environment (`rnaseq_metrics_env`).
```bash
bash 00_setup_conda_envs.sh
```

### Step 1: Download & Prepare References & Reads
Downloads the GRCh38 genome, Ensembl 110 annotations, ERCC-92 spike-ins, GSE83402 RT-qPCR truth table, and prefetches all 10 SRA datasets.
```bash
bash 01_prepare_data.sh
```

### Step 2: Run Alignment Pipelines (STAR / HISAT2)
Builds the full and chr21-masked indices, runs alignments, and quantifies using featureCounts.
```bash
bash 02_run_alignment.sh
```

### Step 3: Run Pseudoalignment Pipelines (Salmon / kallisto)
Builds Salmon (decoy-aware, gentrome-based) and kallisto indices, and runs quantifications.
```bash
bash 03_run_pseudoalignment.sh
```

### Step 4: Parse & Aggregate Benchmark Metrics
Runs the python analysis scripts to parse all aligner logs, abundance files, and execution resource summaries.
```bash
conda activate rnaseq_metrics_env
python 04_aggregate_metrics.py
```

---

## Output Metrics & Summaries

The aggregated results are saved in `results/summaries/`:

1. `aim1_reproducibility.csv`: Correlation analysis of replicates.
2. `aim2_accuracy.csv`: Accuracy metrics (Pearson $r$, Spearman $r$, RMSE) against Everaert qPCR.
3. `aim3_junction_recovery.csv`: Precision, recall, and F1 for chromosome 21 junction recovery.
4. `performance_profile.csv`: Complete resource footprint profile (wall time & peak RAM) for indexing, alignment, and quantification.
