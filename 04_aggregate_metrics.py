#!/usr/bin/env python3
# ==============================================================================
# 04_aggregate_metrics.py
# Purpose : Parse alignment logs and quantification outputs, then produce
#           four benchmark summary CSV files:
#
#   summaries/aim1_reproducibility.csv   — Pearson / Spearman r (tool vs tool)
#   summaries/aim2_accuracy.csv          — log2FC vs qPCR truth (each tool)
#   summaries/aim3_junction_recovery.csv — STAR 2-pass / 1-pass / HISAT2
#                                          precision, recall, F1 on chr21 junctions
#   summaries/performance_profile.csv    — Wall time & peak RAM for every step
#
# Requires: pandas, numpy, scipy  (conda activate rnaseq_metrics_env)
# Usage   : python 04_aggregate_metrics.py
# ==============================================================================

import os
import re
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

# ==============================================================================
# PATHS
# ==============================================================================
BASE_DIR   = Path(__file__).resolve().parent
DATA_DIR   = BASE_DIR / "data"
REF_DIR    = DATA_DIR / "reference"
ALN_DIR    = BASE_DIR / "results" / "alignment"
PSA_DIR    = BASE_DIR / "results" / "pseudoalignment"
LOG_DIR    = BASE_DIR / "results" / "logs"
TRUTH_DIR  = DATA_DIR / "truth_tables"
OUT_DIR    = BASE_DIR / "results" / "summaries"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TX2GENE_PATH       = REF_DIR / "tx2gene.tsv"
CHR21_TRUTH_PATH   = REF_DIR / "chr21_true_junctions.tsv"
QPCR_TRUTH_PATH    = TRUTH_DIR / "GSE83402_qPCR_normalized.txt"
BENCH_TRUTH_PATH   = TRUTH_DIR / "RNAontheBENCH_truth_counts.csv"

# Sample sets
AIM1_SAMPLES = {
    "SRR896663": "SampleA_Rep1",
    "SRR896679": "SampleA_Rep2",
    "SRR896743": "SampleB_Rep1",
    "SRR896759": "SampleB_Rep2",
}

AIM2B_SAMPLES = {
    "SRR896663": "MAQCA",
    "SRR896743": "MAQCB",
}

AIM2A_SAMPLES = {
    "SRR1643426": "WBS_Rep1",
    "SRR1643427": "WBS_Rep2",
    "SRR1643434": "7dup_Rep1",
}

AIM3_SAMPLE = "SRR1293902"


# ==============================================================================
# UTILITY: parse /usr/bin/time -v log
# ==============================================================================
def parse_time_log(filepath: Path) -> dict:
    """Return wall-clock seconds and peak RSS (MB) from a /usr/bin/time -v log."""
    result = {"runtime_s": None, "peak_ram_mb": None}
    if not filepath.exists():
        return result

    with open(filepath) as fh:
        for line in fh:
            # Wall clock — format: h:mm:ss or m:ss
            if "Elapsed (wall clock) time" in line:
                time_str = line.split("): ")[-1].strip()
                parts = [float(x) for x in time_str.split(":")]
                if len(parts) == 3:
                    result["runtime_s"] = parts[0] * 3600 + parts[1] * 60 + parts[2]
                elif len(parts) == 2:
                    result["runtime_s"] = parts[0] * 60 + parts[1]
            # Peak RSS
            elif "Maximum resident set size (kbytes):" in line:
                kb = float(line.split(":")[-1].strip())
                result["peak_ram_mb"] = kb / 1024.0
    return result


# ==============================================================================
# UTILITY: load featureCounts output → Series(gene_id → count)
# ==============================================================================
def load_featurecounts(filepath: Path) -> pd.Series:
    """Parse featureCounts output into a gene→count Series."""
    if not filepath.exists():
        return pd.Series(dtype=float)

    df = pd.read_csv(
        filepath,
        sep="\t",
        comment="#",
        index_col=0,          # Geneid
        usecols=lambda c: c not in ("Chr", "Start", "End", "Strand", "Length"),
    )
    # The last column is the count column (full BAM path as header)
    counts = df.iloc[:, -1].astype(float)
    counts.name = filepath.parent.name  # label by tool directory
    return counts


# ==============================================================================
# UTILITY: load Salmon quant.sf → Series(gene_id → TPM)
# ==============================================================================
def load_salmon(quant_dir: Path, tx2gene: pd.Series) -> pd.Series:
    """Aggregate transcript-level Salmon TPM to gene level."""
    sf_path = quant_dir / "quant.sf"
    if not sf_path.exists():
        return pd.Series(dtype=float)

    df = pd.read_csv(sf_path, sep="\t", index_col="Name")
    df.index = df.index.map(lambda t: t.split("|")[0])  # strip version if present
    df["gene_id"] = df.index.map(tx2gene)
    return df.groupby("gene_id")["TPM"].sum()


# ==============================================================================
# UTILITY: load kallisto abundance.tsv → Series(gene_id → TPM)
# ==============================================================================
def load_kallisto(quant_dir: Path, tx2gene: pd.Series) -> pd.Series:
    """Aggregate transcript-level kallisto TPM to gene level."""
    ab_path = quant_dir / "abundance.tsv"
    if not ab_path.exists():
        return pd.Series(dtype=float)

    df = pd.read_csv(ab_path, sep="\t", index_col="target_id")
    df.index = df.index.map(lambda t: t.split("|")[0])
    df["gene_id"] = df.index.map(tx2gene)
    return df.groupby("gene_id")["tpm"].sum()


# ==============================================================================
# UTILITY: load tx2gene mapping
# ==============================================================================
def load_tx2gene() -> pd.Series:
    if not TX2GENE_PATH.exists():
        print(f"[WARN] tx2gene not found at {TX2GENE_PATH}. Transcript-level TPMs "
              "will not be collapsed to gene level.")
        return pd.Series(dtype=str)
    df = pd.read_csv(TX2GENE_PATH, sep="\t", header=None, names=["tx", "gene"])
    return df.set_index("tx")["gene"]


# ==============================================================================
# BUILD GENE-LEVEL EXPRESSION MATRICES FOR EACH TOOL
# ==============================================================================
def build_expression_matrices(samples: dict, tx2gene: pd.Series) -> dict[str, pd.DataFrame]:
    """
    Returns a dict of DataFrames keyed by tool name.
    Each DataFrame has shape (genes × samples), values = TPM / raw counts.
    """
    matrices = {
        "STAR_featureCounts":  {},
        "HISAT2_featureCounts": {},
        "Salmon":              {},
        "kallisto":            {},
    }

    for srr, label in samples.items():
        # featureCounts (raw counts from aligner BAMs)
        fc_star   = load_featurecounts(ALN_DIR / "star"   / f"{srr}_featureCounts.txt")
        fc_hisat2 = load_featurecounts(ALN_DIR / "hisat2" / f"{srr}_featureCounts.txt")

        if not fc_star.empty:
            matrices["STAR_featureCounts"][label]   = fc_star
        if not fc_hisat2.empty:
            matrices["HISAT2_featureCounts"][label] = fc_hisat2

        # Salmon / kallisto (TPM, already gene-collapsed if tx2gene is available)
        salmon_tpm   = load_salmon(PSA_DIR  / "salmon"  / srr, tx2gene)
        kallisto_tpm = load_kallisto(PSA_DIR / "kallisto" / srr, tx2gene)

        if not salmon_tpm.empty:
            matrices["Salmon"][label]   = salmon_tpm
        if not kallisto_tpm.empty:
            matrices["kallisto"][label] = kallisto_tpm

    # Convert inner dicts to DataFrames
    out = {}
    for tool, cols in matrices.items():
        if cols:
            out[tool] = pd.DataFrame(cols).dropna()
        else:
            out[tool] = pd.DataFrame()
    return out


# ==============================================================================
# SECTION A  —  AIM 1: REPRODUCIBILITY
# ==============================================================================
def compute_aim1_reproducibility(tx2gene: pd.Series) -> pd.DataFrame:
    """
    For each tool, compute Pearson r and Spearman r between:
      - SampleA_Rep1 vs SampleA_Rep2  (within-sample technical replicates)
      - SampleB_Rep1 vs SampleB_Rep2
    Use log2(count+1) or log2(TPM+1) for correlation.
    """
    print("[AIM 1] Computing reproducibility metrics …")

    mats = build_expression_matrices(AIM1_SAMPLES, tx2gene)
    rows = []

    replicate_pairs = [
        ("SampleA_Rep1", "SampleA_Rep2", "SEQC_SampleA"),
        ("SampleB_Rep1", "SampleB_Rep2", "SEQC_SampleB"),
    ]

    for tool, mat in mats.items():
        if mat.empty:
            print(f"  [WARN] No data for {tool} — skipping.")
            continue

        for c1, c2, pair_label in replicate_pairs:
            if c1 not in mat.columns or c2 not in mat.columns:
                continue

            # log2-transform to linearize dynamic range
            x = np.log2(mat[c1].values + 1)
            y = np.log2(mat[c2].values + 1)

            # Filter genes with at least 1 count in either replicate
            mask = (x > 0) | (y > 0)
            x, y = x[mask], y[mask]

            pearson_r,  pearson_p  = stats.pearsonr(x, y)
            spearman_r, spearman_p = stats.spearmanr(x, y)

            rows.append({
                "Tool":            tool,
                "Replicate_Pair":  pair_label,
                "N_Genes":         int(mask.sum()),
                "Pearson_r":       round(pearson_r,  4),
                "Pearson_p":       f"{pearson_p:.2e}",
                "Spearman_r":      round(spearman_r, 4),
                "Spearman_p":      f"{spearman_p:.2e}",
            })

    df = pd.DataFrame(rows)
    out_path = OUT_DIR / "aim1_reproducibility.csv"
    df.to_csv(out_path, index=False)
    print(f"[OK] Aim 1 reproducibility → {out_path}  ({len(df)} rows)")
    return df


# ==============================================================================
# SECTION B  —  AIM 2: ACCURACY AGAINST ORTHOGONAL TRUTH
# ==============================================================================
def load_qpcr_truth() -> pd.DataFrame | None:
    """
    Load Everaert et al. 2017 qPCR normalized Cq values.
    Expected columns: GeneSymbol (or HGNC), NormCt_MAQCA, NormCt_MAQCB
    Lower NormCt = higher expression.  log2FC = NormCt_B - NormCt_A (sign flip).
    Returns DataFrame with index = gene symbol and column log2FC_qPCR.
    """
    if not QPCR_TRUTH_PATH.exists():
        print(f"  [WARN] qPCR truth not found: {QPCR_TRUTH_PATH}")
        return None

    df = pd.read_csv(QPCR_TRUTH_PATH, sep="\t")
    # Auto-detect column names (handle different releases of the file)
    gene_col  = next((c for c in df.columns if "gene" in c.lower()), df.columns[0])
    cta_col   = next((c for c in df.columns if "maqca" in c.lower() or "_a" in c.lower()), None)
    ctb_col   = next((c for c in df.columns if "maqcb" in c.lower() or "_b" in c.lower()), None)

    if cta_col is None or ctb_col is None:
        print("  [WARN] Could not identify MAQCA/MAQCB columns in qPCR truth table.")
        return None

    df = df[[gene_col, cta_col, ctb_col]].dropna()
    df.columns = ["gene", "NormCt_A", "NormCt_B"]
    df["gene"] = df["gene"].astype(str).str.upper()
    df["log2FC_qPCR"] = df["NormCt_B"] - df["NormCt_A"]  # reversed because lower Ct = higher expr
    df = df.set_index("gene")[["log2FC_qPCR"]]
    df = df.groupby(df.index).mean()
    return df


def load_gene_id_to_symbol_map() -> dict[str, str]:
    """Parse GTF to map Ensembl gene_id to gene_name (symbol) quickly."""
    gtf_path = REF_DIR / "GRCh38_ERCC92.gtf"
    mapping = {}
    if not gtf_path.exists():
        print(f"  [WARN] GTF not found at {gtf_path}. Cannot map gene IDs to symbols.")
        return mapping

    gene_id_pat = re.compile(r'gene_id "([^"]+)"')
    gene_name_pat = re.compile(r'gene_name "([^"]+)"')

    with open(gtf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) > 8 and parts[2] in ("gene", "transcript"):
                attrs = parts[8]
                gid_match = gene_id_pat.search(attrs)
                gname_match = gene_name_pat.search(attrs)
                if gid_match:
                    gid = gid_match.group(1)
                    gname = gname_match.group(1) if gname_match else gid
                    mapping[gid] = gname
    return mapping


def compute_aim2_accuracy(tx2gene: pd.Series, gene_map: dict[str, str]) -> pd.DataFrame:
    """
    Compare each tool's log2 fold-change (MAQCA vs MAQCB) against the qPCR truth.
    Metric: Pearson r, Spearman r, RMSE on log2FC over shared genes.
    """
    print("[AIM 2] Computing accuracy metrics against qPCR truth …")

    truth = load_qpcr_truth()
    mats  = build_expression_matrices(AIM2B_SAMPLES, tx2gene)
    rows  = []

    for tool, mat in mats.items():
        if mat.empty:
            print(f"  [WARN] No data for {tool} — skipping.")
            continue

        if "MAQCA" not in mat.columns or "MAQCB" not in mat.columns:
            print(f"  [WARN] MAQCA/MAQCB columns not found in {tool} matrix.")
            continue

        # log2FC from tool (TPM or counts)
        tool_fc = np.log2(mat["MAQCA"] + 1) - np.log2(mat["MAQCB"] + 1)
        tool_fc.name = "log2FC_tool"

        # Convert index from Ensembl gene IDs to Gene Symbols
        if gene_map:
            tool_fc.index = tool_fc.index.map(lambda x: gene_map.get(x, x).upper())
            tool_fc = tool_fc.groupby(tool_fc.index).mean()

        row_base = {"Tool": tool}

        if truth is not None:
            # Align on gene index (gene symbols from qPCR may differ from Ensembl IDs)
            merged = pd.concat([tool_fc, truth["log2FC_qPCR"]], axis=1, join="inner").dropna()
            if merged.shape[0] < 10:
                print(f"  [WARN] Only {merged.shape[0]} overlapping genes for {tool} — check ID mapping.")
                continue

            x = merged["log2FC_tool"].values
            y = merged["log2FC_qPCR"].values

            pearson_r,  pearson_p  = stats.pearsonr(x, y)
            spearman_r, spearman_p = stats.spearmanr(x, y)
            rmse = float(np.sqrt(np.mean((x - y) ** 2)))

            row_base.update({
                "N_Genes_Overlap":  merged.shape[0],
                "Pearson_r":        round(pearson_r,  4),
                "Pearson_p":        f"{pearson_p:.2e}",
                "Spearman_r":       round(spearman_r, 4),
                "Spearman_p":       f"{spearman_p:.2e}",
                "RMSE_log2FC":      round(rmse, 4),
                "Truth_Source":     "Everaert2017_qPCR_GSE83402",
            })
        else:
            # No external truth: at minimum report dynamic range and gene count
            row_base.update({
                "N_Genes":          int(tool_fc.shape[0]),
                "FC_range_log2":    round(float(tool_fc.abs().max()), 4),
                "Truth_Source":     "NOT_AVAILABLE",
            })

        rows.append(row_base)

    df = pd.DataFrame(rows)
    out_path = OUT_DIR / "aim2_accuracy.csv"
    df.to_csv(out_path, index=False)
    print(f"[OK] Aim 2 accuracy → {out_path}  ({len(df)} rows)")
    return df


# ==============================================================================
# SECTION C  —  AIM 3: NOVEL SPLICE-JUNCTION RECOVERY
# ==============================================================================
def load_chr21_truth() -> set[tuple]:
    """
    Load the chr21 junction truth set generated by 01_prepare_data.sh.
    Returns a set of (chr, intron_start, intron_end, strand) tuples.
    Coordinates are 1-based, matching STAR's SJ.out.tab.
    """
    if not CHR21_TRUTH_PATH.exists():
        return set()
    df = pd.read_csv(CHR21_TRUTH_PATH, sep="\t", header=None,
                     names=["chr", "start", "end", "strand"])
    df["chr"] = df["chr"].astype(str)
    return set(zip(df["chr"], df["start"].astype(int),
                   df["end"].astype(int), df["strand"]))


def load_star_sj(sj_path: Path) -> pd.DataFrame:
    """
    Parse STAR SJ.out.tab.
    Columns: chr, start, end, strand_code, motif, annotated, unique, multi, overhang
    strand_code: 1=+, 2=-, 0=undefined.
    Returns DataFrame with normalised (chr, start, end, strand) columns.
    """
    if not sj_path.exists():
        return pd.DataFrame()

    cols = ["chr","start","end","strand_code","motif",
            "annotated","unique_reads","multi_reads","overhang"]
    df = pd.read_csv(sj_path, sep="\t", header=None, names=cols)
    df["chr"] = df["chr"].astype(str)
    df["strand"] = df["strand_code"].map({1: "+", 2: "-", 0: "."})
    return df


def load_hisat2_novel_junctions(jxn_path: Path) -> pd.DataFrame:
    """
    Parse HISAT2 --novel-splicesite-outfile.
    Columns: chr, left_anchor, right_anchor, strand
    Anchors are the last base of the left exon / first base of right exon (0-based).
    We convert to 1-based intron start / end to match STAR SJ.out.tab.
    """
    if not jxn_path.exists():
        return pd.DataFrame()

    df = pd.read_csv(jxn_path, sep="\t", header=None,
                     names=["chr", "left_anchor", "right_anchor", "strand"])
    df["chr"] = df["chr"].astype(str)
    df["start"] = df["left_anchor"] + 2   # convert left exon end → intron start (1-based)
    df["end"]   = df["right_anchor"]       # right_anchor is already intron end (0-based → 1-based: same)
    return df


def precision_recall_f1(found: set, truth: set) -> tuple[float, float, float]:
    tp = len(found & truth)
    fp = len(found - truth)
    fn = len(truth - found)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1        = (2 * precision * recall / (precision + recall)
                 if (precision + recall) > 0 else 0.0)
    return round(precision, 4), round(recall, 4), round(f1, 4)


def compute_aim3_junctions() -> pd.DataFrame:
    """
    Evaluate chr21 novel-junction recovery for:
      STAR 1-pass  (SJ.out.tab, column annotated == 0  → novel)
      STAR 2-pass  (SJ.out.tab, column annotated == 0  → novel)
      HISAT2       (--novel-splicesite-outfile is by definition novel)
    """
    print("[AIM 3] Computing junction recovery metrics …")

    truth_set = load_chr21_truth()
    if not truth_set:
        print(f"  [WARN] Chr21 truth file not found: {CHR21_TRUTH_PATH}. "
              "Run 01_prepare_data.sh first.")

    print(f"  Chr21 truth junctions: {len(truth_set)}")

    rows = []

    def sj_to_set(df: pd.DataFrame, novel_only: bool = False) -> set[tuple]:
        """Convert SJ dataframe to a set of (chr,start,end,strand) tuples."""
        if df.empty:
            return set()
        sub = df[df["annotated"] == 0] if novel_only else df
        return set(zip(sub["chr"].astype(str),
                       sub["start"].astype(int),
                       sub["end"].astype(int),
                       sub["strand"]))

    # ---- STAR 1-pass ----
    sj_1p = load_star_sj(ALN_DIR / "star" / f"{AIM3_SAMPLE}_1pass" / "SJ.out.tab")
    star1_novel = sj_to_set(sj_1p, novel_only=True)
    star1_chr21 = {j for j in star1_novel if j[0] == "21"}
    p, r, f = precision_recall_f1(star1_chr21, truth_set) if truth_set else (0,0,0)
    rows.append({
        "Tool": "STAR_1pass",
        "Total_Junctions":       len(sj_1p),
        "Novel_Junctions":       len(star1_novel),
        "Chr21_Novel_Reported":  len(star1_chr21),
        "Chr21_Truth":           len(truth_set),
        "Precision":             p,
        "Recall":                r,
        "F1":                    f,
    })

    # ---- STAR 2-pass ----
    sj_2p = load_star_sj(ALN_DIR / "star" / f"{AIM3_SAMPLE}_2pass" / "SJ.out.tab")
    star2_novel = sj_to_set(sj_2p, novel_only=True)
    star2_chr21 = {j for j in star2_novel if j[0] == "21"}
    p, r, f = precision_recall_f1(star2_chr21, truth_set) if truth_set else (0,0,0)
    rows.append({
        "Tool": "STAR_2pass",
        "Total_Junctions":       len(sj_2p),
        "Novel_Junctions":       len(star2_novel),
        "Chr21_Novel_Reported":  len(star2_chr21),
        "Chr21_Truth":           len(truth_set),
        "Precision":             p,
        "Recall":                r,
        "F1":                    f,
    })

    # ---- HISAT2 ----
    h2_jxn = load_hisat2_novel_junctions(
        ALN_DIR / "hisat2" / f"{AIM3_SAMPLE}_novel_junctions.txt"
    )
    h2_set = set()
    if not h2_jxn.empty:
        h2_set = set(zip(h2_jxn["chr"].astype(str),
                         h2_jxn["start"].astype(int),
                         h2_jxn["end"].astype(int),
                         h2_jxn["strand"]))
    h2_chr21 = {j for j in h2_set if j[0] == "21"}
    p, r, f = precision_recall_f1(h2_chr21, truth_set) if truth_set else (0,0,0)
    rows.append({
        "Tool": "HISAT2",
        "Total_Junctions":       len(h2_jxn) if not h2_jxn.empty else 0,
        "Novel_Junctions":       len(h2_set),
        "Chr21_Novel_Reported":  len(h2_chr21),
        "Chr21_Truth":           len(truth_set),
        "Precision":             p,
        "Recall":                r,
        "F1":                    f,
    })

    df = pd.DataFrame(rows)
    out_path = OUT_DIR / "aim3_junction_recovery.csv"
    df.to_csv(out_path, index=False)
    print(f"[OK] Aim 3 junctions → {out_path}")
    return df


# ==============================================================================
# SECTION D  —  PERFORMANCE PROFILING
# ==============================================================================
def compute_performance_profile() -> pd.DataFrame:
    """
    Scan all /usr/bin/time -v log files in LOG_DIR and aggregate into a
    performance table: tool, step, sample, wall time (s), peak RAM (MB).
    """
    print("[PERF] Parsing timing logs …")

    # Only parse the primary -v logs (not the .stdout.log tees)
    log_files = sorted(LOG_DIR.glob("*.log"))
    # Exclude stdout mirrors
    log_files = [p for p in log_files if not p.name.endswith(".stdout.log")]

    rows = []
    for log_path in log_files:
        stem = log_path.stem  # e.g. "star_align_SRR1212207"

        tool = "unknown"
        step = "unknown"
        sample = "N/A"

        if stem == "salmon_index":
            tool, step = "Salmon", "index"
        elif stem == "kallisto_index":
            tool, step = "Kallisto", "index"
        elif "index" in stem:
            m = re.match(r"^(star|hisat2)_index_(full|noChr21)$", stem)
            if m:
                tool = m.group(1).capitalize()
                step = "index"
                sample = m.group(2)
        elif "align" in stem:
            if "1pass" in stem:
                tool, step, sample = "Star", "align_1pass", stem.split("_")[-1]
            elif "2pass" in stem:
                tool, step, sample = "Star", "align_2pass", stem.split("_")[-1]
            elif "aim3" in stem:
                tool, step, sample = "Hisat2", "align_aim3", stem.split("_")[-1]
            else:
                m = re.match(r"^(star|hisat2)_align_(SRR\d+)$", stem)
                if m:
                    tool = m.group(1).capitalize()
                    step = "align"
                    sample = m.group(2)
        elif "quant" in stem:
            m = re.match(r"^(salmon|kallisto)_quant_(SRR\d+)$", stem)
            if m:
                tool = m.group(1).capitalize()
                step = "quant"
                sample = m.group(2)
        elif "featureCounts" in stem or "featurecounts" in stem.lower():
            m = re.match(r"^featureCounts_(star|hisat2)_(SRR\d+)$", stem, re.IGNORECASE)
            if m:
                tool = "FeatureCounts_" + m.group(1).capitalize()
                step = "quant"
                sample = m.group(2)

        timing = parse_time_log(log_path)
        rows.append({
            "Log_File":    log_path.name,
            "Tool":        tool,
            "Step":        step,
            "Sample":      sample,
            "Runtime_s":   timing["runtime_s"],
            "Peak_RAM_MB": timing["peak_ram_mb"],
        })

    df = pd.DataFrame(rows)
    # Annotate aims for readability
    aim1_ids = set(AIM1_SAMPLES.keys())
    aim2_ids = set(AIM2A_SAMPLES.keys()) | set(AIM2B_SAMPLES.keys())

    def assign_aim(sample: str) -> str:
        if sample in aim1_ids:       return "Aim1_Reproducibility"
        if sample in aim2_ids:       return "Aim2_Accuracy"
        if AIM3_SAMPLE in sample:    return "Aim3_Junctions"
        if sample in ("N/A", "full", "noChr21"): return "Indexing"
        return "Unknown"

    df["Aim"] = df["Sample"].apply(assign_aim)

    out_path = OUT_DIR / "performance_profile.csv"
    df.to_csv(out_path, index=False)
    print(f"[OK] Performance profile → {out_path}  ({len(df)} log files parsed)")
    return df


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    print("=" * 62)
    print("  RNA-seq Benchmark — Metric Aggregation")
    print("=" * 62)

    # Load shared tx2gene map once
    tx2gene = load_tx2gene()
    if tx2gene.empty:
        print("[WARN] tx2gene not loaded — TPM will not be collapsed to genes.")

    # Section A: Aim 1 Reproducibility
    df_a1 = compute_aim1_reproducibility(tx2gene)
    print(df_a1.to_string(index=False))
    print()

    # Load gene ID to symbol mapping for Aim 2
    gene_map = load_gene_id_to_symbol_map()

    # Section B: Aim 2 Accuracy
    df_a2 = compute_aim2_accuracy(tx2gene, gene_map)
    print(df_a2.to_string(index=False))
    print()

    # Section C: Aim 3 Junction Recovery
    df_a3 = compute_aim3_junctions()
    print(df_a3.to_string(index=False))
    print()

    # Section D: Performance Profile
    df_perf = compute_performance_profile()
    print(df_perf.to_string(index=False))
    print()

    print("=" * 62)
    print(f"  [DONE] All summaries written to: {OUT_DIR}")
    print("=" * 62)


if __name__ == "__main__":
    main()
