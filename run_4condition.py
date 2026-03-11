#!/usr/bin/env python3
"""
4-Condition RNA-seq Analysis Launcher
=====================================
Runs the DESeq2/rMATS filtering pipeline across all four conditions:

  1. MIAT OE vs Control       (overexpression)
  2. QKI-KO vs WT             (full knockout)
  3. polyQKI-KO vs WT         (polyalanine-expansion knockout)
  4. MIAT KD vs Non-Targeting  (knockdown, from miat-kd-rnaseq pipeline)

Cutoffs match the P01 Transcriptome Analysis Workflow:
  DESeq2:
    - baseMean > 100
    - padj < 0.01
    - |log2FC| >= 0.4
  rMATS (dual filter -- both required):
    - FDR < 0.05 AND PValue < 0.01
    - |dPSI| >= 0.1

Usage:
    python run_4condition.py            # Run full analysis
    python run_4condition.py --check    # Validate inputs only
"""

import sys
import os
from datetime import datetime
from pathlib import Path

# Add pipeline directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deseq2_rmats_filter_pipeline as pipeline

# =============================================================================
# DIRECTORIES
# =============================================================================

BASE = "/mnt/f/MIAT OE v QKI-KO v polyKQI-KO"
PREV_OUTPUT = f"{BASE}/pipeline_output"

MIAT_KD_BASE = "/mnt/f/miat-kd-rnaseq/phase2_results"

# Always create a fresh output directory with timestamp
_ts = datetime.now().strftime("%Y_%m_%d_%H%M%S")
OUTPUT_DIR = f"/mnt/f/RNA-seq_4Condition_Analysis_{_ts}/"

# =============================================================================
# CONDITIONS
# =============================================================================

CONDITIONS = [
    {
        "name": "MIAT_OE_vs_Control",
        "label": "MIAT OE vs Control",
        "deseq2_file": f"{PREV_OUTPUT}/MIAT_OE_vs_Control_all_genes.xlsx",
        "rmats_dir": f"{BASE}/rMATs_MIAT.vs.Control/",
        # "vasttools_file": "",  # Future: VAST-tools output for this condition
    },
    {
        "name": "QKI_KO_vs_WT",
        "label": "QKI-KO vs WT",
        "deseq2_file": f"{PREV_OUTPUT}/QKI_KO_vs_WT_all_genes.xlsx",
        "rmats_dir": f"{BASE}/rMATS QKIKOvsWT/",
        # "vasttools_file": "",  # Future: VAST-tools output for this condition
    },
    {
        "name": "polyQKI_KO_vs_WT",
        "label": "polyQKI-KO vs WT",
        "deseq2_file": f"{PREV_OUTPUT}/polyQKI_KO_vs_WT_all_genes.xlsx",
        "rmats_dir": f"{BASE}/rMATS_NKX25-GFP-UD_WT-polyQKIKO_RNAseq/",
        # "vasttools_file": "",  # Future: VAST-tools output for this condition
    },
    {
        "name": "MIAT_KD_vs_NT",
        "label": "MIAT KD vs Non-Targeting",
        "deseq2_file": f"{MIAT_KD_BASE}/deseq2/WT_MIAT_KD_vs_WT_NT.results.tsv",
        "rmats_dir": f"{MIAT_KD_BASE}/rmats/WT_MIAT_KD_vs_WT_NT/",
        # "vasttools_file": "",  # Future: VAST-tools output for this condition
    },
]

# rMATS event types expected in each rmats_dir
RMATS_EVENT_TYPES = ["SE", "A3SS", "A5SS", "RI", "MXE"]
RMATS_FILE_SUFFIX = ".MATS.JCEC.txt"

# =============================================================================
# PRE-PROCESSING: Normalize MIAT KD column names
# =============================================================================


def preprocess_miat_kd(conditions, output_dir):
    """Normalize DESeq2 column names for the MIAT KD condition.

    Standard DESeq2 output uses camelCase (log2FoldChange, baseMean).
    The first 3 conditions' Excel files already have lowercase columns
    (log2fc, basemean) from prior processing. This function converts
    the MIAT KD TSV to match.

    Returns the (possibly modified) conditions list.
    """
    miat_kd = None
    miat_kd_idx = None
    for i, cond in enumerate(conditions):
        if cond["name"] == "MIAT_KD_vs_NT":
            miat_kd = cond
            miat_kd_idx = i
            break

    if miat_kd is None:
        return conditions

    src = miat_kd["deseq2_file"]
    if not os.path.isfile(src):
        print(f"[WARN] MIAT KD DESeq2 file not found: {src}")
        print("       Phase 2 pipeline may not be complete yet.")
        print("       Removing MIAT KD condition from this run.\n")
        conditions.pop(miat_kd_idx)
        return conditions

    # Read the TSV
    import pandas as pd

    print(f"[INFO] Pre-processing MIAT KD DESeq2 file: {src}")
    df = pd.read_csv(src, sep="\t")

    # Rename standard DESeq2 columns to lowercase convention
    rename_map = {
        "log2FoldChange": "log2fc",
        "baseMean": "basemean",
    }
    # Only rename columns that actually exist
    actual_renames = {k: v for k, v in rename_map.items() if k in df.columns}
    if actual_renames:
        df.rename(columns=actual_renames, inplace=True)
        print(f"  Renamed columns: {actual_renames}")
    else:
        print("  Column names already normalized (or different naming scheme).")

    # If there's a gene ID index column (standard DESeq2 uses row names as gene IDs),
    # check if it needs to be promoted to a column
    if df.columns[0] not in ("gene_id", "GeneID", "Geneid"):
        # Check if the index looks like Ensembl IDs
        first_val = str(df.iloc[0, 0]) if len(df) > 0 else ""
        if first_val.startswith("ENSG"):
            print("  First column appears to contain Ensembl gene IDs.")

    # Save normalized file to the output directory
    os.makedirs(output_dir, exist_ok=True)
    norm_path = os.path.join(output_dir, "MIAT_KD_vs_NT_normalized.tsv")
    df.to_csv(norm_path, sep="\t", index=False)
    print(f"  Saved normalized file: {norm_path}\n")

    # Update the condition to point to the normalized file
    conditions[miat_kd_idx]["deseq2_file"] = norm_path
    return conditions


# =============================================================================
# INPUT VALIDATION (--check flag)
# =============================================================================


def check_inputs(conditions):
    """Validate that all input files and directories exist.

    Returns True if all checks pass, False otherwise.
    """
    print("=" * 80)
    print("  INPUT VALIDATION")
    print("=" * 80)

    all_ok = True
    for cond in conditions:
        print(f"\n  [{cond['name']}] {cond['label']}")

        # Check DESeq2 file
        deseq2_path = cond["deseq2_file"]
        if os.path.isfile(deseq2_path):
            size_mb = os.path.getsize(deseq2_path) / (1024 * 1024)
            print(f"    DESeq2 file: OK ({size_mb:.1f} MB)")
        else:
            print(f"    DESeq2 file: MISSING - {deseq2_path}")
            all_ok = False

        # Check rMATS directory
        rmats_dir = cond["rmats_dir"]
        if os.path.isdir(rmats_dir):
            print(f"    rMATS dir:   OK")
            # Check for the 5 event type files
            for event in RMATS_EVENT_TYPES:
                fname = f"{event}{RMATS_FILE_SUFFIX}"
                fpath = os.path.join(rmats_dir, fname)
                if os.path.isfile(fpath):
                    print(f"      {fname}: OK")
                else:
                    print(f"      {fname}: MISSING")
                    all_ok = False
        else:
            print(f"    rMATS dir:   MISSING - {rmats_dir}")
            all_ok = False

    print()
    if all_ok:
        print("  All inputs validated successfully.")
    else:
        print("  WARNING: Some inputs are missing. Pipeline may fail or skip conditions.")
    print("=" * 80)
    return all_ok


# =============================================================================
# PIPELINE CONFIGURATION
# =============================================================================

config = {
    "CONDITIONS": CONDITIONS,
    "OUTPUT_DIR": OUTPUT_DIR,
    # P01 workflow cutoffs
    "LOG2FC_CUTOFF": 0.4,
    "BASEMEAN_CUTOFF": 100,
    "PADJ_CUTOFF": 0.01,
    # rMATS dual filter: BOTH FDR AND PValue required
    "RMATS_FDR_CUTOFF": 0.05,
    "RMATS_PVAL_CUTOFF": 0.01,
    "INCLEVEL_DIFF_CUTOFF": 0.1,
    "USE_FDR": True,
    "RMATS_DUAL_FILTER": True,
    # Figure settings
    "FIG_DPI": 300,
    "FIG_FORMAT": "png",
    "FONT_SIZE": 12,
    "COLOR_UP": "#E69F00",
    "COLOR_DOWN": "#0072B2",
    "COLOR_NS": "#BFBFBF",
    "INTERACTIVE_PLOTS": True,
    "AUTO_BIOTYPE_SPLIT": True,
    "GENE_NAME_LOOKUP": True,
    "SPECIES": "human",
    # Dual ORA mode
    "ORA_METHOD": "both",
    # GSEA settings
    "GSEA_RANKING": "stat",
    "GSEA_MIN_SIZE": 10,
    "GSEA_MAX_SIZE": 500,
    "GSEA_PERMUTATIONS": 1000,
    # RBP annotation
    "RBP_FILE": "/mnt/f/RBP Lists/RBP-E-A-C-Complex.xlsx",
    # Genes of interest — highlighted in volcano plots and tables
    "GENES_OF_INTEREST": ["MIAT", "QKI", "QKI-5", "QKI-6", "QKI-7"],
    # Column mappings (pre-processed files have normalized names)
    "DESEQ2_COLS": {
        "gene_id": "gene_id",
        "gene_name": "gene_name",
        "log2fc": "log2fc",
        "basemean": "basemean",
        "padj": "padj",
        "pvalue": "padj",
        "biotype": "biotype_group",
        "stat": "stat",
        "lfcSE": "lfcSE",
    },
    "RMATS_COLS": {
        "event_id": "ID",
        "gene_id": "GeneID",
        "gene_name": "geneSymbol",
        "pvalue": "PValue",
        "fdr": "FDR",
        "inclevel_diff": "IncLevelDifference",
    },
}

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # --check mode: validate inputs and exit
    if "--check" in sys.argv:
        ok = check_inputs(CONDITIONS)
        sys.exit(0 if ok else 1)

    # Pre-process MIAT KD columns (gracefully skips if file doesn't exist yet)
    config["CONDITIONS"] = preprocess_miat_kd(config["CONDITIONS"], OUTPUT_DIR)

    n_conditions = len(config["CONDITIONS"])
    if n_conditions == 0:
        print("ERROR: No valid conditions remaining after pre-processing.")
        sys.exit(1)

    print("=" * 80)
    print("  4-CONDITION RNA-seq ANALYSIS — PIPELINE RUN")
    print("=" * 80)
    print(f"\nOutput: {config['OUTPUT_DIR']}")
    print(f"DESeq2: padj < {config['PADJ_CUTOFF']}, |log2FC| >= {config['LOG2FC_CUTOFF']}, "
          f"baseMean >= {config['BASEMEAN_CUTOFF']}")
    print(f"rMATS:  FDR < {config['RMATS_FDR_CUTOFF']} AND PValue < {config['RMATS_PVAL_CUTOFF']}, "
          f"|dPSI| >= {config['INCLEVEL_DIFF_CUTOFF']} [DUAL FILTER]")
    print(f"ORA: {config['ORA_METHOD']} (Enrichr + g:Profiler)")
    print(f"GSEA ranking: {config['GSEA_RANKING']} (min_size={config['GSEA_MIN_SIZE']})")
    print(f"RBP annotation: {config['RBP_FILE']}")
    print(f"Genes of interest: {', '.join(config['GENES_OF_INTEREST'])}")
    print(f"Conditions: {n_conditions}")
    for cond in config["CONDITIONS"]:
        print(f"  - {cond['label']}")
    print()

    try:
        pipeline.run_pipeline(config)
        print("\n" + "=" * 80)
        print("  PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        sys.exit(0)
    except Exception as e:
        import traceback
        print("\n" + "=" * 80)
        print(f"  PIPELINE FAILED: {e}")
        traceback.print_exc()
        print("=" * 80)
        sys.exit(1)
