#!/usr/bin/env python3
"""
Pipeline run with P01 Transcriptome Analysis Workflow cutoffs.

Cutoffs (from P01 workflow PDF, with L2FC adjusted to |0.4|):
  DESeq2:
    - baseMean > 100
    - padj < 0.01
    - |log2FC| >= 0.4
  rMATS (dual filter — both required):
    - FDR < 0.05 AND PValue < 0.01
    - |dPSI| >= 0.1
"""

import sys
import os
from datetime import datetime

# Add pipeline directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deseq2_rmats_filter_pipeline as pipeline

BASE = "/mnt/f/MIAT OE v QKI-KO v polyKQI-KO"
PREV_OUTPUT = f"{BASE}/pipeline_output"

# Always create a fresh output directory with timestamp
_ts = datetime.now().strftime("%Y_%m_%d_%H%M%S")
OUTPUT_DIR = f"/mnt/f/RNA-seq_Analysis_{_ts}_P01/"

config = {
    "CONDITIONS": [
        {
            "name": "MIAT_OE_vs_Control",
            "label": "MIAT OE vs Control",
            "deseq2_file": f"{PREV_OUTPUT}/MIAT_OE_vs_Control_all_genes.xlsx",
            "rmats_dir": f"{BASE}/rMATs_MIAT.vs.Control/",
        },
        {
            "name": "QKI_KO_vs_WT",
            "label": "QKI-KO vs WT",
            "deseq2_file": f"{PREV_OUTPUT}/QKI_KO_vs_WT_all_genes.xlsx",
            "rmats_dir": f"{BASE}/rMATS QKIKOvsWT/",
        },
        {
            "name": "polyQKI_KO_vs_WT",
            "label": "polyQKI-KO vs WT",
            "deseq2_file": f"{PREV_OUTPUT}/polyQKI_KO_vs_WT_all_genes.xlsx",
            "rmats_dir": f"{BASE}/rMATS_NKX25-GFP-UD_WT-polyQKIKO_RNAseq/",
        },
    ],
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
    "GSEA_MIN_SIZE": 15,
    "GSEA_MAX_SIZE": 500,
    "GSEA_PERMUTATIONS": 1000,
    # RBP annotation
    "RBP_FILE": "/mnt/f/RBP Lists/RBP-E-A-C-Complex.xlsx",
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

print("=" * 80)
print("  P01 TRANSCRIPTOME ANALYSIS WORKFLOW — PIPELINE RUN")
print("=" * 80)
print(f"\nOutput: {config['OUTPUT_DIR']}")
print(f"DESeq2: padj < {config['PADJ_CUTOFF']}, |log2FC| >= {config['LOG2FC_CUTOFF']}, "
      f"baseMean >= {config['BASEMEAN_CUTOFF']}")
print(f"rMATS:  FDR < {config['RMATS_FDR_CUTOFF']} AND PValue < {config['RMATS_PVAL_CUTOFF']}, "
      f"|dPSI| >= {config['INCLEVEL_DIFF_CUTOFF']} [DUAL FILTER]")
print(f"ORA: {config['ORA_METHOD']} (Enrichr + g:Profiler)")
print(f"GSEA ranking: {config['GSEA_RANKING']}")
print(f"Conditions: {len(config['CONDITIONS'])}\n")

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
