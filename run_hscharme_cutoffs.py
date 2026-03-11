#!/usr/bin/env python3
"""
Pipeline re-run with HSCHARME-aligned cutoffs + dual ORA (Enrichr + g:Profiler).
Reference: Buonaiuto et al. 2025, Nature Comms 16:7880

Cutoff changes from previous run:
  - padj: 0.01 → 0.05
  - |log2FC|: 0.4 → 1.0
  - baseMean: 20 → 10
  - rMATS FDR: 0.01 → 0.05
  - ORA: gprofiler → both (Enrichr + g:Profiler side-by-side)
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
OUTPUT_DIR = f"/mnt/f/RNA-seq_Analysis_{_ts}_HSCHARME/"

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
    # HSCHARME-aligned cutoffs
    "LOG2FC_CUTOFF": 1.0,
    "BASEMEAN_CUTOFF": 10,
    "PADJ_CUTOFF": 0.05,
    "RMATS_FDR_CUTOFF": 0.05,
    "RMATS_PVAL_CUTOFF": 0.05,
    "INCLEVEL_DIFF_CUTOFF": 0.1,
    "USE_FDR": True,
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
print("  HSCHARME-ALIGNED PIPELINE RUN")
print("  Buonaiuto et al. 2025, Nature Comms 16:7880")
print("=" * 80)
print(f"\nOutput: {config['OUTPUT_DIR']}")
print(f"Cutoffs: padj < {config['PADJ_CUTOFF']}, |log2FC| >= {config['LOG2FC_CUTOFF']}, "
      f"baseMean >= {config['BASEMEAN_CUTOFF']}")
print(f"rMATS: FDR < {config['RMATS_FDR_CUTOFF']}, |dPSI| >= {config['INCLEVEL_DIFF_CUTOFF']}")
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
