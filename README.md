# RNA-seq DESeq2 + rMATS Analysis Pipeline

A comprehensive pipeline for differential expression (DESeq2) and alternative splicing (rMATS) analysis with publication-quality figures, GO/GSEA enrichment, and GraphPad Prism export.

## Architecture

Two files:

| File | Purpose |
|------|---------|
| `pipeline_launcher.py` | Tkinter GUI — configure and run without editing code |
| `deseq2_rmats_filter_pipeline.py` | Analysis engine — all filtering, figures, enrichment, export |

## Installation

```bash
git clone git@github.com:SpaceSorcerer/rnaseq-pipeline.git
cd rnaseq-pipeline
pip install -r requirements.txt
```

Requires Python 3.10+.

## Usage

### GUI (recommended)

```bash
python pipeline_launcher.py
```

The launcher provides tabbed configuration for:
- **Conditions**: Add DESeq2 files and rMATS directories per condition
- **Thresholds**: padj, log2FC, baseMean, FDR, dPSI cutoffs
- **Output**: Directory, figure format/DPI, colors
- **Enrichment**: GSEA/ORA database selection, species, genes of interest
- **Columns**: Custom column name mapping for non-standard input files

### Programmatic

```python
import deseq2_rmats_filter_pipeline as pipeline

config = {
    "CONDITIONS": [
        {"name": "Treatment_vs_Control", "label": "Treatment vs Control",
         "deseq2_file": "/path/to/deseq2.xlsx", "rmats_dir": "/path/to/rmats/"},
    ],
    "OUTPUT_DIR": "/path/to/output",
    "LOG2FC_CUTOFF": 0.4,
    "BASEMEAN_CUTOFF": 20,
    "PADJ_CUTOFF": 0.01,
    "RMATS_FDR_CUTOFF": 0.01,
    "INCLEVEL_DIFF_CUTOFF": 0.1,
    "USE_FDR": True,
}

pipeline.run_pipeline(config)
```

All config keys are optional except `CONDITIONS` and `OUTPUT_DIR`. Defaults are applied for any missing keys.

## Outputs

Per condition:
- Volcano plots (static, labeled, interactive HTML)
- MA plots, p-value histograms, expression rank plots
- Biotype distribution and enrichment analysis
- rMATS scatter plots per event type (SE, A3SS, A5SS, RI, MXE)
- GO ORA combined dot plots (up/down, all databases)
- GSEA combined multi-database plots + enrichment score plots
- Excel exports (DESeq2 + rMATS filtered results)

Cross-condition:
- Venn diagrams (3-way + pairwise, DEG + splicing)
- UpSet plots, direction concordance heatmaps
- log2FC heatmap, pairwise scatter plots
- DE vs AS overlap analysis
- Gene overlap summary table
- Summary statistics dashboard

Exports:
- GraphPad Prism .pzfx files (18+ per run)
- PowerPoint report
- Leading edge gene tables (GSEA)

## Dependencies

See `requirements.txt`. Key packages: pandas, matplotlib, seaborn, gseapy, openpyxl, scipy, adjustText, upsetplot, mygene, plotly, python-pptx, matplotlib-venn.

## Author

Brian Amburn, MD-PhD Candidate
University of Texas Medical Branch (UTMB)
