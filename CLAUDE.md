# RNA-seq Explorer

## Description
Comprehensive interactive RNA-seq downstream analysis workstation. Processes DESeq2 differential expression and rMATS alternative splicing results with GSEA/ORA enrichment analysis, publication-quality visualization, and multi-format export (Excel, Prism, PowerPoint).

## Architecture
- **Modular engine** (`rnaseq_explorer/engine/`) -- core analysis logic split into focused modules
- **Plotly/Matplotlib viz** (`rnaseq_explorer/viz/`) -- unified colorblind-safe theming
- **Streamlit UI** (`rnaseq_explorer/ui/`) -- planned interactive web interface (future)
- **Legacy monolith** (`deseq2_rmats_filter_pipeline.py`) -- original 8,497-line engine, kept as reference
- **Tkinter launcher** (`pipeline_launcher.py`) -- original GUI launcher

## Status
Being modularized from a monolithic pipeline. The original monolith is preserved and still functional. The new modular engine provides the same functionality through `rnaseq_explorer.engine.pipeline.run_pipeline(config)`.

## Key Files

### Modular Engine (new)
| File | Purpose |
|------|---------|
| `rnaseq_explorer/engine/deseq2.py` | DESeq2 loading, filtering, gene name lookup, biotype assignment, RBP annotation |
| `rnaseq_explorer/engine/rmats.py` | rMATS event parsing, filtering (single/dual mode), all 5 event types |
| `rnaseq_explorer/engine/gsea.py` | GSEA prerank runner with gseapy, ranked list creation, result normalization |
| `rnaseq_explorer/engine/ora.py` | ORA via Enrichr (gseapy) and g:Profiler, dual-method support |
| `rnaseq_explorer/engine/qc.py` | PCA, sample correlation heatmap, top DEG heatmap (requires counts matrix) |
| `rnaseq_explorer/engine/cross_condition.py` | Venn diagrams, UpSet plots, concordance heatmaps, pairwise scatter |
| `rnaseq_explorer/engine/exports.py` | Excel, Prism .pzfx, PowerPoint, unfiltered merged exports |
| `rnaseq_explorer/engine/pipeline.py` | Backward-compatible `run_pipeline(config)` orchestrator |
| `rnaseq_explorer/viz/theme.py` | Okabe-Ito palette, Plotly/Matplotlib theming, helper functions |

### Original Files (reference)
| File | Purpose |
|------|---------|
| `deseq2_rmats_filter_pipeline.py` | Original 8,497-line monolithic engine (do not delete) |
| `pipeline_launcher.py` | Tkinter GUI launcher |
| `new_functions.py` | Additional functions for event-level matching and heatmaps |
| `run_4condition.py` | Batch script for 4-condition analysis |
| `run_hscharme_cutoffs.py` | Batch script with HSCHARME paper cutoffs |
| `run_p01_cutoffs.py` | Batch script with P01 grant cutoffs |

## Build/Run Instructions

### Install
```bash
pip install -e ".[full]"
```

### Run (modular engine)
```python
from rnaseq_explorer.engine.pipeline import run_pipeline

config = {
    "CONDITIONS": [
        {"name": "Treatment_vs_Control", "label": "Treatment vs Control",
         "deseq2_file": "/path/to/deseq2.xlsx", "rmats_dir": "/path/to/rmats/"},
    ],
    "OUTPUT_DIR": "/path/to/output",
    "LOG2FC_CUTOFF": 0.4,
    "PADJ_CUTOFF": 0.01,
}
run_pipeline(config)
```

### Run (legacy monolith, still works)
```python
import deseq2_rmats_filter_pipeline as pipeline
pipeline.run_pipeline(config)
```

### Run (Tkinter GUI)
```bash
python pipeline_launcher.py
```

## Configuration
All config keys are optional except CONDITIONS and OUTPUT_DIR. See `rnaseq_explorer/engine/pipeline.py` DEFAULT_CONFIG for the full list with defaults.
