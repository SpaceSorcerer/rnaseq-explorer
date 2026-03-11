# RNA-seq Explorer

> Interactive downstream analysis workstation for RNA-seq data.

RNA-seq Explorer is a comprehensive, modular toolkit for exploring differential
expression (DESeq2) and alternative splicing (rMATS) results, running pathway
enrichment, and generating publication-quality figures -- all from a single
Streamlit web interface or programmatic Python API.

---

## Features

- **DESeq2 Analysis** -- Volcano plots, MA plots, p-value distributions,
  log2FC histograms, top gene bar charts, biotype breakdowns, and gene
  filtering with auto-detected column names.
- **Alternative Splicing** -- rMATS event visualization for all 5 event types
  (SE, A3SS, A5SS, RI, MXE), dPSI distributions, event type summaries, and
  gene-level event counts.
- **Pathway Enrichment** -- GSEA prerank analysis across 7+ databases
  (Hallmark, KEGG, Reactome, GO, WikiPathways, etc.) and ORA via g:Profiler
  and Enrichr with dual-method support.
- **QC Dashboards** -- PCA scatter plots, sample-sample Pearson correlation
  heatmaps, and top-DEG z-scored expression heatmaps.
- **Cross-Condition Comparison** -- Direction concordance heatmaps, pairwise
  log2FC scatter plots, gene overlap bar charts, Venn diagrams, and UpSet
  plots for multi-condition experiments.
- **GeneWalk Integration** -- Import GeneWalk functional annotation results;
  visualize with volcano plots, gene-GO networks, similarity heatmaps, domain
  pie charts, and per-gene GO term summaries.
- **Gene Investigator** -- Search any gene and see all evidence across DEG,
  GSEA, ORA, splicing, and GeneWalk in a single evidence card.
- **Publication Export** -- ggplot2-quality static figures via Matplotlib
  (300 DPI), formatted Excel workbooks, GraphPad Prism .pzfx files, and
  PowerPoint slide reports.
- **Colorblind-Safe Theming** -- All plots use the Okabe-Ito palette with
  dark/light mode support.

---

## Installation

### From source (recommended)

```bash
git clone https://github.com/SpaceSorcerer/rnaseq-explorer.git
cd rnaseq-explorer
pip install -e ".[full,dev]"
```

### Minimal install

```bash
pip install -e .
```

### Docker (coming soon)

```bash
docker build -t rnaseq-explorer .
docker run -p 8501:8501 rnaseq-explorer
```

Requires **Python 3.9+**.

---

## Quick Start

### Streamlit UI (interactive)

```bash
streamlit run rnaseq_explorer/ui/app.py
```

Then open `http://localhost:8501` in your browser. Upload your files via the
sidebar and explore your data interactively.

### Programmatic (batch mode)

```python
from rnaseq_explorer.engine.pipeline import run_pipeline

config = {
    "CONDITIONS": [
        {
            "name": "Treatment_vs_Control",
            "label": "Treatment vs Control",
            "deseq2_file": "/path/to/deseq2_results.xlsx",
            "rmats_dir": "/path/to/rmats_output/",
        },
    ],
    "OUTPUT_DIR": "/path/to/output",
    "LOG2FC_CUTOFF": 0.4,
    "PADJ_CUTOFF": 0.01,
}

run_pipeline(config)
```

### Legacy monolith (still works)

```python
import deseq2_rmats_filter_pipeline as pipeline
pipeline.run_pipeline(config)
```

### Tkinter GUI

```bash
python pipeline_launcher.py
```

---

## Input Files

| Data Type | Format | Required |
|-----------|--------|----------|
| DESeq2 results | `.xlsx` or `.csv` with columns: gene name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj | Yes (at least one condition) |
| rMATS output | Directory containing `SE.MATS.JC.txt`, `A3SS.MATS.JC.txt`, etc. | Optional |
| Normalized counts | `.csv` or `.xlsx` matrix (genes x samples) | Optional (for QC) |
| GeneWalk results | `.csv` with hgnc_symbol, go_name, go_id, go_domain, sim, gene_padj | Optional |
| GSEA results | Pre-computed `.csv` with Term, NES, FDR columns | Optional (or run live) |
| ORA results | Pre-computed Enrichr/g:Profiler `.csv` | Optional (or run live) |

Column names are auto-detected. Custom mappings can be provided via the config
or the sidebar.

---

## Pages

| Page | Description |
|------|-------------|
| **Overview** | Summary metrics (total DEGs, up/down counts, splicing events) and quick-look charts. |
| **DESeq2** | All DEG visualizations: volcano, MA, p-value dist, log2FC dist, top genes, biotype breakdown. Supports genes-of-interest highlighting and filterable data table with CSV export. |
| **Splicing** | rMATS visualizations: dPSI volcano (colored by event type), event pie chart, dPSI distributions, top events, gene-by-event-count. Multi-select event type filter. |
| **Enrichment** | Tabbed GSEA/ORA/Comparison views: NES bar charts, enrichment dot plots, leading edge gene tables, ORA combined dot plots, up-vs-down comparison. |
| **QC** | PCA scatter, sample correlation heatmap, top DEG z-scored expression heatmap. Requires normalized counts matrix. |
| **Cross-Condition** | Direction concordance heatmap, pairwise log2FC scatter with condition selectors, gene overlap bar chart. Supports additional condition uploads. |
| **GeneWalk** | Volcano, domain pie, gene summary rankings, per-gene GO bar chart, bipartite gene-GO network, similarity heatmap. |
| **Gene Investigator** | Search any gene; see aggregated evidence from every analysis type with metric cards and visual figures. |

---

## Batch Mode

The engine modules can be used programmatically in existing batch scripts
without the Streamlit UI:

```python
from rnaseq_explorer.engine.deseq2 import load_deseq2, filter_significant
from rnaseq_explorer.engine.rmats import load_rmats_results
from rnaseq_explorer.engine.gsea import run_gsea_prerank
from rnaseq_explorer.viz.deseq2_viz import volcano_plot

# Load and filter
df = load_deseq2("results.xlsx")
sig = filter_significant(df, log2fc_cutoff=1.0, padj_cutoff=0.05)

# Visualize
fig = volcano_plot(sig, genes_of_interest=["MIAT", "QKI"])
fig.write_html("volcano.html")
```

---

## Development

```bash
# Install with dev dependencies
pip install -e ".[full,dev]"

# Run tests
pytest tests/ -v

# Lint
ruff check .
```

---

## Project Structure

```
rnaseq-explorer/
  rnaseq_explorer/
    engine/          # Core analysis logic (DESeq2, rMATS, GSEA, ORA, QC, exports)
    viz/             # Plotly/Matplotlib visualization functions + theme
    ui/              # Streamlit web interface (app, sidebar, pages)
  tests/             # pytest test suite
  configs/           # Configuration presets (planned)
  sample_data/       # Sample datasets for testing (planned)
```

---

## Citation

If you use RNA-seq Explorer in your research, please cite:

```
Amburn, B. (2024-2026). RNA-seq Explorer: Interactive downstream analysis
workstation for RNA-seq data. University of Texas Medical Branch.
https://github.com/SpaceSorcerer/rnaseq-explorer
```

See `CITATION.cff` for machine-readable citation metadata.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

## Author

**Brian Amburn**, MD-PhD Candidate
University of Texas Medical Branch (UTMB)
Biochemistry & Molecular Biology
