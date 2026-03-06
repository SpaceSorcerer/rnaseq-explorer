# RNA-seq Analysis Pipeline

Modular RNA-seq differential expression and alternative splicing analysis pipeline for human embryonic stem cells (hESCs).

## Overview

This pipeline performs comprehensive RNA-seq analysis including:
- DESeq2 differential expression analysis
- rMATS alternative splicing analysis
- Cross-condition concordance analysis
- Publication-quality figure generation
- GSEA and GO enrichment analysis

## Features

- **Modular design**: Separate modules for each analysis step
- **Config-driven**: All parameters in YAML config file
- **User-friendly**: Single CLI entry point with comprehensive options
- **Publication-ready**: 300 DPI figures with color-blind friendly Okabe-Ito palette
- **Flexible**: Works for any RNA-seq experiment, not just MIAT/QKI

## Installation

```bash
# Clone repository
git clone https://github.com/SpaceSorcerer/rnaseq-pipeline.git
cd rnaseq-pipeline

# Install dependencies
pip install pandas numpy matplotlib seaborn openpyxl pyyaml
```

## Quick Start

```bash
# Run complete pipeline with default config
python run_pipeline.py

# Run with custom config
python run_pipeline.py --config my_config.yaml

# Run only specific conditions
python run_pipeline.py --conditions MIAT_OE_vs_Control,QKI_KO_vs_WT

# Run only DESeq2 analysis, skip splicing
python run_pipeline.py --modules deseq2

# Quick filtering check without regenerating figures
python run_pipeline.py --skip-figures
```

## Configuration

Edit `rnaseq_config.yaml` to customize:
- Input file paths
- Filtering thresholds (log2FC, padj, baseMean, dPSI)
- Conditions to analyze
- Figure settings (DPI, format, which figures to generate)
- Analysis modules to run

## Project Structure

```
RNA-seq_MIAT_OE_KCN_QKI-KO/
├── run_pipeline.py          # Single CLI entry point
├── rnaseq_config.yaml       # All parameters
├── modules/
│   ├── deseq2.py            # DESeq2 differential expression
│   ├── splicing.py          # rMATS alternative splicing
│   ├── concordance.py       # Directional overlap analysis
│   ├── figures.py           # Publication figures
│   └── standardize.py       # HGNC symbol mapping, biotype annotation
├── PIPELINE_USAGE.md        # Detailed usage guide
└── CONSOLIDATION_ANALYSIS.md # Technical consolidation notes
```

## Output

The pipeline generates:
- **Excel files**: Filtered DEGs and splicing events with gene symbols and biotypes
- **Figures**: Volcano plots, MA plots, Venn diagrams, concordance scatters, violin plots, heatmaps
- **Summary statistics**: DEG counts, pathway enrichment results
- **Prism files**: CSV files compatible with GraphPad Prism for re-plotting

## Modules

### deseq2.py
- Load and normalize DESeq2 output
- Apply filtering thresholds
- Gene name lookup (MyGene → BioMart → Ensembl fallback)
- Biotype annotation
- Export filtered DEGs to Excel

### splicing.py
- Load rMATS output files (SE, A3SS, A5SS, MXE, RI)
- Filter by dPSI, FDR, and p-value
- Event type summaries
- Directional analysis (included/excluded)
- Export to Excel

### concordance.py
- Cross-condition DEG overlap analysis
- Directional concordance (up/down agreement)
- Scatter plots with Pearson correlation
- Venn diagrams (2-way and 3-way)
- Excel exports with concordance labels

### figures.py
- All publication-quality figure generation
- 300 DPI resolution
- Okabe-Ito color-blind friendly palette
- Customizable layouts

### standardize.py
- Gene ID/name lookup functions
- Biotype normalization
- Column standardization
- Data preprocessing utilities

## Filtering Thresholds

Default thresholds (can be customized in config):
- **DESeq2**: log2FC ≥ 0.4, padj ≤ 0.05, baseMean ≥ 20
- **rMATS**: |dPSI| ≥ 0.1, FDR ≤ 0.05

## Citation

If you use this pipeline, please cite:
- DESeq2: Love, M.I., Huber, W., Anders, S. (2014) [Genome Biology]
- rMATS: Shen, S. et al. (2014) [PNAS]

## Author

Brian Amburn, MD-PhD Candidate  
University of Texas Medical Branch (UTMB)  
Dissertation research: lncRNA MIAT and QKI protein interactions

## License

MIT License - see LICENSE file for details
