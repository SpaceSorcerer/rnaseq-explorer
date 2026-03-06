# RNA-seq Pipeline Usage Guide

## Overview

`run_pipeline.py` is the single entry point for all RNA-seq analysis. It replaces the previous 6+ run scripts with one unified command-line interface.

## Quick Start

### Run the complete pipeline (default)
```bash
python run_pipeline.py
```

This will:
- Load DESeq2 and rMATS data for all conditions in `rnaseq_config.yaml`
- Apply filtering thresholds from config
- Generate all figures (volcano, MA, biotype, splicing, Venns, concordance, etc.)
- Export Excel files with filtered results
- Run cross-condition concordance analysis
- Generate summary statistics

### Use a different config file
```bash
python run_pipeline.py --config my_custom_config.yaml
```

### Override output directory
```bash
python run_pipeline.py --output-dir results_2026_03_06/
```

## Common Use Cases

### Analyze only specific conditions
```bash
# Run only MIAT OE and QKI-KO
python run_pipeline.py --conditions MIAT_OE_vs_Control,QKI_KO_vs_WT
```

### Run only specific modules
```bash
# DESeq2 analysis only (skip splicing)
python run_pipeline.py --modules deseq2

# DESeq2 and splicing only (skip concordance)
python run_pipeline.py --modules deseq2,splicing

# Only concordance (assumes data already processed)
python run_pipeline.py --modules concordance
```

### Skip figure generation
```bash
# Useful for quick filtering checks without regenerating plots
python run_pipeline.py --skip-figures
```

## Output Structure

After running, you'll find:

```
pipeline_output/
├── pipeline.log                    # Detailed execution log
├── summary_statistics.csv          # Summary table of all conditions
├── excel/                          # Filtered results
│   ├── MIAT_OE_vs_Control_DEGs.xlsx
│   ├── MIAT_OE_vs_Control_splicing.xlsx
│   ├── QKI_KO_vs_WT_DEGs.xlsx
│   └── ...
├── figures/                        # All publication-quality figures
│   ├── MIAT_OE_vs_Control_volcano.png
│   ├── MIAT_OE_vs_Control_ma.png
│   ├── MIAT_OE_vs_Control_biotype.png
│   ├── MIAT_OE_vs_Control_splicing_summary.png
│   ├── log2fc_violin.png
│   ├── top_genes_heatmap.png
│   └── ...
└── concordance/                    # Cross-condition comparisons
    ├── MIAT_vs_QKI_deg_venn.png
    ├── MIAT_vs_QKI_deg_overlap.xlsx
    ├── MIAT_vs_QKI_splicing_SE_venn.png
    └── ...
```

## Configuration

All parameters are defined in `rnaseq_config.yaml`:

- **Filtering thresholds**: log2FC, padj, baseMean, dPSI, FDR
- **Conditions**: Input file paths and labels
- **Comparisons**: Which condition pairs to compare
- **Figures**: Which plots to generate, DPI, format
- **Modules**: Which analysis modules to run

## Modules

The pipeline uses modular components from `modules/`:

- **deseq2**: Load, filter, and summarize DESeq2 differential expression data
- **splicing**: Load, filter, and parse rMATS alternative splicing events
- **figures**: Generate all publication-quality plots
- **concordance**: Cross-condition overlap and concordance analysis
- **standardize**: Gene name lookup and biotype normalization

## Logging

The pipeline writes detailed logs to `pipeline_output/pipeline.log`. Check this file if you encounter errors or want to see detailed progress.

Console output shows:
- High-level progress messages
- Counts of significant genes/events per condition
- Summary statistics table
- Final output location

## Examples

### Daily analysis workflow
```bash
# Edit rnaseq_config.yaml with new thresholds
# Then run:
python run_pipeline.py --output-dir results_$(date +%Y%m%d)/
```

### Quick check of new filtering thresholds
```bash
python run_pipeline.py --skip-figures
# Check summary_statistics.csv to see counts
```

### Regenerate only figures (data already filtered)
```bash
# Not directly supported - but you can:
# 1. Comment out data loading in run_pipeline.py
# 2. Or use individual figure scripts
```

### Debug a single condition
```bash
python run_pipeline.py --conditions MIAT_OE_vs_Control --modules deseq2
```

## Troubleshooting

### Pipeline fails with "Module not found"
Make sure you're in the correct directory:
```bash
cd /workspace/extra/drive-f/RNA-seq_MIAT_OE_KCN_QKI-KO/
python run_pipeline.py
```

### "File not found" errors
Check that:
- `base_dir` in config points to the correct data directory
- DESeq2 Excel files exist at the specified paths
- rMATS directories exist and contain JCEC files

### Figures fail to generate
Check the log file for specific errors. Common issues:
- Missing required columns in input data
- Empty dataframes after filtering (thresholds too strict)
- Font or matplotlib backend issues

### Memory issues
If you run out of memory:
- Process one condition at a time: `--conditions MIAT_OE_vs_Control`
- Skip figure generation: `--skip-figures`
- Increase Docker memory allocation

## Migration from Old Scripts

This pipeline replaces:
- `deseq2_rmats_filter_pipeline.py` → Now handled by modules
- `run_analysis.py` → Use `run_pipeline.py` instead
- `run_full_pipeline.py` → Use `run_pipeline.py` instead
- Individual plotting scripts → Integrated into `modules/figures.py`
- Manual Excel exports → Automated with `--export`

All functionality is preserved, just unified under one command.
