# RNA-seq Explorer -- File Map

## Root Files

- **E:\Claude\rnaseq-explorer\CLAUDE.md** -- Project instructions, architecture overview, and build/run instructions for the RNA-seq Explorer workstation.
- **E:\Claude\rnaseq-explorer\file_map.md** -- This file. Full path and 2-sentence description for every file in the project.
- **E:\Claude\rnaseq-explorer\pyproject.toml** -- Python package configuration with dependencies, optional extras (full, dev), and project metadata.
- **E:\Claude\rnaseq-explorer\LICENSE** -- MIT license, Brian Amburn / University of Texas Medical Branch, 2024-2026.
- **E:\Claude\rnaseq-explorer\requirements.txt** -- Legacy pip requirements file listing all direct dependencies.
- **E:\Claude\rnaseq-explorer\README.md** -- Original project README describing the pipeline architecture, usage, and outputs.
- **E:\Claude\rnaseq-explorer\.gitignore** -- Git ignore rules for Python artifacts, data files, figures, and IDE files.
- **E:\Claude\rnaseq-explorer\enhancement_plan.md** -- Planning document for pipeline enhancements and feature additions.
- **E:\Claude\rnaseq-explorer\extraction_plan.md** -- Detailed extraction plan listing all ~30+ visualization functions to modularize, their monolith line numbers, target modules, and parameterized signatures. All tasks now complete.

## Original Monolith (reference, do not delete)

- **E:\Claude\rnaseq-explorer\deseq2_rmats_filter_pipeline.py** -- Original 8,497-line monolithic analysis engine with all DESeq2, rMATS, GSEA, ORA, visualization, and export logic. Kept as reference during modularization.
- **E:\Claude\rnaseq-explorer\pipeline_launcher.py** -- Tkinter GUI launcher that configures and runs the pipeline without editing code. Provides tabbed configuration for conditions, thresholds, output, enrichment, and column mappings.
- **E:\Claude\rnaseq-explorer\new_functions.py** -- Additional functions for event-level coordinate matching, directional Venn diagrams, clustered heatmaps, event pie charts, and pairwise workbook export.
- **E:\Claude\rnaseq-explorer\run_4condition.py** -- Batch script running the 4-condition analysis (MIAT OE, QKI-KO, polyQKI-KO, MIAT KD) with P01 cutoffs.
- **E:\Claude\rnaseq-explorer\run_hscharme_cutoffs.py** -- Batch script running analysis with HSCHARME paper-standard cutoffs (padj<0.05, |log2FC|>1, baseMean>10).
- **E:\Claude\rnaseq-explorer\run_p01_cutoffs.py** -- Batch script running analysis with P01 grant cutoffs (padj<0.01, |log2FC|>=0.4, baseMean>100).

## Modular Engine Package

- **E:\Claude\rnaseq-explorer\rnaseq_explorer\__init__.py** -- Package init with version and author metadata.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\__init__.py** -- Engine subpackage init listing all available modules.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\deseq2.py** -- DESeq2 loading, column auto-detection, gene name lookup via MyGene.info, biotype assignment, RBP annotation, significance filtering, and 17 visualization functions (volcano, MA, biotype charts, p-value histogram, expression rank, lollipop, interactive plots, cross-condition biotype comparison, RBP heatmap/summary). All globals parameterized.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\rmats.py** -- rMATS event parsing for all 5 types (SE, A3SS, A5SS, RI, MXE), column validation, single/dual-filter mode, event key generation, plus 5 visualization functions (scatter, combined volcano, event summary chart, dPSI distribution, PSI scatter). All globals parameterized.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\gsea.py** -- GSEA prerank enrichment runner wrapping gseapy with ranked list creation, result normalization, per-condition analysis, plus 4 visualization/export functions (gsea_combined_plot, gsea_enrichment_plots, export_gsea_leading_edge, gsea_dotplot_legacy) and helpers (_db_short_label, _collect_gsea_rows). All globals parameterized.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\ora.py** -- Over-representation analysis with g:Profiler and Enrichr backends, dual-mode runner, plus 2 visualization/export functions (go_enrichment_combined_plot, export_go_prism). All globals parameterized; uses CATEGORY_COLORS from viz.theme.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\qc.py** -- QC analyses requiring a normalized counts matrix: PCA scatter plot, sample-sample Pearson correlation heatmap, and top-50-DEG z-scored expression heatmap. Supports both pre-computed PCA files and from-counts computation.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\cross_condition.py** -- Multi-condition comparison analyses with ~22 functions: DESeq2 Venn/UpSet/concordance/heatmap/scatter, rMATS cross-condition Venn/UpSet/concordance/scatter/heatmap/pie, combined DESeq2+rMATS Venn and scatter, gene overlap summary, and summary dashboard. Works for both DESeq2 gene sets and rMATS splicing events. All globals parameterized.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\exports.py** -- Export functionality for Excel workbooks (per-condition, combined, pairwise comparison), GraphPad Prism .pzfx files, PowerPoint slide reports, unfiltered merged data, and output validation (validate_outputs). All globals parameterized.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\pipeline.py** -- Full orchestration wrapper calling all ~50+ functions across all submodules in the exact order matching the monolith's main(). Exposes `run_pipeline(config)` with 34-key DEFAULT_CONFIG. Backward-compatible with batch scripts and Tkinter GUI.

## Visualization Package

- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\__init__.py** -- Visualization subpackage init.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\theme.py** -- Unified visualization theme with Okabe-Ito colorblind-safe palette, Plotly/Matplotlib template registration, standard color maps for conditions/events/biotypes, and helper functions (count boxes, Venn styling, grid layout).

## UI Package (planned)

- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\__init__.py** -- UI subpackage init for future Streamlit pages.

## Tests

- **E:\Claude\rnaseq-explorer\tests\__init__.py** -- Test suite init.

## Empty Directories

- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\** -- Planned: Streamlit page modules.
- **E:\Claude\rnaseq-explorer\configs\** -- Planned: YAML/JSON configuration presets.
- **E:\Claude\rnaseq-explorer\sample_data\** -- Planned: Small sample datasets for testing.
