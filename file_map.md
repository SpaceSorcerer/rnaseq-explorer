# RNA-seq Explorer -- File Map

## Root Files

- **E:\Claude\rnaseq-explorer\CLAUDE.md** -- Project instructions, architecture overview, and build/run instructions for the RNA-seq Explorer workstation.
- **E:\Claude\rnaseq-explorer\file_map.md** -- This file. Full path and 2-sentence description for every file in the project.
- **E:\Claude\rnaseq-explorer\pyproject.toml** -- Python package configuration with dependencies, optional extras (full, dev), and project metadata.
- **E:\Claude\rnaseq-explorer\LICENSE** -- MIT license, Brian Amburn / University of Texas Medical Branch, 2024-2026.
- **E:\Claude\rnaseq-explorer\requirements.txt** -- Legacy pip requirements file listing all direct dependencies.
- **E:\Claude\rnaseq-explorer\README.md** -- Comprehensive project README covering features, installation, quick start, input formats, page descriptions, batch mode, and development instructions.
- **E:\Claude\rnaseq-explorer\CITATION.cff** -- Machine-readable citation metadata (CFF format) for academic use and GitHub citation integration.
- **E:\Claude\rnaseq-explorer\.gitignore** -- Git ignore rules for Python artifacts, testing caches, data files, figures, and IDE files.
- **E:\Claude\rnaseq-explorer\.github\workflows\test.yml** -- GitHub Actions CI workflow running pytest across Python 3.9-3.12 with ruff linting on 3.12.
- **E:\Claude\rnaseq-explorer\enhancement_plan.md** -- Planning document for pipeline enhancements and feature additions.
- **E:\Claude\rnaseq-explorer\extraction_plan.md** -- Detailed extraction plan listing all ~30+ visualization functions to modularize, their monolith line numbers, target modules, and parameterized signatures. All tasks now complete.
- **E:\Claude\rnaseq-explorer\review_code_audit.md** -- Deep code audit covering import errors, engine completeness vs monolith, pipeline wrapper verification, exception handling, type issues, viz/UI review, and test coverage gaps. Identifies 1 CRITICAL bug, 2 HIGH issues, 3 MEDIUM issues, and 5 LOW issues.

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
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\exports.py** -- Export functionality for Excel workbooks (per-condition, combined, pairwise comparison), comprehensive GraphPad Prism .pzfx files (11 export types: DEG counts, top DEGs, volcano data, splicing data, GSEA pathways with -log10(FDR), SE/other splicing overlaps, DEG overlaps, summary statistics, concordance scatter from raw data, dPSI distributions), PowerPoint slide reports, unfiltered merged data, and output validation (validate_outputs). All globals parameterized.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\engine\pipeline.py** -- Full orchestration wrapper calling all ~50+ functions across all submodules in the exact order matching the monolith's main(). Exposes `run_pipeline(config)` with 34-key DEFAULT_CONFIG. Backward-compatible with batch scripts and Tkinter GUI.

## Visualization Package

- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\__init__.py** -- Visualization subpackage init listing all available viz modules.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\theme.py** -- Unified visualization theme with Okabe-Ito colorblind-safe palette, Plotly/Matplotlib template registration, standard color maps for conditions/events/biotypes, and helper functions (count boxes, Venn styling, grid layout).
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\deseq2_viz.py** -- Interactive Plotly DEG visualizations: volcano_plot, ma_plot, pvalue_distribution, log2fc_distribution, top_genes_bar, biotype_breakdown. All use theme.py colors and handle empty DataFrames gracefully.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\rmats_viz.py** -- Interactive Plotly splicing visualizations: dpsi_volcano (colored by event type), event_type_pie, dpsi_distribution (faceted), top_splicing_events, genes_by_event_count (stacked bar). Uses EVENT_COLORS from theme.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\gsea_viz.py** -- Interactive Plotly enrichment visualizations: nes_bar_chart, enrichment_dot_plot, leading_edge_table, ora_dot_plot, enrichment_comparison. Auto-detects column names from gseapy and Enrichr outputs.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\genewalk_viz.py** -- GeneWalk Plotly visualizations: gw_volcano, gw_gene_bar, gw_network (bipartite gene-GO layout), gw_heatmap, gw_domain_pie, gw_gene_summary. Supports per-gene exploration and network rendering.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\qc_viz.py** -- Interactive QC visualizations: pca_plot (scatter by condition), correlation_heatmap (annotated Pearson), top_deg_heatmap (z-scored expression). All use theme colors and handle missing data gracefully.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\cross_condition_viz.py** -- Cross-condition Plotly visualizations: direction_concordance_heatmap, log2fc_scatter (with fit line and significance coloring), overlap_bar (stacked unique/shared genes). Supports pairwise condition comparison.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\viz\gene_investigator.py** -- Per-gene evidence aggregation: investigate_gene() collects DEG/GSEA/ORA/splicing/GeneWalk data, gene_evidence_card() renders visual summary with figures and text. Searches across multiple column name conventions.

## UI Package

- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\__init__.py** -- UI subpackage init listing all Streamlit modules and pages.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\app.py** -- Main Streamlit entry point. Multi-tab app with sidebar, data loading, dark/light mode, and dynamic tab visibility based on uploaded data. Run with: `streamlit run rnaseq_explorer/ui/app.py`.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\sidebar.py** -- Shared sidebar: file uploads (DESeq2, rMATS, GeneWalk, counts, GSEA, ORA), threshold sliders (log2FC, padj, dPSI, FDR), display settings (top N, dark mode), and export info. Returns settings dict.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\styles.py** -- CSS styling for Streamlit: colored metric card borders (up/down/total/splicing), sidebar styling, tab styling, dark/light mode support, and utility classes (info/warning/success boxes).
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\__init__.py** -- Pages subpackage init.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\overview.py** -- Overview page: summary metrics (total DEGs, up/down counts, splicing events) and quick-look charts (volcano thumbnail, top enrichment).
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\deseq2_page.py** -- DESeq2 page: all DEG visualizations (volcano, MA, p-value dist, log2FC dist, top genes, biotype), genes-of-interest highlighting, filterable data table with CSV export.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\splicing_page.py** -- Splicing page: all rMATS visualizations (dPSI volcano, event pie, distributions, top events, gene counts), event type multi-selector, data table with CSV export.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\enrichment_page.py** -- Enrichment page: tabbed GSEA/ORA/Comparison views, database selector, NES bars, dot plots, leading edge gene table, up-vs-down comparison chart.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\qc_page.py** -- QC page: PCA scatter (with sklearn), sample correlation heatmap, top DEG z-scored heatmap. Only shown when normalized counts are uploaded.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\cross_condition_page.py** -- Cross-condition page: direction concordance heatmap, pairwise log2FC scatter with condition selectors, gene overlap bar chart. Supports upload of additional condition datasets.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\genewalk_page.py** -- GeneWalk page: volcano, domain pie, gene summaries, per-gene GO bar chart with gene selector, network graph, heatmap, data table. Only shown when GeneWalk data uploaded.
- **E:\Claude\rnaseq-explorer\rnaseq_explorer\ui\pages\gene_investigator_page.py** -- Gene Investigator page: gene search with autocomplete, evidence aggregation across all data types, metric cards, visual evidence figures, and expandable detail tables.

## Tests

- **E:\Claude\rnaseq-explorer\tests\__init__.py** -- Test suite init.
- **E:\Claude\rnaseq-explorer\tests\conftest.py** -- Shared pytest fixtures: sample DataFrames for DESeq2 (~30 rows), rMATS (~20 rows), GeneWalk (~15 rows), GSEA (~10 rows), ORA (~10 rows), and an empty DataFrame.
- **E:\Claude\rnaseq-explorer\tests\test_viz_deseq2.py** -- Tests for deseq2_viz.py: volcano_plot, ma_plot, pvalue_distribution, log2fc_distribution, top_genes_bar, biotype_breakdown. Covers normal data, empty data, single gene, custom cutoffs, and genes-of-interest highlighting.
- **E:\Claude\rnaseq-explorer\tests\test_viz_rmats.py** -- Tests for rmats_viz.py: dpsi_volcano, event_type_pie, dpsi_distribution, top_splicing_events, genes_by_event_count. Covers normal data, empty data, single event type, and missing columns.
- **E:\Claude\rnaseq-explorer\tests\test_viz_gsea.py** -- Tests for gsea_viz.py: nes_bar_chart, enrichment_dot_plot, leading_edge_table, ora_dot_plot, enrichment_comparison. Covers normal data, empty data, strict FDR fallback, and gene extraction from leading edges.
- **E:\Claude\rnaseq-explorer\tests\test_viz_genewalk.py** -- Tests for genewalk_viz.py: gw_volcano, gw_gene_bar, gw_network, gw_heatmap, gw_domain_pie, gw_gene_summary. Covers normal data, empty data, few nodes, single gene, and strict filter edge cases.
- **E:\Claude\rnaseq-explorer\tests\test_gene_investigator.py** -- Tests for gene_investigator.py: investigate_gene with all/partial/no sources, gene_evidence_card output types. Covers gene-found, gene-not-found, and empty evidence scenarios.
- **E:\Claude\rnaseq-explorer\tests\test_theme.py** -- Tests for theme.py: setup_plotly_theme execution, PALETTE keys/values, CONDITION_COLORS, EVENT_COLORS, condition_color_map, and grid_dims.

## Review Documents

- **E:\Claude\rnaseq-explorer\review_engine_completeness.md** -- Comprehensive cross-reference verifying all 111 monolith functions are accounted for in the modular engine. Documents renamed, moved, inlined, and intentionally dropped functions with rationale.

## Empty Directories

- **E:\Claude\rnaseq-explorer\configs\** -- Planned: YAML/JSON configuration presets.
- **E:\Claude\rnaseq-explorer\sample_data\** -- Planned: Small sample datasets for testing.
