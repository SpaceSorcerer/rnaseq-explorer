# Visualization Function Extraction Plan

## Overview
Extract ~30+ visualization functions from the monolith (`deseq2_rmats_filter_pipeline.py`) into the modular engine package. Each function must be parameterized (no global variables) and placed in the correct module.

## Global Variables to Parameterize
All functions in the monolith use these globals that must become function parameters:
- `DESEQ2_COLS` -> `cols: dict[str, str]`
- `RMATS_COLS` -> `rmats_cols: dict[str, str]`
- `LOG2FC_CUTOFF`, `PADJ_CUTOFF`, `BASEMEAN_CUTOFF` -> explicit params or a `cutoffs` dict
- `RMATS_FDR_CUTOFF`, `RMATS_PVAL_CUTOFF`, `INCLEVEL_DIFF_CUTOFF`, `USE_FDR`, `RMATS_DUAL_FILTER`
- `COLOR_UP`, `COLOR_DOWN`, `COLOR_NS` -> import from viz.theme
- `EVENT_COLORS` -> import from viz.theme
- `FIG_FORMAT`, `FIG_DPI` -> `fig_format: str = "png"`, `fig_dpi: int = 300`
- `INTERACTIVE_PLOTS` -> caller decides whether to call interactive functions

## Module Assignments

### 1. deseq2.py — Add Per-Condition DESeq2 Viz (monolith lines 1537-2055)
Functions to add:
- `volcano_plot(df, outdir, cols, padj_cutoff, log2fc_cutoff, basemean_cutoff, label, suffix, fig_format, fig_dpi)`
- `ma_plot(df, outdir, cols, padj_cutoff, log2fc_cutoff, basemean_cutoff, label, suffix, fig_format, fig_dpi)`
- `volcano_plot_interactive(df, outdir, cols, padj_cutoff, log2fc_cutoff, basemean_cutoff, label, suffix)`
- `ma_plot_interactive(df, outdir, cols, padj_cutoff, log2fc_cutoff, basemean_cutoff, label, suffix)`
- `biotype_chart(filtered_df, outdir, cols, label, suffix, fig_format, fig_dpi)`
- `biotype_direction_chart(filtered_df, outdir, cols, label, suffix, fig_format, fig_dpi)`
- `biotype_enrichment_test(filtered_df, all_df, outdir, cols, label, suffix, fig_format, fig_dpi)`
- `biotype_volcano(all_df, outdir, cols, padj_cutoff, log2fc_cutoff, basemean_cutoff, label, suffix, fig_format, fig_dpi)`
- `ecdf_log2fc_by_biotype(all_df, outdir, cols, log2fc_cutoff, basemean_cutoff, label, suffix, fig_format, fig_dpi)`
- `pvalue_histogram(df, outdir, cols, padj_cutoff, label, suffix, fig_format)`
- `top_genes_lollipop(filtered_df, outdir, cols, log2fc_cutoff, label, suffix, fig_format, top_n=20)`
- `expression_rank_plot(df, outdir, cols, padj_cutoff, log2fc_cutoff, basemean_cutoff, label, suffix, fig_format, fig_dpi)`
- `volcano_plot_labeled(df, outdir, cols, padj_cutoff, log2fc_cutoff, basemean_cutoff, label, suffix, fig_format, fig_dpi, genes_of_interest=None)` (monolith line 6412)
- `rbp_heatmap(condition_results, condition_labels, outdir, cols, fig_format, fig_dpi)` (monolith line 1305)
- `rbp_summary_table(condition_results, condition_labels, outdir, cols)` (monolith line 1382)

Also add these helper functions:
- `_assign_biotype_group(series)` (monolith line 1805)
- `_bh_correction(pvals)` (monolith line 1813)

### 2. deseq2.py — Add Cross-Condition Biotype Viz (monolith lines 2057-2197)
- `cross_condition_biotype_comparison(condition_results, condition_labels, outdir, cols, fig_format, fig_dpi)`
- `cross_condition_biotype_direction(condition_results, condition_labels, outdir, cols, fig_format, fig_dpi)`

### 3. rmats.py — Add Per-Condition rMATS Viz (monolith lines 2274-2488)
Functions to add:
- `rmats_scatter(df, event_type, outdir, rmats_cols, fdr_cutoff, pval_cutoff, dpsi_cutoff, use_fdr, dual_filter, fig_format, fig_dpi)`
- `rmats_combined_volcano(all_data, outdir, rmats_cols, fdr_cutoff, pval_cutoff, dpsi_cutoff, use_fdr, fig_format)`
- `rmats_event_summary_chart(filtered_counts, outdir, use_fdr, fdr_cutoff, pval_cutoff, dpsi_cutoff, fig_format)`
- `rmats_dpsi_distribution(all_filtered, outdir, rmats_cols, dpsi_cutoff, fig_format)`
- `rmats_psi_scatter(rmats_raw, rmats_filtered, event_type, outdir, rmats_cols, dpsi_cutoff, fig_format)`

### 4. cross_condition.py — Add All Missing Cross-Condition Functions
Functions to add (grouped by type):

DESeq2 cross-condition:
- `deseq2_de_counts_chart(condition_results, condition_labels, outdir, cols, fig_format)` (line 3658)
- `pairwise_deg_venns(condition_results, condition_labels, outdir, cols, fig_format, fig_dpi)` (line 3411)
- `deseq2_upset_plot(condition_results, condition_labels, outdir, cols, fig_format)` (line 3700)
- `deseq2_venn_diagrams(gene_sets, condition_labels, outdir, fig_format)` (line 2523)
- `deseq2_direction_concordance(condition_results, condition_labels, outdir, cols, fig_format)` (line 2557)
- `deseq2_log2fc_heatmap(condition_results, condition_labels, outdir, cols, fig_format)` (line 2667)
- `pairwise_log2fc_scatter(condition_results, condition_labels, outdir, cols, fig_format)` (line 3756)

rMATS cross-condition:
- `rmats_cross_condition_venn(condition_results, condition_labels, outdir, rmats_cols, fig_format, fig_dpi, match_by)` (line 2753)
- `rmats_direction_concordance(condition_results, condition_labels, outdir, rmats_cols, fig_format)` (line 2849)
- `rmats_directional_venn_diagrams(condition_results, condition_labels, outdir, rmats_cols, dpsi_cutoff, fig_format, fig_dpi, match_by)` (line 2991)
- `pairwise_splicing_venns(condition_results, condition_labels, outdir, rmats_cols, dpsi_cutoff, fig_format, fig_dpi, match_by)` (line 3268)
- `rmats_event_count_comparison(rmats_conditions, condition_labels, outdir, fig_format)` (line 3911)
- `pairwise_dpsi_scatter(rmats_conditions, condition_labels, outdir, rmats_cols, fig_format, match_by)` (line 3948)
- `rmats_upset_plot(rmats_conditions, condition_labels, outdir, rmats_cols, fig_format, match_by)` (line 4099)
- `rmats_event_heatmap(condition_results, condition_labels, event_type, outdir, rmats_cols, dpsi_cutoff, fig_format, fig_dpi)` (line 4211)
- `rmats_event_pie_chart(condition_results, condition_labels, outdir, fig_format, fig_dpi)` (line 4344)

Combined DESeq2+rMATS:
- `deseq2_vs_rmats_venn(condition_results, condition_labels, outdir, cols, rmats_cols, fig_format)` (line 4750)
- `log2fc_vs_dpsi_scatter(condition_results, condition_labels, outdir, cols, rmats_cols, dpsi_cutoff, log2fc_cutoff, fig_format)` (line 4784)

Summary:
- `gene_overlap_summary(condition_results, condition_labels, outdir, cols)` (line 6562)
- `summary_dashboard(condition_results, condition_labels, go_results, gsea_results, outdir, cols, fig_format)` (line 6655)

### 5. exports.py — Add Missing Export Functions
- `export_pairwise_workbook(condition_results, condition_labels, outdir, cols, rmats_cols, dpsi_cutoff)` (line 4421)

### 6. gsea.py — Add GSEA Visualization Functions (monolith lines 5862-6320)
- `_db_short_label(db_name)` (line 5983)
- `_collect_gsea_rows(gsea_results, condition_labels)` (line 5995)
- `gsea_combined_plot(gsea_results, condition_labels, outdir, fig_format, fig_dpi)` (line 6060)
- `gsea_enrichment_plots(gsea_results, condition_labels, outdir, fig_format, fig_dpi)` (line 6202)
- `export_gsea_leading_edge(gsea_results, condition_labels, outdir)` (line 6320)

### 7. ora.py — Add ORA Visualization Functions (monolith lines 5548-5860)
- `go_enrichment_combined_plot(go_results, condition_labels, outdir, fig_format, fig_dpi, filename_suffix)` (line 5548)
- `export_go_prism(go_results, condition_labels, outdir, filename_suffix)` (line 5716)

### 8. exports.py — Add validate_outputs (monolith line 7819)
- `validate_outputs(condition_results, condition_labels, outdir, cols, rmats_cols)` (line 7819)

### 9. pipeline.py — Rewrite to Call All Functions
Rewrite `run_pipeline()` to match the exact call order from `main()` (monolith lines 7991-8413).

## Execution Order
1. deseq2.py viz functions (Task A)
2. rmats.py viz functions (Task B)
3. cross_condition.py all missing functions (Task C)
4. gsea.py viz functions (Task D)
5. ora.py viz functions (Task E)
6. exports.py missing functions (Task F)
7. pipeline.py rewrite (Task G - depends on A-F)
8. Update file_map.md
