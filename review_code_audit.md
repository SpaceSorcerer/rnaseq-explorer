# RNA-seq Explorer -- Deep Code Audit

**Date:** 2026-03-11
**Scope:** All Python files under `rnaseq_explorer/` and `tests/`
**Reference:** Monolith `deseq2_rmats_filter_pipeline.py` (8,497 lines, 62 functions)

---

## 1. Import Errors

### 1.1 CATEGORY_COLORS Shadowing in ora.py
- **File:** `rnaseq_explorer/engine/ora.py`, lines 28 and 50
- **Severity:** MEDIUM
- **Description:** Line 28 imports `CATEGORY_COLORS` from `rnaseq_explorer.viz.theme`, then line 50 redefines a local `CATEGORY_COLORS` dict with the same name. The import on line 28 is effectively dead code -- the local redefinition shadows it entirely. The local values happen to match the theme values today, but this creates a maintenance hazard: changing `CATEGORY_COLORS` in `theme.py` will not propagate to `ora.py`.
- **Fix:** Remove line 28 import (`from rnaseq_explorer.viz.theme import CATEGORY_COLORS`) and keep the local definition, OR remove lines 50-56 and use the theme import exclusively.

### 1.2 parse_gsea_results Defined but Never Called
- **File:** `rnaseq_explorer/engine/gsea.py`, line 365
- **Severity:** LOW
- **Description:** `parse_gsea_results()` is defined in `gsea.py` but never imported or called anywhere in the codebase (not in `pipeline.py`, not in any UI file, not in tests). It appears to be a utility function that was written but never wired into any workflow.
- **Fix:** Either integrate it into the pipeline where appropriate, export it from `__init__.py` for external use, or remove it to reduce dead code.

### 1.3 All Cross-Module Imports Are Valid
- No broken imports found. All `from rnaseq_explorer.engine.X import Y` and `from rnaseq_explorer.viz.theme import Z` resolve correctly. Verified: `deseq2.py`, `rmats.py`, `gsea.py`, `ora.py`, `qc.py`, `cross_condition.py`, `exports.py`, `pipeline.py`, all viz modules, and all UI modules.

---

## 2. Engine Completeness (Monolith vs Modular)

### 2.1 Function Coverage Summary

All 62 monolith functions have been accounted for in the modular engine. The mapping is:

| Monolith Function | Modular Location | Notes |
|---|---|---|
| `setup_style()` | `viz/theme.py` -> `setup_matplotlib_style()` | Renamed |
| `add_count_box()` | `viz/theme.py` -> `add_count_box()` | Same name |
| `load_file()` | `engine/deseq2.py` | Same name |
| `load_counts_matrix()` | `engine/qc.py` | Same name |
| `pca_plot()` | `engine/qc.py` -> `compute_pca()` | Renamed |
| `sample_correlation_heatmap()` | `engine/qc.py` -> `compute_sample_correlation()` | Renamed |
| `top_deg_heatmap()` | `engine/qc.py` -> `compute_top_deg_heatmap()` | Renamed |
| `_fetch_gene_names()` | `engine/deseq2.py` -> `fetch_gene_names()` | Made public |
| `_enrich_with_gene_names()` | `engine/deseq2.py` -> `enrich_with_gene_names()` | Made public |
| `_reassign_biotypes_from_mygene()` | `engine/deseq2.py` -> `reassign_biotypes_from_mygene()` | Made public |
| `validate_columns()` | `engine/deseq2.py` | Same name |
| `_resolve_column()` | `engine/deseq2.py` | Same name (kept private) |
| `_validate_rmats_columns()` | `engine/rmats.py` | Same name |
| `_normalize_gsea_cols()` | `engine/gsea.py` -> `normalize_gsea_cols()` | Made public |
| `_strip_ensembl_version()` | `engine/deseq2.py` | Same name |
| `normalize_deseq2_columns()` | `engine/deseq2.py` | Same name |
| `_best_gene_key()` | `engine/deseq2.py` -> `best_gene_key()` | Made public |
| `load_rbp_annotations()` | `engine/deseq2.py` | Same name |
| `annotate_rbps()` | `engine/deseq2.py` | Same name |
| `rbp_heatmap()` | `engine/deseq2.py` | Same name |
| `rbp_summary_table()` | `engine/deseq2.py` | Same name |
| `normalize_rmats_columns()` | `engine/rmats.py` | Same name |
| `filter_deseq2()` | `engine/deseq2.py` | Same name |
| `volcano_plot()` | `engine/deseq2.py` | Same name |
| `ma_plot()` | `engine/deseq2.py` | Same name |
| `volcano_plot_interactive()` | `engine/deseq2.py` | Same name |
| `ma_plot_interactive()` | `engine/deseq2.py` | Same name |
| `biotype_chart()` | `engine/deseq2.py` | Same name |
| `_assign_biotype_group()` | `engine/deseq2.py` -> `assign_biotype_group()` | Made public |
| `_bh_correction()` | `engine/deseq2.py` | Same name |
| `biotype_direction_chart()` | `engine/deseq2.py` | Same name |
| `biotype_enrichment_test()` | `engine/deseq2.py` | Same name |
| `biotype_volcano()` | `engine/deseq2.py` | Same name |
| `ecdf_log2fc_by_biotype()` | `engine/deseq2.py` | Same name |
| `cross_condition_biotype_comparison()` | `engine/deseq2.py` | Same name |
| `cross_condition_biotype_direction()` | `engine/deseq2.py` | Same name |
| `load_all_rmats()` | `engine/rmats.py` | Same name |
| `filter_rmats()` | `engine/rmats.py` | Same name |
| `rmats_scatter()` | `engine/rmats.py` | Same name |
| `rmats_combined_volcano()` | `engine/rmats.py` | Same name |
| `rmats_event_summary_chart()` | `engine/rmats.py` | Same name |
| `rmats_dpsi_distribution()` | `engine/rmats.py` | Same name |
| `extract_gene_sets()` | `engine/deseq2.py` | Same name |
| `deseq2_venn_diagrams()` | `engine/cross_condition.py` | Same name |
| `deseq2_direction_concordance()` | `engine/cross_condition.py` -> `compute_direction_heatmap()` | Renamed, refactored |
| `deseq2_log2fc_heatmap()` | `engine/cross_condition.py` | Same name |
| `rmats_cross_condition_venn()` | `engine/cross_condition.py` | Same name |
| `rmats_direction_concordance()` | `engine/cross_condition.py` | Same name |
| `_make_event_key()` | `engine/rmats.py` -> `make_event_key()` | Made public |
| `_style_venn()` | `viz/theme.py` -> `style_venn()` | Moved to theme |
| `rmats_directional_venn_diagrams()` | **NOT EXTRACTED** | Intentional (line 729 of pipeline.py) |
| `pairwise_splicing_venns()` | `engine/cross_condition.py` | Same name |
| `pairwise_deg_venns()` | `engine/cross_condition.py` | Same name |
| `_grid_dims()` | `viz/theme.py` -> `grid_dims()` | Made public, moved |
| `_parse_inclevel_mean()` | `engine/rmats.py` -> `parse_inclevel_mean()` | Made public |
| `pvalue_histogram()` | `engine/deseq2.py` | Same name |
| `top_genes_lollipop()` | `engine/deseq2.py` | Same name |
| `expression_rank_plot()` | `engine/deseq2.py` | Same name |
| `deseq2_de_counts_chart()` | `engine/cross_condition.py` | Same name |
| `deseq2_upset_plot()` | `engine/cross_condition.py` | Same name |
| `pairwise_log2fc_scatter()` | `engine/cross_condition.py` | Same name |
| `rmats_psi_scatter()` | `engine/rmats.py` | Same name |
| `rmats_event_count_comparison()` | `engine/cross_condition.py` | Same name |
| `pairwise_dpsi_scatter()` | `engine/cross_condition.py` | Same name |
| `rmats_upset_plot()` | `engine/cross_condition.py` | Same name |
| `rmats_event_heatmap()` | `engine/cross_condition.py` | Same name |
| `rmats_event_pie_chart()` | `engine/cross_condition.py` | Same name |
| `export_pairwise_workbook()` | `engine/exports.py` | Same name |
| `deseq2_vs_rmats_venn()` | `engine/cross_condition.py` | Same name |
| `log2fc_vs_dpsi_scatter()` | `engine/cross_condition.py` | Same name |
| `export_results()` | `engine/exports.py` -> `export_excel()` | Renamed |
| `build_multi_condition_summary()` | Inlined in `export_combined_results()` | Merged |
| `export_combined_results()` | `engine/exports.py` | Same name |
| `run_gsea_enrichment()` | `engine/gsea.py` | Same name |
| `run_go_enrichment()` | `engine/ora.py` -> `run_enrichr_ora()` + `run_gprofiler_ora()` + `run_dual_ora()` | Split into 3 |
| `go_enrichment_combined_plot()` | `engine/ora.py` | Same name |
| `export_go_prism()` | `engine/ora.py` | Same name |
| `_gsea_dotplot_legacy()` | `engine/gsea.py` -> `gsea_dotplot_legacy()` | Made public |
| `_db_short_label()` | `engine/gsea.py` | Same name |
| `_collect_gsea_rows()` | `engine/gsea.py` | Same name |
| `gsea_combined_plot()` | `engine/gsea.py` | Same name |
| `gsea_enrichment_plots()` | `engine/gsea.py` | Same name |
| `export_gsea_leading_edge()` | `engine/gsea.py` | Same name |
| `volcano_plot_labeled()` | `engine/deseq2.py` | Same name |
| `gene_overlap_summary()` | `engine/cross_condition.py` | Same name |
| `summary_dashboard()` | `engine/cross_condition.py` | Same name |
| `export_prism_files()` | `engine/exports.py` -> `export_prism_pzfx()` | Renamed |
| `export_unfiltered_merged()` | `engine/exports.py` | Same name |
| `generate_powerpoint_report()` | `engine/exports.py` -> `export_powerpoint()` | Renamed |
| `validate_outputs()` | `engine/exports.py` | Same name |
| `main()` | `engine/pipeline.py` -> `run_pipeline()` | Entry point |
| `run_pipeline()` | `engine/pipeline.py` | Same name |

### 2.2 Intentionally Omitted Function
- **File:** `engine/pipeline.py`, line 729
- **Severity:** LOW (documented)
- **Description:** `rmats_directional_venn_diagrams()` (~277 lines in monolith, lines 2991-3267) was intentionally not extracted. The pipeline.py comment reads: "rmats_directional_venn_diagrams was intentionally not extracted". The functionality is partially covered by `pairwise_splicing_venns()` with directional filtering.
- **Fix:** None required. Consider documenting this in CLAUDE.md if not already noted.

### 2.3 Extra Function in Modular Engine
- **File:** `engine/gsea.py`, line 365
- **Severity:** LOW
- **Description:** `parse_gsea_results()` exists in the modular `gsea.py` but has no counterpart in the monolith and is never called. It appears to be a new utility that was added but never integrated.
- **Fix:** Either wire it into a workflow or remove it.

---

## 3. Pipeline Wrapper Verification

### 3.1 run_pipeline(config) vs main() Call Order

The `run_pipeline()` in `engine/pipeline.py` (lines 201-931) was verified against the monolith's `main()` (lines 7991-8418). The call order matches:

1. Style setup
2. Counts matrix loading + PCA + correlation (if counts file provided)
3. RBP annotations loading (if RBP file provided)
4. Per-condition DESeq2 loading, normalization, gene name enrichment, biotype reassignment, RBP annotation
5. Per-condition biotype passes (All, Protein Coding, Non-Protein Coding) with filtering and visualization
6. Per-condition rMATS loading, filtering, visualization
7. Per-condition Excel export
8. Cross-condition DESeq2 comparisons (counts chart, Venn, pairwise Venn, UpSet, concordance, heatmap, scatter, biotype comparison)
9. Cross-condition rMATS comparisons (Venn, event count, UpSet, concordance, pairwise Venn, scatter, heatmap, pie)
10. Combined DESeq2+rMATS (Venn, scatter)
11. Combined export
12. Unfiltered merged export
13. GSEA enrichment
14. ORA (Enrichr/g:Profiler/both)
15. GSEA visualizations (combined plot, enrichment plots, leading edge)
16. Gene overlap summary
17. Summary dashboard
18. Prism export
19. PowerPoint generation
20. Validation

**Verdict:** The pipeline wrapper faithfully reproduces the monolith's `main()` execution order.

### 3.2 DEFAULT_CONFIG Coverage
- **Severity:** N/A
- **Description:** The 34-key `DEFAULT_CONFIG` (lines 133-182) covers all parameters that were previously hardcoded globals in the monolith. Verified: all config keys used in `run_pipeline()` have defaults.

---

## 4. Bare Excepts

**No bare `except:` statements found anywhere in the codebase.** All exception handlers use `except Exception` (with or without variable binding). This is acceptable Python practice.

### 4.1 Silent Exception Swallowing
- **File:** `rnaseq_explorer/engine/gsea.py`, lines 501, 580, 858, 910, 912
- **Severity:** MEDIUM
- **Description:** Five `except Exception:` blocks in `gsea.py` silently `continue` or `pass` without logging. These occur in `_collect_gsea_rows()` (line 501), `gsea_dotplot_legacy()` (line 580), and `gsea_enrichment_plots()` (lines 858, 910, 912). Failed CSV reads, file copies, and gseaplot regenerations are silently skipped, making debugging very difficult.
- **Fix:** Add `logging.debug()` or at minimum `print()` in each catch block to surface failures. Example:
  ```python
  except Exception as e:
      print(f"  [WARN] Skipping {report_csv}: {e}")
      continue
  ```

- **File:** `rnaseq_explorer/engine/deseq2.py`, lines 1316, 1710
- **Severity:** LOW
- **Description:** Two silent `except Exception:` blocks. Line 1316 is in `biotype_enrichment_test()` (catches fisher_exact failure) and line 1710 is in `volcano_plot_labeled()` (catches label placement failure). Both are acceptable for non-critical display operations.

- **File:** `rnaseq_explorer/engine/exports.py`, line 342
- **Severity:** LOW
- **Description:** Silent catch in `_save_prism()` for XML prettification failure. Falls through to a direct write, which is a reasonable fallback.

---

## 5. Type Issues

### 5.1 CRITICAL: Operator Precedence Bug in overview.py
- **File:** `rnaseq_explorer/ui/pages/overview.py`, lines 68-71
- **Severity:** CRITICAL
- **Description:** The expression is:
  ```python
  n_splice = int(
      (rmats_df[fdr_col] < settings["fdr_cutoff"])
      & (rmats_df[dpsi_col].abs() >= settings["dpsi_cutoff"])
  ).sum() if fdr_col else 0
  ```
  The `int()` call wraps the entire boolean Series, which converts it to a single integer (raising `TypeError: cannot convert the series to <class 'int'>` if the Series has more than one element, or returning 0/1 for single-element Series). The `.sum()` is called on the result of `int()`, not on the Series.

  **What was intended:**
  ```python
  n_splice = int(
      (
          (rmats_df[fdr_col] < settings["fdr_cutoff"])
          & (rmats_df[dpsi_col].abs() >= settings["dpsi_cutoff"])
      ).sum()
  ) if fdr_col else 0
  ```
- **Fix:** Move `.sum()` inside the `int()` parentheses so it operates on the boolean Series first.

### 5.2 CONDITION_COLORS Used as List with Index Access
- **File:** `rnaseq_explorer/engine/qc.py` (via cross_condition.py and pipeline.py)
- **Severity:** LOW
- **Description:** `CONDITION_COLORS` is defined as a list in `theme.py`. When conditions exceed the list length, `condition_color_map()` in `theme.py` handles wrapping with modular indexing. However, direct index access (e.g., `CONDITION_COLORS[i]`) without bounds checking could raise `IndexError` if code bypasses `condition_color_map()`.
- **Fix:** All current code uses `condition_color_map()` or direct list comprehension with modular indexing. No immediate fix needed, but document that `CONDITION_COLORS` should always be accessed through `condition_color_map()`.

### 5.3 Optional gprofiler/gseapy Import Handling
- **File:** `rnaseq_explorer/engine/ora.py`, `rnaseq_explorer/engine/gsea.py`
- **Severity:** LOW
- **Description:** Both `gseapy` and `gprofiler-official` are optional dependencies. The engine modules import them at the top level without try/except guards. If these packages are not installed, importing the engine modules will fail with `ModuleNotFoundError`. The viz modules (for Streamlit UI) don't have this problem since they don't import these.
- **Fix:** Consider wrapping the gseapy/gprofiler imports in try/except blocks and raising informative errors only when the relevant functions are actually called, similar to how `deseq2.py` handles `plotly.express` (lines 38-41).

---

## 6. Viz Module Issues

### 6.1 Empty DataFrame Handling
- **Severity:** N/A (PASS)
- **Description:** All 8 viz modules were verified for empty DataFrame handling:
  - `deseq2_viz.py`: All 6 functions return empty `go.Figure()` for empty/None input
  - `rmats_viz.py`: All 5 functions return empty `go.Figure()` for empty/None input
  - `gsea_viz.py`: All 5 functions return empty `go.Figure()` or empty `pd.DataFrame` for empty/None input
  - `genewalk_viz.py`: All 6 functions return empty `go.Figure()` for empty/None input
  - `qc_viz.py`: All 3 functions return empty `go.Figure()` for empty/None input
  - `cross_condition_viz.py`: All 3 functions return empty `go.Figure()` for empty/None input
  - `gene_investigator.py`: Both functions handle None/empty inputs gracefully

### 6.2 Theme Usage Consistency
- **Severity:** N/A (PASS)
- **Description:** All viz modules import and use colors from `theme.py`. Verified: `COLOR_UP`, `COLOR_DOWN`, `COLOR_NS`, `EVENT_COLORS`, `CONDITION_COLORS`, `BIOTYPE_COLORS`, `CATEGORY_COLORS` are used consistently. The Okabe-Ito palette is applied through `setup_plotly_theme()`.

---

## 7. UI Issues

### 7.1 session_state Key Safety
- **Severity:** N/A (PASS)
- **Description:** All UI pages use `st.session_state.get(key)` with fallback to `None`, rather than direct `st.session_state[key]` access. The keys used are: `deseq2_data`, `rmats_data`, `genewalk_data`, `counts_data`, `gsea_data`, `ora_data`, `condition_datasets`. All are initialized in `app.py` via `_load_data()` (line 69) which sets missing keys to `None`.

### 7.2 Viz Import Consistency in UI Pages
- **Severity:** N/A (PASS)
- **Description:** All UI pages import from `rnaseq_explorer.viz.*` modules correctly. Each page imports only the viz functions it needs:
  - `deseq2_page.py` -> `deseq2_viz`
  - `splicing_page.py` -> `rmats_viz`
  - `enrichment_page.py` -> `gsea_viz`
  - `genewalk_page.py` -> `genewalk_viz`
  - `qc_page.py` -> `qc_viz`
  - `cross_condition_page.py` -> `cross_condition_viz`
  - `gene_investigator_page.py` -> `gene_investigator`
  - `overview.py` -> `deseq2_viz`, `gsea_viz`

### 7.3 overview.py Bug (also listed as 5.1)
- **File:** `rnaseq_explorer/ui/pages/overview.py`, lines 68-71
- **Severity:** CRITICAL
- **Description:** See item 5.1 above. This will crash at runtime when rMATS data has more than one row.

---

## 8. Test Coverage

### 8.1 Existing Tests (PASS)

| Test File | Module Tested | Functions Tested | Status |
|---|---|---|---|
| `test_viz_deseq2.py` | `viz/deseq2_viz.py` | 6/6 | Full coverage |
| `test_viz_rmats.py` | `viz/rmats_viz.py` | 5/5 | Full coverage |
| `test_viz_gsea.py` | `viz/gsea_viz.py` | 5/5 | Full coverage |
| `test_viz_genewalk.py` | `viz/genewalk_viz.py` | 6/6 | Full coverage |
| `test_gene_investigator.py` | `viz/gene_investigator.py` | 2/2 | Full coverage |
| `test_theme.py` | `viz/theme.py` | 7 test classes | Full coverage |

### 8.2 Missing Tests -- Engine Modules

- **Severity:** HIGH
- **Description:** No engine modules have any tests at all. This is the largest gap in the test suite.

| Module | Functions | Priority |
|---|---|---|
| `engine/deseq2.py` | 31 functions (17 viz + 14 data) | HIGH -- core data pipeline |
| `engine/rmats.py` | 9 functions (5 viz + 4 data) | HIGH -- core data pipeline |
| `engine/gsea.py` | 11 functions | HIGH -- enrichment analysis |
| `engine/ora.py` | 6 functions | HIGH -- enrichment analysis |
| `engine/qc.py` | 4 functions | MEDIUM |
| `engine/cross_condition.py` | 22 functions | MEDIUM |
| `engine/exports.py` | 7 public + 3 private functions | MEDIUM |
| `engine/pipeline.py` | 2 functions | LOW (integration test) |

**Recommended test priority:**
1. `engine/deseq2.py` data functions: `load_file`, `normalize_deseq2_columns`, `filter_deseq2`, `extract_gene_sets`, `best_gene_key`
2. `engine/rmats.py` data functions: `load_all_rmats`, `filter_rmats`, `make_event_key`
3. `engine/gsea.py` core functions: `normalize_gsea_cols`, `create_ranked_list`
4. `engine/ora.py` core functions: `run_enrichr_ora`, `run_gprofiler_ora`
5. `engine/exports.py`: `export_excel`, `export_prism_pzfx`

### 8.3 Missing Tests -- Viz Modules

- **File:** `viz/qc_viz.py` -- 3 functions, 0 tests
- **File:** `viz/cross_condition_viz.py` -- 3 functions, 0 tests
- **Severity:** MEDIUM
- **Description:** These two viz modules have no dedicated test files. The other 5 viz modules all have thorough test suites.
- **Fix:** Create `tests/test_viz_qc.py` and `tests/test_viz_cross_condition.py` following the same pattern as the existing viz test files (test normal data, empty data, edge cases).

### 8.4 Missing Tests -- UI Modules

- **Severity:** LOW
- **Description:** No UI modules have tests. This is common for Streamlit apps since they require special testing infrastructure (e.g., `streamlit.testing`). The UI pages are thin wrappers around viz functions that are already tested.
- **Fix:** Consider adding basic smoke tests using `streamlit.testing.v1.AppTest` for at least `app.py` and `overview.py` to catch the CRITICAL bug (item 5.1).

---

## Summary of All Issues by Severity

### CRITICAL (1)
1. **overview.py:68-71** -- Operator precedence bug: `int()` wraps boolean Series instead of `.sum()` result. Will crash at runtime with multi-row rMATS data.

### HIGH (2)
1. **No engine tests** -- 8 engine modules with ~90+ functions have zero test coverage.
2. **No qc_viz/cross_condition_viz tests** -- 6 viz functions missing tests.

### MEDIUM (3)
1. **ora.py:28,50** -- `CATEGORY_COLORS` imported then shadowed by local redefinition.
2. **gsea.py:501,580,858,910,912** -- Silent exception swallowing in 5 locations makes debugging difficult.
3. **ora.py/gsea.py** -- `gseapy` and `gprofiler-official` imports at top level will crash if packages not installed, even when those functions aren't needed.

### LOW (5)
1. **gsea.py:365** -- `parse_gsea_results()` defined but never called (dead code).
2. **pipeline.py:729** -- `rmats_directional_venn_diagrams` intentionally not extracted (documented).
3. **deseq2.py:1316,1710** -- Silent catches in non-critical display operations (acceptable).
4. **exports.py:342** -- Silent catch in Prism XML prettification (has fallback).
5. **No UI tests** -- Expected for Streamlit; UI is thin wrapper around tested viz functions.
