# Plan: WSF-Lab-Pipeline Compatibility + New Analysis Features

## Context

Brian compared his RNA-seq pipeline (`/mnt/f/rnaseq-pipeline/`) with the WSF-Lab-Pipeline (`/mnt/f/WSF-Lab-Pipeline/`, Nextflow-based). WSF has features Brian's lacks (PCA, heatmaps, better GSEA methodology). Brian wants WSF as ground truth for analytical methods.

**Files to modify:**
- `/mnt/f/rnaseq-pipeline/deseq2_rmats_filter_pipeline.py` (~6044 lines)
- `/mnt/f/rnaseq-pipeline/pipeline_launcher.py` (~714 lines)
- `/mnt/f/rnaseq-pipeline/requirements.txt`

**Reviewed by 5 subagents** (statistical rigor, backwards compatibility, code architecture, WSF compatibility, UX) **+ 3 NanoClaw bots** (QA, Science Review, RNA-seq Pipeline Dev). Changes below incorporate all feedback.

---

## Step 1: New Config Constants

**Insert at ~line 193** (after `GENES_OF_INTEREST`):

```python
# --- Normalized Counts (optional, for PCA/heatmaps) ---
COUNTS_FILE = ""              # Path to normalized_counts.tsv (genes x samples)
SAMPLE_METADATA = {}          # {"sample_name": "condition_name", ...} or auto-detect

# --- GSEA Settings (WSF ground truth) ---
GSEA_RANKING = "stat"         # "stat" (Wald, WSF ground truth) or "log2fc"
                              # Falls back to log2fc if stat column missing
GSEA_MIN_SIZE = 15            # WSF uses 15; more rigorous than 5
GSEA_MAX_SIZE = 500
GSEA_PERMUTATIONS = 1000      # WSF uses 1000 (was 100); better p-value accuracy

# --- ORA Method ---
ORA_METHOD = "gprofiler"      # "gprofiler" (WSF ground truth, g:SCS correction)
                              # or "enrichr" (legacy). Falls back to enrichr if
                              # gprofiler-official not installed.
```

**Add to `DESEQ2_COLS` (~line 107):**
```python
"stat":   "stat",       # Wald test statistic (for GSEA ranking)
"lfcSE":  "lfcSE",      # log2FC standard error (shrinkage detection)
```

**Add to `_DESEQ2_ALIASES` (~line 374):**
```python
"stat":    ["stat", "wald_statistic", "test_stat", "statistic"],
"lfcSE":   ["lfcse", "lfc_se", "std_error", "lfcstderror"],
"gene_id": [...existing..., "X"],  # Add "X" for R row.names export (WSF compat)
```

---

## Step 2: New Function — `load_counts_matrix()`

**Insert after `load_file()` (~line 267)**

- Load normalized counts TSV/CSV (genes x samples)
- **WSF compat**: Auto-detect R's row.names format (unnamed first column → use `index_col=0`)
- Strip Ensembl version suffixes from index
- Identify numeric sample columns
- **Auto-detect SAMPLE_METADATA from CONDITIONS**: If metadata not provided, try to match sample column names against condition names from `CONDITIONS` list (substring/prefix matching)
- Returns `(counts_df, metadata_dict)` or `(None, {})` on failure
- Print `[INFO] No counts file provided — skipping PCA, correlation heatmap, top DEG heatmap` when file not given

---

## Step 3: New Function — `pca_plot()`

**Insert after `load_counts_matrix()`**

Two modes:
- **Pre-computed**: Reads WSF `pca_data.csv` (columns: sample, PC1, PC2, condition, PC1_variance_pct, PC2_variance_pct) — preferred
- **From counts**: log2(counts+1) → StandardScaler → sklearn PCA(n_components=2)
  - **Note in output**: "PCA computed from log2(counts+1); for publication, use DESeq2 VST"

Okabe-Ito palette, variance % on axes, crosshair reference lines. Saves to `qc_plots/pca_plot.{format}`.

---

## Step 4: New Function — `sample_correlation_heatmap()`

**Insert after `pca_plot()`**

Euclidean distance on log2(counts+1), hierarchical clustering, seaborn clustermap with condition color annotations. Saves to `qc_plots/sample_correlation_heatmap.{format}`.

---

## Step 5: New Function — `top_deg_heatmap()`

**Insert after `sample_correlation_heatmap()`**

Top 50 genes by padj across all conditions (deduplicated). log2(counts+1) → z-score per gene (capped ±3). Seaborn clustermap, RdBu_r diverging palette, condition color bar, gene name labels. Saves to `cross_condition/figures/top_deg_heatmap.{format}`.

---

## Step 6: New Function — `run_gprofiler_ora()`

**Insert before `run_go_enrichment()` (~line 3349)**

- Uses `gprofiler-official` Python package
- Queries GO:BP, GO:MF, GO:CC, KEGG, Reactome with g:SCS FDR correction (hierarchy-aware, superior to BH)
- Processes up/down DEGs separately per condition
- Top 10 terms per database
- **CRITICAL**: Output DataFrame columns MUST match Enrichr schema exactly: `Term`, `Adjusted_P_value`, `Overlap_count`, `Category` — so `go_enrichment_combined_plot()` works unchanged
- Falls back to Enrichr if gprofiler-official not installed
- Species mapping: human→hsapiens, mouse→mmusculus, rat→rnorvegicus, zebrafish→drerio, fly→dmelanogaster, worm→celegans

---

## Step 7: Modify `run_gsea_enrichment()` (~line 3214)

**Change A — Ranking metric** (~line 3252):
- If `GSEA_RANKING == "stat"` AND `stat` column exists in deseq2_raw: rank by Wald statistic
- Else: fall back to log2FC with `[INFO]` message
- Update deduplication to use the chosen rank column

**Change B — Parameters** (~line 3276):
- `min_size=5` → `min_size=GSEA_MIN_SIZE`
- `max_size=500` → `max_size=GSEA_MAX_SIZE`
- `permutation_num=100` → `permutation_num=GSEA_PERMUTATIONS`

---

## Step 8: Modify `run_go_enrichment()` (~line 3349)

Add routing at function entry:
```python
if ORA_METHOD == "gprofiler":
    return run_gprofiler_ora(...)
```
Existing Enrichr code remains untouched as fallback.

---

## Step 9: Modify `normalize_deseq2_columns()` (~line 526)

After existing column resolution, detect `stat` and `lfcSE` columns (optional, silent if missing). If `lfcSE` found and median ratio lfcSE/|log2FC| < 0.3, print note about likely LFC shrinkage.

---

## Step 10: Modify `main()` (~line 5718)

**CRITICAL (from architecture review)**: Declare `counts_df = None` and `sample_metadata = {}` at main() entry (~line 5735, alongside `condition_results = {}`) to ensure scope across Phase 1 and Phase 2.

**Phase 1 insertion** (before per-condition loop, ~line 5737):
```python
if COUNTS_FILE:
    counts_df, sample_metadata = load_counts_matrix(COUNTS_FILE, SAMPLE_METADATA)
    if counts_df is not None:
        qc_dir = outdir / "qc_plots"
        pca_plot(counts_df=counts_df, metadata=sample_metadata, outdir=qc_dir)
        sample_correlation_heatmap(counts_df, sample_metadata, qc_dir)
else:
    print("[INFO] No counts file — skipping PCA, correlation heatmap, top DEG heatmap")

# Per-condition PCA files (from WSF output)
for cond in CONDITIONS:
    if cond.get("pca_file"):
        pca_plot(pca_file=cond["pca_file"], outdir=outdir/cond["name"]/"figures")
```

**Phase 2 insertion** (after log2FC heatmap, ~line 5881):
```python
if counts_df is not None:
    top_deg_heatmap(counts_df, condition_results, condition_labels,
                    sample_metadata, comparison_fig_dir)
```

---

## Step 11: Modify `run_pipeline()` (~line 5978)

**Add to existing global declarations** (after line 6013):
```python
global COUNTS_FILE, SAMPLE_METADATA, GSEA_RANKING, GSEA_MIN_SIZE, GSEA_MAX_SIZE
global GSEA_PERMUTATIONS, ORA_METHOD
```

Add `.get()` config injection:
```python
COUNTS_FILE       = str(config.get("COUNTS_FILE", COUNTS_FILE))
SAMPLE_METADATA   = dict(config.get("SAMPLE_METADATA", SAMPLE_METADATA))
GSEA_RANKING      = str(config.get("GSEA_RANKING", GSEA_RANKING))
GSEA_MIN_SIZE     = int(config.get("GSEA_MIN_SIZE", GSEA_MIN_SIZE))
GSEA_MAX_SIZE     = int(config.get("GSEA_MAX_SIZE", GSEA_MAX_SIZE))
GSEA_PERMUTATIONS = int(config.get("GSEA_PERMUTATIONS", GSEA_PERMUTATIONS))
ORA_METHOD        = str(config.get("ORA_METHOD", ORA_METHOD))
```

---

## Step 12: Update `pipeline_launcher.py`

**New variables in `__init__`:**
- `counts_file_var` (StringVar, default "")
- `gsea_ranking_var` (StringVar, default "stat")
- `gsea_min_size_var` (StringVar, default "15")
- `ora_method_var` (StringVar, default "gprofiler")

**Enrichment tab additions:**
- GSEA ranking dropdown: "Wald statistic (WSF method)" / "log2FC"
  - Tooltip: "Accounts for measurement uncertainty. More rigorous for gene ranking."
- GSEA min gene set size spinner (default 15)
- GSEA permutations spinner (default 1000)
- ORA method dropdown: "g:Profiler (WSF method)" / "Enrichr"
  - Tooltip: "Uses hierarchy-aware FDR correction. Matches WSF-Lab methodology."

**Output tab addition:**
- Counts file browse button (optional, with label: "For PCA/heatmaps")

**Update `_collect_config()`** with all new keys.

---

## Step 13: Update `requirements.txt`

Add:
```
gprofiler-official>=1.0
scikit-learn>=1.0
```

---

## Implementation: Subagent Delegation

1. **Subagent A** (pipeline engine): Steps 1-11 — all changes to `deseq2_rmats_filter_pipeline.py`
2. **Subagent B** (GUI launcher): Step 12 — all changes to `pipeline_launcher.py`
3. **Master**: Step 13 (requirements.txt), commit, push, run pipeline to test

Subagents A and B run in parallel.

---

## Verification

Run with WSF-style config:
```python
config = {
    "CONDITIONS": [...],
    "OUTPUT_DIR": "/mnt/f/RNA-seq_Analysis_2026_03_08_v5",
    # New defaults apply: GSEA_RANKING="stat", GSEA_MIN_SIZE=15, ORA_METHOD="gprofiler"
    # No COUNTS_FILE — should skip QC plots with INFO message
}
```

Checks:
- [ ] Pipeline runs without errors (no counts file, no stat column → graceful fallback)
- [ ] GSEA uses Wald stat when `stat` column present, falls back to log2FC with message
- [ ] GSEA uses min_size=15, permutation_num=1000
- [ ] g:Profiler ORA runs and produces results; falls back to Enrichr if package missing
- [ ] g:Profiler output compatible with `go_enrichment_combined_plot()` (same column schema)
- [ ] If `COUNTS_FILE` provided: PCA plot, correlation heatmap, top DEG heatmap generated
- [ ] WSF DESeq2 CSV with row.names (unnamed `X` column) loads correctly
- [ ] WSF `normalized_counts.tsv` parsed with gene IDs as index (not data column)
- [ ] Python PCA output includes note about VST approximation
- [ ] GUI launcher shows new dropdowns with WSF defaults pre-selected
- [ ] Old configs (without new keys) still work via `.get()` defaults
- [ ] `[INFO]` message printed when counts file not provided
- [ ] Commit + push to SpaceSorcerer/rnaseq-pipeline
