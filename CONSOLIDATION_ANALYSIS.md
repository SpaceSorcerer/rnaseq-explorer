# RNA-seq Pipeline: Run Scripts Consolidation Analysis

**Date:** 2026-03-06  
**Directory:** `/workspace/extra/drive-f/RNA-seq_MIAT_OE_KCN_QKI-KO/`

---

## Executive Summary

Six entry-point scripts exist in the pipeline, representing **three generations of design evolution**:

1. **Generation 1 (Monolithic):** `run_analysis.py`, `run_analysis_enhanced.py` — hardcoded configuration, single-file implementations
2. **Generation 2 (CLI/GUI):** `run_pipeline_cli.py`, `pipeline_launcher.py` — config-driven wrappers with YAML support
3. **Generation 3 (Modular):** `run_full_pipeline.py`, `run_epistasis_analysis.py` — modular design with specialized analysis phases

**Status:** **Over-engineered with redundancy.** Multiple entry points serve overlapping purposes. A single, unified launcher is recommended.

---

## Detailed Script Analysis

### 1. run_analysis.py
**Type:** Monolithic monolith  
**LOC:** ~1,130  
**Purpose:** Self-contained DESeq2 + rMATS analysis for 3 conditions (MIAT OE, QKI-KO, polyQKI-KO)

**Key Features:**
- Hardcoded base directory: `/mnt/f/MIAT OE v QKI-KO v polyKQI-KO`
- Fixed cutoffs: log2FC ≥ 0.4, baseMean ≥ 20, padj < 0.01 (rMATS: FDR < 0.05)
- Handles **different column schemas** per DESeq2 file (normalization)
- Gene name lookup via MyGene.info API (with REST fallback)
- Core plots: volcano, MA, biotype distribution, biotype×direction
- rMATS plots: event summary, splicing volcano
- Cross-condition: overlap comparison, concordance scatters, global FC correlation, top-DEG heatmap
- Exports: Excel tables with filtered DEGs, full gene lists, shared gene lists, summary statistics

**Unique Features:**
- Handles column name heterogeneity across three different DESeq2 files
- Gene name lookup with batch API calls
- Biotype normalization with granular mapping

**Config Mechanism:** None — fully hardcoded

**Output:** 20+ PNG figures + 10+ Excel tables, all to `pipeline_output/`

---

### 2. run_analysis_enhanced.py
**Type:** Monolithic evolution  
**LOC:** ~1,735  
**Purpose:** Extended version of run_analysis.py with richer comparative analysis

**Key Features:**
- **All features from run_analysis.py**, plus:
- Enhanced gene name lookup: MyGene.info → BioMart REST → fallback to Ensembl ID
- **Directional Venn diagrams:** 3-panel for concordant UP, DOWN, all DEGs
- **Splicing Venn diagrams:** 3-panel per event type (SE, A3SS, etc.)
- **Biotype-specific volcano plots:** Separate plots for protein-coding vs lncRNA+other
- **Violin plots:** log2FC and dPSI distributions across conditions
- **Directional Excel export:** Multi-sheet workbook with concordant UP/DOWN/discordant/unique genes
- Dual filtering for rMATS: FDR AND p-value (configurable)

**Config Mechanism:** Limited hardcoded flags:
- `USE_PVAL_AND_FDR` = True
- `USE_PVAL_FOR_SPLICING_VOLCANO` = True

**Output:** 30+ PNG figures + 15+ Excel tables

**vs run_analysis.py:**
- 60% larger codebase
- 50% more output files
- Better handling of missing gene names
- More sophisticated comparative visualization

---

### 3. run_full_pipeline.py
**Type:** Modular orchestrator  
**LOC:** ~433  
**Purpose:** Master pipeline runner that chains multiple analysis phases via YAML config

**Config Mechanism:** YAML configuration file (required: `--config <file>`)

**Architecture:**
```
PipelineRunner class:
  ├─ Phase 0: Configuration loading & validation
  ├─ Phase 1: Epistasis analysis (2×2 factorial design)
  ├─ Phase 2: VAST-TOOLS integration (optional)
  ├─ Phase 3: Cross-omics integration (expression + splicing)
  └─ Phase 4: Publication figures
```

**Features:**
- Orchestrates external modules (not self-contained):
  - `config_loader`, `epistasis`, `directional_concordance`
  - `cross_omics_integration`, `publication_figures`
  - Optional: `vasttools_parser`, `rmats_vast_comparison`
- Selective phase execution via CLI: `--phases epistasis cross_omics figures`
- Skip-figures option: `--skip-figures`
- Graceful fallback to mock data if files missing

**Config Format:** YAML with sections:
```yaml
factorial:
  enabled: true
  data_files: {WT: ..., MIAT_KD: ..., QKI_KO: ..., DKO: ...}
  value_col: log2FoldChange
  statistic: median
  n_bootstrap: 10000
  ci_level: 0.95

vast_tools:
  enabled: false
  diff_file: path/to/vast.txt
  
cross_omics:
  enabled: true
  expression: {...}
  splicing: {...}
  
figures:
  output_dir: figures
```

**Command-line Interface:**
```bash
python run_full_pipeline.py --config rnaseq_config.yaml
python run_full_pipeline.py --config rnaseq_config.yaml --phases epistasis cross_omics
python run_full_pipeline.py --config rnaseq_config.yaml --skip-figures
```

**Unique Features:**
- Factorial epistasis design (2×2: WT, MIAT-KD, QKI-KO, DKO)
- VAST-TOOLS optional integration
- Cross-omics (expression×splicing) correlation
- Modular: can run individual phases
- Mock data fallback for testing

**Dependencies:** Requires many external modules not provided in script review

---

### 4. run_pipeline_cli.py
**Type:** CLI wrapper  
**LOC:** ~163  
**Purpose:** Command-line interface to existing monolithic pipeline via YAML config

**Config Mechanism:** YAML configuration file (required: `--config <file>`)

**Architecture:**
- Wraps `deseq2_rmats_filter_pipeline` module
- Loads YAML config via `config_loader.load_config()`
- Injects global variables into pipeline namespace
- Runs pipeline's `main()` function

**Command-line Interface:**
```bash
python run_pipeline_cli.py --config rnaseq_config.yaml
python run_pipeline_cli.py --config my_analysis.yaml --output results/
python run_pipeline_cli.py --config rnaseq_config.yaml --base-dir /path/to/data/
python run_pipeline_cli.py --config rnaseq_config.yaml --log2fc-cutoff 0.5 --padj-cutoff 0.05
```

**Features:**
- Override any cutoff via CLI args (takes precedence over YAML)
- Converts config to global variables matching pipeline expectations
- Dependency: `deseq2_rmats_filter_pipeline.py` must exist in same directory
- Dependency: `config_loader.py` for YAML parsing

**Unique Aspects:**
- Thinnest wrapper (163 LOC)
- Pure CLI approach
- No GUI
- For automation/scripting

---

### 5. pipeline_launcher.py
**Type:** Interactive GUI launcher  
**LOC:** ~790  
**Purpose:** Tkinter GUI to configure and run pipeline without writing code

**Config Mechanism:**
1. **Interactive GUI** with 4 tabs
2. **YAML Load/Save** (Phase 0 support)
3. **Direct variable entry** into GUI fields

**GUI Tabs:**
- **Conditions:** Add/remove conditions (max 5), set DESeq2 files + rMATS dirs, browse file dialogs
- **Thresholds:** DESeq2 cutoffs (log2FC, baseMean, padj), rMATS cutoffs (FDR, p-value, dPSI)
- **Output & Figures:** Output directory, format (PNG/SVG/PDF), DPI, font size, color picker for up/down/NS
- **Column Names (Advanced):** Customize DESeq2 and rMATS column mappings

**Features:**
- Dynamically add/remove conditions (UI updates in real-time)
- File browser dialogs for DESeq2 files and rMATS directories
- Color picker widget for publication colors
- Validation before pipeline run
- Load/save YAML configs
- Thread-safe stdout/stderr redirection to log widget
- Status indicator (Ready/Running/Completed/Error)

**Command-line Interface:**
```bash
python pipeline_launcher.py  # Launches GUI, no CLI args
```

**Unique Features:**
- **Only GUI option in the suite**
- No CLI needed for configuration
- Beautiful Tkinter interface with tabs
- YAML persistence
- Real-time validation with error messages
- Thread-safe logging to text widget

**Dependencies:** Requires `tkinter`, `config_loader.py`, `deseq2_rmats_filter_pipeline.py`

---

### 6. run_epistasis_analysis.py
**Type:** Specialized standalone  
**LOC:** ~383  
**Purpose:** Dedicated pipeline for 2×2 factorial epistasis design (WT, MIAT-KD, QKI-KO, DKO)

**Config Mechanism:** YAML configuration file (required: `--config <file>`)

**Architecture:**
```
run_epistasis_pipeline() function:
  [1/6] Calculate per-gene epistasis
  [2/6] Analyze directional concordance
  [3/6] Merge epistasis + concordance
  [4/6] Summarize concordance by epistasis category
  [5/6] Identify ceRNA candidates
  [6/6] Generate visualizations (volcano, distribution, heatmap)
```

**Features:**
- Loads DESeq2 for all 4 conditions from `factorial.data_files` in config
- Computes epistasis per gene (bootstrap CI optional)
- Directional concordance analysis
- ceRNA candidate identification
- Epistasis volcano plot
- Epistasis distribution + heatmap visualizations
- Export: CSV tables, PNG figures

**Command-line Interface:**
```bash
python run_epistasis_analysis.py --config rnaseq_config.yaml
python run_epistasis_analysis.py --config rnaseq_config.yaml --output-dir my_results/
```

**Unique Features:**
- **Only script for 2×2 factorial design**
- Epistasis + concordance integration
- ceRNA candidate prediction
- Specialized visualizations

**Dependencies:** External modules:
- `config_loader`, `epistasis`, `directional_concordance`, `epistasis_plots`

---

## Comparative Feature Matrix

| Feature | run_analysis | enhanced | full_pipeline | cli | launcher | epistasis |
|---------|:---:|:---:|:---:|:---:|:---:|:---:|
| **Hardcoded config** | ✓ | ✓ | ✗ | ✗ | ✗ | ✗ |
| **YAML config** | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ |
| **GUI** | ✗ | ✗ | ✗ | ✗ | ✓ | ✗ |
| **CLI flags** | ✗ | ✗ | ~partial | ✓ | ✗ | ~partial |
| **Monolithic code** | ✓ | ✓ | ✗ | wrapper | wrapper | ✗ |
| **DESeq2+rMATS** | ✓ | ✓ | ~via modules | ✓ | ✓ | ✗ |
| **Epistasis analysis** | ✗ | ✗ | ✓ | ✗ | ✗ | ✓ |
| **3-panel Venn plots** | ✗ | ✓ | ✗ | ✗ | ✗ | ✗ |
| **Directional Excel export** | ✗ | ✓ | ✗ | ✗ | ✗ | ✗ |
| **Violin plots** | ✗ | ✓ | ✗ | ✗ | ✗ | ✗ |
| **VAST-TOOLS integration** | ✗ | ✗ | ✓ | ✗ | ✗ | ✗ |
| **Cross-omics integration** | ✗ | ✗ | ✓ | ✗ | ✗ | ✗ |
| **Thread-safe output logging** | ✗ | ✗ | ✗ | ✗ | ✓ | ✗ |
| **Mock data fallback** | ✗ | ✗ | ✓ | ✗ | ✗ | ✗ |
| **Color picker widget** | ✗ | ✗ | ✗ | ✗ | ✓ | ✗ |

---

## Functional Overlap & Redundancy

### **Overlapping Functionality:**

1. **Basic DESeq2 + rMATS analysis:**
   - Implemented in: `run_analysis.py`, `run_analysis_enhanced.py`, `run_pipeline_cli.py`, `pipeline_launcher.py`
   - **Problem:** Code duplication; config scattered across scripts

2. **Configuration loading:**
   - Hardcoded: `run_analysis.py`, `run_analysis_enhanced.py`
   - YAML wrapper: `run_pipeline_cli.py`, `pipeline_launcher.py`
   - Modular: `run_full_pipeline.py`, `run_epistasis_analysis.py`
   - **Problem:** Three different approaches; inconsistent

3. **Column name mapping:**
   - Handled in: `run_analysis.py`, `run_analysis_enhanced.py`
   - **Problem:** Not exposed to CLI/GUI; hardcoded condition definitions

4. **Figure generation:**
   - Implemented in: `run_analysis.py`, `run_analysis_enhanced.py`
   - Abstracted to modules in: `run_full_pipeline.py`
   - **Problem:** Monolithic scripts can't reuse modular figures

5. **Gene name lookup:**
   - In: `run_analysis.py` (basic), `run_analysis_enhanced.py` (enhanced with BioMart)
   - **Problem:** Duplicated logic, not available to other scripts

---

## Critical Design Issues

### 1. **Configuration Fragmentation**
- **run_analysis.py**: All cutoffs hardcoded in Python (lines 54–96)
- **run_analysis_enhanced.py**: Same hardcoding + flags for dual-filtering
- **run_pipeline_cli.py**: Config injected as globals (fragile metaprogramming)
- **pipeline_launcher.py**: Collected from GUI, serialized to YAML
- **run_full_pipeline.py**: Proper YAML loading
- **run_epistasis_analysis.py**: Proper YAML loading

**Impact:** No single source of truth. Users must edit different places depending on entry point.

### 2. **Module Dependency Hell**
- `run_full_pipeline.py` requires: `config_loader`, `epistasis`, `directional_concordance`, `cross_omics_integration`, `publication_figures`
- `run_epistasis_analysis.py` requires: `config_loader`, `epistasis`, `directional_concordance`, `epistasis_plots`
- `run_pipeline_cli.py` requires: `config_loader`, `deseq2_rmats_filter_pipeline`
- **Problem:** No shared interface. Each script imports different modules. Tight coupling to external code.

### 3. **Code Duplication**
- Gene name lookup logic: ~100 LOC in `run_analysis.py`, ~120 LOC in `run_analysis_enhanced.py` (enhanced version)
- Biotype normalization: repeated
- Column mapping logic: repeated (though with variations)
- **Problem:** Bug fixes must be propagated manually

### 4. **No Common Entry Point**
- User must choose between:
  - Direct Python run (requires editing code)
  - CLI with YAML
  - GUI
  - Dedicated epistasis script
  - Full modular pipeline

**Impact:** Confusing UX. New users don't know which to use. No obvious upgrade path.

### 5. **Incomplete Feature Parity**
| Analysis Type | Hardcoded | CLI | GUI | Modular | Epistasis |
|---|:---:|:---:|:---:|:---:|:---:|
| Basic DESeq2 | ✓ | ✓ | ✓ | ~ | ✗ |
| rMATS | ✓ | ✓ | ✓ | ~ | ✗ |
| Directional analysis | enhanced | ✗ | ✗ | ✗ | ✗ |
| Epistasis | ✗ | ✗ | ✗ | ✓ | ✓ |
| VAST-TOOLS | ✗ | ✗ | ✗ | ✓ | ✗ |

**Problem:** Feature X works in script Y but not in script Z. No consistent way to access all functionality.

---

## Recommended Consolidation Approach

### **Option A: Single Entry Point (Recommended)**

Create one unified launcher that includes all functionality:

```
unified_launcher.py  (choose mode: gui | cli | api)
├── GUI mode → Tkinter widget (use pipeline_launcher.py as base)
├── CLI mode → argparse interface (use run_pipeline_cli.py as base)
└── API mode → Python function import (expose run_epistasis_pipeline)

Shared core modules:
├── config.py (unified YAML schema)
├── analysis_engine.py (orchestrate phases)
├── plots.py (all visualization functions)
└── data_loader.py (read/normalize DESeq2 + rMATS)
```

**Pros:**
- Single entry point for all workflows
- Unified configuration schema
- One source of truth
- Consistent feature parity

**Cons:**
- Larger refactoring effort (~3-4 hours)
- Breaking change for existing automation

---

### **Option B: Modular Consolidation (Fast Path)**

Refactor existing scripts to use shared modules:

```
Phase 1: Extract reusable components
├── gene_lookup.py (rename from scattered code)
├── column_mapper.py (abstracted from run_analysis.py)
├── biotype_normalizer.py (reusable mapping)
└── common_plots.py (shared visualization library)

Phase 2: Create common interface
├── BaseAnalysisRunner (abstract class)
├── DESeq2_RMATSRunner (extends base)
├── EpistasisRunner (extends base)
└── MultiOmicsRunner (extends base)

Phase 3: Update launchers
├── pipeline_launcher.py (GUI) → uses BaseAnalysisRunner
├── run_pipeline_cli.py → uses BaseAnalysisRunner
└── run_epistasis_analysis.py → uses EpistasisRunner
```

**Pros:**
- Incremental changes
- Preserves existing scripts
- Can implement in parallel

**Cons:**
- Leaves some redundancy
- Still multiple entry points

---

## Modules Required for Consolidation

**Currently Missing (Not in Script Review):**
- `config_loader.py` — YAML config handling
- `epistasis.py` — 2×2 factorial analysis
- `directional_concordance.py` — concordance patterns
- `cross_omics_integration.py` — expression×splicing
- `publication_figures.py` — figure generation
- `epistasis_plots.py` — epistasis visualizations
- `vasttools_parser.py` — VAST-TOOLS parsing (optional)

**To be created for consolidation:**
- `unified_launcher.py` — Single entry point (GUI + CLI + API)
- `base_analysis_engine.py` — Shared orchestration logic
- `shared_plots.py` — Unified figure generation (consolidate from run_analysis*.py)
- `data_normalization.py` — Column mapping, gene lookup, biotype normalization

---

## Quick Wins (No Refactoring)

1. **Document which script to use:**
   - For quick analysis: `pipeline_launcher.py` (GUI)
   - For automation: `run_pipeline_cli.py` (CLI)
   - For epistasis: `run_epistasis_analysis.py`
   - For advanced: `run_full_pipeline.py`

2. **Create README.md:**
   ```
   ## Choosing an Entry Point
   
   **I want a GUI:** python pipeline_launcher.py
   **I want command-line:** python run_pipeline_cli.py --config my.yaml
   **I want full features:** python run_full_pipeline.py --config my.yaml
   **I want epistasis analysis:** python run_epistasis_analysis.py --config my.yaml
   **I want to edit code:** edit run_analysis.py or run_analysis_enhanced.py
   ```

3. **Standardize YAML schema** across all scripts (use `run_full_pipeline.py` as template)

4. **Add validation** to ensure all required files exist before running

---

## Consolidation Implementation Timeline

| Phase | Duration | Work | Deliverable |
|-------|----------|------|-------------|
| **Phase 1: Audit** | 2h | Document findings (this file) | consolidation_analysis.md |
| **Phase 2: Extract** | 4h | Pull out gene_lookup, column_mapper, biotype, plots into modules | shared module library |
| **Phase 3: Unify Config** | 3h | Single YAML schema, validation | `config.py` with schema |
| **Phase 4: Base Runner** | 4h | Abstract class with phases | `base_analysis_engine.py` |
| **Phase 5: Update Launchers** | 5h | Refactor GUI, CLI to use base class | Updated `pipeline_launcher.py`, `run_pipeline_cli.py` |
| **Phase 6: Consolidate Plots** | 6h | Merge `run_analysis` figure code into unified library | `shared_plots.py` |
| **Phase 7: Testing & Docs** | 4h | Test all entry points, update README | Integration tests, documentation |
| **Total** | ~28h | Full consolidation | Single unified launcher (Option A) |

**Fast Path (Option B):** ~12h (just extract Phase 1-3, document, leave scripts as-is)

---

## Conclusion

**Current state:** Three parallel implementation efforts with overlapping functionality, making the codebase harder to maintain and more error-prone.

**Most viable path:** Implement **Option B (Modular Consolidation)** over 2-3 sprints, prioritizing:
1. Extract gene_lookup, biotype, column_mapper → shared modules
2. Unify YAML schema
3. Create integration test suite
4. Document user-facing decision tree

Once modularization is complete, can consider **Option A (Single Entry Point)** as a longer-term refactor.

---

## References

- **Hardcoded analysis:** run_analysis.py (l. 54-96), run_analysis_enhanced.py (l. 54-100)
- **Config-driven modular:** run_full_pipeline.py (PipelineRunner class)
- **GUI launcher:** pipeline_launcher.py (PipelineLauncherApp class)
- **CLI wrapper:** run_pipeline_cli.py (thin wrapper pattern)
- **Specialized analysis:** run_epistasis_analysis.py (factorial design)

