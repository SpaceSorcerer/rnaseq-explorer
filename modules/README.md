# RNA-seq Analysis Modules

Modular components extracted from the main analysis pipeline for reusability and maintainability.

## modules/deseq2.py

**Purpose**: DESeq2 data loading, filtering, and summarization

**Lines**: 434

**Functions**: 10 total (8 public, 2 private helpers)

### Public API

#### Data Loading & Normalization

**`load_and_normalize_deseq2(condition, base_dir, biotype_map=BIOTYPE_MAP)`**
- Loads DESeq2 Excel files with varying column schemas
- Standardizes column names (gene_id, gene_name, log2fc, padj, etc.)
- Looks up missing gene names via MyGene.info
- Normalizes biotype annotations
- Returns: DataFrame with standardized columns

**`lookup_gene_names_enhanced(ensembl_ids, species="human")`**
- Multi-source gene symbol lookup with fallback strategy
- Priority: MyGene.info → Ensembl BioMart → Use Ensembl ID
- Guarantees no NaN values
- Returns: Dict mapping Ensembl IDs to gene symbols

#### Data Filtering

**`filter_deseq2(df, log2fc_cutoff, basemean_cutoff, padj_cutoff, label="")`**
- Applies significance thresholds to DESeq2 results
- Filters by padj, log2FC, and baseMean
- Adds "direction" column (up/down)
- Returns: Tuple of (cleaned_df, filtered_degs)

#### Summary & Export

**`count_degs_by_direction(df_sig)`**
- Counts upregulated and downregulated DEGs
- Returns: Tuple of (n_up, n_down)

**`count_degs_by_biotype(df_sig)`**
- Counts DEGs per biotype category
- Returns: Series with counts

**`get_deg_summary(df_all, df_sig, condition_label)`**
- Comprehensive statistics for one condition
- Returns: Dict with total genes, DEGs, biotype breakdown

**`export_deg_results(df_sig, output_path, columns=None)`**
- Exports filtered DEGs to Excel
- Auto-selects relevant columns if not specified

### Constants

**`BIOTYPE_MAP`** - Dictionary mapping 26 raw biotype names to 5 standardized categories:
- Protein Coding
- lncRNA
- Pseudogene
- Small ncRNA
- Other

### Usage Example

```python
from modules import deseq2
from pathlib import Path

# Define condition with column mappings
condition = {
    "name": "MIAT_OE_vs_Control",
    "label": "MIAT OE vs Control",
    "deseq2_file": "NKX-MIAT3.vs.Control.xlsx",
    "columns": {
        "gene_id": "gene_id",
        "gene_name": "gene_name",
        "log2fc": "log2FoldChange",
        "basemean": "baseMean",
        "padj": "padj",
        "biotype": "biotype",
    }
}

# Load and normalize
base_dir = Path("/mnt/f/MIAT OE v QKI-KO v polyKQI-KO")
df = deseq2.load_and_normalize_deseq2(condition, base_dir)

# Filter by thresholds
df_clean, df_sig = deseq2.filter_deseq2(
    df, 
    log2fc_cutoff=0.4,
    basemean_cutoff=20,
    padj_cutoff=0.01,
    label="MIAT OE"
)

# Get summary statistics
summary = deseq2.get_deg_summary(df_clean, df_sig, "MIAT OE vs Control")
print(f"Found {summary['total_degs']} DEGs:")
print(f"  {summary['degs_up']} upregulated")
print(f"  {summary['degs_down']} downregulated")
print(f"  {summary['protein_coding_degs']} protein-coding")

# Export results
deseq2.export_deg_results(df_sig, Path("output/MIAT_OE_DEGs.xlsx"))
```

## modules/splicing.py

**Purpose**: rMATS alternative splicing analysis, event filtering, and set operations

**Lines**: 651

**Functions**: 16 total (all public)

### Public API

#### File Loading

**`load_rmats_file(filepath, event_type)`**
- Loads a single rMATS .MATS.JCEC.txt file
- Adds event_type column
- Returns: DataFrame with raw rMATS data

**`load_rmats_all_events(rmats_dir, event_types=None)`**
- Loads all event types from a directory
- Default event types: SE, A3SS, A5SS, RI, MXE
- Returns: Combined DataFrame with all events

**`load_rmats(condition, base_dir, fdr_cutoff=0.05, pvalue_cutoff=0.01, inclevel_diff_cutoff=0.1, use_pval_and_fdr=True, event_types=None)`**
- Main pipeline entry point
- Loads and filters rMATS results for one condition
- Supports dual FDR+p-value filtering
- Returns: Filtered significant splicing events

#### Event Filtering

**`filter_rmats_events(df, fdr_cutoff=0.05, pvalue_cutoff=None, inclevel_diff_cutoff=0.1, use_pval_and_fdr=True)`**
- Filters events by FDR, p-value, and dPSI thresholds
- Adds "direction" column (included/excluded)
- Returns: Filtered DataFrame

**`filter_rmats_by_event_type(df, event_type, **filter_kwargs)`**
- Filters for a specific event type (SE, A3SS, etc.)
- Returns: Filtered events of specified type

#### Summary & Statistics

**`summarize_splicing_events(df)`**
- Generates summary statistics by event type and direction
- Returns: Summary DataFrame with counts

**`count_events_by_type(df)`**
- Counts events per event type
- Returns: Dict mapping event_type → count

**`get_event_genes(df, event_type)`**
- Extracts gene symbols for specific event type
- Returns: List of gene symbols

#### Set Operations

**`get_event_identifiers(df, event_type=None)`**
- Extracts unique event identifiers (geneSymbol or ID)
- Returns: Set of event identifiers

**`compare_event_sets(df_a, df_b, event_type=None)`**
- Compares events between two conditions
- Returns: Tuple of (only_a, only_b, shared)

**`get_directional_event_sets(df, event_type=None)`**
- Splits events into included vs excluded sets
- Returns: Tuple of (included, excluded)

**`compare_directional_events(df_a, df_b, event_type=None)`**
- Compares directional splicing between conditions
- Returns: Dict with both_included, both_excluded, discordant, only_a, only_b

#### Export Functions

**`export_splicing_events(df, output_path, include_raw_counts=False)`**
- Exports filtered events to Excel
- Optionally includes raw junction counts
- Returns: None

**`export_splicing_summary(summary_df, output_path)`**
- Exports summary statistics to CSV/Excel
- Returns: None

#### Utilities

**`validate_rmats_dataframe(df)`**
- Validates required rMATS columns
- Returns: bool

**`get_event_type_description(event_type)`**
- Converts abbreviations to full names (SE → Skipped Exon)
- Returns: str

### Constants

**`DEFAULT_FDR_CUTOFF = 0.05`** - FDR significance threshold

**`DEFAULT_PVALUE_CUTOFF = 0.01`** - P-value threshold for dual filtering

**`DEFAULT_INCLEVEL_DIFF_CUTOFF = 0.1`** - Minimum absolute dPSI

**`DEFAULT_EVENT_TYPES = ["SE", "A3SS", "A5SS", "RI", "MXE"]`** - All event types

**`EVENT_COLORS`** - Dict mapping event types to publication-quality colors

### Usage Example

```python
from modules import splicing
from pathlib import Path

# Define condition
condition = {
    "name": "MIAT_OE_vs_Control",
    "label": "MIAT OE vs Control",
    "rmats_dir": "rMATs_MIAT.vs.Control"
}

# Load and filter splicing events
base_dir = Path("/mnt/f/MIAT OE v QKI-KO v polyKQI-KO")
rmats_df = splicing.load_rmats(
    condition,
    base_dir,
    fdr_cutoff=0.05,
    pvalue_cutoff=0.01,
    inclevel_diff_cutoff=0.1,
    use_pval_and_fdr=True
)

# Get summary statistics
summary = splicing.summarize_splicing_events(rmats_df)
print(summary)

# Count events by type
counts = splicing.count_events_by_type(rmats_df)
print(f"SE events: {counts.get('SE', 0)}")

# Compare two conditions
rmats_qki = splicing.load_rmats(qki_condition, base_dir)
only_miat, only_qki, shared = splicing.compare_event_sets(rmats_df, rmats_qki, event_type="SE")
print(f"Shared SE events: {len(shared)}")

# Compare directional overlap
comparison = splicing.compare_directional_events(rmats_df, rmats_qki, event_type="SE")
print(f"Concordant included: {len(comparison['both_included'])}")
print(f"Concordant excluded: {len(comparison['both_excluded'])}")
print(f"Discordant: {len(comparison['discordant'])}")

# Export results
splicing.export_splicing_events(rmats_df, Path("output/MIAT_OE_splicing.xlsx"))
```

## Design Principles

- All filtering logic preserved exactly from original script
- Functions accept parameters (no hardcoded paths or cutoffs)
- Comprehensive docstrings for all public functions
- Helper functions prefixed with `_` for clarity
- Module-level constants for shared configuration
- No external state or side effects (except file I/O)

## Integration

To use in `run_analysis_enhanced.py`:

```python
from modules import deseq2, splicing

# Replace direct function calls with module calls
df = deseq2.load_and_normalize_deseq2(condition, BASE_DIR)
df_clean, df_sig = deseq2.filter_deseq2(df, LOG2FC_CUTOFF, BASEMEAN_CUTOFF, PADJ_CUTOFF, label)

# Load splicing data
rmats_df = splicing.load_rmats(condition, BASE_DIR, FDR_CUTOFF, PVALUE_CUTOFF, INCLEVEL_DIFF_CUTOFF)
```
