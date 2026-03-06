"""
DESeq2 Data Loading and Filtering Module
=========================================

Functions for loading, normalizing, and filtering DESeq2 differential expression data.
Handles multiple file formats with different column schemas and applies configurable
significance cutoffs.

Functions:
----------
- load_and_normalize_deseq2: Load and standardize DESeq2 Excel files
- filter_deseq2: Apply significance cutoffs to DESeq2 data
- lookup_gene_names: Fetch gene symbols from MyGene.info
- lookup_gene_names_enhanced: Multi-source gene name lookup with fallback
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Biotype grouping map
BIOTYPE_MAP = {
    "protein_coding": "Protein Coding",
    "lncrna": "lncRNA", "lincrna": "lncRNA", "sense_intronic": "lncRNA",
    "sense_overlapping": "lncRNA", "antisense": "lncRNA",
    "processed_transcript": "lncRNA", "bidirectional_promoter_lncrna": "lncRNA",
    "macro_lncrna": "lncRNA", "non_coding": "lncRNA",
    "pseudogene": "Pseudogene", "processed_pseudogene": "Pseudogene",
    "unprocessed_pseudogene": "Pseudogene",
    "transcribed_unprocessed_pseudogene": "Pseudogene",
    "transcribed_processed_pseudogene": "Pseudogene",
    "transcribed_unitary_pseudogene": "Pseudogene",
    "polymorphic_pseudogene": "Pseudogene", "unitary_pseudogene": "Pseudogene",
    "mirna": "Small ncRNA", "snrna": "Small ncRNA", "snorna": "Small ncRNA",
    "misc_ncrna": "Small ncRNA", "rrna": "Small ncRNA", "scrna": "Small ncRNA",
    "scarna": "Small ncRNA", "trna": "Small ncRNA",
}


# ─── Gene Name Lookup ────────────────────────────────────────────────────────

def lookup_gene_names(ensembl_ids, species="human"):
    """
    Fetch gene symbols from MyGene.info for Ensembl IDs.

    Args:
        ensembl_ids: List of Ensembl gene IDs (with or without version numbers)
        species: Species name for MyGene.info query (default: "human")

    Returns:
        Dictionary mapping Ensembl IDs to gene symbols
    """
    try:
        import mygene
        mg = mygene.MyGeneInfo()
    except ImportError:
        # Fallback: try REST API
        return _lookup_gene_names_rest(ensembl_ids)

    # Strip version numbers
    clean_ids = [eid.split(".")[0] for eid in ensembl_ids]
    unique_ids = list(set(clean_ids))

    print(f"  Looking up {len(unique_ids)} gene names from MyGene.info...")

    mapping = {}
    batch_size = 1000
    for i in range(0, len(unique_ids), batch_size):
        batch = unique_ids[i:i + batch_size]
        try:
            results = mg.querymany(batch, scopes="ensembl.gene",
                                   fields="symbol", species=species,
                                   returnall=True)
            for hit in results.get("out", []):
                if "symbol" in hit and "query" in hit:
                    mapping[hit["query"]] = hit["symbol"]
        except Exception as e:
            print(f"  Warning: gene lookup batch failed: {e}")

    # Build full mapping with version numbers
    full_mapping = {}
    for eid in ensembl_ids:
        clean = eid.split(".")[0]
        full_mapping[eid] = mapping.get(clean, eid)

    found = sum(1 for v in full_mapping.values() if not v.startswith("ENS"))
    print(f"  Resolved {found}/{len(ensembl_ids)} gene names")
    return full_mapping


def _lookup_gene_names_rest(ensembl_ids):
    """
    Fallback: REST API gene lookup without mygene package.

    Args:
        ensembl_ids: List of Ensembl gene IDs

    Returns:
        Dictionary mapping Ensembl IDs to gene symbols
    """
    import json
    import urllib.request

    clean_ids = list(set(eid.split(".")[0] for eid in ensembl_ids))
    mapping = {}

    print(f"  Looking up {len(clean_ids)} gene names via REST API...")

    batch_size = 500
    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i + batch_size]
        try:
            data = json.dumps({"ids": batch}).encode()
            req = urllib.request.Request(
                "https://mygene.info/v3/gene",
                data=data,
                headers={"Content-Type": "application/json"},
                method="POST"
            )
            with urllib.request.urlopen(req, timeout=30) as resp:
                results = json.loads(resp.read())
                for hit in results:
                    if isinstance(hit, dict) and "symbol" in hit:
                        qid = hit.get("query", hit.get("_id", ""))
                        mapping[qid] = hit["symbol"]
        except Exception as e:
            print(f"  Warning: REST lookup batch failed: {e}")

    full_mapping = {}
    for eid in ensembl_ids:
        clean = eid.split(".")[0]
        full_mapping[eid] = mapping.get(clean, eid)

    found = sum(1 for v in full_mapping.values() if not v.startswith("ENS"))
    print(f"  Resolved {found}/{len(ensembl_ids)} gene names")
    return full_mapping


def _lookup_biomart(ensembl_ids):
    """
    Enhanced fallback: Ensembl BioMart REST API.

    Args:
        ensembl_ids: List of Ensembl gene IDs

    Returns:
        Dictionary mapping Ensembl IDs to gene symbols
    """
    import urllib.request

    clean_ids = [eid.split(".")[0] for eid in ensembl_ids]
    mapping = {}

    print(f"  Trying BioMart for {len(clean_ids)} remaining IDs...")

    # BioMart query - use simple TSV format
    batch_size = 500
    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i + batch_size]
        try:
            # Use Ensembl REST API (simpler than BioMart XML)
            url = "https://rest.ensembl.org/lookup/id"
            data = {"ids": batch}
            import json
            req = urllib.request.Request(
                url,
                data=json.dumps(data).encode(),
                headers={"Content-Type": "application/json"},
                method="POST"
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                results = json.loads(resp.read())
                for eid, info in results.items():
                    if isinstance(info, dict) and "display_name" in info:
                        mapping[eid] = info["display_name"]
        except Exception as e:
            print(f"  BioMart batch failed: {e}")

    # Map back with version numbers
    full_mapping = {}
    for eid in ensembl_ids:
        clean = eid.split(".")[0]
        full_mapping[eid] = mapping.get(clean, eid)

    found = sum(1 for v in full_mapping.values() if not v.startswith("ENS"))
    print(f"  BioMart resolved {found}/{len(ensembl_ids)} additional names")
    return full_mapping


def lookup_gene_names_enhanced(ensembl_ids, species="human"):
    """
    Enhanced gene name lookup with multiple fallback sources.

    Priority:
    1. MyGene.info API (existing)
    2. Ensembl BioMart REST API
    3. Return Ensembl ID (no NaN - guaranteed)

    Args:
        ensembl_ids: List of Ensembl gene IDs
        species: Species name for MyGene.info query (default: "human")

    Returns:
        Dictionary mapping Ensembl IDs to gene symbols (guaranteed non-null)
    """
    # Step 1: Try MyGene.info (existing code)
    mapping = lookup_gene_names(ensembl_ids, species)

    # Step 2: For unmapped IDs, try BioMart
    unmapped = [eid for eid in ensembl_ids
                if mapping.get(eid, "").startswith("ENS")]

    if unmapped and len(unmapped) < len(ensembl_ids) * 0.9:  # Only if <90% failed
        biomart_mapping = _lookup_biomart(unmapped)
        for eid in unmapped:
            if eid in biomart_mapping and not biomart_mapping[eid].startswith("ENS"):
                mapping[eid] = biomart_mapping[eid]

    # Step 3: Final fallback - use Ensembl ID itself (ZERO NaN guarantee)
    for eid in ensembl_ids:
        if eid not in mapping or pd.isna(mapping[eid]) or mapping[eid] == "":
            mapping[eid] = eid  # Use Ensembl ID as name

    # Final count
    found = sum(1 for v in mapping.values() if not v.startswith("ENS"))
    print(f"  FINAL: Resolved {found}/{len(ensembl_ids)} to gene symbols, "
          f"{len(ensembl_ids)-found} using Ensembl IDs")

    return mapping


# ─── DESeq2 Loading & Filtering ─────────────────────────────────────────────

def load_and_normalize_deseq2(condition, base_dir, biotype_map=BIOTYPE_MAP):
    """
    Load DESeq2 Excel file and normalize columns to standard names.

    Handles different column schemas across files:
    - Renames columns according to condition-specific mapping
    - Fills missing gene names via MyGene.info lookup
    - Standardizes biotype annotations
    - Strips Ensembl version numbers for cross-condition matching

    Args:
        condition: Dictionary with file path and column mappings
        base_dir: Base directory path for input files
        biotype_map: Dictionary mapping raw biotypes to grouped categories

    Returns:
        DataFrame with standardized columns:
        - gene_id: Ensembl gene ID
        - gene_name: Gene symbol
        - log2fc: log2 fold change
        - basemean: Base mean expression
        - padj: Adjusted p-value
        - pvalue: Nominal p-value
        - biotype: Raw biotype annotation
        - biotype_group: Standardized biotype category
        - gene_id_base: Ensembl ID without version number
    """
    filepath = Path(base_dir) / condition["deseq2_file"]
    col_map = condition["columns"]
    label = condition["label"]

    print(f"\n{'='*60}")
    print(f"Loading DESeq2: {label}")
    print(f"  File: {filepath.name}")

    df = pd.read_excel(filepath)
    print(f"  Raw rows: {len(df):,}")
    print(f"  Raw columns: {list(df.columns)}")

    # Rename columns to standard names
    rename_map = {}
    for std_name, file_col in col_map.items():
        if file_col is not None and file_col in df.columns:
            rename_map[file_col] = std_name

    df = df.rename(columns=rename_map)
    print(f"  Renamed: {rename_map}")

    # Handle missing gene_name
    if "gene_name" not in df.columns or df["gene_name"].isna().all():
        if "gene_id" in df.columns:
            gene_names = lookup_gene_names_enhanced(df["gene_id"].tolist())
            df["gene_name"] = df["gene_id"].map(gene_names)
    # Also fill any remaining NaN gene_name values
    elif "gene_name" in df.columns and df["gene_name"].isna().any():
        missing_mask = df["gene_name"].isna()
        if "gene_id" in df.columns and missing_mask.sum() > 0:
            missing_ids = df.loc[missing_mask, "gene_id"].tolist()
            gene_names = lookup_gene_names_enhanced(missing_ids)
            df.loc[missing_mask, "gene_name"] = df.loc[missing_mask, "gene_id"].map(gene_names)

    # Handle missing biotype
    if "biotype" not in df.columns:
        df["biotype"] = "unknown"
        print("  Note: No biotype column — set to 'unknown'")

    # Normalize biotype values
    df["biotype_raw"] = df["biotype"].copy()
    df["biotype_group"] = df["biotype"].str.lower().str.strip().map(biotype_map).fillna("Other")

    # Strip Ensembl version for cross-condition matching
    if "gene_id" in df.columns:
        df["gene_id_base"] = df["gene_id"].str.split(".").str[0]

    print(f"  Standardized columns: {[c for c in df.columns if c in ['gene_id','gene_name','log2fc','basemean','padj','biotype','biotype_group']]}")

    return df


def filter_deseq2(df, log2fc_cutoff, basemean_cutoff, padj_cutoff, label=""):
    """
    Apply significance cutoffs to DESeq2 data.

    Filters genes by:
    - Adjusted p-value < padj_cutoff
    - |log2 fold change| >= log2fc_cutoff
    - Base mean expression >= basemean_cutoff

    Args:
        df: DataFrame with DESeq2 results
        log2fc_cutoff: Minimum absolute log2 fold change
        basemean_cutoff: Minimum base mean expression
        padj_cutoff: Maximum adjusted p-value
        label: Condition label for logging

    Returns:
        Tuple of (cleaned_df, filtered_df):
        - cleaned_df: All genes after dropping NaN in critical columns
        - filtered_df: Significant DEGs with "direction" column (up/down)
    """
    # Drop NaN in critical columns
    df_clean = df.dropna(subset=["padj", "log2fc", "basemean"]).copy()

    mask = (
        (df_clean["padj"] < padj_cutoff) &
        (df_clean["log2fc"].abs() >= log2fc_cutoff) &
        (df_clean["basemean"] >= basemean_cutoff)
    )
    filtered = df_clean[mask].copy()
    filtered["direction"] = np.where(filtered["log2fc"] > 0, "up", "down")

    n_up = (filtered["direction"] == "up").sum()
    n_down = (filtered["direction"] == "down").sum()
    print(f"  {label} filtered: {len(filtered):,} DEGs ({n_up:,} up, {n_down:,} down)")
    print(f"    Cutoffs: |log2FC| >= {log2fc_cutoff}, baseMean >= {basemean_cutoff}, padj < {padj_cutoff}")

    return df_clean, filtered


# ─── DEG Summary Functions ──────────────────────────────────────────────────

def count_degs_by_direction(df_sig):
    """
    Count upregulated and downregulated DEGs.

    Args:
        df_sig: DataFrame with filtered DEGs containing "direction" column

    Returns:
        Tuple of (n_up, n_down)
    """
    if "direction" not in df_sig.columns or len(df_sig) == 0:
        return 0, 0

    n_up = (df_sig["direction"] == "up").sum()
    n_down = (df_sig["direction"] == "down").sum()
    return n_up, n_down


def count_degs_by_biotype(df_sig):
    """
    Count DEGs by biotype category.

    Args:
        df_sig: DataFrame with filtered DEGs containing "biotype_group" column

    Returns:
        Series with counts per biotype category
    """
    if "biotype_group" not in df_sig.columns or len(df_sig) == 0:
        return pd.Series(dtype=int)

    return df_sig["biotype_group"].value_counts()


def get_deg_summary(df_all, df_sig, condition_label):
    """
    Generate comprehensive summary statistics for a single condition.

    Args:
        df_all: All genes (after cleanup)
        df_sig: Significant DEGs
        condition_label: Condition name for display

    Returns:
        Dictionary with summary statistics
    """
    n_up, n_down = count_degs_by_direction(df_sig)
    biotype_counts = count_degs_by_biotype(df_sig)

    return {
        "condition": condition_label,
        "total_genes": len(df_all),
        "total_degs": len(df_sig),
        "degs_up": n_up,
        "degs_down": n_down,
        "protein_coding_degs": biotype_counts.get("Protein Coding", 0),
        "lncrna_degs": biotype_counts.get("lncRNA", 0),
        "pseudogene_degs": biotype_counts.get("Pseudogene", 0),
        "small_ncrna_degs": biotype_counts.get("Small ncRNA", 0),
        "other_degs": biotype_counts.get("Other", 0),
    }


def export_deg_results(df_sig, output_path, columns=None):
    """
    Export filtered DEG results to Excel.

    Args:
        df_sig: DataFrame with significant DEGs
        output_path: Path object for output file
        columns: List of columns to export (default: all relevant columns)
    """
    if columns is None:
        columns = ["gene_id", "gene_name", "log2fc", "basemean", "padj",
                  "direction", "biotype_group"]

    # Filter to only existing columns
    export_cols = [c for c in columns if c in df_sig.columns]
    df_sig[export_cols].to_excel(output_path, index=False)
    print(f"  Exported: {output_path.name}")
