#!/usr/bin/env python3
"""
Standardization and Gene Lookup Utilities for RNA-seq Analysis
===============================================================

This module provides functions for:
- Gene ID/symbol lookup via MyGene.info, BioMart, and Ensembl REST APIs
- Biotype normalization and mapping
- Data loading and preprocessing utilities

Extracted from run_analysis_enhanced.py for reusability and maintainability.
"""

import numpy as np
import pandas as pd


# ─── Biotype Mapping ─────────────────────────────────────────────────────────

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


def normalize_biotype(biotype_series):
    """
    Normalize biotype names to standardized groups.

    Args:
        biotype_series: pandas Series containing biotype names

    Returns:
        pandas Series with normalized biotype group names
    """
    return biotype_series.str.lower().str.strip().map(BIOTYPE_MAP).fillna("Other")


# ─── Gene Name Lookup ────────────────────────────────────────────────────────

def lookup_gene_names(ensembl_ids, species="human"):
    """
    Fetch gene symbols from MyGene.info for Ensembl IDs.

    Args:
        ensembl_ids: List of Ensembl gene IDs (with or without version numbers)
        species: Species name (default: "human")

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
    import json

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
        species: Species name (default: "human")

    Returns:
        Dictionary mapping Ensembl IDs to gene symbols (guaranteed no NaN)
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


# ─── Data Loading & Preprocessing ────────────────────────────────────────────

def standardize_column_names(df, column_mapping):
    """
    Rename DataFrame columns according to a standard mapping.

    Args:
        df: pandas DataFrame
        column_mapping: Dictionary mapping standard names to file-specific column names

    Returns:
        Tuple of (renamed DataFrame, dictionary of applied renames)
    """
    rename_map = {}
    for std_name, file_col in column_mapping.items():
        if file_col is not None and file_col in df.columns:
            rename_map[file_col] = std_name

    df_renamed = df.rename(columns=rename_map)
    return df_renamed, rename_map


def add_gene_names(df, gene_id_col="gene_id", gene_name_col="gene_name"):
    """
    Add or fill missing gene names using gene ID lookup.

    Args:
        df: pandas DataFrame containing gene IDs
        gene_id_col: Name of the gene ID column
        gene_name_col: Name of the gene name column (will be created if missing)

    Returns:
        pandas DataFrame with gene names added/filled
    """
    df = df.copy()

    # Handle missing gene_name column
    if gene_name_col not in df.columns or df[gene_name_col].isna().all():
        if gene_id_col in df.columns:
            gene_names = lookup_gene_names_enhanced(df[gene_id_col].tolist())
            df[gene_name_col] = df[gene_id_col].map(gene_names)

    # Also fill any remaining NaN gene_name values
    elif gene_name_col in df.columns and df[gene_name_col].isna().any():
        missing_mask = df[gene_name_col].isna()
        if gene_id_col in df.columns and missing_mask.sum() > 0:
            missing_ids = df.loc[missing_mask, gene_id_col].tolist()
            gene_names = lookup_gene_names_enhanced(missing_ids)
            df.loc[missing_mask, gene_name_col] = df.loc[missing_mask, gene_id_col].map(gene_names)

    return df


def add_biotype_groups(df, biotype_col="biotype"):
    """
    Add normalized biotype groups to DataFrame.

    Args:
        df: pandas DataFrame containing biotype column
        biotype_col: Name of the biotype column

    Returns:
        pandas DataFrame with biotype_group column added
    """
    df = df.copy()

    # Handle missing biotype
    if biotype_col not in df.columns:
        df[biotype_col] = "unknown"
        print("  Note: No biotype column — set to 'unknown'")

    # Normalize biotype values
    df["biotype_raw"] = df[biotype_col].copy()
    df["biotype_group"] = normalize_biotype(df[biotype_col])

    return df


def strip_ensembl_version(df, gene_id_col="gene_id", output_col="gene_id_base"):
    """
    Strip version numbers from Ensembl gene IDs for cross-condition matching.

    Args:
        df: pandas DataFrame containing gene IDs
        gene_id_col: Name of the gene ID column
        output_col: Name for the output column with base IDs

    Returns:
        pandas DataFrame with base gene IDs added
    """
    df = df.copy()

    if gene_id_col in df.columns:
        df[output_col] = df[gene_id_col].str.split(".").str[0]

    return df


def load_and_normalize_deseq2(filepath, column_mapping, label=""):
    """
    Load DESeq2 Excel file and normalize columns to standard names.

    Args:
        filepath: Path to Excel file
        column_mapping: Dictionary mapping standard names to file-specific columns
        label: Human-readable label for logging

    Returns:
        pandas DataFrame with standardized columns
    """
    print(f"\n{'='*60}")
    print(f"Loading DESeq2: {label}")
    print(f"  File: {filepath.name if hasattr(filepath, 'name') else filepath}")

    df = pd.read_excel(filepath)
    print(f"  Raw rows: {len(df):,}")
    print(f"  Raw columns: {list(df.columns)}")

    # Rename columns to standard names
    df, rename_map = standardize_column_names(df, column_mapping)
    print(f"  Renamed: {rename_map}")

    # Add gene names
    df = add_gene_names(df)

    # Add biotype groups
    df = add_biotype_groups(df)

    # Strip Ensembl version for cross-condition matching
    df = strip_ensembl_version(df)

    print(f"  Standardized columns: {[c for c in df.columns if c in ['gene_id','gene_name','log2fc','basemean','padj','biotype','biotype_group']]}")

    return df


if __name__ == "__main__":
    # Test mode - demonstrate functionality
    print("Standardization module loaded successfully")
    print(f"Available biotype mappings: {len(BIOTYPE_MAP)} entries")
    print(f"Standard biotype groups: {set(BIOTYPE_MAP.values())}")
