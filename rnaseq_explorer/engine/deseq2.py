"""DESeq2 loading, filtering, annotation, and visualization.

This module handles all DESeq2-related processing:
- File loading with auto-detection of column names
- Gene name enrichment via MyGene.info
- Biotype assignment and re-assignment
- DESeq2 filtering with configurable cutoffs
- RBP annotation
- Volcano, MA, biotype, and other per-condition plots
- Cross-condition biotype comparisons
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns

from rnaseq_explorer.viz.theme import (
    COLOR_DOWN,
    COLOR_NS,
    COLOR_UP,
    add_count_box,
    grid_dims,
)

try:
    import plotly.express as px
except ImportError:
    px = None

try:
    from scipy.stats import fisher_exact
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False

# ---------------------------------------------------------------------------
# Default column mappings
# ---------------------------------------------------------------------------

DEFAULT_DESEQ2_COLS: dict[str, str] = {
    "gene_id": "gene_id",
    "gene_name": "gene_name",
    "log2fc": "log2FoldChange",
    "basemean": "baseMean",
    "padj": "padj",
    "pvalue": "pvalue",
    "biotype": "biotype",
    "stat": "stat",
    "lfcSE": "lfcSE",
}

# Column-name aliases for auto-detection (checked case-insensitively)
_DESEQ2_ALIASES: dict[str, list[str]] = {
    "gene_id": [
        "gene_id", "ensembl_gene_id", "ensembl_geneid", "ensemblgeneid",
        "ensemblgene", "enzemblgeneid", "ensgene",
        "geneid", "gene", "id", "feature_id", "X",
    ],
    "gene_name": [
        "gene_name", "gene_symbol", "symbol", "hgnc_symbol",
        "name", "genename", "external_gene_name", "mgi_symbol",
    ],
    "log2fc": [
        "log2foldchange", "log2fc", "log2_fold_change", "lfc",
        "logfc", "log2ratio", "log2_ratio",
    ],
    "basemean": [
        "basemean", "base_mean", "aveexpr", "mean_expression",
        "meanexpr", "avgexpr", "averageexpression",
    ],
    "padj": [
        "padj", "p.adj", "adjusted_p_value", "adj.p.val",
        "adj_pval", "fdr", "bh", "bonferroni",
    ],
    "pvalue": ["pvalue", "pval", "p.value", "p_value", "p", "rawp"],
    "biotype": [
        "biotype", "gene_biotype", "gene_type", "transcript_biotype",
        "transcript_type",
    ],
    "stat": ["stat", "wald_statistic", "test_stat", "statistic"],
    "lfcSE": ["lfcse", "lfc_se", "std_error", "lfcstderror"],
}

# ---------------------------------------------------------------------------
# Biotype grouping
# ---------------------------------------------------------------------------

BIOTYPE_GROUPS: dict[str, str] = {
    "protein_coding": "Protein Coding",
    "lncrna": "lncRNA", "lincrna": "lncRNA",
    "sense_intronic": "lncRNA", "sense_overlapping": "lncRNA",
    "antisense": "lncRNA", "processed_transcript": "lncRNA",
    "bidirectional_promoter_lncrna": "lncRNA", "macro_lncrna": "lncRNA",
    "non_coding": "lncRNA",
    "pseudogene": "Pseudogene", "processed_pseudogene": "Pseudogene",
    "unprocessed_pseudogene": "Pseudogene",
    "transcribed_unprocessed_pseudogene": "Pseudogene",
    "transcribed_processed_pseudogene": "Pseudogene",
    "transcribed_unitary_pseudogene": "Pseudogene",
    "polymorphic_pseudogene": "Pseudogene", "unitary_pseudogene": "Pseudogene",
    "ig_pseudogene": "Pseudogene", "ig_c_pseudogene": "Pseudogene",
    "ig_v_pseudogene": "Pseudogene", "tr_v_pseudogene": "Pseudogene",
    "tr_j_pseudogene": "Pseudogene",
    "mirna": "Small ncRNA", "snrna": "Small ncRNA", "snorna": "Small ncRNA",
    "misc_ncrna": "Small ncRNA", "rrna": "Small ncRNA", "scrna": "Small ncRNA",
    "scarna": "Small ncRNA", "pirna": "Small ncRNA", "vault_rna": "Small ncRNA",
    "y_rna": "Small ncRNA", "ribozyme": "Small ncRNA", "srp_rna": "Small ncRNA",
    "trna": "Small ncRNA",
}

BIOTYPE_ORDER: list[str] = [
    "Protein Coding", "lncRNA", "Pseudogene", "Small ncRNA", "Other",
]

BIOTYPE_COLORS: dict[str, str] = {
    "Protein Coding": "#4C72B0",
    "lncRNA": "#DD8452",
    "Pseudogene": "#55A868",
    "Small ncRNA": "#C44E52",
    "Other": "#8172B2",
}

# mygene type_of_gene -> Ensembl-style biotype mapping
_MYGENE_BIOTYPE_MAP: dict[str, str] = {
    "protein-coding": "protein_coding",
    "ncrna": "lncrna",
    "pseudo": "pseudogene",
    "snrna": "snrna",
    "snorna": "snorna",
    "rrna": "rrna",
    "trna": "trna",
    "scrna": "scrna",
    "mirna": "mirna",
}

# Runtime cache for gene name lookups
_GENE_NAME_CACHE: dict[str, Optional[str]] = {}


# ---------------------------------------------------------------------------
# File loading utilities
# ---------------------------------------------------------------------------


def load_file(filepath: str | Path, name: str = "file") -> pd.DataFrame:
    """Load CSV/TSV/XLSX file based on extension.

    Handles quoted fields (e.g., ``"ENSG00000123456.15"``) automatically.

    Parameters
    ----------
    filepath : str or Path
        Path to the input file.
    name : str
        Human-readable label for log messages.

    Returns
    -------
    pd.DataFrame
        Loaded data.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Not found: {filepath}")

    ext = path.suffix.lower()
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    elif ext in (".tsv", ".tab"):
        df = pd.read_csv(path, sep="\t", quotechar='"')
    else:
        df = pd.read_csv(path, sep="\t", quotechar='"')
        if len(df.columns) <= 1:
            df = pd.read_csv(path, sep=",", quotechar='"')

    print(f"  Loaded {name}: {df.shape[0]:,} rows x {df.shape[1]} columns")
    return df


def validate_columns(
    df: pd.DataFrame, required_cols: list[str], name: str = "file"
) -> None:
    """Check that expected columns exist in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to check.
    required_cols : list[str]
        List of required column names.
    name : str
        Human-readable label for error messages.

    Raises
    ------
    KeyError
        If any required columns are missing.
    """
    missing = [c for c in required_cols if c and c not in df.columns]
    if missing:
        print(f"  Available columns: {list(df.columns)}")
        raise KeyError(f"Missing columns in {name}: {missing}")


# ---------------------------------------------------------------------------
# Column resolution
# ---------------------------------------------------------------------------


def _resolve_column(
    df: pd.DataFrame,
    col_key: str,
    configured_name: str,
    alias_map: dict[str, list[str]],
    file_label: str = "file",
) -> Optional[str]:
    """Return the actual column name in *df* for a given column key.

    Tries (in order):
    1. Exact match of configured name.
    2. Case-insensitive exact match of configured name.
    3. Case-insensitive exact match of any known alias.
    4. Prefix match on configured name.
    5. Prefix match on aliases.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to search.
    col_key : str
        Logical column key (e.g. 'gene_id').
    configured_name : str
        The user-configured column name.
    alias_map : dict
        Mapping of col_key to list of known aliases.
    file_label : str
        Label for log messages.

    Returns
    -------
    str or None
        Matched column name as it appears in df, or None if not found.
    """
    if configured_name and configured_name in df.columns:
        return configured_name

    lower_cols = {c.lower(): c for c in df.columns}

    if configured_name and configured_name.lower() in lower_cols:
        found = lower_cols[configured_name.lower()]
        print(f"  Column '{configured_name}' not found -- using '{found}' (case-insensitive match)")
        return found

    for alias in alias_map.get(col_key, []):
        if alias.lower() in lower_cols:
            found = lower_cols[alias.lower()]
            print(f"  Column '{configured_name}' not found -- auto-detected '{found}' as {col_key}")
            return found

    if configured_name:
        prefix = configured_name.lower() + "_"
        candidates = [orig for lc, orig in lower_cols.items() if lc.startswith(prefix)]
        if candidates:
            candidates.sort(key=len)
            found = candidates[0]
            if len(candidates) > 1:
                print(
                    f"  WARNING: multiple '{prefix}*' columns in {file_label}: "
                    f"{candidates} -- using shortest '{found}'."
                )
            else:
                print(
                    f"  Column '{configured_name}' not found -- using '{found}' "
                    f"(prefix match for {col_key})"
                )
            return found

    for alias in alias_map.get(col_key, []):
        prefix = alias.lower() + "_"
        candidates = [orig for lc, orig in lower_cols.items() if lc.startswith(prefix)]
        if candidates:
            candidates.sort(key=len)
            found = candidates[0]
            if len(candidates) > 1:
                print(
                    f"  WARNING: multiple '{prefix}*' columns in {file_label}: "
                    f"{candidates} -- using shortest '{found}'."
                )
            else:
                print(
                    f"  Column '{configured_name}' not found -- using '{found}' "
                    f"(prefix match via alias '{alias}' for {col_key})"
                )
            return found

    return None


def _strip_ensembl_version(series: pd.Series) -> pd.Series:
    """Strip version suffixes from Ensembl IDs: ENSG00000123456.12 -> ENSG00000123456.

    Only strips when the value looks like an Ensembl ID (starts with ENS).
    Leaves non-Ensembl values (gene symbols, etc.) unchanged.

    Parameters
    ----------
    series : pd.Series
        Series of gene identifiers.

    Returns
    -------
    pd.Series
        Series with version suffixes stripped from Ensembl IDs.
    """
    def _strip(val):
        s = str(val)
        if s.upper().startswith("ENS") and "." in s:
            return s.rsplit(".", 1)[0]
        return val

    return series.apply(_strip)


def normalize_deseq2_columns(
    df: pd.DataFrame,
    cols: dict[str, str],
    file_label: str = "DESeq2 file",
) -> pd.DataFrame:
    """Rename df columns to match the names in *cols*, then strip Ensembl versions.

    Parameters
    ----------
    df : pd.DataFrame
        Raw DESeq2 data.
    cols : dict
        Column name mapping (same structure as DEFAULT_DESEQ2_COLS).
    file_label : str
        Label for log messages.

    Returns
    -------
    pd.DataFrame
        DataFrame with normalized column names.
    """
    rename_map = {}
    for key, configured in cols.items():
        if not configured:
            continue
        if configured in df.columns:
            continue
        actual = _resolve_column(df, key, configured, _DESEQ2_ALIASES, file_label)
        if actual and actual != configured:
            rename_map[actual] = configured

    if rename_map:
        df = df.rename(columns=rename_map)

    # Strip Ensembl version numbers from gene_id column
    id_col = cols.get("gene_id", "")
    if id_col and id_col in df.columns:
        sample = df[id_col].dropna().astype(str).head(500)
        ens_frac = sample.str.upper().str.startswith("ENS").sum() / max(len(sample), 1)
        if ens_frac > 0.1:
            df = df.copy()
            df[id_col] = _strip_ensembl_version(df[id_col])

    # Detect stat and lfcSE columns (optional -- silent if missing)
    lfcse_col = cols.get("lfcSE", "lfcSE")
    log2fc_col = cols.get("log2fc", "log2FoldChange")

    if lfcse_col and lfcse_col in df.columns and log2fc_col and log2fc_col in df.columns:
        lfc_abs = df[log2fc_col].abs()
        lfcse_vals = df[lfcse_col]
        valid = (lfc_abs > 0) & lfcse_vals.notna()
        if valid.sum() > 0:
            ratio = (lfcse_vals[valid] / lfc_abs[valid]).median()
            if ratio < 0.3:
                print(
                    f"  [INFO] lfcSE/|log2FC| median ratio = {ratio:.3f} -- "
                    f"likely LFC shrinkage applied (apeglm/ashr)"
                )

    return df


def best_gene_key(
    df: pd.DataFrame, cols: dict[str, str]
) -> tuple[str, str]:
    """Return the column name that gives the most reliable unique gene identifier.

    Prefers gene_id (Ensembl) when it looks like Ensembl IDs, because Ensembl IDs
    are stable unique identifiers. Falls back to gene_name.

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 data.
    cols : dict
        Column name mapping.

    Returns
    -------
    tuple
        (col_name, description_string)
    """
    id_col = cols.get("gene_id", "")
    name_col = cols.get("gene_name", "")

    if id_col and id_col in df.columns:
        sample = df[id_col].dropna().astype(str).head(500)
        ens_frac = sample.str.upper().str.startswith("ENS").sum() / max(len(sample), 1)
        if ens_frac > 0.1:
            return id_col, "Ensembl ID"

    return name_col, "gene name"


# ---------------------------------------------------------------------------
# Gene name lookup (MyGene.info REST API)
# ---------------------------------------------------------------------------


def fetch_gene_names(
    ensembl_ids: list[str], species: str = "human"
) -> dict[str, str]:
    """Query MyGene.info to resolve Ensembl IDs to gene symbols.

    Batches up to 1,000 IDs per POST request.

    Parameters
    ----------
    ensembl_ids : list[str]
        Ensembl gene IDs to resolve.
    species : str
        Species name (e.g. 'human', 'mouse').

    Returns
    -------
    dict
        Mapping of ensembl_id -> gene symbol.
    """
    result: dict[str, str] = {}
    ids = [str(i) for i in ensembl_ids if str(i).upper().startswith("ENS")]
    if not ids:
        return result

    batch_size = 1000
    url = "https://mygene.info/v3/query"
    for start in range(0, len(ids), batch_size):
        batch = ids[start: start + batch_size]
        payload = urllib.parse.urlencode({
            "q": ",".join(batch),
            "scopes": "ensembl.gene",
            "fields": "symbol",
            "species": species,
        }).encode()
        req = urllib.request.Request(
            url, data=payload,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read())
                hits = data if isinstance(data, list) else data.get("hits", [])
            for hit in hits:
                if "symbol" in hit and "query" in hit:
                    result[hit["query"]] = hit["symbol"]
        except Exception as exc:
            print(f"  WARNING: MyGene.info lookup failed for batch starting at {start}: {exc}")
    return result


def enrich_with_gene_names(
    df: pd.DataFrame,
    cols: dict[str, str],
    species: str = "human",
    lookup_enabled: bool = True,
    file_label: str = "",
) -> pd.DataFrame:
    """Add a gene_name column from Ensembl ID lookup when it is absent.

    Uses ``_GENE_NAME_CACHE`` to avoid redundant network requests across conditions.

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 data.
    cols : dict
        Column name mapping.
    species : str
        Species for MyGene.info query.
    lookup_enabled : bool
        If False, skip the lookup entirely.
    file_label : str
        Label for log messages.

    Returns
    -------
    pd.DataFrame
        DataFrame with gene_name column added (if missing and lookup succeeds).
    """
    name_col = cols.get("gene_name", "")
    id_col = cols.get("gene_id", "")

    if not lookup_enabled:
        return df
    if name_col and name_col in df.columns:
        return df
    if not id_col or id_col not in df.columns:
        return df

    sample = df[id_col].dropna().astype(str).head(500)
    ens_frac = sample.str.upper().str.startswith("ENS").sum() / max(len(sample), 1)
    if ens_frac <= 0.1:
        return df

    all_ids = df[id_col].dropna().astype(str).unique().tolist()
    to_fetch = [i for i in all_ids if i not in _GENE_NAME_CACHE]

    if to_fetch:
        print(f"  Fetching gene names for {len(to_fetch):,} Ensembl IDs from MyGene.info (species={species})...")
        new_mappings = fetch_gene_names(to_fetch, species=species)
        _GENE_NAME_CACHE.update(new_mappings)
        resolved = sum(1 for i in to_fetch if i in new_mappings)
        print(f"  Resolved {resolved:,} / {len(to_fetch):,} gene names")
        for _id in to_fetch:
            if _id not in _GENE_NAME_CACHE:
                _GENE_NAME_CACHE[_id] = None

    if not name_col:
        return df

    df = df.copy()
    df[name_col] = df[id_col].map(
        lambda x: _GENE_NAME_CACHE.get(str(x)) or f"{str(x)} (no symbol)"
    )
    print(f"  Added gene_name column ('{name_col}') via Ensembl lookup")
    return df


# ---------------------------------------------------------------------------
# Biotype assignment
# ---------------------------------------------------------------------------


def assign_biotype_group(series: pd.Series) -> pd.Series:
    """Map detailed Ensembl biotypes to 5 broad groups using BIOTYPE_GROUPS.

    Case-insensitive.

    Parameters
    ----------
    series : pd.Series
        Series of raw biotype strings.

    Returns
    -------
    pd.Series
        Series of broad biotype group names.
    """
    return series.map(lambda x: BIOTYPE_GROUPS.get(str(x).lower(), "Other"))


def reassign_biotypes_from_mygene(
    df: pd.DataFrame,
    cols: dict[str, str],
    species: str = "human",
    file_label: str = "",
) -> pd.DataFrame:
    """Re-assign biotype_group when all values are 'Other' using mygene type_of_gene.

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 data.
    cols : dict
        Column name mapping.
    species : str
        Species for MyGene.info query.
    file_label : str
        Label for log messages.

    Returns
    -------
    pd.DataFrame
        DataFrame with biotype column updated.
    """
    bio_col = cols.get("biotype", "")
    id_col = cols.get("gene_id", "")
    if not bio_col or bio_col not in df.columns or not id_col or id_col not in df.columns:
        return df

    unique_biotypes = set(df[bio_col].dropna().str.strip().unique())
    if unique_biotypes - {"Other", "other", ""}:
        return df

    print(f"  [INFO] All biotypes are 'Other' in {file_label} -- fetching from MyGene.info...")
    all_ids = df[id_col].dropna().astype(str).unique().tolist()
    ens_ids = [i for i in all_ids if i.upper().startswith("ENS")]
    if not ens_ids:
        return df

    type_map: dict[str, str] = {}
    batch_size = 1000
    url = "https://mygene.info/v3/query"
    for start in range(0, len(ens_ids), batch_size):
        batch = ens_ids[start: start + batch_size]
        payload = urllib.parse.urlencode({
            "q": ",".join(batch),
            "scopes": "ensembl.gene",
            "fields": "type_of_gene",
            "species": species,
        }).encode()
        req = urllib.request.Request(
            url, data=payload,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read())
                hits = data if isinstance(data, list) else data.get("hits", [])
            for hit in hits:
                if "type_of_gene" in hit and "query" in hit:
                    type_map[hit["query"]] = hit["type_of_gene"]
        except Exception as exc:
            print(f"  WARNING: MyGene.info biotype lookup failed: {exc}")

    if not type_map:
        return df

    df = df.copy()

    def _resolve(gene_id):
        raw = type_map.get(str(gene_id), "")
        ensembl_bt = _MYGENE_BIOTYPE_MAP.get(raw.lower(), raw.lower())
        return BIOTYPE_GROUPS.get(ensembl_bt, "Other")

    df[bio_col] = df[id_col].map(_resolve)
    resolved = (df[bio_col] != "Other").sum()
    print(f"  Resolved {resolved:,} / {len(df):,} gene biotypes via MyGene.info")
    return df


def _bh_correction(pvals: list[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (no external deps required).

    Parameters
    ----------
    pvals : list[float]
        Raw p-values.

    Returns
    -------
    np.ndarray
        FDR-adjusted p-values.
    """
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    adjusted = np.array(pvals, dtype=float)[order] * n / (np.arange(1, n + 1))
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    result = np.empty(n)
    result[order] = np.minimum(adjusted, 1.0)
    return result


# ---------------------------------------------------------------------------
# RBP annotation
# ---------------------------------------------------------------------------


def load_rbp_annotations(rbp_file: str | Path) -> dict[str, dict[str, bool]]:
    """Load RBP annotations from an Excel file.

    Supports two formats:
      1) RBP-E-A-C-Complex.xlsx  (Gene Name, TriSNRP, B Complex, ...)
      2) RBP_yael_MW_lists.xlsx  (Gene Symbol, Yael_RBP?, MW_RBP?)

    Parameters
    ----------
    rbp_file : str or Path
        Path to the RBP annotation Excel file.

    Returns
    -------
    dict
        Mapping of uppercase gene_symbol -> dict of annotation booleans.
    """
    path = Path(rbp_file)
    if not path.exists():
        print(f"  WARNING: RBP file not found: {rbp_file}")
        return {}

    df = pd.read_excel(path, sheet_name=0)
    print(f"  Loaded RBP annotations: {df.shape[0]:,} rows x {df.shape[1]} columns")

    gene_col = None
    for candidate in ["Gene Name", "Gene Symbol", "gene_name", "gene_symbol"]:
        if candidate in df.columns:
            gene_col = candidate
            break
    if gene_col is None:
        print(f"  WARNING: Cannot find gene column in RBP file. Columns: {list(df.columns)}")
        return {}

    complex_cols: dict[str, str] = {}
    has_complexes = False
    for raw_col, clean_key in [
        ("TriSNRP", "TriSNRP"),
        ("B Complex", "B_Complex"),
        ("C Complex", "C_Complex"),
        ("U2-Ecomplex", "U2_Ecomplex"),
    ]:
        if raw_col in df.columns:
            complex_cols[raw_col] = clean_key
            has_complexes = True

    annotations: dict[str, dict[str, bool]] = {}
    for _, row in df.iterrows():
        gene = row[gene_col]
        if pd.isna(gene) or str(gene).strip() == "":
            continue
        gene_upper = str(gene).strip().upper()
        entry = {
            "is_RBP_Yael": str(row.get("Yael_RBP?", "")).strip().upper() == "Y",
            "is_RBP_MW": str(row.get("MW_RBP?", "")).strip().upper() == "Y",
        }
        if has_complexes:
            for raw_col, clean_key in complex_cols.items():
                entry[clean_key] = str(row.get(raw_col, "")).strip().upper() == "Y"
        annotations[gene_upper] = entry

    n_yael = sum(1 for v in annotations.values() if v["is_RBP_Yael"])
    n_mw = sum(1 for v in annotations.values() if v["is_RBP_MW"])
    print(f"  RBP annotations: {len(annotations):,} genes (Yael: {n_yael:,}, MW: {n_mw:,})")
    if has_complexes:
        print(f"  Spliceosome complex columns detected: {list(complex_cols.values())}")
    return annotations


def annotate_rbps(
    deg_df: pd.DataFrame,
    rbp_annotations: dict[str, dict[str, bool]],
    gene_col: str = "gene_name",
) -> pd.DataFrame:
    """Add RBP annotation columns to a DEG DataFrame.

    Adds: is_RBP, is_RBP_Yael, is_RBP_MW, and spliceosome complex columns
    if available. Matching is case-insensitive.

    Parameters
    ----------
    deg_df : pd.DataFrame
        Filtered DEG data.
    rbp_annotations : dict
        Output of load_rbp_annotations().
    gene_col : str
        Column containing gene names/symbols.

    Returns
    -------
    pd.DataFrame
        DataFrame with RBP annotation columns added.
    """
    if not rbp_annotations or gene_col not in deg_df.columns:
        return deg_df

    df = deg_df.copy()
    gene_upper = df[gene_col].astype(str).str.strip().str.upper()

    sample_entry = next(iter(rbp_annotations.values()))
    all_keys = list(sample_entry.keys())

    for key in all_keys:
        df[key] = gene_upper.map(
            lambda g, k=key: rbp_annotations.get(g, {}).get(k, False)
        )

    df["is_RBP"] = df["is_RBP_Yael"] | df["is_RBP_MW"]

    cols = [c for c in df.columns if c not in ["is_RBP"] + all_keys]
    df = df[cols + ["is_RBP"] + all_keys]

    n_rbp = int(df["is_RBP"].sum())
    print(f"  RBP-annotated: {n_rbp:,} / {len(df):,} genes are RBPs")
    return df


# ---------------------------------------------------------------------------
# DESeq2 filtering
# ---------------------------------------------------------------------------


def filter_deseq2(
    df: pd.DataFrame,
    cols: dict[str, str],
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    padj_cutoff: float = 0.05,
    biotype_filter: Optional[str] = None,
    label: str = "All",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply DESeq2 cutoffs and return (biotype-subset, significance-filtered) DataFrames.

    Parameters
    ----------
    df : pd.DataFrame
        Full DESeq2 results.
    cols : dict
        Column name mapping.
    log2fc_cutoff : float
        Absolute log2 fold change threshold.
    basemean_cutoff : float
        Minimum baseMean expression.
    padj_cutoff : float
        Adjusted p-value threshold.
    biotype_filter : str or None
        'protein_coding', 'non_protein_coding', or None for all genes.
    label : str
        Label for log messages.

    Returns
    -------
    tuple
        (biotype_subset_df, significance_filtered_df)
    """
    n_before = len(df)
    df = df.dropna(subset=[cols["padj"], cols["log2fc"], cols["basemean"]])
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(
            f"  NOTE: {n_dropped:,} rows dropped (NA in padj/log2FC/baseMean -- "
            f"normal DESeq2 behaviour for low-count/outlier genes)"
        )

    if biotype_filter in ("protein_coding", "non_protein_coding"):
        _bio_norm = df[cols["biotype"]].fillna("").str.lower().str.replace(" ", "_")
        if biotype_filter == "protein_coding":
            df = df[_bio_norm == "protein_coding"]
        else:
            df = df[_bio_norm != "protein_coding"]

    mask = (
        (df[cols["padj"]] < padj_cutoff)
        & (df[cols["log2fc"]].abs() >= log2fc_cutoff)
        & (df[cols["basemean"]] >= basemean_cutoff)
    )

    filtered = df[mask].copy()
    filtered["direction"] = np.where(filtered[cols["log2fc"]] > 0, "up", "down")

    print(f"\n-- DESeq2 Filtering [{label}] --")
    print(f"  Input rows:  {len(df):,}")
    print(
        f"  |log2FoldChange| >= {log2fc_cutoff}, baseMean >= {basemean_cutoff}, padj < {padj_cutoff}"
    )
    if biotype_filter:
        print(f"  Biotype filter: {biotype_filter}")
    print(
        f"  -> Filtered: {len(filtered):,} genes "
        f"({filtered['direction'].value_counts().to_dict()})"
    )

    return df, filtered


# ---------------------------------------------------------------------------
# Gene set extraction
# ---------------------------------------------------------------------------


def extract_gene_sets(
    condition_results: dict[str, dict],
    cols: dict[str, str],
) -> dict[str, dict[str, set[str]]]:
    """Extract significant gene identifier sets from per-condition DESeq2 results.

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    cols : dict
        Column name mapping.

    Returns
    -------
    dict
        Mapping of condition_name -> {"all": set, "up": set, "down": set}.
    """
    gene_sets: dict[str, dict[str, set[str]]] = {}

    for cond_name, data in condition_results.items():
        deg_df = data.get("deseq2_filtered", {}).get("all_genes", pd.DataFrame())
        if len(deg_df) == 0:
            gene_sets[cond_name] = {"all": set(), "up": set(), "down": set()}
            continue

        key_col, _ = best_gene_key(deg_df, cols)
        if key_col not in deg_df.columns:
            gene_sets[cond_name] = {"all": set(), "up": set(), "down": set()}
            continue

        all_genes = set(deg_df[key_col].dropna().astype(str).unique())
        up_genes = set()
        down_genes = set()
        if "direction" in deg_df.columns:
            up_genes = set(
                deg_df.loc[deg_df["direction"] == "up", key_col]
                .dropna().astype(str).unique()
            )
            down_genes = set(
                deg_df.loc[deg_df["direction"] == "down", key_col]
                .dropna().astype(str).unique()
            )

        gene_sets[cond_name] = {"all": all_genes, "up": up_genes, "down": down_genes}

    return gene_sets


# ---------------------------------------------------------------------------
# Per-condition DESeq2 visualization
# ---------------------------------------------------------------------------


def volcano_plot(
    df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate volcano plot from full DESeq2 results."""
    outdir = Path(outdir)
    data = df.dropna(subset=[cols["padj"], cols["log2fc"]]).copy()
    data["-log10padj"] = -np.log10(data[cols["padj"]].clip(lower=1e-300))

    basemean_ok = (data[cols["basemean"]] >= basemean_cutoff) if cols["basemean"] in data.columns else True
    conditions = [
        (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]] >= log2fc_cutoff) & basemean_ok,
        (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]] <= -log2fc_cutoff) & basemean_ok,
    ]
    data["status"] = np.select(conditions, ["Up", "Down"], default="NS")

    color_map = {"Up": COLOR_UP, "Down": COLOR_DOWN, "NS": COLOR_NS}
    fig, ax = plt.subplots(figsize=(8, 6))
    for status in ["NS", "Down", "Up"]:
        subset = data[data["status"] == status]
        lbl = "NS" if status == "NS" else f"{status} ({len(subset):,})"
        ax.scatter(subset[cols["log2fc"]], subset["-log10padj"],
                   c=color_map[status], s=8, alpha=0.5, edgecolors="none",
                   label=lbl, rasterized=True)

    ax.axhline(-np.log10(padj_cutoff), color="grey", ls="--", lw=0.8)
    ax.axvline(log2fc_cutoff, color="grey", ls="--", lw=0.8)
    ax.axvline(-log2fc_cutoff, color="grey", ls="--", lw=0.8)

    n_up = int((data["status"] == "Up").sum())
    n_down = int((data["status"] == "Down").sum())
    add_count_box(ax, n_up, n_down, n_up + n_down, position="lower left")

    ax.set_xlabel("log$_2$ Fold Change")
    ax.set_ylabel("-log$_{10}$ (adjusted p-value)")
    ax.set_title(f"Volcano Plot - DESeq2 ({label})")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=10, markerscale=2)

    outpath = outdir / f"volcano_plot{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def ma_plot(
    df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate MA plot (baseMean vs log2FC)."""
    outdir = Path(outdir)
    data = df.dropna(subset=[cols["padj"], cols["log2fc"], cols["basemean"]]).copy()

    basemean_ok = data[cols["basemean"]] >= basemean_cutoff
    sig = (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]].abs() >= log2fc_cutoff) & basemean_ok
    data["significant"] = np.where(
        sig & (data[cols["log2fc"]] > 0), "Up",
        np.where(sig & (data[cols["log2fc"]] < 0), "Down", "NS")
    )

    color_map = {"Up": COLOR_UP, "Down": COLOR_DOWN, "NS": COLOR_NS}
    fig, ax = plt.subplots(figsize=(8, 6))
    for status in ["NS", "Down", "Up"]:
        subset = data[data["significant"] == status]
        lbl = "NS" if status == "NS" else f"{status} ({len(subset):,})"
        ax.scatter(np.log10(subset[cols["basemean"]].clip(lower=0.1)),
                   subset[cols["log2fc"]], c=color_map[status], s=8, alpha=0.5,
                   edgecolors="none", label=lbl, rasterized=True)

    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(log2fc_cutoff, color="grey", ls="--", lw=0.8)
    ax.axhline(-log2fc_cutoff, color="grey", ls="--", lw=0.8)

    n_up = int((data["significant"] == "Up").sum())
    n_down = int((data["significant"] == "Down").sum())
    add_count_box(ax, n_up, n_down, n_up + n_down, position="lower left")

    ax.set_xlabel("log$_{10}$ (baseMean)")
    ax.set_ylabel("log$_2$ Fold Change")
    ax.set_title(f"MA Plot - DESeq2 ({label})")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=10, markerscale=2)

    outpath = outdir / f"ma_plot{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def volcano_plot_interactive(
    df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    label: str = "All",
    suffix: str = "",
) -> None:
    """Generate interactive volcano plot with hover tooltips (HTML output)."""
    if px is None:
        return
    outdir = Path(outdir)
    data = df.dropna(subset=[cols["padj"], cols["log2fc"]]).copy()
    data["-log10padj"] = -np.log10(data[cols["padj"]].clip(lower=1e-300))

    basemean_ok = (data[cols["basemean"]] >= basemean_cutoff) if cols["basemean"] in data.columns else True
    conditions = [
        (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]] >= log2fc_cutoff) & basemean_ok,
        (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]] <= -log2fc_cutoff) & basemean_ok,
    ]
    data["Status"] = np.select(conditions, ["Up", "Down"], default="NS")

    n_up = int((data["Status"] == "Up").sum())
    n_down = int((data["Status"] == "Down").sum())
    color_map = {"NS": COLOR_NS, "Down": COLOR_DOWN, "Up": COLOR_UP}

    _hover_name = cols["gene_name"] if cols["gene_name"] in data.columns else cols["gene_id"]
    fig = px.scatter(
        data.sort_values("Status", key=lambda s: s.map({"NS": 0, "Down": 1, "Up": 2})),
        x=cols["log2fc"], y="-log10padj", color="Status",
        color_discrete_map=color_map,
        category_orders={"Status": ["NS", "Down", "Up"]},
        hover_name=_hover_name,
        hover_data={cols["log2fc"]: ":.3f", cols["padj"]: ":.2e",
                    cols["basemean"]: ":.1f", "-log10padj": ":.2f", "Status": False},
        opacity=0.5, title=f"Volcano Plot - DESeq2 ({label})",
    )
    fig.update_traces(marker=dict(size=5))
    fig.add_hline(y=-np.log10(padj_cutoff), line_dash="dash", line_color="grey", line_width=0.8)
    fig.add_vline(x=log2fc_cutoff, line_dash="dash", line_color="grey", line_width=0.8)
    fig.add_vline(x=-log2fc_cutoff, line_dash="dash", line_color="grey", line_width=0.8)
    fig.add_annotation(
        text=f"Up: {n_up:,}<br>Down: {n_down:,}<br>Total: {n_up + n_down:,}",
        xref="paper", yref="paper", x=0.02, y=0.98, showarrow=False,
        bgcolor="rgba(255,255,255,0.85)", bordercolor="grey", borderwidth=1,
        font=dict(size=11), align="left", xanchor="left", yanchor="top",
    )
    fig.update_layout(xaxis_title="log\u2082 Fold Change",
                      yaxis_title="-log\u2081\u2080 (adjusted p-value)",
                      hovermode="closest", template="plotly_white", width=900, height=650)

    outpath = outdir / f"volcano_plot{suffix}_interactive.html"
    fig.write_html(str(outpath))
    print(f"  Saved: {outpath}")


def ma_plot_interactive(
    df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    label: str = "All",
    suffix: str = "",
) -> None:
    """Generate interactive MA plot with hover tooltips (HTML output)."""
    if px is None:
        return
    outdir = Path(outdir)
    data = df.dropna(subset=[cols["padj"], cols["log2fc"], cols["basemean"]]).copy()
    data["log10_basemean"] = np.log10(data[cols["basemean"]].clip(lower=0.1))

    basemean_ok = data[cols["basemean"]] >= basemean_cutoff
    sig = (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]].abs() >= log2fc_cutoff) & basemean_ok
    data["Status"] = np.where(
        sig & (data[cols["log2fc"]] > 0), "Up",
        np.where(sig & (data[cols["log2fc"]] < 0), "Down", "NS"),
    )

    n_up = int((data["Status"] == "Up").sum())
    n_down = int((data["Status"] == "Down").sum())
    color_map = {"NS": COLOR_NS, "Down": COLOR_DOWN, "Up": COLOR_UP}

    _hover_name = cols["gene_name"] if cols["gene_name"] in data.columns else cols["gene_id"]
    fig = px.scatter(
        data.sort_values("Status", key=lambda s: s.map({"NS": 0, "Down": 1, "Up": 2})),
        x="log10_basemean", y=cols["log2fc"], color="Status",
        color_discrete_map=color_map,
        category_orders={"Status": ["NS", "Down", "Up"]},
        hover_name=_hover_name,
        hover_data={cols["log2fc"]: ":.3f", cols["padj"]: ":.2e",
                    cols["basemean"]: ":.1f", "log10_basemean": False, "Status": False},
        opacity=0.5, title=f"MA Plot - DESeq2 ({label})",
    )
    fig.update_traces(marker=dict(size=5))
    fig.add_hline(y=0, line_color="black", line_width=0.8)
    fig.add_hline(y=log2fc_cutoff, line_dash="dash", line_color="grey", line_width=0.8)
    fig.add_hline(y=-log2fc_cutoff, line_dash="dash", line_color="grey", line_width=0.8)
    fig.add_annotation(
        text=f"Up: {n_up:,}<br>Down: {n_down:,}<br>Total: {n_up + n_down:,}",
        xref="paper", yref="paper", x=0.02, y=0.98, showarrow=False,
        bgcolor="rgba(255,255,255,0.85)", bordercolor="grey", borderwidth=1,
        font=dict(size=11), align="left", xanchor="left", yanchor="top",
    )
    fig.update_layout(xaxis_title="log\u2081\u2080 (baseMean)",
                      yaxis_title="log\u2082 Fold Change",
                      hovermode="closest", template="plotly_white", width=900, height=650)

    outpath = outdir / f"ma_plot{suffix}_interactive.html"
    fig.write_html(str(outpath))
    print(f"  Saved: {outpath}")


def biotype_chart(
    filtered_df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate biotype distribution bar chart from filtered DESeq2 results."""
    outdir = Path(outdir)
    if cols["biotype"] not in filtered_df.columns:
        print("  Skipping biotype chart (no biotype column)")
        return

    counts = filtered_df[cols["biotype"]].value_counts()
    if len(counts) == 0:
        return

    if len(counts) > 10:
        top = counts.head(9)
        other = pd.Series({"Other": counts.iloc[9:].sum()})
        counts = pd.concat([top, other])

    palette = sns.color_palette("Set2", n_colors=len(counts))
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax1 = axes[0]
    bars = ax1.barh(counts.index[::-1], counts.values[::-1], color=palette[::-1])
    ax1.set_xlabel("Number of DE Genes")
    ax1.set_title(f"Biotype Distribution ({label})")
    for bar, val in zip(bars, counts.values[::-1]):
        ax1.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
                 f"{val}", va="center", fontsize=9)

    ax2 = axes[1]
    wedges, texts, autotexts = ax2.pie(
        counts.values, labels=counts.index, autopct="%1.1f%%",
        colors=palette, startangle=90, pctdistance=0.8)
    for t in autotexts:
        t.set_fontsize(8)
    ax2.set_title(f"Biotype Proportions ({label})")

    plt.tight_layout()
    outpath = outdir / f"biotype_distribution{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def biotype_direction_chart(
    filtered_df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Diverging horizontal bar chart: Up/Down DE gene counts split by biotype group."""
    outdir = Path(outdir)
    if cols["biotype"] not in filtered_df.columns:
        print("  Skipping biotype direction chart (no biotype column)")
        return
    if len(filtered_df) == 0:
        return

    df = filtered_df.copy()
    df["_group"] = assign_biotype_group(df[cols["biotype"]])

    up_counts = df[df["direction"] == "up"]["_group"].value_counts()
    dn_counts = df[df["direction"] == "down"]["_group"].value_counts()

    groups = [g for g in BIOTYPE_ORDER if g in up_counts.index or g in dn_counts.index]
    if not groups:
        return

    up_vals = [up_counts.get(g, 0) for g in groups]
    dn_vals = [-dn_counts.get(g, 0) for g in groups]

    fig, ax = plt.subplots(figsize=(9, max(3, len(groups) * 1.0 + 1.5)))
    y = range(len(groups))
    ax.barh(list(y), up_vals, color=COLOR_UP, label="Up-regulated", alpha=0.85)
    ax.barh(list(y), dn_vals, color=COLOR_DOWN, label="Down-regulated", alpha=0.85)
    ax.axvline(0, color="black", linewidth=0.8)

    for i, (u, d) in enumerate(zip(up_vals, dn_vals)):
        if u > 0:
            ax.text(u + 0.3, i, str(u), va="center", fontsize=9, color=COLOR_UP)
        if d < 0:
            ax.text(d - 0.3, i, str(-d), va="center", ha="right", fontsize=9, color=COLOR_DOWN)

    ax.set_yticks(list(y))
    ax.set_yticklabels(groups)
    ax.set_xlabel("Number of DE Genes")
    ax.set_title(f"DE Genes by Biotype & Direction ({label})")
    ax.legend(loc="lower right", fontsize=9)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: str(int(abs(x)))))
    plt.tight_layout()

    outpath = outdir / f"biotype_direction_chart{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def biotype_enrichment_test(
    filtered_df: pd.DataFrame,
    all_df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Fisher's exact test: is each biotype group enriched among DE genes vs background?"""
    if not _SCIPY_AVAILABLE:
        print("  Skipping biotype enrichment (scipy not installed)")
        return
    outdir = Path(outdir)
    if cols["biotype"] not in all_df.columns:
        print("  Skipping biotype enrichment (no biotype column)")
        return

    bg = all_df.copy()
    bg["_group"] = assign_biotype_group(bg[cols["biotype"]])
    de = filtered_df.copy()
    if cols["biotype"] in de.columns:
        de["_group"] = assign_biotype_group(de[cols["biotype"]])
    else:
        return

    n_bg = len(bg)
    n_de = len(de)
    if n_de == 0:
        return

    results = []
    for grp in BIOTYPE_ORDER:
        a = int((de["_group"] == grp).sum())
        b = n_de - a
        c = int((bg["_group"] == grp).sum()) - a
        d = n_bg - n_de - c
        if a + c == 0:
            continue
        c = max(c, 0)
        d = max(d, 0)
        try:
            odds, pval = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        except Exception:
            continue
        if odds == 0:
            log2_or = -10.0
        elif not np.isfinite(odds):
            log2_or = 10.0
        else:
            log2_or = float(np.clip(np.log2(odds), -10, 10))
        results.append({"group": grp, "log2OR": log2_or, "pval": pval,
                        "n_de": a, "n_group": a + c})

    if not results:
        return

    res_df = pd.DataFrame(results)
    res_df["fdr"] = _bh_correction(res_df["pval"].tolist())
    res_df = res_df.sort_values("log2OR")

    fig, ax = plt.subplots(figsize=(9, max(3, len(res_df) * 0.9 + 1.5)))
    colors = [("#2ecc71" if f < 0.05 else COLOR_NS) for f in res_df["fdr"]]
    ax.hlines(range(len(res_df)), 0, res_df["log2OR"], color="grey", linewidth=1.5, zorder=1)
    ax.scatter(res_df["log2OR"], range(len(res_df)), color=colors, s=80, zorder=2)
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")

    x_max = max(abs(res_df["log2OR"].max()), abs(res_df["log2OR"].min()), 1.0)
    ax.set_xlim(-x_max * 1.4, x_max * 1.8)
    for i, row in enumerate(res_df.itertuples()):
        ax.text(x_max * 1.05, i, f" {row.n_de}/{row.n_group}", va="center", fontsize=8, ha="left")
    ax.text(x_max * 1.05, len(res_df) - 0.5, "DE/Total", va="bottom", fontsize=7,
            ha="left", color="grey")

    ax.set_yticks(range(len(res_df)))
    ax.set_yticklabels(res_df["group"].tolist())
    ax.set_xlabel("log\u2082(Odds Ratio)")
    ax.set_title(f"Biotype Enrichment Among DE Genes ({label})\nGreen = FDR < 0.05")
    plt.tight_layout()

    outpath = outdir / f"biotype_enrichment{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def biotype_volcano(
    all_df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Volcano plot colored by biotype group instead of direction."""
    outdir = Path(outdir)
    if cols["biotype"] not in all_df.columns:
        print("  Skipping biotype volcano (no biotype column)")
        return

    df = all_df.copy().dropna(subset=[cols["log2fc"], cols["padj"]])
    df["_group"] = assign_biotype_group(df[cols["biotype"]])
    df["_neg_log10p"] = -np.log10(df[cols["padj"]].clip(lower=1e-300))
    _bm_ok = (df[cols["basemean"]] >= basemean_cutoff) if cols["basemean"] in df.columns else True
    df["_sig"] = (df[cols["padj"]] < padj_cutoff) & (df[cols["log2fc"]].abs() >= log2fc_cutoff) & _bm_ok

    fig, ax = plt.subplots(figsize=(9, 7))

    ns = df[~df["_sig"]]
    ax.scatter(ns[cols["log2fc"]], ns["_neg_log10p"], color=COLOR_NS,
               alpha=0.25, s=4, rasterized=True, label=None)

    sig = df[df["_sig"]]
    legend_handles = []
    for grp in BIOTYPE_ORDER:
        grp_sig = sig[sig["_group"] == grp]
        if len(grp_sig) == 0:
            continue
        color = BIOTYPE_COLORS.get(grp, "#999999")
        ax.scatter(grp_sig[cols["log2fc"]], grp_sig["_neg_log10p"],
                   color=color, alpha=0.75, s=10, rasterized=True)
        legend_handles.append(mpatches.Patch(color=color, label=f"{grp} (n={len(grp_sig)})"))

    ax.axhline(-np.log10(padj_cutoff), color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(log2fc_cutoff, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(-log2fc_cutoff, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log\u2082 Fold Change")
    ax.set_ylabel("-log\u2081\u2080(adjusted p-value)")
    ax.set_title(f"Volcano Plot by Biotype ({label})")
    if legend_handles:
        ax.legend(handles=legend_handles, fontsize=8, loc="upper left")
    plt.tight_layout()

    outpath = outdir / f"volcano_biotype{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def ecdf_log2fc_by_biotype(
    all_df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """ECDF of |log2FC| per biotype group."""
    outdir = Path(outdir)
    if cols["biotype"] not in all_df.columns:
        print("  Skipping ECDF biotype plot (no biotype column)")
        return

    df = all_df.copy().dropna(subset=[cols["log2fc"]])
    if cols["basemean"] in df.columns:
        df = df[df[cols["basemean"]] >= basemean_cutoff]
    df["_group"] = assign_biotype_group(df[cols["biotype"]])
    df["_abs_lfc"] = df[cols["log2fc"]].abs()

    min_genes = 10
    plotted = 0
    fig, ax = plt.subplots(figsize=(8, 5))

    for grp in BIOTYPE_ORDER:
        sub = df[df["_group"] == grp]["_abs_lfc"].dropna().sort_values()
        if len(sub) < min_genes:
            continue
        ecdf_y = np.arange(1, len(sub) + 1) / len(sub)
        color = BIOTYPE_COLORS.get(grp, "#999999")
        ax.plot(sub.values, ecdf_y, label=f"{grp} (n={len(sub)})", color=color, linewidth=2)
        plotted += 1

    if plotted < 2:
        plt.close(fig)
        print("  Skipping ECDF biotype (fewer than 2 groups with >=10 genes)")
        return

    ax.axvline(log2fc_cutoff, color="grey", linestyle="--", linewidth=0.8,
               label=f"cutoff ({log2fc_cutoff})")
    ax.set_xlabel("|log\u2082 Fold Change|")
    ax.set_ylabel("Cumulative Fraction")
    ax.set_title(f"ECDF of |log\u2082FC| by Biotype Group ({label})")
    ax.legend(fontsize=8)
    ax.set_xlim(left=0)
    plt.tight_layout()

    outpath = outdir / f"ecdf_log2fc_biotype{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def pvalue_histogram(
    df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    padj_cutoff: float = 0.05,
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
) -> None:
    """Raw p-value distribution histogram -- QC diagnostic."""
    outdir = Path(outdir)
    pval_col = cols["pvalue"]
    if pval_col not in df.columns:
        print(f"  Skipping p-value histogram: column '{pval_col}' not found")
        return

    data = df[pval_col].dropna()
    if len(data) == 0:
        print("  Skipping p-value histogram: no data")
        return

    fig, ax = plt.subplots(figsize=(7, 5))
    _, edges, patches = ax.hist(data, bins=50, range=(0, 1), edgecolor="white", lw=0.3)
    for patch, left in zip(patches, edges[:-1]):
        patch.set_facecolor(COLOR_UP if left < padj_cutoff else COLOR_NS)

    ax.axvline(padj_cutoff, color="black", ls="--", lw=1.0, label=f"padj cutoff ({padj_cutoff})")
    n_sig = int((data < padj_cutoff).sum())
    ax.text(0.97, 0.97, f"p < {padj_cutoff}: {n_sig:,}\nTotal: {len(data):,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.85))
    ax.set_xlabel("Raw p-value")
    ax.set_ylabel("Gene count")
    ax.set_title(f"P-value Distribution -- {label}\n(spike near 0 = true DE signal; flat = no enrichment)")
    ax.legend(fontsize=9)

    outpath = outdir / f"pvalue_histogram{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def top_genes_lollipop(
    filtered_df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    log2fc_cutoff: float = 1.0,
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    top_n: int = 20,
) -> None:
    """Horizontal lollipop chart of top N up + top N down genes by log2FC."""
    outdir = Path(outdir)
    if len(filtered_df) == 0:
        print(f"  Skipping lollipop: no significant genes [{label}]")
        return

    fc_col = cols["log2fc"]
    name_col = cols["gene_name"]
    if name_col not in filtered_df.columns:
        print(f"  Skipping lollipop: gene_name column '{name_col}' not found")
        return

    df = filtered_df.dropna(subset=[fc_col, name_col]).copy()
    up = df[df["direction"] == "up"].nlargest(top_n, fc_col)
    down = df[df["direction"] == "down"].nsmallest(top_n, fc_col)
    plot_df = pd.concat([down, up], ignore_index=True)

    fig_h = max(6, len(plot_df) * 0.35)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    colors = [COLOR_UP if d == "up" else COLOR_DOWN for d in plot_df["direction"]]
    y_pos = range(len(plot_df))

    for i, (fc, col) in enumerate(zip(plot_df[fc_col], colors)):
        ax.hlines(i, 0, fc, color=col, lw=1.8, alpha=0.8)
        ax.plot(fc, i, "o", color=col, ms=7, zorder=3)

    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(log2fc_cutoff, color="grey", ls="--", lw=0.7)
    ax.axvline(-log2fc_cutoff, color="grey", ls="--", lw=0.7)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(plot_df[name_col].tolist(), fontsize=max(6, 9 - len(plot_df) // 15))
    ax.set_xlabel("log$_2$ Fold Change")
    ax.set_title(f"Top DE Genes -- {label}\n(top {top_n} up + {top_n} down by |log2FC|)")
    ax.legend(handles=[mpatches.Patch(facecolor=COLOR_UP, label="Up"),
                       mpatches.Patch(facecolor=COLOR_DOWN, label="Down")],
              loc="lower right", fontsize=9)

    outpath = outdir / f"top_genes_lollipop{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def expression_rank_plot(
    df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Waterfall rank plot: all genes sorted by log2FC, significant genes highlighted."""
    outdir = Path(outdir)
    data = df.dropna(subset=[cols["padj"], cols["log2fc"], cols["basemean"]]).copy()
    if len(data) == 0:
        print(f"  Skipping rank plot: no data [{label}]")
        return

    data = data.sort_values(cols["log2fc"]).reset_index(drop=True)
    rank = np.arange(len(data))

    _bm_ok = (data[cols["basemean"]] >= basemean_cutoff) if cols["basemean"] in data.columns else True
    conds = [
        (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]] >= log2fc_cutoff) & _bm_ok,
        (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]] <= -log2fc_cutoff) & _bm_ok,
    ]
    data["status"] = np.select(conds, ["Up", "Down"], default="NS")

    fig, ax = plt.subplots(figsize=(9, 5))
    for status, color, size, alpha, z in [
        ("NS", COLOR_NS, 3, 0.25, 1),
        ("Down", COLOR_DOWN, 6, 0.75, 2),
        ("Up", COLOR_UP, 6, 0.75, 2),
    ]:
        mask = data["status"] == status
        lbl = "NS" if status == "NS" else f"{status} ({mask.sum():,})"
        ax.scatter(rank[mask], data.loc[mask, cols["log2fc"]],
                   c=color, s=size, alpha=alpha, edgecolors="none", rasterized=True, zorder=z,
                   label=lbl)

    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(log2fc_cutoff, color="grey", ls="--", lw=0.7)
    ax.axhline(-log2fc_cutoff, color="grey", ls="--", lw=0.7)
    n_up = int((data["status"] == "Up").sum())
    n_down = int((data["status"] == "Down").sum())
    add_count_box(ax, n_up, n_down, n_up + n_down, position="lower left")
    ax.set_xlabel("Gene Rank (sorted by log$_2$ FC)")
    ax.set_ylabel("log$_2$ Fold Change")
    ax.set_title(f"Expression Rank Plot -- {label}")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9, markerscale=2)

    outpath = outdir / f"expression_rank_plot{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def volcano_plot_labeled(
    df: pd.DataFrame,
    outdir: str | Path,
    cols: dict[str, str],
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    basemean_cutoff: float = 10.0,
    label: str = "All",
    suffix: str = "",
    fig_format: str = "png",
    fig_dpi: int = 300,
    genes_of_interest: list[str] | None = None,
    top_n_label: int = 10,
) -> None:
    """Volcano plot with gene name labels on top DE genes and genes of interest."""
    outdir = Path(outdir)
    data = df.dropna(subset=[cols["padj"], cols["log2fc"]]).copy()
    data["-log10padj"] = -np.log10(data[cols["padj"]].clip(lower=1e-300))

    basemean_ok = (data[cols["basemean"]] >= basemean_cutoff) if cols["basemean"] in data.columns else True
    conditions = [
        (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]] >= log2fc_cutoff) & basemean_ok,
        (data[cols["padj"]] < padj_cutoff) & (data[cols["log2fc"]] <= -log2fc_cutoff) & basemean_ok,
    ]
    data["status"] = np.select(conditions, ["Up", "Down"], default="NS")

    color_map = {"Up": COLOR_UP, "Down": COLOR_DOWN, "NS": COLOR_NS}
    fig, ax = plt.subplots(figsize=(10, 8))
    for status in ["NS", "Down", "Up"]:
        subset = data[data["status"] == status]
        lbl_str = "NS" if status == "NS" else f"{status} ({len(subset):,})"
        ax.scatter(subset[cols["log2fc"]], subset["-log10padj"],
                   c=color_map[status], s=8, alpha=0.5, edgecolors="none",
                   label=lbl_str, rasterized=True)

    ax.axhline(-np.log10(padj_cutoff), color="grey", ls="--", lw=0.8)
    ax.axvline(log2fc_cutoff, color="grey", ls="--", lw=0.8)
    ax.axvline(-log2fc_cutoff, color="grey", ls="--", lw=0.8)

    # Label top genes
    name_col = cols.get("gene_name", "")
    if name_col and name_col in data.columns:
        sig_data = data[data["status"].isin(["Up", "Down"])].copy()
        sig_data["_abs_lfc"] = sig_data[cols["log2fc"]].abs()
        top_genes = sig_data.nlargest(top_n_label, "_abs_lfc")

        texts = []
        try:
            from adjustText import adjust_text
            _ADJUST_TEXT = True
        except ImportError:
            _ADJUST_TEXT = False

        for _, row in top_genes.iterrows():
            gene_name = row[name_col]
            if pd.isna(gene_name):
                continue
            txt = ax.annotate(
                str(gene_name),
                (row[cols["log2fc"]], row["-log10padj"]),
                fontsize=7, alpha=0.9,
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.5) if _ADJUST_TEXT else None,
            )
            texts.append(txt)

        # Label genes of interest
        if genes_of_interest:
            goi_upper = {g.upper() for g in genes_of_interest}
            goi_data = data[data[name_col].astype(str).str.upper().isin(goi_upper)]
            for _, row in goi_data.iterrows():
                gene_name = row[name_col]
                if pd.isna(gene_name):
                    continue
                txt = ax.annotate(
                    str(gene_name),
                    (row[cols["log2fc"]], row["-log10padj"]),
                    fontsize=8, fontweight="bold", color="#D55E00",
                    arrowprops=dict(arrowstyle="-", color="#D55E00", lw=0.8) if _ADJUST_TEXT else None,
                )
                texts.append(txt)

        if _ADJUST_TEXT and texts:
            try:
                adjust_text(texts, ax=ax)
            except Exception:
                pass

    n_up = int((data["status"] == "Up").sum())
    n_down = int((data["status"] == "Down").sum())
    add_count_box(ax, n_up, n_down, n_up + n_down, position="lower left")

    ax.set_xlabel("log$_2$ Fold Change")
    ax.set_ylabel("-log$_{10}$ (adjusted p-value)")
    ax.set_title(f"Volcano Plot (Labeled) - DESeq2 ({label})")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=10, markerscale=2)

    outpath = outdir / f"volcano_plot_labeled{suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# RBP visualization
# ---------------------------------------------------------------------------


def rbp_heatmap(
    condition_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str],
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Heatmap of log2FC for RBP genes across conditions."""
    outdir = Path(outdir)
    fc_col = cols["log2fc"]
    name_col = cols.get("gene_name", "")

    fc_data = {}
    for name, data in condition_results.items():
        lbl = condition_labels[name]
        filt = data["deseq2_filtered"].get("all_genes", pd.DataFrame())
        if "is_RBP" not in filt.columns or len(filt) == 0:
            continue
        rbps = filt[filt["is_RBP"] == True]  # noqa: E712
        if len(rbps) == 0 or name_col not in rbps.columns:
            continue
        series = rbps.set_index(name_col)[fc_col]
        series = series[~series.index.duplicated(keep="first")]
        fc_data[lbl] = series

    if len(fc_data) < 1:
        print("  No RBP data for heatmap")
        return

    matrix = pd.DataFrame(fc_data).dropna(how="all")
    if len(matrix) == 0:
        print("  No overlapping RBPs for heatmap")
        return

    if len(matrix) > 80:
        matrix["_max"] = matrix.abs().max(axis=1)
        matrix = matrix.nlargest(80, "_max").drop(columns="_max")

    g = sns.clustermap(matrix.fillna(0), cmap="RdBu_r", center=0,
                       figsize=(8, max(6, len(matrix) * 0.18)),
                       row_cluster=True, col_cluster=False,
                       yticklabels=True, linewidths=0.3, linecolor="white")
    g.fig.suptitle("RBP log2FC Heatmap", y=1.02, fontsize=12, fontweight="bold")

    outpath = outdir / f"rbp_heatmap.{fig_format}"
    g.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(g.fig)
    print(f"  Saved: {outpath}")


def rbp_summary_table(
    condition_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str],
) -> None:
    """Export a summary Excel table of RBP genes across conditions."""
    outdir = Path(outdir)
    fc_col = cols["log2fc"]
    padj_col = cols["padj"]
    name_col = cols.get("gene_name", "")

    all_rows = []
    for name, data in condition_results.items():
        lbl = condition_labels[name]
        filt = data["deseq2_filtered"].get("all_genes", pd.DataFrame())
        if "is_RBP" not in filt.columns or len(filt) == 0:
            continue
        rbps = filt[filt["is_RBP"] == True]  # noqa: E712
        if len(rbps) == 0:
            continue
        for _, row in rbps.iterrows():
            entry = {"Condition": lbl}
            if name_col and name_col in rbps.columns:
                entry["Gene"] = row[name_col]
            entry["log2FC"] = row.get(fc_col, np.nan)
            entry["padj"] = row.get(padj_col, np.nan)
            entry["direction"] = row.get("direction", "")
            for rbp_key in ["is_RBP_Yael", "is_RBP_MW"]:
                if rbp_key in rbps.columns:
                    entry[rbp_key] = row.get(rbp_key, False)
            all_rows.append(entry)

    if not all_rows:
        print("  No RBP data for summary table")
        return

    summary_df = pd.DataFrame(all_rows)
    xlsx_path = outdir / "rbp_summary.xlsx"
    summary_df.to_excel(xlsx_path, index=False)
    print(f"  Saved: {xlsx_path}")


# ---------------------------------------------------------------------------
# Cross-condition biotype visualization
# ---------------------------------------------------------------------------


def cross_condition_biotype_comparison(
    condition_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str],
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Grouped bar + stacked % chart of DE gene counts by biotype group across conditions."""
    outdir = Path(outdir)
    count_data = {}
    for name, res in condition_results.items():
        filt = res["deseq2_filtered"].get("all_genes")
        if filt is None or len(filt) == 0:
            continue
        if cols["biotype"] not in filt.columns:
            continue
        groups = assign_biotype_group(filt[cols["biotype"]])
        count_data[condition_labels[name]] = groups.value_counts()

    if len(count_data) < 1:
        print("  Skipping cross-condition biotype comparison (no biotype data)")
        return

    all_groups = BIOTYPE_ORDER
    df = pd.DataFrame(count_data, index=all_groups).fillna(0).astype(int)
    df = df.loc[(df.sum(axis=1) > 0)]

    fig, axes = plt.subplots(1, 2, figsize=(14, max(5, len(df) * 0.6 + 2)))

    ax = axes[0]
    x = np.arange(len(df))
    n_conds = len(df.columns)
    width = 0.8 / max(n_conds, 1)
    cond_colors = sns.color_palette("tab10", n_conds)
    for i, cond_lbl in enumerate(df.columns):
        vals = df[cond_lbl].values
        ax.bar(x + i * width - (n_conds - 1) * width / 2,
               vals, width * 0.9, label=cond_lbl, color=cond_colors[i], alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=30, ha="right")
    ax.set_ylabel("DE Gene Count")
    ax.set_title("DE Genes per Biotype Group")
    ax.legend(fontsize=8)

    ax2 = axes[1]
    col_sums = df.sum(axis=0).replace(0, np.nan)
    pct = df.div(col_sums, axis=1) * 100
    bottom = np.zeros(n_conds)
    x2 = np.arange(n_conds)
    for grp in df.index:
        vals = pct.loc[grp].values
        color = BIOTYPE_COLORS.get(grp, "#999999")
        ax2.bar(x2, vals, bottom=bottom, color=color, alpha=0.85, label=grp)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 5:
                ax2.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                         fontsize=8, color="white" if v > 10 else "black")
        bottom += vals
    ax2.set_xticks(x2)
    ax2.set_xticklabels(df.columns, rotation=20, ha="right")
    ax2.set_ylabel("Percentage of DE Genes (%)")
    ax2.set_title("Biotype Composition per Condition")
    ax2.set_ylim(0, 105)
    ax2.legend(fontsize=8, loc="upper right")

    plt.tight_layout()
    outpath = outdir / f"cross_condition_biotype_comparison.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def cross_condition_biotype_direction(
    condition_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str],
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Faceted diverging bar chart: Up/Down DE gene counts by biotype group per condition."""
    outdir = Path(outdir)
    n = len(condition_results)
    if n == 0:
        return

    nrows, ncols = grid_dims(n)
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 5, nrows * max(3, len(BIOTYPE_ORDER) * 0.8 + 1)),
                              squeeze=False)
    axes_flat = [axes[r][c] for r in range(nrows) for c in range(ncols)]

    max_abs = 0
    panel_data = []
    for name, res in condition_results.items():
        filt = res["deseq2_filtered"].get("all_genes")
        cond_lbl = condition_labels[name]
        if filt is None or len(filt) == 0 or cols["biotype"] not in filt.columns:
            panel_data.append((cond_lbl, None, None))
            continue
        filt = filt.copy()
        filt["_group"] = assign_biotype_group(filt[cols["biotype"]])
        up_c = filt[filt["direction"] == "up"]["_group"].value_counts()
        dn_c = filt[filt["direction"] == "down"]["_group"].value_counts()
        max_abs = max(max_abs, up_c.max() if len(up_c) else 0,
                      dn_c.max() if len(dn_c) else 0)
        panel_data.append((cond_lbl, up_c, dn_c))

    xlim = max_abs * 1.15 if max_abs > 0 else 10

    for ax_idx, (cond_lbl, up_c, dn_c) in enumerate(panel_data):
        ax = axes_flat[ax_idx]
        if up_c is None and dn_c is None:
            ax.set_visible(False)
            continue

        groups = [g for g in BIOTYPE_ORDER
                  if (up_c is not None and g in up_c.index) or
                     (dn_c is not None and g in dn_c.index)]
        if not groups:
            ax.set_visible(False)
            continue

        y = range(len(groups))
        up_vals = [up_c.get(g, 0) if up_c is not None else 0 for g in groups]
        dn_vals = [-dn_c.get(g, 0) if dn_c is not None else 0 for g in groups]

        ax.barh(list(y), up_vals, color=COLOR_UP, alpha=0.85, label="Up")
        ax.barh(list(y), dn_vals, color=COLOR_DOWN, alpha=0.85, label="Down")
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlim(-xlim, xlim)
        ax.set_yticks(list(y))
        ax.set_yticklabels(groups, fontsize=9)
        ax.set_title(cond_lbl, fontsize=10)
        ax.set_xlabel("DE Gene Count", fontsize=8)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: str(int(abs(x)))))
        if ax_idx == 0:
            ax.legend(fontsize=8, loc="lower right")

    for ax_idx in range(len(panel_data), len(axes_flat)):
        axes_flat[ax_idx].set_visible(False)

    fig.suptitle("DE Gene Direction by Biotype Group (per Condition)", fontsize=12, y=1.01)
    plt.tight_layout()
    outpath = outdir / f"cross_condition_biotype_direction.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")
