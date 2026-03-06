#!/usr/bin/env python3
"""
Alternative Splicing Analysis Module (rMATS)
============================================

This module handles rMATS alternative splicing analysis:
- Loading rMATS output files (JCEC format)
- Filtering splicing events by significance (FDR, p-value, dPSI)
- Parsing event types (SE, A3SS, A5SS, MXE, RI)
- Generating splicing-specific visualizations
- Exporting filtered splicing events

All filtering thresholds and event types are configurable via parameters
or can use global defaults from the main pipeline.

Author: Data - RNA-seq Pipeline Engineer
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional


# ─── Constants ───────────────────────────────────────────────────────────────

# Default rMATS cutoffs (can be overridden via function parameters)
DEFAULT_FDR_CUTOFF = 0.05
DEFAULT_PVALUE_CUTOFF = 0.01
DEFAULT_INCLEVEL_DIFF_CUTOFF = 0.1
DEFAULT_EVENT_TYPES = ["SE", "A3SS", "A5SS", "RI", "MXE"]

# Event type colors for plotting
EVENT_COLORS = {
    "SE": "#E64B35",
    "A3SS": "#4DBBD5",
    "A5SS": "#00A087",
    "RI": "#3C5488",
    "MXE": "#F39B7F"
}


# ─── rMATS File Loading ──────────────────────────────────────────────────────

def load_rmats_file(filepath: Path, event_type: str) -> pd.DataFrame:
    """
    Load a single rMATS output file.

    Parameters
    ----------
    filepath : Path
        Path to rMATS .MATS.JCEC.txt file
    event_type : str
        Event type (SE, A3SS, A5SS, RI, MXE)

    Returns
    -------
    pd.DataFrame
        Loaded rMATS data with event_type column added
    """
    if not filepath.exists():
        return pd.DataFrame()

    df = pd.read_csv(filepath, sep="\t")
    df["event_type"] = event_type

    return df


def load_rmats_all_events(
    rmats_dir: Path,
    event_types: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Load all rMATS event types from a directory.

    Parameters
    ----------
    rmats_dir : Path
        Directory containing rMATS output files
    event_types : List[str], optional
        List of event types to load. If None, uses DEFAULT_EVENT_TYPES

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with all events from all types
    """
    if event_types is None:
        event_types = DEFAULT_EVENT_TYPES

    if not rmats_dir.exists():
        print(f"  Warning: rMATS directory not found: {rmats_dir}")
        return pd.DataFrame()

    all_events = []
    for event_type in event_types:
        filepath = rmats_dir / f"{event_type}.MATS.JCEC.txt"
        df = load_rmats_file(filepath, event_type)

        if len(df) > 0:
            all_events.append(df)
            print(f"  Loaded {event_type}: {len(df):,} events")
        else:
            print(f"  Warning: {event_type} file not found or empty")

    if all_events:
        combined = pd.concat(all_events, ignore_index=True)
        return combined

    return pd.DataFrame()


# ─── Splicing Event Filtering ────────────────────────────────────────────────

def filter_rmats_events(
    df: pd.DataFrame,
    fdr_cutoff: float = DEFAULT_FDR_CUTOFF,
    pvalue_cutoff: Optional[float] = None,
    inclevel_diff_cutoff: float = DEFAULT_INCLEVEL_DIFF_CUTOFF,
    use_pval_and_fdr: bool = True
) -> pd.DataFrame:
    """
    Filter rMATS events by significance thresholds.

    Parameters
    ----------
    df : pd.DataFrame
        Raw rMATS data
    fdr_cutoff : float, default=0.05
        FDR threshold for significance
    pvalue_cutoff : float, optional
        P-value threshold (only used if use_pval_and_fdr=True)
    inclevel_diff_cutoff : float, default=0.1
        Minimum absolute IncLevelDifference (dPSI)
    use_pval_and_fdr : bool, default=True
        If True, require BOTH FDR and p-value to pass cutoffs

    Returns
    -------
    pd.DataFrame
        Filtered significant splicing events with direction column
    """
    if len(df) == 0:
        return pd.DataFrame()

    # Build filter mask
    if use_pval_and_fdr and pvalue_cutoff is not None:
        # Dual filtering: FDR + p-value
        mask = (
            (df["FDR"] < fdr_cutoff) &
            (df["PValue"] < pvalue_cutoff) &
            (df["IncLevelDifference"].abs() >= inclevel_diff_cutoff)
        )
    else:
        # FDR-only filtering
        mask = (
            (df["FDR"] < fdr_cutoff) &
            (df["IncLevelDifference"].abs() >= inclevel_diff_cutoff)
        )

    filtered = df[mask].copy()

    # Add direction annotation
    filtered["direction"] = np.where(
        filtered["IncLevelDifference"] > 0,
        "included",
        "excluded"
    )

    return filtered


def filter_rmats_by_event_type(
    df: pd.DataFrame,
    event_type: str,
    **filter_kwargs
) -> pd.DataFrame:
    """
    Filter rMATS data for a specific event type.

    Parameters
    ----------
    df : pd.DataFrame
        Combined rMATS data (all event types)
    event_type : str
        Event type to filter (SE, A3SS, A5SS, RI, MXE)
    **filter_kwargs
        Additional keyword arguments passed to filter_rmats_events()

    Returns
    -------
    pd.DataFrame
        Filtered events of specified type
    """
    if len(df) == 0:
        return pd.DataFrame()

    subset = df[df["event_type"] == event_type].copy()
    return filter_rmats_events(subset, **filter_kwargs)


# ─── rMATS Pipeline Wrapper ──────────────────────────────────────────────────

def load_rmats(
    condition: Dict,
    base_dir: Path,
    fdr_cutoff: float = DEFAULT_FDR_CUTOFF,
    pvalue_cutoff: float = DEFAULT_PVALUE_CUTOFF,
    inclevel_diff_cutoff: float = DEFAULT_INCLEVEL_DIFF_CUTOFF,
    use_pval_and_fdr: bool = True,
    event_types: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Load and filter rMATS results for a single condition.

    This is the main entry point function used by the analysis pipeline.
    It combines loading all event types and filtering by significance.

    Parameters
    ----------
    condition : Dict
        Condition dictionary containing 'rmats_dir' and 'label' keys
    base_dir : Path
        Base directory containing rMATS output folders
    fdr_cutoff : float, default=0.05
        FDR threshold
    pvalue_cutoff : float, default=0.01
        P-value threshold
    inclevel_diff_cutoff : float, default=0.1
        Minimum |dPSI| threshold
    use_pval_and_fdr : bool, default=True
        Require both FDR and p-value to pass cutoffs
    event_types : List[str], optional
        List of event types to load

    Returns
    -------
    pd.DataFrame
        Filtered significant splicing events across all event types
    """
    if event_types is None:
        event_types = DEFAULT_EVENT_TYPES

    rmats_dir = base_dir / condition["rmats_dir"]
    label = condition["label"]

    if not rmats_dir.exists():
        print(f"  Warning: rMATS directory not found: {rmats_dir}")
        return pd.DataFrame()

    all_events = []
    for event_type in event_types:
        filepath = rmats_dir / f"{event_type}.MATS.JCEC.txt"
        if not filepath.exists():
            print(f"  Warning: {event_type} file not found")
            continue

        # Load raw data
        df = pd.read_csv(filepath, sep="\t")
        df["event_type"] = event_type

        # Filter by significance
        if use_pval_and_fdr:
            sig = df[
                (df["FDR"] < fdr_cutoff) &
                (df["PValue"] < pvalue_cutoff) &
                (df["IncLevelDifference"].abs() >= inclevel_diff_cutoff)
            ].copy()
        else:
            sig = df[
                (df["FDR"] < fdr_cutoff) &
                (df["IncLevelDifference"].abs() >= inclevel_diff_cutoff)
            ].copy()

        # Add direction
        sig["direction"] = np.where(sig["IncLevelDifference"] > 0, "included", "excluded")
        all_events.append(sig)

        print(f"  {label} {event_type}: {len(sig):,} significant events "
              f"(of {len(df):,} total)")

    if all_events:
        combined = pd.concat(all_events, ignore_index=True)
        print(f"  Total significant splicing events: {len(combined):,}")
        return combined

    return pd.DataFrame()


# ─── Splicing Summary Statistics ─────────────────────────────────────────────

def summarize_splicing_events(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics for splicing events.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered rMATS events

    Returns
    -------
    pd.DataFrame
        Summary table with counts per event type and direction
    """
    if len(df) == 0:
        return pd.DataFrame()

    summary = df.groupby(["event_type", "direction"]).size().reset_index(name="count")
    return summary


def count_events_by_type(df: pd.DataFrame) -> Dict[str, int]:
    """
    Count splicing events by event type.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered rMATS events

    Returns
    -------
    Dict[str, int]
        Dictionary mapping event_type -> count
    """
    if len(df) == 0:
        return {}

    counts = df["event_type"].value_counts().to_dict()
    return counts


def get_event_genes(df: pd.DataFrame, event_type: str) -> List[str]:
    """
    Extract gene symbols for a specific event type.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered rMATS events
    event_type : str
        Event type to extract (SE, A3SS, etc.)

    Returns
    -------
    List[str]
        List of gene symbols (geneSymbol column)
    """
    if len(df) == 0:
        return []

    subset = df[df["event_type"] == event_type]

    if "geneSymbol" in subset.columns:
        return subset["geneSymbol"].dropna().unique().tolist()
    elif "GeneID" in subset.columns:
        return subset["GeneID"].dropna().unique().tolist()

    return []


# ─── Event Set Operations ────────────────────────────────────────────────────

def get_event_identifiers(df: pd.DataFrame, event_type: Optional[str] = None) -> set:
    """
    Extract unique event identifiers from rMATS data.

    Uses geneSymbol as the primary identifier. Falls back to ID column
    if geneSymbol is not available.

    Parameters
    ----------
    df : pd.DataFrame
        rMATS events
    event_type : str, optional
        If provided, filter to specific event type first

    Returns
    -------
    set
        Set of unique event identifiers
    """
    if len(df) == 0:
        return set()

    # Filter by event type if specified
    if event_type is not None:
        df = df[df["event_type"] == event_type]

    # Choose identifier column
    id_col = "geneSymbol" if "geneSymbol" in df.columns else "ID"

    return set(df[id_col].dropna())


def compare_event_sets(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    event_type: Optional[str] = None
) -> Tuple[set, set, set]:
    """
    Compare splicing events between two conditions.

    Parameters
    ----------
    df_a : pd.DataFrame
        rMATS events from condition A
    df_b : pd.DataFrame
        rMATS events from condition B
    event_type : str, optional
        If provided, compare only this event type

    Returns
    -------
    Tuple[set, set, set]
        (events in A only, events in B only, shared events)
    """
    events_a = get_event_identifiers(df_a, event_type)
    events_b = get_event_identifiers(df_b, event_type)

    only_a = events_a - events_b
    only_b = events_b - events_a
    shared = events_a & events_b

    return only_a, only_b, shared


def get_directional_event_sets(
    df: pd.DataFrame,
    event_type: Optional[str] = None
) -> Tuple[set, set]:
    """
    Split events into included vs excluded sets.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered rMATS events with direction column
    event_type : str, optional
        If provided, filter to specific event type

    Returns
    -------
    Tuple[set, set]
        (included events, excluded events)
    """
    if len(df) == 0:
        return set(), set()

    # Filter by event type if specified
    if event_type is not None:
        df = df[df["event_type"] == event_type]

    # Choose identifier column
    id_col = "geneSymbol" if "geneSymbol" in df.columns else "ID"

    included = set(df[df["direction"] == "included"][id_col].dropna())
    excluded = set(df[df["direction"] == "excluded"][id_col].dropna())

    return included, excluded


def compare_directional_events(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    event_type: Optional[str] = None
) -> Dict[str, set]:
    """
    Compare directional splicing events between two conditions.

    Returns sets of events that are:
    - Both included (concordant)
    - Both excluded (concordant)
    - Discordant (different directions)

    Parameters
    ----------
    df_a : pd.DataFrame
        rMATS events from condition A
    df_b : pd.DataFrame
        rMATS events from condition B
    event_type : str, optional
        If provided, compare only this event type

    Returns
    -------
    Dict[str, set]
        Dictionary with keys:
        - 'both_included': Events included in both
        - 'both_excluded': Events excluded in both
        - 'discordant': Events with opposite directions
        - 'only_a': Events only in condition A
        - 'only_b': Events only in condition B
    """
    # Filter by event type if specified
    if event_type is not None:
        df_a = df_a[df_a["event_type"] == event_type]
        df_b = df_b[df_b["event_type"] == event_type]

    # Get directional sets for each condition
    included_a, excluded_a = get_directional_event_sets(df_a)
    included_b, excluded_b = get_directional_event_sets(df_b)

    # All events in each condition
    all_a = included_a | excluded_a
    all_b = included_b | excluded_b

    # Concordant events
    both_included = included_a & included_b
    both_excluded = excluded_a & excluded_b

    # Discordant events
    discordant = (included_a & excluded_b) | (excluded_a & included_b)

    # Unique to each condition
    only_a = all_a - all_b
    only_b = all_b - all_a

    return {
        'both_included': both_included,
        'both_excluded': both_excluded,
        'discordant': discordant,
        'only_a': only_a,
        'only_b': only_b
    }


# ─── Export Functions ─────────────────────────────────────────────────────────

def export_splicing_events(
    df: pd.DataFrame,
    output_path: Path,
    include_raw_counts: bool = False
) -> None:
    """
    Export filtered splicing events to Excel.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered rMATS events
    output_path : Path
        Output Excel file path
    include_raw_counts : bool, default=False
        If True, include raw junction/read counts in export
    """
    if len(df) == 0:
        print(f"  No events to export")
        return

    # Select columns to export
    base_cols = [
        "event_type", "geneSymbol", "chr", "strand",
        "IncLevelDifference", "PValue", "FDR", "direction"
    ]

    # Add raw count columns if requested
    if include_raw_counts:
        count_cols = [
            "IJC_SAMPLE_1", "SJC_SAMPLE_1",
            "IJC_SAMPLE_2", "SJC_SAMPLE_2",
            "IncLevel1", "IncLevel2"
        ]
        export_cols = base_cols + count_cols
    else:
        export_cols = base_cols

    # Filter to available columns
    export_cols = [c for c in export_cols if c in df.columns]

    df[export_cols].to_excel(output_path, index=False)
    print(f"  Exported: {output_path.name}")


def export_splicing_summary(
    summary_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Export splicing summary statistics to CSV/Excel.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary table from summarize_splicing_events()
    output_path : Path
        Output file path (.csv or .xlsx)
    """
    if len(summary_df) == 0:
        print(f"  No summary data to export")
        return

    if output_path.suffix == ".xlsx":
        summary_df.to_excel(output_path, index=False)
    else:
        summary_df.to_csv(output_path, index=False)

    print(f"  Exported summary: {output_path.name}")


# ─── Utility Functions ────────────────────────────────────────────────────────

def validate_rmats_dataframe(df: pd.DataFrame) -> bool:
    """
    Validate that DataFrame contains required rMATS columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate

    Returns
    -------
    bool
        True if valid rMATS data, False otherwise
    """
    required_cols = ["FDR", "PValue", "IncLevelDifference"]

    for col in required_cols:
        if col not in df.columns:
            print(f"  Error: Missing required column '{col}'")
            return False

    return True


def get_event_type_description(event_type: str) -> str:
    """
    Get human-readable description of rMATS event type.

    Parameters
    ----------
    event_type : str
        Event type abbreviation (SE, A3SS, A5SS, RI, MXE)

    Returns
    -------
    str
        Full description of event type
    """
    descriptions = {
        "SE": "Skipped Exon",
        "A3SS": "Alternative 3' Splice Site",
        "A5SS": "Alternative 5' Splice Site",
        "RI": "Retained Intron",
        "MXE": "Mutually Exclusive Exons"
    }

    return descriptions.get(event_type, event_type)
