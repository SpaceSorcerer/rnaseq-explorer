"""Comprehensive tests for rnaseq_explorer.engine.rmats module."""

from __future__ import annotations


import numpy as np
import pandas as pd
import pytest

from rnaseq_explorer.engine.rmats import (
    COORD_COLS,
    DEFAULT_RMATS_COLS,
    RMATS_EVENT_TYPES,
    _validate_rmats_columns,
    filter_rmats,
    load_all_rmats,
    make_event_key,
    normalize_rmats_columns,
    parse_inclevel_mean,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def se_df():
    """Minimal SE (Skipped Exon) rMATS DataFrame."""
    n = 15
    np.random.seed(50)
    return pd.DataFrame({
        "ID": list(range(1, n + 1)),
        "GeneID": [f"ENSG{i:011d}" for i in range(1, n + 1)],
        "geneSymbol": [f"GENE{i}" for i in range(1, n + 1)],
        "chr": ["chr1"] * n,
        "strand": ["+"] * n,
        "exonStart_0base": np.arange(1000, 1000 + n * 100, 100),
        "exonEnd": np.arange(1050, 1050 + n * 100, 100),
        "upstreamES": np.arange(800, 800 + n * 100, 100),
        "upstreamEE": np.arange(900, 900 + n * 100, 100),
        "downstreamES": np.arange(1200, 1200 + n * 100, 100),
        "downstreamEE": np.arange(1300, 1300 + n * 100, 100),
        "PValue": np.concatenate([
            np.random.uniform(1e-8, 0.01, 8),
            np.random.uniform(0.1, 1.0, 7),
        ]),
        "FDR": np.concatenate([
            np.random.uniform(1e-6, 0.01, 8),
            np.random.uniform(0.1, 1.0, 7),
        ]),
        "IncLevelDifference": np.concatenate([
            np.random.uniform(0.15, 0.5, 5),
            np.random.uniform(-0.5, -0.15, 3),
            np.random.uniform(-0.05, 0.05, 7),
        ]),
        "IncLevel1": ["0.3,0.4,0.35"] * n,
        "IncLevel2": ["0.5,0.6,0.55"] * n,
    })


@pytest.fixture
def mxe_df():
    """Minimal MXE (Mutually Exclusive Exon) DataFrame."""
    n = 5
    return pd.DataFrame({
        "ID": list(range(1, n + 1)),
        "GeneID": [f"ENSG{i:011d}" for i in range(100, 100 + n)],
        "geneSymbol": [f"MXE_GENE{i}" for i in range(1, n + 1)],
        "chr": ["chr2"] * n,
        "strand": ["-"] * n,
        "1stExonStart_0base": [2000, 3000, 4000, 5000, 6000],
        "1stExonEnd": [2100, 3100, 4100, 5100, 6100],
        "2ndExonStart_0base": [2200, 3200, 4200, 5200, 6200],
        "2ndExonEnd": [2300, 3300, 4300, 5300, 6300],
        "upstreamES": [1800, 2800, 3800, 4800, 5800],
        "upstreamEE": [1900, 2900, 3900, 4900, 5900],
        "downstreamES": [2500, 3500, 4500, 5500, 6500],
        "downstreamEE": [2600, 3600, 4600, 5600, 6600],
        "PValue": [0.001, 0.01, 0.05, 0.1, 0.5],
        "FDR": [0.005, 0.03, 0.08, 0.2, 0.8],
        "IncLevelDifference": [0.3, -0.2, 0.15, -0.05, 0.02],
    })


# ---------------------------------------------------------------------------
# _validate_rmats_columns tests
# ---------------------------------------------------------------------------


class TestValidateRmatsColumns:
    """Tests for _validate_rmats_columns()."""

    def test_no_renaming_needed(self, se_df):
        result = _validate_rmats_columns(se_df, DEFAULT_RMATS_COLS)
        assert "FDR" in result.columns
        assert "PValue" in result.columns

    def test_alias_renaming(self):
        df = pd.DataFrame({
            "id": [1],
            "gene_id": ["ENSG1"],
            "gene_name": ["A"],
            "pvalue": [0.01],
            "fdr": [0.05],
            "IncLevelDifference": [0.2],
        })
        result = _validate_rmats_columns(df, DEFAULT_RMATS_COLS)
        assert "FDR" in result.columns
        assert "PValue" in result.columns

    def test_missing_required_col_raises(self):
        df = pd.DataFrame({"random_col": [1]})
        with pytest.raises(ValueError, match="Required rMATS column"):
            _validate_rmats_columns(df, DEFAULT_RMATS_COLS)


# ---------------------------------------------------------------------------
# normalize_rmats_columns tests
# ---------------------------------------------------------------------------


class TestNormalizeRmatsColumns:
    """Tests for normalize_rmats_columns()."""

    def test_strips_ensembl_versions(self):
        df = pd.DataFrame({
            "GeneID": ["ENSG00000123456.12", "ENSG00000654321.3"],
            "geneSymbol": ["A", "B"],
            "ID": [1, 2],
            "PValue": [0.01, 0.02],
            "FDR": [0.05, 0.1],
            "IncLevelDifference": [0.2, -0.3],
        })
        result = normalize_rmats_columns(df, DEFAULT_RMATS_COLS)
        assert result["GeneID"].iloc[0] == "ENSG00000123456"


# ---------------------------------------------------------------------------
# filter_rmats tests
# ---------------------------------------------------------------------------


class TestFilterRmats:
    """Tests for filter_rmats()."""

    def test_basic_fdr_filtering(self, se_df):
        raw, filtered = filter_rmats(se_df, fdr_cutoff=0.05, dpsi_cutoff=0.1)
        assert len(filtered) <= len(raw)
        assert (filtered["FDR"] < 0.05).all()
        assert (filtered["IncLevelDifference"].abs() >= 0.1).all()

    def test_pvalue_filtering(self, se_df):
        raw, filtered = filter_rmats(
            se_df, pval_cutoff=0.01, dpsi_cutoff=0.1, use_fdr=False,
        )
        assert (filtered["PValue"] < 0.01).all()

    def test_dual_filter_mode(self, se_df):
        raw, filtered = filter_rmats(
            se_df,
            fdr_cutoff=0.05,
            pval_cutoff=0.01,
            dpsi_cutoff=0.1,
            dual_filter=True,
        )
        if len(filtered) > 0:
            assert (filtered["FDR"] < 0.05).all()
            assert (filtered["PValue"] < 0.01).all()
            assert (filtered["IncLevelDifference"].abs() >= 0.1).all()

    def test_strict_cutoffs_reduce_results(self, se_df):
        _, filtered_loose = filter_rmats(se_df, fdr_cutoff=0.1, dpsi_cutoff=0.05)
        _, filtered_strict = filter_rmats(se_df, fdr_cutoff=0.01, dpsi_cutoff=0.2)
        assert len(filtered_strict) <= len(filtered_loose)

    def test_empty_dataframe(self):
        df = pd.DataFrame({
            "ID": pd.Series(dtype=int),
            "GeneID": pd.Series(dtype=str),
            "geneSymbol": pd.Series(dtype=str),
            "PValue": pd.Series(dtype=float),
            "FDR": pd.Series(dtype=float),
            "IncLevelDifference": pd.Series(dtype=float),
        })
        raw, filtered = filter_rmats(df)
        assert len(filtered) == 0

    def test_all_filtered_out(self, se_df):
        raw, filtered = filter_rmats(se_df, fdr_cutoff=1e-20, dpsi_cutoff=0.99)
        assert len(filtered) == 0


# ---------------------------------------------------------------------------
# make_event_key tests
# ---------------------------------------------------------------------------


class TestMakeEventKey:
    """Tests for make_event_key()."""

    def test_se_event_key(self, se_df):
        keys = make_event_key(se_df, "SE")
        assert len(keys) == len(se_df)
        assert keys.iloc[0] != ""
        # Should contain colon-separated coordinate fields
        assert ":" in keys.iloc[0]

    def test_mxe_event_key(self, mxe_df):
        keys = make_event_key(mxe_df, "MXE")
        assert len(keys) == len(mxe_df)
        assert keys.iloc[0] != ""

    def test_all_event_types_have_coord_cols(self):
        for et in RMATS_EVENT_TYPES:
            assert et in COORD_COLS
            assert len(COORD_COLS[et]) > 0

    def test_unknown_event_type(self, se_df):
        keys = make_event_key(se_df, "UNKNOWN")
        assert (keys == "").all()

    def test_missing_coordinate_columns(self):
        df = pd.DataFrame({"chr": ["chr1"], "strand": ["+"]})
        keys = make_event_key(df, "SE")
        assert (keys == "").all()

    def test_unique_keys_for_different_events(self, se_df):
        keys = make_event_key(se_df, "SE")
        # Each row should have a unique key
        assert keys.nunique() == len(se_df)


# ---------------------------------------------------------------------------
# parse_inclevel_mean tests
# ---------------------------------------------------------------------------


class TestParseInclevelMean:
    """Tests for parse_inclevel_mean()."""

    def test_basic_parsing(self):
        s = pd.Series(["0.3,0.4,0.35"])
        result = parse_inclevel_mean(s)
        assert result.iloc[0] == pytest.approx(0.35, abs=1e-6)

    def test_handles_na(self):
        s = pd.Series([np.nan, "0.5,0.6", "NA,0.3"])
        result = parse_inclevel_mean(s)
        assert pd.isna(result.iloc[0])
        assert result.iloc[1] == pytest.approx(0.55, abs=1e-6)
        assert result.iloc[2] == pytest.approx(0.3, abs=1e-6)

    def test_single_value(self):
        s = pd.Series(["0.75"])
        result = parse_inclevel_mean(s)
        assert result.iloc[0] == pytest.approx(0.75)

    def test_empty_string(self):
        s = pd.Series([""])
        result = parse_inclevel_mean(s)
        assert pd.isna(result.iloc[0])


# ---------------------------------------------------------------------------
# load_all_rmats tests
# ---------------------------------------------------------------------------


class TestLoadAllRmats:
    """Tests for load_all_rmats()."""

    def test_loads_available_files(self, tmp_path, se_df):
        se_path = tmp_path / "SE.MATS.JCEC.txt"
        se_df.to_csv(se_path, sep="\t", index=False)

        result = load_all_rmats(tmp_path, event_types=["SE"])
        assert "SE" in result
        assert len(result["SE"]) == len(se_df)

    def test_skips_missing_files(self, tmp_path):
        result = load_all_rmats(tmp_path, event_types=["SE"])
        assert "SE" not in result

    def test_empty_directory(self, tmp_path):
        result = load_all_rmats(tmp_path)
        assert len(result) == 0
