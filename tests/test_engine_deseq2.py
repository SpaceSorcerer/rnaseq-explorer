"""Comprehensive tests for rnaseq_explorer.engine.deseq2 module."""

from __future__ import annotations

from unittest.mock import patch, MagicMock

import numpy as np
import pandas as pd
import pytest

from rnaseq_explorer.engine.deseq2 import (
    DEFAULT_DESEQ2_COLS,
    _DESEQ2_ALIASES,
    _resolve_column,
    _strip_ensembl_version,
    assign_biotype_group,
    best_gene_key,
    filter_deseq2,
    load_file,
    normalize_deseq2_columns,
    validate_columns,
    extract_gene_sets,
    fetch_gene_names,
    _bh_correction,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def deseq2_df():
    """DataFrame with realistic DESeq2 columns (~30 rows)."""
    np.random.seed(42)
    n = 30
    genes = [f"GENE{i}" for i in range(1, n + 1)]
    genes[0] = "MIAT"
    genes[1] = "QKI"
    genes[2] = "TP53"
    genes[3] = "BRCA1"

    gene_ids = [f"ENSG{i:011d}" for i in range(1, n + 1)]

    log2fc = np.concatenate([
        np.random.uniform(1.5, 4.0, 8),
        np.random.uniform(-4.0, -1.5, 7),
        np.random.uniform(-0.3, 0.3, 15),
    ])
    padj = np.concatenate([
        np.random.uniform(1e-10, 0.01, 8),
        np.random.uniform(1e-10, 0.01, 7),
        np.random.uniform(0.1, 1.0, 15),
    ])
    basemean = np.random.uniform(50, 5000, n)
    pvalue = padj * np.random.uniform(0.5, 1.0, n)
    lfcSE = np.abs(log2fc) * np.random.uniform(0.1, 0.3, n)
    stat = log2fc / lfcSE

    biotypes = (
        ["protein_coding"] * 10
        + ["lncrna"] * 5
        + ["pseudogene"] * 5
        + ["mirna"] * 5
        + ["unknown_type"] * 5
    )

    return pd.DataFrame({
        "gene_id": gene_ids,
        "gene_name": genes,
        "baseMean": basemean,
        "log2FoldChange": log2fc,
        "lfcSE": lfcSE,
        "stat": stat,
        "pvalue": pvalue,
        "padj": padj,
        "biotype": biotypes,
    })


# ---------------------------------------------------------------------------
# best_gene_key tests
# ---------------------------------------------------------------------------


class TestBestGeneKey:
    """Tests for best_gene_key()."""

    def test_prefers_ensembl_ids(self, deseq2_df):
        col, desc = best_gene_key(deseq2_df, DEFAULT_DESEQ2_COLS)
        assert col == "gene_id"
        assert "Ensembl" in desc

    def test_falls_back_to_gene_name(self):
        df = pd.DataFrame({
            "gene_id": ["BRCA1", "TP53", "MIAT"],
            "gene_name": ["BRCA1", "TP53", "MIAT"],
        })
        col, desc = best_gene_key(df, DEFAULT_DESEQ2_COLS)
        assert col == "gene_name"
        assert "name" in desc.lower()

    def test_missing_gene_id_col(self):
        df = pd.DataFrame({"gene_name": ["A", "B"]})
        col, _ = best_gene_key(df, DEFAULT_DESEQ2_COLS)
        assert col == "gene_name"

    def test_empty_dataframe(self):
        df = pd.DataFrame({"gene_id": pd.Series(dtype=str), "gene_name": pd.Series(dtype=str)})
        col, _ = best_gene_key(df, DEFAULT_DESEQ2_COLS)
        assert col == "gene_name"


# ---------------------------------------------------------------------------
# _strip_ensembl_version tests
# ---------------------------------------------------------------------------


class TestStripEnsemblVersion:
    """Tests for _strip_ensembl_version()."""

    def test_strips_version_suffix(self):
        s = pd.Series(["ENSG00000123456.12", "ENSG00000654321.3"])
        result = _strip_ensembl_version(s)
        assert result.iloc[0] == "ENSG00000123456"
        assert result.iloc[1] == "ENSG00000654321"

    def test_preserves_non_ensembl(self):
        s = pd.Series(["BRCA1", "TP53", "MIAT"])
        result = _strip_ensembl_version(s)
        assert result.tolist() == ["BRCA1", "TP53", "MIAT"]

    def test_handles_mixed_ids(self):
        s = pd.Series(["ENSG00000123456.5", "BRCA1", "ENSMUSG00000000001.10"])
        result = _strip_ensembl_version(s)
        assert result.iloc[0] == "ENSG00000123456"
        assert result.iloc[1] == "BRCA1"
        assert result.iloc[2] == "ENSMUSG00000000001"

    def test_no_version_ensembl(self):
        s = pd.Series(["ENSG00000123456"])
        result = _strip_ensembl_version(s)
        assert result.iloc[0] == "ENSG00000123456"

    def test_empty_series(self):
        s = pd.Series(dtype=str)
        result = _strip_ensembl_version(s)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# _resolve_column tests
# ---------------------------------------------------------------------------


class TestResolveColumn:
    """Tests for _resolve_column()."""

    def test_exact_match(self):
        df = pd.DataFrame({"log2FoldChange": [1.0]})
        result = _resolve_column(df, "log2fc", "log2FoldChange", _DESEQ2_ALIASES)
        assert result == "log2FoldChange"

    def test_case_insensitive_match(self):
        df = pd.DataFrame({"LOG2FOLDCHANGE": [1.0]})
        result = _resolve_column(df, "log2fc", "log2FoldChange", _DESEQ2_ALIASES)
        assert result == "LOG2FOLDCHANGE"

    def test_alias_match(self):
        df = pd.DataFrame({"logFC": [1.0]})
        result = _resolve_column(df, "log2fc", "log2FoldChange", _DESEQ2_ALIASES)
        assert result == "logFC"

    def test_no_match_returns_none(self):
        df = pd.DataFrame({"totally_random": [1.0]})
        result = _resolve_column(df, "log2fc", "log2FoldChange", _DESEQ2_ALIASES)
        assert result is None

    def test_padj_alias_detection(self):
        df = pd.DataFrame({"adj.P.Val": [0.05]})
        result = _resolve_column(df, "padj", "padj", _DESEQ2_ALIASES)
        assert result == "adj.P.Val"

    def test_gene_id_alias(self):
        df = pd.DataFrame({"ensembl_gene_id": ["ENSG123"]})
        result = _resolve_column(df, "gene_id", "gene_id", _DESEQ2_ALIASES)
        assert result == "ensembl_gene_id"


# ---------------------------------------------------------------------------
# normalize_deseq2_columns tests
# ---------------------------------------------------------------------------


class TestNormalizeDeseq2Columns:
    """Tests for normalize_deseq2_columns()."""

    def test_renames_aliased_columns(self):
        df = pd.DataFrame({
            "ensembl_gene_id": ["ENSG00000123456.5"],
            "symbol": ["BRCA1"],
            "logFC": [2.0],
            "AveExpr": [100.0],
            "adj.P.Val": [0.01],
            "P.Value": [0.001],
        })
        result = normalize_deseq2_columns(df, DEFAULT_DESEQ2_COLS)
        assert "gene_id" in result.columns
        assert "log2FoldChange" in result.columns

    def test_strips_ensembl_versions(self):
        df = pd.DataFrame({
            "gene_id": ["ENSG00000123456.12", "ENSG00000654321.3"],
            "gene_name": ["A", "B"],
            "log2FoldChange": [1.0, -1.0],
            "baseMean": [100, 200],
            "padj": [0.01, 0.02],
            "pvalue": [0.001, 0.002],
        })
        result = normalize_deseq2_columns(df, DEFAULT_DESEQ2_COLS)
        assert result["gene_id"].iloc[0] == "ENSG00000123456"


# ---------------------------------------------------------------------------
# load_file tests
# ---------------------------------------------------------------------------


class TestLoadFile:
    """Tests for load_file()."""

    def test_load_csv(self, tmp_path):
        csv_path = tmp_path / "test.csv"
        df = pd.DataFrame({"gene": ["A", "B"], "value": [1, 2]})
        df.to_csv(csv_path, index=False)
        result = load_file(csv_path, name="test CSV")
        assert len(result) == 2
        assert "gene" in result.columns

    def test_load_tsv(self, tmp_path):
        tsv_path = tmp_path / "test.tsv"
        df = pd.DataFrame({"gene": ["A", "B"], "value": [1, 2]})
        df.to_csv(tsv_path, sep="\t", index=False)
        result = load_file(tsv_path, name="test TSV")
        assert len(result) == 2

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_file(tmp_path / "nonexistent.csv")

    def test_load_xlsx(self, tmp_path):
        xlsx_path = tmp_path / "test.xlsx"
        df = pd.DataFrame({"gene": ["A", "B"], "value": [1, 2]})
        df.to_excel(xlsx_path, index=False)
        result = load_file(xlsx_path, name="test Excel")
        assert len(result) == 2


# ---------------------------------------------------------------------------
# validate_columns tests
# ---------------------------------------------------------------------------


class TestValidateColumns:
    """Tests for validate_columns()."""

    def test_passes_when_all_present(self, deseq2_df):
        validate_columns(deseq2_df, ["gene_name", "padj", "log2FoldChange"])

    def test_raises_on_missing(self, deseq2_df):
        with pytest.raises(KeyError, match="Missing columns"):
            validate_columns(deseq2_df, ["nonexistent_column"])

    def test_ignores_empty_string_cols(self, deseq2_df):
        validate_columns(deseq2_df, ["gene_name", "", "padj"])


# ---------------------------------------------------------------------------
# assign_biotype_group tests
# ---------------------------------------------------------------------------


class TestAssignBiotypeGroup:
    """Tests for assign_biotype_group()."""

    def test_protein_coding(self):
        s = pd.Series(["protein_coding", "Protein_Coding"])
        result = assign_biotype_group(s)
        assert result.iloc[0] == "Protein Coding"

    def test_lncrna_variants(self):
        s = pd.Series(["lncRNA", "lincRNA", "antisense", "sense_intronic"])
        result = assign_biotype_group(s)
        assert (result == "lncRNA").all()

    def test_pseudogene_variants(self):
        s = pd.Series(["pseudogene", "processed_pseudogene", "unprocessed_pseudogene"])
        result = assign_biotype_group(s)
        assert (result == "Pseudogene").all()

    def test_small_ncrna(self):
        s = pd.Series(["miRNA", "snRNA", "snoRNA", "tRNA"])
        result = assign_biotype_group(s)
        assert (result == "Small ncRNA").all()

    def test_unknown_becomes_other(self):
        s = pd.Series(["totally_unknown_biotype", "some_weird_type"])
        result = assign_biotype_group(s)
        assert (result == "Other").all()

    def test_empty_series(self):
        s = pd.Series(dtype=str)
        result = assign_biotype_group(s)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# filter_deseq2 tests
# ---------------------------------------------------------------------------


class TestFilterDeseq2:
    """Tests for filter_deseq2()."""

    def test_basic_filtering(self, deseq2_df):
        df_bio = deseq2_df.copy()
        df_bio["biotype_group"] = assign_biotype_group(df_bio["biotype"])
        raw, filtered = filter_deseq2(
            df_bio, DEFAULT_DESEQ2_COLS,
            log2fc_cutoff=1.0, basemean_cutoff=10.0, padj_cutoff=0.05,
        )
        assert len(filtered) <= len(raw)
        assert "direction" in filtered.columns
        assert set(filtered["direction"].unique()).issubset({"up", "down"})

    def test_strict_cutoffs_reduce_results(self, deseq2_df):
        _, filtered_loose = filter_deseq2(
            deseq2_df, DEFAULT_DESEQ2_COLS,
            log2fc_cutoff=0.5, padj_cutoff=0.1,
        )
        _, filtered_strict = filter_deseq2(
            deseq2_df, DEFAULT_DESEQ2_COLS,
            log2fc_cutoff=2.0, padj_cutoff=0.001,
        )
        assert len(filtered_strict) <= len(filtered_loose)

    def test_nan_padj_dropped(self):
        df = pd.DataFrame({
            "gene_id": ["A", "B", "C"],
            "gene_name": ["A", "B", "C"],
            "log2FoldChange": [2.0, -2.0, 1.5],
            "baseMean": [100, 200, 150],
            "padj": [0.01, np.nan, 0.02],
            "pvalue": [0.001, np.nan, 0.002],
            "biotype": ["protein_coding"] * 3,
        })
        raw, filtered = filter_deseq2(df, DEFAULT_DESEQ2_COLS)
        assert len(raw) == 2  # B was dropped for NaN padj

    def test_all_nan_padj(self):
        df = pd.DataFrame({
            "gene_id": ["A", "B"],
            "gene_name": ["A", "B"],
            "log2FoldChange": [2.0, -2.0],
            "baseMean": [100, 200],
            "padj": [np.nan, np.nan],
            "pvalue": [np.nan, np.nan],
            "biotype": ["protein_coding"] * 2,
        })
        raw, filtered = filter_deseq2(df, DEFAULT_DESEQ2_COLS)
        assert len(raw) == 0
        assert len(filtered) == 0

    def test_empty_dataframe(self):
        df = pd.DataFrame({
            "gene_id": pd.Series(dtype=str),
            "gene_name": pd.Series(dtype=str),
            "log2FoldChange": pd.Series(dtype=float),
            "baseMean": pd.Series(dtype=float),
            "padj": pd.Series(dtype=float),
            "pvalue": pd.Series(dtype=float),
            "biotype": pd.Series(dtype=str),
        })
        raw, filtered = filter_deseq2(df, DEFAULT_DESEQ2_COLS)
        assert len(filtered) == 0

    def test_biotype_filter_protein_coding(self, deseq2_df):
        df = deseq2_df.copy()
        df["biotype_group"] = assign_biotype_group(df["biotype"])
        raw, filtered = filter_deseq2(
            df, DEFAULT_DESEQ2_COLS,
            biotype_filter="protein_coding",
        )
        # All remaining rows should have protein_coding biotype
        biotypes_lower = raw["biotype"].str.lower().str.replace(" ", "_")
        assert (biotypes_lower == "protein_coding").all()

    def test_direction_assignment(self, deseq2_df):
        _, filtered = filter_deseq2(deseq2_df, DEFAULT_DESEQ2_COLS)
        if len(filtered) > 0:
            up = filtered[filtered["direction"] == "up"]
            down = filtered[filtered["direction"] == "down"]
            assert (up["log2FoldChange"] > 0).all()
            assert (down["log2FoldChange"] < 0).all()

    def test_single_gene_passing(self):
        df = pd.DataFrame({
            "gene_id": ["ENSG00000000001"],
            "gene_name": ["MIAT"],
            "log2FoldChange": [3.0],
            "baseMean": [500.0],
            "padj": [0.001],
            "pvalue": [0.0001],
            "biotype": ["lncrna"],
        })
        _, filtered = filter_deseq2(df, DEFAULT_DESEQ2_COLS)
        assert len(filtered) == 1
        assert filtered["direction"].iloc[0] == "up"


# ---------------------------------------------------------------------------
# extract_gene_sets tests
# ---------------------------------------------------------------------------


class TestExtractGeneSets:
    """Tests for extract_gene_sets()."""

    def test_extracts_up_down_all(self):
        deg_df = pd.DataFrame({
            "gene_name": ["A", "B", "C", "D"],
            "direction": ["up", "up", "down", "down"],
        })
        condition_results = {
            "cond1": {
                "deseq2_filtered": {"all_genes": deg_df},
                "deseq2_raw": deg_df,
            }
        }
        gene_sets = extract_gene_sets(condition_results, DEFAULT_DESEQ2_COLS)
        assert "cond1" in gene_sets
        assert len(gene_sets["cond1"]["all"]) == 4
        assert len(gene_sets["cond1"]["up"]) == 2
        assert len(gene_sets["cond1"]["down"]) == 2

    def test_empty_condition(self):
        condition_results = {
            "cond1": {
                "deseq2_filtered": {"all_genes": pd.DataFrame()},
                "deseq2_raw": pd.DataFrame(),
            }
        }
        gene_sets = extract_gene_sets(condition_results, DEFAULT_DESEQ2_COLS)
        assert gene_sets["cond1"]["all"] == set()


# ---------------------------------------------------------------------------
# _bh_correction tests
# ---------------------------------------------------------------------------


class TestBhCorrection:
    """Tests for _bh_correction()."""

    def test_empty_input(self):
        result = _bh_correction([])
        assert len(result) == 0

    def test_single_value(self):
        result = _bh_correction([0.05])
        assert len(result) == 1
        assert result[0] == pytest.approx(0.05)

    def test_preserves_order(self):
        pvals = [0.01, 0.05, 0.10]
        result = _bh_correction(pvals)
        assert len(result) == 3
        # Adjusted p-values should be >= raw
        for raw, adj in zip(pvals, result):
            assert adj >= raw - 1e-10

    def test_adjusted_bounded_by_1(self):
        pvals = [0.5, 0.8, 0.99]
        result = _bh_correction(pvals)
        assert all(v <= 1.0 for v in result)


# ---------------------------------------------------------------------------
# fetch_gene_names tests (mocked)
# ---------------------------------------------------------------------------


class TestFetchGeneNames:
    """Tests for fetch_gene_names() with mocked API."""

    @patch("rnaseq_explorer.engine.deseq2.urllib.request.urlopen")
    def test_basic_lookup(self, mock_urlopen):
        mock_response = MagicMock()
        mock_response.read.return_value = b'[{"query": "ENSG00000123456", "symbol": "BRCA1"}]'
        mock_response.__enter__ = lambda self: self
        mock_response.__exit__ = lambda self, *args: None
        mock_urlopen.return_value = mock_response

        result = fetch_gene_names(["ENSG00000123456"])
        assert result == {"ENSG00000123456": "BRCA1"}

    def test_skips_non_ensembl_ids(self):
        result = fetch_gene_names(["BRCA1", "TP53"])
        assert result == {}


# ---------------------------------------------------------------------------
# Duplicate gene name edge case
# ---------------------------------------------------------------------------


class TestDuplicateGeneNames:
    """Tests for duplicate gene name handling in deseq2."""

    def test_filter_with_duplicate_gene_names(self):
        df = pd.DataFrame({
            "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            "gene_name": ["MIAT", "MIAT", "QKI"],
            "log2FoldChange": [2.0, -1.5, 3.0],
            "baseMean": [200, 300, 400],
            "padj": [0.001, 0.003, 0.01],
            "pvalue": [0.0001, 0.0003, 0.001],
            "biotype": ["lncrna", "lncrna", "protein_coding"],
        })
        raw, filtered = filter_deseq2(df, DEFAULT_DESEQ2_COLS)
        # Both MIAT entries should be filtered separately
        assert len(filtered) >= 2
