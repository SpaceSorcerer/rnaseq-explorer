"""Comprehensive tests for rnaseq_explorer.engine.gsea module."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rnaseq_explorer.engine.gsea import (
    _db_short_label,
    create_ranked_list,
    normalize_gsea_cols,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def deseq2_full_df():
    """Full unfiltered DESeq2 DataFrame for ranking."""
    np.random.seed(60)
    n = 50
    genes = [f"GENE{i}" for i in range(1, n + 1)]
    genes[0] = "MIAT"
    genes[1] = "QKI"

    stat = np.random.normal(0, 3, n)
    log2fc = stat * np.random.uniform(0.3, 0.6, n)

    return pd.DataFrame({
        "gene_name": genes,
        "stat": stat,
        "log2FoldChange": log2fc,
        "padj": np.random.uniform(0, 1, n),
    })


# ---------------------------------------------------------------------------
# create_ranked_list tests
# ---------------------------------------------------------------------------


class TestCreateRankedList:
    """Tests for create_ranked_list()."""

    def test_basic_ranking_by_stat(self, deseq2_full_df):
        ranked = create_ranked_list(deseq2_full_df, "gene_name")
        assert isinstance(ranked, pd.Series)
        assert len(ranked) <= len(deseq2_full_df)
        # Should be sorted descending
        assert ranked.iloc[0] >= ranked.iloc[-1]

    def test_ranking_by_log2fc(self, deseq2_full_df):
        ranked = create_ranked_list(
            deseq2_full_df, "gene_name", ranking_method="log2fc",
        )
        assert isinstance(ranked, pd.Series)
        assert ranked.iloc[0] >= ranked.iloc[-1]

    def test_falls_back_to_log2fc_when_stat_missing(self):
        df = pd.DataFrame({
            "gene_name": ["A", "B", "C"],
            "log2FoldChange": [2.0, -1.0, 0.5],
        })
        ranked = create_ranked_list(df, "gene_name", ranking_method="stat")
        assert len(ranked) == 3
        assert ranked.iloc[0] >= ranked.iloc[-1]

    def test_handles_duplicate_genes(self):
        df = pd.DataFrame({
            "gene_name": ["A", "A", "B"],
            "stat": [3.0, 1.0, -2.0],
            "log2FoldChange": [2.0, 0.5, -1.0],
        })
        ranked = create_ranked_list(df, "gene_name")
        assert len(ranked) == 2  # A appears once (largest abs value kept)
        assert "A" in ranked.index

    def test_nan_values_dropped(self):
        df = pd.DataFrame({
            "gene_name": ["A", "B", "C"],
            "stat": [3.0, np.nan, -2.0],
            "log2FoldChange": [2.0, np.nan, -1.0],
        })
        ranked = create_ranked_list(df, "gene_name")
        assert len(ranked) == 2
        assert "B" not in ranked.index

    def test_empty_gene_list(self):
        df = pd.DataFrame({
            "gene_name": pd.Series(dtype=str),
            "stat": pd.Series(dtype=float),
            "log2FoldChange": pd.Series(dtype=float),
        })
        ranked = create_ranked_list(df, "gene_name")
        assert len(ranked) == 0

    def test_single_gene(self):
        df = pd.DataFrame({
            "gene_name": ["MIAT"],
            "stat": [5.0],
            "log2FoldChange": [3.0],
        })
        ranked = create_ranked_list(df, "gene_name")
        assert len(ranked) == 1
        assert ranked.index[0] == "MIAT"


# ---------------------------------------------------------------------------
# normalize_gsea_cols tests
# ---------------------------------------------------------------------------


class TestNormalizeGseaCols:
    """Tests for normalize_gsea_cols()."""

    def test_normalizes_fdr_column(self):
        df = pd.DataFrame({"FDR q-val": [0.01], "NES": [2.0], "Term": ["path1"]})
        result = normalize_gsea_cols(df)
        assert "fdr" in result.columns
        assert "nes" in result.columns
        assert "Term" in result.columns

    def test_handles_different_fdr_names(self):
        for col_name in ["FDR", "fdr_bh", "padj", "adj_p_value"]:
            df = pd.DataFrame({col_name: [0.05]})
            result = normalize_gsea_cols(df)
            assert "fdr" in result.columns

    def test_derives_geneset_size_from_tag_pct(self):
        df = pd.DataFrame({
            "Tag_%": ["15/200"],
            "NES": [2.0],
        })
        result = normalize_gsea_cols(df)
        assert "geneset_size" in result.columns
        assert result["geneset_size"].iloc[0] == 200


# ---------------------------------------------------------------------------
# _db_short_label tests
# ---------------------------------------------------------------------------


class TestDbShortLabel:
    """Tests for _db_short_label()."""

    def test_known_databases(self):
        assert _db_short_label("GO_Biological_Process_2023") == "BP"
        assert _db_short_label("GO_Cellular_Component_2023") == "CC"
        assert _db_short_label("GO_Molecular_Function_2023") == "MF"
        assert _db_short_label("KEGG_2021_Human") == "KEGG"
        assert _db_short_label("Reactome_2022") == "Reactome"
        assert _db_short_label("MSigDB_Hallmark_2020") == "Hallmark"
        assert _db_short_label("WikiPathway_2021_Human") == "WikiPath"

    def test_unknown_database_truncated(self):
        result = _db_short_label("some_unknown_very_long_database")
        assert len(result) <= 8

    def test_colon_replaced_variant(self):
        result = _db_short_label("KEGG_2021_Human")
        assert result == "KEGG"


# ---------------------------------------------------------------------------
# _collect_gsea_rows tests
# ---------------------------------------------------------------------------


class TestCollectGseaRows:
    """Tests for _collect_gsea_rows()."""

    @staticmethod
    def _get_func():
        from rnaseq_explorer.engine.gsea import _collect_gsea_rows
        return _collect_gsea_rows

    def test_collects_from_memory(self, tmp_path):
        _collect_gsea_rows = self._get_func()
        gsea_results = {
            "cond1": {
                "KEGG_2021_Human": pd.DataFrame({
                    "Term": ["pathway1", "pathway2"],
                    "nes": [2.0, -1.5],
                    "fdr": [0.01, 0.05],
                    "geneset_size": [100, 200],
                }),
            }
        }
        rows = _collect_gsea_rows(gsea_results, "cond1", tmp_path, top_n=5)
        assert len(rows) == 2
        assert rows[0]["Term"] == "pathway1"
        assert "DatabaseShort" in rows[0]

    def test_empty_results(self, tmp_path):
        _collect_gsea_rows = self._get_func()
        rows = _collect_gsea_rows(None, "cond1", tmp_path)
        assert len(rows) == 0

    def test_missing_condition(self, tmp_path):
        _collect_gsea_rows = self._get_func()
        gsea_results = {"other_cond": {}}
        rows = _collect_gsea_rows(gsea_results, "cond1", tmp_path)
        assert len(rows) == 0
