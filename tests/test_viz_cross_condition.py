"""Comprehensive tests for rnaseq_explorer.viz.cross_condition_viz module."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from rnaseq_explorer.viz.cross_condition_viz import (
    direction_concordance_heatmap,
    log2fc_scatter,
    overlap_bar,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def concordance_matrix():
    """Pairwise concordance matrix."""
    data = np.array([
        [1.0, 0.8, 0.6],
        [0.8, 1.0, 0.7],
        [0.6, 0.7, 1.0],
    ])
    labels = ["Cond1", "Cond2", "Cond3"]
    return pd.DataFrame(data, index=labels, columns=labels)


@pytest.fixture
def deseq2_df_pair():
    """Pair of DESeq2 DataFrames for scatter comparison."""
    np.random.seed(90)
    n = 50
    genes = [f"GENE{i}" for i in range(1, n + 1)]
    df1 = pd.DataFrame({
        "gene_name": genes,
        "log2FoldChange": np.random.normal(0, 2, n),
        "padj": np.random.uniform(0, 1, n),
    })
    df2 = pd.DataFrame({
        "gene_name": genes,
        "log2FoldChange": df1["log2FoldChange"] * 0.8 + np.random.normal(0, 0.5, n),
        "padj": np.random.uniform(0, 1, n),
    })
    return df1, df2


@pytest.fixture
def overlap_data_long():
    """Overlap data in long format."""
    return pd.DataFrame({
        "condition": ["Cond1", "Cond2", "Cond3"],
        "unique": [150, 200, 100],
        "shared": [50, 50, 50],
    })


@pytest.fixture
def overlap_data_matrix():
    """Overlap data in matrix format."""
    return pd.DataFrame(
        {"Cond1 only": [100, 0], "Cond2 only": [0, 80], "Shared": [50, 50]},
        index=["Cond1", "Cond2"],
    )


# ---------------------------------------------------------------------------
# direction_concordance_heatmap tests
# ---------------------------------------------------------------------------


class TestDirectionConcordanceHeatmap:
    """Tests for direction_concordance_heatmap()."""

    def test_returns_figure(self, concordance_matrix):
        fig = direction_concordance_heatmap(concordance_matrix)
        assert isinstance(fig, go.Figure)

    def test_empty_matrix(self):
        fig = direction_concordance_heatmap(pd.DataFrame())
        assert isinstance(fig, go.Figure)
        assert len(fig.layout.annotations) > 0

    def test_custom_labels(self, concordance_matrix):
        fig = direction_concordance_heatmap(
            concordance_matrix,
            condition_labels=["A", "B", "C"],
        )
        assert isinstance(fig, go.Figure)

    def test_single_condition(self):
        matrix = pd.DataFrame([[1.0]], index=["C1"], columns=["C1"])
        fig = direction_concordance_heatmap(matrix)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# log2fc_scatter tests
# ---------------------------------------------------------------------------


class TestLog2fcScatter:
    """Tests for log2fc_scatter()."""

    def test_returns_figure(self, deseq2_df_pair):
        df1, df2 = deseq2_df_pair
        fig = log2fc_scatter(df1, df2, cond1_name="Ctrl", cond2_name="Treat")
        assert isinstance(fig, go.Figure)
        # Should have at least data traces + fit line + identity line
        assert len(fig.data) >= 2

    def test_empty_df1(self, deseq2_df_pair):
        _, df2 = deseq2_df_pair
        fig = log2fc_scatter(pd.DataFrame(), df2)
        assert isinstance(fig, go.Figure)
        assert len(fig.layout.annotations) > 0

    def test_empty_df2(self, deseq2_df_pair):
        df1, _ = deseq2_df_pair
        fig = log2fc_scatter(df1, pd.DataFrame())
        assert isinstance(fig, go.Figure)

    def test_no_overlapping_genes(self):
        df1 = pd.DataFrame({
            "gene_name": ["A", "B"],
            "log2FoldChange": [1.0, -1.0],
            "padj": [0.01, 0.02],
        })
        df2 = pd.DataFrame({
            "gene_name": ["C", "D"],
            "log2FoldChange": [2.0, -2.0],
            "padj": [0.01, 0.02],
        })
        fig = log2fc_scatter(df1, df2)
        assert isinstance(fig, go.Figure)

    def test_missing_gene_col(self):
        df1 = pd.DataFrame({"value": [1.0]})
        df2 = pd.DataFrame({"value": [2.0]})
        fig = log2fc_scatter(df1, df2)
        assert isinstance(fig, go.Figure)

    def test_without_padj_col(self, deseq2_df_pair):
        df1, df2 = deseq2_df_pair
        fig = log2fc_scatter(df1, df2, padj_col=None)
        assert isinstance(fig, go.Figure)

    def test_with_significance_coloring(self, deseq2_df_pair):
        df1, df2 = deseq2_df_pair
        fig = log2fc_scatter(df1, df2, padj_col="padj", padj_cutoff=0.05)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# overlap_bar tests
# ---------------------------------------------------------------------------


class TestOverlapBar:
    """Tests for overlap_bar()."""

    def test_long_format(self, overlap_data_long):
        fig = overlap_bar(overlap_data_long)
        assert isinstance(fig, go.Figure)

    def test_matrix_format(self, overlap_data_matrix):
        fig = overlap_bar(overlap_data_matrix)
        assert isinstance(fig, go.Figure)

    def test_empty_data(self):
        fig = overlap_bar(pd.DataFrame())
        assert isinstance(fig, go.Figure)
        assert len(fig.layout.annotations) > 0

    def test_custom_labels(self, overlap_data_long):
        fig = overlap_bar(overlap_data_long, condition_labels=["A", "B", "C"])
        assert isinstance(fig, go.Figure)

    def test_long_format_unique_only(self):
        df = pd.DataFrame({
            "condition": ["A", "B"],
            "unique": [100, 200],
        })
        fig = overlap_bar(df)
        assert isinstance(fig, go.Figure)

    def test_long_format_shared_only(self):
        df = pd.DataFrame({
            "condition": ["A", "B"],
            "shared": [50, 50],
        })
        fig = overlap_bar(df)
        assert isinstance(fig, go.Figure)
