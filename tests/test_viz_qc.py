"""Comprehensive tests for rnaseq_explorer.viz.qc_viz module."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from rnaseq_explorer.viz.qc_viz import (
    correlation_heatmap,
    pca_plot,
    top_deg_heatmap,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def pca_df():
    """PCA data with sample coordinates."""
    return pd.DataFrame({
        "sample": ["s1", "s2", "s3", "s4", "s5", "s6"],
        "PC1": [2.1, 1.8, 2.3, -1.5, -2.0, -1.8],
        "PC2": [0.5, 0.3, 0.7, -0.2, -0.5, -0.3],
        "condition": ["Ctrl", "Ctrl", "Ctrl", "Treat", "Treat", "Treat"],
    })


@pytest.fixture
def corr_matrix():
    """Sample correlation matrix."""
    np.random.seed(80)
    n = 4
    data = np.random.uniform(0.85, 1.0, (n, n))
    data = (data + data.T) / 2
    np.fill_diagonal(data, 1.0)
    labels = ["s1", "s2", "s3", "s4"]
    return pd.DataFrame(data, index=labels, columns=labels)


@pytest.fixture
def expr_matrix():
    """Expression matrix for heatmap."""
    np.random.seed(81)
    genes = ["MIAT", "QKI", "TP53", "BRCA1", "MYC"]
    samples = ["s1", "s2", "s3", "s4"]
    data = np.random.poisson(50, (5, 4)).astype(float)
    return pd.DataFrame(data, index=genes, columns=samples)


# ---------------------------------------------------------------------------
# pca_plot tests
# ---------------------------------------------------------------------------


class TestPcaPlot:
    """Tests for pca_plot()."""

    def test_returns_figure(self, pca_df):
        fig = pca_plot(pca_df, var1=45.3, var2=20.1)
        assert isinstance(fig, go.Figure)

    def test_empty_dataframe(self):
        fig = pca_plot(pd.DataFrame())
        assert isinstance(fig, go.Figure)
        # Should have an annotation about no data
        assert len(fig.layout.annotations) > 0

    def test_variance_in_axis_labels(self, pca_df):
        fig = pca_plot(pca_df, var1=45.3, var2=20.1)
        assert "45.3" in fig.layout.xaxis.title.text
        assert "20.1" in fig.layout.yaxis.title.text

    def test_zero_variance(self, pca_df):
        fig = pca_plot(pca_df, var1=0.0, var2=0.0)
        assert isinstance(fig, go.Figure)
        assert fig.layout.xaxis.title.text == "PC1"

    def test_traces_per_condition(self, pca_df):
        fig = pca_plot(pca_df)
        # Should have one trace per condition (Ctrl + Treat)
        assert len(fig.data) == 2

    def test_single_condition(self):
        df = pd.DataFrame({
            "sample": ["s1", "s2"],
            "PC1": [1.0, -1.0],
            "PC2": [0.5, -0.5],
            "condition": ["A", "A"],
        })
        fig = pca_plot(df)
        assert isinstance(fig, go.Figure)
        assert len(fig.data) == 1


# ---------------------------------------------------------------------------
# correlation_heatmap tests
# ---------------------------------------------------------------------------


class TestCorrelationHeatmap:
    """Tests for correlation_heatmap()."""

    def test_returns_figure(self, corr_matrix):
        fig = correlation_heatmap(corr_matrix)
        assert isinstance(fig, go.Figure)

    def test_empty_matrix(self):
        fig = correlation_heatmap(pd.DataFrame())
        assert isinstance(fig, go.Figure)
        assert len(fig.layout.annotations) > 0

    def test_custom_labels(self, corr_matrix):
        fig = correlation_heatmap(corr_matrix, sample_labels=["A", "B", "C", "D"])
        assert isinstance(fig, go.Figure)

    def test_single_sample(self):
        corr = pd.DataFrame([[1.0]], index=["s1"], columns=["s1"])
        fig = correlation_heatmap(corr)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# top_deg_heatmap tests
# ---------------------------------------------------------------------------


class TestTopDegHeatmap:
    """Tests for top_deg_heatmap()."""

    def test_returns_figure(self, expr_matrix):
        fig = top_deg_heatmap(expr_matrix, ["MIAT", "QKI", "TP53"])
        assert isinstance(fig, go.Figure)

    def test_empty_matrix(self):
        fig = top_deg_heatmap(pd.DataFrame(), ["A", "B"])
        assert isinstance(fig, go.Figure)
        assert len(fig.layout.annotations) > 0

    def test_empty_gene_list(self, expr_matrix):
        fig = top_deg_heatmap(expr_matrix, [])
        assert isinstance(fig, go.Figure)
        assert len(fig.layout.annotations) > 0

    def test_no_matching_genes(self, expr_matrix):
        fig = top_deg_heatmap(expr_matrix, ["NONEXISTENT1", "NONEXISTENT2"])
        assert isinstance(fig, go.Figure)
        assert len(fig.layout.annotations) > 0

    def test_partial_gene_match(self, expr_matrix):
        fig = top_deg_heatmap(expr_matrix, ["MIAT", "NONEXISTENT"])
        assert isinstance(fig, go.Figure)
        # Should still create a heatmap with the matched gene

    def test_custom_sample_labels(self, expr_matrix):
        fig = top_deg_heatmap(
            expr_matrix, ["MIAT", "QKI"],
            sample_labels=["S1", "S2", "S3", "S4"],
        )
        assert isinstance(fig, go.Figure)
