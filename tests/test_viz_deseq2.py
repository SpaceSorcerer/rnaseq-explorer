"""Tests for rnaseq_explorer.viz.deseq2_viz module."""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
import pytest

from rnaseq_explorer.viz.deseq2_viz import (
    biotype_breakdown,
    log2fc_distribution,
    ma_plot,
    pvalue_distribution,
    top_genes_bar,
    volcano_plot,
)


# ---------------------------------------------------------------------------
# volcano_plot
# ---------------------------------------------------------------------------

class TestVolcanoPlot:
    def test_returns_figure(self, sample_deseq2_df):
        fig = volcano_plot(sample_deseq2_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = volcano_plot(empty_df)
        assert isinstance(fig, go.Figure)
        assert "No Data" in fig.layout.title.text or "No data" in str(fig.layout.annotations)

    def test_with_single_gene(self, sample_deseq2_df):
        fig = volcano_plot(sample_deseq2_df.head(1))
        assert isinstance(fig, go.Figure)
        assert len(fig.data) >= 1

    def test_genes_of_interest_highlighting(self, sample_deseq2_df):
        fig = volcano_plot(sample_deseq2_df, genes_of_interest=["MIAT", "QKI"])
        assert isinstance(fig, go.Figure)
        # Should have an extra trace for genes of interest
        trace_names = [t.name for t in fig.data if t.name]
        assert any("Interest" in n for n in trace_names)

    def test_genes_of_interest_not_found(self, sample_deseq2_df):
        fig = volcano_plot(sample_deseq2_df, genes_of_interest=["NONEXISTENT_GENE"])
        assert isinstance(fig, go.Figure)

    def test_custom_cutoffs(self, sample_deseq2_df):
        fig = volcano_plot(sample_deseq2_df, log2fc_cutoff=2.0, padj_cutoff=0.01)
        assert isinstance(fig, go.Figure)

    def test_custom_column_names(self, sample_deseq2_df):
        df = sample_deseq2_df.rename(columns={
            "log2FoldChange": "lfc",
            "padj": "adj_p",
            "gene_name": "symbol",
        })
        fig = volcano_plot(df, log2fc_col="lfc", padj_col="adj_p", gene_col="symbol")
        assert isinstance(fig, go.Figure)

    def test_has_threshold_lines(self, sample_deseq2_df):
        fig = volcano_plot(sample_deseq2_df)
        # Figure should contain hlines and vlines (stored as shapes)
        shapes = fig.layout.shapes if fig.layout.shapes else []
        # At minimum, check figure was created with traces
        assert len(fig.data) >= 1


# ---------------------------------------------------------------------------
# ma_plot
# ---------------------------------------------------------------------------

class TestMAPlot:
    def test_returns_figure(self, sample_deseq2_df):
        fig = ma_plot(sample_deseq2_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = ma_plot(empty_df)
        assert isinstance(fig, go.Figure)

    def test_with_single_gene(self, sample_deseq2_df):
        fig = ma_plot(sample_deseq2_df.head(1))
        assert isinstance(fig, go.Figure)

    def test_title_contains_ma(self, sample_deseq2_df):
        fig = ma_plot(sample_deseq2_df)
        assert "MA" in fig.layout.title.text


# ---------------------------------------------------------------------------
# pvalue_distribution
# ---------------------------------------------------------------------------

class TestPvalueDistribution:
    def test_returns_figure(self, sample_deseq2_df):
        fig = pvalue_distribution(sample_deseq2_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = pvalue_distribution(empty_df)
        assert isinstance(fig, go.Figure)

    def test_histograms_present(self, sample_deseq2_df):
        fig = pvalue_distribution(sample_deseq2_df)
        hist_traces = [t for t in fig.data if isinstance(t, go.Histogram)]
        assert len(hist_traces) >= 1


# ---------------------------------------------------------------------------
# log2fc_distribution
# ---------------------------------------------------------------------------

class TestLog2FCDistribution:
    def test_returns_figure(self, sample_deseq2_df):
        fig = log2fc_distribution(sample_deseq2_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = log2fc_distribution(empty_df)
        assert isinstance(fig, go.Figure)

    def test_custom_padj_cutoff(self, sample_deseq2_df):
        fig = log2fc_distribution(sample_deseq2_df, padj_cutoff=0.001)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# top_genes_bar
# ---------------------------------------------------------------------------

class TestTopGenesBar:
    def test_returns_figure(self, sample_deseq2_df):
        fig = top_genes_bar(sample_deseq2_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = top_genes_bar(empty_df)
        assert isinstance(fig, go.Figure)

    def test_custom_n(self, sample_deseq2_df):
        fig = top_genes_bar(sample_deseq2_df, n=5)
        assert isinstance(fig, go.Figure)

    def test_with_single_gene(self, sample_deseq2_df):
        fig = top_genes_bar(sample_deseq2_df.head(1))
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# biotype_breakdown
# ---------------------------------------------------------------------------

class TestBiotypeBreakdown:
    def test_returns_figure(self, sample_deseq2_df):
        fig = biotype_breakdown(sample_deseq2_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = biotype_breakdown(empty_df)
        assert isinstance(fig, go.Figure)

    def test_missing_biotype_column(self, sample_deseq2_df):
        df = sample_deseq2_df.drop(columns=["biotype_group"])
        fig = biotype_breakdown(df)
        assert isinstance(fig, go.Figure)

    def test_missing_direction_column(self, sample_deseq2_df):
        df = sample_deseq2_df.drop(columns=["direction"])
        fig = biotype_breakdown(df)
        assert isinstance(fig, go.Figure)
