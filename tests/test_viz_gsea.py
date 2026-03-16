"""Tests for rnaseq_explorer.viz.gsea_viz module."""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

from rnaseq_explorer.viz.gsea_viz import (
    enrichment_comparison,
    enrichment_dot_plot,
    leading_edge_table,
    nes_bar_chart,
    ora_dot_plot,
)


# ---------------------------------------------------------------------------
# nes_bar_chart
# ---------------------------------------------------------------------------

class TestNesBarChart:
    def test_returns_figure(self, sample_gsea_results):
        fig = nes_bar_chart(sample_gsea_results)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = nes_bar_chart(empty_df)
        assert isinstance(fig, go.Figure)

    def test_custom_n(self, sample_gsea_results):
        fig = nes_bar_chart(sample_gsea_results, n=5)
        assert isinstance(fig, go.Figure)

    def test_strict_fdr_cutoff(self, sample_gsea_results):
        fig = nes_bar_chart(sample_gsea_results, fdr_cutoff=0.001)
        assert isinstance(fig, go.Figure)

    def test_no_passing_fdr_falls_back(self, sample_gsea_results):
        fig = nes_bar_chart(sample_gsea_results, fdr_cutoff=1e-20)
        assert isinstance(fig, go.Figure)
        # Should fall back to all results
        assert len(fig.data) >= 1


# ---------------------------------------------------------------------------
# enrichment_dot_plot
# ---------------------------------------------------------------------------

class TestEnrichmentDotPlot:
    def test_returns_figure(self, sample_gsea_results):
        fig = enrichment_dot_plot(sample_gsea_results)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = enrichment_dot_plot(empty_df)
        assert isinstance(fig, go.Figure)

    def test_custom_n(self, sample_gsea_results):
        fig = enrichment_dot_plot(sample_gsea_results, n=5)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# leading_edge_table
# ---------------------------------------------------------------------------

class TestLeadingEdgeTable:
    def test_returns_dataframe(self, sample_gsea_results):
        result = leading_edge_table(
            sample_gsea_results,
            pathway_name="HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        )
        assert isinstance(result, pd.DataFrame)
        assert "gene" in result.columns

    def test_with_empty_df(self, empty_df):
        result = leading_edge_table(empty_df, pathway_name="anything")
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_pathway_not_found(self, sample_gsea_results):
        result = leading_edge_table(
            sample_gsea_results,
            pathway_name="NONEXISTENT_PATHWAY",
        )
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_extracts_genes(self, sample_gsea_results):
        result = leading_edge_table(
            sample_gsea_results,
            pathway_name="HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        )
        assert len(result) == 3  # "MIAT;QKI;TP53" -> 3 genes
        assert "MIAT" in result["gene"].values


# ---------------------------------------------------------------------------
# ora_dot_plot
# ---------------------------------------------------------------------------

class TestOraDotPlot:
    def test_returns_figure(self, sample_ora_results):
        fig = ora_dot_plot(sample_ora_results)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = ora_dot_plot(empty_df)
        assert isinstance(fig, go.Figure)

    def test_strict_fdr_falls_back(self, sample_ora_results):
        fig = ora_dot_plot(sample_ora_results, fdr_cutoff=1e-20)
        assert isinstance(fig, go.Figure)

    def test_custom_n(self, sample_ora_results):
        fig = ora_dot_plot(sample_ora_results, n=3)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# enrichment_comparison
# ---------------------------------------------------------------------------

class TestEnrichmentComparison:
    def test_returns_figure(self, sample_gsea_results):
        fig = enrichment_comparison(sample_gsea_results, sample_gsea_results)
        assert isinstance(fig, go.Figure)

    def test_with_both_empty(self, empty_df):
        fig = enrichment_comparison(empty_df, empty_df)
        assert isinstance(fig, go.Figure)

    def test_with_one_empty(self, sample_gsea_results, empty_df):
        fig = enrichment_comparison(sample_gsea_results, empty_df)
        assert isinstance(fig, go.Figure)

    def test_with_other_empty(self, sample_gsea_results, empty_df):
        fig = enrichment_comparison(empty_df, sample_gsea_results)
        assert isinstance(fig, go.Figure)
