"""Tests for rnaseq_explorer.viz.genewalk_viz module."""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
import pytest

from rnaseq_explorer.viz.genewalk_viz import (
    gw_domain_pie,
    gw_gene_bar,
    gw_gene_summary,
    gw_heatmap,
    gw_network,
    gw_volcano,
)


# ---------------------------------------------------------------------------
# gw_volcano
# ---------------------------------------------------------------------------

class TestGwVolcano:
    def test_returns_figure(self, sample_genewalk_df):
        fig = gw_volcano(sample_genewalk_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = gw_volcano(empty_df)
        assert isinstance(fig, go.Figure)

    def test_custom_padj_cutoff(self, sample_genewalk_df):
        fig = gw_volcano(sample_genewalk_df, padj_cutoff=0.01)
        assert isinstance(fig, go.Figure)

    def test_title_contains_genewalk(self, sample_genewalk_df):
        fig = gw_volcano(sample_genewalk_df)
        assert "GeneWalk" in fig.layout.title.text


# ---------------------------------------------------------------------------
# gw_gene_bar
# ---------------------------------------------------------------------------

class TestGwGeneBar:
    def test_returns_figure(self, sample_genewalk_df):
        fig = gw_gene_bar(sample_genewalk_df, gene="MIAT")
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = gw_gene_bar(empty_df, gene="MIAT")
        assert isinstance(fig, go.Figure)

    def test_gene_not_found(self, sample_genewalk_df):
        fig = gw_gene_bar(sample_genewalk_df, gene="NONEXISTENT")
        assert isinstance(fig, go.Figure)
        # Should show "not found" message
        annotations = [a.text for a in (fig.layout.annotations or [])]
        assert any("not found" in a.lower() or "Not Found" in a for a in annotations) or "Not Found" in (fig.layout.title.text or "")

    def test_custom_n(self, sample_genewalk_df):
        fig = gw_gene_bar(sample_genewalk_df, gene="MIAT", n=3)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# gw_network
# ---------------------------------------------------------------------------

class TestGwNetwork:
    def test_returns_figure(self, sample_genewalk_df):
        fig = gw_network(sample_genewalk_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = gw_network(empty_df)
        assert isinstance(fig, go.Figure)

    def test_with_few_nodes(self, sample_genewalk_df):
        # Use only rows for one gene
        one_gene = sample_genewalk_df[sample_genewalk_df["hgnc_symbol"] == "MIAT"]
        fig = gw_network(one_gene, padj_cutoff=1.0, min_sim=0.0)
        assert isinstance(fig, go.Figure)

    def test_strict_filter_no_results(self, sample_genewalk_df):
        fig = gw_network(sample_genewalk_df, padj_cutoff=1e-20, min_sim=0.99)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# gw_heatmap
# ---------------------------------------------------------------------------

class TestGwHeatmap:
    def test_returns_figure(self, sample_genewalk_df):
        fig = gw_heatmap(sample_genewalk_df, padj_cutoff=1.0)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = gw_heatmap(empty_df)
        assert isinstance(fig, go.Figure)

    def test_with_single_gene(self, sample_genewalk_df):
        one_gene = sample_genewalk_df[sample_genewalk_df["hgnc_symbol"] == "MIAT"]
        fig = gw_heatmap(one_gene, padj_cutoff=1.0)
        assert isinstance(fig, go.Figure)

    def test_strict_filter_empty(self, sample_genewalk_df):
        fig = gw_heatmap(sample_genewalk_df, padj_cutoff=1e-20)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# gw_domain_pie
# ---------------------------------------------------------------------------

class TestGwDomainPie:
    def test_returns_figure(self, sample_genewalk_df):
        fig = gw_domain_pie(sample_genewalk_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = gw_domain_pie(empty_df)
        assert isinstance(fig, go.Figure)

    def test_missing_domain_col(self, sample_genewalk_df):
        df = sample_genewalk_df.drop(columns=["go_domain"])
        fig = gw_domain_pie(df)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# gw_gene_summary
# ---------------------------------------------------------------------------

class TestGwGeneSummary:
    def test_returns_figure_count(self, sample_genewalk_df):
        fig = gw_gene_summary(sample_genewalk_df, padj_cutoff=1.0, metric="count")
        assert isinstance(fig, go.Figure)

    def test_returns_figure_mean_sim(self, sample_genewalk_df):
        fig = gw_gene_summary(sample_genewalk_df, padj_cutoff=1.0, metric="mean_sim")
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = gw_gene_summary(empty_df)
        assert isinstance(fig, go.Figure)

    def test_no_significant_results(self, sample_genewalk_df):
        fig = gw_gene_summary(sample_genewalk_df, padj_cutoff=1e-20)
        assert isinstance(fig, go.Figure)
