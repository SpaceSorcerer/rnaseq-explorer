"""Tests for rnaseq_explorer.viz.rmats_viz module."""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

from rnaseq_explorer.viz.rmats_viz import (
    dpsi_distribution,
    dpsi_volcano,
    event_type_pie,
    genes_by_event_count,
    top_splicing_events,
)


# ---------------------------------------------------------------------------
# dpsi_volcano
# ---------------------------------------------------------------------------

class TestDpsiVolcano:
    def test_returns_figure(self, sample_rmats_df):
        fig = dpsi_volcano(sample_rmats_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = dpsi_volcano(empty_df)
        assert isinstance(fig, go.Figure)

    def test_with_single_event(self, sample_rmats_df):
        fig = dpsi_volcano(sample_rmats_df.head(1))
        assert isinstance(fig, go.Figure)

    def test_custom_cutoffs(self, sample_rmats_df):
        fig = dpsi_volcano(sample_rmats_df, dpsi_cutoff=0.2, fdr_cutoff=0.01)
        assert isinstance(fig, go.Figure)

    def test_single_event_type(self, sample_rmats_df):
        se_only = sample_rmats_df[sample_rmats_df["event_type"] == "SE"]
        fig = dpsi_volcano(se_only)
        assert isinstance(fig, go.Figure)
        assert len(fig.data) >= 1

    def test_title_contains_volcano(self, sample_rmats_df):
        fig = dpsi_volcano(sample_rmats_df)
        assert "Volcano" in fig.layout.title.text


# ---------------------------------------------------------------------------
# event_type_pie
# ---------------------------------------------------------------------------

class TestEventTypePie:
    def test_returns_figure(self, sample_rmats_df):
        fig = event_type_pie(sample_rmats_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = event_type_pie(empty_df)
        assert isinstance(fig, go.Figure)

    def test_missing_event_type_col(self, sample_rmats_df):
        df = sample_rmats_df.drop(columns=["event_type"])
        fig = event_type_pie(df)
        assert isinstance(fig, go.Figure)

    def test_single_event_type_pie(self):
        df = pd.DataFrame({
            "GeneID": ["A", "B", "C"],
            "event_type": ["SE", "SE", "SE"],
            "IncLevelDifference": [0.1, 0.2, 0.3],
            "FDR": [0.01, 0.02, 0.03],
        })
        fig = event_type_pie(df)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# dpsi_distribution
# ---------------------------------------------------------------------------

class TestDpsiDistribution:
    def test_returns_figure(self, sample_rmats_df):
        fig = dpsi_distribution(sample_rmats_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = dpsi_distribution(empty_df)
        assert isinstance(fig, go.Figure)

    def test_single_event_type_facet(self, sample_rmats_df):
        se_only = sample_rmats_df[sample_rmats_df["event_type"] == "SE"]
        fig = dpsi_distribution(se_only)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# top_splicing_events
# ---------------------------------------------------------------------------

class TestTopSplicingEvents:
    def test_returns_figure(self, sample_rmats_df):
        fig = top_splicing_events(sample_rmats_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = top_splicing_events(empty_df)
        assert isinstance(fig, go.Figure)

    def test_custom_n(self, sample_rmats_df):
        fig = top_splicing_events(sample_rmats_df, n=5)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# genes_by_event_count
# ---------------------------------------------------------------------------

class TestGenesByEventCount:
    def test_returns_figure(self, sample_rmats_df):
        fig = genes_by_event_count(sample_rmats_df)
        assert isinstance(fig, go.Figure)

    def test_with_empty_df(self, empty_df):
        fig = genes_by_event_count(empty_df)
        assert isinstance(fig, go.Figure)

    def test_missing_gene_col(self, sample_rmats_df):
        df = sample_rmats_df.drop(columns=["GeneID"])
        fig = genes_by_event_count(df)
        assert isinstance(fig, go.Figure)

    def test_custom_n(self, sample_rmats_df):
        fig = genes_by_event_count(sample_rmats_df, n=5)
        assert isinstance(fig, go.Figure)
