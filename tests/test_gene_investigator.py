"""Tests for rnaseq_explorer.viz.gene_investigator module."""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
import pytest

from rnaseq_explorer.viz.gene_investigator import (
    gene_evidence_card,
    investigate_gene,
)


# ---------------------------------------------------------------------------
# investigate_gene
# ---------------------------------------------------------------------------

class TestInvestigateGene:
    def test_with_all_sources(
        self,
        sample_deseq2_df,
        sample_gsea_results,
        sample_ora_results,
        sample_rmats_df,
        sample_genewalk_df,
    ):
        evidence = investigate_gene(
            "MIAT",
            deseq2_results=sample_deseq2_df,
            gsea_results=sample_gsea_results,
            ora_results=sample_ora_results,
            rmats_results=sample_rmats_df,
            genewalk_results=sample_genewalk_df,
        )
        assert isinstance(evidence, dict)
        assert evidence["gene"] == "MIAT"
        assert "deg" in evidence
        assert "gsea" in evidence
        assert "ora" in evidence
        assert "splicing" in evidence
        assert "genewalk" in evidence

    def test_deg_data_found(self, sample_deseq2_df):
        evidence = investigate_gene("MIAT", deseq2_results=sample_deseq2_df)
        assert evidence["deg"]  # non-empty dict
        assert "log2fc" in evidence["deg"]

    def test_with_only_deseq2(self, sample_deseq2_df):
        evidence = investigate_gene("MIAT", deseq2_results=sample_deseq2_df)
        assert isinstance(evidence, dict)
        assert evidence["deg"]
        assert evidence["gsea"] == []
        assert evidence["splicing"] == []

    def test_with_partial_sources_no_genewalk(
        self,
        sample_deseq2_df,
        sample_rmats_df,
    ):
        evidence = investigate_gene(
            "MIAT",
            deseq2_results=sample_deseq2_df,
            rmats_results=sample_rmats_df,
        )
        assert isinstance(evidence, dict)
        assert evidence["genewalk"] == []

    def test_gene_not_found_in_any(self, sample_deseq2_df, sample_rmats_df):
        evidence = investigate_gene(
            "ZZZZZ_FAKE_GENE",
            deseq2_results=sample_deseq2_df,
            rmats_results=sample_rmats_df,
        )
        assert evidence["deg"] == {}
        assert evidence["splicing"] == []

    def test_with_no_sources(self):
        evidence = investigate_gene("MIAT")
        assert evidence["gene"] == "MIAT"
        assert evidence["deg"] == {}
        assert evidence["gsea"] == []
        assert evidence["ora"] == []
        assert evidence["splicing"] == []
        assert evidence["genewalk"] == []

    def test_splicing_found(self, sample_rmats_df):
        evidence = investigate_gene("MIAT", rmats_results=sample_rmats_df)
        assert len(evidence["splicing"]) >= 1
        assert "event_type" in evidence["splicing"][0]

    def test_genewalk_found(self, sample_genewalk_df):
        evidence = investigate_gene("MIAT", genewalk_results=sample_genewalk_df)
        assert len(evidence["genewalk"]) >= 1
        assert "go_term" in evidence["genewalk"][0]


# ---------------------------------------------------------------------------
# gene_evidence_card
# ---------------------------------------------------------------------------

class TestGeneEvidenceCard:
    def test_returns_figures_and_summary(
        self,
        sample_deseq2_df,
        sample_gsea_results,
        sample_rmats_df,
        sample_genewalk_df,
    ):
        evidence = investigate_gene(
            "MIAT",
            deseq2_results=sample_deseq2_df,
            gsea_results=sample_gsea_results,
            rmats_results=sample_rmats_df,
            genewalk_results=sample_genewalk_df,
        )
        figures, summary = gene_evidence_card(evidence)
        assert isinstance(figures, list)
        assert isinstance(summary, str)
        assert "MIAT" in summary
        for fig in figures:
            assert isinstance(fig, go.Figure)

    def test_with_empty_evidence(self):
        evidence = investigate_gene("FAKE_GENE")
        figures, summary = gene_evidence_card(evidence)
        assert isinstance(figures, list)
        assert len(figures) == 0
        assert "FAKE_GENE" in summary
