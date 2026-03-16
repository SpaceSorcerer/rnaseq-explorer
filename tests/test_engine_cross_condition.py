"""Comprehensive tests for rnaseq_explorer.engine.cross_condition module."""

from __future__ import annotations

import pandas as pd
import pytest

from rnaseq_explorer.engine.cross_condition import (
    compute_concordance,
    compute_venn_data,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def two_condition_results():
    """Pipeline condition_results with two conditions sharing some genes."""
    cond1_deg = pd.DataFrame({
        "gene_name": ["A", "B", "C", "D", "E"],
        "log2FoldChange": [2.0, -1.5, 3.0, -2.0, 1.0],
        "padj": [0.001, 0.002, 0.003, 0.004, 0.005],
        "direction": ["up", "down", "up", "down", "up"],
    })
    cond2_deg = pd.DataFrame({
        "gene_name": ["A", "B", "F", "G", "C"],
        "log2FoldChange": [1.5, -2.0, 2.5, -1.0, -1.5],
        "padj": [0.01, 0.02, 0.001, 0.01, 0.03],
        "direction": ["up", "down", "up", "down", "down"],
    })
    return {
        "cond1": {
            "deseq2_raw": cond1_deg,
            "deseq2_filtered": {"all_genes": cond1_deg},
        },
        "cond2": {
            "deseq2_raw": cond2_deg,
            "deseq2_filtered": {"all_genes": cond2_deg},
        },
    }


@pytest.fixture
def three_condition_gene_sets():
    """Gene sets for three conditions."""
    return {
        "cond1": {
            "all": {"A", "B", "C", "D"},
            "up": {"A", "C"},
            "down": {"B", "D"},
        },
        "cond2": {
            "all": {"A", "B", "E", "F"},
            "up": {"A", "E"},
            "down": {"B", "F"},
        },
        "cond3": {
            "all": {"C", "D", "E", "G"},
            "up": {"C", "E"},
            "down": {"D", "G"},
        },
    }


# ---------------------------------------------------------------------------
# compute_venn_data tests
# ---------------------------------------------------------------------------


class TestComputeVennData:
    """Tests for compute_venn_data()."""

    def test_basic_venn_data(self, three_condition_gene_sets):
        labels = {"cond1": "C1", "cond2": "C2", "cond3": "C3"}
        venn = compute_venn_data(three_condition_gene_sets, labels)

        assert "all" in venn
        assert "up" in venn
        assert "down" in venn

        # Check label mapping
        assert "C1" in venn["all"]
        assert "C2" in venn["all"]
        assert "C3" in venn["all"]

    def test_correct_set_contents(self, three_condition_gene_sets):
        labels = {"cond1": "C1", "cond2": "C2", "cond3": "C3"}
        venn = compute_venn_data(three_condition_gene_sets, labels)

        assert venn["all"]["C1"] == {"A", "B", "C", "D"}
        assert venn["up"]["C2"] == {"A", "E"}

    def test_single_condition(self):
        gene_sets = {"cond1": {"all": {"A", "B"}, "up": {"A"}, "down": {"B"}}}
        labels = {"cond1": "C1"}
        venn = compute_venn_data(gene_sets, labels)
        assert len(venn["all"]["C1"]) == 2

    def test_no_overlap(self):
        gene_sets = {
            "cond1": {"all": {"A", "B"}, "up": {"A"}, "down": {"B"}},
            "cond2": {"all": {"C", "D"}, "up": {"C"}, "down": {"D"}},
        }
        labels = {"cond1": "C1", "cond2": "C2"}
        venn = compute_venn_data(gene_sets, labels)
        overlap = venn["all"]["C1"] & venn["all"]["C2"]
        assert len(overlap) == 0

    def test_complete_overlap(self):
        gene_sets = {
            "cond1": {"all": {"A", "B"}, "up": {"A"}, "down": {"B"}},
            "cond2": {"all": {"A", "B"}, "up": {"A"}, "down": {"B"}},
        }
        labels = {"cond1": "C1", "cond2": "C2"}
        venn = compute_venn_data(gene_sets, labels)
        overlap = venn["all"]["C1"] & venn["all"]["C2"]
        assert overlap == {"A", "B"}

    def test_empty_gene_sets(self):
        gene_sets = {
            "cond1": {"all": set(), "up": set(), "down": set()},
            "cond2": {"all": set(), "up": set(), "down": set()},
        }
        labels = {"cond1": "C1", "cond2": "C2"}
        venn = compute_venn_data(gene_sets, labels)
        assert len(venn["all"]["C1"]) == 0


# ---------------------------------------------------------------------------
# compute_concordance tests
# ---------------------------------------------------------------------------


class TestComputeConcordance:
    """Tests for compute_concordance()."""

    def test_basic_concordance(self, two_condition_results):
        labels = {"cond1": "C1", "cond2": "C2"}
        result = compute_concordance(two_condition_results, labels)

        assert result.shape == (2, 2)
        assert result.loc["C1", "C1"] == 1.0
        assert result.loc["C2", "C2"] == 1.0
        # Off-diagonal should be a concordance rate between 0 and 1
        assert 0.0 <= result.loc["C1", "C2"] <= 1.0

    def test_concordance_symmetry(self, two_condition_results):
        labels = {"cond1": "C1", "cond2": "C2"}
        result = compute_concordance(two_condition_results, labels)
        assert result.loc["C1", "C2"] == result.loc["C2", "C1"]

    def test_perfect_concordance(self):
        deg = pd.DataFrame({
            "gene_name": ["A", "B", "C"],
            "direction": ["up", "down", "up"],
        })
        results = {
            "cond1": {"deseq2_filtered": {"all_genes": deg.copy()}, "deseq2_raw": deg.copy()},
            "cond2": {"deseq2_filtered": {"all_genes": deg.copy()}, "deseq2_raw": deg.copy()},
        }
        labels = {"cond1": "C1", "cond2": "C2"}
        concordance = compute_concordance(results, labels)
        assert concordance.loc["C1", "C2"] == 1.0

    def test_no_shared_genes(self):
        deg1 = pd.DataFrame({"gene_name": ["A", "B"], "direction": ["up", "down"]})
        deg2 = pd.DataFrame({"gene_name": ["C", "D"], "direction": ["up", "down"]})
        results = {
            "cond1": {"deseq2_filtered": {"all_genes": deg1}, "deseq2_raw": deg1},
            "cond2": {"deseq2_filtered": {"all_genes": deg2}, "deseq2_raw": deg2},
        }
        labels = {"cond1": "C1", "cond2": "C2"}
        concordance = compute_concordance(results, labels)
        # No shared genes means 0 concordance
        assert concordance.loc["C1", "C2"] == 0.0

    def test_single_condition(self):
        deg = pd.DataFrame({"gene_name": ["A"], "direction": ["up"]})
        results = {"cond1": {"deseq2_filtered": {"all_genes": deg}, "deseq2_raw": deg}}
        labels = {"cond1": "C1"}
        concordance = compute_concordance(results, labels)
        assert concordance.shape == (1, 1)
        assert concordance.loc["C1", "C1"] == 1.0

    def test_empty_conditions(self):
        results = {
            "cond1": {"deseq2_filtered": {"all_genes": pd.DataFrame()}, "deseq2_raw": pd.DataFrame()},
            "cond2": {"deseq2_filtered": {"all_genes": pd.DataFrame()}, "deseq2_raw": pd.DataFrame()},
        }
        labels = {"cond1": "C1", "cond2": "C2"}
        concordance = compute_concordance(results, labels)
        assert concordance.shape == (2, 2)
