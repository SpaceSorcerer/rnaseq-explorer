"""Comprehensive tests for rnaseq_explorer.engine.qc module."""

from __future__ import annotations


import numpy as np
import pandas as pd
import pytest

from rnaseq_explorer.engine.qc import (
    compute_pca,
    compute_sample_correlation,
    load_counts_matrix,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def counts_df():
    """Normalized counts matrix: 100 genes x 6 samples."""
    np.random.seed(70)
    n_genes = 100
    n_samples = 6
    data = np.random.poisson(lam=50, size=(n_genes, n_samples)).astype(float)
    genes = [f"ENSG{i:011d}" for i in range(1, n_genes + 1)]
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    return pd.DataFrame(data, index=genes, columns=samples)


@pytest.fixture
def sample_metadata():
    """Sample metadata mapping."""
    return {
        "ctrl_1": "Control",
        "ctrl_2": "Control",
        "ctrl_3": "Control",
        "treat_1": "Treatment",
        "treat_2": "Treatment",
        "treat_3": "Treatment",
    }


# ---------------------------------------------------------------------------
# load_counts_matrix tests
# ---------------------------------------------------------------------------


class TestLoadCountsMatrix:
    """Tests for load_counts_matrix()."""

    def test_loads_csv(self, tmp_path, counts_df):
        path = tmp_path / "counts.csv"
        counts_df.to_csv(path)
        result_df, meta = load_counts_matrix(str(path))
        assert result_df is not None
        assert result_df.shape[0] == 100
        assert result_df.shape[1] == 6

    def test_loads_tsv(self, tmp_path, counts_df):
        path = tmp_path / "counts.tsv"
        counts_df.to_csv(path, sep="\t")
        result_df, meta = load_counts_matrix(str(path))
        assert result_df is not None

    def test_returns_none_for_empty_path(self):
        result_df, meta = load_counts_matrix("")
        assert result_df is None
        assert meta == {}

    def test_returns_none_for_nonexistent_file(self):
        result_df, meta = load_counts_matrix("/nonexistent/path.csv")
        assert result_df is None

    def test_autodetects_metadata_from_conditions(self, tmp_path, counts_df):
        path = tmp_path / "counts.csv"
        counts_df.to_csv(path)
        conditions = [
            {"name": "ctrl", "label": "Control"},
            {"name": "treat", "label": "Treatment"},
        ]
        result_df, meta = load_counts_matrix(str(path), conditions=conditions)
        assert result_df is not None
        # Should auto-detect some sample-condition mappings
        assert len(meta) > 0


# ---------------------------------------------------------------------------
# compute_pca tests
# ---------------------------------------------------------------------------


class TestComputePca:
    """Tests for compute_pca()."""

    def test_generates_pca_plot(self, tmp_path, counts_df, sample_metadata):
        try:
            import sklearn  # noqa: F401
            has_sklearn = True
        except ImportError:
            has_sklearn = False

        compute_pca(counts_df, metadata=sample_metadata, outdir=tmp_path)
        pca_file = tmp_path / "pca_plot.png"
        if has_sklearn:
            assert pca_file.exists()
            assert pca_file.stat().st_size > 0
        else:
            # Gracefully skipped when sklearn not installed
            assert not pca_file.exists()

    def test_no_outdir_returns_early(self, counts_df, sample_metadata):
        # Should not raise
        compute_pca(counts_df, metadata=sample_metadata, outdir=None)

    def test_insufficient_samples(self, tmp_path):
        df = pd.DataFrame({"sample1": [1, 2, 3]}, index=["g1", "g2", "g3"])
        compute_pca(df, outdir=tmp_path)
        # Should handle gracefully (no plot or info message)

    def test_precomputed_pca_file(self, tmp_path):
        pca_data = pd.DataFrame({
            "sample": ["s1", "s2", "s3"],
            "PC1": [1.0, -1.0, 0.5],
            "PC2": [0.5, -0.5, 0.0],
            "condition": ["A", "A", "B"],
        })
        pca_path = tmp_path / "pca_data.csv"
        pca_data.to_csv(pca_path, index=False)

        compute_pca(None, outdir=tmp_path, pca_file=str(pca_path))
        assert (tmp_path / "pca_plot.png").exists()


# ---------------------------------------------------------------------------
# compute_sample_correlation tests
# ---------------------------------------------------------------------------


class TestComputeSampleCorrelation:
    """Tests for compute_sample_correlation()."""

    def test_generates_correlation_heatmap(self, tmp_path, counts_df, sample_metadata):
        compute_sample_correlation(counts_df, sample_metadata, outdir=tmp_path)
        heatmap_file = tmp_path / "sample_correlation_heatmap.png"
        assert heatmap_file.exists()
        assert heatmap_file.stat().st_size > 0

    def test_insufficient_samples(self, tmp_path):
        df = pd.DataFrame({"sample1": [1, 2, 3]}, index=["g1", "g2", "g3"])
        compute_sample_correlation(df, {}, outdir=tmp_path)
        # Should handle gracefully

    def test_none_counts_df(self, tmp_path):
        compute_sample_correlation(None, {}, outdir=tmp_path)
        # Should handle gracefully
