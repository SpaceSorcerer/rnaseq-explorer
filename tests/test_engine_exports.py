"""Comprehensive tests for rnaseq_explorer.engine.exports module."""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pandas as pd
import pytest

from rnaseq_explorer.engine.exports import (
    _add_table,
    _create_prism_xml,
    _save_prism,
    export_excel,
    validate_outputs,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def deseq2_sets():
    """DESeq2 filtered data sets for export."""
    all_genes = pd.DataFrame({
        "gene_name": ["A", "B", "C", "D"],
        "log2FoldChange": [2.0, -1.5, 3.0, -2.0],
        "padj": [0.001, 0.002, 0.003, 0.004],
        "baseMean": [100, 200, 150, 300],
        "direction": ["up", "down", "up", "down"],
    })
    return {"all_genes": all_genes}


@pytest.fixture
def rmats_filtered():
    """rMATS filtered data for export."""
    se_df = pd.DataFrame({
        "ID": [1, 2, 3],
        "geneSymbol": ["X", "Y", "Z"],
        "IncLevelDifference": [0.3, -0.2, 0.15],
        "FDR": [0.001, 0.01, 0.02],
        "event_type": ["SE", "SE", "SE"],
    })
    return {"SE": se_df}


# ---------------------------------------------------------------------------
# _create_prism_xml tests
# ---------------------------------------------------------------------------


class TestCreatePrismXml:
    """Tests for _create_prism_xml()."""

    def test_creates_valid_root(self):
        root = _create_prism_xml()
        assert root.tag == "GraphPadPrismFile"
        assert root.get("PrismXMLVersion") == "5.00"

    def test_has_table_sequence(self):
        root = _create_prism_xml()
        table_seq = root.find("TableSequence")
        assert table_seq is not None

    def test_has_info_section(self):
        root = _create_prism_xml()
        info = root.find("Info")
        assert info is not None


# ---------------------------------------------------------------------------
# _add_table tests
# ---------------------------------------------------------------------------


class TestAddTable:
    """Tests for _add_table()."""

    def test_adds_table_with_data(self):
        root = _create_prism_xml()
        _add_table(root, "Test Table", [
            ("Gene", ["A", "B", "C"]),
            ("Value", [1.0, 2.0, 3.0]),
        ])
        tables = root.findall("Table")
        assert len(tables) == 1
        assert tables[0].find("Title").text == "Test Table"

    def test_adds_multiple_tables(self):
        root = _create_prism_xml()
        _add_table(root, "Table 1", [("Col", [1, 2])])
        _add_table(root, "Table 2", [("Col", [3, 4])])
        tables = root.findall("Table")
        assert len(tables) == 2

    def test_numeric_only_columns(self):
        root = _create_prism_xml()
        _add_table(root, "Numeric", [
            ("Values", [1.0, 2.0, 3.0]),
        ])
        tables = root.findall("Table")
        assert len(tables) == 1

    def test_row_titles_detection(self):
        root = _create_prism_xml()
        _add_table(root, "With Titles", [
            ("Gene", ["MIAT", "QKI", "TP53"]),
            ("FC", [2.0, -1.5, 3.0]),
        ])
        tables = root.findall("Table")
        assert len(tables) == 1
        # Should have row titles column
        row_titles = tables[0].find("RowTitlesColumn")
        assert row_titles is not None


# ---------------------------------------------------------------------------
# _save_prism tests
# ---------------------------------------------------------------------------


class TestSavePrism:
    """Tests for _save_prism()."""

    def test_saves_valid_xml(self, tmp_path):
        root = _create_prism_xml()
        _add_table(root, "Test", [("Col", [1, 2, 3])])
        filepath = tmp_path / "test.pzfx"
        _save_prism(root, filepath)

        assert filepath.exists()
        assert filepath.stat().st_size > 0

        # Should be parseable XML
        tree = ET.parse(filepath)
        assert tree.getroot().tag == "GraphPadPrismFile"

    def test_removes_internal_counter(self, tmp_path):
        root = _create_prism_xml()
        filepath = tmp_path / "test.pzfx"
        _save_prism(root, filepath)

        tree = ET.parse(filepath)
        assert tree.getroot().get("_table_count") is None


# ---------------------------------------------------------------------------
# export_excel tests
# ---------------------------------------------------------------------------


class TestExportExcel:
    """Tests for export_excel()."""

    def test_creates_deseq2_xlsx(self, tmp_path, deseq2_sets, rmats_filtered):
        export_excel(deseq2_sets, rmats_filtered, tmp_path)
        deseq2_xlsx = tmp_path / "deseq2_results.xlsx"
        assert deseq2_xlsx.exists()
        assert deseq2_xlsx.stat().st_size > 0

    def test_creates_rmats_xlsx(self, tmp_path, deseq2_sets, rmats_filtered):
        export_excel(deseq2_sets, rmats_filtered, tmp_path)
        rmats_xlsx = tmp_path / "rmats_results.xlsx"
        assert rmats_xlsx.exists()

    def test_deseq2_xlsx_readable(self, tmp_path, deseq2_sets, rmats_filtered):
        export_excel(deseq2_sets, rmats_filtered, tmp_path)
        deseq2_xlsx = tmp_path / "deseq2_results.xlsx"
        sheets = pd.ExcelFile(deseq2_xlsx).sheet_names
        assert "Summary" in sheets

    def test_empty_rmats(self, tmp_path, deseq2_sets):
        export_excel(deseq2_sets, {}, tmp_path)
        deseq2_xlsx = tmp_path / "deseq2_results.xlsx"
        assert deseq2_xlsx.exists()
        # rmats xlsx should not be created
        rmats_xlsx = tmp_path / "rmats_results.xlsx"
        assert not rmats_xlsx.exists()


# ---------------------------------------------------------------------------
# validate_outputs tests
# ---------------------------------------------------------------------------


class TestValidateOutputs:
    """Tests for validate_outputs()."""

    def test_missing_files_return_false(self, tmp_path):
        results = {
            "cond1": {
                "deseq2_raw": pd.DataFrame(),
                "deseq2_filtered": {"all_genes": pd.DataFrame()},
                "rmats_raw": {},
                "rmats_filtered": {},
            }
        }
        labels = {"cond1": "C1"}
        passed = validate_outputs(results, labels, tmp_path)
        assert passed is False

    def test_with_expected_files(self, tmp_path):
        # Create expected output structure
        cross_dir = tmp_path / "cross_condition"
        cross_dir.mkdir()
        (cross_dir / "multi_condition_results.xlsx").write_bytes(b"x" * 1000)
        (tmp_path / "RNA-seq_Analysis_Report.pptx").write_bytes(b"x" * 1000)

        prism_dir = tmp_path / "prism_files"
        prism_dir.mkdir()
        for i in range(16):
            (prism_dir / f"file{i}.pzfx").write_bytes(b"x" * 20000)

        deg_df = pd.DataFrame({
            "gene_name": ["A", "B"],
            "direction": ["up", "down"],
        })
        results = {
            "cond1": {
                "deseq2_raw": deg_df,
                "deseq2_filtered": {"all_genes": deg_df},
                "rmats_raw": {},
                "rmats_filtered": {},
            }
        }
        labels = {"cond1": "C1"}
        passed = validate_outputs(results, labels, tmp_path)
        # Should pass since expected files exist
        assert passed is True
