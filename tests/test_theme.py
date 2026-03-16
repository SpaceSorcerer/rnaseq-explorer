"""Tests for rnaseq_explorer.viz.theme module."""

from __future__ import annotations


from rnaseq_explorer.viz.theme import (
    BIOTYPE_COLORS,
    CATEGORY_COLORS,
    CONDITION_COLORS,
    EVENT_COLORS,
    PALETTE,
    condition_color_map,
    grid_dims,
    setup_plotly_theme,
)


class TestSetupPlotlyTheme:
    def test_runs_without_error(self):
        setup_plotly_theme()

    def test_dark_mode(self):
        setup_plotly_theme(dark_mode=True)
        # Reset to light
        setup_plotly_theme(dark_mode=False)


class TestPalette:
    def test_has_expected_keys(self):
        expected = {"up", "down", "neutral", "highlight", "accent1", "accent2", "accent3", "accent4"}
        assert expected.issubset(set(PALETTE.keys()))

    def test_values_are_hex(self):
        for key, value in PALETTE.items():
            assert value.startswith("#"), f"PALETTE['{key}'] is not a hex color: {value}"


class TestConditionColors:
    def test_is_non_empty(self):
        assert len(CONDITION_COLORS) > 0

    def test_all_hex(self):
        for c in CONDITION_COLORS:
            assert c.startswith("#")


class TestConditionColorMap:
    def test_maps_labels(self):
        labels = ["Control", "Treatment"]
        cmap = condition_color_map(labels)
        assert "Control" in cmap
        assert "Treatment" in cmap
        assert cmap["Control"] != cmap["Treatment"]


class TestGridDims:
    def test_one_panel(self):
        assert grid_dims(1) == (1, 1)

    def test_four_panels(self):
        assert grid_dims(4) == (2, 2)

    def test_six_panels(self):
        nrows, ncols = grid_dims(6)
        assert nrows * ncols >= 6


class TestColorDicts:
    def test_event_colors_has_five_types(self):
        assert len(EVENT_COLORS) == 5
        for key in ["SE", "A3SS", "A5SS", "RI", "MXE"]:
            assert key in EVENT_COLORS

    def test_biotype_colors_non_empty(self):
        assert len(BIOTYPE_COLORS) > 0

    def test_category_colors_non_empty(self):
        assert len(CATEGORY_COLORS) > 0
