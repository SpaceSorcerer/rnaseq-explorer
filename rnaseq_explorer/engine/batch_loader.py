"""Phase 3 pipeline output directory loader.

Scans a Phase 3 output directory structure and auto-discovers
all conditions, results files, and figures for interactive exploration.
"""

from pathlib import Path
import pandas as pd
from typing import Optional


def load_phase3_output(directory: str) -> dict:
    """Scan Phase 3 output directory and return structured results.

    Parameters
    ----------
    directory : str
        Path to Phase 3 output directory (e.g., Phase3_4condition_2026_03_15_155549/)

    Returns
    -------
    dict with keys:
        conditions: list of condition names
        condition_data: {cond_name: {deseq2: DataFrame, rmats: {etype: DataFrame}, figures: [paths]}}
        cross_condition: {figures: [paths], excel_files: [paths]}
        gsea_results: {cond_name: {db: DataFrame}}
        ora_results: {method: {cond_name: {direction: DataFrame}}}
        splicing_enrichment: {cond_name: DataFrame}
        qc_plots: [paths]
        prism_files: [paths]
        pptx_path: str or None
        shortlist: DataFrame or None
        summary: dict with counts
    """
    root = Path(directory)
    if not root.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")

    result = {
        "conditions": [],
        "condition_data": {},
        "cross_condition": {"figures": [], "excel_files": []},
        "gsea_results": {},
        "ora_results": {"enrichr": {}, "gprofiler": {}},
        "splicing_enrichment": {},
        "qc_plots": [],
        "prism_files": [],
        "pptx_path": None,
        "shortlist": None,
        "summary": {},
    }

    # Known non-condition directories
    _SKIP_DIRS = {
        "cross_condition", "gsea_results", "go_ora_enrichr", "go_ora_gprofiler",
        "splicing_enrichment", "prism_files", "qc_plots",
    }

    # Discover conditions (subdirectories with deseq2_results.xlsx or rmats_results.xlsx)
    for subdir in sorted(root.iterdir()):
        if not subdir.is_dir() or subdir.name in _SKIP_DIRS:
            continue
        deseq2_file = subdir / "deseq2_results.xlsx"
        rmats_file = subdir / "rmats_results.xlsx"
        if deseq2_file.exists() or rmats_file.exists():
            cond_name = subdir.name
            result["conditions"].append(cond_name)

            cond_data = {"deseq2": None, "rmats": {}, "figures": [], "vasttools": {}}

            # Load DESeq2
            if deseq2_file.exists():
                try:
                    cond_data["deseq2"] = pd.read_excel(deseq2_file)
                except Exception as e:
                    print(f"  Warning: Could not load {deseq2_file}: {e}")

            # Load rMATS
            if rmats_file.exists():
                try:
                    xl = pd.ExcelFile(rmats_file)
                    for sheet in xl.sheet_names:
                        cond_data["rmats"][sheet] = pd.read_excel(xl, sheet_name=sheet)
                except Exception as e:
                    print(f"  Warning: Could not load {rmats_file}: {e}")

            # Collect figures
            fig_dir = subdir / "figures"
            if fig_dir.exists():
                cond_data["figures"] = sorted(fig_dir.glob("*.png")) + sorted(fig_dir.glob("*.html"))

            result["condition_data"][cond_name] = cond_data

    # Cross-condition data
    cross_dir = root / "cross_condition"
    if cross_dir.exists():
        fig_dir = cross_dir / "figures"
        if fig_dir.exists():
            result["cross_condition"]["figures"] = sorted(fig_dir.glob("*.png"))
            # Also check event_level subfolder
            event_dir = fig_dir / "event_level"
            if event_dir.exists():
                result["cross_condition"]["figures"].extend(sorted(event_dir.glob("*.png")))
        result["cross_condition"]["excel_files"] = sorted(cross_dir.glob("*.xlsx"))

        # DE-splicing shortlist
        shortlist_path = cross_dir / "de_splicing_shortlist.xlsx"
        if shortlist_path.exists():
            try:
                result["shortlist"] = pd.read_excel(shortlist_path, sheet_name=None)
            except Exception:
                pass

    # GSEA results
    gsea_dir = root / "gsea_results"
    if gsea_dir.exists():
        for cond_dir in sorted(gsea_dir.iterdir()):
            if cond_dir.is_dir() and cond_dir.name in result["conditions"]:
                cond_gsea = {}
                for csv_file in sorted(cond_dir.glob("*.csv")):
                    try:
                        cond_gsea[csv_file.stem] = pd.read_csv(csv_file)
                    except Exception:
                        pass
                if cond_gsea:
                    result["gsea_results"][cond_dir.name] = cond_gsea
        # Cross-condition GSEA summary
        summary_xlsx = gsea_dir / "gsea_cross_condition_summary.xlsx"
        if not summary_xlsx.exists():
            summary_xlsx = cross_dir / "gsea_cross_condition_summary.xlsx" if cross_dir.exists() else None

    # ORA results
    for method in ["enrichr", "gprofiler"]:
        ora_dir = root / f"go_ora_{method}"
        if ora_dir.exists():
            for f in sorted(ora_dir.glob("*.xlsx")) + sorted(ora_dir.glob("*.csv")):
                try:
                    result["ora_results"][method][f.stem] = pd.read_csv(f) if f.suffix == ".csv" else pd.read_excel(f)
                except Exception:
                    pass

    # Splicing enrichment
    spl_dir = root / "splicing_enrichment"
    if spl_dir.exists():
        for f in sorted(spl_dir.glob("**/*.xlsx")) + sorted(spl_dir.glob("**/*.csv")):
            try:
                result["splicing_enrichment"][f.stem] = pd.read_csv(f) if f.suffix == ".csv" else pd.read_excel(f)
            except Exception:
                pass

    # QC plots
    qc_dir = root / "qc_plots"
    if qc_dir.exists():
        result["qc_plots"] = sorted(qc_dir.glob("*.png"))

    # Prism files
    prism_dir = root / "prism_files"
    if prism_dir.exists():
        result["prism_files"] = sorted(prism_dir.glob("*.pzfx"))

    # PowerPoint
    for pptx in root.glob("*.pptx"):
        result["pptx_path"] = str(pptx)
        break

    # Summary
    result["summary"] = {
        "n_conditions": len(result["conditions"]),
        "n_figures": sum(len(cd["figures"]) for cd in result["condition_data"].values()) + len(result["cross_condition"]["figures"]),
        "n_prism_files": len(result["prism_files"]),
        "has_gsea": bool(result["gsea_results"]),
        "has_ora": any(result["ora_results"].values()),
        "has_shortlist": result["shortlist"] is not None,
        "has_pptx": result["pptx_path"] is not None,
        "directory": str(root),
    }

    # DEG counts per condition
    for cond_name, cd in result["condition_data"].items():
        if cd["deseq2"] is not None:
            df = cd["deseq2"]
            result["summary"][f"deg_count_{cond_name}"] = len(df)

    return result
