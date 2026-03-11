"""Per-gene evidence aggregation and visual summary.

Collects evidence for a specific gene across DESeq2, GSEA, ORA, rMATS,
and GeneWalk results and produces a visual evidence card.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from rnaseq_explorer.viz.theme import (
    PALETTE,
    CONDITION_COLORS,
    FONT_SIZE_TITLE,
    FONT_SIZE_ANNOTATION,
    setup_plotly_theme,
)


def _find_gene(
    df: pd.DataFrame,
    gene_name: str,
    gene_col_candidates: Sequence[str],
) -> pd.DataFrame:
    """Search for a gene across multiple possible column names.

    Parameters
    ----------
    df : pd.DataFrame
        Data to search.
    gene_name : str
        Gene name/symbol to find (case-insensitive).
    gene_col_candidates : Sequence[str]
        Candidate column names to check.

    Returns
    -------
    pd.DataFrame
        Matching rows.
    """
    if df is None or df.empty:
        return pd.DataFrame()

    gene_upper = gene_name.upper()
    for col in gene_col_candidates:
        if col in df.columns:
            matches = df[df[col].astype(str).str.upper() == gene_upper]
            if not matches.empty:
                return matches
    return pd.DataFrame()


def investigate_gene(
    gene_name: str,
    deseq2_results: Optional[pd.DataFrame] = None,
    gsea_results: Optional[pd.DataFrame] = None,
    ora_results: Optional[pd.DataFrame] = None,
    rmats_results: Optional[pd.DataFrame] = None,
    genewalk_results: Optional[pd.DataFrame] = None,
) -> dict:
    """Aggregate all evidence for a single gene across analysis types.

    Parameters
    ----------
    gene_name : str
        Gene symbol to investigate.
    deseq2_results : pd.DataFrame, optional
        DESeq2 differential expression results.
    gsea_results : pd.DataFrame, optional
        GSEA enrichment results.
    ora_results : pd.DataFrame, optional
        ORA enrichment results.
    rmats_results : pd.DataFrame, optional
        rMATS alternative splicing results.
    genewalk_results : pd.DataFrame, optional
        GeneWalk functional annotation results.

    Returns
    -------
    dict
        Evidence dictionary with keys: "deg", "gsea", "ora", "splicing", "genewalk".
        Each value contains relevant data for the gene.
    """
    evidence: dict = {
        "gene": gene_name,
        "deg": {},
        "gsea": [],
        "ora": [],
        "splicing": [],
        "genewalk": [],
    }

    # DESeq2 evidence
    if deseq2_results is not None and not deseq2_results.empty:
        deg_match = _find_gene(
            deseq2_results, gene_name,
            ["gene_name", "Gene", "gene", "hgnc_symbol", "GeneID", "gene_id"],
        )
        if not deg_match.empty:
            row = deg_match.iloc[0]
            evidence["deg"] = {
                "log2fc": row.get("log2FoldChange", row.get("log2fc", None)),
                "padj": row.get("padj", row.get("pvalue", None)),
                "basemean": row.get("baseMean", row.get("basemean", None)),
                "direction": row.get("direction", None),
                "biotype": row.get("biotype_group", row.get("biotype", None)),
            }

    # GSEA evidence — pathways containing this gene
    if gsea_results is not None and not gsea_results.empty:
        lead_cols = ["Lead_genes", "lead_genes", "genes", "Gene"]
        for col in lead_cols:
            if col in gsea_results.columns:
                mask = gsea_results[col].astype(str).str.contains(
                    gene_name, case=False, na=False
                )
                matches = gsea_results[mask]
                for _, row in matches.head(10).iterrows():
                    term_col = next(
                        (c for c in ["Term", "term", "Name", "pathway"] if c in row.index), None
                    )
                    nes_col = next((c for c in ["NES", "nes"] if c in row.index), None)
                    evidence["gsea"].append({
                        "pathway": row[term_col] if term_col else "Unknown",
                        "nes": row[nes_col] if nes_col else None,
                        "fdr": row.get("FDR q-val", row.get("fdr", row.get("FDR", None))),
                    })
                break

    # ORA evidence
    if ora_results is not None and not ora_results.empty:
        gene_cols = ["Genes", "genes", "Overlap"]
        for col in gene_cols:
            if col in ora_results.columns:
                mask = ora_results[col].astype(str).str.contains(
                    gene_name, case=False, na=False
                )
                matches = ora_results[mask]
                for _, row in matches.head(10).iterrows():
                    term_col = next(
                        (c for c in ["Term", "term", "Name"] if c in row.index), None
                    )
                    evidence["ora"].append({
                        "term": row[term_col] if term_col else "Unknown",
                        "padj": row.get("Adjusted P-value", row.get("padj", None)),
                        "score": row.get("Combined Score", row.get("combined_score", None)),
                    })
                break

    # Splicing evidence
    if rmats_results is not None and not rmats_results.empty:
        splice_match = _find_gene(
            rmats_results, gene_name,
            ["GeneID", "geneSymbol", "gene_name", "Gene", "gene"],
        )
        for _, row in splice_match.head(10).iterrows():
            evidence["splicing"].append({
                "event_type": row.get("event_type", "Unknown"),
                "dpsi": row.get("IncLevelDifference", None),
                "fdr": row.get("FDR", None),
            })

    # GeneWalk evidence
    if genewalk_results is not None and not genewalk_results.empty:
        gw_match = _find_gene(
            genewalk_results, gene_name,
            ["hgnc_symbol", "gene", "gene_name", "GeneID"],
        )
        for _, row in gw_match.head(15).iterrows():
            evidence["genewalk"].append({
                "go_term": row.get("go_name", row.get("GO_name", "Unknown")),
                "similarity": row.get("sim", None),
                "padj": row.get("gene_padj", None),
                "domain": row.get("go_domain", None),
            })

    return evidence


def gene_evidence_card(
    evidence_dict: dict,
) -> tuple[list[go.Figure], str]:
    """Generate a visual evidence card for a single gene.

    Parameters
    ----------
    evidence_dict : dict
        Output from investigate_gene().

    Returns
    -------
    tuple[list[go.Figure], str]
        List of Plotly figures and a plain-text summary string.
    """
    setup_plotly_theme()

    gene = evidence_dict.get("gene", "Unknown")
    figures: list[go.Figure] = []
    summary_parts: list[str] = [f"Evidence summary for {gene}:"]

    # 1. DEG summary
    deg = evidence_dict.get("deg", {})
    if deg:
        log2fc = deg.get("log2fc")
        padj = deg.get("padj")
        direction = deg.get("direction", "")
        basemean = deg.get("basemean")
        biotype = deg.get("biotype", "")

        parts = []
        if log2fc is not None:
            parts.append(f"log2FC={log2fc:.3f}")
        if padj is not None:
            parts.append(f"padj={padj:.2e}")
        if direction:
            parts.append(f"direction={direction}")
        if basemean is not None:
            parts.append(f"baseMean={basemean:.1f}")
        if biotype:
            parts.append(f"biotype={biotype}")
        summary_parts.append(f"  DEG: {', '.join(parts)}")
    else:
        summary_parts.append("  DEG: No differential expression data found.")

    # 2. GSEA pathways bar
    gsea_list = evidence_dict.get("gsea", [])
    if gsea_list:
        pathways = [g["pathway"][:50] for g in gsea_list]
        nes_vals = [g.get("nes", 0) or 0 for g in gsea_list]
        colors = [PALETTE["up"] if n > 0 else PALETTE["down"] for n in nes_vals]

        fig_gsea = go.Figure(
            go.Bar(
                y=pathways,
                x=nes_vals,
                orientation="h",
                marker_color=colors,
                hovertemplate="<b>%{y}</b><br>NES: %{x:.3f}<extra></extra>",
            )
        )
        fig_gsea.update_layout(
            title=f"{gene} — GSEA Pathways ({len(gsea_list)})",
            xaxis_title="NES",
            height=max(250, len(gsea_list) * 28),
        )
        figures.append(fig_gsea)
        summary_parts.append(f"  GSEA: Found in {len(gsea_list)} pathway(s).")
    else:
        summary_parts.append("  GSEA: Not found in any pathway leading edges.")

    # 3. Splicing events
    splice_list = evidence_dict.get("splicing", [])
    if splice_list:
        event_types = [s.get("event_type", "?") for s in splice_list]
        dpsi_vals = [s.get("dpsi", 0) or 0 for s in splice_list]
        labels = [f"{et} ({d:.3f})" for et, d in zip(event_types, dpsi_vals)]

        from rnaseq_explorer.viz.theme import EVENT_COLORS

        bar_colors = [EVENT_COLORS.get(et, PALETTE["neutral"]) for et in event_types]

        fig_splice = go.Figure(
            go.Bar(
                y=labels,
                x=dpsi_vals,
                orientation="h",
                marker_color=bar_colors,
                hovertemplate="<b>%{y}</b><br>dPSI: %{x:.3f}<extra></extra>",
            )
        )
        fig_splice.update_layout(
            title=f"{gene} — Splicing Events ({len(splice_list)})",
            xaxis_title="ΔPSI",
            height=max(200, len(splice_list) * 30),
        )
        figures.append(fig_splice)
        summary_parts.append(f"  Splicing: {len(splice_list)} event(s) found.")
    else:
        summary_parts.append("  Splicing: No alternative splicing events found.")

    # 4. GeneWalk GO terms
    gw_list = evidence_dict.get("genewalk", [])
    if gw_list:
        go_terms = [g.get("go_term", "?")[:45] for g in gw_list]
        sim_vals = [g.get("similarity", 0) or 0 for g in gw_list]

        fig_gw = go.Figure(
            go.Bar(
                y=go_terms,
                x=sim_vals,
                orientation="h",
                marker_color=PALETTE["accent2"],
                hovertemplate="<b>%{y}</b><br>Similarity: %{x:.3f}<extra></extra>",
            )
        )
        fig_gw.update_layout(
            title=f"{gene} — GeneWalk GO Terms ({len(gw_list)})",
            xaxis_title="Similarity Score",
            height=max(250, len(gw_list) * 22),
        )
        figures.append(fig_gw)
        summary_parts.append(f"  GeneWalk: {len(gw_list)} GO term association(s).")
    else:
        summary_parts.append("  GeneWalk: No functional annotations found.")

    # 5. ORA terms
    ora_list = evidence_dict.get("ora", [])
    if ora_list:
        summary_parts.append(f"  ORA: Found in {len(ora_list)} enriched term(s).")
    else:
        summary_parts.append("  ORA: Not found in enriched terms.")

    summary_text = "\n".join(summary_parts)
    return figures, summary_text
