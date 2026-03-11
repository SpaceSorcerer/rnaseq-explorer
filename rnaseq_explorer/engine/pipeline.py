"""Backward-compatible orchestration wrapper.

Imports from all submodules and exposes a ``run_pipeline(config)`` function
that works identically to the original monolith's ``run_pipeline()``.
This ensures existing batch scripts (``run_4condition.py``, etc.) still work.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from rnaseq_explorer.engine.deseq2 import (
    DEFAULT_DESEQ2_COLS,
    BIOTYPE_GROUPS,
    BIOTYPE_ORDER,
    BIOTYPE_COLORS,
    load_file,
    validate_columns,
    normalize_deseq2_columns,
    best_gene_key,
    enrich_with_gene_names,
    reassign_biotypes_from_mygene,
    filter_deseq2,
    extract_gene_sets,
    load_rbp_annotations,
    annotate_rbps,
    # Per-condition DESeq2 visualization functions
    volcano_plot,
    ma_plot,
    volcano_plot_interactive,
    ma_plot_interactive,
    biotype_chart,
    biotype_direction_chart,
    biotype_enrichment_test,
    biotype_volcano,
    ecdf_log2fc_by_biotype,
    pvalue_histogram,
    top_genes_lollipop,
    expression_rank_plot,
    volcano_plot_labeled,
    rbp_heatmap,
    rbp_summary_table,
    # Cross-condition biotype visualization
    cross_condition_biotype_comparison,
    cross_condition_biotype_direction,
)
from rnaseq_explorer.engine.rmats import (
    RMATS_EVENT_TYPES,
    DEFAULT_RMATS_COLS,
    load_all_rmats,
    filter_rmats,
    make_event_key,
    # Per-condition rMATS visualization functions
    rmats_scatter,
    rmats_combined_volcano,
    rmats_event_summary_chart,
    rmats_dpsi_distribution,
    rmats_psi_scatter,
)
from rnaseq_explorer.engine.gsea import (
    run_gsea_enrichment,
    normalize_gsea_cols,
    create_ranked_list,
    run_gsea_prerank,
    # GSEA visualization functions
    gsea_combined_plot,
    gsea_enrichment_plots,
    export_gsea_leading_edge,
    gsea_dotplot_legacy,
)
from rnaseq_explorer.engine.ora import (
    run_gprofiler_ora,
    run_enrichr_ora,
    run_dual_ora,
    # ORA visualization functions
    go_enrichment_combined_plot,
    export_go_prism,
)
from rnaseq_explorer.engine.qc import (
    load_counts_matrix,
    compute_pca,
    compute_sample_correlation,
    compute_top_deg_heatmap,
)
from rnaseq_explorer.engine.cross_condition import (
    compute_venn_data,
    deseq2_venn_diagrams,
    compute_upset_data,
    deseq2_upset_plot,
    compute_concordance,
    compute_direction_heatmap,
    deseq2_log2fc_heatmap,
    pairwise_log2fc_scatter,
    # Additional cross-condition functions
    deseq2_de_counts_chart,
    pairwise_deg_venns,
    rmats_cross_condition_venn,
    rmats_direction_concordance,
    pairwise_splicing_venns,
    rmats_event_count_comparison,
    pairwise_dpsi_scatter,
    rmats_upset_plot,
    rmats_event_heatmap,
    rmats_event_pie_chart,
    deseq2_vs_rmats_venn,
    log2fc_vs_dpsi_scatter,
    gene_overlap_summary,
    summary_dashboard,
)
from rnaseq_explorer.engine.exports import (
    export_excel,
    export_combined_results,
    export_prism_pzfx,
    export_powerpoint,
    export_unfiltered_merged,
    export_pairwise_workbook,
    validate_outputs,
)
from rnaseq_explorer.viz.theme import (
    setup_matplotlib_style,
    COLOR_UP,
    COLOR_DOWN,
    COLOR_NS,
    EVENT_COLORS,
)


# ---------------------------------------------------------------------------
# Default configuration
# ---------------------------------------------------------------------------

DEFAULT_CONFIG: dict[str, Any] = {
    "CONDITIONS": [],
    "OUTPUT_DIR": "./output",
    "LOG2FC_CUTOFF": 1.0,
    "BASEMEAN_CUTOFF": 10,
    "PADJ_CUTOFF": 0.05,
    "AUTO_BIOTYPE_SPLIT": True,
    "GENE_NAME_LOOKUP": True,
    "SPECIES": "human",
    "RMATS_FDR_CUTOFF": 0.05,
    "RMATS_PVAL_CUTOFF": 0.05,
    "INCLEVEL_DIFF_CUTOFF": 0.1,
    "USE_FDR": True,
    "RMATS_DUAL_FILTER": False,
    "FIG_DPI": 300,
    "FIG_FORMAT": "png",
    "FONT_SIZE": 12,
    "COLOR_UP": COLOR_UP,
    "COLOR_DOWN": COLOR_DOWN,
    "COLOR_NS": COLOR_NS,
    "INTERACTIVE_PLOTS": True,
    "DESEQ2_COLS": dict(DEFAULT_DESEQ2_COLS),
    "RMATS_COLS": dict(DEFAULT_RMATS_COLS),
    "GSEA_DATABASES": [
        "GO_Biological_Process_2023",
        "GO_Cellular_Component_2023",
        "GO_Molecular_Function_2023",
        "KEGG_2021_Human",
        "Reactome_2022",
        "MSigDB_Hallmark_2020",
        "WikiPathway_2021_Human",
    ],
    "ORA_DATABASES": [
        "GO_Biological_Process_2023",
        "GO_Cellular_Component_2023",
        "GO_Molecular_Function_2023",
        "KEGG_2021_Human",
        "Reactome_2022",
    ],
    "GENES_OF_INTEREST": [],
    "COUNTS_FILE": "",
    "SAMPLE_METADATA": {},
    "GSEA_RANKING": "stat",
    "GSEA_MIN_SIZE": 15,
    "GSEA_MAX_SIZE": 500,
    "GSEA_PERMUTATIONS": 1000,
    "ORA_METHOD": "both",
    "RBP_FILE": "",
    "RMATS_FILE_SUFFIX": ".MATS.JCEC.txt",
}


def _merge_config(user_config: dict) -> dict:
    """Merge user config with defaults, applying type coercions."""
    cfg = dict(DEFAULT_CONFIG)
    for key, val in user_config.items():
        if key in cfg:
            cfg[key] = val
        else:
            cfg[key] = val
    return cfg


# ---------------------------------------------------------------------------
# Main pipeline entry point
# ---------------------------------------------------------------------------


def run_pipeline(config: dict) -> None:
    """Run the full RNA-seq analysis pipeline.

    This is the backward-compatible entry point that works identically to
    the original monolith's ``run_pipeline(config)`` function.  The call
    order matches the monolith's ``main()`` exactly.

    Parameters
    ----------
    config : dict
        Pipeline configuration. Required keys: CONDITIONS, OUTPUT_DIR.
        All other keys are optional and have sensible defaults.
    """
    cfg = _merge_config(config)

    # Extract config values
    conditions = cfg["CONDITIONS"]
    output_dir = cfg["OUTPUT_DIR"]
    log2fc_cutoff = float(cfg["LOG2FC_CUTOFF"])
    basemean_cutoff = float(cfg["BASEMEAN_CUTOFF"])
    padj_cutoff = float(cfg["PADJ_CUTOFF"])
    auto_biotype_split = bool(cfg["AUTO_BIOTYPE_SPLIT"])
    gene_name_lookup = bool(cfg["GENE_NAME_LOOKUP"])
    species = str(cfg["SPECIES"])
    fdr_cutoff = float(cfg["RMATS_FDR_CUTOFF"])
    pval_cutoff = float(cfg["RMATS_PVAL_CUTOFF"])
    dpsi_cutoff = float(cfg["INCLEVEL_DIFF_CUTOFF"])
    use_fdr = bool(cfg["USE_FDR"])
    dual_filter = bool(cfg["RMATS_DUAL_FILTER"])
    fig_dpi = int(cfg["FIG_DPI"])
    fig_format = str(cfg["FIG_FORMAT"])
    font_size = int(cfg["FONT_SIZE"])
    interactive_plots = bool(cfg["INTERACTIVE_PLOTS"])
    deseq2_cols = dict(cfg["DESEQ2_COLS"])
    rmats_cols = dict(cfg["RMATS_COLS"])
    gsea_databases = list(cfg["GSEA_DATABASES"])
    gsea_ranking = str(cfg["GSEA_RANKING"])
    gsea_min_size = int(cfg["GSEA_MIN_SIZE"])
    gsea_max_size = int(cfg["GSEA_MAX_SIZE"])
    gsea_permutations = int(cfg["GSEA_PERMUTATIONS"])
    ora_method = str(cfg["ORA_METHOD"])
    counts_file = str(cfg["COUNTS_FILE"])
    sample_metadata_cfg = dict(cfg["SAMPLE_METADATA"])
    rbp_file = str(cfg["RBP_FILE"])
    rmats_file_suffix = str(cfg["RMATS_FILE_SUFFIX"])
    genes_of_interest = list(cfg.get("GENES_OF_INTEREST", []))

    if not conditions:
        raise ValueError("No conditions defined in CONDITIONS list")

    # Setup style
    setup_matplotlib_style(dpi=fig_dpi, font_size=font_size)

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  RNA-seq Explorer Pipeline (modular engine)")
    print(f"  {len(conditions)} conditions to process")
    print("=" * 60)

    # ===================================================================
    # PHASE 1: Per-condition analysis
    # ===================================================================
    condition_results: dict[str, dict] = {}
    condition_labels: dict[str, str] = {}
    counts_df = None
    sample_metadata: dict[str, str] = {}

    # Load normalized counts matrix if provided (for PCA, heatmaps)
    if counts_file:
        counts_df, sample_metadata = load_counts_matrix(
            counts_file, sample_metadata_cfg or None, conditions
        )
        if counts_df is not None:
            qc_dir = outdir / "qc_plots"
            qc_dir.mkdir(parents=True, exist_ok=True)
            print("\n-- QC: PCA Plot --")
            compute_pca(
                counts_df=counts_df, metadata=sample_metadata,
                outdir=qc_dir, fig_format=fig_format, fig_dpi=fig_dpi,
            )
            print("\n-- QC: Sample Correlation Heatmap --")
            compute_sample_correlation(
                counts_df, sample_metadata, qc_dir,
                fig_format=fig_format, fig_dpi=fig_dpi,
            )
    else:
        print("[INFO] No counts file -- skipping PCA, correlation heatmap, top DEG heatmap")

    # Load RBP annotations if provided
    rbp_annotations: dict = {}
    if rbp_file:
        print("\n-- Loading RBP Annotations --")
        try:
            rbp_annotations = load_rbp_annotations(rbp_file)
        except Exception as e:
            print(f"  WARNING: Failed to load RBP annotations: {e}")
    else:
        print("[INFO] No RBP_FILE -- skipping RBP annotation")

    # Per-condition PCA files (from WSF output, if provided)
    for cond in conditions:
        if cond.get("pca_file"):
            cond_fig = outdir / cond["name"] / "figures"
            cond_fig.mkdir(parents=True, exist_ok=True)
            compute_pca(pca_file=cond["pca_file"], outdir=cond_fig,
                        fig_format=fig_format, fig_dpi=fig_dpi)

    for cond in conditions:
        cond_name = cond["name"]
        cond_label = cond["label"]
        condition_labels[cond_name] = cond_label

        print(f"\n{'=' * 60}")
        print(f"  Processing condition: {cond_label}")
        print(f"{'=' * 60}")

        cond_outdir = outdir / cond_name
        cond_outdir.mkdir(parents=True, exist_ok=True)
        cond_fig_dir = cond_outdir / "figures"
        cond_fig_dir.mkdir(exist_ok=True)

        # --- DESeq2 ---
        print("\n-- Loading DESeq2 --")
        deseq2_raw = load_file(cond["deseq2_file"], f"DESeq2 ({cond_label})")
        deseq2_raw = normalize_deseq2_columns(
            deseq2_raw, deseq2_cols, f"DESeq2 ({cond_label})"
        )

        optional_keys = {"biotype", "gene_name", "stat", "lfcSE"}
        required_deseq2 = [
            v for k, v in deseq2_cols.items()
            if v is not None and k not in optional_keys
        ]
        validate_columns(deseq2_raw, required_deseq2, f"DESeq2 ({cond_label})")

        deseq2_raw = enrich_with_gene_names(
            deseq2_raw, deseq2_cols, species=species,
            lookup_enabled=gene_name_lookup, file_label=f"DESeq2 ({cond_label})",
        )
        deseq2_raw = reassign_biotypes_from_mygene(
            deseq2_raw, deseq2_cols, species=species,
            file_label=f"DESeq2 ({cond_label})",
        )

        # RBP annotation on raw data (so columns propagate to all exports)
        if rbp_annotations:
            gene_col_rbp = deseq2_cols.get("gene_name", "")
            if gene_col_rbp and gene_col_rbp in deseq2_raw.columns:
                print(f"\n-- RBP Annotation [{cond_label}] --")
                deseq2_raw = annotate_rbps(
                    deseq2_raw, rbp_annotations, gene_col=gene_col_rbp
                )

        # Biotype passes
        biotype_passes = [("All Genes", None, "")]
        if auto_biotype_split and deseq2_cols["biotype"] in deseq2_raw.columns:
            biotype_passes += [
                ("Protein Coding", "protein_coding", "_protein_coding"),
                ("Non-Protein Coding", "non_protein_coding", "_non_protein_coding"),
            ]

        deseq2_filtered_sets: dict[str, Any] = {}

        for label, bio_filter, suffix in biotype_passes:
            full_label = f"{cond_label} - {label}"
            deseq2_all, deseq2_filt = filter_deseq2(
                deseq2_raw, deseq2_cols,
                log2fc_cutoff=log2fc_cutoff,
                basemean_cutoff=basemean_cutoff,
                padj_cutoff=padj_cutoff,
                biotype_filter=bio_filter,
                label=full_label,
            )
            export_key = label.lower().replace(" ", "_").replace("-", "_")
            deseq2_filtered_sets[export_key] = deseq2_filt

            # -- Per-condition DESeq2 visualization --
            print(f"\n-- Generating DESeq2 Figures [{full_label}] --")
            pvalue_histogram(
                deseq2_all, cond_fig_dir, cols=deseq2_cols,
                padj_cutoff=padj_cutoff, label=full_label, suffix=suffix,
                fig_format=fig_format,
            )
            volcano_plot(
                deseq2_all, cond_fig_dir, cols=deseq2_cols,
                padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff,
                basemean_cutoff=basemean_cutoff, label=full_label,
                suffix=suffix, fig_format=fig_format, fig_dpi=fig_dpi,
            )
            volcano_plot_labeled(
                deseq2_all, cond_fig_dir, cols=deseq2_cols,
                padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff,
                basemean_cutoff=basemean_cutoff, label=full_label,
                suffix=suffix, fig_format=fig_format, fig_dpi=fig_dpi,
                genes_of_interest=genes_of_interest or None,
            )
            ma_plot(
                deseq2_all, cond_fig_dir, cols=deseq2_cols,
                padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff,
                basemean_cutoff=basemean_cutoff, label=full_label,
                suffix=suffix, fig_format=fig_format, fig_dpi=fig_dpi,
            )
            expression_rank_plot(
                deseq2_all, cond_fig_dir, cols=deseq2_cols,
                padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff,
                basemean_cutoff=basemean_cutoff, label=full_label,
                suffix=suffix, fig_format=fig_format, fig_dpi=fig_dpi,
            )
            if interactive_plots:
                volcano_plot_interactive(
                    deseq2_all, cond_fig_dir, cols=deseq2_cols,
                    padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff,
                    basemean_cutoff=basemean_cutoff, label=full_label,
                    suffix=suffix,
                )
                ma_plot_interactive(
                    deseq2_all, cond_fig_dir, cols=deseq2_cols,
                    padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff,
                    basemean_cutoff=basemean_cutoff, label=full_label,
                    suffix=suffix,
                )
            if len(deseq2_filt) > 0:
                biotype_chart(
                    deseq2_filt, cond_fig_dir, cols=deseq2_cols,
                    label=full_label, suffix=suffix,
                    fig_format=fig_format, fig_dpi=fig_dpi,
                )
                if bio_filter is None:
                    top_genes_lollipop(
                        deseq2_filt, cond_fig_dir, cols=deseq2_cols,
                        log2fc_cutoff=log2fc_cutoff, label=full_label,
                        suffix=suffix, fig_format=fig_format,
                    )
                    biotype_direction_chart(
                        deseq2_filt, cond_fig_dir, cols=deseq2_cols,
                        label=full_label, suffix=suffix,
                        fig_format=fig_format, fig_dpi=fig_dpi,
                    )
                    biotype_enrichment_test(
                        deseq2_filt, deseq2_all, cond_fig_dir, cols=deseq2_cols,
                        label=full_label, suffix=suffix,
                        fig_format=fig_format, fig_dpi=fig_dpi,
                    )
            else:
                print(f"  No genes passed filter -- skipping biotype chart for {full_label}")
            if bio_filter is None:
                biotype_volcano(
                    deseq2_all, cond_fig_dir, cols=deseq2_cols,
                    padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff,
                    basemean_cutoff=basemean_cutoff, label=full_label,
                    suffix=suffix, fig_format=fig_format, fig_dpi=fig_dpi,
                )
                ecdf_log2fc_by_biotype(
                    deseq2_all, cond_fig_dir, cols=deseq2_cols,
                    log2fc_cutoff=log2fc_cutoff, basemean_cutoff=basemean_cutoff,
                    label=full_label, suffix=suffix,
                    fig_format=fig_format, fig_dpi=fig_dpi,
                )

        # --- rMATS (if directory exists) ---
        rmats_raw: dict = {}
        rmats_filtered: dict = {}

        if cond.get("rmats_dir"):
            print(f"\n-- Loading rMATS for {cond_label} --")
            rmats_raw = load_all_rmats(
                cond["rmats_dir"], cols=rmats_cols, file_suffix=rmats_file_suffix
            )

            print(f"\n-- rMATS Filtering --")
            filtered_counts: dict[str, int] = {}
            for event_type, df in rmats_raw.items():
                raw, filt = filter_rmats(
                    df, cols=rmats_cols,
                    fdr_cutoff=fdr_cutoff,
                    pval_cutoff=pval_cutoff,
                    dpsi_cutoff=dpsi_cutoff,
                    use_fdr=use_fdr,
                    dual_filter=dual_filter,
                    event_type=event_type,
                )
                rmats_raw[event_type] = raw
                rmats_filtered[event_type] = filt
                filtered_counts[event_type] = len(filt)

            total_sig = sum(filtered_counts.values())
            print(f"  TOTAL significant events: {total_sig:,}")

            # RBP annotation on filtered rMATS data
            if rbp_annotations:
                print(f"\n-- RBP Annotation on rMATS [{cond_label}] --")
                for et in rmats_filtered:
                    if not rmats_filtered[et].empty:
                        rmats_filtered[et] = annotate_rbps(
                            rmats_filtered[et], rbp_annotations,
                            gene_col=rmats_cols["gene_name"],
                        )

            # -- Per-condition rMATS visualization --
            print(f"\n-- Generating rMATS Figures [{cond_label}] --")
            for event_type, df in rmats_raw.items():
                rmats_scatter(
                    df, event_type, cond_fig_dir, rmats_cols=rmats_cols,
                    fdr_cutoff=fdr_cutoff, pval_cutoff=pval_cutoff,
                    dpsi_cutoff=dpsi_cutoff, use_fdr=use_fdr,
                    dual_filter=dual_filter,
                    fig_format=fig_format, fig_dpi=fig_dpi,
                )
                rmats_psi_scatter(
                    rmats_raw, rmats_filtered, event_type, cond_fig_dir,
                    rmats_cols=rmats_cols, dpsi_cutoff=dpsi_cutoff,
                    fig_format=fig_format, fig_dpi=fig_dpi,
                )
            rmats_combined_volcano(
                rmats_raw, cond_fig_dir, rmats_cols=rmats_cols,
                fdr_cutoff=fdr_cutoff, pval_cutoff=pval_cutoff,
                dpsi_cutoff=dpsi_cutoff, use_fdr=use_fdr,
                fig_format=fig_format, fig_dpi=fig_dpi,
            )
            rmats_event_summary_chart(
                filtered_counts, cond_fig_dir,
                use_fdr=use_fdr, fdr_cutoff=fdr_cutoff,
                pval_cutoff=pval_cutoff, dpsi_cutoff=dpsi_cutoff,
                fig_format=fig_format, fig_dpi=fig_dpi,
            )
            rmats_dpsi_distribution(
                rmats_filtered, cond_fig_dir, rmats_cols=rmats_cols,
                dpsi_cutoff=dpsi_cutoff,
                fig_format=fig_format, fig_dpi=fig_dpi,
            )
        else:
            print(f"\n  No rMATS data for {cond_label}, skipping rMATS analysis")

        # Per-condition export
        export_excel(
            deseq2_filtered_sets, rmats_filtered, cond_outdir,
            log2fc_cutoff=log2fc_cutoff, basemean_cutoff=basemean_cutoff,
            padj_cutoff=padj_cutoff, use_fdr=use_fdr,
            fdr_cutoff=fdr_cutoff, pval_cutoff=pval_cutoff,
            dpsi_cutoff=dpsi_cutoff,
        )

        # Store for cross-condition comparisons
        condition_results[cond_name] = {
            "deseq2_raw": deseq2_raw,
            "deseq2_filtered": deseq2_filtered_sets,
            "rmats_raw": rmats_raw,
            "rmats_filtered": rmats_filtered,
        }

    # ===================================================================
    # PHASE 2: Cross-condition comparisons
    # ===================================================================
    print(f"\n{'=' * 60}")
    print("  Cross-Condition Comparisons")
    print(f"{'=' * 60}")

    comparison_dir = outdir / "cross_condition"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    comparison_fig_dir = comparison_dir / "figures"
    comparison_fig_dir.mkdir(exist_ok=True)

    cross_data: dict = {}

    # DESeq2 comparisons (all conditions)
    print("\n-- DESeq2 DE Counts Overview --")
    deseq2_de_counts_chart(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
    )

    print("\n-- DESeq2 Venn Diagrams (3-way) --")
    gene_sets = extract_gene_sets(condition_results, deseq2_cols)
    deseq2_venn_diagrams(
        gene_sets, condition_labels, comparison_fig_dir,
        fig_format=fig_format, fig_dpi=fig_dpi,
    )

    print("\n-- Pairwise DEG Venn Diagrams --")
    pairwise_deg_venns(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
    )

    print("\n-- DESeq2 UpSet Plots --")
    deseq2_upset_plot(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
    )

    print("\n-- DESeq2 Direction Concordance --")
    concordance_df = compute_direction_heatmap(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
    )
    cross_data["concordance_matrix"] = concordance_df

    print("\n-- DESeq2 Log2FC Heatmap --")
    log2fc_df = deseq2_log2fc_heatmap(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
    )
    cross_data["log2fc_matrix"] = log2fc_df

    # Top DEG heatmap (requires counts matrix)
    if counts_df is not None:
        print("\n-- Top DEG Expression Heatmap --")
        compute_top_deg_heatmap(
            counts_df, condition_results, condition_labels,
            sample_metadata, comparison_fig_dir,
            cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
        )

    print("\n-- Pairwise log2FC Scatter --")
    pairwise_log2fc_scatter(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
    )

    print("\n-- Cross-Condition Biotype Comparison --")
    cross_condition_biotype_comparison(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
    )
    cross_condition_biotype_direction(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
    )

    # RBP cross-condition analysis (only if RBP annotations were loaded)
    if rbp_annotations:
        print("\n-- RBP Cross-Condition Heatmap --")
        try:
            rbp_heatmap(
                condition_results, condition_labels, comparison_fig_dir,
                cols=deseq2_cols, fig_format=fig_format, fig_dpi=fig_dpi,
            )
        except Exception as e:
            print(f"  WARNING: RBP heatmap failed: {e}")

        print("\n-- RBP Summary Table --")
        try:
            rbp_summary_table(
                condition_results, condition_labels, comparison_dir,
                cols=deseq2_cols,
            )
        except Exception as e:
            print(f"  WARNING: RBP summary table failed: {e}")

    # rMATS comparisons (only conditions with rMATS data)
    rmats_conditions = {
        name: res for name, res in condition_results.items()
        if res["rmats_filtered"]
    }
    if len(rmats_conditions) >= 2:
        print("\n-- rMATS Cross-Condition Venn Diagrams --")
        rmats_cross_condition_venn(
            rmats_conditions, condition_labels, comparison_fig_dir,
            rmats_cols=rmats_cols, fig_format=fig_format, fig_dpi=fig_dpi,
        )
        try:
            rmats_cross_condition_venn(
                rmats_conditions, condition_labels, comparison_fig_dir,
                rmats_cols=rmats_cols, fig_format=fig_format, fig_dpi=fig_dpi,
                match_by="gene",
            )
        except Exception as e:
            print(f"  [WARN] Gene-level cross-condition Venns failed: {e}")

        print("\n-- rMATS Event Count Comparison --")
        rmats_event_count_comparison(
            rmats_conditions, condition_labels, comparison_fig_dir,
            fig_format=fig_format, fig_dpi=fig_dpi,
        )

        print("\n-- rMATS UpSet Plots --")
        rmats_upset_plot(
            rmats_conditions, condition_labels, comparison_fig_dir,
            rmats_cols=rmats_cols, fig_format=fig_format, fig_dpi=fig_dpi,
        )
        try:
            rmats_upset_plot(
                rmats_conditions, condition_labels, comparison_fig_dir,
                rmats_cols=rmats_cols, fig_format=fig_format, fig_dpi=fig_dpi,
                match_by="gene",
            )
        except Exception as e:
            print(f"  [WARN] Gene-level UpSet plots failed: {e}")

        print("\n-- rMATS Direction Concordance --")
        rmats_conc = rmats_direction_concordance(
            rmats_conditions, condition_labels, comparison_fig_dir,
            rmats_cols=rmats_cols, fig_format=fig_format, fig_dpi=fig_dpi,
        )
        cross_data["rmats_concordance"] = rmats_conc

        print("\n-- Pairwise Splicing Venn Diagrams --")
        pairwise_splicing_venns(
            rmats_conditions, condition_labels, comparison_fig_dir,
            rmats_cols=rmats_cols, dpsi_cutoff=dpsi_cutoff,
            fig_format=fig_format, fig_dpi=fig_dpi,
        )
        try:
            pairwise_splicing_venns(
                rmats_conditions, condition_labels, comparison_fig_dir,
                rmats_cols=rmats_cols, dpsi_cutoff=dpsi_cutoff,
                fig_format=fig_format, fig_dpi=fig_dpi,
                match_by="gene",
            )
        except Exception as e:
            print(f"  [WARN] Gene-level splicing Venns failed: {e}")

        print("\n-- Pairwise dPSI Scatter --")
        pairwise_dpsi_scatter(
            rmats_conditions, condition_labels, comparison_fig_dir,
            rmats_cols=rmats_cols, fig_format=fig_format, fig_dpi=fig_dpi,
        )
        try:
            pairwise_dpsi_scatter(
                rmats_conditions, condition_labels, comparison_fig_dir,
                rmats_cols=rmats_cols, fig_format=fig_format, fig_dpi=fig_dpi,
                match_by="gene",
            )
        except Exception as e:
            print(f"  [WARN] Gene-level dPSI scatter failed: {e}")

        # Directional Venn Diagrams -- REMOVED: redundant with pairwise Venns
        # (rmats_directional_venn_diagrams was intentionally not extracted)

        print("\n-- rMATS Event Heatmaps --")
        for _et in ["SE", "RI"]:
            try:
                rmats_event_heatmap(
                    rmats_conditions, condition_labels, _et, comparison_fig_dir,
                    rmats_cols=rmats_cols, dpsi_cutoff=dpsi_cutoff,
                    fig_format=fig_format, fig_dpi=fig_dpi,
                )
            except Exception as e:
                print(f"  WARNING: {_et} event heatmap failed: {e}")

        print("\n-- rMATS Event Type Pie Chart --")
        try:
            rmats_event_pie_chart(
                rmats_conditions, condition_labels, comparison_fig_dir,
                fig_format=fig_format, fig_dpi=fig_dpi,
            )
        except Exception as e:
            print(f"  WARNING: Event type pie chart failed: {e}")

        print("\n-- Pairwise Comparison Workbooks --")
        try:
            export_pairwise_workbook(
                rmats_conditions, condition_labels, comparison_fig_dir,
                cols=deseq2_cols, rmats_cols=rmats_cols,
                dpsi_cutoff=dpsi_cutoff,
            )
        except Exception as e:
            print(f"  WARNING: Pairwise workbook export failed: {e}")
    else:
        print("\n  Fewer than 2 conditions have rMATS data, skipping rMATS comparisons")

    # Combined DESeq2 + rMATS
    print("\n-- Combined DESeq2 + rMATS Analyses --")
    deseq2_vs_rmats_venn(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, rmats_cols=rmats_cols,
        fig_format=fig_format, fig_dpi=fig_dpi,
    )
    log2fc_vs_dpsi_scatter(
        condition_results, condition_labels, comparison_fig_dir,
        cols=deseq2_cols, rmats_cols=rmats_cols,
        dpsi_cutoff=dpsi_cutoff, log2fc_cutoff=log2fc_cutoff,
        fig_format=fig_format, fig_dpi=fig_dpi,
    )

    # Master combined export
    print("\n-- Exporting Combined Multi-Condition Results --")
    export_combined_results(
        condition_results, cross_data, comparison_dir,
        cols=deseq2_cols, log2fc_cutoff=log2fc_cutoff,
        basemean_cutoff=basemean_cutoff, padj_cutoff=padj_cutoff,
    )

    # Unfiltered merged overlap Excel (all genes, all conditions)
    export_unfiltered_merged(
        condition_results, condition_labels, outdir,
        deseq2_cols=deseq2_cols, rmats_cols=rmats_cols,
    )

    # ===================================================================
    # PHASE 3: GSEA, GO ORA, Prism, PowerPoint, and Validation
    # ===================================================================
    print(f"\n{'=' * 60}")
    print("  GSEA, GO ORA, Prism Export, PowerPoint, and Validation")
    print(f"{'=' * 60}")

    # GSEA enrichment
    gsea_results = run_gsea_enrichment(
        condition_results, condition_labels, outdir,
        cols=deseq2_cols, databases=gsea_databases,
        ranking_method=gsea_ranking,
        min_size=gsea_min_size, max_size=gsea_max_size,
        permutations=gsea_permutations,
    )

    # GO Over-Representation Analysis
    if ora_method == "both":
        # Run both Enrichr and g:Profiler side-by-side for comparison
        print("\n-- Running Enrichr ORA --")
        try:
            go_results_enrichr = run_enrichr_ora(
                condition_results, condition_labels, outdir,
                cols=deseq2_cols,
            )
        except Exception as e:
            print(f"  WARNING: Enrichr ORA failed: {e}")
            go_results_enrichr = {}
        if go_results_enrichr:
            go_enrichment_combined_plot(
                go_results_enrichr, condition_labels,
                comparison_fig_dir,
                fig_format=fig_format, fig_dpi=fig_dpi,
                filename_suffix="_enrichr",
            )
            export_go_prism(
                go_results_enrichr, condition_labels,
                outdir / "prism_files", filename_suffix="_enrichr",
            )

        print("\n-- Running g:Profiler ORA --")
        try:
            go_results_gprofiler = run_gprofiler_ora(
                condition_results, condition_labels, outdir,
                cols=deseq2_cols, species=species,
            )
        except Exception as e:
            print(f"  WARNING: g:Profiler ORA failed: {e}")
            go_results_gprofiler = {}
        if go_results_gprofiler:
            go_enrichment_combined_plot(
                go_results_gprofiler, condition_labels,
                comparison_fig_dir,
                fig_format=fig_format, fig_dpi=fig_dpi,
                filename_suffix="_gprofiler",
            )
            export_go_prism(
                go_results_gprofiler, condition_labels,
                outdir / "prism_files", filename_suffix="_gprofiler",
            )

        # Use Enrichr as primary for summary dashboard (HSCHARME paper standard)
        go_results = go_results_enrichr if go_results_enrichr else go_results_gprofiler
    elif ora_method == "gprofiler":
        go_results = run_gprofiler_ora(
            condition_results, condition_labels, outdir,
            cols=deseq2_cols, species=species,
        )
        go_enrichment_combined_plot(
            go_results, condition_labels, comparison_fig_dir,
            fig_format=fig_format, fig_dpi=fig_dpi,
        )
        export_go_prism(go_results, condition_labels, outdir / "prism_files")
    else:
        go_results = run_enrichr_ora(
            condition_results, condition_labels, outdir,
            cols=deseq2_cols,
        )
        go_enrichment_combined_plot(
            go_results, condition_labels, comparison_fig_dir,
            fig_format=fig_format, fig_dpi=fig_dpi,
        )
        export_go_prism(go_results, condition_labels, outdir / "prism_files")

    # Combined GSEA dot plots (replaces legacy gsea_dotplot)
    gsea_combined_plot(
        gsea_results, condition_labels, outdir,
        fig_format=fig_format, fig_dpi=fig_dpi,
    )

    # GSEA enrichment plots and leading edge export
    gsea_enrichment_plots(
        gsea_results, condition_labels, outdir,
        fig_format=fig_format, fig_dpi=fig_dpi,
    )
    export_gsea_leading_edge(gsea_results, condition_labels, outdir)

    # Gene overlap summary
    gene_overlap_summary(
        condition_results, condition_labels,
        comparison_dir, cols=deseq2_cols, rmats_cols=rmats_cols,
    )

    # Summary dashboard
    summary_dashboard(
        condition_results, condition_labels, go_results, gsea_results,
        comparison_fig_dir, cols=deseq2_cols, rmats_cols=rmats_cols,
        fig_format=fig_format, fig_dpi=fig_dpi,
    )

    # Prism export
    export_prism_pzfx(
        condition_results, condition_labels, outdir,
        cols=deseq2_cols, gsea_results=gsea_results,
    )

    # PowerPoint generation
    export_powerpoint(
        condition_results, condition_labels, outdir,
        padj_cutoff=padj_cutoff, log2fc_cutoff=log2fc_cutoff,
        basemean_cutoff=basemean_cutoff, fdr_cutoff=fdr_cutoff,
        pval_cutoff=pval_cutoff, dpsi_cutoff=dpsi_cutoff,
        use_fdr=use_fdr, dual_filter=dual_filter,
        fig_format=fig_format,
    )

    # Validation
    validation_passed = validate_outputs(
        condition_results, condition_labels, outdir,
        cols=deseq2_cols, rmats_cols=rmats_cols,
    )

    print(f"\n{'=' * 60}")
    print(f"  Done! All outputs in: {outdir.resolve()}")
    if validation_passed:
        print("  [PASS] All validation checks passed")
    else:
        print("  [WARN] Some validation checks failed (see report above)")
    print(f"{'=' * 60}")
