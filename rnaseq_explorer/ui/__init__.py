"""UI subpackage: Streamlit pages and components for RNA-seq Explorer.

Modules:
    app      - Main Streamlit entry point (multi-page app)
    sidebar  - Shared sidebar configuration (uploads, thresholds, settings)
    styles   - CSS styling for metric cards, tabs, dark/light mode

Pages (rnaseq_explorer.ui.pages):
    overview              - Summary metrics and quick-look charts
    deseq2_page           - DEG visualizations and data table
    splicing_page         - rMATS splicing visualizations
    enrichment_page       - GSEA and ORA enrichment visualizations
    qc_page               - PCA, correlation, and DEG heatmaps
    cross_condition_page  - Multi-condition comparison
    genewalk_page         - GeneWalk functional annotations
    gene_investigator_page - Per-gene evidence aggregation
"""
