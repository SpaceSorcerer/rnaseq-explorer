"""Engine subpackage: core analysis modules for RNA-seq Explorer.

Modules:
    deseq2          - DESeq2 loading, filtering, annotation
    rmats           - rMATS event parsing and filtering
    gsea            - GSEA prerank enrichment runner
    ora             - Over-representation analysis (Enrichr + g:Profiler)
    qc              - QC analyses (PCA, correlation, heatmaps)
    cross_condition - Multi-condition comparisons (Venn, UpSet, concordance)
    exports         - Export to Excel, Prism .pzfx, PowerPoint
    pipeline        - Backward-compatible orchestration wrapper
"""
