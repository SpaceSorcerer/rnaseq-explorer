"""
RNA-seq Analysis Modules
========================

Modular components for RNA-seq differential expression and splicing analysis.

Modules:
--------
- deseq2: DESeq2 data loading, filtering, and summarization functions
- standardize: Gene ID/name lookup, biotype normalization, and data preprocessing
- splicing: rMATS alternative splicing analysis, filtering, and event parsing
"""

from . import deseq2
from . import standardize
from . import splicing

__all__ = ["deseq2", "standardize", "splicing"]
