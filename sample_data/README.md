# Sample Data for RNA-seq Explorer

These files let you try the app immediately without your own data.

## Files

### `sample_deseq2_results.csv`
Simulated DESeq2 differential expression output (50 genes). Contains a mix of significantly up-regulated genes (e.g., MYC, VEGFA, EGFR), significantly down-regulated genes (e.g., TP53, PTEN, MTOR), and non-significant housekeeping genes (e.g., GAPDH, ACTB).

**Columns:** gene_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

### `sample_rmats_SE.txt`
Simulated rMATS Skipped Exon (SE) output (20 events). Contains splicing events in RNA-binding proteins and splicing factors with a range of significance levels and delta-PSI values.

**Columns:** Standard rMATS SE format (ID, GeneID, geneSymbol, chr, strand, coordinates, junction counts, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference)

### `sample_genewalk_results.csv`
Simulated GeneWalk functional annotation output (30 gene-GO term associations). Includes genes with significant GO associations across biological process, molecular function, and cellular component domains.

**Columns:** hgnc_symbol, go_name, go_id, go_domain, sim, gene_padj, global_padj

## Usage

1. Launch the app: `streamlit run rnaseq_explorer/app.py`
2. Use the sidebar file uploaders to load sample files
3. For DESeq2: upload `sample_deseq2_results.csv`
4. For Splicing: upload `sample_rmats_SE.txt`
5. For GeneWalk: upload `sample_genewalk_results.csv`
