"""Shared pytest fixtures for RNA-seq Explorer test suite."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def sample_deseq2_df() -> pd.DataFrame:
    """DataFrame mimicking DESeq2 results: ~30 rows with a mix of sig up/down/NS."""
    np.random.seed(42)
    n = 30
    genes = [f"GENE{i}" for i in range(1, n + 1)]
    genes[0] = "MIAT"
    genes[1] = "QKI"
    genes[2] = "TP53"
    genes[3] = "BRCA1"

    log2fc = np.concatenate([
        np.random.uniform(1.5, 4.0, 8),     # significant up
        np.random.uniform(-4.0, -1.5, 7),   # significant down
        np.random.uniform(-0.3, 0.3, 15),   # not significant
    ])
    padj = np.concatenate([
        np.random.uniform(1e-10, 0.01, 8),
        np.random.uniform(1e-10, 0.01, 7),
        np.random.uniform(0.1, 1.0, 15),
    ])
    basemean = np.random.uniform(50, 5000, n)
    pvalue = padj * np.random.uniform(0.5, 1.0, n)
    lfcSE = np.abs(log2fc) * np.random.uniform(0.1, 0.3, n)
    stat = log2fc / lfcSE

    biotypes = (
        ["Protein Coding"] * 10
        + ["lncRNA"] * 5
        + ["Pseudogene"] * 5
        + ["Small ncRNA"] * 5
        + ["Other"] * 5
    )
    directions = []
    for lfc, p in zip(log2fc, padj):
        if p < 0.05 and lfc > 1.0:
            directions.append("Up")
        elif p < 0.05 and lfc < -1.0:
            directions.append("Down")
        else:
            directions.append("NS")

    return pd.DataFrame({
        "gene_name": genes,
        "baseMean": basemean,
        "log2FoldChange": log2fc,
        "lfcSE": lfcSE,
        "stat": stat,
        "pvalue": pvalue,
        "padj": padj,
        "biotype_group": biotypes,
        "direction": directions,
    })


@pytest.fixture
def sample_rmats_df() -> pd.DataFrame:
    """DataFrame mimicking rMATS output: ~20 rows with multiple event types."""
    np.random.seed(43)
    n = 20
    event_types = ["SE"] * 6 + ["MXE"] * 4 + ["A3SS"] * 4 + ["A5SS"] * 3 + ["RI"] * 3
    gene_ids = [f"GENE{i}" for i in range(1, n + 1)]
    gene_ids[0] = "MIAT"
    gene_ids[1] = "QKI"

    dpsi = np.random.uniform(-0.5, 0.5, n)
    fdr = np.concatenate([
        np.random.uniform(1e-6, 0.01, 10),
        np.random.uniform(0.05, 1.0, 10),
    ])
    pvalue = fdr * np.random.uniform(0.5, 1.0, n)
    inc1 = np.random.uniform(0.2, 0.8, n)
    inc2 = inc1 + dpsi

    return pd.DataFrame({
        "GeneID": gene_ids,
        "event_type": event_types,
        "IncLevelDifference": dpsi,
        "FDR": fdr,
        "PValue": pvalue,
        "IncLevel1": [f"{v:.3f}" for v in inc1],
        "IncLevel2": [f"{v:.3f}" for v in inc2],
    })


@pytest.fixture
def sample_genewalk_df() -> pd.DataFrame:
    """DataFrame mimicking genewalk_results.csv: ~15 rows."""
    np.random.seed(44)
    n = 15
    genes = ["MIAT"] * 5 + ["QKI"] * 5 + ["TP53"] * 5
    go_names = [
        "RNA binding", "nucleus", "transcription", "mRNA splicing", "chromatin",
        "RNA processing", "cytoplasm", "mRNA transport", "translation", "ribosome",
        "apoptosis", "cell cycle", "DNA repair", "tumor suppression", "p53 pathway",
    ]
    go_ids = [f"GO:{i:07d}" for i in range(1, n + 1)]
    domains = (
        ["molecular_function"] * 3 + ["cellular_component"] * 2
        + ["biological_process"] * 3 + ["cellular_component"] * 2
        + ["biological_process"] * 5
    )
    sim = np.random.uniform(0.1, 0.9, n)
    gene_padj = np.concatenate([
        np.random.uniform(0.001, 0.05, 8),
        np.random.uniform(0.1, 0.5, 7),
    ])
    global_padj = gene_padj * np.random.uniform(0.8, 1.2, n)

    return pd.DataFrame({
        "hgnc_symbol": genes,
        "go_name": go_names,
        "go_id": go_ids,
        "go_domain": domains,
        "sim": sim,
        "gene_padj": gene_padj,
        "global_padj": global_padj,
    })


@pytest.fixture
def sample_gsea_results() -> pd.DataFrame:
    """DataFrame mimicking GSEA output: ~10 rows."""
    np.random.seed(45)
    n = 10
    terms = [
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_E2F_TARGETS",
        "KEGG_SPLICEOSOME",
        "REACTOME_MRNA_PROCESSING",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "KEGG_RIBOSOME",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "REACTOME_CELL_CYCLE",
        "HALLMARK_APOPTOSIS",
    ]
    nes = np.array([2.1, 1.8, 1.5, -1.9, -1.6, 1.3, -2.0, 1.1, -1.2, -0.9])
    fdr = np.array([0.001, 0.005, 0.01, 0.002, 0.008, 0.03, 0.001, 0.05, 0.1, 0.25])
    gene_pct = np.random.uniform(5, 40, n)
    lead_genes = [
        "MIAT;QKI;TP53", "MYC;MAX;CDK4", "E2F1;RB1;CDK2",
        "SRSF1;U2AF2;PRPF8", "HNRNPA1;SRSF3;QKI",
        "NFKB1;RELA;TNF", "RPL5;RPS6;RPL11",
        "STAT1;IRF1;GBP1", "CDK1;CCNB1;CDC20",
        "BAX;BCL2;CASP3",
    ]

    return pd.DataFrame({
        "Term": terms,
        "NES": nes,
        "FDR q-val": fdr,
        "Gene %": gene_pct,
        "Lead_genes": lead_genes,
    })


@pytest.fixture
def sample_ora_results() -> pd.DataFrame:
    """DataFrame mimicking Enrichr/g:Profiler ORA output: ~10 rows."""
    np.random.seed(46)
    n = 10
    terms = [
        "GO:0003723 RNA binding",
        "GO:0005634 nucleus",
        "KEGG:03040 Spliceosome",
        "GO:0006397 mRNA processing",
        "GO:0008380 RNA splicing",
        "REAC:R-HSA-72163 mRNA Splicing",
        "GO:0005829 cytosol",
        "GO:0006915 apoptotic process",
        "KEGG:03010 Ribosome",
        "GO:0006412 translation",
    ]
    overlap = [
        "15/200", "12/150", "8/120", "10/180", "9/160",
        "7/140", "11/190", "6/130", "8/110", "5/100",
    ]
    pvalue = np.array([1e-8, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 0.01, 0.02, 0.05])
    adj_pvalue = pvalue * 5
    adj_pvalue = np.clip(adj_pvalue, 0, 1)
    combined_score = np.random.uniform(10, 500, n)

    return pd.DataFrame({
        "Term": terms,
        "Overlap": overlap,
        "P-value": pvalue,
        "Adjusted P-value": adj_pvalue,
        "Combined Score": combined_score,
    })


@pytest.fixture
def empty_df() -> pd.DataFrame:
    """Empty DataFrame."""
    return pd.DataFrame()
