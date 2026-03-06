#!/usr/bin/env python3
"""
RNA-seq Analysis Pipeline - Main Entry Point
============================================

Single command-line interface for running the complete RNA-seq analysis pipeline:
- DESeq2 differential expression analysis
- rMATS alternative splicing analysis
- Cross-condition concordance analysis
- Publication-quality figure generation
- Excel export of filtered results

This script replaces all previous run_*.py scripts with a unified, configurable CLI.

Usage:
    python run_pipeline.py
    python run_pipeline.py --config my_config.yaml
    python run_pipeline.py --conditions MIAT_OE_vs_Control,QKI_KO_vs_WT
    python run_pipeline.py --modules deseq2,splicing --skip-figures
    python run_pipeline.py --output-dir custom_output/

Author: Data - RNA-seq Pipeline Engineer
Date: 2026-03-06
"""

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings

import pandas as pd
import numpy as np
import yaml

# Import pipeline modules
from modules import deseq2, splicing, figures, concordance, standardize


# ─── Configuration & Setup ───────────────────────────────────────────────────

def setup_logging(output_dir: Path, log_level: str = "INFO", log_file: str = "pipeline.log"):
    """
    Configure logging to both console and file.

    Args:
        output_dir: Directory for log file
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Name of log file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / log_file

    # Create formatters
    file_formatter = logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_formatter = logging.Formatter('%(levelname)-8s | %(message)s')

    # File handler
    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, log_level.upper()))
    console_handler.setFormatter(console_formatter)

    # Configure root logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers = []  # Clear existing handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def load_config(config_path: Path) -> dict:
    """
    Load YAML configuration file.

    Args:
        config_path: Path to rnaseq_config.yaml

    Returns:
        Configuration dictionary
    """
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    logging.info(f"Loaded configuration from {config_path}")
    return config


def validate_config(config: dict) -> None:
    """
    Validate that all required configuration fields are present.

    Args:
        config: Configuration dictionary

    Raises:
        ValueError: If required fields are missing
    """
    required_fields = ['base_dir', 'output_dir', 'conditions', 'filtering']
    for field in required_fields:
        if field not in config:
            raise ValueError(f"Missing required field in config: {field}")

    if not config['conditions']:
        raise ValueError("No conditions defined in config")

    logging.info("Configuration validated successfully")


# ─── Data Loading ────────────────────────────────────────────────────────────

def load_condition_data(
    condition: dict,
    config: dict,
    load_deseq: bool = True,
    load_splicing: bool = True
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Load DESeq2 and rMATS data for a single condition.

    Args:
        condition: Condition dictionary from config
        config: Full configuration dictionary
        load_deseq: Whether to load DESeq2 data
        load_splicing: Whether to load splicing data

    Returns:
        Tuple of (deseq2_df, rmats_df) - either can be None if not loaded
    """
    condition_name = condition['name']
    base_dir = Path(config['base_dir'])

    deseq_df = None
    rmats_df = None

    # Load DESeq2 data
    if load_deseq:
        try:
            logging.info(f"Loading DESeq2 data for {condition_name}...")
            deseq_df = deseq2.load_and_normalize_deseq2(
                condition=condition,
                base_dir=base_dir
            )
            logging.info(f"  → Loaded {len(deseq_df)} genes")
        except Exception as e:
            logging.error(f"Failed to load DESeq2 data for {condition_name}: {e}")
            raise

    # Load rMATS splicing data
    if load_splicing:
        try:
            logging.info(f"Loading rMATS data for {condition_name}...")
            rmats_dir = base_dir / condition['rmats_dir']
            rmats_df = splicing.load_rmats(
                rmats_dir=rmats_dir,
                event_types=['SE', 'A3SS', 'A5SS', 'MXE', 'RI']
            )
            logging.info(f"  → Loaded {len(rmats_df)} splicing events")
        except Exception as e:
            logging.warning(f"Failed to load rMATS data for {condition_name}: {e}")
            # Don't raise - splicing data is optional
            rmats_df = None

    return deseq_df, rmats_df


def load_all_conditions(
    config: dict,
    condition_subset: Optional[List[str]] = None,
    load_deseq: bool = True,
    load_splicing: bool = True
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """
    Load data for all conditions specified in config.

    Args:
        config: Configuration dictionary
        condition_subset: Optional list of condition names to load (None = all)
        load_deseq: Whether to load DESeq2 data
        load_splicing: Whether to load splicing data

    Returns:
        Tuple of (deseq_dict, rmats_dict) mapping condition names to dataframes
    """
    deseq_dict = {}
    rmats_dict = {}

    conditions = config['conditions']
    if condition_subset:
        conditions = [c for c in conditions if c['name'] in condition_subset]
        if not conditions:
            raise ValueError(f"No conditions matched subset: {condition_subset}")

    logging.info(f"Loading data for {len(conditions)} condition(s)...")

    for condition in conditions:
        condition_name = condition['name']
        deseq_df, rmats_df = load_condition_data(
            condition, config, load_deseq, load_splicing
        )

        if deseq_df is not None:
            deseq_dict[condition_name] = deseq_df
        if rmats_df is not None:
            rmats_dict[condition_name] = rmats_df

    return deseq_dict, rmats_dict


# ─── Filtering ───────────────────────────────────────────────────────────────

def filter_all_conditions(
    deseq_dict: Dict[str, pd.DataFrame],
    rmats_dict: Dict[str, pd.DataFrame],
    config: dict
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """
    Apply filtering thresholds to all conditions.

    Args:
        deseq_dict: DESeq2 dataframes by condition name
        rmats_dict: rMATS dataframes by condition name
        config: Configuration dictionary

    Returns:
        Tuple of (filtered_deseq_dict, filtered_rmats_dict)
    """
    deseq_cutoffs = config['filtering']['deseq2']
    splicing_cutoffs = config['filtering']['splicing']

    filtered_deseq = {}
    filtered_rmats = {}

    # Filter DESeq2 data
    logging.info("Applying DESeq2 filters...")
    for condition_name, df in deseq_dict.items():
        df_sig = deseq2.filter_deseq2(
            df=df,
            log2fc_cutoff=deseq_cutoffs['log2fc_cutoff'],
            basemean_cutoff=deseq_cutoffs['basemean_cutoff'],
            padj_cutoff=deseq_cutoffs['padj_cutoff'],
            label=condition_name
        )
        filtered_deseq[condition_name] = df_sig
        logging.info(f"  {condition_name}: {len(df_sig)} significant DEGs")

    # Filter rMATS data
    if rmats_dict:
        logging.info("Applying splicing filters...")
        for condition_name, df in rmats_dict.items():
            df_sig = splicing.filter_rmats_events(
                df=df,
                fdr_cutoff=splicing_cutoffs['fdr_cutoff'],
                pval_cutoff=splicing_cutoffs['pval_cutoff'],
                inclevel_diff_cutoff=splicing_cutoffs['dpsi_cutoff'],
                require_both=splicing_cutoffs.get('use_pval_and_fdr', True)
            )
            filtered_rmats[condition_name] = df_sig
            logging.info(f"  {condition_name}: {len(df_sig)} significant splicing events")

    return filtered_deseq, filtered_rmats


# ─── Summary Statistics ──────────────────────────────────────────────────────

def generate_summary_stats(
    deseq_all: Dict[str, pd.DataFrame],
    deseq_filtered: Dict[str, pd.DataFrame],
    rmats_all: Dict[str, pd.DataFrame],
    rmats_filtered: Dict[str, pd.DataFrame],
    config: dict
) -> pd.DataFrame:
    """
    Generate summary statistics table for all conditions.

    Args:
        deseq_all: All DESeq2 data by condition
        deseq_filtered: Filtered DESeq2 data by condition
        rmats_all: All rMATS data by condition
        rmats_filtered: Filtered rMATS data by condition
        config: Configuration dictionary

    Returns:
        DataFrame with summary statistics
    """
    rows = []

    for condition in config['conditions']:
        condition_name = condition['name']
        condition_label = condition['label']

        row = {
            'Condition': condition_label,
            'Total_Genes': len(deseq_all.get(condition_name, [])),
            'Significant_DEGs': len(deseq_filtered.get(condition_name, [])),
        }

        # Count up/down-regulated
        if condition_name in deseq_filtered:
            df_sig = deseq_filtered[condition_name]
            row['DEGs_Upregulated'] = (df_sig['log2FoldChange'] > 0).sum()
            row['DEGs_Downregulated'] = (df_sig['log2FoldChange'] < 0).sum()

        # Splicing stats
        if condition_name in rmats_all:
            row['Total_Splicing_Events'] = len(rmats_all[condition_name])
            row['Significant_Splicing'] = len(rmats_filtered.get(condition_name, []))

        rows.append(row)

    return pd.DataFrame(rows)


# ─── Figure Generation ───────────────────────────────────────────────────────

def generate_all_figures(
    deseq_all: Dict[str, pd.DataFrame],
    deseq_filtered: Dict[str, pd.DataFrame],
    rmats_all: Dict[str, pd.DataFrame],
    rmats_filtered: Dict[str, pd.DataFrame],
    config: dict,
    output_dir: Path
) -> None:
    """
    Generate all publication-quality figures.

    Args:
        deseq_all: All DESeq2 data by condition
        deseq_filtered: Filtered DESeq2 data by condition
        rmats_all: All rMATS data by condition
        rmats_filtered: Filtered rMATS data by condition
        config: Configuration dictionary
        output_dir: Output directory for figures
    """
    fig_config = config.get('figures', {})
    fig_generate = fig_config.get('generate', {})
    deseq_cutoffs = config['filtering']['deseq2']
    splicing_cutoffs = config['filtering']['splicing']

    figures_dir = output_dir / 'figures'
    figures_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Generating figures...")

    # Per-condition figures
    for condition in config['conditions']:
        condition_name = condition['name']
        condition_label = condition['label']

        df_all = deseq_all.get(condition_name)
        df_sig = deseq_filtered.get(condition_name)

        if df_all is not None and df_sig is not None:
            # Volcano plot
            if fig_generate.get('volcano', True):
                try:
                    figures.plot_volcano(
                        df_all=df_all,
                        df_sig=df_sig,
                        label=condition_label,
                        output_path=figures_dir / f'{condition_name}_volcano.png',
                        log2fc_cutoff=deseq_cutoffs['log2fc_cutoff'],
                        padj_cutoff=deseq_cutoffs['padj_cutoff']
                    )
                except Exception as e:
                    logging.error(f"Failed to generate volcano plot for {condition_name}: {e}")

            # MA plot
            if fig_generate.get('ma_plot', True):
                try:
                    figures.plot_ma(
                        df_all=df_all,
                        df_sig=df_sig,
                        label=condition_label,
                        output_path=figures_dir / f'{condition_name}_ma.png',
                        log2fc_cutoff=deseq_cutoffs['log2fc_cutoff']
                    )
                except Exception as e:
                    logging.error(f"Failed to generate MA plot for {condition_name}: {e}")

            # Biotype distribution
            if fig_generate.get('biotype', True):
                try:
                    figures.plot_biotype_distribution(
                        df_sig=df_sig,
                        label=condition_label,
                        output_path=figures_dir / f'{condition_name}_biotype.png'
                    )
                except Exception as e:
                    logging.error(f"Failed to generate biotype plot for {condition_name}: {e}")

        # Splicing figures
        rmats_all_cond = rmats_all.get(condition_name)
        rmats_sig = rmats_filtered.get(condition_name)

        if rmats_all_cond is not None and rmats_sig is not None:
            # Splicing summary
            if fig_generate.get('splicing_summary', True):
                try:
                    figures.plot_splicing_summary(
                        rmats_df=rmats_sig,
                        label=condition_label,
                        output_path=figures_dir / f'{condition_name}_splicing_summary.png'
                    )
                except Exception as e:
                    logging.error(f"Failed to generate splicing summary for {condition_name}: {e}")

            # Splicing volcano
            if fig_generate.get('splicing_volcano', True):
                try:
                    figures.plot_splicing_volcano(
                        rmats_df=rmats_all_cond,
                        label=condition_label,
                        output_path=figures_dir / f'{condition_name}_splicing_volcano.png',
                        use_pval=splicing_cutoffs.get('use_pval_for_volcano', True),
                        dpsi_cutoff=splicing_cutoffs['dpsi_cutoff'],
                        fdr_cutoff=splicing_cutoffs['fdr_cutoff']
                    )
                except Exception as e:
                    logging.error(f"Failed to generate splicing volcano for {condition_name}: {e}")

    # Cross-condition figures
    if len(config['conditions']) > 1:
        # Log2FC violin plots
        if fig_generate.get('violin_plots', True):
            try:
                figures.plot_log2fc_violin(
                    filtered_dict=deseq_filtered,
                    output_path=figures_dir / 'log2fc_violin.png',
                    conditions=config['conditions']
                )
            except Exception as e:
                logging.error(f"Failed to generate log2FC violin plot: {e}")

        # dPSI violin plots
        if rmats_filtered and fig_generate.get('violin_plots', True):
            try:
                figures.plot_dpsi_violin(
                    rmats_dict=rmats_filtered,
                    output_path=figures_dir / 'dpsi_violin.png',
                    conditions=config['conditions']
                )
            except Exception as e:
                logging.error(f"Failed to generate dPSI violin plot: {e}")

        # Concordance scatter plots
        if fig_generate.get('concordance_scatter', True):
            try:
                figures.plot_concordance_all_pairs(
                    all_data=deseq_all,
                    filtered_dict=deseq_filtered,
                    output_path=figures_dir,
                    conditions=config['conditions']
                )
            except Exception as e:
                logging.error(f"Failed to generate concordance scatter plots: {e}")

        # Heatmap
        if fig_generate.get('heatmaps', True):
            try:
                figures.plot_top_genes_heatmap(
                    all_data=deseq_all,
                    filtered_dict=deseq_filtered,
                    output_path=figures_dir / 'top_genes_heatmap.png',
                    conditions=config['conditions'],
                    n_top=50
                )
            except Exception as e:
                logging.error(f"Failed to generate heatmap: {e}")

    logging.info(f"Figures saved to {figures_dir}/")


# ─── Cross-Condition Concordance ─────────────────────────────────────────────

def run_concordance_analysis(
    deseq_filtered: Dict[str, pd.DataFrame],
    rmats_filtered: Dict[str, pd.DataFrame],
    config: dict,
    output_dir: Path
) -> None:
    """
    Run cross-condition concordance analysis for all comparison pairs.

    Args:
        deseq_filtered: Filtered DESeq2 data by condition
        rmats_filtered: Filtered rMATS data by condition
        config: Configuration dictionary
        output_dir: Output directory for results
    """
    if 'comparisons' not in config or not config['comparisons']:
        logging.info("No cross-condition comparisons defined in config - skipping concordance")
        return

    concordance_dir = output_dir / 'concordance'
    concordance_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Running cross-condition concordance analysis...")

    for comparison in config['comparisons']:
        comp_name = comparison['name']
        cond1_name = comparison['condition1']
        cond2_name = comparison['condition2']

        # Get condition labels
        cond1_label = next(c['label'] for c in config['conditions'] if c['name'] == cond1_name)
        cond2_label = next(c['label'] for c in config['conditions'] if c['name'] == cond2_name)

        logging.info(f"  Comparing {cond1_label} vs {cond2_label}...")

        # DEG concordance
        if cond1_name in deseq_filtered and cond2_name in deseq_filtered:
            df_a = deseq_filtered[cond1_name]
            df_b = deseq_filtered[cond2_name]

            # Directional Venn diagrams
            try:
                concordance.plot_directional_venn_3panel(
                    df_a=df_a,
                    df_b=df_b,
                    label_a=cond1_label,
                    label_b=cond2_label,
                    output_path=concordance_dir / f'{comp_name}_deg_venn.png'
                )
            except Exception as e:
                logging.error(f"Failed to generate DEG Venn for {comp_name}: {e}")

            # Export overlap Excel
            try:
                concordance.export_directional_overlap_excel(
                    df_a=df_a,
                    df_b=df_b,
                    cond_a=cond1_label,
                    cond_b=cond2_label,
                    output_path=concordance_dir / f'{comp_name}_deg_overlap.xlsx'
                )
            except Exception as e:
                logging.error(f"Failed to export DEG overlap for {comp_name}: {e}")

        # Splicing concordance
        if cond1_name in rmats_filtered and cond2_name in rmats_filtered:
            rmats_a = rmats_filtered[cond1_name]
            rmats_b = rmats_filtered[cond2_name]

            # Venn for each event type
            for event_type in ['SE', 'A3SS', 'A5SS', 'MXE', 'RI']:
                try:
                    concordance.plot_splicing_venn_3panel(
                        rmats_a=rmats_a,
                        rmats_b=rmats_b,
                        label_a=cond1_label,
                        label_b=cond2_label,
                        event_type=event_type,
                        output_path=concordance_dir / f'{comp_name}_splicing_{event_type}_venn.png'
                    )
                except Exception as e:
                    logging.warning(f"Failed to generate {event_type} Venn for {comp_name}: {e}")

    logging.info(f"Concordance analysis saved to {concordance_dir}/")


# ─── Excel Export ────────────────────────────────────────────────────────────

def export_excel_results(
    deseq_filtered: Dict[str, pd.DataFrame],
    rmats_filtered: Dict[str, pd.DataFrame],
    config: dict,
    output_dir: Path
) -> None:
    """
    Export filtered results to Excel files.

    Args:
        deseq_filtered: Filtered DESeq2 data by condition
        rmats_filtered: Filtered rMATS data by condition
        config: Configuration dictionary
        output_dir: Output directory for Excel files
    """
    excel_dir = output_dir / 'excel'
    excel_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Exporting Excel files...")

    for condition in config['conditions']:
        condition_name = condition['name']
        condition_label = condition['label']

        # Export DESeq2 results
        if condition_name in deseq_filtered:
            try:
                output_path = excel_dir / f'{condition_name}_DEGs.xlsx'
                deseq2.export_deg_results(
                    df_sig=deseq_filtered[condition_name],
                    output_path=output_path
                )
                logging.info(f"  Exported {condition_label} DEGs to {output_path.name}")
            except Exception as e:
                logging.error(f"Failed to export DEGs for {condition_name}: {e}")

        # Export splicing results
        if condition_name in rmats_filtered:
            try:
                output_path = excel_dir / f'{condition_name}_splicing.xlsx'
                splicing.export_splicing_events(
                    df=rmats_filtered[condition_name],
                    output_path=output_path,
                    by_event_type=True
                )
                logging.info(f"  Exported {condition_label} splicing events to {output_path.name}")
            except Exception as e:
                logging.error(f"Failed to export splicing for {condition_name}: {e}")


# ─── Main Pipeline ───────────────────────────────────────────────────────────

def run_pipeline(
    config_path: Path,
    output_dir: Optional[Path] = None,
    condition_subset: Optional[List[str]] = None,
    skip_figures: bool = False,
    modules: Optional[List[str]] = None
) -> None:
    """
    Execute the complete RNA-seq analysis pipeline.

    Args:
        config_path: Path to rnaseq_config.yaml
        output_dir: Override output directory from config
        condition_subset: Optional list of condition names to analyze
        skip_figures: Skip figure generation if True
        modules: List of modules to run (deseq2, splicing, concordance)
    """
    # Load and validate config
    config = load_config(config_path)
    validate_config(config)

    # Determine output directory
    if output_dir is None:
        output_dir = Path(config['output_dir'])
    else:
        output_dir = Path(output_dir)

    # Setup logging
    log_config = config.get('logging', {})
    setup_logging(
        output_dir=output_dir,
        log_level=log_config.get('level', 'INFO'),
        log_file=log_config.get('file', 'pipeline.log')
    )

    # Determine which modules to run
    if modules is None:
        module_config = config.get('modules', {})
        run_deseq = module_config.get('deseq2', True)
        run_splicing = module_config.get('splicing', True)
        run_concordance = module_config.get('concordance', True)
    else:
        run_deseq = 'deseq2' in modules
        run_splicing = 'splicing' in modules
        run_concordance = 'concordance' in modules

    # Pipeline start
    start_time = time.time()
    logging.info("="*80)
    logging.info("RNA-seq Analysis Pipeline - Starting")
    logging.info("="*80)
    logging.info(f"Config: {config_path}")
    logging.info(f"Output: {output_dir}")
    logging.info(f"Modules: DESeq2={run_deseq}, Splicing={run_splicing}, Concordance={run_concordance}")
    logging.info("")

    # Load data
    deseq_all, rmats_all = load_all_conditions(
        config=config,
        condition_subset=condition_subset,
        load_deseq=run_deseq,
        load_splicing=run_splicing
    )

    # Filter data
    deseq_filtered, rmats_filtered = filter_all_conditions(
        deseq_dict=deseq_all,
        rmats_dict=rmats_all,
        config=config
    )

    # Generate summary statistics
    if config.get('export', {}).get('summary_stats', True):
        logging.info("Generating summary statistics...")
        summary_df = generate_summary_stats(
            deseq_all=deseq_all,
            deseq_filtered=deseq_filtered,
            rmats_all=rmats_all,
            rmats_filtered=rmats_filtered,
            config=config
        )
        summary_path = output_dir / 'summary_statistics.csv'
        summary_df.to_csv(summary_path, index=False)
        logging.info(f"Summary statistics saved to {summary_path}")
        logging.info("")
        logging.info(summary_df.to_string(index=False))
        logging.info("")

    # Export Excel files
    if config.get('export', {}).get('excel', True):
        export_excel_results(
            deseq_filtered=deseq_filtered,
            rmats_filtered=rmats_filtered,
            config=config,
            output_dir=output_dir
        )

    # Generate figures
    if not skip_figures and config.get('export', {}).get('figures', True):
        generate_all_figures(
            deseq_all=deseq_all,
            deseq_filtered=deseq_filtered,
            rmats_all=rmats_all,
            rmats_filtered=rmats_filtered,
            config=config,
            output_dir=output_dir
        )

    # Cross-condition concordance
    if run_concordance:
        run_concordance_analysis(
            deseq_filtered=deseq_filtered,
            rmats_filtered=rmats_filtered,
            config=config,
            output_dir=output_dir
        )

    # Pipeline complete
    elapsed_time = time.time() - start_time
    logging.info("")
    logging.info("="*80)
    logging.info(f"Pipeline completed successfully in {elapsed_time:.1f} seconds")
    logging.info("="*80)
    logging.info(f"Results written to: {output_dir.resolve()}")


# ─── Command-Line Interface ──────────────────────────────────────────────────

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='RNA-seq Analysis Pipeline - Unified CLI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline with default config
  python run_pipeline.py

  # Use custom config file
  python run_pipeline.py --config my_config.yaml

  # Run only specific conditions
  python run_pipeline.py --conditions MIAT_OE_vs_Control,QKI_KO_vs_WT

  # Run only DESeq2 analysis, skip splicing
  python run_pipeline.py --modules deseq2

  # Skip figure generation
  python run_pipeline.py --skip-figures

  # Override output directory
  python run_pipeline.py --output-dir results_2026_03_06/
        """
    )

    parser.add_argument(
        '--config',
        type=Path,
        default=Path('rnaseq_config.yaml'),
        help='Path to configuration YAML file (default: rnaseq_config.yaml)'
    )

    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        help='Override output directory from config'
    )

    parser.add_argument(
        '--conditions',
        type=str,
        default=None,
        help='Comma-separated list of condition names to analyze (default: all)'
    )

    parser.add_argument(
        '--skip-figures',
        action='store_true',
        help='Skip figure generation'
    )

    parser.add_argument(
        '--modules',
        type=str,
        default=None,
        help='Comma-separated list of modules to run: deseq2,splicing,concordance (default: all)'
    )

    return parser.parse_args()


def main():
    """Main entry point."""
    # Suppress warnings
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=UserWarning)

    # Parse arguments
    args = parse_args()

    # Parse condition subset
    condition_subset = None
    if args.conditions:
        condition_subset = [c.strip() for c in args.conditions.split(',')]

    # Parse modules
    modules = None
    if args.modules:
        modules = [m.strip() for m in args.modules.split(',')]

    # Run pipeline
    try:
        run_pipeline(
            config_path=args.config,
            output_dir=args.output_dir,
            condition_subset=condition_subset,
            skip_figures=args.skip_figures,
            modules=modules
        )
    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
