"""
Command-line interface for fly-stocker-v2 pipelines.

Provides entry points for:
    - get-allele-refs: Pipeline 1 (Gene → Stocks → References)
    - split-stocks: Pipeline 2 (Organize using JSON configuration + Ref++-scoped validation)

Usage:
    # Pipeline 1: Get allele references
    python -m fly_stocker_v2.cli get-allele-refs /path/to/input --config /path/to/config.json
    
    # Pipeline 2: Split stocks using JSON config
    python -m fly_stocker_v2.cli split-stocks /path/to/input --config /path/to/config.json
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

try:
    from .config import Settings
    from .pipeline_references import AlleleReferencesPipeline
    from .pipeline_split import StockSplitterPipeline, load_split_config
except ImportError:
    # Support flat-repo execution (for example `python -m cli ...`).
    from config import Settings
    from pipeline_references import AlleleReferencesPipeline
    from pipeline_split import StockSplitterPipeline, load_split_config


def create_parser() -> argparse.ArgumentParser:
    """Create the main argument parser."""
    parser = argparse.ArgumentParser(
        prog="fly-stocker-v2",
        description="Modular Drosophila Stock Processing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Pipeline 1: Get allele references for genes
    python -m fly_stocker_v2.cli get-allele-refs ./gene_lists --config ./my_config.json
    
    # Pipeline 2: Organize stocks using JSON config
    python -m fly_stocker_v2.cli split-stocks ./gene_lists/Stocks --config ./my_config.json

    # Full pipeline: Pipeline 1 then Pipeline 2
    python -m fly_stocker_v2.cli run-full-pipeline ./gene_lists --config ./my_config.json
    
    # Pipeline 2: Use default config
    python -m fly_stocker_v2.cli split-stocks ./gene_lists/Stocks
    
    # Run Pipeline 1 validation in soft-run mode (count queries without calling OpenAI)
    python -m fly_stocker_v2.cli get-allele-refs ./gene_lists --run-validation --soft-run
"""
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # =========================================================================
    # Pipeline 1: Get Allele References
    # =========================================================================
    refs_parser = subparsers.add_parser(
        "get-allele-refs",
        help="Pipeline 1: Map genes to BDSC stocks and get references",
        description="""
Get BDSC stock-references chain for genes.

This pipeline:
1. Converts gene symbols to FlyBase gene IDs (FBgn)
2. Maps genes to BDSC stocks via stockgenes.csv
3. Gets references for each stock's alleles
4. Uses split JSON relevantSearchTerms for keyword matching
5. Outputs Excel file with stocks and references
        """
    )
    
    refs_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing CSV files with gene IDs"
    )
    refs_parser.add_argument(
        "--gene-col",
        default="flybase_gene_id",
        help="Column name for FlyBase gene IDs (after conversion). Default: flybase_gene_id"
    )
    refs_parser.add_argument(
        "--input-gene-col",
        default="ext_gene",
        help="Column name for gene symbols in input (for conversion). Default: ext_gene"
    )
    refs_parser.add_argument(
        "--config", "-c",
        type=Path,
        default=None,
        help="Path to JSON configuration file used for keywords. "
             "If not provided, uses default split config."
    )
    refs_parser.add_argument(
        "--batch-size", "-b",
        type=int,
        default=50,
        help="Number of genes to process per batch. Default: 50"
    )
    refs_parser.add_argument(
        "--skip-fbgnid-conversion",
        action="store_true",
        help="Skip FBgnID conversion step (use if input already has flybase_gene_id)"
    )
    refs_parser.add_argument(
        "--run-validation",
        action="store_true",
        help="Run GPT functional validation during Pipeline 1. "
             "By default, validation is skipped here and performed in split-stocks."
    )
    refs_parser.add_argument(
        "--soft-run",
        action="store_true",
        help="Run up to but not including OpenAI API calls. Prints expected query count. "
             "Requires validation to be enabled for the command."
    )
    refs_parser.add_argument(
        "--test-log",
        action="store_true",
        help="Log all GPT queries to data/logs/GPT-Queries/ directory"
    )
    refs_parser.add_argument(
        "--max-gpt-calls-per-stock",
        type=int,
        default=5,
        help="Maximum OpenAI GPT API calls per stock. "
             "Only actual GPT calls count against this limit - references without accessible "
             "full text or without stock/allele pattern matches are marked ambiguous and don't count. "
             "References are processed most recent first. "
             "Validation short-circuits when 'Functionally validated' is found. Default: 5"
    )
    
    # =========================================================================
    # Pipeline 2: Split Stocks using JSON Configuration
    # =========================================================================
    split_parser = subparsers.add_parser(
        "split-stocks",
        help="Pipeline 2: Organize stocks using JSON configuration into formatted Excel sheets",
        description="""
Organize stocks using JSON-based filter configurations.

This pipeline:
1. Loads stock data from Excel files (output of Pipeline 1)
2. Computes derived columns (Balancers via structural genotype parsing, multiple_insertions via component data, ALLELE_PAPER_RELEVANCE_SCORE via title/abstract keyword hits)
3. Applies configurable filters from JSON to create sheets
4. Creates aggregated Excel output with Contents, Sheet1..N, References sheets

The JSON configuration defines:
- relevantSearchTerms: Keywords for determining Ref++ status
- maxStocksPerGene/maxStocksPerAllele: Limits for stock selection
- filters: Named filter specifications (column, type, value)
- combinations: Lists of filter names to AND together for each sheet
        """
    )
    
    split_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing Excel files from Pipeline 1"
    )
    split_parser.add_argument(
        "--config", "-c",
        type=Path,
        default=None,
        help="Path to JSON configuration file. If not provided, uses the default config at data/JSON/stock_split_config_example.json"
    )
    split_parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress verbose output"
    )
    split_parser.add_argument(
        "--soft-run",
        action="store_true",
        help="Run up to but not including OpenAI API calls. Prints expected query count."
    )
    split_parser.add_argument(
        "--test-log",
        action="store_true",
        help="Log all GPT queries to data/logs/GPT-Queries/ directory"
    )
    split_parser.add_argument(
        "--max-gpt-calls-per-stock",
        type=int,
        default=5,
        help="Maximum OpenAI GPT API calls per stock. "
             "Only actual GPT calls count against this limit - references without accessible "
             "full text or without stock/allele pattern matches are marked ambiguous and don't count. "
             "References are processed most recent first. "
             "Validation short-circuits when 'Functionally validated' is found. Default: 5"
    )

    # =========================================================================
    # Full Pipeline: Pipeline 1 then Pipeline 2
    # =========================================================================
    full_parser = subparsers.add_parser(
        "run-full-pipeline",
        help="Run Pipeline 1 then Pipeline 2 end-to-end",
        description="""
Run both pipelines in one command.

This command:
1. Runs get-allele-refs on an input gene-list directory
2. Uses the generated Stocks directory as input to split-stocks
3. Produces final organized outputs in Organized Stocks
        """
    )
    full_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing CSV files with gene IDs"
    )
    full_parser.add_argument(
        "--config", "-c",
        type=Path,
        default=None,
        help="Path to JSON configuration file shared by both pipeline steps. "
             "If not provided, uses the default split config."
    )
    full_parser.add_argument(
        "--gene-col",
        default="flybase_gene_id",
        help="Column name for FlyBase gene IDs (after conversion). Default: flybase_gene_id"
    )
    full_parser.add_argument(
        "--input-gene-col",
        default="ext_gene",
        help="Column name for gene symbols in input (for conversion). Default: ext_gene"
    )
    full_parser.add_argument(
        "--batch-size", "-b",
        type=int,
        default=50,
        help="Number of genes to process per batch for Pipeline 1. Default: 50"
    )
    full_parser.add_argument(
        "--skip-fbgnid-conversion",
        action="store_true",
        help="Skip FBgnID conversion in Pipeline 1 (use if input already has flybase_gene_id)"
    )
    full_parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress verbose output for Pipeline 2"
    )
    full_parser.add_argument(
        "--soft-run",
        action="store_true",
        help="Run up to but not including OpenAI API calls where validation would occur."
    )
    full_parser.add_argument(
        "--test-log",
        action="store_true",
        help="Log all GPT queries to data/logs/GPT-Queries/ directory"
    )
    full_parser.add_argument(
        "--max-gpt-calls-per-stock",
        type=int,
        default=5,
        help="Maximum OpenAI GPT API calls per stock during validation. Default: 5"
    )
    
    return parser


def run_get_allele_refs(args) -> int:
    """Run Pipeline 1: Get Allele References."""
    print(f"\n{'='*70}")
    print("FLY-STOCKER v2: Get Allele References Pipeline")
    print(f"{'='*70}")
    
    # Create settings
    settings = Settings(
        batch_size=args.batch_size,
        soft_run=args.soft_run,
        skip_fbgnid_conversion=args.skip_fbgnid_conversion,
        enable_gpt_logging=args.test_log,
        max_gpt_calls_per_stock=args.max_gpt_calls_per_stock,
    )
    
    # Create and run pipeline
    pipeline = AlleleReferencesPipeline(settings)
    
    try:
        # Load keywords from split configuration
        config_path = args.config if args.config else settings.split_config_path
        config = load_split_config(config_path)
        keywords = config.get("settings", {}).get("relevantSearchTerms", [])
        print(f"Using split config for keywords: {config_path}")
        if keywords:
            print(f"  relevantSearchTerms: {', '.join(keywords)}")
        else:
            print("  Warning: relevantSearchTerms is empty; no title/abstract keyword hits will be detected.")

        if not args.run_validation and (args.soft_run or args.test_log):
            print("  Note: --soft-run/--test-log only apply when --run-validation is enabled.")

        output_path = pipeline.run(
            input_dir=args.input_dir,
            keywords=keywords,
            gene_col=args.gene_col,
            input_gene_col=args.input_gene_col,
            skip_fbgnid_conversion=args.skip_fbgnid_conversion,
            run_functional_validation=args.run_validation,
        )
        
        if output_path:
            print(f"\nSuccess! Output: {output_path}")
            return 0
        else:
            print("\nPipeline completed with errors.")
            return 1
            
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return 1


def run_split_stocks(args) -> int:
    """Run Pipeline 2: Split Stocks using JSON Configuration."""
    # Create settings
    settings = Settings(
        soft_run=args.soft_run,
        enable_gpt_logging=args.test_log,
        max_gpt_calls_per_stock=args.max_gpt_calls_per_stock,
    )
    
    # Override config path if provided
    if args.config:
        settings.split_config_path = args.config
    
    # Create and run pipeline
    pipeline = StockSplitterPipeline(settings)
    
    try:
        output_dir = pipeline.run(
            input_dir=args.input_dir,
            config_path=args.config,
            verbose=not args.quiet,
        )
        
        if output_dir:
            print(f"\nSuccess! Output directory: {output_dir}")
            return 0
        else:
            print("\nPipeline completed with errors.")
            return 1
            
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return 1


def run_full_pipeline(args) -> int:
    """Run Pipeline 1 then Pipeline 2 in one command."""
    print(f"\n{'='*70}")
    print("FLY-STOCKER v2: Full Pipeline (Pipeline 1 -> Pipeline 2)")
    print(f"{'='*70}")

    # Step 1 settings (Pipeline 1)
    refs_settings = Settings(
        batch_size=args.batch_size,
        soft_run=args.soft_run,
        skip_fbgnid_conversion=args.skip_fbgnid_conversion,
        enable_gpt_logging=args.test_log,
        max_gpt_calls_per_stock=args.max_gpt_calls_per_stock,
    )

    # Step 2 settings (Pipeline 2)
    split_settings = Settings(
        soft_run=args.soft_run,
        enable_gpt_logging=args.test_log,
        max_gpt_calls_per_stock=args.max_gpt_calls_per_stock,
    )

    if args.config:
        split_settings.split_config_path = args.config

    refs_pipeline = AlleleReferencesPipeline(refs_settings)
    split_pipeline = StockSplitterPipeline(split_settings)

    try:
        # Load keywords from split configuration so both stages use same config.
        config_path = args.config if args.config else refs_settings.split_config_path
        config = load_split_config(config_path)
        keywords = config.get("settings", {}).get("relevantSearchTerms", [])
        print(f"Using shared config: {config_path}")
        if keywords:
            print(f"  relevantSearchTerms: {', '.join(keywords)}")
        else:
            print("  Warning: relevantSearchTerms is empty; no title/abstract keyword hits will be detected.")

        # Pipeline 1 (default behavior: no GPT validation here)
        print(f"\n{'-'*70}")
        print("Step 1/2: Running get-allele-refs")
        print(f"{'-'*70}")
        refs_output_path = refs_pipeline.run(
            input_dir=args.input_dir,
            keywords=keywords,
            gene_col=args.gene_col,
            input_gene_col=args.input_gene_col,
            skip_fbgnid_conversion=args.skip_fbgnid_conversion,
            run_functional_validation=False,
        )
        if not refs_output_path:
            print("\nPipeline 1 failed. Stopping full pipeline.")
            return 1

        split_input_dir = Path(refs_output_path).parent

        # Pipeline 2
        print(f"\n{'-'*70}")
        print("Step 2/2: Running split-stocks")
        print(f"{'-'*70}")
        split_output_dir = split_pipeline.run(
            input_dir=split_input_dir,
            config_path=args.config,
            verbose=not args.quiet,
        )
        if not split_output_dir:
            print("\nPipeline 2 failed.")
            return 1

        print("\nSuccess! Full pipeline completed.")
        print(f"  Pipeline 1 output: {refs_output_path}")
        print(f"  Pipeline 2 output directory: {split_output_dir}")
        return 0

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return 1


def main(argv: Optional[list] = None) -> int:
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args(argv)
    
    if args.command is None:
        parser.print_help()
        return 0
    
    if args.command == "get-allele-refs":
        return run_get_allele_refs(args)
    elif args.command == "split-stocks":
        return run_split_stocks(args)
    elif args.command == "run-full-pipeline":
        return run_full_pipeline(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
