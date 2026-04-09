"""
Command-line interface for the ``fl_ai_reagent_stocker`` package.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from .config import Settings
from .pipelines.stock_finding import StockFindingPipeline
from .pipelines.stock_splitting import StockSplittingPipeline, load_split_config


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fl-ai-reagent-stocker",
        description="Drosophila stock finding, splitting, and validation pipelines",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m fl_ai_reagent_stocker find-stocks ./gene_lists --config ./my_config.json
    python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_config.json
    python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_config.json
    python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_config.json
""",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    find_parser = subparsers.add_parser(
        "find-stocks",
        help="Stage 1: map genes to FlyBase stocks and references",
    )
    find_parser.add_argument("input_dir", type=Path, help="Directory containing CSV gene lists")
    find_parser.add_argument("--gene-col", default="flybase_gene_id")
    find_parser.add_argument("--input-gene-col", default="ext_gene")
    find_parser.add_argument("--config", "-c", type=Path, default=None)
    find_parser.add_argument("--batch-size", "-b", type=int, default=50)
    find_parser.add_argument("--skip-fbgnid-conversion", action="store_true")

    split_parser = subparsers.add_parser(
        "split-stocks",
        help="Stage 2: organize stocks using JSON filters only",
    )
    split_parser.add_argument("input_dir", type=Path, help="Directory containing Stage 1 Excel files")
    split_parser.add_argument("--config", "-c", type=Path, default=None)
    split_parser.add_argument("--quiet", "-q", action="store_true")
    split_parser.add_argument(
        "--soft-run",
        action="store_true",
        help="Produce the soft-run phenotype sheet instead of normal downstream summary sheets.",
    )
    split_parser.add_argument(
        "--OAI-embedding",
        dest="oai_embedding",
        action="store_true",
        help="When combined with --soft-run, add OpenAI embedding cosine similarity to the Stock Phenotype Sheet.",
    )
    split_parser.add_argument(
        "--simple-buckets",
        action="store_true",
        help="When combined with --soft-run --OAI-embedding, replace cosine-threshold tier tabs with rule-based combination tabs in the phenotype similarity workbook.",
    )

    validate_parser = subparsers.add_parser(
        "validate-stocks",
        help="Stage 3: append GPT validation results to Ref++ output sheets",
    )
    validate_parser.add_argument("input_dir", type=Path, help="Directory containing Stage 1 Excel files")
    validate_parser.add_argument("--config", "-c", type=Path, default=None)
    validate_parser.add_argument("--quiet", "-q", action="store_true")
    validate_parser.add_argument("--soft-run", action="store_true")
    validate_parser.add_argument(
        "--OAI-embedding",
        dest="oai_embedding",
        action="store_true",
        help="When combined with --soft-run, add OpenAI embedding cosine similarity to the Stock Phenotype Sheet.",
    )
    validate_parser.add_argument(
        "--simple-buckets",
        action="store_true",
        help="When combined with --soft-run --OAI-embedding, replace cosine-threshold tier tabs with rule-based combination tabs in the phenotype similarity workbook.",
    )
    validate_parser.add_argument("--test-log", action="store_true")
    validate_parser.add_argument("--max-gpt-calls-per-stock", type=int, default=5)

    full_parser = subparsers.add_parser(
        "run-full-pipeline",
        help="Run stock finding, splitting, and validation end-to-end",
    )
    full_parser.add_argument("input_dir", type=Path, help="Directory containing CSV gene lists")
    full_parser.add_argument("--config", "-c", type=Path, default=None)
    full_parser.add_argument("--gene-col", default="flybase_gene_id")
    full_parser.add_argument("--input-gene-col", default="ext_gene")
    full_parser.add_argument("--batch-size", "-b", type=int, default=50)
    full_parser.add_argument("--skip-fbgnid-conversion", action="store_true")
    full_parser.add_argument("--quiet", "-q", action="store_true")
    full_parser.add_argument("--soft-run", action="store_true")
    full_parser.add_argument(
        "--OAI-embedding",
        dest="oai_embedding",
        action="store_true",
        help="When combined with --soft-run, add OpenAI embedding cosine similarity to the Stock Phenotype Sheet.",
    )
    full_parser.add_argument(
        "--simple-buckets",
        action="store_true",
        help="When combined with --soft-run --OAI-embedding, replace cosine-threshold tier tabs with rule-based combination tabs in the phenotype similarity workbook.",
    )
    full_parser.add_argument("--test-log", action="store_true")
    full_parser.add_argument("--max-gpt-calls-per-stock", type=int, default=5)

    return parser


def _load_keywords(config_path: Path) -> list[str]:
    config = load_split_config(config_path)
    return config.get("settings", {}).get("relevantSearchTerms", [])


def run_find_stocks(args) -> int:
    settings = Settings(batch_size=args.batch_size, skip_fbgnid_conversion=args.skip_fbgnid_conversion)
    pipeline = StockFindingPipeline(settings)
    config_path = args.config if args.config else settings.split_config_path
    keywords = _load_keywords(config_path)

    output_path = pipeline.run(
        input_dir=args.input_dir,
        keywords=keywords,
        gene_col=args.gene_col,
        input_gene_col=args.input_gene_col,
        skip_fbgnid_conversion=args.skip_fbgnid_conversion,
        run_functional_validation=False,
    )
    return 0 if output_path else 1


def run_split_stocks(args) -> int:
    settings = Settings(
        soft_run=args.soft_run,
        enable_oai_embedding=args.oai_embedding,
        simple_buckets=args.simple_buckets,
    )
    if args.config:
        settings.split_config_path = args.config
    pipeline = StockSplittingPipeline(settings)
    output_dir = pipeline.run(
        input_dir=args.input_dir,
        config_path=args.config,
        verbose=not args.quiet,
    )
    return 0 if output_dir else 1


def run_validate_stocks(args) -> int:
    settings = Settings(
        soft_run=args.soft_run,
        enable_oai_embedding=args.oai_embedding,
        simple_buckets=args.simple_buckets,
        enable_gpt_logging=args.test_log,
        max_gpt_calls_per_stock=args.max_gpt_calls_per_stock,
    )
    if args.config:
        settings.split_config_path = args.config
    pipeline = StockSplittingPipeline(settings)
    output_dir = pipeline.run(
        input_dir=args.input_dir,
        config_path=args.config,
        verbose=not args.quiet,
        run_validation=True,
    )
    return 0 if output_dir else 1


def run_full_pipeline(args) -> int:
    refs_settings = Settings(
        batch_size=args.batch_size,
        skip_fbgnid_conversion=args.skip_fbgnid_conversion,
    )
    split_settings = Settings(
        soft_run=args.soft_run,
        enable_oai_embedding=args.oai_embedding,
        simple_buckets=args.simple_buckets,
    )
    validate_settings = Settings(
        soft_run=args.soft_run,
        enable_oai_embedding=args.oai_embedding,
        simple_buckets=args.simple_buckets,
        enable_gpt_logging=args.test_log,
        max_gpt_calls_per_stock=args.max_gpt_calls_per_stock,
    )

    if args.config:
        split_settings.split_config_path = args.config
        validate_settings.split_config_path = args.config

    config_path = args.config if args.config else refs_settings.split_config_path
    keywords = _load_keywords(config_path)

    refs_pipeline = StockFindingPipeline(refs_settings)
    split_pipeline = StockSplittingPipeline(split_settings)
    validate_pipeline = StockSplittingPipeline(validate_settings)

    print(f"\n{'=' * 70}")
    print("FL.AI REAGENT STOCKER: FULL PIPELINE")
    print(f"{'=' * 70}")

    print(f"\n{'-' * 70}")
    print("Step 1/3: Running find-stocks")
    print(f"{'-' * 70}")
    refs_output_path = refs_pipeline.run(
        input_dir=args.input_dir,
        keywords=keywords,
        gene_col=args.gene_col,
        input_gene_col=args.input_gene_col,
        skip_fbgnid_conversion=args.skip_fbgnid_conversion,
        run_functional_validation=False,
    )
    if not refs_output_path:
        return 1

    split_input_dir = Path(refs_output_path).parent

    print(f"\n{'-' * 70}")
    print("Step 2/3: Running split-stocks")
    print(f"{'-' * 70}")
    split_output_dir = split_pipeline.run(
        input_dir=split_input_dir,
        config_path=args.config,
        verbose=not args.quiet,
    )
    if not split_output_dir:
        return 1

    print(f"\n{'-' * 70}")
    print("Step 3/3: Running validate-stocks")
    print(f"{'-' * 70}")
    validated_output_dir = validate_pipeline.run(
        input_dir=split_input_dir,
        config_path=args.config,
        verbose=not args.quiet,
        run_validation=True,
    )
    if not validated_output_dir:
        return 1

    print("\nSuccess! Full pipeline completed.")
    print(f"  Stage 1 output: {refs_output_path}")
    print(f"  Stage 2 output directory: {split_output_dir}")
    print(f"  Stage 3 output directory: {validated_output_dir}")
    return 0


def main(argv: Optional[list] = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "find-stocks":
        return run_find_stocks(args)
    if args.command == "split-stocks":
        return run_split_stocks(args)
    if args.command == "validate-stocks":
        return run_validate_stocks(args)
    if args.command == "run-full-pipeline":
        return run_full_pipeline(args)

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
