"""
Command-line interface for the ``fl_ai_reagent_stocker`` package.

Public commands:

- ``find-stocks``: Stage 1 stock and reference discovery; persists the
  ``Gene Reagent Index`` sheet alongside ``Stocks`` and ``References``.
- ``split-stocks``: Stage 2 JSON stock organization only. Does not produce
  the phenotype-sheet workbook.
- ``phenotype-sheet``: gene-first phenotype reagent flow over Stage 1
  workbook(s). Selects the optional similarity sidecar layout via
  ``--similarity``.
- ``validate-stocks``: Stage 3 GPT validation only. Does not produce the
  phenotype-sheet workbook.
- ``run-full-pipeline``: end-to-end wrapper from raw gene-list CSVs.
  ``--mode normal`` runs Stage 1 + ``split-stocks`` + ``validate-stocks``.
  ``--mode phenotype`` runs Stage 1 + ``phenotype-sheet`` only and never
  invokes GPT validation.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd

from .config import Settings
from .pipelines.stock_finding import StockFindingPipeline
from .pipelines.stock_splitting import StockSplittingPipeline, load_split_config


SIMILARITY_CHOICES = ("none", "tiers", "simple-buckets", "keyword-buckets")
PIPELINE_MODE_CHOICES = ("normal", "phenotype")


###############################################################################
# Argument helpers
###############################################################################


def _add_config_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        "-c",
        type=Path,
        default=None,
        help="Path to a JSON split-config (defaults to data/config/stock_split_config_example.json).",
    )


def _add_quiet_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress per-file progress logging.",
    )


def _add_stage1_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--gene-col",
        default="flybase_gene_id",
        help="Name of the FBgn ID column in input CSVs (default: flybase_gene_id).",
    )
    parser.add_argument(
        "--input-gene-col",
        default="ext_gene",
        help="Name of the gene-symbol column in input CSVs (default: ext_gene).",
    )
    parser.add_argument(
        "--batch-size",
        "-b",
        type=int,
        default=50,
        help="PubMed/full-text batch size (default: 50).",
    )
    parser.add_argument(
        "--skip-fbgnid-conversion",
        action="store_true",
        help="Trust the FBgn IDs in the input CSV without rebuilding them.",
    )


def _add_validation_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--test-log",
        action="store_true",
        help="Write per-call GPT logs under data/logs/gpt_queries/.",
    )
    parser.add_argument(
        "--max-gpt-calls-per-stock",
        type=int,
        default=5,
        help="Cap actual GPT calls per stock during validation (default: 5).",
    )


def _add_similarity_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--similarity",
        choices=SIMILARITY_CHOICES,
        default="none",
        help=(
            "Phenotype similarity sidecar layout. "
            "'none' writes only the Stock Phenotype Sheet. "
            "'tiers' adds the cosine-threshold tier workbook. "
            "'simple-buckets' adds the collection/UAS/sleep-circ/balancer workbook. "
            "'keyword-buckets' adds the Keyword Hits / No Keyword Hits workbook. "
            "(default: none)"
        ),
    )


###############################################################################
# Parser construction
###############################################################################


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fl-ai-reagent-stocker",
        description="Drosophila stock finding, splitting, phenotype-sheet, and validation pipelines",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m fl_ai_reagent_stocker find-stocks ./gene_lists --config ./my_config.json
    python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_config.json
    python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_config.json --similarity tiers
    python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_config.json --test-log
    python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_config.json --mode normal
    python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_config.json --mode phenotype --similarity tiers
""",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # find-stocks
    find_parser = subparsers.add_parser(
        "find-stocks",
        help="Stage 1: map genes to FlyBase stocks, references, and the Gene Reagent Index",
    )
    find_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing CSV gene lists.",
    )
    _add_config_arg(find_parser)
    _add_stage1_args(find_parser)

    # split-stocks (organization-only)
    split_parser = subparsers.add_parser(
        "split-stocks",
        help="Stage 2: organize stocks using JSON filters only (no phenotype sheet)",
    )
    split_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing Stage 1 Excel files (e.g., ./gene_lists/Stocks).",
    )
    _add_config_arg(split_parser)
    _add_quiet_arg(split_parser)

    # phenotype-sheet (gene-first phenotype reagent flow)
    phenotype_parser = subparsers.add_parser(
        "phenotype-sheet",
        help="Build the gene-first Stock Phenotype Sheet from Stage 1 workbook(s)",
        description=(
            "Build the gene-first Stock Phenotype Sheet from Stage 1 workbook(s). "
            "Consumes a Stage 1 workbook directory such as ./gene_lists/Stocks "
            "and expects a 'Gene Reagent Index' sheet for full-scope output. "
            "Workbooks without the index still run via a legacy fallback that is "
            "best-effort and not true full scope; re-run find-stocks to refresh "
            "the index."
        ),
    )
    phenotype_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing Stage 1 Excel files (e.g., ./gene_lists/Stocks).",
    )
    _add_config_arg(phenotype_parser)
    _add_quiet_arg(phenotype_parser)
    _add_similarity_arg(phenotype_parser)

    # validate-stocks (GPT validation only)
    validate_parser = subparsers.add_parser(
        "validate-stocks",
        help="Stage 3: append GPT validation results to Ref++ output sheets (no phenotype sheet)",
    )
    validate_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing Stage 1 Excel files (e.g., ./gene_lists/Stocks).",
    )
    _add_config_arg(validate_parser)
    _add_quiet_arg(validate_parser)
    _add_validation_args(validate_parser)

    # run-full-pipeline (end-to-end wrapper)
    full_parser = subparsers.add_parser(
        "run-full-pipeline",
        help="End-to-end wrapper from raw gene-list CSVs",
        description=(
            "End-to-end wrapper from raw gene-list CSVs. "
            "--mode normal runs Stage 1 + split-stocks + validate-stocks. "
            "--mode phenotype runs Stage 1 + phenotype-sheet only and never "
            "invokes GPT validation."
        ),
    )
    full_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing CSV gene lists.",
    )
    _add_config_arg(full_parser)
    _add_quiet_arg(full_parser)
    _add_stage1_args(full_parser)
    full_parser.add_argument(
        "--mode",
        choices=PIPELINE_MODE_CHOICES,
        default="normal",
        help=(
            "Pipeline mode. 'normal' runs Stage 1 + split-stocks + validate-stocks. "
            "'phenotype' runs Stage 1 + phenotype-sheet only. (default: normal)"
        ),
    )
    _add_similarity_arg(full_parser)
    _add_validation_args(full_parser)

    return parser


###############################################################################
# Settings helpers
###############################################################################


def _similarity_to_settings_kwargs(similarity: str) -> dict:
    """Map ``--similarity`` choice to ``Settings`` feature flags."""
    if similarity == "none":
        return {
            "soft_run": True,
            "enable_oai_embedding": False,
            "simple_buckets": False,
            "keyword_bucketing": False,
        }
    if similarity == "tiers":
        return {
            "soft_run": True,
            "enable_oai_embedding": True,
            "simple_buckets": False,
            "keyword_bucketing": False,
        }
    if similarity == "simple-buckets":
        return {
            "soft_run": True,
            "enable_oai_embedding": True,
            "simple_buckets": True,
            "keyword_bucketing": False,
        }
    if similarity == "keyword-buckets":
        return {
            "soft_run": True,
            "enable_oai_embedding": True,
            "simple_buckets": False,
            "keyword_bucketing": True,
        }
    raise ValueError(f"Unknown --similarity choice: {similarity!r}")


def _load_keywords(config_path: Path) -> list[str]:
    config = load_split_config(config_path)
    return config.get("settings", {}).get("relevantSearchTerms", [])


###############################################################################
# Reagent-index pre-flight warning
###############################################################################


def _iter_stage1_workbooks(input_dir: Path) -> Iterable[Path]:
    """Yield Stage 1 Excel workbooks under ``input_dir``.

    Mirrors the discovery logic used by ``StockSplittingPipeline.run`` so the
    warning surfaces for the same files Stage 2 will read.
    """
    if not input_dir.exists() or not input_dir.is_dir():
        return

    def _candidates(directory: Path) -> List[Path]:
        return [
            f
            for f in directory.glob("*.xlsx")
            if "Organized Stocks" not in str(f)
            and "Organized Stock Sheets" not in str(f)
            and "Uncategorized" not in str(f)
            and not f.name.startswith("~$")
            and not f.name.startswith(".")
        ]

    stocks_dir = input_dir / "Stocks"
    direct = _candidates(input_dir)
    nested = _candidates(stocks_dir) if stocks_dir.exists() else []
    yield from (nested if nested else direct)


def _workbook_has_reagent_index(path: Path) -> bool:
    """Return True iff ``path`` has a ``Gene Reagent Index`` sheet."""
    try:
        excel = pd.ExcelFile(path)
    except Exception:
        return False
    return "Gene Reagent Index" in excel.sheet_names


def _warn_if_missing_reagent_index(input_dir: Path) -> None:
    """Emit a stderr warning for Stage 1 workbooks missing the index.

    The warning bubbles even when ``--quiet`` is set so users running
    ``phenotype-sheet`` against an old workbook are told to re-run
    ``find-stocks`` to regenerate the ``Gene Reagent Index`` sheet.
    """
    missing: List[Path] = []
    for workbook in _iter_stage1_workbooks(input_dir):
        if not _workbook_has_reagent_index(workbook):
            missing.append(workbook)

    if not missing:
        return

    print(
        "WARNING: 'Gene Reagent Index' sheet not found in the following Stage 1 "
        "workbook(s); phenotype-sheet output will use the legacy stocks_df "
        "fallback (best-effort, not true full scope).",
        file=sys.stderr,
    )
    for workbook in missing:
        print(f"  - {workbook}", file=sys.stderr)
    print(
        "Re-run `python -m fl_ai_reagent_stocker find-stocks ...` to regenerate "
        "the Gene Reagent Index sheet for full-scope output.",
        file=sys.stderr,
    )


###############################################################################
# Subcommand entry points
###############################################################################


def run_find_stocks(args) -> int:
    settings = Settings(
        batch_size=args.batch_size,
        skip_fbgnid_conversion=args.skip_fbgnid_conversion,
    )
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
    settings = Settings()
    if args.config:
        settings.split_config_path = args.config
    pipeline = StockSplittingPipeline(settings)
    output_dir = pipeline.run(
        input_dir=args.input_dir,
        config_path=args.config,
        verbose=not args.quiet,
    )
    return 0 if output_dir else 1


def run_phenotype_sheet(args) -> int:
    _warn_if_missing_reagent_index(args.input_dir)

    settings = Settings(**_similarity_to_settings_kwargs(args.similarity))
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


def _run_full_pipeline_normal(args, *, refs_output_path: Path) -> int:
    split_input_dir = Path(refs_output_path).parent

    print(f"\n{'-' * 70}")
    print("Step 2/3: Running split-stocks")
    print(f"{'-' * 70}")
    split_settings = Settings()
    if args.config:
        split_settings.split_config_path = args.config
    split_pipeline = StockSplittingPipeline(split_settings)
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
    validate_settings = Settings(
        enable_gpt_logging=args.test_log,
        max_gpt_calls_per_stock=args.max_gpt_calls_per_stock,
    )
    if args.config:
        validate_settings.split_config_path = args.config
    validate_pipeline = StockSplittingPipeline(validate_settings)
    validated_output_dir = validate_pipeline.run(
        input_dir=split_input_dir,
        config_path=args.config,
        verbose=not args.quiet,
        run_validation=True,
    )
    if not validated_output_dir:
        return 1

    print("\nSuccess! Full pipeline completed (mode: normal).")
    print(f"  Stage 1 output: {refs_output_path}")
    print(f"  Stage 2 output directory: {split_output_dir}")
    print(f"  Stage 3 output directory: {validated_output_dir}")
    return 0


def _run_full_pipeline_phenotype(args, *, refs_output_path: Path) -> int:
    split_input_dir = Path(refs_output_path).parent

    print(f"\n{'-' * 70}")
    print(f"Step 2/2: Running phenotype-sheet (--similarity {args.similarity})")
    print(f"{'-' * 70}")
    # Stage 1 just produced the workbook; the index is guaranteed present, but
    # keep the same safety net other phenotype-sheet runs use.
    _warn_if_missing_reagent_index(split_input_dir)

    pheno_settings = Settings(**_similarity_to_settings_kwargs(args.similarity))
    if args.config:
        pheno_settings.split_config_path = args.config
    pheno_pipeline = StockSplittingPipeline(pheno_settings)
    pheno_output_dir = pheno_pipeline.run(
        input_dir=split_input_dir,
        config_path=args.config,
        verbose=not args.quiet,
    )
    if not pheno_output_dir:
        return 1

    print("\nSuccess! Full pipeline completed (mode: phenotype).")
    print(f"  Stage 1 output: {refs_output_path}")
    print(f"  phenotype-sheet output directory: {pheno_output_dir}")
    return 0


def run_full_pipeline(args) -> int:
    refs_settings = Settings(
        batch_size=args.batch_size,
        skip_fbgnid_conversion=args.skip_fbgnid_conversion,
    )

    config_path = args.config if args.config else refs_settings.split_config_path
    keywords = _load_keywords(config_path)

    refs_pipeline = StockFindingPipeline(refs_settings)

    print(f"\n{'=' * 70}")
    print(f"FL.AI REAGENT STOCKER: FULL PIPELINE (mode: {args.mode})")
    print(f"{'=' * 70}")

    print(f"\n{'-' * 70}")
    if args.mode == "normal":
        print("Step 1/3: Running find-stocks")
    else:
        print("Step 1/2: Running find-stocks")
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

    if args.mode == "normal":
        return _run_full_pipeline_normal(args, refs_output_path=refs_output_path)
    if args.mode == "phenotype":
        return _run_full_pipeline_phenotype(args, refs_output_path=refs_output_path)

    print(f"Unknown --mode: {args.mode!r}", file=sys.stderr)
    return 1


###############################################################################
# Dispatcher
###############################################################################


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
    if args.command == "phenotype-sheet":
        return run_phenotype_sheet(args)
    if args.command == "validate-stocks":
        return run_validate_stocks(args)
    if args.command == "run-full-pipeline":
        return run_full_pipeline(args)

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
