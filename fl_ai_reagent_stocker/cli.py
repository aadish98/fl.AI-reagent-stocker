"""
Command-line interface for the ``fl_ai_reagent_stocker`` package.

Primary command:

- ``run``: end-to-end pipeline from raw gene-list CSVs. It builds the Stage 1
  stock/reference workbook, organizes stocks with the JSON config (including
  any phenotype-relevance filters), and validates ``Ref++`` output stocks.
  The JSON config controls whether phenotype embeddings are computed.

Advanced per-stage entry points (``find-stocks``, ``split-stocks``,
``phenotype-sheet``, ``validate-stocks``) remain available for power users and
tests but are not the primary interface.
"""

from __future__ import annotations

import argparse
import gc
import shutil
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import pandas as pd

from .config import Settings
from .pipelines.stock_finding import StockFindingPipeline
from .pipelines.stock_splitting import StockSplittingPipeline, load_split_config


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


###############################################################################
# Parser construction
###############################################################################


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fl-ai-reagent-stocker",
        description="Drosophila reagent-stocker pipeline: gene lists to organized, validated stock sheets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m fl_ai_reagent_stocker run ./gene_lists --config ./my_config.json
""",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # run (primary, config-driven, end-to-end)
    run_parser = subparsers.add_parser(
        "run",
        help="Run the end-to-end pipeline from gene-list CSVs (config-driven)",
        description=(
            "End-to-end pipeline from gene-list CSVs. Builds the Stage 1 "
            "stock/reference workbook, organizes stocks with the JSON config "
            "(including any phenotype-relevance filters), validates Ref++ "
            "output stocks when enabled by the JSON config, and computes "
            "phenotype-embedding cosine columns when enabled by the JSON config. "
            "Outputs are written cleanly: each gene set's final "
            "<gene set>_aggregated.xlsx lives directly in its Stocks folder, "
            "and a combined combination-counts summary is written at the input "
            "root. The similarity sidecar workbook and t-SNE/plot folders are "
            "not emitted by the primary run."
        ),
    )
    run_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing CSV gene lists.",
    )
    _add_config_arg(run_parser)
    _add_quiet_arg(run_parser)

    # ----- Advanced per-stage entry points (not the primary interface) -----

    # find-stocks
    find_parser = subparsers.add_parser(
        "find-stocks",
        help="Advanced: Stage 1 mapping of genes to FlyBase stocks, references, and the Gene Reagent Index",
    )
    find_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing CSV gene lists.",
    )
    _add_config_arg(find_parser)

    # split-stocks (organization-only)
    split_parser = subparsers.add_parser(
        "split-stocks",
        help="Advanced: organize stocks using JSON filters only",
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
        help="Advanced: build the gene-first All Phenotypic Stocks Sheet from Stage 1 workbook(s)",
        description=(
            "Build the gene-first All Phenotypic Stocks Sheet from Stage 1 workbook(s). "
            "Consumes a Stage 1 workbook directory such as ./gene_lists/Stocks "
            "and expects a 'Gene Reagent Index' sheet for full-scope output. "
            "OpenAI phenotype-embedding analysis is controlled by the JSON config."
        ),
    )
    phenotype_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing Stage 1 Excel files (e.g., ./gene_lists/Stocks).",
    )
    _add_config_arg(phenotype_parser)
    _add_quiet_arg(phenotype_parser)

    # validate-stocks (GPT validation only)
    validate_parser = subparsers.add_parser(
        "validate-stocks",
        help="Advanced: append GPT validation results to Ref++ output sheets",
    )
    validate_parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing Stage 1 Excel files (e.g., ./gene_lists/Stocks).",
    )
    _add_config_arg(validate_parser)
    _add_quiet_arg(validate_parser)

    return parser


###############################################################################
# Settings helpers
###############################################################################


def _resolve_config_path(config_path: Optional[Path]) -> Path:
    """Return the explicit config path, or the Settings default config path."""
    if config_path:
        return config_path
    return Settings().split_config_path


def _load_config_settings(config_path: Path) -> dict:
    """Load validated JSON settings."""
    return load_split_config(config_path)["settings"]


def _load_input_policy(config_path: Path) -> dict:
    return _load_config_settings(config_path)["input"]


def _load_pubmed_policy(config_path: Path) -> dict:
    return _load_config_settings(config_path)["pubmed"]


def _load_embeddings_policy(config_path: Path) -> dict:
    return _load_config_settings(config_path)["embeddings"]


def _load_output_policy(config_path: Path) -> dict:
    return _load_config_settings(config_path)["output"]


def _resolve_phenotype_settings_kwargs(embeddings: bool) -> dict:
    """Map the phenotype-embedding toggle to ``Settings`` feature flags.

    The phenotype-sheet flow always builds the gene-first Stock Phenotype
    Sheet. Enabling embeddings additionally computes the OpenAI cosine analysis
    and writes the similarity-tier sidecar workbook and plots.
    """
    enable_embeddings = bool(embeddings)
    return {
        "soft_run": True,
        "enable_oai_embedding": enable_embeddings,
        "phenotype_similarity_sidecar": enable_embeddings,
    }


def _build_phenotype_settings(config_path: Path) -> Settings:
    embeddings_enabled = _load_embeddings_policy(config_path)["enabled"]
    pubmed_policy = _load_pubmed_policy(config_path)
    settings = Settings(
        **_resolve_phenotype_settings_kwargs(embeddings_enabled),
        batch_size=pubmed_policy["batchSize"],
    )
    settings.split_config_path = config_path
    return settings


def _load_keywords(config_path: Path) -> list[str]:
    config = load_split_config(config_path)
    return config.get("settings", {}).get("relevantSearchTerms", [])


def _load_validation_policy(config_path: Path) -> dict:
    """Return config-driven validation settings."""
    config = load_split_config(config_path)
    return config["settings"]["validation"]


def _build_validation_settings(config_path: Path, **overrides) -> Settings:
    """Build Settings populated from the config validation and PubMed policy."""
    validation_policy = _load_validation_policy(config_path)
    pubmed_policy = _load_pubmed_policy(config_path)
    settings = Settings(
        batch_size=pubmed_policy["batchSize"],
        enable_gpt_logging=validation_policy["enableGptLogging"],
        max_gpt_calls_per_stock=validation_policy["maxGptCallsPerStock"],
        min_fulltext_chars=validation_policy["minFullTextChars"],
        gpt_call_delay_seconds=validation_policy["gptCallDelaySeconds"],
        short_circuit_on_functional_validation=(
            validation_policy["shortCircuitOnFunctionalValidation"]
        ),
        **overrides,
    )
    settings.split_config_path = config_path
    return settings


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
            and not f.stem.endswith("_aggregated")  # generated outputs
            and "_similarity" not in f.name          # similarity sidecars
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
    config_path = _resolve_config_path(args.config)
    input_policy = _load_input_policy(config_path)
    pubmed_policy = _load_pubmed_policy(config_path)
    settings = Settings(
        batch_size=pubmed_policy["batchSize"],
        skip_fbgnid_conversion=input_policy["skipFbgnidConversion"],
    )
    settings.split_config_path = config_path
    pipeline = StockFindingPipeline(settings)
    keywords = _load_keywords(config_path)

    output_path = pipeline.run(
        input_dir=args.input_dir,
        keywords=keywords,
        gene_col=input_policy["geneCol"],
        input_gene_col=input_policy["inputGeneCol"],
        skip_fbgnid_conversion=input_policy["skipFbgnidConversion"],
        run_functional_validation=False,
    )
    return 0 if output_path else 1


def run_split_stocks(args) -> int:
    config_path = _resolve_config_path(args.config)
    pubmed_policy = _load_pubmed_policy(config_path)
    settings = Settings(batch_size=pubmed_policy["batchSize"])
    settings.split_config_path = config_path
    pipeline = StockSplittingPipeline(settings)
    output_dir = pipeline.run(
        input_dir=args.input_dir,
        config_path=config_path,
        verbose=not args.quiet,
    )
    return 0 if output_dir else 1


def run_phenotype_sheet(args) -> int:
    _warn_if_missing_reagent_index(args.input_dir)

    config_path = _resolve_config_path(args.config)
    settings = _build_phenotype_settings(config_path)
    pipeline = StockSplittingPipeline(settings)
    output_dir = pipeline.run(
        input_dir=args.input_dir,
        config_path=config_path,
        verbose=not args.quiet,
    )
    return 0 if output_dir else 1


def run_validate_stocks(args) -> int:
    config_path = _resolve_config_path(args.config)
    validation_policy = _load_validation_policy(config_path)
    settings = _build_validation_settings(config_path)
    if not _should_run_validation(
        settings.openai_api_key,
        validation_enabled=validation_policy["enabled"],
    ):
        if not validation_policy["enabled"]:
            print("Config settings.validation.enabled is false; skipping validation.")
        else:
            print("No OpenAI API key configured; skipping validation.")
        return 0
    pipeline = StockSplittingPipeline(settings)
    output_dir = pipeline.run(
        input_dir=args.input_dir,
        config_path=config_path,
        verbose=not args.quiet,
        run_validation=True,
    )
    return 0 if output_dir else 1


def _run_embedding_analysis(
    *,
    split_input_dir: Path,
    config_path: Path,
    quiet: bool,
) -> int:
    """Run the phenotype-sheet flow with OpenAI embeddings for analysis output."""
    print(f"\n{'-' * 70}")
    print("Step 4/4: Running phenotype-embedding analysis")
    print(f"{'-' * 70}")
    _warn_if_missing_reagent_index(split_input_dir)

    pubmed_policy = _load_pubmed_policy(config_path)
    pheno_settings = Settings(**_resolve_phenotype_settings_kwargs(embeddings=True))
    pheno_settings.batch_size = pubmed_policy["batchSize"]
    pheno_settings.split_config_path = config_path
    # Keep cosine columns in the aggregated workbook, but do not emit the
    # similarity sidecar workbook, the t-SNE/plot folders, or .txt reports.
    pheno_settings.phenotype_similarity_sidecar = False
    pheno_settings.phenotype_similarity_plots = False
    pheno_settings.emit_no_pmid_report = False
    # Write the aggregated workbook directly into the gene set's Stocks folder.
    pheno_settings.organized_output_dir = split_input_dir
    pheno_pipeline = StockSplittingPipeline(pheno_settings)
    pheno_output_dir = pheno_pipeline.run(
        input_dir=split_input_dir,
        config_path=config_path,
        verbose=not quiet,
    )
    return 0 if pheno_output_dir else 1


def _run_pipeline_for_gene_set(
    args,
    *,
    csv_path: Path,
    run_dir: Path,
    config_path: Path,
    keywords: Optional[List[str]],
    input_policy: dict,
    pubmed_policy: dict,
    embeddings_enabled: bool,
    output_policy: dict,
    refs_pipeline: StockFindingPipeline,
    validation_policy: dict,
) -> Optional[Path]:
    """Run Stage 1 -> organize -> validate (+ embeddings) for a single CSV.

    The gene set is processed inside its own isolated ``run_dir`` so that
    stock matching, gene-count accounting, and the organized/embedding outputs
    are scoped to exactly one input CSV. Returns the organized output
    directory, or ``None`` if any stage failed for this gene set.
    """
    stem = csv_path.stem

    # Stage the single CSV in a clean, isolated working directory so the
    # directory-level globbing in Stage 1/2 only ever sees this one gene set.
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(csv_path, run_dir / csv_path.name)

    print(f"\n{'-' * 70}")
    print(f"[{stem}] Step 1/3: Running find-stocks")
    print(f"{'-' * 70}")
    refs_output_path = refs_pipeline.run(
        input_dir=run_dir,
        keywords=keywords,
        gene_col=input_policy["geneCol"],
        input_gene_col=input_policy["inputGeneCol"],
        skip_fbgnid_conversion=input_policy["skipFbgnidConversion"],
        run_functional_validation=False,
        output_name=f"{stem}.xlsx",
    )
    if not refs_output_path:
        return None

    split_input_dir = Path(refs_output_path).parent

    print(f"\n{'-' * 70}")
    print(f"[{stem}] Step 2/3: Organizing stocks")
    print(f"{'-' * 70}")
    split_settings = Settings(
        batch_size=pubmed_policy["batchSize"],
        emit_no_pmid_report=False,
        organized_output_dir=split_input_dir,
    )
    split_settings.split_config_path = config_path
    split_output_dir = StockSplittingPipeline(split_settings).run(
        input_dir=split_input_dir,
        config_path=config_path,
        verbose=not args.quiet,
    )
    if not split_output_dir:
        return None

    print(f"\n{'-' * 70}")
    print(f"[{stem}] Step 3/3: Validating Ref++ stocks")
    print(f"{'-' * 70}")
    validate_settings = _build_validation_settings(
        config_path,
        emit_no_pmid_report=False,
        organized_output_dir=split_input_dir,
    )
    if _should_run_validation(
        validate_settings.openai_api_key,
        validation_enabled=validation_policy["enabled"],
    ):
        validated_output_dir = StockSplittingPipeline(validate_settings).run(
            input_dir=split_input_dir,
            config_path=config_path,
            verbose=not args.quiet,
            run_validation=True,
        )
        if not validated_output_dir:
            return None
    elif not validation_policy["enabled"]:
        print("  Config settings.validation.enabled is false; skipping Ref++ validation pass.")
    else:
        print("  No OpenAI API key configured; skipping Ref++ validation pass.")

    if embeddings_enabled:
        if _run_embedding_analysis(
            split_input_dir=split_input_dir,
            config_path=config_path,
            quiet=args.quiet,
        ) != 0:
            return None

    # Enforce a clean output layout: final workbook in Stocks, no unsplit
    # workbook (unless config output policy preserves it), no .txt or _similarity artifacts, no
    # Organized Stocks wrapper directory.
    return _finalize_gene_set_run_outputs(
        run_dir=run_dir,
        stem=stem,
        log_mode=output_policy["preserveUnsplitWorkbook"],
    )


def _should_run_validation(
    openai_api_key: Optional[str],
    *,
    validation_enabled: bool = True,
) -> bool:
    """Return whether the GPT validation pass can run with current settings."""
    return validation_enabled and bool(str(openai_api_key or "").strip())


def _finalize_gene_set_run_outputs(
    run_dir: Path,
    stem: str,
    *,
    log_mode: bool = False,
) -> Optional[Path]:
    """Enforce a clean output layout for a single gene-set run.

    Guarantees the invariant the pipeline is designed to produce (the artifacts
    are also prevented at the source):

    - the final ``<stem>_aggregated.xlsx`` lives directly in ``Stocks``;
    - no ``Organized Stocks`` wrapper directory remains;
    - no ``.txt`` sidecar reports remain;
    - no ``*_similarity`` files/folders remain;
    - the unsplit Stage 1 ``<stem>.xlsx`` is removed unless ``log_mode``.

    Returns the gene set's ``Stocks`` directory (where the final workbook
    lives), or ``None`` if the expected final workbook cannot be located.
    """
    stocks_dir = run_dir / "Stocks"
    final_name = f"{stem}_aggregated.xlsx"
    final_path = stocks_dir / final_name

    # If, for any reason, the aggregated workbook still sits under a nested
    # "Organized Stocks" directory, move it up into Stocks.
    if not final_path.exists():
        for nested in run_dir.rglob(final_name):
            if "Organized Stocks" in nested.parts and nested.is_file():
                stocks_dir.mkdir(parents=True, exist_ok=True)
                shutil.move(str(nested), str(final_path))
                break

    if not final_path.exists():
        print(
            f"ERROR: expected final workbook not found: {final_path}",
            file=sys.stderr,
        )
        return None

    # Remove the unsplit Stage 1 workbook unless preserving for audit.
    unsplit_path = stocks_dir / f"{stem}.xlsx"
    if unsplit_path.exists() and not log_mode:
        try:
            unsplit_path.unlink()
        except OSError:
            pass
    elif unsplit_path.exists() and log_mode:
        print(f"  [log-mode] preserved unsplit workbook: {unsplit_path}")

    # Remove any leftover .txt sidecars and *_similarity artifacts, and any
    # Organized Stocks wrapper directories under this run.
    for txt_file in run_dir.rglob("*.txt"):
        try:
            txt_file.unlink()
        except OSError:
            pass
    for similarity_path in sorted(
        run_dir.rglob("*_similarity*"), key=lambda p: len(p.parts), reverse=True
    ):
        try:
            if similarity_path.is_dir():
                shutil.rmtree(similarity_path, ignore_errors=True)
            else:
                similarity_path.unlink()
        except OSError:
            pass
    for organized_dir in run_dir.rglob("Organized Stocks"):
        if organized_dir.is_dir():
            shutil.rmtree(organized_dir, ignore_errors=True)

    return stocks_dir


# Directory names produced by the pipeline itself. Discovery skips anything that
# lives under one of these so recursive runs never reprocess prior outputs or
# the per-run staged input copies.
_GENERATED_DIR_NAMES = {
    "Per Gene Set Runs",
    "Stocks",
    "Organized Stocks",
    "Organized Stock Sheets",
    "Uncategorized",
}
_GENERATED_CSV_NAMES = {
    "combination_counts_summary.csv",
}


def _discover_input_csvs(input_dir: Path) -> List[Path]:
    """Recursively find gene-set CSVs under ``input_dir``.

    Generated/output trees (``Per Gene Set Runs``, ``Stocks``, ``Organized
    Stocks``, ``*_similarity`` plot folders, etc.) and hidden files are excluded
    so re-runs never reprocess prior outputs or staged copies.
    """
    discovered: List[Path] = []
    for path in sorted(input_dir.rglob("*.csv")):
        if path.name.startswith(".") or path.name in _GENERATED_CSV_NAMES:
            continue
        rel_parts = path.relative_to(input_dir).parts[:-1]
        if any(part in _GENERATED_DIR_NAMES for part in rel_parts):
            continue
        if any(part.endswith("_similarity") for part in rel_parts):
            continue
        discovered.append(path)
    return discovered


def _validate_gene_columns(
    csv_files: List[Path], gene_col: str, input_gene_col: str
) -> List[str]:
    """Return error messages for any CSV missing the required gene columns.

    A single config-provided gene column pair applies to the whole
    recursive run, so any discovered CSV that lacks either column is a fatal
    error (the caller aborts before heavy work).
    """
    required = [gene_col, input_gene_col]
    errors: List[str] = []
    for path in csv_files:
        try:
            header = pd.read_csv(path, nrows=0)
        except Exception as exc:  # noqa: BLE001 - report unreadable headers
            errors.append(f"{path}: could not read header ({exc})")
            continue
        missing = [col for col in required if col not in header.columns]
        if missing:
            errors.append(f"{path}: missing required column(s): {', '.join(missing)}")
    return errors


def run_pipeline(args) -> int:
    """Primary end-to-end command: Stage 1 -> organize -> validate (+ embeddings).

    Each gene-set CSV in ``input_dir`` is processed independently and produces
    its own organized workbook (named after the CSV), rather than being merged
    into a single aggregated workbook. Behavior is config-driven: the JSON
    config determines the organized output sheets (including any
    phenotype-relevance filters) and whether to run phenotype embeddings.
    """
    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        print(f"ERROR: input path is not a directory: {input_dir}", file=sys.stderr)
        return 1

    csv_files = _discover_input_csvs(input_dir)
    if not csv_files:
        print(f"ERROR: no gene-list CSV files found in {input_dir}", file=sys.stderr)
        return 1

    config_path = _resolve_config_path(args.config)
    input_policy = _load_input_policy(config_path)
    pubmed_policy = _load_pubmed_policy(config_path)
    embeddings_policy = _load_embeddings_policy(config_path)
    output_policy = _load_output_policy(config_path)
    validation_policy = _load_validation_policy(config_path)

    column_errors = _validate_gene_columns(
        csv_files, input_policy["geneCol"], input_policy["inputGeneCol"]
    )
    if column_errors:
        print(
            "ERROR: the following gene-set CSV(s) are missing required columns "
            f"(settings.input.geneCol '{input_policy['geneCol']}', "
            f"settings.input.inputGeneCol '{input_policy['inputGeneCol']}'):",
            file=sys.stderr,
        )
        for message in column_errors:
            print(f"  - {message}", file=sys.stderr)
        return 1

    refs_settings = Settings(
        batch_size=pubmed_policy["batchSize"],
        skip_fbgnid_conversion=input_policy["skipFbgnidConversion"],
        emit_no_pmid_report=False,
    )
    refs_settings.split_config_path = config_path
    keywords = _load_keywords(config_path)

    print(f"\n{'=' * 70}")
    print("FL.AI REAGENT STOCKER: RUN")
    print(f"{'=' * 70}")
    print(f"  Input directory: {input_dir}")
    print(f"  Gene-set CSVs: {len(csv_files)} (one organized workbook each)")

    # Reuse a single Stage 1 pipeline so the heavy FlyBase tables load once and
    # are cached across gene sets.
    refs_pipeline = StockFindingPipeline(refs_settings)
    runs_root = input_dir / "Per Gene Set Runs"
    runs_root.mkdir(exist_ok=True)

    results: List[Tuple[str, Optional[Path]]] = []
    for csv_path in csv_files:
        rel_key = csv_path.relative_to(input_dir).with_suffix("")
        print(f"\n{'#' * 70}")
        print(f"# Gene set: {rel_key}")
        print(f"{'#' * 70}")
        organized_dir = _run_pipeline_for_gene_set(
            args,
            csv_path=csv_path,
            run_dir=runs_root / rel_key,
            config_path=config_path,
            keywords=keywords,
            input_policy=input_policy,
            pubmed_policy=pubmed_policy,
            embeddings_enabled=embeddings_policy["enabled"],
            output_policy=output_policy,
            refs_pipeline=refs_pipeline,
            validation_policy=validation_policy,
        )
        results.append((str(csv_path.relative_to(input_dir)), organized_dir))
        del organized_dir, rel_key
        gc.collect()

    succeeded = [(name, path) for name, path in results if path is not None]
    failed = [name for name, path in results if path is None]

    print(f"\n{'=' * 70}")
    print("RUN SUMMARY")
    print(f"{'=' * 70}")
    for name, path in results:
        stem = Path(name).stem
        if path is None:
            print(f"  FAILED  {name}")
        else:
            workbook = path / f"{stem}_aggregated.xlsx"
            print(f"  OK      {name} -> {workbook}")
    print(f"\n  {len(succeeded)}/{len(results)} gene set(s) completed successfully.")

    # Automatically generate the combined combination-counts summary at the
    # input root so the user does not need to run the script separately.
    if succeeded:
        print(f"\n{'-' * 70}")
        print("Generating combination-counts summary")
        print(f"{'-' * 70}")
        try:
            _generate_combination_summary(input_dir, config_path=config_path)
        except Exception as exc:  # noqa: BLE001 - summary is best-effort
            print(f"WARNING: failed to generate summary: {exc}", file=sys.stderr)

    if not succeeded:
        return 1
    return 1 if failed else 0


def _generate_combination_summary(
    input_dir: Path, *, config_path: Optional[Path]
) -> None:
    """Load the summarizer script module and write the combined summary."""
    import importlib.util

    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "summarize_combination_counts.py"
    )
    if not script_path.exists():
        print(
            f"WARNING: summary script not found at {script_path}; skipping.",
            file=sys.stderr,
        )
        return

    spec = importlib.util.spec_from_file_location(
        "summarize_combination_counts", script_path
    )
    if spec is None or spec.loader is None:
        print("WARNING: could not load summary script; skipping.", file=sys.stderr)
        return
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.summarize_combination_counts(input_dir, config_path=config_path)


###############################################################################
# Dispatcher
###############################################################################


def main(argv: Optional[list] = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "run":
        return run_pipeline(args)
    if args.command == "find-stocks":
        return run_find_stocks(args)
    if args.command == "split-stocks":
        return run_split_stocks(args)
    if args.command == "phenotype-sheet":
        return run_phenotype_sheet(args)
    if args.command == "validate-stocks":
        return run_validate_stocks(args)

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
