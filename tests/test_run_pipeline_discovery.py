from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from openpyxl import Workbook


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker import cli  # noqa: E402


def _write_csv(path: Path, header: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + "\nval\n", encoding="utf-8")


def _minimal_validation_policy(**overrides) -> dict:
    policy = {
        "enabled": False,
        "maxGptCallsPerStock": 5,
        "minFullTextChars": 500,
        "gptCallDelaySeconds": 0.5,
        "shortCircuitOnFunctionalValidation": True,
        "enableGptLogging": False,
    }
    policy.update(overrides)
    return policy


def _write_minimal_config(path: Path, validation: dict | None) -> None:
    payload = {
        "settings": {
            "input": {
                "geneCol": "flybase_gene_id",
                "inputGeneCol": "ext_gene",
                "skipFbgnidConversion": False,
            },
            "pubmed": {
                "batchSize": 50,
            },
            "embeddings": {
                "enabled": False,
            },
            "output": {
                "preserveUnsplitWorkbook": False,
            },
            "relevantSearchTerms": ["sleep"],
        },
        "filters": {
            "all": {
                "column": "collection",
                "type": "contains",
                "value": "",
            }
        },
        "combinations": [["all"]],
    }
    if validation is not None:
        payload["settings"]["validation"] = validation
    path.write_text(json.dumps(payload), encoding="utf-8")


def _write_summary_fixture_workbook(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "Contents"
    for row in (
        ["Prioritized sheet breakdown"],
        ["Sheet criteria", "# Stocks", "# Alleles", "# Genes", "Sheet Name"],
        ["all", 2, 1, 1, "All"],
        ["Total", 2, 1, 1, ""],
    ):
        worksheet.append(row)
    workbook.save(path)


class TestRecursiveDiscovery(unittest.TestCase):
    def test_discovers_nested_csvs_and_excludes_generated_dirs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cols = "flybase_gene_id,ext_gene"

            # Real gene-set inputs, nested.
            _write_csv(root / "CSW FC0.5" / "a.csv", cols)
            _write_csv(root / "CSW FC1" / "b.csv", cols)
            _write_csv(root / "top.csv", cols)

            # Generated / output trees that must be excluded.
            _write_csv(
                root / "Per Gene Set Runs" / "a" / "a.csv", cols
            )
            _write_csv(
                root / "CSW FC0.5" / "Stocks" / "leftover.csv", cols
            )
            _write_csv(
                root / "x" / "Organized Stocks" / "org.csv", cols
            )
            _write_csv(
                root / "y_aggregated_similarity" / "plot.csv", cols
            )
            # Hidden file.
            _write_csv(root / ".hidden.csv", cols)

            found = {
                p.relative_to(root).as_posix()
                for p in cli._discover_input_csvs(root)
            }
            self.assertEqual(
                found,
                {"CSW FC0.5/a.csv", "CSW FC1/b.csv", "top.csv"},
            )


class TestGeneColumnValidation(unittest.TestCase):
    def test_reports_missing_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            good = root / "good.csv"
            bad = root / "bad.csv"
            _write_csv(good, "flybase_gene_id,ext_gene")
            _write_csv(bad, "FlyBase_ID,Fly_Gene_symbol")

            errors = cli._validate_gene_columns(
                [good, bad], "flybase_gene_id", "ext_gene"
            )
            self.assertEqual(len(errors), 1)
            self.assertIn("bad.csv", errors[0])
            self.assertIn("flybase_gene_id", errors[0])
            self.assertIn("ext_gene", errors[0])

    def test_no_errors_when_all_columns_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            csv_path = root / "ok.csv"
            _write_csv(csv_path, "flybase_gene_id,ext_gene,extra")
            self.assertEqual(
                cli._validate_gene_columns(
                    [csv_path], "flybase_gene_id", "ext_gene"
                ),
                [],
            )


class TestRunPipelineFatalOnMissingColumn(unittest.TestCase):
    def test_run_pipeline_aborts_before_heavy_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_csv(root / "CSW FC1" / "b.csv", "FlyBase_ID,Fly_Gene_symbol")

            args = SimpleNamespace(
                input_dir=root,
                gene_col="flybase_gene_id",
                input_gene_col="ext_gene",
                config=None,
            )
            # Validation fails before any Settings/pipeline construction.
            self.assertEqual(cli.run_pipeline(args), 1)


class TestGenerateCombinationSummary(unittest.TestCase):
    def test_dynamic_import_registers_dataclass_module(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_summary_fixture_workbook(
                root
                / "Per Gene Set Runs"
                / "Folder"
                / "SetA"
                / "Stocks"
                / "SetA_aggregated.xlsx"
            )
            config_path = root / "config.json"
            _write_minimal_config(
                config_path,
                _minimal_validation_policy(maxGptCallsPerStock=1),
            )

            cli._generate_combination_summary(root, config_path=config_path)

            self.assertTrue((root / "combination_counts_summary.xlsx").exists())
            self.assertTrue((root / "combination_counts_summary.csv").exists())


class TestValidationPolicy(unittest.TestCase):
    def test_loads_validation_policy_from_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            config_path = Path(tmp) / "config.json"
            _write_minimal_config(
                config_path,
                _minimal_validation_policy(maxGptCallsPerStock=2),
            )

            validation_policy = cli._load_validation_policy(config_path)
            self.assertEqual(validation_policy["enabled"], False)
            self.assertEqual(validation_policy["maxGptCallsPerStock"], 2)
            self.assertEqual(validation_policy["minFullTextChars"], 500)

    def test_validation_policy_requires_validation_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            config_path = Path(tmp) / "config.json"
            _write_minimal_config(config_path, None)

            with self.assertRaisesRegex(ValueError, "settings.validation is required"):
                cli._load_validation_policy(config_path)

    def test_validation_policy_requires_explicit_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            config_path = Path(tmp) / "config.json"
            _write_minimal_config(config_path, {})

            with self.assertRaisesRegex(ValueError, "settings.validation.enabled is required"):
                cli._load_validation_policy(config_path)

    def test_removed_pipeline_policy_flags_are_rejected(self):
        parser = cli.create_parser()
        removed_flags = (
            "--gene-col",
            "--input-gene-col",
            "--skip-fbgnid-conversion",
            "--batch-size",
            "-b",
            "--embeddings",
            "--no-embeddings",
            "--log-mode",
            "--test-log",
        )
        for flag in removed_flags:
            with self.subTest(flag=flag):
                with self.assertRaises(SystemExit):
                    parser.parse_args(["run", "./gene_lists", flag])


def _build_dirty_run_dir(root: Path, stem: str) -> Path:
    """Create a gene-set run dir with all the artifacts finalization removes."""
    run_dir = root
    stocks = run_dir / "Stocks"
    stocks.mkdir(parents=True, exist_ok=True)
    (stocks / f"{stem}_aggregated.xlsx").write_text("final", encoding="utf-8")
    (stocks / f"{stem}.xlsx").write_text("unsplit", encoding="utf-8")
    (stocks / "references_without_pmid_fbrf.txt").write_text("txt", encoding="utf-8")
    (stocks / f"{stem}_aggregated_similarity_tiers.xlsx").write_text("s", encoding="utf-8")
    plots = stocks / f"{stem}_aggregated_similarity"
    plots.mkdir(parents=True, exist_ok=True)
    (plots / "tsne.png").write_text("png", encoding="utf-8")
    return run_dir


class TestFinalizeGeneSetRunOutputs(unittest.TestCase):
    def test_default_cleans_artifacts_and_removes_unsplit(self):
        with tempfile.TemporaryDirectory() as tmp:
            stem = "FC1_Wake_frequency_2_n=77genes"
            run_dir = _build_dirty_run_dir(Path(tmp), stem)

            result = cli._finalize_gene_set_run_outputs(run_dir, stem, log_mode=False)

            stocks = run_dir / "Stocks"
            self.assertEqual(result, stocks)
            self.assertTrue((stocks / f"{stem}_aggregated.xlsx").exists())
            self.assertFalse((stocks / f"{stem}.xlsx").exists())
            self.assertEqual(list(run_dir.rglob("*.txt")), [])
            self.assertEqual(list(run_dir.rglob("*_similarity*")), [])
            self.assertEqual(list(run_dir.rglob("Organized Stocks")), [])

    def test_log_mode_preserves_unsplit_workbook(self):
        with tempfile.TemporaryDirectory() as tmp:
            stem = "FC1_Wake_frequency_2_n=77genes"
            run_dir = _build_dirty_run_dir(Path(tmp), stem)

            cli._finalize_gene_set_run_outputs(run_dir, stem, log_mode=True)

            stocks = run_dir / "Stocks"
            self.assertTrue((stocks / f"{stem}_aggregated.xlsx").exists())
            self.assertTrue((stocks / f"{stem}.xlsx").exists())
            # Other artifacts are still removed even in log mode.
            self.assertEqual(list(run_dir.rglob("*.txt")), [])
            self.assertEqual(list(run_dir.rglob("*_similarity*")), [])

    def test_moves_aggregated_from_legacy_organized_stocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            stem = "SetA"
            run_dir = Path(tmp)
            organized = run_dir / "Stocks" / "Organized Stocks"
            organized.mkdir(parents=True, exist_ok=True)
            (organized / f"{stem}_aggregated.xlsx").write_text("final", encoding="utf-8")

            result = cli._finalize_gene_set_run_outputs(run_dir, stem, log_mode=False)

            stocks = run_dir / "Stocks"
            self.assertEqual(result, stocks)
            self.assertTrue((stocks / f"{stem}_aggregated.xlsx").exists())
            self.assertEqual(list(run_dir.rglob("Organized Stocks")), [])

    def test_returns_none_when_final_workbook_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "Stocks").mkdir(parents=True, exist_ok=True)
            self.assertIsNone(
                cli._finalize_gene_set_run_outputs(run_dir, "Missing", log_mode=False)
            )


if __name__ == "__main__":
    unittest.main()
