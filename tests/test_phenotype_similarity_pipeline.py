from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker.config import Settings  # noqa: E402
from fl_ai_reagent_stocker.cli import create_parser  # noqa: E402
from fl_ai_reagent_stocker.integrations.phenotype_similarity import (  # noqa: E402
    EmbeddingSimilarityScorer,
)
from fl_ai_reagent_stocker.pipelines.stock_splitting import (  # noqa: E402
    PHENOTYPE_SIMILARITY_EMBEDDING_MODEL,
    SIMILARITY_TIER_SHEET_COUNT,
    StockSplittingPipeline,
    _build_simple_bucket_workbook_entries,
    _build_similarity_tier_sheets,
    _compute_max_cosine_similarity,
    load_split_config,
)


TEST_FIXTURE_DIR = REPO_ROOT / "data" / "gene_sets" / "TEST"
CONFIG_PATH = REPO_ROOT / "data" / "config" / "stock_split_config_no_bloomington.json"


def _deterministic_embedding_map(texts):
    mapping = {}
    for text in texts:
        normalized = str(text or "").strip().lower()
        if not normalized:
            continue
        mapping[text] = [
            float(len(normalized)),
            float(normalized.count("sleep")),
            float(normalized.count("circadian")),
            float(normalized.count("rhythm")),
            float(sum(ord(ch) for ch in normalized) % 101),
        ]
    return mapping


def _fake_ensure_embeddings(self, texts, cache):
    mapping = _deterministic_embedding_map(texts)
    for text, embedding in mapping.items():
        cache.set(self.model, text, embedding)
    cache.save()
    return mapping


class TestPhenotypeSimilarityPipeline(unittest.TestCase):
    def test_cli_parses_simple_buckets_flag(self):
        parser = create_parser()
        args = parser.parse_args(
            [
                "split-stocks",
                "./gene_lists/Stocks",
                "--soft-run",
                "--OAI-embedding",
                "--simple-buckets",
            ]
        )
        self.assertTrue(args.soft_run)
        self.assertTrue(args.oai_embedding)
        self.assertTrue(args.simple_buckets)

    def test_config_requires_explicit_phenotype_similarity_targets(self):
        self.assertTrue(CONFIG_PATH.exists(), f"Missing config: {CONFIG_PATH}")

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_config = Path(tmp_dir) / "missing_targets.json"
            payload = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
            payload["settings"].pop("phenotypeSimilarityTargets", None)
            tmp_config.write_text(json.dumps(payload), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "phenotypeSimilarityTargets"):
                load_split_config(tmp_config)

    def test_soft_run_outputs_similarity_columns_and_plots(self):
        self.assertTrue(TEST_FIXTURE_DIR.exists(), f"Missing test fixture: {TEST_FIXTURE_DIR}")
        self.assertTrue(CONFIG_PATH.exists(), f"Missing config: {CONFIG_PATH}")

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            fixture_copy = tmp_root / "TEST"
            shutil.copytree(TEST_FIXTURE_DIR, fixture_copy)

            settings = Settings(
                soft_run=True,
                enable_oai_embedding=True,
                openai_api_key="test-key",
                phenotype_embedding_cache_path=tmp_root / "cache" / "phenotype_embeddings.csv",
                phenotype_target_embedding_cache_path=tmp_root / "cache" / "target_embeddings.csv",
            )
            pipeline = StockSplittingPipeline(settings)

            with patch.object(
                EmbeddingSimilarityScorer,
                "_ensure_embeddings",
                new=_fake_ensure_embeddings,
            ):
                output_dir = pipeline.run(
                    input_dir=fixture_copy,
                    config_path=CONFIG_PATH,
                    verbose=False,
                )

            self.assertIsNotNone(output_dir)
            workbook_path = Path(output_dir) / "aggregated_stock_refs_aggregated.xlsx"
            self.assertTrue(workbook_path.exists(), workbook_path)
            tier_workbook_path = (
                Path(output_dir) / "aggregated_stock_refs_aggregated_similarity_tiers.xlsx"
            )
            self.assertTrue(tier_workbook_path.exists(), tier_workbook_path)

            phenotype_df = pd.read_excel(workbook_path, sheet_name="Stock Phenotype Sheet")
            expected_columns = [
                "Phenotype",
                "Qualifier",
                "PMID",
                "PMCID",
                "Cosine Similarity (sleep)",
                "Cosine Similarity (circadian)",
            ]
            for column in expected_columns:
                self.assertIn(column, phenotype_df.columns)
            self.assertNotIn("_reference_url", phenotype_df.columns)

            tier_workbook = pd.ExcelFile(tier_workbook_path)
            try:
                expected_tiers = _build_similarity_tier_sheets(phenotype_df)
                self.assertEqual(
                    tier_workbook.sheet_names,
                    ["Contents", "Gene Set", "Stock Phenotype Sheet", *[sheet_name for sheet_name, _, _ in expected_tiers]],
                )
                self.assertLessEqual(len(expected_tiers), SIMILARITY_TIER_SHEET_COUNT)
                self.assertTrue(expected_tiers)
                self.assertTrue(
                    all(not tier_df.empty for _, tier_df, _ in expected_tiers)
                )
                expected_gene_set_df = pd.read_csv(fixture_copy / "vGAT_genes.csv", dtype=str)
                tier_gene_set_df = pd.read_excel(
                    tier_workbook_path,
                    sheet_name="Gene Set",
                )
                pd.testing.assert_frame_equal(
                    tier_gene_set_df,
                    expected_gene_set_df,
                    check_dtype=False,
                )
                tier_first_df = pd.read_excel(
                    tier_workbook_path,
                    sheet_name="Stock Phenotype Sheet",
                )
                self.assertEqual(
                    tier_first_df.columns.tolist(),
                    phenotype_df.columns.tolist(),
                )
                pd.testing.assert_frame_equal(
                    tier_first_df,
                    phenotype_df,
                    check_dtype=False,
                )
                tier_contents_df = pd.read_excel(
                    tier_workbook_path,
                    sheet_name="Contents",
                    header=None,
                )
                contents_text = "\n".join(
                    tier_contents_df.fillna("").astype(str).agg(" ".join, axis=1).tolist()
                )
                self.assertIn("Tier Workbook Contents", contents_text)
                self.assertIn("Max Cosine Similarity", contents_text)
                self.assertIn("empty buckets are skipped", contents_text)
                self.assertIn("Gene Set", contents_text)
                self.assertIn("Stock Phenotype Sheet", contents_text)

                scored_rows = int(_compute_max_cosine_similarity(phenotype_df).notna().sum())
                tier_row_total = 0
                previous_max_score = None
                for sheet_name, expected_df, metadata in expected_tiers:
                    actual_df = pd.read_excel(tier_workbook_path, sheet_name=sheet_name)
                    self.assertEqual(len(actual_df), len(expected_df))
                    pd.testing.assert_frame_equal(
                        actual_df,
                        expected_df,
                        check_dtype=False,
                    )
                    if not actual_df.empty:
                        self.assertIn("Similarity Range", actual_df.columns)
                        self.assertIn("Max Cosine Similarity", actual_df.columns)
                        actual_scores = pd.to_numeric(
                            actual_df["Max Cosine Similarity"], errors="coerce"
                        )
                        self.assertTrue(
                            actual_df["Similarity Range"].eq(metadata["range_label"]).all()
                        )
                        if metadata["lower_bound"] is None:
                            self.assertTrue((actual_scores < metadata["upper_bound"]).all())
                        elif metadata["upper_bound"] >= 1.0:
                            self.assertTrue((actual_scores >= metadata["lower_bound"]).all())
                        else:
                            self.assertTrue(
                                actual_scores.between(
                                    metadata["lower_bound"],
                                    metadata["upper_bound"],
                                    inclusive="left",
                                ).all()
                            )
                        current_max_score = float(actual_scores.max())
                        if previous_max_score is not None:
                            self.assertLessEqual(current_max_score, previous_max_score)
                        previous_max_score = current_max_score
                    tier_row_total += len(actual_df)
                self.assertEqual(tier_row_total, scored_rows)
            finally:
                tier_workbook.close()

            phenotype_cache_df = pd.read_csv(settings.phenotype_embedding_cache_path)
            target_cache_df = pd.read_csv(settings.phenotype_target_embedding_cache_path)
            self.assertEqual(
                set(phenotype_cache_df["model"].astype(str)),
                {PHENOTYPE_SIMILARITY_EMBEDDING_MODEL},
            )
            self.assertEqual(
                set(target_cache_df["model"].astype(str)),
                {PHENOTYPE_SIMILARITY_EMBEDDING_MODEL},
            )

            similarity_dir = Path(output_dir) / "aggregated_stock_refs_aggregated_similarity"
            self.assertTrue(similarity_dir.exists(), similarity_dir)
            expected_artifacts = [
                "sleep_similarity_density.png",
                "circadian_similarity_density.png",
                "similarity_vs_frequency.png",
                "tsne_cosine_sleep.png",
                "tsne_cosine_circadian.png",
            ]
            for filename in expected_artifacts:
                self.assertTrue((similarity_dir / filename).exists(), filename)

    def test_simple_buckets_build_ordered_combinations_without_double_counting(self):
        stock_df = pd.DataFrame(
            {
                "FBst": ["FBst0001", "FBst0002", "FBst0003"],
                "stock_number": ["101", "202", "303"],
                "collection": ["BDSC", "VDRC", "NIG"],
                "UAS": [True, True, False],
                "num_Balancers": [0, 0, 1],
                "relevant_gene_symbols": ["gene-a", "gene-a;gene-b", "gene-c"],
                "AlleleSymbol": ["allele-a", "allele-a;allele-b", "allele-c"],
            }
        )
        phenotype_df = pd.DataFrame(
            {
                "Source/ Stock #": ["BDSC (101)", "VDRC (202)", "NIG (303)"],
                "Phenotype": ["sleep phenotype", "circadian phenotype", "other phenotype"],
                "Reference": ["ref-1", "ref-2", "ref-3"],
            }
        )

        entries = _build_simple_bucket_workbook_entries(
            phenotype_sheet_df=phenotype_df,
            combination_outputs=[(["dummy"], stock_df, {})],
            csv_input_genes={"gene-a", "gene-b", "gene-c"},
        )

        self.assertEqual(len(entries), 24)
        self.assertEqual(len({sheet_name for sheet_name, _, _ in entries}), 24)

        metadata_by_combo = {
            metadata["combination"]: metadata
            for _sheet_name, _bucket_df, metadata in entries
        }

        self.assertEqual(
            entries[0][2]["combination"],
            "BDSC > UAS (true) > sleep/ circ (true) > has balancer (false)",
        )
        self.assertEqual(
            metadata_by_combo["BDSC > UAS (true) > sleep/ circ (true) > has balancer (false)"]["stock_count"],
            1,
        )
        self.assertEqual(
            metadata_by_combo["BDSC > UAS (true) > sleep/ circ (true) > has balancer (false)"]["allele_count"],
            1,
        )
        self.assertEqual(
            metadata_by_combo["BDSC > UAS (true) > sleep/ circ (true) > has balancer (false)"]["gene_count"],
            1,
        )

        self.assertEqual(
            metadata_by_combo["VDRC > UAS (true) > sleep/ circ (true) > has balancer (false)"]["stock_count"],
            1,
        )
        self.assertEqual(
            metadata_by_combo["VDRC > UAS (true) > sleep/ circ (true) > has balancer (false)"]["allele_count"],
            1,
        )
        self.assertEqual(
            metadata_by_combo["VDRC > UAS (true) > sleep/ circ (true) > has balancer (false)"]["gene_count"],
            1,
        )

        self.assertEqual(
            metadata_by_combo["NIG > UAS (false) > sleep/ circ (false) > has balancer (true)"]["stock_count"],
            1,
        )
        self.assertEqual(
            metadata_by_combo["NIG > UAS (false) > sleep/ circ (false) > has balancer (true)"]["allele_count"],
            1,
        )
        self.assertEqual(
            metadata_by_combo["NIG > UAS (false) > sleep/ circ (false) > has balancer (true)"]["gene_count"],
            1,
        )

        self.assertEqual(sum(metadata["stock_count"] for _n, _d, metadata in entries), 3)
        self.assertEqual(sum(metadata["allele_count"] for _n, _d, metadata in entries), 3)
        self.assertEqual(sum(metadata["gene_count"] for _n, _d, metadata in entries), 3)
        self.assertTrue(
            any(metadata["stock_count"] == 0 for _sheet_name, _bucket_df, metadata in entries)
        )

    def test_soft_run_simple_buckets_outputs_combination_workbook(self):
        self.assertTrue(TEST_FIXTURE_DIR.exists(), f"Missing test fixture: {TEST_FIXTURE_DIR}")
        self.assertTrue(CONFIG_PATH.exists(), f"Missing config: {CONFIG_PATH}")

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            fixture_copy = tmp_root / "TEST"
            shutil.copytree(TEST_FIXTURE_DIR, fixture_copy)

            settings = Settings(
                soft_run=True,
                enable_oai_embedding=True,
                simple_buckets=True,
                openai_api_key="test-key",
                phenotype_embedding_cache_path=tmp_root / "cache" / "phenotype_embeddings.csv",
                phenotype_target_embedding_cache_path=tmp_root / "cache" / "target_embeddings.csv",
            )
            pipeline = StockSplittingPipeline(settings)

            with patch.object(
                EmbeddingSimilarityScorer,
                "_ensure_embeddings",
                new=_fake_ensure_embeddings,
            ):
                output_dir = pipeline.run(
                    input_dir=fixture_copy,
                    config_path=CONFIG_PATH,
                    verbose=False,
                )

            self.assertIsNotNone(output_dir)
            tier_workbook_path = (
                Path(output_dir) / "aggregated_stock_refs_aggregated_similarity_tiers.xlsx"
            )
            self.assertTrue(tier_workbook_path.exists(), tier_workbook_path)

            tier_workbook = pd.ExcelFile(tier_workbook_path)
            try:
                self.assertGreater(len(tier_workbook.sheet_names), 3)
                self.assertEqual(
                    tier_workbook.sheet_names[:3],
                    ["Contents", "Gene Set", "Stock Phenotype Sheet"],
                )

                expected_gene_set_df = pd.read_csv(fixture_copy / "vGAT_genes.csv", dtype=str)
                tier_gene_set_df = pd.read_excel(tier_workbook_path, sheet_name="Gene Set")
                pd.testing.assert_frame_equal(
                    tier_gene_set_df,
                    expected_gene_set_df,
                    check_dtype=False,
                )

                contents_df = pd.read_excel(
                    tier_workbook_path,
                    sheet_name="Contents",
                    header=None,
                )
                contents_text = "\n".join(
                    contents_df.fillna("").astype(str).agg(" ".join, axis=1).tolist()
                )
                self.assertIn("Simple Bucket Workbook Contents", contents_text)
                self.assertIn("Sheet name", contents_text)
                self.assertIn("has balancer", contents_text)
                self.assertNotIn("Max Cosine Similarity", contents_text)

                header_row_idx = next(
                    idx
                    for idx, row in contents_df.iterrows()
                    if row.fillna("").astype(str).tolist()[:3] == ["Sheet name", "Combination", "Collection"]
                )
                data_rows = contents_df.iloc[header_row_idx + 1 :].copy()
                data_rows = data_rows[data_rows.iloc[:, 0].fillna("").astype(str).str.strip().ne("")]
                sheet_names_from_contents = data_rows.iloc[:, 0].astype(str).tolist()
                self.assertEqual(tier_workbook.sheet_names[3:], sheet_names_from_contents)
                self.assertTrue(data_rows.iloc[:, 6].fillna(0).astype(int).ge(0).all())
                self.assertTrue(data_rows.iloc[:, 7].fillna(0).astype(int).ge(0).all())
                self.assertTrue(data_rows.iloc[:, 8].fillna(0).astype(int).ge(0).all())
            finally:
                tier_workbook.close()

    def test_similarity_tiers_follow_fixed_threshold_bins(self):
        phenotype_df = pd.DataFrame(
            {
                "Source/ Stock #": [f"stock-{idx}" for idx in range(12)],
                "Phenotype": [f"phenotype-{idx}" for idx in range(12)],
                "Reference": [f"reference-{idx}" for idx in range(12)],
                "Cosine Similarity (sleep)": [
                    0.95,
                    0.95,
                    0.82,
                    0.82,
                    0.71,
                    0.71,
                    0.58,
                    0.58,
                    0.41,
                    0.41,
                    0.05,
                    -0.20,
                ],
                "_reference_url": [f"https://example.org/{idx}" for idx in range(12)],
            }
        )

        tiers = _build_similarity_tier_sheets(phenotype_df)
        self.assertLessEqual(len(tiers), SIMILARITY_TIER_SHEET_COUNT)

        total_rows = 0
        expected_sheet_names = [
            "0.95-1",
            "0.8-0.85",
            "0.7-0.75",
            "0.55-0.6",
            "0.4-0.45",
            "0.05-0.1",
            "<0.05",
        ]
        expected_lengths = [2, 2, 2, 2, 2, 1, 1]
        expected_labels = [
            "0.95-1",
            "0.8-0.85",
            "0.7-0.75",
            "0.55-0.6",
            "0.4-0.45",
            "0.05-0.1",
            "<0.05",
        ]
        for idx, (sheet_name, tier_df, metadata) in enumerate(tiers):
            self.assertEqual(sheet_name, expected_sheet_names[idx])
            self.assertEqual(len(tier_df), expected_lengths[idx])
            self.assertEqual(metadata["range_label"], expected_labels[idx])
            if not tier_df.empty:
                self.assertTrue(tier_df["Similarity Range"].eq(expected_labels[idx]).all())
                tier_scores = pd.to_numeric(tier_df["Max Cosine Similarity"], errors="coerce")
                if metadata["lower_bound"] is None:
                    self.assertTrue((tier_scores < metadata["upper_bound"]).all())
                elif metadata["upper_bound"] >= 1.0:
                    self.assertTrue((tier_scores >= metadata["lower_bound"]).all())
                else:
                    self.assertTrue(
                        tier_scores.between(
                            metadata["lower_bound"], metadata["upper_bound"], inclusive="left"
                        ).all()
                    )
            total_rows += len(tier_df)

        self.assertEqual(total_rows, len(phenotype_df))

    def test_similarity_tiers_group_by_gene_and_sort_reagents_by_max_cosine(self):
        phenotype_df = pd.DataFrame(
            {
                "Gene": ["gene-b", "gene-a", "gene-b", "gene-a", "gene-a", "gene-c"],
                "Source/ Stock #": [
                    "stock-b1",
                    "stock-a1",
                    "stock-b2",
                    "stock-a2",
                    "stock-a3",
                    "stock-c1",
                ],
                "Phenotype": ["phb1", "pha1", "phb2", "pha2", "pha3", "phc1"],
                "Reference": ["ref-b1", "ref-a1", "ref-b2", "ref-a2", "ref-a3", "ref-c1"],
                "Cosine Similarity (sleep)": [0.58, 0.57, 0.56, 0.59, 0.55, 0.82],
            }
        )

        tiers = _build_similarity_tier_sheets(phenotype_df)
        self.assertEqual(
            [sheet_name for sheet_name, _, _ in tiers],
            ["0.8-0.85", "0.55-0.6"],
        )
        self.assertEqual(
            tiers[1][1]["Gene"].tolist(),
            ["gene-a", "gene-a", "gene-a", "gene-b", "gene-b"],
        )
        self.assertEqual(
            tiers[1][1]["Source/ Stock #"].tolist(),
            ["stock-a2", "stock-a1", "stock-a3", "stock-b1", "stock-b2"],
        )

    def test_similarity_tiers_handle_sparse_or_missing_scores(self):
        one_score_df = pd.DataFrame(
            {
                "Phenotype": ["a", "b", "c"],
                "Cosine Similarity (sleep)": [0.75, 0.75, 0.75],
            }
        )
        tiers = _build_similarity_tier_sheets(one_score_df)
        self.assertEqual(
            [(sheet_name, len(tier_df)) for sheet_name, tier_df, _ in tiers],
            [("0.75-0.8", 3)],
        )

        low_score_df = pd.DataFrame(
            {
                "Phenotype": ["low-a", "low-b"],
                "Cosine Similarity (sleep)": [0.05, -0.15],
            }
        )
        low_score_tiers = _build_similarity_tier_sheets(low_score_df)
        self.assertEqual(
            [(sheet_name, len(tier_df)) for sheet_name, tier_df, _ in low_score_tiers],
            [("0.05-0.1", 1), ("<0.05", 1)],
        )

        no_score_df = pd.DataFrame(
            {
                "Phenotype": ["a", "b"],
                "Reference": ["ref-a", "ref-b"],
            }
        )
        self.assertEqual(_build_similarity_tier_sheets(no_score_df), [])


if __name__ == "__main__":
    unittest.main()
