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
    CO_REAGENT_FBIDS_COLUMN,
    CO_REAGENT_SYMBOLS_COLUMN,
    PHENOTYPE_SIMILARITY_EMBEDDING_MODEL,
    PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN,
    PARTNER_DRIVER_SYMBOLS_COLUMN,
    SIMILARITY_TIER_SHEET_COUNT,
    StockSplittingPipeline,
    _build_keyword_bucketing_sheets,
    _build_simple_bucket_workbook_entries,
    _build_similarity_tier_sheets,
    _compute_max_cosine_similarity,
    _reorder_to_masterlist_columns,
    load_split_config,
    write_aggregated_excel,
)
from fl_ai_reagent_stocker.utils import (  # noqa: E402
    REAGENT_BUCKET_COLUMNS,
    RNAI_TYPE_COLUMN,
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
    def _phenotype_settings_kwargs(self, similarity: str) -> dict:
        from fl_ai_reagent_stocker.cli import _similarity_to_settings_kwargs

        return _similarity_to_settings_kwargs(similarity)

    def test_contents_sheet_is_researcher_facing(self):
        config = {
            "settings": {
                "relevantSearchTerms": ["sleep", "circadian"],
                "phenotypeSimilarityTargets": [
                    {"keyword": "sleep", "embedding_text": "sleep"},
                ],
                "maxStocksPerGene": 5,
                "maxStocksPerAllele": 3,
            },
            "filterDescriptions": {
                "Bloomington": "Stocks from the Bloomington collection in FlyBase.",
                "RNAi": "Unified RNAi proxy.",
                "No sgRNA": "Genotype does not include sgRNA.",
                "no_other_insertion": "Stock has only one unique transgenic construct.",
                "Ref-": "Allele has no publication citations.",
                "UAS": "UAS (Upstream Activating Sequence) driver / GAL4-responsive construct.",
            },
            "filters": {
                "Bloomington": {
                    "column": "collection",
                    "type": "contains",
                    "value": "Bloomington",
                },
                "RNAi": {"column": "RNAi", "type": "equals", "value": True},
                "No sgRNA": {"column": "sgRNA", "type": "equals", "value": False},
                "no_other_insertion": {
                    "column": "multiple_insertions",
                    "type": "equals",
                    "value": False,
                },
                "Ref-": {
                    "column": "ALLELE_PAPER_RELEVANCE_SCORE",
                    "type": "equals",
                    "value": 0,
                },
                "UAS": {"column": "UAS", "type": "equals", "value": True},
            },
            "combinations": [
                ["Bloomington", "RNAi", "No sgRNA", "no_other_insertion", "Ref-"],
                ["UAS", "No sgRNA"],
            ],
        }
        included_df = pd.DataFrame(
            [
                {
                    "stock_number": "100",
                    "relevant_gene_symbols": "geneA",
                    "relevant_flybase_gene_ids": "FBgn0000001",
                    "AlleleSymbol": "alleleA",
                },
                {
                    "stock_number": "101",
                    "relevant_gene_symbols": pd.NA,
                    "relevant_flybase_gene_ids": "FBgn0000002",
                    "AlleleSymbol": "alleleB",
                },
                {
                    "stock_number": "102",
                    "relevant_gene_symbols": "sharedSymbol",
                    "relevant_flybase_gene_ids": "FBgn0000010",
                    "AlleleSymbol": "alleleC",
                },
                {
                    "stock_number": "103",
                    "relevant_gene_symbols": "sharedSymbol",
                    "relevant_flybase_gene_ids": "FBgn0000011",
                    "AlleleSymbol": "alleleD",
                },
            ]
        )
        empty_df = included_df.iloc[0:0].copy()
        combination_outputs = [
            (
                ["Bloomington", "RNAi", "No sgRNA", "no_other_insertion", "Ref-"],
                included_df,
                {
                    "Category": "Bloomington >> RNAi >> No sgRNA >> no_other_insertion >> Ref-",
                    "# Stocks": 4,
                    "# Alleles": 4,
                    "# Genes": 4,
                },
            ),
            (
                ["UAS", "No sgRNA"],
                empty_df,
                {
                    "Category": "UAS >> No sgRNA",
                    "# Stocks": 0,
                    "# Alleles": 0,
                    "# Genes": 0,
                },
            ),
        ]

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "researcher_contents.xlsx"
            write_aggregated_excel(
                output_path=output_path,
                source_workbook_path=Path(tmp_dir) / "source.xlsx",
                config=config,
                combination_outputs=combination_outputs,
                all_stocks_df=pd.DataFrame(),
                references_df=None,
                verbose=False,
                all_input_genes={
                    "FBgn0000001",
                    "FBgn0000002",
                    "FBgn0000010",
                    "FBgn0000011",
                },
                genes_no_stocks={"geneD"},
                csv_input_genes={
                    "geneA",
                    "geneB",
                    "geneC",
                    "geneD",
                    "geneAA",
                    "geneBB",
                },
                csv_input_fbgns={
                    "FBgn0000001",
                    "FBgn0000002",
                    "FBgn0000010",
                    "FBgn0000011",
                    "FBgn0000099",
                    "FBgn0000100",
                },
                n_input_genes=6,
                current_to_input_map={
                    "FBgn0000001": "geneA",
                    "FBgn0000002": "geneB",
                    "FBgn0000010": "geneAA",
                    "FBgn0000011": "geneBB",
                    "FBgn0000099": "geneC",
                    "FBgn0000100": "geneD",
                },
            )

            contents_df = pd.read_excel(output_path, sheet_name="Contents", header=None)

        first_col = contents_df.iloc[:, 0].fillna("").astype(str).tolist()
        contents_text = "\n".join(
            contents_df.fillna("").astype(str).agg(" ".join, axis=1).tolist()
        )

        def row_index(label: str) -> int:
            return next(i for i, value in enumerate(first_col) if value == label)

        section_order = [
            "What this workbook is",
            "How sheets were selected",
            "What each filter means",
            "Prioritized sheet breakdown",
            "Genes that need follow-up",
            "Per-sheet gene rosters",
            "References and Stock Sheet by Gene definitions",
        ]
        section_indices = [row_index(label) for label in section_order]
        self.assertEqual(section_indices, sorted(section_indices))

        self.assertIn("A stock appears in only one sheet", contents_text)
        self.assertIn("Sheet criteria", contents_text)
        self.assertIn("Total (unique stocks / alleles / genes across sheets)", contents_text)
        self.assertIn("Total Input Genes (across all sets): 6", contents_text)
        self.assertIn("Total Genes Included in output sheets: 4", contents_text)
        # FBgn0000010 and FBgn0000011 share the display symbol "sharedSymbol"
        # but must each be counted as a distinct input gene.
        self.assertIn("Sheet1 (4 stocks, 4 unique genes)", contents_text)
        self.assertIn("geneA", contents_text)
        self.assertIn("geneAA", contents_text)
        self.assertIn("geneBB", contents_text)
        self.assertNotIn("Keyword reference count column", contents_text)
        self.assertNotIn("multiple_insertions equals false", contents_text)
        self.assertNotIn("ALLELE_PAPER_RELEVANCE_SCORE equals 0", contents_text)

        matched_idx = next(
            i
            for i, value in enumerate(first_col)
            if value.startswith("Genes with matched stocks but not included")
        )
        no_stocks_idx = next(
            i
            for i, value in enumerate(first_col)
            if value.startswith("Input Genes with 0 matched stocks")
        )
        roster_idx = row_index("Per-sheet gene rosters")
        matched_block = set(first_col[matched_idx + 1 : no_stocks_idx])
        no_stocks_block = set(first_col[no_stocks_idx + 1 : roster_idx])
        # geneC has no matched stocks (its FBgn doesn't appear in the
        # combination_outputs) and so must appear in the no-stocks block.
        self.assertIn("geneC", no_stocks_block)
        self.assertIn("geneD", no_stocks_block)
        # geneC is also not in the "matched but not in any sheet" list.
        self.assertNotIn("geneC", matched_block)
        self.assertNotIn("geneD", matched_block)

    def test_input_gene_detection_collapses_secondary_fbgns(self):
        """Two input FBgns aliased to the same primary FBgn collapse to one
        canonical gene, with stock rows recognized regardless of whether they
        carry the secondary or the primary FBgn."""
        stocks_df = pd.DataFrame(
            [
                {
                    "stock_number": "200",
                    "relevant_gene_symbols": "oldName",
                    "relevant_flybase_gene_ids": "FBgn0001000",
                },
                {
                    "stock_number": "201",
                    "relevant_gene_symbols": "primaryName",
                    "relevant_flybase_gene_ids": "FBgn0099000",
                },
                {
                    "stock_number": "202",
                    "relevant_gene_symbols": "other",
                    "relevant_flybase_gene_ids": "FBgn0002000",
                },
            ]
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_dir = Path(tmp_dir) / "GeneSets" / "Stocks"
            input_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"flybase_gene_id": "FBgn0001000", "Fly Symbol": "oldName"},
                    {"flybase_gene_id": "FBgn0099000", "Fly Symbol": "primaryName"},
                    {"flybase_gene_id": "FBgn0002000", "Fly Symbol": "other"},
                ]
            ).to_csv(input_dir.parent / "input.csv", index=False)

            with patch(
                "fl_ai_reagent_stocker.pipelines.stock_splitting._load_secondary_fbgn_to_primary",
                return_value={"FBgn0001000": "FBgn0099000"},
            ):
                pipeline = StockSplittingPipeline(Settings())
                (
                    genes_no_stocks,
                    csv_input_genes,
                    current_to_input_map,
                    n_input_genes,
                    _gene_to_datasets,
                    csv_input_fbgns,
                    secondary_to_primary,
                ) = pipeline._find_genes_without_stocks(
                    input_dir,
                    stocks_df,
                    verbose=False,
                )

        self.assertEqual(secondary_to_primary, {"FBgn0001000": "FBgn0099000"})
        # The two input rows (oldName / primaryName) collapse to one
        # canonical gene, plus the unrelated "other" gene -> 2 total.
        self.assertEqual(csv_input_fbgns, {"FBgn0099000", "FBgn0002000"})
        self.assertEqual(n_input_genes, 2)
        self.assertEqual(genes_no_stocks, set())
        # Both the secondary FBgn and the primary FBgn must resolve to an
        # input-side display label so stock rows with either ID render
        # correctly.
        self.assertEqual(
            current_to_input_map.get("FBgn0001000"),
            current_to_input_map.get("FBgn0099000"),
        )

    def test_input_gene_detection_uses_fly_symbol_and_fbgn_fallbacks(self):
        stocks_df = pd.DataFrame(
            [
                {
                    "relevant_gene_symbols": pd.NA,
                    "relevant_flybase_gene_ids": "FBgn0004364",
                },
                {
                    "relevant_gene_symbols": "known",
                    "relevant_flybase_gene_ids": "FBgn0000001",
                },
            ]
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_dir = Path(tmp_dir) / "GeneSets" / "Stocks"
            input_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"flybase_gene_id": "FBgn0004364", "Fly Symbol": "18w"},
                    {"flybase_gene_id": "FBgn0000001", "Fly Symbol": "known"},
                ]
            ).to_csv(input_dir.parent / "input.csv", index=False)

            with patch(
                "fl_ai_reagent_stocker.pipelines.stock_splitting._load_secondary_fbgn_to_primary",
                return_value={},
            ):
                pipeline = StockSplittingPipeline(Settings())
                (
                    genes_no_stocks,
                    csv_input_genes,
                    current_to_input_map,
                    n_input_genes,
                    _gene_to_datasets,
                    csv_input_fbgns,
                    secondary_to_primary,
                ) = pipeline._find_genes_without_stocks(
                    input_dir,
                    stocks_df,
                    verbose=False,
                )

        self.assertEqual(genes_no_stocks, set())
        self.assertEqual(n_input_genes, 2)
        self.assertIn("18w", csv_input_genes)
        self.assertIn("known", csv_input_genes)
        self.assertEqual(current_to_input_map["FBgn0004364"], "18w")
        self.assertEqual(csv_input_fbgns, {"FBgn0004364", "FBgn0000001"})
        self.assertIsInstance(secondary_to_primary, dict)

    def test_cli_parses_phenotype_sheet_similarity_none(self):
        parser = create_parser()
        args = parser.parse_args(
            [
                "phenotype-sheet",
                "./gene_lists/Stocks",
                "--similarity",
                "none",
            ]
        )
        self.assertEqual(args.command, "phenotype-sheet")
        self.assertEqual(args.similarity, "none")
        kwargs = self._phenotype_settings_kwargs(args.similarity)
        self.assertEqual(
            kwargs,
            {
                "soft_run": True,
                "enable_oai_embedding": False,
                "simple_buckets": False,
                "keyword_bucketing": False,
            },
        )

    def test_cli_parses_phenotype_sheet_similarity_tiers(self):
        parser = create_parser()
        args = parser.parse_args(
            [
                "phenotype-sheet",
                "./gene_lists/Stocks",
                "--similarity",
                "tiers",
            ]
        )
        self.assertEqual(args.command, "phenotype-sheet")
        self.assertEqual(args.similarity, "tiers")
        kwargs = self._phenotype_settings_kwargs(args.similarity)
        self.assertEqual(
            kwargs,
            {
                "soft_run": True,
                "enable_oai_embedding": True,
                "simple_buckets": False,
                "keyword_bucketing": False,
            },
        )

    def test_cli_parses_phenotype_sheet_similarity_simple_buckets(self):
        parser = create_parser()
        args = parser.parse_args(
            [
                "phenotype-sheet",
                "./gene_lists/Stocks",
                "--similarity",
                "simple-buckets",
            ]
        )
        self.assertEqual(args.command, "phenotype-sheet")
        self.assertEqual(args.similarity, "simple-buckets")
        kwargs = self._phenotype_settings_kwargs(args.similarity)
        self.assertEqual(
            kwargs,
            {
                "soft_run": True,
                "enable_oai_embedding": True,
                "simple_buckets": True,
                "keyword_bucketing": False,
            },
        )

    def test_cli_parses_phenotype_sheet_similarity_keyword_buckets(self):
        parser = create_parser()
        args = parser.parse_args(
            [
                "phenotype-sheet",
                "./gene_lists/Stocks",
                "--similarity",
                "keyword-buckets",
            ]
        )
        self.assertEqual(args.command, "phenotype-sheet")
        self.assertEqual(args.similarity, "keyword-buckets")
        kwargs = self._phenotype_settings_kwargs(args.similarity)
        self.assertEqual(
            kwargs,
            {
                "soft_run": True,
                "enable_oai_embedding": True,
                "simple_buckets": False,
                "keyword_bucketing": True,
            },
        )

    def test_cli_split_stocks_rejects_removed_phenotype_flags(self):
        parser = create_parser()
        for flag in (
            "--soft-run",
            "--OAI-embedding",
            "--simple-buckets",
            "--keyword-bucketing",
        ):
            with self.subTest(flag=flag):
                with self.assertRaises(SystemExit):
                    parser.parse_args(
                        [
                            "split-stocks",
                            "./gene_lists/Stocks",
                            flag,
                        ]
                    )

    def test_cli_validate_stocks_rejects_removed_phenotype_flags(self):
        parser = create_parser()
        for flag in (
            "--soft-run",
            "--OAI-embedding",
            "--simple-buckets",
            "--keyword-bucketing",
        ):
            with self.subTest(flag=flag):
                with self.assertRaises(SystemExit):
                    parser.parse_args(
                        [
                            "validate-stocks",
                            "./gene_lists/Stocks",
                            flag,
                        ]
                    )

    def test_cli_run_full_pipeline_modes(self):
        parser = create_parser()
        args_normal = parser.parse_args(
            [
                "run-full-pipeline",
                "./gene_lists",
            ]
        )
        self.assertEqual(args_normal.command, "run-full-pipeline")
        self.assertEqual(args_normal.mode, "normal")
        # --similarity defaults to "none" but is only honored in --mode phenotype.
        self.assertEqual(args_normal.similarity, "none")

        args_phenotype = parser.parse_args(
            [
                "run-full-pipeline",
                "./gene_lists",
                "--mode",
                "phenotype",
                "--similarity",
                "tiers",
            ]
        )
        self.assertEqual(args_phenotype.mode, "phenotype")
        self.assertEqual(args_phenotype.similarity, "tiers")

    def test_load_stocks_from_excel_backfills_rnai_type_from_flybase(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            flybase_root = tmp_root / "flybase"
            constructs_dir = flybase_root / "transgenic_constructs"
            constructs_dir.mkdir(parents=True)
            (constructs_dir / "transgenic_construct_descriptions_fb_test.tsv").write_text(
                "\n".join(
                    [
                        "\t".join(
                            [
                                "Component Allele (symbol)",
                                "Component Allele (id)",
                                "Transgenic Construct (symbol)",
                                "Transgenic Construct (id)",
                                "Transgenic Product class (term)",
                                "Transgenic Product class (id)",
                                "Regulatory region (symbol)",
                                "Regulatory region (id)",
                                "Encoded product/tool (symbol)",
                                "Encoded product/tool (id)",
                                "Tagged with (symbol)",
                                "Tagged with (id)",
                                "Also carries (symbol)",
                                "Also carries (id)",
                                "Description (text)",
                                "Description (supporting reference)",
                                "Stocks (number)",
                            ]
                        ),
                        "\t".join(
                            [
                                "dpp[RNAi.UAS.cNa]",
                                "FBal0000001",
                                "P{UAS-dpp.RNAi.N}",
                                "FBtp0000001",
                                "RNAi_reagent",
                                "FBcv:0000000",
                                "",
                                "",
                                "dpp",
                                "FBgn0000001",
                                "",
                                "",
                                "",
                                "",
                                "UAS regulatory sequences drive expression of an artificial microRNA (shRNA) targeting dpp.",
                                "FBrf0000001",
                                "1",
                            ]
                        ),
                    ]
                ),
                encoding="utf-8",
            )

            workbook_path = tmp_root / "stocks.xlsx"
            stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst0000001",
                        "stock_number": "12345",
                        "collection": "BDSC",
                        "genotype": "P{UAS-dpp.RNAi.N}",
                        "relevant_fbal_ids": "FBal0000001",
                        "relevant_fbal_symbols": "dpp[RNAi.UAS.cNa]",
                        "relevant_fbtp_ids": "FBtp0000001",
                        "relevant_fbtp_symbols": "P{UAS-dpp.RNAi.N}",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "RNAi": True,
                    }
                ]
            )
            with pd.ExcelWriter(workbook_path) as writer:
                stocks_df.to_excel(writer, sheet_name="Stocks", index=False)

            pipeline = StockSplittingPipeline(
                Settings(
                    flybase_data_path=flybase_root,
                )
            )
            (
                loaded_stocks_df,
                references_df,
                keywords,
                reagent_index_df,
            ) = pipeline._load_stocks_from_excel(
                workbook_path,
                verbose=False,
            )

            self.assertIsNone(references_df)
            self.assertIsNone(reagent_index_df)
            self.assertEqual(keywords, [])
            self.assertIn(RNAI_TYPE_COLUMN, loaded_stocks_df.columns)
            self.assertEqual(loaded_stocks_df[RNAI_TYPE_COLUMN].iloc[0], "shRNA")

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
            phenotype_masterlist_df = _reorder_to_masterlist_columns(phenotype_df).fillna("")
            expected_columns = [
                "Balancers",
                "matched_component_types",
                *REAGENT_BUCKET_COLUMNS,
                "allele_class_terms",
                "transgenic_product_class_terms",
                RNAI_TYPE_COLUMN,
                "Source",
                "Stock #",
                "Genotype",
                CO_REAGENT_FBIDS_COLUMN,
                CO_REAGENT_SYMBOLS_COLUMN,
                PARTNER_DRIVER_SYMBOLS_COLUMN,
                PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN,
                "Phenotype",
                "Qualifier",
                "PMID",
                "PMCID",
                "Cosine Similarity (sleep)",
                "Cosine Similarity (circadian)",
            ]
            for column in expected_columns:
                self.assertIn(column, phenotype_df.columns)
            self.assertNotIn("Source/ Stock #", phenotype_df.columns)
            source_idx = phenotype_df.columns.get_loc("Source")
            self.assertEqual(
                phenotype_df.columns[source_idx:source_idx + 2].tolist(),
                ["Source", "Stock #"],
            )
            self.assertNotIn("_reference_url", phenotype_df.columns)
            self.assertTrue(
                phenotype_df[REAGENT_BUCKET_COLUMNS].fillna(False).astype(bool).sum(axis=1).eq(1).all()
            )

            aggregated_contents_df = pd.read_excel(
                workbook_path,
                sheet_name="Contents",
                header=None,
            )
            aggregated_contents_text = "\n".join(
                aggregated_contents_df.fillna("").astype(str).agg(" ".join, axis=1).tolist()
            )
            self.assertIn("Stock Phenotype Sheet columns", aggregated_contents_text)
            self.assertIn("Source", aggregated_contents_text)
            self.assertIn("Stock #", aggregated_contents_text)
            self.assertIn("mutually exclusive one-hot set", aggregated_contents_text)
            self.assertIn("matched_component_types", aggregated_contents_text)
            self.assertIn(RNAI_TYPE_COLUMN, aggregated_contents_text)
            self.assertIn(CO_REAGENT_FBIDS_COLUMN, aggregated_contents_text)
            self.assertIn(PARTNER_DRIVER_SYMBOLS_COLUMN, aggregated_contents_text)
            self.assertIn(PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN, aggregated_contents_text)
            self.assertIn("GAL4 / mutant", aggregated_contents_text)
            self.assertIn("Vienna-style knockdown-family signals", aggregated_contents_text)

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
                ).fillna("")
                self.assertEqual(
                    tier_first_df.columns.tolist(),
                    phenotype_masterlist_df.columns.tolist(),
                )
                pd.testing.assert_frame_equal(
                    tier_first_df,
                    phenotype_masterlist_df,
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
                self.assertIn("Stock Phenotype Sheet columns", contents_text)
                self.assertIn("mutually exclusive one-hot set", contents_text)
                self.assertIn("mutant/UAS", contents_text)
                self.assertIn("rnai_reagent", contents_text)

                scored_rows = int(_compute_max_cosine_similarity(phenotype_df).notna().sum())
                tier_row_total = 0
                previous_max_score = None
                for sheet_name, expected_df, metadata in expected_tiers:
                    actual_df = pd.read_excel(tier_workbook_path, sheet_name=sheet_name)
                    self.assertEqual(len(actual_df), len(expected_df))
                    pd.testing.assert_frame_equal(
                        actual_df.fillna(""),
                        _reorder_to_masterlist_columns(expected_df).fillna(""),
                        check_dtype=False,
                    )
                    if not actual_df.empty:
                        self.assertIn("Similarity Range", actual_df.columns)
                        self.assertIn(
                            "Circadian/Sleep Relevance (embedding max score)",
                            actual_df.columns,
                        )
                        expected_scores = pd.to_numeric(
                            expected_df["Max Cosine Similarity"], errors="coerce"
                        )
                        self.assertTrue(
                            actual_df["Similarity Range"].eq(metadata["range_label"]).all()
                        )
                        if metadata["lower_bound"] is None:
                            self.assertTrue((expected_scores < metadata["upper_bound"]).all())
                        elif metadata["upper_bound"] >= 1.0:
                            self.assertTrue((expected_scores >= metadata["lower_bound"]).all())
                        else:
                            self.assertTrue(
                                expected_scores.between(
                                    metadata["lower_bound"],
                                    metadata["upper_bound"],
                                    inclusive="left",
                                ).all()
                            )
                        current_max_score = float(expected_scores.max())
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

    def test_keyword_bucketing_keeps_only_matching_rows_in_keyword_hits(self):
        phenotype_df = pd.DataFrame(
            {
                "Gene": ["gene-a", "gene-a", "gene-a", "gene-a", "gene-a", "gene-b"],
                "Source": ["BDSC", "BDSC", "VDRC", "VDRC", "NIG", "TRiP"],
                "Stock #": ["101", "101", "202", "202", "303", "404"],
                "Phenotype": [
                    "sleep phenotype",
                    "viability phenotype",
                    "unrelated phenotype",
                    "another unrelated phenotype",
                    "plain phenotype",
                    "no circadian change observed",
                ],
                "Reference": ["ref-1", "ref-2", "ref-3", "ref-4", "ref-5", "ref-6"],
                "Cosine Similarity (sleep)": [0.20, 0.10, 0.91, 0.88, 0.30, 0.95],
            }
        )

        entries = _build_keyword_bucketing_sheets(
            phenotype_df,
            keywords=["sleep", "circadian"],
        )

        self.assertEqual(
            [sheet_name for sheet_name, _, _ in entries],
            ["Keyword Hits", "No Keyword Hits"],
        )

        keyword_hits_df = entries[0][1]
        no_keyword_hits_df = entries[1][1]

        self.assertEqual(keyword_hits_df["Stock #"].tolist(), ["101", "404"])
        self.assertEqual(
            keyword_hits_df["Phenotype"].tolist(),
            ["sleep phenotype", "no circadian change observed"],
        )
        self.assertIn("Max Cosine Similarity", keyword_hits_df.columns)
        self.assertEqual(entries[0][2]["stock_count"], 2)
        self.assertTrue(
            keyword_hits_df["Phenotype"].str.contains("sleep|circadian", case=False, regex=True).all()
        )

        self.assertEqual(
            no_keyword_hits_df["Gene"].tolist(),
            ["gene-a", "gene-a", "gene-a", "gene-a"],
        )
        self.assertEqual(
            no_keyword_hits_df["Stock #"].tolist(),
            ["202", "202", "303", "101"],
        )
        self.assertIn("Max Cosine Similarity", no_keyword_hits_df.columns)
        self.assertEqual(entries[1][2]["stock_count"], 3)

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
                self.assertIn("Stock Phenotype Sheet columns", contents_text)
                self.assertIn("mutually exclusive one-hot set", contents_text)
                self.assertIn("GAL4 / mutant", contents_text)

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

    def test_soft_run_keyword_bucketing_outputs_keyword_workbook(self):
        self.assertTrue(TEST_FIXTURE_DIR.exists(), f"Missing test fixture: {TEST_FIXTURE_DIR}")
        self.assertTrue(CONFIG_PATH.exists(), f"Missing config: {CONFIG_PATH}")

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            fixture_copy = tmp_root / "TEST"
            shutil.copytree(TEST_FIXTURE_DIR, fixture_copy)

            settings = Settings(
                soft_run=True,
                enable_oai_embedding=True,
                keyword_bucketing=True,
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
            tier_workbook_path = (
                Path(output_dir) / "aggregated_stock_refs_aggregated_similarity_tiers.xlsx"
            )
            self.assertTrue(workbook_path.exists(), workbook_path)
            self.assertTrue(tier_workbook_path.exists(), tier_workbook_path)

            phenotype_df = pd.read_excel(workbook_path, sheet_name="Stock Phenotype Sheet")
            configured_keywords = load_split_config(CONFIG_PATH)["settings"]["relevantSearchTerms"]
            expected_entries = _build_keyword_bucketing_sheets(
                phenotype_df,
                keywords=configured_keywords,
            )

            tier_workbook = pd.ExcelFile(tier_workbook_path)
            try:
                self.assertEqual(
                    tier_workbook.sheet_names,
                    [
                        "Contents",
                        "Gene Set",
                        "Stock Phenotype Sheet",
                        "Keyword Hits",
                        "No Keyword Hits",
                    ],
                )

                contents_df = pd.read_excel(
                    tier_workbook_path,
                    sheet_name="Contents",
                    header=None,
                )
                contents_text = "\n".join(
                    contents_df.fillna("").astype(str).agg(" ".join, axis=1).tolist()
                )
                self.assertIn("Keyword Bucket Workbook Contents", contents_text)
                self.assertIn("Keyword Hits", contents_text)
                self.assertIn("No Keyword Hits", contents_text)
                self.assertIn("relevantSearchTerms", contents_text)
                self.assertIn("sleep, circadian, locomotor, rhythm", contents_text)
                self.assertIn("Stock Phenotype Sheet columns", contents_text)

                for sheet_name, expected_df, _metadata in expected_entries:
                    actual_df = pd.read_excel(tier_workbook_path, sheet_name=sheet_name)
                    pd.testing.assert_frame_equal(
                        actual_df.fillna(""),
                        _reorder_to_masterlist_columns(expected_df).fillna(""),
                        check_dtype=False,
                    )
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
                "Source": [
                    "VDRC",
                    "BDSC",
                    "VDRC",
                    "BDSC",
                    "BDSC",
                    "NIG",
                ],
                "Stock #": [
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
            tiers[1][1]["Stock #"].tolist(),
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
