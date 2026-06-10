from __future__ import annotations

import gzip
import importlib.util
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
from fl_ai_reagent_stocker.pipelines.stock_finding import (  # noqa: E402
    StockFindingPipeline,
    _derive_one_hot_reagent_buckets,
)
from fl_ai_reagent_stocker.pipelines.stock_splitting import (  # noqa: E402
    CO_REAGENT_FBIDS_COLUMN,
    CO_REAGENT_SYMBOLS_COLUMN,
    MASTERLIST_TEMPLATE_COLUMNS,
    PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN,
    PARTNER_DRIVER_SYMBOLS_COLUMN,
    _build_stock_phenotype_sheet,
    _reorder_to_masterlist_columns,
)
from fl_ai_reagent_stocker.utils import (  # noqa: E402
    REAGENT_BUCKET_COLUMNS,
    RNAI_TYPE_COLUMN,
    infer_rnai_type_from_text,
)


def _load_fbsf_script_module():
    path = REPO_ROOT / "scripts" / "build_fbsf_to_fbgn_mapping.py"
    spec = importlib.util.spec_from_file_location("build_fbsf_to_fbgn_mapping_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _settings() -> Settings:
    return Settings(
        openai_api_key="test-key",
        ncbi_api_key="test-key",
        unpaywall_token="test-token",
    )


def _empty_insertions_df() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "Allele (id)",
            "Insertion (id)",
            "Insertion (symbol)",
            "Allele Class (term)",
            "Description (supporting reference)",
        ]
    )


def _empty_fbtp_to_fbti_df() -> pd.DataFrame:
    return pd.DataFrame(columns=["FBtp", "FBti"])


class TestBuildFbsfToFbgnMappingScript(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.script = _load_fbsf_script_module()

    def test_extract_fbsf_to_fbgn_pairs_from_sql_dump(self):
        sql_fixture = """\
COPY public.feature (feature_id, uniquename, type_id) FROM stdin;
1\tFBsf0000001\t1
2\tFBgn0000001\t2
3\tFBgn0000002\t2
4\tFBal0000001\t3
\\.
COPY public.feature_relationship (feature_relationship_id, subject_id, object_id, type_id) FROM stdin;
10\t1\t2\t99
11\t1\t3\t99
12\t4\t2\t99
\\.
"""

        with tempfile.TemporaryDirectory() as tmp_dir:
            sql_path = Path(tmp_dir) / "fixture.sql.gz"
            with gzip.open(sql_path, "wt", encoding="utf-8") as handle:
                handle.write(sql_fixture)

            pairs, summary = self.script.extract_fbsf_to_fbgn_pairs(sql_path)

        self.assertEqual(
            pairs,
            {
                ("FBsf0000001", "FBgn0000001"),
                ("FBsf0000001", "FBgn0000002"),
            },
        )
        self.assertEqual(summary["feature_rows"], 4)
        self.assertEqual(summary["feature_relationship_rows"], 3)


class TestPhenotypeReagentCapture(unittest.TestCase):
    def setUp(self):
        self.pipeline = StockFindingPipeline(_settings())

    def test_derive_one_hot_reagent_buckets_classifies_representative_cases(self):
        test_cases = [
            (
                "pure_uas",
                {
                    "genotype": "w[*]; P{UAS-gene}1",
                    "relevant_fbal_symbols": "gene[RNAi.UAS]",
                    "transgenic_product_class_terms": "rnai_reagent",
                    "RNAi": True,
                },
                "UAS",
            ),
            (
                "pure_gal4",
                {
                    "genotype": "Scer\\GAL4[repo]",
                    "relevant_fbal_symbols": "Scer\\GAL4[repo]",
                    "transgenic_product_class_terms": "driver",
                    "RNAi": False,
                },
                "GAL4",
            ),
            (
                "pure_mutant",
                {
                    "genotype": "Ddc[1]",
                    "relevant_fbal_symbols": "Ddc[1]",
                    "match_provenance": "direct_allele",
                    "RNAi": False,
                },
                "mutant",
            ),
            (
                "mutant_uas",
                {
                    "genotype": "Ddc[1]; P{UAS-gene}1",
                    "relevant_fbal_symbols": "Ddc[1]; gene[RNAi.UAS]",
                    "match_provenance": "direct_allele",
                    "transgenic_product_class_terms": "rnai_reagent",
                    "RNAi": True,
                },
                "mutant/UAS",
            ),
            (
                "gal4_mutant",
                {
                    "genotype": "Ddc[1] Scer\\GAL4[repo]",
                    "relevant_fbal_symbols": "Ddc[1]; Scer\\GAL4[repo]",
                    "match_provenance": "direct_allele; construct_regulatory_region",
                    "transgenic_product_class_terms": "driver",
                    "RNAi": False,
                },
                "GAL4 / mutant",
            ),
            (
                "gal4_mutant_uas",
                {
                    "genotype": "Ddc[1] Scer\\GAL4[repo] P{UAS-gene}1",
                    "relevant_fbal_symbols": "Ddc[1]; Scer\\GAL4[repo]; gene[RNAi.UAS]",
                    "match_provenance": "direct_allele; construct_regulatory_region",
                    "transgenic_product_class_terms": "driver; rnai_reagent",
                    "RNAi": True,
                },
                "GAL4 / mutant",
            ),
            (
                "gal4_uas_only",
                {
                    "genotype": "Scer\\GAL4[repo] P{UAS-gene}1",
                    "relevant_fbal_symbols": "Scer\\GAL4[repo]; gene[RNAi.UAS]",
                    "transgenic_product_class_terms": "driver; rnai_reagent",
                    "RNAi": True,
                },
                "Other",
            ),
        ]

        for case_name, payload, expected_bucket in test_cases:
            with self.subTest(case_name=case_name):
                result = _derive_one_hot_reagent_buckets(pd.Series(payload))
                self.assertEqual(
                    [bucket for bucket in REAGENT_BUCKET_COLUMNS if result[bucket]],
                    [expected_bucket],
                )

    def test_infer_rnai_type_from_text_distinguishes_shrna_and_dsrna(self):
        self.assertEqual(
            infer_rnai_type_from_text(
                "UAS regulatory sequences drive expression of an artificial microRNA (shRNA) targeting dpp."
            ),
            "shRNA",
        )
        self.assertEqual(
            infer_rnai_type_from_text(
                "UAS regulatory sequences drive expression of two copies of cv-c sequence, arranged in an inverted repeat."
            ),
            "dsRNA",
        )

    def test_build_stock_mapping_keeps_direct_alleles_with_provenance(self):
        fbal_to_fbgn_df = pd.DataFrame(
            [
                {
                    "GeneID": "FBgn0000422",
                    "GeneSymbol": "Ddc",
                    "AlleleID": "FBal0000001",
                    "AlleleSymbol": "Ddc[1]",
                }
            ]
        )
        derived_components_df = pd.DataFrame(
            [
                {
                    "FBst": "FBst1000001",
                    "stock_number": "1001",
                    "collection": "BDSC",
                    "FB_genotype": "Ddc[1]",
                    "derived_stock_component": "FBal0000001",
                    "embedded_type": "FBal",
                    "object_symbol": "Ddc[1]",
                    "GeneSymbol": "Ddc",
                }
            ]
        )

        result = self.pipeline._build_stock_mapping(
            input_genes=["FBgn0000422"],
            derived_components_df=derived_components_df,
            fbal_to_fbgn_df=fbal_to_fbgn_df,
            fbsf_to_fbgn_df=pd.DataFrame(columns=["FBsf", "FBgn"]),
            transgenic_constructs_df=pd.DataFrame(),
            insertion_alleles_df=_empty_insertions_df(),
            fbtp_to_fbti_df=_empty_fbtp_to_fbti_df(),
            fb_stocks_df=None,
        )

        self.assertEqual(result["FBst"].tolist(), ["FBst1000001"])
        self.assertEqual(result["relevant_flybase_gene_ids"].tolist(), ["FBgn0000422"])
        self.assertEqual(result["match_provenance"].tolist(), ["direct_allele"])

    def test_build_gene_component_tables_captures_ddc_regulatory_region_match(self):
        fbal_to_fbgn_df = pd.DataFrame(
            [
                {
                    "GeneID": "FBgn0000422",
                    "GeneSymbol": "Ddc",
                    "AlleleID": "FBal0000001",
                    "AlleleSymbol": "Ddc[1]",
                },
                {
                    "GeneID": "FBgn0014445",
                    "GeneSymbol": "Scer\\GAL4",
                    "AlleleID": "FBal0104725",
                    "AlleleSymbol": "Scer\\GAL4[Ddc.PL]",
                },
            ]
        )
        transgenic_constructs_df = pd.DataFrame(
            [
                {
                    "Component Allele (id)": "FBal0104725",
                    "Component Allele (symbol)": "Scer\\GAL4[Ddc.PL]",
                    "Regulatory region (id)": "FBgn0000422",
                    "Encoded product/tool (id)": "FBto0000001",
                    "Transgenic Construct (id)": "FBtp0104725",
                    "Transgenic Construct (symbol)": "P{GawB}Ddc",
                    "Transgenic Product class (term)": "driver",
                    "Description (supporting reference)": "FBrf0000001",
                }
            ]
        )
        derived_components_df = pd.DataFrame(
            [
                {
                    "FBst": "FBst1000001",
                    "stock_number": "1001",
                    "collection": "BDSC",
                    "FB_genotype": "Ddc[1]",
                    "derived_stock_component": "FBal0000001",
                    "embedded_type": "FBal",
                    "object_symbol": "Ddc[1]",
                    "GeneSymbol": "Ddc",
                },
                {
                    "FBst": "FBst1000002",
                    "stock_number": "1002",
                    "collection": "BDSC",
                    "FB_genotype": "Scer\\GAL4[Ddc.PL]",
                    "derived_stock_component": "FBal0104725",
                    "embedded_type": "FBal",
                    "object_symbol": "Scer\\GAL4[Ddc.PL]",
                    "GeneSymbol": "Scer\\GAL4",
                },
            ]
        )

        result = self.pipeline._build_stock_mapping(
            input_genes=["FBgn0000422"],
            derived_components_df=derived_components_df,
            fbal_to_fbgn_df=fbal_to_fbgn_df,
            fbsf_to_fbgn_df=pd.DataFrame(columns=["FBsf", "FBgn"]),
            transgenic_constructs_df=transgenic_constructs_df,
            insertion_alleles_df=_empty_insertions_df(),
            fbtp_to_fbti_df=_empty_fbtp_to_fbti_df(),
            fb_stocks_df=None,
        )

        self.assertEqual(set(result["FBst"]), {"FBst1000001", "FBst1000002"})
        construct_row = result[result["FBst"] == "FBst1000002"].iloc[0]
        self.assertEqual(construct_row["relevant_flybase_gene_ids"], "FBgn0000422")
        self.assertEqual(construct_row["relevant_gene_symbols"], "Ddc")
        self.assertEqual(
            construct_row["match_provenance"], "construct_regulatory_region"
        )
        self.assertEqual(construct_row["relevant_fbtp_ids"], "FBtp0104725")

    def test_build_gene_component_tables_resolves_clk_fbsf_match(self):
        fbal_to_fbgn_df = pd.DataFrame(
            [
                {
                    "GeneID": "FBgn0023076",
                    "GeneSymbol": "Clk",
                    "AlleleID": "FBal0000002",
                    "AlleleSymbol": "Clk[1]",
                },
                {
                    "GeneID": "FBgn0014445",
                    "GeneSymbol": "Scer\\GAL4",
                    "AlleleID": "FBal0264191",
                    "AlleleSymbol": "Scer\\GAL4[Clk.PX]",
                },
            ]
        )
        fbsf_to_fbgn_df = pd.DataFrame(
            [{"FBsf": "FBsf0000872949", "FBgn": "FBgn0023076"}]
        )
        transgenic_constructs_df = pd.DataFrame(
            [
                {
                    "Component Allele (id)": "FBal0264191",
                    "Component Allele (symbol)": "Scer\\GAL4[Clk.PX]",
                    "Regulatory region (id)": "FBsf0000872949",
                    "Encoded product/tool (id)": "FBto0000001",
                    "Transgenic Construct (id)": "FBtp0264191",
                    "Transgenic Construct (symbol)": "P{GawB}Clk",
                    "Transgenic Product class (term)": "driver",
                    "Description (supporting reference)": "FBrf0000002",
                }
            ]
        )
        derived_components_df = pd.DataFrame(
            [
                {
                    "FBst": "FBst2000001",
                    "stock_number": "2001",
                    "collection": "BDSC",
                    "FB_genotype": "Scer\\GAL4[Clk.PX]",
                    "derived_stock_component": "FBal0264191",
                    "embedded_type": "FBal",
                    "object_symbol": "Scer\\GAL4[Clk.PX]",
                    "GeneSymbol": "Scer\\GAL4",
                }
            ]
        )

        result = self.pipeline._build_stock_mapping(
            input_genes=["FBgn0023076"],
            derived_components_df=derived_components_df,
            fbal_to_fbgn_df=fbal_to_fbgn_df,
            fbsf_to_fbgn_df=fbsf_to_fbgn_df,
            transgenic_constructs_df=transgenic_constructs_df,
            insertion_alleles_df=_empty_insertions_df(),
            fbtp_to_fbti_df=_empty_fbtp_to_fbti_df(),
            fb_stocks_df=None,
        )

        self.assertEqual(result["FBst"].tolist(), ["FBst2000001"])
        self.assertEqual(result["relevant_flybase_gene_ids"].tolist(), ["FBgn0023076"])
        self.assertEqual(result["relevant_gene_symbols"].tolist(), ["Clk"])
        self.assertEqual(
            result["match_provenance"].tolist(), ["construct_regulatory_region"]
        )

    def test_stock_phenotype_sheet_uses_all_stocks_df(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "9001",
                        "collection": "BDSC",
                        "FB_genotype": "Scer\\GAL4[Ddc.PL]",
                        "derived_stock_component": "FBal0999999",
                        "embedded_type": "FBal",
                        "object_symbol": "Scer\\GAL4[Ddc.PL]",
                        "GeneSymbol": "Scer\\GAL4",
                    }
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            phenotype_tsv = (
                "reference\tphenotype_name\tphenotype_id\tqualifier_names\t"
                "genotype_FBids\tgenotype_symbols\n"
                "FBrf0000001\tabnormal locomotor rhythm\tFBcv0000001\tabnormal\t"
                "FBal0999999\tScer\\\\GAL4[Ddc.PL]\n"
            )
            (alleles_dir / "genotype_phenotype_data_fb_test.tsv").write_text(
                phenotype_tsv,
                encoding="utf-8",
            )
            refs_tsv = (
                "FBrf\tPMID\tPMCID\tDOI\tminiref\n"
                "FBrf0000001\t\t\t\tTest et al., 2024, Journal\n"
            )
            (refs_dir / "fbrf_pmid_pmcid_doi_fb_test.tsv").write_text(
                refs_tsv,
                encoding="utf-8",
            )

            all_stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "9001",
                        "collection": "BDSC",
                        "relevant_gene_symbols": "Ddc",
                        "relevant_component_ids": "FBal0999999",
                        "relevant_fbal_ids": "FBal0999999",
                        "relevant_fbal_symbols": "Scer\\GAL4[Ddc.PL]",
                        "Balancers": "-",
                        "matched_component_types": "FBal",
                        "UAS": False,
                        "GAL4": True,
                        "mutant/UAS": False,
                        "mutant": False,
                        "GAL4 / mutant": False,
                        "Other": False,
                        "allele_class_terms": "gain_of_function_allele",
                        "transgenic_product_class_terms": "gal4",
                        "custom_stock": False,
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=all_stocks_df,
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
                gene_to_datasets={"Ddc": {"dataset_alpha"}},
            )

        self.assertEqual(len(result), 1)
        self.assertIn("Phenotype", result.columns)
        self.assertIn("PMID", result.columns)
        self.assertIn("PMCID", result.columns)
        self.assertIn("Balancers", result.columns)
        self.assertIn("matched_component_types", result.columns)
        for column in REAGENT_BUCKET_COLUMNS:
            self.assertIn(column, result.columns)
        self.assertIn("allele_class_terms", result.columns)
        self.assertIn("transgenic_product_class_terms", result.columns)
        self.assertIn(RNAI_TYPE_COLUMN, result.columns)
        self.assertIn(CO_REAGENT_FBIDS_COLUMN, result.columns)
        self.assertIn(CO_REAGENT_SYMBOLS_COLUMN, result.columns)
        self.assertIn(PARTNER_DRIVER_SYMBOLS_COLUMN, result.columns)
        self.assertIn(PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN, result.columns)
        self.assertEqual(result["Phenotype"].iloc[0], "abnormal locomotor rhythm")
        self.assertEqual(result["PMID"].iloc[0], "")
        self.assertEqual(result["PMCID"].iloc[0], "")
        self.assertEqual(result["Balancers"].iloc[0], "-")
        self.assertEqual(result["matched_component_types"].iloc[0], "FBal")
        self.assertTrue(result["GAL4"].iloc[0])
        self.assertEqual(result["Data Set"].iloc[0], "dataset_alpha")
        self.assertEqual(
            int(result.loc[0, REAGENT_BUCKET_COLUMNS].astype(bool).sum()),
            1,
        )
        self.assertEqual(result["allele_class_terms"].iloc[0], "gain_of_function_allele")
        self.assertEqual(result["transgenic_product_class_terms"].iloc[0], "gal4")
        self.assertEqual(result[RNAI_TYPE_COLUMN].iloc[0], "")
        self.assertEqual(result[CO_REAGENT_FBIDS_COLUMN].iloc[0], "")
        self.assertEqual(result[CO_REAGENT_SYMBOLS_COLUMN].iloc[0], "")
        self.assertEqual(result[PARTNER_DRIVER_SYMBOLS_COLUMN].iloc[0], "")
        self.assertEqual(result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0], "")

    def test_stock_phenotype_sheet_recovers_partner_driver_candidates(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "FB_genotype": "P{TRiP.HMS00001}",
                        "derived_stock_component": "FBal0999999",
                        "embedded_type": "FBal",
                        "object_symbol": "P{TRiP.HMS00001}",
                        "GeneSymbol": "CG31400",
                    },
                    {
                        "FBst": "FBst9000002",
                        "stock_number": "7126",
                        "collection": "BDSC",
                        "FB_genotype": "tim-GAL4",
                        "derived_stock_component": "FBal0888888",
                        "embedded_type": "FBal",
                        "object_symbol": "tim-GAL4",
                        "GeneSymbol": "Scer\\GAL4",
                    },
                    {
                        "FBst": "FBst9000003",
                        "stock_number": "7127",
                        "collection": "BDSC",
                        "FB_genotype": "tim-GAL4",
                        "derived_stock_component": "FBal0888888",
                        "embedded_type": "FBal",
                        "object_symbol": "tim-GAL4",
                        "GeneSymbol": "Scer\\GAL4",
                    },
                    {
                        "FBst": "FBst9000004",
                        "stock_number": "7128",
                        "collection": "BDSC",
                        "FB_genotype": "tim-GAL4 / UAS-GFP",
                        "derived_stock_component": "FBal0888888",
                        "embedded_type": "FBal",
                        "object_symbol": "tim-GAL4",
                        "GeneSymbol": "Scer\\GAL4",
                    },
                    {
                        "FBst": "FBst9000004",
                        "stock_number": "7128",
                        "collection": "BDSC",
                        "FB_genotype": "tim-GAL4 / UAS-GFP",
                        "derived_stock_component": "FBal0777777",
                        "embedded_type": "FBal",
                        "object_symbol": "UAS-GFP",
                        "GeneSymbol": "GFP",
                    },
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            phenotype_tsv = (
                "reference\tphenotype_name\tphenotype_id\tqualifier_names\t"
                "genotype_FBids\tgenotype_symbols\n"
                "FBrf0000001\tabnormal circadian rhythm\tFBcv0000002\tabnormal\t"
                "FBal0999999/FBal0888888\tP{TRiP.HMS00001}/tim-GAL4\n"
            )
            (alleles_dir / "genotype_phenotype_data_fb_test.tsv").write_text(
                phenotype_tsv,
                encoding="utf-8",
            )
            refs_tsv = (
                "FBrf\tPMID\tPMCID\tDOI\tminiref\n"
                "FBrf0000001\t\t\t\tDriver et al., 2025, Journal\n"
            )
            (refs_dir / "fbrf_pmid_pmcid_doi_fb_test.tsv").write_text(
                refs_tsv,
                encoding="utf-8",
            )

            all_stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "relevant_gene_symbols": "tim",
                        "relevant_component_ids": "FBal0999999",
                        "relevant_fbal_ids": "FBal0999999",
                        "relevant_fbal_symbols": "P{TRiP.HMS00001}",
                        "Balancers": "-",
                        "matched_component_types": "FBal",
                        "UAS": True,
                        "GAL4": False,
                        "mutant/UAS": False,
                        "mutant": False,
                        "GAL4 / mutant": False,
                        "Other": False,
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "custom_stock": False,
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=all_stocks_df,
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
            )

        self.assertEqual(len(result), 1)
        self.assertEqual(result["Phenotype"].iloc[0], "abnormal circadian rhythm")
        self.assertEqual(result[RNAI_TYPE_COLUMN].iloc[0], "shRNA")
        self.assertEqual(result[CO_REAGENT_FBIDS_COLUMN].iloc[0], "FBal0888888")
        self.assertEqual(result[CO_REAGENT_SYMBOLS_COLUMN].iloc[0], "tim-GAL4")
        self.assertEqual(result[PARTNER_DRIVER_SYMBOLS_COLUMN].iloc[0], "tim-GAL4")
        self.assertIn(
            "(7126, BDSC)",
            result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0],
        )
        self.assertIn(
            "(7127, BDSC)",
            result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0],
        )
        self.assertNotIn(
            "(7128, BDSC)",
            result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0],
        )

    def test_stock_phenotype_sheet_filters_partner_driver_symbols_when_only_mixed_stock_exists(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000100",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "FB_genotype": "P{TRiP.HMS00001}",
                        "derived_stock_component": "FBal4999999",
                        "embedded_type": "FBal",
                        "object_symbol": "P{TRiP.HMS00001}",
                        "GeneSymbol": "CG31400",
                    },
                    {
                        "FBst": "FBst9000101",
                        "stock_number": "8123",
                        "collection": "BDSC",
                        "FB_genotype": "tim-GAL4 / UAS-GFP",
                        "derived_stock_component": "FBal4888888",
                        "embedded_type": "FBal",
                        "object_symbol": "tim-GAL4",
                        "GeneSymbol": "Scer\\GAL4",
                    },
                    {
                        "FBst": "FBst9000101",
                        "stock_number": "8123",
                        "collection": "BDSC",
                        "FB_genotype": "tim-GAL4 / UAS-GFP",
                        "derived_stock_component": "FBal4777777",
                        "embedded_type": "FBal",
                        "object_symbol": "UAS-GFP",
                        "GeneSymbol": "GFP",
                    },
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            phenotype_tsv = (
                "reference\tphenotype_name\tphenotype_id\tqualifier_names\t"
                "genotype_FBids\tgenotype_symbols\n"
                "FBrf0000012\tabnormal circadian rhythm\tFBcv0000002\tabnormal\t"
                "FBal4999999/FBal4888888\tP{TRiP.HMS00001}/tim-GAL4\n"
            )
            (alleles_dir / "genotype_phenotype_data_fb_test.tsv").write_text(
                phenotype_tsv,
                encoding="utf-8",
            )
            refs_tsv = (
                "FBrf\tPMID\tPMCID\tDOI\tminiref\n"
                "FBrf0000012\t\t\t\tFiltered et al., 2025, Journal\n"
            )
            (refs_dir / "fbrf_pmid_pmcid_doi_fb_test.tsv").write_text(
                refs_tsv,
                encoding="utf-8",
            )

            all_stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000100",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "relevant_gene_symbols": "tim",
                        "relevant_component_ids": "FBal4999999",
                        "relevant_fbal_ids": "FBal4999999",
                        "relevant_fbal_symbols": "P{TRiP.HMS00001}",
                        "Balancers": "-",
                        "matched_component_types": "FBal",
                        "UAS": True,
                        "GAL4": False,
                        "mutant/UAS": False,
                        "mutant": False,
                        "GAL4 / mutant": False,
                        "Other": False,
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "custom_stock": False,
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=all_stocks_df,
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
            )

        self.assertEqual(result[PARTNER_DRIVER_SYMBOLS_COLUMN].iloc[0], "")
        self.assertEqual(result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0], "")

    def test_stock_phenotype_sheet_prefers_exact_genotype_driver_symbol(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000010",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "FB_genotype": "P{TRiP.HMS00001}",
                        "derived_stock_component": "FBal1999999",
                        "embedded_type": "FBal",
                        "object_symbol": "P{TRiP.HMS00001}",
                        "GeneSymbol": "CG31400",
                    },
                    {
                        "FBst": "FBst9000011",
                        "stock_number": "8001",
                        "collection": "BDSC",
                        "FB_genotype": "tim-GAL4",
                        "derived_stock_component": "FBal1888888",
                        "embedded_type": "FBal",
                        "object_symbol": "Scer\\GAL4[tim]",
                        "GeneSymbol": "Scer\\GAL4",
                    },
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            phenotype_tsv = (
                "reference\tphenotype_name\tphenotype_id\tqualifier_names\t"
                "genotype_FBids\tgenotype_symbols\n"
                "FBrf0000009\tabnormal circadian rhythm\tFBcv0000002\tabnormal\t"
                "FBal1999999/FBal1888888\tP{TRiP.HMS00001}/tim-GAL4\n"
            )
            (alleles_dir / "genotype_phenotype_data_fb_test.tsv").write_text(
                phenotype_tsv,
                encoding="utf-8",
            )
            refs_tsv = (
                "FBrf\tPMID\tPMCID\tDOI\tminiref\n"
                "FBrf0000009\t\t\t\tDriver et al., 2025, Journal\n"
            )
            (refs_dir / "fbrf_pmid_pmcid_doi_fb_test.tsv").write_text(
                refs_tsv,
                encoding="utf-8",
            )

            all_stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000010",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "relevant_gene_symbols": "tim",
                        "relevant_component_ids": "FBal1999999",
                        "relevant_fbal_ids": "FBal1999999",
                        "relevant_fbal_symbols": "P{TRiP.HMS00001}",
                        "Balancers": "-",
                        "matched_component_types": "FBal",
                        "UAS": True,
                        "GAL4": False,
                        "mutant/UAS": False,
                        "mutant": False,
                        "GAL4 / mutant": False,
                        "Other": False,
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "custom_stock": False,
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=all_stocks_df,
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
            )

        self.assertEqual(result[CO_REAGENT_SYMBOLS_COLUMN].iloc[0], "tim-GAL4")
        self.assertEqual(result[PARTNER_DRIVER_SYMBOLS_COLUMN].iloc[0], "tim-GAL4")
        self.assertEqual(
            result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0],
            "(8001, BDSC)",
        )

    def test_stock_phenotype_sheet_recovers_gal4_partner_when_symbol_pairing_is_imperfect(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000030",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "FB_genotype": "P{TRiP.HMS00001}",
                        "derived_stock_component": "FBal3999999",
                        "embedded_type": "FBal",
                        "object_symbol": "P{TRiP.HMS00001}",
                        "GeneSymbol": "CG31400",
                    },
                    {
                        "FBst": "FBst9000031",
                        "stock_number": "8765",
                        "collection": "BDSC",
                        "FB_genotype": "tim-GAL4",
                        "derived_stock_component": "FBal3888888",
                        "embedded_type": "FBal",
                        "object_symbol": "",
                        "GeneSymbol": "",
                    },
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            phenotype_tsv = (
                "reference\tphenotype_name\tphenotype_id\tqualifier_names\t"
                "genotype_FBids\tgenotype_symbols\n"
                "FBrf0000011\tabnormal circadian rhythm\tFBcv0000002\tabnormal\t"
                "FBal3999999/FBal3888888\tP{TRiP.HMS00001}; tim-GAL4; extra_marker\n"
            )
            (alleles_dir / "genotype_phenotype_data_fb_test.tsv").write_text(
                phenotype_tsv,
                encoding="utf-8",
            )
            refs_tsv = (
                "FBrf\tPMID\tPMCID\tDOI\tminiref\n"
                "FBrf0000011\t\t\t\tImperfect et al., 2025, Journal\n"
            )
            (refs_dir / "fbrf_pmid_pmcid_doi_fb_test.tsv").write_text(
                refs_tsv,
                encoding="utf-8",
            )

            all_stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000030",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "relevant_gene_symbols": "tim",
                        "relevant_component_ids": "FBal3999999",
                        "relevant_fbal_ids": "FBal3999999",
                        "relevant_fbal_symbols": "P{TRiP.HMS00001}",
                        "Balancers": "-",
                        "matched_component_types": "FBal",
                        "UAS": True,
                        "GAL4": False,
                        "mutant/UAS": False,
                        "mutant": False,
                        "GAL4 / mutant": False,
                        "Other": False,
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "custom_stock": False,
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=all_stocks_df,
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
            )

        self.assertEqual(result[PARTNER_DRIVER_SYMBOLS_COLUMN].iloc[0], "tim-GAL4")
        self.assertEqual(
            result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0],
            "(8765, BDSC)",
        )

    def test_stock_phenotype_sheet_excludes_non_gal4_partner_reagents(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000020",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "FB_genotype": "P{TRiP.HMS00001}",
                        "derived_stock_component": "FBal2999999",
                        "embedded_type": "FBal",
                        "object_symbol": "P{TRiP.HMS00001}",
                        "GeneSymbol": "CG31400",
                    },
                    {
                        "FBst": "FBst9000021",
                        "stock_number": "9002",
                        "collection": "BDSC",
                        "FB_genotype": "LexAop-GFP",
                        "derived_stock_component": "FBal2888888",
                        "embedded_type": "FBal",
                        "object_symbol": "LexAop-GFP",
                        "GeneSymbol": "GFP",
                    },
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            phenotype_tsv = (
                "reference\tphenotype_name\tphenotype_id\tqualifier_names\t"
                "genotype_FBids\tgenotype_symbols\n"
                "FBrf0000010\tabnormal circadian rhythm\tFBcv0000002\tabnormal\t"
                "FBal2999999/FBal2888888\tP{TRiP.HMS00001}/LexAop-GFP\n"
            )
            (alleles_dir / "genotype_phenotype_data_fb_test.tsv").write_text(
                phenotype_tsv,
                encoding="utf-8",
            )
            refs_tsv = (
                "FBrf\tPMID\tPMCID\tDOI\tminiref\n"
                "FBrf0000010\t\t\t\tNonGal4 et al., 2025, Journal\n"
            )
            (refs_dir / "fbrf_pmid_pmcid_doi_fb_test.tsv").write_text(
                refs_tsv,
                encoding="utf-8",
            )

            all_stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000020",
                        "stock_number": "v9241",
                        "collection": "Vienna",
                        "relevant_gene_symbols": "tim",
                        "relevant_component_ids": "FBal2999999",
                        "relevant_fbal_ids": "FBal2999999",
                        "relevant_fbal_symbols": "P{TRiP.HMS00001}",
                        "Balancers": "-",
                        "matched_component_types": "FBal",
                        "UAS": True,
                        "GAL4": False,
                        "mutant/UAS": False,
                        "mutant": False,
                        "GAL4 / mutant": False,
                        "Other": False,
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "custom_stock": False,
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=all_stocks_df,
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
            )

        self.assertEqual(result[CO_REAGENT_FBIDS_COLUMN].iloc[0], "")
        self.assertEqual(result[CO_REAGENT_SYMBOLS_COLUMN].iloc[0], "")
        self.assertEqual(result[PARTNER_DRIVER_SYMBOLS_COLUMN].iloc[0], "")
        self.assertEqual(result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0], "")

    def test_reorder_to_masterlist_columns_matches_template_order(self):
        sheet_df = pd.DataFrame(
            [
                {
                    "Gene": "tim",
                    "Reagent Type or Allele Symbol": "HMS00001",
                    "Balancers": "-",
                    "matched_component_types": "FBal",
                    "UAS": True,
                    "GAL4": False,
                    "mutant/UAS": False,
                    "mutant": False,
                    "GAL4 / mutant": False,
                    "Other": False,
                    "allele_class_terms": "RNAi_reagent",
                    "transgenic_product_class_terms": "rnai_reagent",
                    RNAI_TYPE_COLUMN: "shRNA",
                    "Source/ Stock #": "BDSC (12345)",
                    "Genotype": "P{TRiP.HMS00001}",
                    CO_REAGENT_FBIDS_COLUMN: "FBal0001",
                    CO_REAGENT_SYMBOLS_COLUMN: "tim-GAL4",
                    PARTNER_DRIVER_SYMBOLS_COLUMN: "tim-GAL4",
                    PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN: "(7126, BDSC)",
                    "Data Set": "dataset_alpha",
                    "Phenotype": "abnormal circadian rhythm",
                    "Qualifier": "abnormal",
                    "PMID": "12345678",
                    "PMCID": "PMC123456",
                    "Reference": "Title of paper",
                    "Authors": "Author A; Author B",
                    "Journal": "Cell Reports",
                    "Year of Publication": "2025",
                    "Cosine Similarity (sleep)": 0.42,
                    "Cosine Similarity (circadian)": 0.87,
                    "Similarity Range": "0.85-0.90",
                }
            ]
        )

        result = _reorder_to_masterlist_columns(sheet_df)

        self.assertEqual(
            result.columns[: len(MASTERLIST_TEMPLATE_COLUMNS)].tolist(),
            MASTERLIST_TEMPLATE_COLUMNS,
        )
        self.assertEqual(result["Stock Source"].iloc[0], "BDSC")
        self.assertEqual(result["ID #"].iloc[0], "12345")
        self.assertEqual(result["allele shorthand"].iloc[0], "HMS00001")
        self.assertEqual(
            result["Published Gal4/ Positive control"].iloc[0],
            "tim-GAL4",
        )
        self.assertEqual(result["Reference"].iloc[0], "12345678")
        self.assertEqual(result["Column 31"].iloc[0], "PMC123456")
        self.assertEqual(result["Column 30"].iloc[0], "Title of paper")
        self.assertEqual(
            result["Circadian/Sleep Relevance (embedding max score)"].iloc[0],
            0.87,
        )
        self.assertEqual(result["Data Set"].iloc[0], "dataset_alpha")
        self.assertEqual(result[RNAI_TYPE_COLUMN].iloc[0], "shRNA")
        self.assertIn("Similarity Range", result.columns)
        self.assertIn("matched_component_types", result.columns)
        self.assertIn("Cosine Similarity (sleep)", result.columns)
        self.assertIn("Cosine Similarity (circadian)", result.columns)


def _write_fbrf_refs(refs_dir: Path) -> None:
    """Write a minimal fbrf_pmid_pmcid_doi fixture used across tests."""
    refs_tsv = (
        "FBrf\tPMID\tPMCID\tDOI\tminiref\n"
        "FBrf0000001\t\t\t\tTest et al., 2024, Journal\n"
    )
    (refs_dir / "fbrf_pmid_pmcid_doi_fb_test.tsv").write_text(
        refs_tsv,
        encoding="utf-8",
    )


def _write_phenotype_tsv(
    alleles_dir: Path,
    phenotype_rows: list[str],
) -> None:
    """Write a minimal genotype_phenotype_data fixture from row strings."""
    header = (
        "reference\tphenotype_name\tphenotype_id\tqualifier_names\t"
        "genotype_FBids\tgenotype_symbols"
    )
    body = "\n".join([header, *phenotype_rows]) + "\n"
    (alleles_dir / "genotype_phenotype_data_fb_test.tsv").write_text(
        body,
        encoding="utf-8",
    )


def _empty_all_stocks_df() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "FBst",
            "stock_number",
            "collection",
            "relevant_gene_symbols",
            "relevant_component_ids",
            "relevant_fbal_ids",
            "relevant_fbal_symbols",
            "Balancers",
            "matched_component_types",
            "UAS",
            "GAL4",
            "mutant/UAS",
            "mutant",
            "GAL4 / mutant",
            "Other",
            "allele_class_terms",
            "transgenic_product_class_terms",
            "custom_stock",
        ]
    )


class TestPhenotypeSheetGeneFirstFlow(unittest.TestCase):
    """Regression tests for the gene-first All Phenotypic Stocks Sheet flow."""

    def test_find_stocks_writes_reagent_index_when_no_stocks_found(self):
        """Stage 1 should persist Gene Reagent Index even with zero stock rows."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_dir = Path(tmp_dir) / "genes"
            input_dir.mkdir()
            pd.DataFrame(
                [
                    {
                        "flybase_gene_id": "FBgn0000422",
                        "ext_gene": "Ddc",
                    }
                ]
            ).to_csv(input_dir / "genes.csv", index=False)

            pipeline = StockFindingPipeline(_settings())
            tables = {
                "derived_stock_components": pd.DataFrame(
                    columns=[
                        "FBst",
                        "stock_number",
                        "collection",
                        "FB_genotype",
                        "derived_stock_component",
                        "embedded_type",
                        "object_symbol",
                        "GeneSymbol",
                    ]
                ),
                "fbal_to_fbgn": pd.DataFrame(
                    [
                        {
                            "GeneID": "FBgn0000422",
                            "GeneSymbol": "Ddc",
                            "AlleleID": "FBal0999999",
                            "AlleleSymbol": "Ddc[1]",
                        }
                    ]
                ),
                "fbsf_to_fbgn": pd.DataFrame(columns=["FBsf", "FBgn"]),
                "transgenic_construct_descriptions": pd.DataFrame(),
                "insertion_allele_descriptions": _empty_insertions_df(),
                "fbtp_to_fbti": _empty_fbtp_to_fbti_df(),
                "fb_stocks": pd.DataFrame(),
                "genotype_phenotype_data": pd.DataFrame(),
                "entity_publication": pd.DataFrame(),
                "ref_metadata": pd.DataFrame(),
            }

            with patch.object(pipeline, "_load_tables", return_value=tables):
                output_path = pipeline.run(
                    input_dir=input_dir,
                    keywords=[],
                    skip_fbgnid_conversion=True,
                )

            self.assertIsNotNone(output_path)
            output_path = Path(output_path)
            self.assertTrue(output_path.exists())
            workbook = pd.ExcelFile(output_path)
            self.assertIn("Stocks", workbook.sheet_names)
            self.assertIn("Gene Reagent Index", workbook.sheet_names)
            stock_sheet = pd.read_excel(output_path, sheet_name="Stocks")
            reagent_index = pd.read_excel(
                output_path,
                sheet_name="Gene Reagent Index",
            )

        self.assertEqual(len(stock_sheet), 0)
        self.assertEqual(reagent_index["component_id"].tolist(), ["FBal0999999"])
        self.assertEqual(reagent_index["input_gene_symbol"].tolist(), ["Ddc"])

    def test_reagent_index_surfaces_phenotype_reagent_absent_from_all_stocks_df(self):
        """Plan step 8 case 1: reagent index covers a reagent not in stocks_df."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            # Global lookup includes the FBst, but stocks_df is empty.
            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "9001",
                        "collection": "BDSC",
                        "FB_genotype": "Ddc[1]",
                        "derived_stock_component": "FBal0999999",
                        "embedded_type": "FBal",
                        "object_symbol": "Ddc[1]",
                        "GeneSymbol": "Ddc",
                    }
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            _write_phenotype_tsv(
                alleles_dir,
                [
                    "FBrf0000001\tabnormal locomotor rhythm\tFBcv0000001\tabnormal\t"
                    "FBal0999999\tDdc[1]",
                ],
            )
            _write_fbrf_refs(refs_dir)

            reagent_index_df = pd.DataFrame(
                [
                    {
                        "input_flybase_gene_id": "FBgn0000422",
                        "input_gene_symbol": "Ddc",
                        "component_id": "FBal0999999",
                        "component_type": "FBal",
                        "component_symbol": "Ddc[1]",
                        "source_allele_id": "FBal0999999",
                        "source_allele_symbol": "Ddc[1]",
                        "allele_class_terms": "loss_of_function_allele",
                        "transgenic_product_class_terms": "",
                        "rnai_type": "",
                        "match_provenance": "direct_allele",
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=_empty_all_stocks_df(),
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
                gene_to_datasets={"Ddc": {"dataset_alpha"}},
                reagent_index_df=reagent_index_df,
            )

        self.assertEqual(len(result), 1)
        self.assertEqual(result["Phenotype"].iloc[0], "abnormal locomotor rhythm")
        self.assertEqual(result["Gene"].iloc[0], "Ddc")
        self.assertIn("9001", result["Source/ Stock #"].iloc[0])
        self.assertEqual(result["matched_component_types"].iloc[0], "FBal")
        self.assertEqual(result["allele_class_terms"].iloc[0], "loss_of_function_allele")
        self.assertEqual(
            int(result.loc[0, REAGENT_BUCKET_COLUMNS].astype(bool).sum()),
            1,
        )

    def test_reagent_index_surfaces_all_global_fbsts(self):
        """Plan step 8 case 2: every FBst FlyBase associates with the component."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            # Three FBsts back the same FBal in the global lookup; only the
            # first one made it into the Stage 1 stocks_df. The new flow must
            # surface all three.
            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "9001",
                        "collection": "BDSC",
                        "FB_genotype": "P{TRiP.HMS00001}",
                        "derived_stock_component": "FBal0999999",
                        "embedded_type": "FBal",
                        "object_symbol": "Ddc[RNAi.HMS00001]",
                        "GeneSymbol": "Ddc",
                    },
                    {
                        "FBst": "FBst9000002",
                        "stock_number": "9002",
                        "collection": "BDSC",
                        "FB_genotype": "P{TRiP.HMS00001}",
                        "derived_stock_component": "FBal0999999",
                        "embedded_type": "FBal",
                        "object_symbol": "Ddc[RNAi.HMS00001]",
                        "GeneSymbol": "Ddc",
                    },
                    {
                        "FBst": "FBst9000003",
                        "stock_number": "9003",
                        "collection": "VDRC",
                        "FB_genotype": "P{TRiP.HMS00001}",
                        "derived_stock_component": "FBal0999999",
                        "embedded_type": "FBal",
                        "object_symbol": "Ddc[RNAi.HMS00001]",
                        "GeneSymbol": "Ddc",
                    },
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            _write_phenotype_tsv(
                alleles_dir,
                [
                    "FBrf0000001\tabnormal locomotor rhythm\tFBcv0000001\tabnormal\t"
                    "FBal0999999\tDdc[RNAi.HMS00001]",
                ],
            )
            _write_fbrf_refs(refs_dir)

            reagent_index_df = pd.DataFrame(
                [
                    {
                        "input_flybase_gene_id": "FBgn0000422",
                        "input_gene_symbol": "Ddc",
                        "component_id": "FBal0999999",
                        "component_type": "FBal",
                        "component_symbol": "Ddc[RNAi.HMS00001]",
                        "source_allele_id": "FBal0999999",
                        "source_allele_symbol": "Ddc[RNAi.HMS00001]",
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "rnai_type": "shRNA",
                        "match_provenance": "direct_allele",
                    }
                ]
            )

            # Stage 1 stocks_df only contains one of the three.
            all_stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "9001",
                        "collection": "BDSC",
                        "relevant_gene_symbols": "Ddc",
                        "relevant_component_ids": "FBal0999999",
                        "relevant_fbal_ids": "FBal0999999",
                        "relevant_fbal_symbols": "Ddc[RNAi.HMS00001]",
                        "Balancers": "-",
                        "matched_component_types": "FBal",
                        "UAS": True,
                        "GAL4": False,
                        "mutant/UAS": False,
                        "mutant": False,
                        "GAL4 / mutant": False,
                        "Other": False,
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "custom_stock": False,
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=all_stocks_df,
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
                gene_to_datasets={"Ddc": {"dataset_alpha"}},
                reagent_index_df=reagent_index_df,
            )

        source_stocks = set(result["Source/ Stock #"].tolist())
        # All three stocks must appear in the sheet.
        self.assertEqual(len(result), 3)
        self.assertTrue(any("9001" in label for label in source_stocks))
        self.assertTrue(any("9002" in label for label in source_stocks))
        self.assertTrue(any("9003" in label for label in source_stocks))
        # All rows should resolve to a one-hot reagent bucket.
        for _, row in result.iterrows():
            n_buckets = int(sum(bool(row[col]) for col in REAGENT_BUCKET_COLUMNS))
            self.assertEqual(n_buckets, 1)

    def test_phenotype_fbal_attaches_stock_keyed_by_linked_fbti(self):
        """Phenotype FBal rows should find stocks stored under linked FBti IDs."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            # The phenotype row references the FBal, but FlyBase's global
            # stock-component lookup records the stock under the linked FBti.
            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000101",
                        "stock_number": "9101",
                        "collection": "BDSC",
                        "FB_genotype": "P{UAS-Ddc.RNAi}attP2",
                        "derived_stock_component": "FBti0777777",
                        "embedded_type": "FBti",
                        "object_symbol": "P{UAS-Ddc.RNAi}attP2",
                        "GeneSymbol": "Ddc",
                    }
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            _write_phenotype_tsv(
                alleles_dir,
                [
                    "FBrf0000001\tabnormal locomotor rhythm\tFBcv0000001\tabnormal\t"
                    "FBal0999999\tDdc[RNAi.HMS00001]",
                ],
            )
            _write_fbrf_refs(refs_dir)

            reagent_index_df = pd.DataFrame(
                [
                    {
                        "input_flybase_gene_id": "FBgn0000422",
                        "input_gene_symbol": "Ddc",
                        "component_id": "FBal0999999",
                        "component_type": "FBal",
                        "component_symbol": "Ddc[RNAi.HMS00001]",
                        "source_allele_id": "FBal0999999",
                        "source_allele_symbol": "Ddc[RNAi.HMS00001]",
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "rnai_type": "shRNA",
                        "match_provenance": "direct_allele",
                    },
                    {
                        "input_flybase_gene_id": "FBgn0000422",
                        "input_gene_symbol": "Ddc",
                        "component_id": "FBti0777777",
                        "component_type": "FBti",
                        "component_symbol": "P{UAS-Ddc.RNAi}attP2",
                        "source_allele_id": "FBal0999999",
                        "source_allele_symbol": "Ddc[RNAi.HMS00001]",
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "rnai_type": "shRNA",
                        "match_provenance": "direct_allele",
                    },
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=_empty_all_stocks_df(),
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
                gene_to_datasets={"Ddc": {"dataset_alpha"}},
                reagent_index_df=reagent_index_df,
            )

        self.assertEqual(len(result), 1)
        self.assertIn("9101", result["Source/ Stock #"].iloc[0])
        self.assertNotIn("No-stock phenotype reagent", result["Source/ Stock #"].iloc[0])
        self.assertEqual(result["Reagent Type or Allele Symbol"].iloc[0], "Ddc[RNAi.HMS00001]")
        self.assertEqual(result["matched_component_types"].iloc[0], "FBal")

    def test_reagent_index_no_stock_phenotype_row(self):
        """Plan step 8 case 3: phenotype reagent with no FBst."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            # No row in the global lookup for FBal0888888.
            pd.DataFrame(
                columns=[
                    "FBst",
                    "stock_number",
                    "collection",
                    "FB_genotype",
                    "derived_stock_component",
                    "embedded_type",
                    "object_symbol",
                    "GeneSymbol",
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            _write_phenotype_tsv(
                alleles_dir,
                [
                    "FBrf0000001\tabnormal locomotor rhythm\tFBcv0000001\tabnormal\t"
                    "FBal0888888\tDdc[lab.custom]",
                ],
            )
            _write_fbrf_refs(refs_dir)

            reagent_index_df = pd.DataFrame(
                [
                    {
                        "input_flybase_gene_id": "FBgn0000422",
                        "input_gene_symbol": "Ddc",
                        "component_id": "FBal0888888",
                        "component_type": "FBal",
                        "component_symbol": "Ddc[lab.custom]",
                        "source_allele_id": "FBal0888888",
                        "source_allele_symbol": "Ddc[lab.custom]",
                        "allele_class_terms": "loss_of_function_allele",
                        "transgenic_product_class_terms": "",
                        "rnai_type": "",
                        "match_provenance": "direct_allele",
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=_empty_all_stocks_df(),
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
                gene_to_datasets={"Ddc": {"dataset_alpha"}},
                reagent_index_df=reagent_index_df,
            )

        self.assertEqual(len(result), 1)
        self.assertEqual(result["Phenotype"].iloc[0], "abnormal locomotor rhythm")
        self.assertEqual(result["Gene"].iloc[0], "Ddc")
        # No-stock phenotype reagent rows use the component label as Source/ Stock #.
        self.assertIn("Ddc[lab.custom]", result["Source/ Stock #"].iloc[0])
        self.assertEqual(result["allele_class_terms"].iloc[0], "loss_of_function_allele")
        # Reagent bucket must still be one-hot.
        self.assertEqual(
            int(result.loc[0, REAGENT_BUCKET_COLUMNS].astype(bool).sum()),
            1,
        )

    def test_co_reagent_not_in_reagent_index_is_informational_only(self):
        """Plan step 8 case 4: co-reagent stocks must not emit focal rows."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            # Two FBsts: one for the focal input-gene reagent (FBal0999999)
            # and one for a co-reagent driver (FBal0888777) that is NOT in
            # the input-gene reagent index. The driver stock must not be
            # emitted as a focal Source/ Stock # row.
            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "9001",
                        "collection": "BDSC",
                        "FB_genotype": "Ddc[RNAi.HMS00001]",
                        "derived_stock_component": "FBal0999999",
                        "embedded_type": "FBal",
                        "object_symbol": "Ddc[RNAi.HMS00001]",
                        "GeneSymbol": "Ddc",
                    },
                    {
                        "FBst": "FBst9000099",
                        "stock_number": "9099",
                        "collection": "BDSC",
                        "FB_genotype": "Scer\\GAL4[tim]",
                        "derived_stock_component": "FBal0888777",
                        "embedded_type": "FBal",
                        "object_symbol": "Scer\\GAL4[tim]",
                        "GeneSymbol": "Scer\\GAL4",
                    },
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            _write_phenotype_tsv(
                alleles_dir,
                [
                    "FBrf0000001\tabnormal locomotor rhythm\tFBcv0000001\tabnormal\t"
                    "FBal0999999/FBal0888777\tDdc[RNAi.HMS00001]/Scer\\\\GAL4[tim]",
                ],
            )
            _write_fbrf_refs(refs_dir)

            # Reagent index only contains the focal input-gene reagent.
            reagent_index_df = pd.DataFrame(
                [
                    {
                        "input_flybase_gene_id": "FBgn0000422",
                        "input_gene_symbol": "Ddc",
                        "component_id": "FBal0999999",
                        "component_type": "FBal",
                        "component_symbol": "Ddc[RNAi.HMS00001]",
                        "source_allele_id": "FBal0999999",
                        "source_allele_symbol": "Ddc[RNAi.HMS00001]",
                        "allele_class_terms": "RNAi_reagent",
                        "transgenic_product_class_terms": "rnai_reagent",
                        "rnai_type": "shRNA",
                        "match_provenance": "direct_allele",
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=_empty_all_stocks_df(),
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
                gene_to_datasets={"Ddc": {"dataset_alpha"}},
                reagent_index_df=reagent_index_df,
            )

        # Only the focal stock 9001 should appear as a Source/ Stock #.
        source_stocks = set(result["Source/ Stock #"].tolist())
        self.assertTrue(any("9001" in label for label in source_stocks))
        self.assertFalse(any("9099" in label for label in source_stocks))
        self.assertEqual(len(result), 1)
        # The driver should still appear as a partner-driver hint.
        partner_drivers = str(
            result[PARTNER_DRIVER_SYMBOLS_COLUMN].iloc[0] or ""
        ).lower()
        self.assertIn("gal4", partner_drivers)

    def test_legacy_flow_without_reagent_index_still_works(self):
        """Older Stage 1 workbooks without `Gene Reagent Index` keep working."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            flybase_root = Path(tmp_dir) / "flybase"
            alleles_dir = flybase_root / "alleles_and_stocks"
            refs_dir = flybase_root / "references"
            alleles_dir.mkdir(parents=True)
            refs_dir.mkdir(parents=True)

            pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "9001",
                        "collection": "BDSC",
                        "FB_genotype": "Scer\\GAL4[Ddc.PL]",
                        "derived_stock_component": "FBal0999999",
                        "embedded_type": "FBal",
                        "object_symbol": "Scer\\GAL4[Ddc.PL]",
                        "GeneSymbol": "Scer\\GAL4",
                    }
                ]
            ).to_csv(alleles_dir / "fbst_to_derived_stock_component.csv", index=False)

            _write_phenotype_tsv(
                alleles_dir,
                [
                    "FBrf0000001\tabnormal locomotor rhythm\tFBcv0000001\tabnormal\t"
                    "FBal0999999\tScer\\\\GAL4[Ddc.PL]",
                ],
            )
            _write_fbrf_refs(refs_dir)

            all_stocks_df = pd.DataFrame(
                [
                    {
                        "FBst": "FBst9000001",
                        "stock_number": "9001",
                        "collection": "BDSC",
                        "relevant_gene_symbols": "Ddc",
                        "relevant_component_ids": "FBal0999999",
                        "relevant_fbal_ids": "FBal0999999",
                        "relevant_fbal_symbols": "Scer\\GAL4[Ddc.PL]",
                        "Balancers": "-",
                        "matched_component_types": "FBal",
                        "UAS": False,
                        "GAL4": True,
                        "mutant/UAS": False,
                        "mutant": False,
                        "GAL4 / mutant": False,
                        "Other": False,
                        "allele_class_terms": "gain_of_function_allele",
                        "transgenic_product_class_terms": "gal4",
                        "custom_stock": False,
                    }
                ]
            )

            result = _build_stock_phenotype_sheet(
                all_stocks_df=all_stocks_df,
                flybase_data_path=flybase_root,
                references_df=None,
                unfiltered_references_df=None,
                pubmed_cache_path=None,
                pubmed_client=None,
                similarity_targets=None,
                embedding_scorer=None,
                verbose=False,
                gene_to_datasets={"Ddc": {"dataset_alpha"}},
                # No reagent_index_df: legacy fallback path.
            )

        # Same behavior as before the rework: one row for the matched stock.
        self.assertEqual(len(result), 1)
        self.assertEqual(result["Phenotype"].iloc[0], "abnormal locomotor rhythm")
        self.assertEqual(result["Balancers"].iloc[0], "-")
        self.assertTrue(result["GAL4"].iloc[0])


if __name__ == "__main__":
    unittest.main()
