from __future__ import annotations

import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

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
    PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN,
    PARTNER_DRIVER_SYMBOLS_COLUMN,
    _build_stock_phenotype_sheet,
)
from fl_ai_reagent_stocker.utils import REAGENT_BUCKET_COLUMNS  # noqa: E402


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
        self.assertEqual(
            int(result.loc[0, REAGENT_BUCKET_COLUMNS].astype(bool).sum()),
            1,
        )
        self.assertEqual(result["allele_class_terms"].iloc[0], "gain_of_function_allele")
        self.assertEqual(result["transgenic_product_class_terms"].iloc[0], "gal4")
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
        self.assertEqual(result[CO_REAGENT_FBIDS_COLUMN].iloc[0], "FBal0888888")
        self.assertEqual(result[CO_REAGENT_SYMBOLS_COLUMN].iloc[0], "tim-GAL4")
        self.assertEqual(result[PARTNER_DRIVER_SYMBOLS_COLUMN].iloc[0], "tim-GAL4")
        self.assertIn(
            "BDSC (7126) [FBst9000002]",
            result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0],
        )
        self.assertIn(
            "BDSC (7127) [FBst9000003]",
            result[PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN].iloc[0],
        )

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
            "BDSC (8001) [FBst9000011]",
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


if __name__ == "__main__":
    unittest.main()
