from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker.cli import (  # noqa: E402
    _discover_input_csvs,
    _validate_gene_columns,
)
from fl_ai_reagent_stocker.tx_omics_revision.categories import (  # noqa: E402
    audit_category_identities,
    build_csw_4plus,
    build_hlh,
    build_homeostatic,
    build_mechanistic,
    load_master_tables,
)
from fl_ai_reagent_stocker.tx_omics_revision.constants import (  # noqa: E402
    BREAKDOWN_DIR,
    CSW_4PLUS_SOURCES,
    GO_ANALYSIS_DIR,
    GO_SCRIPT_BLOBS,
    HLH_PUBLICATION_GENES,
    MASTER_TABLES,
    MECHANISTIC_IDENTITIES,
    SET_NAMES,
    STOCKER_DIR,
    TRACHEALESS_FBGN,
    TRHN_FBGN,
    TRHN_SYMBOL,
    TX_INPUT_DIR,
    audit_paths,
)
from fl_ai_reagent_stocker.tx_omics_revision.flybase import (  # noqa: E402
    case_insensitive_current_candidates,
    maps_from_dataframe,
    resolve_symbol,
)
from fl_ai_reagent_stocker.tx_omics_revision.go_analysis import (  # noqa: E402
    parse_chart_records,
    run_david_go,
    verify_go_script_hashes,
)
from fl_ai_reagent_stocker.tx_omics_revision.io_utils import git_blob_sha1, sha256_file  # noqa: E402
from fl_ai_reagent_stocker.tx_omics_revision.keywords import (  # noqa: E402
    match_buckets,
    parse_david_term,
    proposed_decision,
    split_david_genes,
)
from fl_ai_reagent_stocker.tx_omics_revision.overlap import (  # noqa: E402
    exact_intersections,
    membership_table,
    overlapping_genes,
    pairwise_overlap,
)
from fl_ai_reagent_stocker.tx_omics_revision.pathways import (  # noqa: E402
    approve_proposed_terms,
    build_pathway_tables,
    build_term_review,
    reconstruct_summary_ok,
)
from fl_ai_reagent_stocker.tx_omics_revision.plots import plot_all, plot_membership_matrix  # noqa: E402


def _fixture_maps():
    df = pd.DataFrame(
        [
            {
                "primary_FBid": "FBgn0035187",
                "current_symbol": "Trhn",
                "symbol_synonym(s)": "Trh|trh|TRH",
            },
            {
                "primary_FBid": "FBgn0262139",
                "current_symbol": "trh",
                "symbol_synonym(s)": "Trh|trachealess",
            },
            {
                "primary_FBid": "FBgn0001208",
                "current_symbol": "Hn",
                "symbol_synonym(s)": "Trh|TRH",
            },
            {
                "primary_FBid": "FBgn0039509",
                "current_symbol": "bigmax",
                "symbol_synonym(s)": "",
            },
            {
                "primary_FBid": "FBgn0011276",
                "current_symbol": "HLH3B",
                "symbol_synonym(s)": "",
            },
            {
                "primary_FBid": "FBgn0002609",
                "current_symbol": "E(spl)m3-HLH",
                "symbol_synonym(s)": "",
            },
            {
                "primary_FBid": "FBgn0002733",
                "current_symbol": "E(spl)mbeta-HLH",
                "symbol_synonym(s)": "",
            },
            {
                "primary_FBid": "FBgn0000001",
                "current_symbol": "alpha",
                "symbol_synonym(s)": "uniqueSyn",
            },
            {
                "primary_FBid": "FBgn0000002",
                "current_symbol": "beta",
                "symbol_synonym(s)": "sharedSyn",
            },
            {
                "primary_FBid": "FBgn0000003",
                "current_symbol": "gamma",
                "symbol_synonym(s)": "sharedSyn",
            },
            {
                "primary_FBid": "FBgn0261608",
                "current_symbol": "RpL37A",
                "symbol_synonym(s)": "RpL37a",
            },
            {
                "primary_FBid": "FBgn0030616",
                "current_symbol": "RpL37-1",
                "symbol_synonym(s)": "RpL37a",
            },
        ]
    )
    return maps_from_dataframe(df)


class KeywordTests(unittest.TestCase):
    def test_parse_david_term(self):
        self.assertEqual(parse_david_term("GO:0006412~translation"), ("GO:0006412", "translation"))
        self.assertEqual(parse_david_term("dme03010:Ribosome"), ("dme03010", "Ribosome"))

    def test_split_and_match(self):
        self.assertEqual(split_david_genes("FBgn0000001, FBgn0000002"), ["FBgn0000001", "FBgn0000002"])
        self.assertEqual(split_david_genes("FBGN0000001;FBgn0000002"), ["FBgn0000001", "FBgn0000002"])
        buckets, _ = match_buckets("cytoplasmic translation")
        self.assertEqual(buckets, ["Ribosomal/translation"])
        buckets, _ = match_buckets("mitochondrial translation")
        self.assertEqual(buckets, ["Ribosomal/translation", "Mitochondrial/metabolism"])
        self.assertEqual(proposed_decision(buckets), "conflict")
        buckets, _ = match_buckets("defense response to bacterium")
        self.assertEqual(buckets, ["Immune"])
        buckets, _ = match_buckets("unrelated DNA binding")
        self.assertEqual(buckets, [])
        self.assertEqual(proposed_decision(buckets), "exclude")


class FlyBaseResolutionTests(unittest.TestCase):
    def test_trh_family(self):
        maps = _fixture_maps()
        trhn = resolve_symbol("Trhn", maps)
        self.assertEqual(trhn["exact_current_fbgn"], "FBgn0035187")
        self.assertEqual(trhn["match_type"], "exact_current_symbol")
        trh = resolve_symbol("trh", maps)
        self.assertEqual(trh["exact_current_fbgn"], "FBgn0262139")
        bare = resolve_symbol("Trh", maps)
        self.assertEqual(bare["match_type"], "ambiguous_synonym")
        self.assertGreaterEqual(bare["candidate_count"], 3)

    def test_unique_and_unmapped(self):
        maps = _fixture_maps()
        unique = resolve_symbol("uniqueSyn", maps)
        self.assertEqual(unique["match_type"], "unique_synonym")
        self.assertEqual(unique["candidate_fbgn_ids"], ["FBgn0000001"])
        missing = resolve_symbol("notAGene", maps)
        self.assertEqual(missing["match_type"], "unmapped")
        shared = resolve_symbol("sharedSyn", maps)
        self.assertEqual(shared["match_type"], "ambiguous_synonym")
        self.assertEqual(shared["candidate_count"], 2)

    def test_case_insensitive_current_among_synonyms(self):
        maps = _fixture_maps()
        resolved = resolve_symbol("RpL37a", maps)
        self.assertEqual(resolved["match_type"], "ambiguous_synonym")
        hits = case_insensitive_current_candidates(
            "RpL37a", maps, list(resolved["candidate_fbgn_ids"])
        )
        self.assertEqual(hits, ["FBgn0261608"])
        frames = {
            "CSW 4+": pd.DataFrame(
                {"ext_gene": ["RpL37a"], "flybase_gene_id": ["FBgn0261608"]}
            )
        }
        audit, unresolved = audit_category_identities(frames, maps)
        self.assertFalse(unresolved)
        self.assertEqual(
            audit.iloc[0]["match_type"], "case_insensitive_current_among_synonyms"
        )


@unittest.skipUnless(BREAKDOWN_DIR.exists() and TX_INPUT_DIR.exists(), "source gene sets not present")
class SourceCategoryTests(unittest.TestCase):
    def test_master_invariants(self):
        observations = load_master_tables(TX_INPUT_DIR)
        self.assertEqual(len(observations), 902)
        self.assertEqual(observations["flybase_gene_id"].nunique(), 791)
        key = observations[["flybase_gene_id", "threshold", "direction"]]
        self.assertFalse(key.duplicated().any())
        self.assertFalse(key.isna().any().any())
        for meta in MASTER_TABLES.values():
            path = TX_INPUT_DIR / meta["filename"]
            self.assertTrue(path.exists(), path)
            self.assertEqual(len(sha256_file(path)), 64)
        for name in CSW_4PLUS_SOURCES:
            path = BREAKDOWN_DIR / name
            self.assertTrue(path.exists(), path)
            self.assertEqual(len(sha256_file(path)), 64)

    def test_mechanistic_and_homeostatic_and_csw_and_hlh(self):
        mechanistic = build_mechanistic(BREAKDOWN_DIR)
        self.assertEqual(list(zip(mechanistic["ext_gene"], mechanistic["flybase_gene_id"])), MECHANISTIC_IDENTITIES)
        self.assertTrue((mechanistic["gene_set"] == "Mechanistic").all())
        self.assertTrue(mechanistic["gene_set_definition"].str.startswith("Consistent or correlated").all())
        self.assertNotIn("Mechanistic Category", mechanistic.columns)
        self.assertFalse(mechanistic["mechanistic_subcategory"].str.startswith("Homeostatic").any())
        homeostatic = build_homeostatic(BREAKDOWN_DIR)
        self.assertEqual(len(homeostatic), 20)
        self.assertTrue((homeostatic["gene_set"] == "Homeostatic genes").all())
        self.assertTrue(homeostatic["gene_set_definition"].str.contains("opposite direction").all())
        self.assertNotIn(TRACHEALESS_FBGN, set(homeostatic["flybase_gene_id"]))
        self.assertFalse(
            homeostatic.astype(str).apply(lambda col: col.str.contains(TRACHEALESS_FBGN)).any().any()
        )
        self.assertIn(TRHN_FBGN, set(homeostatic["flybase_gene_id"]))
        self.assertEqual(
            homeostatic.loc[homeostatic["flybase_gene_id"] == TRHN_FBGN, "ext_gene"].iloc[0],
            TRHN_SYMBOL,
        )
        self.assertEqual(
            homeostatic.loc[homeostatic["flybase_gene_id"] == TRHN_FBGN, "identity_status"].iloc[0],
            "publication_resolved",
        )
        csw = build_csw_4plus(BREAKDOWN_DIR)
        self.assertEqual(len(csw), 97)
        self.assertTrue((csw["gene_set"] == "CSW 4+ genes").all())
        self.assertTrue(csw["gene_set_definition"].str.contains("four or more").all())
        self.assertTrue((csw["qualifies_FC0.5_4plus"] == "TRUE").all())
        self.assertEqual((csw["qualifies_FC1_4plus"] == "TRUE").sum(), 7)
        self.assertNotIn("frequency", csw.columns)
        self.assertNotIn("wake_corr_exps", csw.columns)
        self.assertIn("frequency_FC0.5_wake", csw.columns)
        self.assertIn("frequency_FC1_wake", csw.columns)
        self.assertIn("wake_corr_exps_FC0.5_wake", csw.columns)
        self.assertIn("wake_corr_exps_FC1_wake", csw.columns)
        fc1 = pd.read_csv(BREAKDOWN_DIR / "FC1_Wake_frequency_4_n=7genes.csv")
        self.assertTrue(set(fc1["flybase_gene_id"]).issubset(set(csw["flybase_gene_id"])))
        hlh = build_hlh(_fixture_maps())
        self.assertEqual(list(zip(hlh["ext_gene"], hlh["flybase_gene_id"])), HLH_PUBLICATION_GENES)
        self.assertTrue((hlh["gene_set"] == "HLH genes").all())
        self.assertTrue(hlh["gene_set_definition"].str.contains("upstream regulators").all())
        self.assertEqual(set(hlh.columns) & {"frequency", "wake_corr_exps"}, set())
        self.assertTrue((hlh["publication_section"] == "Results").all())
        self.assertTrue((hlh["publication_table"] == "Table S1").all())
        self.assertTrue((hlh["motif"] == "CAGCTG E-box").all())
        for frame in (mechanistic, homeostatic, csw, hlh):
            self.assertNotIn(TRACHEALESS_FBGN, set(frame["flybase_gene_id"]))

    def test_category_audit_accepts_publication_resolved_trhn(self):
        maps = _fixture_maps()
        homeostatic = build_homeostatic(BREAKDOWN_DIR)
        trhn_row = homeostatic[homeostatic["flybase_gene_id"] == TRHN_FBGN]
        frames = {
            "Homeostatic History/Rebound": trhn_row,
            "HLH Upstream Regulators": build_hlh(maps),
        }
        audit, unresolved = audit_category_identities(frames, maps)
        self.assertFalse(unresolved)
        trhn = audit[audit["ext_gene"] == "Trhn"].iloc[0]
        self.assertEqual(trhn["flybase_gene_id"], TRHN_FBGN)


class GoParseTests(unittest.TestCase):
    def test_token_count_and_fdr(self):
        go_df = pd.DataFrame(
            [
                {
                    "Category": "GOTERM_BP_DIRECT",
                    "Term": "GO:0006412~translation",
                    "Count": 2,
                    "%": "10",
                    "Pvalue": 0.001,
                    "Genes": "FBgn0000001,FBgn0000002",
                    "List Total": 10,
                    "Pop Hits": 20,
                    "Pop Total": 100,
                    "Fold Enrichment": 5,
                    "Bonferroni": 0.01,
                    "Benjamini": 0.01,
                    "FDR": 1.5,
                }
            ]
        )
        terms, genes = parse_chart_records(go_df, "Wake FC0.5", "wake.xlsx")
        self.assertEqual(len(terms), 1)
        self.assertEqual(len(genes), 2)
        self.assertAlmostEqual(float(terms.iloc[0]["fdr_q"]), 0.015)
        with self.assertRaises(ValueError):
            parse_chart_records(
                go_df.assign(Count=3),
                "Wake FC0.5",
                "wake.xlsx",
            )
        empty_terms, empty_genes = parse_chart_records(
            go_df.assign(FDR=20.0),
            "Sleep FC1",
            "sleep.xlsx",
        )
        self.assertTrue(empty_terms.empty)
        self.assertTrue(empty_genes.empty)
        normalized = go_df.copy()
        normalized.loc[0, "Genes"] = "FBGN0000001,FBgn0000002"
        terms, genes = parse_chart_records(normalized, "Wake FC0.5", "wake.xlsx")
        self.assertEqual(list(genes["flybase_gene_id"]), ["FBgn0000001", "FBgn0000002"])
        self.assertTrue((terms["fdr_percent"] <= 10).all())


class GoHashTests(unittest.TestCase):
    def test_stale_or_wrong_blob_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in GO_SCRIPT_BLOBS:
                (root / name).write_text("not the pinned script\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_go_script_hashes(root)

    @unittest.skipUnless(GO_ANALYSIS_DIR.exists(), "GO_Analysis scripts not present")
    def test_pinned_local_blobs(self):
        checked = verify_go_script_hashes(GO_ANALYSIS_DIR)
        self.assertEqual(checked, GO_SCRIPT_BLOBS)
        for name, expected in GO_SCRIPT_BLOBS.items():
            self.assertEqual(git_blob_sha1(GO_ANALYSIS_DIR / name), expected)


class OfficialGoAnalysisTests(unittest.TestCase):
    def test_wrapper_does_not_reimplement_david_soap(self):
        source = (
            REPO_ROOT / "fl_ai_reagent_stocker" / "tx_omics_revision" / "go_analysis.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("getChartReport", source)
        self.assertNotIn("import ProcessGOresults", source)
        self.assertIn("init_david_client", source)
        self.assertIn("fetch_go_report", source)

    def test_run_david_go_calls_generate_go_chart_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "input"
            raw_dir = root / "raw"
            processed_dir = root / "processed"
            evidence_dir = root / "evidence"
            input_dir.mkdir()
            evidence_dir.mkdir()
            submitted = {}
            table_ids = {
                "Wake FC0.5": "FBgn0000001",
                "Wake FC1": "FBgn0000002",
                "Sleep FC0.5": "FBgn0000003",
                "Sleep FC1": "FBgn0000004",
            }
            for name, meta in MASTER_TABLES.items():
                fbgn = table_ids[name]
                pd.DataFrame({"ext_gene": ["g"], "flybase_gene_id": [fbgn]}).to_csv(
                    input_dir / meta["filename"], index=False
                )
                submitted[name] = {fbgn}

            fake = MagicMock()
            fake.init_david_client.return_value = MagicMock(name="david_client")

            def load_input_table(path):
                return pd.read_csv(path)["flybase_gene_id"].to_numpy()

            def fetch_go_report(_client, flybase_ids, list_name):
                if "Sleep gene list FC1" in list_name:
                    return pd.DataFrame()
                gene = str(flybase_ids[0])
                return pd.DataFrame(
                    [
                        {
                            "Category": "GOTERM_BP_DIRECT",
                            "Term": "GO:0006412~translation",
                            "Count": 1,
                            "%": "10",
                            "Pvalue": 0.001,
                            "Genes": gene,
                            "List Total": 10,
                            "Pop Hits": 20,
                            "Pop Total": 100,
                            "Fold Enrichment": 5,
                            "Bonferroni": 0.01,
                            "Benjamini": 0.01,
                            "FDR": 1.0,
                        }
                    ]
                )

            fake.load_input_table.side_effect = load_input_table
            fake.fetch_go_report.side_effect = fetch_go_report
            with patch(
                "fl_ai_reagent_stocker.tx_omics_revision.go_analysis.verify_go_script_hashes",
                return_value=GO_SCRIPT_BLOBS,
            ), patch(
                "fl_ai_reagent_stocker.tx_omics_revision.go_analysis._load_go_module",
                return_value=fake,
            ):
                payload = run_david_go(
                    input_dir=input_dir,
                    raw_dir=raw_dir,
                    processed_dir=processed_dir,
                    evidence_dir=evidence_dir,
                    submitted_ids=submitted,
                    go_dir=root,
                )
            fake.init_david_client.assert_called_once()
            self.assertEqual(fake.load_input_table.call_count, 4)
            self.assertEqual(fake.fetch_go_report.call_count, 4)
            self.assertEqual(payload["go_analysis_repo"], "https://github.com/aadish98/GO_Analysis")
            self.assertEqual(
                payload["skipped_go_analysis_step"], "ProcessGOresults.process_csv_files"
            )
            statuses = {row["source_table"]: row["status"] for row in payload["outcomes"]}
            self.assertEqual(statuses["Sleep FC1"], "no_fdr_passing_terms")
            self.assertEqual(statuses["Wake FC0.5"], "success")


class PathwayAndOverlapTests(unittest.TestCase):
    def test_review_conflict_and_reconstruction(self):
        terms = pd.DataFrame(
            [
                {
                    "source_table": "Wake FC0.5",
                    "source_workbook": "w.xlsx",
                    "category": "GOTERM_BP_DIRECT",
                    "term_id": "GO:1",
                    "term_name": "cytoplasmic translation",
                    "count": 1,
                    "percent": "1",
                    "pvalue": 0.01,
                    "list_total": 10,
                    "pop_hits": 10,
                    "pop_total": 100,
                    "fold_enrichment": 2,
                    "bonferroni": 0.1,
                    "benjamini": 0.1,
                    "fdr_percent": 1.0,
                    "fdr_q": 0.01,
                    "genes": "FBgn0000001",
                },
                {
                    "source_table": "Wake FC0.5",
                    "source_workbook": "w.xlsx",
                    "category": "GOTERM_BP_DIRECT",
                    "term_id": "GO:2",
                    "term_name": "mitochondrial translation",
                    "count": 1,
                    "percent": "1",
                    "pvalue": 0.01,
                    "list_total": 10,
                    "pop_hits": 10,
                    "pop_total": 100,
                    "fold_enrichment": 2,
                    "bonferroni": 0.1,
                    "benjamini": 0.1,
                    "fdr_percent": 2.0,
                    "fdr_q": 0.02,
                    "genes": "FBgn0000002",
                },
            ]
        )
        genes = pd.DataFrame(
            [
                {
                    "source_table": "Wake FC0.5",
                    "category": "GOTERM_BP_DIRECT",
                    "term_id": "GO:1",
                    "term_name": "cytoplasmic translation",
                    "flybase_gene_id": "FBgn0000001",
                },
                {
                    "source_table": "Wake FC0.5",
                    "category": "GOTERM_BP_DIRECT",
                    "term_id": "GO:2",
                    "term_name": "mitochondrial translation",
                    "flybase_gene_id": "FBgn0000002",
                },
            ]
        )
        review = build_term_review(terms, genes)
        self.assertIn("conflict", set(review["proposed_decision"]))
        approved = approve_proposed_terms(review, "include_all_matched_buckets")
        observations = pd.DataFrame(
            [
                {
                    "flybase_gene_id": "FBgn0000001",
                    "ext_gene": "geneA",
                    "source_table": "Wake FC0.5",
                    "threshold": "FC0.5",
                    "direction": "wake",
                    "frequency": "4",
                    "is_cycling": "FALSE",
                    "ZT_cyclers": "",
                    "sleep_corr_exps": "",
                    "wake_corr_exps": "MechSD6",
                },
                {
                    "flybase_gene_id": "FBgn0000002",
                    "ext_gene": "geneB",
                    "source_table": "Wake FC0.5",
                    "threshold": "FC0.5",
                    "direction": "wake",
                    "frequency": "5",
                    "is_cycling": "FALSE",
                    "ZT_cyclers": "",
                    "sleep_corr_exps": "",
                    "wake_corr_exps": "MechSD3",
                },
            ]
        )
        category_ids = {
            "CSW 4+": {"FBgn0000001"},
            "Mechanistic": set(),
            "Homeostatic History/Rebound": set(),
            "HLH Upstream Regulators": set(),
        }
        tables = build_pathway_tables(
            approved, genes, observations, category_ids, {"FBgn0000001": "FC0.5"}
        )
        self.assertTrue(
            reconstruct_summary_ok(
                tables["Ribosomal/translation"], genes, approved
            )
        )
        self.assertIn("FBgn0000002", set(tables["Ribosomal/translation"]["flybase_gene_id"]))
        self.assertIn("FBgn0000002", set(tables["Mitochondrial/metabolism"]["flybase_gene_id"]))

    def test_membership_invariants(self):
        sets = {
            "Mechanistic": pd.DataFrame({"ext_gene": ["a"], "flybase_gene_id": ["FBgn0000001"]}),
            "Homeostatic History/Rebound": pd.DataFrame(
                {"ext_gene": ["a", "b"], "flybase_gene_id": ["FBgn0000001", "FBgn0000002"]}
            ),
            "CSW 4+": pd.DataFrame({"ext_gene": ["c"], "flybase_gene_id": ["FBgn0000003"]}),
            "HLH Upstream Regulators": pd.DataFrame(
                {"ext_gene": ["d"], "flybase_gene_id": ["FBgn0000004"]}
            ),
            "CSW Ribosomal/Translation": pd.DataFrame(
                {"ext_gene": ["a"], "flybase_gene_id": ["FBgn0000001"]}
            ),
            "CSW Mitochondrial/Metabolism": pd.DataFrame(
                {"ext_gene": ["c"], "flybase_gene_id": ["FBgn0000003"]}
            ),
            "CSW Immune": pd.DataFrame(columns=["ext_gene", "flybase_gene_id"]),
        }
        membership = membership_table(sets)
        self.assertEqual(len(membership), 4)
        overlapping = overlapping_genes(membership)
        self.assertEqual(len(overlapping), 2)
        self.assertEqual(set(overlapping["primary_symbol"]), {"a", "c"})
        exact, exact_genes = exact_intersections(membership)
        self.assertEqual(int(exact["intersection_size"].sum()), 4)
        self.assertEqual(len(exact_genes), 4)
        pairs, _ = pairwise_overlap(membership)
        ab = pairs[(pairs["set_a"] == "Mechanistic") & (pairs["set_b"] == "Homeostatic History/Rebound")]
        ba = pairs[(pairs["set_a"] == "Homeostatic History/Rebound") & (pairs["set_b"] == "Mechanistic")]
        self.assertEqual(int(ab.iloc[0]["n_overlapping_genes"]), int(ba.iloc[0]["n_overlapping_genes"]))
        reconstructed = int(
            ((membership["Mechanistic"] == "TRUE") & (membership["Homeostatic History/Rebound"] == "TRUE")).sum()
        )
        self.assertEqual(int(ab.iloc[0]["n_overlapping_genes"]), reconstructed)
        self.assertTrue((overlapping["membership_count"] >= 2).all())

    def test_pending_conflict_without_all_bucket_policy(self):
        review = pd.DataFrame(
            [
                {
                    "proposed_decision": "conflict",
                    "candidate_buckets": "Ribosomal/translation;Mitochondrial/metabolism",
                }
            ]
        )
        pending = approve_proposed_terms(review, "wait")
        self.assertEqual(pending.iloc[0]["approval_status"], "pending_conflict")

    def test_empty_plots(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            empty = pd.DataFrame(columns=["ext_gene", "flybase_gene_id", "directions"])
            membership = membership_table({name: empty.copy() for name in SET_NAMES})
            plot_all(
                {
                    "Ribosomal/translation": empty,
                    "Mitochondrial/metabolism": empty,
                    "Immune": empty,
                },
                pd.DataFrame(),
                membership,
                overlapping_genes(membership),
                pairwise_overlap(membership)[0],
                root / "figures",
                root / "figure_data",
            )
            self.assertTrue((root / "figures" / "pathway_overlap_venn.png").exists())
            self.assertTrue((root / "figure_data" / "pathway_gene_counts.csv").exists())

    def test_singleton_identical_disjoint_plots_and_labels(self):
        empty = pd.DataFrame(columns=["ext_gene", "flybase_gene_id", "directions"])
        gene_a = pd.DataFrame({"ext_gene": ["Trhn"], "flybase_gene_id": ["FBgn0035187"], "directions": ["wake"]})
        gene_b = pd.DataFrame({"ext_gene": ["RpL23"], "flybase_gene_id": ["FBgn0010078"], "directions": ["wake"]})
        fixtures = {
            "singleton": {
                "Mechanistic": gene_a,
                "Homeostatic History/Rebound": empty,
                "CSW 4+": empty,
                "HLH Upstream Regulators": empty,
                "CSW Ribosomal/Translation": empty,
                "CSW Mitochondrial/Metabolism": empty,
                "CSW Immune": empty,
            },
            "identical": {
                "Mechanistic": gene_a,
                "Homeostatic History/Rebound": gene_a,
                "CSW 4+": empty,
                "HLH Upstream Regulators": empty,
                "CSW Ribosomal/Translation": gene_a,
                "CSW Mitochondrial/Metabolism": gene_a,
                "CSW Immune": empty,
            },
            "disjoint": {
                "Mechanistic": gene_a,
                "Homeostatic History/Rebound": gene_b,
                "CSW 4+": empty,
                "HLH Upstream Regulators": empty,
                "CSW Ribosomal/Translation": gene_a,
                "CSW Mitochondrial/Metabolism": gene_b,
                "CSW Immune": empty,
            },
        }
        for name, sets in fixtures.items():
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                membership = membership_table(sets)
                overlapping = overlapping_genes(membership)
                plot_all(
                    {
                        "Ribosomal/translation": sets["CSW Ribosomal/Translation"],
                        "Mitochondrial/metabolism": sets["CSW Mitochondrial/Metabolism"],
                        "Immune": sets["CSW Immune"],
                    },
                    pd.DataFrame(),
                    membership,
                    overlapping,
                    pairwise_overlap(membership)[0],
                    root / "figures",
                    root / "figure_data",
                )
                self.assertTrue((root / "figures" / "seven_set_upset.png").exists(), name)
                if not overlapping.empty:
                    plot_membership_matrix(
                        overlapping, root / "figures" / "matrix.png", root / "figure_data"
                    )
                    labels = pd.read_csv(root / "figure_data" / "overlapping_gene_membership_matrix.csv")
                    self.assertEqual(
                        set(labels["primary_symbol"]),
                        set(overlapping["primary_symbol"]),
                    )


@unittest.skipUnless(STOCKER_DIR.exists(), "generated stocker gene sets not present")
class GeneratedOutputTests(unittest.TestCase):
    def test_stocker_discovery_and_columns(self):
        found = _discover_input_csvs(STOCKER_DIR)
        self.assertEqual(len(found), 7)
        self.assertEqual(
            [path.name for path in found],
            [
                "01_Mechanistic_n=6genes.csv",
                "02_Homeostatic_HistoryxRebound_n=20genes.csv",
                "03_CSW_4plus_n=97genes.csv",
                "04_HLH_upstream_regulators_n=4genes.csv",
                "05_CSW_Ribosomal_Translation_n=99genes.csv",
                "06_CSW_Mitochondrial_Metabolism_n=167genes.csv",
                "07_CSW_Immune_n=5genes.csv",
            ],
        )
        self.assertEqual(_validate_gene_columns(found, "flybase_gene_id", "ext_gene"), [])
        for path in found:
            df = pd.read_csv(path)
            self.assertEqual(df["flybase_gene_id"].nunique(), len(df))
            self.assertFalse(df["flybase_gene_id"].isna().any())
            self.assertNotIn(TRACHEALESS_FBGN, set(df["flybase_gene_id"]))

    def test_pathway_non_id_metadata_restored(self):
        paths = audit_paths()
        identity_cols = {"flybase_gene_id", "primary_symbol", "corrected_ext_gene"}
        for pre, approved in zip(
            sorted(paths["pathways_preaudit"].glob("*.csv")),
            sorted(paths["pathways_approved"].glob("*.csv")),
        ):
            self.assertEqual(pre.name, approved.name)
            before = pd.read_csv(pre, dtype=str, keep_default_na=False)
            after = pd.read_csv(approved, dtype=str, keep_default_na=False)
            self.assertEqual(len(before), len(after))
            for col in before.columns:
                if col in identity_cols:
                    continue
                self.assertEqual(list(before[col]), list(after[col]), f"{pre.name}:{col}")

    def test_overlap_report_names_priority_genes(self):
        paths = audit_paths()
        report = paths["overlap_report"].read_text(encoding="utf-8")
        for token in ("Trhn", "FBgn0035187", "AstA-R2", "RpL23", "bigmax"):
            self.assertIn(token, report)
        membership = pd.read_csv(paths["evidence"] / "Set_Membership.csv")
        overlapping = pd.read_csv(paths["evidence"] / "Overlapping_Genes.csv")
        self.assertEqual(len(overlapping), int((membership["membership_count"] >= 2).sum()))
        figure_labels = pd.read_csv(paths["figure_data"] / "overlapping_gene_membership_matrix.csv")
        self.assertEqual(set(figure_labels["primary_symbol"]), set(overlapping["primary_symbol"]))


if __name__ == "__main__":
    unittest.main()
