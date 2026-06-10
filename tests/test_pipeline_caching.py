from __future__ import annotations

import os
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker import cli  # noqa: E402
from fl_ai_reagent_stocker.integrations import pubmed  # noqa: E402
from fl_ai_reagent_stocker.integrations.pubmed import PubMedCache  # noqa: E402
from fl_ai_reagent_stocker.pipelines import stock_splitting as ss  # noqa: E402


def _write_alias_fixture(root: Path) -> None:
    (root / "fbgn_annotation_ID_test.tsv").write_text(
        "\n".join(
            [
                "##primary_FBgn#\tsecondary_FBgn#(s)\torganism_abbreviation",
                "FBgn0000001\tFBgn0000002,FBgn0000003\tDmel",
            ]
        ),
        encoding="utf-8",
    )


def _write_synonym_fixture(root: Path) -> None:
    genes_dir = root / "genes"
    genes_dir.mkdir(parents=True, exist_ok=True)
    (genes_dir / "fb_synonym_test.tsv").write_text(
        "\n".join(
            [
                "organism_abbreviation\tcurrent_symbol\tsymbol_synonym(s)\tprimary_FBid",
                "Dmel\tGeneA\tOldGeneA|AliasA\tFBgn0000001",
            ]
        ),
        encoding="utf-8",
    )


def _write_rnai_fixture(root: Path) -> None:
    constructs_dir = root / "transgenic_constructs"
    constructs_dir.mkdir(parents=True, exist_ok=True)
    (constructs_dir / "transgenic_construct_descriptions_test.tsv").write_text(
        "\n".join(
            [
                "Component Allele (id)\tTransgenic Construct (symbol)\t"
                "Transgenic Construct (id)\tDescription (text)",
                "FBal0001\tUAS-GeneA-RNAi\tFBtp0001\tdouble stranded RNA construct",
            ]
        ),
        encoding="utf-8",
    )


class TestFlyBaseAuxCaching(unittest.TestCase):
    def setUp(self) -> None:
        ss._load_secondary_fbgn_to_primary_cached.cache_clear()
        ss._load_gene_synonyms_map_cached.cache_clear()
        ss._load_rnai_type_lookup_cached.cache_clear()

    def test_secondary_fbgn_alias_loader_reads_once_for_same_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_alias_fixture(root)

            with patch("builtins.open", wraps=open) as mocked_open:
                first = ss._load_secondary_fbgn_to_primary(root)
                second = ss._load_secondary_fbgn_to_primary(root)

            self.assertEqual(first, {"FBgn0000002": "FBgn0000001", "FBgn0000003": "FBgn0000001"})
            self.assertEqual(second, first)
            self.assertEqual(mocked_open.call_count, 1)

    def test_gene_synonym_loader_reads_once_for_same_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_synonym_fixture(root)

            real_loader = ss.load_flybase_tsv
            with patch.object(ss, "load_flybase_tsv", wraps=real_loader) as mocked_loader:
                first = ss._load_gene_synonyms_map(root)
                second = ss._load_gene_synonyms_map(root)

            self.assertEqual(first["GeneA"], "AliasA; OldGeneA")
            self.assertEqual(second, first)
            self.assertEqual(mocked_loader.call_count, 1)

    def test_rnai_type_loader_reads_once_for_same_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_rnai_fixture(root)

            real_loader = ss.load_flybase_tsv
            with patch.object(ss, "load_flybase_tsv", wraps=real_loader) as mocked_loader:
                first = ss._load_rnai_type_lookup(root)
                second = ss._load_rnai_type_lookup(root)

            self.assertIn("FBal0001", first[0])
            self.assertIn("FBtp0001", first[1])
            self.assertEqual(second, first)
            self.assertEqual(mocked_loader.call_count, 1)


class TestPubMedCacheReuse(unittest.TestCase):
    def setUp(self) -> None:
        pubmed._PARSED_PUBMED_CACHE_BY_PATH_MTIME.clear()

    def test_cache_parse_reused_until_file_mtime_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "pubmed_cache.csv"
            cache_path.write_text(
                "pmid,title,abstract\n123,Title one,Abstract one\n",
                encoding="utf-8",
            )

            with patch.object(pubmed.pd, "read_csv", wraps=pd.read_csv) as mocked_read:
                first = PubMedCache(cache_path)
                first.load()
                first._cache["999"] = {"title": "local only"}

                second = PubMedCache(cache_path)
                second.load()

                self.assertEqual(mocked_read.call_count, 1)
                self.assertIn("123", second._cache)
                self.assertNotIn("999", second._cache)

                cache_path.write_text(
                    "pmid,title,abstract\n123,Title one,Abstract one\n456,Title two,Abstract two\n",
                    encoding="utf-8",
                )
                future = time.time() + 2
                os.utime(cache_path, (future, future))

                third = PubMedCache(cache_path)
                third.load()

                self.assertEqual(mocked_read.call_count, 2)
                self.assertIn("456", third._cache)


class TestValidateSkipGuard(unittest.TestCase):
    def test_validation_runs_only_when_openai_key_present(self):
        self.assertFalse(cli._should_run_validation(None))
        self.assertFalse(cli._should_run_validation(""))
        self.assertFalse(cli._should_run_validation("   "))
        self.assertTrue(cli._should_run_validation("sk-test"))
        self.assertFalse(
            cli._should_run_validation("sk-test", validation_enabled=False)
        )


if __name__ == "__main__":
    unittest.main()
