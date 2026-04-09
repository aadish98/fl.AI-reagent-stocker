"""Golden checks for Stock Phenotype / Chado audit (plan: optional-golden-tests)."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _load_build_fbst_module():
    path = REPO_ROOT / "scripts" / "build_fbst_derived_stock_components.py"
    spec = importlib.util.spec_from_file_location("fbst_build_test", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


from fl_ai_reagent_stocker.pipelines.stock_splitting import _extract_flybase_ids  # noqa: E402
from scripts.audit_stock_phenotype_pipeline import CANONICAL_FB_ID_RE  # noqa: E402


class TestChadoTokenPattern(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.build = _load_build_fbst_module()

    def test_token_pattern_extracts_fbal(self):
        g = "@FBal0018186:w<up>1118</up>@"
        matches = self.build.TOKEN_PATTERN.findall(g)
        self.assertEqual([m[0] for m in matches], ["FBal0018186"])

    def test_token_pattern_multiple(self):
        g = (
            "@FBal0018607:y<up>1</up>@ w<up>*</up>; "
            "@FBti0210151:PBac{UAS-hFBXL20.HA}VK00033@"
        )
        ids = [m[0] for m in self.build.TOKEN_PATTERN.findall(g)]
        self.assertEqual(ids, ["FBal0018607", "FBti0210151"])


class TestCanonicalFlyBaseIdRegex(unittest.TestCase):
    def test_ignores_human_fbxo_style(self):
        s = "P{GawB}FBXO11<up>NP2786</up>"
        self.assertEqual(CANONICAL_FB_ID_RE.findall(s), [])

    def test_finds_fbti_not_human_symbol_in_construct_label(self):
        s = "@FBti0039590:P{EPgy2}FBXO11<up>EY09314</up>@"
        self.assertEqual(CANONICAL_FB_ID_RE.findall(s), ["FBti0039590"])


class TestExtractFlybaseIds(unittest.TestCase):
    def test_genotype_phenotype_style(self):
        self.assertEqual(
            _extract_flybase_ids("106y[106y] FBal0151008"),
            ["FBal0151008"],
        )

    def test_dedup_order(self):
        out = _extract_flybase_ids("FBal0151008 / FBal0151008")
        self.assertEqual(out, ["FBal0151008"])


if __name__ == "__main__":
    unittest.main()
