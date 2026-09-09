"""Guard against the default split-config pointing at a file that does not ship.

Only ``data/config/stock_split_config_priority_example.json`` is git-tracked
(see ``.gitignore``); the other example configs are local-only. If
``Settings.split_config_path`` (or its underlying ``DEFAULT_SPLIT_CONFIG_PATH``)
ever drifts back to one of those gitignored files, a fresh clone running
``python -m fl_ai_reagent_stocker run`` without ``--config`` would fail with a
missing-file error. This test catches that regression.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker.config import DEFAULT_SPLIT_CONFIG_PATH, Settings  # noqa: E402
from fl_ai_reagent_stocker.pipelines.stock_splitting import load_split_config  # noqa: E402


class DefaultSplitConfigPathTests(unittest.TestCase):
    def test_default_split_config_path_exists(self) -> None:
        self.assertTrue(
            DEFAULT_SPLIT_CONFIG_PATH.exists(),
            f"DEFAULT_SPLIT_CONFIG_PATH does not exist on disk: "
            f"{DEFAULT_SPLIT_CONFIG_PATH}. Only "
            f"stock_split_config_priority_example.json is git-tracked; "
            f"the default must point at a file that ships with a fresh clone.",
        )

    def test_default_split_config_path_is_the_tracked_priority_config(self) -> None:
        self.assertEqual(
            DEFAULT_SPLIT_CONFIG_PATH.name,
            "stock_split_config_priority_example.json",
        )

    def test_settings_default_split_config_path_exists(self) -> None:
        settings = Settings()
        self.assertTrue(
            settings.split_config_path.exists(),
            f"Settings().split_config_path does not exist on disk: "
            f"{settings.split_config_path}",
        )

    def test_public_priority_config_is_a_complete_documented_catalog(self) -> None:
        config = load_split_config(DEFAULT_SPLIT_CONFIG_PATH)
        filters = config["filters"]

        self.assertEqual(len(config["combination_names"]), len(config["combinations"]))
        self.assertTrue(
            all(filter_name in filters for combo in config["combinations"] for filter_name in combo)
        )
        self.assertEqual(set(filters), set(config["filterDescriptions"]))
        self.assertIn(
            ["Bloomington", "RNAi_reagent", "Ref++"],
            config["combinations"],
        )
        self.assertIn(
            ["Bloomington", "RNAi_reagent", "Ref+"],
            config["combinations"],
        )


if __name__ == "__main__":
    unittest.main()
