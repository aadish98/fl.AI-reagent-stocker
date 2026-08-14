#!/usr/bin/env python3
"""Regenerate Tx-Omics Revision figures from saved evidence tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker.tx_omics_revision.constants import PATHWAY_BUCKETS, audit_paths
from fl_ai_reagent_stocker.tx_omics_revision.io_utils import read_csv
from fl_ai_reagent_stocker.tx_omics_revision.plots import plot_all

_STEMS = {
    "Ribosomal/translation": "05_CSW_Ribosomal_Translation",
    "Mitochondrial/metabolism": "06_CSW_Mitochondrial_Metabolism",
    "Immune": "07_CSW_Immune",
}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-root", type=Path, default=None)
    args = parser.parse_args(argv)
    paths = audit_paths(args.audit_root)
    membership = read_csv(paths["evidence"] / "Set_Membership.csv")
    overlapping = read_csv(paths["evidence"] / "Overlapping_Genes.csv")
    pairs = read_csv(paths["evidence"] / "Pairwise_Overlap.csv")
    approved = read_csv(paths["go_review"] / "GO_Term_Bucket_Review.approved.csv")
    pathway_sets = {}
    for bucket, stem in _STEMS.items():
        matches = list(paths["pathways_approved"].glob(f"{stem}_*.csv"))
        pathway_sets[bucket] = read_csv(matches[0]) if matches else pd.DataFrame()
    plot_all(
        pathway_sets,
        approved,
        membership,
        overlapping,
        pairs,
        paths["figures"],
        paths["figure_data"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
