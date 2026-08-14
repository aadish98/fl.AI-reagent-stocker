#!/usr/bin/env python3
"""Build Tx-Omics Revision gene sets, audits, overlap tables, and figures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker.tx_omics_revision.pipeline import STAGES, run_pipeline


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--through", default="all", choices=["all", *STAGES])
    parser.add_argument(
        "--approve-proposed-terms",
        action="store_true",
        help="Approve proposed keyword includes; conflict terms go into every matched bucket.",
    )
    parser.add_argument(
        "--promote-if-audit-clean",
        action="store_true",
        help="Promote pathway CSVs and copy the seven stocker inputs if the FlyBase audit is clean.",
    )
    parser.add_argument("--audit-root", type=Path, default=None)
    args = parser.parse_args(argv)
    run_pipeline(
        through=args.through,
        approve_proposed_terms=args.approve_proposed_terms,
        promote_if_audit_clean=args.promote_if_audit_clean,
        audit_root=args.audit_root,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
