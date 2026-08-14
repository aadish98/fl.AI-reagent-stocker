#!/usr/bin/env python3
"""Run hash-pinned aadish98/GO_Analysis GenerateGOChartReport DAVID enrichment."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker.tx_omics_revision.pipeline import run_pipeline


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-root", type=Path, default=None)
    args = parser.parse_args(argv)
    run_pipeline(through="go", audit_root=args.audit_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
