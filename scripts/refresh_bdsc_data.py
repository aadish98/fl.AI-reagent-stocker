#!/usr/bin/env python3
"""Download the native BDSC datasets used by the GAL4 workbook exporter."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path
from typing import Dict, List, Optional
from urllib.request import Request, urlopen


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "bdsc"
BDSC_DOWNLOAD_BASE = "https://bdsc.indiana.edu/pdf"
DATASETS: Dict[str, List[str]] = {
    "bloomington.csv": ["Stk #", "Genotype", "A.K.A"],
    "stockgenes.csv": [
        "stknum",
        "genotype",
        "component_symbol",
        "gene_symbol",
        "bdsc_symbol_id",
    ],
}
REQUEST_HEADERS = {
    "User-Agent": "fl.AI-reagent-stocker BDSC data refresher",
    "Accept": "text/csv,*/*",
}


def _validate_csv(path: Path, required_columns: List[str]) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        columns = next(csv.reader(handle), [])
    missing = [column for column in required_columns if column not in columns]
    if missing:
        raise ValueError(
            f"{path.name} is missing required column(s): {', '.join(missing)}"
        )


def _download_dataset(filename: str, destination_dir: Path) -> Path:
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / filename
    temporary = destination.with_suffix(f"{destination.suffix}.tmp")
    request = Request(
        f"{BDSC_DOWNLOAD_BASE}/{filename}",
        headers=REQUEST_HEADERS,
    )
    try:
        with urlopen(request, timeout=120) as response, temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        if temporary.stat().st_size == 0:
            raise ValueError(f"Downloaded an empty file for {filename}")
        _validate_csv(temporary, DATASETS[filename])
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def refresh_bdsc_data(destination_dir: Path = DEFAULT_OUTPUT_DIR) -> List[Path]:
    return [
        _download_dataset(filename, destination_dir)
        for filename in DATASETS
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download current native BDSC stock datasets.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination directory (default: data/bdsc).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    paths = refresh_bdsc_data(args.output_dir.resolve())
    for path in paths:
        print(f"Downloaded {path} ({path.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
