#!/usr/bin/env python3
"""
Extract FBtp -> FBti mappings from FlyBase chado_FBti XML.

The XML is streamed line-by-line so the full file is never loaded into memory.
Each FBti feature can reference one or more FBtp construct IDs via
feature_relationship/object_id/feature/uniquename.

Usage:
    python scripts/build_fbtp_to_fbti_mapping.py
    python scripts/build_fbtp_to_fbti_mapping.py --dry-run
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Iterator, Set, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
TRANSGENIC_INSERTIONS_DIR = REPO_ROOT / "data" / "flybase" / "transgenic_insertions"
DEFAULT_INPUT_PATH = TRANSGENIC_INSERTIONS_DIR / "chado_FBti.xml.gz"
DEFAULT_OUTPUT_PATH = TRANSGENIC_INSERTIONS_DIR / "fbtp_to_fbti.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract FBtp-to-FBti mappings from chado_FBti XML.")
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help=f"Path to chado_FBti XML(.gz). Default: {DEFAULT_INPUT_PATH}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"CSV output path. Default: {DEFAULT_OUTPUT_PATH}",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse and report counts without writing the CSV.",
    )
    return parser.parse_args()


def _extract_xml_text(line: str, tag: str) -> str:
    start_token = f"<{tag}>"
    end_token = f"</{tag}>"
    start_idx = line.find(start_token)
    end_idx = line.find(end_token)
    if start_idx == -1 or end_idx == -1 or end_idx < start_idx:
        return ""
    return line[start_idx + len(start_token):end_idx].strip()


def iter_fbtp_fbti_pairs(xml_path: Path) -> Iterator[Tuple[str, str]]:
    current_fbti = ""
    feature_relationship_depth = 0
    object_id_depth = 0
    relationship_fbtp = ""

    with gzip.open(xml_path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if "<accession>FBti" in line:
                current_fbti = _extract_xml_text(line, "accession")
                continue

            if "<feature_relationship>" in line:
                feature_relationship_depth += 1
                if feature_relationship_depth == 1:
                    relationship_fbtp = ""
                continue

            if "</feature_relationship>" in line:
                if feature_relationship_depth == 1 and current_fbti and relationship_fbtp:
                    yield relationship_fbtp, current_fbti
                if feature_relationship_depth > 0:
                    feature_relationship_depth -= 1
                if feature_relationship_depth == 0:
                    object_id_depth = 0
                    relationship_fbtp = ""
                continue

            if feature_relationship_depth > 0 and "<object_id>" in line:
                object_id_depth += 1
                continue

            if feature_relationship_depth > 0 and "</object_id>" in line:
                if object_id_depth > 0:
                    object_id_depth -= 1
                continue

            if (
                feature_relationship_depth == 1
                and object_id_depth == 1
                and "<uniquename>FBtp" in line
            ):
                relationship_fbtp = _extract_xml_text(line, "uniquename")


def write_pairs(output_path: Path, pairs: Set[Tuple[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["FBtp", "FBti"])
        for fbtp_id, fbti_id in sorted(pairs):
            writer.writerow([fbtp_id, fbti_id])


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input XML not found: {input_path}")

    print("Parsing chado_FBti XML")
    print(f"  Input: {input_path}")
    if args.dry_run:
        print("  Mode: dry-run (no CSV will be written)")
    else:
        print(f"  Output: {output_path}")

    pairs: Set[Tuple[str, str]] = set()
    for idx, pair in enumerate(iter_fbtp_fbti_pairs(input_path), start=1):
        pairs.add(pair)
        if idx % 100000 == 0:
            print(f"  Parsed {idx:,} relationship rows...")

    unique_fbtp = {fbtp_id for fbtp_id, _fbti_id in pairs}
    unique_fbti = {fbti_id for _fbtp_id, fbti_id in pairs}

    print("\nSummary:")
    print(f"  Unique FBtp -> FBti pairs: {len(pairs):,}")
    print(f"  Unique FBtp IDs: {len(unique_fbtp):,}")
    print(f"  Unique FBti IDs: {len(unique_fbti):,}")

    if args.dry_run:
        return 0

    write_pairs(output_path, pairs)
    print(f"  Wrote CSV: {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
