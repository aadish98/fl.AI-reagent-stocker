#!/usr/bin/env python3
"""
Build an FBsf -> FBgn mapping CSV from a FlyBase PostgreSQL dump.

The parser streams the dump and reads only the ``feature`` and
``feature_relationship`` COPY sections. It keeps only FBsf / FBgn feature IDs in
memory, then joins feature relationships where the subject is an FBsf and the
object is an FBgn.

Usage:
    python scripts/build_fbsf_to_fbgn_mapping.py
    python scripts/build_fbsf_to_fbgn_mapping.py --sql-dump data/flybase/FB2026_01.sql.gz
    python scripts/build_fbsf_to_fbgn_mapping.py --dry-run
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Set, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
FLYBASE_DIR = REPO_ROOT / "data" / "flybase"
SEQUENCE_FEATURES_DIR = FLYBASE_DIR / "sequence_features"
DEFAULT_OUTPUT_PATH = SEQUENCE_FEATURES_DIR / "fbsf_to_fbgn.csv"


def _default_sql_dump() -> Path:
    matches = sorted(FLYBASE_DIR.glob("FB*.sql.gz"), reverse=True)
    if matches:
        return matches[0]
    return FLYBASE_DIR / "FB2026_01.sql.gz"


DEFAULT_SQL_DUMP = _default_sql_dump()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an FBsf-to-FBgn mapping CSV from a FlyBase SQL dump."
    )
    parser.add_argument(
        "--sql-dump",
        type=Path,
        default=DEFAULT_SQL_DUMP,
        help=(
            "Path to a FlyBase PostgreSQL dump (.sql or .sql.gz) containing "
            "the feature and feature_relationship COPY sections."
        ),
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


def open_text(path: Path):
    """Open a plain-text or gzip-compressed SQL dump."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _parse_copy_columns(line: str, table_name: str) -> List[str]:
    prefix = f"COPY public.{table_name} ("
    suffix = ") FROM stdin;"
    if not line.startswith(prefix) or not line.rstrip().endswith(suffix):
        return []
    columns_blob = line[len(prefix) : line.rfind(suffix)]
    return [column.strip() for column in columns_blob.split(",") if column.strip()]


def _row_dict(columns: Sequence[str], raw_line: str) -> Dict[str, str]:
    values = raw_line.rstrip("\n").split("\t")
    padded_values = values + [""] * max(0, len(columns) - len(values))
    row = dict(zip(columns, padded_values))
    return {
        key: ("" if value == r"\N" else value.strip())
        for key, value in row.items()
    }


def extract_fbsf_to_fbgn_pairs(sql_dump_path: Path) -> Tuple[Set[Tuple[str, str]], Dict[str, int]]:
    """
    Stream a FlyBase SQL dump and extract FBsf -> FBgn relationships.

    Returns:
        (pairs, summary_counts)
    """
    feature_lookup: Dict[str, str] = {}
    relationship_pairs: Set[Tuple[str, str]] = set()
    current_table = ""
    current_columns: List[str] = []

    feature_rows = 0
    feature_relationship_rows = 0
    kept_features = 0
    kept_relationships = 0
    saw_feature = False
    saw_feature_relationship = False

    with open_text(sql_dump_path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if not current_table:
                feature_columns = _parse_copy_columns(line, "feature")
                if feature_columns:
                    current_table = "feature"
                    current_columns = feature_columns
                    saw_feature = True
                    continue

                feature_relationship_columns = _parse_copy_columns(
                    line, "feature_relationship"
                )
                if feature_relationship_columns:
                    current_table = "feature_relationship"
                    current_columns = feature_relationship_columns
                    saw_feature_relationship = True
                    continue

                continue

            if line == r"\.":
                if saw_feature and saw_feature_relationship and current_table == "feature_relationship":
                    break
                current_table = ""
                current_columns = []
                continue

            row = _row_dict(current_columns, raw_line)
            if current_table == "feature":
                feature_rows += 1
                feature_id = row.get("feature_id", "")
                uniquename = row.get("uniquename", "")
                if feature_id and (
                    uniquename.startswith("FBsf") or uniquename.startswith("FBgn")
                ):
                    feature_lookup[feature_id] = uniquename
                    kept_features += 1
            elif current_table == "feature_relationship":
                feature_relationship_rows += 1
                subject_id = row.get("subject_id", "")
                object_id = row.get("object_id", "")
                subject_unique = feature_lookup.get(subject_id, "")
                object_unique = feature_lookup.get(object_id, "")
                if subject_unique.startswith("FBsf") and object_unique.startswith("FBgn"):
                    relationship_pairs.add((subject_unique, object_unique))
                    kept_relationships += 1

    if not saw_feature:
        raise ValueError(f"Did not find COPY public.feature in {sql_dump_path}")
    if not saw_feature_relationship:
        raise ValueError(
            f"Did not find COPY public.feature_relationship in {sql_dump_path}"
        )

    summary = {
        "feature_rows": feature_rows,
        "feature_relationship_rows": feature_relationship_rows,
        "kept_features": kept_features,
        "kept_relationships": kept_relationships,
    }
    return relationship_pairs, summary


def write_pairs(output_path: Path, pairs: Iterable[Tuple[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["FBsf", "FBgn"])
        for fbsf_id, fbgn_id in sorted(set(pairs)):
            writer.writerow([fbsf_id, fbgn_id])


def main() -> int:
    args = parse_args()
    sql_dump_path = Path(args.sql_dump)
    output_path = Path(args.output)

    if not sql_dump_path.exists():
        raise FileNotFoundError(f"SQL dump not found: {sql_dump_path}")

    print("Building FBsf -> FBgn mapping CSV")
    print(f"  SQL dump: {sql_dump_path}")
    if args.dry_run:
        print("  Mode: dry-run (no CSV will be written)")
    else:
        print(f"  Output: {output_path}")

    pairs, summary = extract_fbsf_to_fbgn_pairs(sql_dump_path)
    unique_fbsf = {fbsf_id for fbsf_id, _fbgn_id in pairs}
    unique_fbgn = {fbgn_id for _fbsf_id, fbgn_id in pairs}

    print("\nSummary:")
    print(f"  feature rows scanned: {summary['feature_rows']:,}")
    print(
        "  feature_relationship rows scanned: "
        f"{summary['feature_relationship_rows']:,}"
    )
    print(f"  FBsf/FBgn features kept in lookup: {summary['kept_features']:,}")
    print(
        "  qualifying feature_relationship rows: "
        f"{summary['kept_relationships']:,}"
    )
    print(f"  Unique FBsf -> FBgn pairs: {len(pairs):,}")
    print(f"  Unique FBsf IDs: {len(unique_fbsf):,}")
    print(f"  Unique FBgn IDs: {len(unique_fbgn):,}")

    if args.dry_run:
        return 0

    write_pairs(output_path, pairs)
    print(f"  Wrote CSV: {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
