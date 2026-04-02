#!/usr/bin/env python3
"""
Build a CSV of FBst -> derived stock component mappings from FlyBase Chado XML.

The Chado XML is used only for the raw FBst -> FB... association extraction.
All other metadata is filled from TSV sources:
  - stocks_FB*.tsv for stock metadata
  - fbal_to_fbgn*.tsv(.gz) for FBal symbol/gene metadata
  - entity_publication*.tsv(.gz) for non-FBal object symbols

Usage:
    python scripts/build_fbst_derived_stock_components.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
ALLELES_AND_STOCKS_DIR = REPO_ROOT / "data" / "flybase" / "alleles_and_stocks"
REFERENCES_DIR = REPO_ROOT / "data" / "flybase" / "references"

DEFAULT_CHADO_PATH = ALLELES_AND_STOCKS_DIR / "chado_FBst.xml.gz"
DEFAULT_OUTPUT_PATH = ALLELES_AND_STOCKS_DIR / "fbst_to_derived_stock_component.csv"

TOKEN_PATTERN = re.compile(r"@(FB[A-Za-z]{2}\d+):([^@]+)@")


def open_text(path: Path):
    """Open plain or gzip-compressed text files."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def find_latest_tsv(directory: Path, pattern: str) -> Path:
    """Return the latest TSV or TSV.GZ matching a prefix."""
    gz_matches = sorted(directory.glob(f"{pattern}*.tsv.gz"), reverse=True)
    if gz_matches:
        return gz_matches[0]
    tsv_matches = sorted(directory.glob(f"{pattern}*.tsv"), reverse=True)
    if tsv_matches:
        return tsv_matches[0]
    raise FileNotFoundError(f"No TSV file matching '{pattern}' found in {directory}")


DEFAULT_STOCKS_PATH = find_latest_tsv(ALLELES_AND_STOCKS_DIR, "stocks_FB")
DEFAULT_FBAL_PATH = find_latest_tsv(ALLELES_AND_STOCKS_DIR, "fbal_to_fbgn_fb_")


def iter_flybase_tsv_rows(path: Path) -> Iterator[Dict[str, str]]:
    """
    Iterate FlyBase TSV rows while handling commented headers.

    Supports files where the header line is prefixed by '#'.
    """
    with open_text(path) as handle:
        header: List[str] | None = None
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#"):
                candidate = line[1:].lstrip("\t")
                if "\t" in candidate and header is None:
                    header = [part.strip() for part in candidate.split("\t") if part.strip()]
                continue

            parts = line.split("\t")
            if header is None:
                header = [part.strip() for part in parts]
                continue

            if len(parts) != len(header):
                continue

            yield dict(zip(header, parts))


def load_stock_metadata(path: Path) -> Dict[str, Dict[str, str]]:
    """Load stock-level metadata keyed by FBst."""
    stock_lookup: Dict[str, Dict[str, str]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            fbst = (row.get("FBst") or "").strip()
            if not fbst:
                continue
            stock_lookup[fbst] = {
                "stock_number": (row.get("stock_number") or "").strip(),
                "collection": (
                    (row.get("collection") or row.get("collection_short_name") or "").strip()
                ),
                "FB_genotype": (row.get("FB_genotype") or "").strip(),
            }
    return stock_lookup


def load_fbal_lookup(path: Path) -> Dict[str, Dict[str, str]]:
    """Load FBal -> symbol/gene metadata."""
    fbal_lookup: Dict[str, Dict[str, str]] = {}
    for row in iter_flybase_tsv_rows(path):
        fbal = (row.get("AlleleID") or "").strip()
        if not fbal:
            continue
        fbal_lookup[fbal] = {
            "AlleleSymbol": (row.get("AlleleSymbol") or "").strip(),
            "FBgnID": (row.get("GeneID") or "").strip(),
            "GeneSymbol": (row.get("GeneSymbol") or "").strip(),
        }
    return fbal_lookup


def load_entity_name_lookup(path: Path) -> Tuple[Dict[str, str], Counter, List[Tuple[str, str, str]]]:
    """
    Load a generic entity_id -> entity_name lookup.

    Returns the lookup, per-type counts, and a few examples of conflicting names.
    """
    entity_lookup: Dict[str, str] = {}
    per_type_counts: Counter = Counter()
    conflicting_examples: List[Tuple[str, str, str]] = []

    for row in iter_flybase_tsv_rows(path):
        entity_id = (row.get("entity_id") or "").strip()
        entity_name = (row.get("entity_name") or "").strip()
        if not entity_id or not entity_name:
            continue

        per_type_counts[entity_id[:4]] += 1
        previous_name = entity_lookup.get(entity_id)
        if previous_name is None:
            entity_lookup[entity_id] = entity_name
        elif previous_name != entity_name and len(conflicting_examples) < 10:
            conflicting_examples.append((entity_id, previous_name, entity_name))

    return entity_lookup, per_type_counts, conflicting_examples


def extract_stock_component_rows(
    chado_path: Path,
) -> Tuple[List[Dict[str, str]], Dict[str, object]]:
    """Stream the Chado XML and extract raw FBst -> FB... mappings."""
    extracted_rows: List[Dict[str, str]] = []
    token_counts: Counter = Counter()
    stocks_seen = 0
    stocks_with_genotype = 0
    stocks_with_tokens = 0
    unexpected_embedded_fbst: List[Tuple[str, str, str]] = []

    with gzip.open(chado_path, "rt", encoding="utf-8") as handle:
        for _, elem in ET.iterparse(handle, events=("end",)):
            if elem.tag != "stock":
                continue

            stocks_seen += 1
            fbst = (elem.findtext("./dbxref_id/dbxref/accession") or "").strip()
            if not fbst:
                fbst = (elem.findtext("./uniquename") or "").strip()

            genotype_name = (
                elem.findtext("./stock_genotype/genotype_id/genotype/name") or ""
            ).strip()

            if genotype_name:
                stocks_with_genotype += 1
                matches = TOKEN_PATTERN.findall(genotype_name)
                if matches:
                    stocks_with_tokens += 1

                for derived_stock_component, raw_token_text in matches:
                    embedded_type = derived_stock_component[:4]
                    token_counts[embedded_type] += 1
                    if embedded_type == "FBst" and len(unexpected_embedded_fbst) < 20:
                        unexpected_embedded_fbst.append(
                            (fbst, derived_stock_component, raw_token_text)
                        )
                    extracted_rows.append(
                        {
                            "FBst": fbst,
                            "derived_stock_component": derived_stock_component,
                            "embedded_type": embedded_type,
                            "_xml_token_raw": raw_token_text,
                        }
                    )

            elem.clear()

    summary = {
        "stocks_seen": stocks_seen,
        "stocks_with_genotype": stocks_with_genotype,
        "stocks_with_tokens": stocks_with_tokens,
        "token_counts": token_counts,
        "total_tokens": sum(token_counts.values()),
        "unexpected_embedded_fbst": unexpected_embedded_fbst,
    }
    return extracted_rows, summary


def format_count(value: int) -> str:
    return f"{value:,}"


def print_counts_by_type(title: str, counts: Iterable[Tuple[str, int]]) -> None:
    print(title)
    for token_type, count in counts:
        print(f"  {token_type}: {format_count(count)}")


def build_dataframe(
    extracted_rows: List[Dict[str, str]],
    stock_lookup: Dict[str, Dict[str, str]],
    fbal_lookup: Dict[str, Dict[str, str]],
    entity_lookup: Dict[str, str],
) -> Tuple[pd.DataFrame, Dict[str, int | Counter]]:
    """Combine extracted XML rows with TSV-based metadata."""
    df = pd.DataFrame(extracted_rows)
    if df.empty:
        empty_df = pd.DataFrame(
            columns=[
                "FBst",
                "stock_number",
                "collection",
                "FB_genotype",
                "derived_stock_component",
                "embedded_type",
                "object_symbol",
                "AlleleSymbol",
                "FBgnID",
                "GeneSymbol",
            ]
        )
        return empty_df, {
            "rows_pre_dedup": 0,
            "rows_post_dedup": 0,
            "duplicates_removed": 0,
            "stock_metadata_matched": 0,
            "stock_metadata_unmatched": 0,
            "fbal_matched": 0,
            "fbal_unmatched": 0,
            "non_fbal_matched": Counter(),
            "non_fbal_unmatched": Counter(),
        }

    rows_pre_dedup = len(df)
    df = df.drop_duplicates(subset=["FBst", "derived_stock_component"]).copy()
    rows_post_dedup = len(df)

    df["stock_number"] = df["FBst"].map(
        lambda value: stock_lookup.get(value, {}).get("stock_number", "")
    )
    df["collection"] = df["FBst"].map(
        lambda value: stock_lookup.get(value, {}).get("collection", "")
    )
    df["FB_genotype"] = df["FBst"].map(
        lambda value: stock_lookup.get(value, {}).get("FB_genotype", "")
    )

    stock_metadata_matched = int(df["FB_genotype"].fillna("").astype(bool).sum())
    stock_metadata_unmatched = int(len(df) - stock_metadata_matched)

    df["AlleleSymbol"] = ""
    df["FBgnID"] = ""
    df["GeneSymbol"] = ""
    df["object_symbol"] = ""

    fbal_mask = df["embedded_type"] == "FBal"
    non_fbal_mask = ~fbal_mask

    if fbal_mask.any():
        df.loc[fbal_mask, "AlleleSymbol"] = df.loc[fbal_mask, "derived_stock_component"].map(
            lambda value: fbal_lookup.get(value, {}).get("AlleleSymbol", "")
        )
        df.loc[fbal_mask, "FBgnID"] = df.loc[fbal_mask, "derived_stock_component"].map(
            lambda value: fbal_lookup.get(value, {}).get("FBgnID", "")
        )
        df.loc[fbal_mask, "GeneSymbol"] = df.loc[fbal_mask, "derived_stock_component"].map(
            lambda value: fbal_lookup.get(value, {}).get("GeneSymbol", "")
        )
        df.loc[fbal_mask, "object_symbol"] = df.loc[fbal_mask, "AlleleSymbol"]

    if non_fbal_mask.any():
        df.loc[non_fbal_mask, "object_symbol"] = df.loc[
            non_fbal_mask, "derived_stock_component"
        ].map(lambda value: entity_lookup.get(value, ""))

    fbal_matched = int(df.loc[fbal_mask, "AlleleSymbol"].fillna("").astype(bool).sum())
    fbal_unmatched = int(fbal_mask.sum() - fbal_matched)

    non_fbal_matched: Counter = Counter()
    non_fbal_unmatched: Counter = Counter()
    if non_fbal_mask.any():
        for token_type, group in df.loc[non_fbal_mask].groupby("embedded_type"):
            matched = int(group["object_symbol"].fillna("").astype(bool).sum())
            unmatched = int(len(group) - matched)
            non_fbal_matched[token_type] = matched
            non_fbal_unmatched[token_type] = unmatched

    df = df[
        [
            "FBst",
            "stock_number",
            "collection",
            "FB_genotype",
            "derived_stock_component",
            "embedded_type",
            "object_symbol",
            "AlleleSymbol",
            "FBgnID",
            "GeneSymbol",
        ]
    ].sort_values(
        by=["FBst", "embedded_type", "derived_stock_component"], kind="stable"
    )

    summary = {
        "rows_pre_dedup": rows_pre_dedup,
        "rows_post_dedup": rows_post_dedup,
        "duplicates_removed": rows_pre_dedup - rows_post_dedup,
        "stock_metadata_matched": stock_metadata_matched,
        "stock_metadata_unmatched": stock_metadata_unmatched,
        "fbal_matched": fbal_matched,
        "fbal_unmatched": fbal_unmatched,
        "non_fbal_matched": non_fbal_matched,
        "non_fbal_unmatched": non_fbal_unmatched,
    }
    return df, summary


def verify_examples(df: pd.DataFrame) -> None:
    """Print a few concrete stock examples for sanity-checking."""
    print("\nSanity-check examples:")
    for stock_id in ("FBst0000002", "FBst0000006"):
        subset = df[df["FBst"] == stock_id]
        print(f"  {stock_id}: {format_count(len(subset))} derived stock component row(s)")
        if subset.empty:
            continue
        preview = subset[
            ["derived_stock_component", "embedded_type", "object_symbol"]
        ].head(10)
        for _, row in preview.iterrows():
            print(
                "    "
                f"{row['derived_stock_component']} "
                f"({row['embedded_type']}) -> {row['object_symbol'] or '[no symbol match]'}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a CSV of FBst -> derived stock component mappings."
    )
    parser.add_argument(
        "--chado-path",
        type=Path,
        default=DEFAULT_CHADO_PATH,
        help="Path to chado_FBst.xml.gz",
    )
    parser.add_argument(
        "--stocks-path",
        type=Path,
        default=DEFAULT_STOCKS_PATH,
        help="Path to stocks_FB*.tsv",
    )
    parser.add_argument(
        "--fbal-path",
        type=Path,
        default=DEFAULT_FBAL_PATH,
        help="Path to fbal_to_fbgn*.tsv(.gz)",
    )
    parser.add_argument(
        "--entity-lookup-path",
        type=Path,
        default=find_latest_tsv(REFERENCES_DIR, "entity_publication"),
        help="Path to entity_publication*.tsv(.gz) used for non-FBal symbols",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="Output CSV path",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    print("Building FBst -> derived_stock_component CSV")
    print(f"  Chado XML: {args.chado_path}")
    print(f"  Stocks TSV: {args.stocks_path}")
    print(f"  FBal TSV: {args.fbal_path}")
    print(f"  Entity lookup TSV: {args.entity_lookup_path}")
    print(f"  Output CSV: {args.output}")

    stock_lookup = load_stock_metadata(args.stocks_path)
    fbal_lookup = load_fbal_lookup(args.fbal_path)
    entity_lookup, entity_lookup_counts, entity_conflicts = load_entity_name_lookup(
        args.entity_lookup_path
    )

    print("\nLoaded TSV lookups:")
    print(f"  Stock metadata rows: {format_count(len(stock_lookup))}")
    print(f"  FBal lookup rows: {format_count(len(fbal_lookup))}")
    print(f"  Generic entity lookup rows: {format_count(len(entity_lookup))}")
    if entity_conflicts:
        print(
            "  Warning: found entity IDs with conflicting names in entity lookup "
            f"(showing {len(entity_conflicts)} example(s))"
        )
        for entity_id, first_name, second_name in entity_conflicts:
            print(f"    {entity_id}: '{first_name}' vs '{second_name}'")

    extracted_rows, extraction_summary = extract_stock_component_rows(args.chado_path)

    print("\nXML extraction summary:")
    print(
        f"  Identified {format_count(extraction_summary['stocks_seen'])} "
        "FBst record(s) in the XML"
    )
    print(
        f"  Identified {format_count(extraction_summary['stocks_with_genotype'])} "
        "FBst record(s) with a genotype block"
    )
    print(
        f"  Identified {format_count(extraction_summary['stocks_with_tokens'])} "
        "FBst record(s) with at least one embedded FB... token"
    )
    print(
        f"  Identified {format_count(extraction_summary['total_tokens'])} "
        "embedded FB... token(s) across all stocks"
    )
    print_counts_by_type(
        "  Embedded token counts by type:",
        extraction_summary["token_counts"].most_common(),
    )

    unexpected_embedded_fbst = extraction_summary["unexpected_embedded_fbst"]
    if unexpected_embedded_fbst:
        print(
            "  Warning: found unexpected embedded FBst token(s): "
            f"{format_count(len(unexpected_embedded_fbst))} example(s) collected"
        )
        for fbst, embedded_id, raw_text in unexpected_embedded_fbst:
            print(f"    Parent {fbst}: {embedded_id} -> {raw_text}")
    else:
        print("  No unexpected embedded FBst tokens were found")

    df, build_summary = build_dataframe(
        extracted_rows=extracted_rows,
        stock_lookup=stock_lookup,
        fbal_lookup=fbal_lookup,
        entity_lookup=entity_lookup,
    )

    print("\nMetadata join summary:")
    print(
        f"  Stock metadata matched: {format_count(build_summary['stock_metadata_matched'])}"
    )
    print(
        f"  Stock metadata unmatched: {format_count(build_summary['stock_metadata_unmatched'])}"
    )
    print(f"  FBal rows matched: {format_count(build_summary['fbal_matched'])}")
    print(f"  FBal rows unmatched: {format_count(build_summary['fbal_unmatched'])}")

    non_fbal_matched: Counter = build_summary["non_fbal_matched"]  # type: ignore[assignment]
    non_fbal_unmatched: Counter = build_summary["non_fbal_unmatched"]  # type: ignore[assignment]
    if non_fbal_matched or non_fbal_unmatched:
        print("  Non-FBal symbol lookup by type:")
        all_types = sorted(set(non_fbal_matched) | set(non_fbal_unmatched))
        for token_type in all_types:
            matched = non_fbal_matched.get(token_type, 0)
            unmatched = non_fbal_unmatched.get(token_type, 0)
            print(
                f"    {token_type}: matched {format_count(matched)}, "
                f"unmatched {format_count(unmatched)}"
            )
    else:
        print("  No non-FBal rows were extracted")

    print("\nOutput summary:")
    print(f"  Rows before deduplication: {format_count(build_summary['rows_pre_dedup'])}")
    print(f"  Rows after deduplication: {format_count(build_summary['rows_post_dedup'])}")
    print(f"  Duplicate rows removed: {format_count(build_summary['duplicates_removed'])}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"  Wrote CSV: {args.output}")

    verify_examples(df)
    return 0


if __name__ == "__main__":
    sys.exit(main())
