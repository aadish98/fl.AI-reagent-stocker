#!/usr/bin/env python3
"""
Stock Phenotype Sheet / Chado curation audit (plan implementation).

Writes machine-readable JSON and a short Markdown summary under audit_outputs/.
Run from repo root:

    python scripts/audit_stock_phenotype_pipeline.py
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fl_ai_reagent_stocker.pipelines.stock_splitting import (  # noqa: E402
    _build_stock_phenotype_sheet,
    _extract_flybase_ids,
)
from fl_ai_reagent_stocker.utils import clean_id, load_flybase_tsv, unique_join  # noqa: E402

ALLELES = REPO_ROOT / "data" / "flybase" / "alleles_and_stocks"
# Chado / stocks genotype strings may contain human symbols like "FBXO11" that match a loose
# "FB + letters + digits" heuristic but are not FlyBase database IDs.
CANONICAL_FB_ID_RE = re.compile(
    r"\b(FB(?:st|al|ti|tp|ba|ab|gn|sn)\d+)\b",
    re.IGNORECASE,
)
REFS = REPO_ROOT / "data" / "flybase" / "references"
OUT_DIR = REPO_ROOT / "audit_outputs"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_build_module():
    path = REPO_ROOT / "scripts" / "build_fbst_derived_stock_components.py"
    spec = importlib.util.spec_from_file_location("fbst_build", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def freeze_artifacts() -> Dict[str, Any]:
    paths = {
        "chado_FBst_xml_gz": ALLELES / "chado_FBst.xml.gz",
        "fbst_to_derived_stock_component_csv": ALLELES / "fbst_to_derived_stock_component.csv",
        "stocks_FB": ALLELES / "stocks_FB2026_01.tsv.gz",
        "fbal_to_fbgn": ALLELES / "fbal_to_fbgn_fb_2026_01.tsv.gz",
        "genotype_phenotype_data": ALLELES / "genotype_phenotype_data_fb_2026_01.tsv",
        "entity_publication": REFS / "entity_publication_fb_2026_01.tsv.gz",
        "fbrf_pmid_pmcid_doi": REFS / "fbrf_pmid_pmcid_doi_fb_2026_01.tsv.gz",
    }
    out: Dict[str, Any] = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": {},
    }
    for key, p in paths.items():
        if not p.exists():
            out["files"][key] = {"path": str(p), "exists": False}
            continue
        st = p.stat()
        out["files"][key] = {
            "path": str(p.relative_to(REPO_ROOT)),
            "exists": True,
            "size_bytes": st.st_size,
            "mtime_epoch": int(st.st_mtime),
            "sha256": _sha256(p),
        }
    # Header lines from genotype_phenotype (FlyBase snapshot)
    gp = paths["genotype_phenotype_data"]
    if gp.exists():
        with open(gp, "r", encoding="utf-8", errors="replace") as f:
            out["genotype_phenotype_header_comments"] = [next(f).rstrip() for _ in range(5)]
    out["reporting_snapshot_note"] = "Alleles/references TSVs and phenotype report use fb_2026_01 in filenames; verify against header comments."
    # Plan: freeze-downstream-artifacts (Stage 1 workbook path when present)
    candidates = sorted(REPO_ROOT.glob("gene_lists/**/aggregated_stock_refs.xlsx"))
    out["stage1_workbook"] = {
        "canonical_relative_path": "gene_lists/Stocks/aggregated_stock_refs.xlsx",
        "matches_in_repo": [str(p.relative_to(REPO_ROOT)) for p in candidates],
        "workbook_found": bool(candidates),
    }
    return out


def _serialize_counter(c: Counter) -> Dict[str, int]:
    return dict(sorted(c.items(), key=lambda x: (-x[1], x[0])))


def format_chado_xml_extraction_audit(summary: Dict[str, Any], extracted_row_count: int) -> Dict[str, Any]:
    """Serialize extract_stock_component_rows summary (plan: audit-chado-extraction)."""
    token_counts = summary.get("token_counts") or Counter()
    unexpected = summary.get("unexpected_embedded_fbst") or []
    return {
        "stocks_seen": summary.get("stocks_seen"),
        "stocks_with_genotype": summary.get("stocks_with_genotype"),
        "stocks_with_tokens": summary.get("stocks_with_tokens"),
        "stocks_with_genotype_but_zero_tokens": (
            (summary.get("stocks_with_genotype") or 0)
            - (summary.get("stocks_with_tokens") or 0)
        ),
        "total_embedded_tokens": summary.get("total_tokens"),
        "embedded_token_counts_by_type": _serialize_counter(token_counts),
        "unexpected_embedded_FBst_examples": [
            {"parent_FBst": a, "embedded_FBst": b, "raw_token_excerpt": (c or "")[:120]}
            for a, b, c in unexpected[:25]
        ],
        "unexpected_embedded_FBst_count": len(unexpected),
        "extracted_row_count": extracted_row_count,
        "TOKEN_PATTERN_doc": "Only @(FBxxdigits):label@ substrings become derived_stock_component rows.",
    }


def audit_tsv_enrichment_and_dedup(build_mod, extracted_rows: List[Dict[str, str]]) -> Dict[str, Any]:
    """Quantify join stats from build_dataframe; verify (FBst, derived_stock_component) dedup."""
    stock_lookup = build_mod.load_stock_metadata(build_mod.DEFAULT_STOCKS_PATH)
    fbal_lookup = build_mod.load_fbal_lookup(build_mod.DEFAULT_FBAL_PATH)
    ep_path = build_mod.find_latest_tsv(build_mod.REFERENCES_DIR, "entity_publication")
    entity_lookup, _per_entity_prefix, entity_conflicts = build_mod.load_entity_name_lookup(ep_path)

    df, build_summary = build_mod.build_dataframe(
        extracted_rows=extracted_rows,
        stock_lookup=stock_lookup,
        fbal_lookup=fbal_lookup,
        entity_lookup=entity_lookup,
    )

    # Duplicate (FBst, component) rows with differing raw token text (pre-build dedup)
    raw_by_pair: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for r in extracted_rows:
        fbst = (r.get("FBst") or "").strip()
        comp = (r.get("derived_stock_component") or "").strip()
        raw_by_pair[(fbst, comp)].add((r.get("_xml_token_raw") or "").strip())

    pairs_multi_raw = sum(1 for s in raw_by_pair.values() if len(s) > 1)

    bs = dict(build_summary)
    if isinstance(bs.get("non_fbal_matched"), Counter):
        bs["non_fbal_matched"] = _serialize_counter(bs["non_fbal_matched"])
    if isinstance(bs.get("non_fbal_unmatched"), Counter):
        bs["non_fbal_unmatched"] = _serialize_counter(bs["non_fbal_unmatched"])

    return {
        "build_dataframe_summary": bs,
        "output_csv_row_count": len(df),
        "entity_publication_conflict_examples": [
            {"entity_id": a, "first_name": b, "second_name": c} for a, b, c in entity_conflicts
        ],
        "dedup_verification": {
            "unique_fbst_component_pairs_pre_dedup": len(raw_by_pair),
            "pairs_with_multiple_distinct_xml_token_raw": pairs_multi_raw,
            "note": "build_dataframe drops duplicate (FBst, derived_stock_component); multiple raw tokens for same pair would collapse to one row.",
        },
    }


def audit_chado_token_gaps(build_mod) -> Dict[str, Any]:
    """Find canonical FlyBase IDs in Chado genotype names that are not captured by TOKEN_PATTERN."""
    TOKEN_PATTERN = build_mod.TOKEN_PATTERN

    stocks_with_genotype = 0
    stocks_with_any_canonical_id_in_name = 0
    stocks_with_uncaptured_ids = 0
    uncaptured_examples: List[Dict[str, str]] = []
    genotype_no_tokens = 0

    chado_path = build_mod.DEFAULT_CHADO_PATH
    with gzip.open(chado_path, "rt", encoding="utf-8") as handle:
        for _, elem in ET.iterparse(handle, events=("end",)):
            if elem.tag != "stock":
                continue
            fbst = (elem.findtext("./dbxref_id/dbxref/accession") or "").strip()
            if not fbst:
                fbst = (elem.findtext("./uniquename") or "").strip()
            genotype_name = (
                elem.findtext("./stock_genotype/genotype_id/genotype/name") or ""
            ).strip()
            if genotype_name:
                stocks_with_genotype += 1
                token_ids = {clean_id(t[0]) for t in TOKEN_PATTERN.findall(genotype_name)}
                canonical_ids = {clean_id(m.group(1)) for m in CANONICAL_FB_ID_RE.finditer(genotype_name)}
                if canonical_ids:
                    stocks_with_any_canonical_id_in_name += 1
                if not TOKEN_PATTERN.findall(genotype_name):
                    genotype_no_tokens += 1
                uncaptured = canonical_ids - token_ids
                if uncaptured:
                    stocks_with_uncaptured_ids += 1
                    if len(uncaptured_examples) < 15:
                        uncaptured_examples.append(
                            {
                                "FBst": fbst,
                                "uncaptured_canonical_ids": sorted(uncaptured)[:8],
                                "genotype_name_excerpt": genotype_name[:240],
                            }
                        )
            elem.clear()

    return {
        "stocks_with_genotype": stocks_with_genotype,
        "stocks_with_any_canonical_fb_id_in_name": stocks_with_any_canonical_id_in_name,
        "stocks_with_canonical_ids_not_in_token_pattern": stocks_with_uncaptured_ids,
        "stocks_with_genotype_but_zero_tokens": genotype_no_tokens,
        "uncaptured_id_examples": uncaptured_examples,
        "note": "Canonical ID regex is FBst|FBal|FBti|FBtp|FBba|FBab|FBgn|FBsn + digits; excludes human symbols like FBXO11.",
    }


def audit_derived_vs_stocks_tsv(
    derived_path: Path, stocks_path: Path, max_compare: int = 200_000
) -> Dict[str, Any]:
    """Compare derived_stock_component IDs per FBst to IDs appearing in stocks TSV FB_genotype."""
    derived = pd.read_csv(derived_path, dtype=str).fillna("")
    derived["FBst"] = derived["FBst"].apply(clean_id)
    derived["derived_stock_component"] = derived["derived_stock_component"].apply(clean_id)
    by_fst: Dict[str, Set[str]] = defaultdict(set)
    for _, row in derived.iterrows():
        f, c = row["FBst"], row["derived_stock_component"]
        if f and c:
            by_fst[f].add(c)

    stocks = pd.read_csv(stocks_path, sep="\t", compression="gzip", dtype=str).fillna("")
    if "FBst" not in stocks.columns:
        return {"error": "stocks TSV missing FBst"}
    stocks["FBst"] = stocks["FBst"].apply(clean_id)
    gt_col = "FB_genotype" if "FB_genotype" in stocks.columns else None
    if not gt_col:
        return {"error": "stocks TSV missing FB_genotype"}

    only_in_derived = 0
    only_in_tsv = 0
    both_nonempty_diff = 0
    examples: List[Dict[str, Any]] = []

    checked = 0
    for _, srow in stocks.iterrows():
        fbst = srow["FBst"]
        if not fbst:
            continue
        dset = by_fst.get(fbst, set())
        tset = {clean_id(m.group(1)) for m in CANONICAL_FB_ID_RE.finditer(str(srow.get(gt_col, "")))}
        checked += 1
        if checked > max_compare:
            break
        if not dset and not tset:
            continue
        if dset and not tset:
            only_in_derived += 1
        elif tset and not dset:
            only_in_tsv += 1
        elif dset != tset:
            both_nonempty_diff += 1
            if len(examples) < 12:
                examples.append(
                    {
                        "FBst": fbst,
                        "only_in_derived": sorted(dset - tset)[:6],
                        "only_in_tsv_genotype": sorted(tset - dset)[:6],
                    }
                )

    return {
        "stocks_rows_scanned": checked,
        "fbst_only_in_derived_csv_not_in_tsv_fb_genotype_ids": only_in_derived,
        "fbst_only_in_tsv_fb_genotype_not_in_derived_csv": only_in_tsv,
        "fbst_both_nonempty_but_set_mismatch": both_nonempty_diff,
        "mismatch_examples": examples,
        "note": "TSV FB_genotype is parsed with the same canonical FB* prefixes as Chado audit; human gene text in braces is ignored.",
    }


def audit_phenotype_edge_cases(phenotype_path: Path) -> Dict[str, Any]:
    df = load_flybase_tsv(phenotype_path)
    df["genotype_FBids"] = df.get("genotype_FBids", "").fillna("").astype(str)
    df["phenotype_name"] = df.get("phenotype_name", "").fillna("").astype(str)
    df["reference"] = df.get("reference", "").fillna("").astype(str)

    loose = df["genotype_FBids"].apply(lambda s: _extract_flybase_ids(s))
    empty_ids = int((loose.apply(len) == 0).sum())

    def would_drop_sheet_logic(row) -> bool:
        pn = str(row["phenotype_name"] or "").strip()
        qual = str(row.get("qualifier_names", "") or "").strip()
        if pn and qual:
            qs = [q.strip() for q in qual.split("|") if q.strip()]
            if qs:
                pn = f"{pn} ({', '.join(qs)})"
        raw_fbrfs = [clean_id(v) for v in str(row.get("reference", "")).split("|") if clean_id(v)]
        return not pn and not raw_fbrfs

    drop_both = int(df.apply(would_drop_sheet_logic, axis=1).sum())

    # genotype_FBids containing substrings that might not match _extract (none expected if standard)
    weird = df[df["genotype_FBids"].str.contains(r"FB[a-zA-Z]{5,}[0-9]", regex=True, na=False)]
    long_class = int(len(weird))

    return {
        "total_rows": len(df),
        "rows_empty_after_id_extraction": empty_ids,
        "rows_empty_phenotype_and_reference_would_drop_in_sheet_builder": drop_both,
        "rows_with_unusual_id_pattern_FB_plus_5plus_letters": long_class,
    }


def audit_entity_lookup_conflicts(build_mod) -> Dict[str, Any]:
    """entity_publication first-wins conflicts (same as build script)."""
    ep_path = build_mod.find_latest_tsv(build_mod.REFERENCES_DIR, "entity_publication")
    _lookup, per_type, conflicts = build_mod.load_entity_name_lookup(ep_path)
    return {
        "entity_publication_path": str(ep_path.relative_to(REPO_ROOT)),
        "entity_lookup_rows": len(_lookup),
        "conflict_examples_count_returned": len(conflicts),
        "conflict_examples": [{"entity_id": a, "first": b, "second": c} for a, b, c in conflicts],
        "per_type_prefix_counts_top_12": per_type.most_common(12),
    }


def validate_phenotype_ids_vs_derived_for_stocks(
    derived_path: Path,
    phenotype_path: Path,
    sample_fbsts: List[str],
    max_examples: int = 25,
) -> Dict[str, Any]:
    """
    For phenotype rows that intersect a stock's Chado-derived components, flag rows whose
    genotype_FBids also mention IDs not on that stock (compound genotypes; plan: validate-derived-vs-flybase).
    """
    derived = pd.read_csv(derived_path, dtype=str).fillna("")
    derived["FBst"] = derived["FBst"].apply(clean_id)
    derived["derived_stock_component"] = derived["derived_stock_component"].apply(clean_id)
    by_fst: Dict[str, Set[str]] = defaultdict(set)
    for _, row in derived.iterrows():
        f, c = row["FBst"], row["derived_stock_component"]
        if f and c:
            by_fst[f].add(c)

    pheno = load_flybase_tsv(phenotype_path)
    pheno = pheno.copy()
    pheno["genotype_FBids"] = pheno.get("genotype_FBids", "").fillna("").astype(str)
    pheno["_G"] = pheno["genotype_FBids"].apply(lambda s: set(_extract_flybase_ids(s)))

    flagged_rows = 0
    examples: List[Dict[str, Any]] = []
    rows_checked = 0

    for fbst in sample_fbsts:
        cset = by_fst.get(fbst, set())
        if not cset:
            continue
        hit_mask = pheno["_G"].apply(lambda g: bool(g & cset))
        sub = pheno[hit_mask]
        rows_checked += len(sub)
        for _, row in sub.iterrows():
            g = row["_G"]
            extra = g - cset
            if not extra:
                continue
            flagged_rows += 1
            if len(examples) < max_examples:
                inter = g & cset
                examples.append(
                    {
                        "FBst": fbst,
                        "extra_ids_not_on_this_stock_chado": sorted(extra)[:10],
                        "intersection_with_stock": sorted(inter)[:10],
                        "genotype_FBids_excerpt": str(row["genotype_FBids"])[:160],
                    }
                )

    return {
        "sample_FBsts": sample_fbsts,
        "phenotype_rows_checked_across_samples": rows_checked,
        "linked_rows_with_some_FBid_not_on_stock_derived_set": flagged_rows,
        "examples": examples,
        "note": (
            "Expected for multi-allele genotype_phenotype rows: sheet still links via intersection. "
            "relevant_component_ids can add IDs not present in Chado-derived CSV (Stage 1 expansion)."
        ),
    }


def n_fb_n_rel_sample(
    derived_path: Path,
    phenotype_path: Path,
    sample_fbsts: List[str],
) -> List[Dict[str, Any]]:
    derived = pd.read_csv(derived_path, dtype=str).fillna("")
    derived["FBst"] = derived["FBst"].apply(clean_id)
    derived["derived_stock_component"] = derived["derived_stock_component"].apply(clean_id)

    pheno = load_flybase_tsv(phenotype_path)
    pheno["genotype_FBids"] = pheno.get("genotype_FBids", "").fillna("").astype(str)

    def count_hits(component_ids: Set[str]) -> int:
        if not component_ids:
            return 0

        def hit(s):
            return bool(set(_extract_flybase_ids(s)) & component_ids)

        return int(pheno["genotype_FBids"].apply(hit).sum())

    rows_out: List[Dict[str, Any]] = []
    for fbst in sample_fbsts:
        sub = derived[derived["FBst"] == fbst]
        full_ids = set(sub["derived_stock_component"].tolist())
        fbal_only = {
            clean_id(x)
            for x in sub["derived_stock_component"].tolist()
            if str(x).startswith("FBal")
        }
        n_full = count_hits(full_ids)
        n_fbal = count_hits(fbal_only) if fbal_only else 0
        rows_out.append(
            {
                "FBst": fbst,
                "n_derived_components": len(full_ids),
                "N_phenotype_rows_any_derived_component": n_full,
                "N_phenotype_rows_if_only_FBal_in_all_relevant": n_fbal,
                "gap_full_minus_fbal_restricted": n_full - n_fbal,
            }
        )
    return rows_out


def spot_check_stocks(
    derived_path: Path,
    phenotype_path: Path,
    fbsts: List[str],
) -> List[Dict[str, Any]]:
    derived = pd.read_csv(derived_path, dtype=str).fillna("")
    derived["FBst"] = derived["FBst"].apply(clean_id)
    pheno = load_flybase_tsv(phenotype_path)
    pheno["genotype_FBids"] = pheno.get("genotype_FBids", "").fillna("").astype(str)

    out = []
    for fbst in fbsts:
        sub = derived[derived["FBst"] == fbst]
        ids = set(sub["derived_stock_component"].apply(clean_id))

        def hit(s):
            return bool(set(_extract_flybase_ids(s)) & ids)

        subp = pheno[pheno["genotype_FBids"].apply(hit)]
        out.append(
            {
                "FBst": fbst,
                "derived_components": sorted(ids),
                "embedded_type_counts": sub["embedded_type"].value_counts().to_dict(),
                "n_matching_phenotype_rows": len(subp),
                "sample_phenotypes": subp["phenotype_name"].drop_duplicates().head(8).tolist(),
            }
        )
    return out


def relevant_component_ablation(flybase_data_path: Path) -> List[Dict[str, Any]]:
    """Reconstruct join: phenotype sheet row counts when relevant_component_ids is narrowed."""
    fbst = "FBst0000006"
    full_ids = "FBal0018186; FBal0028921; FBba0000025; FBti0002044"
    narrow_ids = "FBal0018186"
    rows: List[Dict[str, Any]] = []
    for label, rc in [("full_relevant_components", full_ids), ("fbal_only_FBal0018186", narrow_ids)]:
        included = pd.DataFrame(
            [
                {
                    "FBst": fbst,
                    "stock_number": "25040",
                    "collection": "Bloomington",
                    "relevant_gene_symbols": "w; Adh",
                    "relevant_component_ids": rc,
                    "relevant_fbal_ids": "",
                    "relevant_fbal_symbols": "",
                    "custom_stock": False,
                }
            ]
        )
        df = _build_stock_phenotype_sheet(
            included,
            flybase_data_path,
            references_df=None,
            unfiltered_references_df=None,
            pubmed_cache_path=None,
            pubmed_client=None,
            verbose=False,
        )
        rows.append({"scenario": label, "sheet_row_count": len(df)})
    if len(rows) == 2:
        rows.append(
            {
                "scenario": "delta_full_minus_narrow",
                "sheet_row_count": rows[0]["sheet_row_count"] - rows[1]["sheet_row_count"],
            }
        )
    return rows


def synthetic_sheet_smoke(flybase_data_path: Path) -> Dict[str, Any]:
    """Call _build_stock_phenotype_sheet with minimal included_df; ensure non-empty or explain."""
    fbst = "FBst0000006"
    included = pd.DataFrame(
        [
            {
                "FBst": fbst,
                "stock_number": "25040",
                "collection": "Bloomington",
                "relevant_gene_symbols": "w; Adh",
                "relevant_component_ids": unique_join(
                    [
                        "FBal0018186",
                        "FBal0028921",
                        "FBba0000025",
                        "FBti0002044",
                    ]
                ),
                "relevant_fbal_ids": "FBal0018186; FBal0028921",
                "relevant_fbal_symbols": "w[1118]; Adv[1]",
                "custom_stock": False,
            }
        ]
    )
    df = _build_stock_phenotype_sheet(
        included,
        flybase_data_path,
        references_df=None,
        unfiltered_references_df=None,
        pubmed_cache_path=None,
        pubmed_client=None,
        verbose=False,
    )
    return {
        "synthetic_FBst": fbst,
        "phenotype_sheet_rows": len(df),
        "columns": list(df.columns) if len(df) else [],
    }


def write_markdown(path: Path, payload: Dict[str, Any]) -> None:
    lines = [
        "# Stock Phenotype / Chado curation audit",
        "",
        f"Generated (UTC): `{payload['generated_at_utc']}`",
        "",
        "## Frozen artifacts",
        "",
        "SHA256 and paths recorded in `audit_summary.json` (includes `stage1_workbook` search).",
        "",
        "## Chado XML extraction (build script parity)",
        "",
        "```json",
        json.dumps(payload.get("audit_chado_extraction", {}), indent=2)[:10000],
        "```",
        "",
        "## TSV enrichment + dedup (build_dataframe)",
        "",
        "```json",
        json.dumps(payload.get("audit_tsv_enrichment", {}), indent=2)[:12000],
        "```",
        "",
        "## Canonical ID vs TOKEN_PATTERN (second XML pass)",
        "",
        "```json",
        json.dumps(payload.get("chado_token_gaps", {}), indent=2)[:8000],
        "```",
        "",
        "## Derived CSV vs stocks TSV FB_genotype",
        "",
        "```json",
        json.dumps(payload.get("derived_vs_stocks", {}), indent=2)[:8000],
        "```",
        "",
        "## genotype_phenotype_data edge cases",
        "",
        "```json",
        json.dumps(payload.get("phenotype_edges", {}), indent=2),
        "```",
        "",
        "## entity_publication lookup",
        "",
        "```json",
        json.dumps(payload.get("entity_lookup_audit", {}), indent=2)[:6000],
        "```",
        "",
        "## Phenotype genotype_FBids vs derived components (reachability)",
        "",
        "```json",
        json.dumps(payload.get("phenotype_derived_reachability", {}), indent=2)[:12000],
        "```",
        "",
        "## N_fb vs restricted (FBal-only) phenotype row counts (sample FBsts)",
        "",
        "```json",
        json.dumps(payload.get("n_fb_n_rel_sample", []), indent=2),
        "```",
        "",
        "## Spot-check stocks",
        "",
        "```json",
        json.dumps(payload.get("spot_check", []), indent=2)[:12000],
        "```",
        "",
        "## Synthetic _build_stock_phenotype_sheet smoke",
        "",
        "```json",
        json.dumps(payload.get("synthetic_sheet", {}), indent=2),
        "```",
        "",
        "## Relevant-component ablation (reconstruct join)",
        "",
        "```json",
        json.dumps(payload.get("relevant_component_ablation", []), indent=2),
        "```",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    build_mod = _load_build_module()

    payload: Dict[str, Any] = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    payload["frozen_artifacts"] = freeze_artifacts()
    _blog = OUT_DIR / "build_fbst_derived_console.log"
    payload["build_fbst_derived_console_log"] = (
        str(_blog.relative_to(REPO_ROOT)) if _blog.exists() else None
    )

    print("    Chado XML: extract_stock_component_rows + build_dataframe audit...")
    extracted_rows, ex_summary = build_mod.extract_stock_component_rows(build_mod.DEFAULT_CHADO_PATH)
    payload["audit_chado_extraction"] = format_chado_xml_extraction_audit(
        ex_summary, len(extracted_rows)
    )
    payload["audit_tsv_enrichment"] = audit_tsv_enrichment_and_dedup(build_mod, extracted_rows)
    extracted_rows.clear()

    print("    Chado XML: canonical ID vs TOKEN_PATTERN pass...")
    payload["chado_token_gaps"] = audit_chado_token_gaps(build_mod)
    payload["derived_vs_stocks"] = audit_derived_vs_stocks_tsv(
        ALLELES / "fbst_to_derived_stock_component.csv",
        ALLELES / "stocks_FB2026_01.tsv.gz",
    )
    payload["phenotype_edges"] = audit_phenotype_edge_cases(
        ALLELES / "genotype_phenotype_data_fb_2026_01.tsv"
    )
    payload["entity_lookup_audit"] = audit_entity_lookup_conflicts(build_mod)

    # Pick sample FBsts: plan examples + high component count
    derived = pd.read_csv(ALLELES / "fbst_to_derived_stock_component.csv", dtype=str)
    derived["FBst"] = derived["FBst"].apply(clean_id)
    counts = derived.groupby("FBst").size().sort_values(ascending=False)
    raw_sample = [
        "FBst0000002",
        "FBst0000006",
        "FBst0051299",
        counts.index[0] if len(counts) else "",
        counts.index[min(100, len(counts) - 1)] if len(counts) > 100 else "",
    ]
    seen: Set[str] = set()
    sample = []
    for s in raw_sample:
        if s and s not in seen:
            seen.add(s)
            sample.append(s)

    payload["phenotype_derived_reachability"] = validate_phenotype_ids_vs_derived_for_stocks(
        ALLELES / "fbst_to_derived_stock_component.csv",
        ALLELES / "genotype_phenotype_data_fb_2026_01.tsv",
        sample,
    )

    payload["n_fb_n_rel_sample"] = n_fb_n_rel_sample(
        ALLELES / "fbst_to_derived_stock_component.csv",
        ALLELES / "genotype_phenotype_data_fb_2026_01.tsv",
        sample,
    )
    payload["spot_check"] = spot_check_stocks(
        ALLELES / "fbst_to_derived_stock_component.csv",
        ALLELES / "genotype_phenotype_data_fb_2026_01.tsv",
        sample[:6],
    )
    payload["synthetic_sheet"] = synthetic_sheet_smoke(REPO_ROOT / "data" / "flybase")
    payload["relevant_component_ablation"] = relevant_component_ablation(
        REPO_ROOT / "data" / "flybase"
    )

    json_path = OUT_DIR / "audit_summary.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_markdown(OUT_DIR / "STOCK_PHENOTYPE_AUDIT.md", payload)
    print(f"Wrote {json_path} and {OUT_DIR / 'STOCK_PHENOTYPE_AUDIT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
