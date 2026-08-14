"""FlyBase identity lookup that enumerates all synonym candidates."""

from __future__ import annotations

import importlib.util
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from fl_ai_reagent_stocker.config import REPO_ROOT


def _load_flybase_tsv(filepath, **kwargs):
    script = REPO_ROOT / "scripts" / "fetch_fbgn_ids.py"
    spec = importlib.util.spec_from_file_location("fetch_fbgn_ids", script)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.load_flybase_tsv(filepath, **kwargs)


@dataclass(frozen=True)
class FlyBaseMaps:
    current_symbol_to_fbgn: dict[str, str]
    fbgn_to_current_symbol: dict[str, str]
    lookup_to_fbgns: dict[str, tuple[str, ...]]
    source_path: Path


def load_flybase_maps(synonym_path: Path) -> FlyBaseMaps:
    mappings = _load_flybase_tsv(
        synonym_path,
        dtype={
            "current_symbol": str,
            "current_fullname": str,
            "fullname_synonym(s)": str,
            "symbol_synonym(s)": str,
            "primary_FBid": str,
            "organism_abbreviation": str,
        },
        keep_default_na=False,
    )
    mappings = mappings[mappings["organism_abbreviation"] == "Dmel"]
    mappings = mappings[mappings["primary_FBid"].str.startswith("FBgn", na=False)]

    current_symbol_to_fbgn = dict(
        zip(
            mappings["current_symbol"].astype(str).str.strip(),
            mappings["primary_FBid"].astype(str).str.strip(),
        )
    )
    fbgn_to_current_symbol = dict(
        zip(
            mappings["primary_FBid"].astype(str).str.strip(),
            mappings["current_symbol"].astype(str).str.strip(),
        )
    )
    lookup: dict[str, set[str]] = defaultdict(set)
    for _, row in mappings.iterrows():
        fbgn = str(row["primary_FBid"]).strip()
        values = [str(row["current_symbol"]).strip(), str(row["current_fullname"]).strip()]
        for column in ("fullname_synonym(s)", "symbol_synonym(s)"):
            raw = str(row.get(column, "") or "")
            if raw and raw != "-":
                values.extend(part.strip() for part in raw.split("|"))
        for value in values:
            if value and value != "-":
                lookup[value].add(fbgn)
    lookup_to_fbgns = {key: tuple(sorted(ids)) for key, ids in lookup.items()}
    return FlyBaseMaps(
        current_symbol_to_fbgn=current_symbol_to_fbgn,
        fbgn_to_current_symbol=fbgn_to_current_symbol,
        lookup_to_fbgns=lookup_to_fbgns,
        source_path=synonym_path,
    )


def maps_from_dataframe(mappings: pd.DataFrame, source_path: Path | None = None) -> FlyBaseMaps:
    """Build maps from an already-filtered Dmel FBgn table. Used by tests."""
    current_symbol_to_fbgn = dict(
        zip(
            mappings["current_symbol"].astype(str).str.strip(),
            mappings["primary_FBid"].astype(str).str.strip(),
        )
    )
    fbgn_to_current_symbol = dict(
        zip(
            mappings["primary_FBid"].astype(str).str.strip(),
            mappings["current_symbol"].astype(str).str.strip(),
        )
    )
    lookup: dict[str, set[str]] = defaultdict(set)
    for _, row in mappings.iterrows():
        fbgn = str(row["primary_FBid"]).strip()
        values = [str(row["current_symbol"]).strip()]
        raw = str(row.get("symbol_synonym(s)", "") or "")
        if raw and raw != "-":
            values.extend(part.strip() for part in raw.split("|"))
        for value in values:
            if value:
                lookup[value].add(fbgn)
    return FlyBaseMaps(
        current_symbol_to_fbgn=current_symbol_to_fbgn,
        fbgn_to_current_symbol=fbgn_to_current_symbol,
        lookup_to_fbgns={key: tuple(sorted(ids)) for key, ids in lookup.items()},
        source_path=source_path or Path("fixture"),
    )


def resolve_symbol(symbol: str, maps: FlyBaseMaps) -> dict[str, object]:
    exact_fbgn = maps.current_symbol_to_fbgn.get(symbol)
    candidates = list(maps.lookup_to_fbgns.get(symbol, ()))
    if exact_fbgn and exact_fbgn not in candidates:
        candidates.insert(0, exact_fbgn)
    candidate_symbols = [
        maps.fbgn_to_current_symbol.get(fbgn, "") for fbgn in candidates
    ]
    match_type = "unmapped"
    if exact_fbgn:
        match_type = "exact_current_symbol"
    elif len(candidates) == 1:
        match_type = "unique_synonym"
    elif len(candidates) > 1:
        match_type = "ambiguous_synonym"
    return {
        "symbol": symbol,
        "exact_current_symbol_match": bool(exact_fbgn),
        "exact_current_fbgn": exact_fbgn or "",
        "candidate_fbgn_ids": candidates,
        "candidate_primary_symbols": candidate_symbols,
        "candidate_count": len(candidates),
        "match_type": match_type,
    }


def case_insensitive_current_candidates(symbol: str, maps: FlyBaseMaps, candidates: list[str]) -> list[str]:
    needle = symbol.lower()
    return [
        cand
        for cand in candidates
        if maps.fbgn_to_current_symbol.get(cand, "").lower() == needle
    ]


def trh_family_flag(symbol: str) -> str:
    if symbol == "Trhn":
        return "Trhn"
    if symbol == "trh":
        return "trh"
    if symbol == "Trh":
        return "Trh"
    return ""
