"""Build the four publication-defined category CSVs and CSW observations."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .constants import (
    BREAKDOWN_DIR,
    CSW_4PLUS_DEFINITION,
    CSW_4PLUS_SOURCES,
    EXPERIMENT_COLUMNS,
    EXPECTED_MASTER_ROWS,
    EXPECTED_UNIQUE_FBGN,
    HLH_DEFINITION,
    HLH_PUBLICATION_GENES,
    HOMEOSTATIC_DEFINITION,
    MASTER_TABLES,
    MECHANISTIC_DEFINITION,
    MECHANISTIC_IDENTITIES,
    MIN_FREQUENCY,
    PUBLICATION_RESOLUTION_BASIS,
    TRACHEALESS_FBGN,
    TRACHEALESS_SYMBOL,
    TRHN_FBGN,
    TRHN_SYMBOL,
    TX_INPUT_DIR,
)
from .flybase import FlyBaseMaps, case_insensitive_current_candidates, resolve_symbol, trh_family_flag
from .io_utils import join_sorted, read_csv, write_csv


def _blank(value: object) -> str:
    if value is None:
        return ""
    text = str(value)
    if text in {"nan", "None", "<NA>"}:
        return ""
    return text


def load_master_tables(input_dir: Path | None = None) -> pd.DataFrame:
    root = Path(input_dir) if input_dir is not None else TX_INPUT_DIR
    frames = []
    for source_table, meta in MASTER_TABLES.items():
        path = root / meta["filename"]
        df = read_csv(path)
        missing = {"ext_gene", "gene_id", "flybase_gene_id", "frequency"} - set(df.columns)
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        if len(df) != meta["expected_rows"]:
            raise ValueError(
                f"{path} has {len(df)} rows; expected {meta['expected_rows']}"
            )
        freq = pd.to_numeric(df["frequency"], errors="coerce")
        if freq.min() < MIN_FREQUENCY:
            raise ValueError(f"{path} contains frequency < {MIN_FREQUENCY}")
        if not (df["gene_id"] == df["flybase_gene_id"]).all():
            raise ValueError(f"{path} gene_id != flybase_gene_id")
        df = df.copy()
        df["source_table"] = source_table
        df["threshold"] = meta["threshold"]
        df["direction"] = meta["direction"]
        frames.append(df)
    observations = pd.concat(frames, ignore_index=True)
    if len(observations) != EXPECTED_MASTER_ROWS:
        raise ValueError(
            f"CSW observations have {len(observations)} rows; expected {EXPECTED_MASTER_ROWS}"
        )
    unique_fbgn = observations["flybase_gene_id"].nunique()
    if unique_fbgn != EXPECTED_UNIQUE_FBGN:
        raise ValueError(
            f"CSW observations have {unique_fbgn} unique FBgn IDs; expected {EXPECTED_UNIQUE_FBGN}"
        )
    key = observations[["flybase_gene_id", "threshold", "direction"]]
    if key.duplicated().any() or key.isna().any().any() or (key == "").any().any():
        raise ValueError("Duplicate or missing (flybase_gene_id, threshold, direction)")
    return observations


def write_csw_observations(observations: pd.DataFrame, dest: Path) -> Path:
    keep = [
        "flybase_gene_id",
        "ext_gene",
        "gene_id",
        "source_table",
        "threshold",
        "direction",
        "frequency",
        "is_cycling",
        "ZT_cyclers",
        "sleep_corr_exps",
        "wake_corr_exps",
        *EXPERIMENT_COLUMNS,
    ]
    present = [col for col in keep if col in observations.columns]
    out = observations[present].copy()
    write_csv(out, dest)
    return dest


def _mechanistic_evidence_type(evidence: str) -> str:
    text = evidence.lower()
    if "hist+" in text or "history correlate" in text or "pos sleep history" in text:
        return "correlated"
    if "consistent" in text or "freq=" in text or "present in" in text:
        return "consistent"
    return "curated"


def build_mechanistic(breakdown: Path | None = None) -> pd.DataFrame:
    root = Path(breakdown) if breakdown is not None else BREAKDOWN_DIR
    core = read_csv(root / "Tx-Omics_mechanistic_genes_n=6genes.csv")
    screen = read_csv(root / "Tx-Omics_Mechanistic_Screen_Candidates_nsyb_elav.csv")
    screen = screen.rename(columns={"Gene": "ext_gene"})
    subcategory_column = (
        "Mechanistic Subcategory"
        if "Mechanistic Subcategory" in screen.columns
        else "Mechanistic Category"
    )
    screen["Priority Rank"] = pd.to_numeric(screen["Priority Rank"], errors="coerce")
    screen = screen[screen["Priority Rank"] <= 6]
    merged = core.merge(screen, on=["ext_gene", "flybase_gene_id"], how="left")
    expected = list(MECHANISTIC_IDENTITIES)
    got = list(zip(merged["ext_gene"], merged["flybase_gene_id"]))
    if got != expected:
        raise ValueError(f"Mechanistic identities {got} != {expected}")
    merged["mechanistic_evidence_type"] = merged[
        "Tx-Omics Evidence (Rosensweig-Shah 2026)"
    ].map(_mechanistic_evidence_type)
    merged["gene_set"] = "Mechanistic"
    merged["gene_set_definition"] = MECHANISTIC_DEFINITION
    merged["mechanistic_subcategory"] = (
        merged[subcategory_column]
        .astype(str)
        .str.replace(r"^Homeostatic\s*/\s*", "", regex=True)
    )
    merged["identity_status"] = "source"
    columns = [
        "ext_gene",
        "flybase_gene_id",
        "gene_set",
        "gene_set_definition",
        "mechanistic_evidence_type",
        "mechanistic_subcategory",
        "Literature Status",
        "Tx-Omics Evidence (Rosensweig-Shah 2026)",
        "Proposed Mechanism",
        "fl.ai Category",
        "fl.ai Confidence",
        "fl.ai Rationale",
        "Key References (PMCIDs)",
        "identity_status",
    ]
    return merged[columns]


def build_homeostatic(breakdown: Path | None = None) -> pd.DataFrame:
    root = Path(breakdown) if breakdown is not None else BREAKDOWN_DIR
    pos = read_csv(root / "History(pos)_x_Rebound(neg)_overlap_n=15genes.csv")
    neg = read_csv(root / "History(neg)_x_Rebound(pos)_overlap_n=5genes.csv")
    if len(pos) != 15 or len(neg) != 5:
        raise ValueError(f"Homeostatic source counts were {len(pos)} and {len(neg)}")
    pos = pos.copy()
    pos["history_direction"] = "positive"
    pos["rebound_direction"] = "negative"
    pos["source_file"] = "History(pos)_x_Rebound(neg)_overlap_n=15genes.csv"
    neg = neg.copy()
    neg["history_direction"] = "negative"
    neg["rebound_direction"] = "positive"
    neg["source_file"] = "History(neg)_x_Rebound(pos)_overlap_n=5genes.csv"
    combined = pd.concat([pos, neg], ignore_index=True)
    combined["gene_set"] = "Homeostatic genes"
    combined["gene_set_definition"] = HOMEOSTATIC_DEFINITION
    combined["identity_status"] = "source"
    combined["identity_resolution_basis"] = ""
    mask = (combined["ext_gene"] == TRACHEALESS_SYMBOL) & (
        combined["flybase_gene_id"] == TRACHEALESS_FBGN
    )
    if int(mask.sum()) != 1:
        raise ValueError(
            "Expected exactly one trachealess source row to publication-resolve to Trhn"
        )
    combined.loc[mask, "ext_gene"] = TRHN_SYMBOL
    combined.loc[mask, "flybase_gene_id"] = TRHN_FBGN
    combined.loc[mask, "identity_status"] = "publication_resolved"
    combined.loc[mask, "identity_resolution_basis"] = PUBLICATION_RESOLUTION_BASIS
    if (combined["flybase_gene_id"] == TRACHEALESS_FBGN).any():
        raise ValueError("Generated homeostatic CSV still contains trachealess FBgn")
    if (combined["ext_gene"] == TRACHEALESS_SYMBOL).any():
        raise ValueError("Generated homeostatic CSV still contains trachealess symbol")
    keep = [
        "ext_gene",
        "flybase_gene_id",
        "gene_set",
        "gene_set_definition",
        "history_direction",
        "rebound_direction",
        "source_file",
        "identity_status",
        "identity_resolution_basis",
    ]
    extra = [
        col
        for col in combined.columns
        if col not in keep and col not in {"flybase_gene_id", "ext_gene"}
        and not col.startswith("History_")
        and not col.startswith("Rebound_")
    ]
    # Keep literature annotations without identity-like extras that could reintroduce the wrong ID.
    lit_cols = [
        col
        for col in combined.columns
        if col.startswith("History_") or col.startswith("Rebound_")
    ]
    return combined[keep + lit_cols + extra]


def _scope_frame(df: pd.DataFrame, threshold: str) -> pd.DataFrame:
    scoped = df.copy()
    rename = {
        "frequency": f"frequency_{threshold}_wake",
        "is_cycling": f"is_cycling_{threshold}_wake",
        "ZT_cyclers": f"ZT_cyclers_{threshold}_wake",
        "sleep_corr_exps": f"sleep_corr_exps_{threshold}_wake",
        "wake_corr_exps": f"wake_corr_exps_{threshold}_wake",
    }
    for column in EXPERIMENT_COLUMNS:
        rename[column] = f"{column}_{threshold}_wake"
    scoped = scoped.rename(columns=rename)
    keep = ["ext_gene", "flybase_gene_id", "source_file", *rename.values()]
    return scoped[[col for col in keep if col in scoped.columns]]


def build_csw_4plus(breakdown: Path | None = None) -> pd.DataFrame:
    root = Path(breakdown) if breakdown is not None else BREAKDOWN_DIR
    frames_fc05 = []
    frames_fc1 = []
    for name in CSW_4PLUS_SOURCES:
        df = read_csv(root / name)
        df = df.copy()
        df["source_file"] = name
        threshold = "FC0.5" if name.startswith("FC0.5") else "FC1"
        scoped = _scope_frame(df, threshold)
        if threshold == "FC0.5":
            frames_fc05.append(scoped)
        else:
            frames_fc1.append(scoped)
    fc05 = pd.concat(frames_fc05, ignore_index=True)
    fc1 = pd.concat(frames_fc1, ignore_index=True)
    if fc05["flybase_gene_id"].nunique() != 97:
        raise ValueError(
            f"FC0.5 4+ unique genes were {fc05['flybase_gene_id'].nunique()}, expected 97"
        )
    if set(fc1["flybase_gene_id"]) - set(fc05["flybase_gene_id"]):
        raise ValueError("FC1 frequency-4 genes are not a subset of FC0.5 4+")
    if len(fc1) != 7:
        raise ValueError(f"FC1 frequency-4 count was {len(fc1)}, expected 7")

    merged = fc05.merge(fc1, on=["ext_gene", "flybase_gene_id"], how="outer", suffixes=("", "_fc1"))
    if "source_file_fc1" in merged.columns:
        merged["source_files"] = [
            join_sorted([a, b])
            for a, b in zip(merged["source_file"].fillna(""), merged["source_file_fc1"].fillna(""))
        ]
        merged = merged.drop(columns=["source_file", "source_file_fc1"])
    else:
        merged["source_files"] = merged["source_file"]
        merged = merged.drop(columns=["source_file"])

    freq05 = pd.to_numeric(merged.get("frequency_FC0.5_wake"), errors="coerce")
    freq1 = pd.to_numeric(merged.get("frequency_FC1_wake"), errors="coerce")
    merged["qualifies_FC0.5_4plus"] = (freq05 >= 4).fillna(False).map(lambda v: "TRUE" if v else "FALSE")
    merged["qualifies_FC1_4plus"] = (freq1 >= 4).fillna(False).map(lambda v: "TRUE" if v else "FALSE")
    merged["qualifying_thresholds"] = [
        join_sorted(
            (["FC0.5"] if q05 == "TRUE" else []) + (["FC1"] if q1 == "TRUE" else [])
        )
        for q05, q1 in zip(merged["qualifies_FC0.5_4plus"], merged["qualifies_FC1_4plus"])
    ]
    if len(merged) != 97:
        raise ValueError(f"CSW 4+ union was {len(merged)}, expected 97")
    merged["gene_set"] = "CSW 4+ genes"
    merged["gene_set_definition"] = CSW_4PLUS_DEFINITION
    sleep_files = list(root.glob("FC*_Sleep_frequency_*.csv"))
    for path in sleep_files:
        sleep = read_csv(path)
        if pd.to_numeric(sleep["frequency"], errors="coerce").max() >= 4:
            raise ValueError(f"Sleep file {path.name} unexpectedly reaches frequency 4")
    front = [
        "ext_gene",
        "flybase_gene_id",
        "gene_set",
        "gene_set_definition",
        "frequency_FC0.5_wake",
        "frequency_FC1_wake",
        "qualifies_FC0.5_4plus",
        "qualifies_FC1_4plus",
        "qualifying_thresholds",
        "source_files",
    ]
    rest = [col for col in merged.columns if col not in front]
    return merged[front + rest].sort_values(["ext_gene", "flybase_gene_id"]).reset_index(drop=True)


def build_hlh(maps: FlyBaseMaps) -> pd.DataFrame:
    rows = []
    for symbol, expected_fbgn in HLH_PUBLICATION_GENES:
        resolved = resolve_symbol(symbol, maps)
        exact = str(resolved["exact_current_fbgn"] or "")
        if exact != expected_fbgn:
            raise ValueError(
                f"HLH {symbol} resolved to {exact!r}, expected {expected_fbgn}"
            )
        rows.append(
            {
                "ext_gene": symbol,
                "flybase_gene_id": expected_fbgn,
                "gene_set": "HLH genes",
                "gene_set_definition": HLH_DEFINITION,
                "publication_section": "Results",
                "publication_table": "Table S1",
                "motif": "CAGCTG E-box",
                "publication_evidence": (
                    "HOMER identified the CAGCTG E-box motif upstream of wake-induced "
                    "genes; the paper names this bHLH factor as a potential direct "
                    "upstream regulator of wake-induced gene expression."
                ),
                "identity_status": "publication_named",
            }
        )
    return pd.DataFrame(rows)


def audit_category_identities(
    frames: dict[str, pd.DataFrame], maps: FlyBaseMaps
) -> tuple[pd.DataFrame, list[str]]:
    rows = []
    unresolved: list[str] = []
    for set_name, df in frames.items():
        for _, rec in df.iterrows():
            symbol = str(rec["ext_gene"])
            fbgn = str(rec["flybase_gene_id"])
            resolved = resolve_symbol(symbol, maps)
            exact = str(resolved["exact_current_fbgn"] or "")
            candidates = list(resolved["candidate_fbgn_ids"])
            agrees = bool(exact and fbgn == exact)
            if not agrees and resolved["match_type"] == "unique_synonym":
                agrees = candidates == [fbgn]
            family = trh_family_flag(symbol)
            publication_conflict = family == "Trh"
            identity_status = str(rec.get("identity_status", "source"))
            proposed = "none"
            case_insensitive_current = case_insensitive_current_candidates(
                symbol, maps, candidates
            )
            if family == "Trh":
                proposed = "block_bare_Trh"
                unresolved.append(f"{set_name}:{symbol}:{fbgn}:bare_Trh")
            elif not fbgn.startswith("FBgn"):
                proposed = "block_unmapped"
                unresolved.append(f"{set_name}:{symbol}:{fbgn}:unmapped")
            elif identity_status == "publication_resolved" and symbol == TRHN_SYMBOL and fbgn == TRHN_FBGN:
                proposed = "none"
                agrees = True
            elif exact and fbgn != exact:
                proposed = "block_id_mismatch"
                unresolved.append(f"{set_name}:{symbol}:{fbgn}:mismatch_exact_{exact}")
            elif resolved["match_type"] == "ambiguous_synonym" and not exact:
                if case_insensitive_current == [fbgn]:
                    proposed = "none"
                    agrees = True
                    resolved["match_type"] = "case_insensitive_current_among_synonyms"
                else:
                    proposed = "block_ambiguous_synonym"
                    unresolved.append(f"{set_name}:{symbol}:{fbgn}:ambiguous")
            elif resolved["match_type"] == "unmapped":
                proposed = "block_unmapped_symbol"
                unresolved.append(f"{set_name}:{symbol}:{fbgn}:unmapped_symbol")
            rows.append(
                {
                    "set_name": set_name,
                    "ext_gene": symbol,
                    "flybase_gene_id": fbgn,
                    "exact_current_symbol_match": "TRUE" if resolved["exact_current_symbol_match"] else "FALSE",
                    "exact_current_fbgn": exact,
                    "candidate_fbgn_ids": ";".join(candidates),
                    "candidate_primary_symbols": ";".join(resolved["candidate_primary_symbols"]),
                    "candidate_count": resolved["candidate_count"],
                    "match_type": resolved["match_type"],
                    "final_id_agrees_with_exact_current": "TRUE" if agrees else "FALSE",
                    "trh_family_flag": family,
                    "publication_context_conflict": "TRUE" if publication_conflict else "FALSE",
                    "identity_status": identity_status,
                    "proposed_action": proposed,
                }
            )
    return pd.DataFrame(rows), unresolved


def write_named_csv(df: pd.DataFrame, directory: Path, stem: str) -> Path:
    n = len(df)
    path = directory / f"{stem}_n={n}genes.csv"
    write_csv(df, path)
    return path
