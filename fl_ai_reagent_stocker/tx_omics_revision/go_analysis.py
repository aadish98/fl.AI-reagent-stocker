"""Run DAVID GO via pinned aadish98/GO_Analysis GenerateGOChartReport functions.

The enrichment call is ``init_david_client``, ``load_input_table``, and
``fetch_go_report`` from the hash-pinned script. ``ProcessGOresults.process_csv_files``
is not invoked: it requires a ``primary_FBid`` mapping column that the bundled
GO category workbook does not have, and it crashes at the end of
``GenerateGOChartReport.main``.
"""

from __future__ import annotations

import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

from .constants import (
    DAVID_CATEGORIES,
    DAVID_EASE,
    DAVID_EMAIL_DEFAULT,
    DAVID_ENDPOINT,
    DAVID_FDR_PERCENT_MAX,
    DAVID_MIN_COUNT,
    DAVID_SPECIES,
    DAVID_WSDL,
    FBGN_RE,
    GO_ANALYSIS_DIR,
    GO_SCRIPT_BLOBS,
    MASTER_TABLES,
)
from .io_utils import git_blob_sha1, write_csv
from .keywords import parse_david_term, split_david_genes


def verify_go_script_hashes(go_dir: Path | None = None) -> dict[str, str]:
    root = Path(go_dir) if go_dir is not None else GO_ANALYSIS_DIR
    checked = {}
    for name, expected in GO_SCRIPT_BLOBS.items():
        path = root / name
        if not path.exists():
            raise FileNotFoundError(f"Pinned GO_Analysis script missing: {path}")
        digest = git_blob_sha1(path)
        if digest != expected:
            raise ValueError(
                f"{name} blob hash {digest} != pinned {expected}. Aborting."
            )
        checked[name] = digest
    return checked


GO_ANALYSIS_ENTRYPOINTS = (
    "init_david_client",
    "load_input_table",
    "fetch_go_report",
)


def _load_go_module(go_dir: Path):
    path = go_dir / "GenerateGOChartReport.py"
    go_dir_str = str(go_dir)
    if go_dir_str not in sys.path:
        sys.path.insert(0, go_dir_str)
    spec = importlib.util.spec_from_file_location("GenerateGOChartReport", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    missing = [
        name for name in GO_ANALYSIS_ENTRYPOINTS if not callable(getattr(module, name, None))
    ]
    if missing:
        raise ImportError(f"{path} is missing GO_Analysis entrypoints: {missing}")
    return module


def _fetch_go_report(module, client, flybase_ids, list_name) -> tuple[pd.DataFrame, object]:
    """Call GenerateGOChartReport.fetch_go_report, capturing addList coverage if possible."""
    add_raw = None
    original = None
    try:
        original = client.service.addList

        def _add_list(*args, **kwargs):
            nonlocal add_raw
            add_raw = original(*args, **kwargs)
            return add_raw

        client.service.addList = _add_list
    except Exception:  # noqa: BLE001 - suds clients are not always patchable
        original = None
    try:
        go_df = module.fetch_go_report(client, flybase_ids, list_name)
    except TypeError as exc:
        # Official loop is `for record in chart_report` and crashes if DAVID returns None.
        message = str(exc)
        if "NoneType" in message and "iter" in message:
            go_df = pd.DataFrame()
        else:
            raise
    finally:
        if original is not None:
            try:
                client.service.addList = original
            except Exception:  # noqa: BLE001
                pass
    if go_df is None:
        go_df = pd.DataFrame()
    return go_df, add_raw


def _add_list_coverage(raw) -> dict[str, str]:
    text = "" if raw is None else str(raw)
    mapped = ""
    submitted = ""
    match = re.search(r"(\d+)\s+out of\s+(\d+)", text)
    if match:
        mapped, submitted = match.group(1), match.group(2)
    return {
        "addList_raw": text,
        "mapped_count": mapped,
        "submitted_count": submitted,
    }


def parse_chart_records(go_df: pd.DataFrame, source_table: str, source_workbook: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    term_rows = []
    gene_rows = []
    for _, rec in go_df.iterrows():
        fdr_percent = float(rec["FDR"])
        if fdr_percent > DAVID_FDR_PERCENT_MAX:
            continue
        term_id, term_name = parse_david_term(rec["Term"])
        genes = split_david_genes(rec["Genes"])
        count = int(rec["Count"])
        if len(genes) != count:
            raise ValueError(
                f"{source_table} term {term_id} Count={count} but {len(genes)} gene tokens"
            )
        for gene in genes:
            if not re.match(FBGN_RE, gene):
                raise ValueError(f"{source_table} term {term_id} has non-FBgn token {gene!r}")
            gene_rows.append(
                {
                    "source_table": source_table,
                    "source_workbook": source_workbook,
                    "category": rec["Category"],
                    "term_id": term_id,
                    "term_name": term_name,
                    "flybase_gene_id": gene,
                }
            )
        term_rows.append(
            {
                "source_table": source_table,
                "source_workbook": source_workbook,
                "category": rec["Category"],
                "term_id": term_id,
                "term_name": term_name,
                "term_raw": rec["Term"],
                "count": count,
                "percent": rec["%"],
                "pvalue": rec["Pvalue"],
                "list_total": rec["List Total"],
                "pop_hits": rec["Pop Hits"],
                "pop_total": rec["Pop Total"],
                "fold_enrichment": rec["Fold Enrichment"],
                "bonferroni": rec["Bonferroni"],
                "benjamini": rec["Benjamini"],
                "fdr_percent": fdr_percent,
                "fdr_q": fdr_percent / 100.0,
                "genes": ",".join(genes),
            }
        )
    return pd.DataFrame(term_rows), pd.DataFrame(gene_rows)


def run_david_go(
    input_dir: Path,
    raw_dir: Path,
    processed_dir: Path,
    evidence_dir: Path,
    submitted_ids: dict[str, set[str]],
    email: str = DAVID_EMAIL_DEFAULT,
    go_dir: Path | None = None,
) -> dict[str, object]:
    go_root = Path(go_dir) if go_dir is not None else GO_ANALYSIS_DIR
    verify_go_script_hashes(go_root)
    module = _load_go_module(go_root)
    if raw_dir.exists():
        shutil.rmtree(raw_dir)
    if processed_dir.exists():
        shutil.rmtree(processed_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    client = module.init_david_client(email)

    outcomes = []
    term_frames = []
    gene_frames = []
    input_files = sorted(input_dir.glob("*.csv"))
    expected = [meta["filename"] for meta in MASTER_TABLES.values()]
    found = [path.name for path in input_files]
    if found != sorted(expected):
        raise ValueError(f"GO input files {found} != {sorted(expected)}")

    for path in input_files:
        source_table = next(
            name for name, meta in MASTER_TABLES.items() if meta["filename"] == path.name
        )
        flybase_ids = module.load_input_table(path)
        outcome = {
            "source_table": source_table,
            "input_file": path.name,
            "n_ids": int(len(flybase_ids)),
            "status": "service_failure",
            "go_analysis_function": "GenerateGOChartReport.fetch_go_report",
        }
        try:
            go_df, add_raw = _fetch_go_report(module, client, flybase_ids, path.stem)
            outcome.update(_add_list_coverage(add_raw))
            raw_path = raw_dir / f"{path.stem}_GO_Analysis.xlsx"
            processed_path = processed_dir / f"{path.stem}_GO_Analysis.xlsx"
            if go_df.empty:
                outcome["status"] = "no_fdr_passing_terms"
                pd.DataFrame(
                    columns=[
                        "Category",
                        "Term",
                        "Count",
                        "%",
                        "Pvalue",
                        "Genes",
                        "List Total",
                        "Pop Hits",
                        "Pop Total",
                        "Fold Enrichment",
                        "Bonferroni",
                        "Benjamini",
                        "FDR",
                    ]
                ).to_excel(raw_path, index=False)
                outcome["raw_path"] = str(raw_path)
            else:
                go_df.to_excel(raw_path, index=False)
                saved = pd.read_excel(raw_path)
                workbook = path.stem + "_GO_Analysis.xlsx"
                terms, genes = parse_chart_records(saved, source_table, workbook)
                allowed = submitted_ids[source_table]
                extra = set(genes["flybase_gene_id"]) - allowed
                if extra:
                    raise ValueError(f"{source_table} GO hits not in submitted IDs: {sorted(extra)[:5]}")
                term_frames.append(terms)
                gene_frames.append(genes)
                with pd.ExcelWriter(processed_path, engine="openpyxl") as writer:
                    for category, group in terms.groupby("category"):
                        group.to_excel(writer, sheet_name=str(category)[:31], index=False)
                outcome["status"] = "success"
                outcome["n_terms"] = int(len(terms))
                outcome["raw_path"] = str(raw_path)
                outcome["processed_path"] = str(processed_path)
        except Exception as exc:  # noqa: BLE001 - record and abort after writing status
            outcome["status"] = "service_failure"
            outcome["error"] = str(exc)
            outcomes.append(outcome)
            (raw_dir / "outcomes.json").write_text(
                json.dumps(outcomes, indent=2), encoding="utf-8"
            )
            raise
        outcomes.append(outcome)

    terms_all = pd.concat(term_frames, ignore_index=True) if term_frames else pd.DataFrame()
    genes_all = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
    write_csv(terms_all, evidence_dir / "GO_Term_Results.csv")
    write_csv(genes_all, evidence_dir / "GO_Term_Genes.csv")
    payload = {
        "go_analysis_repo": "https://github.com/aadish98/GO_Analysis",
        "go_analysis_entrypoints": list(GO_ANALYSIS_ENTRYPOINTS),
        "skipped_go_analysis_step": "ProcessGOresults.process_csv_files",
        "authenticate": "GenerateGOChartReport.init_david_client",
        "species": DAVID_SPECIES,
        "ease": DAVID_EASE,
        "min_count": DAVID_MIN_COUNT,
        "fdr_percent_max": DAVID_FDR_PERCENT_MAX,
        "categories": DAVID_CATEGORIES,
        "wsdl": DAVID_WSDL,
        "endpoint": DAVID_ENDPOINT,
        "email": email,
        "outcomes": outcomes,
    }
    (raw_dir / "outcomes.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    failed = [row for row in outcomes if row["status"] not in {"success", "no_fdr_passing_terms"}]
    if failed:
        raise RuntimeError(f"DAVID outcomes failed: {failed}")
    return payload
