"""Orchestrate Tx-Omics Revision construction, GO, audit, overlap, and plots."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import pandas as pd

from fl_ai_reagent_stocker.config import REPO_ROOT
from fl_ai_reagent_stocker.cli import _discover_input_csvs, _validate_gene_columns

from .categories import (
    audit_category_identities,
    build_csw_4plus,
    build_hlh,
    build_homeostatic,
    build_mechanistic,
    load_master_tables,
    write_csw_observations,
    write_named_csv,
)
from .constants import (
    AUDIT_ROOT,
    BREAKDOWN_DIR,
    BUCKET_TO_SET_NAME,
    BUCKET_TO_STEM,
    DAVID_CATEGORIES,
    DAVID_EASE,
    DAVID_EMAIL_DEFAULT,
    DAVID_ENDPOINT,
    DAVID_FDR_PERCENT_MAX,
    DAVID_MIN_COUNT,
    DAVID_SPECIES,
    DAVID_WSDL,
    FETCH_FBGN_SCRIPT,
    FLYBASE_SYNONYM_PATH,
    GO_ANALYSIS_DIR,
    GO_ANALYSIS_REMOTE_COMMIT,
    MASTER_TABLES,
    SET_NAMES,
    STOCKER_DIR,
    TX_INPUT_DIR,
    audit_paths,
)
from .flybase import FlyBaseMaps, case_insensitive_current_candidates, load_flybase_maps, resolve_symbol, trh_family_flag
from .go_analysis import run_david_go, verify_go_script_hashes
from .io_utils import atomic_write_json, ensure_dirs, git_status, read_csv, sha256_file, write_csv
from .overlap import membership_table, overlap_report_markdown, write_overlap_tables
from .pathways import approve_proposed_terms, build_pathway_tables, build_term_review, reconstruct_summary_ok
from .plots import plot_all


STAGES = [
    "preflight",
    "categories",
    "category_audit",
    "go",
    "review_terms",
    "pathways",
    "pathway_audit",
    "overlap",
    "plots",
    "docs",
]


def _package_versions() -> dict[str, str]:
    names = ["pandas", "openpyxl", "matplotlib", "matplotlib-venn", "UpSetPlot", "suds"]
    versions = {"python": sys.version.split()[0]}
    for name in names:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = "missing"
    return versions


def _input_hashes() -> dict[str, str]:
    hashes = {}
    for path in sorted(BREAKDOWN_DIR.glob("*.csv")):
        hashes[f"breakdown/{path.name}"] = sha256_file(path)
    for meta in MASTER_TABLES.values():
        path = TX_INPUT_DIR / meta["filename"]
        hashes[f"masters/{path.name}"] = sha256_file(path)
    hashes[f"flybase/{FLYBASE_SYNONYM_PATH.name}"] = sha256_file(FLYBASE_SYNONYM_PATH)
    return hashes


def run_preflight(paths: dict[str, Path], manifest: dict) -> None:
    ensure_dirs(
        [
            paths["staging_categories"],
            paths["go_input"],
            paths["go_raw"],
            paths["go_processed"],
            paths["go_review"],
            paths["pathways_preaudit"],
            paths["pathways_auditrun"],
            paths["pathways_approved"],
            paths["evidence"],
            paths["audits"],
            paths["figures"],
            paths["figure_data"],
            STOCKER_DIR,
        ]
    )
    go_hashes = verify_go_script_hashes(GO_ANALYSIS_DIR)
    observations = load_master_tables(TX_INPUT_DIR)
    write_csw_observations(observations, paths["evidence"] / "CSW_Observations.csv")
    for meta in MASTER_TABLES.values():
        shutil.copy2(TX_INPUT_DIR / meta["filename"], paths["go_input"] / meta["filename"])
    manifest["git"] = git_status(REPO_ROOT)
    manifest["go_script_blobs"] = go_hashes
    manifest["go_analysis_remote_commit_pin"] = GO_ANALYSIS_REMOTE_COMMIT
    manifest["flybase_mapping"] = {
        "path": str(FLYBASE_SYNONYM_PATH),
        "sha256": sha256_file(FLYBASE_SYNONYM_PATH),
        "release": "fb_2026_01",
    }
    manifest["david"] = {
        "email": DAVID_EMAIL_DEFAULT,
        "wsdl": DAVID_WSDL,
        "endpoint": DAVID_ENDPOINT,
        "species": DAVID_SPECIES,
        "categories": DAVID_CATEGORIES,
        "ease": DAVID_EASE,
        "min_count": DAVID_MIN_COUNT,
        "fdr_percent_max": DAVID_FDR_PERCENT_MAX,
    }
    manifest["package_versions"] = _package_versions()
    manifest["input_sha256"] = _input_hashes()
    manifest["csw_observations"] = {
        "n_rows": int(len(observations)),
        "n_unique_fbgn": int(observations["flybase_gene_id"].nunique()),
    }


def run_categories(paths: dict[str, Path], maps: FlyBaseMaps, manifest: dict) -> dict[str, Path]:
    mechanistic = build_mechanistic(BREAKDOWN_DIR)
    homeostatic = build_homeostatic(BREAKDOWN_DIR)
    csw4 = build_csw_4plus(BREAKDOWN_DIR)
    hlh = build_hlh(maps)
    written = {
        "Mechanistic": write_named_csv(mechanistic, paths["staging_categories"], "01_Mechanistic"),
        "Homeostatic History/Rebound": write_named_csv(
            homeostatic, paths["staging_categories"], "02_Homeostatic_HistoryxRebound"
        ),
        "CSW 4+": write_named_csv(csw4, paths["staging_categories"], "03_CSW_4plus"),
        "HLH Upstream Regulators": write_named_csv(
            hlh, paths["staging_categories"], "04_HLH_upstream_regulators"
        ),
    }
    manifest["category_files"] = {name: str(path) for name, path in written.items()}
    manifest["category_counts"] = {name: int(len(read_csv(path))) for name, path in written.items()}
    return written


def run_category_audit(
    paths: dict[str, Path],
    maps: FlyBaseMaps,
    category_files: dict[str, Path],
    manifest: dict,
) -> None:
    frames = {name: read_csv(path) for name, path in category_files.items()}
    audit, unresolved = audit_category_identities(frames, maps)
    dest = paths["audits"] / "category_identity_audit.csv"
    write_csv(audit, dest)
    manifest["category_identity_audit"] = str(dest)
    manifest["category_unresolved"] = unresolved
    if unresolved:
        raise RuntimeError(
            "Unresolved category identity discrepancies:\n" + "\n".join(unresolved)
        )


def run_go(paths: dict[str, Path], manifest: dict) -> None:
    observations = read_csv(paths["evidence"] / "CSW_Observations.csv")
    submitted = {
        table: set(group["flybase_gene_id"])
        for table, group in observations.groupby("source_table")
    }
    payload = run_david_go(
        input_dir=paths["go_input"],
        raw_dir=paths["go_raw"],
        processed_dir=paths["go_processed"],
        evidence_dir=paths["evidence"],
        submitted_ids=submitted,
        email=DAVID_EMAIL_DEFAULT,
        go_dir=GO_ANALYSIS_DIR,
    )
    manifest["david_run"] = payload


def run_review_terms(paths: dict[str, Path], manifest: dict, approve_proposed: bool) -> Path:
    terms = read_csv(paths["evidence"] / "GO_Term_Results.csv")
    genes = read_csv(paths["evidence"] / "GO_Term_Genes.csv")
    review = build_term_review(terms, genes)
    review_path = paths["go_review"] / "GO_Term_Bucket_Review.csv"
    write_csv(review, review_path)
    manifest["go_term_review"] = str(review_path)
    if not approve_proposed:
        raise RuntimeError(
            f"GO term review written to {review_path}. Re-run with "
            "--approve-proposed-terms after inspecting the file."
        )
    approved = approve_proposed_terms(review, "include_all_matched_buckets")
    pending = approved[approved["approval_status"] == "pending_conflict"]
    if not pending.empty:
        raise RuntimeError("Conflict terms remain unapproved.")
    approved_path = paths["go_review"] / "GO_Term_Bucket_Review.approved.csv"
    write_csv(approved, approved_path)
    manifest["go_term_review_approved"] = str(approved_path)
    manifest["go_term_review_approved_sha256"] = sha256_file(approved_path)
    manifest["term_approval_policy"] = "proposed_includes_plus_conflict_all_matched_buckets"
    return approved_path


def _category_id_sets(category_files: dict[str, Path]) -> tuple[dict[str, set[str]], dict[str, str]]:
    ids = {}
    csw4_thresholds = {}
    for name, path in category_files.items():
        df = read_csv(path)
        ids[name] = set(df["flybase_gene_id"])
        if name == "CSW 4+" and "qualifying_thresholds" in df.columns:
            csw4_thresholds = dict(zip(df["flybase_gene_id"], df["qualifying_thresholds"]))
    return ids, csw4_thresholds


def run_pathways(paths: dict[str, Path], category_files: dict[str, Path], manifest: dict) -> dict[str, Path]:
    approved = read_csv(paths["go_review"] / "GO_Term_Bucket_Review.approved.csv")
    genes = read_csv(paths["evidence"] / "GO_Term_Genes.csv")
    observations = read_csv(paths["evidence"] / "CSW_Observations.csv")
    ids, csw4_thresholds = _category_id_sets(category_files)
    tables = build_pathway_tables(approved, genes, observations, ids, csw4_thresholds)
    written = {}
    for bucket, df in tables.items():
        if not reconstruct_summary_ok(df, genes, approved):
            raise ValueError(f"Pathway summary for {bucket} does not reconstruct from evidence")
        stem = BUCKET_TO_STEM[bucket]
        path = write_named_csv(df, paths["pathways_preaudit"], stem)
        written[BUCKET_TO_SET_NAME[bucket]] = path
    manifest["pathway_preaudit_files"] = {name: str(path) for name, path in written.items()}
    return written


def _snapshot_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def run_pathway_audit(
    paths: dict[str, Path],
    maps: FlyBaseMaps,
    pathway_files: dict[str, Path],
    manifest: dict,
    promote_if_clean: bool,
) -> dict[str, Path]:
    audit_dir = paths["pathways_auditrun"]
    if audit_dir.exists():
        shutil.rmtree(audit_dir)
    audit_dir.mkdir(parents=True)
    snapshots = {}
    for src in pathway_files.values():
        dest = audit_dir / src.name
        shutil.copy2(src, dest)
        snapshots[dest.name] = _snapshot_csv(src)
    present = sorted(p.name for p in audit_dir.glob("*.csv"))
    expected = sorted(src.name for src in pathway_files.values())
    if present != expected:
        raise RuntimeError(f"AuditRun CSVs {present} != {expected}")
    cmd = [sys.executable, str(FETCH_FBGN_SCRIPT), str(audit_dir), "ext_gene"]
    subprocess.run(cmd, check=True)
    audit_rows = []
    blockers = []
    restored_paths = {}
    for name in expected:
        before = snapshots[name]
        after = _snapshot_csv(audit_dir / name)
        if len(before) != len(after):
            blockers.append(f"{name}: row count {len(before)} -> {len(after)}")
        identity_cols = {"flybase_gene_id", "primary_symbol", "corrected_ext_gene"}
        restored = after.copy()
        for col in before.columns:
            if col not in identity_cols and col in restored.columns:
                restored[col] = before[col].values
            elif col not in identity_cols:
                restored[col] = before[col].values
        # Keep converter identity columns, restore everything else from snapshot.
        ordered = [
            col
            for col in [
                "ext_gene",
                "corrected_ext_gene",
                "primary_symbol",
                "source_flybase_gene_id",
                "flybase_gene_id",
            ]
            if col in restored.columns
        ]
        rest = [col for col in before.columns if col not in ordered]
        extra = [col for col in restored.columns if col not in ordered + rest]
        restored = restored[ordered + rest + extra]
        restored_path = paths["pathways_approved"] / name
        if len(before) == len(after):
            for idx in range(len(before)):
                symbol = str(before.iloc[idx]["ext_gene"])
                source_id = str(before.iloc[idx]["flybase_gene_id"])
                selected = str(after.iloc[idx]["flybase_gene_id"])
                resolved = resolve_symbol(symbol, maps)
                metadata_equal = True
                for col in before.columns:
                    if col in identity_cols:
                        continue
                    if col in after.columns and str(after.iloc[idx][col]) != str(before.iloc[idx][col]):
                        # converter fillna may have changed blanks; restored file will fix this
                        if str(before.iloc[idx][col]) == "" and str(after.iloc[idx][col]) == "-":
                            continue
                        metadata_equal = False
                status = "agree"
                if selected in {"", "-"}:
                    status = "unmapped"
                    blockers.append(f"{name}:{symbol}:unmapped")
                elif selected != source_id:
                    status = "disagree"
                    blockers.append(f"{name}:{symbol}:{source_id}->{selected}")
                elif resolved["match_type"] == "ambiguous_synonym" and not resolved["exact_current_symbol_match"]:
                    case_hits = case_insensitive_current_candidates(
                        symbol, maps, list(resolved["candidate_fbgn_ids"])
                    )
                    if case_hits == [source_id] and selected == source_id:
                        status = "agree_case_insensitive_current"
                    else:
                        status = "ambiguous_non_exact"
                        blockers.append(f"{name}:{symbol}:ambiguous")
                family = trh_family_flag(symbol)
                if family == "Trh":
                    blockers.append(f"{name}:{symbol}:bare_Trh")
                audit_rows.append(
                    {
                        "file": name,
                        "ext_gene": symbol,
                        "match_type": resolved["match_type"],
                        "matched_value": symbol,
                        "candidate_fbgn_ids": ";".join(resolved["candidate_fbgn_ids"]),
                        "candidate_primary_symbols": ";".join(resolved["candidate_primary_symbols"]),
                        "candidate_count": resolved["candidate_count"],
                        "converter_selected_id": selected,
                        "source_id": source_id,
                        "status": status,
                        "trh_family_flag": family,
                        "non_id_metadata_equal": "TRUE" if metadata_equal else "FALSE",
                    }
                )
        write_csv(restored, restored_path)
        restored_paths[name] = restored_path
    audit_path = paths["audits"] / "pathway_fbgn_resolution_audit.csv"
    write_csv(pd.DataFrame(audit_rows), audit_path)
    manifest["pathway_fbgn_audit"] = str(audit_path)
    manifest["pathway_audit_blockers"] = blockers
    if blockers:
        raise RuntimeError(
            "Pathway FlyBase audit blockers:\n" + "\n".join(blockers)
        )
    if not promote_if_clean:
        raise RuntimeError("Pathway audit is clean. Re-run with --promote-if-audit-clean to promote.")
    approved = {}
    for set_name, src in pathway_files.items():
        approved[set_name] = restored_paths[src.name]
    manifest["pathway_approved_files"] = {name: str(path) for name, path in approved.items()}
    return approved


def promote_stocker_inputs(
    category_files: dict[str, Path], pathway_files: dict[str, Path]
) -> list[Path]:
    if STOCKER_DIR.exists():
        for path in STOCKER_DIR.glob("*.csv"):
            path.unlink()
    STOCKER_DIR.mkdir(parents=True, exist_ok=True)
    copied = []
    for src in list(category_files.values()) + list(pathway_files.values()):
        dest = STOCKER_DIR / src.name
        shutil.copy2(src, dest)
        copied.append(dest)
    discovered = _discover_input_csvs(STOCKER_DIR)
    if sorted(p.name for p in discovered) != sorted(p.name for p in copied):
        raise RuntimeError("Stocker discovery did not find exactly the seven promoted CSVs")
    errors = _validate_gene_columns(discovered, "flybase_gene_id", "ext_gene")
    if errors:
        raise RuntimeError("Stocker column validation failed: " + "; ".join(errors))
    return copied


def run_overlap(paths: dict[str, Path], sets: dict[str, pd.DataFrame], manifest: dict) -> pd.DataFrame:
    membership = membership_table(sets)
    write_overlap_tables(membership, paths["evidence"])
    report = overlap_report_markdown(membership)
    paths["overlap_report"].write_text(report, encoding="utf-8")
    manifest["overlap_unique_genes"] = int(len(membership))
    return membership


def write_readme(paths: dict[str, Path], manifest: dict) -> None:
    text = f"""# Tx-Omics Revision

Publication-defined gene sets for stocker input, plus FDR-filtered DAVID pathway-hit lists.

## Stocker inputs

`data/gene_sets/Tx-Omics_Revision/` contains exactly seven CSVs. Each has `ext_gene` and `flybase_gene_id`. Analysis, audit, GO, overlap, and figure files are kept out of that directory so stocker recursive discovery cannot ingest them.

| File | Definition |
| --- | --- |
| `01_Mechanistic_n=6genes.csv` | Six publication genes with in-vivo homeostatic support |
| `02_Homeostatic_HistoryxRebound_n=20genes.csv` | 15 History+/Rebound− plus 5 History−/Rebound+ |
| `03_CSW_4plus_n=97genes.csv` | Wake genes at frequency ≥4 (FC0.5 and/or FC1) |
| `04_HLH_upstream_regulators_n=4genes.csv` | Publication-named bHLH upstream regulators |
| `05_CSW_Ribosomal_Translation_n=99genes.csv` | CSW genes hit by approved ribosomal/translation terms |
| `06_CSW_Mitochondrial_Metabolism_n=167genes.csv` | CSW genes hit by approved mitochondrial/metabolism terms |
| `07_CSW_Immune_n=5genes.csv` | CSW genes hit by approved immune terms |

Run later with:

```bash
python -m fl_ai_reagent_stocker run ./data/gene_sets/Tx-Omics_Revision/ \\
  --config ./data/config/stock_split_config_priority_example.json
```

The stocker pipeline was not run as part of this revision. Compatibility was checked with the same recursive CSV discovery and `ext_gene`/`flybase_gene_id` column validation the CLI uses.

## Sources

Do not modify these inputs. They were hashed into `run_manifest.json`.

- CSW masters (`private/lab-materials/Tx-Omics_FollowUp/Data/Transcriptomics_Input/`): Wake FC0.5 647, Wake FC1 98, Sleep FC0.5 141, Sleep FC1 16; 902 observations and 791 unique FBgns; every frequency ≥2.
- Category sources (`data/gene_sets/Tx-Omics-FollowUp_v3/Breakdown/`): homeostatic 6-gene list, mechanistic screen ranks 1–6, history×rebound overlap tables, and CSW frequency-4/5/6 wake tables.
- HLH membership is from Rosensweig–Shah 2026 Results / Table S1, not from CSW-frequency tables.
- FlyBase synonym table: `data/flybase/genes/fb_synonym_fb_2026_01.tsv.gz` (release `fb_2026_01`).
- GO analysis uses [`aadish98/GO_Analysis`](https://github.com/aadish98/GO_Analysis) `GenerateGOChartReport.py` functions `init_david_client`, `load_input_table`, and `fetch_go_report` (local copy in `private/lab-materials/Tx-Omics_FollowUp/`, pinned to git blobs `GenerateGOChartReport.py` `8430bf1fb0209ce2ff4c956750672d10981b3e34`, `ProcessGOresults.py` `e8287d99e6163f111cdfa2dbfaf7832f520c9b80`, `VisualizeGOResults.py` `56ea2c11817fa31aa589e76712ee27ed34065288`, remote commit `95501af480473d33f906d0c72714ef7613ace722`). `GenerateGOChartReport.main` is not used because it always calls `ProcessGOresults.process_csv_files`, which crashes on the bundled mapping workbook (no `primary_FBid` column). DAVID parameters and FDR ≤ 10 filtering therefore come from `fetch_go_report`, not a reimplemented SOAP client.

## Category definitions

- Mechanistic: unc79, SIFa, rumpel, AstA-R2, Trhn, RpL23, joined to ranks 1–6 of the mechanistic screen CSV.
- Homeostatic history/rebound: 15 History+/Rebound− plus 5 History−/Rebound+. The source row labeled `trh`/`FBgn0262139` was normalized to `Trhn`/`FBgn0035187` from Rosensweig–Shah 2026 Results, Discussion, and Figure 14 (tryptophan hydroxylase, not trachealess). The generated CSV does not retain the wrong FBgn.
- CSW 4+: FC0.5 wake frequency 4/5/6 union FC1 wake frequency 4. All 7 FC1 frequency-4 genes are already in the FC0.5 4+ set. Sleep never reaches frequency 4. Frequency, correlation, and cycling fields are stored per threshold (`frequency_FC0.5_wake`, `frequency_FC1_wake`, `wake_corr_exps_FC0.5_wake`, …), never as a collapsed maximum.
- HLH upstream regulators: bigmax, HLH3B, E(spl)m3-HLH, E(spl)mbeta-HLH. Columns are publication provenance only (`publication_section=Results`, `publication_table=Table S1`, `motif=CAGCTG E-box`).

History/rebound correlation coefficients and time bins are not in this repository and were not fabricated.

## Pathway lists

Pathway genes are CSW genes hit by FDR-filtered DAVID terms whose names match the approved keywords. They are not expanded to all annotated genes.

DAVID parameters: email `aadishms@umich.edu`, species 7227, categories {DAVID_CATEGORIES}, EASE {DAVID_EASE}, minimum count {DAVID_MIN_COUNT}, adaptive FDR ≤ {DAVID_FDR_PERCENT_MAX}%. `FBGN########` tokens are normalized to `FBgn`. Sleep FC1 produced no FDR-passing terms; that is an explicit recorded outcome, not a missing file.

### Term-approval policy

`--approve-proposed-terms` accepted proposed single-bucket includes. The three conflict terms (mitochondrial translation; mitochondrial large ribosomal subunit; mitochondrial small ribosomal subunit) were assigned to every matched bucket (`include_all_matched_buckets`). Adjacent unmatched terms remain in the review file as excluded/reviewable rows. The approved review hash is stored in `run_manifest.json`.

## Identity audit

Category and pathway FlyBase audits used the pinned 2026_01 synonym table. Exact `Trhn` resolves to `FBgn0035187`; exact `trh` is trachealess (`FBgn0262139`); bare `Trh` is ambiguous and is never auto-approved.

Accepted without changing IDs:

- Publication-resolved Homeostatic `Trhn` / `FBgn0035187`.
- `RpL37a` / `FBgn0261608` and `arg` / `FBgn0023535`: synonym-ambiguous, but exactly one candidate has a case-insensitive current symbol match that equals the source ID (`agree_case_insensitive_current`). Current symbols are `RpL37A` and `Arg`.
- Unique-synonym CSW 4+ symbols (several `CG*` genes and `Rpb11`) whose single synonym candidate equals the source ID.

No generated identity field contains `FBgn0262139`. Converter `fillna("-")` was restored from pre-audit snapshots for non-ID columns.

## Overlap

Membership is by FBgn ID only. See `Overlap_Report.md` for named genes and set membership. Figures and `FigureData/` CSVs reconstruct from `Evidence/Set_Membership.csv`.

## Caveats and local-only artifacts

- Generated data under `data/gene_sets/`, `audit_outputs/`, and sources under `private/` are gitignored.
- Working files live in `audit_outputs/Tx-Omics_Revision/` (`Staging/`, `GO/`, `Pathways/`, `Evidence/`, `Audits/`, `Figures/`, `FigureData/`, `Overlap_Report.md`, this README, `run_manifest.json`).
- Re-running DAVID requires network access to the DAVID web service.

## Rerun

```bash
python scripts/build_tx_omics_revision_sets.py --through all \\
  --approve-proposed-terms --promote-if-audit-clean
python scripts/run_tx_omics_revision_go.py
python scripts/plot_tx_omics_revision.py
```
"""
    paths["readme"].write_text(text, encoding="utf-8")
    manifest["readme"] = str(paths["readme"])


def run_pipeline(
    through: str = "all",
    approve_proposed_terms: bool = False,
    promote_if_audit_clean: bool = False,
    audit_root: Path | None = None,
) -> dict:
    if through != "all" and through not in STAGES:
        raise ValueError(f"Unknown stage {through}")
    stop_at = STAGES[-1] if through == "all" else through
    paths = audit_paths(audit_root)
    manifest = {
        "started_at": datetime.now(timezone.utc).isoformat(),
        "through": stop_at,
        "commands": list(sys.argv),
    }
    maps = None
    category_files = None
    pathway_files = None
    approved_pathways = None
    membership = None
    sets = None

    def _done(stage: str) -> bool:
        return STAGES.index(stage) <= STAGES.index(stop_at)

    if _done("preflight"):
        run_preflight(paths, manifest)
    if _done("categories"):
        maps = load_flybase_maps(FLYBASE_SYNONYM_PATH)
        category_files = run_categories(paths, maps, manifest)
    if _done("category_audit"):
        if maps is None:
            maps = load_flybase_maps(FLYBASE_SYNONYM_PATH)
        if category_files is None:
            category_files = {
                "Mechanistic": next(paths["staging_categories"].glob("01_Mechanistic_*.csv")),
                "Homeostatic History/Rebound": next(
                    paths["staging_categories"].glob("02_Homeostatic_*.csv")
                ),
                "CSW 4+": next(paths["staging_categories"].glob("03_CSW_4plus_*.csv")),
                "HLH Upstream Regulators": next(
                    paths["staging_categories"].glob("04_HLH_*.csv")
                ),
            }
        run_category_audit(paths, maps, category_files, manifest)
    if _done("go"):
        run_go(paths, manifest)
    if _done("review_terms"):
        run_review_terms(paths, manifest, approve_proposed_terms)
    if _done("pathways"):
        if category_files is None:
            category_files = {
                "Mechanistic": next(paths["staging_categories"].glob("01_Mechanistic_*.csv")),
                "Homeostatic History/Rebound": next(
                    paths["staging_categories"].glob("02_Homeostatic_*.csv")
                ),
                "CSW 4+": next(paths["staging_categories"].glob("03_CSW_4plus_*.csv")),
                "HLH Upstream Regulators": next(
                    paths["staging_categories"].glob("04_HLH_*.csv")
                ),
            }
        pathway_files = run_pathways(paths, category_files, manifest)
    if _done("pathway_audit"):
        if maps is None:
            maps = load_flybase_maps(FLYBASE_SYNONYM_PATH)
        if pathway_files is None:
            pathway_files = {
                BUCKET_TO_SET_NAME[bucket]: next(
                    paths["pathways_preaudit"].glob(f"{stem}_*.csv")
                )
                for bucket, stem in BUCKET_TO_STEM.items()
            }
        approved_pathways = run_pathway_audit(
            paths, maps, pathway_files, manifest, promote_if_audit_clean
        )
        promote_stocker_inputs(category_files, approved_pathways)
    if _done("overlap"):
        if category_files is None or approved_pathways is None:
            raise RuntimeError("overlap requires promoted category and pathway files")
        sets = {name: read_csv(path) for name, path in category_files.items()}
        sets.update({name: read_csv(path) for name, path in approved_pathways.items()})
        membership = run_overlap(paths, sets, manifest)
    if _done("plots"):
        if membership is None or sets is None:
            raise RuntimeError("plots require overlap tables")
        approved = read_csv(paths["go_review"] / "GO_Term_Bucket_Review.approved.csv")
        overlapping = read_csv(paths["evidence"] / "Overlapping_Genes.csv")
        pairs = read_csv(paths["evidence"] / "Pairwise_Overlap.csv")
        pathway_sets = {name: sets[BUCKET_TO_SET_NAME[name]] for name in BUCKET_TO_STEM}
        plot_all(
            pathway_sets,
            approved,
            membership,
            overlapping,
            pairs,
            paths["figures"],
            paths["figure_data"],
        )
        manifest["figures"] = [str(path) for path in sorted(paths["figures"].glob("*.png"))]
    if _done("docs"):
        write_readme(paths, manifest)
        manifest["finished_at"] = datetime.now(timezone.utc).isoformat()
        atomic_write_json(paths["manifest"], manifest)
    else:
        atomic_write_json(paths["manifest"], manifest)
    return manifest
