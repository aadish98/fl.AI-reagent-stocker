#!/usr/bin/env python3
"""
Download the latest FlyBase TSV/TSV.GZ report families used by this repo.

The script discovers files from the official FlyBase current-release index under:
    https://s3ftp.flybase.org/releases/current/precomputed_files/

It downloads only the report families this repository uses, saves them into the
expected local data directories, and removes older versions of the same prefix
only after the latest file is present and non-empty. After all manifest families
have been processed, it also removes orphan TSV/TSV.GZ files from the managed
FlyBase directories so stale helper files do not accumulate.

Usage:
    python scripts/refresh_flybase_data.py
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen


REPO_ROOT = Path(__file__).resolve().parents[1]
FLYBASE_ROOT = REPO_ROOT / "data" / "flybase"
ALLELES_AND_STOCKS_DIR = FLYBASE_ROOT / "alleles_and_stocks"
REFERENCES_DIR = FLYBASE_ROOT / "references"
GENES_DIR = FLYBASE_ROOT / "genes"
TRANSGENIC_CONSTRUCTS_DIR = FLYBASE_ROOT / "transgenic_constructs"

CURRENT_PRECOMPUTED_BASE = "https://s3ftp.flybase.org/releases/current/precomputed_files/"


@dataclass(frozen=True)
class ReportFamily:
    name: str
    prefix: str
    remote_subdir: str
    destination_dir: Path
    required: bool = True
    preferred_extensions: tuple[str, ...] = (".tsv.gz", ".tsv")


MANIFEST: List[ReportFamily] = [
    ReportFamily(
        name="stocks",
        prefix="stocks",
        remote_subdir="stocks",
        destination_dir=ALLELES_AND_STOCKS_DIR,
    ),
    ReportFamily(
        name="fbal_to_fbgn",
        prefix="fbal_to_fbgn",
        remote_subdir="alleles",
        destination_dir=ALLELES_AND_STOCKS_DIR,
    ),
    ReportFamily(
        name="dmel_classical_and_insertion_allele_descriptions",
        prefix="dmel_classical_and_insertion_allele_descriptions",
        remote_subdir="alleles",
        destination_dir=ALLELES_AND_STOCKS_DIR,
    ),
    ReportFamily(
        name="genotype_phenotype_data",
        prefix="genotype_phenotype_data",
        remote_subdir="alleles",
        destination_dir=ALLELES_AND_STOCKS_DIR,
    ),
    ReportFamily(
        name="transgenic_construct_descriptions",
        prefix="transgenic_construct_descriptions",
        remote_subdir="transposons",
        destination_dir=TRANSGENIC_CONSTRUCTS_DIR,
    ),
    ReportFamily(
        name="entity_publication",
        prefix="entity_publication",
        remote_subdir="references",
        destination_dir=REFERENCES_DIR,
    ),
    ReportFamily(
        name="fbrf_pmid_pmcid_doi",
        prefix="fbrf_pmid_pmcid_doi",
        remote_subdir="references",
        destination_dir=REFERENCES_DIR,
    ),
    ReportFamily(
        name="fb_synonym",
        prefix="fb_synonym",
        remote_subdir="synonyms",
        destination_dir=GENES_DIR,
    ),
]


USER_AGENT = "fl.AI-reagent-stocker FlyBase refresher"
REQUEST_HEADERS = {
    "User-Agent": USER_AGENT,
    "Accept": "*/*",
    "Accept-Language": "en-US,en;q=0.9",
}
MD5SUM_URL = urljoin(CURRENT_PRECOMPUTED_BASE, "md5sum.txt")
_CURRENT_FILE_LISTING: Optional[List[str]] = None


class FlyBaseRequestBlocked(RuntimeError):
    """Raised when FlyBase responds with a WAF challenge instead of content."""


def _check_for_waf_challenge(response, url: str) -> None:
    status = getattr(response, "status", response.getcode())
    waf_action = response.headers.get("x-amzn-waf-action", "")
    if status == 202 or waf_action.lower() == "challenge":
        raise FlyBaseRequestBlocked(
            f"FlyBase blocked urllib request for {url} with status {status}"
        )


def _curl_path() -> str:
    curl = shutil.which("curl")
    if not curl:
        raise RuntimeError(
            "FlyBase rejected urllib requests and no 'curl' executable was found "
            "for fallback downloads."
        )
    return curl


def _fetch_bytes_with_curl(url: str) -> bytes:
    result = subprocess.run(
        [_curl_path(), "--fail", "--silent", "--show-error", "--location", url],
        check=True,
        capture_output=True,
    )
    return result.stdout


def _download_with_curl(url: str, destination_path: Path) -> None:
    subprocess.run(
        [
            _curl_path(),
            "--fail",
            "--silent",
            "--show-error",
            "--location",
            "--output",
            str(destination_path),
            url,
        ],
        check=True,
    )


def fetch_text(url: str) -> str:
    req = Request(url, headers=REQUEST_HEADERS)
    try:
        with urlopen(req, timeout=60) as response:
            _check_for_waf_challenge(response, url)
            payload = response.read()
    except FlyBaseRequestBlocked:
        payload = _fetch_bytes_with_curl(url)
    return payload.decode("utf-8", errors="replace")


def basename_from_url(url: str) -> str:
    return Path(urlparse(url).path).name


def current_precomputed_listing() -> List[str]:
    global _CURRENT_FILE_LISTING
    if _CURRENT_FILE_LISTING is None:
        entries: List[str] = []
        for raw_line in fetch_text(MD5SUM_URL).splitlines():
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            entries.append(parts[1].strip())
        _CURRENT_FILE_LISTING = entries
    return _CURRENT_FILE_LISTING


def discover_latest_url(family: ReportFamily) -> str:
    matches = [
        relative_path
        for relative_path in current_precomputed_listing()
        if relative_path.startswith(f"{family.remote_subdir}/")
        and basename_from_url(relative_path).startswith(family.prefix)
        and basename_from_url(relative_path).endswith(family.preferred_extensions)
    ]
    if not matches:
        if family.required:
            raise FileNotFoundError(
                f"Could not discover a download URL for prefix '{family.prefix}' from {MD5SUM_URL}"
            )
        return ""

    def sort_key(relative_path: str) -> tuple[int, str]:
        name = basename_from_url(relative_path)
        ext_rank = next(
            (
                idx
                for idx, ext in enumerate(family.preferred_extensions)
                if name.endswith(ext)
            ),
            len(family.preferred_extensions),
        )
        return (ext_rank, name)

    matches.sort(key=sort_key)
    return urljoin(CURRENT_PRECOMPUTED_BASE, matches[0])


def atomic_download(url: str, destination_path: Path) -> int:
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination_path.with_name(destination_path.name + ".tmp")
    if temp_path.exists():
        temp_path.unlink()

    req = Request(url, headers=REQUEST_HEADERS)
    try:
        with urlopen(req, timeout=60) as response, temp_path.open("wb") as out_handle:
            _check_for_waf_challenge(response, url)
            shutil.copyfileobj(response, out_handle)
    except FlyBaseRequestBlocked:
        temp_path.unlink(missing_ok=True)
        _download_with_curl(url, temp_path)

    size = temp_path.stat().st_size
    if size <= 0:
        temp_path.unlink(missing_ok=True)
        raise ValueError(f"Downloaded file is empty: {destination_path.name}")

    temp_path.replace(destination_path)
    return size


def _list_family_candidates(family: ReportFamily) -> List[Path]:
    candidates = list(family.destination_dir.glob(f"{family.prefix}*.tsv")) + list(
        family.destination_dir.glob(f"{family.prefix}*.tsv.gz")
    )
    return sorted(set(candidates))


def prune_older_versions(family: ReportFamily, keep_path: Path) -> List[Path]:
    old_files: List[Path] = []
    for candidate in _list_family_candidates(family):
        if candidate == keep_path:
            continue
        candidate.unlink(missing_ok=True)
        old_files.append(candidate)
    return old_files


def list_older_versions(family: ReportFamily, keep_path: Path) -> List[Path]:
    return [path for path in _list_family_candidates(family) if path != keep_path]


def _manifest_prefixes_by_dir(manifest: Iterable[ReportFamily]) -> Dict[Path, Set[str]]:
    prefixes_by_dir: Dict[Path, Set[str]] = {}
    for family in manifest:
        prefixes_by_dir.setdefault(family.destination_dir, set()).add(family.prefix)
    return prefixes_by_dir


def _list_orphan_tsvs(manifest: Iterable[ReportFamily]) -> List[Path]:
    orphan_paths: List[Path] = []
    for directory, prefixes in _manifest_prefixes_by_dir(manifest).items():
        if not directory.exists():
            continue
        candidates = sorted(set(directory.glob("*.tsv")) | set(directory.glob("*.tsv.gz")))
        for candidate in candidates:
            if any(candidate.name.startswith(prefix) for prefix in prefixes):
                continue
            orphan_paths.append(candidate)
    return orphan_paths


def cleanup_orphan_tsvs(manifest: Iterable[ReportFamily]) -> List[Path]:
    orphan_paths = _list_orphan_tsvs(manifest)
    for path in orphan_paths:
        path.unlink(missing_ok=True)
    return orphan_paths


def find_latest_tsv_like_repo(directory: Path, prefix: str) -> Path:
    gz_files = sorted(directory.glob(f"{prefix}*.tsv.gz"), reverse=True)
    if gz_files:
        return gz_files[0]

    tsv_files = sorted(directory.glob(f"{prefix}*.tsv"), reverse=True)
    if tsv_files:
        return tsv_files[0]

    raise FileNotFoundError(f"No local files found for prefix '{prefix}' in {directory}")


def format_size(size: int) -> str:
    units = ["B", "KB", "MB", "GB"]
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{size} B"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Refresh FlyBase TSV report files.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Discover and report actions without downloading or deleting files.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    print("Refreshing FlyBase TSV report families")
    print(f"  Source base: {CURRENT_PRECOMPUTED_BASE}")
    print(f"  Repo root: {REPO_ROOT}")
    if args.dry_run:
        print("  Mode: dry-run (no downloads or deletions)")

    refreshed = 0
    skipped = 0
    removed_total = 0

    for family in MANIFEST:
        print(f"\n[{family.name}]")
        latest_url = discover_latest_url(family)
        latest_name = basename_from_url(latest_url)
        destination_path = family.destination_dir / latest_name

        print(f"  Latest URL: {latest_url}")
        print(f"  Destination: {destination_path}")

        status = "download"
        if destination_path.exists() and destination_path.stat().st_size > 0:
            status = "already current"
            size = destination_path.stat().st_size
        elif args.dry_run:
            size = 0
        else:
            size = atomic_download(latest_url, destination_path)
            refreshed += 1

        if status == "already current":
            skipped += 1
            print(f"  Status: already current ({format_size(size)})")
        elif args.dry_run:
            print("  Status: would download")
        else:
            print(f"  Status: downloaded ({format_size(size)})")

        removed: List[Path] = []
        if args.dry_run:
            removed = list_older_versions(family, destination_path)
        elif destination_path.exists() and destination_path.stat().st_size > 0:
            removed = prune_older_versions(family, destination_path)
            removed_total += len(removed)

        if removed:
            if args.dry_run:
                print("  Old versions that would be removed:")
            else:
                print("  Removed old versions:")
            for path in removed:
                print(f"    - {path.name}")
        else:
            print("  Removed old versions: none")

        try:
            selected = find_latest_tsv_like_repo(family.destination_dir, family.prefix)
        except FileNotFoundError:
            selected = None

        if args.dry_run and not destination_path.exists():
            if selected is None:
                print("  Active latest file now: none")
            else:
                print(f"  Active latest file now: {selected.name}")
            print(f"  Active latest file after download: {destination_path.name}")
            print("  Loader compatibility after download: OK")
        else:
            if selected is None:
                raise RuntimeError(
                    f"Compatibility check failed for {family.name}: no local file exists "
                    f"for prefix {family.prefix} after refresh"
                )
            print(f"  Active latest file: {selected.name}")
            if selected != destination_path:
                raise RuntimeError(
                    f"Compatibility check failed for {family.name}: expected {destination_path.name}, "
                    f"but repo-style latest selection chose {selected.name}"
                )
            print("  Loader compatibility: OK")

    orphan_paths = _list_orphan_tsvs(MANIFEST) if args.dry_run else cleanup_orphan_tsvs(MANIFEST)
    if orphan_paths:
        if args.dry_run:
            print("\nOrphan TSV/TSV.GZ files that would be removed:")
        else:
            print("\nRemoved orphan TSV/TSV.GZ files:")
        for path in orphan_paths:
            print(f"  - {path.relative_to(REPO_ROOT)}")
        if not args.dry_run:
            removed_total += len(orphan_paths)
    else:
        print("\nOrphan TSV/TSV.GZ files removed: none")

    print("\nSummary:")
    print(f"  Families processed: {len(MANIFEST)}")
    print(f"  Newly downloaded: {refreshed}")
    print(f"  Already current: {skipped}")
    print(f"  Older files removed: {removed_total if not args.dry_run else 'dry-run'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
