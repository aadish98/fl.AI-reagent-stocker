#!/usr/bin/env python3
"""
Download the latest FlyBase report families used by this repo.

The script discovers TSV/TSV.GZ files from the official FlyBase
current-release precomputed index under:
    https://s3ftp.flybase.org/releases/current/precomputed_files/

When ``--with-xml`` is passed, it also discovers chado XML/XML.GZ files from:
    https://s3ftp.flybase.org/releases/current/chado-xml/

It downloads only the report families this repository uses, saves them into the
expected local data directories, and removes older versions of the same prefix
only after the latest file is present and non-empty. After all manifest families
have been processed, it also removes orphan files from the managed FlyBase
directories so stale helper files do not accumulate. Orphan cleanup is scoped
per file kind (TSV vs XML), so omitting ``--with-xml`` never deletes existing
local XML files.

Usage:
    python scripts/refresh_flybase_data.py
    python scripts/refresh_flybase_data.py --with-xml
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Set, Tuple
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen


REPO_ROOT = Path(__file__).resolve().parents[1]
FLYBASE_ROOT = REPO_ROOT / "data" / "flybase"
ALLELES_AND_STOCKS_DIR = FLYBASE_ROOT / "alleles_and_stocks"
REFERENCES_DIR = FLYBASE_ROOT / "references"
GENES_DIR = FLYBASE_ROOT / "genes"
TRANSGENIC_CONSTRUCTS_DIR = FLYBASE_ROOT / "transgenic_constructs"
TRANSGENIC_INSERTIONS_DIR = FLYBASE_ROOT / "transgenic_insertions"

CURRENT_PRECOMPUTED_BASE = "https://s3ftp.flybase.org/releases/current/precomputed_files/"
CURRENT_CHADO_XML_BASE = "https://s3ftp.flybase.org/releases/current/chado-xml/"


def _parse_md5sum_listing(payload: str) -> List[str]:
    """Parse a FlyBase ``md5sum.txt`` payload into relative paths."""
    entries: List[str] = []
    for raw_line in payload.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        entries.append(parts[1].strip())
    return entries


_CHADO_XML_HREF_RE = re.compile(
    r"href=['\"][^'\"]*?/chado-xml/([^'\"/]+\.xml(?:\.gz)?)['\"]",
    re.IGNORECASE,
)


def _parse_chado_xml_listing(payload: str) -> List[str]:
    """Parse the FlyBase ``chado-xml/`` HTML index into relative filenames."""
    seen: Dict[str, None] = {}
    for match in _CHADO_XML_HREF_RE.finditer(payload):
        seen.setdefault(match.group(1), None)
    return list(seen.keys())


@dataclass(frozen=True)
class RemoteSource:
    """Describes a remote FlyBase listing/download source.

    ``base_url`` is joined with the relative paths returned by ``listing_parser``
    when constructing download URLs. ``listing_url`` is the listing endpoint
    fetched once per run and parsed by ``listing_parser``.
    """

    name: str
    base_url: str
    listing_url: str
    listing_parser: Callable[[str], List[str]]


PRECOMPUTED_SOURCE = RemoteSource(
    name="precomputed",
    base_url=CURRENT_PRECOMPUTED_BASE,
    listing_url=urljoin(CURRENT_PRECOMPUTED_BASE, "md5sum.txt"),
    listing_parser=_parse_md5sum_listing,
)

CHADO_XML_SOURCE = RemoteSource(
    name="chado_xml",
    base_url=CURRENT_CHADO_XML_BASE,
    listing_url=urljoin(CURRENT_CHADO_XML_BASE, "index.html"),
    listing_parser=_parse_chado_xml_listing,
)


@dataclass(frozen=True)
class ReportFamily:
    name: str
    prefix: str
    remote_subdir: str
    destination_dir: Path
    source: RemoteSource = PRECOMPUTED_SOURCE
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
    ReportFamily(
        name="fbgn_annotation_ID",
        prefix="fbgn_annotation_ID",
        remote_subdir="genes",
        destination_dir=GENES_DIR,
    ),
]


# Optional XML families. Files in chado-xml/ live directly under the source
# base URL, so ``remote_subdir`` is empty. Only the families used by this repo
# are listed here so the rest of the chado-xml/ directory is never touched.
XML_MANIFEST: List[ReportFamily] = [
    ReportFamily(
        name="chado_FBst",
        prefix="chado_FBst",
        remote_subdir="",
        destination_dir=ALLELES_AND_STOCKS_DIR,
        source=CHADO_XML_SOURCE,
        preferred_extensions=(".xml.gz", ".xml"),
    ),
    ReportFamily(
        name="chado_FBti",
        prefix="chado_FBti",
        remote_subdir="",
        destination_dir=TRANSGENIC_INSERTIONS_DIR,
        source=CHADO_XML_SOURCE,
        preferred_extensions=(".xml.gz", ".xml"),
    ),
]


USER_AGENT = "fl.AI-reagent-stocker FlyBase refresher"
REQUEST_HEADERS = {
    "User-Agent": USER_AGENT,
    "Accept": "*/*",
    "Accept-Language": "en-US,en;q=0.9",
}
_LISTING_CACHE: Dict[str, List[str]] = {}


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


def current_listing(source: RemoteSource) -> List[str]:
    """Return cached listing entries for ``source`` (relative paths)."""
    cached = _LISTING_CACHE.get(source.name)
    if cached is None:
        cached = source.listing_parser(fetch_text(source.listing_url))
        _LISTING_CACHE[source.name] = cached
    return cached


def discover_latest_url(family: ReportFamily) -> str:
    listing = current_listing(family.source)
    subdir_prefix = f"{family.remote_subdir}/" if family.remote_subdir else ""

    matches: List[str] = []
    for relative_path in listing:
        if subdir_prefix and not relative_path.startswith(subdir_prefix):
            continue
        if not subdir_prefix and "/" in relative_path:
            continue
        name = basename_from_url(relative_path)
        if not name.startswith(family.prefix):
            continue
        if not name.endswith(family.preferred_extensions):
            continue
        matches.append(relative_path)

    if not matches:
        if family.required:
            raise FileNotFoundError(
                f"Could not discover a download URL for prefix '{family.prefix}' "
                f"from {family.source.listing_url}"
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
    return urljoin(family.source.base_url, matches[0])


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
    candidates: Set[Path] = set()
    for extension in family.preferred_extensions:
        candidates.update(family.destination_dir.glob(f"{family.prefix}*{extension}"))
    return sorted(candidates)


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


@dataclass(frozen=True)
class _DirExtensionScope:
    directory: Path
    extensions: Tuple[str, ...]


def _scope_prefixes(
    manifest: Iterable[ReportFamily],
) -> Dict[_DirExtensionScope, Set[str]]:
    """Group manifest prefixes by destination directory + file-extension family.

    Orphan cleanup must respect the file kind a manifest manages: a TSV-only
    refresh should never remove sibling XML files (and vice versa) that live
    in the same destination directory.
    """
    scopes: Dict[_DirExtensionScope, Set[str]] = {}
    for family in manifest:
        scope = _DirExtensionScope(
            directory=family.destination_dir,
            extensions=tuple(family.preferred_extensions),
        )
        scopes.setdefault(scope, set()).add(family.prefix)
    return scopes


def _list_orphan_files(manifest: Iterable[ReportFamily]) -> List[Path]:
    orphan_paths: List[Path] = []
    for scope, prefixes in _scope_prefixes(manifest).items():
        if not scope.directory.exists():
            continue
        candidates: Set[Path] = set()
        for extension in scope.extensions:
            candidates.update(scope.directory.glob(f"*{extension}"))
        for candidate in sorted(candidates):
            if any(candidate.name.startswith(prefix) for prefix in prefixes):
                continue
            orphan_paths.append(candidate)
    return orphan_paths


def cleanup_orphan_files(manifest: Iterable[ReportFamily]) -> List[Path]:
    orphan_paths = _list_orphan_files(manifest)
    for path in orphan_paths:
        path.unlink(missing_ok=True)
    return orphan_paths


def find_latest_local_for_family(family: ReportFamily) -> Path:
    """Return the latest local file for ``family`` honoring its preferred extensions."""
    for extension in family.preferred_extensions:
        matches = sorted(
            family.destination_dir.glob(f"{family.prefix}*{extension}"),
            reverse=True,
        )
        if matches:
            return matches[0]

    raise FileNotFoundError(
        f"No local files found for prefix '{family.prefix}' in {family.destination_dir}"
    )


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
    parser = argparse.ArgumentParser(description="Refresh FlyBase report files.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Discover and report actions without downloading or deleting files.",
    )
    parser.add_argument(
        "--with-xml",
        action="store_true",
        help=(
            "Also refresh chado XML report families used by the helper scripts. "
            "These files (e.g. chado_FBst, chado_FBti) are large and live under "
            "https://s3ftp.flybase.org/releases/current/chado-xml/."
        ),
    )
    return parser.parse_args()


def _process_family(family: ReportFamily, dry_run: bool) -> Tuple[bool, bool, int]:
    """Refresh a single family and return (refreshed, skipped, removed_count)."""
    print(f"\n[{family.name}]")
    latest_url = discover_latest_url(family)
    latest_name = basename_from_url(latest_url)
    destination_path = family.destination_dir / latest_name

    print(f"  Source: {family.source.name} ({family.source.base_url})")
    print(f"  Latest URL: {latest_url}")
    print(f"  Destination: {destination_path}")

    refreshed = False
    skipped = False
    if destination_path.exists() and destination_path.stat().st_size > 0:
        skipped = True
        size = destination_path.stat().st_size
        print(f"  Status: already current ({format_size(size)})")
    elif dry_run:
        print("  Status: would download")
    else:
        size = atomic_download(latest_url, destination_path)
        refreshed = True
        print(f"  Status: downloaded ({format_size(size)})")

    removed: List[Path] = []
    if dry_run:
        removed = list_older_versions(family, destination_path)
    elif destination_path.exists() and destination_path.stat().st_size > 0:
        removed = prune_older_versions(family, destination_path)

    if removed:
        if dry_run:
            print("  Old versions that would be removed:")
        else:
            print("  Removed old versions:")
        for path in removed:
            print(f"    - {path.name}")
    else:
        print("  Removed old versions: none")

    try:
        selected = find_latest_local_for_family(family)
    except FileNotFoundError:
        selected = None

    if dry_run and not destination_path.exists():
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

    return refreshed, skipped, len(removed) if not dry_run else 0


def main() -> int:
    args = parse_args()

    manifest: List[ReportFamily] = list(MANIFEST)
    if args.with_xml:
        manifest.extend(XML_MANIFEST)

    print("Refreshing FlyBase report families")
    print(f"  Precomputed source: {CURRENT_PRECOMPUTED_BASE}")
    if args.with_xml:
        print(f"  Chado XML source:   {CURRENT_CHADO_XML_BASE}")
    else:
        print("  Chado XML source:   skipped (use --with-xml to include)")
    print(f"  Repo root: {REPO_ROOT}")
    if args.dry_run:
        print("  Mode: dry-run (no downloads or deletions)")

    refreshed = 0
    skipped = 0
    removed_total = 0

    for family in manifest:
        family_refreshed, family_skipped, family_removed = _process_family(
            family, dry_run=args.dry_run
        )
        if family_refreshed:
            refreshed += 1
        if family_skipped:
            skipped += 1
        removed_total += family_removed

    orphan_paths = (
        _list_orphan_files(manifest) if args.dry_run else cleanup_orphan_files(manifest)
    )
    if orphan_paths:
        if args.dry_run:
            print("\nOrphan files that would be removed:")
        else:
            print("\nRemoved orphan files:")
        for path in orphan_paths:
            print(f"  - {path.relative_to(REPO_ROOT)}")
        if not args.dry_run:
            removed_total += len(orphan_paths)
    else:
        print("\nOrphan files removed: none")

    print("\nSummary:")
    print(f"  Families processed: {len(manifest)}")
    print(f"  Newly downloaded: {refreshed}")
    print(f"  Already current: {skipped}")
    print(f"  Older files removed: {removed_total if not args.dry_run else 'dry-run'}")
    if not args.with_xml:
        print(
            "  Chado XML refresh: skipped "
            "(re-run with --with-xml to refresh chado_FBst / chado_FBti)"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
