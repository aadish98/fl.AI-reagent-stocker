"""
PubMed/Entrez client for fetching metadata and managing cache.

This module provides:
    - PubMedCache: Persistent cache for PubMed metadata (title, abstract)
    - PubMedClient: Client for fetching metadata from PubMed via Entrez
"""

import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import requests

from ..utils import clean_id

# Try to import Biopython
try:
    from Bio import Entrez
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


class PubMedCache:
    """
    Persistent cache for PubMed metadata.

    Cache schema is compatible with the richer fl.AI-CLI metadata payload while
    remaining backward-compatible with legacy title/abstract-only rows.
    """

    CACHE_COLUMNS = [
        "pmid",
        "pmcid",
        "title",
        "abstract",
        "year",
        "journal",
        "authors",
        "doi",
        "source",
        "updated_at",
    ]
    COLUMN_ALIASES = {
        "pmid": ["PMID", "Pmid"],
        "pmcid": ["PMCID", "Pmcid"],
        "title": ["Title"],
        "abstract": ["Abstract"],
        "year": ["Year", "publication_year", "PublicationYear", "date_published"],
        "journal": ["Journal"],
        "authors": ["Authors", "all_authors", "AllAuthors"],
        "doi": ["DOI", "Doi"],
        "source": ["Source"],
        "updated_at": ["UpdatedAt", "updated", "Updated"],
    }

    def __init__(self, cache_path: Path):
        """
        Initialize the cache.
        
        Args:
            cache_path: Path to the CSV cache file
        """
        self.cache_path = Path(cache_path)
        self._cache: Dict[str, dict] = {}
        self._loaded = False

    @staticmethod
    def _normalize_authors(raw_authors: Any) -> List[str]:
        """Normalize authors from list or semicolon-delimited string."""
        if raw_authors is None:
            return []
        if isinstance(raw_authors, str):
            return [a.strip() for a in raw_authors.split(";") if a and a.strip()]
        try:
            return [str(a).strip() for a in list(raw_authors) if str(a).strip()]
        except Exception:
            return []

    def _normalize_entry(self, metadata: Optional[dict]) -> dict:
        """Normalize metadata to canonical cache schema."""
        metadata = metadata or {}

        def _get_value(field: str) -> Any:
            keys = [field] + self.COLUMN_ALIASES.get(field, [])
            for key in keys:
                if key in metadata:
                    value = metadata.get(key)
                    if value is not None and str(value).strip() != "":
                        return value
            for key in keys:
                if key in metadata:
                    return metadata.get(key)
            return ""

        return {
            "title": str(_get_value("title") or ""),
            "abstract": str(_get_value("abstract") or ""),
            "year": str(_get_value("year") or ""),
            "journal": str(_get_value("journal") or ""),
            "authors": self._normalize_authors(_get_value("authors")),
            "doi": str(_get_value("doi") or ""),
            "pmcid": str(_get_value("pmcid") or ""),
            "source": str(_get_value("source") or ""),
            "updated_at": str(_get_value("updated_at") or ""),
        }

    def _merged_entry(self, pmid: str, metadata: Optional[dict]) -> dict:
        """Merge metadata with existing cached values, preferring non-empty new values."""
        existing = self._cache.get(pmid, {})
        normalized = self._normalize_entry(metadata)
        return {
            "title": normalized["title"] or str(existing.get("title", "") or ""),
            "abstract": normalized["abstract"] or str(existing.get("abstract", "") or ""),
            "year": normalized["year"] or str(existing.get("year", "") or ""),
            "journal": normalized["journal"] or str(existing.get("journal", "") or ""),
            "authors": normalized["authors"] or self._normalize_authors(existing.get("authors", [])),
            "doi": normalized["doi"] or str(existing.get("doi", "") or ""),
            "pmcid": normalized["pmcid"] or str(existing.get("pmcid", "") or ""),
            "source": normalized["source"] or str(existing.get("source", "") or ""),
            "updated_at": normalized["updated_at"] or str(existing.get("updated_at", "") or ""),
        }

    def _save_all(self):
        """Persist the full in-memory cache to disk in canonical schema."""
        rows = []
        for pmid, meta in sorted(self._cache.items()):
            if not str(pmid).isdigit():
                continue
            rows.append(
                {
                    "pmid": str(pmid),
                    "pmcid": str(meta.get("pmcid", "") or ""),
                    "title": str(meta.get("title", "") or ""),
                    "abstract": str(meta.get("abstract", "") or ""),
                    "year": str(meta.get("year", "") or ""),
                    "journal": str(meta.get("journal", "") or ""),
                    "authors": "; ".join(self._normalize_authors(meta.get("authors", []))),
                    "doi": str(meta.get("doi", "") or ""),
                    "source": str(meta.get("source", "") or ""),
                    "updated_at": str(meta.get("updated_at", "") or ""),
                }
            )
        out_df = pd.DataFrame(rows, columns=self.CACHE_COLUMNS)
        try:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            out_df.to_csv(self.cache_path, index=False)
        except OSError as e:
            print(f"    Warning: Could not persist PubMed cache to {self.cache_path}: {e}")
    
    def load(self) -> Dict[str, dict]:
        """
        Load the cache from disk.
        
        Returns:
            Dict mapping PMID to canonical metadata dict.
        """
        if self._loaded:
            return self._cache

        if self.cache_path.exists():
            try:
                cache_df = pd.read_csv(self.cache_path, dtype=str, keep_default_na=False)
                drop_cols = [c for c in cache_df.columns if str(c).startswith("Unnamed")]
                if drop_cols:
                    cache_df = cache_df.drop(columns=drop_cols)
                normalized_columns = {str(col).strip().lower(): col for col in cache_df.columns}
                rename_map = {}
                for canonical in self.CACHE_COLUMNS:
                    if canonical in cache_df.columns:
                        continue
                    aliases = [canonical] + self.COLUMN_ALIASES.get(canonical, [])
                    for alias in aliases:
                        matched_col = normalized_columns.get(str(alias).strip().lower())
                        if matched_col is not None:
                            rename_map[matched_col] = canonical
                            break
                if rename_map:
                    cache_df = cache_df.rename(columns=rename_map)
                if "pmid" not in cache_df.columns:
                    cache_df = pd.DataFrame(columns=self.CACHE_COLUMNS)
                for col in self.CACHE_COLUMNS:
                    if col not in cache_df.columns:
                        cache_df[col] = ""
                cache_df["pmid"] = cache_df["pmid"].apply(clean_id)
                cache_df = cache_df[cache_df["pmid"].str.isdigit()].copy()
                if len(cache_df) > 0:
                    for _, row in cache_df.iterrows():
                        pmid = str(row.get("pmid", "")).strip()
                        if not pmid:
                            continue
                        self._cache[pmid] = self._normalize_entry(row.to_dict())
                print(f"    Loaded {len(self._cache)} entries from PubMed cache")
            except Exception as e:
                print(f"    Warning: Could not load PubMed cache: {e}")
        else:
            print(f"    PubMed cache file not found, will create: {self.cache_path}")
        
        self._loaded = True
        return self._cache
    
    def get(self, pmid: str) -> Optional[dict]:
        """
        Get cached metadata for a PMID.
        
        Args:
            pmid: The PMID to look up
        
        Returns:
            Metadata dict, or None if not cached
        """
        self.load()
        pmid_clean = clean_id(pmid)
        return self._cache.get(pmid_clean)
    
    def set(self, pmid: str, metadata: dict):
        """
        Add a single entry to the cache (in memory only).
        
        Args:
            pmid: The PMID
            metadata: Dict with metadata fields
        """
        self.load()
        pmid_clean = clean_id(pmid)
        if pmid_clean:
            self._cache[pmid_clean] = self._merged_entry(pmid_clean, metadata)
    
    def save_entries(self, new_entries: Dict[str, dict]):
        """
        Save new/updated entries into the cache file on disk.
        
        Args:
            new_entries: Dict mapping PMID to metadata dict
        """
        if not new_entries:
            return

        self.load()
        for pmid, metadata in new_entries.items():
            pmid_clean = clean_id(pmid)
            if not pmid_clean:
                continue
            self._cache[pmid_clean] = self._merged_entry(pmid_clean, metadata)
        self._save_all()
        print(f"    Saved {len(new_entries)} PubMed metadata updates")
    
    def __contains__(self, pmid: str) -> bool:
        """Check if a PMID is in the cache."""
        self.load()
        return clean_id(pmid) in self._cache
    
    def __len__(self) -> int:
        """Return the number of cached entries."""
        self.load()
        return len(self._cache)


class PubMedClient:
    """
    Client for fetching metadata from PubMed via NCBI Entrez.
    
    Uses Biopython's Entrez module for API calls.
    Integrates with PubMedCache for persistent caching.
    """
    
    def __init__(
        self,
        cache: Optional[PubMedCache] = None,
        email: str = "flybase_scraper@example.com",
        api_key: Optional[str] = None
    ):
        """
        Initialize the PubMed client.
        
        Args:
            cache: Optional PubMedCache for persistent caching
            email: Email for NCBI Entrez (required by NCBI)
            api_key: Optional NCBI API key for higher rate limits
        """
        self.cache = cache
        self.email = email
        self.api_key = api_key
        
        if not BIOPYTHON_AVAILABLE:
            print("    Warning: Biopython not installed. PubMed fetching will be limited.")
            print("    Install with: pip install biopython")
    
    def fetch_metadata(self, pmids: List[str]) -> Dict[str, dict]:
        """
        Batch fetch metadata for PMIDs from PubMed.

        Uses cache-first retrieval and only queries PubMed for missing metadata.
        By default, partial cached metadata is accepted. Richer author/date/journal
        backfilling is only attempted when explicitly requested.
        
        Args:
            pmids: List of PMIDs to fetch
        
        Returns:
            Dict mapping PMID to metadata dict:
            {
              'title','abstract','year','journal','authors','doi','pmcid',
              'source','updated_at'
            }
        """
        return self._fetch_metadata(pmids, enrich_incomplete=False)

    def _fetch_metadata(self, pmids: List[str], *, enrich_incomplete: bool) -> Dict[str, dict]:
        """Internal metadata fetch with optional author/date enrichment."""
        results: Dict[str, dict] = {}
        pmids_to_fetch = []

        def _has_keyword_text(meta: dict) -> bool:
            if not isinstance(meta, dict):
                return False
            return bool(
                str(meta.get("title", "")).strip()
                or str(meta.get("abstract", "")).strip()
            )

        def _needs_enrichment(meta: dict) -> bool:
            already_enriched = bool(
                str(meta.get("source", "")).strip()
                or str(meta.get("updated_at", "")).strip()
            )
            if not _has_keyword_text(meta):
                return not already_enriched
            if not enrich_incomplete:
                return False
            has_authors = bool(meta.get("authors", []))
            has_year = bool(str(meta.get("year", "")).strip())
            has_journal = bool(str(meta.get("journal", "")).strip())
            if has_authors and has_year and has_journal:
                return False
            return not already_enriched

        for pmid in pmids:
            pmid_str = clean_id(pmid)
            if not pmid_str or not pmid_str.isdigit():
                continue

            if self.cache is not None and pmid_str in self.cache:
                cached = self.cache.get(pmid_str)
                if cached:
                    results[pmid_str] = cached
                    if _needs_enrichment(cached):
                        pmids_to_fetch.append(pmid_str)
                else:
                    pmids_to_fetch.append(pmid_str)
            else:
                pmids_to_fetch.append(pmid_str)

        pmids_to_fetch = sorted(set(pmids_to_fetch))
        cache_hits = len(results) - len([p for p in pmids_to_fetch if p in results])
        if pmids_to_fetch or cache_hits:
            print(f"    PubMed cache: {cache_hits} hits, {len(pmids_to_fetch)} to fetch")
        if not pmids_to_fetch:
            return results

        if not BIOPYTHON_AVAILABLE:
            print("    Warning: Biopython not installed. Cannot fetch PubMed metadata.")
            return results

        print(f"    Fetching metadata for {len(pmids_to_fetch)} PMIDs from PubMed...")

        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key

        def _normalize_authors(raw_authors: Any) -> List[str]:
            if raw_authors is None:
                return []
            if isinstance(raw_authors, str):
                return [a.strip() for a in raw_authors.split(";") if a and a.strip()]
            try:
                return [str(a).strip() for a in list(raw_authors) if str(a).strip()]
            except Exception:
                return []

        def _extract_article_ids(pubmed_data: dict) -> Dict[str, str]:
            out = {"doi": "", "pmcid": ""}
            id_list = pubmed_data.get("ArticleIdList", []) or []
            for aid in id_list:
                try:
                    id_type = str(aid.attributes.get("IdType", "")).lower()
                except Exception:
                    id_type = ""
                value = str(aid).strip()
                if not value:
                    continue
                if id_type == "doi" and not out["doi"]:
                    out["doi"] = value
                elif id_type == "pmc" and not out["pmcid"]:
                    out["pmcid"] = value if value.upper().startswith("PMC") else f"PMC{value}"
            return out

        def _extract_authors(article_data: dict) -> List[str]:
            author_list = article_data.get("AuthorList", []) or []
            authors: List[str] = []
            for author_obj in author_list:
                if not isinstance(author_obj, dict):
                    continue
                ln = str(author_obj.get("LastName", "")).strip()
                fn = str(author_obj.get("ForeName", "")).strip()
                if ln:
                    authors.append(f"{ln}, {fn}" if fn else ln)
            return _normalize_authors(authors)

        def _extract_year(article_data: dict) -> str:
            try:
                pub_date = article_data.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                year = str(pub_date.get("Year", "")).strip()
                if year:
                    return year
                medline_date = str(pub_date.get("MedlineDate", "")).strip()
                digits = "".join(ch for ch in medline_date if ch.isdigit())
                if len(digits) >= 4:
                    return digits[:4]
            except Exception:
                pass
            return ""

        batch_size = 200  # NCBI limit
        total_fetched = 0

        for i in range(0, len(pmids_to_fetch), batch_size):
            batch = pmids_to_fetch[i:i + batch_size]
            batch_entries: Dict[str, dict] = {}

            try:
                handle = Entrez.efetch(
                    db="pubmed",
                    id=",".join(batch),
                    rettype="xml",
                    retmode="xml"
                )
                records = Entrez.read(handle)
                handle.close()

                for article in records.get('PubmedArticle', []):
                    try:
                        medline = article['MedlineCitation']
                        pmid = str(medline['PMID'])

                        article_data = medline.get('Article', {})
                        title = str(article_data.get('ArticleTitle', '') or '')
                        abstract_data = article_data.get('Abstract', {})
                        abstract_parts = abstract_data.get('AbstractText', [])
                        if isinstance(abstract_parts, list):
                            abstract = ' '.join(str(p) for p in abstract_parts)
                        else:
                            abstract = str(abstract_parts or '')
                        journal_info = article_data.get('Journal', {})
                        journal = str(
                            journal_info.get('Title', '')
                            or journal_info.get('ISOAbbreviation', '')
                            or ''
                        )
                        ids = _extract_article_ids(article.get('PubmedData', {}))
                        metadata = {
                            'title': title,
                            'abstract': abstract,
                            'year': _extract_year(article_data),
                            'journal': journal,
                            'authors': _extract_authors(article_data),
                            'doi': ids.get("doi", ""),
                            'pmcid': ids.get("pmcid", ""),
                            'source': "entrez",
                            'updated_at': datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        }
                        results[pmid] = metadata
                        batch_entries[pmid] = metadata

                        if self.cache is not None:
                            self.cache.set(pmid, metadata)

                    except (KeyError, TypeError):
                        continue

                # Mark PMIDs that returned no result as already-attempted so
                # they are not re-fetched on subsequent runs.
                now_ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
                for pmid_str in batch:
                    if pmid_str not in batch_entries and pmid_str not in results:
                        empty_marker = {
                            'title': '', 'abstract': '', 'year': '',
                            'journal': '', 'authors': [], 'doi': '',
                            'pmcid': '', 'source': 'entrez',
                            'updated_at': now_ts,
                        }
                        results[pmid_str] = empty_marker
                        batch_entries[pmid_str] = empty_marker
                        if self.cache is not None:
                            self.cache.set(pmid_str, empty_marker)

                total_fetched += len(batch_entries)

                if self.cache is not None and batch_entries:
                    self.cache.save_entries(batch_entries)

                time.sleep(0.4)

            except Exception as e:
                print(f"    Warning: PubMed fetch failed for batch: {e}")
                continue

        if self.cache is not None:
            for pmid in pmids:
                pmid_str = clean_id(pmid)
                if not pmid_str:
                    continue
                if pmid_str not in results:
                    cached = self.cache.get(pmid_str)
                    if cached:
                        results[pmid_str] = cached

        return results

    def fetch_author_date_metadata(self, pmids: List[str]) -> Dict[str, Dict]:
        """
        Return author/date/journal metadata for PMIDs.

        This method reuses the unified PubMed metadata cache and only relies on
        Entrez via :meth:`fetch_metadata` for missing metadata.

        Args:
            pmids: List of PMIDs to fetch

        Returns:
            Dict mapping PMID to {
                'first_author', 'corresponding_author', 'all_authors',
                'date_published', 'title', 'journal'
            }
        """
        results: Dict[str, Dict[str, str]] = {}
        meta_by_pmid = self._fetch_metadata(pmids, enrich_incomplete=True)
        for pmid in pmids:
            pmid_str = clean_id(pmid)
            if not pmid_str:
                continue
            meta = meta_by_pmid.get(pmid_str, {})
            authors = meta.get("authors", [])
            if isinstance(authors, str):
                authors = [a.strip() for a in authors.split(";") if a.strip()]
            first_author = str(authors[0]).strip() if authors else ""
            results[pmid_str] = {
                "first_author": first_author,
                "corresponding_author": "",
                "all_authors": "; ".join([str(a).strip() for a in authors if str(a).strip()]),
                "date_published": str(meta.get("year", "") or ""),
                "title": str(meta.get("title", "") or ""),
                "journal": str(meta.get("journal", "") or ""),
            }
        return results
    
    def check_keywords(
        self,
        pmid: str,
        keywords: List[str],
        metadata: Optional[Dict[str, dict]] = None
    ) -> bool:
        """
        Check if a PMID's title/abstract contains any of the keywords.
        
        Args:
            pmid: The PMID to check
            keywords: List of lowercase keywords to search for
            metadata: Optional pre-fetched metadata dict
        
        Returns:
            True if any keyword is found in title or abstract
        """
        pmid_str = clean_id(pmid)
        if not pmid_str:
            return False
        
        # Get metadata from provided dict, cache, or fetch
        meta = None
        if metadata and pmid_str in metadata:
            meta = metadata[pmid_str]
        elif self.cache is not None and pmid_str in self.cache:
            meta = self.cache.get(pmid_str)
        
        if not meta:
            return False
        
        text = (meta.get('title', '') + ' ' + meta.get('abstract', '')).lower()
        return any(kw.lower() in text for kw in keywords)


def convert_pmid_pmcid(ids: List[str], id_type: str = "pmid") -> Dict[str, str]:
    """
    Convert between PMID and PMCID using NCBI's ID Converter API.
    
    Args:
        ids: List of IDs to convert
        id_type: "pmid" to convert PMIDs to PMCIDs, or "pmcid" to convert PMCIDs to PMIDs
    
    Returns:
        Dict mapping input ID to converted ID (only successful conversions)
    """
    if not ids:
        return {}
    
    base_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
    results = {}
    batch_size = 200
    
    clean_ids = [clean_id(id_val) for id_val in ids if clean_id(id_val)]
    
    if not clean_ids:
        return {}
    
    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i + batch_size]
        ids_str = ",".join(batch)
        
        try:
            params = {
                "ids": ids_str,
                "format": "json",
                "tool": "flybase_stock_scraper",
                "email": "flybase_scraper@example.com"
            }
            
            response = requests.get(base_url, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                records = data.get("records", [])
                
                for record in records:
                    pmid = record.get("pmid", "")
                    pmcid = record.get("pmcid", "")
                    
                    if id_type == "pmid" and pmid and pmcid:
                        results[str(pmid)] = pmcid
                    elif id_type == "pmcid" and pmid and pmcid:
                        pmcid_key = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
                        results[pmcid_key] = str(pmid)
            
            time.sleep(0.5)
            
        except Exception as e:
            print(f"    Warning: NCBI ID conversion failed for batch: {e}")
            continue
    
    return results
