"""
Full-text retrieval and functional validation with OpenAI.

This module provides:
    - FullTextFetcher: Fetches full paper text from PMC, Unpaywall, etc.
    - FunctionalValidator: Uses OpenAI to validate stock function from paper text
"""

import io
import os
import random
import re
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote

import pandas as pd
import requests

from ..config import ValidationStatus, Settings
from ..utils import clean_id

# Try to import PyPDF2 for PDF extraction
try:
    from PyPDF2 import PdfReader
    PYPDF2_AVAILABLE = True
except ImportError:
    PYPDF2_AVAILABLE = False

# Try to import OpenAI
try:
    from openai import OpenAI
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False

# Try to import Pydantic for structured responses
try:
    from pydantic import BaseModel
    PYDANTIC_AVAILABLE = True
except ImportError:
    PYDANTIC_AVAILABLE = False

# Try to import metapub
try:
    from metapub import PubMedFetcher
    METAPUB_AVAILABLE = True
except ImportError:
    METAPUB_AVAILABLE = False

# Try to import Biopython
try:
    from Bio import Entrez
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


###############################################################################
# Full Text Method Cache
###############################################################################

class FullTextMethodCache:
    """
    Persistent cache for PMID -> full text retrieval method mappings.
    
    Stores successful retrieval methods in a CSV file so subsequent runs
    can skip the cascade and use the known-good method directly.
    """
    
    def __init__(self, cache_path: Optional[Path] = None):
        """
        Initialize the cache.
        
        Args:
            cache_path: Path to the CSV cache file. If None, cache is disabled.
        """
        self.cache_path = Path(cache_path) if cache_path else None
        self._cache: Dict[str, str] = {}  # PMID -> method label (or "" if none)
        self._loaded = False
        self._pending_entries: Dict[str, str] = {}  # New entries to save
    
    def load(self) -> Dict[str, str]:
        """
        Load the cache from disk.
        
        Returns:
            Dict mapping PMID to successful retrieval method
        """
        if self._loaded:
            return self._cache
        
        if self.cache_path and self.cache_path.exists():
            try:
                cache_df = pd.read_csv(self.cache_path, dtype={'pmid': str, 'method': str})
                cache_df['pmid'] = cache_df['pmid'].apply(lambda x: str(x).strip() if pd.notna(x) else '')
                cache_df = cache_df[cache_df['pmid'] != ''].copy()
                
                if len(cache_df) > 0:
                    self._cache = {
                        row['pmid']: str(row['method']) if pd.notna(row['method']) else ''
                        for _, row in cache_df.iterrows()
                    }
                print(f"    Loaded {len(self._cache)} entries from full-text method cache")
            except Exception as e:
                print(f"    Warning: Could not load full-text method cache: {e}")
        elif self.cache_path:
            print(f"    Full-text method cache file not found, will create: {self.cache_path}")
        
        self._loaded = True
        return self._cache
    
    def get(self, pmid: str) -> Optional[str]:
        """
        Get cached method for a PMID.
        
        Args:
            pmid: The PMID to look up
        
        Returns:
            Method label string, or None if not cached
        """
        self.load()
        pmid_clean = str(pmid).strip()
        return self._cache.get(pmid_clean)
    
    def set(self, pmid: str, method: str):
        """
        Add a single entry to the cache (in memory and pending save).
        
        Args:
            pmid: The PMID
            method: The successful retrieval method label
        """
        self.load()
        pmid_clean = str(pmid).strip()
        if pmid_clean and method:
            self._cache[pmid_clean] = method
            self._pending_entries[pmid_clean] = method
    
    def save_pending(self):
        """
        Append pending entries to the cache file on disk.
        """
        if not self._pending_entries or not self.cache_path:
            return
        
        # Prepare DataFrame for new entries
        rows = [{'pmid': pmid, 'method': method} for pmid, method in self._pending_entries.items()]
        new_df = pd.DataFrame(rows)
        
        # Ensure parent directory exists
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Append to existing file or create new one
        if self.cache_path.exists():
            new_df.to_csv(self.cache_path, mode='a', header=False, index=False)
        else:
            new_df.to_csv(self.cache_path, index=False)
        
        print(f"    Appended {len(self._pending_entries)} new entries to full-text method cache")
        self._pending_entries.clear()
    
    def __contains__(self, pmid: str) -> bool:
        """Check if a PMID is in the cache."""
        self.load()
        return str(pmid).strip() in self._cache
    
    def __len__(self) -> int:
        """Return the number of cached entries."""
        self.load()
        return len(self._cache)


###############################################################################
# Full Text Fetcher
###############################################################################

class FullTextFetcher:
    """
    Fetches full paper text from multiple sources.
    
    Strategies (in priority order):
    1. PMC Open Access PDF
    2. PMC Open Access XML
    3. Europe PMC API
    4. PMC HTML scraping
    5. Unpaywall (if DOI available)
    
    Uses a persistent method cache to remember which strategy worked for each PMID,
    avoiding repeated cascade attempts across runs.
    """
    
    def __init__(
        self,
        unpaywall_token: Optional[str] = None,
        ncbi_api_key: Optional[str] = None,
        method_cache_path: Optional[Path] = None
    ):
        """
        Initialize the full-text fetcher.
        
        Args:
            unpaywall_token: Email/token for Unpaywall API
            ncbi_api_key: Optional NCBI API key
            method_cache_path: Path to persistent cache CSV for PMID -> method mapping
        """
        self.unpaywall_token = unpaywall_token or "aadish98@gmail.com"
        self.ncbi_api_key = ncbi_api_key
        self._method_cache = FullTextMethodCache(method_cache_path)
        self._session_cache: Dict[str, str] = {}  # PMID -> failure reason for this session
        self._session = requests.Session()
        self._default_headers = {
            "User-Agent": "fly-stocker-v2/1.0 (+https://github.com/Allada-Lab)",
            "Accept": "*/*",
        }
    
    def clear_cache(self):
        """Clear the in-memory session cache (not the persistent cache)."""
        self._session_cache.clear()
    
    def save_cache(self):
        """Save pending entries to the persistent cache file."""
        self._method_cache.save_pending()
    
    def _http_get(
        self,
        url: str,
        headers: dict = None,
        params: dict = None,
        timeout: int = 15,
        max_retries: int = 3,
        backoff: float = 0.75
    ):
        """HTTP GET with retry on transient failures; no exceptions on failure."""
        merged_headers = dict(self._default_headers)
        if headers:
            merged_headers.update(headers)
        
        for attempt in range(max_retries):
            try:
                r = self._session.get(
                    url,
                    headers=merged_headers,
                    params=params or {},
                    timeout=timeout,
                    allow_redirects=True,
                )
                if r.status_code == 200:
                    return r
                
                # Retry on transient server/rate-limit errors only.
                if r.status_code in (429, 500, 502, 503, 504):
                    wait_s = backoff * (2 ** attempt) + random.uniform(0, backoff)
                    time.sleep(wait_s)
                    continue
                return None
            except Exception:
                wait_s = backoff * (2 ** attempt) + random.uniform(0, backoff)
                time.sleep(wait_s)
        return None
    
    def _strip_html(self, text: str) -> str:
        """Remove HTML tags from text."""
        return re.sub(r"<[^>]+>", " ", text or "")
    
    def pmid_to_pmcid(self, pmid: str) -> Optional[str]:
        """Convert a single PMID to PMCID using NCBI's ID converter."""
        pmcid, _doi = self._lookup_ids_from_ncbi(pmid)
        return pmcid
    
    def _lookup_ids_from_ncbi(self, pmid: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Resolve PMCID and DOI for a PMID using NCBI ID Converter.
        
        Returns:
            (pmcid, doi) where either value can be None.
        """
        pmid_clean = clean_id(pmid)
        if not pmid_clean:
            return None, None
        
        try:
            url = (
                "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
                f"?tool=flybase_scraper&email=flybase_scraper@example.com&ids={pmid_clean}&format=json"
            )
            if self.ncbi_api_key:
                url += f"&api_key={self.ncbi_api_key}"
            r = self._http_get(url, timeout=15, max_retries=3)
            if r is None:
                return None, None
            data = r.json()
            recs = data.get("records", [])
            if not recs:
                return None, None
            
            rec0 = recs[0]
            pmcid = clean_id(rec0.get("pmcid", ""))
            doi = clean_id(rec0.get("doi", ""))
            return (pmcid or None), (doi or None)
        except:
            return None, None
    
    def _lookup_ids_from_europepmc(
        self,
        pmid: Optional[str] = None,
        doi: Optional[str] = None
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Resolve missing IDs via Europe PMC search endpoint.
        
        Returns:
            (pmcid, doi) where either value can be None.
        """
        pmid_clean = clean_id(pmid) if pmid else ""
        doi_clean = clean_id(doi) if doi else ""
        query = ""
        if pmid_clean:
            query = f"EXT_ID:{pmid_clean} AND SRC:MED"
        elif doi_clean:
            query = f'DOI:"{doi_clean}"'
        else:
            return None, None
        
        r = self._http_get(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={"query": query, "format": "json", "pageSize": 1},
            timeout=20,
            max_retries=3,
        )
        if r is None:
            return None, None
        
        try:
            data = r.json()
            result_list = data.get("resultList", {}).get("result", [])
            if not result_list:
                return None, None
            result = result_list[0]
            pmcid = clean_id(result.get("pmcid", ""))
            doi_val = clean_id(result.get("doi", ""))
            return (pmcid or None), (doi_val or None)
        except Exception:
            return None, None
    
    def _is_pdf_response(self, response) -> bool:
        """Heuristic check for PDF content."""
        if response is None:
            return False
        content_type = str(response.headers.get("content-type", "")).lower()
        if "application/pdf" in content_type or "pdf" in content_type:
            return True
        return bool(response.content[:5] == b"%PDF-")
    
    def _extract_pdf_text(self, content: bytes) -> str:
        """Extract text from raw PDF bytes."""
        if not content or len(content) <= 500 or not PYPDF2_AVAILABLE:
            return ""
        try:
            reader = PdfReader(io.BytesIO(content))
            text = ""
            for page in reader.pages:
                text += page.extract_text() or ""
                text += "\n"
            return text
        except Exception:
            return ""
    
    def _extract_html_text(self, html_text: str, min_chars: int = 300) -> str:
        """Extract text from HTML with basic cleanup."""
        if not html_text or len(html_text) <= 300:
            return ""
        txt = re.sub(r"<script[\s\S]*?</script>", " ", html_text, flags=re.IGNORECASE)
        txt = re.sub(r"<style[\s\S]*?</style>", " ", txt, flags=re.IGNORECASE)
        txt = re.sub(r"\s+", " ", self._strip_html(txt)).strip()
        if len(txt) >= min_chars:
            return txt
        return ""
    
    def _fetch_pmcoa_pdf(self, pmcid: str) -> str:
        """Attempt to download PDF from NCBI PMC and extract text."""
        if not PYPDF2_AVAILABLE:
            return ""
        pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf"
        r = self._http_get(pdf_url, timeout=25)
        if r is not None and r.status_code == 200:
            text = self._extract_pdf_text(r.content)
            if len(text) > 300:
                return text
        return ""
    
    def _fetch_pmcoa_xml(self, pmcid: str) -> str:
        """Attempt to get full text XML from NCBI PMC OAI."""
        xml_url = f"https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:{pmcid.replace('PMC', '')}&metadataPrefix=pmc"
        r = self._http_get(xml_url, headers={"User-Agent": "Mozilla/5.0"}, timeout=25)
        if r is not None and r.status_code == 200 and len(r.content) > 500:
            text = re.sub(r"<[^>]+>", " ", r.text)
            return re.sub(r"\s+", " ", text).strip()
        return ""
    
    def _fetch_europepmc(self, identifier: str) -> str:
        """Try Europe PMC full text XML API."""
        if identifier.lower().startswith("pmc"):
            id_param = "PMC" + identifier.replace("PMC", "")
        else:
            id_param = identifier
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{id_param}/fullTextXML"
        r = self._http_get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=25)
        if r is not None and r.status_code == 200 and len(r.content) > 500:
            text = re.sub(r"<[^>]+>", " ", r.text)
            return re.sub(r"\s+", " ", text).strip()
        return ""
    
    def _fetch_pmc_html(self, pmcid: str) -> str:
        """Fallback to scraping article HTML from PMC page."""
        page_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
        r = self._http_get(page_url, timeout=25)
        if r is not None and r.status_code == 200:
            text = self._extract_html_text(r.text, min_chars=1000)
            if text:
                return text
        return ""
    
    def _fetch_unpaywall(self, doi: str) -> Tuple[str, str]:
        """Return (text, source_label) using Unpaywall OA locations."""
        if not doi or not self.unpaywall_token:
            return "", ""
        
        r = self._http_get(
            f"https://api.unpaywall.org/v2/{quote(doi, safe='')}",
            params={"email": self.unpaywall_token},
            timeout=20
        )
        if r is None:
            return "", ""
        
        try:
            data = r.json()
            locs = data.get("oa_locations", []) or []
            if data.get("best_oa_location"):
                locs = [data["best_oa_location"]] + [l for l in locs if l != data["best_oa_location"]]
            
            for loc in locs:
                url_for_pdf = loc.get("url_for_pdf") or ""
                url = loc.get("url") or ""
                
                # Try PDF first
                if url_for_pdf:
                    pr = self._http_get(url_for_pdf, timeout=25)
                    if pr is not None and pr.status_code == 200 and self._is_pdf_response(pr):
                        text = self._extract_pdf_text(pr.content)
                        if len(text) > 300:
                            return text, "Unpaywall PDF"
                
                # Then try HTML
                if url:
                    hr = self._http_get(url, timeout=20)
                    if hr is not None and hr.status_code == 200:
                        if self._is_pdf_response(hr):
                            txt_pdf = self._extract_pdf_text(hr.content)
                            if len(txt_pdf) > 300:
                                return txt_pdf, "Unpaywall PDF"
                        txt_html = self._extract_html_text(hr.text, min_chars=300)
                        if txt_html:
                            return txt_html, "Unpaywall HTML"
        except Exception:
            return "", ""
        
        return "", ""
    
    def _fetch_openalex(self, doi: str) -> Tuple[str, str]:
        """Try OpenAlex OA links for DOI-based retrieval."""
        if not doi:
            return "", ""
        
        doi_clean = clean_id(doi).lower()
        if doi_clean.startswith("https://doi.org/"):
            doi_clean = doi_clean.replace("https://doi.org/", "", 1)
        elif doi_clean.startswith("http://doi.org/"):
            doi_clean = doi_clean.replace("http://doi.org/", "", 1)
        elif doi_clean.startswith("doi:"):
            doi_clean = doi_clean.replace("doi:", "", 1)
        
        r = self._http_get(
            "https://api.openalex.org/works",
            params={"filter": f"doi:{doi_clean}", "per-page": 1},
            timeout=20,
            max_retries=3,
        )
        if r is None:
            return "", ""
        
        try:
            data = r.json()
            results = data.get("results", []) or []
            if not results:
                return "", ""
            work = results[0]
            oa_url = (
                (work.get("open_access") or {}).get("oa_url")
                or ((work.get("primary_location") or {}).get("pdf_url"))
                or ((work.get("primary_location") or {}).get("landing_page_url"))
            )
            if not oa_url:
                return "", ""
            
            hr = self._http_get(oa_url, timeout=25)
            if hr is None:
                return "", ""
            if self._is_pdf_response(hr):
                text = self._extract_pdf_text(hr.content)
                if len(text) > 300:
                    return text, "OpenAlex PDF"
            txt = self._extract_html_text(hr.text, min_chars=300)
            if txt:
                return txt, "OpenAlex HTML"
        except Exception:
            return "", ""
        
        return "", ""
    
    def _fetch_crossref(self, doi: str) -> Tuple[str, str]:
        """Try DOI links exposed by Crossref metadata."""
        if not doi:
            return "", ""
        doi_encoded = quote(clean_id(doi), safe="")
        r = self._http_get(
            f"https://api.crossref.org/works/{doi_encoded}",
            timeout=20,
            max_retries=3,
        )
        if r is None:
            return "", ""
        
        try:
            message = r.json().get("message", {})
            links = message.get("link", []) or []
            for link in links:
                url = str(link.get("URL", "")).strip()
                if not url:
                    continue
                lr = self._http_get(url, timeout=25)
                if lr is None:
                    continue
                if self._is_pdf_response(lr):
                    text = self._extract_pdf_text(lr.content)
                    if len(text) > 300:
                        return text, "Crossref PDF"
                txt = self._extract_html_text(lr.text, min_chars=300)
                if txt:
                    return txt, "Crossref HTML"
        except Exception:
            return "", ""
        
        return "", ""
    
    def _fetch_doi_resolver(self, doi: str) -> Tuple[str, str]:
        """Last-resort DOI landing page retrieval."""
        if not doi:
            return "", ""
        doi_norm = clean_id(doi)
        if doi_norm.startswith("https://doi.org/") or doi_norm.startswith("http://doi.org/"):
            doi_url = doi_norm
        elif doi_norm.startswith("doi:"):
            doi_url = f"https://doi.org/{doi_norm.replace('doi:', '', 1)}"
        else:
            doi_url = f"https://doi.org/{doi_norm}"
        
        r = self._http_get(doi_url, timeout=25, max_retries=2)
        if r is None:
            return "", ""
        if self._is_pdf_response(r):
            text = self._extract_pdf_text(r.content)
            if len(text) > 300:
                return text, "DOI Resolver PDF"
        txt = self._extract_html_text(r.text, min_chars=300)
        if txt:
            return txt, "DOI Resolver HTML"
        return "", ""
    
    def fetch(
        self,
        pmid: str,
        pmcid: Optional[str] = None,
        doi: Optional[str] = None
    ) -> Tuple[str, str]:
        """
        Fetch full text for a paper.
        
        Uses cascade of strategies, with persistent caching to remember successful methods.
        Failures are cached only in-memory for the session to avoid repeated attempts.
        
        Args:
            pmid: PubMed ID
            pmcid: Optional PMC ID
            doi: Optional DOI
        
        Returns:
            Tuple of (text, source_label). source_label is empty if not found.
        """
        text, label, _reason = self.fetch_with_diagnostics(pmid, pmcid, doi)
        return text, label
    
    def fetch_with_diagnostics(
        self,
        pmid: str,
        pmcid: Optional[str] = None,
        doi: Optional[str] = None
    ) -> Tuple[str, str, str]:
        """
        Fetch full text with source and failure diagnostics.
        
        Returns:
            (text, source_label, failure_reason). failure_reason is empty on success.
        """
        pmid_clean = clean_id(pmid)
        if not pmid_clean:
            return "", "", "invalid_pmid"
        
        # Check session cache for known failures (not persisted)
        if pmid_clean in self._session_cache:
            reason = self._session_cache.get(pmid_clean, "") or "session_cached_failure"
            return "", "", reason
        
        pmcid_clean = clean_id(pmcid)
        doi_clean = clean_id(doi)
        
        # Enrich missing IDs from NCBI and EuropePMC.
        if not pmcid_clean or not doi_clean:
            ncbi_pmcid, ncbi_doi = self._lookup_ids_from_ncbi(pmid_clean)
            if not pmcid_clean and ncbi_pmcid:
                pmcid_clean = ncbi_pmcid
            if not doi_clean and ncbi_doi:
                doi_clean = ncbi_doi
        
        if not pmcid_clean or not doi_clean:
            epmc_pmcid, epmc_doi = self._lookup_ids_from_europepmc(
                pmid=pmid_clean,
                doi=doi_clean
            )
            if not pmcid_clean and epmc_pmcid:
                pmcid_clean = epmc_pmcid
            if not doi_clean and epmc_doi:
                doi_clean = epmc_doi
        
        # Check persistent cache for previously successful method
        cached_method = self._method_cache.get(pmid_clean)
        if cached_method:
            if cached_method == "PMC OA PDF" and pmcid_clean:
                text = self._fetch_pmcoa_pdf(pmcid_clean)
                if len(text) > 300:
                    self._session_cache.pop(pmid_clean, None)
                    return text, cached_method, ""
            elif cached_method == "PMC OA XML" and pmcid_clean:
                text = self._fetch_pmcoa_xml(pmcid_clean)
                if len(text) > 300:
                    self._session_cache.pop(pmid_clean, None)
                    return text, cached_method, ""
            elif cached_method == "Europe PMC" and pmcid_clean:
                text = self._fetch_europepmc(pmcid_clean)
                if len(text) > 300:
                    self._session_cache.pop(pmid_clean, None)
                    return text, cached_method, ""
            elif cached_method == "PMC HTML" and pmcid_clean:
                text = self._fetch_pmc_html(pmcid_clean)
                if len(text) > 300:
                    self._session_cache.pop(pmid_clean, None)
                    return text, cached_method, ""
            elif "Unpaywall" in cached_method and doi_clean:
                text, label = self._fetch_unpaywall(doi_clean)
                if text and len(text) > 300:
                    self._session_cache.pop(pmid_clean, None)
                    return text, label, ""
            elif "OpenAlex" in cached_method and doi_clean:
                text, label = self._fetch_openalex(doi_clean)
                if text and len(text) > 300:
                    self._session_cache.pop(pmid_clean, None)
                    return text, label, ""
            elif "Crossref" in cached_method and doi_clean:
                text, label = self._fetch_crossref(doi_clean)
                if text and len(text) > 300:
                    self._session_cache.pop(pmid_clean, None)
                    return text, label, ""
            elif "DOI Resolver" in cached_method and doi_clean:
                text, label = self._fetch_doi_resolver(doi_clean)
                if text and len(text) > 300:
                    self._session_cache.pop(pmid_clean, None)
                    return text, label, ""
        
        # Full cascade
        pmcid_strategies = [
            ("PMC OA PDF", self._fetch_pmcoa_pdf),
            ("PMC OA XML", self._fetch_pmcoa_xml),
            ("Europe PMC", self._fetch_europepmc),
            ("PMC HTML", self._fetch_pmc_html),
        ]
        if pmcid_clean:
            for label, method in pmcid_strategies:
                try:
                    text = method(pmcid_clean)
                    if len(text) > 300:
                        self._method_cache.set(pmid_clean, label)
                        self._session_cache.pop(pmid_clean, None)
                        return text, label, ""
                except Exception:
                    continue
        
        # DOI-based fallbacks
        if doi_clean:
            doi_strategies = [
                self._fetch_unpaywall,
                self._fetch_openalex,
                self._fetch_crossref,
                self._fetch_doi_resolver,
            ]
            for method in doi_strategies:
                try:
                    text, label = method(doi_clean)
                    if text and len(text) > 300:
                        self._method_cache.set(pmid_clean, label)
                        self._session_cache.pop(pmid_clean, None)
                        return text, label, ""
                except Exception:
                    continue
        
        if not pmcid_clean and not doi_clean:
            reason = "missing_pmcid_and_doi"
        elif not pmcid_clean:
            reason = "missing_pmcid"
        elif not doi_clean:
            reason = "missing_doi"
        else:
            reason = "all_sources_failed"
        
        # Cache failure in session only (not persistent)
        self._session_cache[pmid_clean] = reason
        return "", "", reason
    
    def get_last_failure_reason(self, pmid: str) -> str:
        """Get most recent session failure reason for a PMID."""
        return self._session_cache.get(clean_id(pmid), "")
    
    def fetch_title_abstract(self, pmcid: str) -> Tuple[str, str]:
        """
        Fetch title and abstract for a PMCID using metapub or Entrez.
        
        Returns:
            Tuple of (title, abstract)
        """
        if METAPUB_AVAILABLE:
            try:
                fetcher = PubMedFetcher(
                    email="flybase_scraper@example.com",
                    api_key=self.ncbi_api_key
                )
                art = fetcher.article_by_pmcid(pmcid)
                return (art.title or ""), (art.abstract or "")
            except:
                pass
        
        if BIOPYTHON_AVAILABLE:
            try:
                # Convert PMCID to PMID
                url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=flybase_scraper&email=flybase_scraper@example.com&ids={pmcid}&format=json"
                r = requests.get(url, timeout=15)
                if r.status_code == 200:
                    data = r.json()
                    recs = data.get("records", [])
                    if recs and "pmid" in recs[0]:
                        pmid = recs[0]["pmid"]
                        Entrez.email = "flybase_scraper@example.com"
                        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="xml")
                        records = Entrez.read(handle)
                        handle.close()
                        for article in records.get('PubmedArticle', []):
                            medline = article['MedlineCitation']
                            title = medline['Article'].get('ArticleTitle', '')
                            abstract_data = medline['Article'].get('Abstract', {})
                            abstract_parts = abstract_data.get('AbstractText', [])
                            if isinstance(abstract_parts, list):
                                abstract = ' '.join(str(p) for p in abstract_parts)
                            else:
                                abstract = str(abstract_parts)
                            return str(title), abstract
            except:
                pass
        
        return "", ""


###############################################################################
# Functional Validator
###############################################################################

if PYDANTIC_AVAILABLE:
    class FunctionalValidationResult(BaseModel):
        """Pydantic model for functional validation query results."""
        functional_validity: str
        phenotypes: Optional[str]
        confidence: int
        rationale: str


class FunctionalValidator:
    """
    Uses OpenAI to determine if a stock was functionally validated in a paper.
    
    The validator:
    1. Searches paper text for stock/allele mentions
    2. Extracts relevant paragraphs around mentions
    3. Queries OpenAI with the context
    4. Returns validation result
    """
    
    def __init__(
        self,
        openai_api_key: Optional[str] = None,
        model: str = "gpt-5-mini",
        log_dir: Optional[Path] = None
    ):
        """
        Initialize the validator.
        
        Args:
            openai_api_key: OpenAI API key
            model: OpenAI model to use
            log_dir: Optional directory for logging GPT queries
        """
        self.model = model
        self.log_dir = log_dir
        self._log_counter = 0
        self._system_prompt_saved = False
        
        self._client = None
        if OPENAI_AVAILABLE and openai_api_key:
            self._client = OpenAI(api_key=openai_api_key)
    
    @property
    def is_available(self) -> bool:
        """Check if OpenAI client is configured."""
        return self._client is not None
    
    def _log_query(
        self,
        system_prompt: str,
        user_prompt: str,
        context: Dict[str, str] = None
    ):
        """Log a GPT query to the log directory."""
        if not self.log_dir:
            return
        
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        # Save system prompt once
        if not self._system_prompt_saved:
            system_path = self.log_dir / "system_prompt.txt"
            with open(system_path, 'w', encoding='utf-8') as f:
                f.write(system_prompt)
            self._system_prompt_saved = True
        
        # Save user prompt
        self._log_counter += 1
        query_path = self.log_dir / f"query_{self._log_counter:04d}.txt"
        
        with open(query_path, 'w', encoding='utf-8') as f:
            if context:
                f.write("=" * 60 + "\n")
                f.write("CONTEXT\n")
                f.write("=" * 60 + "\n")
                for key, value in context.items():
                    f.write(f"{key}: {value}\n")
                f.write("=" * 60 + "\n\n")
            f.write(user_prompt)
    
    def _build_search_patterns(
        self,
        stock_number: str,
        allele_symbol: str,
        is_custom_stock: bool = False
    ) -> List[str]:
        """Build search patterns for finding stock/allele mentions in text."""
        patterns = []
        
        # Stock number patterns (skip for custom stocks)
        stock_clean = clean_id(stock_number)
        if stock_clean and stock_clean.upper() != "CUSTOM":
            stock_patterns = [
                stock_clean,
                f"#{stock_clean}",
                f"bl{stock_clean}",
                f"bdsc{stock_clean}",
                f"bdsc_{stock_clean}",
                f"rrid:bdsc_{stock_clean}",
                f"stock {stock_clean}",
                f"bloomington {stock_clean}",
            ]
            patterns.extend(p.lower() for p in stock_patterns)
        
        # Allele symbol patterns
        allele_symbols_raw = str(allele_symbol).strip() if allele_symbol else ""
        allele_list = [s.strip() for s in allele_symbols_raw.split(';') if s.strip()]
        
        for allele_clean in allele_list:
            if not allele_clean:
                continue
            
            allele_patterns = [
                allele_clean,
                re.sub(r'[\[\]]', '', allele_clean),
                re.sub(r'\[([^\]]+)\]', r'-\1', allele_clean),
                re.sub(r'\[([^\]]+)\]', r'^\1', allele_clean),
                re.sub(r'\[([^\]]+)\]', r'(\1)', allele_clean),
            ]
            patterns.extend(p.lower() for p in allele_patterns if p)
            # NOTE: Gene name patterns intentionally NOT included
            # Only stock numbers and allele symbols are searched
        
        return list(set(p for p in patterns if p and len(p) > 1))
    
    def _calculate_context_window(self, mention_count: int) -> int:
        """Calculate context window size based on mention count."""
        if mention_count <= 2:
            return 2000
        elif mention_count <= 5:
            return 1200
        elif mention_count <= 10:
            return 800
        elif mention_count <= 20:
            return 500
        else:
            return 300
    
    def _extract_context_from_patterns(
        self,
        full_text: str,
        patterns: List[str],
        context_chars: int,
        prefer_after_position: Optional[int] = None,
        used_ranges: Optional[List[Tuple[int, int]]] = None
    ) -> Tuple[List[str], List[Tuple[int, int]], int]:
        """
        Extract non-overlapping text windows around pattern matches.
        
        Returns:
            (excerpts, updated_used_ranges, mention_count)
        """
        if not full_text or not patterns:
            return [], used_ranges or [], 0
        
        text_lower = full_text.lower()
        found_positions: List[Tuple[int, int, str]] = []
        
        for pattern in patterns:
            start = 0
            while True:
                pos = text_lower.find(pattern, start)
                if pos == -1:
                    break
                found_positions.append((pos, len(pattern), pattern))
                start = pos + 1
        
        if not found_positions:
            return [], used_ranges or [], 0
        
        found_positions.sort(key=lambda x: x[0])
        mention_count = len(found_positions)
        
        # Optionally prioritize matches that occur in/after a likely "Results" section.
        if prefer_after_position is not None:
            preferred = [p for p in found_positions if p[0] >= prefer_after_position]
            non_preferred = [p for p in found_positions if p[0] < prefer_after_position]
            found_positions = preferred + non_preferred
        
        excerpts: List[str] = []
        ranges = list(used_ranges or [])
        
        for pos, _match_len, _pattern in found_positions:
            is_covered = any(start <= pos <= end for start, end in ranges)
            if is_covered:
                continue
            
            para_start = max(0, pos - context_chars)
            para_end = min(len(full_text), pos + context_chars)
            excerpt = full_text[para_start:para_end].strip()
            if excerpt and len(excerpt) > 50:
                excerpts.append(excerpt)
                ranges.append((para_start, para_end))
        
        return excerpts, ranges, mention_count
    
    def _normalize_gene_terms(self, gene_terms: Optional[List[str]]) -> List[str]:
        """Normalize and deduplicate gene symbols/synonyms for text searching."""
        if not gene_terms:
            return []
        
        normalized: List[str] = []
        seen = set()
        for term in gene_terms:
            t = str(term).strip()
            if not t:
                continue
            variants = [t, re.sub(r"\s+", " ", t)]
            for v in variants:
                k = v.lower()
                if len(k) < 2 or k in seen:
                    continue
                seen.add(k)
                normalized.append(k)
        return normalized
    
    def _extract_relevant_paragraphs(
        self,
        full_text: str,
        stock_number: str,
        allele_symbol: str,
        is_custom_stock: bool = False,
        gene_terms: Optional[List[str]] = None,
        abstract: str = "",
    ) -> Tuple[str, bool, int, int]:
        """
        Extract paragraphs containing stock/allele mentions from full text.
        
        Returns:
            Tuple of
            (extracted_text, stock_found, stock_mention_count, gene_mention_count)
        """
        if not full_text or not full_text.strip():
            return "", False, 0, 0
        
        search_patterns = self._build_search_patterns(stock_number, allele_symbol, is_custom_stock)
        stock_context_chars = self._calculate_context_window(2)
        
        # 1) Extract stock/allele-centric excerpts first.
        stock_excerpts, used_ranges, stock_mention_count = self._extract_context_from_patterns(
            full_text=full_text,
            patterns=search_patterns,
            context_chars=stock_context_chars,
            used_ranges=None,
        )
        stock_found = len(stock_excerpts) > 0
        
        # 2) If stock context is sparse, append gene-centric excerpts.
        gene_terms_norm = self._normalize_gene_terms(gene_terms)
        gene_excerpts: List[str] = []
        gene_mention_count = 0
        
        stock_excerpt_chars = sum(len(x) for x in stock_excerpts)
        if stock_excerpt_chars == 0:
            target_gene_chars = 12000
        elif stock_excerpt_chars < 6000:
            target_gene_chars = 4000
        elif stock_excerpt_chars < 12000:
            target_gene_chars = 2000
        else:
            target_gene_chars = 0
        
        if target_gene_chars > 0 and gene_terms_norm:
            results_pos = full_text.lower().find("results")
            gene_context_chars = 800 if stock_found else 1200
            extracted_gene_excerpts, used_ranges, gene_mention_count = self._extract_context_from_patterns(
                full_text=full_text,
                patterns=gene_terms_norm,
                context_chars=gene_context_chars,
                prefer_after_position=results_pos if results_pos >= 0 else None,
                used_ranges=used_ranges,
            )
            
            # Prefer adding abstract context when gene term appears there.
            abstract_text = str(abstract or "").strip()
            abstract_lower = abstract_text.lower()
            if abstract_text and any(t in abstract_lower for t in gene_terms_norm):
                snippet = abstract_text if len(abstract_text) <= 2500 else abstract_text[:2500] + "..."
                gene_excerpts.append("[---ABSTRACT CONTEXT---]\n\n" + snippet)
            
            for excerpt in extracted_gene_excerpts:
                if sum(len(x) for x in gene_excerpts) >= target_gene_chars:
                    break
                gene_excerpts.append(excerpt)
        
        if not stock_excerpts and not gene_excerpts:
            return "", False, 0, 0
        
        sections: List[str] = []
        if stock_excerpts:
            sections.append(
                "[---STOCK/ALLELE EXCERPTS---]\n\n" +
                "\n\n[---EXCERPT---]\n\n".join(stock_excerpts)
            )
        if gene_excerpts:
            sections.append(
                "[---GENE-CONTEXT EXCERPTS---]\n\n" +
                "\n\n[---EXCERPT---]\n\n".join(gene_excerpts)
            )
        
        extracted_text = "\n\n".join(sections)
        
        # Truncate if too long
        max_chars = 100000
        if len(extracted_text) > max_chars:
            extracted_text = extracted_text[:max_chars] + "\n\n[...truncated...]"
        
        return extracted_text, stock_found, stock_mention_count, gene_mention_count
    
    def check_patterns_exist(
        self,
        full_text: str,
        stock_number: str,
        allele_symbol: str,
        is_custom_stock: bool = False,
        gene_terms: Optional[List[str]] = None
    ) -> bool:
        """
        Check if stock/allele or fallback gene patterns are found in full text.
        
        This is a lightweight check to determine if GPT validation should be called.
        If no stock/allele/gene patterns are found, we can skip the GPT call and
        mark as ambiguous.
        
        Args:
            full_text: Full paper text
            stock_number: Stock number (or 'CUSTOM')
            allele_symbol: Allele symbol
            is_custom_stock: Whether this is a custom/lab stock
        
        Returns:
            True if any stock/allele or gene patterns are found in the text
        """
        if not full_text or not full_text.strip():
            return False
        
        text_lower = full_text.lower()
        search_patterns = self._build_search_patterns(stock_number, allele_symbol, is_custom_stock)
        if gene_terms:
            search_patterns.extend(self._normalize_gene_terms(gene_terms))
        
        for pattern in search_patterns:
            if pattern in text_lower:
                return True
        
        return False
    
    def validate(
        self,
        stock_number: str,
        allele_symbol: str,
        full_text: str,
        title: str = "",
        abstract: str = "",
        is_custom_stock: bool = False,
        gene_terms: Optional[List[str]] = None,
    ) -> Dict[str, any]:
        """
        Validate if a stock was functionally tested in a paper.
        
        Args:
            stock_number: Stock number (or 'CUSTOM')
            allele_symbol: Allele symbol
            full_text: Full paper text
            title: Paper title
            abstract: Paper abstract
            is_custom_stock: Whether this is a custom/lab stock
            gene_terms: Primary gene symbols and synonyms for fallback context
        
        Returns:
            Dict with 'functional_validity', 'phenotypes', 'confidence', 'rationale'
        """
        if not self.is_available:
            return {
                "functional_validity": ValidationStatus.AMBIGUOUS,
                "phenotypes": None,
                "confidence": 0,
                "rationale": "OpenAI not available",
                "from_gpt": False,
            }
        
        if not full_text or not full_text.strip():
            return {
                "functional_validity": ValidationStatus.AMBIGUOUS,
                "phenotypes": None,
                "confidence": 0,
                "rationale": "No text provided",
                "from_gpt": False,
            }
        
        # Extract relevant paragraphs
        extracted_text, stock_found, stock_mention_count, gene_mention_count = self._extract_relevant_paragraphs(
            full_text,
            stock_number,
            allele_symbol,
            is_custom_stock,
            gene_terms=gene_terms,
            abstract=abstract,
        )
        
        if not stock_found and gene_mention_count == 0:
            if is_custom_stock:
                rationale = f"Allele symbol '{allele_symbol}' not found in paper text."
            else:
                rationale = f"Stock number '{stock_number}' and allele symbol '{allele_symbol}' not found in paper text."
            return {
                "functional_validity": ValidationStatus.AMBIGUOUS,
                "phenotypes": None,
                "confidence": 95,
                "rationale": rationale,
                "from_gpt": False,
            }
        
        # Build prompts
        sys_prompt = (
            "You are an expert Drosophila geneticist analyzing scientific papers. "
            "Your task is to determine if a specific fly stock or allele was functionally validated "
            "(i.e., showed a phenotype in experimental results) in the provided paper text."
        )
        
        if is_custom_stock:
            stock_info = f"Allele symbol: {allele_symbol}\n(Custom/lab-generated stock)"
        else:
            stock_info = f"Stock number: {stock_number}\nAllele symbol: {allele_symbol}"
        
        user_prompt = f"""Analyze if the following fly stock/allele was functionally validated:

{stock_info}

TITLE: {title or 'Not provided'}
ABSTRACT: {abstract or 'Not provided'}

Paper excerpts (stock/allele mentions: {stock_mention_count}; gene-context mentions: {gene_mention_count}):

{extracted_text}

Return JSON with:
- functional_validity: One of: "Functionally validated", "Tested, but no phenotype", "Ambiguous functional validity"
- phenotypes: brief description of observed phenotypes (null if none)
- confidence: 0-100 confidence score
- rationale: 1-2 sentence explanation

Classification guidelines:
- "Functionally validated": Clear experimental evidence of a phenotype for THIS specific stock/allele.
- "Tested, but no phenotype": Clear evidence of NO phenotype (negative result).
- "Ambiguous functional validity": Stock appears in text but no phenotype directly attributed. DEFAULT if uncertain.
- If only gene-context excerpts are present (without clear stock/allele evidence), classify as "Ambiguous functional validity".
"""
        
        # Log query
        self._log_query(
            sys_prompt, user_prompt,
            {"stock_number": stock_number, "allele_symbol": allele_symbol, "title": title}
        )
        
        try:
            if PYDANTIC_AVAILABLE:
                completion = self._client.beta.chat.completions.parse(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": sys_prompt},
                        {"role": "user", "content": user_prompt}
                    ],
                    response_format=FunctionalValidationResult,
                    max_completion_tokens=1000,
                    reasoning_effort="medium",
                )
                result = completion.choices[0].message.parsed
                if result:
                    validity = ValidationStatus.normalize(result.functional_validity)
                    return {
                        "functional_validity": validity,
                        "phenotypes": result.phenotypes if result.phenotypes else None,
                        "confidence": max(0, min(100, int(result.confidence))) if result.confidence else 0,
                        "rationale": str(result.rationale) if result.rationale else "",
                        "from_gpt": True,
                    }
            else:
                # Fallback without Pydantic
                completion = self._client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": sys_prompt},
                        {"role": "user", "content": user_prompt}
                    ],
                    max_completion_tokens=1000,
                    reasoning_effort="medium",
                )
                # Parse JSON from response
                import json
                content = completion.choices[0].message.content
                result = json.loads(content)
                validity = ValidationStatus.normalize(result.get("functional_validity", ""))
                return {
                    "functional_validity": validity,
                    "phenotypes": result.get("phenotypes"),
                    "confidence": result.get("confidence", 0),
                    "rationale": result.get("rationale", ""),
                    "from_gpt": True,
                }
        
        except Exception as e:
            return {
                "functional_validity": ValidationStatus.AMBIGUOUS,
                "phenotypes": None,
                "confidence": 0,
                "rationale": f"Classification error: {e}",
                "from_gpt": False,
            }
        
        return {
            "functional_validity": ValidationStatus.AMBIGUOUS,
            "phenotypes": None,
            "confidence": 0,
            "rationale": "No result from API",
            "from_gpt": False,
        }
