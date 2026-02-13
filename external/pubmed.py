"""
PubMed/Entrez client for fetching metadata and managing cache.

This module provides:
    - PubMedCache: Persistent cache for PubMed metadata (title, abstract)
    - PubMedClient: Client for fetching metadata from PubMed via Entrez
"""

import time
from pathlib import Path
from typing import Dict, List, Optional

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
    Persistent cache for PubMed metadata (title, abstract).
    
    Stores metadata in a CSV file for persistence across sessions.
    New entries are appended incrementally to avoid data loss.
    """
    
    def __init__(self, cache_path: Path):
        """
        Initialize the cache.
        
        Args:
            cache_path: Path to the CSV cache file
        """
        self.cache_path = Path(cache_path)
        self._cache: Dict[str, dict] = {}
        self._loaded = False
    
    def load(self) -> Dict[str, dict]:
        """
        Load the cache from disk.
        
        Returns:
            Dict mapping PMID to {'title': str, 'abstract': str}
        """
        if self._loaded:
            return self._cache
        
        if self.cache_path.exists():
            try:
                cache_df = pd.read_csv(self.cache_path, dtype={'pmid': str})
                cache_df['pmid_clean'] = cache_df['pmid'].apply(clean_id)
                cache_df = cache_df[cache_df['pmid_clean'] != ''].copy()
                
                if len(cache_df) > 0:
                    cache_df['title_str'] = cache_df['title'].apply(lambda x: str(x) if pd.notna(x) else '')
                    cache_df['abstract_str'] = cache_df['abstract'].apply(lambda x: str(x) if pd.notna(x) else '')
                    self._cache = {
                        pmid: {'title': title, 'abstract': abstract}
                        for pmid, title, abstract in zip(
                            cache_df['pmid_clean'],
                            cache_df['title_str'],
                            cache_df['abstract_str']
                        )
                    }
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
            Dict with 'title' and 'abstract', or None if not cached
        """
        self.load()
        pmid_clean = clean_id(pmid)
        return self._cache.get(pmid_clean)
    
    def set(self, pmid: str, metadata: dict):
        """
        Add a single entry to the cache (in memory only).
        
        Args:
            pmid: The PMID
            metadata: Dict with 'title' and 'abstract'
        """
        self.load()
        pmid_clean = clean_id(pmid)
        if pmid_clean:
            self._cache[pmid_clean] = metadata
    
    def save_entries(self, new_entries: Dict[str, dict]):
        """
        Append new entries to the cache file on disk.
        
        Args:
            new_entries: Dict mapping PMID to {'title': str, 'abstract': str}
        """
        if not new_entries:
            return
        
        # Update in-memory cache
        for pmid, metadata in new_entries.items():
            self._cache[clean_id(pmid)] = metadata
        
        # Prepare DataFrame for new entries
        rows = []
        for pmid, metadata in new_entries.items():
            rows.append({
                'pmid': pmid,
                'title': metadata.get('title', ''),
                'abstract': metadata.get('abstract', '')
            })
        
        new_df = pd.DataFrame(rows)
        
        # Ensure parent directory exists
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Check if existing file has an index column for compatibility
        if self.cache_path.exists():
            try:
                existing_cols = pd.read_csv(self.cache_path, nrows=0).columns.tolist()
                if 'Unnamed: 0' in existing_cols or (existing_cols and existing_cols[0].startswith('Unnamed')):
                    new_df.insert(0, 'Unnamed: 0', '')
            except Exception:
                pass
        
        # Append to existing file or create new one
        if self.cache_path.exists():
            new_df.to_csv(self.cache_path, mode='a', header=False, index=False)
        else:
            new_df.to_csv(self.cache_path, index=False)
        
        print(f"    Appended {len(new_entries)} new entries to PubMed cache")
    
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
        Batch fetch title and abstract for PMIDs from PubMed.
        
        Results are cached to avoid re-fetching.
        
        Args:
            pmids: List of PMIDs to fetch
        
        Returns:
            Dict mapping PMID to {'title': str, 'abstract': str}
        """
        results = {}
        pmids_to_fetch = []
        
        # Clean PMIDs and check cache
        for pmid in pmids:
            pmid_str = clean_id(pmid)
            if not pmid_str or not pmid_str.isdigit():
                continue
            
            if self.cache and pmid_str in self.cache:
                cached = self.cache.get(pmid_str)
                if cached:
                    results[pmid_str] = cached
            else:
                pmids_to_fetch.append(pmid_str)
        
        if not pmids_to_fetch:
            return results
        
        if not BIOPYTHON_AVAILABLE:
            print("    Warning: Biopython not installed. Cannot fetch PubMed metadata.")
            return results
        
        print(f"    Fetching metadata for {len(pmids_to_fetch)} PMIDs from PubMed...")
        
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        
        batch_size = 200  # NCBI limit
        new_entries = {}
        
        for i in range(0, len(pmids_to_fetch), batch_size):
            batch = pmids_to_fetch[i:i + batch_size]
            
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
                        
                        # Get title
                        title = medline['Article'].get('ArticleTitle', '')
                        if title:
                            title = str(title)
                        
                        # Get abstract (may be in parts)
                        abstract_data = medline['Article'].get('Abstract', {})
                        abstract_parts = abstract_data.get('AbstractText', [])
                        if isinstance(abstract_parts, list):
                            abstract = ' '.join(str(p) for p in abstract_parts)
                        else:
                            abstract = str(abstract_parts)
                        
                        metadata = {'title': title, 'abstract': abstract}
                        results[pmid] = metadata
                        new_entries[pmid] = metadata
                        
                        # Update cache in memory
                        if self.cache:
                            self.cache.set(pmid, metadata)
                        
                    except (KeyError, TypeError):
                        continue
                
                # Rate limiting - NCBI allows 3 requests/sec with API key, 1/sec without
                time.sleep(0.4)
                
            except Exception as e:
                print(f"    Warning: PubMed fetch failed for batch: {e}")
                continue
        
        # Save newly fetched entries to cache
        if self.cache and new_entries:
            self.cache.save_entries(new_entries)
        
        return results
    
    def fetch_author_date_metadata(self, pmids: List[str]) -> Dict[str, Dict]:
        """
        Fetch author, date, and journal metadata for PMIDs from PubMed.
        
        Args:
            pmids: List of PMIDs to fetch
        
        Returns:
            Dict mapping PMID to {'first_author': str, 'corresponding_author': str, 
                                  'date_published': str, 'title': str, 'journal': str}
        """
        if not BIOPYTHON_AVAILABLE:
            return {}
        
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        
        results = {}
        
        # Clean PMIDs
        pmids_to_fetch = []
        for pmid in pmids:
            pmid_str = str(pmid).strip()
            if pmid_str and pmid_str.isdigit():
                pmids_to_fetch.append(pmid_str)
        
        if not pmids_to_fetch:
            return results
        
        batch_size = 200
        for i in range(0, len(pmids_to_fetch), batch_size):
            batch = pmids_to_fetch[i:i + batch_size]
            
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
                        
                        # Get title
                        title = str(medline['Article'].get('ArticleTitle', ''))
                        
                        # Get journal name
                        journal = ""
                        try:
                            journal_info = medline['Article'].get('Journal', {})
                            journal = str(journal_info.get('Title', ''))
                            if not journal:
                                # Try ISOAbbreviation as fallback
                                journal = str(journal_info.get('ISOAbbreviation', ''))
                        except (KeyError, TypeError):
                            pass
                        
                        # Get authors
                        author_list = medline['Article'].get('AuthorList', [])
                        first_author = ""
                        corresponding_author = ""
                        all_authors = []
                        
                        if author_list:
                            # First author
                            first_author_obj = author_list[0]
                            last_name = first_author_obj.get('LastName', '')
                            first_name = first_author_obj.get('ForeName', '')
                            if last_name:
                                first_author = f"{last_name}, {first_name}" if first_name else last_name
                            
                            # Collect all authors
                            for author_obj in author_list:
                                if isinstance(author_obj, dict):
                                    ln = author_obj.get('LastName', '')
                                    fn = author_obj.get('ForeName', '')
                                    if ln:
                                        all_authors.append(f"{ln}, {fn}" if fn else ln)
                                    
                                    # Look for corresponding author with email
                                    if not corresponding_author:
                                        affiliation_info = author_obj.get('AffiliationInfo', [])
                                        if not isinstance(affiliation_info, list):
                                            affiliation_info = [affiliation_info] if affiliation_info else []
                                        
                                        for aff in affiliation_info:
                                            if isinstance(aff, dict):
                                                affiliation = aff.get('Affiliation', '')
                                                if isinstance(affiliation, str) and '@' in affiliation:
                                                    corresponding_author = f"{ln}, {fn}" if fn else ln
                                                    break
                        
                        # Get publication date
                        date_published = ""
                        try:
                            journal_issue = medline['Article'].get('Journal', {}).get('JournalIssue', {})
                            pub_date = journal_issue.get('PubDate', {})
                            if pub_date:
                                year = pub_date.get('Year', '')
                                month = pub_date.get('Month', '')
                                day = pub_date.get('Day', '')
                                if year:
                                    date_parts = [str(year)]
                                    if month:
                                        date_parts.append(str(month))
                                    if day:
                                        date_parts.append(str(day))
                                    date_published = "-".join(date_parts)
                        except (KeyError, TypeError):
                            pass
                        
                        results[pmid] = {
                            'first_author': first_author,
                            'corresponding_author': corresponding_author,
                            'all_authors': '; '.join(all_authors) if all_authors else '',
                            'date_published': date_published,
                            'title': title,
                            'journal': journal
                        }
                        
                    except (KeyError, TypeError, IndexError):
                        continue
                
                time.sleep(0.4)
                
            except Exception as e:
                print(f"    Warning: PubMed fetch failed for batch: {e}")
                continue
        
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
        elif self.cache and pmid_str in self.cache:
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
