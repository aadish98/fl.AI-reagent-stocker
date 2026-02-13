"""
External API clients for fly-stocker-v2.

This module contains classes and functions for interacting with external services:
    - PubMed/Entrez for metadata fetching
    - PMC for full-text retrieval
    - Unpaywall for open access content
    - OpenAI for functional validation queries
"""

from .pubmed import PubMedClient, PubMedCache
from .fulltext import FullTextFetcher, FunctionalValidator, FullTextMethodCache

__all__ = [
    "PubMedClient",
    "PubMedCache", 
    "FullTextFetcher",
    "FullTextMethodCache",
    "FunctionalValidator",
]
