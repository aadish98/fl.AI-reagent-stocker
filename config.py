"""
Configuration, constants, and settings for fly-stocker-v2.

This module contains:
    - ValidationStatus: Validation status constants and priorities
    - Settings: Environment and path configuration
    - Path constants for data files
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from dotenv import load_dotenv


###############################################################################
# Validation Status Constants
###############################################################################

class ValidationStatus:
    """
    Manages validation status constants and priorities.
    
    Simplified priority order (highest to lowest):
    1. Functionally validated - Clear experimental evidence of phenotype
    2. Tested, but no phenotype - Tested but showed no phenotype (negative result)
    3. Ambiguous functional validity - Default for: mentioned but not validated,
       text not accessible, no PMID, non-keyword refs, etc.
    4. No references - Stock has no associated references
    """
    FUNCTIONALLY_VALIDATED = "Functionally validated"
    TESTED_NO_PHENOTYPE = "Tested, but no phenotype"
    UNCERTAIN = "Mentioned but not validated (uncertain)"  # Legacy - maps to AMBIGUOUS
    NOT_FOUND = "Not found in paper text"  # Legacy - maps to AMBIGUOUS
    NO_TEXT = "Full text not accessible"  # Legacy - maps to AMBIGUOUS
    AMBIGUOUS = "Ambiguous functional validity"
    NO_REFS = "No references"
    
    PRIORITIES = {
        FUNCTIONALLY_VALIDATED: 1,
        TESTED_NO_PHENOTYPE: 2,
        AMBIGUOUS: 3,
        UNCERTAIN: 3,  # Maps to same priority as AMBIGUOUS
        NOT_FOUND: 3,  # Maps to same priority as AMBIGUOUS
        NO_TEXT: 3,    # Maps to same priority as AMBIGUOUS
        NO_REFS: 4,
    }
    
    VALID_VALUES = [FUNCTIONALLY_VALIDATED, TESTED_NO_PHENOTYPE, AMBIGUOUS, NO_REFS]
    
    # Map legacy/detailed statuses to output categories
    CATEGORY_MAP = {
        FUNCTIONALLY_VALIDATED: "Functionally Validated",
        TESTED_NO_PHENOTYPE: "Tested But No Phenotype",
        UNCERTAIN: "Ambiguous Functional Validity",
        NOT_FOUND: "Ambiguous Functional Validity",
        NO_TEXT: "Ambiguous Functional Validity",
        AMBIGUOUS: "Ambiguous Functional Validity",
        NO_REFS: "No references",
    }
    
    # Output categories for SplitStocksByValidation
    OUTPUT_CATEGORIES = [
        "Functionally Validated",
        "Ambiguous Functional Validity",
        "Tested But No Phenotype",
        "No references",
    ]
    
    # Priority values for validation categories (lower is better/more validated)
    CATEGORY_PRIORITY = {
        "Functionally Validated": 0,
        "Ambiguous Functional Validity": 1,
        "Tested But No Phenotype": 2,
        "No references": 3,
    }
    
    @classmethod
    def normalize(cls, validity: str) -> str:
        """Normalize a validity string to a standard value."""
        validity = str(validity).strip()
        if validity in cls.VALID_VALUES:
            return validity
        
        # Try to map common variations
        validity_lower = validity.lower()
        if "functionally validated" in validity_lower and "not" not in validity_lower:
            return cls.FUNCTIONALLY_VALIDATED
        elif "tested, but no phenotype" in validity_lower or "no phenotype" in validity_lower:
            return cls.TESTED_NO_PHENOTYPE
        elif "no references" in validity_lower or validity_lower == "no refs":
            return cls.NO_REFS
        else:
            return cls.AMBIGUOUS
    
    @classmethod
    def categorize(cls, validation_status: str) -> str:
        """Map a validation status to one of the 4 output categories."""
        if pd.isna(validation_status) or not str(validation_status).strip():
            return "Ambiguous Functional Validity"
        
        status = str(validation_status).strip()
        return cls.CATEGORY_MAP.get(status, "Ambiguous Functional Validity")
    
    @classmethod
    def get_highest(cls, validity_str: str) -> str:
        """
        Extract the highest level of validation from an aggregated string.
        
        Parses strings like:
        - "PMID:123=Functionally validated; PMID:456=Ambiguous functional validity"
        - "PMID:123: Functionally validated; PMID:456: Ambiguous functional validity"
        and returns the highest priority validation status.
        """
        if not validity_str:
            return ""
        
        best_priority = 999
        best_validity = ""
        
        for part in validity_str.split("; "):
            part = str(part).strip()
            if not part:
                continue
            validity = ""
            if "=" in part:
                validity = part.split("=", 1)[1].strip()
            elif ":" in part:
                # Handle "PMID:<id>: <status>" and similar colon-delimited forms.
                validity = part.rsplit(":", 1)[1].strip()
            if not validity:
                continue
            priority = cls.PRIORITIES.get(validity, 999)
            if priority < best_priority:
                best_priority = priority
                best_validity = validity
        
        return best_validity


###############################################################################
# Path Configuration
###############################################################################

def _find_project_root() -> Path:
    """Find the project root directory (FlyStocker)."""
    current = Path(__file__).resolve()
    # Walk up to find FlyStocker directory
    for parent in current.parents:
        if parent.name == "FlyStocker":
            return parent
    # Fallback to parent of Core
    return Path(__file__).parent.parent.parent


def _get_package_dir() -> Path:
    """Get the directory containing this package (fly_stocker_v2)."""
    return Path(__file__).parent.resolve()


# Package-local data directory (for portable installations)
PACKAGE_DIR = _get_package_dir()
LOCAL_DATA_DIR = PACKAGE_DIR / "data"

# Default paths relative to project structure
PROJECT_ROOT = _find_project_root()
EXTERNAL_DATA_ROOT = PROJECT_ROOT.parent / "Data"  # ../Data from FlyStocker (i.e., Fly_Stock_Scrapers/Data)


def _get_data_root() -> Path:
    """
    Get the data root directory, preferring local data if available.
    
    This enables portable installations - if the local data directory exists
    and contains the required files, use it. Otherwise, fall back to the
    external data directory.
    """
    local_flybase = LOCAL_DATA_DIR / "FlyBase"
    local_bdsc = LOCAL_DATA_DIR / "BDSC" / "stockgenes.csv"
    
    # Check if local data is available
    if LOCAL_DATA_DIR.exists() and local_flybase.exists() and local_bdsc.exists():
        return LOCAL_DATA_DIR
    
    # Fall back to external data
    return EXTERNAL_DATA_ROOT


# Get the active data root
DATA_ROOT = _get_data_root()

# FlyBase data paths
FLYBASE_DATA = DATA_ROOT / "FlyBase"
FLYBASE_ALLELES_STOCKS = FLYBASE_DATA / "Alleles_And_Stocks"
FLYBASE_REFERENCES = FLYBASE_DATA / "FlyBase_References"

# BDSC data paths
BDSC_DATA = DATA_ROOT / "BDSC"
BDSC_STOCKGENES_PATH = BDSC_DATA / "stockgenes.csv"
BDSC_BLOOMINGTON_PATH = BDSC_DATA / "bloomington.csv"
BDSC_STOCKCOMPS_PATH = BDSC_DATA / "stockcomps_map_comments.csv"
BDSC_BALANCERS_PATH = BDSC_DATA / "balancers.csv"

# JSON config paths
JSON_CONFIG_DIR = PACKAGE_DIR / "data" / "JSON"
DEFAULT_SPLIT_CONFIG_PATH = JSON_CONFIG_DIR / "stock_split_config_example.json"

# Network BDSC data paths (fallback if local not available)
NETWORK_BDSC_DATA = Path("/Volumes/umms-rallada/UM Lab Users/Aadish/Data/BDSC")
NETWORK_BDSC_BLOOMINGTON_PATH = NETWORK_BDSC_DATA / "bloomington.csv"
NETWORK_BDSC_STOCKCOMPS_PATH = NETWORK_BDSC_DATA / "stockcomps_map_comments.csv"

# PubMed cache path (on network drive)
# Always read from network location as requested
PUBMED_CACHE_PATH = Path("/Volumes/umms-rallada/UM Lab Users/Aadish/Data/PubMed Cache/pmid_to_title_abstract.csv")

# Full-text method cache path (in same folder as PubMed cache)
FULLTEXT_METHOD_CACHE_PATH = Path("/Volumes/umms-rallada/UM Lab Users/Aadish/Data/PubMed Cache/pmid_to_fulltext_method.csv")

# Helper scripts
HELPER_SCRIPTS = PROJECT_ROOT / "HelperScripts"
GET_FBGN_IDS_SCRIPT = HELPER_SCRIPTS / "GetFBgnIDs.py"

# Logs directory (repo-local)
LOGS_DIR = PACKAGE_DIR / "data" / "logs"
GPT_LOGS_DIR = LOGS_DIR / "GPT-Queries"


def is_portable_mode() -> bool:
    """Check if the package is running in portable mode (using local data)."""
    return DATA_ROOT == LOCAL_DATA_DIR


###############################################################################
# Settings Dataclass
###############################################################################

@dataclass
class Settings:
    """
    Central configuration for fly-stocker-v2 pipelines.
    
    Loads environment variables and provides access to paths and API keys.
    Can be instantiated with custom values or defaults from environment.
    
    Example:
        settings = Settings()
        settings = Settings(batch_size=100, openai_model="gpt-5-mini")
    """
    # API Keys (loaded from environment if not provided)
    ncbi_api_key: Optional[str] = None
    openai_api_key: Optional[str] = None
    unpaywall_token: Optional[str] = None
    
    # OpenAI configuration
    openai_model: str = "gpt-5-mini"
    
    # Processing configuration
    batch_size: int = 50
    
    # Path overrides (use defaults if None)
    flybase_data_path: Optional[Path] = None
    bdsc_stockgenes_path: Optional[Path] = None
    bdsc_bloomington_path: Optional[Path] = None
    bdsc_stockcomps_path: Optional[Path] = None
    bdsc_balancers_path: Optional[Path] = None
    pubmed_cache_path: Optional[Path] = None
    fulltext_method_cache_path: Optional[Path] = None
    split_config_path: Optional[Path] = None
    
    # Logging configuration
    enable_gpt_logging: bool = False
    gpt_log_dir: Optional[Path] = None
    
    # Feature flags
    soft_run: bool = False  # If True, stop before OpenAI calls
    skip_fbgnid_conversion: bool = False
    
    # Functional validation limits
    max_gpt_calls_per_stock: Optional[int] = 5  # Max OpenAI GPT calls per stock (None = no limit)
    
    def __post_init__(self):
        """Load environment variables if API keys not provided."""
        # Load .env file from package directory
        load_dotenv(dotenv_path=PACKAGE_DIR / ".env")
        
        # Load API keys from environment if not provided
        if self.ncbi_api_key is None:
            self.ncbi_api_key = os.getenv("NCBI_API_KEY")
        
        if self.openai_api_key is None:
            self.openai_api_key = os.getenv("OPENAI_API_KEY")
        
        if self.unpaywall_token is None:
            self.unpaywall_token = os.getenv("UNPAYWALL_TOKEN", "aadish98@gmail.com")
        
        # Override model from environment if set
        env_model = os.getenv("OPENAI_MODEL")
        if env_model:
            self.openai_model = env_model
        
        # Set default paths if not provided
        if self.flybase_data_path is None:
            self.flybase_data_path = FLYBASE_DATA
        
        if self.bdsc_stockgenes_path is None:
            self.bdsc_stockgenes_path = BDSC_STOCKGENES_PATH
        
        # Set bloomington.csv path (try local first, then network)
        if self.bdsc_bloomington_path is None:
            if BDSC_BLOOMINGTON_PATH.exists():
                self.bdsc_bloomington_path = BDSC_BLOOMINGTON_PATH
            elif NETWORK_BDSC_BLOOMINGTON_PATH.exists():
                self.bdsc_bloomington_path = NETWORK_BDSC_BLOOMINGTON_PATH
            else:
                self.bdsc_bloomington_path = BDSC_BLOOMINGTON_PATH  # Will fail gracefully
        
        # Set stockcomps_map_comments.csv path (try local first, then network)
        if self.bdsc_stockcomps_path is None:
            if BDSC_STOCKCOMPS_PATH.exists():
                self.bdsc_stockcomps_path = BDSC_STOCKCOMPS_PATH
            elif NETWORK_BDSC_STOCKCOMPS_PATH.exists():
                self.bdsc_stockcomps_path = NETWORK_BDSC_STOCKCOMPS_PATH
            else:
                self.bdsc_stockcomps_path = BDSC_STOCKCOMPS_PATH  # Will fail gracefully
        
        # Set balancers.csv path
        if self.bdsc_balancers_path is None:
            self.bdsc_balancers_path = BDSC_BALANCERS_PATH
        
        # Set default split config path
        if self.split_config_path is None:
            self.split_config_path = DEFAULT_SPLIT_CONFIG_PATH
        
        if self.pubmed_cache_path is None:
            self.pubmed_cache_path = PUBMED_CACHE_PATH
        
        if self.fulltext_method_cache_path is None:
            self.fulltext_method_cache_path = FULLTEXT_METHOD_CACHE_PATH
        
        if self.gpt_log_dir is None:
            self.gpt_log_dir = GPT_LOGS_DIR
    
    @property
    def flybase_alleles_stocks_path(self) -> Path:
        """Path to FlyBase Alleles_And_Stocks directory."""
        return self.flybase_data_path / "Alleles_And_Stocks"
    
    @property
    def flybase_references_path(self) -> Path:
        """Path to FlyBase FlyBase_References directory."""
        return self.flybase_data_path / "FlyBase_References"
    
    def validate(self) -> List[str]:
        """
        Validate configuration and return list of warnings/errors.
        
        Returns:
            List of warning/error messages (empty if all OK)
        """
        issues = []
        
        # Check required paths exist
        if not self.flybase_data_path.exists():
            issues.append(f"FlyBase data directory not found: {self.flybase_data_path}")
        
        if not self.bdsc_stockgenes_path.exists():
            issues.append(f"BDSC stockgenes.csv not found: {self.bdsc_stockgenes_path}")
        
        # Check API keys
        if not self.openai_api_key:
            issues.append("OPENAI_API_KEY not set - functional validation will be disabled")
        
        if not self.ncbi_api_key:
            issues.append("NCBI_API_KEY not set - PubMed queries may be rate-limited")
        
        return issues
    
    def has_openai(self) -> bool:
        """Check if OpenAI API is configured."""
        return bool(self.openai_api_key)
    
    def has_ncbi_api_key(self) -> bool:
        """Check if NCBI API key is configured."""
        return bool(self.ncbi_api_key)
