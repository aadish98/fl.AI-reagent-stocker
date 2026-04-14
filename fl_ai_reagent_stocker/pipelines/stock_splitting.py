"""
Stage 2/3 shared stock organization pipeline.

This module provides the StockSplittingPipeline class that:
1. Loads stock data from Excel files (output of Stage 1)
2. Computes derived columns (Balancers, multiple_insertions, ALLELE_PAPER_RELEVANCE_SCORE)
3. Applies configurable filters from JSON to create sheets
4. Writes organized Excel output
5. Optionally appends Ref++ validation columns when `run_validation=True`

Usage:
    from fl_ai_reagent_stocker.pipelines.stock_splitting import StockSplittingPipeline
    from fl_ai_reagent_stocker.config import Settings

    settings = Settings()
    pipeline = StockSplittingPipeline(settings)
    pipeline.run(input_dir, config_path="path/to/config.json")
"""

import json
import re
import shutil
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd

from ..config import (
    DEFAULT_SPLIT_CONFIG_PATH,
    Settings,
    ValidationStatus,
    normalize_phenotype_similarity_targets,
)
from ..integrations.pubmed import PubMedClient, PubMedCache
from ..integrations.fulltext import FullTextFetcher, FunctionalValidator
from ..integrations.phenotype_similarity import (
    EmbeddingSimilarityScorer,
    PhenotypeSimilarityTarget,
    build_similarity_targets,
    normalize_phenotype_text,
    normalize_qualifier_text,
    plot_similarity_outputs,
)
from .stock_finding import _derive_one_hot_reagent_buckets
from ..utils import (
    clean_id,
    parse_semicolon_list,
    unique_join,
    find_keyword_column,
    find_keyword_title_abstract_column,
    extract_keywords_from_column,
    extract_keywords_from_title_abstract_column,
    has_keyword_validation,
    generate_keyword_column_name,
    get_gpt_derived_columns,
    apply_experimental_prefix,
    EXPERIMENTAL_PREFIX,
    find_latest_tsv,
    load_flybase_tsv,
    REAGENT_BUCKET_COLUMNS,
)
from ..validation_runner import run_functional_validation


# Human-readable filter type descriptions for Contents sheet
FILTER_TYPE_PHRASES = {
    "contains": "contains",
    "doesnt_contain": "doesn't contain",
    "equals": "equals",
    "doesnt_equal": "doesn't equal",
    "insertion_site": "has single insertion at",
}

PHENOTYPE_SIMILARITY_EMBEDDING_MODEL = "text-embedding-3-large"
SIMILARITY_TIER_BIN_WIDTH = 0.05
SIMILARITY_TIER_SHEET_COUNT = int(round(1.0 / SIMILARITY_TIER_BIN_WIDTH))
SOURCE_STOCK_COLUMN = "Source/ Stock #"
SOURCE_COLUMN = "Source"
STOCK_NUMBER_COLUMN = "Stock #"
CO_REAGENT_FBIDS_COLUMN = "Co-reagent FBids"
CO_REAGENT_SYMBOLS_COLUMN = "Co-reagent symbols"
PARTNER_DRIVER_SYMBOLS_COLUMN = "Partner driver symbols (best-effort)"
PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN = "Partner driver stock candidates"
MASTERLIST_TEMPLATE_COLUMNS = [
    "Screening Group",
    "Actual cross set date",
    "Finished (=27+B)",
    "Projected Start Date",
    "Projected End Date",
    "Notes/Function",
    "allele shorthand",
    "Gene",
    "Stock Source",
    "ID #",
    "RNAi Shorthand",
    "Full Stock Genotype",
    "balancers in stock?",
    "ordered? (date)",
    "location of stock",
    "Published Gal4/ Positive control",
    "Published Gal4 source",
    "Published phenotype",
    "Reference",
    "Column 31",
    "Column 30",
    "Column 29",
    "Column 1",
    "Column 32",
    "Circadian/Sleep Relevance (embedding max score)",
    "Experimental Driver",
    "Gal4 control",
    "RNAi control",
    "Data Set",
    "Sleep promoting?",
    "Wake promoting?",
    "link to graph     ",
    "   link to data",
]
MASTERLIST_TEMPLATE_SOURCE_COLUMNS = {
    "allele shorthand": "Reagent Type or Allele Symbol",
    "Gene": "Gene",
    "Stock Source": SOURCE_COLUMN,
    "ID #": STOCK_NUMBER_COLUMN,
    "Full Stock Genotype": "Genotype",
    "balancers in stock?": "Balancers",
    "Published Gal4/ Positive control": PARTNER_DRIVER_SYMBOLS_COLUMN,
    "Published Gal4 source": PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN,
    "Published phenotype": "Phenotype",
    "Reference": "PMID",
    "Column 31": "PMCID",
    "Column 30": "Reference",
    "Column 29": "Authors",
    "Column 1": "Journal",
    "Column 32": "Year of Publication",
    "Data Set": "Data Set",
}


###############################################################################
# JSON Configuration Loading and Validation
###############################################################################

def load_split_config(config_path: Path) -> Dict[str, Any]:
    """
    Load and validate a stock split configuration JSON file.
    
    Args:
        config_path: Path to the JSON configuration file
        
    Returns:
        Parsed configuration dictionary
        
    Raises:
        ValueError: If configuration is invalid
    """
    with open(config_path, 'r', encoding='utf-8') as f:
        config = json.load(f)
    if not isinstance(config, dict):
        raise ValueError("Config root must be a JSON object")
    
    # Validate required keys
    required_keys = ['filters', 'combinations']
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Config missing required key: {key}")

    filters = config.get('filters')
    combinations = config.get('combinations')
    if not isinstance(filters, dict):
        raise ValueError("Config 'filters' must be an object mapping names to filter specs")
    if not isinstance(combinations, list):
        raise ValueError("Config 'combinations' must be a list of filter-name lists")
    
    # Set defaults for optional settings
    if 'settings' not in config:
        config['settings'] = {}
    if not isinstance(config['settings'], dict):
        raise ValueError("Config 'settings' must be an object")
    if 'filterDescriptions' in config and not isinstance(config['filterDescriptions'], dict):
        raise ValueError("Config 'filterDescriptions' must be an object when provided")
    
    settings = config['settings']
    settings.setdefault('relevantSearchTerms', [])
    if 'phenotypeSimilarityTargets' not in settings:
        raise ValueError("Config 'settings.phenotypeSimilarityTargets' is required")
    settings['phenotypeSimilarityTargets'] = normalize_phenotype_similarity_targets(
        settings.get('phenotypeSimilarityTargets')
    )
    settings.setdefault('maxStocksPerGene', None)
    settings.setdefault('maxStocksPerAllele', None)
    if not isinstance(settings['relevantSearchTerms'], list):
        raise ValueError("settings.relevantSearchTerms must be a list")

    valid_filter_types = set(FILTER_TYPE_PHRASES)
    for filter_name, spec in filters.items():
        if not isinstance(spec, dict):
            raise ValueError(f"Filter '{filter_name}' must be an object")
        if not str(spec.get("column", "")).strip():
            raise ValueError(f"Filter '{filter_name}' is missing a non-empty 'column'")
        filter_type = spec.get("type")
        if filter_type not in valid_filter_types:
            raise ValueError(
                f"Filter '{filter_name}' has unsupported type '{filter_type}'. "
                f"Expected one of: {', '.join(sorted(valid_filter_types))}"
            )
        if filter_type in {"contains", "doesnt_contain", "equals", "doesnt_equal", "insertion_site"} and "value" not in spec:
            raise ValueError(f"Filter '{filter_name}' is missing required key 'value'")

    for idx, combo in enumerate(combinations, start=1):
        if not isinstance(combo, list) or not combo:
            raise ValueError(f"Combination #{idx} must be a non-empty list of filter names")
        undefined = [name for name in combo if name not in filters]
        if undefined:
            raise ValueError(
                f"Combination #{idx} references undefined filter(s): {', '.join(undefined)}"
            )
    
    return config


def filter_spec_to_string(spec: Dict[str, Any]) -> str:
    """Convert a filter spec dict to a human-readable string."""
    col = spec.get("column", "?")
    typ = spec.get("type", "?")
    val = spec.get("value")
    
    if val is True:
        val_str = "true"
    elif val is False:
        val_str = "false"
    elif isinstance(val, str):
        val_str = f'"{val}"'
    else:
        val_str = str(val)
    
    phrase = FILTER_TYPE_PHRASES.get(typ, typ)
    return f"{col} {phrase} {val_str}"


###############################################################################
# Column Computation Functions
###############################################################################

def compute_balancers_column(df: pd.DataFrame) -> pd.Series:
    """
    Compute the Balancers column.

    Balancers must come from the FlyBase-derived Stage 1 workbook, where they
    are backed by `FBba` component associations. No genotype parsing fallback
    is allowed.
    """
    if "Balancers" not in df.columns:
        raise ValueError(
            "Input stocks data must include a 'Balancers' column derived from "
            "FlyBase FBba components."
        )
    existing = df["Balancers"].fillna("").astype(str).str.strip()
    existing = existing.mask(existing.eq(""), "-")
    existing = existing.mask(existing.str.lower().eq("nan"), "-")
    return existing


def normalize_num_balancers_column(df: pd.DataFrame) -> pd.Series:
    """Normalize the precomputed FBba-derived balancer count column."""
    if "num_Balancers" not in df.columns:
        raise ValueError(
            "Input stocks data must include a 'num_Balancers' column derived "
            "from FlyBase FBba components."
        )
    values = pd.to_numeric(df["num_Balancers"], errors="coerce").fillna(0)
    return values.astype(int)


def compute_multiple_insertions_column(df: pd.DataFrame) -> pd.Series:
    """
    Compute the multiple_insertions column using stock-level construct data.

    Prefers the `all_stock_constructs` column from Pipeline 1, which is now
    populated from FlyBase-derived stock components. If component data is not
    available, falls back to counting distinct construct-like elements in the
    genotype string via regex.
    
    Returns True if the stock has more than one unique construct, False otherwise.
    """
    # Prefer component-based detection from Pipeline 1 stock-component summaries
    constructs_col = None
    for col in ['all_stock_constructs']:
        if col in df.columns:
            constructs_col = col
            break
    
    if constructs_col is not None:
        def _count_from_components(constructs_str) -> bool:
            if pd.isna(constructs_str) or not str(constructs_str).strip():
                return False
            constructs = [c.strip() for c in str(constructs_str).split(';') if c.strip()]
            return len(set(constructs)) > 1
        
        return df[constructs_col].apply(_count_from_components)
    
    # Fallback: parse genotype for distinct P-element / PBac / Mi / TI constructs
    genotype_col = None
    for col in ['genotype', 'stock_genotype', 'Genotype', 'FB_genotype']:
        if col in df.columns:
            genotype_col = col
            break
    
    if genotype_col is None:
        return pd.Series([False] * len(df), index=df.index)
    
    def _count_from_genotype(genotype: str) -> bool:
        if pd.isna(genotype) or not str(genotype).strip():
            return False
        # Match distinct construct elements: P{...}Name, PBac{...}Name, Mi{...}Name, TI{...}Name
        elements = re.findall(r'(?:P|PBac|Mi|TI|M)\{[^}]+\}\S*', str(genotype))
        return len(set(elements)) > 1
    
    return df[genotype_col].apply(_count_from_genotype)


def compute_allele_paper_relevance_score(
    df: pd.DataFrame,
    keywords: List[str]
) -> pd.Series:
    """
    Compute ALLELE_PAPER_RELEVANCE_SCORE based on keyword hits in publication
    title / abstract.

    Score meanings (mapped to filter labels in the split config):
        0  (Ref-):  Allele has **no** surviving references.
        1  (Ref+):  Allele has paper references, but **none** contain any of
                     the search keywords in their title or abstract.
        2  (Ref++): Allele has at least one paper reference whose title or
                     abstract contains one or more of the search keywords.

    The determination is based entirely on the title/abstract keyword-hit
    column (``"<kw> in title/ abstract of publication?"``) and/or the
    ``keyword_ref_count`` column produced by ``pipeline_references``.  It does
    **not** depend on GPT / functional-validity output.
    """
    # ── Locate the title/abstract keyword boolean column ──────────────
    kw_title_col = find_keyword_title_abstract_column(df)

    # ── Locate keyword_ref_count (integer, from pipeline_references) ──
    # pipeline_references renames this column to "<kw> references count",
    # so use the keywords to construct the exact renamed column name.
    kw_count_col = None
    if 'keyword_ref_count' in df.columns:
        kw_count_col = 'keyword_ref_count'
    elif keywords:
        renamed_col = f"{' OR '.join(keywords)} references count"
        if renamed_col in df.columns:
            kw_count_col = renamed_col

    # ── Locate PMID column ────────────────────────────────────────────
    pmid_col = None
    for col in ['PMID', 'pmid', 'PMIDs']:
        if col in df.columns:
            pmid_col = col
            break
    
    # Optional counts for references lacking PMIDs.
    no_pmid_count_col = 'references_without_pmid_count' if 'references_without_pmid_count' in df.columns else None
    total_refs_col = 'total_refs' if 'total_refs' in df.columns else None

    def calculate_score(row) -> int:
        # --- Score 2: at least one keyword hit in title/abstract ------
        has_keyword_ref = False

        # Primary: boolean column "<kw> in title/ abstract of publication?"
        if kw_title_col and kw_title_col in row.index:
            val = row[kw_title_col]
            if val is True or str(val).lower() in ('true', 'yes'):
                has_keyword_ref = True

        # Secondary: keyword_ref_count > 0
        if not has_keyword_ref and kw_count_col and kw_count_col in row.index:
            try:
                if int(row[kw_count_col]) > 0:
                    has_keyword_ref = True
            except (ValueError, TypeError):
                pass

        if has_keyword_ref:
            return 2  # Ref++

        # --- Score 1: has references but no keyword hit ---------------
        if total_refs_col and total_refs_col in row.index:
            try:
                if int(float(row[total_refs_col])) > 0:
                    return 1  # Ref+
            except (ValueError, TypeError):
                pass
        if no_pmid_count_col and no_pmid_count_col in row.index:
            try:
                if int(float(row[no_pmid_count_col])) > 0:
                    return 1  # Ref+
            except (ValueError, TypeError):
                pass
        if pmid_col and pmid_col in row.index:
            pmid_val = row[pmid_col]
            if pd.notna(pmid_val):
                pmid_str = str(pmid_val).strip()
                if pmid_str and pmid_str != "-" and pmid_str.lower() != "nan":
                    return 1  # Ref+

        return 0  # Ref-

    return df.apply(calculate_score, axis=1)


###############################################################################
# Filter Application Logic
###############################################################################

def apply_filter(
    df: pd.DataFrame,
    filter_spec: Dict[str, Any]
) -> pd.Series:
    """
    Apply a single filter specification to a DataFrame.
    
    Returns a boolean Series indicating which rows match the filter.
    """
    column = filter_spec.get("column")
    filter_type = filter_spec.get("type")
    value = filter_spec.get("value")
    
    if column not in df.columns:
        # Column doesn't exist - return all False (no matches)
        print(f"      Warning: Filter column '{column}' not found in data")
        return pd.Series([False] * len(df), index=df.index)
    
    col_data = df[column]
    
    if filter_type == "contains":
        return col_data.astype(str).str.contains(str(value), case=True, na=False)
    
    elif filter_type == "doesnt_contain":
        return ~col_data.astype(str).str.contains(str(value), case=True, na=False)
    
    elif filter_type == "equals":
        # Handle boolean and numeric comparisons
        if isinstance(value, bool):
            return col_data == value
        elif isinstance(value, (int, float)):
            return col_data == value
        else:
            return col_data.astype(str) == str(value)
    
    elif filter_type == "doesnt_equal":
        if isinstance(value, bool):
            return col_data != value
        elif isinstance(value, (int, float)):
            return col_data != value
        else:
            return col_data.astype(str) != str(value)
    
    elif filter_type == "insertion_site":
        # Single-construct stock at specific site (contains value AND only one construct)
        contains_site = col_data.astype(str).str.contains(str(value), case=True, na=False)
        single_insertion = ~df.get('multiple_insertions', pd.Series([False] * len(df), index=df.index))
        return contains_site & single_insertion
    
    else:
        print(f"      Warning: Unknown filter type '{filter_type}'")
        return pd.Series([False] * len(df), index=df.index)


def apply_filter_combination(
    df: pd.DataFrame,
    combination: List[str],
    filters_config: Dict[str, Dict[str, Any]]
) -> pd.DataFrame:
    """
    Apply a combination of filters (ANDed together) to get matching stocks.
    
    Args:
        df: DataFrame to filter
        combination: List of filter names to apply
        filters_config: Dict mapping filter names to specifications
        
    Returns:
        Filtered DataFrame containing only matching rows
    """
    mask = pd.Series([True] * len(df), index=df.index)
    
    for filter_name in combination:
        if filter_name not in filters_config:
            print(f"      Warning: Filter '{filter_name}' not defined in config")
            continue
        
        filter_mask = apply_filter(df, filters_config[filter_name])
        mask = mask & filter_mask
    
    return df[mask].copy()


###############################################################################
# Stock Limiting Logic
###############################################################################

def apply_stock_limits(
    df: pd.DataFrame,
    max_stocks_per_gene: Optional[int],
    max_stocks_per_allele: Optional[int],
    gene_counters: Dict[str, int],
    allele_counters: Dict[str, int],
    seen_stocks: Set[str],
    verbose: bool = True
) -> pd.DataFrame:
    """
    Apply per-gene and per-allele stock limits.
    
    Modifies counters and seen_stocks in place.
    
    Args:
        df: DataFrame of stocks to limit
        max_stocks_per_gene: Maximum stocks per gene (None = no limit)
        max_stocks_per_allele: Maximum stocks per allele (None = no limit)
        gene_counters: Dict tracking stocks per gene (modified in place)
        allele_counters: Dict tracking stocks per allele (modified in place)
        seen_stocks: Set of already-seen stock IDs (modified in place)
        verbose: Print progress
        
    Returns:
        Limited DataFrame
    """
    
    # Find relevant columns
    stock_col = _find_stock_id_column(df)
    gene_col = _find_gene_id_column(df)
    allele_col = _find_allele_column(df)
    
    selected_indices = []
    
    for idx, row in df.iterrows():
        stock_id = str(row.get(stock_col, '')) if stock_col else str(idx)
        gene_id = str(row.get(gene_col, '')) if gene_col else ''
        allele_id = str(row.get(allele_col, '')) if allele_col else ''
        
        # Skip if already seen
        if stock_id in seen_stocks:
            continue
        
        # Check gene limit
        if max_stocks_per_gene is not None and gene_id:
            if gene_counters.get(gene_id, 0) >= max_stocks_per_gene:
                continue
        
        # Check allele limit
        if max_stocks_per_allele is not None and allele_id:
            if allele_counters.get(allele_id, 0) >= max_stocks_per_allele:
                continue
        
        # Accept this stock
        selected_indices.append(idx)
        seen_stocks.add(stock_id)
        
        if gene_id:
            gene_counters[gene_id] = gene_counters.get(gene_id, 0) + 1
        if allele_id:
            allele_counters[allele_id] = allele_counters.get(allele_id, 0) + 1
    
    return df.loc[selected_indices].copy()


def _find_stock_id_column(df: pd.DataFrame) -> Optional[str]:
    """Find the stock ID column."""
    for col in ['stock_number', 'StockID', 'stock_id', 'FBst']:
        if col in df.columns:
            return col
    return None


def _find_gene_id_column(df: pd.DataFrame) -> Optional[str]:
    """Find the gene ID column."""
    for col in ['flybase_gene_id', 'relevant_gene_symbols', 'FBgn', 'gene_id', 'relevant_flybase_gene_ids', 'all_gene_symbols']:
        if col in df.columns:
            return col
    return None


def _find_allele_column(df: pd.DataFrame) -> Optional[str]:
    """Find the allele column."""
    for col in ['AlleleSymbol', 'flybase_allele_id', 'AlleleID', 'allele_symbol', 'relevant_component_symbols']:
        if col in df.columns:
            return col
    return None


###############################################################################
# Excel Output Functions
###############################################################################

def _apply_grey_fill_xlsxwriter(
    worksheet,
    df: pd.DataFrame,
    grey_format,
    dark_header_format=None,
    dark_header_columns: Optional[Set[str]] = None,
) -> None:
    """
    Apply light grey background to GPT-derived column **headers** via
    xlsxwriter.  Only the header row (row 0) is shaded.

    GPT-derived columns are detected by the '[EXPERIMENTAL] ' prefix that
    :func:`apply_experimental_prefix` adds to their header names.
    """
    dark_header_columns = dark_header_columns or set()
    for col_idx, col_name in enumerate(df.columns):
        col_name_str = str(col_name)
        if dark_header_format is not None and col_name_str in dark_header_columns:
            worksheet.write(0, col_idx, col_name, dark_header_format)
        elif col_name_str.startswith(EXPERIMENTAL_PREFIX):
            worksheet.write(0, col_idx, col_name, grey_format)


def _find_column_with_experimental_fallback(
    columns: List[str],
    base_name: str
) -> Optional[str]:
    """Find a column by exact name or [EXPERIMENTAL]-prefixed variant."""
    if base_name in columns:
        return base_name
    exp_name = f"{EXPERIMENTAL_PREFIX}{base_name}"
    if exp_name in columns:
        return exp_name
    return None


def _find_keyword_ref_count_column(df: pd.DataFrame) -> Optional[str]:
    """Locate per-stock keyword-hit reference count column."""
    if 'keyword_ref_count' in df.columns:
        return 'keyword_ref_count'
    for col in df.columns:
        if str(col).endswith(" references count"):
            return col
    return None


def _find_keyword_ref_pmids_column(df: pd.DataFrame) -> Optional[str]:
    """Locate per-stock keyword-hit PMID list column."""
    if 'keyword_ref_pmids' in df.columns:
        return 'keyword_ref_pmids'
    for col in df.columns:
        if re.match(r"^Stock \(.+\) .+ references$", str(col)):
            return col
    return None


def _normalize_pmid_list_string(value: Any) -> str:
    """Normalize semicolon-delimited PMID-like strings for sheet display."""
    items = [clean_id(v) for v in str(value or "").split(";")]
    deduped: List[str] = []
    seen: Set[str] = set()
    for item in items:
        if not item or item in seen:
            continue
        seen.add(item)
        deduped.append(item)
    return "; ".join(deduped) if deduped else "-"


def _extract_row_pmid_entry(cell_value: Any, row_pmid: str) -> Any:
    """
    Keep only the `PMID: <id>: ...` entry matching the current row PMID.
    Returns original value when no PMID-tagged structure is present.
    """
    text = str(cell_value or "").strip()
    if not text or "PMID:" not in text:
        return cell_value
    pmid_clean = clean_id(row_pmid)
    if not pmid_clean:
        return cell_value

    pattern = re.compile(
        r"PMID:\s*([0-9]+)\s*:\s*(.*?)(?=(?:\s*[|;]\s*)?PMID:\s*[0-9]+\s*:|$)",
        re.IGNORECASE | re.DOTALL,
    )
    matches = pattern.findall(text)
    if not matches:
        return cell_value

    for pmid_found, content in matches:
        if clean_id(pmid_found) == pmid_clean:
            entry = str(content).strip()
            return entry if entry else f"PMID: {pmid_clean}"
    return ""


def _is_keyword_validation_pmid_column(col_name: str) -> bool:
    """Return True for dynamic PMID-of-keyword-valid-references column names."""
    col = str(col_name or "").strip()
    if col.startswith(EXPERIMENTAL_PREFIX):
        col = col[len(EXPERIMENTAL_PREFIX):].strip()
    return bool(
        re.match(
            r"^PMID of .+ references that showed stocks functional validity$",
            col
        )
    )


def _derive_stock_keyword_refs_header(keyword_pmids_col: Optional[str]) -> str:
    """
    Build a Stock Sheet by Gene keyword-reference column header.

    If source column already encodes explicit keywords (for example
    "Stock (allele) sleep OR circadian references"), preserve that wording.
    """
    base = "Stock (FBal / FBtp / FBti) keyword-hit references"
    raw = str(keyword_pmids_col or "").strip()
    if raw:
        if raw.startswith("Stock (") and "references" in raw:
            base = raw.replace("(all for stock)", "").strip()
        elif raw == "keyword_ref_pmids":
            base = "Stock (FBal / FBtp / FBti) keyword-hit references"
    return f"{base} (all for stock)"


def _load_gene_synonyms_map(
    flybase_data_path: Optional[Path],
    verbose: bool = False
) -> Dict[str, str]:
    """
    Load FlyBase gene symbol/synonym table and build symbol -> synonyms display.

    Maps both current symbols and known synonyms to a normalized semicolon list.
    """
    def _find_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
        lower_map = {str(c).strip().lower(): c for c in df.columns}
        for cand in candidates:
            col = lower_map.get(str(cand).strip().lower())
            if col is not None:
                return col
        return None

    def _find_synonym_tsv(base_path: Path) -> Optional[Path]:
        genes_dir = base_path / "genes"
        if genes_dir.exists():
            try:
                return find_latest_tsv(genes_dir, "fb_synonym")
            except FileNotFoundError:
                pass
            except Exception:
                pass

        # Fallback: some local data layouts do not include a top-level "Genes"
        # folder, so search within FlyBase for any fb_synonym TSV/GZ.
        candidates = sorted(base_path.rglob("fb_synonym*.tsv*"))
        if candidates:
            return max(candidates, key=lambda p: p.stat().st_mtime)
        return None

    if flybase_data_path is None:
        return {}
    try:
        synonym_path = _find_synonym_tsv(Path(flybase_data_path))
        if synonym_path is None:
            if verbose:
                print(
                    f"    Warning: Could not find fb_synonym TSV under {flybase_data_path}"
                )
            return {}
        syn_df = load_flybase_tsv(synonym_path, keep_default_na=False)
    except Exception as e:
        if verbose:
            print(f"    Warning: Could not load FlyBase gene synonyms: {e}")
        return {}

    organism_col = _find_col(syn_df, ["organism_abbreviation"])
    current_symbol_col = _find_col(syn_df, ["current_symbol"])
    synonyms_col = _find_col(syn_df, ["symbol_synonym(s)", "symbol_synonyms"])
    primary_fbid_col = _find_col(syn_df, ["primary_FBid", "primary_fbid"])

    if organism_col:
        syn_df = syn_df[
            syn_df[organism_col]
            .astype(str)
            .str.strip()
            .str.lower()
            .eq("dmel")
        ].copy()
    if len(syn_df) == 0 or current_symbol_col is None:
        return {}

    fbgn_to_names: Dict[str, Set[str]] = {}
    for _, row in syn_df.iterrows():
        fbgn = str(row.get(primary_fbid_col, "") or "").strip() if primary_fbid_col else ""
        if not fbgn:
            continue
        names = fbgn_to_names.setdefault(fbgn, set())
        current_symbol = str(row.get(current_symbol_col, "") or "").strip()
        if current_symbol:
            names.add(current_symbol)
        syn_field = str(row.get(synonyms_col, "") or "").strip() if synonyms_col else ""
        if syn_field:
            for syn in syn_field.split("|"):
                syn = syn.strip()
                if syn:
                    names.add(syn)

    name_to_synonyms: Dict[str, str] = {}
    for names in fbgn_to_names.values():
        normalized = sorted({n for n in names if n})
        if not normalized:
            continue
        for name in normalized:
            others = [n for n in normalized if n != name]
            name_to_synonyms[name] = "; ".join(others) if others else ""

    # Add case-insensitive lookup keys so stock-sheet symbols still resolve even
    # when casing differs from FlyBase's synonym table.
    for name, syns in list(name_to_synonyms.items()):
        key = str(name).strip().casefold()
        if key and key not in name_to_synonyms:
            name_to_synonyms[key] = syns
    if verbose:
        print(f"    Loaded FlyBase gene synonyms for {len(name_to_synonyms)} symbols")
    return name_to_synonyms


def _lookup_gene_synonyms(
    gene_label: str,
    gene_synonyms_map: Optional[Dict[str, str]],
) -> str:
    """Resolve synonyms for a gene symbol with robust normalization."""
    if not gene_synonyms_map:
        return ""
    label = str(gene_label or "").strip()
    if not label:
        return ""
    direct = str(gene_synonyms_map.get(label, "") or "").strip()
    if direct:
        return direct
    # Fallback to case-insensitive key.
    return str(gene_synonyms_map.get(label.casefold(), "") or "").strip()


def _normalize_stock_sheet_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Drop requested columns and normalize empty keyword-reference PMID lists."""
    if df is None or len(df) == 0:
        return df
    out = df.copy()
    drop_cols = ['FBrf without PMID', 'references_without_pmid_count']
    out = out.drop(columns=[c for c in drop_cols if c in out.columns], errors='ignore')

    kw_pmids_col = _find_keyword_ref_pmids_column(out)
    if kw_pmids_col and kw_pmids_col in out.columns:
        out[kw_pmids_col] = out[kw_pmids_col].apply(
            lambda v: _normalize_pmid_list_string(v) if str(v or "").strip() else "-"
        )
    return out


def _build_phenotype_similarity_context(
    config: Dict[str, Any],
    pipeline_settings: Optional[Settings],
    verbose: bool = False,
) -> Tuple[
    List[PhenotypeSimilarityTarget],
    Optional[EmbeddingSimilarityScorer],
]:
    """Instantiate similarity targets and scorers for the soft-run sheet."""
    config_settings = config.get("settings", {})
    raw_targets = config_settings.get("phenotypeSimilarityTargets")
    targets = build_similarity_targets(
        normalize_phenotype_similarity_targets(raw_targets)
    )

    embedding_scorer: Optional[EmbeddingSimilarityScorer] = None

    if pipeline_settings is None:
        return targets, embedding_scorer

    if pipeline_settings.soft_run and pipeline_settings.enable_oai_embedding:
        embedding_scorer = EmbeddingSimilarityScorer(
            openai_api_key=pipeline_settings.openai_api_key,
            model=PHENOTYPE_SIMILARITY_EMBEDDING_MODEL,
            phenotype_cache_path=pipeline_settings.phenotype_embedding_cache_path,
            target_cache_path=pipeline_settings.phenotype_target_embedding_cache_path,
        )
        if verbose and not embedding_scorer.is_available:
            print(
                "    Warning: OpenAI embedding similarity requested, but the embedding "
                "client is unavailable or OPENAI_API_KEY is not set."
            )

    return targets, embedding_scorer


def _extract_flybase_ids(value: Any) -> List[str]:
    """Extract FlyBase IDs from slash/space-delimited genotype fields."""
    matches = re.findall(r"FB[a-zA-Z]{2,4}[0-9]+", str(value or ""))
    deduped: List[str] = []
    seen: Set[str] = set()
    for match in matches:
        match_clean = clean_id(match)
        if not match_clean or match_clean in seen:
            continue
        seen.add(match_clean)
        deduped.append(match_clean)
    return deduped


def _normalize_symbol_token(value: Any) -> str:
    """Return a comparison-friendly normalization for genotype symbol fragments."""
    return re.sub(r"\s+", " ", str(value or "").strip()).lower()


def _split_genotype_symbol_tokens(value: Any) -> List[str]:
    """Split genotype text into readable component-like tokens.

    FlyBase genotype_symbols are space-delimited allele components
    (e.g. ``"A23[A23] Scer\\GAL4[Act5C.PI]"``).  Chromosomal arms are
    separated by ``";"``, homologous alleles by ``"/"``.
    """
    tokens: List[str] = []
    seen: Set[str] = set()
    raw_parts = re.split(r"[;|,\n]+", str(value or ""))
    for raw_part in raw_parts:
        for slash_part in str(raw_part).split("/"):
            for candidate in re.split(r"\s+", str(slash_part).strip()):
                token = candidate.strip()
                if not token or token in {"+", "-"}:
                    continue
                normalized = _normalize_symbol_token(token)
                if normalized in seen:
                    continue
                seen.add(normalized)
                tokens.append(token)
    return tokens


def _extract_genotype_id_symbol_pairs(
    genotype_fbids: Any,
    genotype_symbols: Any,
) -> List[Tuple[str, str]]:
    """Best-effort alignment between genotype FBids and symbol tokens."""
    ids = _extract_flybase_ids(genotype_fbids)
    if not ids:
        return []

    symbol_text = str(genotype_symbols or "").strip()
    if not symbol_text:
        return []

    candidate_lists: List[List[str]] = []

    slash_split = [
        token
        for token in (part.strip() for part in re.split(r"[/|;\n,]+", symbol_text))
        if token and token not in {"+", "-"}
    ]
    if slash_split:
        candidate_lists.append(slash_split)

    whitespace_split = [
        token
        for token in (part.strip() for part in re.split(r"\s+", symbol_text))
        if token and token not in {"+", "-"}
    ]
    if whitespace_split and whitespace_split != slash_split:
        candidate_lists.append(whitespace_split)

    if len(ids) == 1:
        return [(ids[0], symbol_text)]

    for tokens in candidate_lists:
        if len(tokens) == len(ids):
            return [
                (clean_id(component_id), token)
                for component_id, token in zip(ids, tokens)
                if clean_id(component_id) and token
            ]

    return []


def _looks_like_gal4_symbol(symbol: Any, gene_symbol: Any = "") -> bool:
    """Return True when the symbol or linked gene clearly indicates GAL4."""
    for raw_value in (symbol, gene_symbol):
        normalized = _normalize_symbol_token(raw_value)
        if not normalized:
            continue
        if "gal4" in normalized:
            return True
    return False


def _extract_unmatched_gal4_tokens(
    genotype_label: Any,
    matched_symbol_norms: Set[str],
) -> List[str]:
    """Return GAL4-containing genotype tokens that are not part of the focal stock."""
    gal4_tokens: List[str] = []
    seen: Set[str] = set()
    for token in _split_genotype_symbol_tokens(genotype_label):
        normalized_token = _normalize_symbol_token(token)
        if normalized_token in matched_symbol_norms or "gal4" not in normalized_token:
            continue
        if normalized_token in seen:
            continue
        seen.add(normalized_token)
        gal4_tokens.append(token)
    return gal4_tokens


def _format_stock_candidate_label(
    fbst: Any,
    collection: Any,
    stock_number: Any,
) -> str:
    """Format a stock candidate label for best-effort partner resolution."""
    stock_number_str = clean_id(stock_number)
    collection_str = str(collection or "").strip()
    if stock_number_str and collection_str:
        return f"({stock_number_str}, {collection_str})"
    if stock_number_str:
        return f"({stock_number_str})"
    if collection_str:
        return f"({collection_str})"
    return clean_id(fbst)


def _build_stock_phenotype_sheet(
    all_stocks_df: pd.DataFrame,
    flybase_data_path: Optional[Path],
    references_df: Optional[pd.DataFrame] = None,
    unfiltered_references_df: Optional[pd.DataFrame] = None,
    pubmed_cache_path: Optional[Path] = None,
    pubmed_client: Optional[PubMedClient] = None,
    similarity_targets: Optional[List[PhenotypeSimilarityTarget]] = None,
    embedding_scorer: Optional[EmbeddingSimilarityScorer] = None,
    verbose: bool = False,
    gene_to_datasets: Optional[Dict[str, Set[str]]] = None,
) -> pd.DataFrame:
    """Build a soft-run replacement sheet from genotype_phenotype_data."""
    if all_stocks_df is None or len(all_stocks_df) == 0 or flybase_data_path is None:
        return pd.DataFrame()

    included_df = all_stocks_df.copy()
    if 'FBst' not in included_df.columns:
        return pd.DataFrame()

    def _normalize_authors(value: Any) -> str:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return ''
        if isinstance(value, str):
            return '; '.join([a.strip() for a in value.split(';') if a.strip()])
        try:
            return '; '.join([str(a).strip() for a in list(value) if str(a).strip()])
        except Exception:
            return str(value).strip()

    def _empty_reference_meta() -> Dict[str, str]:
        return {
            'title': '',
            'authors': '',
            'journal': '',
            'publication_date': '',
            'miniref': '',
        }

    def _merge_reference_meta(base: Dict[str, str], incoming: Dict[str, Any]) -> Dict[str, str]:
        merged = dict(base)
        normalized = {
            'title': str(incoming.get('title', '') or '').strip(),
            'authors': _normalize_authors(incoming.get('authors', '')),
            'journal': str(incoming.get('journal', '') or '').strip(),
            'publication_date': str(
                incoming.get('publication_date', '')
                or incoming.get('publication date', '')
                or incoming.get('date of publication', '')
                or incoming.get('year', '')
                or ''
            ).strip(),
            'miniref': str(incoming.get('miniref', '') or '').strip(),
        }
        for key, value in normalized.items():
            if value and not merged.get(key, '').strip():
                merged[key] = value
        return merged

    def _parse_miniref(value: Any) -> Dict[str, str]:
        text = str(value or '').strip()
        if not text:
            return _empty_reference_meta()

        authors = ''
        publication_date = ''
        journal = ''
        parts = [part.strip() for part in text.split(',') if part.strip()]
        if parts:
            authors = parts[0]
            year_idx = None
            for idx, part in enumerate(parts[1:], start=1):
                year_match = re.search(r'\b(19|20)\d{2}\b', part)
                if year_match:
                    publication_date = year_match.group(0)
                    year_idx = idx
                    break

            if year_idx is not None:
                journal = ', '.join(parts[year_idx + 1:]).strip()
            elif len(parts) > 1:
                journal = ', '.join(parts[1:]).strip()

        return {
            'title': '',
            'authors': authors,
            'journal': journal,
            'publication_date': publication_date,
            'miniref': text,
        }

    def _needs_reference_enrichment(meta: Dict[str, str]) -> bool:
        return not all(
            str(meta.get(key, '') or '').strip()
            for key in ('title', 'authors', 'journal', 'publication_date')
        )

    if 'FBst' in included_df.columns:
        included_df['FBst'] = included_df['FBst'].apply(clean_id)
    else:
        included_df['FBst'] = ''
    included_df['stock_number'] = included_df.get('stock_number', '').apply(clean_id)
    included_df['collection'] = included_df.get('collection', '').fillna('').astype(str).str.strip()
    included_df['relevant_gene_symbols'] = included_df.get('relevant_gene_symbols', '').fillna('').astype(str)

    def _is_custom_row(row):
        v = row.get('custom_stock', False)
        return v is True or str(v).strip().lower() in ('true', '1', 'yes')

    custom_included_df = pd.DataFrame()
    if 'custom_stock' in included_df.columns:
        custom_mask = included_df.apply(_is_custom_row, axis=1)
        custom_included_df = included_df[custom_mask].copy()

    metadata_text_columns = (
        'relevant_gene_symbols',
        'Balancers',
        'matched_component_types',
        'allele_class_terms',
        'transgenic_product_class_terms',
    )
    for df in (included_df, custom_included_df):
        for column in metadata_text_columns:
            if column not in df.columns:
                df[column] = ''
            df[column] = df[column].fillna('').astype(str)
    included_df = _normalize_reagent_bucket_columns(included_df)
    custom_included_df = _normalize_reagent_bucket_columns(custom_included_df)
    _assert_one_hot_reagent_buckets(included_df, "Included stock rows")
    _assert_one_hot_reagent_buckets(custom_included_df, "Custom phenotype rows")

    included_df = included_df[included_df['FBst'].astype(bool)].copy()
    if len(included_df) == 0 and len(custom_included_df) == 0:
        return pd.DataFrame()

    flybase_data_path = Path(flybase_data_path)
    alleles_dir = flybase_data_path / 'alleles_and_stocks'
    derived_components_path = alleles_dir / 'fbst_to_derived_stock_component.csv'

    lookup_derived_df = pd.DataFrame()
    derived_df = pd.DataFrame()
    if not derived_components_path.exists():
        if verbose and len(included_df) > 0:
            print(f"    Warning: Derived stock component CSV not found: {derived_components_path}")
    else:
        lookup_derived_df = pd.read_csv(derived_components_path, dtype=str).fillna('')

    try:
        phenotype_path = find_latest_tsv(alleles_dir, 'genotype_phenotype_data')
        phenotype_df = load_flybase_tsv(phenotype_path)
    except FileNotFoundError:
        if verbose:
            print(f"    Warning: Could not find genotype_phenotype_data under {alleles_dir}")
        return pd.DataFrame()

    if len(phenotype_df) == 0:
        return pd.DataFrame()

    if len(lookup_derived_df) > 0:
        if 'collection_short_name' in lookup_derived_df.columns and 'collection' not in lookup_derived_df.columns:
            lookup_derived_df = lookup_derived_df.rename(columns={'collection_short_name': 'collection'})
        for col in (
            'FBst',
            'stock_number',
            'collection',
            'derived_stock_component',
            'object_symbol',
            'GeneSymbol',
        ):
            if col not in lookup_derived_df.columns:
                lookup_derived_df[col] = ''
        lookup_derived_df['FBst'] = lookup_derived_df['FBst'].apply(clean_id)
        lookup_derived_df['stock_number'] = lookup_derived_df['stock_number'].apply(clean_id)
        lookup_derived_df['derived_stock_component'] = lookup_derived_df['derived_stock_component'].apply(clean_id)
        lookup_derived_df['collection'] = lookup_derived_df['collection'].fillna('').astype(str).str.strip()

    if len(lookup_derived_df) > 0 and len(included_df) > 0:
        derived_df = lookup_derived_df[lookup_derived_df['FBst'].isin(set(included_df['FBst']))].copy()
        relevant_component_pairs: Set[Tuple[str, str]] = set()
        if 'relevant_component_ids' in included_df.columns:
            for _, included_row in included_df[['FBst', 'relevant_component_ids']].iterrows():
                fbst = clean_id(included_row.get('FBst', ''))
                if not fbst:
                    continue
                for component_id in parse_semicolon_list(included_row.get('relevant_component_ids', '')):
                    component_id = clean_id(component_id)
                    if component_id:
                        relevant_component_pairs.add((fbst, component_id))
        if relevant_component_pairs:
            derived_pairs = pd.Series(
                list(zip(derived_df['FBst'], derived_df['derived_stock_component'])),
                index=derived_df.index,
            )
            derived_df = derived_df[derived_pairs.isin(relevant_component_pairs)].copy()

    if len(derived_df) == 0 and len(custom_included_df) == 0:
        return pd.DataFrame()

    phenotype_df['reference'] = phenotype_df.get('reference', '').fillna('').astype(str)
    phenotype_df['phenotype_name'] = phenotype_df.get('phenotype_name', '').fillna('').astype(str)
    phenotype_df['phenotype_id'] = phenotype_df.get('phenotype_id', '').fillna('').astype(str)
    phenotype_df['qualifier_names'] = phenotype_df.get('qualifier_names', '').fillna('').astype(str)
    phenotype_df['genotype_FBids'] = phenotype_df.get('genotype_FBids', '').fillna('').astype(str)
    phenotype_df['genotype_symbols'] = phenotype_df.get('genotype_symbols', '').fillna('').astype(str)

    all_relevant_ids = set()
    if len(derived_df) > 0:
        all_relevant_ids.update(derived_df['derived_stock_component'].astype(str).str.strip())
    if 'relevant_component_ids' in included_df.columns:
        for ids_str in included_df['relevant_component_ids'].fillna('').astype(str):
            all_relevant_ids.update(cid for cid in parse_semicolon_list(ids_str) if cid)
    if len(custom_included_df) > 0 and 'relevant_component_ids' in custom_included_df.columns:
        for ids_str in custom_included_df['relevant_component_ids'].fillna('').astype(str):
            all_relevant_ids.update(cid for cid in parse_semicolon_list(ids_str) if cid)
    all_relevant_ids.discard('')

    def _row_has_relevant_component(fbids_str):
        return bool(set(_extract_flybase_ids(fbids_str)) & all_relevant_ids)

    relevant_mask = phenotype_df['genotype_FBids'].apply(_row_has_relevant_component)
    phenotype_df = phenotype_df[relevant_mask].copy()

    if len(phenotype_df) == 0:
        return pd.DataFrame()

    similarity_targets = similarity_targets or []
    cosine_similarity_columns = [
        target.cosine_similarity_column for target in similarity_targets
    ]
    similarity_columns: List[str] = list(cosine_similarity_columns)

    phenotype_df['phenotype_text'] = phenotype_df['phenotype_name'].apply(
        normalize_phenotype_text
    )
    phenotype_df['qualifier_text'] = phenotype_df['qualifier_names'].apply(
        normalize_qualifier_text
    )
    phenotype_df['phenotype_similarity_text'] = phenotype_df.apply(
        lambda row: normalize_phenotype_text(row.get('phenotype_name', '')),
        axis=1,
    )

    embedding_similarity_by_text: Dict[str, Dict[str, Optional[float]]] = {}
    if embedding_scorer is not None:
        embedding_similarity_by_text = embedding_scorer.score_texts(
            phenotype_df['phenotype_similarity_text'].tolist(),
            similarity_targets,
        )

    # ------------------------------------------------------------------
    # FBrf → title / URL resolution
    # ------------------------------------------------------------------
    refs_dir = flybase_data_path / 'references'
    fbrf_to_pmid: Dict[str, str] = {}
    fbrf_to_pmcid: Dict[str, str] = {}
    fbrf_to_doi: Dict[str, str] = {}
    fbrf_to_miniref: Dict[str, str] = {}
    reference_meta_by_pmid: Dict[str, Dict[str, str]] = {}
    pmid_to_resolved_doi: Dict[str, str] = {}
    try:
        refs_path = find_latest_tsv(refs_dir, 'fbrf_pmid_pmcid_doi')
        refs_tsv = load_flybase_tsv(refs_path)
        for c in ('FBrf', 'PMID', 'PMCID', 'DOI', 'miniref'):
            if c not in refs_tsv.columns:
                refs_tsv[c] = ''
        refs_tsv['FBrf'] = refs_tsv['FBrf'].apply(clean_id)
        refs_tsv['PMID'] = refs_tsv['PMID'].apply(clean_id)
        refs_tsv['PMCID'] = refs_tsv['PMCID'].fillna('').astype(str).str.strip()
        refs_tsv['DOI'] = refs_tsv['DOI'].fillna('').astype(str).str.strip()
        refs_tsv['miniref'] = refs_tsv['miniref'].fillna('').astype(str).str.strip()
        fbrf_to_pmid = dict(zip(refs_tsv['FBrf'], refs_tsv['PMID']))
        fbrf_to_pmcid = dict(zip(refs_tsv['FBrf'], refs_tsv['PMCID']))
        fbrf_to_doi = dict(zip(refs_tsv['FBrf'], refs_tsv['DOI']))
        fbrf_to_miniref = dict(zip(refs_tsv['FBrf'], refs_tsv['miniref']))
        pmid_to_resolved_doi.update({
            pmid: doi
            for pmid, doi in zip(refs_tsv['PMID'], refs_tsv['DOI'])
            if pmid and doi
        })
    except FileNotFoundError:
        if verbose:
            print(f"    Warning: Could not find fbrf_pmid_pmcid_doi under {refs_dir}")

    reference_sources: List[pd.DataFrame] = []
    if unfiltered_references_df is not None and len(unfiltered_references_df) > 0:
        reference_sources.append(unfiltered_references_df)
    elif references_df is not None and len(references_df) > 0:
        reference_sources.append(references_df)

    for reference_source_df in reference_sources:
        if 'PMID' not in reference_source_df.columns:
            continue
        for _, ref_row in reference_source_df.iterrows():
            pmid = clean_id(ref_row.get('PMID', ''))
            if not pmid:
                continue
            doi = str(
                ref_row.get('DOI', ref_row.get('doi', '')) or ''
            ).strip()
            reference_meta_by_pmid[pmid] = _merge_reference_meta(
                reference_meta_by_pmid.get(pmid, _empty_reference_meta()),
                {
                    'title': ref_row.get('title', ''),
                    'authors': ref_row.get('author(s)', ref_row.get('authors', '')),
                    'journal': ref_row.get('journal', ''),
                    'publication_date': ref_row.get(
                        'publication_date',
                        ref_row.get('publication date', ref_row.get('date of publication', '')),
                    ),
                },
            )
            if doi:
                pmid_to_resolved_doi[pmid] = doi

    cache_to_read = (
        pubmed_client.cache if pubmed_client is not None and pubmed_client.cache is not None
        else PubMedCache(pubmed_cache_path) if pubmed_cache_path
        else None
    )
    if cache_to_read is not None:
        try:
            cache_data = cache_to_read.load()
            for pmid, meta in cache_data.items():
                doi = str(meta.get('doi', '') or '').strip()
                reference_meta_by_pmid[pmid] = _merge_reference_meta(
                    reference_meta_by_pmid.get(pmid, _empty_reference_meta()),
                    {
                        'title': meta.get('title', ''),
                        'authors': meta.get('authors', ''),
                        'journal': meta.get('journal', ''),
                        'publication_date': meta.get('year', ''),
                    },
                )
                if doi:
                    pmid_to_resolved_doi[pmid] = doi
        except Exception as exc:
            print(f"    Warning: Could not load PubMed cache: {exc}")

    phenotype_fbrfs: Set[str] = set()
    for raw_reference_ids in phenotype_df['reference'].dropna().astype(str):
        for value in raw_reference_ids.split('|'):
            fbrf = clean_id(value)
            if fbrf:
                phenotype_fbrfs.add(fbrf)

    phenotype_pmids = sorted({
        pmid
        for fbrf in phenotype_fbrfs
        for pmid in [clean_id(fbrf_to_pmid.get(fbrf, ''))]
        if pmid
    })
    pmids_to_fetch = [
        pmid
        for pmid in phenotype_pmids
        if _needs_reference_enrichment(reference_meta_by_pmid.get(pmid, _empty_reference_meta()))
    ]
    if pmids_to_fetch and pubmed_client is not None:
        try:
            fetched_metadata = pubmed_client.fetch_metadata(pmids_to_fetch)
            for pmid, meta in fetched_metadata.items():
                doi = str(meta.get('doi', '') or '').strip()
                reference_meta_by_pmid[pmid] = _merge_reference_meta(
                    reference_meta_by_pmid.get(pmid, _empty_reference_meta()),
                    {
                        'title': meta.get('title', ''),
                        'authors': meta.get('authors', ''),
                        'journal': meta.get('journal', ''),
                        'publication_date': meta.get('year', ''),
                    },
                )
                if doi:
                    pmid_to_resolved_doi[pmid] = doi
        except Exception as exc:
            print(f"    Warning: Could not fetch PubMed metadata for phenotype references: {exc}")

    def _resolve_fbrf(fbrf: str) -> Dict[str, str]:
        """Return reference display + identifier metadata for an FBrf ID."""
        pmid = fbrf_to_pmid.get(fbrf, '')
        pmcid = fbrf_to_pmcid.get(fbrf, '')
        doi = fbrf_to_doi.get(fbrf, '')
        meta = _empty_reference_meta()
        if pmid:
            meta = dict(reference_meta_by_pmid.get(pmid, _empty_reference_meta()))
            doi = pmid_to_resolved_doi.get(pmid, '') or doi
        miniref = fbrf_to_miniref.get(fbrf, '')
        if miniref:
            meta = _merge_reference_meta(meta, _parse_miniref(miniref))
        if doi:
            url = f"https://doi.org/{doi}" if not doi.startswith('http') else doi
        elif pmid:
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        elif pmcid:
            normalized_pmcid = pmcid if pmcid.startswith('PMC') else f"PMC{pmcid}"
            url = f"https://pmc.ncbi.nlm.nih.gov/articles/{normalized_pmcid}/"
        else:
            url = f"https://flybase.org/reports/{fbrf}"
        display = meta.get('title', '') or meta.get('miniref', '') or doi or fbrf
        return {
            'display': display,
            'url': url,
            'pmid': pmid,
            'pmcid': pmcid,
            'authors': meta.get('authors', ''),
            'journal': meta.get('journal', ''),
            'publication_date': meta.get('publication_date', ''),
        }

    # ------------------------------------------------------------------

    stock_meta_by_fbst: Dict[str, Dict] = {}
    if len(included_df) > 0:
        stock_meta = (
            included_df.groupby('FBst', sort=False)
            .agg(
                stock_number=('stock_number', lambda s: unique_join(s.tolist())),
                collection=('collection', lambda s: unique_join(s.tolist())),
                gene_symbol=('relevant_gene_symbols', lambda s: unique_join([g for raw in s.tolist() for g in parse_semicolon_list(raw)])),
                Balancers=('Balancers', lambda s: unique_join(s.tolist())),
                matched_component_types=('matched_component_types', lambda s: unique_join(s.tolist())),
                allele_class_terms=('allele_class_terms', lambda s: unique_join(s.tolist())),
                transgenic_product_class_terms=(
                    'transgenic_product_class_terms',
                    lambda s: unique_join(s.tolist()),
                ),
                **{
                    column: (column, lambda s: any(_is_truthy(v) for v in s.tolist()))
                    for column in REAGENT_BUCKET_COLUMNS
                },
            )
            .reset_index()
        )
        stock_meta_by_fbst = stock_meta.set_index('FBst').to_dict('index')

    # Build metadata for custom_stock rows keyed by stock_number (allele symbol)
    custom_meta_by_stock_num: Dict[str, Dict] = {}
    if len(custom_included_df) > 0:
        for _, crow in custom_included_df.iterrows():
            sn = str(crow.get('stock_number', '') or '').strip()
            if not sn:
                continue
            custom_meta_by_stock_num[sn] = {
                'stock_number': sn,
                'collection': str(crow.get('collection', '') or '').strip(),
                'gene_symbol': unique_join(parse_semicolon_list(str(crow.get('relevant_gene_symbols', '')))),
                'Balancers': str(crow.get('Balancers', '') or '').strip(),
                'matched_component_types': str(
                    crow.get('matched_component_types', '') or ''
                ).strip(),
                'allele_class_terms': str(crow.get('allele_class_terms', '') or '').strip(),
                'transgenic_product_class_terms': str(
                    crow.get('transgenic_product_class_terms', '') or ''
                ).strip(),
                **{
                    column: _is_truthy(crow.get(column, False))
                    for column in REAGENT_BUCKET_COLUMNS
                },
            }

    # Build reverse index: component_id -> set of stock keys (FBst for
    # stock-backed rows, stock_number for custom rows), using both
    # derived_stock_component (Chado-level) and relevant_component_ids
    # (Stage 1 allele-level) so phenotype rows referencing FBal IDs can
    # still resolve to stocks whose Chado entries only have FBti/FBtp.
    component_to_stock_keys: Dict[str, Set[str]] = {}
    if len(derived_df) > 0:
        for _, row in derived_df.iterrows():
            cid = clean_id(row.get('derived_stock_component', ''))
            fbst = clean_id(row.get('FBst', ''))
            if cid and fbst:
                component_to_stock_keys.setdefault(cid, set()).add(fbst)
    if 'relevant_component_ids' in included_df.columns:
        for _, row in included_df.iterrows():
            fbst = clean_id(row.get('FBst', ''))
            if not fbst:
                continue
            for cid in parse_semicolon_list(str(row.get('relevant_component_ids', ''))):
                cid = clean_id(cid)
                if cid:
                    component_to_stock_keys.setdefault(cid, set()).add(fbst)
    if len(custom_included_df) > 0 and 'relevant_component_ids' in custom_included_df.columns:
        for _, crow in custom_included_df.iterrows():
            sn = str(crow.get('stock_number', '') or '').strip()
            if not sn:
                continue
            for cid in parse_semicolon_list(str(crow.get('relevant_component_ids', ''))):
                cid = clean_id(cid)
                if cid:
                    component_to_stock_keys.setdefault(cid, set()).add(sn)

    # Build component_id -> symbol lookup from derived_df and included_df
    component_id_to_symbol: Dict[str, str] = {}
    component_id_to_gene_symbol: Dict[str, str] = {}
    component_id_to_stock_candidates: Dict[str, Set[str]] = defaultdict(set)
    component_id_to_stock_candidate_details: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    gal4_only_fbst: Set[str] = set()
    if len(lookup_derived_df) > 0:
        component_rows_by_fbst: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
        for _, row in lookup_derived_df.iterrows():
            cid = clean_id(row.get('derived_stock_component', ''))
            sym = str(row.get('object_symbol', '') or '').strip()
            gene_symbol = str(row.get('GeneSymbol', '') or '').strip()
            fallback_symbol = str(row.get('FB_genotype', '') or '').strip()
            fbst = clean_id(row.get('FBst', ''))
            candidate_label = _format_stock_candidate_label(
                fbst,
                row.get('collection', ''),
                row.get('stock_number', ''),
            )
            if cid and sym:
                component_id_to_symbol.setdefault(cid, sym)
            if cid and gene_symbol:
                component_id_to_gene_symbol.setdefault(cid, gene_symbol)
            if cid and fbst and candidate_label:
                component_id_to_stock_candidates[cid].add(candidate_label)
                component_id_to_stock_candidate_details[cid].add((fbst, candidate_label))
            if fbst:
                component_rows_by_fbst[fbst].append((sym or fallback_symbol, gene_symbol))
        gal4_only_fbst = {
            fbst
            for fbst, component_rows in component_rows_by_fbst.items()
            if component_rows
            and all(
                _looks_like_gal4_symbol(symbol, gene_symbol)
                for symbol, gene_symbol in component_rows
            )
        }
    for src_df in [included_df, custom_included_df]:
        if len(src_df) > 0 and 'relevant_fbal_ids' in src_df.columns and 'relevant_fbal_symbols' in src_df.columns:
            for _, row in src_df.iterrows():
                for cid, sym in zip(
                    parse_semicolon_list(str(row.get('relevant_fbal_ids', ''))),
                    parse_semicolon_list(str(row.get('relevant_fbal_symbols', ''))),
                ):
                    cid = clean_id(cid)
                    sym = str(sym).strip()
                    if cid and sym and cid not in component_id_to_symbol:
                        component_id_to_symbol[cid] = sym
                    gene_values = parse_semicolon_list(str(row.get('relevant_gene_symbols', '') or ''))
                    if cid and gene_values and cid not in component_id_to_gene_symbol:
                        component_id_to_gene_symbol[cid] = unique_join(gene_values)

    def _filter_gal4_only_candidate_labels(
        candidate_details: Set[Tuple[str, str]],
    ) -> Set[str]:
        return {
            label
            for fbst, label in candidate_details
            if fbst in gal4_only_fbst and label
        }

    # ----------------------------------------------------------------
    # FBal → FBtp → FBti chain for resolving partner GAL4 alleles to
    # stock candidates.  Phenotype genotype data references alleles
    # (FBal), but the derived stock CSV indexes GAL4 stocks via their
    # FBti insertion IDs.  We bridge the gap by loading construct
    # descriptions (FBal → FBtp) and the insertion map (FBtp → FBti),
    # then forwarding stock candidates from any matching FBti entries.
    # ----------------------------------------------------------------
    constructs_dir = flybase_data_path / 'transgenic_constructs'
    insertions_dir = flybase_data_path / 'transgenic_insertions'
    fbal_to_fbtps: Dict[str, Set[str]] = {}
    fbtp_to_fbtis: Dict[str, Set[str]] = {}
    try:
        construct_path = find_latest_tsv(constructs_dir, 'transgenic_construct_descriptions')
        construct_df = load_flybase_tsv(construct_path)
        fbal_col = 'Component Allele (id)' if 'Component Allele (id)' in construct_df.columns else None
        fbtp_col = 'Transgenic Construct (id)' if 'Transgenic Construct (id)' in construct_df.columns else None
        if fbal_col is None or fbtp_col is None:
            raw_cols = list(construct_df.columns)
            if len(raw_cols) >= 4:
                fbal_col, fbtp_col = raw_cols[1], raw_cols[3]
        if fbal_col and fbtp_col:
            for _, c_row in construct_df[[fbal_col, fbtp_col]].iterrows():
                fbal_id = clean_id(str(c_row[fbal_col] or ''))
                if not fbal_id:
                    continue
                for fbtp_raw in str(c_row[fbtp_col] or '').split('|'):
                    fbtp_id = clean_id(fbtp_raw)
                    if fbtp_id:
                        fbal_to_fbtps.setdefault(fbal_id, set()).add(fbtp_id)
    except FileNotFoundError:
        pass
    fbtp_to_fbti_path = insertions_dir / 'fbtp_to_fbti.csv'
    if fbtp_to_fbti_path.exists():
        fbtp_fbti_df = pd.read_csv(fbtp_to_fbti_path, dtype=str).fillna('')
        for _, m_row in fbtp_fbti_df.iterrows():
            fbtp_id = clean_id(str(m_row.get('FBtp', '') or ''))
            fbti_id = clean_id(str(m_row.get('FBti', '') or ''))
            if fbtp_id and fbti_id:
                fbtp_to_fbtis.setdefault(fbtp_id, set()).add(fbti_id)

    def _resolve_fbal_stock_candidates(
        fbal_id: str,
        gal4_only: bool = False,
    ) -> Set[str]:
        """Chain FBal → FBtp → FBti → stock candidates."""
        candidates: Set[str] = set()
        for fbtp_id in fbal_to_fbtps.get(fbal_id, set()):
            for fbti_id in fbtp_to_fbtis.get(fbtp_id, set()):
                if gal4_only:
                    candidates.update(
                        _filter_gal4_only_candidate_labels(
                            component_id_to_stock_candidate_details.get(fbti_id, set())
                        )
                    )
                else:
                    candidates.update(component_id_to_stock_candidates.get(fbti_id, set()))
        return candidates

    phenotype_rows: List[Dict[str, str]] = []
    for _, pheno_row in phenotype_df.iterrows():
        component_ids = _extract_flybase_ids(pheno_row.get('genotype_FBids', ''))
        if not component_ids:
            continue

        matched_keys: Dict[str, Set[str]] = {}
        for cid in component_ids:
            for key in component_to_stock_keys.get(cid, set()):
                matched_keys.setdefault(key, set()).add(cid)
        if not matched_keys:
            continue

        phenotype_name = str(
            pheno_row.get('phenotype_text', '')
            or normalize_phenotype_text(pheno_row.get('phenotype_name', ''))
        ).strip()
        qualifier_text = str(
            pheno_row.get('qualifier_text', '')
            or normalize_qualifier_text(pheno_row.get('qualifier_names', ''))
        ).strip()
        genotype_label = str(pheno_row.get('genotype_symbols', '') or '').strip()
        genotype_id_symbol_pairs = _extract_genotype_id_symbol_pairs(
            pheno_row.get('genotype_FBids', ''),
            genotype_label,
        )
        genotype_symbols_by_id: Dict[str, List[str]] = defaultdict(list)
        for component_id, symbol in genotype_id_symbol_pairs:
            if component_id and symbol:
                genotype_symbols_by_id[component_id].append(symbol)
        raw_fbrfs = [clean_id(v) for v in str(pheno_row.get('reference', '') or '').split('|') if clean_id(v)]
        if not phenotype_name and not raw_fbrfs:
            continue

        cosine_scores = embedding_similarity_by_text.get(phenotype_name, {})

        for stock_key, matched_cids in matched_keys.items():
            stock_info = stock_meta_by_fbst.get(stock_key) or custom_meta_by_stock_num.get(stock_key)
            if not stock_info:
                continue
            matched_component_symbols = [
                unique_join(
                    genotype_symbols_by_id.get(cid, []) or [component_id_to_symbol.get(cid, cid)]
                )
                for cid in matched_cids
            ]
            component_symbols = unique_join(matched_component_symbols)
            matched_symbol_norms = {
                _normalize_symbol_token(symbol)
                for symbol in matched_component_symbols
                if str(symbol).strip()
            }
            co_reagent_ids = [cid for cid in component_ids if cid not in matched_cids]
            unmatched_gal4_tokens = _extract_unmatched_gal4_tokens(
                genotype_label,
                matched_symbol_norms,
            ) if co_reagent_ids else []
            partner_gal4_ids: List[str] = []
            partner_gal4_symbol_values: List[str] = []
            partner_symbol_values_by_id: Dict[str, List[str]] = defaultdict(list)
            for cid in co_reagent_ids:
                paired_symbol = unique_join(genotype_symbols_by_id.get(cid, []))
                symbol = paired_symbol or component_id_to_symbol.get(cid, '')
                gene_symbol = component_id_to_gene_symbol.get(cid, '')
                if _looks_like_gal4_symbol(symbol, gene_symbol):
                    partner_gal4_ids.append(cid)
                    if symbol:
                        partner_gal4_symbol_values.append(symbol)
                        partner_symbol_values_by_id[cid].append(symbol)
                    elif gene_symbol:
                        partner_gal4_symbol_values.append(gene_symbol)
                        partner_symbol_values_by_id[cid].append(gene_symbol)
                    else:
                        partner_gal4_symbol_values.append(cid)
                        partner_symbol_values_by_id[cid].append(cid)
            if unmatched_gal4_tokens and not partner_gal4_ids and len(co_reagent_ids) == 1:
                cid = co_reagent_ids[0]
                partner_gal4_ids.append(cid)
                partner_gal4_symbol_values.extend(unmatched_gal4_tokens)
                partner_symbol_values_by_id[cid].extend(unmatched_gal4_tokens)
            co_reagent_fbids = unique_join(partner_gal4_ids)
            co_reagent_symbols = unique_join(partner_gal4_symbol_values)
            partner_symbol_norms = {
                _normalize_symbol_token(s) for s in partner_gal4_symbol_values if s
            }
            extra_gal4_tokens = [
                t for t in unmatched_gal4_tokens
                if _normalize_symbol_token(t) not in partner_symbol_norms
            ]
            partner_stock_candidate_set: Set[str] = set()
            filtered_partner_driver_values: List[str] = []
            for cid in partner_gal4_ids:
                direct = _filter_gal4_only_candidate_labels(
                    component_id_to_stock_candidate_details.get(cid, set())
                )
                if not direct:
                    direct = _resolve_fbal_stock_candidates(cid, gal4_only=True)
                if not direct:
                    continue
                partner_stock_candidate_set.update(direct)
                filtered_partner_driver_values.extend(
                    partner_symbol_values_by_id.get(cid, [cid])
                )
            filtered_partner_driver_values.extend(
                token
                for token in extra_gal4_tokens
                if filtered_partner_driver_values
                and _normalize_symbol_token(token)
                not in {
                    _normalize_symbol_token(value)
                    for value in filtered_partner_driver_values
                    if value
                }
            )
            partner_driver_symbols = unique_join(filtered_partner_driver_values)
            partner_driver_stock_candidates = unique_join(
                sorted(partner_stock_candidate_set)
            )
            component_gene_symbols = ''
            if len(derived_df) > 0:
                matched_derived = derived_df[
                    (derived_df['FBst'] == stock_key)
                    & derived_df['derived_stock_component'].isin(matched_cids)
                ]
                if len(matched_derived) > 0:
                    component_gene_symbols = unique_join(
                        matched_derived.get('GeneSymbol', pd.Series(dtype=str)).tolist()
                    )
            if not component_gene_symbols:
                component_gene_symbols = str(stock_info.get('gene_symbol', '') or '').strip()
            dataset_label = _get_dataset_label_for_gene_value(
                unique_join([
                    component_gene_symbols,
                    str(stock_info.get('gene_symbol', '') or '').strip(),
                ]),
                gene_to_datasets,
            )

            refs_to_emit = raw_fbrfs if raw_fbrfs else ['']
            for fbrf in refs_to_emit:
                ref_details = _resolve_fbrf(fbrf) if fbrf else {
                    'display': '',
                    'url': '',
                    'pmid': '',
                    'pmcid': '',
                    'authors': '',
                    'journal': '',
                    'publication_date': '',
                }
                stock_num = str(stock_info.get('stock_number', '') or '').strip()
                collection = str(stock_info.get('collection', '') or '').strip()
                source_stock = _format_source_stock_label(collection, stock_num)
                phenotype_row = {
                    'FBst': stock_key if stock_key.startswith('FBst') else '',
                    'Gene': component_gene_symbols,
                    'Reagent Type or Allele Symbol': component_symbols,
                    'Balancers': str(stock_info.get('Balancers', '') or '').strip(),
                    'matched_component_types': str(
                        stock_info.get('matched_component_types', '') or ''
                    ).strip(),
                    **{
                        column: _is_truthy(stock_info.get(column, False))
                        for column in REAGENT_BUCKET_COLUMNS
                    },
                    'allele_class_terms': str(stock_info.get('allele_class_terms', '') or '').strip(),
                    'transgenic_product_class_terms': str(
                        stock_info.get('transgenic_product_class_terms', '') or ''
                    ).strip(),
                    'Source/ Stock #': source_stock,
                    'Genotype': genotype_label,
                    CO_REAGENT_FBIDS_COLUMN: co_reagent_fbids,
                    CO_REAGENT_SYMBOLS_COLUMN: co_reagent_symbols,
                    PARTNER_DRIVER_SYMBOLS_COLUMN: partner_driver_symbols,
                    PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN: partner_driver_stock_candidates,
                    'Data Set': dataset_label,
                    'Phenotype': phenotype_name,
                    'Qualifier': qualifier_text,
                    'PMID': ref_details['pmid'],
                    'PMCID': ref_details['pmcid'],
                    'Reference': ref_details['display'],
                    'Authors': ref_details['authors'],
                    'Journal': ref_details['journal'],
                    'Year of Publication': ref_details['publication_date'],
                    '_reference_url': ref_details['url'],
                }
                for target in similarity_targets:
                    phenotype_row[target.cosine_similarity_column] = cosine_scores.get(
                        target.cosine_similarity_column
                    )
                phenotype_rows.append(phenotype_row)

    if not phenotype_rows:
        return pd.DataFrame()

    phenotype_sheet = pd.DataFrame(phenotype_rows)
    phenotype_sheet = (
        phenotype_sheet.groupby(
            [
                'Source/ Stock #',
                'Genotype',
                'Phenotype',
                'Qualifier',
                'PMID',
                'PMCID',
                'Reference',
                'Authors',
                'Journal',
                'Year of Publication',
                '_reference_url',
            ],
            sort=False,
            as_index=False,
        )
        .agg({
            'Gene': lambda s: unique_join(s.tolist()),
            'Reagent Type or Allele Symbol': lambda s: unique_join(s.tolist()),
            'Balancers': lambda s: unique_join(s.tolist()),
            'matched_component_types': lambda s: unique_join(s.tolist()),
            **{
                column: (lambda s: any(_is_truthy(v) for v in s.tolist()))
                for column in REAGENT_BUCKET_COLUMNS
            },
            'allele_class_terms': lambda s: unique_join(s.tolist()),
            'transgenic_product_class_terms': lambda s: unique_join(s.tolist()),
            CO_REAGENT_FBIDS_COLUMN: lambda s: unique_join(s.tolist()),
            CO_REAGENT_SYMBOLS_COLUMN: lambda s: unique_join(s.tolist()),
            PARTNER_DRIVER_SYMBOLS_COLUMN: lambda s: unique_join(s.tolist()),
            PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN: lambda s: unique_join(s.tolist()),
            'Data Set': lambda s: unique_join(s.tolist()),
            **{
                col: 'max'
                for col in similarity_columns
                if col in phenotype_sheet.columns
            },
        })
    )
    # Order rows by:
    # 1) Gene blocks (genes with more unique reagents first),
    # 2) Reagent blocks within each gene (reagents with more unique phenotypes first),
    # while keeping rows for the same reagent adjacent.
    def _split_sheet_values(value: Any) -> List[str]:
        return [v for v in parse_semicolon_list(str(value or '')) if v and v != '-']

    def _phenotype_key(row: pd.Series) -> str:
        phenotype = str(row.get('Phenotype', '') or '').strip()
        qualifier = str(row.get('Qualifier', '') or '').strip()
        return f"{phenotype} ({qualifier})" if phenotype and qualifier else phenotype or qualifier

    gene_to_reagents: Dict[str, Set[str]] = defaultdict(set)
    gene_reagent_to_phenotypes: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for _, row in phenotype_sheet.iterrows():
        genes = _split_sheet_values(row.get('Gene', ''))
        if not genes:
            genes = ['']
        reagent = str(row.get('Source/ Stock #', '') or '').strip()
        phenotypes = [_phenotype_key(row)]
        for gene in genes:
            if reagent:
                gene_to_reagents[gene].add(reagent)
            gene_reagent_to_phenotypes[(gene, reagent)].update(phenotypes)

    gene_sort_order = sorted(
        gene_to_reagents.keys(),
        key=lambda gene: (
            -len(gene_to_reagents.get(gene, set())),
            str(gene).lower(),
            str(gene),
        ),
    )
    gene_rank = {gene: idx for idx, gene in enumerate(gene_sort_order)}
    gene_reagent_count = {
        gene: len(reagents)
        for gene, reagents in gene_to_reagents.items()
    }

    reagent_rank_by_gene: Dict[str, Dict[str, int]] = {}
    for gene in gene_sort_order:
        ordered_reagents = sorted(
            gene_to_reagents.get(gene, set()),
            key=lambda reagent: (
                -len(gene_reagent_to_phenotypes.get((gene, reagent), set())),
                reagent.lower(),
                reagent,
            ),
        )
        reagent_rank_by_gene[gene] = {
            reagent: idx for idx, reagent in enumerate(ordered_reagents)
        }

    phenotype_sheet['_gene_values'] = phenotype_sheet['Gene'].apply(_split_sheet_values)
    phenotype_sheet['_primary_gene'] = phenotype_sheet['_gene_values'].apply(
        lambda genes: min(
            genes if genes else [''],
            key=lambda gene: (
                gene_rank.get(gene, float('inf')),
                str(gene).lower(),
                str(gene),
            ),
        )
    )
    phenotype_sheet['_primary_reagent'] = (
        phenotype_sheet['Source/ Stock #'].fillna('').astype(str).str.strip()
    )
    phenotype_sheet['_gene_reagent_count'] = phenotype_sheet['_primary_gene'].map(
        lambda gene: gene_reagent_count.get(gene, 0)
    )
    phenotype_sheet['_reagent_rank_within_gene'] = phenotype_sheet.apply(
        lambda row: reagent_rank_by_gene.get(row['_primary_gene'], {}).get(
            row['_primary_reagent'],
            float('inf'),
        ),
        axis=1,
    )

    phenotype_sheet = phenotype_sheet.sort_values(
        by=[
            '_gene_reagent_count',
            '_primary_gene',
            '_reagent_rank_within_gene',
            '_primary_reagent',
            'Genotype',
            'Phenotype',
            'PMID',
            'PMCID',
            'Reference',
            'Authors',
            'Journal',
            'Year of Publication',
        ],
        ascending=[False, True, True, True, True, True, True, True, True, True, True, True],
    ).reset_index(drop=True)
    _assert_one_hot_reagent_buckets(
        phenotype_sheet,
        "Stock Phenotype Sheet rows",
    )
    visible_columns = _get_stock_phenotype_sheet_output_columns(
        [col for col in similarity_columns if col in phenotype_sheet.columns]
    )
    return phenotype_sheet[
        visible_columns + ['_reference_url']
    ]


def _summarize_stock_phenotype_sheet(phenotype_sheet_df: pd.DataFrame) -> Dict[str, int]:
    """Return unique gene, stock, and reference counts for the phenotype sheet."""
    if phenotype_sheet_df is None or len(phenotype_sheet_df) == 0:
        return {
            'genes': 0,
            'stocks': 0,
            'references': 0,
        }

    unique_genes: Set[str] = set()
    if 'Gene' in phenotype_sheet_df.columns:
        for raw in phenotype_sheet_df['Gene'].dropna().astype(str):
            for gene in parse_semicolon_list(raw):
                gene = gene.strip()
                if gene and gene != '-':
                    unique_genes.add(gene)

    unique_stocks: Set[str] = set()
    source_stocks = _get_source_stock_series(phenotype_sheet_df)
    for source_stock in source_stocks:
        source_stock = source_stock.strip()
        if source_stock and source_stock != '-':
            unique_stocks.add(source_stock)

    unique_references: Set[Tuple[str, str]] = set()
    reference_labels = (
        phenotype_sheet_df['Reference'].fillna('').astype(str)
        if 'Reference' in phenotype_sheet_df.columns
        else pd.Series(dtype=str)
    )
    reference_urls = (
        phenotype_sheet_df['_reference_url'].fillna('').astype(str)
        if '_reference_url' in phenotype_sheet_df.columns
        else pd.Series(dtype=str)
    )
    for reference_label, reference_url in zip(reference_labels, reference_urls):
        reference_label = reference_label.strip()
        reference_url = reference_url.strip()
        if reference_label or reference_url:
            unique_references.add((reference_label, reference_url))

    return {
        'genes': len(unique_genes),
        'stocks': len(unique_stocks),
        'references': len(unique_references),
    }


def _get_cosine_similarity_columns(phenotype_sheet_df: pd.DataFrame) -> List[str]:
    """Return the phenotype-sheet cosine similarity columns."""
    if phenotype_sheet_df is None or len(phenotype_sheet_df.columns) == 0:
        return []
    return [
        str(col)
        for col in phenotype_sheet_df.columns
        if str(col).startswith("Cosine Similarity (")
    ]


def _compute_max_cosine_similarity(phenotype_sheet_df: pd.DataFrame) -> pd.Series:
    """Return the per-row maximum cosine similarity across all target columns."""
    if phenotype_sheet_df is None or len(phenotype_sheet_df) == 0:
        return pd.Series(dtype=float)

    cosine_columns = _get_cosine_similarity_columns(phenotype_sheet_df)
    if not cosine_columns:
        return pd.Series(index=phenotype_sheet_df.index, dtype=float)

    cosine_df = phenotype_sheet_df[cosine_columns].apply(pd.to_numeric, errors="coerce")
    return cosine_df.max(axis=1, skipna=True)


def _format_input_dataset_name(csv_path: Path) -> str:
    """Return the workbook-facing dataset label for an input CSV."""
    return csv_path.stem.strip()


def _get_dataset_label_for_gene_value(
    gene_value: Any,
    gene_to_datasets: Optional[Dict[str, Set[str]]],
) -> str:
    """Return a semicolon-joined dataset label for a phenotype-sheet gene cell."""
    if not gene_to_datasets:
        return ""
    dataset_names: Set[str] = set()
    for gene in parse_semicolon_list(str(gene_value or "")):
        gene = gene.strip()
        if not gene or gene == "-":
            continue
        dataset_names.update(gene_to_datasets.get(gene, set()))
    return unique_join(sorted(dataset_names))


def _get_masterlist_embedding_score_series(sheet_df: pd.DataFrame) -> pd.Series:
    """Return the template's embedding-score column from the phenotype sheet."""
    if sheet_df is None:
        return pd.Series(dtype=float)
    index = sheet_df.index
    if "Max Cosine Similarity" in sheet_df.columns:
        return pd.to_numeric(
            sheet_df["Max Cosine Similarity"],
            errors="coerce",
        ).round(6)
    return _compute_max_cosine_similarity(sheet_df).reindex(index).round(6)


def _reorder_to_masterlist_columns(sheet_df: pd.DataFrame) -> pd.DataFrame:
    """Return a template-ordered phenotype sheet with extra columns appended."""
    if sheet_df is None:
        return pd.DataFrame()
    sheet_out = _export_source_stock_columns(sheet_df)
    sheet_out = sheet_out.drop(columns=["_reference_url"], errors="ignore")
    if len(sheet_out.columns) == 0:
        return sheet_out

    used_source_columns: Set[str] = set()
    formatted_columns: Dict[str, pd.Series] = {}
    blank_series = pd.Series("", index=sheet_out.index, dtype=object)
    for template_column in MASTERLIST_TEMPLATE_COLUMNS:
        source_column = MASTERLIST_TEMPLATE_SOURCE_COLUMNS.get(template_column)
        if template_column == "Circadian/Sleep Relevance (embedding max score)":
            formatted_columns[template_column] = _get_masterlist_embedding_score_series(
                sheet_out
            )
            used_source_columns.add("Max Cosine Similarity")
            continue
        if (
            source_column
            and source_column != template_column
            and source_column in sheet_out.columns
        ):
            formatted_columns[template_column] = sheet_out[source_column]
            used_source_columns.add(source_column)
            continue
        if template_column in sheet_out.columns:
            formatted_columns[template_column] = sheet_out[template_column]
            used_source_columns.add(template_column)
            continue
        if source_column and source_column in sheet_out.columns:
            formatted_columns[template_column] = sheet_out[source_column]
            used_source_columns.add(source_column)
            continue
        formatted_columns[template_column] = blank_series.copy()

    masterlist_df = pd.DataFrame(formatted_columns, index=sheet_out.index)
    extra_columns = [
        column
        for column in sheet_out.columns
        if column not in used_source_columns and column not in MASTERLIST_TEMPLATE_COLUMNS
    ]
    return masterlist_df.join(sheet_out[extra_columns]) if extra_columns else masterlist_df


def _split_phenotype_sheet_values(value: Any) -> List[str]:
    """Split semicolon-delimited phenotype sheet cells into non-empty values."""
    return [v for v in parse_semicolon_list(str(value or '')) if v and v != '-']


def _sort_similarity_tier_rows(tier_df: pd.DataFrame) -> pd.DataFrame:
    """Group tier rows by gene, then sort reagents by tier-local max cosine."""
    if tier_df is None or len(tier_df) == 0:
        return pd.DataFrame(columns=tier_df.columns if tier_df is not None else None)

    sorted_df = tier_df.copy().reset_index(drop=True)
    sorted_df["_tier_original_order"] = range(len(sorted_df))
    sorted_df["_gene_values"] = (
        sorted_df.get("Gene", pd.Series("", index=sorted_df.index))
        .apply(_split_phenotype_sheet_values)
        .apply(lambda genes: genes if genes else [""])
    )
    sorted_df["_primary_reagent"] = _get_source_stock_series(sorted_df)
    sorted_df["_max_cosine_similarity"] = pd.to_numeric(
        sorted_df.get("Max Cosine Similarity", pd.Series(index=sorted_df.index, dtype=float)),
        errors="coerce",
    )

    gene_to_reagents: Dict[str, Set[str]] = defaultdict(set)
    reagent_score_by_gene: Dict[Tuple[str, str], float] = {}
    gene_values_list = sorted_df["_gene_values"].tolist()
    reagent_list = sorted_df["_primary_reagent"].tolist()
    score_list = sorted_df["_max_cosine_similarity"].tolist()

    for genes, reagent, score in zip(gene_values_list, reagent_list, score_list):
        for gene in genes:
            gene_to_reagents.setdefault(gene, set())
            if reagent:
                gene_to_reagents[gene].add(reagent)
                current_score = reagent_score_by_gene.get((gene, reagent), float("-inf"))
                if pd.notna(score):
                    reagent_score_by_gene[(gene, reagent)] = max(current_score, float(score))
                else:
                    reagent_score_by_gene.setdefault((gene, reagent), current_score)

    gene_sort_order = sorted(
        gene_to_reagents.keys(),
        key=lambda gene: (
            -len(gene_to_reagents.get(gene, set())),
            str(gene).lower(),
            str(gene),
        ),
    )
    gene_rank = {gene: idx for idx, gene in enumerate(gene_sort_order)}

    sorted_df["_primary_gene"] = sorted_df["_gene_values"].apply(
        lambda genes: min(
            genes if genes else [""],
            key=lambda gene: (
                gene_rank.get(gene, float("inf")),
                str(gene).lower(),
                str(gene),
            ),
        )
    )
    sorted_df["_gene_reagent_count"] = sorted_df["_primary_gene"].map(
        lambda gene: len(gene_to_reagents.get(gene, set()))
    )
    sorted_df["_reagent_max_cosine_within_gene"] = sorted_df.apply(
        lambda row: reagent_score_by_gene.get(
            (row["_primary_gene"], row["_primary_reagent"]),
            float("-inf"),
        ),
        axis=1,
    )

    sorted_df = sorted_df.sort_values(
        by=[
            "_gene_reagent_count",
            "_primary_gene",
            "_reagent_max_cosine_within_gene",
            "_primary_reagent",
            "_tier_original_order",
        ],
        ascending=[False, True, False, True, True],
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)

    helper_columns = {
        "_tier_original_order",
        "_gene_values",
        "_primary_reagent",
        "_max_cosine_similarity",
        "_primary_gene",
        "_gene_reagent_count",
        "_reagent_max_cosine_within_gene",
    }
    return sorted_df.drop(columns=[col for col in helper_columns if col in sorted_df.columns])


def _format_similarity_score_for_sheet_name(value: float) -> str:
    """Format a similarity bound compactly for Excel sheet names."""
    return f"{float(value):.3f}".rstrip("0").rstrip(".")


def _format_similarity_tier_sheet_name(
    lower_bound: Optional[float],
    upper_bound: Optional[float],
) -> str:
    """Build a compact, Excel-safe tier sheet name."""
    if lower_bound is None:
        return f"<{_format_similarity_score_for_sheet_name(float(upper_bound))}"[:31]
    if upper_bound is None:
        return _format_similarity_score_for_sheet_name(float(lower_bound))[:31]
    lower_label = _format_similarity_score_for_sheet_name(lower_bound)
    upper_label = _format_similarity_score_for_sheet_name(upper_bound)
    return f"{lower_label}-{upper_label}"[:31]


def _build_similarity_threshold_bins() -> List[Tuple[Optional[float], float]]:
    """Return descending similarity bins using the configured bin width."""
    bins: List[Tuple[Optional[float], float]] = []
    upper_bound = 1.0
    step = float(SIMILARITY_TIER_BIN_WIDTH)
    for _ in range(max(SIMILARITY_TIER_SHEET_COUNT - 1, 0)):
        lower_bound = round(upper_bound - step, 10)
        bins.append((lower_bound, upper_bound))
        upper_bound = lower_bound
    bins.append((None, upper_bound))
    return bins


def _build_similarity_tier_sheets(
    phenotype_sheet_df: pd.DataFrame,
) -> List[Tuple[str, pd.DataFrame, Dict[str, Any]]]:
    """Partition scored phenotype rows into fixed max-cosine bins.

    Tier sheets preserve the original Stock Phenotype Sheet row order and skip
    empty bins entirely so the workbook only contains reagent-bearing tabs.
    """
    if phenotype_sheet_df is None or len(phenotype_sheet_df) == 0:
        return []

    max_cosine_similarity = _compute_max_cosine_similarity(phenotype_sheet_df)
    scored_mask = max_cosine_similarity.notna()
    if not scored_mask.any():
        return []

    scored_df = phenotype_sheet_df.loc[scored_mask].copy()
    scored_df.insert(
        min(len(scored_df.columns), 8),
        "Max Cosine Similarity",
        max_cosine_similarity.loc[scored_df.index].round(6),
    )

    tiers: List[Tuple[str, pd.DataFrame, Dict[str, Any]]] = []
    threshold_bins = _build_similarity_threshold_bins()
    scores = pd.to_numeric(scored_df["Max Cosine Similarity"], errors="coerce")

    for lower_bound, upper_bound in threshold_bins:
        if lower_bound is None:
            tier_mask = scores < upper_bound
            range_label = f"<{_format_similarity_score_for_sheet_name(upper_bound)}"
        elif upper_bound >= 1.0:
            tier_mask = scores >= lower_bound
            range_label = (
                f"{_format_similarity_score_for_sheet_name(lower_bound)}-"
                f"{_format_similarity_score_for_sheet_name(upper_bound)}"
            )
        else:
            tier_mask = (scores >= lower_bound) & (scores < upper_bound)
            range_label = (
                f"{_format_similarity_score_for_sheet_name(lower_bound)}-"
                f"{_format_similarity_score_for_sheet_name(upper_bound)}"
            )

        tier_df = _sort_similarity_tier_rows(scored_df.loc[tier_mask].copy())
        if tier_df.empty:
            continue

        tier_metadata = {
            "lower_bound": lower_bound,
            "upper_bound": upper_bound,
            "row_count": len(tier_df),
            "range_label": range_label,
        }
        tier_df.insert(
            0,
            "Similarity Range",
            range_label,
        )
        tiers.append(
            (
                _format_similarity_tier_sheet_name(
                    lower_bound=lower_bound,
                    upper_bound=upper_bound,
                ),
                tier_df,
                tier_metadata,
            )
        )

    return tiers


def _write_reference_hyperlinks(
    worksheet,
    workbook,
    sheet_df: pd.DataFrame,
    reference_urls: Optional[pd.Series],
) -> None:
    """Apply hyperlink formatting to the Reference column when URLs are available."""
    if reference_urls is None or 'Reference' not in sheet_df.columns:
        return

    ref_col_idx = sheet_df.columns.get_loc('Reference')
    url_fmt = workbook.add_format({'font_color': 'blue', 'underline': 1})
    for excel_row, url in enumerate(reference_urls, start=1):
        url = str(url or '').strip()
        if not url:
            continue
        display = str(sheet_df.iloc[excel_row - 1]['Reference'] or '').strip()
        if display:
            worksheet.write_url(excel_row, ref_col_idx, url, url_fmt, display)


def _write_phenotype_similarity_sheet(
    writer: pd.ExcelWriter,
    _workbook,
    sheet_name: str,
    phenotype_sheet_df: pd.DataFrame,
    format_as_masterlist: bool = False,
) -> None:
    """Write a phenotype-derived sheet without emitting workbook hyperlinks."""
    sheet_out = phenotype_sheet_df.copy()
    if format_as_masterlist:
        sheet_out = _reorder_to_masterlist_columns(sheet_out)
    else:
        sheet_out = _export_source_stock_columns(sheet_out)
        sheet_out = sheet_out.drop(columns=['_reference_url'], errors='ignore')
    sheet_out.to_excel(writer, sheet_name=sheet_name, index=False)


def _format_source_stock_label(collection: Any, stock_number: Any) -> str:
    """Format the stock label used on the Stock Phenotype sheet."""
    collection_str = str(collection or "").strip()
    stock_number_str = str(stock_number or "").strip()
    if collection_str and collection_str.lower() == "nan":
        collection_str = ""
    if stock_number_str and stock_number_str.lower() == "nan":
        stock_number_str = ""
    if collection_str and stock_number_str:
        return f"{collection_str} ({stock_number_str})"
    return stock_number_str or collection_str


def _normalize_source_stock_text(value: Any) -> str:
    """Normalize display text used for source and stock workbook columns."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _normalize_source_stock_series(
    values: Any,
    index: Optional[pd.Index] = None,
) -> pd.Series:
    """Return a normalized string series for source/stock display values."""
    if values is None:
        return pd.Series("", index=index, dtype=str) if index is not None else pd.Series(dtype=str)
    series = values if isinstance(values, pd.Series) else pd.Series(values, index=index)
    if index is not None:
        series = series.reindex(index)
    series = series.fillna("").astype(str).str.strip()
    return series.mask(series.str.lower().eq("nan"), "")


def _split_source_stock_label(source_stock: Any) -> Tuple[str, str]:
    """Split a combined source/stock label into workbook display columns."""
    source_stock_str = _normalize_source_stock_text(source_stock)
    if not source_stock_str:
        return "", ""
    match = re.fullmatch(r"(.+?)\s*\(([^()]*)\)\s*", source_stock_str)
    if match:
        return (
            _normalize_source_stock_text(match.group(1)),
            _normalize_source_stock_text(match.group(2)),
        )
    return "", source_stock_str


def _get_exported_source_stock_columns(
    sheet_df: pd.DataFrame,
) -> Tuple[pd.Series, pd.Series]:
    """Return normalized Source and Stock # columns for workbook output."""
    if sheet_df is None:
        return pd.Series(dtype=str), pd.Series(dtype=str)
    index = sheet_df.index
    if SOURCE_COLUMN in sheet_df.columns or STOCK_NUMBER_COLUMN in sheet_df.columns:
        return (
            _normalize_source_stock_series(sheet_df.get(SOURCE_COLUMN), index=index),
            _normalize_source_stock_series(sheet_df.get(STOCK_NUMBER_COLUMN), index=index),
        )
    parsed = _normalize_source_stock_series(
        sheet_df.get(SOURCE_STOCK_COLUMN, pd.Series("", index=index, dtype=str)),
        index=index,
    ).apply(_split_source_stock_label)
    return (
        parsed.map(lambda parts: parts[0]).astype(str),
        parsed.map(lambda parts: parts[1]).astype(str),
    )


def _get_source_stock_series(sheet_df: pd.DataFrame) -> pd.Series:
    """Return the normalized combined reagent label from either column layout."""
    if sheet_df is None:
        return pd.Series(dtype=str)
    index = sheet_df.index
    if SOURCE_STOCK_COLUMN in sheet_df.columns:
        return _normalize_source_stock_series(sheet_df[SOURCE_STOCK_COLUMN], index=index)
    source_series, stock_series = _get_exported_source_stock_columns(sheet_df)
    return pd.Series(
        [
            _format_source_stock_label(source, stock)
            for source, stock in zip(source_series.tolist(), stock_series.tolist())
        ],
        index=index,
        dtype=str,
    )


def _export_source_stock_columns(sheet_df: pd.DataFrame) -> pd.DataFrame:
    """Replace the combined source/stock column with Source and Stock # in place."""
    if sheet_df is None or len(sheet_df.columns) == 0 or SOURCE_STOCK_COLUMN not in sheet_df.columns:
        return sheet_df.copy() if sheet_df is not None else pd.DataFrame()
    source_series, stock_series = _get_exported_source_stock_columns(sheet_df)
    sheet_out = sheet_df.copy()
    ordered_columns: List[str] = []
    for column in sheet_out.columns:
        if column == SOURCE_STOCK_COLUMN:
            ordered_columns.extend([SOURCE_COLUMN, STOCK_NUMBER_COLUMN])
        elif column not in {SOURCE_COLUMN, STOCK_NUMBER_COLUMN}:
            ordered_columns.append(column)
    sheet_out = sheet_out.drop(
        columns=[SOURCE_STOCK_COLUMN, SOURCE_COLUMN, STOCK_NUMBER_COLUMN],
        errors='ignore',
    )
    sheet_out[SOURCE_COLUMN] = source_series
    sheet_out[STOCK_NUMBER_COLUMN] = stock_series
    return sheet_out[ordered_columns]


def _has_source_stock_columns(sheet_df: Optional[pd.DataFrame]) -> bool:
    """Return True when either the combined or split source/stock columns exist."""
    return (
        sheet_df is not None
        and any(
            column in sheet_df.columns
            for column in (SOURCE_STOCK_COLUMN, SOURCE_COLUMN, STOCK_NUMBER_COLUMN)
        )
    )


def _is_truthy(value: Any) -> bool:
    """Return True for common truthy encodings used in workbook data."""
    if value is True:
        return True
    if value is False or value is None:
        return False
    if isinstance(value, (int, float)) and not pd.isna(value):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _normalize_reagent_bucket_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure one-hot reagent bucket columns exist and are normalized to bools."""
    df = df.copy()
    for column in REAGENT_BUCKET_COLUMNS:
        if column not in df.columns:
            df[column] = False
        df[column] = df[column].map(_is_truthy)
    return df


def _recompute_reagent_bucket_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Recompute the one-hot reagent buckets from stock metadata."""
    if df is None or len(df) == 0:
        return df.copy()
    df = df.copy()
    bucket_df = df.apply(
        lambda row: pd.Series(_derive_one_hot_reagent_buckets(row)),
        axis=1,
    )
    for column in REAGENT_BUCKET_COLUMNS:
        df[column] = bucket_df[column].map(_is_truthy)
    _assert_one_hot_reagent_buckets(df, "Recomputed stock rows")
    return df


def _assert_one_hot_reagent_buckets(df: pd.DataFrame, context: str) -> None:
    """Raise if rows fail the exactly-one-true reagent bucket invariant."""
    if df is None or len(df) == 0:
        return
    if any(column not in df.columns for column in REAGENT_BUCKET_COLUMNS):
        return
    bucket_df = df[REAGENT_BUCKET_COLUMNS].apply(lambda col: col.map(_is_truthy))
    invalid_mask = bucket_df.sum(axis=1) != 1
    if invalid_mask.any():
        raise ValueError(
            f"{context} must have exactly one true reagent bucket per row."
        )


def _get_stock_phenotype_sheet_output_columns(
    similarity_columns: List[str],
) -> List[str]:
    """Return the visible Stock Phenotype Sheet column order."""
    return [
        'Gene',
        'Reagent Type or Allele Symbol',
        'Balancers',
        'matched_component_types',
        *REAGENT_BUCKET_COLUMNS,
        'allele_class_terms',
        'transgenic_product_class_terms',
        SOURCE_STOCK_COLUMN,
        'Genotype',
        CO_REAGENT_FBIDS_COLUMN,
        CO_REAGENT_SYMBOLS_COLUMN,
        PARTNER_DRIVER_SYMBOLS_COLUMN,
        PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN,
        'Data Set',
        'Phenotype',
        'Qualifier',
        *similarity_columns,
        'PMID',
        'PMCID',
        'Reference',
        'Authors',
        'Journal',
        'Year of Publication',
    ]


def _describe_stock_phenotype_sheet_column(column: str) -> str:
    """Return the workbook-contents definition for a phenotype-sheet column."""
    masterlist_definitions = {
        'Screening Group': "Template planning column kept blank by the pipeline.",
        'Actual cross set date': "Template planning column kept blank by the pipeline.",
        'Finished (=27+B)': "Template planning column kept blank by the pipeline.",
        'Projected Start Date': "Template planning column kept blank by the pipeline.",
        'Projected End Date': "Template planning column kept blank by the pipeline.",
        'Notes/Function': "Template notes column kept blank by the pipeline.",
        'allele shorthand': "Gene-relevant allele or reagent symbols extracted from the matched FlyBase components.",
        'Stock Source': "Collection or source label for the reagent when FlyBase or stock-center metadata provides one.",
        'ID #': "Stock-center stock number or custom reagent label used to identify the reagent.",
        'RNAi Shorthand': "Template column with no stable standalone source in the phenotype sheet; left blank.",
        'Full Stock Genotype': "FlyBase genotype label copied from the matched genotype_phenotype_data row.",
        'balancers in stock?': "FlyBase balancer symbols carried by the source stock. '-' indicates no balancer.",
        'ordered? (date)': "Template inventory column kept blank by the pipeline.",
        'location of stock': "Template inventory column kept blank by the pipeline.",
        'Published Gal4/ Positive control': "Best-effort partner GAL4 driver symbols, filtered to partner stocks whose derived components are GAL4-only.",
        'Published Gal4 source': "Best-effort FlyBase stock candidates for partner GAL4 drivers, filtered to GAL4-only partner stocks and formatted as '(<stock #>, <collection>)'.",
        'Published phenotype': "Normalized phenotype term derived from FlyBase genotype_phenotype_data.",
        'Reference': "PubMed ID resolved for the phenotype-supporting reference when available.",
        'Column 31': "PubMed Central ID resolved for the phenotype-supporting reference when available.",
        'Column 30': "Human-readable phenotype-supporting reference label, preferring title or mini-reference text.",
        'Column 29': "Authors for the phenotype-supporting reference when available.",
        'Column 1': "Journal for the phenotype-supporting reference when available.",
        'Column 32': "Publication year for the phenotype-supporting reference when available.",
        'Circadian/Sleep Relevance (embedding max score)': "Maximum cosine similarity across the configured phenotype similarity target columns.",
        'Experimental Driver': "Template experiment column kept blank by the pipeline.",
        'Gal4 control': "Template experiment column kept blank by the pipeline.",
        'RNAi control': "Template experiment column kept blank by the pipeline.",
        'Data Set': "Input CSV filename stem(s) for the gene symbols on this phenotype row.",
        'Sleep promoting?': "Template result column kept blank by the pipeline.",
        'Wake promoting?': "Template result column kept blank by the pipeline.",
        'link to graph     ': "Template results-link column kept blank by the pipeline.",
        '   link to data': "Template results-link column kept blank by the pipeline.",
    }
    if column in masterlist_definitions:
        return masterlist_definitions[column]
    definitions = {
        'Screening Group': (
            "Template planning column kept blank by the pipeline."
        ),
        'Actual cross set date': (
            "Template planning column kept blank by the pipeline."
        ),
        'Finished (=27+B)': (
            "Template planning column kept blank by the pipeline."
        ),
        'Projected Start Date': (
            "Template planning column kept blank by the pipeline."
        ),
        'Projected End Date': (
            "Template planning column kept blank by the pipeline."
        ),
        'Notes/Function': (
            "Template notes column kept blank by the pipeline."
        ),
        'allele shorthand': (
            "Gene-relevant allele or reagent symbols extracted from the matched FlyBase components."
        ),
        'Gene': (
            "Unique gene symbols linked to the matched stock components on this phenotype row."
        ),
        'Stock Source': (
            "Collection or source label for the reagent when FlyBase or stock-center metadata provides one."
        ),
        'ID #': (
            "Stock-center stock number or custom reagent label used to identify the reagent."
        ),
        'RNAi Shorthand': (
            "Template column with no stable standalone source in the phenotype sheet; left blank."
        ),
        'Full Stock Genotype': (
            "FlyBase genotype label copied from the matched genotype_phenotype_data row."
        ),
        'balancers in stock?': (
            "FlyBase balancer symbols carried by the source stock. '-' indicates no balancer."
        ),
        'ordered? (date)': (
            "Template inventory column kept blank by the pipeline."
        ),
        'location of stock': (
            "Template inventory column kept blank by the pipeline."
        ),
        'Published Gal4/ Positive control': (
            "Best-effort partner GAL4 driver symbols, filtered to partner stocks whose derived components are GAL4-only."
        ),
        'Published Gal4 source': (
            "Best-effort FlyBase stock candidates for partner GAL4 drivers, filtered to GAL4-only partner stocks and formatted as '(<stock #>, <collection>)'."
        ),
        'Published phenotype': (
            "Normalized phenotype term derived from FlyBase genotype_phenotype_data."
        ),
        'Reference': (
            "PubMed ID resolved for the phenotype-supporting reference when available."
        ),
        'Column 31': (
            "PubMed Central ID resolved for the phenotype-supporting reference when available."
        ),
        'Column 30': (
            "Human-readable phenotype-supporting reference label, preferring title or mini-reference text."
        ),
        'Column 29': (
            "Authors for the phenotype-supporting reference when available."
        ),
        'Column 1': (
            "Journal for the phenotype-supporting reference when available."
        ),
        'Column 32': (
            "Publication year for the phenotype-supporting reference when available."
        ),
        'Circadian/Sleep Relevance (embedding max score)': (
            "Maximum cosine similarity across the configured phenotype similarity target columns."
        ),
        'Experimental Driver': (
            "Template experiment column kept blank by the pipeline."
        ),
        'Gal4 control': (
            "Template experiment column kept blank by the pipeline."
        ),
        'RNAi control': (
            "Template experiment column kept blank by the pipeline."
        ),
        'Data Set': (
            "Input CSV filename stem(s) for the gene symbols on this phenotype row."
        ),
        'Sleep promoting?': (
            "Template result column kept blank by the pipeline."
        ),
        'Wake promoting?': (
            "Template result column kept blank by the pipeline."
        ),
        'link to graph     ': (
            "Template results-link column kept blank by the pipeline."
        ),
        '   link to data': (
            "Template results-link column kept blank by the pipeline."
        ),
        'Reagent Type or Allele Symbol': (
            "Gene-relevant allele or reagent symbols extracted from the matched FlyBase components."
        ),
        'Balancers': (
            "FlyBase balancer symbols carried by the source stock. '-' indicates no balancer."
        ),
        'matched_component_types': (
            "FlyBase component ID types that linked the reagent to the input gene(s), aggregated across matches (for example FBal, FBtp, FBti)."
        ),
        'UAS': (
            "One-hot reagent bucket. True only when the reagent resolves to the UAS bucket using the widened UAS proxy, which includes UAS or RNAi signals, rnai_reagent class terms, and Vienna-style knockdown-family signals."
        ),
        'GAL4': (
            "One-hot reagent bucket. True only when GAL4 or driver signals are present and the reagent does not qualify for a higher-precedence mixed bucket."
        ),
        'mutant/UAS': (
            "One-hot reagent bucket. True when mutant and UAS signals co-occur, unless the reagent is promoted to 'GAL4 / mutant' by the precedence rule."
        ),
        'mutant': (
            "One-hot reagent bucket. True for direct non-transgenic mutant-like reagents when neither mixed bucket wins."
        ),
        'GAL4 / mutant': (
            "One-hot reagent bucket. True whenever GAL4 and mutant signals co-occur, including GAL4 + mutant + UAS reagents under the selected precedence."
        ),
        'Other': (
            "One-hot reagent bucket. True for reagents that do not land in the named buckets, including UAS-GAL4 style reagents without mutant evidence."
        ),
        'allele_class_terms': (
            "FlyBase allele class terms aggregated from the gene-relevant allele(s) associated with the reagent."
        ),
        'transgenic_product_class_terms': (
            "FlyBase transgenic product class terms aggregated from the gene-relevant construct-linked reagent(s)."
        ),
        SOURCE_COLUMN: (
            "Collection or source label for the reagent when FlyBase or stock-center metadata provides one."
        ),
        STOCK_NUMBER_COLUMN: (
            "Stock-center stock number or custom reagent label used to identify the reagent."
        ),
        SOURCE_STOCK_COLUMN: (
            "Display label for the reagent source, typically 'Collection (stock_number)' for stock-center lines or the custom reagent label when no FBst stock exists."
        ),
        'Genotype': (
            "FlyBase genotype label copied from the matched genotype_phenotype_data row."
        ),
        CO_REAGENT_FBIDS_COLUMN: (
            "FlyBase IDs for partner genotype components that are not on the focal matched stock and that resolve specifically to GAL4-linked co-reagents."
        ),
        CO_REAGENT_SYMBOLS_COLUMN: (
            "Best-effort symbols for partner GAL4 co-reagents recovered from aligned genotype text or FlyBase stock-component lookups."
        ),
        PARTNER_DRIVER_SYMBOLS_COLUMN: (
            "Best-effort partner GAL4 driver symbols inferred from non-focal genotype components and GAL4-containing genotype text, filtered to drivers with GAL4-only stock candidates."
        ),
        PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN: (
            "Best-effort FlyBase stock candidates for partner GAL4 driver components, filtered to GAL4-only partner stocks and formatted as '(<stock #>, <collection>)'."
        ),
        'Phenotype': (
            "Normalized phenotype term derived from FlyBase genotype_phenotype_data."
        ),
        'Qualifier': (
            "Normalized FlyBase phenotype qualifier text for the matched phenotype record."
        ),
        'PMID': (
            "PubMed ID resolved for the phenotype-supporting reference when available."
        ),
        'PMCID': (
            "PubMed Central ID resolved for the phenotype-supporting reference when available."
        ),
        'Reference': (
            "Human-readable phenotype-supporting reference label, preferring title or mini-reference text."
        ),
        'Authors': (
            "Authors for the phenotype-supporting reference when available."
        ),
        'Journal': (
            "Journal for the phenotype-supporting reference when available."
        ),
        'Year of Publication': (
            "Publication year for the phenotype-supporting reference when available."
        ),
    }
    if column in definitions:
        return definitions[column]
    if column.startswith('Cosine Similarity (') and column.endswith(')'):
        target = column[len('Cosine Similarity ('):-1]
        return (
            f"Cosine similarity between the phenotype text and the configured '{target}' target concept."
        )
    return "Phenotype-sheet column carried through from the stock and reference aggregation pipeline."


def _get_stock_phenotype_sheet_column_definitions(
    phenotype_sheet_df: pd.DataFrame,
    format_as_masterlist: bool = False,
) -> List[Tuple[str, str]]:
    """Return visible Stock Phenotype Sheet columns with workbook definitions."""
    if phenotype_sheet_df is None or len(phenotype_sheet_df.columns) == 0:
        return []
    phenotype_sheet_df = (
        _reorder_to_masterlist_columns(phenotype_sheet_df)
        if format_as_masterlist
        else _export_source_stock_columns(phenotype_sheet_df)
    )
    visible_columns = [
        column for column in phenotype_sheet_df.columns
        if column != '_reference_url'
    ]
    return [
        (column, _describe_stock_phenotype_sheet_column(column))
        for column in visible_columns
    ]


def _get_reagent_bucket_one_hot_note() -> str:
    """Return the shared workbook note for the mutually exclusive bucket set."""
    return (
        "The reagent-bucket columns "
        "`UAS`, `GAL4`, `mutant/UAS`, `mutant`, `GAL4 / mutant`, and `Other` "
        "form a mutually exclusive one-hot set. Exactly one of these columns "
        "should be true for each reagent row."
    )


def _normalize_simple_bucket_collection(row: pd.Series) -> str:
    """Normalize the collection label used by simple-bucket summaries."""
    collection = str(row.get("collection", "") or "").strip()
    if not collection or collection.lower() == "nan":
        collection = "Unknown"
    return collection


def _build_simple_bucket_reagent_key(row: pd.Series) -> str:
    """Return a stable reagent key across stock-center and custom reagents."""
    fbst = clean_id(row.get("FBst", ""))
    if fbst:
        return fbst
    stock_number = str(row.get("stock_number", "") or "").strip()
    collection = _normalize_simple_bucket_collection(row)
    component_ids = unique_join(parse_semicolon_list(str(row.get("relevant_component_ids", "") or "")))
    return f"CUSTOM::{collection}::{stock_number}::{component_ids}"


def _sanitize_simple_bucket_sheet_component(value: str, max_len: int = 12) -> str:
    """Build a compact Excel-safe sheet-name token."""
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "").strip()).strip("_")
    if not token:
        token = "bucket"
    return token[:max_len]


def _format_simple_bucket_sheet_name(
    index: int,
    collection: str,
    uas: bool,
    sleep_circ: bool,
    has_balancer: bool,
) -> str:
    """Build a deterministic, Excel-safe simple-bucket sheet name."""
    collection_token = _sanitize_simple_bucket_sheet_component(collection, max_len=12)
    base_name = (
        f"{index:02d}_{collection_token}_u{int(uas)}_s{int(sleep_circ)}_b{int(has_balancer)}"
    )
    return base_name[:31]


def _format_simple_bucket_combination_label(
    collection: str,
    uas: bool,
    sleep_circ: bool,
    has_balancer: bool,
) -> str:
    """Build the user-facing combination label for simple buckets."""
    return (
        f"{collection} > UAS ({str(uas).lower()}) > sleep/ circ ({str(sleep_circ).lower()}) "
        f"> has balancer ({str(has_balancer).lower()})"
    )


def _read_gene_set_sheet_for_similarity_workbook(source_workbook_path: Path) -> Optional[pd.DataFrame]:
    """Load the input gene set, preferring the original CSV gene-list inputs."""
    workbook_path = Path(source_workbook_path)
    csv_dirs: List[Path] = []
    for candidate in (workbook_path.parent, workbook_path.parent.parent):
        if candidate not in csv_dirs:
            csv_dirs.append(candidate)

    csv_frames: List[pd.DataFrame] = []
    for csv_dir in csv_dirs:
        for csv_path in sorted(csv_dir.glob("*.csv")):
            try:
                csv_frames.append(pd.read_csv(csv_path, dtype=str))
            except Exception:
                continue
        if csv_frames:
            break

    if csv_frames:
        if len(csv_frames) == 1:
            return csv_frames[0]
        return pd.concat(csv_frames, ignore_index=True, sort=False)

    source_workbook = pd.ExcelFile(workbook_path)
    if not source_workbook.sheet_names:
        return None
    return pd.read_excel(workbook_path, sheet_name=source_workbook.sheet_names[0], header=None)


def _build_simple_bucket_workbook_entries(
    phenotype_sheet_df: pd.DataFrame,
    combination_outputs: List[Tuple[List[str], pd.DataFrame, Dict[str, Any]]],
    csv_input_genes: Optional[Set[str]] = None,
) -> List[Tuple[str, pd.DataFrame, Dict[str, Any]]]:
    """Build ordered simple-bucket workbook entries and counts."""
    if phenotype_sheet_df is None or len(phenotype_sheet_df) == 0:
        return []

    phenotype_sources = _get_source_stock_series(phenotype_sheet_df)
    phenotype_source_set = {value for value in phenotype_sources if value}
    if not phenotype_source_set:
        return []

    sleep_circ_sources: Set[str] = set()
    if "Phenotype" in phenotype_sheet_df.columns:
        phenotype_sheet_with_source = phenotype_sheet_df.assign(
            _source_stock_key=phenotype_sources
        )
        for source_stock, phenotype_group in phenotype_sheet_with_source.groupby(
            "_source_stock_key",
            sort=False,
        ):
            source_stock_str = str(source_stock or "").strip()
            if not source_stock_str:
                continue
            phenotypes = phenotype_group["Phenotype"].fillna("").astype(str)
            if phenotypes.str.contains(r"sleep|circadian", case=False, na=False).any():
                sleep_circ_sources.add(source_stock_str)

    reagent_records: List[Dict[str, Any]] = []
    seen_reagents: Set[str] = set()
    allele_col_cache: Dict[int, Optional[str]] = {}
    collection_order: List[str] = []

    for _combo, stock_df, _summary in combination_outputs:
        if stock_df is None or len(stock_df) == 0:
            continue
        allele_col = allele_col_cache.setdefault(id(stock_df), _find_allele_column(stock_df))
        gene_col = "relevant_gene_symbols" if "relevant_gene_symbols" in stock_df.columns else _find_gene_id_column(stock_df)
        for _, row in stock_df.iterrows():
            source_stock = _format_source_stock_label(
                row.get("collection", ""),
                row.get("stock_number", ""),
            )
            if not source_stock or source_stock not in phenotype_source_set:
                continue

            reagent_key = _build_simple_bucket_reagent_key(row)
            if reagent_key in seen_reagents:
                continue
            seen_reagents.add(reagent_key)

            collection = _normalize_simple_bucket_collection(row)
            if collection not in collection_order:
                collection_order.append(collection)

            genes = {
                gene.strip()
                for gene in parse_semicolon_list(str(row.get(gene_col, "") or ""))
                if gene.strip() and gene.strip() != "-"
            }
            if csv_input_genes is not None:
                genes &= csv_input_genes

            alleles: Set[str] = set()
            if allele_col:
                for allele in parse_semicolon_list(str(row.get(allele_col, "") or "")):
                    allele = allele.strip()
                    if allele and allele != "-":
                        alleles.add(allele)
            if csv_input_genes is not None and not genes:
                alleles.clear()

            has_balancer = False
            try:
                has_balancer = int(float(row.get("num_Balancers", 0) or 0)) > 0
            except (TypeError, ValueError):
                balancers_value = str(row.get("Balancers", "") or "").strip()
                has_balancer = bool(balancers_value and balancers_value not in {"-", "nan"})

            reagent_records.append(
                {
                    "reagent_key": reagent_key,
                    "source_stock": source_stock,
                    "collection": collection,
                    "uas": _is_truthy(row.get("UAS", False)),
                    "sleep_circ": source_stock in sleep_circ_sources,
                    "has_balancer": has_balancer,
                    "genes": genes,
                    "alleles": alleles,
                }
            )

    if not reagent_records:
        return []

    records_by_bucket: Dict[Tuple[str, bool, bool, bool], List[Dict[str, Any]]] = defaultdict(list)
    for record in reagent_records:
        bucket_key = (
            record["collection"],
            record["uas"],
            record["sleep_circ"],
            record["has_balancer"],
        )
        records_by_bucket[bucket_key].append(record)

    global_seen_genes: Set[str] = set()
    global_seen_alleles: Set[str] = set()
    entries: List[Tuple[str, pd.DataFrame, Dict[str, Any]]] = []
    entry_index = 1
    for collection in collection_order:
        for uas in (True, False):
            for sleep_circ in (True, False):
                for has_balancer in (False, True):
                    bucket_key = (collection, uas, sleep_circ, has_balancer)
                    bucket_records = records_by_bucket.get(bucket_key, [])
                    source_stocks = {
                        record["source_stock"]
                        for record in bucket_records
                        if record["source_stock"]
                    }
                    bucket_df = phenotype_sheet_df[
                        phenotype_sources.isin(source_stocks)
                    ].copy()

                    combo_genes: Set[str] = set()
                    combo_alleles: Set[str] = set()
                    for record in bucket_records:
                        combo_genes.update(record["genes"])
                        combo_alleles.update(record["alleles"])

                    owned_genes = combo_genes - global_seen_genes
                    owned_alleles = combo_alleles - global_seen_alleles
                    global_seen_genes.update(owned_genes)
                    global_seen_alleles.update(owned_alleles)

                    sheet_name = _format_simple_bucket_sheet_name(
                        entry_index,
                        collection=collection,
                        uas=uas,
                        sleep_circ=sleep_circ,
                        has_balancer=has_balancer,
                    )
                    entry_index += 1

                    metadata = {
                        "sheet_name": sheet_name,
                        "combination": _format_simple_bucket_combination_label(
                            collection=collection,
                            uas=uas,
                            sleep_circ=sleep_circ,
                            has_balancer=has_balancer,
                        ),
                        "collection": collection,
                        "uas": uas,
                        "sleep_circ": sleep_circ,
                        "has_balancer": has_balancer,
                        "stock_count": len(source_stocks),
                        "allele_count": len(owned_alleles),
                        "gene_count": len(owned_genes),
                    }
                    entries.append((sheet_name, bucket_df, metadata))

    return entries


def _write_similarity_tier_contents_sheet(
    workbook,
    similarity_tiers: List[Tuple[str, pd.DataFrame, Dict[str, Any]]],
    phenotype_sheet_df: Optional[pd.DataFrame] = None,
    format_as_masterlist: bool = False,
) -> None:
    """Write a concise contents sheet for the similarity tier workbook."""
    worksheet = workbook.add_worksheet("Contents")
    title_fmt = workbook.add_format({"bold": True, "font_size": 14})
    body_fmt = workbook.add_format({"font_size": 11, "text_wrap": True, "valign": "top"})
    header_fmt = workbook.add_format({"bold": True, "bottom": 1, "font_size": 11})

    worksheet.set_column(0, 0, 26)
    worksheet.set_column(1, 1, 70)
    worksheet.set_column(2, 2, 16)

    row = 0
    worksheet.write(row, 0, "Tier Workbook Contents", title_fmt)
    row += 2
    worksheet.write(
        row,
        0,
        "Bucket assignment",
        header_fmt,
    )
    worksheet.write(
        row,
        1,
        "Rows are assigned by Max Cosine Similarity, defined as the per-row maximum across all `Cosine Similarity (...)` columns.",
        body_fmt,
    )
    row += 1
    worksheet.write(row, 0, "Thresholds", header_fmt)
    worksheet.write(
        row,
        1,
        "Fixed 0.05 cosine bins are evaluated from 0.95-1.0 downward to <0.05; empty buckets are skipped.",
        body_fmt,
    )
    row += 1
    worksheet.write(row, 0, "Ordering within tiers", header_fmt)
    worksheet.write(
        row,
        1,
        "Rows are grouped by gene using only reagents present in that tier. Gene groups with more unique reagents appear first, and reagents within each gene are sorted by max cosine similarity.",
        body_fmt,
    )
    row += 2

    worksheet.write_row(row, 0, ["Sheet", "Meaning", "Rows"], header_fmt)
    row += 1
    worksheet.write_row(row, 0, ["Gene Set", "Copied from the input gene-list CSV data when available.", ""])
    row += 1
    worksheet.write_row(row, 0, ["Stock Phenotype Sheet", "Full phenotype table used to build all similarity tiers.", ""])
    row += 1
    for sheet_name, tier_df, metadata in similarity_tiers:
        worksheet.write_row(
            row,
            0,
            [
                sheet_name,
                f"Tier for similarity range {metadata['range_label']}.",
                int(len(tier_df)),
            ],
        )
        row += 1
    column_definitions = _get_stock_phenotype_sheet_column_definitions(
        phenotype_sheet_df if phenotype_sheet_df is not None else pd.DataFrame(),
        format_as_masterlist=format_as_masterlist,
    )
    if column_definitions:
        row += 2
        worksheet.write(row, 0, "Stock Phenotype Sheet columns", header_fmt)
        worksheet.write(row, 1, _get_reagent_bucket_one_hot_note(), body_fmt)
        row += 1
        worksheet.write_row(row, 0, ["Column", "Definition"], header_fmt)
        row += 1
        for column, definition in column_definitions:
            worksheet.write(row, 0, column, body_fmt)
            worksheet.write(row, 1, definition, body_fmt)
            row += 1


def _format_bucket_combo_header(uas: bool, sleep_circ: bool, has_balancer: bool) -> str:
    """Short human-readable column header for one boolean permutation."""
    parts: List[str] = []
    parts.append("UAS" if uas else "non-UAS")
    parts.append("slp/ circ" if sleep_circ else "Non slp/ circ")
    parts.append("Has bal" if has_balancer else "No bal")
    return " | ".join(parts)


def _write_simple_bucket_contents_sheet(
    workbook,
    simple_bucket_entries: List[Tuple[str, pd.DataFrame, Dict[str, Any]]],
    phenotype_sheet_df: Optional[pd.DataFrame] = None,
    format_as_masterlist: bool = False,
) -> None:
    """Write a pivot-matrix contents sheet for simple-bucket workbooks.

    Renders three stacked matrices (Stocks, Alleles, Genes) with rows =
    collections and columns = boolean permutations that have at least one
    non-zero cell.  Each matrix includes row totals and a column-total
    footer row.  A sheet-name legend follows the matrices.
    """
    worksheet = workbook.add_worksheet("Contents")
    title_fmt = workbook.add_format({"bold": True, "font_size": 14})
    body_fmt = workbook.add_format({"font_size": 11, "text_wrap": True, "valign": "top"})
    header_fmt = workbook.add_format({"bold": True, "bottom": 1, "font_size": 11})
    num_fmt = workbook.add_format({"font_size": 11, "align": "center"})
    num_zero_fmt = workbook.add_format({"font_size": 11, "align": "center", "font_color": "#BFBFBF"})
    total_fmt = workbook.add_format({"bold": True, "font_size": 11, "top": 1, "align": "center"})
    total_label_fmt = workbook.add_format({"bold": True, "font_size": 11, "top": 1})
    row_total_hdr_fmt = workbook.add_format({"bold": True, "bottom": 1, "font_size": 11, "align": "center"})
    section_fmt = workbook.add_format({"bold": True, "font_size": 12})
    legend_hdr_fmt = workbook.add_format({"bold": True, "bottom": 1, "font_size": 11})
    legend_fmt = workbook.add_format({"font_size": 11})

    # -- Collect per-entry data keyed by (collection, combo_tuple) ----------
    collections_ordered: List[str] = []
    all_combos_ordered: List[Tuple[bool, bool, bool]] = []
    data: Dict[Tuple[str, Tuple[bool, bool, bool]], Dict[str, int]] = {}
    legend_rows: List[Tuple[str, Dict[str, Any]]] = []

    for sheet_name, _bucket_df, metadata in simple_bucket_entries:
        coll: str = metadata["collection"]
        combo = (metadata["uas"], metadata["sleep_circ"], metadata["has_balancer"])
        if coll not in collections_ordered:
            collections_ordered.append(coll)
        if combo not in all_combos_ordered:
            all_combos_ordered.append(combo)
        data[(coll, combo)] = {
            "stock_count": metadata["stock_count"],
            "allele_count": metadata["allele_count"],
            "gene_count": metadata["gene_count"],
        }
        legend_rows.append((sheet_name, metadata))

    # Keep only combos that have at least one non-zero stock across all
    # collections (avoids entirely-empty columns).
    active_combos = [
        combo for combo in all_combos_ordered
        if any(data.get((c, combo), {}).get("stock_count", 0) for c in collections_ordered)
    ]
    if not active_combos:
        active_combos = all_combos_ordered

    combo_headers = [_format_bucket_combo_header(*c) for c in active_combos]
    n_combos = len(active_combos)

    # Column widths: col 0 = collection label, cols 1..n = combo columns, col n+1 = total
    worksheet.set_column(0, 0, 18)
    for ci in range(1, n_combos + 2):
        worksheet.set_column(ci, ci, max(14, max((len(h) for h in combo_headers), default=10) + 2))

    row = 0
    worksheet.write(row, 0, "Simple Bucket Workbook Contents", title_fmt)
    row += 2

    # -- Preamble ----------------------------------------------------------
    worksheet.write(row, 0, "Bucket assignment", header_fmt)
    worksheet.write(
        row, 1,
        "Each reagent is assigned to exactly one collection / UAS / sleep-circ / balancer combination. "
        "The matrices below show per-collection counts for every non-empty combination.",
        body_fmt,
    )
    row += 1
    worksheet.write(row, 0, "Counting invariant", header_fmt)
    worksheet.write(
        row, 1,
        "Genes, alleles, and reagents are never double-counted within or across combinations. "
        "Gene and allele counts are owned by the first matching combination in workbook order.",
        body_fmt,
    )
    row += 2

    # -- Column header definitions -----------------------------------------
    worksheet.write(row, 0, "Column definitions", section_fmt)
    row += 1
    worksheet.write(row, 0, "UAS / non-UAS", header_fmt)
    worksheet.write(
        row, 1,
        "Whether the stock genotype contains a UAS (Upstream Activation Sequence) construct. "
        "UAS stocks carry transgenes driven by the GAL4-UAS binary expression system.",
        body_fmt,
    )
    row += 1
    worksheet.write(row, 0, "slp/ circ / Non slp/ circ", header_fmt)
    worksheet.write(
        row, 1,
        "Whether any phenotype row linked to the stock mentions 'sleep' or 'circadian'. "
        "Derived from the FlyBase genotype_phenotype_data phenotype descriptions.",
        body_fmt,
    )
    row += 1
    worksheet.write(row, 0, "No bal / Has bal", header_fmt)
    worksheet.write(
        row, 1,
        "Whether the stock carries at least one balancer chromosome (e.g. CyO, TM3, TM6B, FM7). "
        "Derived from the FBba (balancer) component count in the Stage 1 workbook.",
        body_fmt,
    )
    row += 2

    # -- Helper to write one pivot matrix ----------------------------------
    def _write_matrix(start_row: int, metric_label: str, metric_key: str) -> int:
        r = start_row
        worksheet.write(r, 0, metric_label, section_fmt)
        r += 1
        # Header row: Collection | combo1 | combo2 | ... | Total
        worksheet.write(r, 0, "Collection", header_fmt)
        for ci, hdr in enumerate(combo_headers):
            worksheet.write(r, 1 + ci, hdr, header_fmt)
        worksheet.write(r, 1 + n_combos, "Total", row_total_hdr_fmt)
        r += 1

        col_totals = [0] * n_combos
        grand_total = 0
        for coll in collections_ordered:
            worksheet.write(r, 0, coll, legend_fmt)
            row_sum = 0
            for ci, combo in enumerate(active_combos):
                val = data.get((coll, combo), {}).get(metric_key, 0)
                fmt = num_zero_fmt if val == 0 else num_fmt
                worksheet.write_number(r, 1 + ci, val, fmt)
                col_totals[ci] += val
                row_sum += val
            worksheet.write_number(r, 1 + n_combos, row_sum, num_fmt)
            grand_total += row_sum
            r += 1

        # Column totals footer
        worksheet.write(r, 0, "Total", total_label_fmt)
        for ci, ct in enumerate(col_totals):
            worksheet.write_number(r, 1 + ci, ct, total_fmt)
        worksheet.write_number(r, 1 + n_combos, grand_total, total_fmt)
        r += 1
        return r

    # -- Write three matrices ----------------------------------------------
    row = _write_matrix(row, "Stocks per bucket", "stock_count")
    row += 1
    row = _write_matrix(row, "Alleles per bucket", "allele_count")
    row += 1
    row = _write_matrix(row, "Genes per bucket", "gene_count")
    row += 2

    # -- Sheet-name legend -------------------------------------------------
    worksheet.write(row, 0, "Sheet legend", section_fmt)
    row += 1
    legend_headers = [
        "Sheet name", "Combination", "Collection",
        "UAS", "sleep/ circ", "has balancer",
        "# Stocks", "# Alleles", "# Genes",
    ]
    for ci, hdr in enumerate(legend_headers):
        worksheet.write(row, ci, hdr, legend_hdr_fmt)
    row += 1
    for sn, meta in legend_rows:
        worksheet.write(row, 0, sn, legend_fmt)
        worksheet.write(row, 1, meta["combination"], legend_fmt)
        worksheet.write(row, 2, meta["collection"], legend_fmt)
        worksheet.write(row, 3, str(meta["uas"]).lower(), legend_fmt)
        worksheet.write(row, 4, str(meta["sleep_circ"]).lower(), legend_fmt)
        worksheet.write(row, 5, str(meta["has_balancer"]).lower(), legend_fmt)
        worksheet.write_number(row, 6, meta["stock_count"], num_fmt)
        worksheet.write_number(row, 7, meta["allele_count"], num_fmt)
        worksheet.write_number(row, 8, meta["gene_count"], num_fmt)
        row += 1
    column_definitions = _get_stock_phenotype_sheet_column_definitions(
        phenotype_sheet_df if phenotype_sheet_df is not None else pd.DataFrame(),
        format_as_masterlist=format_as_masterlist,
    )
    if column_definitions:
        definition_col = n_combos + 3
        worksheet.set_column(definition_col, definition_col, 28)
        worksheet.set_column(definition_col + 1, definition_col + 1, 80)
        row += 2
        worksheet.write(row, definition_col, "Stock Phenotype Sheet columns", section_fmt)
        row += 1
        worksheet.write(row, definition_col, "Column", header_fmt)
        worksheet.write(row, definition_col + 1, "Definition", header_fmt)
        row += 1
        worksheet.write(row, definition_col, "Note", header_fmt)
        worksheet.write(
            row,
            definition_col + 1,
            _get_reagent_bucket_one_hot_note(),
            body_fmt,
        )
        row += 1
        for column, definition in column_definitions:
            worksheet.write(row, definition_col, column, legend_fmt)
            worksheet.write(row, definition_col + 1, definition, body_fmt)
            row += 1


def _write_similarity_tier_workbook(
    output_path: Path,
    source_workbook_path: Path,
    stock_phenotype_sheet_df: pd.DataFrame,
    combination_outputs: List[Tuple[List[str], pd.DataFrame, Dict[str, Any]]],
    csv_input_genes: Optional[Set[str]] = None,
    simple_buckets: bool = False,
    verbose: bool = False,
    format_as_masterlist: bool = False,
) -> Optional[Path]:
    """Write the phenotype similarity workbook alongside the aggregated workbook."""
    if stock_phenotype_sheet_df is None or len(stock_phenotype_sheet_df) == 0:
        return None

    similarity_tiers = _build_similarity_tier_sheets(stock_phenotype_sheet_df)
    simple_bucket_entries = _build_simple_bucket_workbook_entries(
        phenotype_sheet_df=stock_phenotype_sheet_df,
        combination_outputs=combination_outputs,
        csv_input_genes=csv_input_genes,
    ) if simple_buckets else []
    gene_set_df = _read_gene_set_sheet_for_similarity_workbook(source_workbook_path)
    tier_workbook_path = output_path.parent / f"{output_path.stem}_similarity_tiers.xlsx"
    with pd.ExcelWriter(tier_workbook_path, engine='xlsxwriter') as writer:
        workbook = writer.book
        if simple_buckets:
            _write_simple_bucket_contents_sheet(
                workbook,
                simple_bucket_entries,
                phenotype_sheet_df=stock_phenotype_sheet_df,
                format_as_masterlist=format_as_masterlist,
            )
        else:
            _write_similarity_tier_contents_sheet(
                workbook,
                similarity_tiers,
                phenotype_sheet_df=stock_phenotype_sheet_df,
                format_as_masterlist=format_as_masterlist,
            )
        if gene_set_df is not None and len(gene_set_df.columns) > 0:
            write_header = not all(isinstance(col, int) for col in gene_set_df.columns)
            gene_set_df.to_excel(writer, sheet_name="Gene Set", index=False, header=write_header)
        _write_phenotype_similarity_sheet(
            writer,
            workbook,
            "Stock Phenotype Sheet",
            stock_phenotype_sheet_df,
            format_as_masterlist=format_as_masterlist,
        )
        if simple_buckets:
            for sheet_name, bucket_df, _metadata in simple_bucket_entries:
                _write_phenotype_similarity_sheet(
                    writer,
                    workbook,
                    sheet_name,
                    bucket_df,
                    format_as_masterlist=format_as_masterlist,
                )
        else:
            for sheet_name, tier_df, _metadata in similarity_tiers:
                _write_phenotype_similarity_sheet(
                    writer,
                    workbook,
                    sheet_name,
                    tier_df,
                    format_as_masterlist=format_as_masterlist,
                )

    if verbose:
        sheet_count = len(simple_bucket_entries) if simple_buckets else len(similarity_tiers)
        sheet_kind = "simple bucket" if simple_buckets else "similarity tier"
        print(f"    Saved: {tier_workbook_path.name} ({sheet_count} {sheet_kind} sheet(s))")
    return tier_workbook_path


def _build_stock_sheet_by_gene(
    combination_outputs: List[Tuple[List[str], pd.DataFrame, Dict]],
    references_df: Optional[pd.DataFrame],
    csv_input_genes: Optional[Set[str]] = None,
    gene_synonyms_map: Optional[Dict[str, str]] = None,
    current_to_input_map: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Build a grouped sheet with unique (stock, PMID) keyword-hit rows.
    
    Inclusion rules:
    - Stocks must come from defined output categories (Sheet1..N dataframes).
    - Stocks must be Ref++ (>=1 keyword-hit title/abstract reference).
    - Rows are unique (stock, PMID) pairs for keyword-hit PMIDs only.
    """
    # Collect all stocks that made it into defined category sheets.
    included_rows: List[pd.Series] = []
    for _combo, out_df, _summary in combination_outputs:
        if out_df is None or len(out_df) == 0:
            continue
        included_rows.extend([row for _, row in out_df.iterrows()])
    
    if not included_rows:
        return pd.DataFrame()
    
    # Build PMID -> reference metadata map and keyword-hit set from References sheet.
    ref_meta: Dict[str, Dict[str, str]] = {}
    keyword_hit_pmids: Set[str] = set()
    if references_df is not None and len(references_df) > 0 and 'PMID' in references_df.columns:
        kw_ref_col = find_keyword_title_abstract_column(references_df)
        for _, r in references_df.iterrows():
            pmid = clean_id(r.get('PMID', ''))
            if not pmid:
                continue
            ref_meta[pmid] = {
                'title': str(r.get('title', '') or ''),
                'journal': str(r.get('journal', '') or ''),
                'publication_date': str(r.get('publication_date', '') or ''),
                'authors': str(r.get('author(s)', '') or ''),
            }
            if kw_ref_col:
                kw_val = r.get(kw_ref_col, False)
                if kw_val is True or str(kw_val).strip().lower() in ('true', '1', 'yes'):
                    keyword_hit_pmids.add(pmid)
    
    # If keyword-hit PMIDs are unavailable from references, we can still use stock-level
    # keyword PMID lists where present.
    sample_df = pd.DataFrame(included_rows)
    stock_col = _find_stock_id_column(sample_df)
    if stock_col is None:
        return pd.DataFrame()
    
    pmid_col = None
    for c in ['PMID', 'pmid', 'PMIDs']:
        if c in sample_df.columns:
            pmid_col = c
            break
    if pmid_col is None:
        return pd.DataFrame()
    
    kw_count_col = _find_keyword_ref_count_column(sample_df)
    kw_pmids_col = _find_keyword_ref_pmids_column(sample_df)
    keyword_refs_header = _derive_stock_keyword_refs_header(kw_pmids_col)
    has_keyword_signal = bool(kw_count_col or kw_pmids_col or keyword_hit_pmids)
    if not has_keyword_signal:
        return pd.DataFrame()
    validity_col = _find_column_with_experimental_fallback(
        list(sample_df.columns), "Functionally Valid Stock?"
    )
    
    gpt_cols = []
    seen_gpt = set()
    for c in get_gpt_derived_columns(sample_df):
        if c not in seen_gpt:
            gpt_cols.append(c)
            seen_gpt.add(c)
    for c in sample_df.columns:
        if str(c).startswith(EXPERIMENTAL_PREFIX) and c not in seen_gpt:
            gpt_cols.append(c)
            seen_gpt.add(c)
    
    rows: List[Dict[str, Any]] = []
    seen_pairs: Set[Tuple[str, str]] = set()
    
    for row in included_rows:
        stock_num = clean_id(row.get(stock_col, ''))
        if not stock_num:
            continue
        
        # Pick a single gene label per stock (prefer genes from input list).
        genes = parse_semicolon_list(row.get('relevant_gene_symbols', ''))
        if not genes:
            # Gene curation for this sheet is exclusively from relevant_gene_symbols.
            continue
        if csv_input_genes:
            genes_filtered = [g for g in genes if g in csv_input_genes]
            if not genes_filtered:
                continue
            gene_label = genes_filtered[0]
        else:
            gene_label = genes[0]
        
        if current_to_input_map:
            gene_label = current_to_input_map.get(gene_label, gene_label)
        
        # Determine keyword-hit PMIDs for this stock.
        pmids_all = [clean_id(p) for p in str(row.get(pmid_col, '')).split(';') if clean_id(p)]
        
        pmids_hit: List[str] = []
        if kw_pmids_col and str(row.get(kw_pmids_col, '')).strip():
            pmids_hit = [clean_id(p) for p in str(row.get(kw_pmids_col, '')).split(';') if clean_id(p)]
        elif keyword_hit_pmids:
            pmids_hit = [p for p in pmids_all if p in keyword_hit_pmids]
        else:
            pmids_hit = pmids_all
        
        # Ref++ stock gate: must have at least one keyword-hit PMID.
        kw_ref_count = 0
        if kw_count_col:
            try:
                kw_ref_count = int(row.get(kw_count_col, 0))
            except (ValueError, TypeError):
                kw_ref_count = 0
        if kw_ref_count <= 0:
            kw_ref_count = len(set(pmids_hit))
        if kw_ref_count <= 0:
            continue
        
        # Determine stock-level functional validity ordering.
        is_functionally_valid = False
        if validity_col:
            validity_val = str(row.get(validity_col, '') or '').strip().lower()
            is_functionally_valid = validity_val == ValidationStatus.FUNCTIONALLY_VALIDATED.lower()
        
        for pmid in sorted(set(pmids_hit)):
            pair = (stock_num, pmid)
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            
            meta = ref_meta.get(pmid, {})
            kw_refs_for_stock = (
                _normalize_pmid_list_string(row.get(kw_pmids_col, ""))
                if kw_pmids_col and kw_pmids_col in row.index
                else _normalize_pmid_list_string(";".join(pmids_hit))
            )
            all_refs_for_stock = _normalize_pmid_list_string(row.get(pmid_col, ""))

            out_row: Dict[str, Any] = {
                'gene': gene_label,
                'gene synonyms': _lookup_gene_synonyms(gene_label, gene_synonyms_map),
                'stock #': stock_num,
                'Balancers': str(row.get('Balancers', '') or '').strip(),
                'matched_component_types': str(
                    row.get('matched_component_types', '') or ''
                ).strip(),
                'pmid': pmid,
                'title': meta.get('title', ''),
                'journal': meta.get('journal', ''),
                'publication date': meta.get('publication_date', ''),
                'authors': meta.get('authors', ''),
                keyword_refs_header: kw_refs_for_stock,
                'Stock (FBal / FBtp / FBti) all references (all for stock)': all_refs_for_stock,
                '_is_functionally_valid': 1 if is_functionally_valid else 0,
                '_keyword_ref_count': kw_ref_count,
            }
            for column in REAGENT_BUCKET_COLUMNS:
                out_row[column] = _is_truthy(row.get(column, False))
            for gc in gpt_cols:
                cell_val = row.get(gc, '')
                out_row[gc] = _extract_row_pmid_entry(cell_val, pmid)
            rows.append(out_row)
    
    if not rows:
        return pd.DataFrame()
    
    out_df = pd.DataFrame(rows)
    
    # Group order: genes with more distinct Ref++ stocks first.
    gene_stock_counts = (
        out_df.groupby('gene')['stock #']
        .nunique()
        .sort_values(ascending=False)
        .to_dict()
    )
    out_df['_gene_stock_count'] = out_df['gene'].map(lambda g: gene_stock_counts.get(g, 0))
    
    # Parse publication date for deterministic within-stock PMID sorting.
    out_df['_pub_date_sort'] = pd.to_datetime(out_df['publication date'], errors='coerce')
    
    out_df = out_df.sort_values(
        by=[
            '_gene_stock_count',    # gene groups with most Ref++ stocks first
            'gene',
            '_is_functionally_valid',  # functionally valid stocks first
            '_keyword_ref_count',      # then more keyword refs first
            'stock #',                 # keep same stock adjacent
            '_pub_date_sort',
            'pmid',
        ],
        ascending=[False, True, False, False, True, False, True],
    ).reset_index(drop=True)
    
    # Final column order: core columns, then non-aggregate GPT columns.
    keyword_validation_col = None
    for gc in gpt_cols:
        if _is_keyword_validation_pmid_column(gc):
            keyword_validation_col = gc
            break

    ordered_cols = [
        'gene', 'gene synonyms', 'stock #', 'Balancers', 'matched_component_types',
        *REAGENT_BUCKET_COLUMNS,
        'pmid', 'title', 'journal',
        'publication date', 'authors',
    ]
    for gc in gpt_cols:
        if gc == keyword_validation_col:
            continue
        if gc not in ordered_cols and gc in out_df.columns:
            ordered_cols.append(gc)

    # Stock-level aggregate reference columns should appear last, in this order.
    aggregate_cols = [
        keyword_refs_header,
        'Stock (FBal / FBtp / FBti) all references (all for stock)',
    ]
    for col in aggregate_cols:
        if col in out_df.columns and col not in ordered_cols:
            ordered_cols.append(col)
    if keyword_validation_col and keyword_validation_col in out_df.columns:
        ordered_cols.append(keyword_validation_col)

    return out_df[ordered_cols]


def write_aggregated_excel(
    output_path: Path,
    source_workbook_path: Path,
    config: Dict[str, Any],
    combination_outputs: List[Tuple[List[str], pd.DataFrame, Dict]],
    all_stocks_df: Optional[pd.DataFrame],
    references_df: Optional[pd.DataFrame],
    unfiltered_references_df: Optional[pd.DataFrame] = None,
    verbose: bool = True,
    all_input_genes: Optional[Set[str]] = None,
    genes_no_stocks: Optional[Set[str]] = None,
    csv_input_genes: Optional[Set[str]] = None,
    n_input_genes: Optional[int] = None,
    gene_synonyms_map: Optional[Dict[str, str]] = None,
    soft_run: bool = False,
    flybase_data_path: Optional[Path] = None,
    pubmed_cache_path: Optional[Path] = None,
    pubmed_client: Optional[PubMedClient] = None,
    current_to_input_map: Optional[Dict[str, str]] = None,
    gene_to_datasets: Optional[Dict[str, Set[str]]] = None,
    pipeline_settings: Optional[Settings] = None,
) -> None:
    """
    Write aggregated Excel file with Contents, Sheet1..N, References.
    
    Args:
        output_path: Path to output Excel file
        source_workbook_path: Original Stage 1 workbook used to build the aggregate
        config: Configuration dictionary
        combination_outputs: List of (combination, limited_df, summary_dict) tuples
        references_df: References DataFrame (may be None)
        unfiltered_references_df: Full References sheet prior to output-sheet PMID trimming.
        verbose: Print progress
        all_input_genes: Set of gene symbols from stocks_df (genes that have
            matched stocks). Used to identify genes not in any split sheet.
        genes_no_stocks: Set of gene symbols from the original CSV input that
            had 0 matched stocks found by Pipeline 1. Displayed separately from
            genes that have stocks but were filtered out by combinations.
        csv_input_genes: Set of unique gene symbols from the original CSV input
            (augmented with FlyBase current symbols for filtering purposes).
        n_input_genes: True count of unique input genes from CSVs (before
            augmentation with remapped symbols). Used for display.
        gene_synonyms_map: Optional map of gene symbol -> semicolon-delimited synonyms.
        pubmed_cache_path: Path to PubMed cache CSV (for resolving reference titles).
        pubmed_client: PubMed client used for targeted phenotype-reference enrichment.
    """
    filters_config = config.get('filters', {})
    filter_descriptions = config.get('filterDescriptions', {})
    settings = config.get('settings', {})
    max_per_gene = settings.get('maxStocksPerGene')
    max_per_allele = settings.get('maxStocksPerAllele')
    
    # Build summary data, assigning sheet names based on stock presence
    summary_rows = []
    sheet_name_idx = 0
    for combo, limited_df, summary_dict in combination_outputs:
        has_stocks = len(limited_df) > 0 if limited_df is not None else False
        if has_stocks:
            sheet_name_idx += 1
            summary_dict["Sheet Name"] = f"Sheet{sheet_name_idx}"
        else:
            summary_dict["Sheet Name"] = "-"
        summary_rows.append(summary_dict)
    summary_df = pd.DataFrame(summary_rows)
    # Reorder columns to place Sheet Name last
    cols = summary_df.columns.tolist()
    if "Sheet Name" in cols:
        cols.remove("Sheet Name")
        cols.append("Sheet Name")
        summary_df = summary_df[cols]
    
    # Get used filter names
    used_filters = set()
    for combo, _, _ in combination_outputs:
        used_filters.update(combo)
    
    # Build filter definitions with human-readable descriptions
    filter_def_rows = []
    for name in sorted(used_filters):
        if name not in filters_config:
            continue
        spec = filters_config[name]
        meaning = filter_descriptions.get(name, "Filter applied (see config for technical details).")
        exact_str = filter_spec_to_string(spec)
        filter_def_rows.append((name, f"{meaning} [{exact_str}]"))
    
    # Determine which combinations have stocks
    combo_has_stocks = [
        (combo, limited_df, len(limited_df) > 0 if limited_df is not None else False)
        for combo, limited_df, _ in combination_outputs
    ]
    
    stock_sheet_by_gene_df = pd.DataFrame()
    stock_phenotype_sheet_df = pd.DataFrame()
    similarity_targets: List[PhenotypeSimilarityTarget] = []
    embedding_scorer: Optional[EmbeddingSimilarityScorer] = None
    stock_phenotype_sheet_counts = {
        'genes': 0,
        'stocks': 0,
        'references': 0,
    }
    if soft_run:
        similarity_targets, embedding_scorer = _build_phenotype_similarity_context(
            config,
            pipeline_settings,
            verbose=verbose,
        )
        stock_phenotype_sheet_df = _build_stock_phenotype_sheet(
            all_stocks_df if all_stocks_df is not None else pd.DataFrame(),
            flybase_data_path,
            references_df=references_df,
            unfiltered_references_df=unfiltered_references_df,
            pubmed_cache_path=pubmed_cache_path,
            pubmed_client=pubmed_client,
            similarity_targets=similarity_targets,
            embedding_scorer=embedding_scorer,
            verbose=verbose,
            gene_to_datasets=gene_to_datasets,
        )
        stock_phenotype_sheet_counts = _summarize_stock_phenotype_sheet(
            stock_phenotype_sheet_df
        )

        keywords_for_pheno = [
            kw.lower()
            for kw in settings.get('relevantSearchTerms', [])
            if kw
        ]
        if (
            keywords_for_pheno
            and len(stock_phenotype_sheet_df) > 0
            and 'Phenotype' in stock_phenotype_sheet_df.columns
            and _has_source_stock_columns(stock_phenotype_sheet_df)
        ):
            source_stock_series = _get_source_stock_series(stock_phenotype_sheet_df)
            pheno_lower = stock_phenotype_sheet_df['Phenotype'].fillna('').str.lower()
            kw_mask = pheno_lower.apply(
                lambda p: any(kw in p for kw in keywords_for_pheno)
            )
            matched_reagents = (
                source_stock_series.loc[kw_mask]
                .dropna()
                .loc[lambda s: s.str.strip().ne('')]
            )
            n_unique = matched_reagents.nunique()
            kw_display = ' or '.join(f"'{kw}'" for kw in keywords_for_pheno)
            print(
                f"\n    Phenotype keyword hits: {n_unique} unique reagent(s) "
                f"have {kw_display} in their Phenotype "
                f"(across {int(kw_mask.sum())} phenotype rows in the "
                f"Stock Phenotype Sheet)."
            )
            for kw in keywords_for_pheno:
                per_kw_mask = pheno_lower.str.contains(kw, na=False)
                per_kw_reagents = (
                    source_stock_series.loc[per_kw_mask]
                    .dropna()
                    .loc[lambda s: s.str.strip().ne('')]
                )
                print(
                    f"      - '{kw}': {per_kw_reagents.nunique()} unique reagent(s) "
                    f"across {int(per_kw_mask.sum())} rows"
                )
    else:
        stock_sheet_by_gene_df = _build_stock_sheet_by_gene(
            combination_outputs,
            references_df,
            csv_input_genes=csv_input_genes,
            gene_synonyms_map=gene_synonyms_map,
            current_to_input_map=current_to_input_map,
        )
    
    # Pre-compute gene symbols that appear in at least one sheet
    genes_in_sheets: Set[str] = set()
    for _combo, _ldf, _has in combo_has_stocks:
        if _ldf is not None and len(_ldf) > 0:
            if 'relevant_gene_symbols' in _ldf.columns:
                for syms in _ldf['relevant_gene_symbols'].dropna():
                    for s in str(syms).split(';'):
                        s = s.strip()
                        if s and s != '-':
                            genes_in_sheets.add(s)
            else:
                _g_col = _find_gene_id_column(_ldf)
                if _g_col:
                    for val in _ldf[_g_col].dropna().astype(str):
                        for v in val.split(';'):
                            v = v.strip()
                            if v and v != '-':
                                genes_in_sheets.add(v)
    
    # Write Excel
    with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
        workbook = writer.book
        
        # Formats
        fmt_13 = workbook.add_format({'font_size': 13, 'align': 'left'})
        fmt_13_wrap = workbook.add_format({'font_size': 13, 'align': 'left', 'text_wrap': True, 'valign': 'top'})
        fmt_13_bold = workbook.add_format({'font_size': 13, 'bold': True, 'align': 'left'})
        bold_bottom = workbook.add_format({'bold': True, 'bottom': 2, 'font_size': 13, 'align': 'left'})
        fmt_faint_bottom = workbook.add_format({'font_size': 13, 'bottom': 2, 'bottom_color': '#808080', 'align': 'left'})
        fmt_grey = workbook.add_format({'bg_color': '#D9D9D9'})
        fmt_dark_grey_white = workbook.add_format(
            {'bg_color': '#595959', 'font_color': '#FFFFFF'}
        )
        
        # Contents sheet
        contents_ws = workbook.add_worksheet("Contents")
        col_widths = defaultdict(int)
        
        def write_cell(r, c, val, fmt=None, skip_width=False):
            f = fmt if fmt else fmt_13
            if isinstance(val, bool):
                s = str(val)
                contents_ws.write_string(r, c, s, f)
                if not skip_width:
                    col_widths[c] = max(col_widths[c], len(s))
            elif isinstance(val, (int, float)) and not pd.isna(val):
                contents_ws.write_number(r, c, val, f)
                if not skip_width:
                    col_widths[c] = max(col_widths[c], len(str(val)))
            else:
                s = str(val) if val is not None and not (isinstance(val, float) and pd.isna(val)) else ""
                contents_ws.write_string(r, c, s, f)
                if not skip_width:
                    col_widths[c] = max(col_widths[c], len(s))
        
        def write_row(r, c, vals, fmt=None):
            for i, v in enumerate(vals):
                write_cell(r, c + i, v, fmt)
        
        row = 0
        
        # Section 1: Counts table
        write_cell(row, 0, "Stock / Allele / Gene counts per combination")
        row += 1
        
        # Format limits (treat very large numbers as unlimited)
        gene_limit_str = 'unlimited' if (max_per_gene is None or max_per_gene >= 1000000) else str(max_per_gene)
        allele_limit_str = 'unlimited' if (max_per_allele is None or max_per_allele >= 1000000) else str(max_per_allele)
        limits_str = f"Limits: up to {gene_limit_str} stocks per gene, up to {allele_limit_str} stocks per allele"
        write_cell(row, 0, limits_str)
        row += 1
        
        # Show gene counts (intersect with input genes to avoid counting hitchhiker genes)
        n_in_sheets = len(genes_in_sheets & csv_input_genes) if csv_input_genes else len(genes_in_sheets)
        n_no_stocks = len(genes_no_stocks) if genes_no_stocks else 0
        if csv_input_genes is not None:
            display_input_count = n_input_genes if n_input_genes else len(csv_input_genes)
            counts_str = (
                f"Total Input Genes (across all sets): {display_input_count}\n"
                f"Total Genes Included in below Categories: {n_in_sheets}\n"
                f"Input Genes with 0 matched stocks: {n_no_stocks}"
            )
        else:
            counts_str = (
                f"Total Input Genes (across all sets): unknown (no CSV gene lists found)\n"
                f"Total Genes Included in below Categories: {n_in_sheets}"
            )
        write_cell(row, 0, counts_str)
        row += 1
        
        # Show search terms used for Ref++ determination
        keywords_list = settings.get('relevantSearchTerms', [])
        if keywords_list:
            write_cell(row, 0, f"Search terms (for Ref++ status, case-insensitive): {', '.join(keywords_list)}")
            row += 1
            kw_ref_col_name = f"{' OR '.join(keywords_list)} references count"
            write_cell(row, 0, f"Keyword reference count column: \"{kw_ref_col_name}\"")
            row += 1
        
        if len(summary_df) > 0:
            write_row(row, 0, summary_df.columns.tolist(), bold_bottom)
            row += 1
            for i, (_, r) in enumerate(summary_df.iterrows()):
                row_fmt = fmt_faint_bottom if (i + 1) % 3 == 0 else fmt_13
                write_row(row, 0, ["" if pd.isna(x) else x for x in r.tolist()], row_fmt)
                row += 1
            
            # Total row with unique counts across all combinations
            all_stock_ids = set()
            all_gene_ids = set()
            all_allele_ids = set()
            for combo, limited_df, _ in combination_outputs:
                if limited_df is not None and len(limited_df) > 0:
                    s_col = _find_stock_id_column(limited_df)
                    # Prefer relevant_gene_symbols for gene counting (consistent
                    # with csv_input_genes which contains symbols, not FBgn IDs)
                    g_col = 'relevant_gene_symbols' if 'relevant_gene_symbols' in limited_df.columns else _find_gene_id_column(limited_df)
                    a_col = _find_allele_column(limited_df)
                    if s_col:
                        all_stock_ids.update(
                            v for v in limited_df[s_col].dropna().astype(str).unique()
                            if v and v != '-'
                        )
                    if g_col:
                        for cell in limited_df[g_col].dropna().astype(str).unique():
                            for v in cell.split(';'):
                                v = v.strip()
                                if v and v != '-':
                                    all_gene_ids.add(v)
                    if a_col:
                        # Only count alleles from rows whose gene(s) are in the input set
                        for _, _row in limited_df.iterrows():
                            if csv_input_genes is not None and g_col:
                                row_genes = set()
                                raw_g = str(_row.get(g_col, ''))
                                for g in raw_g.split(';'):
                                    g = g.strip()
                                    if g and g != '-' and g.lower() != 'nan':
                                        row_genes.add(g)
                                if not (row_genes & csv_input_genes):
                                    continue
                            raw_a = str(_row.get(a_col, ''))
                            for v in raw_a.split(';'):
                                v = v.strip()
                                if v and v != '-' and v.lower() != 'nan':
                                    all_allele_ids.add(v)
            # Filter gene IDs to only input genes (exclude hitchhiker genes)
            if csv_input_genes is not None:
                all_gene_ids &= csv_input_genes
            fmt_total = workbook.add_format({'bold': True, 'font_size': 13, 'top': 2, 'align': 'left'})
            total_row = []
            for col in summary_df.columns.tolist():
                if col == "Category":
                    total_row.append("Total (unique across sheets)")
                elif col == "# Stocks":
                    total_row.append(len(all_stock_ids))
                elif col == "# Alleles":
                    total_row.append(len(all_allele_ids))
                elif col == "# Genes":
                    total_row.append(len(all_gene_ids))
                else:
                    total_row.append("")
            write_row(row, 0, total_row, fmt_total)
            row += 1
        row += 2
        
        # Section 2: Filter definitions
        write_cell(row, 0, "What each filter means")
        row += 1
        
        if len(filter_def_rows) > 0:
            write_row(row, 0, ["Filter / Meaning"], bold_bottom)
            row += 1
            for name, meaning in filter_def_rows:
                combined_text = f"{name}: {meaning}"
                contents_ws.write_rich_string(
                    row, 0,
                    fmt_13_bold, name,
                    fmt_13, f": {meaning}",
                )
                col_widths[0] = max(col_widths[0], len(combined_text))
                row += 1
        row += 2
        
        # Section 3: Sheet / combination and gene lists
        write_cell(row, 0, "Unique gene symbols per defined sheet")
        row += 1
        
        sheet_index = 0
        wrote_any_sheet = False
        for combo, limited_df, has_stocks in combo_has_stocks:
            if not has_stocks or limited_df is None or len(limited_df) == 0:
                continue
            
            combo_label = " >> ".join(combo)
            sheet_index += 1
            label = f"Sheet{sheet_index}: {combo_label}"
            
            write_cell(row, 0, label, fmt_13_bold)
            row += 1
            wrote_any_sheet = True
            
            gene_col = _find_gene_id_column(limited_df)
            if gene_col:
                gene_ids = limited_df[gene_col].dropna().astype(str).unique()
                # Try to get gene symbols from relevant_gene_symbols column
                symbols = set()
                if 'relevant_gene_symbols' in limited_df.columns:
                    for syms in limited_df['relevant_gene_symbols'].dropna():
                        for s in str(syms).split(';'):
                            s = s.strip()
                            if s and s != '-':
                                symbols.add(s)
                if not symbols:
                    symbols = set(gene_ids)
                
                for sym in sorted(symbols):
                    write_cell(row, 0, sym)
                    row += 1
            else:
                write_cell(row, 0, "(no gene column found)")
                row += 1
            row += 1
        
        if not wrote_any_sheet:
            write_cell(row, 0, "(no sheets with stocks)")
            row += 2

        # Section 4: What the references-focused sheets include.
        sheet_label = "Stock Phenotype Sheet" if soft_run else "Stock Sheet by Gene"
        write_cell(row, 0, f"References and {sheet_label} inclusion criteria", fmt_13_bold)
        row += 1
        write_cell(
            row,
            0,
            "References sheet includes PMIDs cited by stocks that are present in output sheets (Sheet1..N).",
            fmt_13_wrap,
            skip_width=True,
        )
        row += 1
        if soft_run:
            write_cell(
                row,
                0,
                "Stock Phenotype Sheet includes unique (source/stock, genotype, reference) rows for stocks present in output sheets, based on shared gene-relevant stock-component IDs found in FlyBase genotype_phenotype_data. Genes, reagent type or allele symbols, Balancers, matched_component_types, one-hot reagent buckets, allele_class_terms, transgenic_product_class_terms, and phenotype terms are aggregated within each row.",
                fmt_13_wrap,
                skip_width=True,
            )
            row += 1
            write_cell(
                row,
                0,
                "Stock Phenotype Sheet totals: "
                f"{stock_phenotype_sheet_counts['genes']} unique genes, "
                f"{stock_phenotype_sheet_counts['stocks']} unique stocks, "
                f"{stock_phenotype_sheet_counts['references']} unique references.",
                fmt_13_wrap,
                skip_width=True,
            )
            row += 2
            column_definitions = _get_stock_phenotype_sheet_column_definitions(
                stock_phenotype_sheet_df
            )
            if column_definitions:
                write_cell(row, 0, "Stock Phenotype Sheet columns", fmt_13_bold)
                write_cell(
                    row,
                    1,
                    _get_reagent_bucket_one_hot_note(),
                    fmt_13_wrap,
                    skip_width=True,
                )
                row += 1
                write_row(row, 0, ["Column", "Definition"], bold_bottom)
                row += 1
                for column, definition in column_definitions:
                    write_cell(row, 0, column)
                    write_cell(row, 1, definition, fmt_13_wrap, skip_width=True)
                    row += 1
        else:
            write_cell(
                row,
                0,
                "Stock Sheet by Gene includes unique (stock, keyword-hit PMID) rows for stocks present in output sheets.",
                fmt_13_wrap,
                skip_width=True,
            )
        row += 2
        
        # Section 5a: Genes with stocks but not included in any split sheet
        if all_input_genes is not None:
            genes_not_in_sheets = sorted(all_input_genes - genes_in_sheets)
            
            write_cell(
                row, 0,
                f"Genes with matched stocks but not included in any split sheet "
                f"({len(genes_not_in_sheets)} of {len(all_input_genes)} genes with stocks)",
                fmt_13_bold,
            )
            row += 1
            
            if genes_not_in_sheets:
                _c2i = current_to_input_map or {}
                for sym in genes_not_in_sheets:
                    write_cell(row, 0, _c2i.get(sym, sym))
                    row += 1
            else:
                write_cell(row, 0, "(all genes with stocks appear in at least one sheet)")
                row += 1
            row += 1
        
        # Section 5b: Input genes with 0 matched stocks
        if csv_input_genes is not None:
            genes_no_stocks_sorted = sorted(genes_no_stocks) if genes_no_stocks else []
            write_cell(
                row, 0,
                f"Input Genes with 0 matched stocks "
                f"({len(genes_no_stocks_sorted)} of {n_input_genes or len(csv_input_genes)} input genes)",
                fmt_13_bold,
            )
            row += 1
            
            if genes_no_stocks_sorted:
                for sym in genes_no_stocks_sorted:
                    write_cell(row, 0, sym)
                    row += 1
            else:
                write_cell(row, 0, "(none)")
                row += 1
            row += 1
        
        # Auto-resize Contents columns
        for c, w in col_widths.items():
            contents_ws.set_column(c, c, min(w + 2, 255))
        if soft_run and len(stock_phenotype_sheet_df) > 0:
            contents_ws.set_column(1, 1, 96)
        
        # Data sheets (Sheet1, Sheet2, ...)
        sheet_index = 0
        for combo, limited_df, has_stocks in combo_has_stocks:
            if has_stocks and limited_df is not None and len(limited_df) > 0:
                sheet_index += 1
                sheet_name = f"Sheet{sheet_index}"
                # Prefix GPT-derived column headers with [EXPERIMENTAL]
                cleaned_df = _normalize_stock_sheet_columns(limited_df)
                limited_df_out = apply_experimental_prefix(cleaned_df)
                limited_df_out.to_excel(writer, sheet_name=sheet_name, index=False)
                # Apply light grey fill to GPT-derived column headers
                _apply_grey_fill_xlsxwriter(
                    writer.sheets[sheet_name], limited_df_out, fmt_grey
                )
        
        # References sheet
        if references_df is not None and len(references_df) > 0:
            references_df.to_excel(writer, sheet_name="References", index=False)
            kw_title_col = find_keyword_title_abstract_column(references_df)
            if kw_title_col and kw_title_col in references_df.columns:
                refs_ws = writer.sheets["References"]
                true_fmt = workbook.add_format({'bg_color': '#CFE8CF'})
                false_fmt = workbook.add_format({'bg_color': '#EBC7C7'})
                kw_col_idx = references_df.columns.get_loc(kw_title_col)
                for excel_row_idx, value in enumerate(references_df[kw_title_col].tolist(), start=1):
                    val_norm = str(value).strip().lower()
                    if value is True or val_norm in ("true", "1", "yes"):
                        refs_ws.write_boolean(excel_row_idx, kw_col_idx, True, true_fmt)
                    elif value is False or val_norm in ("false", "0", "no"):
                        refs_ws.write_boolean(excel_row_idx, kw_col_idx, False, false_fmt)
        
        if soft_run:
            if stock_phenotype_sheet_df is not None and len(stock_phenotype_sheet_df) > 0:
                _write_phenotype_similarity_sheet(
                    writer,
                    workbook,
                    "Stock Phenotype Sheet",
                    stock_phenotype_sheet_df,
                )
        else:
            if stock_sheet_by_gene_df is not None and len(stock_sheet_by_gene_df) > 0:
                stock_sheet_by_gene_out = apply_experimental_prefix(stock_sheet_by_gene_df)
                stock_sheet_by_gene_out.to_excel(
                    writer, sheet_name="Stock Sheet by Gene", index=False
                )
                stock_aggregate_columns: Set[str] = set()
                for col_name in stock_sheet_by_gene_out.columns:
                    col_name_str = str(col_name)
                    if (
                        col_name_str.startswith("Stock (")
                        and col_name_str.endswith("(all for stock)")
                    ) or _is_keyword_validation_pmid_column(col_name_str):
                        stock_aggregate_columns.add(col_name_str)
                _apply_grey_fill_xlsxwriter(
                    writer.sheets["Stock Sheet by Gene"],
                    stock_sheet_by_gene_out,
                    fmt_grey,
                    dark_header_format=fmt_dark_grey_white,
                    dark_header_columns=stock_aggregate_columns,
                )

    if soft_run and stock_phenotype_sheet_df is not None and len(stock_phenotype_sheet_df) > 0:
        if pipeline_settings is not None and pipeline_settings.enable_oai_embedding:
            _write_similarity_tier_workbook(
                output_path=output_path,
                source_workbook_path=source_workbook_path,
                stock_phenotype_sheet_df=stock_phenotype_sheet_df,
                combination_outputs=combination_outputs,
                csv_input_genes=csv_input_genes,
                simple_buckets=bool(pipeline_settings and pipeline_settings.simple_buckets),
                verbose=verbose,
                format_as_masterlist=True,
            )
        similarity_output_dir = output_path.parent / f"{output_path.stem}_similarity"
        written_visuals = plot_similarity_outputs(
            phenotype_sheet_df=stock_phenotype_sheet_df,
            targets=similarity_targets,
            output_dir=similarity_output_dir,
            embedding_scorer=embedding_scorer,
        )
        if verbose and written_visuals:
            print(
                f"    Wrote {len(written_visuals)} phenotype similarity plot(s) "
                f"to {similarity_output_dir}"
            )
    
    if verbose:
        print(f"    Saved: {output_path.name}")


###############################################################################
# Main Pipeline Class
###############################################################################

class StockSplittingPipeline:
    """
    Pipeline for organizing stocks using JSON-based filter configurations.
    
    This is Pipeline 2 in the fly-stocker workflow:
    Input: Excel files from Stage 1 (`find-stocks`) output
    Output: Aggregated Excel file with Contents, Sheet1..N, References sheets
    """
    
    def __init__(self, settings: Optional[Settings] = None):
        """
        Initialize the pipeline.
        
        Args:
            settings: Configuration settings (uses defaults if None)
        """
        self.settings = settings or Settings()
        self._pubmed_cache = PubMedCache(self.settings.pubmed_cache_path)
        self._pubmed_client = PubMedClient(
            cache=self._pubmed_cache,
            api_key=self.settings.ncbi_api_key,
        )
        self._fulltext_fetcher = FullTextFetcher(
            unpaywall_token=self.settings.unpaywall_token,
            ncbi_api_key=self.settings.ncbi_api_key,
            method_cache_path=self.settings.fulltext_method_cache_path,
        )
        self._validator = FunctionalValidator(
            openai_api_key=self.settings.openai_api_key,
            model=self.settings.openai_model,
            log_dir=self.settings.gpt_log_dir / datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            if self.settings.enable_gpt_logging
            else None,
        )
    
    def _find_genes_without_stocks(
        self,
        input_dir: Path,
        stocks_df: pd.DataFrame,
        verbose: bool = True,
    ) -> Tuple[Set[str], Optional[Set[str]], Dict[str, str], int, Dict[str, Set[str]]]:
        """
        Identify input genes that have 0 matched stocks.

        Uses FBgn IDs as ground truth so that renamed/synonym gene symbols
        do not cause false negatives.

        Returns
        -------
        genes_no_stocks : set[str]
            Gene symbols from the input that had zero matched stocks.
        csv_input_genes : set[str] | None
            All unique gene symbols from the input CSVs **augmented** with
            FlyBase current symbols for input genes, or ``None`` if no
            suitable CSVs were found.
        current_to_input_map : dict[str, str]
            Mapping of FlyBase current gene symbols to the user's original
            input symbols (via FBgn bridge).  Empty when unavailable.
        n_input_genes : int
            The true count of unique input genes (FBgn IDs) from the CSVs,
            before augmentation with remapped symbols.
        gene_to_datasets : dict[str, set[str]]
            Mapping of gene symbols to the input CSV filename stem(s) where
            they appeared.
        """
        # ── 1. Collect FBgn IDs + gene symbols from input CSVs ──────
        csv_dirs = [input_dir]
        if input_dir.parent != input_dir:
            csv_dirs.append(input_dir.parent)

        fbgn_to_symbol: Dict[str, str] = {}   # FBgn  → display symbol
        symbol_to_fbgn: Dict[str, str] = {}   # symbol → FBgn (reverse)
        fbgn_to_datasets: Dict[str, Set[str]] = defaultdict(set)

        for d in csv_dirs:
            for csv_path in sorted(d.glob("*.csv")):
                try:
                    df = pd.read_csv(csv_path, dtype=str)
                    if 'flybase_gene_id' not in df.columns:
                        continue
                    dataset_name = _format_input_dataset_name(csv_path)

                    sym_col = next(
                        (c for c in ('ext_gene', 'gene', 'gene_symbol',
                                     'Gene', 'Gene_Symbol')
                         if c in df.columns),
                        None,
                    )

                    for _, row in df.iterrows():
                        fbgn = str(row.get('flybase_gene_id', '')).strip()
                        if not fbgn.startswith('FBgn'):
                            continue
                        sym = fbgn
                        if sym_col:
                            s = str(row.get(sym_col, '')).strip()
                            if s and s.lower() != 'nan':
                                sym = s
                        fbgn_to_symbol[fbgn] = sym
                        symbol_to_fbgn[sym] = fbgn
                        if dataset_name:
                            fbgn_to_datasets[fbgn].add(dataset_name)
                except Exception:
                    continue

        if not fbgn_to_symbol:
            if verbose:
                print("    Note: No CSV gene lists with flybase_gene_id "
                      "found; cannot identify genes with 0 stocks")
            return set(), None, {}, 0, {}

        input_fbgns = set(fbgn_to_symbol.keys())
        csv_input_genes = set(fbgn_to_symbol.values())
        gene_to_datasets: Dict[str, Set[str]] = defaultdict(set)
        for fbgn, symbol in fbgn_to_symbol.items():
            if symbol:
                gene_to_datasets[symbol].update(fbgn_to_datasets.get(fbgn, set()))

        # ── 2. Determine stock FBgn IDs ──────────────────────────────
        if 'relevant_gene_symbols' not in stocks_df.columns:
            if verbose:
                print("    Warning: 'relevant_gene_symbols' column missing "
                      "from stocks data; cannot compute 0-stock genes")
            return set(), csv_input_genes, {}, len(input_fbgns), dict(gene_to_datasets)

        current_to_input_map: Dict[str, str] = {}
        stock_fbgns: Set[str] = set()

        has_fbgn_col = 'relevant_flybase_gene_ids' in stocks_df.columns

        if has_fbgn_col:
            # Preferred path: read FBgn IDs directly from the stocks sheet,
            # no symbol indirection needed.
            for _, row in stocks_df[['relevant_gene_symbols', 'relevant_flybase_gene_ids']].iterrows():
                symbols = [s.strip() for s in str(row.get('relevant_gene_symbols', '')).split(';') if s.strip() and s.strip().lower() != 'nan']
                fbgns = [f.strip() for f in str(row.get('relevant_flybase_gene_ids', '')).split(';') if f.strip() and f.strip().lower() != 'nan']
                stock_fbgns.update(f for f in fbgns if f.startswith('FBgn'))
                for sym, fbgn in zip(symbols, fbgns):
                    if not fbgn.startswith('FBgn'):
                        continue
                    input_sym = fbgn_to_symbol.get(fbgn)
                    if input_sym and input_sym != sym:
                        current_to_input_map[sym] = input_sym
                        csv_input_genes.add(sym)
                        gene_to_datasets[sym].update(fbgn_to_datasets.get(fbgn, set()))
        else:
            # Fallback for old Stage 1 outputs without the FBgn column.
            stock_symbols: Set[str] = set()
            for val in stocks_df['relevant_gene_symbols'].dropna().astype(str):
                for s in val.split(';'):
                    s = s.strip()
                    if s and s.lower() != 'nan':
                        stock_symbols.add(s)
            for sym in stock_symbols:
                fbgn = symbol_to_fbgn.get(sym)
                if fbgn:
                    stock_fbgns.add(fbgn)

        # ── 3. Compare the two FBgn sets ─────────────────────────────
        missing_fbgns = input_fbgns - stock_fbgns
        genes_no_stocks = {fbgn_to_symbol[fbgn] for fbgn in missing_fbgns}

        if verbose:
            print(f"    Input genes (from CSVs): {len(input_fbgns)}")
            print(f"    Genes covered in stocks: "
                  f"{len(input_fbgns - missing_fbgns)}")
            print(f"    Genes with 0 matched stocks: {len(genes_no_stocks)}")
            if current_to_input_map:
                print(f"    Gene symbols remapped (FlyBase current -> input): "
                      f"{len(current_to_input_map)}")

        return (
            genes_no_stocks,
            csv_input_genes,
            current_to_input_map,
            len(input_fbgns),
            dict(gene_to_datasets),
        )
    
    def _compute_derived_columns(
        self,
        df: pd.DataFrame,
        keywords: List[str],
        verbose: bool = True
    ) -> pd.DataFrame:
        """
        Compute derived columns needed for filtering.
        
        Adds: Balancers, multiple_insertions, ALLELE_PAPER_RELEVANCE_SCORE
        """
        if verbose:
            print("    Computing derived columns...")
        
        df = df.copy()
        
        # Compute Balancers
        df['num_Balancers'] = normalize_num_balancers_column(df)
        df['Balancers'] = compute_balancers_column(df)
        if verbose:
            has_balancers = (df['num_Balancers'] > 0).sum()
            print(f"      Balancers: {has_balancers} stocks with balancers")
        
        # Compute multiple_insertions
        df['multiple_insertions'] = compute_multiple_insertions_column(df)
        if verbose:
            multi_count = df['multiple_insertions'].sum()
            print(f"      Multiple insertions: {multi_count} stocks")
        
        # Compute ALLELE_PAPER_RELEVANCE_SCORE
        df['ALLELE_PAPER_RELEVANCE_SCORE'] = compute_allele_paper_relevance_score(df, keywords)
        if verbose:
            score_counts = df['ALLELE_PAPER_RELEVANCE_SCORE'].value_counts().to_dict()
            print(f"      Relevance scores: Ref++ (2)={score_counts.get(2, 0)}, Ref+ (1)={score_counts.get(1, 0)}, Ref- (0)={score_counts.get(0, 0)}")

        # Recompute one-hot reagent bucket columns from stock metadata so
        # downstream outputs stay consistent even for older Stage 1 workbooks.
        df = _recompute_reagent_bucket_columns(df)
        if verbose:
            bucket_counts = ", ".join(
                f"{column}={int(df[column].sum())}"
                for column in REAGENT_BUCKET_COLUMNS
                if column in df.columns
            )
            print(f"      Reagent buckets: {bucket_counts}")
        
        # Normalize PMID column for filtering
        pmid_col = None
        for col in ['PMID', 'pmid', 'PMIDs']:
            if col in df.columns:
                pmid_col = col
                break
        
        if pmid_col:
            # Replace empty/NaN with "-" for consistent filtering
            df[pmid_col] = df[pmid_col].fillna("-").astype(str)
            df.loc[df[pmid_col].str.strip() == "", pmid_col] = "-"
        
        return df
    
    def _load_stocks_from_excel(
        self,
        excel_path: Path,
        verbose: bool = True
    ) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], List[str]]:
        """Load stocks and references from an Excel file."""
        if verbose:
            print(f"\n  Loading: {excel_path.name}")
        
        # Read Stocks sheet
        stocks_df = pd.read_excel(excel_path, sheet_name='Stocks')
        
        # Read References sheet if it exists
        references_df = None
        try:
            references_df = pd.read_excel(excel_path, sheet_name='References')
        except ValueError:
            pass
        
        # Extract keywords from the title/abstract keyword-hit column
        keywords = []
        kw_ta_col = find_keyword_title_abstract_column(stocks_df)
        if kw_ta_col:
            keywords = extract_keywords_from_title_abstract_column(kw_ta_col)
        else:
            # Fallback: try the legacy GPT-derived column (for older files)
            kw_col = find_keyword_column(stocks_df)
            if kw_col:
                keywords = extract_keywords_from_column(kw_col)
        
        if verbose:
            print(f"    Loaded {len(stocks_df)} stocks")
            if references_df is not None:
                print(f"    Loaded {len(references_df)} references")
            if keywords:
                print(f"    Keywords: {', '.join(keywords)}")
        
        return stocks_df, references_df, keywords
    
    def _sort_stocks_for_priority(
        self,
        df: pd.DataFrame,
        keywords: List[str]
    ) -> pd.DataFrame:
        """
        Sort stocks for priority selection (best stocks first).
        
        Priority order:
        1. Higher ALLELE_PAPER_RELEVANCE_SCORE (2 > 1 > 0)
        2. More keyword-matching references (keyword_ref_count)
        3. More total references
        """
        df = df.copy()
        
        # pipeline_references renames keyword_ref_count to "<kw> references count"
        kw_count_col = None
        if 'keyword_ref_count' in df.columns:
            kw_count_col = 'keyword_ref_count'
        elif keywords:
            renamed_col = f"{' OR '.join(keywords)} references count"
            if renamed_col in df.columns:
                kw_count_col = renamed_col
        
        def get_sort_key(row):
            score = row.get('ALLELE_PAPER_RELEVANCE_SCORE', 0)
            
            # Count keyword refs (title/abstract hits from pipeline_references)
            kw_ref_count = 0
            if kw_count_col:
                try:
                    kw_ref_count = int(row.get(kw_count_col, 0))
                except (ValueError, TypeError):
                    kw_ref_count = 0
            
            # Count total refs
            total_refs = 0
            pmid_val = row.get('PMID', '')
            if pd.notna(pmid_val) and str(pmid_val).strip() and str(pmid_val).strip() != "-":
                total_refs = len([p for p in str(pmid_val).split(";") if p.strip()])
            
            return (-score, -kw_ref_count, -total_refs)
        
        df['_sort_key'] = df.apply(get_sort_key, axis=1)
        df = df.sort_values('_sort_key')
        df = df.drop(columns=['_sort_key'], errors='ignore')
        
        return df

    def _build_refpp_validation_tasks(
        self,
        combination_outputs: List[Tuple[List[str], pd.DataFrame, Dict]],
        references_df: Optional[pd.DataFrame],
        verbose: bool = True,
    ) -> Tuple[Set[str], Dict[str, List[Tuple[str, str]]]]:
        """
        Build (stock -> [(allele, pmid), ...]) tasks only from output sheets whose
        combinations include Ref++.
        """
        keyword_hit_pmids: Set[str] = set()
        if references_df is not None and len(references_df) > 0 and "PMID" in references_df.columns:
            kw_ref_col = find_keyword_title_abstract_column(references_df)
            if kw_ref_col:
                for _, r in references_df.iterrows():
                    pmid = clean_id(r.get("PMID", ""))
                    if not pmid:
                        continue
                    kw_val = r.get(kw_ref_col, False)
                    if kw_val is True or str(kw_val).strip().lower() in ("true", "1", "yes"):
                        keyword_hit_pmids.add(pmid)

        selected_stock_ids: Set[str] = set()
        stock_tasks: Dict[str, List[Tuple[str, str]]] = defaultdict(list)

        for combo, out_df, _summary in combination_outputs:
            if "Ref++" not in combo:
                continue
            if out_df is None or len(out_df) == 0:
                continue

            stock_col = _find_stock_id_column(out_df)
            allele_col = _find_allele_column(out_df)
            kw_pmids_col = _find_keyword_ref_pmids_column(out_df)
            pmid_col = None
            for c in ["PMID", "pmid", "PMIDs"]:
                if c in out_df.columns:
                    pmid_col = c
                    break

            if stock_col is None or allele_col is None:
                continue

            for _, row in out_df.iterrows():
                stock_num = clean_id(row.get(stock_col, ""))
                if not stock_num:
                    continue
                selected_stock_ids.add(stock_num)

                alleles = parse_semicolon_list(row.get(allele_col, ""))
                if not alleles:
                    continue

                pmids_hit: List[str] = []
                if kw_pmids_col and str(row.get(kw_pmids_col, "")).strip():
                    pmids_hit = [
                        clean_id(p) for p in str(row.get(kw_pmids_col, "")).split(";") if clean_id(p)
                    ]
                elif pmid_col and keyword_hit_pmids:
                    pmids_all = [clean_id(p) for p in str(row.get(pmid_col, "")).split(";") if clean_id(p)]
                    pmids_hit = [p for p in pmids_all if p in keyword_hit_pmids]

                if not pmids_hit:
                    continue

                for allele in alleles:
                    allele_clean = str(allele).strip()
                    if not allele_clean:
                        continue
                    for pmid in sorted(set(pmids_hit)):
                        stock_tasks[stock_num].append((allele_clean, pmid))

        # De-duplicate tasks per stock
        deduped_tasks: Dict[str, List[Tuple[str, str]]] = {}
        for stock_num, tasks in stock_tasks.items():
            seen: Set[Tuple[str, str]] = set()
            ordered: List[Tuple[str, str]] = []
            for allele, pmid in tasks:
                key = (allele, pmid)
                if key in seen:
                    continue
                seen.add(key)
                ordered.append(key)
            deduped_tasks[stock_num] = ordered

        if verbose:
            num_task_stocks = len([s for s, t in deduped_tasks.items() if t])
            num_tasks = sum(len(t) for t in deduped_tasks.values())
            print(f"\n    Ref++ output-sheet stocks selected for validation: {len(selected_stock_ids)}")
            print(f"    Ref++ stocks with keyword-hit validation tasks: {num_task_stocks}")
            print(f"    Ref++ (stock, allele, PMID) tasks: {num_tasks}")

        return selected_stock_ids, deduped_tasks

    def _run_refpp_functional_validation(
        self,
        stocks_df: pd.DataFrame,
        references_df: Optional[pd.DataFrame],
        combination_outputs: List[Tuple[List[str], pd.DataFrame, Dict]],
        keywords: List[str],
        verbose: bool = True,
    ) -> Tuple[pd.DataFrame, List[Tuple[List[str], pd.DataFrame, Dict]]]:
        """
        Run GPT validation only for stocks appearing in Ref++ output sheets.
        """
        selected_stock_ids, stock_tasks = self._build_refpp_validation_tasks(
            combination_outputs, references_df, verbose=verbose
        )
        if not selected_stock_ids:
            if verbose:
                print("    No Ref++ output-sheet stocks found; skipping GPT validation.")
            return stocks_df, combination_outputs

        stock_col = _find_stock_id_column(stocks_df)
        if stock_col is None:
            if verbose:
                print("    Warning: stock ID column not found in stocks_df; skipping GPT validation.")
            return stocks_df, combination_outputs

        stock_clean_series = stocks_df[stock_col].apply(clean_id)
        selected_mask = stock_clean_series.isin(selected_stock_ids)
        selected_stocks_df = stocks_df.loc[selected_mask].copy()

        if len(selected_stocks_df) == 0:
            if verbose:
                print("    No selected Ref++ stocks found in source stocks dataframe; skipping GPT validation.")
            return stocks_df, combination_outputs

        selected_stocks_df = run_functional_validation(
            stocks_df=selected_stocks_df,
            references_df=references_df,
            keywords=keywords,
            soft_run=self.settings.soft_run,
            validator=self._validator,
            fulltext_fetcher=self._fulltext_fetcher,
            pubmed_client=self._pubmed_client,
            pubmed_cache=self._pubmed_cache,
            max_gpt_calls_per_stock=self.settings.max_gpt_calls_per_stock,
            stock_tasks_override=stock_tasks,
        )

        if self.settings.soft_run:
            return stocks_df, combination_outputs

        gpt_cols = get_gpt_derived_columns(selected_stocks_df)
        if not gpt_cols:
            return stocks_df, combination_outputs

        validated_values = selected_stocks_df[[stock_col] + gpt_cols].copy()
        validated_values["_stock_clean"] = validated_values[stock_col].apply(clean_id)
        validated_values = validated_values.drop(columns=[stock_col]).drop_duplicates("_stock_clean")

        # Merge GPT results back onto full stocks dataframe (non-selected stocks remain empty).
        stocks_df = stocks_df.copy()
        stocks_df = stocks_df.drop(columns=[c for c in gpt_cols if c in stocks_df.columns], errors="ignore")
        stocks_df["_stock_clean"] = stocks_df[stock_col].apply(clean_id)
        stocks_df = stocks_df.merge(validated_values, on="_stock_clean", how="left")
        for col in gpt_cols:
            if col not in stocks_df.columns:
                stocks_df[col] = ""
            stocks_df[col] = stocks_df[col].fillna("")
        stocks_df = stocks_df.drop(columns=["_stock_clean"])

        # Merge GPT columns into each output sheet dataframe.
        merged_outputs: List[Tuple[List[str], pd.DataFrame, Dict]] = []
        for combo, out_df, summary in combination_outputs:
            if out_df is None or len(out_df) == 0:
                merged_outputs.append((combo, out_df, summary))
                continue

            out_stock_col = _find_stock_id_column(out_df)
            if out_stock_col is None:
                merged_outputs.append((combo, out_df, summary))
                continue

            out_df2 = out_df.copy()
            out_df2 = out_df2.drop(columns=[c for c in gpt_cols if c in out_df2.columns], errors="ignore")
            out_df2["_stock_clean"] = out_df2[out_stock_col].apply(clean_id)
            out_df2 = out_df2.merge(validated_values, on="_stock_clean", how="left")
            out_df2 = out_df2.drop(columns=["_stock_clean"])
            for col in gpt_cols:
                if col in out_df2.columns:
                    out_df2[col] = out_df2[col].fillna("")
            merged_outputs.append((combo, out_df2, summary))

        return stocks_df, merged_outputs
    
    def run(
        self,
        input_dir: Path,
        config_path: Optional[Path] = None,
        verbose: bool = True,
        run_validation: bool = False,
    ) -> Optional[Path]:
        """
        Run the stock splitting pipeline.
        
        Args:
            input_dir: Directory containing Excel files from Pipeline 1
            config_path: Path to JSON configuration file (uses default if None)
            verbose: Print progress information
            run_validation: If True, append GPT-derived validation columns to
                Ref++ output-sheet stocks. If False, perform split/filtering only.
        
        Returns:
            Path to output directory, or None if failed
        """
        input_dir = Path(input_dir)
        if not input_dir.is_dir():
            raise ValueError(f"Input path is not a directory: {input_dir}")
        
        # Load configuration
        if config_path is None:
            config_path = self.settings.split_config_path
        
        if not config_path.exists():
            raise ValueError(f"Config file not found: {config_path}")
        
        print(f"\n{'='*70}")
        if run_validation:
            print("STOCK VALIDATION PIPELINE (JSON Configuration)")
        else:
            print("STOCK SPLITTER PIPELINE (JSON Configuration)")
        print(f"{'='*70}")
        print(f"  Config: {config_path}")
        
        config = load_split_config(config_path)
        settings = config.get('settings', {})
        filters_config = config.get('filters', {})
        combinations = config.get('combinations', [])
        
        keywords = settings.get('relevantSearchTerms', [])
        max_per_gene = settings.get('maxStocksPerGene')
        max_per_allele = settings.get('maxStocksPerAllele')
        
        print(f"  Keywords: {', '.join(keywords) if keywords else '(none)'}")
        gene_limit_display = 'unlimited' if (max_per_gene is None or max_per_gene >= 1000000) else str(max_per_gene)
        allele_limit_display = 'unlimited' if (max_per_allele is None or max_per_allele >= 1000000) else str(max_per_allele)
        print(f"  Max per gene: {gene_limit_display}")
        print(f"  Max per allele: {allele_limit_display}")
        print(f"  Combinations: {len(combinations)}")

        def _filtered_excel_files(directory: Path) -> List[Path]:
            return [
                f for f in directory.glob("*.xlsx")
                if "Organized Stocks" not in str(f)
                and "Organized Stock Sheets" not in str(f)
                and "Uncategorized" not in str(f)
                and not f.name.startswith("~$")  # Excel temp files
                and not f.name.startswith(".")   # Hidden files
            ]

        stocks_dir = input_dir / "Stocks"
        direct_excel_files = _filtered_excel_files(input_dir)
        stocks_excel_files = _filtered_excel_files(stocks_dir) if stocks_dir.exists() else []

        # Prefer the canonical Stocks/ workspace when it contains Stage 1 workbooks.
        workbook_dir = stocks_dir if stocks_excel_files else input_dir
        excel_files = stocks_excel_files if stocks_excel_files else direct_excel_files

        if verbose and workbook_dir != input_dir:
            print(f"  Resolved workbook directory: {workbook_dir}")

        # Filter out files in output directories and temp files
        excel_files = [
            f for f in excel_files
            if "Organized Stocks" not in str(f)
            and "Organized Stock Sheets" not in str(f)
            and "Uncategorized" not in str(f)
            and not f.name.startswith("~$")  # Excel temp files
            and not f.name.startswith(".")   # Hidden files
        ]

        if not excel_files:
            print(f"No Excel files found in {input_dir}")
            return None

        print(f"\n  Found {len(excel_files)} Excel file(s)")

        # Create output directory under the resolved workbook directory.
        output_dir = workbook_dir / "Organized Stocks"
        output_dir.mkdir(exist_ok=True)

        # Carry forward no-PMID FBrf report from pipeline 1, if present.
        # Keep exactly one report file in the final Organized Stocks output.
        no_pmid_report_name = "references_without_pmid_fbrf.txt"
        no_pmid_sources = [
            workbook_dir / no_pmid_report_name,
            input_dir / no_pmid_report_name,
        ]
        existing_sources = []
        for path in no_pmid_sources:
            if path.exists() and path not in existing_sources:
                existing_sources.append(path)
        if existing_sources:
            dst = output_dir / no_pmid_report_name
            # Prefer the report that lives alongside the workbook(s) we resolved.
            preferred_src = workbook_dir / no_pmid_report_name
            src = preferred_src if preferred_src.exists() else existing_sources[0]
            shutil.move(str(src), str(dst))
            print(f"  Moved no-PMID FBrf report to: {dst}")
            
            # Remove any leftover duplicate source reports.
            for extra in existing_sources:
                if extra != src and extra.exists():
                    try:
                        extra.unlink()
                        print(f"  Removed duplicate no-PMID report: {extra}")
                    except Exception as e:
                        print(f"  Warning: Could not remove duplicate report {extra}: {e}")
        else:
            print("  No no-PMID FBrf report found to copy")

        gene_synonyms_map = _load_gene_synonyms_map(
            self.settings.flybase_data_path,
            verbose=verbose,
        )
        
        # Process each file
        for excel_path in excel_files:
            # Load data
            stocks_df, references_df, file_keywords = self._load_stocks_from_excel(excel_path, verbose)
            
            # Identify genes from original CSV input with 0 matched stocks
            (
                genes_no_stocks,
                csv_input_genes,
                current_to_input_map,
                n_input_genes,
                gene_to_datasets,
            ) = self._find_genes_without_stocks(workbook_dir, stocks_df, verbose)
            
            # Use keywords from config, fallback to file keywords
            active_keywords = keywords if keywords else file_keywords
            
            # Compute derived columns
            stocks_df = self._compute_derived_columns(stocks_df, active_keywords, verbose)
            
            # Sort for priority selection
            stocks_df = self._sort_stocks_for_priority(stocks_df, active_keywords)
            
            # Extract gene symbols from stocks_df (genes that have matched stocks)
            # Filter to input genes only so counts exclude hitchhiker genes
            all_input_genes: Set[str] = set()
            if 'relevant_gene_symbols' in stocks_df.columns:
                for syms in stocks_df['relevant_gene_symbols'].dropna():
                    for s in str(syms).split(';'):
                        s = s.strip()
                        if s and s != '-':
                            all_input_genes.add(s)
            else:
                _input_gene_col = _find_gene_id_column(stocks_df)
                if _input_gene_col:
                    for val in stocks_df[_input_gene_col].dropna().astype(str):
                        for v in val.split(';'):
                            v = v.strip()
                            if v and v != '-':
                                all_input_genes.add(v)
            if csv_input_genes is not None:
                all_input_genes &= csv_input_genes
            if verbose:
                print(f"    Genes with matched stocks: {len(all_input_genes)}")
            
            # Initialize counters for stock limits
            gene_counters: Dict[str, int] = {}
            allele_counters: Dict[str, int] = {}
            seen_stocks: Set[str] = set()
            
            # Process each combination
            combination_outputs = []
            
            # Debug: Show filter match counts for first combination
            if verbose and combinations:
                print(f"\n    Filter match summary (before combinations):")
                for filter_name, filter_spec in filters_config.items():
                    col = filter_spec.get('column', '')
                    if col in stocks_df.columns:
                        mask = apply_filter(stocks_df, filter_spec)
                        print(f"      {filter_name}: {mask.sum()}/{len(stocks_df)} stocks match")
                    else:
                        print(f"      {filter_name}: column '{col}' not found")
                print()
            
            for combo in combinations:
                # Apply filters
                filtered_df = apply_filter_combination(stocks_df, combo, filters_config)
                
                if verbose:
                    print(f"    Combination: {' >> '.join(combo)}")
                    print(f"      Matched: {len(filtered_df)} stocks")
                
                # Apply stock limits
                limited_df = apply_stock_limits(
                    filtered_df,
                    max_per_gene,
                    max_per_allele,
                    gene_counters,
                    allele_counters,
                    seen_stocks,
                    verbose
                )
                
                if verbose:
                    print(f"      After limits: {len(limited_df)} stocks")
                
                # Build summary
                stock_col = _find_stock_id_column(limited_df)
                gene_col = _find_gene_id_column(limited_df)
                allele_col = _find_allele_column(limited_df)
                
                num_stocks = limited_df[stock_col].nunique() if stock_col and len(limited_df) > 0 else 0
                # Split semicolon-separated values before counting unique genes/alleles
                # Prefer relevant_gene_symbols (symbols) so counts are consistent
                # with csv_input_genes; filter to input genes only.
                _gene_sym_col = 'relevant_gene_symbols' if 'relevant_gene_symbols' in limited_df.columns else gene_col
                if _gene_sym_col and len(limited_df) > 0:
                    _genes = set()
                    for cell in limited_df[_gene_sym_col].dropna().astype(str):
                        for v in cell.split(';'):
                            v = v.strip()
                            if v and v != '-':
                                _genes.add(v)
                    if csv_input_genes is not None:
                        _genes &= csv_input_genes
                    num_genes = len(_genes)
                else:
                    num_genes = 0
                if allele_col and len(limited_df) > 0:
                    _alleles = set()
                    # Only count alleles from rows whose gene(s) are in the input set
                    for _, _row in limited_df.iterrows():
                        if csv_input_genes is not None and _gene_sym_col:
                            row_genes = set()
                            raw = str(_row.get(_gene_sym_col, ''))
                            for g in raw.split(';'):
                                g = g.strip()
                                if g and g != '-' and g.lower() != 'nan':
                                    row_genes.add(g)
                            if not (row_genes & csv_input_genes):
                                continue
                        raw_a = str(_row.get(allele_col, ''))
                        for v in raw_a.split(';'):
                            v = v.strip()
                            if v and v != '-' and v.lower() != 'nan':
                                _alleles.add(v)
                    num_alleles = len(_alleles)
                else:
                    num_alleles = 0
                
                summary_dict = {
                    "Category": " >> ".join(combo),
                    "# Stocks": num_stocks,
                    "# Alleles": num_alleles,
                    "# Genes": num_genes,
                }
                
                # Remove internal columns before output
                internal_cols = ['multiple_insertions', 'ALLELE_PAPER_RELEVANCE_SCORE', '_stock_type', '_stock_collection', '_sort_key', '_priority_key']
                output_df = limited_df.drop(columns=[c for c in internal_cols if c in limited_df.columns], errors='ignore')
                
                combination_outputs.append((combo, output_df, summary_dict))

            if run_validation:
                # Run GPT validation only for stocks that actually made it into
                # Ref++ output sheets (post-filters and post-limits).
                stocks_df, combination_outputs = self._run_refpp_functional_validation(
                    stocks_df=stocks_df,
                    references_df=references_df,
                    combination_outputs=combination_outputs,
                    keywords=active_keywords,
                    verbose=verbose,
                )
            elif verbose:
                print("    Skipping GPT validation in split-stocks.")
            
            unfiltered_references_df = references_df.copy() if references_df is not None else None

            # Filter references to only those cited by stocks in output sheets
            if references_df is not None and len(references_df) > 0:
                included_pmids: Set[str] = set()
                for _combo, _out_df, _ in combination_outputs:
                    if _out_df is not None and len(_out_df) > 0:
                        _pmid_col = None
                        for col in ['PMID', 'pmid', 'PMIDs']:
                            if col in _out_df.columns:
                                _pmid_col = col
                                break
                        if _pmid_col:
                            for pmids_str in _out_df[_pmid_col].dropna().astype(str):
                                for p in pmids_str.split(';'):
                                    p = clean_id(p)
                                    if p and p != '-':
                                        included_pmids.add(p)
                
                if 'PMID' in references_df.columns and included_pmids:
                    original_ref_count = len(references_df)
                    references_df = references_df[
                        references_df['PMID'].apply(
                            lambda x: clean_id(x) in included_pmids
                            if pd.notna(x) else False
                        )
                    ].copy()
                    if verbose:
                        print(f"    References filtered: {original_ref_count} → {len(references_df)} "
                              f"(kept only references cited by stocks in output sheets)")
            
            # Create aggregated Excel
            output_name = f"{excel_path.stem}_aggregated.xlsx"
            output_path = output_dir / output_name
            
            write_aggregated_excel(
                output_path=output_path,
                source_workbook_path=excel_path,
                config=config,
                combination_outputs=combination_outputs,
                all_stocks_df=stocks_df,
                references_df=references_df,
                unfiltered_references_df=unfiltered_references_df,
                verbose=verbose,
                all_input_genes=all_input_genes,
                genes_no_stocks=genes_no_stocks,
                csv_input_genes=csv_input_genes,
                n_input_genes=n_input_genes,
                gene_synonyms_map=gene_synonyms_map,
                soft_run=self.settings.soft_run,
                flybase_data_path=self.settings.flybase_data_path,
                pubmed_cache_path=self.settings.pubmed_cache_path,
                pubmed_client=self._pubmed_client,
                current_to_input_map=current_to_input_map,
                gene_to_datasets=gene_to_datasets,
                pipeline_settings=self.settings,
            )
        
        # Print summary
        self._fulltext_fetcher.save_cache()
        self._fulltext_fetcher.clear_cache()

        print(f"\n{'='*70}")
        print("PROCESSING COMPLETE")
        print(f"  Processed: {len(excel_files)} file(s)")
        print(f"  Output directory: {output_dir}")
        
        return output_dir
