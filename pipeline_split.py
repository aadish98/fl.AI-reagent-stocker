"""
Pipeline 2: Split and Organize Stocks using JSON Configuration

This module provides the StockSplitterPipeline class that:
1. Loads stock data from Excel files (output of Pipeline 1)
2. Computes derived columns (Balancers via structural genotype parsing, multiple_insertions via component data, ALLELE_PAPER_RELEVANCE_SCORE via title/abstract keyword hits)
3. Applies configurable filters from JSON to create sheets
4. Creates aggregated Excel output with Contents, Sheet1..N, References sheets

Usage:
    from fly_stocker_v2 import StockSplitterPipeline
    from fly_stocker_v2.config import Settings
    
    settings = Settings()
    pipeline = StockSplitterPipeline(settings)
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

from .config import Settings, ValidationStatus, DEFAULT_SPLIT_CONFIG_PATH
from .external.pubmed import PubMedClient, PubMedCache
from .external.fulltext import FullTextFetcher, FunctionalValidator
from .utils import (
    clean_id,
    parse_semicolon_list,
    find_keyword_column,
    find_keyword_title_abstract_column,
    extract_keywords_from_column,
    extract_keywords_from_title_abstract_column,
    has_keyword_validation,
    generate_keyword_column_name,
    get_gpt_derived_columns,
    apply_experimental_prefix,
    EXPERIMENTAL_PREFIX,
)
from .validation_runner import run_functional_validation


# Human-readable filter type descriptions for Contents sheet
FILTER_TYPE_PHRASES = {
    "contains": "contains",
    "doesnt_contain": "doesn't contain",
    "equals": "equals",
    "doesnt_equal": "doesn't equal",
    "insertion_site": "has single insertion at",
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
    
    # Validate required keys
    required_keys = ['filters', 'combinations']
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Config missing required key: {key}")
    
    # Set defaults for optional settings
    if 'settings' not in config:
        config['settings'] = {}
    
    settings = config['settings']
    settings.setdefault('relevantSearchTerms', [])
    settings.setdefault('maxStocksPerGene', None)
    settings.setdefault('maxStocksPerAllele', None)
    
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

def load_balancers(balancers_path: Path) -> List[str]:
    """Load balancer names from CSV file."""
    if not balancers_path.exists():
        print(f"    Warning: Balancers file not found: {balancers_path}")
        return []
    
    df = pd.read_csv(balancers_path)
    if 'balancer_name' in df.columns:
        return df['balancer_name'].dropna().tolist()
    return []


def compute_balancers_column(df: pd.DataFrame, balancer_names: List[str]) -> pd.Series:
    """
    Compute the Balancers column from genotype using structural parsing.
    
    Identifies balancers by requiring a '/' in the genotype (indicating a
    balanced stock), then searching for balancer names on *both* sides of
    '/' (since Drosophila convention allows balancers on either side).
    
    Before matching, construct contents (P{...}, PBac{...}, Mi{...}, TI{...},
    M{...}) are masked out to prevent false positives from gene names embedded
    in construct symbols (e.g. 'hCMTM6' inside P{UAS-hCMTM6.HA}).
    
    Uses longest-match-first deduplication to avoid substring overlap
    (e.g. 'TM6B' also matching 'TM6' and 'M6').
    
    Returns comma-separated list of balancers found, or "-" if none.
    """
    # Sort balancer names longest-first so longer names are matched (and
    # consumed) before shorter substrings (e.g. TM6B before TM6 before M6)
    sorted_balancers = sorted(balancer_names, key=len, reverse=True)
    
    # Regex to mask construct brace contents only: P{...}, PBac{...}, etc.
    # We intentionally keep the post-'}' name (insertion site / chromosome)
    # because it can itself be a balancer (e.g. P{...}CyO means the construct
    # is inserted on a CyO balancer chromosome).
    _construct_re = re.compile(r'(?:P|PBac|Mi|TI|M)\{[^}]*\}')
    
    def find_balancers(genotype: str) -> str:
        if pd.isna(genotype) or not str(genotype).strip():
            return "-"
        genotype_str = str(genotype)
        
        # Only look for balancers in genotypes with '/' (indicates balanced stock)
        if '/' not in genotype_str:
            return "-"
        
        # Mask out construct contents to avoid false positives from gene names
        # embedded in construct symbols (e.g. P{UAS-hCMTM6.HA} matching TM6)
        masked = _construct_re.sub('', genotype_str)
        
        # Search the masked genotype for balancer names using longest-match-first
        found = set()
        temp = masked
        for bal in sorted_balancers:
            if bal in temp:
                found.add(bal)
                temp = temp.replace(bal, '', 1)
        
        return ", ".join(sorted(found)) if found else "-"
    
    genotype_col = None
    for col in ['genotype', 'stock_genotype', 'Genotype', 'FB_genotype']:
        if col in df.columns:
            genotype_col = col
            break
    
    if genotype_col is None:
        return pd.Series(["-"] * len(df), index=df.index)
    
    return df[genotype_col].apply(find_balancers)


def compute_multiple_insertions_column(df: pd.DataFrame) -> pd.Series:
    """
    Compute the multiple_insertions column using component data from stockgenes.
    
    Uses the 'all_stock_constructs' column (semicolon-separated construct symbols
    curated by BDSC) to count unique transgenic constructs per stock. If the
    component data is not available, falls back to counting distinct P-element
    constructs in the genotype string via regex.
    
    Returns True if the stock has more than one unique construct, False otherwise.
    """
    # Prefer component-based detection (from stockgenes.csv via pipeline_references)
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


def merge_stockcomps_data(
    df: pd.DataFrame,
    stockcomps_path: Optional[Path]
) -> pd.DataFrame:
    """
    Merge comment1 column from stockcomps_map_comments.csv for BDSC stocks.
    
    Only merges for stocks with collection containing 'BDSC' or 'Bloomington'.
    """
    if stockcomps_path is None or not stockcomps_path.exists():
        df['comment1'] = ""
        return df
    
    try:
        stockcomps = pd.read_csv(stockcomps_path, dtype=str)
        stockcomps = stockcomps.rename(columns={'Stk #': 'stk_num'})
        
        # Ensure we have the required columns
        if 'stk_num' not in stockcomps.columns or 'comment1' not in stockcomps.columns:
            df['comment1'] = ""
            return df
        
        # Find stock number column
        stock_num_col = None
        for col in ['stock_number', 'StockID', 'stock_id']:
            if col in df.columns:
                stock_num_col = col
                break
        
        if stock_num_col is None:
            df['comment1'] = ""
            return df
        
        # Prepare merge
        df = df.copy()
        df['_merge_stock_num'] = df[stock_num_col].astype(str)
        stockcomps['stk_num'] = stockcomps['stk_num'].astype(str)
        
        # Only merge for BDSC stocks
        collection_col = None
        for col in ['collection', 'Collection']:
            if col in df.columns:
                collection_col = col
                break
        
        if collection_col:
            is_bdsc = df[collection_col].str.contains('BDSC|Bloomington', case=False, na=False)
        else:
            is_bdsc = pd.Series([True] * len(df), index=df.index)
        
        # Perform merge
        merged = df.merge(
            stockcomps[['stk_num', 'comment1']],
            left_on='_merge_stock_num',
            right_on='stk_num',
            how='left',
            suffixes=('', '_stockcomps')
        )
        
        # Fill comment1 for non-BDSC stocks with empty string
        if 'comment1' not in merged.columns:
            merged['comment1'] = ""
        merged.loc[~is_bdsc, 'comment1'] = ""
        merged['comment1'] = merged['comment1'].fillna("")
        
        # Clean up merge columns
        merged = merged.drop(columns=['_merge_stock_num', 'stk_num'], errors='ignore')
        
        return merged
    
    except Exception as e:
        print(f"    Warning: Could not merge stockcomps data: {e}")
        df['comment1'] = ""
        return df


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
    worksheet, df: pd.DataFrame, grey_format
) -> None:
    """
    Apply light grey background to GPT-derived column **headers** via
    xlsxwriter.  Only the header row (row 0) is shaded.

    GPT-derived columns are detected by the '[EXPERIMENTAL] ' prefix that
    :func:`apply_experimental_prefix` adds to their header names.
    """
    for col_idx, col_name in enumerate(df.columns):
        if str(col_name).startswith(EXPERIMENTAL_PREFIX):
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
        if re.match(r"^Stock \(allele\) .+ references$", str(col)):
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

    pattern = re.compile(r"PMID:\s*([0-9]+)\s*:\s*(.*?)(?=(?:\s*\|\s*)?PMID:\s*[0-9]+\s*:|$)", re.IGNORECASE | re.DOTALL)
    matches = pattern.findall(text)
    if not matches:
        return cell_value

    for pmid_found, content in matches:
        if clean_id(pmid_found) == pmid_clean:
            entry = str(content).strip()
            return entry if entry else f"PMID: {pmid_clean}"
    return ""


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


def _build_stock_sheet_by_gene(
    combination_outputs: List[Tuple[List[str], pd.DataFrame, Dict]],
    references_df: Optional[pd.DataFrame],
    csv_input_genes: Optional[Set[str]] = None
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
                'stock #': stock_num,
                'pmid': pmid,
                'title': meta.get('title', ''),
                'journal': meta.get('journal', ''),
                'publication date': meta.get('publication_date', ''),
                'authors': meta.get('authors', ''),
                'Stock (allele) keyword-hit references (all for stock)': kw_refs_for_stock,
                'Stock (allele) all references (all for stock)': all_refs_for_stock,
                '_is_functionally_valid': 1 if is_functionally_valid else 0,
                '_keyword_ref_count': kw_ref_count,
            }
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
    
    # Final column order: core columns, then GPT columns.
    ordered_cols = [
        'gene', 'stock #', 'pmid', 'title', 'journal',
        'publication date', 'authors',
        'Stock (allele) keyword-hit references (all for stock)',
        'Stock (allele) all references (all for stock)',
    ]
    for gc in gpt_cols:
        if gc not in ordered_cols and gc in out_df.columns:
            ordered_cols.append(gc)
    
    return out_df[ordered_cols]


def write_aggregated_excel(
    output_path: Path,
    config: Dict[str, Any],
    combination_outputs: List[Tuple[List[str], pd.DataFrame, Dict]],
    references_df: Optional[pd.DataFrame],
    verbose: bool = True,
    all_input_genes: Optional[Set[str]] = None,
    genes_no_stocks: Optional[Set[str]] = None,
    csv_input_genes: Optional[Set[str]] = None,
) -> None:
    """
    Write aggregated Excel file with Contents, Sheet1..N, References.
    
    Args:
        output_path: Path to output Excel file
        config: Configuration dictionary
        combination_outputs: List of (combination, limited_df, summary_dict) tuples
        references_df: References DataFrame (may be None)
        verbose: Print progress
        all_input_genes: Set of gene symbols from stocks_df (genes that have
            BDSC stocks). Used to identify genes not in any split sheet.
        genes_no_stocks: Set of gene symbols from the original CSV input that
            had 0 BDSC stocks found by Pipeline 1. Displayed separately from
            genes that have stocks but were filtered out by combinations.
        csv_input_genes: Set of unique gene symbols from the original CSV input.
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
    
    stock_sheet_by_gene_df = _build_stock_sheet_by_gene(
        combination_outputs, references_df, csv_input_genes=csv_input_genes
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
        fmt_13_bold = workbook.add_format({'font_size': 13, 'bold': True, 'align': 'left'})
        bold_bottom = workbook.add_format({'bold': True, 'bottom': 2, 'font_size': 13, 'align': 'left'})
        fmt_faint_bottom = workbook.add_format({'font_size': 13, 'bottom': 2, 'bottom_color': '#808080', 'align': 'left'})
        fmt_grey = workbook.add_format({'bg_color': '#D9D9D9'})
        
        # Contents sheet
        contents_ws = workbook.add_worksheet("Contents")
        col_widths = defaultdict(int)
        
        def write_cell(r, c, val, fmt=None):
            f = fmt if fmt else fmt_13
            if isinstance(val, bool):
                s = str(val)
                contents_ws.write_string(r, c, s, f)
                col_widths[c] = max(col_widths[c], len(s))
            elif isinstance(val, (int, float)) and not pd.isna(val):
                contents_ws.write_number(r, c, val, f)
                col_widths[c] = max(col_widths[c], len(str(val)))
            else:
                s = str(val) if val is not None and not (isinstance(val, float) and pd.isna(val)) else ""
                contents_ws.write_string(r, c, s, f)
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
            n_input_csv = len(csv_input_genes)
            input_genes_str = str(n_input_csv)
            counts_str = (
                f"Total Input Genes (across all sets): {input_genes_str}\n"
                f"Total Genes Included in below Categories: {n_in_sheets}\n"
                f"Input Genes with 0 BDSC stocks: {n_no_stocks}"
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
        write_cell(row, 0, "References and Stock Sheet by Gene inclusion criteria", fmt_13_bold)
        row += 1
        write_cell(
            row,
            0,
            "References sheet includes PMIDs cited by stocks that are present in output sheets (Sheet1..N).",
        )
        row += 1
        write_cell(
            row,
            0,
            "Stock Sheet by Gene includes unique (stock, keyword-hit PMID) rows for stocks present in output sheets.",
        )
        row += 2
        
        # Section 5a: Genes with stocks but not included in any split sheet
        if all_input_genes is not None:
            genes_not_in_sheets = sorted(all_input_genes - genes_in_sheets)
            
            write_cell(
                row, 0,
                f"Genes with BDSC stocks but not included in any split sheet "
                f"({len(genes_not_in_sheets)} of {len(all_input_genes)} genes with stocks)",
                fmt_13_bold,
            )
            row += 1
            
            if genes_not_in_sheets:
                for sym in genes_not_in_sheets:
                    write_cell(row, 0, sym)
                    row += 1
            else:
                write_cell(row, 0, "(all genes with stocks appear in at least one sheet)")
                row += 1
            row += 1
        
        # Section 5b: Input genes with 0 BDSC stocks
        if csv_input_genes is not None:
            genes_no_stocks_sorted = sorted(genes_no_stocks) if genes_no_stocks else []
            write_cell(
                row, 0,
                f"Input Genes with 0 BDSC stocks "
                f"({len(genes_no_stocks_sorted)} of {len(csv_input_genes)} input genes)",
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
        
        # Stock Sheet by Gene
        if stock_sheet_by_gene_df is not None and len(stock_sheet_by_gene_df) > 0:
            stock_sheet_by_gene_out = apply_experimental_prefix(stock_sheet_by_gene_df)
            stock_sheet_by_gene_out.to_excel(
                writer, sheet_name="Stock Sheet by Gene", index=False
            )
            _apply_grey_fill_xlsxwriter(
                writer.sheets["Stock Sheet by Gene"],
                stock_sheet_by_gene_out,
                fmt_grey,
            )
    
    if verbose:
        print(f"    Saved: {output_path.name}")


###############################################################################
# Main Pipeline Class
###############################################################################

class StockSplitterPipeline:
    """
    Pipeline for organizing stocks using JSON-based filter configurations.
    
    This is Pipeline 2 in the fly-stocker workflow:
    Input: Excel files from Pipeline 1 (get-allele-refs output)
    Output: Aggregated Excel file with Contents, Sheet1..N, References sheets
    """
    
    def __init__(self, settings: Optional[Settings] = None):
        """
        Initialize the pipeline.
        
        Args:
            settings: Configuration settings (uses defaults if None)
        """
        self.settings = settings or Settings()
        self._balancers: Optional[List[str]] = None
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
    
    def _load_balancers(self) -> List[str]:
        """Load balancer names (cached)."""
        if self._balancers is None:
            self._balancers = load_balancers(self.settings.bdsc_balancers_path)
            print(f"    Loaded {len(self._balancers)} balancer names")
        return self._balancers
    
    def _find_genes_without_stocks(
        self,
        input_dir: Path,
        stocks_df: pd.DataFrame,
        verbose: bool = True,
    ) -> Tuple[Set[str], Optional[Set[str]]]:
        """
        Identify input genes that have 0 BDSC stocks.

        Algorithm
        ---------
        1. Read every input CSV's ``flybase_gene_id`` column and deduplicate
           to get the full set of input FBgn IDs.
        2. Read the stocks sheet's ``relevant_gene_symbols`` column (genes
           are semicolon-separated), deduplicate, and convert back to FBgn
           IDs via the input-CSV symbol→FBgn mapping.
        3. Missing = input FBgn set − stock FBgn set.

        Returns
        -------
        genes_no_stocks : set[str]
            Gene symbols from the input that had zero BDSC stocks.
        csv_input_genes : set[str] | None
            All unique gene symbols from the input CSVs, or ``None`` if no
            suitable CSVs were found.
        """
        # ── 1. Collect FBgn IDs + gene symbols from input CSVs ──────
        csv_dirs = [input_dir]
        if input_dir.parent != input_dir:
            csv_dirs.append(input_dir.parent)

        fbgn_to_symbol: Dict[str, str] = {}   # FBgn  → display symbol
        symbol_to_fbgn: Dict[str, str] = {}   # symbol → FBgn (reverse)

        for d in csv_dirs:
            for csv_path in sorted(d.glob("*.csv")):
                try:
                    df = pd.read_csv(csv_path, dtype=str)
                    if 'flybase_gene_id' not in df.columns:
                        continue

                    # Locate a gene-symbol column for display / mapping
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
                        sym = fbgn                       # fallback display name
                        if sym_col:
                            s = str(row.get(sym_col, '')).strip()
                            if s and s.lower() != 'nan':
                                sym = s
                        fbgn_to_symbol[fbgn] = sym
                        symbol_to_fbgn[sym] = fbgn
                except Exception:
                    continue

        if not fbgn_to_symbol:
            if verbose:
                print("    Note: No CSV gene lists with flybase_gene_id "
                      "found; cannot identify genes with 0 stocks")
            return set(), None

        input_fbgns = set(fbgn_to_symbol.keys())
        csv_input_genes = set(fbgn_to_symbol.values())

        # ── 2. Read relevant_gene_symbols from stocks ────────────────
        if 'relevant_gene_symbols' not in stocks_df.columns:
            if verbose:
                print("    Warning: 'relevant_gene_symbols' column missing "
                      "from stocks data; cannot compute 0-stock genes")
            return set(), csv_input_genes

        stock_symbols: Set[str] = set()
        for val in stocks_df['relevant_gene_symbols'].dropna().astype(str):
            for s in val.split(';'):
                s = s.strip()
                if s and s.lower() != 'nan':
                    stock_symbols.add(s)

        # Convert stock gene symbols → FBgn IDs via input-CSV mapping
        stock_fbgns: Set[str] = set()
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
            print(f"    Genes with 0 BDSC stocks: {len(genes_no_stocks)}")

        return genes_no_stocks, csv_input_genes
    
    def _compute_derived_columns(
        self,
        df: pd.DataFrame,
        keywords: List[str],
        verbose: bool = True
    ) -> pd.DataFrame:
        """
        Compute derived columns needed for filtering.
        
        Adds: Balancers, multiple_insertions, ALLELE_PAPER_RELEVANCE_SCORE, comment1
        """
        if verbose:
            print("    Computing derived columns...")
        
        df = df.copy()
        
        # Compute Balancers
        balancers = self._load_balancers()
        df['Balancers'] = compute_balancers_column(df, balancers)
        if verbose:
            has_balancers = (df['Balancers'] != "-").sum()
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
        
        # Merge comment1 from stockcomps
        df = merge_stockcomps_data(df, self.settings.bdsc_stockcomps_path)
        if verbose:
            has_comment = (df['comment1'].astype(str).str.strip() != "").sum()
            print(f"      comment1: {has_comment} stocks with comments")
        
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
        verbose: bool = True
    ) -> Optional[Path]:
        """
        Run the stock splitting pipeline.
        
        Args:
            input_dir: Directory containing Excel files from Pipeline 1
            config_path: Path to JSON configuration file (uses default if None)
            verbose: Print progress information
        
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
        
        # Find Excel files
        excel_files = list(input_dir.glob("*.xlsx"))
        
        # Also check Stocks subdirectory
        stocks_dir = input_dir / "Stocks"
        if stocks_dir.exists():
            excel_files.extend(stocks_dir.glob("*.xlsx"))
        
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
        
        # Create output directory
        output_dir = input_dir / "Organized Stocks"
        output_dir.mkdir(exist_ok=True)
        
        # Carry forward no-PMID FBrf report from pipeline 1, if present.
        # Keep exactly one report file in the final Organized Stocks output.
        no_pmid_report_name = "references_without_pmid_fbrf.txt"
        no_pmid_sources = [
            input_dir / no_pmid_report_name,
            stocks_dir / no_pmid_report_name,
        ]
        existing_sources = [p for p in no_pmid_sources if p.exists()]
        if existing_sources:
            dst = output_dir / no_pmid_report_name
            # Prefer Stocks/ report when both exist.
            preferred_src = stocks_dir / no_pmid_report_name
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
        
        # Process each file
        for excel_path in excel_files:
            # Load data
            stocks_df, references_df, file_keywords = self._load_stocks_from_excel(excel_path, verbose)
            
            # Identify genes from original CSV input with 0 BDSC stocks
            # IMPORTANT: Must be called BEFORE _compute_derived_columns which may overwrite comment1
            genes_no_stocks, csv_input_genes = self._find_genes_without_stocks(
                input_dir, stocks_df, verbose
            )
            
            # Use keywords from config, fallback to file keywords
            active_keywords = keywords if keywords else file_keywords
            
            # Compute derived columns
            stocks_df = self._compute_derived_columns(stocks_df, active_keywords, verbose)
            
            # Sort for priority selection
            stocks_df = self._sort_stocks_for_priority(stocks_df, active_keywords)
            
            # Extract gene symbols from stocks_df (genes that have BDSC stocks)
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
                print(f"    Genes with BDSC stocks: {len(all_input_genes)}")
            
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
                internal_cols = ['Balancers', 'multiple_insertions', 'ALLELE_PAPER_RELEVANCE_SCORE', 'comment1', '_stock_type', '_stock_collection', '_sort_key', '_priority_key']
                output_df = limited_df.drop(columns=[c for c in internal_cols if c in limited_df.columns], errors='ignore')
                
                combination_outputs.append((combo, output_df, summary_dict))

            # Run GPT validation only for stocks that actually made it into
            # Ref++ output sheets (post-filters and post-limits).
            stocks_df, combination_outputs = self._run_refpp_functional_validation(
                stocks_df=stocks_df,
                references_df=references_df,
                combination_outputs=combination_outputs,
                keywords=active_keywords,
                verbose=verbose,
            )
            
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
                                    p = p.strip()
                                    if p and p != '-' and p.lower() != 'nan':
                                        included_pmids.add(p)
                
                if 'PMID' in references_df.columns and included_pmids:
                    original_ref_count = len(references_df)
                    references_df = references_df[
                        references_df['PMID'].apply(
                            lambda x: str(x).strip() in included_pmids
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
                output_path,
                config,
                combination_outputs,
                references_df,
                verbose,
                all_input_genes=all_input_genes,
                genes_no_stocks=genes_no_stocks,
                csv_input_genes=csv_input_genes,
            )
        
        # Print summary
        self._fulltext_fetcher.save_cache()
        self._fulltext_fetcher.clear_cache()

        print(f"\n{'='*70}")
        print("PROCESSING COMPLETE")
        print(f"  Processed: {len(excel_files)} file(s)")
        print(f"  Output directory: {output_dir}")
        
        return output_dir
