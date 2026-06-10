#!/usr/bin/env python3
"""
Map input gene symbols to FlyBase FBgn IDs in CSV files.

Usage:
    python scripts/fetch_fbgn_ids.py <input_directory> [input_gene_col]
"""

import argparse
import gzip
import os
import sys
from glob import glob
from itertools import combinations
from pathlib import Path

import pandas as pd

symbol_to_name = {
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",  # Lowercase delta
    "ε": "epsilon",
    "ζ": "zeta",
    "η": "eta",
    "θ": "theta",
    "ι": "iota",
    "κ": "kappa",
    "λ": "lambda",
    "μ": "mu",
    "ν": "nu",
    "ξ": "xi",
    "ο": "omicron",
    "π": "pi",
    "ρ": "rho",
    "σ": "sigma",
    "τ": "tau",
    "υ": "upsilon",
    "φ": "phi",
    "χ": "chi",
    "ψ": "psi",
    "ω": "omega",
    "Α": "Alpha",  # Uppercase Alpha
    "Β": "Beta",
    "Γ": "Gamma",
    "Δ": "Delta",  # Uppercase Delta
    "Ε": "Epsilon",
    "Ζ": "Zeta",
    "Η": "Eta",
    "Θ": "Theta",
    "Ι": "Iota",
    "Κ": "Kappa",
    "Λ": "Lambda",
    "Μ": "Mu",
    "Ν": "Nu",
    "Ξ": "Xi",
    "Ο": "Omicron",
    "Π": "Pi",
    "Ρ": "Rho",
    "Σ": "Sigma",
    "Τ": "Tau",
    "Υ": "Upsilon",
    "Φ": "Phi",
    "Χ": "Chi",
    "Ψ": "Psi",
    "Ω": "Omega",
}


_SCRIPT_DIR = Path(__file__).parent.resolve()
_REPO_ROOT = _SCRIPT_DIR.parent
DEFAULT_FLYBASE_DATA = _REPO_ROOT / "data" / "flybase"


def replace_symbol(gene_series):
    for symbol, name in symbol_to_name.items():
        gene_series = gene_series.str.replace(symbol, name, regex=False)
    return gene_series


def create_expanded_mappings(mappings_df):
    """
    Create expanded mappings for each relevant column to ensure all entries are processed
    (e.g., pipe-separated values in 'fullname_synonym(s)' and 'symbol_synonym(s)').
    """

    def expand_synonyms(column_name):
        """
        Expand a column containing pipe-separated synonyms into individual mappings.
        """
        if column_name not in mappings_df.columns:
            raise ValueError(f"Column '{column_name}' not found in the DataFrame.")

        # Select relevant rows where the column is not null
        relevant_rows = mappings_df.dropna(subset=[column_name]).copy()

        # Split pipe-separated values into lists
        relevant_rows[column_name] = relevant_rows[column_name].str.split("|")

        # Explode to create individual rows for each synonym
        expanded = relevant_rows.explode(column_name).rename(columns={column_name: "synonym"})
        expanded["synonym"] = expanded["synonym"].str.strip()  # Remove whitespace

        return expanded[["synonym", "primary_FBid"]]

    # Create mappings for each relevant column
    columns_to_expand = [
        "current_symbol",
        "current_fullname",
        "fullname_synonym(s)",
        "symbol_synonym(s)",
    ]
    all_mappings = pd.DataFrame()

    for column in columns_to_expand:
        print(f"Processing column: {column}")
        expanded_mapping = expand_synonyms(column)
        all_mappings = pd.concat([all_mappings, expanded_mapping], ignore_index=True)

    # Drop duplicates to ensure unique mappings
    all_mappings = all_mappings.drop_duplicates()

    # Convert to a dictionary for efficient lookups
    synonym_to_fbgnid_map = dict(zip(all_mappings["synonym"], all_mappings["primary_FBid"]))

    return synonym_to_fbgnid_map


def map_gene_ids(df, gene_to_fbgnid_main, gene_to_fbgnid_synonym, gene_col, corrected_col):
    """
    Efficiently maps ext_gene in df to flybase_gene_id using main and synonym mappings,
    while applying prioritized step combinations in a vectorized manner.

    `corrected_col` records the exact (possibly spelling-/case-normalized) form that
    produced each successful match, leaving the caller's original column untouched.
    """
    import re

    def clean_gene_vectorized(series, steps):
        """
        Apply a sequence of cleaning steps to a pandas Series in a vectorized manner.
        """
        for step in steps:
            if step == "lowercase":
                series = series.str.lower()
            elif step == "remove_zeros":
                series = (
                    series.str.rstrip("0")
                    .str.replace("-0", "-", regex=False)
                    .str.replace("-00", "-", regex=False)
                )
            elif step == "capitalize":
                series = series.str.capitalize()
            elif step == "remove_hyphen":
                series = series.str.replace("-", "", regex=False)
            elif step == "number_to_end":
                series = series.str.replace(
                    r"^(\d+)-(\D+)$", r"\2\1", regex=True
                )  # Handles cases like "1-Sep" -> "Sep1"
            elif step == "number_to_start":
                series = series.str.replace(
                    r"^(\D+)-(\d+)$", r"\2\1", regex=True
                )  # Handles cases like "Sep-1" -> "1Sep"
            elif step == "add_cr":
                series = "CR" + series
            elif step == "uppercase":
                series = series.str.upper()
        return series

    # Step 1: Initial mapping. The corrected form starts as the symbol-normalized,
    # stripped input and is overwritten whenever a cleaning combination resolves a gene.
    df[gene_col] = df[gene_col].str.strip()
    df[corrected_col] = df[gene_col]
    df["flybase_gene_id"] = df[gene_col].map(gene_to_fbgnid_main)

    # Step 2: Synonym mapping for initially unmapped genes
    unmapped_genes_mask = df["flybase_gene_id"].isna()
    print(f"TOTAL GENES: {df.shape[0]}")
    print(f"Initially missed genes: {unmapped_genes_mask.sum()}")
    df.loc[unmapped_genes_mask, "flybase_gene_id"] = df.loc[unmapped_genes_mask, gene_col].map(
        gene_to_fbgnid_synonym
    )
    print(f'Unmatched synonyms: {df["flybase_gene_id"].isna().sum()}')
    print(f'{df[df["flybase_gene_id"].isna()]}')
    # Step 3: Prioritized step combinations
    steps = ["lowercase", "remove_zeros", "remove_hyphen", "number_to_end", "uppercase"]  # 'add_cr'
    all_step_combinations = [combo for r in range(1, len(steps) + 1) for combo in combinations(steps, r)]
    # for r in range(1, len(steps) + 1):
    #     all_step_combinations.extend(permutations(steps, r))

    prev_unmapped_count = df["flybase_gene_id"].isna().sum()

    # Step 4: Apply combinations iteratively
    for step_combo in all_step_combinations:
        unmapped_idx = df.index[df["flybase_gene_id"].isna()]
        if len(unmapped_idx) == 0:
            break  # Exit if all genes are mapped

        # Clean only the still-unmapped genes for this combination.
        cleaned_genes = clean_gene_vectorized(df.loc[unmapped_idx, gene_col], step_combo)

        # Try the main symbol mapping, recording the corrected spelling on hits.
        main_hits = cleaned_genes.map(gene_to_fbgnid_main).dropna()
        df.loc[main_hits.index, "flybase_gene_id"] = main_hits
        df.loc[main_hits.index, corrected_col] = cleaned_genes.loc[main_hits.index]

        # Fall back to the synonym mapping for whatever is still unmapped.
        remaining_idx = unmapped_idx.difference(main_hits.index)
        syn_hits = cleaned_genes.loc[remaining_idx].map(gene_to_fbgnid_synonym).dropna()
        df.loc[syn_hits.index, "flybase_gene_id"] = syn_hits
        df.loc[syn_hits.index, corrected_col] = cleaned_genes.loc[syn_hits.index]

        current_unmapped_count = df["flybase_gene_id"].isna().sum()
        if current_unmapped_count == prev_unmapped_count:
            continue  # Skip redundant combinations
        else:
            print(
                f"{step_combo} was succesful in resolving {prev_unmapped_count - current_unmapped_count} genes"
            )
        prev_unmapped_count = current_unmapped_count

    # Step 5: Log remaining unmapped genes
    remaining_unmapped_mask = df["flybase_gene_id"].isna()
    print(f"Remaining unmapped genes after all steps: {remaining_unmapped_mask.sum()}")
    if remaining_unmapped_mask.sum() > 0:
        print("Remaining unmapped genes:")
        print(df.loc[remaining_unmapped_mask, gene_col])

    return df


def find_latest_tsv(directory: Path, pattern: str) -> Path:
    """Find latest FlyBase TSV, preferring .tsv.gz over .tsv."""
    gz_files = sorted(glob(str(directory / f"{pattern}*.tsv.gz")), reverse=True)
    if gz_files:
        return Path(gz_files[0])

    tsv_files = sorted(glob(str(directory / f"{pattern}*.tsv")), reverse=True)
    if tsv_files:
        return Path(tsv_files[0])

    raise FileNotFoundError(f"No TSV file matching '{pattern}' found in {directory}")


def load_flybase_tsv(filepath, **kwargs) -> pd.DataFrame:
    """
    Load a FlyBase TSV file, handling metadata/comments and header variants.
    Supports both .tsv and .tsv.gz files.
    """
    filepath = Path(filepath)

    if filepath.suffix == ".gz" or str(filepath).endswith(".tsv.gz"):
        opener = lambda f: gzip.open(f, "rt", encoding="utf-8")
    else:
        opener = lambda f: open(f, "r", encoding="utf-8")

    skip_rows = 0
    header_line = None
    use_custom_header = False
    header_found = False

    with opener(filepath) as f:
        for i, line in enumerate(f):
            stripped = line.strip()

            if not stripped:
                skip_rows = i + 1
                continue

            if header_found:
                if stripped.startswith("#"):
                    skip_rows = i + 1
                    continue
                skip_rows = i
                break

            if stripped.startswith("#"):
                if "\t" not in stripped:
                    skip_rows = i + 1
                    continue

                # Some FlyBase files use "##col1\tcol2..." headers.
                header_content = stripped.lstrip("#")
                cols = header_content.split("\t")
                col_names = [c.strip() for c in cols if c.strip()]

                if len(col_names) >= 2:
                    header_line = col_names
                    use_custom_header = True
                    header_found = True
                    skip_rows = i + 1
                    continue
                else:
                    skip_rows = i + 1
                    continue
            else:
                cols = stripped.split("\t")
                first_field = cols[0].strip() if cols else ""
                if first_field.startswith("FB") and len(first_field) > 4:
                    break
                else:
                    break

    if use_custom_header and header_line:
        df = pd.read_csv(
            filepath,
            sep="\t",
            skiprows=skip_rows,
            names=header_line,
            low_memory=False,
            on_bad_lines="warn",
            **kwargs,
        )
    else:
        df = pd.read_csv(
            filepath,
            sep="\t",
            skiprows=skip_rows,
            low_memory=False,
            on_bad_lines="warn",
            **kwargs,
        )

    return df


def load_mappings(flybase_data_dir: Path):
    columns_as_strings = {
        "current_symbol": str,
        "current_fullname": str,
        "fullname_synonym(s)": str,
        "symbol_synonym(s)": str,
        "primary_FBid": str,
    }

    synonym_path = find_latest_tsv(flybase_data_dir / "Genes", "fb_synonym")
    print(f"Loading FlyBase synonym table: {synonym_path}")
    mappings_df = load_flybase_tsv(
        synonym_path,
        dtype=columns_as_strings,
        keep_default_na=False,
    )
    mappings_df = mappings_df[mappings_df.organism_abbreviation == "Dmel"]
    mappings_df = mappings_df[mappings_df["primary_FBid"].str.startswith("FBgn", na=False)]

    synonym_to_fbgnid_map = create_expanded_mappings(mappings_df)
    gene_to_fbgnid_main = dict(
        zip(
            mappings_df["current_symbol"].astype(str).str.strip(),
            mappings_df["primary_FBid"].astype(str).str.strip(),
        )
    )
    fbgnid_to_symbol_map = dict(
        zip(
            mappings_df["primary_FBid"].astype(str).str.strip(),
            mappings_df["current_symbol"].astype(str).str.strip(),
        )
    )
    return gene_to_fbgnid_main, synonym_to_fbgnid_map, fbgnid_to_symbol_map


def process_csv_file(
    csv_file, gene_column, gene_to_fbgnid_main, synonym_to_fbgnid_map, fbgnid_to_symbol_map
):
    print(f"Processing file: {csv_file}")

    header_row = 0
    try:
        temp_df = pd.read_csv(csv_file, nrows=10, header=None, encoding="utf-8-sig")
        for i, row in temp_df.iterrows():
            if gene_column in row.astype(str).values:
                header_row = i
                print(f"Detected header at row {i}")
                break
    except Exception as e:
        print(f"Header detection skipped/failed (using default 0): {e}")

    custom_na_values = [
        "",
        "#N/A",
        "#N/A N/A",
        "#NA",
        "-1.#IND",
        "-1.#QNAN",
        "-NaN",
        "-nan",
        "1.#IND",
        "1.#QNAN",
        "<NA>",
        "N/A",
        "NA",
        "NULL",
        "NaN",
        "n/a",
        "nan",
        "null",
    ]
    df = pd.read_csv(
        csv_file,
        encoding="utf-8-sig",
        header=header_row,
        na_values=custom_na_values,
        keep_default_na=False,
        low_memory=False,
    )

    if gene_column not in df.columns:
        if "GO" in csv_file and "Genes" in df.columns:
            df["Genes"] = df["Genes"].astype(str).str.split(", ")
            df = df.explode("Genes").reset_index(drop=True)
            df[gene_column] = df["Genes"]
        else:
            print(f"'{gene_column}' column not found in {csv_file}. Skipping file.")
            return

    working_col = gene_column + "_new"
    corrected_col = "corrected_" + gene_column

    for existing_col in ("flybase_gene_id", "primary_symbol", corrected_col):
        if existing_col in df.columns:
            df.drop(columns=[existing_col], inplace=True)

    # Keep the user's original `gene_column` untouched. All spelling/symbol/
    # capitalization normalization happens on a throwaway working column; the
    # form that actually matched is captured in `corrected_col`.
    df[gene_column] = df[gene_column].astype(str)
    df[working_col] = replace_symbol(df[gene_column])
    df = df.dropna(subset=[working_col])
    df = map_gene_ids(
        df, gene_to_fbgnid_main, synonym_to_fbgnid_map, working_col, corrected_col
    )
    df["primary_symbol"] = df["flybase_gene_id"].map(fbgnid_to_symbol_map)

    # Only surface `corrected_col` when it actually differs from the user's
    # original input; otherwise mark it as "-" to signal no spelling correction.
    no_correction_mask = df[corrected_col].astype(str) == df[gene_column].astype(str)
    df.loc[no_correction_mask, corrected_col] = "-"

    df.fillna("-", inplace=True)
    df.drop(columns=[working_col], inplace=True)

    # Order the derived columns up front, preserving any other original columns after.
    ordered = [gene_column, corrected_col, "primary_symbol", "flybase_gene_id"]
    remaining = [c for c in df.columns if c not in ordered]
    df = df[ordered + remaining]
    df.to_csv(csv_file, index=False, encoding="utf-8-sig")

    print(f"File saved to {csv_file} with {df.shape[0]} rows.")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Map input fly gene symbols to FlyBase FBgn IDs in CSV files."
    )
    parser.add_argument("input_directory", help="Directory containing input CSV files")
    parser.add_argument(
        "input_gene_col",
        nargs="?",
        default="ext_gene",
        help="Input gene symbol column (default: ext_gene)",
    )
    parser.add_argument(
        "--flybase-data-dir",
        default=os.environ.get("FLYBASE_DATA_DIR", str(DEFAULT_FLYBASE_DATA)),
        help="Path to FlyBase data directory (default: Data/FlyBase path)",
    )
    args = parser.parse_args(argv)

    input_dir = Path(args.input_directory)
    if not input_dir.is_dir():
        print(f"Error: '{input_dir}' is not a valid directory")
        return 1

    flybase_data_dir = Path(args.flybase_data_dir)
    if not flybase_data_dir.exists():
        print(f"Error: FlyBase data directory not found at: {flybase_data_dir}")
        return 1

    csv_files = sorted(glob(str(input_dir / "*.csv")))
    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        return 0

    try:
        gene_to_fbgnid_main, synonym_to_fbgnid_map, fbgnid_to_symbol_map = load_mappings(
            flybase_data_dir
        )
    except Exception as e:
        print(f"Error loading FlyBase synonym mappings: {e}")
        return 1

    for csv_file in csv_files:
        process_csv_file(
            csv_file,
            args.input_gene_col,
            gene_to_fbgnid_main,
            synonym_to_fbgnid_map,
            fbgnid_to_symbol_map,
        )

    print("All files processed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
