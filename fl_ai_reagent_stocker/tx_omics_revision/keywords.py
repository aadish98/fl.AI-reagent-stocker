"""Keyword matching and DAVID term parsing."""

from __future__ import annotations

import re

from .constants import KEYWORD_RULES, PATHWAY_BUCKETS


def parse_david_term(term: str) -> tuple[str, str]:
    text = str(term).strip()
    if "~" in text:
        term_id, term_name = text.split("~", 1)
        return term_id.strip(), term_name.strip()
    if ":" in text and re.match(r"^[A-Za-z]{2,}\d+", text):
        term_id, term_name = text.split(":", 1)
        return term_id.strip(), term_name.strip()
    return text, text


def split_david_genes(genes: str) -> list[str]:
    tokens = [part.strip() for part in str(genes).replace(";", ",").split(",")]
    normalized = []
    for token in tokens:
        if not token:
            continue
        if token.upper().startswith("FBGN") and len(token) == 11:
            token = "FBgn" + token[4:]
        normalized.append(token)
    return normalized


def match_buckets(term_name: str) -> tuple[list[str], list[str]]:
    text = str(term_name).lower()
    buckets: list[str] = []
    rules: list[str] = []
    for bucket in PATHWAY_BUCKETS:
        for rule_name, pattern in KEYWORD_RULES[bucket]:
            if re.search(pattern, text, flags=re.IGNORECASE):
                if bucket not in buckets:
                    buckets.append(bucket)
                rules.append(f"{bucket}:{rule_name}")
                break
    return buckets, rules


def proposed_decision(buckets: list[str]) -> str:
    if not buckets:
        return "exclude"
    if len(buckets) == 1:
        return "include"
    return "conflict"


def adjacent_unmatched(term_name: str, matched_names: list[str]) -> bool:
    tokens = {tok for tok in re.findall(r"[A-Za-z]{6,}", str(term_name).lower())}
    if not tokens:
        return False
    matched_tokens = {
        tok
        for name in matched_names
        for tok in re.findall(r"[A-Za-z]{6,}", str(name).lower())
    }
    return bool(tokens & matched_tokens)
