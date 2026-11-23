"""Functional categories for genes (SPEC section 8.1).

Priority: Pharokka ``function`` qualifier, then PHROG ID lookup (M3), then keyword rules on the
product name. Keyword categories are approximate and the UI labels them as such.
"""

from __future__ import annotations

import re

import pandas as pd

from phageviz.models import CATEGORIES, DEFAULT_PRODUCT, RNA_CATEGORY, RNA_TYPES

# Tried top to bottom; the first match wins. "connector" comes first so that
# "head-tail connector" isn't caught by the "tail" rule.
KEYWORD_RULES: tuple[tuple[str, re.Pattern[str]], ...] = tuple(
    (category, re.compile(pattern, re.IGNORECASE))
    for category, pattern in (
        ("connector", r"head.tail|connector|neck|adaptor"),
        ("lysis", r"holin|endolysin|lysozyme|spanin|lysis"),
        ("head and packaging", r"terminase|portal|capsid|scaffold|packag|prohead"),
        ("tail", r"tail|fib(?:er|re)|baseplate|sheath|tape.measure"),
        ("integration and excision", r"integrase|excisionase|recombinase|transposase"),
        (
            "DNA, RNA and nucleotide metabolism",
            r"polymerase|helicase|primase|nuclease|ligase|reductase|kinase",
        ),
        ("transcription regulation", r"repressor|regulator|antiterminat|sigma|anti-sigma"),
        (
            "moron, auxiliary metabolic gene and host takeover",
            r"superinfection|exclusion|anti-restriction|toxin",
        ),
        ("unknown function", r"hypothetical|unknown|uncharacteri[sz]ed|DUF"),
    )
)
NO_MATCH_CATEGORY = "other"

_CATEGORY_BY_KEY = {category.casefold(): category for category in CATEGORIES}


def categorize(product: str) -> str:
    """Assign a category to one product name using the keyword rules.

    Args:
        product: Product name, e.g. ``"portal protein"``. Blank means hypothetical.

    Returns:
        One of the 10 PHROG categories; ``"other"`` when no rule matches.
    """
    # Plain-Python path over the same KEYWORD_RULES; building a Series for one value costs
    # ~25 ms. tests/test_annotation.py checks it agrees with the vectorised path.
    text = product if product and product.strip() else DEFAULT_PRODUCT
    for category, pattern in KEYWORD_RULES:
        if pattern.search(text):
            return category
    return NO_MATCH_CATEGORY


def assign_categories(types: pd.Series, products: pd.Series, functions: pd.Series) -> pd.DataFrame:
    """Assign a category and its source to every feature.

    CDS: a Pharokka ``function`` that names a PHROG category (case and spacing ignored) gives
    source ``"phrog"``; otherwise the keyword rules on ``product`` give source ``"keyword"``.
    RNA genes are always ``"other"`` / ``"keyword"`` (SPEC 4.2). PHROG ID lookup is applied
    later by ``apply_categories``, which only upgrades ``"keyword"`` rows.

    Args:
        types: Feature types (``CDS``, ``tRNA``, ``tmRNA``).
        products: Product names; missing or blank means ``"hypothetical protein"``.
        functions: Pharokka ``function`` qualifier values; missing when absent.

    Returns:
        DataFrame with ``category`` and ``category_source`` columns (plain strings), sharing
        the index of ``types``.

    Raises:
        ValueError: The three inputs have different lengths.
    """
    if not len(types) == len(products) == len(functions):
        raise ValueError("types, products and functions must have the same length")
    index = types.index
    types = pd.Series(types.to_numpy(), index=index, dtype="string")
    products = pd.Series(products.to_numpy(), index=index, dtype="string")
    functions = pd.Series(functions.to_numpy(), index=index, dtype="string")

    category = _keyword_categories(products)
    source = pd.Series("keyword", index=index, dtype="string")

    from_function = _normalise(functions).map(_CATEGORY_BY_KEY).astype("string")
    use_function = from_function.notna()  # RNA rows are overridden just below
    category = category.mask(use_function, from_function)
    source = source.mask(use_function, "phrog")

    is_rna = types.isin(RNA_TYPES).fillna(False)
    category = category.mask(is_rna, RNA_CATEGORY)
    source = source.mask(is_rna, "keyword")
    return pd.DataFrame({"category": category, "category_source": source}, index=index)


def _keyword_categories(products: pd.Series) -> pd.Series:
    """Vectorised keyword rules: each rule fills only rows no earlier rule matched."""
    blank = products.isna() | (products.str.strip() == "")
    text = products.mask(blank, DEFAULT_PRODUCT)
    result = pd.Series(pd.NA, index=products.index, dtype="string")
    for category, pattern in KEYWORD_RULES:
        hit = result.isna() & text.str.contains(pattern, regex=True).fillna(False)
        result = result.mask(hit, category)
    return result.fillna(NO_MATCH_CATEGORY)


def _normalise(values: pd.Series) -> pd.Series:
    return values.str.strip().str.replace(r"\s+", " ", regex=True).str.casefold()
