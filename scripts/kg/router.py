#!/usr/bin/env python3
"""
Intent Detection Router for Gene Knowledge Graph

Routes natural language queries to appropriate functions based on detected intent.
Supports size, tier, API, category, detail, relationship, compare, list, and search intents.

Usage:
    from scripts.kg.router import detect_intent, route_to_function

    # Detect intent from query
    intent, params = detect_intent("show me the 10 largest databases")
    # Returns: ('size', {'limit': 10})

    # Get function name for intent
    func_name = route_to_function(intent)
    # Returns: 'find_large_datasets'
"""

import re
from typing import Dict, List, Tuple, Optional, Any


# Category keyword mapping to full category names
CATEGORY_MAP: Dict[str, str] = {
    'genetics': 'Genetics & Genomics',
    'genomics': 'Genetics & Genomics',
    'variant': 'Genetics & Genomics',
    'drug': 'Compounds & Molecules',
    'compound': 'Compounds & Molecules',
    'chemical': 'Compounds & Molecules',
    'disease': 'Diseases & Phenotypes',
    'phenotype': 'Diseases & Phenotypes',
    'pathway': 'Pathways & Networks',
    'network': 'Pathways & Networks',
    'tcm': 'Traditional Medicine',
    'traditional': 'Traditional Medicine',
    'herbal': 'Traditional Medicine',
    'ayurveda': 'Traditional Medicine',
    'nutrition': 'Nutrition & Food',
    'food': 'Nutrition & Food',
    'protein': 'Proteins & Molecular Biology',
    'literature': 'Literature & Knowledge',
    'pubmed': 'Literature & Knowledge',
    'microbiome': 'Microbiome',
    'gut': 'Microbiome',
}

# Intent to function name mapping
INTENT_FUNCTION_MAP: Dict[str, str] = {
    'size': 'find_large_datasets',
    'tier': 'find_by_tier',
    'api': 'find_with_api',
    'category': 'find_by_category',
    'detail': 'get_source',
    'relationship': 'get_relationships',
    'compare': 'compare_sources',
    'list': 'list_sources',
    'search': 'semantic_search',
}


def detect_intent(query: str) -> Tuple[str, Dict[str, Any]]:
    """
    Detect the intent from a natural language query.

    Analyzes the query text to determine what type of information
    the user is seeking and extracts relevant parameters.

    Args:
        query: Natural language query string

    Returns:
        Tuple of (intent_name, parameters_dict)
        - intent_name: One of 'size', 'tier', 'api', 'category', 'detail',
                       'relationship', 'compare', 'list', 'search'
        - parameters_dict: Dict with intent-specific parameters

    Examples:
        >>> detect_intent("10 largest databases")
        ('size', {'limit': 10})

        >>> detect_intent("production-ready sources")
        ('tier', {'tier': 1})

        >>> detect_intent("databases with REST API")
        ('api', {})

        >>> detect_intent("drug databases")
        ('category', {'category': 'Compounds & Molecules'})

        >>> detect_intent("tell me about ClinVar")
        ('detail', {'source_id': 'ClinVar'})

        >>> detect_intent("compare DrugBank and ChEMBL")
        ('compare', {'sources': ['DrugBank', 'ChEMBL']})
    """
    query_lower = query.lower().strip()

    # Check for size intent (largest, biggest, N largest)
    size_result = _detect_size_intent(query_lower)
    if size_result:
        return size_result

    # Check for tier intent (production-ready, tier 1/2/3, stable)
    tier_result = _detect_tier_intent(query_lower)
    if tier_result:
        return tier_result

    # Check for API intent
    api_result = _detect_api_intent(query_lower)
    if api_result:
        return api_result

    # Check for compare intent (X vs Y, compare X and Y)
    compare_result = _detect_compare_intent(query, query_lower)
    if compare_result:
        return compare_result

    # Check for detail intent (about X, tell me about X)
    detail_result = _detect_detail_intent(query, query_lower)
    if detail_result:
        return detail_result

    # Check for relationship intent (connects to, links to)
    relationship_result = _detect_relationship_intent(query, query_lower)
    if relationship_result:
        return relationship_result

    # Check for category intent (domain keywords)
    category_result = _detect_category_intent(query_lower)
    if category_result:
        return category_result

    # Check for list intent (list all, show databases)
    list_result = _detect_list_intent(query_lower)
    if list_result:
        return list_result

    # Default: semantic search fallback
    return ('search', {'query': query})


def _detect_size_intent(query_lower: str) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Detect size-related queries (largest, biggest, top N)."""
    # Patterns: "10 largest", "biggest databases", "top 5"
    size_patterns = [
        r'(\d+)\s*(?:largest|biggest|top|major)',
        r'(?:largest|biggest|top|major)\s*(\d+)?',
        r'(?:show|find|list)\s*(?:the\s*)?(\d+)?\s*(?:largest|biggest)',
    ]

    for pattern in size_patterns:
        match = re.search(pattern, query_lower)
        if match:
            # Extract limit from match groups
            limit = None
            for group in match.groups():
                if group and group.isdigit():
                    limit = int(group)
                    break
            return ('size', {'limit': limit or 10})

    return None


def _detect_tier_intent(query_lower: str) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Detect tier-related queries (production-ready, stable, tier N)."""
    # Tier 1: production-ready, stable, reference
    if any(kw in query_lower for kw in [
        'production-ready', 'production ready', 'stable',
        'reference', 'gold standard', 'tier 1', 'tier1'
    ]):
        return ('tier', {'tier': 1})

    # Tier 2: established, well-maintained
    if any(kw in query_lower for kw in [
        'established', 'well-maintained', 'well maintained',
        'mature', 'tier 2', 'tier2'
    ]):
        return ('tier', {'tier': 2})

    # Tier 3: emerging, experimental
    if any(kw in query_lower for kw in [
        'emerging', 'experimental', 'new', 'beta',
        'tier 3', 'tier3'
    ]):
        return ('tier', {'tier': 3})

    return None


def _detect_api_intent(query_lower: str) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Detect API access queries."""
    api_keywords = [
        'with api', 'rest api', 'api access', 'programmatic',
        'programmatically', 'web service', 'webservice',
        'endpoint', 'api endpoint'
    ]

    if any(kw in query_lower for kw in api_keywords):
        return ('api', {})

    return None


def _detect_compare_intent(
    query: str,
    query_lower: str
) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Detect comparison queries (X vs Y, compare X and Y)."""
    # Pattern: "compare X and Y" or "X vs Y" or "X versus Y"
    compare_patterns = [
        r'compare\s+(\w+)\s+(?:and|with|to|vs\.?|versus)\s+(\w+)',
        r'(\w+)\s+(?:vs\.?|versus)\s+(\w+)',
        r'difference(?:s)?\s+between\s+(\w+)\s+and\s+(\w+)',
    ]

    for pattern in compare_patterns:
        match = re.search(pattern, query_lower)
        if match:
            source1, source2 = match.groups()
            # Try to find original case from query
            sources = _extract_source_names(query, [source1, source2])
            return ('compare', {'sources': sources})

    return None


def _detect_detail_intent(
    query: str,
    query_lower: str
) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Detect detail/information queries (about X, tell me about X)."""
    # Patterns: "about X", "tell me about X", "what is X", "details on X"
    detail_patterns = [
        r'(?:tell\s+me\s+)?about\s+(\w+)',
        r'what\s+is\s+(\w+)',
        r'(?:details?|info(?:rmation)?)\s+(?:on|about|for)\s+(\w+)',
        r'describe\s+(\w+)',
        r'explain\s+(\w+)',
    ]

    for pattern in detail_patterns:
        match = re.search(pattern, query_lower)
        if match:
            source_lower = match.group(1)
            # Try to find original case from query
            source_id = _extract_source_name(query, source_lower)
            return ('detail', {'source_id': source_id})

    return None


def _detect_relationship_intent(
    query: str,
    query_lower: str
) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Detect relationship queries (connects to, links to, related to)."""
    # Patterns: "connects to X", "links to X", "related to X"
    rel_patterns = [
        r'(?:connect(?:s|ed)?|link(?:s|ed)?|relate[sd]?)\s+to\s+(\w+)',
        r'what\s+(?:connect(?:s|ed)?|link(?:s|ed)?)\s+(?:to\s+)?(\w+)',
        r'(?:cross-?references?|xrefs?)\s+(?:for|of|from)\s+(\w+)',
        r'(\w+)\s+(?:connect(?:s|ed)?|link(?:s|ed)?)\s+to',
    ]

    for pattern in rel_patterns:
        match = re.search(pattern, query_lower)
        if match:
            source_lower = match.group(1)
            source_id = _extract_source_name(query, source_lower)
            return ('relationship', {'source_id': source_id})

    return None


def _detect_category_intent(
    query_lower: str
) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Detect category queries based on domain keywords."""
    # Check each category keyword
    for keyword, category in CATEGORY_MAP.items():
        if keyword in query_lower:
            return ('category', {'category': category})

    return None


def _detect_list_intent(query_lower: str) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Detect list/show all queries."""
    list_patterns = [
        r'list\s+(?:all\s+)?(?:data)?(?:bases?|sources?)?',
        r'show\s+(?:all\s+)?(?:data)?(?:bases?|sources?)',
        r'(?:all|available)\s+(?:data)?(?:bases?|sources?)',
        r'what\s+(?:data)?(?:bases?|sources?)\s+(?:are\s+)?(?:there|available)',
    ]

    for pattern in list_patterns:
        if re.search(pattern, query_lower):
            # Extract limit if specified
            limit_match = re.search(r'(\d+)', query_lower)
            limit = int(limit_match.group(1)) if limit_match else 20
            return ('list', {'limit': limit})

    return None


def _extract_source_name(query: str, source_lower: str) -> str:
    """
    Extract the original-case source name from the query.

    Attempts to find the source name in its original case by
    searching the original query string.
    """
    # Try to find the word in original case
    pattern = re.compile(rf'\b({re.escape(source_lower)})\b', re.IGNORECASE)
    match = pattern.search(query)
    if match:
        return match.group(1)
    return source_lower


def _extract_source_names(query: str, sources_lower: List[str]) -> List[str]:
    """Extract original-case source names from the query."""
    return [_extract_source_name(query, s) for s in sources_lower]


def route_to_function(intent: str) -> str:
    """
    Get the function name to call based on detected intent.

    Args:
        intent: Intent name from detect_intent()

    Returns:
        Function name string that handles this intent

    Examples:
        >>> route_to_function('size')
        'find_large_datasets'

        >>> route_to_function('tier')
        'find_by_tier'

        >>> route_to_function('search')
        'semantic_search'
    """
    return INTENT_FUNCTION_MAP.get(intent, 'semantic_search')


def get_all_intents() -> List[str]:
    """
    Get list of all supported intents.

    Returns:
        List of intent names
    """
    return list(INTENT_FUNCTION_MAP.keys())


def get_intent_description(intent: str) -> str:
    """
    Get a human-readable description of an intent.

    Args:
        intent: Intent name

    Returns:
        Description string
    """
    descriptions = {
        'size': 'Find largest datasets by size',
        'tier': 'Filter by production tier (1=stable, 2=established, 3=emerging)',
        'api': 'Find databases with REST API access',
        'category': 'Filter by domain category',
        'detail': 'Get detailed information about a specific source',
        'relationship': 'Find cross-references and connections',
        'compare': 'Compare two or more data sources',
        'list': 'List all available data sources',
        'search': 'Semantic search across all sources',
    }
    return descriptions.get(intent, 'Unknown intent')


# ─────────────────────────────────────────────────────────────────────
# CLI Interface
# ─────────────────────────────────────────────────────────────────────

def main() -> None:
    """CLI interface for intent detection testing."""
    import sys
    import json

    if len(sys.argv) < 2:
        print("Gene Knowledge Graph - Intent Router")
        print("=" * 50)
        print()
        print("Supported intents:")
        for intent in get_all_intents():
            func = route_to_function(intent)
            desc = get_intent_description(intent)
            print(f"  {intent:12} -> {func:20} | {desc}")
        print()
        print("Usage: python3 -m scripts.kg.router '<query>'")
        print()
        print("Examples:")
        print("  python3 -m scripts.kg.router '10 largest databases'")
        print("  python3 -m scripts.kg.router 'production-ready sources'")
        print("  python3 -m scripts.kg.router 'drug databases'")
        print("  python3 -m scripts.kg.router 'tell me about ClinVar'")
        print("  python3 -m scripts.kg.router 'compare DrugBank and ChEMBL'")
        return

    query = " ".join(sys.argv[1:])
    intent, params = detect_intent(query)
    func_name = route_to_function(intent)

    result = {
        'query': query,
        'intent': intent,
        'function': func_name,
        'params': params,
        'description': get_intent_description(intent)
    }

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
