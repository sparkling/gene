"""
Gene Knowledge Graph - Hybrid Query System

Combines SPARQL (Fuseki) for structured queries with RuVector for semantic search.

API Module:
    Provides clean, MCP-ready functions that return user-friendly dicts.
    See `scripts.kg.api` for the clean API.
"""

from .hybrid import HybridKnowledgeService
from .sparql_client import SparqlClient
from .api import (
    search_sources,
    list_sources,
    get_source,
    find_connections,
    compare_sources,
    find_by_size,
    find_with_api,
    get_stats,
    clean_output,
)
from .router import (
    detect_intent,
    route_to_function,
    get_all_intents,
    get_intent_description,
    CATEGORY_MAP,
    INTENT_FUNCTION_MAP,
)

__all__ = [
    # Services
    'HybridKnowledgeService',
    'SparqlClient',
    # Clean API functions
    'search_sources',
    'list_sources',
    'get_source',
    'find_connections',
    'compare_sources',
    'find_by_size',
    'find_with_api',
    'get_stats',
    'clean_output',
    # Router functions
    'detect_intent',
    'route_to_function',
    'get_all_intents',
    'get_intent_description',
    'CATEGORY_MAP',
    'INTENT_FUNCTION_MAP',
]
