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
]
