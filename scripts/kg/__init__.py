"""
Gene Knowledge Graph - Hybrid Query System

Combines SPARQL (Fuseki) for structured queries with RuVector for semantic search.
"""

from .hybrid import HybridKnowledgeService
from .sparql_client import SparqlClient

__all__ = ['HybridKnowledgeService', 'SparqlClient']
