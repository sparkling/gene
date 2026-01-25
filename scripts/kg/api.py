#!/usr/bin/env python3
"""
Clean MCP-Ready API for Gene Knowledge Graph

Provides user-friendly functions that return clean dicts without URIs or SPARQL syntax.
Designed for use from:
- Skills (via Bash/Python)
- MCP server (direct tool calls)

All functions return clean data structures suitable for end users.
"""

import asyncio
from typing import Dict, List, Any, Optional

from .sparql_client import SparqlClient
from .hybrid import HybridKnowledgeService


# Default service instances
_sparql_client: Optional[SparqlClient] = None
_hybrid_service: Optional[HybridKnowledgeService] = None


def _get_sparql() -> SparqlClient:
    """Get or create SPARQL client singleton."""
    global _sparql_client
    if _sparql_client is None:
        _sparql_client = SparqlClient()
    return _sparql_client


def _get_hybrid() -> HybridKnowledgeService:
    """Get or create hybrid service singleton."""
    global _hybrid_service
    if _hybrid_service is None:
        _hybrid_service = HybridKnowledgeService()
    return _hybrid_service


def _extract_tier(tier_value: Optional[str]) -> Optional[int]:
    """
    Extract tier number from URI or string.

    Examples:
        "https://gene.ai/taxonomy/Tier1" -> 1
        "tax:Tier2" -> 2
        "Tier3" -> 3
        None -> None
    """
    if not tier_value:
        return None
    # Extract the last part and find the digit
    tier_str = tier_value.split("/")[-1].replace("Tier", "")
    try:
        return int(tier_str)
    except (ValueError, TypeError):
        return None


def _extract_category(category_value: Optional[str]) -> Optional[str]:
    """
    Extract clean category name from URI or string.

    Examples:
        "https://gene.ai/taxonomy/SequenceDatabases" -> "Sequence Databases"
        "tax:ClinicalVariantDBs" -> "Clinical Variant DBs"
    """
    if not category_value:
        return None
    # If it's a URI, extract the last part
    name = category_value.split("/")[-1].split("#")[-1]
    # Convert CamelCase to spaces
    import re
    name = re.sub(r'([a-z])([A-Z])', r'\1 \2', name)
    name = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1 \2', name)
    return name


def _extract_id(uri: str) -> str:
    """Extract clean ID from URI."""
    return uri.split("/")[-1].split("#")[-1]


def clean_output(sources: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Clean and normalize source data for user-friendly output.

    Strips URIs, extracts tier numbers, cleans category names,
    and returns only user-relevant fields.

    Args:
        sources: Raw source data from SPARQL queries

    Returns:
        List of clean dicts with user-friendly field names and values
    """
    cleaned = []
    for source in sources:
        clean_item = {}

        # Extract name/label
        for key in ["label", "name", "rdfs:label"]:
            if key in source and source[key]:
                clean_item["name"] = source[key]
                break

        # Extract source ID
        for key in ["source", "uri", "id"]:
            if key in source and source[key]:
                clean_item["id"] = _extract_id(source[key])
                break

        # Extract and clean tier
        if "tier" in source:
            tier = _extract_tier(source["tier"])
            if tier:
                clean_item["tier"] = tier

        # Extract and clean category
        for key in ["category", "cat", "catLabel"]:
            if key in source and source[key]:
                clean_item["category"] = _extract_category(source[key])
                break

        # Copy other user-relevant fields directly
        if "description" in source and source["description"]:
            clean_item["description"] = source["description"]

        if "maintainer" in source and source["maintainer"]:
            clean_item["maintainer"] = source["maintainer"]

        # Handle API URL (from access methods)
        for key in ["url", "api_url", "baseUrl"]:
            if key in source and source[key]:
                clean_item["api_url"] = source[key]
                break

        if "license" in source and source["license"]:
            clean_item["license"] = source["license"]

        if "size" in source and source["size"]:
            clean_item["size"] = source["size"]

        if "method" in source and source["method"]:
            clean_item["access_method"] = source["method"]

        cleaned.append(clean_item)

    return cleaned


def search_sources(query: str, limit: int = 10) -> List[Dict[str, Any]]:
    """
    Semantic search for sources matching query.

    Searches data source names, descriptions, and categories
    using both keyword matching and semantic similarity.

    Args:
        query: Search query (e.g., "clinical variant databases")
        limit: Maximum number of results to return

    Returns:
        List of matching sources with clean, user-friendly data.
        Returns empty list on failure.

    Example:
        >>> results = search_sources("genomic variation", limit=5)
        >>> for r in results:
        ...     print(f"{r['name']} (Tier {r.get('tier', '?')})")
    """
    try:
        sparql = _get_sparql()
        results = sparql.search_labels(query, limit)

        # Enrich with more details
        enriched = []
        for r in results:
            if "DataSource" in r.get("type", ""):
                source_id = _extract_id(r.get("entity", ""))
                try:
                    details = sparql.get_source(source_id)
                    props = details.get("properties", {})
                    enriched.append({
                        "source": r.get("entity", ""),
                        "label": r.get("label", ""),
                        "tier": props.get("tier", ""),
                        "category": props.get("category", ""),
                        "description": props.get("description", ""),
                        "maintainer": props.get("maintainer", "")
                    })
                except Exception:
                    enriched.append({
                        "source": r.get("entity", ""),
                        "label": r.get("label", "")
                    })
            else:
                enriched.append({
                    "source": r.get("entity", ""),
                    "label": r.get("label", "")
                })

        return clean_output(enriched)
    except Exception:
        return []


def list_sources(
    category: Optional[str] = None,
    tier: Optional[int] = None,
    limit: int = 20
) -> List[Dict[str, Any]]:
    """
    List sources with optional filters.

    Args:
        category: Filter by category name (e.g., "Sequence", "Clinical")
        tier: Filter by tier level (1, 2, or 3)
        limit: Maximum number of results to return

    Returns:
        List of sources with clean, user-friendly data.
        Returns empty list on failure.

    Example:
        >>> tier1 = list_sources(tier=1)
        >>> clinical = list_sources(category="Clinical")
    """
    try:
        sparql = _get_sparql()

        if tier is not None:
            results = sparql.find_by_tier(tier)
        elif category is not None:
            results = sparql.find_by_category(category)
        else:
            results = sparql.list_sources(limit)

        # Limit results
        results = results[:limit]

        return clean_output(results)
    except Exception:
        return []


def get_source(source_id: str) -> Dict[str, Any]:
    """
    Get full details for a single source.

    Args:
        source_id: The source identifier (e.g., "ClinVar", "NCBI-Gene")

    Returns:
        Dict with source details including name, tier, category,
        description, maintainer, access methods, etc.
        Returns empty dict on failure.

    Example:
        >>> clinvar = get_source("ClinVar")
        >>> print(clinvar['name'], clinvar.get('tier'))
    """
    try:
        sparql = _get_sparql()
        raw = sparql.get_source(source_id)

        if not raw or "properties" not in raw:
            return {}

        props = raw.get("properties", {})

        result = {
            "id": source_id,
            "name": props.get("label", source_id),
        }

        # Add tier if present
        tier = _extract_tier(props.get("tier"))
        if tier:
            result["tier"] = tier

        # Add category if present
        if "category" in props:
            result["category"] = _extract_category(props["category"])

        # Add optional fields
        if props.get("description"):
            result["description"] = props["description"]

        if props.get("maintainer"):
            result["maintainer"] = props["maintainer"]

        if props.get("license"):
            result["license"] = props["license"]

        if props.get("homepage"):
            result["homepage"] = props["homepage"]

        # Add access methods
        access_methods = raw.get("access_methods", [])
        if access_methods:
            result["access"] = []
            for am in access_methods:
                access_entry = {}
                if am.get("method"):
                    access_entry["type"] = am["method"]
                if am.get("url"):
                    access_entry["url"] = am["url"]
                if am.get("format"):
                    access_entry["format"] = am["format"]
                if access_entry:
                    result["access"].append(access_entry)

        # Add version info
        versions = raw.get("versions", [])
        if versions:
            # Get most recent version
            latest = versions[0]
            if latest.get("version"):
                result["version"] = latest["version"]
            if latest.get("date"):
                result["release_date"] = latest["date"]
            if latest.get("size"):
                result["size"] = latest["size"]

        # Add cross-references
        xrefs = raw.get("cross_references", [])
        if xrefs:
            result["connections"] = []
            for xref in xrefs:
                conn = {}
                if xref.get("target"):
                    conn["database"] = _extract_id(xref["target"])
                if xref.get("idType"):
                    conn["identifier_type"] = xref["idType"]
                if conn:
                    result["connections"].append(conn)

        return result
    except Exception:
        return {}


def find_connections(source_id: str) -> List[Dict[str, Any]]:
    """
    Find sources that connect to this one via shared identifiers.

    Args:
        source_id: The source identifier to find connections for

    Returns:
        List of connected sources with connection details.
        Returns empty list on failure.

    Example:
        >>> connected = find_connections("ClinVar")
        >>> for c in connected:
        ...     print(f"{c['name']} via {c.get('connection_type', 'unknown')}")
    """
    try:
        sparql = _get_sparql()

        # Find outgoing connections (this source references others)
        outgoing = sparql.query_simple(f"""
            SELECT ?target ?targetLabel ?idType WHERE {{
                ?xref a ds:CrossReference ;
                      ds:forSource data:{source_id} ;
                      ds:targetDatabase ?target .
                ?target rdfs:label ?targetLabel .
                OPTIONAL {{ ?xref ds:identifierType ?idType }}
            }}
        """)

        # Find incoming connections (others reference this source)
        incoming = sparql.query_simple(f"""
            SELECT ?source ?sourceLabel ?idType WHERE {{
                ?xref a ds:CrossReference ;
                      ds:forSource ?source ;
                      ds:targetDatabase data:{source_id} .
                ?source rdfs:label ?sourceLabel .
                OPTIONAL {{ ?xref ds:identifierType ?idType }}
            }}
        """)

        connections = []

        for conn in outgoing:
            connections.append({
                "name": conn.get("targetLabel", ""),
                "id": _extract_id(conn.get("target", "")),
                "direction": "outgoing",
                "connection_type": conn.get("idType", "shared identifier")
            })

        for conn in incoming:
            connections.append({
                "name": conn.get("sourceLabel", ""),
                "id": _extract_id(conn.get("source", "")),
                "direction": "incoming",
                "connection_type": conn.get("idType", "shared identifier")
            })

        return connections
    except Exception:
        return []


def compare_sources(source_a: str, source_b: str) -> Dict[str, Any]:
    """
    Compare two sources side by side.

    Args:
        source_a: First source identifier
        source_b: Second source identifier

    Returns:
        Dict with comparison data for both sources.
        Returns empty dict on failure.

    Example:
        >>> cmp = compare_sources("ClinVar", "dbSNP")
        >>> print(cmp['source_a']['name'], "vs", cmp['source_b']['name'])
    """
    try:
        a_details = get_source(source_a)
        b_details = get_source(source_b)

        if not a_details or not b_details:
            return {}

        # Find shared connections
        a_connections = set(c.get("id", "") for c in find_connections(source_a))
        b_connections = set(c.get("id", "") for c in find_connections(source_b))
        shared = a_connections & b_connections

        return {
            "source_a": a_details,
            "source_b": b_details,
            "shared_connections": list(shared),
            "comparison": {
                "same_tier": a_details.get("tier") == b_details.get("tier"),
                "same_category": a_details.get("category") == b_details.get("category"),
                "shared_connection_count": len(shared)
            }
        }
    except Exception:
        return {}


def find_by_size(limit: int = 10) -> List[Dict[str, Any]]:
    """
    Find largest databases by size.

    Args:
        limit: Maximum number of results to return

    Returns:
        List of sources ordered by size (largest first).
        Returns empty list on failure.

    Example:
        >>> large = find_by_size(5)
        >>> for s in large:
        ...     print(f"{s['name']}: {s.get('size', 'unknown')}")
    """
    try:
        sparql = _get_sparql()
        results = sparql.find_large_datasets(min_size_gb=0)
        results = results[:limit]
        return clean_output(results)
    except Exception:
        return []


def find_with_api() -> List[Dict[str, Any]]:
    """
    Find sources with REST API access.

    Returns:
        List of sources that have REST API endpoints available.
        Returns empty list on failure.

    Example:
        >>> api_sources = find_with_api()
        >>> for s in api_sources:
        ...     print(f"{s['name']}: {s.get('api_url', 'N/A')}")
    """
    try:
        sparql = _get_sparql()
        results = sparql.find_with_api()
        return clean_output(results)
    except Exception:
        return []


def get_stats() -> Dict[str, Any]:
    """
    Get knowledge graph statistics.

    Returns:
        Dict with counts of entities, triples, etc.
        Returns empty dict on failure.
    """
    try:
        sparql = _get_sparql()
        counts = sparql.get_entity_counts()
        total_triples = sparql.count_triples()

        return {
            "total_triples": total_triples,
            "data_sources": counts.get("DataSource", 0),
            "access_methods": counts.get("AccessMethod", 0),
            "versions": counts.get("Version", 0),
            "cross_references": counts.get("CrossReference", 0),
            "categories": counts.get("Category", 0)
        }
    except Exception:
        return {}


# Async variants for use with HybridKnowledgeService
async def search_sources_async(query: str, limit: int = 10) -> List[Dict[str, Any]]:
    """
    Async semantic search for sources (uses hybrid service).

    Same as search_sources but uses async hybrid search with
    RuVector embeddings for better semantic matching.
    """
    try:
        hybrid = _get_hybrid()
        results = await hybrid.hybrid_search(query, k=limit)

        cleaned = []
        for r in results:
            item = {
                "id": r.id,
                "name": r.label,
                "score": r.score
            }

            if r.description:
                item["description"] = r.description

            if r.sparql_data:
                tier = _extract_tier(r.sparql_data.get("tier"))
                if tier:
                    item["tier"] = tier
                if r.sparql_data.get("category"):
                    item["category"] = r.sparql_data["category"]
                if r.sparql_data.get("maintainer"):
                    item["maintainer"] = r.sparql_data["maintainer"]

            cleaned.append(item)

        return cleaned
    except Exception:
        return []


if __name__ == "__main__":
    import json

    print("Gene Knowledge Graph - Clean API Test")
    print("=" * 50)

    # Test get_stats
    print("\n1. Knowledge Graph Stats:")
    stats = get_stats()
    if stats:
        print(f"   Triples: {stats.get('total_triples', 0):,}")
        print(f"   Data Sources: {stats.get('data_sources', 0)}")
    else:
        print("   (Could not connect to SPARQL endpoint)")

    # Test list_sources
    print("\n2. List Sources (limit=5):")
    sources = list_sources(limit=5)
    if sources:
        for s in sources:
            tier_str = f"Tier {s['tier']}" if s.get('tier') else "No tier"
            print(f"   - {s.get('name', 'Unknown')} ({tier_str})")
    else:
        print("   (No sources found or connection error)")

    # Test search_sources
    print("\n3. Search for 'clinical':")
    results = search_sources("clinical", limit=3)
    if results:
        for r in results:
            print(f"   - {r.get('name', 'Unknown')}")
    else:
        print("   (No results)")

    # Test find_with_api
    print("\n4. Sources with API access:")
    api_sources = find_with_api()
    if api_sources:
        for s in api_sources[:3]:
            print(f"   - {s.get('name', 'Unknown')}: {s.get('api_url', 'N/A')}")
    else:
        print("   (No API sources found)")

    print("\n" + "=" * 50)
    print("API ready for MCP integration")
