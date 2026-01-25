#!/usr/bin/env python3
"""
SPARQL Query Helper for Gene Knowledge Graph

Usage:
    python3 scripts/sparql.py "SELECT * WHERE { ?s ?p ?o } LIMIT 10"

Or as a module:
    from scripts.sparql import query, list_sources, get_source_details
    results = query("SELECT * WHERE { ?s ?p ?o } LIMIT 10")
"""

import json
import urllib.request
import urllib.parse
from typing import Dict, List, Any, Optional

ENDPOINT = "http://localhost:3030/gene/sparql"

# Common prefixes
PREFIXES = """
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX data: <https://gene.ai/data/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
"""


def query(sparql: str, add_prefixes: bool = True) -> Dict[str, Any]:
    """Execute a SPARQL query and return JSON results."""
    if add_prefixes and not sparql.strip().upper().startswith("PREFIX"):
        sparql = PREFIXES + sparql

    data = urllib.parse.urlencode({"query": sparql}).encode()
    req = urllib.request.Request(
        ENDPOINT,
        data=data,
        headers={"Accept": "application/json"}
    )

    with urllib.request.urlopen(req) as response:
        return json.loads(response.read().decode())


def query_simple(sparql: str) -> List[Dict[str, str]]:
    """Execute query and return simplified list of dicts."""
    result = query(sparql)
    return [
        {k: v.get("value", "") for k, v in binding.items()}
        for binding in result.get("results", {}).get("bindings", [])
    ]


def count_triples() -> int:
    """Count total triples in the graph."""
    result = query("SELECT (COUNT(*) AS ?c) WHERE { ?s ?p ?o }")
    return int(result["results"]["bindings"][0]["c"]["value"])


def list_sources(limit: int = 100) -> List[Dict[str, str]]:
    """List all data sources."""
    return query_simple(f"""
        SELECT ?source ?label ?tier ?maintainer WHERE {{
            ?source a ds:DataSource ;
                    rdfs:label ?label .
            OPTIONAL {{ ?source ds:tier ?tier }}
            OPTIONAL {{ ?source ds:maintainer ?maintainer }}
        }}
        ORDER BY ?label
        LIMIT {limit}
    """)


def get_source_details(source_id: str) -> Dict[str, Any]:
    """Get full details for a data source."""
    # Get main source info
    source = query_simple(f"""
        SELECT * WHERE {{
            data:{source_id} ?p ?o .
        }}
    """)

    # Get access methods
    access = query_simple(f"""
        SELECT ?method ?url ?format WHERE {{
            ?access a ds:AccessMethod ;
                    ds:forSource data:{source_id} ;
                    ds:methodType ?method .
            OPTIONAL {{ ?access ds:baseUrl ?url }}
            OPTIONAL {{ ?access ds:format ?format }}
        }}
    """)

    # Get version/size info
    versions = query_simple(f"""
        SELECT ?version ?date ?size WHERE {{
            ?v a ds:Version ;
               ds:forSource data:{source_id} .
            OPTIONAL {{ ?v ds:versionName ?version }}
            OPTIONAL {{ ?v ds:releaseDate ?date }}
            OPTIONAL {{ ?v ds:totalSize ?size }}
        }}
    """)

    # Get cross-references
    xrefs = query_simple(f"""
        SELECT ?target ?idType WHERE {{
            ?xref a ds:CrossReference ;
                  ds:forSource data:{source_id} ;
                  ds:targetDatabase ?target .
            OPTIONAL {{ ?xref ds:identifierType ?idType }}
        }}
    """)

    return {
        "properties": source,
        "access_methods": access,
        "versions": versions,
        "cross_references": xrefs
    }


def find_sources_by_size(min_size_gb: float = 10.0) -> List[Dict[str, str]]:
    """Find data sources with large datasets."""
    return query_simple(f"""
        SELECT ?source ?label ?size WHERE {{
            ?v a ds:Version ;
               ds:forSource ?source ;
               ds:totalSize ?size .
            ?source rdfs:label ?label .
        }}
        ORDER BY DESC(?size)
    """)


def find_sources_with_api() -> List[Dict[str, str]]:
    """Find sources with REST API access."""
    return query_simple("""
        SELECT DISTINCT ?source ?label ?url WHERE {
            ?access a ds:AccessMethod ;
                    ds:forSource ?source ;
                    ds:methodType ?method ;
                    ds:baseUrl ?url .
            ?source rdfs:label ?label .
            FILTER (CONTAINS(LCASE(?method), "api") || CONTAINS(LCASE(?method), "rest"))
        }
    """)


def get_entity_counts() -> Dict[str, int]:
    """Get counts of all entity types."""
    results = query_simple("""
        SELECT ?type (COUNT(?s) AS ?count) WHERE {
            ?s a ?type .
            FILTER(STRSTARTS(STR(?type), "https://gene.ai"))
        }
        GROUP BY ?type
        ORDER BY DESC(?count)
    """)
    return {r["type"].split("#")[-1]: int(r["count"]) for r in results}


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Gene Knowledge Graph - SPARQL Query Helper")
        print("=" * 50)
        print(f"Endpoint: {ENDPOINT}")
        print(f"Triples: {count_triples():,}")
        print()
        print("Entity counts:")
        for entity, count in get_entity_counts().items():
            print(f"  {entity}: {count}")
        print()
        print("Usage: python3 sparql.py 'SPARQL QUERY'")
        sys.exit(0)

    sparql_query = sys.argv[1]
    result = query(sparql_query)
    print(json.dumps(result, indent=2))
