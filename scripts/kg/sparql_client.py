#!/usr/bin/env python3
"""
SPARQL Client for Gene Knowledge Graph

Provides a clean interface to the Fuseki SPARQL endpoint.
"""

import json
import urllib.request
import urllib.parse
from typing import Dict, List, Any, Optional

# Default endpoint
DEFAULT_ENDPOINT = "http://localhost:3030/gene/sparql"

# Common prefixes for the gene.ai namespace
PREFIXES = """
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX data: <https://gene.ai/data/>
PREFIX tax: <https://gene.ai/taxonomy/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
"""


class SparqlClient:
    """Client for executing SPARQL queries against Fuseki."""

    def __init__(self, endpoint: str = DEFAULT_ENDPOINT):
        self.endpoint = endpoint

    def query(self, sparql: str, add_prefixes: bool = True) -> Dict[str, Any]:
        """Execute a SPARQL query and return raw JSON results."""
        if add_prefixes and not sparql.strip().upper().startswith("PREFIX"):
            sparql = PREFIXES + sparql

        data = urllib.parse.urlencode({"query": sparql}).encode()
        req = urllib.request.Request(
            self.endpoint,
            data=data,
            headers={"Accept": "application/json"}
        )

        with urllib.request.urlopen(req) as response:
            return json.loads(response.read().decode())

    def query_simple(self, sparql: str) -> List[Dict[str, str]]:
        """Execute query and return simplified list of dicts."""
        result = self.query(sparql)
        return [
            {k: v.get("value", "") for k, v in binding.items()}
            for binding in result.get("results", {}).get("bindings", [])
        ]

    def count_triples(self) -> int:
        """Count total triples in the graph."""
        result = self.query("SELECT (COUNT(*) AS ?c) WHERE { ?s ?p ?o }")
        return int(result["results"]["bindings"][0]["c"]["value"])

    def list_sources(self, limit: int = 100) -> List[Dict[str, str]]:
        """List all data sources."""
        return self.query_simple(f"""
            SELECT ?source ?label ?tier ?maintainer WHERE {{
                ?source a ds:DataSource ;
                        rdfs:label ?label .
                OPTIONAL {{ ?source ds:tier ?tier }}
                OPTIONAL {{ ?source ds:maintainer ?maintainer }}
            }}
            ORDER BY ?label
            LIMIT {limit}
        """)

    def get_source(self, source_id: str) -> Dict[str, Any]:
        """Get full details for a data source by ID."""
        # Get main properties
        props = self.query_simple(f"""
            SELECT ?p ?o WHERE {{
                data:{source_id} ?p ?o .
            }}
        """)

        # Get access methods
        access = self.query_simple(f"""
            SELECT ?method ?url ?format WHERE {{
                ?access a ds:AccessMethod ;
                        ds:forSource data:{source_id} ;
                        ds:methodType ?method .
                OPTIONAL {{ ?access ds:baseUrl ?url }}
                OPTIONAL {{ ?access ds:format ?format }}
            }}
        """)

        # Get version info
        versions = self.query_simple(f"""
            SELECT ?version ?date ?size WHERE {{
                ?v a ds:Version ;
                   ds:forSource data:{source_id} .
                OPTIONAL {{ ?v ds:versionName ?version }}
                OPTIONAL {{ ?v ds:releaseDate ?date }}
                OPTIONAL {{ ?v ds:totalSize ?size }}
            }}
        """)

        # Get cross-references
        xrefs = self.query_simple(f"""
            SELECT ?target ?idType WHERE {{
                ?xref a ds:CrossReference ;
                      ds:forSource data:{source_id} ;
                      ds:targetDatabase ?target .
                OPTIONAL {{ ?xref ds:identifierType ?idType }}
            }}
        """)

        return {
            "id": source_id,
            "uri": f"https://gene.ai/data/{source_id}",
            "properties": {p["p"].split("#")[-1]: p["o"] for p in props},
            "access_methods": access,
            "versions": versions,
            "cross_references": xrefs
        }

    def find_by_category(self, category: str) -> List[Dict[str, str]]:
        """Find sources by category name."""
        return self.query_simple(f"""
            SELECT ?source ?label WHERE {{
                ?source a ds:DataSource ;
                        rdfs:label ?label ;
                        ds:category ?cat .
                ?cat rdfs:label ?catLabel .
                FILTER (CONTAINS(LCASE(?catLabel), LCASE("{category}")))
            }}
            ORDER BY ?label
        """)

    def find_by_tier(self, tier: int) -> List[Dict[str, str]]:
        """Find sources by tier (1, 2, or 3)."""
        return self.query_simple(f"""
            SELECT ?source ?label WHERE {{
                ?source a ds:DataSource ;
                        rdfs:label ?label ;
                        ds:tier tax:Tier{tier} .
            }}
            ORDER BY ?label
        """)

    def find_with_api(self) -> List[Dict[str, str]]:
        """Find sources with REST API access."""
        return self.query_simple("""
            SELECT DISTINCT ?source ?label ?url WHERE {
                ?access a ds:AccessMethod ;
                        ds:forSource ?source ;
                        ds:methodType ?method ;
                        ds:baseUrl ?url .
                ?source rdfs:label ?label .
                FILTER (CONTAINS(LCASE(?method), "api") || CONTAINS(LCASE(?method), "rest"))
            }
        """)

    def find_large_datasets(self, min_size_gb: float = 1.0) -> List[Dict[str, str]]:
        """Find sources with large datasets."""
        return self.query_simple("""
            SELECT ?source ?label ?size WHERE {
                ?v a ds:Version ;
                   ds:forSource ?source ;
                   ds:totalSize ?size .
                ?source rdfs:label ?label .
            }
            ORDER BY DESC(?size)
        """)

    def get_entity_counts(self) -> Dict[str, int]:
        """Get counts of all entity types."""
        results = self.query_simple("""
            SELECT ?type (COUNT(?s) AS ?count) WHERE {
                ?s a ?type .
                FILTER(STRSTARTS(STR(?type), "https://gene.ai"))
            }
            GROUP BY ?type
            ORDER BY DESC(?count)
        """)
        return {r["type"].split("#")[-1]: int(r["count"]) for r in results}

    def search_labels(self, term: str, limit: int = 20) -> List[Dict[str, str]]:
        """Search for entities by label."""
        return self.query_simple(f"""
            SELECT ?entity ?label ?type WHERE {{
                ?entity rdfs:label ?label ;
                        a ?type .
                FILTER (CONTAINS(LCASE(?label), LCASE("{term}")))
                FILTER (STRSTARTS(STR(?type), "https://gene.ai"))
            }}
            ORDER BY ?label
            LIMIT {limit}
        """)


# Convenience functions for CLI usage
def create_client(endpoint: str = DEFAULT_ENDPOINT) -> SparqlClient:
    """Create a SPARQL client instance."""
    return SparqlClient(endpoint)


if __name__ == "__main__":
    import sys

    client = SparqlClient()

    if len(sys.argv) < 2:
        print("Gene Knowledge Graph - SPARQL Client")
        print("=" * 50)
        print(f"Endpoint: {client.endpoint}")
        print(f"Triples: {client.count_triples():,}")
        print()
        print("Entity counts:")
        for entity, count in client.get_entity_counts().items():
            print(f"  {entity}: {count}")
        print()
        print("Usage: python3 sparql_client.py 'SPARQL QUERY'")
        sys.exit(0)

    sparql_query = sys.argv[1]
    result = client.query(sparql_query)
    print(json.dumps(result, indent=2))
