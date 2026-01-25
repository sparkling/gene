#!/usr/bin/env python3
"""
Gene Knowledge Graph Query CLI

A simple CLI for querying the Gene Knowledge Graph.

Usage:
    ./scripts/kg-query.py                    # Show stats
    ./scripts/kg-query.py list               # List data sources
    ./scripts/kg-query.py search <query>     # Search by label
    ./scripts/kg-query.py source <id>        # Get source details
    ./scripts/kg-query.py sparql '<query>'   # Raw SPARQL query
    ./scripts/kg-query.py tier1              # List Tier 1 sources
    ./scripts/kg-query.py api                # List sources with APIs
    ./scripts/kg-query.py category <name>    # Find by category
    ./scripts/kg-query.py query '<question>' # Natural language query
"""

import sys
import json
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.kg.sparql_client import SparqlClient
from scripts.kg.router import detect_intent, route_to_function
from scripts.kg.api import (
    search_sources, list_sources, get_source,
    find_connections, compare_sources, find_by_size, find_with_api
)


def print_json(data):
    """Pretty print JSON data."""
    print(json.dumps(data, indent=2))


def print_table(rows, headers=None):
    """Print data as a simple table."""
    if not rows:
        print("No results found.")
        return

    if headers:
        print(" | ".join(headers))
        print("-" * 60)

    for row in rows:
        if isinstance(row, dict):
            values = [str(v)[:40] for v in row.values()]
        else:
            values = [str(row)]
        print(" | ".join(values))


def main():
    client = SparqlClient()

    if len(sys.argv) < 2:
        # Show stats
        print("Gene Knowledge Graph")
        print("=" * 50)
        print(f"Endpoint: {client.endpoint}")
        print(f"Namespace: https://gene.ai/")
        print()

        try:
            print(f"Total triples: {client.count_triples():,}")
            print()
            print("Entity counts:")
            for entity, count in client.get_entity_counts().items():
                print(f"  {entity}: {count}")
        except Exception as e:
            print(f"Error connecting to SPARQL endpoint: {e}")
            print()
            print("Start the server with: ./scripts/start-fuseki.sh")
            return 1

        print()
        print("Commands:")
        print("  list              - List all data sources")
        print("  search <query>    - Search by label")
        print("  source <id>       - Get source details")
        print("  sparql '<query>'  - Execute SPARQL query")
        print("  tier1             - List Tier 1 sources")
        print("  api               - List sources with REST APIs")
        print("  category <name>   - Find by category")
        print("  query '<question>'- Natural language query")
        return 0

    command = sys.argv[1].lower()

    try:
        if command == "list":
            limit = int(sys.argv[2]) if len(sys.argv) > 2 else 50
            sources = client.list_sources(limit)
            print(f"Data Sources ({len(sources)} of {limit} max):")
            print("-" * 60)
            for s in sources:
                tier = s.get('tier', '').split('/')[-1] if s.get('tier') else '-'
                print(f"  [{tier}] {s['label']}")

        elif command == "search":
            if len(sys.argv) < 3:
                print("Usage: kg-query.py search <term>")
                return 1
            term = " ".join(sys.argv[2:])
            results = client.search_labels(term)
            print(f"Search results for '{term}':")
            print("-" * 60)
            for r in results:
                entity_type = r['type'].split('#')[-1]
                print(f"  [{entity_type}] {r['label']}")
                print(f"    URI: {r['entity']}")

        elif command == "source":
            if len(sys.argv) < 3:
                print("Usage: kg-query.py source <id>")
                return 1
            source_id = sys.argv[2]
            result = client.get_source(source_id)
            print_json(result)

        elif command == "sparql":
            if len(sys.argv) < 3:
                print("Usage: kg-query.py sparql '<SPARQL query>'")
                return 1
            query = sys.argv[2]
            results = client.query_simple(query)
            print_json(results)

        elif command == "tier1":
            sources = client.find_by_tier(1)
            print(f"Tier 1 Sources ({len(sources)}):")
            print("-" * 60)
            for s in sources:
                print(f"  {s['label']}")

        elif command == "tier2":
            sources = client.find_by_tier(2)
            print(f"Tier 2 Sources ({len(sources)}):")
            print("-" * 60)
            for s in sources:
                print(f"  {s['label']}")

        elif command == "tier3":
            sources = client.find_by_tier(3)
            print(f"Tier 3 Sources ({len(sources)}):")
            print("-" * 60)
            for s in sources:
                print(f"  {s['label']}")

        elif command == "api":
            sources = client.find_with_api()
            print(f"Sources with REST API ({len(sources)}):")
            print("-" * 60)
            for s in sources:
                url = s.get('url', 'N/A')[:50]
                print(f"  {s['label']}")
                print(f"    URL: {url}")

        elif command == "category":
            if len(sys.argv) < 3:
                print("Usage: kg-query.py category <name>")
                print()
                print("Categories: genetics, compounds, diseases, pathways,")
                print("            traditional, nutrition, proteins, literature, microbiome")
                return 1
            category = " ".join(sys.argv[2:])
            sources = client.find_by_category(category)
            print(f"Sources in '{category}' ({len(sources)}):")
            print("-" * 60)
            for s in sources:
                print(f"  {s['label']}")

        elif command == "counts":
            counts = client.get_entity_counts()
            print("Entity Counts:")
            print("-" * 40)
            for entity, count in sorted(counts.items(), key=lambda x: -x[1]):
                print(f"  {entity}: {count}")

        elif command == "prefixes":
            print("Namespace Prefixes:")
            print("-" * 60)
            print("  ds:    https://gene.ai/ontology/datasource#")
            print("  data:  https://gene.ai/data/")
            print("  tax:   https://gene.ai/taxonomy/")
            print("  rdfs:  http://www.w3.org/2000/01/rdf-schema#")
            print("  skos:  http://www.w3.org/2004/02/skos/core#")
            print("  dct:   http://purl.org/dc/terms/")

        elif command == "query":
            if len(sys.argv) < 3:
                print("Usage: kg-query.py query '<natural language question>'")
                return 1

            question = " ".join(sys.argv[2:])
            intent, params = detect_intent(question)

            # Route to appropriate function
            if intent == 'size':
                results = find_by_size(params.get('limit', 10))
            elif intent == 'tier':
                results = list_sources(tier=params.get('tier', 1))
            elif intent == 'api':
                results = find_with_api()
            elif intent == 'category':
                results = list_sources(category=params.get('category'))
            elif intent == 'detail':
                results = get_source(params.get('source_id', ''))
            elif intent == 'relationship':
                results = find_connections(params.get('source_id', ''))
            elif intent == 'compare':
                sources = params.get('sources', [])
                if len(sources) >= 2:
                    results = compare_sources(sources[0], sources[1])
                else:
                    results = {"error": "Need two sources to compare"}
            elif intent == 'list':
                results = list_sources(limit=params.get('limit', 20))
            else:  # 'search' - semantic search
                results = search_sources(params.get('query', question))

            # Output clean JSON for Claude to format
            print(json.dumps({
                "intent": intent,
                "params": params,
                "results": results
            }, indent=2))

        else:
            print(f"Unknown command: {command}")
            print("Use: list, search, source, sparql, tier1, api, category, query")
            return 1

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
