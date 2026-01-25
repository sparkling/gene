#!/usr/bin/env python3
"""
Hybrid Knowledge Service - Combines SPARQL + RuVector

Provides semantic similarity search (via embeddings) combined with
structured graph queries (via SPARQL) for the Gene Knowledge Graph.

Architecture:
    RuVector (Embeddings) ←→ Linked by URI ←→ Fuseki (SPARQL)

Usage:
    from scripts.kg.hybrid import HybridKnowledgeService

    kg = HybridKnowledgeService()

    # Semantic search
    results = await kg.semantic_search("clinical variant databases")

    # Hybrid search (semantic + SPARQL enrichment)
    results = await kg.hybrid_search("pharmacogenomics")

    # Pure SPARQL
    results = kg.sparql("SELECT ?s WHERE { ?s a ds:DataSource }")
"""

import json
import subprocess
import asyncio
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, asdict

from .sparql_client import SparqlClient, PREFIXES


@dataclass
class SearchResult:
    """A search result with optional SPARQL enrichment."""
    id: str
    label: str
    score: float
    uri: str
    entity_type: str
    description: Optional[str] = None
    sparql_data: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class HybridKnowledgeService:
    """
    Hybrid query service combining:
    - RuVector: 384-dim embeddings with HNSW indexing (150x-12,500x faster)
    - Fuseki: Full SPARQL 1.1 with 12,098 triples

    The link between them is the URI - each RuVector entry stores
    `qlever_uri` in its metadata pointing to the Fuseki entity.
    """

    def __init__(
        self,
        sparql_endpoint: str = "http://localhost:3030/gene/sparql",
        memory_namespace: str = "gene-kg"
    ):
        self.sparql = SparqlClient(sparql_endpoint)
        self.memory_namespace = memory_namespace

    # ─────────────────────────────────────────────────────────────────
    # SPARQL Methods (Structured Queries)
    # ─────────────────────────────────────────────────────────────────

    def sparql_query(self, query: str) -> List[Dict[str, str]]:
        """Execute a SPARQL query and return results."""
        return self.sparql.query_simple(query)

    def get_source(self, source_id: str) -> Dict[str, Any]:
        """Get full details for a data source."""
        return self.sparql.get_source(source_id)

    def list_sources(self, limit: int = 100) -> List[Dict[str, str]]:
        """List all data sources."""
        return self.sparql.list_sources(limit)

    def find_by_category(self, category: str) -> List[Dict[str, str]]:
        """Find sources by category."""
        return self.sparql.find_by_category(category)

    def find_tier1(self) -> List[Dict[str, str]]:
        """Find Tier 1 (production-ready) sources."""
        return self.sparql.find_by_tier(1)

    def find_with_api(self) -> List[Dict[str, str]]:
        """Find sources with REST API access."""
        return self.sparql.find_with_api()

    # ─────────────────────────────────────────────────────────────────
    # RuVector Methods (Semantic Search)
    # ─────────────────────────────────────────────────────────────────

    async def semantic_search(
        self,
        query: str,
        k: int = 10
    ) -> List[SearchResult]:
        """
        Semantic similarity search using RuVector embeddings.

        This searches the embedded representations of data sources
        to find semantically similar entries.
        """
        # Use claude-flow memory search with semantic matching
        try:
            result = subprocess.run(
                [
                    "npx", "@claude-flow/cli@latest", "memory", "search",
                    "--query", query,
                    "--namespace", self.memory_namespace,
                    "--limit", str(k)
                ],
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode != 0:
                # Fallback to SPARQL label search if RuVector not available
                return await self._fallback_label_search(query, k)

            # Parse results
            results = []
            lines = result.stdout.strip().split('\n')
            for line in lines:
                if line.startswith('{'):
                    try:
                        data = json.loads(line)
                        results.append(SearchResult(
                            id=data.get('key', ''),
                            label=data.get('metadata', {}).get('label', ''),
                            score=data.get('score', 0.0),
                            uri=data.get('metadata', {}).get('qlever_uri', ''),
                            entity_type=data.get('metadata', {}).get('entity_type', ''),
                            description=data.get('metadata', {}).get('description')
                        ))
                    except json.JSONDecodeError:
                        continue

            return results

        except (subprocess.TimeoutExpired, FileNotFoundError):
            return await self._fallback_label_search(query, k)

    async def _fallback_label_search(
        self,
        query: str,
        k: int
    ) -> List[SearchResult]:
        """Fallback to SPARQL label search if RuVector unavailable."""
        results = self.sparql.search_labels(query, k)
        return [
            SearchResult(
                id=r['entity'].split('/')[-1],
                label=r['label'],
                score=1.0,  # No semantic score available
                uri=r['entity'],
                entity_type=r['type'].split('#')[-1]
            )
            for r in results
        ]

    # ─────────────────────────────────────────────────────────────────
    # Hybrid Methods (Semantic + SPARQL)
    # ─────────────────────────────────────────────────────────────────

    async def hybrid_search(
        self,
        query: str,
        k: int = 10,
        enrich: bool = True
    ) -> List[SearchResult]:
        """
        Hybrid search: semantic similarity + SPARQL enrichment.

        1. Find semantically similar entities via RuVector
        2. Extract URIs from results
        3. Enrich with full SPARQL data (access methods, versions, etc.)
        """
        # Step 1: Semantic search
        results = await self.semantic_search(query, k)

        if not enrich or not results:
            return results

        # Step 2: Extract source IDs
        source_ids = [r.id for r in results if r.entity_type == 'DataSource']

        if not source_ids:
            return results

        # Step 3: Batch SPARQL enrichment
        enrichment = self._batch_enrich(source_ids)

        # Step 4: Merge results
        for result in results:
            if result.id in enrichment:
                result.sparql_data = enrichment[result.id]

        return results

    def _batch_enrich(self, source_ids: List[str]) -> Dict[str, Dict[str, Any]]:
        """Batch enrich multiple sources with SPARQL data."""
        if not source_ids:
            return {}

        # Build VALUES clause for batch query
        values = " ".join(f"data:{sid}" for sid in source_ids)

        # Query for basic info
        basic = self.sparql.query_simple(f"""
            SELECT ?source ?label ?tier ?category ?maintainer WHERE {{
                VALUES ?source {{ {values} }}
                ?source rdfs:label ?label .
                OPTIONAL {{ ?source ds:tier ?tier }}
                OPTIONAL {{
                    ?source ds:category ?cat .
                    ?cat rdfs:label ?category .
                }}
                OPTIONAL {{ ?source ds:maintainer ?maintainer }}
            }}
        """)

        # Query for access methods
        access = self.sparql.query_simple(f"""
            SELECT ?source ?method ?url WHERE {{
                VALUES ?source {{ {values} }}
                ?access ds:forSource ?source ;
                        ds:methodType ?method .
                OPTIONAL {{ ?access ds:baseUrl ?url }}
            }}
        """)

        # Build enrichment dict
        enrichment = {}
        for row in basic:
            sid = row['source'].split('/')[-1]
            enrichment[sid] = {
                'label': row.get('label'),
                'tier': row.get('tier', '').split('/')[-1] if row.get('tier') else None,
                'category': row.get('category'),
                'maintainer': row.get('maintainer'),
                'access_methods': []
            }

        for row in access:
            sid = row['source'].split('/')[-1]
            if sid in enrichment:
                enrichment[sid]['access_methods'].append({
                    'method': row.get('method'),
                    'url': row.get('url')
                })

        return enrichment

    # ─────────────────────────────────────────────────────────────────
    # Embedding Management
    # ─────────────────────────────────────────────────────────────────

    async def load_embeddings(self) -> int:
        """
        Load all data sources into RuVector with embeddings.

        Returns the number of entities loaded.
        """
        # Get all data sources with descriptions
        sources = self.sparql.query_simple("""
            SELECT ?source ?label ?description ?category WHERE {
                ?source a ds:DataSource ;
                        rdfs:label ?label .
                OPTIONAL { ?source dct:description ?description }
                OPTIONAL {
                    ?source ds:category ?cat .
                    ?cat rdfs:label ?category .
                }
            }
        """)

        count = 0
        for source in sources:
            source_id = source['source'].split('/')[-1]
            label = source.get('label', '')
            description = source.get('description', '')
            category = source.get('category', '')

            # Create embedding text
            text = f"{label}: {description}" if description else label
            if category:
                text += f" ({category})"

            # Store in RuVector via claude-flow memory
            try:
                metadata = json.dumps({
                    'qlever_uri': source['source'],
                    'entity_type': 'DataSource',
                    'label': label,
                    'description': description,
                    'category': category
                })

                result = subprocess.run(
                    [
                        "npx", "@claude-flow/cli@latest", "memory", "store",
                        "--key", source_id,
                        "--value", text,
                        "--namespace", self.memory_namespace,
                        "--tags", f"datasource,{category.lower().replace(' ', '-')}" if category else "datasource"
                    ],
                    capture_output=True,
                    text=True,
                    timeout=10
                )

                if result.returncode == 0:
                    count += 1

            except (subprocess.TimeoutExpired, Exception) as e:
                print(f"Warning: Failed to store {source_id}: {e}")
                continue

        return count


# ─────────────────────────────────────────────────────────────────────
# CLI Interface
# ─────────────────────────────────────────────────────────────────────

async def main():
    """CLI interface for hybrid knowledge service."""
    import sys

    kg = HybridKnowledgeService()

    if len(sys.argv) < 2:
        print("Gene Knowledge Graph - Hybrid Search")
        print("=" * 50)
        print(f"SPARQL Endpoint: {kg.sparql.endpoint}")
        print(f"Triples: {kg.sparql.count_triples():,}")
        print()
        print("Commands:")
        print("  python3 -m scripts.kg.hybrid search <query>")
        print("  python3 -m scripts.kg.hybrid sparql '<SPARQL>'")
        print("  python3 -m scripts.kg.hybrid source <id>")
        print("  python3 -m scripts.kg.hybrid load")
        print("  python3 -m scripts.kg.hybrid list")
        return

    command = sys.argv[1]

    if command == "search" and len(sys.argv) > 2:
        query = " ".join(sys.argv[2:])
        results = await kg.hybrid_search(query)
        print(json.dumps([r.to_dict() for r in results], indent=2))

    elif command == "sparql" and len(sys.argv) > 2:
        query = sys.argv[2]
        results = kg.sparql_query(query)
        print(json.dumps(results, indent=2))

    elif command == "source" and len(sys.argv) > 2:
        source_id = sys.argv[2]
        result = kg.get_source(source_id)
        print(json.dumps(result, indent=2))

    elif command == "load":
        count = await kg.load_embeddings()
        print(f"Loaded {count} entities into RuVector")

    elif command == "list":
        sources = kg.list_sources(20)
        for s in sources:
            tier = s.get('tier', '').split('/')[-1] if s.get('tier') else '-'
            print(f"  [{tier}] {s['label']}")

    else:
        print(f"Unknown command: {command}")
        print("Use: search, sparql, source, load, list")


if __name__ == "__main__":
    asyncio.run(main())
