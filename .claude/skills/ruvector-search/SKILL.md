---
name: "RuVector Semantic Search"
description: "Perform semantic search over vector embeddings using natural language. Use when searching for similar concepts, finding related entities, or when structured queries need semantic enrichment. Triggers: 'semantic search', 'similar to', 'find related', 'embedding search', 'vector search'."
---

# RuVector Semantic Search

## What This Skill Does

Provides semantic search capabilities over vector embeddings stored in RuVector/AgentDB. This enables:
1. Finding semantically similar entities even with different terminology
2. Enriching structured queries with semantic context
3. Ranking results by semantic relevance
4. Cross-referencing between graph data and embeddings

## Prerequisites

- RuVector/AgentDB with HNSW indexing enabled
- Embeddings generated for target corpus
- claude-flow CLI with embeddings support

## Quick Start

### Basic Semantic Search

```
Find concepts similar to "cancer treatment databases"
```

This skill will:
1. Generate an embedding for your query
2. Search the vector index for similar embeddings
3. Return ranked results with similarity scores

### Hybrid Search (Graph + Semantic)

Combine with structured queries for best results:
1. Semantic search finds relevant entities
2. Graph queries filter and enrich results
3. Results are merged and deduplicated

---

## How It Works

### 1. Query Embedding Generation

The skill converts natural language queries into vector embeddings using the configured embedding model.

### 2. HNSW Index Search

Searches the Hierarchical Navigable Small World (HNSW) index for approximate nearest neighbors. Performance: 150x-12,500x faster than brute force.

### 3. Result Ranking

Results are ranked by cosine similarity. Default threshold: 0.7 (adjustable).

### 4. Entity Resolution

Maps embedding results back to knowledge graph entities using stored identifiers.

---

## Integration Patterns

### For Skill Developers

This skill provides semantic search as a building block. To use in your skill:

1. **Intent Detection**: Use embeddings to classify user intent
2. **Entity Disambiguation**: Find the most relevant entity when names are ambiguous
3. **Query Expansion**: Add semantically related terms to improve recall
4. **Result Reranking**: Reorder graph query results by semantic relevance

### API Interface

The skill uses claude-flow CLI commands internally:

```bash
# Search embeddings (internal - not exposed to users)
npx @claude-flow/cli@latest embeddings search --query "..." --top-k 10 --threshold 0.7

# Generate embedding for comparison
npx @claude-flow/cli@latest embeddings generate --text "..."

# Compare two texts semantically
npx @claude-flow/cli@latest embeddings compare --text1 "..." --text2 "..."
```

---

## Configuration

### Embedding Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| top_k | 10 | Number of results to return |
| threshold | 0.7 | Minimum similarity score |
| namespace | default | Embedding namespace to search |

### Supported Namespaces

- `datasources` - Data source descriptions and metadata
- `categories` - Category and subcategory concepts
- `documentation` - Full documentation text
- `patterns` - Learned patterns from previous queries

---

## Output Format

Results are returned as structured data:

```
Found 5 semantically similar results:

1. **ClinVar** (similarity: 0.92)
   Category: Genetics and Genomics > Variant Repositories

2. **dbSNP** (similarity: 0.87)
   Category: Genetics and Genomics > Variant Repositories

3. **COSMIC** (similarity: 0.84)
   Category: Genetics and Genomics > Cancer Genomics
```

---

## Troubleshooting

### No Results Found

- Lower the similarity threshold (try 0.5)
- Check that embeddings exist for the target namespace
- Verify the embedding index is initialized

### Irrelevant Results

- Increase the similarity threshold (try 0.8+)
- Add more specific terms to your query
- Use hybrid search with graph constraints

### Performance Issues

- Ensure HNSW index is built (not linear search)
- Check index parameters (ef_search, M values)
- Consider reducing top_k for faster results

---

## Related Skills

- [sparql](../sparql/) - Structured graph queries
- [qlever](../qlever/) - High-performance SPARQL engine
- [datasource](../datasource/) - Unified data source queries
