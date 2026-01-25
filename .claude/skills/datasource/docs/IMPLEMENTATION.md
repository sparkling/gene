# Data Source Skill Implementation Guide

This document describes how the datasource skill works internally and how it orchestrates other skills.

## Architecture Overview

```
User Query (Natural Language)
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│                    datasource skill                          │
│  ┌─────────────┐  ┌──────────────┐  ┌─────────────────────┐ │
│  │   Intent    │  │    Query     │  │      Result         │ │
│  │ Classifier  │──▶│  Generator   │──▶│    Formatter       │ │
│  └─────────────┘  └──────────────┘  └─────────────────────┘ │
│         │                │                     │             │
└─────────┼────────────────┼─────────────────────┼─────────────┘
          │                │                     │
          ▼                ▼                     ▼
   ┌────────────┐  ┌─────────────┐      ┌─────────────┐
   │ ruvector-  │  │   sparql    │      │    owl/     │
   │  search    │  │   + qlever  │      │ skos/shacl  │
   └────────────┘  └─────────────┘      └─────────────┘
   (semantic)       (structured)         (schema)
```

## Execution Flow

### Phase 1: Intent Classification

The LLM analyzes the user query to determine:

1. **Query Type**: browse, search, details, compare, stats, related
2. **Entities**: Category names, source IDs, properties mentioned
3. **Filters**: Tier, license, access method, format
4. **Output Mode**: List, detail, comparison table, statistics

### Phase 2: Context Loading

Based on the intent, load relevant context:

1. **Ontology Context** (`resources/ontology-context.md`)
   - Prefixes, classes, properties
   - Taxonomy concepts for filtering

2. **Query Examples** (`resources/query-guide.md`)
   - Similar query patterns as few-shot examples

3. **Semantic Enrichment** (via ruvector-search)
   - If query is ambiguous, use embeddings to clarify
   - Find semantically related terms to expand search

### Phase 3: SPARQL Generation

The LLM generates a SPARQL query using:

1. **Schema Knowledge**: OWL class/property definitions
2. **Vocabulary**: SKOS concept labels and hierarchy
3. **Constraints**: SHACL shape requirements
4. **Examples**: Few-shot patterns from query guide

The query is generated dynamically, NOT from templates.

### Phase 4: Query Execution

Execute the generated SPARQL via qlever:

```bash
# Internal execution (not exposed to user)
curl -s "http://localhost:7001/api?query=$(urlencode "$SPARQL")"
```

Or via the qlever skill for managed execution.

### Phase 5: Result Formatting

Transform raw results into user-friendly output:

1. **Remove technical details**: No URIs, prefixes, or SPARQL
2. **Apply structure**: Tables, lists, or detailed cards
3. **Sort appropriately**: By tier, alphabetically, or by relevance
4. **Limit output**: Pagination for large result sets
5. **Add context**: Category breadcrumbs, related info

## Skill Invocation Patterns

### Using the sparql Skill

When detailed SPARQL optimization is needed:

```
/sparql optimize query against qlever for performance
```

The sparql skill provides:
- Query validation
- Performance optimization hints
- QLever-specific dialect adjustments

### Using the qlever Skill

For query execution and endpoint management:

```
/qlever execute query against datasource endpoint
```

The qlever skill handles:
- Endpoint configuration
- Connection management
- Result streaming
- Error handling

### Using the ruvector-search Skill

For semantic search integration:

```
/ruvector-search find similar concepts to "cancer treatment"
```

Use cases:
- Disambiguate user queries
- Expand search with related terms
- Rank results by semantic relevance
- Find entities when exact match fails

### Using owl/skos/shacl Skills

For schema understanding:

```
/owl explain DataSource class properties
/skos show taxonomy hierarchy for categories
/shacl validate query constraints
```

These skills provide:
- Schema documentation
- Property cardinality info
- Valid value ranges
- Hierarchy navigation

## Hybrid Query Strategy

For complex queries, combine structured and semantic approaches:

```
1. User: "Find sources similar to UniProt for protein analysis"

2. Semantic Phase:
   - Generate embedding for "similar to UniProt for protein analysis"
   - Search ruvector for nearest neighbors
   - Get candidate sources: [PDB, InterPro, PFAM, ...]

3. Structured Phase:
   - Generate SPARQL to get full details for candidates
   - Filter by additional criteria (tier, license, etc.)

4. Merge Results:
   - Rank by semantic similarity + tier priority
   - Deduplicate
   - Format for user
```

## Error Handling

### Query Generation Errors

If SPARQL generation fails:
1. Fall back to broader query
2. Use semantic search as backup
3. Ask user for clarification

### Execution Errors

If qlever execution fails:
1. Check endpoint availability
2. Simplify query (remove OPTIONALs)
3. Return partial results with explanation

### No Results

If query returns empty:
1. Relax filters
2. Suggest alternative queries
3. Use semantic search to find similar

## Performance Optimization

### Caching

- Cache taxonomy lookups (categories, tiers, etc.)
- Cache frequent query results
- Use semantic cache for similar queries

### Query Optimization

- Limit results early in query
- Use indexed properties for filtering
- Avoid unnecessary JOINs

### Context Management

- Load minimal context for simple queries
- Progressive disclosure of schema info
- Batch multiple queries when possible

## Configuration

### Endpoints

| Service | Default | Environment Variable |
|---------|---------|---------------------|
| QLever | localhost:7001 | `QLEVER_ENDPOINT` |
| RuVector | localhost:7002 | `RUVECTOR_ENDPOINT` |

### Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| max_results | 50 | Maximum results to return |
| semantic_threshold | 0.7 | Minimum similarity score |
| query_timeout | 30s | SPARQL execution timeout |

## Testing

### Unit Tests

Test individual components:
- Intent classification accuracy
- SPARQL generation correctness
- Result formatting consistency

### Integration Tests

Test full flow:
- Natural language → formatted results
- Error handling scenarios
- Edge cases (no results, timeouts)

### Regression Tests

Ensure consistent behavior:
- Same query → same results
- Performance benchmarks
- Semantic search quality

## Debugging

### Verbose Mode

Enable verbose output for troubleshooting:
- Show classified intent
- Show generated SPARQL (for debugging only)
- Show raw results before formatting

### Logging

Log key events:
- Query received
- Intent classified
- SPARQL generated
- Results returned
- Errors encountered

## Security Considerations

1. **No SPARQL injection**: LLM-generated queries are validated
2. **No credential exposure**: API keys handled by qlever skill
3. **Rate limiting**: Respect endpoint rate limits
4. **Input sanitization**: User input is escaped in queries
