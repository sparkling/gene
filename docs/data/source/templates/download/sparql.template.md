## SPARQL

### Endpoint

```
{{SPARQL_ENDPOINT}}
```

### Prefixes

```sparql
PREFIX {{prefix1}}: <{{uri1}}>
PREFIX {{prefix2}}: <{{uri2}}>
PREFIX {{prefix3}}: <{{uri3}}>
```

### Limits

| Limit | Value |
|-------|-------|
| Timeout | {{N}} seconds |
| Max Results | {{N}} rows |
| Max Query Size | {{N}} KB |

### Example Query

```sparql
{{EXAMPLE_QUERY}}
```

### Common Patterns

#### Pattern 1: {{PATTERN_NAME}}

```sparql
{{PATTERN_QUERY}}
```

#### Pattern 2: {{PATTERN_NAME}}

```sparql
{{PATTERN_QUERY}}
```

### Usage via curl

```bash
curl -X POST "{{SPARQL_ENDPOINT}}" \
  -H "Accept: application/sparql-results+json" \
  --data-urlencode "query={{QUERY}}"
```
