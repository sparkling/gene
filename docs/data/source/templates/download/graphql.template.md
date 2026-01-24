## GraphQL API

### Endpoint

```
{{GRAPHQL_ENDPOINT}}
```

### Authentication

| Requirement | Details |
|-------------|---------|
| API Key | {{Required/Optional/None}} |
| Header | `Authorization: {{Bearer}} {{TOKEN}}` |

### Rate Limits

| Limit | Value |
|-------|-------|
| Requests/minute | {{N}} |
| Query complexity | {{N}} |

### Schema Explorer

{{SCHEMA_EXPLORER_URL}}

### Example Query

```graphql
{{EXAMPLE_QUERY}}
```

### Example Variables

```json
{{EXAMPLE_VARIABLES}}
```

### Usage via curl

```bash
curl -X POST "{{GRAPHQL_ENDPOINT}}" \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer {{TOKEN}}" \
  -d '{"query": "{{QUERY}}", "variables": {{VARIABLES}}}'
```

### Common Queries

#### Query 1: {{QUERY_NAME}}

```graphql
{{COMMON_QUERY_1}}
```

#### Query 2: {{QUERY_NAME}}

```graphql
{{COMMON_QUERY_2}}
```
