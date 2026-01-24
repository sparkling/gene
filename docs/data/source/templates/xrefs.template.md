# {{SOURCE_NAME}} Cross-References

## Identifier Mappings

### Primary Identifiers

| This Source | External Database | External ID | Coverage |
|-------------|-------------------|-------------|----------|
| {{LOCAL_ID_1}} | {{EXTERNAL_DB_1}} | {{EXTERNAL_ID_1}} | {{HIGH/MEDIUM/LOW}} |
| {{LOCAL_ID_2}} | {{EXTERNAL_DB_2}} | {{EXTERNAL_ID_2}} | {{HIGH/MEDIUM/LOW}} |
| {{LOCAL_ID_3}} | {{EXTERNAL_DB_3}} | {{EXTERNAL_ID_3}} | {{HIGH/MEDIUM/LOW}} |

### Coverage by Database

| External DB | Records Mapped | Coverage % | Notes |
|-------------|----------------|------------|-------|
| {{DB_1}} | {{COUNT}} | {{PERCENT}} | {{NOTES}} |
| {{DB_2}} | {{COUNT}} | {{PERCENT}} | {{NOTES}} |

### Mapping Files

| Mapping | File/Endpoint | Format |
|---------|---------------|--------|
| {{MAPPING_1}} | {{LOCATION_1}} | {{FORMAT_1}} |
| {{MAPPING_2}} | {{LOCATION_2}} | {{FORMAT_2}} |

---

## Integration Pairs

| Pair With | Purpose | Join Key | Recommended |
|-----------|---------|----------|-------------|
| {{SOURCE_1}} | {{PURPOSE}} | {{KEY}} | {{Yes/No}} |
| {{SOURCE_2}} | {{PURPOSE}} | {{KEY}} | {{Yes/No}} |

---

## Related Data Sources

### Same Category

| Source | Relationship | Use Case |
|--------|--------------|----------|
| [{{SOURCE_1}}]({{PATH_1}}) | {{RELATIONSHIP}} | {{USE_CASE}} |
| [{{SOURCE_2}}]({{PATH_2}}) | {{RELATIONSHIP}} | {{USE_CASE}} |

### Complementary Sources

| Source | What It Adds | Integration Key |
|--------|--------------|-----------------|
| {{COMP_SOURCE_1}} | {{VALUE_ADD_1}} | {{KEY_1}} |
| {{COMP_SOURCE_2}} | {{VALUE_ADD_2}} | {{KEY_2}} |

---

## Integration Workflow

1. {{STEP_1: e.g., "Query source by primary ID"}}
2. {{STEP_2: e.g., "Map to external ID using xref table"}}
3. {{STEP_3: e.g., "Fetch enrichment data from external source"}}

---

## Integration Examples

### Join with {{RELATED_SOURCE}}

```sql
-- Example join query
{{JOIN_QUERY_EXAMPLE}}
```

### API Cross-Reference

```bash
# Fetch related data from {{EXTERNAL_API}}
{{API_XREF_EXAMPLE}}
```

---

## Notes

{{INTEGRATION_NOTES}}
