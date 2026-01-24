# CLAUDE.md Template Guide

## Purpose

The CLAUDE.md template creates **LLM-optimized** data source references that:
- Load efficiently into AI context windows
- Provide actionable technical information
- Enable quick API/query construction
- Support data integration decisions

## Design Principles

1. **Concise**: 50-100 lines max (vs 500+ lines in full docs)
2. **Structured**: Tables, lists, code blocks - not prose
3. **Actionable**: API endpoints, query examples, sample data
4. **Cross-ref focused**: Identifiers for integration

## Template Sections

### Quick Reference Table (Required)
The most important metadata in scannable format.

```markdown
| Field | Value |
|-------|-------|
| **URL** | Primary access point |
| **Maintainer** | Who maintains it |
| **License** | License name + type |
| **Commercial OK** | Yes/No/Conditional |
| **Update Freq** | Daily/Weekly/Monthly/Quarterly |
| **Version** | Current version |
```

### Content Summary (Required)
2-3 sentences maximum describing:
- What the database contains
- Primary entity type
- Key differentiator

### Record Counts (Required)
Bullet list of primary entities with counts:
```markdown
- Primary entities: X
- Secondary entities: Y
- Tertiary entities: Z
```

### Key Identifiers (Required)
Table of ID types for cross-referencing:
```markdown
| ID Type | Format | Example |
|---------|--------|---------|
| Primary ID | format pattern | concrete example |
```

Plus comma-separated list of cross-reference databases.

### Core Schema (Required)
Only the 5-7 most important fields of the primary entity.
Include relationship diagram using ASCII:
```
EntityA --[relation]--> EntityB
```

### Access Methods (Required)
- API endpoint with example URL
- Rate limits
- Bulk download URL
- File formats and sizes

### Query Examples (Required)
2 practical examples showing common use cases.
Use appropriate query language (SQL, GraphQL, REST, etc.)

### Sample Record (Required)
One representative JSON record showing key fields.
Keep to ~15 lines max.

### Integration Notes (Required)
Four bullet points:
- **Primary use:** Main use case
- **Best for:** Strengths
- **Limitations:** Known issues/gaps
- **Pairs with:** Complementary sources

## What to Exclude

- Historical background/context
- Detailed schema (all fields)
- Full glossary
- Change logs
- Methodology explanations
- License full text
- Detailed API documentation

## Naming Convention

```
{source-name}-CLAUDE.md
```

Examples:
- `chembl-CLAUDE.md`
- `gnomad-CLAUDE.md`
- `uniprot-idmapping-CLAUDE.md`

## Quality Checklist

- [ ] Under 100 lines
- [ ] All required sections present
- [ ] Working example queries
- [ ] Valid sample JSON
- [ ] Cross-references listed
- [ ] License clearly stated
- [ ] API endpoint included
- [ ] Download URL included

## Example Files

See `/docs/data/source.new/resource/examples/`:
- `chembl-CLAUDE.md` - Bioactivity database
- `gnomad-CLAUDE.md` - Population genetics
- `clinvar-CLAUDE.md` - Clinical variants
