# Migration Templates - Quick Reference

**Purpose:** Copy-paste templates for creating data source documentation

---

## Template 1: _index.md (Primary File)

```markdown
---
id: {source.name}
title: "{Full Source Name}"
type: data-source
category: {genetics|traditional|nutrition|shared}
subcategory: {tcm|ayurveda|kampo|western-herbal|global}
parent: ../_index.md
tier: {1|2|3}
last_updated: YYYY-MM-DD
status: {active|planned|deprecated}
tags: [{tag1}, {tag2}, {tag3}]
---

# {Source Name}

**Category:** [{Category}](../../_index.md) > [{Subcategory}](../_index.md)

## Overview

{2-3 paragraph description}

## Key Statistics

| Metric | Value |
|--------|-------|
| Records | |
| Last Update | |
| Coverage | |
| Storage | |

## Primary Use Cases

1. {Use case 1}
2. {Use case 2}
3. {Use case 3}

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| | | |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| | | |

## License

| Aspect | Value |
|--------|-------|
| License | |
| Commercial Use | {Yes/No/Restricted} |
| Attribution | {Required/Optional} |

## See Also

- [Schema Documentation](./schema.md)
- [Download Instructions](./download.md)
```

---

## Template 2: schema.md (Technical Schema)

```markdown
---
id: schema-{source.name}
title: "{Source Name} Schema Documentation"
type: schema
parent: _index.md
last_updated: YYYY-MM-DD
status: {draft|migrated|final}
tags: [schema, database, {domain-tags}]
---

# {Source Name} Schema Documentation

**Document ID:** SCHEMA-{SOURCE-ID}
**Version:** {version}
**Source Version:** {upstream version/date}

---

## TL;DR

{2-3 sentence summary}

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| | | |

---

## Entity Relationship Overview

```
{ASCII diagram}
```

---

## Core Tables/Entities

### {Entity 1}

**Description:** {purpose}

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| | | | |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| | | |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | |
| Alternative | |
| Encoding | UTF-8 |

---

## Sample Record

```json
{
  "example": "record"
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| | |

---

## References

1. {URL}
```

---

## Template 3: download.md (Acquisition Guide)

```markdown
---
id: download-{source.name}
title: "{Source Name} Download Instructions"
type: download
parent: _index.md
last_updated: YYYY-MM-DD
---

# {Source Name} Download Instructions

## Quick Start

```bash
# Minimal download command
```

## Prerequisites

- {Requirement}

## Download Methods

### Primary: {Method Name}

```bash
# Commands
```

### Alternative: {Method Name}

```bash
# Commands
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| | | |

## Post-Download Processing

```bash
# Processing commands
```

## Verification

```bash
# Validation commands
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| | |
```

---

## Template 4: xrefs.md (Cross-References)

```markdown
---
id: xrefs-{source.name}
title: "{Source Name} Cross-References"
type: xrefs
parent: _index.md
last_updated: YYYY-MM-DD
---

# {Source Name} Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `{path}` |
| Secondary | Symlink | `{path}` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| | | |

## Integration Notes

{Notes on how this source integrates with others}
```

---

## Template 5: Subcategory _index.md

```markdown
---
id: {subcategory.id}
title: "{Subcategory Name}"
type: subcategory
parent: ../_index.md
last_updated: YYYY-MM-DD
---

# {Subcategory Name}

**Parent:** [{Category Name}](../_index.md)

## Overview

{Description}

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [{Name}](./{folder}/_index.md) | | |

## Integration Notes

{Integration guidance}
```

---

## Naming Quick Reference

### Folder Names
- Lowercase only
- Use dots (`.`) as separators
- No hyphens, underscores, or spaces
- Examples: `dbsnp`, `gwas.catalog`, `open.food.facts`

### File Names
- Standard: `_index.md`, `schema.md`, `download.md`, `xrefs.md`
- Custom: `{descriptor}.{type}.md` (e.g., `api.endpoints.md`)

### IDs in Frontmatter
- Data source: `{folder.name}` (e.g., `dbsnp`)
- Schema: `schema-{folder.name}` (e.g., `schema-dbsnp`)
- Download: `download-{folder.name}`
- Cross-refs: `xrefs-{folder.name}`

---

## Checklist Per Source

- [ ] Create folder with correct naming
- [ ] Create `_index.md` with overview
- [ ] Create `schema.md` with technical details
- [ ] Create `download.md` if non-trivial acquisition
- [ ] Create `xrefs.md` if polyhierarchical (has symlinks)
- [ ] Verify all frontmatter fields
- [ ] Verify all internal links
- [ ] Update subcategory `_index.md` to list source
