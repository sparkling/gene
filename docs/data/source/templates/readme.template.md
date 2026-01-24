---
id: {{source-id}}
title: "{{Source Name}}"
type: source
parent: ../README.md
category: {{category}}
subcategory: {{subcategory}}
tier: {{1|2|3}}
status: {{draft|active|deprecated}}
last_updated: {{YYYY-MM-DD}}
tags:
  - {{tag1}}
  - {{tag2}}
---

# {{Source Name}}

{{2-3 sentence description of what this data source is, who maintains it, and its primary purpose.}}

## Overview

| Property | Value |
|----------|-------|
| **Maintainer** | {{Organization}} |
| **Website** | {{https://example.com}} |
| **Update Frequency** | {{Daily/Weekly/Monthly/Quarterly/Annually}} |
| **Records** | {{X,XXX,XXX}} |
| **Latest Release** | {{version}} ({{YYYY-MM-DD}}) |

## Primary Use Cases

1. {{Use case 1: Brief description}}
2. {{Use case 2: Brief description}}
3. {{Use case 3: Brief description}}

## Key Identifiers

| Identifier | Format | Example |
|------------|--------|---------|
| {{Primary ID}} | `{{pattern}}` | `{{example}}` |
| {{Secondary ID}} | `{{pattern}}` | `{{example}}` |

## Limitations

- {{Limitation 1: e.g., "Limited to human variants only"}}
- {{Limitation 2: e.g., "No frequency data for rare populations"}}
- {{Limitation 3: e.g., "Updates may lag behind primary literature"}}

## Data Quality Notes

{{Brief notes about data quality, curation process, known issues, or coverage gaps.}}

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [Example Records](./sample.json) - Sample data for code generation
- [Download Guide](./download.md) - Access methods and API configuration
- [License Terms](./license.md) - Usage rights and restrictions
- [Schema Mapping](./mapping.xslt) - XSLT 3.0 transformation to unified schema
