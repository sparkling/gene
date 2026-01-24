---
id: resource-claude-md
title: CLAUDE.md Templates
description: LLM-optimized data source reference templates
last_updated: 2026-01-23
---

# CLAUDE.md Templates

LLM-optimized reference files for data sources, designed for efficient context loading.

## Files

| File | Description |
|------|-------------|
| [CLAUDE.md.template](./CLAUDE.md.template) | Template with placeholder variables |
| [TEMPLATE-GUIDE.md](./TEMPLATE-GUIDE.md) | Instructions for filling out templates |

## Examples

| File | Data Source | Category |
|------|-------------|----------|
| [chembl-CLAUDE.md](./examples/chembl-CLAUDE.md) | ChEMBL | Bioactivity/Drugs |
| [gnomad-CLAUDE.md](./examples/gnomad-CLAUDE.md) | gnomAD | Population Genetics |
| [clinvar-CLAUDE.md](./examples/clinvar-CLAUDE.md) | ClinVar | Clinical Variants |

## Purpose

These concise (~50-100 line) reference files provide:

1. **Quick metadata** - URL, license, update frequency
2. **Key identifiers** - ID formats for cross-referencing
3. **Core schema** - Essential fields only
4. **Access methods** - API endpoints, downloads
5. **Query examples** - Practical usage patterns
6. **Integration notes** - When to use, what pairs with

## Design Goals

- **Token-efficient**: Minimal prose, maximum structure
- **Actionable**: Working queries and endpoints
- **Integration-focused**: Cross-reference identifiers prominent
- **Scannable**: Tables and lists over paragraphs
