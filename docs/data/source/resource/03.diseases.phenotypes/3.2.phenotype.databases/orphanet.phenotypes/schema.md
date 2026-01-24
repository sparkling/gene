---
id: schema-orphanet-phenotypes
title: "Orphanet Phenotype Annotations Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: draft
tags: [schema, phenotypes, rare-diseases, hpo, clinical-features, orphanet]
---

# Orphanet Phenotype Annotations Schema Documentation

**Document ID:** SCHEMA-ORPHANET-PHENOTYPES

---

## TL;DR

Orphanet Phenotype Annotations provides comprehensive mappings between rare diseases and their associated clinical features using HPO terms. Each association includes frequency information (obligate to very rare) and validation status, curated by clinical experts across 41 countries.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Annotated Diseases | 4,337+ | Orphanet |
| HPO Annotations | 100,000+ | Orphanet |
| Frequency Categories | 6 | Orphanet |
| Update Frequency | Monthly | Orphanet |

---

## Data Format

| Aspect | Value |
|--------|-------|
| Primary Format | XML, JSON, RDF |
| API | Yes (REST API) |

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Orphanet ID | `Orphanet:[0-9]+` | Orphanet:558 |
| ORPHA Code | `[0-9]+` | 558 (Marfan) |
| HPO Term | `HP:[0-9]{7}` | HP:0001166 |
| Frequency | HPO frequency terms | HP:0040281 (Very frequent) |

---

## Frequency Categories

| Category | HPO Term | Percentage |
|----------|----------|------------|
| Obligate | HP:0040280 | 100% |
| Very frequent | HP:0040281 | 80-99% |
| Frequent | HP:0040282 | 30-79% |
| Occasional | HP:0040283 | 5-29% |
| Very rare | HP:0040284 | 1-4% |
| Excluded | HP:0040285 | 0% |

---

## References

See [Overview](./_index.md) for full details.
