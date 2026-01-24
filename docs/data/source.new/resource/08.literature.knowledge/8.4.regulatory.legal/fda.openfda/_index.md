---
id: fda.openfda
title: "FDA openFDA"
type: data-source
category: literature
subcategory: regulatory.legal
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [fda, regulatory, drugs, devices, adverse-events, recalls]
---

# FDA openFDA

**Category:** [Literature](../../../_index.md) > [Regulatory & Legal](../_index.md)

## Overview

openFDA provides open access to FDA public data through a modern, searchable API. It includes data on drugs (adverse events, labeling, NDC codes), devices (adverse events, recalls, 510(k)s, PMA), and foods (adverse events, recalls).

The API provides Elasticsearch-powered search across millions of records, enabling sophisticated queries for pharmacovigilance, regulatory research, and drug safety analysis. Data is harmonized and enriched with standard identifiers.

openFDA is essential for adverse event analysis, drug safety research, and regulatory intelligence.

## Key Statistics

| Metric | Value |
|--------|-------|
| Drug Adverse Events | 20M+ |
| Device Adverse Events | 15M+ |
| Drug Labels | 150,000+ |
| Food Recalls | 25,000+ |
| Last Update | Weekly |

## Primary Use Cases

1. Adverse event signal detection
2. Drug safety analysis
3. Recall monitoring
4. Label information retrieval
5. Regulatory timeline research

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| NDC Code | `[0-9]{4,5}-[0-9]{3,4}-[0-9]{1,2}` | 0002-1200-01 |
| Application Number | `[AN]DA[0-9]+` | NDA019501 |
| Device ID | Various | K123456 (510(k)) |
| RxNorm | Numeric | 310965 |
| UNII | Alphanumeric | 362O9ITL9D |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| API | https://api.fda.gov | REST + Elasticsearch |
| Web Interface | https://open.fda.gov | Documentation site |
| Downloads | https://open.fda.gov/apis/downloads | Bulk data |
| Interactive | https://open.fda.gov/apis/drug | Try queries |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (US Government) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Documentation](./schema.md)
- [ClinicalTrials.gov](../clinicaltrials.gov/_index.md) - Trial data
- [DailyMed](../dailymed/) - Drug labeling (symlink)
