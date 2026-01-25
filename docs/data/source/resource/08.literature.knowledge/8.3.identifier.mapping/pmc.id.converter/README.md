---
id: pmc.id.converter
title: "PMC ID Converter"
type: data-source
category: literature
subcategory: identifier.mapping
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [identifiers, pmid, pmcid, doi, conversion]
---

# PMC ID Converter

**Category:** [Literature](../../../README.md) > [Identifier Mapping](../README.md)

## Overview

The PMC ID Converter is a web service provided by NCBI that converts between different publication identifier types: PubMed IDs (PMID), PubMed Central IDs (PMCID), and DOIs. It supports batch conversion of up to 200 IDs per request.

The service is essential for literature data integration when records from different sources use different identifier types. It also provides publication metadata including release dates and versioning information.

PMC ID Converter is widely used in systematic reviews, bibliography management, and literature database integration.

## Key Statistics

| Metric | Value |
|--------|-------|
| ID Types | PMID, PMCID, DOI |
| Batch Size | Up to 200 IDs |
| Response Format | XML, JSON, CSV |
| Coverage | PMC articles |
| Last Update | Real-time |

## Primary Use Cases

1. Converting between PMID/PMCID/DOI
2. Bibliography standardization
3. Literature database integration
4. Systematic review workflows
5. Reference management

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| PMID | Numeric | 12345678 |
| PMCID | `PMC[0-9]+` | PMC1234567 |
| DOI | Standard | 10.1000/xyz123 |
| Manuscript ID | `MID:NIHMS[0-9]+` | MID:NIHMS123456 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0 | Interactive |
| REST API | https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/ | Programmatic |
| Documentation | https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api | Full guide |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (US Government) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Documentation](./schema.md)
- [PubMed](../../8.1.scientific.literature/pubmed/README.md) - Abstract database
- [PubMed Central](../../8.1.scientific.literature/pubmed.central/README.md) - Full-text archive
