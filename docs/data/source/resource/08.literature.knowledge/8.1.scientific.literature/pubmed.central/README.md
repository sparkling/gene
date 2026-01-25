---
id: pubmed.central
title: "PubMed Central (PMC)"
type: data-source
category: literature
subcategory: scientific.literature
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [literature, full-text, open-access, ncbi, archive]
---

# PubMed Central (PMC)

**Category:** [Literature](../../../README.md) > [Scientific Literature](../README.md)

## Overview

PubMed Central (PMC) is a free full-text archive of biomedical and life sciences journal literature at the U.S. National Institutes of Health's National Library of Medicine (NIH/NLM). It serves as a permanent repository for NIH-funded research and participating publishers' content.

PMC provides XML-formatted full text enabling sophisticated text mining and data extraction. It includes supplementary data, figures, and tables not available in abstracts. PMC also hosts author manuscripts submitted in compliance with NIH public access policy.

PMC is essential for text mining research, systematic reviews requiring full text, and accessing open access biomedical literature.

## Key Statistics

| Metric | Value |
|--------|-------|
| Full-Text Articles | 9M+ |
| Journals | 3,000+ |
| Author Manuscripts | 1M+ |
| New Articles/Year | 500,000+ |
| Last Update | Daily |

## Primary Use Cases

1. Full-text literature mining
2. Systematic review data extraction
3. Figure and table extraction
4. Supplementary data access
5. Open access research

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| PMCID | `PMC[0-9]+` | PMC1234567 |
| PMID | Numeric | 12345678 |
| DOI | Standard | 10.1000/xyz123 |
| Manuscript ID | `NIHMS[0-9]+` | NIHMS123456 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.ncbi.nlm.nih.gov/pmc | Search and read |
| OA Subset | https://www.ncbi.nlm.nih.gov/pmc/tools/openftlist | Commercial mining OK |
| FTP | https://ftp.ncbi.nlm.nih.gov/pub/pmc | Bulk downloads |
| E-utilities | https://eutils.ncbi.nlm.nih.gov | Programmatic access |
| PMC ID Converter | https://www.ncbi.nlm.nih.gov/pmc/utils/idconv | ID mapping |

## License

| Aspect | Value |
|--------|-------|
| License | Varies by article |
| Commercial Use | OA Subset only |
| Attribution | Required (per article license) |

## See Also

- [Schema Documentation](./schema.md)
- [PubMed](../pubmed/README.md) - Abstract database
- [PMC ID Converter](../../8.3.identifier.mapping/pmc.id.converter/README.md) - ID conversion
