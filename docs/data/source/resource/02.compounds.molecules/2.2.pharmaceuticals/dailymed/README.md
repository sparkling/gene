---
id: dailymed
title: "DailyMed - FDA Drug Labeling"
type: source
parent: ../README.md
tier: 2
status: active
category: compounds.molecules
subcategory: pharmaceuticals
tags:
  - fda
  - drug-labels
  - prescribing-information
  - spl
  - regulatory
---

# DailyMed - FDA Drug Labeling

## Overview

DailyMed provides high-quality information about marketed drugs in the United States, including FDA-approved prescription and over-the-counter drug labeling (package inserts). The service is provided by the National Library of Medicine (NLM) in coordination with the FDA.

DailyMed uses Structured Product Labeling (SPL) format, an XML-based standard that enables structured representation of drug label content. This makes it possible to extract and analyze specific sections of drug labels programmatically, supporting clinical decision support systems and drug information databases.

The database includes drug package inserts, patient information, and medication guides, making it the authoritative source for current FDA-approved drug labeling in the US market.

## Key Statistics

| Metric | Value |
|--------|-------|
| Drug Labels | 140,000+ |
| Active Labels | 100,000+ |
| Drug Products | 45,000+ |
| NDC Numbers | 200,000+ |
| Updates | Daily |

## Primary Use Cases

1. **Drug Label Lookup** - Access official FDA-approved labeling
2. **Clinical Reference** - Prescribing information for healthcare providers
3. **Regulatory Compliance** - Verify current approved indications
4. **Data Integration** - SPL-formatted data for health IT systems
5. **Label History** - Track changes in drug labeling over time

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Set ID (SPL) | UUID | 0000a193-4512-468c-93b2-5aae1c... |
| NDC | 10-11 digits | 0002-3228-30 |
| Drug Name | Text | TYLENOL |
| Application Number | NDA/ANDA + digits | NDA020503 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://dailymed.nlm.nih.gov | Search and browse |
| REST API | https://dailymed.nlm.nih.gov/dailymed/services | Programmatic access |
| FTP Download | https://dailymed.nlm.nih.gov/dailymed/spl-resources | Bulk SPL files |
| RSS Feeds | Available | New label alerts |

## SPL Sections

| Section | Description |
|---------|-------------|
| Indications | FDA-approved uses |
| Dosage | Recommended dosing |
| Contraindications | When not to use |
| Warnings | Safety information |
| Adverse Reactions | Side effects |
| Drug Interactions | Interaction data |

## Limitations

- US market drugs only; international products not included
- SPL parsing can be complex for programmatic access
- Historical labels may not reflect current approved use
- OTC products less consistently formatted

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [RxNorm](../rxnorm/_index.md) - Drug nomenclature
- [DrugBank](../drugbank/_index.md) - Drug database
- [Orange Book](../orange.book/_index.md) - FDA approval data

## References

1. National Library of Medicine. DailyMed. https://dailymed.nlm.nih.gov/
