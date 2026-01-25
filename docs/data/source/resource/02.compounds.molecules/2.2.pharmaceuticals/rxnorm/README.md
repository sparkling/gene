---
id: rxnorm
title: "RxNorm - Normalized Drug Nomenclature"
type: source
parent: ../README.md
tier: 1
status: active
category: compounds.molecules
subcategory: pharmaceuticals
tags:
  - terminology
  - drug-names
  - normalization
  - interoperability
  - nlm
---

# RxNorm - Normalized Drug Nomenclature

## Overview

RxNorm is a standardized nomenclature for clinical drugs produced by the National Library of Medicine (NLM). It provides normalized names for clinical drugs and links drug names from various vocabularies to a common identifier, enabling interoperability between pharmacy systems, electronic health records, and drug information sources.

RxNorm serves as a translation layer between the many different drug vocabularies used in pharmacy management, drug interaction software, and electronic health records. By providing a unique concept identifier (RXCUI) for each drug concept, RxNorm enables consistent representation of medications across diverse healthcare systems.

The terminology includes ingredients, drug forms, strengths, and branded/generic products, organized in a hierarchical structure that supports both clinical and administrative use cases.

## Key Statistics

| Metric | Value |
|--------|-------|
| Drug Concepts (RXCUIs) | 300,000+ |
| Source Vocabularies | 12+ |
| Ingredient Concepts | 14,000+ |
| Branded Products | 45,000+ |
| Updates | Monthly |

## Primary Use Cases

1. **Drug Name Normalization** - Standardize drug names across systems
2. **EHR Integration** - Consistent medication representation
3. **Interaction Checking** - Enable drug-drug interaction databases
4. **Analytics** - Aggregate medication data across vocabularies
5. **Clinical Decision Support** - Foundation for drug alerts

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| RXCUI | Integer | 161 (Acetaminophen) |
| RxNorm Name | Standardized | Acetaminophen 500 MG Oral Tablet |
| NDC | 10-11 digits | Cross-referenced |
| ATC Code | 7 characters | Cross-referenced |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| RxNav | https://rxnav.nlm.nih.gov | Interactive browser |
| REST API | https://rxnav.nlm.nih.gov/REST | Programmatic access |
| UMLS Download | https://www.nlm.nih.gov/research/umls/ | Full terminology files |
| RxMix | https://mor.nlm.nih.gov/RxMix/ | Batch processing |

## Concept Types (TTY)

| Type | Description |
|------|-------------|
| IN | Ingredient |
| SCD | Semantic Clinical Drug |
| SBD | Semantic Branded Drug |
| GPCK | Generic Pack |
| BPCK | Branded Pack |
| DF | Dose Form |

## Limitations

- Primarily US drug products; international coverage limited
- UMLS license required for full download
- Veterinary drugs not included
- Non-drug products (devices, supplements) not covered

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [DailyMed](../dailymed/README.md) - Drug labeling
- [DrugBank](../drugbank/README.md) - Drug information
- [ChEMBL](../chembl/README.md) - Bioactivity data

## References

1. Nelson SJ, et al. (2011) "Normalized names for clinical drugs: RxNorm at 6 years." J Am Med Inform Assoc. 18(4):441-448.
