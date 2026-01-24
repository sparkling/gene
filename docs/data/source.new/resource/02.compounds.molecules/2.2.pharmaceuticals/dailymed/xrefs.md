---
id: xrefs-dailymed
title: "DailyMed Cross-References"
type: xrefs
parent: _index.md
last_updated: 2026-01-23
---

# DailyMed Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `02.compounds.molecules/2.2.pharmaceuticals/dailymed` |
| Secondary | Symlink | `08.literature.knowledge/8.4.regulatory.legal/dailymed` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| NDC | ndc | High |
| RxNorm | rxcui | High |
| UNII | unii | High |
| SPL | setid | High |
| NDA/ANDA | application_number | High |
| FDA Orange Book | orange_book_id | Medium |
| MedlinePlus | medlineplus | Medium |

## Integration Notes

DailyMed provides official FDA-approved drug labeling (SPL documents), bridging pharmaceutical compound data with regulatory documentation.

**Primary use (Pharmaceuticals):** Drug formulations, active ingredients, dosage forms, strengths, manufacturers.

**Secondary use (Regulatory/Legal):** FDA-approved labeling, indications, contraindications, warnings, boxed warnings, drug interactions, prescribing information.

DailyMed is the authoritative source for current US drug labeling, essential for regulatory compliance, drug information systems, and clinical decision support.
