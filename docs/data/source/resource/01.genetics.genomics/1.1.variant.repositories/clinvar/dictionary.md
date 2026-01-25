# ClinVar - Data Dictionary

## Overview

This data dictionary documents the schema for ClinVar clinical variant interpretation database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | clinvar |
| **Name** | ClinVar |
| **Parent** | 1.1.variant.repositories |
| **Total Fields** | 12 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### VariationArchive (VCV)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| VariationID | integer | 1:1 | Yes | Internal variant identifier | 17661, 12345 |
| Name | string | 1:1 | Yes | HGVS expression | NM_000546.6(TP53):c.743G>A (p.Arg248Gln) |
| ClassifiedRecord.Classification | string | 1:1 | Yes | Aggregated clinical interpretation | Pathogenic, Likely pathogenic |
| ReviewStatus | string | 1:1 | Yes | Evidence level (0-4 stars) | criteria provided, multiple submitters |
| ConditionList | array | 1:N | No | Associated diseases | Li-Fraumeni syndrome |

### Submission (SCV)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| SubmitterName | string | 1:1 | Yes | Contributing organization | GeneDx, Invitae |
| DateSubmitted | date | 1:1 | Yes | Submission date | 2024-01-15 |
| ClinicalSignificance | string | 1:1 | Yes | Submitter interpretation | Pathogenic |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| VCV | VCV + 9 digits + version | VCV000012345.3 | Variation accession |
| RCV | RCV + 9 digits + version | RCV000012345.6 | Variant-condition accession |
| Variation ID | Integer | 12345 | Internal database ID |
| rsID | rs + digits | rs28934576 | dbSNP cross-reference |

---

## Enumerations

### Classification

| Value | Meaning |
|-------|---------|
| Pathogenic | Disease-causing variant |
| Likely pathogenic | >90% certainty pathogenic |
| Uncertain significance | Insufficient evidence (VUS) |
| Likely benign | >90% certainty benign |
| Benign | Not disease-causing |
| Conflicting interpretations | Submitter disagreement |

### Review Status (Stars)

| Value | Stars | Meaning |
|-------|-------|---------|
| practice guideline | 4 | Expert panel guideline |
| reviewed by expert panel | 3 | Expert review |
| criteria provided, multiple submitters, no conflicts | 2 | Multi-submitter consensus |
| criteria provided, single submitter | 1 | Single lab with criteria |
| no assertion criteria provided | 0 | No criteria given |

---

## Entity Relationships

### VCV to RCV
- **Cardinality:** 1:N
- **Description:** One VCV (aggregate) contains multiple RCVs (variant+condition pairs)
- **Key Fields:** VariationID, RCV

### RCV to SCV
- **Cardinality:** 1:N
- **Description:** One RCV contains multiple SCVs (individual submissions)
- **Key Fields:** RCV, SubmitterID

### Variant to Gene
- **Cardinality:** N:1
- **Description:** Many variants map to one gene
- **Key Fields:** VariationID, GeneID

### Variant to Condition
- **Cardinality:** N:M
- **Description:** Variants associate with multiple conditions
- **Key Fields:** VariationID, MedGenCUI

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| VCV | Variation Archive | Aggregate record |
| RCV | Reference Condition Variation | Variant-condition pair |
| SCV | Submission Condition Variation | Individual submission |
| VUS | Variant of Uncertain Significance | Classification category |
| ACMG | American College of Medical Genetics | Classification guidelines |
| AMP | Association for Molecular Pathology | Classification guidelines |
| HGVS | Human Genome Variation Society | Nomenclature standard |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant identifier |
| OMIM | MIM number | Disease reference |
| MedGen | CUI | Condition ontology |
| Gene | GeneID | Gene reference |

---

## Data Quality Notes

1. **Cardinality:** VCV→RCV→SCV forms hierarchical structure
2. **Conflicts:** Same variant may have conflicting classifications from different submitters
3. **Review Status:** Higher stars indicate stronger evidence consensus
4. **Versioning:** VCV and RCV IDs include version numbers for tracking updates
