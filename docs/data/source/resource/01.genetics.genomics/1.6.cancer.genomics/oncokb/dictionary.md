# OncoKB - Data Dictionary

## Overview

This data dictionary documents the schema for OncoKB precision oncology knowledge base.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | oncokb |
| **Name** | OncoKB |
| **Parent** | 1.6.cancer.genomics |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant Annotation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| hugoSymbol | string | 1:1 | Yes | Gene symbol | BRAF |
| entrezGeneId | integer | 1:1 | Yes | NCBI Gene ID | 673 |
| alteration | string | 1:1 | Yes | Variant name | V600E |
| name | string | 1:1 | Yes | Full alteration name | BRAF V600E |
| consequence | string | 1:1 | Yes | Variant consequence | missense_variant |
| proteinStart | integer | 1:1 | No | Protein position start | 600 |
| proteinEnd | integer | 1:1 | No | Protein position end | 600 |
| refResidues | string | 1:1 | No | Reference amino acid | V |
| variantResidues | string | 1:1 | No | Variant amino acid | E |

### Oncogenicity

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| oncogenic | string | 1:1 | Yes | Oncogenicity status | Oncogenic |
| mutationEffect | string | 1:1 | No | Functional effect | Gain-of-function |
| mutationEffectPmids | array | 1:N | No | Supporting PMIDs | [12345678] |
| mutationEffectDescription | string | 1:1 | No | Effect description | Activates MAPK signaling |

### Treatment Implication

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| level | string | 1:1 | Yes | Evidence level | 1 |
| alterations | array | 1:N | Yes | Target alterations | ["V600E", "V600K"] |
| tumorTypes | array | 1:N | Yes | Cancer types | ["Melanoma"] |
| drugs | array | 1:N | Yes | Associated drugs | ["Vemurafenib", "Dabrafenib"] |
| fdaLevel | string | 1:1 | No | FDA level | Fda2 |
| description | string | 1:1 | No | Treatment description | FDA-approved for... |
| pmids | array | 1:N | No | Reference PMIDs | [12345678, 23456789] |

### Diagnostic/Prognostic

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| level | string | 1:1 | Yes | Diagnostic level | Dx1 |
| tumorTypes | array | 1:N | Yes | Cancer types | ["AML"] |
| description | string | 1:1 | No | Clinical significance | Defines subtype... |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Hugo Symbol | HGNC symbol | BRAF | Gene identifier |
| Entrez Gene ID | Integer | 673 | NCBI Gene ID |
| Alteration | AA change | V600E | Variant notation |
| PubMed ID | Integer | 12345678 | Publication reference |
| OncoTree | Code | MEL | Cancer type code |

---

## Enumerations

### Oncogenicity

| Status | Description |
|--------|-------------|
| Oncogenic | Confirmed driver mutation |
| Likely Oncogenic | Strong evidence for driver |
| Predicted Oncogenic | Predicted driver |
| Likely Neutral | Probably passenger |
| Inconclusive | Insufficient evidence |
| Unknown | Not characterized |

### Mutation Effect

| Effect | Description |
|--------|-------------|
| Gain-of-function | Activating mutation |
| Loss-of-function | Inactivating mutation |
| Likely Gain-of-function | Probably activating |
| Likely Loss-of-function | Probably inactivating |
| Switch-of-function | Altered function |
| Likely Switch-of-function | Probably altered |
| Neutral | No functional impact |
| Inconclusive | Unknown impact |

### Therapeutic Levels

| Level | Description | Evidence |
|-------|-------------|----------|
| 1 | FDA-approved therapy | Same indication |
| 2 | Standard care | Strong evidence |
| 3A | Compelling clinical evidence | Same cancer |
| 3B | Standard care or investigation | Other cancer |
| 4 | Compelling biological evidence | Preclinical |
| R1 | Standard care resistance | Confirmed |
| R2 | Compelling resistance evidence | Investigation |

### FDA Levels

| Level | Description |
|-------|-------------|
| Fda1 | FDA-recognized biomarker |
| Fda2 | FDA-approved companion diagnostic |
| Fda3 | FDA-approved therapy |

### Diagnostic Levels

| Level | Description |
|-------|-------------|
| Dx1 | FDA and/or professional guideline |
| Dx2 | FDA and/or guideline in another tumor |
| Dx3 | Clinical evidence |

### Prognostic Levels

| Level | Description |
|-------|-------------|
| Px1 | FDA and/or professional guideline |
| Px2 | FDA and/or guideline in another tumor |
| Px3 | Clinical evidence |

---

## Entity Relationships

### Gene to Variant
- **Cardinality:** 1:N
- **Description:** One gene has multiple variants
- **Key Fields:** hugoSymbol, alteration

### Variant to Treatment
- **Cardinality:** 1:N
- **Description:** One variant has multiple treatment implications
- **Key Fields:** alteration, tumorType, drugs

### Variant to Oncogenicity
- **Cardinality:** 1:1
- **Description:** One oncogenicity status per variant
- **Key Fields:** hugoSymbol, alteration

### Treatment to Cancer Type
- **Cardinality:** 1:N
- **Description:** Treatments apply to specific cancer types
- **Key Fields:** alteration, tumorTypes

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| OncoKB | Oncology Knowledge Base | Database name |
| FDA | Food and Drug Administration | Regulatory agency |
| NCCN | National Comprehensive Cancer Network | Guidelines |
| CDx | Companion Diagnostic | Biomarker test |
| PMID | PubMed Identifier | Reference ID |
| MSK | Memorial Sloan Kettering | Maintaining institution |
| LOF | Loss of Function | Mutation effect |
| GOF | Gain of Function | Mutation effect |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| OncoTree | Cancer Type | Tumor classification |
| PubMed | PMID | Publication |
| HGNC | Symbol | Gene naming |
| dbSNP | rsID | Variant identifier |
| ClinVar | VCV | Clinical significance |
| CIViC | Variant | Evidence comparison |
| DrugBank | DB ID | Drug information |

---

## OncoTree Cancer Types (Selected)

| Code | Name |
|------|------|
| MEL | Melanoma |
| NSCLC | Non-Small Cell Lung Cancer |
| CRC | Colorectal Cancer |
| BRCA | Breast Cancer |
| AML | Acute Myeloid Leukemia |
| GIST | Gastrointestinal Stromal Tumor |
| THCA | Thyroid Cancer |
| GBM | Glioblastoma |

---

## Data Quality Notes

1. **Cardinality:** One oncogenicity per gene-variant pair
2. **Curation:** Expert-curated by MSK oncologists
3. **Evidence Hierarchy:** Level 1 highest, Level 4 lowest
4. **Updates:** Regular updates with new FDA approvals
5. **Access:** Free for research, license for commercial
6. **Integration:** Available in cBioPortal, clinical reports
