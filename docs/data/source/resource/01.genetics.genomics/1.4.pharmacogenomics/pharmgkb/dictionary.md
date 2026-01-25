# PharmGKB - Data Dictionary

## Overview

This data dictionary documents the schema for PharmGKB pharmacogenomics knowledge base.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pharmgkb |
| **Name** | PharmGKB |
| **Parent** | 1.4.pharmacogenomics |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Gene

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| PharmGKB_Accession | string | 1:1 | Yes | PA + integer | PA124 |
| Symbol | string | 1:1 | Yes | HGNC symbol | CYP2D6 |
| Name | string | 1:1 | Yes | Full gene name | cytochrome P450 family 2 subfamily D member 6 |
| Alternative_Names | array | 1:N | No | Aliases | debrisoquine 4-hydroxylase |
| Alternative_Symbols | array | 1:N | No | Symbol aliases | CYP2D |
| Is_VIP | boolean | 1:1 | No | Very Important Pharmacogene | true |
| Has_Variant_Annotation | boolean | 1:1 | No | Has curated variants | true |
| Cross_References | object | 1:1 | No | External IDs | {"HGNC": "2625"} |

### Drug

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| PharmGKB_Accession | string | 1:1 | Yes | PA + integer | PA449053 |
| Name | string | 1:1 | Yes | Generic name | codeine |
| Trade_Names | array | 1:N | No | Brand names | ["Tylenol with Codeine"] |
| Type | string | 1:1 | Yes | Drug type | Small molecule |
| ATC_Codes | array | 1:N | No | ATC classification | ["R05DA04"] |
| RxNorm_Identifiers | array | 1:N | No | RxNorm CUIs | ["2670"] |
| DrugBank_ID | string | 1:1 | No | DrugBank identifier | DB00318 |
| Has_PGx_On_Label | boolean | 1:1 | No | FDA label PGx | true |

### Clinical Annotation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Annotation_ID | integer | 1:1 | Yes | Unique identifier | 1449309937 |
| Location | string | 1:1 | Yes | Variant or haplotype | CYP2D6, rs3892097 |
| Gene | string | 1:1 | Yes | Gene symbol | CYP2D6 |
| Level_of_Evidence | string | 1:1 | Yes | Evidence level | 1A, 1B, 2A, 2B, 3, 4 |
| Phenotype_Category | string | 1:1 | Yes | Clinical category | Efficacy, Toxicity, Dosage |
| Drug | string | 1:1 | Yes | Drug name | codeine |
| Phenotype | string | 1:1 | Yes | Clinical phenotype | Reduced morphine formation |
| URL | string | 1:1 | Yes | PharmGKB URL | https://www.pharmgkb.org/... |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| PA Accession | PA + digits | PA124 | PharmGKB internal ID |
| Gene Symbol | HGNC symbol | CYP2D6 | Gene identifier |
| rsID | rs + digits | rs3892097 | dbSNP identifier |
| Star Allele | Gene*number | CYP2D6*4 | Haplotype designation |
| ATC Code | Letter + digits | R05DA04 | Drug classification |
| DrugBank | DB + 5 digits | DB00318 | Drug database ID |

---

## Enumerations

### Evidence Levels

| Level | Description | Criteria |
|-------|-------------|----------|
| 1A | Guideline annotation | CPIC/DPWG guideline |
| 1B | Label annotation | FDA/EMA label |
| 2A | Multiple studies | Replicated association |
| 2B | Moderate evidence | Single study, strong |
| 3 | Low evidence | Limited evidence |
| 4 | Case report | Case reports only |

### Phenotype Categories

| Category | Description |
|----------|-------------|
| Efficacy | Drug effectiveness |
| Toxicity | Adverse reactions |
| Dosage | Dose modifications |
| Metabolism/PK | Pharmacokinetics |
| Other | Other phenotypes |

### Pathway Types

| Type | Description |
|------|-------------|
| PK | Pharmacokinetics - drug metabolism |
| PD | Pharmacodynamics - drug action |

### Drug Label Sources

| Source | Description |
|--------|-------------|
| FDA | US Food and Drug Administration |
| EMA | European Medicines Agency |
| PMDA | Japan Pharmaceuticals Agency |
| HCSC | Health Canada |

### PGx Testing Levels

| Level | Meaning |
|-------|---------|
| Testing required | Mandatory before prescribing |
| Testing recommended | Strongly suggested |
| Actionable PGx | Information in label |
| Informative PGx | Background information |

---

## Entity Relationships

### Gene to Drug
- **Cardinality:** N:M
- **Description:** Many genes affect many drugs
- **Key Fields:** gene_accession, drug_accession

### Gene to Variant
- **Cardinality:** 1:N
- **Description:** Multiple variants per gene
- **Key Fields:** gene_accession, variant_id

### Variant to Haplotype
- **Cardinality:** N:M
- **Description:** Variants define haplotypes
- **Key Fields:** variant_id, haplotype_name

### Drug to Clinical Annotation
- **Cardinality:** 1:N
- **Description:** Multiple annotations per drug
- **Key Fields:** drug_accession, annotation_id

### Guideline to Gene-Drug
- **Cardinality:** N:N
- **Description:** Guidelines cover multiple pairs
- **Key Fields:** guideline_id, gene, drug

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PharmGKB | Pharmacogenomics Knowledge Base | Database name |
| ClinPGx | Clinical Pharmacogenomics | Rebrand name |
| VIP | Very Important Pharmacogene | Key genes |
| CPIC | Clinical Pharmacogenetics Implementation Consortium | Guideline source |
| DPWG | Dutch Pharmacogenetics Working Group | Guideline source |
| PGx | Pharmacogenomics | Field abbreviation |
| PK | Pharmacokinetics | Drug metabolism |
| PD | Pharmacodynamics | Drug action |
| CYP | Cytochrome P450 | Enzyme family |
| HGVS | Human Genome Variation Society | Nomenclature |
| rsID | Reference SNP identifier | dbSNP format |
| ATC | Anatomical Therapeutic Chemical | Drug classification |
| FDA | Food and Drug Administration | Regulatory agency |
| EMA | European Medicines Agency | Regulatory agency |
| CC BY-SA | Creative Commons Attribution-ShareAlike | License type |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| CPIC | Guideline | Clinical guidelines |
| DPWG | Guideline | Clinical guidelines |
| PharmVar | Allele | Star allele definitions |
| DrugBank | DB ID | Drug information |
| RxNorm | CUI | Drug ontology |
| HGNC | Symbol | Gene naming |
| dbSNP | rsID | Variant identifier |
| ClinVar | VCV | Clinical significance |

---

## Data Quality Notes

1. **Cardinality:** One annotation per gene-drug-phenotype-evidence combination
2. **Curation:** Expert-curated from literature
3. **Updates:** Continuously updated with new publications
4. **Coverage:** 800+ drugs, 150+ VIP genes
5. **Evidence Hierarchy:** Level 1A/1B highest confidence
6. **Guideline Integration:** Links to CPIC/DPWG recommendations
