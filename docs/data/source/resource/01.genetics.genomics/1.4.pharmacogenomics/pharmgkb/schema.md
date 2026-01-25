---
id: schema-pharmgkb
title: "PharmGKB Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: migrated
tags: [schema, database, pharmacogenomics, drug-response, clinical, pgx]
---

# PharmGKB Schema Documentation

**Document ID:** SCHEMA-PHARMGKB
**Version:** 1.0
**Source Version:** Current (continuously updated)

---

## TL;DR

PharmGKB provides curated pharmacogenomics knowledge linking genetic variants to drug response. Data includes clinical annotations with evidence levels, variant-drug associations, gene-drug relationships, pathway diagrams, and curated drug label information from regulatory agencies (FDA, EMA, PMDA).

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Clinical Annotations | 20,000+ | Evidence-based |
| Variant Annotations | 170,000+ | Literature curated |
| Drugs | 800+ | With PGx data |
| Genes | 150+ | VIP genes |
| Pathways | 200+ | Drug pathways |
| Drug Labels | 400+ | FDA/EMA/PMDA |
| CPIC Guidelines | 33 | Linked |

---

## Entity Relationship Overview

```
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│    Gene       │────▶│  Relationship │────▶│    Drug       │
├───────────────┤     ├───────────────┤     ├───────────────┤
│ PA accession  │     │ evidence_level│     │ PA accession  │
│ symbol        │     │ PK/PD         │     │ generic_name  │
│ VIP status    │     │ variant_annot │     │ trade_names   │
└───────────────┘     └───────────────┘     └───────────────┘
        │                     │
        ▼                     ▼
┌───────────────┐     ┌───────────────┐
│   Variant     │     │Clinical Annot │
├───────────────┤     ├───────────────┤
│ rs_id         │     │ annotation_id │
│ haplotype     │     │ level         │
│ function      │     │ phenotypes    │
└───────────────┘     └───────────────┘
```

---

## Core Tables/Entities

### Gene

**Description:** Pharmacogenomics gene information

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| PharmGKB_Accession | string | Yes | PA + integer (PA124) |
| Symbol | string | Yes | HGNC symbol |
| Name | string | Yes | Full gene name |
| Alternative_Names | array | No | Aliases |
| Alternative_Symbols | array | No | Symbol aliases |
| Is_VIP | boolean | No | Very Important Pharmacogene |
| Has_Variant_Annotation | boolean | No | Has curated variants |
| Cross_References | object | No | External IDs |

### Drug

**Description:** Drug with pharmacogenomic associations

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| PharmGKB_Accession | string | Yes | PA + integer |
| Name | string | Yes | Generic name |
| Trade_Names | array | No | Brand names |
| Type | string | Yes | Drug type |
| ATC_Codes | array | No | ATC classification |
| RxNorm_Identifiers | array | No | RxNorm CUIs |
| DrugBank_ID | string | No | DrugBank identifier |
| Has_PGx_On_Label | boolean | No | FDA label PGx |

### Clinical Annotation

**Description:** Evidence-based clinical annotation

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Annotation_ID | integer | Yes | Unique identifier |
| Location | string | Yes | Variant or haplotype |
| Gene | string | Yes | Gene symbol |
| Level_of_Evidence | string | Yes | 1A, 1B, 2A, 2B, 3, 4 |
| Phenotype_Category | string | Yes | Efficacy, Toxicity, etc. |
| Drug | string | Yes | Drug name |
| Phenotype | string | Yes | Clinical phenotype |
| URL | string | Yes | PharmGKB URL |

### Evidence Levels

| Level | Description | Criteria |
|-------|-------------|----------|
| 1A | Guideline annotation | CPIC/DPWG guideline |
| 1B | Label annotation | FDA/EMA label |
| 2A | Multiple studies | Replicated association |
| 2B | Moderate evidence | Single study, strong |
| 3 | Low evidence | Limited evidence |
| 4 | Case report | Case reports only |

### Variant Annotation

**Description:** Literature-curated variant information

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Variant_ID | string | Yes | Variant identifier |
| Gene | string | Yes | Gene symbol |
| RS_ID | string | No | dbSNP RS ID |
| Location | string | No | Genomic location |
| Variant_Haplotypes | array | No | Associated haplotypes |
| Variant_Summary | string | No | Functional summary |
| Biogeographical_Groups | object | No | Population frequencies |

### Drug Label

**Description:** Regulatory agency drug label PGx

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Label_ID | string | Yes | Label identifier |
| Drug | string | Yes | Drug name |
| Source | string | Yes | FDA, EMA, PMDA, HCSC |
| Genes | array | Yes | PGx genes mentioned |
| PGx_Level | string | Yes | Testing required/recommended |
| Summary | string | Yes | Label content summary |

### Pathway

**Description:** Drug mechanism pathway

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Pathway_ID | string | Yes | PA + integer |
| Name | string | Yes | Pathway name |
| Drug | string | Yes | Associated drug |
| Genes | array | Yes | Pathway genes |
| Type | string | Yes | PK (pharmacokinetics), PD (pharmacodynamics) |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | TSV (tab-separated) |
| API | JSON |
| Alternative | Excel (some files) |
| Encoding | UTF-8 |

---

## Sample Record

```json
{
  "gene": {
    "PharmGKB_Accession": "PA124",
    "Symbol": "CYP2D6",
    "Name": "cytochrome P450 family 2 subfamily D member 6",
    "Is_VIP": true
  },
  "drug": {
    "PharmGKB_Accession": "PA449053",
    "Name": "codeine",
    "Trade_Names": ["Tylenol with Codeine"],
    "Has_PGx_On_Label": true
  },
  "clinical_annotation": {
    "Annotation_ID": 1449309937,
    "Location": "CYP2D6",
    "Gene": "CYP2D6",
    "Level_of_Evidence": "1A",
    "Phenotype_Category": "Efficacy",
    "Drug": "codeine",
    "Phenotype": "CYP2D6 poor metabolizers have reduced conversion of codeine to morphine"
  },
  "variant": {
    "RS_ID": "rs3892097",
    "Gene": "CYP2D6",
    "Location": "chr22:42128174",
    "Haplotype": "*4"
  },
  "drug_label": {
    "Drug": "codeine",
    "Source": "FDA",
    "PGx_Level": "Genetic testing recommended",
    "Genes": ["CYP2D6"]
  }
}
```

---

## API Endpoints

**API Base:** `https://api.pharmgkb.org/v1/` (redirects to ClinPGx)
**S3 Downloads:** `https://s3.pgkb.org/`
**Rate Limit:** 2 requests/second

### Core Data Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/data/clinicalAnnotation` | GET | Clinical annotations |
| `/data/gene` | GET | Gene information |
| `/data/chemical` | GET | Drug/chemical data |
| `/data/variant` | GET | Genetic variants |
| `/data/label` | GET | Drug label annotations |
| `/data/guideline` | GET | Dosing guidelines |
| `/data/pathway` | GET | Pharmacokinetic pathways |
| `/data/automaticAnnotation` | GET | Automated literature annotations |

### Reference Data

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/data/haplotype` | GET | Star allele definitions |
| `/data/diplotype` | GET | Diplotype-phenotype mappings |
| `/data/phenotype` | GET | Phenotype definitions |

---

## S3 Download Files

| Resource | Format | Size | URL |
|----------|--------|------|-----|
| Clinical Annotations | TSV | ~20 MB | https://s3.pgkb.org/data/clinical_annotations.tsv.zip |
| Drug Labels | TSV | ~10 MB | https://s3.pgkb.org/data/drug_labels.tsv.zip |
| Variant-Drug Annotations | TSV | ~15 MB | https://s3.pgkb.org/data/var_drug_ann.tsv.zip |
| Variant-Phenotype Annotations | TSV | ~8 MB | https://s3.pgkb.org/data/var_pheno_ann.tsv.zip |
| Pathways Data | TSV | ~5 MB | https://s3.pgkb.org/data/pathways.tsv.zip |
| Guidelines | JSON/PDF | ~50 MB | https://cpicpgx.org/guidelines/ |
| REST API | JSON | - | https://api.pharmgkb.org/v1/ |

---

## Glossary

| Term | Definition |
|------|------------|
| PA Accession | PharmGKB accession identifier |
| VIP | Very Important Pharmacogene |
| CPIC | Clinical Pharmacogenetics Implementation Consortium |
| DPWG | Dutch Pharmacogenetics Working Group |
| PK | Pharmacokinetics - drug metabolism |
| PD | Pharmacodynamics - drug action |
| Level 1A | Highest evidence - guideline annotation |
| Star Allele | Haplotype designation (e.g., *1, *2) |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PharmGKB | Pharmacogenomics Knowledge Base | Database name |
| ClinPGx | Clinical Pharmacogenomics | Rebrand name |
| CPIC | Clinical Pharmacogenetics Implementation Consortium | Guideline source |
| DPWG | Dutch Pharmacogenetics Working Group | Guideline source |
| PGx | Pharmacogenomics | Field abbreviation |
| CYP | Cytochrome P450 | Enzyme family |
| VKORC1 | Vitamin K Epoxide Reductase Complex 1 | Warfarin target |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| rsID | Reference SNP identifier | dbSNP format |
| HGVS | Human Genome Variation Society | Nomenclature standard |
| CC BY-SA | Creative Commons Attribution-ShareAlike | License type |
| FDA | Food and Drug Administration | Drug label source |
| EMA | European Medicines Agency | Drug label source |

---

## References

1. https://www.pharmgkb.org/
2. Whirl-Carrillo et al. (2021) Clin Pharmacol Ther. DOI: 10.1002/cpt.2122
3. https://www.pharmgkb.org/downloads
