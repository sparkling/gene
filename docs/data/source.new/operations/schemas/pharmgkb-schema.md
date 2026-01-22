---
id: schemas-pharmgkb-schema
title: "PharmGKB Schema Documentation"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, pharmgkb, pharmacogenomics, cpic, drug-response, clinical-annotation]
---

**Parent:** [Schema Documentation](./_index.md)

# PharmGKB Schema Documentation

**Document ID:** SCHEMA-PHARMGKB
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** PharmGKB/ClinPGx API, CPIC API, S3 Downloads

---

## Overview

PharmGKB (Pharmacogenomics Knowledge Base) provides curated pharmacogenomic data including clinical annotations, drug labels, dosing guidelines, and variant-drug associations. The database has rebranded to ClinPGx but maintains the PharmGKB data model.

**API Base:** `https://api.pharmgkb.org/v1/` (redirects to ClinPGx)
**S3 Downloads:** `https://s3.pgkb.org/`
**Rate Limit:** 2 requests/second
**License:** CC BY-SA 4.0

---

## Database Statistics

Based on CPIC API and available documentation:

| Metric | Value | Source |
|--------|-------|--------|
| **CPIC Guidelines** | 33 guidelines | CPIC API |
| **Unique Genes** | ~25 pharmacogenes | CPIC guidelines |
| **Multi-gene Guidelines** | 8 guidelines | 2-6 genes each |
| **S3 File Categories** | 1000+ files | Various formats |

---

## API Endpoints

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

## Core Data Schemas

### Clinical Annotation

Clinical annotations link genetic variants to drug responses with evidence levels.

```json
{
  "id": "1449309937",
  "objCls": "ClinicalAnnotation",
  "name": "Annotation of CPIC Guideline for warfarin and CYP2C9, VKORC1, CYP4F2",
  "location": {
    "id": "1449309937",
    "_url": "/clinicalAnnotation/1449309937"
  },
  "relatedGenes": [
    {
      "id": "PA126",
      "name": "CYP2C9",
      "symbol": "CYP2C9"
    }
  ],
  "relatedChemicals": [
    {
      "id": "PA451906",
      "name": "warfarin"
    }
  ],
  "relatedVariants": [
    {
      "id": "PA166155091",
      "name": "rs1799853"
    }
  ],
  "evidenceLevel": "1A",
  "phenotypeCategories": [
    {
      "id": "PA166182818",
      "name": "Dosage"
    }
  ]
}
```

### Evidence Levels

| Level | Description |
|-------|-------------|
| 1A | CPIC guideline, variant-drug combination with prescribing action |
| 1B | CPIC guideline, high clinical evidence |
| 2A | PharmGKB annotation, moderate clinical evidence |
| 2B | PharmGKB annotation, moderate evidence |
| 3 | Low evidence, annotation only |
| 4 | Case reports only |

---

## Gene Schema

```json
{
  "id": "PA126",
  "objCls": "Gene",
  "name": "cytochrome P450 family 2 subfamily C member 9",
  "symbol": "CYP2C9",
  "chromosome": "10",
  "chromosomeLocation": "10q23.33",
  "ncbiGeneId": "1559",
  "hgncId": "HGNC:2623",
  "ensemblGeneId": "ENSG00000138109",
  "uniprotId": "P11712",
  "hasPharmgkbStar": true,
  "hasCpicGuideline": true
}
```

---

## Variant Schema

```json
{
  "id": "PA166155091",
  "objCls": "Variant",
  "name": "rs1799853",
  "chromosome": "10",
  "position": 94942290,
  "referenceAllele": "C",
  "variantAlleles": ["T"],
  "hgvs": ["NC_000010.11:g.94942290C>T"],
  "clinicalSignificance": "CYP2C9*2 allele",
  "relatedGenes": [
    {
      "id": "PA126",
      "symbol": "CYP2C9"
    }
  ]
}
```

---

## Haplotype/Star Allele Schema

```json
{
  "id": "PA166128353",
  "objCls": "Haplotype",
  "name": "CYP2C9*2",
  "gene": {
    "id": "PA126",
    "symbol": "CYP2C9"
  },
  "activityValue": "Decreased function",
  "variants": [
    {
      "id": "PA166155091",
      "name": "rs1799853"
    }
  ],
  "frequency": {
    "European": 0.13,
    "African": 0.03,
    "Asian": 0.01
  }
}
```

---

## CPIC Guideline Schema

```json
{
  "id": "PA166161537",
  "objCls": "Guideline",
  "name": "CPIC Guideline for warfarin and CYP2C9, CYP4F2, VKORC1",
  "source": "CPIC",
  "version": "3.0",
  "publicationYear": 2023,
  "genes": ["CYP2C9", "CYP4F2", "VKORC1"],
  "drugs": ["warfarin"],
  "recommendations": [
    {
      "phenotype": "CYP2C9 Poor Metabolizer",
      "recommendation": "Decrease initial dose by 40-50%",
      "strength": "Strong"
    }
  ],
  "pmid": "28198005"
}
```

---

## S3 Download Files

| File Category | Format | Description |
|---------------|--------|-------------|
| `clinical_annotations.tsv` | TSV | All clinical annotations |
| `automated_annotations.tsv` | TSV | Literature-mined annotations |
| `drug_labels.tsv` | TSV | FDA drug label annotations |
| `var_drug_ann.tsv` | TSV | Variant-drug annotations |
| `var_pheno_ann.tsv` | TSV | Variant-phenotype annotations |
| `pathways.tsv` | TSV | Pharmacokinetic pathways |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `PA######` | PharmGKB accession for genes, chemicals, or variants | PA126 (CYP2C9) |
| `evidenceLevel` | Clinical annotation evidence strength (1A-4) | 1A = CPIC guideline |
| `haplotype` | Combination of alleles on a chromosome | CYP2C9*2 |
| `diplotype` | Pair of haplotypes inherited from parents | CYP2C9*1/*2 |
| `phenotype` | Observable trait resulting from genotype | Poor Metabolizer |
| `activityValue` | Functional impact of allele | Decreased function |
| `star_allele` | Named genetic variant using star notation | *2, *3 |
| `guideline` | Evidence-based prescribing recommendation | CPIC Guideline |
| `frequency` | Population allele frequency | European: 0.13 |
| `recommendation_strength` | Guideline recommendation strength | Strong/Moderate |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Pharmacogenomics | Study of genetic variation affecting drug response | PharmGKB scope |
| Clinical Annotation | Curated gene-drug-phenotype association with evidence | Core data model |
| Pharmacogene | Gene whose variants affect drug response | CYP2C9, CYP2D6, VKORC1 |
| Star Allele | Nomenclature for named haplotypes in pharmacogenes | *1 = wild-type |
| Metabolizer Status | Phenotype classification based on enzyme activity | Poor/Intermediate/Normal/Ultra-rapid |
| Dosing Guideline | Evidence-based recommendation for drug dosing | CPIC/DPWG |
| Drug Label | FDA-approved prescribing information with PGx content | FDA annotations |
| Clinical Significance | Relevance of variant-drug combination for prescribing | 1A = actionable |

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

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
