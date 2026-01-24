---
id: schema-cosmic
title: "COSMIC Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: migrated
tags: [schema, database, cancer, somatic, mutation, signature]
---

# COSMIC Schema Documentation

**Document ID:** SCHEMA-COSMIC
**Version:** 1.0
**Source Version:** v99+ (current)

---

## TL;DR

COSMIC catalogs somatic mutations in cancer including coding variants, gene fusions, copy number alterations, and mutational signatures. Data structures include mutation records (COSM/COSV IDs), Cancer Gene Census classifications, and mutational signature profiles (SBS/DBS/ID) derived from tumor sequencing.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Coding Mutations | 17,000,000+ | Expert curated |
| Samples | 1,500,000+ | Tumor samples |
| Tumor Types | 500+ | Primary sites |
| Cancer Gene Census | 736 | Driver genes |
| Mutational Signatures | 100+ | SBS, DBS, ID |
| Gene Fusions | 25,000+ | Fusion pairs |
| Copy Number Variants | 2,000,000+ | CNVs |

---

## Entity Relationship Overview

```
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│    Sample     │────▶│   Mutation    │────▶│    Gene       │
├───────────────┤     ├───────────────┤     ├───────────────┤
│ sample_id     │     │ COSV_id       │     │ gene_symbol   │
│ tumour_type   │     │ COSM_id       │     │ census_tier   │
│ site          │     │ consequence   │     │ role          │
└───────────────┘     └───────────────┘     └───────────────┘
        │                     │
        ▼                     │
┌───────────────┐            │
│   Signature   │◀───────────┘
├───────────────┤
│ SBS/DBS/ID    │
│ contribution  │
│ etiology      │
└───────────────┘
```

---

## Core Tables/Entities

### Mutation

**Description:** Somatic mutation record with genomic details

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| COSV_ID | string | Yes | Genomic mutation ID (COSV + integer) |
| COSM_ID | string | No | Legacy mutation ID (COSM + integer) |
| Gene | string | Yes | HGNC gene symbol |
| Transcript | string | Yes | Ensembl transcript ID |
| CDS_Mutation | string | No | CDS-level HGVS |
| AA_Mutation | string | No | Protein-level HGVS |
| Mutation_Type | string | Yes | Substitution, Insertion, etc. |
| Mutation_Description | string | Yes | Functional consequence |
| GRCh38_Position | string | Yes | chr:start-end |
| GRCh37_Position | string | No | Legacy coordinates |
| Strand | string | Yes | + or - |

### Sample

**Description:** Tumor sample with clinical metadata

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Sample_ID | integer | Yes | Internal sample ID |
| Sample_Name | string | Yes | Sample identifier |
| Tumour_ID | integer | Yes | Tumor record ID |
| Primary_Site | string | Yes | Anatomical site |
| Site_Subtype1 | string | No | Subtype level 1 |
| Primary_Histology | string | Yes | Histological type |
| Histology_Subtype1 | string | No | Histology subtype |
| Sample_Type | string | Yes | cell line, tumour, etc. |
| Age | integer | No | Patient age |
| Pubmed_PMID | integer | No | Publication reference |

### Cancer Gene Census

**Description:** Curated cancer driver gene list

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Gene_Symbol | string | Yes | HGNC symbol |
| Entrez_GeneId | integer | Yes | NCBI Gene ID |
| Tier | integer | Yes | 1 or 2 |
| Hallmark | string | No | Cancer hallmarks |
| Role_in_Cancer | string | Yes | oncogene, TSG, fusion |
| Mutation_Types | string | Yes | Types of mutations |
| Translocation_Partner | string | No | Fusion partners |
| Tumour_Types_Somatic | string | No | Associated cancers |
| Tumour_Types_Germline | string | No | Hereditary syndromes |
| Molecular_Genetics | string | No | Dom/Rec inheritance |

### Mutational Signature

**Description:** Mutational signature profiles

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Signature | string | Yes | SBS/DBS/ID + number |
| Type | string | Yes | SBS, DBS, ID |
| Contribution | float | Yes | Proportion (0-1) |
| Proposed_Etiology | string | No | Mutagenic process |
| Associated_Cancers | string | No | Cancer types |
| Features | array | No | Trinucleotide context |

### Signature Types

| Type | Description | Example |
|------|-------------|---------|
| SBS | Single Base Substitution | SBS1, SBS2, SBS13 |
| DBS | Doublet Base Substitution | DBS1, DBS2 |
| ID | Small Insertion/Deletion | ID1, ID2, ID8 |

### Gene Fusion

**Description:** Gene fusion events

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| Fusion_ID | integer | Yes | Fusion record ID |
| Gene_5prime | string | Yes | 5' partner gene |
| Gene_3prime | string | Yes | 3' partner gene |
| Breakpoint_5prime | string | No | 5' breakpoint |
| Breakpoint_3prime | string | No | 3' breakpoint |
| Fusion_Type | string | Yes | Inframe, etc. |
| Sample_ID | integer | Yes | Source sample |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | VCF, TSV |
| Alternative | MySQL dump |
| Compression | gzip |
| Encoding | UTF-8 |
| Reference | GRCh38, GRCh37 |

---

## Sample Record

```json
{
  "mutation": {
    "COSV_ID": "COSV54736898",
    "COSM_ID": "COSM476",
    "Gene": "BRAF",
    "Transcript": "ENST00000288602.11",
    "CDS_Mutation": "c.1799T>A",
    "AA_Mutation": "p.V600E",
    "Mutation_Type": "Substitution - Missense",
    "GRCh38_Position": "7:140753336-140753336"
  },
  "sample": {
    "Sample_ID": 123456,
    "Primary_Site": "skin",
    "Primary_Histology": "malignant_melanoma",
    "Sample_Type": "tumour"
  },
  "census_gene": {
    "Gene_Symbol": "BRAF",
    "Tier": 1,
    "Role_in_Cancer": "oncogene",
    "Hallmark": "proliferative signalling",
    "Mutation_Types": "Mis, F"
  },
  "signature": {
    "Signature": "SBS7a",
    "Type": "SBS",
    "Proposed_Etiology": "Ultraviolet light exposure",
    "Contribution": 0.85
  }
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| COSV | Genomic mutation identifier (GRCh38-based) |
| COSM | Legacy mutation identifier |
| SBS | Single Base Substitution signature |
| DBS | Doublet Base Substitution signature |
| ID | Insertion/Deletion signature |
| CGC | Cancer Gene Census |
| TSG | Tumor Suppressor Gene |
| Tier 1 | Strong evidence for cancer causation |
| Tier 2 | Evidence consistent with cancer role |

---

## References

1. https://cancer.sanger.ac.uk/cosmic
2. Tate et al. (2019) Nucleic Acids Res. DOI: 10.1093/nar/gky1015
3. https://cancer.sanger.ac.uk/signatures/
