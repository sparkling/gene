---
id: schema-npatlas
title: "NPAtlas Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, natural-products, microbial, marine]
---

# NPAtlas - Natural Products Atlas Schema

**Document ID:** SCHEMA-NPATLAS
**Version:** 2.0
**Source Version:** 2024

---

## TL;DR

NPAtlas is a curated database of microbial natural products focusing on bacteria and fungi. The schema provides validated chemical structures linked to producing organisms (with NCBI taxonomy) and original literature citations, enabling microbial natural product dereplication and biosynthetic analysis.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Natural Products | 33,000+ | Compound collection |
| Bacterial Compounds | 15,000+ | Actinomycetes, etc. |
| Fungal Compounds | 18,000+ | Ascomycetes, basidiomycetes |
| Marine-Derived | 5,000+ | Marine microorganisms |
| Literature References | 25,000+ | Original publications |

---

## Entity Relationship Overview

```
Compounds (1) ←→ (1) Structures
     ↓                  ↓
  NPAtlas ID       SMILES/InChI

Compounds (1) ←→ (many) Compound_Origins (many) ←→ (1) Organisms
                              ↓
                     NCBI Taxonomy, isolation source

Compounds (1) ←→ (many) Compound_References (many) ←→ (1) References
                              ↓
                         DOI, PubMed ID
```

---

## Core Tables/Entities

### compounds

**Description:** Microbial natural product entries with validated structures.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| npa_id | string | Yes | NPAtlas ID (NPA+digits) |
| name | string | Yes | Compound name |
| canonical_smiles | string | Yes | Canonical SMILES |
| inchi | string | Yes | InChI identifier |
| inchi_key | string | Yes | InChI Key |
| molecular_formula | string | Yes | Molecular formula |
| molecular_weight | decimal | Yes | Molecular weight |
| exact_mass | decimal | No | Exact monoisotopic mass |
| compound_class | string | No | NP class (polyketide, etc.) |
| cluster_type | string | No | BGC type if known |
| origin_type | string | No | bacterial, fungal, marine |

### organisms

**Description:** Source organisms producing natural products.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| organism_id | integer | Yes | Primary identifier |
| scientific_name | string | Yes | Species name |
| ncbi_taxon_id | integer | No | NCBI Taxonomy ID |
| superkingdom | string | No | Bacteria, Eukaryota |
| phylum | string | No | Taxonomic phylum |
| class | string | No | Taxonomic class |
| order | string | No | Taxonomic order |
| family | string | No | Taxonomic family |
| genus | string | No | Genus name |

### compound_origins

**Description:** Links compounds to their source organisms.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| npa_id | string | Yes | Foreign key to compounds |
| organism_id | integer | Yes | Foreign key to organisms |
| isolation_source | string | No | Environment (marine, soil, etc.) |
| geographic_origin | string | No | Collection location |

### references

**Description:** Literature citations for compound isolation/characterization.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| reference_id | integer | Yes | Primary identifier |
| doi | string | No | Digital Object Identifier |
| pmid | integer | No | PubMed ID |
| title | string | No | Article title |
| authors | string | No | Author list |
| journal | string | No | Journal name |
| year | integer | No | Publication year |

### compound_references

**Description:** Links compounds to literature references.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| npa_id | string | Yes | Foreign key to compounds |
| reference_id | integer | Yes | Foreign key to references |
| reference_type | string | No | isolation, synthesis, activity |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /api/v1/compound/{npa_id} | GET | Get compound by ID |
| /api/v1/search | POST | Search compounds |
| /api/v1/organism/{taxon_id}/compounds | GET | Compounds by organism |
| /api/v1/substructure | POST | Substructure search |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | JSON (API), SDF (download) |
| Alternative | TSV, CSV |
| Encoding | UTF-8 |
| Structure | SMILES, InChI, MOL |

---

## Sample Record

```json
{
  "npa_id": "NPA012345",
  "name": "Streptomycin",
  "canonical_smiles": "CC1OC(OC2C(O)C(O)C(NC(N)=N)C(O)C2O)C(OC2OC(CO)C(O)C(N)C2O)C(O)C1NC(N)=N",
  "inchi_key": "UCSJYZPVAKXKNQ-HZYVHMACSA-N",
  "molecular_formula": "C21H39N7O12",
  "molecular_weight": 581.57,
  "compound_class": "Aminoglycoside",
  "origin_type": "bacterial",
  "organisms": [
    {
      "name": "Streptomyces griseus",
      "ncbi_taxon_id": 1911,
      "phylum": "Actinobacteria"
    }
  ],
  "references": [
    {
      "doi": "10.1038/xxx",
      "pmid": 12345678,
      "year": 1944
    }
  ]
}
```

---

## Compound Classes

| Class | Description | Examples |
|-------|-------------|----------|
| Polyketide | PKS-derived | Erythromycin |
| Peptide/NRP | NRPS-derived | Vancomycin |
| Terpene | Terpene synthases | Taxol |
| Alkaloid | Nitrogen-containing | Staurosporine |
| Aminoglycoside | Sugar-amino | Streptomycin |

---

## Glossary

| Term | Definition |
|------|------------|
| NPA ID | NPAtlas identifier (NPA + digits) |
| BGC | Biosynthetic Gene Cluster |
| Dereplication | Identifying known compounds |
| NCBI Taxon ID | NCBI Taxonomy identifier |

---

## References

1. van Santen JA, et al. (2022) Nucleic Acids Res. 50(D1):D1317-D1323
2. NPAtlas: https://www.npatlas.org
3. API Documentation: https://www.npatlas.org/api/v1/docs
