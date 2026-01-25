# NPAtlas - Data Dictionary

## Overview

This data dictionary documents the schema for NPAtlas (Natural Products Atlas) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | npatlas |
| **Name** | NPAtlas |
| **Parent** | 2.1.natural.products |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| npa_id | string | 1:1 | Yes | NPAtlas identifier | NPA012345 |
| name | string | 1:1 | Yes | Compound name | Streptomycin |
| canonical_smiles | string | 1:1 | Yes | Canonical SMILES | CC1OC(OC2C(O)... |
| inchi | string | 1:1 | Yes | InChI identifier | InChI=1S/C21H39N7O12/... |
| inchi_key | string | 1:1 | Yes | InChI Key | UCSJYZPVAKXKNQ-HZYVHMACSA-N |
| molecular_formula | string | 1:1 | Yes | Molecular formula | C21H39N7O12 |
| molecular_weight | decimal | 1:1 | Yes | Molecular weight | 581.57 |
| exact_mass | decimal | 1:1 | No | Monoisotopic mass | 581.2657 |
| compound_class | string | 1:1 | No | NP class | Aminoglycoside |
| cluster_type | string | 1:1 | No | BGC type if known | NRPS, PKS |
| origin_type | string | 1:1 | No | Source category | bacterial, fungal, marine |

### Organism Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| organism_id | integer | 1:1 | Yes | Primary identifier | 101 |
| scientific_name | string | 1:1 | Yes | Species name | Streptomyces griseus |
| ncbi_taxon_id | integer | 1:1 | No | NCBI Taxonomy ID | 1911 |
| superkingdom | string | 1:1 | No | Domain | Bacteria, Eukaryota |
| phylum | string | 1:1 | No | Taxonomic phylum | Actinobacteria |
| class | string | 1:1 | No | Taxonomic class | Actinobacteria |
| order | string | 1:1 | No | Taxonomic order | Streptomycetales |
| family | string | 1:1 | No | Taxonomic family | Streptomycetaceae |
| genus | string | 1:1 | No | Genus name | Streptomyces |

### Compound Origin Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| npa_id | string | 1:1 | Yes | Foreign key to compounds | NPA012345 |
| organism_id | integer | 1:1 | Yes | Foreign key to organisms | 101 |
| isolation_source | string | 1:1 | No | Environment | marine, soil, endophyte |
| geographic_origin | string | 1:1 | No | Collection location | Pacific Ocean |

### Reference Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| reference_id | integer | 1:1 | Yes | Primary identifier | 201 |
| doi | string | 1:1 | No | Digital Object Identifier | 10.1038/xxx |
| pmid | integer | 1:1 | No | PubMed ID | 12345678 |
| title | string | 1:1 | No | Article title | Discovery of... |
| authors | string | 1:1 | No | Author list | Smith J, et al. |
| journal | string | 1:1 | No | Journal name | J Nat Prod |
| year | integer | 1:1 | No | Publication year | 2020 |
| reference_type | string | 1:1 | No | Citation type | isolation, synthesis, activity |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| NPA ID | NPA + digits | NPA012345 | Primary identifier |
| InChI Key | 27 characters | UCSJYZPVAKXKNQ-HZYVHMACSA-N | Structure hash |
| NCBI Taxon ID | Integer | 1911 | Organism identifier |
| DOI | 10.xxxx/... | 10.1038/s41586 | Publication identifier |
| PMID | Integer | 12345678 | PubMed article ID |

---

## Enumerations

### Compound Classes

| Class | Description | Examples |
|-------|-------------|----------|
| Polyketide | PKS-derived | Erythromycin |
| Peptide/NRP | NRPS-derived | Vancomycin |
| Terpene | Terpene synthases | Taxol |
| Alkaloid | Nitrogen-containing | Staurosporine |
| Aminoglycoside | Sugar-amino | Streptomycin |
| Hybrid | Multiple pathways | Epothilone |

### BGC (Biosynthetic Gene Cluster) Types

| Type | Description |
|------|-------------|
| NRPS | Nonribosomal peptide synthetase |
| PKS | Polyketide synthase |
| Terpene | Terpene cyclase |
| RiPP | Ribosomally synthesized peptide |
| Hybrid | Multiple BGC types |

### Origin Types

| Type | Description |
|------|-------------|
| bacterial | Bacterial metabolites |
| fungal | Fungal metabolites |
| marine | Marine-derived microbes |
| terrestrial | Land-derived microbes |

### Isolation Sources

| Source | Description |
|--------|-------------|
| marine | Ocean samples |
| soil | Terrestrial soil |
| endophyte | Plant-associated |
| symbiont | Animal-associated |
| freshwater | Lake/river samples |

### Reference Types

| Type | Description |
|------|-------------|
| isolation | Original isolation paper |
| synthesis | Total synthesis report |
| activity | Bioactivity study |
| structure | Structure elucidation |

---

## Entity Relationships

### Compound to Organisms
- **Cardinality:** N:M
- **Description:** Compounds from multiple organisms; organisms produce multiple compounds
- **Key Fields:** npa_id, organism_id

### Compound to References
- **Cardinality:** N:M
- **Description:** Compounds cited in multiple papers; papers describe multiple compounds
- **Key Fields:** npa_id, reference_id

### Organism to Taxonomy
- **Cardinality:** 1:1
- **Description:** Each organism has one taxonomic classification
- **Key Fields:** organism_id, ncbi_taxon_id

### Compound-Organism-Reference Triplet
- **Cardinality:** Core model
- **Description:** Documented occurrence requires organism and literature reference
- **Key Fields:** npa_id, organism_id, reference_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NPAtlas | Natural Products Atlas | Database name |
| NPA | NPAtlas | Identifier prefix |
| BGC | Biosynthetic Gene Cluster | Gene region |
| NRPS | Nonribosomal Peptide Synthetase | Biosynthesis |
| PKS | Polyketide Synthase | Biosynthesis |
| RiPP | Ribosomally Synthesized and Post-translationally modified Peptide | Biosynthesis |
| SMILES | Simplified Molecular Input Line Entry System | Structure |
| InChI | International Chemical Identifier | IUPAC standard |
| NCBI | National Center for Biotechnology Information | Taxonomy |
| DOI | Digital Object Identifier | Publication ID |
| PMID | PubMed Identifier | Literature reference |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PubChem | InChI Key | Chemical data |
| ChEMBL | InChI Key | Bioactivity |
| NCBI Taxonomy | Taxon ID | Organism classification |
| PubMed | PMID | Literature |
| CrossRef | DOI | Publication metadata |
| antiSMASH | BGC | Biosynthetic clusters |
| MIBiG | BGC ID | Validated BGCs |

---

## Data Quality Notes

1. **Microbial Focus:** Exclusively bacterial and fungal natural products
2. **Structure Validation:** All structures manually curated and validated
3. **Literature Linked:** Every compound linked to original publication
4. **Taxonomy Verified:** NCBI Taxonomy cross-references for organisms
5. **Marine Emphasis:** Strong coverage of marine microbial NPs
6. **Dereplication Tool:** Designed for known compound identification
7. **Update Frequency:** Regular updates with newly published compounds
