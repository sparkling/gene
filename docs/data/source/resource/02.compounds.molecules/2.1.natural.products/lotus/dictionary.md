# LOTUS - Data Dictionary

## Overview

This data dictionary documents the schema for LOTUS (Linked Open Total Unified Structure-organism) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | lotus |
| **Name** | LOTUS |
| **Parent** | 2.1.natural.products |
| **Total Fields** | 35+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| wikidataId | string | 1:1 | Yes | Wikidata Q-identifier | Q312879 |
| canonicalSmiles | string | 1:1 | Yes | Canonical SMILES notation | CC(C)CCCC(C)C |
| isomericSmiles | string | 1:1 | No | SMILES with stereochemistry | CC[C@H](CC[C@@H](C)... |
| inchi | string | 1:1 | Yes | Standard InChI string | InChI=1S/C29H50O/... |
| inchiKey | string | 1:1 | Yes | 27-character InChI hash | KZJWDPNRJALLNS-VJSFXXLFSA-N |
| molecular_formula | string | 1:1 | No | Molecular formula | C29H50O |
| molecular_weight | float | 1:1 | No | Molecular weight | 414.71 |
| iupac_name | string | 1:1 | No | IUPAC systematic name | beta-sitosterol |
| traditional_name | string | 1:1 | No | Common/traditional name | Sitosterol |

### Organism/Taxon Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| taxonWikidataId | string | 1:1 | Yes | Wikidata Q-ID for organism | Q158695 |
| name | string | 1:1 | Yes | Scientific name | Arabidopsis thaliana |
| rank | string | 1:1 | No | Taxonomic rank | species |
| ncbiId | integer | 1:1 | No | NCBI Taxonomy ID | 3702 |
| ottId | integer | 1:1 | No | Open Tree of Life ID | 1054861 |

### Reference Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| referenceWikidataId | string | 1:1 | Yes | Wikidata Q-ID for reference | Q56636231 |
| doi | string | 1:1 | No | Digital Object Identifier | 10.1016/j.phytochem.2018.01.001 |
| pmid | integer | 1:1 | No | PubMed ID | 29371045 |
| pmcid | string | 1:1 | No | PubMed Central ID | PMC1234567 |
| title | string | 1:1 | No | Publication title | Article title |

### LNPN Extended Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| lotus_id | string | 1:1 | Yes | LNPN internal identifier | LOTUS0001234 |
| inchi2D | string | 1:1 | No | 2D InChI (no stereochemistry) | InChI=1S/... |
| inchikey2D | string | 1:1 | No | 2D InChI key | XXXX-XXXXXX-X |
| smiles2D | string | 1:1 | No | 2D SMILES | Flat structure |
| heavy_atom_number | integer | 1:1 | No | Heavy atom count | 30 |
| npl_score | float | 1:1 | No | Natural product likeness | 0.85 |
| fsp3 | float | 1:1 | No | Fraction sp3 carbons | 0.65 |
| lipinskiRuleOf5Failures | integer | 1:1 | No | Lipinski violations | 0 |
| deep_smiles | string | 1:1 | No | DeepSMILES representation | Alternative encoding |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Wikidata QID | Q + digits | Q312879 | Universal identifier |
| InChI Key | 27 characters | KZJWDPNRJALLNS-VJSFXXLFSA-N | Structure hash |
| NCBI Taxon ID | Integer | 3702 | Organism identifier |
| OTT ID | Integer | 1054861 | Open Tree taxonomy |
| DOI | 10.xxxx/... | 10.1016/j.phytochem.2018.01.001 | Publication identifier |
| PMID | Integer | 29371045 | PubMed article ID |

---

## Enumerations

### Wikidata Properties

| Property | Name | Description |
|----------|------|-------------|
| P235 | InChIKey | Chemical structure identifier |
| P234 | InChI | Full InChI string |
| P233 | Canonical SMILES | Canonical structure |
| P2017 | Isomeric SMILES | Stereochemistry-aware |
| P703 | found in taxon | Compound-organism link |
| P248 | stated in | Bibliographic reference |
| P31 | instance of | Type classification |
| P171 | parent taxon | Taxonomic parent |
| P105 | taxon rank | Species, genus, etc. |
| P225 | taxon name | Scientific name |
| P685 | NCBI taxon ID | NCBI Taxonomy |
| P9157 | Open Tree of Life ID | OTT identifier |
| P356 | DOI | Digital Object Identifier |
| P698 | PubMed ID | PubMed article ID |

### Taxonomic Ranks

| Rank | Description |
|------|-------------|
| species | Lowest standard rank |
| genus | Group of species |
| family | Group of genera |
| order | Group of families |
| class | Group of orders |
| phylum | Group of classes |
| kingdom | Highest standard rank |

### Mapping Predicates

| Predicate | Meaning |
|-----------|---------|
| exactMatch | 1:1 equivalence |
| closeMatch | Very similar |
| narrowMatch | More specific |
| broadMatch | More general |

---

## Entity Relationships

### Compound to Organism (P703)
- **Cardinality:** N:M
- **Description:** Compounds found in multiple organisms; organisms produce multiple compounds
- **Key Fields:** compoundWikidataId, taxonWikidataId

### Compound to Reference (P248)
- **Cardinality:** N:M
- **Description:** Compound-organism pairs documented in literature
- **Key Fields:** compoundWikidataId, referenceWikidataId

### Organism to Parent Taxon (P171)
- **Cardinality:** N:1
- **Description:** Hierarchical taxonomic relationship
- **Key Fields:** taxonWikidataId, parentWikidataId

### Structure-Organism-Reference Triplet
- **Cardinality:** Core model
- **Description:** Every documented occurrence requires all three elements
- **Key Fields:** compound, taxon, reference

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| LOTUS | Linked Open Total Unified Structure-organism | Database name |
| LNPN | LOTUS Natural Products Network | Web interface |
| SSOT | Single Source of Truth | PostgreSQL canonical |
| CC0 | Creative Commons Zero | License (public domain) |
| SMILES | Simplified Molecular Input Line Entry System | Structure notation |
| InChI | International Chemical Identifier | IUPAC standard |
| RDF | Resource Description Framework | Wikidata format |
| SPARQL | SPARQL Protocol and RDF Query Language | Query language |
| IDSM | Integrated Database for Structure and Metabolism | Substructure search |
| OTT | Open Tree of Life Taxonomy | Taxonomy system |
| NCBI | National Center for Biotechnology Information | Taxonomy source |
| DOI | Digital Object Identifier | Publication identifier |
| PMID | PubMed Identifier | Literature reference |
| QID | Q Identifier | Wikidata entity ID |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PubChem | InChI Key | Chemical data |
| ChEMBL | InChI Key, ChEMBL ID | Bioactivity |
| COCONUT | InChI Key | Natural products |
| UniProt | NCBI Taxon ID | Protein targets |
| CrossRef | DOI | Publication metadata |
| PubMed | PMID | Literature |
| Open Tree of Life | OTT ID | Taxonomy |

---

## Data Quality Notes

1. **Validation Pipeline:** 97% true positive rate after manual validation
2. **Retention Rate:** ~30% after validation (2.5M to 750K)
3. **Cardinality:** Every compound-organism pair requires literature reference
4. **Wikidata Integration:** Live data with continuous community curation
5. **License:** CC0 (Public Domain) - completely unrestricted
6. **SPARQL Access:** Full query capability via Wikidata endpoint
