# Data Dictionary: Category 07 - Proteins & Molecular Biology

## Overview

This data dictionary documents all fields across the Proteins & Molecular Biology category, which integrates data from protein sequence databases, structural repositories, and molecular interaction sources.

**Category ID:** 07
**Extraction Date:** 2026-01-24
**Subcategories:** 3

## Subcategories

| ID | Name | Data Sources |
|----|------|--------------|
| 7.1 | Protein Sequences & Annotations | UniProt, RefSeq |
| 7.2 | Protein Structures | PDB, AlphaFold DB, SWISS-MODEL |
| 7.3 | Molecular Interactions | IntAct, Reactome, STRING (via Category 04) |

---

## Global Cross-Category Mapping Fields

These fields serve as primary keys for linking data across subcategories and external databases.

| Field | Description | Used By | Semantic Definition |
|-------|-------------|---------|---------------------|
| `uniprot_accession` | UniProt accession serves as primary key linking sequences to structures | UniProt, RefSeq, PDB, AlphaFold DB, SWISS-MODEL | Universal protein identifier enabling cross-database integration |
| `gene_symbol` | Gene symbols link proteins to genomic context | UniProt, RefSeq, AlphaFold DB | Official gene nomenclature (HGNC for human) enabling gene-protein mapping |
| `pdb_id` | PDB identifier linking to 3D structural data | PDB, UniProt, RefSeq, SWISS-MODEL | 4-character alphanumeric PDB identifier for experimental structures |
| `taxonomy_id` | NCBI taxonomy ID for organism identification | UniProt, RefSeq, PDB, AlphaFold DB | Numeric identifier from NCBI taxonomy enabling organism-based queries |
| `sequence` | Amino acid sequence enabling sequence-based matching | UniProt, RefSeq, PDB, AlphaFold DB, SWISS-MODEL | One-letter amino acid sequence for alignment and identity calculations |

---

## Quality Metrics Summary

### Experimental Structures (PDB)

| Metric | Description |
|--------|-------------|
| Resolution | Diffraction/reconstruction resolution in Angstroms |
| R-factor | Crystallographic R-factor measuring model-data agreement |
| R-free | Cross-validation R-factor using test set reflections |
| Ramachandran outliers | Percentage of residues with unusual backbone angles |
| Clashscore | Steric clash score for model validation |

### Predicted Structures - AI (AlphaFold DB)

| Metric | Description |
|--------|-------------|
| pLDDT | Per-residue confidence score (0-100) |
| PAE | Predicted Aligned Error matrix between residue pairs |
| pTM | Global fold quality metric predicting TM-score |

### Predicted Structures - Homology (SWISS-MODEL)

| Metric | Description |
|--------|-------------|
| QMEAN | Global model quality Z-score (-4 to 0) |
| QMEANDisCo | Distance constraint-based quality score (0-1) |
| Sequence Identity | Fraction identical residues to template |
| Coverage | Fraction of target sequence modeled |

### Sequence Evidence (UniProt)

| Metric | Description |
|--------|-------------|
| Protein Existence Level | Evidence level 1-5 (1=protein level, 5=uncertain) |
| Reviewed Status | Swiss-Prot (reviewed) vs TrEMBL (unreviewed) |

---

## Cross-Reference Mappings

### 7.1 Protein Sequences & Annotations

| From Source | To Source | From Field | To Field | Mapping Type | Notes |
|-------------|-----------|------------|----------|--------------|-------|
| UniProt | RefSeq | accession | accession | 1:N | One UniProt entry may map to multiple RefSeq proteins (isoforms) |
| UniProt | RefSeq | dbReferences[database=RefSeq].id | accession | 1:N | Direct cross-reference via UniProt dbReferences |

### 7.2 Protein Structures

| From Source | To Source | From Field | To Field | Mapping Type | Notes |
|-------------|-----------|------------|----------|--------------|-------|
| PDB | UniProt | entity_poly.pdbx_seq_one_letter_code | sequence.value | N:1 | Multiple PDB entries may represent same UniProt protein |
| AlphaFold DB | UniProt | uniprotAccession | accession | 1:1 | Direct mapping via UniProt accession |
| SWISS-MODEL | UniProt | uniprot_ac | accession | N:1 | Multiple models may exist per UniProt entry |
| SWISS-MODEL | PDB | template | entry.id | N:1 | Models reference PDB templates |
| AlphaFold DB | PDB | none | none | comparison | AlphaFold predictions can be compared to experimental PDB structures |

---

## Subcategory Data Dictionaries

For detailed field documentation, see:

- [7.1 Protein Sequences & Annotations Data Dictionary](./7.1.protein.sequences.annotations/_data-dictionary.md)
- [7.2 Protein Structures Data Dictionary](./7.2.protein.structures/_data-dictionary.md)
- [7.3 Molecular Interactions Data Dictionary](./7.3.molecular.interactions/_data-dictionary.md)

---

## Related Categories

| Category | Relationship |
|----------|-------------|
| 04. Pathways & Networks | Molecular interactions (IntAct, Reactome, STRING) |
| 01. Gene Resources | Gene-protein mappings via gene symbols and IDs |
| 06. Ontologies | GO annotations for proteins |
