# UniProt - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | uniprot |
| **Name** | Universal Protein Resource |
| **Total Fields** | 100+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| primaryAccession | string | 1:1 | Yes | UniProt accession | `P04637` |
| uniProtkbId | string | 1:1 | Yes | Entry name | `P53_HUMAN` |
| entryType | string | 1:1 | Yes | Review status | `UniProtKB reviewed` |
| proteinExistence | string | 1:1 | No | Evidence level | `1: Evidence at protein level` |

### Sequence

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| sequence.value | string | 1:1 | Yes | Amino acid sequence | `MEEPQSDPSV...` |
| sequence.length | integer | 1:1 | Yes | Sequence length | `393` |
| sequence.molWeight | integer | 1:1 | No | Molecular weight (Da) | `43653` |
| sequence.crc64 | string | 1:1 | No | CRC64 checksum | `AD5C149FD8106131` |

### Organism

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| organism.scientificName | string | 1:1 | Yes | Scientific name | `Homo sapiens` |
| organism.taxonId | integer | 1:1 | Yes | NCBI Taxonomy ID | `9606` |
| organism.lineage | array | 1:N | No | Taxonomic lineage | `Eukaryota, Metazoa...` |

---

## Enumerations

### Entry Type

| Value | Description |
|-------|-------------|
| UniProtKB reviewed (Swiss-Prot) | Manually reviewed and annotated |
| UniProtKB unreviewed (TrEMBL) | Automatically annotated |

### Protein Existence Level

| Level | Description | Evidence |
|-------|-------------|----------|
| 1 | Evidence at protein level | MS, X-ray, NMR |
| 2 | Evidence at transcript level | mRNA, EST |
| 3 | Inferred from homology | Sequence similarity |
| 4 | Predicted | Gene prediction |
| 5 | Uncertain | Dubious sequences |

### Feature Types

| Type | Description |
|------|-------------|
| Signal | Signal peptide |
| Chain | Mature protein |
| Domain | Structural domain |
| Binding site | Ligand binding |
| Active site | Catalytic residue |
| Modified residue | PTM site |
| Natural variant | Polymorphism/mutation |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| UniProt | Universal Protein Resource | Consortium |
| Swiss-Prot | Swiss Protein | Reviewed section |
| TrEMBL | Translated EMBL | Unreviewed section |
| PE | Protein Existence | Evidence level |
| PTM | Post-Translational Modification | Chemical change |
| ECO | Evidence & Conclusion Ontology | Evidence codes |
| GO | Gene Ontology | Functional annotation |

---

## Data Quality Notes

1. **Swiss-Prot vs TrEMBL**: Swiss-Prot is manually curated, TrEMBL automatic
2. **Evidence codes**: All annotations include ECO evidence
3. **180+ cross-references**: Links to external databases
4. **Weekly updates**: New release every week

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
