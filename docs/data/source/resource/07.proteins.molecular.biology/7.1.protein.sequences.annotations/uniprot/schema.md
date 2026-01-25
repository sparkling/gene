---
id: schema-uniprot
title: "UniProt Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: final
tags: [schema, database, proteins, sequences, annotations]
---

# UniProt Schema Documentation

**Document ID:** SCHEMA-UNIPROT
**Version:** 2026.01
**Source Version:** Release 2026_01

---

## TL;DR

UniProt provides the most comprehensive protein sequence and annotation database with two main sections: Swiss-Prot (manually curated, 570K+ entries) and TrEMBL (automatically annotated, 250M+ entries). Data is available in multiple formats (FASTA, XML, RDF, TSV) with extensive cross-references to 180+ external databases.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Swiss-Prot entries | 570,000+ | UniProt Release Notes |
| TrEMBL entries | 250,000,000+ | UniProt Release Notes |
| Organisms | 20,000+ | UniProt Statistics |
| Human proteins (reviewed) | 20,400 | UniProt Human Proteome |
| Cross-reference databases | 180+ | UniProt Documentation |

---

## Entity Relationship Overview

```
UniProtKB Entry
  ├── Protein Names (recommended, alternative, submitted)
  ├── Gene Names (primary, synonyms, ordered locus, ORF)
  ├── Organism (taxonomy, lineage)
  ├── Protein Existence (PE level 1-5)
  ├── Sequence (canonical, isoforms)
  │     ├── Features (domains, sites, variants)
  │     └── Processing (signal, propeptide, chain)
  ├── Function
  │     ├── GO Terms (molecular function, biological process, cellular component)
  │     ├── Catalytic Activity (EC numbers, Rhea reactions)
  │     └── Pathways (Reactome, KEGG)
  ├── Subcellular Location
  ├── Disease Associations
  ├── PTM/Processing
  ├── Interaction Partners
  └── Cross-References (180+ databases)
```

---

## Core Tables/Entities

### Entry (Primary Record)

**Description:** Core protein entry containing all annotation data.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| accession | string | Yes | Primary accession (e.g., P04637) |
| id | string | Yes | Entry name (e.g., P53_HUMAN) |
| reviewed | boolean | Yes | Swiss-Prot (true) or TrEMBL (false) |
| proteinExistence | enum | Yes | Evidence level (1-5) |
| sequence | object | Yes | Amino acid sequence and metadata |
| organism | object | Yes | Source organism information |
| proteinDescription | object | Yes | Recommended and alternative names |
| genes | array | No | Gene name information |
| comments | array | No | Functional annotations |
| features | array | No | Sequence features |
| dbReferences | array | No | Cross-references |
| keywords | array | No | UniProt keywords |
| references | array | No | Literature citations |

### Sequence Object

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| value | string | Yes | Amino acid sequence |
| length | integer | Yes | Sequence length |
| molWeight | integer | Yes | Molecular weight (Da) |
| crc64 | string | Yes | CRC64 checksum |
| md5 | string | Yes | MD5 checksum |

### Feature Object

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| type | string | Yes | Feature type (domain, binding site, variant, etc.) |
| location | object | Yes | Start/end positions |
| description | string | No | Feature description |
| evidences | array | No | Evidence codes (ECO) |
| featureId | string | No | Feature identifier |

### Cross-Reference Object

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| database | string | Yes | Database name |
| id | string | Yes | External identifier |
| properties | array | No | Additional properties |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/uniprotkb/{accession}` | GET | Retrieve single entry |
| `/uniprotkb/search` | GET | Search entries |
| `/uniprotkb/stream` | GET | Bulk download |
| `/uniref/{id}` | GET | UniRef cluster |
| `/uniparc/{id}` | GET | UniParc archive entry |
| `/idmapping/run` | POST | Submit ID mapping job |
| `/idmapping/results/{jobId}` | GET | Get mapping results |

### Search Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| query | Search query | `gene:TP53 AND organism_id:9606` |
| fields | Return fields | `accession,id,protein_name,length` |
| format | Output format | `json`, `tsv`, `fasta`, `xml` |
| size | Results per page | `25` (max 500) |
| sort | Sort field | `accession asc` |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | UniProtKB flat file (.dat) |
| Alternative | XML, JSON, RDF/XML, FASTA, GFF, TSV |
| Compression | gzip |
| Encoding | UTF-8 |

---

## Sample Record

```json
{
  "primaryAccession": "P04637",
  "uniProtkbId": "P53_HUMAN",
  "entryType": "UniProtKB reviewed (Swiss-Prot)",
  "proteinExistence": "1: Evidence at protein level",
  "organism": {
    "scientificName": "Homo sapiens",
    "taxonId": 9606,
    "lineage": ["Eukaryota", "Metazoa", "Chordata", "Mammalia", "Primates", "Hominidae", "Homo"]
  },
  "proteinDescription": {
    "recommendedName": {
      "fullName": {"value": "Cellular tumor antigen p53"}
    },
    "alternativeNames": [
      {"fullName": {"value": "Tumor suppressor p53"}},
      {"fullName": {"value": "Phosphoprotein p53"}}
    ]
  },
  "genes": [
    {
      "geneName": {"value": "TP53"},
      "synonyms": [{"value": "P53"}]
    }
  ],
  "sequence": {
    "value": "MEEPQSDPSVEPPLSQETFSDLWKLL...",
    "length": 393,
    "molWeight": 43653,
    "crc64": "AD5C149FD8106131"
  },
  "features": [
    {
      "type": "Domain",
      "location": {"start": {"value": 102}, "end": {"value": 292}},
      "description": "p53 DNA-binding"
    }
  ],
  "comments": [
    {
      "commentType": "FUNCTION",
      "texts": [{"value": "Acts as a tumor suppressor in many tumor types..."}]
    }
  ],
  "uniProtKBCrossReferences": [
    {"database": "PDB", "id": "1TUP"},
    {"database": "RefSeq", "id": "NP_000537.3"},
    {"database": "Ensembl", "id": "ENSG00000141510"}
  ]
}
```

---

## Feature Types

| Type | Description | Example |
|------|-------------|---------|
| Signal | Signal peptide | 1-22 |
| Chain | Mature protein | 23-393 |
| Domain | Structural domain | DBD domain |
| Binding site | Ligand binding | DNA binding |
| Active site | Catalytic residue | Ser-315 |
| Metal binding | Metal coordination | Zinc binding |
| Modified residue | PTM sites | Phosphoserine |
| Natural variant | Polymorphism/mutation | p.R175H |
| Mutagenesis | Experimental mutation | Site-directed |
| Disulfide bond | Cysteine linkage | C176-C238 |

---

## Protein Existence Levels

| Level | Description | Evidence |
|-------|-------------|----------|
| 1 | Evidence at protein level | MS, X-ray, NMR |
| 2 | Evidence at transcript level | mRNA, EST |
| 3 | Inferred from homology | Sequence similarity |
| 4 | Predicted | Gene prediction |
| 5 | Uncertain | Dubious sequences |

---

## Cross-Reference Categories

| Category | Databases | Examples |
|----------|-----------|----------|
| Sequence | EMBL, RefSeq, CCDS | NM_000546 |
| 3D Structure | PDB, AlphaFold | 1TUP |
| Protein Families | Pfam, InterPro, SMART | PF00870 |
| PTM | PhosphoSitePlus, iPTMnet | Modification sites |
| Protein-Protein | STRING, IntAct, BioGRID | Interactions |
| Pathways | Reactome, KEGG, BioCyc | R-HSA-69206 |
| Disease | OMIM, DisGeNET, Orphanet | 191170 |
| Expression | Expression Atlas, Bgee | Tissue specificity |
| Organism-specific | HGNC, MGI, FlyBase | Gene nomenclature |

---

## Glossary

| Term | Definition |
|------|------------|
| Accession | Stable identifier for an entry (e.g., P04637) |
| Entry name | Mnemonic identifier (e.g., P53_HUMAN) |
| Swiss-Prot | Manually reviewed and annotated section |
| TrEMBL | Automatically annotated section |
| PE level | Protein existence evidence level (1-5) |
| Isoform | Alternative protein form from same gene |
| Evidence code | ECO ontology term for annotation source |

---

## References

1. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Research. https://doi.org/10.1093/nar/gkac1052
2. UniProt REST API Documentation: https://www.uniprot.org/help/api
3. UniProt User Manual: https://www.uniprot.org/help/
