---
id: downloads-pharmaceuticals
title: "Pharmaceutical and Drug Database Bulk Downloads"
type: download-guide
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [downloads, bulk-data, pharmaceuticals, drugs]
---

**Parent:** [Download Guides](./_index.md)

# Pharmaceutical and Drug Database Bulk Downloads

This document provides comprehensive information on bulk download methods for major pharmaceutical and drug databases, including direct URLs, file sizes, licensing requirements, and recommended processing approaches.

---

## Table of Contents

1. [DrugBank](#1-drugbank)
2. [ChEMBL](#2-chembl)
3. [PubChem](#3-pubchem)
4. [PharmGKB](#4-pharmgkb)
5. [DGIdb](#5-dgidb)
6. [Open Targets](#6-open-targets)
7. [BindingDB](#7-bindingdb)
8. [OpenFDA](#8-openfda)
9. [DailyMed](#9-dailymed)
10. [RxNorm](#10-rxnorm)

---

[Full content from pharmaceuticals.md continues...]

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **DrugBank** | Download | `https://go.drugbank.com/releases/latest` (registration) |
| **ChEMBL** | FTP | `ftp://ftp.ebi.ac.uk/pub/databases/chembl/` |
| **PubChem** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/pubchem/` |
| **PharmGKB** | Download | `https://www.pharmgkb.org/downloads` |
| **DGIdb** | Download | `https://www.dgidb.org/downloads` |
| **Open Targets** | Download | `https://platform.opentargets.org/downloads` |
| **BindingDB** | Download | `https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFdownload.jsp` |
| **OpenFDA** | API | `https://open.fda.gov/apis/` |
| **DailyMed** | Download | `https://dailymed.nlm.nih.gov/dailymed/spl-resources.cfm` |

**Access Requirements:** Most are freely accessible; DrugBank full data requires academic license.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | SDF, TSV, JSON |
| Alternative | CSV, XML, SQLite |
| Chemical structures | SMILES, InChI, MOL |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `drug_id` | string | Drug identifier | "DB00945" |
| `name` | string | Generic name | "Aspirin" |
| `targets` | array | Target proteins | ["PTGS1", "PTGS2"] |
| `interactions` | array | Drug interactions | ["warfarin"] |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `targets` | Protein | N:M |
| `interacts_with` | Drug | N:M |

## Sample Data

### Example Drug Record
```json
{
  "drugbank_id": "DB00945",
  "name": "Aspirin",
  "chembl_id": "CHEMBL25",
  "targets": [
    {"uniprot": "P23219", "gene": "PTGS1", "action": "inhibitor"}
  ],
  "indication": "Pain and inflammation"
}
```

### Sample Query Result
| drug | target | action | Ki_nM |
|------|--------|--------|-------|
| Aspirin | PTGS1 | inhibitor | 1.67 |
| Ibuprofen | PTGS2 | inhibitor | 6000 |

## License

| Source | License | Commercial Use |
|--------|---------|----------------|
| DrugBank | Academic | License required |
| ChEMBL | CC BY-SA 3.0 | Yes |
| PubChem | Public domain | Yes |
| PharmGKB | CC BY-SA 4.0 | Yes |
| OpenFDA | Public domain | Yes |

## Data Set Size

| Metric | Value |
|--------|-------|
| ChEMBL compounds | 2.4M+ |
| PubChem compounds | 115M+ |
| DrugBank drugs | 15K+ |
| PharmGKB annotations | 20K+ |
| Total storage estimate | ~50-100 GB combined |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Bulk Download | Retrieving entire datasets via FTP/HTTP | Downloading ChEMBL SQLite |
| Drug Database | Structured collection of pharmaceutical compound information | DrugBank, ChEMBL |
| Licensing | Terms governing data access and usage | Academic vs commercial |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| DrugBank | Comprehensive drug and target database | Drug information |
| ChEMBL | Bioactivity database of drug-like compounds | EMBL-EBI |
| PubChem | NIH chemical compound database | Compound structures |
| PharmGKB | Pharmacogenomics knowledge base | Drug-gene relationships |
| DGIdb | Drug-Gene Interaction Database | Drug targets |
| Open Targets | Platform for drug target identification | Target validation |
| BindingDB | Database of binding affinities | Drug-target interactions |
| OpenFDA | FDA open data API | Adverse events |
| DailyMed | FDA drug labeling database | Package inserts |
| RxNorm | Normalized drug naming system | Drug nomenclature |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST endpoints |
| EMBL | European Molecular Biology Laboratory | Research organization |
| FDA | Food and Drug Administration | US regulatory agency |
| FTP | File Transfer Protocol | Download method |
| NIH | National Institutes of Health | US agency |
| SDF | Structure Data File | Chemical structures |
| SQL | Structured Query Language | Database queries |

---

*Last Updated: January 2026*
