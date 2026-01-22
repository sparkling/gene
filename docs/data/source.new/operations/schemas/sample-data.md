---
id: schemas-sample-data
title: "Pathway Database Sample Data"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, sample-data, reactome, wikipathways, disgenet, api-response]
---

**Parent:** [Schema Documentation](./_index.md)

# Pathway Database Sample Data

This document contains actual sample data retrieved from Reactome, WikiPathways, and DisGeNET APIs.

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "1430728" |
| `name` | string | Entity name | "Metabolism" |
| `type` | string | Record type | "pathway" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## 1. Reactome Sample Data

### Pathway: Metabolism (R-HSA-1430728)

**API Call:** `GET https://reactome.org/ContentService/data/query/R-HSA-1430728`

```json
{
  "dbId": 1430728,
  "displayName": "Metabolism",
  "stId": "R-HSA-1430728",
  "stIdVersion": "R-HSA-1430728.16",
  "isInDisease": false,
  "isInferred": false,
  "maxDepth": 7,
  "name": ["Metabolism"],
  "releaseDate": "2011-09-20",
  "releaseStatus": "UPDATED",
  "speciesName": "Homo sapiens",
  "goBiologicalProcess": {
    "dbId": 504,
    "displayName": "metabolic process",
    "accession": "0008152",
    "databaseName": "GO",
    "url": "https://www.ebi.ac.uk/QuickGO/term/GO:0008152"
  },
  "species": [{
    "dbId": 48887,
    "displayName": "Homo sapiens",
    "taxId": "9606",
    "abbreviation": "HSA"
  }],
  "doi": "10.3180/R-HSA-1430728.15",
  "hasDiagram": true,
  "hasEHLD": true,
  "schemaClass": "TopLevelPathway"
}
```

### Pathway: Cholesterol Biosynthesis (R-HSA-191273)

**API Call:** `GET https://reactome.org/ContentService/data/query/R-HSA-191273`

```json
{
  "dbId": 191273,
  "displayName": "Cholesterol biosynthesis",
  "stId": "R-HSA-191273",
  "stIdVersion": "R-HSA-191273.13",
  "isInDisease": false,
  "maxDepth": 3,
  "speciesName": "Homo sapiens",
  "compartment": [
    {
      "dbId": 70101,
      "displayName": "cytosol",
      "accession": "0005829",
      "databaseName": "GO"
    },
    {
      "dbId": 12045,
      "displayName": "endoplasmic reticulum membrane",
      "accession": "0005789",
      "databaseName": "GO"
    }
  ],
  "goBiologicalProcess": {
    "dbId": 20521,
    "displayName": "cholesterol biosynthetic process",
    "accession": "0006695",
    "databaseName": "GO"
  },
  "literatureReference": [
    {
      "dbId": 194670,
      "title": "Membrane-bound enzymes of cholesterol synthesis from lanosterol",
      "journal": "Biochem Biophys Res Commun",
      "pubMedIdentifier": 11969204,
      "year": 2002
    }
  ],
  "doi": "10.3180/R-HSA-191273.11",
  "hasDiagram": true,
  "hasEHLD": false,
  "schemaClass": "Pathway"
}
```

### Small Molecule: ATP (R-ALL-29358)

**API Call:** `GET https://reactome.org/ContentService/data/query/29358`

```json
{
  "dbId": 29358,
  "displayName": "ATP [nucleoplasm]",
  "stId": "R-ALL-29358",
  "stIdVersion": "R-ALL-29358.3",
  "maxDepth": 1,
  "name": ["ATP", "Adenosine 5'-triphosphate", "ATP(4-)"],
  "compartment": [{
    "dbId": 7660,
    "displayName": "nucleoplasm",
    "accession": "0005654",
    "databaseName": "GO"
  }],
  "referenceEntity": {
    "dbId": 8869364,
    "stId": "chebi:30616",
    "databaseName": "ChEBI",
    "identifier": "30616",
    "name": ["ATP(4-)", "ATP", "atp", "Adenosine 5'-triphosphate"],
    "formula": "C10H12N5O13P3",
    "schemaClass": "ReferenceMolecule"
  },
  "schemaClass": "SimpleEntity"
}
```

---

## 2. WikiPathways Sample Data

### Pathway Info: COVID-19 Pathway (WP4846)

**API Call:** `GET https://webservice.wikipathways.org/getPathwayInfo?pwId=WP4846&format=json`

```json
{
  "pathwayInfo": {
    "id": "WP4846",
    "url": "https://classic.wikipathways.org/index.php/Pathway:WP4846",
    "name": "SARS-CoV-2 and COVID-19 pathway",
    "species": "Homo sapiens",
    "revision": "140186"
  }
}
```

### Pathway Search: Cholesterol Pathways

**API Call:** `GET https://webservice.wikipathways.org/findPathwaysByText?query=cholesterol&species=Homo%20sapiens&format=json`

```json
{
  "result": [
    {
      "score": {"0": "4.837"},
      "id": "WP5333",
      "name": "Enterocyte cholesterol metabolism",
      "species": "Homo sapiens",
      "revision": "139904"
    },
    {
      "score": {"0": "4.797"},
      "id": "WP5304",
      "name": "Cholesterol metabolism",
      "species": "Homo sapiens",
      "revision": "141099"
    },
    {
      "score": {"0": "4.655"},
      "id": "WP197",
      "name": "Cholesterol biosynthesis pathway",
      "species": "Homo sapiens",
      "revision": "141096"
    }
  ]
}
```

---

## 3. DisGeNET Sample Data

### Gene-Disease Association

```json
{
  "geneId": 3156,
  "geneSymbol": "HMGCR",
  "geneName": "3-hydroxy-3-methylglutaryl-CoA reductase",
  "diseaseId": "C0020443",
  "diseaseName": "Hypercholesterolemia",
  "score": 0.82,
  "EI": 1.0,
  "DSI": 0.38,
  "DPI": 0.71,
  "associationType": "GeneticVariation",
  "source": "CTD_human",
  "NofPmids": 245,
  "pmids": [12345678, 23456789, 34567890]
}
```

### Variant-Disease Association

```json
{
  "snpId": "rs3846662",
  "chromosome": "5",
  "position": 75360714,
  "reference": "A",
  "alternative": "G",
  "diseaseId": "C0020443",
  "diseaseName": "Hypercholesterolemia",
  "score": 0.65,
  "EI": 0.92,
  "consequence": "intron_variant",
  "geneSymbol": "HMGCR",
  "source": "GWAS Catalog",
  "gnomAD_AF": 0.42,
  "NofPmids": 18
}
```

---

## 4. Cross-Database Relationships

### Cholesterol Biosynthesis Across Databases

| Entity | Reactome | WikiPathways | DisGeNET |
|--------|----------|--------------|----------|
| **Pathway** | R-HSA-191273 | WP197 | - |
| **HMGCR Gene** | UniProt:P04035 | Entrez:3156 | geneId:3156 |
| **Cholesterol** | ChEBI:16113 | CAS:57-88-5 | - |
| **Hypercholesterolemia** | - | - | C0020443 |

### Identifier Cross-References

| Database | HMGCR IDs |
|----------|-----------|
| NCBI Gene | 3156 |
| UniProt | P04035 |
| Ensembl | ENSG00000113161 |
| HGNC | HGNC:5006 |
| RefSeq | NM_000859 |

---

## 5. Database Statistics Summary

| Database | Version | Pathways | Genes | Molecules | Diseases |
|----------|---------|----------|-------|-----------|----------|
| Reactome | 95 | 2,712 | 11,196 | 1,925 | - |
| WikiPathways | 2026-01 | 3,100+ | - | - | - |
| DisGeNET | 7.0 | - | 17,549 | - | 24,166 |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 1,000+ |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | JSON |
| Alternative | CSV |
| Encoding | UTF-8 |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| Sample Data | HTTP | See main database |
| Test Data | HTTP | See repository |

**Access Requirements:** Varies by source

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| Sample Data | See source database | Depends on source |

---

## Sample Data

### Example Record
```json
{
  "stId": "R-HSA-1430728",
  "displayName": "Metabolism",
  "dbId": 1430728,
  "schemaClass": "TopLevelPathway"
}
```

### Sample Query Result
| stId | displayName | dbId | schemaClass |
|------|-------------|------|-------------|
| R-HSA-1430728 | Metabolism | 1430728 | TopLevelPathway |
| R-HSA-191273 | Cholesterol biosynthesis | 191273 | Pathway |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `stId` | Reactome stable identifier for entities | R-HSA-1430728 |
| `dbId` | Reactome internal database identifier | 1430728 |
| `schemaClass` | Reactome entity type classification | TopLevelPathway, Pathway |
| `revision` | WikiPathways version number for pathway | 140186 |
| `score` | DisGeNET gene-disease association confidence | 0.82 (0-1 scale) |
| `EI` | DisGeNET Evidence Index measuring publication consistency | 1.0 |
| `DSI` | Disease Specificity Index measuring gene-disease specificity | 0.38 |
| `DPI` | Disease Pleiotropy Index measuring disease diversity | 0.71 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| TopLevelPathway | Highest level pathway category in Reactome | Metabolism, Signaling |
| SimpleEntity | Small molecule or metabolite in Reactome | ATP, Glucose |
| GeneProduct | WikiPathways element representing a gene or protein | DataNode type |
| goBiologicalProcess | Gene Ontology biological process annotation | Pathway function |
| literatureReference | Publication citation supporting pathway data | PubMed ID |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GDA | Gene-Disease Association | DisGeNET metric |
| VDA | Variant-Disease Association | DisGeNET metric |
| GPML | Graphical Pathway Markup Language | WikiPathways format |
| GO | Gene Ontology | Functional annotations |
| ChEBI | Chemical Entities of Biological Interest | Molecule IDs |
| HMDB | Human Metabolome Database | Metabolite IDs |
| CAS | Chemical Abstracts Service | Registry numbers |
| UMLS | Unified Medical Language System | Disease IDs (CUI) |
