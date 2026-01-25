---
id: schema-dgidb
title: "DGIdb Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, drug-gene-interactions, druggable-genome, graphql]
---

# DGIdb - Drug Gene Interaction Database Schema

**Document ID:** SCHEMA-DGIDB
**Version:** 5.0
**Source Version:** 2024

---

## TL;DR

DGIdb aggregates drug-gene interactions from 40+ sources into a unified queryable resource. The schema links drugs to genes with interaction types and source evidence, plus gene druggability categories. Accessible via GraphQL API for flexible queries supporting drug repurposing and target validation.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Drug-Gene Interactions | 90,000+ | Interaction records |
| Genes | 13,000+ | Gene entries |
| Drugs | 10,000+ | Drug entries |
| Source Databases | 40+ | Integrated sources |
| Gene Categories | 40+ | Druggability classes |
| Interaction Types | 30+ | Type annotations |

---

## Entity Relationship Overview

```
Drugs (1) ←→ (many) Interactions (many) ←→ (1) Genes
  ↓                      ↓                     ↓
Concept ID         Interaction type       Gene symbol

Interactions (many) ←→ (many) Sources
                           ↓
                    ChEMBL, DrugBank, etc.

Genes (1) ←→ (many) Gene_Categories
                         ↓
               KINASE, DRUGGABLE GENOME
```

---

## Core Tables/Entities

### genes

**Description:** Gene entries with identifiers and druggability info.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| concept_id | string | Yes | Normalized gene ID |
| name | string | Yes | Gene symbol (HGNC) |
| long_name | string | No | Full gene name |
| gene_categories | array | No | Druggability categories |
| entrez_id | integer | No | NCBI Entrez Gene ID |
| ensembl_id | string | No | Ensembl gene ID |

### drugs

**Description:** Drug/compound entries with identifiers.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| concept_id | string | Yes | Normalized drug ID |
| name | string | Yes | Drug name |
| approved | boolean | No | FDA approval status |
| chembl_id | string | No | ChEMBL identifier |
| drugbank_id | string | No | DrugBank identifier |
| aliases | array | No | Alternative names |

### interactions

**Description:** Drug-gene interaction records.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| id | string | Yes | Interaction identifier |
| gene | object | Yes | Gene reference |
| drug | object | Yes | Drug reference |
| interaction_types | array | No | Interaction type labels |
| interaction_score | decimal | No | Confidence score |
| sources | array | Yes | Source databases |
| pmids | array | No | PubMed references |

### gene_categories

**Description:** Gene druggability classifications.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| name | string | Yes | Category name |
| source_db_name | string | Yes | Source of category |
| genes | array | No | Genes in category |

### sources

**Description:** Source databases integrated in DGIdb.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| source_db_name | string | Yes | Database name |
| source_db_version | string | No | Version |
| citation | string | No | Reference citation |
| base_url | string | No | Database URL |

---

## Interaction Types

| Type | Description |
|------|-------------|
| inhibitor | Drug inhibits gene product |
| activator | Drug activates gene product |
| antagonist | Drug blocks receptor |
| agonist | Drug activates receptor |
| antibody | Antibody therapeutic |
| antisense | Antisense oligonucleotide |
| binder | Drug binds to target |
| modulator | Drug modulates activity |
| blocker | Drug blocks channel/transporter |
| vaccine | Vaccine target |

---

## Gene Categories

| Category | Description |
|----------|-------------|
| KINASE | Protein kinase genes |
| DRUGGABLE GENOME | Potentially tractable targets |
| CLINICALLY ACTIONABLE | Clinical relevance |
| TUMOR SUPPRESSOR | Cancer suppressor genes |
| TRANSCRIPTION FACTOR | TF genes |
| G PROTEIN COUPLED RECEPTOR | GPCR genes |
| ION CHANNEL | Ion channel genes |
| TRANSPORTER | Transporter genes |
| PROTEASE | Protease genes |

---

## GraphQL API

**Endpoint**: `https://dgidb.org/api/graphql`
**Version**: v5

### Core GraphQL Types

#### Gene
```graphql
type Gene {
  id: ID!
  name: String!
  conceptId: String!
  longName: String
  description: String
  geneCategories: [GeneCategory!]!
  geneAliases: [GeneAlias!]!
  geneClaims: [GeneClaim!]!
  interactions: [Interaction!]!
}
```

#### Drug
```graphql
type Drug {
  id: ID!
  name: String!
  conceptId: String!
  approved: Boolean!
  immunotherapy: Boolean!
  antineoplastic: Boolean!
  drugAliases: [DrugAlias!]!
  drugAttributes: [DrugAttribute!]!
  drugApplications: [DrugApplication!]!
  drugApprovalRatings: [DrugApprovalRating!]!
  interactions: [Interaction!]!
}
```

#### Interaction
```graphql
type Interaction {
  id: ID!
  gene: Gene!
  drug: Drug!
  interactionScore: Float
  interactionTypes: [InteractionType!]!
  interactionClaims: [InteractionClaim!]!
  publications: [Publication!]!
  sources: [Source!]!
}
```

#### InteractionType
```graphql
type InteractionType {
  id: ID!
  type: String!
  directionality: String
}
```

#### Source
```graphql
type Source {
  id: ID!
  fullName: String!
  sourceDbName: String!
  sourceDbVersion: String
  citationShort: String
  citationFull: String
  pmid: String
  pmcid: String
  doi: String
  sourceUrl: String
  license: String
  licenseLink: String
}
```

#### Publication
```graphql
type Publication {
  id: ID!
  pmid: Int
  citation: String
}
```

### Query Examples

```graphql
# Get interactions for a gene
query {
  gene(name: "EGFR") {
    name
    longName
    geneCategories {
      name
    }
    interactions {
      drug {
        name
        approved
      }
      interactionTypes {
        type
      }
      sources {
        sourceDbName
      }
    }
  }
}

# Search drugs
query {
  drugs(name: "gefitinib") {
    name
    chemblId
    approved
    interactions {
      gene {
        name
      }
      interactionTypes {
        type
      }
    }
  }
}

# Get druggable genes
query {
  genes(geneCategories: ["KINASE"]) {
    name
    interactionCount
  }
}

# Search for gene interactions (extended)
query {
  genes(names: ["EGFR"]) {
    nodes {
      name
      longName
      interactions {
        drug {
          name
          approved
        }
        interactionTypes {
          type
          directionality
        }
        interactionScore
      }
    }
  }
}

# Search by interaction type
query {
  interactions(
    interactionTypes: ["inhibitor"]
    geneName: "BRAF"
  ) {
    nodes {
      drug {
        name
        approved
      }
      gene {
        name
      }
      interactionTypes {
        type
      }
      publications {
        pmid
        citation
      }
    }
  }
}
```

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | GraphQL API |
| Alternative | TSV downloads |
| Encoding | UTF-8 |
| API Response | JSON |

---

## Sample Record

```json
{
  "gene": {
    "name": "EGFR",
    "longName": "Epidermal growth factor receptor",
    "conceptId": "hgnc:3236",
    "geneCategories": ["KINASE", "DRUGGABLE GENOME"]
  },
  "drug": {
    "name": "Gefitinib",
    "conceptId": "chembl:CHEMBL939",
    "approved": true,
    "chemblId": "CHEMBL939"
  },
  "interactionTypes": [
    {"type": "inhibitor", "directionality": "inhibits"}
  ],
  "sources": [
    {"sourceDbName": "ChEMBL"},
    {"sourceDbName": "DrugBank"}
  ],
  "pmids": [12345678, 23456789]
}
```

---

## Source Databases

| Source | Type | Coverage |
|--------|------|----------|
| ChEMBL | Bioactivity | Comprehensive |
| DrugBank | Drug info | Approved drugs |
| PharmGKB | PGx | Clinical |
| TTD | Targets | Therapeutic |
| DTC | Clinical trials | Recent |
| CIViC | Cancer | Somatic |
| OncoKB | Cancer | Actionable |
| Guide to Pharmacology | Pharmacology | Curated |

---

## Glossary

| Term | Definition |
|------|------------|
| Concept ID | Normalized identifier |
| Druggable genome | Genes amenable to targeting |
| Interaction type | Nature of drug-gene interaction |
| Source evidence | Supporting database |
| HGNC | HUGO Gene Nomenclature Committee |

---

## References

1. DGIdb: https://dgidb.org
2. Freshour SL, et al. (2021) Nucleic Acids Res. 49(D1):D1144-D1151
3. GraphQL API: https://dgidb.org/api/graphql
