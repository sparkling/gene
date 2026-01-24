---
id: reference-sample-data
title: "Pathway Database Sample Data"
category: reference
parent: _index.md
last_updated: 2026-01-23
status: active
migrated_from: operations/schemas/sample-data.md
tags: [schema, sample-data, reactome, wikipathways, disgenet, api-response]
---

**Parent:** [Architecture Documentation](./_index.md)

# Pathway Database Sample Data

This document contains actual sample data retrieved from Reactome, WikiPathways, and DisGeNET APIs.

---

## 1. Reactome Sample Data

### Pathway: Metabolism (R-HSA-1430728)

**API Call:** \`GET https://reactome.org/ContentService/data/query/R-HSA-1430728\`

\`\`\`json
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
  "speciesName": "Homo sapiens",
  "goBiologicalProcess": {
    "dbId": 504,
    "displayName": "metabolic process",
    "accession": "0008152",
    "databaseName": "GO"
  },
  "doi": "10.3180/R-HSA-1430728.15",
  "hasDiagram": true,
  "schemaClass": "TopLevelPathway"
}
\`\`\`

### Pathway: Cholesterol Biosynthesis (R-HSA-191273)

\`\`\`json
{
  "dbId": 191273,
  "displayName": "Cholesterol biosynthesis",
  "stId": "R-HSA-191273",
  "speciesName": "Homo sapiens",
  "compartment": [
    {"displayName": "cytosol", "accession": "0005829"},
    {"displayName": "endoplasmic reticulum membrane", "accession": "0005789"}
  ],
  "goBiologicalProcess": {
    "displayName": "cholesterol biosynthetic process",
    "accession": "0006695"
  },
  "hasDiagram": true,
  "schemaClass": "Pathway"
}
\`\`\`

---

## 2. WikiPathways Sample Data

### Pathway Search: Cholesterol Pathways

**API Call:** \`GET https://webservice.wikipathways.org/findPathwaysByText?query=cholesterol&species=Homo%20sapiens&format=json\`

\`\`\`json
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
\`\`\`

---

## 3. DisGeNET Sample Data

### Gene-Disease Association

\`\`\`json
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
  "NofPmids": 245
}
\`\`\`

### Variant-Disease Association

\`\`\`json
{
  "snpId": "rs3846662",
  "chromosome": "5",
  "position": 75360714,
  "diseaseId": "C0020443",
  "diseaseName": "Hypercholesterolemia",
  "score": 0.65,
  "consequence": "intron_variant",
  "geneSymbol": "HMGCR",
  "source": "GWAS Catalog",
  "gnomAD_AF": 0.42
}
\`\`\`

---

## 4. Cross-Database Relationships

### Cholesterol Biosynthesis Across Databases

| Entity | Reactome | WikiPathways | DisGeNET |
|--------|----------|--------------|----------|
| Pathway | R-HSA-191273 | WP197 | - |
| HMGCR Gene | UniProt:P04035 | Entrez:3156 | geneId:3156 |
| Cholesterol | ChEBI:16113 | CAS:57-88-5 | - |
| Hypercholesterolemia | - | - | C0020443 |

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

| Database | Version | Pathways | Genes | Diseases |
|----------|---------|----------|-------|----------|
| Reactome | 95 | 2,712 | 11,196 | - |
| WikiPathways | 2026-01 | 3,100+ | - | - |
| DisGeNET | 7.0 | - | 17,549 | 24,166 |

---

## Glossary

| Term | Definition |
|------|------------|
| stId | Reactome stable identifier |
| dbId | Reactome internal database ID |
| score | DisGeNET gene-disease association confidence (0-1) |
| EI | Evidence Index - publication consistency |
| DSI | Disease Specificity Index |
| DPI | Disease Pleiotropy Index |

---

*Full content preserved from original source at operations/schemas/sample-data.md*
