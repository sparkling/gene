---
id: schema-orphanet
title: "Orphanet Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: migrated
tags: [schema, database, rare-diseases, orphan-drugs, ordo, clinical-resources]
---

# Orphanet Schema Documentation

**Document ID:** SCHEMA-ORPHANET
**Version:** 1.0
**Source Version:** Orphadata January 2026

---

## TL;DR

Orphanet provides structured rare disease data with 6,500+ diseases, gene associations, phenotypes (HPO), and epidemiology. Data available as XML products via Orphadata or the ORDO ontology in OWL format. ORPHA codes are the primary disease identifiers with extensive cross-references to OMIM, ICD, and other terminologies.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Rare Diseases | 6,528 | Product 1 |
| Genes with Associations | 4,512 | Product 6 |
| Diagnostic Tests | 36,595 | Product 7 |
| Expert Centres | 8,722 | Product 5 |
| Phenotype Annotations | 100,000+ | Product 4 |
| Patient Organizations | 3,500+ | Product 5 |

---

## Entity Relationship Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                      RARE DISEASE                                │
│  ORPHAcode, Name, ORPHAgroup, SynonymList                       │
│  Classification, DisorderType, AverageAgeOfOnset                │
└─────────────────────────────────────────────────────────────────┘
         │                    │                    │
         │ Gene Assoc.       │ Phenotypes         │ Epidemiology
         ▼                    ▼                    ▼
┌──────────────────┐ ┌──────────────────┐ ┌──────────────────────┐
│      GENE        │ │    HPO TERM      │ │   EPIDEMIOLOGY       │
│  Symbol, HGNC    │ │  HP:XXXXXXX      │ │  Prevalence          │
│  GeneLocus       │ │  Frequency       │ │  Incidence           │
│  AssocType       │ │                  │ │  Geographic          │
└──────────────────┘ └──────────────────┘ └──────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                    EXTERNAL REFERENCES                           │
│  OMIM, ICD-10, ICD-11, UMLS, MeSH, MedDRA, SNOMED CT           │
└─────────────────────────────────────────────────────────────────┘
```

---

## Core Tables/Entities

### Disorder (Disease)

**Description:** Rare disease entity with core attributes

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ORPHAcode | integer | Yes | Primary Orphanet identifier |
| Name | string | Yes | Preferred disease name |
| ORPHAgroup | enum | Yes | Disorder, Group, Subtype |
| DisorderType | object | Yes | Disease/malformation/etc. |
| SynonymList | list | No | Alternative names |
| ExternalReferenceList | list | No | Cross-references |
| TextualInformation | object | No | Definition, description |
| AverageAgeOfOnset | list | No | Onset age categories |
| AverageAgeOfDeath | list | No | Death age categories |
| TypeOfInheritanceList | list | No | Inheritance patterns |

### Gene Association

**Description:** Gene linked to a rare disease

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ORPHAcode | integer | Yes | Disease identifier |
| Gene | object | Yes | Gene information |
| Symbol | string | Yes | HGNC symbol |
| Name | string | Yes | Gene name |
| GeneLocus | string | No | Chromosomal location |
| HGNC | integer | No | HGNC ID |
| DisorderGeneAssociationType | object | Yes | Association type |
| DisorderGeneAssociationStatus | object | Yes | Validated/Not validated |

### HPO Phenotype Annotation

**Description:** Clinical feature linked to disease

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ORPHAcode | integer | Yes | Disease identifier |
| HPO | object | Yes | HPO term information |
| HPOId | string | Yes | HP:XXXXXXX identifier |
| HPOTerm | string | Yes | Phenotype name |
| HPOFrequency | object | No | Frequency classification |
| Diagnostic | boolean | No | Diagnostic criterion |

### Epidemiological Data

**Description:** Disease prevalence and incidence

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| ORPHAcode | integer | Yes | Disease identifier |
| PrevalenceList | list | No | Prevalence records |
| PrevalenceType | string | No | Point/Birth/Lifetime |
| PrevalenceQualification | string | No | Mean/Range/Class |
| PrevalenceClass | string | No | Categorical prevalence |
| ValMoy | float | No | Mean prevalence value |
| PrevalenceGeographic | string | No | Geographic scope |
| PrevalenceValidationStatus | string | No | Validated status |

---

## Data Products (Orphadata)

| Product | File | Description |
|---------|------|-------------|
| Product 1 | en_product1.xml | Disease nomenclature |
| Product 3 | en_product3.xml | Disease classifications |
| Product 4 | en_product4.xml | HPO phenotype annotations |
| Product 6 | en_product6.xml | Gene-disease associations |
| Product 7 | en_product7.xml | Diagnostic tests |
| Product 9 | en_product9.xml | Epidemiological data |
| ORDO | ordo.owl | Ontology (OWL format) |

---

## Identifier Format

### ORPHA Code

```
Orphanet:558
ORPHA:558
```

### ORDO URI

```
http://www.orpha.net/ORDO/Orphanet_558
```

### Pattern

```
Orphanet:[0-9]+
```

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /rd-api/disorder/{orphacode} | GET | Disease details |
| /rd-api/disorder/search | GET | Search diseases |
| /rd-api/disorder/{orphacode}/genes | GET | Gene associations |
| /rd-api/disorder/{orphacode}/phenotypes | GET | HPO phenotypes |
| /rd-api/gene/{hgnc}/disorders | GET | Diseases for gene |

### API Example

```bash
curl "https://api.orphanet.org/rd-api/disorder/558" \
  -H "Accept: application/json"
```

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | XML (Orphadata products) |
| Alternative | OWL (ORDO ontology) |
| API | JSON responses |
| Encoding | UTF-8 |

---

## Sample Records

### Disease (XML)

```xml
<Disorder id="558">
  <OrphaCode>558</OrphaCode>
  <Name lang="en">Marfan syndrome</Name>
  <ORPHAcode>558</ORPHAcode>
  <DisorderType>
    <Name lang="en">Disease</Name>
  </DisorderType>
  <SynonymList count="3">
    <Synonym lang="en">MFS</Synonym>
    <Synonym lang="en">Marfan's syndrome</Synonym>
  </SynonymList>
  <ExternalReferenceList count="4">
    <ExternalReference>
      <Source>OMIM</Source>
      <Reference>154700</Reference>
    </ExternalReference>
    <ExternalReference>
      <Source>ICD-10</Source>
      <Reference>Q87.4</Reference>
    </ExternalReference>
  </ExternalReferenceList>
</Disorder>
```

### Gene Association (XML)

```xml
<DisorderGeneAssociation>
  <SourceOfValidation>PMID:1852208</SourceOfValidation>
  <Gene id="5714">
    <Symbol>FBN1</Symbol>
    <Name lang="en">fibrillin 1</Name>
    <HGNC>3603</HGNC>
    <GeneLocus>15q21.1</GeneLocus>
  </Gene>
  <DisorderGeneAssociationType id="1">
    <Name lang="en">Disease-causing germline mutation(s) in</Name>
  </DisorderGeneAssociationType>
  <DisorderGeneAssociationStatus id="1">
    <Name lang="en">Assessed</Name>
  </DisorderGeneAssociationStatus>
</DisorderGeneAssociation>
```

### HPO Annotation (XML)

```xml
<HPODisorderAssociation>
  <HPO id="14213">
    <HPOId>HP:0001166</HPOId>
    <HPOTerm>Arachnodactyly</HPOTerm>
  </HPO>
  <HPOFrequency id="28440">
    <Name lang="en">Very frequent (99-80%)</Name>
  </HPOFrequency>
</HPODisorderAssociation>
```

### ORDO Term (OWL/Turtle)

```turtle
ordo:Orphanet_558 a owl:Class ;
    rdfs:label "Marfan syndrome"@en ;
    rdfs:subClassOf ordo:Orphanet_98249 ;
    oboInOwl:hasDbXref "OMIM:154700" ;
    oboInOwl:hasDbXref "ICD-10:Q87.4" ;
    obo:IAO_0000115 "A hereditary connective tissue disorder characterized by..." .
```

---

## Association Types

### Gene-Disease Association Types

| ID | Name | Description |
|----|------|-------------|
| 1 | Disease-causing germline mutation(s) | Direct causative |
| 2 | Disease-causing germline mutation(s) (loss of function) | LOF mechanism |
| 3 | Disease-causing germline mutation(s) (gain of function) | GOF mechanism |
| 4 | Major susceptibility factor | Strong risk factor |
| 5 | Role in phenotype | Modifier gene |
| 6 | Biomarker | Diagnostic marker |

### Inheritance Types

| ID | Name |
|----|------|
| 1 | Autosomal dominant |
| 2 | Autosomal recessive |
| 3 | X-linked dominant |
| 4 | X-linked recessive |
| 5 | Mitochondrial |
| 6 | Multigenic/multifactorial |

---

## Glossary

| Term | Definition |
|------|------------|
| ORPHAcode | Unique Orphanet disease identifier |
| ORDO | Orphanet Rare Disease Ontology |
| ORPHAgroup | Category: Disorder, Group of disorders, Subtype |
| Disorder | Single clinical entity |
| Group | Collection of related disorders |
| Subtype | Variant of a parent disorder |
| Orphadata | Scientific data download portal |
| HOOM | HPO-Orphanet Ontological Module |

---

## Cross-Reference Coverage

| External DB | Coverage | Notes |
|-------------|----------|-------|
| OMIM | High | Gene and phenotype |
| ICD-10 | High | WHO classification |
| ICD-11 | Growing | New WHO version |
| UMLS | High | CUI mappings |
| MeSH | Medium | NLM vocabulary |
| SNOMED CT | Medium | Clinical terms |
| MedDRA | Medium | Drug safety |
| GARD | High | NIH rare diseases |

---

## References

1. https://www.orpha.net/
2. https://www.orphadata.com/
3. https://www.ebi.ac.uk/ols4/ontologies/ordo
4. Rath et al. (2012) "Representation of rare diseases in health information systems" Orphanet J Rare Dis

---

## See Also

- [Orphanet Overview](./_index.md)
- [Orphanet Download Instructions](./download.md)
- [OMIM](../../3.2.phenotype.databases/omim/_index.md) - Mendelian inheritance database
- [HPO](../../3.2.phenotype.databases/hpo/_index.md) - Phenotype ontology
