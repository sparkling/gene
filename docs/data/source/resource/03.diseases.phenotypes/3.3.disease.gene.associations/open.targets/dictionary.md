# Open Targets - Data Dictionary

## Overview

This data dictionary documents the schema for Open Targets Platform.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | open.targets |
| **Name** | Open Targets |
| **Parent** | 3.3.disease.gene.associations |
| **Total Fields** | 35+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Target (Gene)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | Ensembl gene ID | ENSG00000146648 |
| approvedSymbol | string | 1:1 | Yes | HGNC approved symbol | EGFR |
| approvedName | string | 1:1 | Yes | Full gene name | epidermal growth factor receptor |
| biotype | string | 1:1 | Yes | Gene biotype | protein_coding |
| transcriptIds | array | 1:N | No | Ensembl transcript IDs | [ENST00000275493] |
| proteinIds | object | 1:1 | No | UniProt mappings | {id: P00533} |
| tractability | object | 1:1 | No | Druggability assessment | {smallmolecule: {...}} |
| safety | object | 1:1 | No | Safety liabilities | {liabilities: [...]} |
| pathways | array | 1:N | No | Reactome pathways | [{id: R-HSA-1234}] |
| go | array | 1:N | No | GO annotations | [{id: GO:0005524}] |

### Disease

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | EFO ID | EFO_0000311 |
| name | string | 1:1 | Yes | Disease name | cancer |
| description | string | 1:1 | No | Definition | A disease characterized by... |
| therapeuticAreas | array | 1:N | No | Therapeutic categories | [EFO_0000651] |
| synonyms | array | 1:N | No | Alternative names | [neoplasm, tumor] |
| dbXRefs | object | 1:1 | No | Cross-references | {MONDO: [...], OMIM: [...]} |
| ancestors | array | 1:N | No | Parent EFO terms | [EFO_0000408] |
| descendants | array | 1:N | No | Child terms | [EFO_0000178] |

### Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| targetId | string | 1:1 | Yes | Ensembl gene ID | ENSG00000146648 |
| diseaseId | string | 1:1 | Yes | EFO disease ID | EFO_0000311 |
| score | float | 1:1 | Yes | Overall score (0-1) | 0.87 |
| datatypeScores | array | 1:N | No | Score per data type | [{id: known_drug, score: 1.0}] |
| evidenceCount | integer | 1:1 | Yes | Total evidence count | 1523 |

### Evidence

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | Evidence ID | d4ef8a7c-bc12... |
| targetId | string | 1:1 | Yes | Ensembl gene ID | ENSG00000146648 |
| diseaseId | string | 1:1 | Yes | EFO disease ID | EFO_0000311 |
| datasourceId | string | 1:1 | Yes | Source identifier | chembl |
| datatypeId | string | 1:1 | Yes | Data type category | known_drug |
| score | float | 1:1 | Yes | Evidence score (0-1) | 1.0 |
| literature | array | 1:N | No | PubMed references | [12345678] |

### Drug

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | ChEMBL ID | CHEMBL941 |
| name | string | 1:1 | Yes | Drug name | GEFITINIB |
| drugType | string | 1:1 | Yes | Drug modality | Small molecule |
| mechanismOfAction | string | 1:1 | No | MOA description | EGFR inhibitor |
| indications | array | 1:N | No | Approved indications | [EFO_0000311] |
| linkedTargets | array | 1:N | No | Drug targets | [ENSG00000146648] |
| phase | integer | 1:1 | No | Max clinical phase | 4 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Ensembl Gene | ENSG[0-9]{11} | ENSG00000146648 | Target identifier |
| EFO ID | EFO_[0-9]{7} | EFO_0000311 | Disease identifier |
| ChEMBL ID | CHEMBL[0-9]+ | CHEMBL941 | Drug identifier |
| Reactome | R-HSA-[0-9]+ | R-HSA-1234 | Pathway identifier |
| GO | GO:[0-9]{7} | GO:0005524 | Gene ontology |
| MONDO | MONDO:[0-9]{7} | MONDO:0005015 | Disease cross-ref |
| OMIM | OMIM:[0-9]{6} | OMIM:154700 | Disease cross-ref |

---

## Enumerations

### Data Types

| Data Type | Description | Sources |
|-----------|-------------|---------|
| genetic_association | Genetic evidence | GWAS Catalog, ClinVar, G2P |
| somatic_mutation | Cancer mutations | Cancer Gene Census, IntOGen |
| known_drug | Approved drugs | ChEMBL |
| affected_pathway | Pathway involvement | Reactome, SLAPenrich |
| rna_expression | Expression changes | Expression Atlas |
| literature | Text mining | Europe PMC |
| animal_model | Animal phenotypes | IMPC, MGI |

### Clinical Phases

| Phase | Description |
|-------|-------------|
| 0 | Exploratory IND |
| 1 | Safety trials |
| 2 | Efficacy trials |
| 3 | Pivotal trials |
| 4 | Post-marketing |

### Drug Types

| Type | Description |
|------|-------------|
| Small molecule | Chemical compounds |
| Antibody | Monoclonal antibodies |
| Protein | Therapeutic proteins |
| Oligonucleotide | ASOs, siRNAs |
| Cell therapy | CAR-T, stem cells |
| Gene therapy | Gene delivery |

### Tractability Categories

| Category | Description |
|----------|-------------|
| Clinical_Precedence_sm | Small molecule in clinic |
| Discovery_Precedence_sm | Small molecule in discovery |
| Predicted_Tractable_sm | Predicted druggable (SM) |
| Clinical_Precedence_ab | Antibody in clinic |
| Predicted_Tractable_ab | Predicted druggable (Ab) |

### Gene Biotypes

| Biotype | Description |
|---------|-------------|
| protein_coding | Protein-coding gene |
| lncRNA | Long non-coding RNA |
| miRNA | MicroRNA |
| pseudogene | Pseudogene |
| rRNA | Ribosomal RNA |

### Score Calculation

| Method | Description |
|--------|-------------|
| Harmonic sum | 1/(1/s1 + 1/s2 + ... + 1/sn) |
| Direct | Exact disease term match |
| Indirect | Including disease descendants |

---

## Entity Relationships

### Target to Diseases
- **Cardinality:** N:M
- **Description:** Genes associated with multiple diseases
- **Key Fields:** targetId, diseaseId

### Association to Evidence
- **Cardinality:** 1:N
- **Description:** Associations supported by multiple evidence
- **Key Fields:** targetId, diseaseId, evidenceId

### Drug to Targets
- **Cardinality:** N:M
- **Description:** Drugs target multiple genes
- **Key Fields:** drugId, targetId

### Disease to Ancestors
- **Cardinality:** N:M
- **Description:** EFO hierarchy
- **Key Fields:** diseaseId, ancestors

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| EFO | Experimental Factor Ontology | Disease IDs |
| ENSG | Ensembl Gene | Target IDs |
| MOA | Mechanism of Action | Drug action |
| GWAS | Genome-Wide Association Study | Evidence source |
| G2P | Gene2Phenotype | Evidence source |
| CGC | Cancer Gene Census | Evidence source |
| IMPC | International Mouse Phenotyping Consortium | Evidence source |
| MGI | Mouse Genome Informatics | Evidence source |
| PMC | PubMed Central | Literature source |
| SM | Small Molecule | Drug type |
| Ab | Antibody | Drug type |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| Ensembl | Gene ID | Target annotation |
| UniProt | Accession | Protein data |
| ChEMBL | ChEMBL ID | Drug data |
| EFO | EFO ID | Disease ontology |
| MONDO | MONDO ID | Disease mapping |
| OMIM | MIM number | Mendelian diseases |
| Reactome | Pathway ID | Pathway data |
| GO | GO ID | Gene ontology |

---

## Data Quality Notes

1. **Association Coverage:** 14M+ target-disease associations
2. **Target Coverage:** 62,000+ genes
3. **Disease Coverage:** 20,000+ EFO terms
4. **Evidence Strings:** 17M+ evidence records
5. **Data Sources:** 20+ integrated sources
6. **GraphQL API:** Recommended access method
7. **Parquet Downloads:** Bulk data access
8. **Quarterly Updates:** Regular data releases

