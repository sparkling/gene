# DGIdb - Data Dictionary

## Overview

This data dictionary documents the schema for DGIdb (Drug Gene Interaction Database).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | dgidb |
| **Name** | DGIdb |
| **Parent** | 2.7.compound.target.interactions |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Gene Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| concept_id | string | 1:1 | Yes | Normalized gene ID | hgnc:3236 |
| name | string | 1:1 | Yes | Gene symbol (HGNC) | EGFR |
| long_name | string | 1:1 | No | Full gene name | Epidermal growth factor receptor |
| gene_categories | array | 1:N | No | Druggability categories | [KINASE, DRUGGABLE GENOME] |
| entrez_id | integer | 1:1 | No | NCBI Entrez Gene ID | 1956 |
| ensembl_id | string | 1:1 | No | Ensembl gene ID | ENSG00000146648 |

### Drug Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| concept_id | string | 1:1 | Yes | Normalized drug ID | chembl:CHEMBL939 |
| name | string | 1:1 | Yes | Drug name | Gefitinib |
| approved | boolean | 1:1 | No | FDA approval status | true |
| immunotherapy | boolean | 1:1 | No | Immunotherapy drug | false |
| antineoplastic | boolean | 1:1 | No | Anticancer drug | true |
| chembl_id | string | 1:1 | No | ChEMBL identifier | CHEMBL939 |
| drugbank_id | string | 1:1 | No | DrugBank identifier | DB00317 |
| aliases | array | 1:N | No | Alternative names | [Iressa, ZD1839] |

### Interaction Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | string | 1:1 | Yes | Interaction identifier | int123 |
| gene | object | 1:1 | Yes | Gene reference | {name: EGFR} |
| drug | object | 1:1 | Yes | Drug reference | {name: Gefitinib} |
| interaction_types | array | 1:N | No | Interaction type labels | [inhibitor] |
| interaction_score | decimal | 1:1 | No | Confidence score | 5.2 |
| sources | array | 1:N | Yes | Source databases | [ChEMBL, DrugBank] |
| pmids | array | 1:N | No | PubMed references | [12345678, 23456789] |

### Gene Category Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| name | string | 1:1 | Yes | Category name | KINASE |
| source_db_name | string | 1:1 | Yes | Source of category | dGene |

### Source Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| source_db_name | string | 1:1 | Yes | Database name | ChEMBL |
| source_db_version | string | 1:1 | No | Version | 33 |
| citation | string | 1:1 | No | Reference citation | - |
| base_url | string | 1:1 | No | Database URL | https://www.ebi.ac.uk/chembl |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Gene Concept ID | hgnc: + ID | hgnc:3236 | Normalized gene ID |
| Drug Concept ID | chembl: + ID | chembl:CHEMBL939 | Normalized drug ID |
| HGNC ID | Integer | 3236 | HUGO Gene nomenclature |
| ChEMBL ID | CHEMBL + digits | CHEMBL939 | Chemical database |
| DrugBank ID | DB + 5 digits | DB00317 | Drug database |
| Entrez ID | Integer | 1956 | NCBI gene ID |
| Ensembl ID | ENSG + digits | ENSG00000146648 | Ensembl gene ID |

---

## Enumerations

### Interaction Types

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
| substrate | Drug is substrate |
| inducer | Drug induces expression |
| suppressor | Drug suppresses activity |

### Gene Categories

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
| PHOSPHATASE | Phosphatase genes |
| DNA REPAIR | DNA repair genes |
| CELL SURFACE | Cell surface proteins |

### Directionality

| Direction | Description |
|-----------|-------------|
| inhibits | Decreases target activity |
| activates | Increases target activity |
| unknown | Direction not specified |

### Source Databases

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

## Entity Relationships

### Drug to Interactions
- **Cardinality:** 1:N
- **Description:** One drug can have interactions with multiple genes
- **Key Fields:** drug.concept_id, interaction.id

### Gene to Interactions
- **Cardinality:** 1:N
- **Description:** One gene can have interactions with multiple drugs
- **Key Fields:** gene.concept_id, interaction.id

### Gene to Categories
- **Cardinality:** N:M
- **Description:** Genes belong to multiple druggability categories
- **Key Fields:** gene.concept_id, gene_categories

### Interaction to Sources
- **Cardinality:** N:M
- **Description:** Interactions supported by multiple databases
- **Key Fields:** interaction.id, sources

### Interaction to Publications
- **Cardinality:** 1:N
- **Description:** Interactions linked to literature
- **Key Fields:** interaction.id, pmids

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DGIdb | Drug Gene Interaction Database | Primary resource |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| GraphQL | Graph Query Language | API format |
| PGx | Pharmacogenomics | Drug-gene relationships |
| GPCR | G Protein-Coupled Receptor | Gene category |
| TF | Transcription Factor | Gene category |
| DDI | Drug-Drug Interaction | Related concept |
| CAR-T | Chimeric Antigen Receptor T-cell | Therapy type |
| mAb | Monoclonal Antibody | Drug type |
| EGFR | Epidermal Growth Factor Receptor | Example target |
| BRAF | B-Raf Proto-Oncogene | Example target |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| ChEMBL | ChEMBL ID | Bioactivity data |
| DrugBank | DrugBank ID | Drug information |
| PharmGKB | Accession | Pharmacogenomics |
| TTD | TTD ID | Therapeutic targets |
| CIViC | Evidence ID | Clinical evidence |
| OncoKB | Alteration | Oncology |
| HGNC | HGNC ID | Gene nomenclature |
| Ensembl | Ensembl ID | Gene annotation |

---

## Data Quality Notes

1. **Aggregated Data:** 40+ source databases integrated
2. **Interaction Coverage:** 90,000+ drug-gene interactions
3. **Gene Coverage:** 13,000+ genes with druggability info
4. **Drug Coverage:** 10,000+ drugs and compounds
5. **GraphQL API:** Flexible query interface
6. **Gene Categories:** 40+ druggability classifications
7. **Interaction Types:** 30+ type annotations
8. **Open Access:** Free for research use

