# TTD - Data Dictionary

## Overview

This data dictionary documents the schema for TTD (Therapeutic Target Database).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | ttd |
| **Name** | TTD |
| **Parent** | 2.7.compound.target.interactions |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Target Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ttd_target_id | string | 1:1 | Yes | TTD identifier | TTDT00001 |
| name | string | 1:1 | Yes | Target name | Epidermal growth factor receptor |
| uniprot_id | string | 1:1 | No | UniProt accession | P00533 |
| gene_name | string | 1:1 | No | Gene symbol | EGFR |
| target_type | string | 1:1 | Yes | Target molecule type | Protein |
| biochemical_class | string | 1:1 | No | Biochemical classification | Kinase |
| ecl_classification | string | 1:1 | No | EC number if enzyme | 2.7.10.1 |
| target_validation | string | 1:1 | No | Validation status | Successful target |
| organism | string | 1:1 | No | Species | Homo sapiens |

### Drug Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ttd_drug_id | string | 1:1 | Yes | TTD drug ID | D0A9YA |
| name | string | 1:1 | Yes | Drug name | Erlotinib |
| cas_number | string | 1:1 | No | CAS registry | 183321-74-6 |
| pubchem_cid | integer | 1:1 | No | PubChem compound | 176870 |
| chembl_id | string | 1:1 | No | ChEMBL ID | CHEMBL553 |
| drugbank_id | string | 1:1 | No | DrugBank ID | DB00530 |
| drug_class | string | 1:1 | No | Drug modality | Small Molecule |
| highest_status | string | 1:1 | No | Development status | Approved |
| company | string | 1:1 | No | Developer company | Genentech |

### Target-Drug Relationship

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ttd_target_id | string | 1:1 | Yes | Foreign key to targets | TTDT00001 |
| ttd_drug_id | string | 1:1 | Yes | Foreign key to drugs | D0A9YA |
| activity_type | string | 1:1 | No | Interaction type | Inhibitor |
| mechanism | string | 1:1 | No | Mechanism of action | ATP-competitive |
| development_status | string | 1:1 | Yes | Clinical status | Approved |
| indication | string | 1:1 | No | Therapeutic indication | NSCLC |

### Disease Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| disease_id | string | 1:1 | Yes | TTD disease ID | D0001 |
| name | string | 1:1 | Yes | Disease name | Non-small cell lung cancer |
| icd_11_code | string | 1:1 | No | ICD-11 classification | 2C25 |
| therapeutic_area | string | 1:1 | No | Therapeutic area | Oncology |

### Pathway Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| pathway_id | string | 1:1 | Yes | Pathway identifier | hsa05200 |
| pathway_name | string | 1:1 | Yes | Pathway name | Pathways in cancer |
| source | string | 1:1 | Yes | Pathway source | KEGG |
| external_id | string | 1:1 | No | Source pathway ID | hsa05200 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| TTD Target ID | TTDT + 5 digits | TTDT00001 | Target identifier |
| TTD Drug ID | 6 characters | D0A9YA | Drug identifier |
| UniProt ID | Accession | P00533 | Protein identifier |
| ChEMBL ID | CHEMBL + digits | CHEMBL553 | Chemical database |
| DrugBank ID | DB + 5 digits | DB00530 | Drug database |
| PubChem CID | Integer | 176870 | Compound identifier |
| CAS Number | varies | 183321-74-6 | Chemical registry |
| ICD-11 Code | varies | 2C25 | Disease classification |
| KEGG ID | hsa + digits | hsa05200 | Pathway identifier |

---

## Enumerations

### Target Validation Status

| Status | Count | Description |
|--------|-------|-------------|
| Successful target | 426 | Has approved drug |
| Clinical trial target | 1,014 | Phase I-III drugs |
| Preclinical/Patented | 212 | Early development |
| Literature-reported | 1,479 | Research targets |

### Drug Development Status

| Status | Description |
|--------|-------------|
| Approved | FDA/EMA approved |
| Phase III | Late-stage trials |
| Phase II | Efficacy trials |
| Phase I | Safety trials |
| Preclinical | Pre-IND |
| Discontinued | Development stopped |
| Withdrawn | Removed from market |

### Target Types

| Type | Description |
|------|-------------|
| Protein | Protein target |
| DNA | DNA target |
| RNA | RNA target |

### Biochemical Classes

| Class | Description |
|-------|-------------|
| Kinase | Protein kinase |
| Receptor | Membrane receptor |
| Enzyme | Non-kinase enzyme |
| Ion channel | Ion channel |
| Transporter | Membrane transporter |
| Nuclear receptor | NHR |
| Transcription factor | TF |

### Drug Classes

| Type | Description | Examples |
|------|-------------|----------|
| Small Molecule | Chemical compounds | Imatinib |
| Antibody | Monoclonal antibodies | Trastuzumab |
| ADC | Antibody-drug conjugate | T-DM1 |
| Peptide | Therapeutic peptides | Exenatide |
| ASO | Antisense oligonucleotide | Nusinersen |
| siRNA | Small interfering RNA | Patisiran |
| mRNA | Messenger RNA | COVID vaccines |
| Cell Therapy | CAR-T, stem cells | Axicabtagene |
| Gene Therapy | Gene delivery | Zolgensma |

### Activity Types

| Type | Description |
|------|-------------|
| Inhibitor | Inhibits target |
| Agonist | Activates target |
| Antagonist | Blocks target |
| Modulator | Modulates target |
| Antibody | Binds target (mAb) |

### Therapeutic Areas

| Area | Examples |
|------|----------|
| Oncology | Cancer indications |
| CNS | Neurological diseases |
| Cardiovascular | Heart disease |
| Infectious | Viral, bacterial |
| Metabolic | Diabetes, obesity |
| Immunology | Autoimmune diseases |

---

## Entity Relationships

### Target to Drugs
- **Cardinality:** N:M
- **Description:** Targets have multiple drugs; drugs target multiple proteins
- **Key Fields:** ttd_target_id, ttd_drug_id

### Target to Diseases
- **Cardinality:** N:M
- **Description:** Targets associated with multiple diseases
- **Key Fields:** ttd_target_id, disease_id

### Target to Pathways
- **Cardinality:** N:M
- **Description:** Targets participate in multiple pathways
- **Key Fields:** ttd_target_id, pathway_id

### Drug to Diseases
- **Cardinality:** N:M
- **Description:** Drugs treat multiple diseases
- **Key Fields:** ttd_drug_id, disease_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| TTD | Therapeutic Target Database | Database name |
| ICD-11 | International Classification of Diseases 11 | Disease coding |
| NSCLC | Non-Small Cell Lung Cancer | Example indication |
| IND | Investigational New Drug | Regulatory milestone |
| mAb | Monoclonal Antibody | Drug modality |
| ADC | Antibody-Drug Conjugate | Drug modality |
| ASO | Antisense Oligonucleotide | Drug modality |
| siRNA | Small Interfering RNA | Drug modality |
| CAR-T | Chimeric Antigen Receptor T-cell | Therapy type |
| FDA | Food and Drug Administration | US regulator |
| EMA | European Medicines Agency | EU regulator |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| UniProt | Accession | Target proteins |
| ChEMBL | ChEMBL ID | Bioactivity data |
| DrugBank | DrugBank ID | Drug information |
| PubChem | CID | Chemical data |
| KEGG | Pathway ID | Pathway context |
| Reactome | Pathway ID | Reactions |
| WikiPathways | WP ID | Community pathways |

---

## Data Quality Notes

1. **Target Coverage:** 3,131 total targets
2. **Successful Targets:** 426 with approved drugs
3. **Drug Coverage:** 39,862 total drugs
4. **Approved Drugs:** 2,895 marketed drugs
5. **Disease Coding:** ICD-11 standardization
6. **Pathway Integration:** KEGG, Reactome, WikiPathways
7. **Development Tracking:** Phase I through approval
8. **Annual Updates:** Regular data releases

