# GMrepo - Data Dictionary

## Overview

This data dictionary documents GMrepo gut microbiome sample records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gmrepo |
| **Name** | GMrepo |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| run_id | string | Yes | Sequencing run accession (SRR/ERR) |
| project_id | string | Yes | BioProject accession (PRJNA/PRJEB) |
| sample_id | string | No | BioSample accession |
| phenotype | string | No | Disease/health phenotype |
| phenotype_mesh_id | string | No | MeSH identifier for phenotype |
| sequencing_type | string | Yes | 16S or WGS |
| relative_abundance | number | Yes | Proportion of taxon (0-1) |

---

## Phenotype Categories

| Category | Example Conditions |
|----------|-------------------|
| Metabolic Disorders | Diabetes, Obesity, NAFLD |
| Gastrointestinal | IBD, IBS, Colorectal cancer |
| Neurological | Parkinson's, Alzheimer's |
| Cardiovascular | Atherosclerosis, Hypertension |
| Autoimmune | RA, Multiple sclerosis |
| Healthy Controls | Health baseline samples |

---

## Sequencing Types

| Type | Description | Sample Count (v3) |
|------|-------------|-------------------|
| 16S | 16S rRNA amplicon | 87,048 |
| WGS | Whole-genome shotgun | 31,917 |

---

## Marker Taxa Fields

| Field | Description |
|-------|-------------|
| direction | increased or decreased in disease |
| effect_size | Magnitude of association |
| p_value | Statistical significance |
| study_count | Number of supporting studies |
| sample_count | Total samples analyzed |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| GMrepo | Gut Microbiome Repository |
| MeSH | Medical Subject Headings |
| SRA | Sequence Read Archive |
| ENA | European Nucleotide Archive |
| WGS | Whole-Genome Shotgun |
| 16S | 16S Ribosomal RNA |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
