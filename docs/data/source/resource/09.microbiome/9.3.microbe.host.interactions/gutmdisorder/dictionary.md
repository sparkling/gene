# gutMDisorder - Data Dictionary

## Overview

This data dictionary documents gutMDisorder microbe-disease association entries.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gutmdisorder |
| **Name** | gutMDisorder |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| association_id | integer | Yes | Unique association identifier |
| microbe.taxid | integer | Yes | NCBI Taxonomy ID |
| microbe.name | string | Yes | Scientific name |
| disease.name | string | Yes | Disease name |
| association.direction | string | Yes | increased, decreased, altered |
| evidence.pmid | integer | Yes | PubMed literature reference |

---

## Direction Values

| Value | Description | Interpretation |
|-------|-------------|----------------|
| increased | Higher in disease | Potential pathogenic role |
| decreased | Lower in disease | Potential protective role |
| altered | Changed (inconsistent) | Context-dependent |

---

## Disease Categories

| Category | Examples |
|----------|----------|
| Metabolic | Diabetes, Obesity, NAFLD |
| Gastrointestinal | IBD, IBS, Colorectal cancer |
| Neurological | Parkinson's, Alzheimer's, Autism |
| Cardiovascular | Atherosclerosis, Hypertension |
| Autoimmune | Rheumatoid arthritis, MS |
| Psychiatric | Depression, Anxiety |

---

## Study Design Types

| Design | Evidence Level |
|--------|----------------|
| Meta-analysis | Highest |
| RCT | Highest |
| Cohort | High |
| Longitudinal | High |
| Case-control | Moderate |
| Cross-sectional | Lower |

---

## Taxonomic Distribution

| Phylum | Typical Pattern |
|--------|-----------------|
| Firmicutes | Often decreased in disease |
| Bacteroidetes | Variable by disease |
| Proteobacteria | Often increased (dysbiosis marker) |
| Actinobacteria | Variable |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| IBD | Inflammatory Bowel Disease |
| IBS | Irritable Bowel Syndrome |
| NAFLD | Non-Alcoholic Fatty Liver Disease |
| MS | Multiple Sclerosis |
| RCT | Randomized Controlled Trial |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
