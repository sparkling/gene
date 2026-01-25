# gutMGene - Data Dictionary

## Overview

This data dictionary documents gutMGene microbe-host gene expression associations.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gutmgene |
| **Name** | gutMGene |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| association_id | integer | Yes | Unique association identifier |
| microbe.taxid | integer | Yes | NCBI Taxonomy ID |
| microbe.name | string | Yes | Scientific species name |
| gene.gene_id | integer | Yes | NCBI Entrez Gene ID |
| gene.symbol | string | Yes | Official gene symbol |
| effect.direction | string | Yes | up, down, altered |
| evidence.pmid | integer | Yes | PubMed literature reference |

---

## Effect Direction Values

| Value | Description | Interpretation |
|-------|-------------|----------------|
| up | Increased expression | Microbe activates gene |
| down | Decreased expression | Microbe suppresses gene |
| altered | Changed expression | Direction varies by context |

---

## Experiment Types

| Type | Description |
|------|-------------|
| Germ-free colonization | Mono-association of GF animals |
| Antibiotic treatment | Microbiome depletion |
| Probiotic supplementation | Addition of specific strains |
| FMT | Fecal microbiota transplant |
| Co-culture | In vitro cell-microbe culture |
| Correlation | Observational association |

---

## Model Systems

| Model | Description |
|-------|-------------|
| Germ-free mice | Axenic animals for causation studies |
| SPF mice | Specific pathogen-free baseline |
| Cell lines | Caco-2, HT-29 intestinal cells |
| Organoids | 3D intestinal models |
| Human cohorts | Observational studies |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| GF | Germ-free |
| SPF | Specific Pathogen-Free |
| FMT | Fecal Microbiota Transplant |
| qPCR | Quantitative Polymerase Chain Reaction |
| RNA-seq | RNA Sequencing |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
