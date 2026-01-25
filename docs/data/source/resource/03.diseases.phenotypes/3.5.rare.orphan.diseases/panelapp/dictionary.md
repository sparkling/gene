# PanelApp - Data Dictionary

## Overview

This data dictionary documents the schema for Genomics England PanelApp.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | panelapp |
| **Name** | PanelApp |
| **Parent** | 3.5.rare.orphan.diseases |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Gene Panel

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | integer | 1:1 | Yes | Panel identifier | 245 |
| name | string | 1:1 | Yes | Panel name | Intellectual disability |
| version | string | 1:1 | Yes | Panel version | 4.12 |
| relevant_disorders | array | 1:N | No | Associated disorders | [Intellectual disability] |
| signed_off | boolean | 1:1 | No | Signed off for clinical use | true |
| status | string | 1:1 | Yes | Panel status | public |
| types | array | 1:N | No | Panel types | [Rare Disease] |
| stats | object | 1:1 | No | Gene statistics | {green: 150, amber: 20} |

### Gene Entry

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_symbol | string | 1:1 | Yes | HGNC gene symbol | MECP2 |
| hgnc_id | string | 1:1 | Yes | HGNC identifier | HGNC:6990 |
| gene_name | string | 1:1 | Yes | Full gene name | methyl-CpG binding protein 2 |
| confidence_level | integer | 1:1 | Yes | Evidence level (0-3) | 3 |
| mode_of_inheritance | string | 1:1 | No | Inheritance pattern | X-LINKED: hemizygous |
| penetrance | string | 1:1 | No | Penetrance | Complete |
| publications | array | 1:N | No | PMID references | [10508514, 10521577] |
| phenotypes | array | 1:N | No | Associated phenotypes | [Rett syndrome] |
| evidence | array | 1:N | No | Evidence sources | [Expert Review] |

### Review

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| review_id | integer | 1:1 | Yes | Review identifier | 12345 |
| gene_symbol | string | 1:1 | Yes | Gene reviewed | MECP2 |
| reviewer | string | 1:1 | Yes | Reviewer name | Expert Curator |
| rating | integer | 1:1 | Yes | Evidence rating (1-4) | 4 |
| comment | string | 1:1 | No | Review comments | Well-established... |
| date | datetime | 1:1 | Yes | Review date | 2024-01-15 |

### STR (Short Tandem Repeat)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| name | string | 1:1 | Yes | STR name | FMR1 CGG repeat |
| gene_symbol | string | 1:1 | Yes | Associated gene | FMR1 |
| chromosome | string | 1:1 | Yes | Chromosome | X |
| position | object | 1:1 | Yes | Genomic position | {start: 147912110} |
| normal_range | string | 1:1 | No | Normal repeat range | 5-44 |
| pathogenic_range | string | 1:1 | No | Pathogenic range | >200 |
| confidence_level | integer | 1:1 | Yes | Evidence level | 3 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Panel ID | Integer | 245 | Panel identifier |
| HGNC ID | HGNC:[0-9]+ | HGNC:6990 | Gene identifier |
| Gene Symbol | Text | MECP2 | HGNC symbol |
| PMID | Integer | 10508514 | PubMed reference |
| Ensembl | ENSG[0-9]{11} | ENSG00000169057 | Gene annotation |
| OMIM | OMIM:[0-9]{6} | OMIM:300005 | Disease reference |

---

## Enumerations

### Evidence Levels

| Level | Color | Description | Clinical Use |
|-------|-------|-------------|--------------|
| 3 | Green | High evidence, diagnostic grade | Yes |
| 2 | Amber | Moderate evidence, emerging | With caution |
| 1 | Red | Low evidence | No |
| 0 | Gray | No evidence yet | No |

### Panel Status

| Status | Description |
|--------|-------------|
| public | Publicly visible |
| internal | Internal use only |
| superpanel | Collection of panels |
| retired | No longer maintained |

### Panel Types

| Type | Description |
|------|-------------|
| Rare Disease | Rare disease panels |
| Cancer | Cancer gene panels |
| Pharmacogenomics | Drug response panels |
| GMS | NHS Genomic Medicine Service |
| Research | Research-only panels |

### Mode of Inheritance

| Mode | Description |
|------|-------------|
| MONOALLELIC | Autosomal dominant |
| BIALLELIC | Autosomal recessive |
| BOTH | Both mono and biallelic |
| X-LINKED: hemizygous | X-linked recessive |
| X-LINKED: monoallelic | X-linked dominant |
| MITOCHONDRIAL | Mitochondrial inheritance |
| Unknown | Not determined |

### Penetrance

| Value | Description |
|-------|-------------|
| Complete | Full penetrance |
| Incomplete | Variable penetrance |
| Unknown | Not determined |

### Evidence Sources

| Source | Description |
|--------|-------------|
| Expert Review | Expert panel review |
| Literature | Published evidence |
| NHS GMS | NHS testing evidence |
| Radboud | Radboud panels |
| UKGTN | UK Genetic Testing Network |
| Illumina | Illumina panels |
| OMIM | OMIM curation |

### Gene Rating Scale

| Rating | Description |
|--------|-------------|
| 1 | Refuted evidence |
| 2 | Limited evidence |
| 3 | Moderate evidence |
| 4 | Strong evidence |

---

## Entity Relationships

### Panel to Genes
- **Cardinality:** 1:N
- **Description:** Panels contain multiple genes
- **Key Fields:** panel_id, gene_symbol

### Gene to Reviews
- **Cardinality:** 1:N
- **Description:** Genes have multiple reviews
- **Key Fields:** gene_symbol, review_id

### Panel to Disorders
- **Cardinality:** N:M
- **Description:** Panels relate to disorders
- **Key Fields:** panel_id, relevant_disorders

### Gene to Publications
- **Cardinality:** 1:N
- **Description:** Literature evidence
- **Key Fields:** gene_symbol, publications

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NHS | National Health Service | UK healthcare |
| GMS | Genomic Medicine Service | NHS program |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| STR | Short Tandem Repeat | Repeat expansion |
| UKGTN | UK Genetic Testing Network | Testing network |
| MOI | Mode of Inheritance | Inheritance pattern |
| VUS | Variant of Uncertain Significance | Variant class |
| OMIM | Online Mendelian Inheritance in Man | Disease database |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| HGNC | HGNC ID | Gene nomenclature |
| Ensembl | Gene ID | Gene annotation |
| OMIM | MIM number | Disease associations |
| PubMed | PMID | Literature evidence |
| HPO | HP ID | Phenotype terms |
| Orphanet | ORPHA ID | Rare diseases |

---

## Data Quality Notes

1. **Gene Panels:** 350+ curated panels
2. **Green Genes:** 4,500+ diagnostic-grade genes
3. **Expert Reviews:** 50,000+ community reviews
4. **NHS Integration:** Used in NHS GMS
5. **International Use:** 20+ countries
6. **REST API:** Free, no registration
7. **Version Control:** Full version history
8. **Weekly Updates:** Regular panel reviews

