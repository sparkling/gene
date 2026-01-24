---
id: schema-mbodymap
title: "mBodyMap Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: final
tags: [schema, database, microbiome, body-sites, atlas, visualization]
---

# mBodyMap Schema Documentation

**Document ID:** SCHEMA-MBODYMAP
**Version:** 2023.01
**Source Version:** mBodyMap v2.0

---

## TL;DR

mBodyMap is a comprehensive atlas of human microbiome data organized by body site, integrating 500+ studies and 100,000+ samples across 50+ anatomical locations. Provides standardized taxonomic profiles, diversity metrics, and visualization tools for comparing microbiomes across body habitats. Essential for understanding microbiome biogeography.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Body sites | 50+ | mBodyMap Stats |
| Samples | 100,000+ | mBodyMap Stats |
| Studies integrated | 500+ | mBodyMap Stats |
| Species cataloged | 5,000+ | mBodyMap Stats |
| Data points | 10M+ | mBodyMap Stats |

---

## Entity Relationship Overview

```
mBodyMap Data
  ├── Body Site
  │     ├── Site code (hierarchical)
  │     ├── Anatomical location
  │     └── Site category
  ├── Sample
  │     ├── Sample ID
  │     ├── Study source
  │     └── Metadata
  ├── Taxonomic Profile
  │     ├── Species abundances
  │     ├── Alpha diversity
  │     └── Beta diversity
  └── Reference Data
        ├── Site-specific taxa
        ├── Core microbiome
        └── Diversity benchmarks
```

---

## Core Tables/Entities

### Body Site

**Description:** Anatomical location classification.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| site_code | string | Yes | Hierarchical site code |
| site_name | string | Yes | Common name |
| site_category | string | Yes | Major body region |
| anatomy_term | string | No | Anatomical ontology term |
| sample_count | integer | Yes | Number of samples |
| study_count | integer | No | Number of studies |

### Sample Entry

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| sample_id | string | Yes | Unique sample identifier |
| site_code | string | Yes | Body site code |
| study_id | string | Yes | Source study |
| subject_id | string | No | Subject identifier |
| age | integer | No | Subject age |
| sex | string | No | male/female |
| health_status | string | No | Healthy/Disease |
| sequencing_method | string | Yes | 16S/WGS |

### Taxonomic Profile

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| sample_id | string | Yes | Sample identifier |
| taxon_id | integer | Yes | NCBI Taxonomy ID |
| taxon_name | string | Yes | Scientific name |
| rank | string | Yes | Taxonomic rank |
| relative_abundance | float | Yes | Proportion (0-1) |
| absolute_count | integer | No | Read count |

### Diversity Metrics

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| sample_id | string | Yes | Sample identifier |
| shannon | float | No | Shannon diversity |
| simpson | float | No | Simpson diversity |
| chao1 | float | No | Chao1 richness |
| observed_species | integer | No | Species count |
| evenness | float | No | Pielou evenness |

---

## Body Site Hierarchy

### Site Code Format

```
{Region}_{Subregion}_{Specific_Site}

Examples:
GI_LI_COLON     - Gastrointestinal > Large Intestine > Colon
ORAL_BUCCAL     - Oral > Buccal mucosa
SKIN_ARM_FOSSA  - Skin > Arm > Antecubital fossa
URO_VAG         - Urogenital > Vagina
RESP_NASAL      - Respiratory > Nasal cavity
```

### Major Body Regions

| Code | Region | Subsites |
|------|--------|----------|
| GI | Gastrointestinal | Stomach, Small intestine, Large intestine, Rectum |
| ORAL | Oral cavity | Tongue, Buccal, Gingiva, Palate, Saliva |
| SKIN | Skin | Various body locations |
| URO | Urogenital | Vagina, Urethra, Bladder |
| RESP | Respiratory | Nasal, Pharynx, Lung |
| OTHER | Other | Blood, Eye, Ear |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | TSV (abundance matrices) |
| Alternative | BIOM, JSON |
| Visualization | Interactive web |
| Encoding | UTF-8 |

---

## Sample Record

### Body Site JSON

```json
{
  "site_code": "GI_LI_COLON",
  "site_name": "Colon",
  "site_category": "Gastrointestinal",
  "anatomy_term": "UBERON:0001155",
  "parent_site": "GI_LI",
  "sample_count": 45000,
  "study_count": 150,
  "dominant_phyla": ["Firmicutes", "Bacteroidetes", "Proteobacteria"],
  "typical_shannon": 3.5,
  "typical_richness": 200
}
```

### Sample Profile

```json
{
  "sample_id": "MBMAP_S0001234",
  "site_code": "GI_LI_COLON",
  "study_id": "PRJNA123456",
  "subject_id": "SUB001",
  "age": 45,
  "sex": "female",
  "health_status": "healthy",
  "sequencing_method": "16S_V4",
  "taxonomic_profile": [
    {"taxon_name": "Bacteroides vulgatus", "taxon_id": 821, "abundance": 0.15},
    {"taxon_name": "Faecalibacterium prausnitzii", "taxon_id": 853, "abundance": 0.12},
    {"taxon_name": "Prevotella copri", "taxon_id": 165179, "abundance": 0.08}
  ],
  "diversity": {
    "shannon": 3.8,
    "simpson": 0.92,
    "observed_species": 185
  }
}
```

### Abundance Matrix (TSV)

```
sample_id	Bacteroides	Prevotella	Faecalibacterium	Ruminococcus
MBMAP_S0001234	0.15	0.08	0.12	0.05
MBMAP_S0001235	0.22	0.02	0.18	0.08
MBMAP_S0001236	0.08	0.25	0.05	0.03
```

---

## Site-Specific Reference Profiles

### Gut Microbiome Reference

| Taxon | Typical Abundance | Variation |
|-------|-------------------|-----------|
| Bacteroides | 10-30% | High individual |
| Firmicutes (total) | 40-60% | Moderate |
| Faecalibacterium | 5-15% | Moderate |
| Prevotella | 0-30% | Diet-dependent |

### Oral Microbiome Reference

| Taxon | Typical Abundance | Site Preference |
|-------|-------------------|-----------------|
| Streptococcus | 20-40% | Throughout oral cavity |
| Veillonella | 10-20% | Tongue, saliva |
| Haemophilus | 5-15% | Pharynx |
| Neisseria | 5-10% | Tongue, mucosa |

### Skin Microbiome Reference

| Taxon | Typical Abundance | Site Preference |
|-------|-------------------|-----------------|
| Cutibacterium | 20-50% | Sebaceous sites |
| Staphylococcus | 10-30% | Moist sites |
| Corynebacterium | 5-20% | Moist sites |
| Malassezia | Variable | Sebaceous sites |

---

## Cross-Site Comparisons

| Metric | Gut | Oral | Skin | Vaginal |
|--------|-----|------|------|---------|
| Diversity (Shannon) | 3-4 | 2-3 | 1-3 | 1-2 |
| Dominant phylum | Firmicutes | Firmicutes | Actinobacteria | Firmicutes |
| Species richness | 100-300 | 50-200 | 50-150 | 20-100 |

---

## Cross-References

| Database | ID Type | Usage |
|----------|---------|-------|
| NCBI Taxonomy | TaxID | Species identification |
| UBERON | Anatomy terms | Body site ontology |
| HMP | Sample IDs | Reference data |
| GMrepo | Project IDs | Gut samples |
| HOMD | HOMT IDs | Oral taxa |

---

## Glossary

| Term | Definition |
|------|------------|
| Body site | Anatomical location for sampling |
| Alpha diversity | Within-sample diversity |
| Beta diversity | Between-sample diversity |
| Core microbiome | Taxa present in most samples |
| Biogeography | Spatial distribution of microbes |
| Relative abundance | Proportion of total community |

---

## References

1. mBodyMap Website: https://mbodymap.microbiome.cloud
2. Human Microbiome Project Consortium. (2012). Structure, function and diversity of the healthy human microbiome. Nature.
3. Lloyd-Price J, et al. (2017). Strains, functions and dynamics in the expanded Human Microbiome Project. Nature.
