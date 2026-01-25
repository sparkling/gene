---
id: symmap
title: "SymMap 2.0"
type: source
parent: ../README.md
tier: 1
status: active
category: traditional.medicine
subcategory: traditional.chinese.medicine
tags:
  - tcm
  - traditional-chinese-medicine
  - symptom-mapping
  - phenotype
  - network-medicine
---

# SymMap 2.0

**Category:** [Traditional Medicine](../../_index.md) > [Traditional Chinese Medicine](../_index.md)

## Overview

SymMap 2.0 is a unique integrative database that bridges Traditional Chinese Medicine symptoms to modern medical phenotypes and molecular mechanisms. It maps the relationship between TCM symptoms (Zheng), herbs, targets, and diseases, enabling systematic analysis of how traditional symptom descriptions relate to contemporary disease understanding.

The database addresses a critical challenge in TCM research: translating traditional diagnostic categories and symptom descriptions into modern biomedical terminology. By connecting TCM symptoms to Human Phenotype Ontology (HPO) terms and disease associations, SymMap enables computational analysis of traditional medicine from a phenotype-driven perspective.

SymMap 2.0 significantly expanded the symptom-target associations and improved the integration with molecular interaction networks.

## Key Statistics

| Metric | Value |
|--------|-------|
| TCM Symptoms | 1,717 |
| Modern Phenotypes (HPO) | 5,000+ |
| TCM Herbs | 961 |
| Compounds | 19,595 |
| Target Genes | 4,302 |
| Diseases | 5,235 |
| Symptom-Herb Associations | 28,212 |
| Symptom-Gene Associations | 322,073 |

## Primary Use Cases

1. Mapping TCM symptoms to modern disease phenotypes
2. Understanding molecular basis of TCM symptom patterns
3. Identifying herbs for specific symptom profiles
4. Network analysis of symptom-gene-disease relationships
5. Phenotype-driven drug repurposing from TCM

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Symptom ID | SMSY internal | SMSY00001 |
| Herb ID | SymMap internal | SMHB00001 |
| Compound ID | SymMap / PubChem | SMCP00001 / 5280343 |
| Target ID | Gene Symbol | TP53 |
| Disease ID | MESH / OMIM | D001234 / 123456 |
| HPO ID | HP:XXXXXXX | HP:0001250 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://www.symmap.org/ | Interactive browsing |
| Download | Bulk data available | TSV files |
| Network View | Cytoscape integration | Visual analysis |

## Data Model

```
TCM Symptoms (1,717) <----> Modern Phenotypes (HPO)
        |                           |
        v                           v
    Herbs (961)              Diseases (5,235)
        |                           |
        v                           v
  Compounds (19,595) --------> Targets (4,302)
```

## Symptom Mapping Features

| Feature | Description |
|---------|-------------|
| TCM Zheng | Traditional symptom patterns |
| HPO Mapping | Modern phenotype ontology links |
| MESH Disease | Disease concept associations |
| Pathway Links | KEGG pathway annotations |

## Unique Value Proposition

SymMap is the only major database focusing on symptom-level integration:

| Integration | Description |
|-------------|-------------|
| TCM Theory | Traditional symptom classification |
| Modern Phenotyping | HPO standardization |
| Molecular | Gene and pathway connections |
| Therapeutic | Herb recommendations per symptom |

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

## Limitations

- Symptom-phenotype mappings may be approximate
- Not all TCM concepts map cleanly to HPO
- Manual curation limits coverage expansion

## See Also

- [BATMAN-TCM](../batman.tcm/_index.md)
- [HERB Database](../herb/_index.md)
- [ETCM](../etcm/_index.md)
