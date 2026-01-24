---
id: hit.2.0
title: "HIT 2.0 - Herbal Ingredients' Targets"
type: data-source
category: traditional-medicine
subcategory: multi-system-integration
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [herbal-medicine, validated-targets, multi-system, experimental, drug-targets]
---

# HIT 2.0 - Herbal Ingredients' Targets

**Category:** [Traditional Medicine](../../_index.md) > [Multi-System Integration](../_index.md)

## Overview

HIT 2.0 (Herbal Ingredients' Targets) is a comprehensive database of experimentally validated interactions between herbal ingredients and their molecular targets. Unlike prediction-based databases, HIT 2.0 focuses exclusively on interactions supported by experimental evidence from literature, making it the gold standard for validated herbal compound-target relationships.

The database integrates data from multiple traditional medicine systems including Traditional Chinese Medicine, Ayurveda, and Western herbal medicine, providing a unified resource for validated targets. Each interaction includes the experimental method, binding affinity (when available), and literature citations.

HIT 2.0 serves as a critical validation resource for network pharmacology studies, helping researchers distinguish between predicted and confirmed interactions.

## Key Statistics

| Metric | Value |
|--------|-------|
| Herbal Ingredients | 10,000+ |
| Validated Targets | 5,000+ |
| Herb-Target Interactions | 100,000+ |
| Literature References | 50,000+ |
| Traditional Systems | Multiple |
| Binding Affinity Data | 30,000+ |

## Primary Use Cases

1. Validating predicted compound-target interactions
2. Finding experimentally confirmed targets for herbal compounds
3. Cross-referencing multiple TCM databases
4. Building validated interaction networks
5. Identifying binding affinity data for natural products

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Ingredient ID | HIT internal | HIT001234 |
| Target ID | UniProt / Gene Symbol | P12345 / TP53 |
| Interaction ID | HIT internal | HIT_INT_001234 |
| PubChem CID | Integer | 969516 |
| PubMed ID | PMID | 12345678 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://hit2.badd-cao.net/ | Interactive search |
| Download | Bulk data available | TSV/Excel files |
| Search | Multi-field search | Ingredient, target, herb |

## Data Model

```
Herbal Sources (Multiple Systems)
├── TCM Herbs
├── Ayurvedic Plants
└── Western Herbs
            |
            v
   Ingredients (10,000+)
            |
            v (experimental validation)
   Target Proteins (5,000+)
            |
            +-- Binding Affinity (30,000+)
            +-- Literature (50,000+)
            +-- Experimental Method
```

## Validation Evidence Types

| Evidence Type | Description | Count |
|---------------|-------------|-------|
| Binding Assay | Direct binding measurement | High |
| Enzyme Inhibition | IC50, Ki values | High |
| Cell-Based Assay | Functional validation | Medium |
| In Vivo | Animal model confirmation | Medium |
| Clinical | Human studies | Limited |

## Binding Affinity Data

| Metric | Description |
|--------|-------------|
| IC50 | Half-maximal inhibitory concentration |
| Ki | Inhibition constant |
| Kd | Dissociation constant |
| EC50 | Half-maximal effective concentration |
| Activity Type | Inhibitor, agonist, antagonist |

## Multi-System Coverage

| System | Coverage |
|--------|----------|
| Traditional Chinese Medicine | Comprehensive |
| Ayurveda | Major herbs |
| Japanese Kampo | Selected |
| Western Herbal | Common herbs |
| African Traditional | Limited |

## Integration Value

HIT 2.0 serves as validation for other databases:

| Database | Integration Purpose |
|----------|---------------------|
| BATMAN-TCM | Confirm predicted TTIs |
| KampoDB | Validate docking predictions |
| IMPPAT | Cross-reference STITCH predictions |
| SymMap | Ground symptom-target links |

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

## Key Strengths

- **Experimental Validation**: No predictions, only confirmed
- **Multi-System**: Spans multiple traditional medicine systems
- **Quantitative**: Binding affinity data included
- **Literature-Linked**: Full citation support

## Limitations

- Literature curation lag
- Not all interactions have affinity data
- Coverage biased toward studied compounds
- English literature bias

## See Also

- [BATMAN-TCM](../../5.1.traditional.chinese.medicine/batman.tcm/_index.md)
- [KampoDB](../../5.2.south.east.asian.systems/kampodb/_index.md)
- [IMPPAT](../../5.2.south.east.asian.systems/imppat/_index.md)
- [SymMap](../../5.1.traditional.chinese.medicine/symmap/_index.md)
