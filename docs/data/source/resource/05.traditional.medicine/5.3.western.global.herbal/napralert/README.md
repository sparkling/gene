---
id: napralert
title: "NAPRALERT - Natural Products Alert"
type: source
parent: ../README.md
tier: 1
status: active
category: traditional.medicine
subcategory: western.global.herbal
tags:
  - natural-products
  - ethnobotany
  - pharmacognosy
  - literature
  - bioactivity
---

# NAPRALERT - Natural Products Alert

**Category:** [Traditional Medicine](../../README.md) > [Western & Global Herbal](../README.md)

## Overview

NAPRALERT (Natural Products Alert) is the world's largest relational database on natural products, maintained by the University of Illinois at Chicago since 1975. It contains comprehensive literature data on the ethnomedical, pharmacological, and biochemical properties of natural products from plant, microbial, and animal sources.

The database represents over 50 years of systematic literature curation, covering ethnobotanical uses, biological activities, chemical constituents, and pharmacological studies. NAPRALERT is particularly valuable for its historical depth and global coverage of traditional medicine literature.

Unlike compound-centric databases, NAPRALERT excels at connecting traditional uses with scientific studies, enabling researchers to trace the evidence base for ethnobotanical claims.

## Key Statistics

| Metric | Value |
|--------|-------|
| Literature Records | 200,000+ |
| Plant Species | 60,000+ |
| Natural Products | 180,000+ |
| Ethnobotanical Records | 100,000+ |
| Pharmacological Records | 150,000+ |
| Years of Coverage | 1975-present |

## Primary Use Cases

1. Literature review for natural product research
2. Ethnobotanical documentation and validation
3. Identifying bioactive compounds from traditional sources
4. Pharmacognosy research support
5. Prior art searches for natural product patents

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Organism ID | NAPRALERT internal | NAP_ORG_001234 |
| Compound ID | NAPRALERT internal | NAP_CPD_001234 |
| Literature ID | PubMed/internal | PMID:12345678 |
| CAS Number | CAS format | 458-37-7 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://napralert.pharmacy.uic.edu/ | Subscription required |
| Subscription | UIC licensing | Academic/commercial tiers |
| Custom Searches | By request | Comprehensive reports |

## Data Model

```
Literature Sources (200,000+)
            |
            v
    +-------+-------+
    |               |
    v               v
Organisms      Natural Products
(60,000+)       (180,000+)
    |               |
    v               v
+---+---+      +----+----+
|       |      |         |
v       v      v         v
Ethnobotany  Bioactivity  Chemistry
(100,000+)   (150,000+)   (180,000+)
```

## Data Categories

| Category | Description |
|----------|-------------|
| Ethnobotany | Traditional uses, preparation methods |
| Pharmacology | In vitro/vivo activity, clinical studies |
| Biochemistry | Compound isolation, structure elucidation |
| Toxicology | Safety data, adverse effects |
| Chemistry | Structural data, synthesis |

## Organism Coverage

| Kingdom | Examples |
|---------|----------|
| Plantae | Medicinal plants worldwide |
| Fungi | Medicinal mushrooms |
| Bacteria | Antibiotic producers |
| Animalia | Marine organisms, insects |

## Ethnobotanical Data

| Field | Description |
|-------|-------------|
| Traditional Use | Documented therapeutic application |
| Geographic Origin | Region/country of use |
| Preparation | Method of preparation |
| Administration | Route and dosage |
| Source | Literature reference |

## License

| Aspect | Value |
|--------|-------|
| License | Subscription-based |
| Academic Use | Institutional subscription |
| Commercial Use | Commercial license required |
| Individual | Limited free searches |

## Key Strengths

- **Historical Depth**: 50+ years of literature
- **Global Coverage**: Worldwide ethnobotany
- **Literature Integration**: Primary source citations
- **Expert Curation**: Pharmacognosy specialists

## Limitations

- Subscription required for full access
- Interface may be dated
- Focus on literature (not primary data)
- Updates lag behind publication

## Subscription Tiers

| Tier | Access Level |
|------|--------------|
| Free | Limited queries |
| Academic | Institutional unlimited |
| Commercial | Full access + reports |

## See Also

- [EMA Herbal](../ema.herbal/README.md)
- [Dr. Duke's](../../02.compounds.molecules/2.1.natural.products/dr.dukes/README.md)
- [IMPPAT](../../5.2.south.east.asian.systems/imppat/README.md)
