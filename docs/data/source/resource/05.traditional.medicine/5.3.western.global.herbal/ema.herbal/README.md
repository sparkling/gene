---
id: ema.herbal
title: "EMA Herbal Medicines - European Medicines Agency"
type: source
parent: ../README.md
tier: 1
status: active
category: traditional.medicine
subcategory: western.global.herbal
tags:
  - herbal-medicine
  - european
  - regulatory
  - monographs
  - traditional-use
---

# EMA Herbal Medicines - European Medicines Agency

**Category:** [Traditional Medicine](../../README.md) > [Western & Global Herbal](../README.md)

## Overview

The European Medicines Agency (EMA) maintains an authoritative database of herbal monographs through its Committee on Herbal Medicinal Products (HMPC). These monographs provide evidence-based assessments of herbal substances, defining conditions for well-established use and traditional use registrations in the European Union.

EMA herbal monographs represent the gold standard for regulatory-grade herbal medicine documentation in Europe. Each monograph includes comprehensive information on the herbal substance, pharmaceutical preparations, therapeutic indications, posology, contraindications, and safety data derived from systematic literature review.

The database distinguishes between "well-established use" (supported by clinical evidence) and "traditional use" (based on long-standing use history), providing clear regulatory pathways for herbal products.

## Key Statistics

| Metric | Value |
|--------|-------|
| Herbal Monographs | 200+ |
| Community Lists | 80+ |
| Assessment Reports | 200+ |
| Herbal Substances | 250+ |
| Traditional Use Registrations | 150+ |
| Well-Established Use | 50+ |

## Primary Use Cases

1. Regulatory documentation for herbal product registration
2. Evidence-based safety and efficacy assessment
3. Standardizing herbal medicine quality requirements
4. Defining acceptable therapeutic indications
5. Supporting pharmacovigilance for herbal products

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Monograph ID | EMA reference | EMA/HMPC/xxxxxx |
| Herbal Substance | Latin binomial | Valeriana officinalis L. |
| Plant Part | Standardized term | radix (root) |
| Preparation | Pharmaceutical form | Dry extract, tincture |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.ema.europa.eu/en/human-regulatory-overview/herbal-medicinal-products | Official portal |
| Monograph Search | EMA website | PDF downloads |
| Community Lists | Published lists | Herbal substances |

## Data Model

```
Herbal Substance (Latin name + plant part)
            |
            v
    EMA Monograph
    ├── Well-Established Use
    │   └── Clinical evidence
    └── Traditional Use
        └── Long-standing use
            |
            v
    ├── Indications
    ├── Posology
    ├── Contraindications
    ├── Interactions
    └── Safety Profile
```

## Monograph Structure

| Section | Content |
|---------|---------|
| Introduction | Botanical identity, quality standards |
| Pharmaceutical | Preparations, extracts, standardization |
| Clinical | Indications, dosage, duration |
| Safety | Contraindications, warnings, interactions |
| Pharmacology | Mechanism (when known) |
| Non-clinical | Toxicology data |

## Regulatory Categories

| Category | Evidence Level | Example |
|----------|----------------|---------|
| Well-Established Use | Clinical trials, systematic reviews | Ginkgo biloba for cognitive function |
| Traditional Use | >=30 years of use (>=15 in EU) | Valerian for sleep disorders |
| Community List | Harmonized across EU | Entry substances |

## Quality Standards

| Standard | Description |
|----------|-------------|
| European Pharmacopoeia | Reference quality standards |
| Herbal Drug | Whole or processed plant material |
| Herbal Preparation | Extract, tincture, powder |
| Marker Compounds | Standardization references |

## License

| Aspect | Value |
|--------|-------|
| License | Public domain (regulatory documents) |
| Commercial Use | Reference freely available |
| Attribution | Citation recommended |

## Key Strengths

- **Regulatory Authority**: Official EU assessment
- **Evidence-Based**: Systematic literature review
- **Safety Focus**: Comprehensive risk evaluation
- **Quality Standards**: Pharmacopoeia-linked

## Limitations

- Focused on European market
- Limited to assessed substances
- Traditional use != efficacy proof
- Updates follow regulatory timelines

## See Also

- [NAPRALERT](../napralert/README.md)
- [Dr. Duke's Phytochemical Database](../../02.compounds.molecules/2.1.natural.products/dr.dukes/README.md)
