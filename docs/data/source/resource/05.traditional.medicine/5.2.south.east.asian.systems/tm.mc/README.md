---
id: tm.mc
title: "TM-MC - Traditional Medicine - Molecular Correlates"
type: source
parent: ../README.md
tier: 2
status: active
category: traditional.medicine
subcategory: south.east.asian.systems
tags:
  - traditional-medicine
  - multi-system
  - molecular-targets
  - ayurveda
  - tcm
---

# TM-MC - Traditional Medicine - Molecular Correlates

**Category:** [Traditional Medicine](../../_index.md) > [South East Asian Systems](../_index.md)

## Overview

TM-MC (Traditional Medicine - Molecular Correlates) is an integrative database that connects traditional medicine knowledge from multiple Asian systems to modern molecular understanding. It bridges Ayurveda, Traditional Chinese Medicine, and other regional systems through shared compounds, targets, and therapeutic indications.

The database emphasizes the molecular basis of traditional therapeutic claims, correlating ethnobotanical usage patterns with protein targets and pathway mechanisms. TM-MC is particularly valuable for comparative analysis across different traditional medicine systems.

By identifying compounds and targets shared between multiple traditional medicine systems, TM-MC supports hypothesis generation for drug discovery based on convergent traditional knowledge.

## Key Statistics

| Metric | Value |
|--------|-------|
| Medicinal Plants | 2,500+ |
| Traditional Systems | 5+ |
| Compounds | 8,000+ |
| Molecular Targets | 1,500+ |
| Therapeutic Indications | 500+ |
| Plant-Compound Links | 25,000+ |

## Primary Use Cases

1. Cross-system comparison of traditional medicine
2. Identifying plants with shared therapeutic uses across cultures
3. Molecular validation of traditional indications
4. Discovering convergent ethnobotanical knowledge
5. Target-based analysis of traditional remedies

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Plant ID | TM-MC internal | TMMC_P_001234 |
| Compound ID | TM-MC internal | TMMC_C_001234 |
| System Code | Abbreviation | AYU, TCM, UNA |
| PubChem CID | Integer | 5280343 |
| UniProt ID | Accession | P12345 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://tm-mc.org/ | Interactive search |
| Download | Via web interface | CSV/Excel |
| Browse | By system, plant, compound | Structured navigation |

## Data Model

```
Traditional Systems (Ayurveda, TCM, Unani, Siddha, etc.)
                    |
                    v
          Medicinal Plants (2,500+)
                    |
         +----------+----------+
         |                     |
         v                     v
  Compounds (8,000+)    Therapeutic Uses (500+)
         |
         v
  Molecular Targets (1,500+)
         |
         v
  Pathways & Mechanisms
```

## Traditional Systems Covered

| System | Region | Coverage |
|--------|--------|----------|
| Ayurveda | India | Comprehensive |
| TCM | China | Major herbs |
| Unani | Middle East/India | Selected |
| Siddha | South India | Selected |
| Kampo | Japan | Selected |
| Jamu | Indonesia | Limited |

## Cross-System Features

| Feature | Description |
|---------|-------------|
| Shared Plants | Species used across multiple systems |
| Common Compounds | Chemicals linking different traditions |
| Convergent Uses | Similar indications across cultures |
| Molecular Overlap | Shared target proteins |

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

## Unique Value

TM-MC provides unique cross-cultural perspective:

| Comparison | Insight |
|------------|---------|
| Ayurveda vs TCM | Identify overlapping herbs and uses |
| Multi-system convergence | Validate traditional claims |
| Molecular mechanism | Link tradition to modern biology |

## Limitations

- Coverage varies by traditional system
- Not all plants have molecular data
- Therapeutic term standardization challenges
- Updates depend on literature curation

## See Also

- [IMPPAT](../imppat/_index.md)
- [BATMAN-TCM](../../5.1.traditional.chinese.medicine/batman.tcm/_index.md)
- [KampoDB](../kampodb/_index.md)
