---
id: traditional.chinese.medicine
title: "Traditional Chinese Medicine"
type: subcategory
parent: ../README.md
last_updated: 2026-01-23
status: active
tags: [tcm, chinese-medicine, herbs, formulas, acupuncture]
---

# Traditional Chinese Medicine

**Parent:** [Traditional Medicine](../README.md)

## Overview

Traditional Chinese Medicine (TCM) databases document the herbs, formulas, and therapeutic principles of this ancient medical system. These resources catalog thousands of medicinal plants, classical formulations, and their constituent bioactive compounds.

Key resources include HERB (high-throughput experiment-based database), ETCM (Encyclopedia of TCM), TCMBANK, TCMSID, BATMAN-TCM (target prediction), and SymMap (symptom mapping). Together they provide comprehensive TCM data integration.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [HERB](./herb/README.md) | 1 | High-throughput TCM database |
| [ETCM](./etcm/README.md) | 1 | Encyclopedia of TCM |
| [TCMBANK](./tcmbank/README.md) | 2 | TCM compound and target database |
| [TCMSID](./tcmsid/README.md) | 2 | TCM ingredient database |
| [BATMAN-TCM](./batman.tcm/README.md) | 2 | TCM target prediction |
| [SymMap](./symmap/README.md) | 2 | Symptom-herb-compound mapping |

## Integration Notes

HERB provides the most comprehensive integration of TCM ingredients with modern targets. ETCM offers detailed herb and formula annotations. BATMAN-TCM predicts targets for TCM compounds. SymMap links TCM symptoms to modern phenotypes. Cross-reference compounds with PubChem for standardized structures.

---

## Database Selection Guide

### TL;DR

Traditional Chinese Medicine has the most extensive collection of specialized databases among all intervention categories, with 30+ databases covering herbs, formulas, compounds, targets, diseases, symptoms, and gene expression data. Priority integration targets BATMAN-TCM 2.0 (API available, 2.3M predicted interactions), HERB 2.0 (transcriptomics), and TCMSID (ADME properties, CC BY 4.0 license), with TCMBank and SymMap 2.0 as secondary priorities for comprehensive coverage.

### Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Primary API source | BATMAN-TCM 2.0 | Only TCM database with REST API + bulk download |
| Gene expression source | HERB 2.0 | Unique transcriptomic data (2,231 experiments) |
| Commercial-friendly source | TCMSID (CC BY 4.0) | ADME properties + open license |
| Symptom mapping | SymMap 2.0 | Bridges TCM-modern medicine phenotypes |
| Formula coverage | ETCM v2.0 + BATMAN-TCM 2.0 | 48K + 55K formulas respectively |
| Clinical evidence | CMAUP 2024 + HERB 2.0 | 691 + 8,558 clinical trials |

### Comparative Statistics

| Database | Herbs | Ingredients | Targets | Formulas | License |
|----------|-------|-------------|---------|----------|---------|
| BATMAN-TCM 2.0 | 8,404 | 39,171 | 2.3M pairs | 54,832 | CC BY-NC 4.0 |
| TCMBank | 9,192 | 61,966 | 15,179 | - | CC BY 4.0 |
| HERB 2.0 | 7,263 | 49,258 | 12,933 | 9 | Academic |
| TCM-Mesh | 6,235 | 383,840 | 14,298 | - | Academic |
| ETCM v2.0 | 2,079 | 38,298 | 1,040 | 48,442 | Academic |
| TCMSID | 499 | 20,015 | 3,270 | - | CC BY 4.0 |
| SymMap 2.0 | 499 | 19,595 | 4,302 | - | Academic |
| CMAUP 2024 | 7,865 | 60,222 | 758 | - | Academic |

### Unique Features by Database

| Feature | Best Database(s) |
|---------|-----------------|
| API Access | BATMAN-TCM 2.0 (REST), TCMNP (R package) |
| ADME Properties | TCMSID (14 ADME/T), YaTCM (50 properties) |
| Gene Expression | HERB 2.0 (2,231 experiments) |
| Target Prediction | BATMAN-TCM 2.0 (2.3M predictions, AUC=0.97) |
| Clinical Trials | CMAUP 2024 (691 trials), HERB 2.0 (8,558 trials) |
| Symptom Mapping | SymMap 2.0 (1,717 TCM symptoms mapped) |
| Blood Constituents | DCABM-TCM (actual bioavailable compounds) |
| Toxicity | TCM-Mesh (163,221 side effect records) |

### Tier Priorities

| Tier | Databases | Priority Reason |
|------|-----------|-----------------|
| **1 (Critical)** | BATMAN-TCM 2.0, HERB 2.0, TCMSID, TCMBank, CMAUP 2024 | API access; gene expression; commercial license; clinical trials |
| **2 (High Value)** | SymMap 2.0, DCABM-TCM, HIT 2.0, TM-MC 2.0, ETCM v2.0 | Symptom mapping; bioavailability; validated interactions |
| **3 (Enrichment)** | TCM-Mesh, TCMID 2.0, TCMNP, SuperTCM, TCMPG 2.0 | Toxicity data; genomics; programmatic tools |

### Licensing Summary

| License Type | Databases | Commercial Use |
|--------------|-----------|----------------|
| CC BY 4.0 | TCMBank, TCMSID, TCM Database@Taiwan, TCMPG 2.0 | Yes (with attribution) |
| CC BY-NC 4.0 | BATMAN-TCM 2.0, DCABM-TCM | Contact for commercial |
| CC BY | TCMM | Yes |
| Academic Use | HERB, SymMap, ETCM, YaTCM, HIT, CMAUP | Research only |
