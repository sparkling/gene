---
id: downloads-traditional-medicine
title: "Traditional Medicine Database Bulk Downloads"
type: download-guide
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [downloads, bulk-data, traditional-medicine, tcm, ayurveda]
---

**Parent:** [Download Guides](./_index.md)

# Traditional Medicine Database Bulk Downloads

This document provides comprehensive information on bulk download methods for Traditional Chinese Medicine (TCM), Ayurveda, and related natural products databases.

---

## Table of Contents

1. [BATMAN-TCM 2.0](#1-batman-tcm-20)
2. [TCMBank](#2-tcmbank)
3. [HERB 2.0](#3-herb-20)
4. [IMPPAT 2.0](#4-imppat-20)
5. [KampoDB](#5-kampodb)
6. [Dr. Duke's Phytochemical Database](#6-dr-dukes-phytochemical-database)
7. [LOTUS](#7-lotus)
8. [KNApSAcK](#8-knapsack)
9. [Processing Recommendations](#9-processing-recommendations)

---

[Full content from traditional-medicine.md continues...]

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **BATMAN-TCM 2.0** | API | `http://bionet.ncpsb.org.cn/batman-tcm/` |
| **TCMBank** | Web | `https://tcmbank.cn/` |
| **HERB 2.0** | Download | `http://herb.ac.cn/` |
| **IMPPAT 2.0** | Web | `https://cb.imsc.res.in/imppat/` |
| **KampoDB** | Web | `http://wakanmoview.inm.u-toyama.ac.jp/kampo/` |
| **Dr. Duke's** | CSV | `https://phytochem.nal.usda.gov/` |
| **LOTUS** | SPARQL | `https://lotus.naturalproducts.net/` |
| **KNApSAcK** | Web | `http://www.knapsackfamily.com/` |

**Access Requirements:** Most are freely accessible for academic use; Dr. Duke's is CC0.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | CSV, TSV, JSON |
| Alternative | SDF, XML |
| Chemical structures | SMILES, InChI |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `compound_id` | string | Compound identifier | "IMPPAT001234" |
| `herb_name` | string | Plant source | "Withania somnifera" |
| `traditional_use` | string | Traditional indication | "Rasayana" |
| `targets` | array | Predicted targets | ["ESR1", "AR"] |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `found_in` | Herb | N:M |
| `targets` | Protein | N:M |

## Sample Data

### Example Compound Record
```json
{
  "compound": "Withaferin A",
  "smiles": "CC1(C)CCC2C...",
  "herb": "Ashwagandha",
  "traditional_system": "Ayurveda",
  "predicted_targets": ["HSP90AA1", "NFE2L2"]
}
```

### Sample Query Result
| compound | herb | system | targets |
|----------|------|--------|---------|
| Withaferin A | Ashwagandha | Ayurveda | 45 |
| Ginsenoside Rg1 | Ginseng | TCM | 78 |

## License

| Source | License | Commercial Use |
|--------|---------|----------------|
| Dr. Duke's | CC0 | Yes |
| LOTUS | CC0 | Yes |
| BATMAN-TCM | Academic | Research only |
| KampoDB | CC BY-SA 4.0 | Yes |
| IMPPAT | Academic | Research only |

## Data Set Size

| Metric | Value |
|--------|-------|
| IMPPAT compounds | 17K+ phytochemicals |
| BATMAN-TCM herbs | 700+ herbs |
| LOTUS compounds | 750K+ structure-organism pairs |
| Dr. Duke's records | 60K+ plant-compound links |
| Total storage estimate | ~5-10 GB combined |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Bulk Download | Retrieving entire datasets via FTP/HTTP rather than API queries | Downloading full BATMAN-TCM database |
| Phytochemical | Bioactive chemical compound produced by plants | Curcumin, ginsenosides |
| Traditional Medicine | Healthcare practices developed before modern medicine | TCM, Ayurveda, Kampo |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| BATMAN-TCM | Bioinformatics Analysis Tool for Molecular mechANism of TCM | TCM target prediction |
| TCMBank | Database of TCM herbal compounds and targets | TCM compounds |
| HERB 2.0 | High-throughput Experiment- and Reference-guided database of TCM | Gene expression |
| IMPPAT | Indian Medicinal Plants, Phytochemistry And Therapeutics | Ayurveda database |
| KampoDB | Database of Kampo medicines and ingredients | Japanese herbalism |
| LOTUS | Natural Products Online database | Natural products |
| KNApSAcK | Comprehensive Species-Metabolite Relationship database | Metabolite data |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST endpoints |
| CC BY | Creative Commons Attribution | Open license |
| CC BY-NC | Creative Commons Attribution Non-Commercial | Restricted license |
| FTP | File Transfer Protocol | Download method |
| HERB | High-throughput Experiment- and Reference-guided database | TCM database |
| TCM | Traditional Chinese Medicine | Herbal system |

---

*Document created: 2026-01-18*
*Last updated: 2026-01-18*
