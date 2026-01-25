# HOMD - Data Dictionary

## Overview

This data dictionary documents Human Oral Microbiome Database taxon records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | homd |
| **Name** | Human Oral Microbiome Database |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| homt_id | string | Yes | HOMD taxon ID (HMT-###) |
| ncbi_taxid | integer | No | NCBI Taxonomy ID |
| genus | string | Yes | Genus name |
| species | string | Yes | Species epithet or provisional |
| full_name | string | No | Complete scientific name |
| cultivation_status | string | Yes | Named, Unnamed, Provisional |

---

## HOMT ID Format

```
Format: HMT-###
Where ### is a 3-digit number (001-999)

Examples:
HMT-096: Streptococcus mutans
HMT-707: Porphyromonas gingivalis
HMT-274: Aggregatibacter actinomycetemcomitans
```

---

## Cultivation Status

| Status | Description | Naming Convention |
|--------|-------------|-------------------|
| Named | Validly published name | Full binomial (e.g., Streptococcus mutans) |
| Unnamed | Cultivated but not named | Genus sp. HMT-### |
| Provisional | Uncultivated, sequence-only | Genus [G-#] sp. HMT-### |

---

## Disease Associations

| Disease | Key Taxa (HOMT IDs) |
|---------|---------------------|
| Dental caries | HMT-096, HMT-686, HMT-416 |
| Periodontitis | HMT-707, HMT-274, HMT-690 |
| Gingivitis | HMT-582, HMT-626 |
| Halitosis | HMT-707, HMT-879 |

---

## Oral Site Categories

| Category | Specific Sites |
|----------|----------------|
| Teeth | Enamel, Dentin, Root |
| Gingiva | Supragingival, Subgingival |
| Mucosal | Buccal, Tongue, Palate |
| Saliva | Whole saliva |
| Tonsils | Palatine, Lingual |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| HOMD | Human Oral Microbiome Database |
| HOMT | Human Oral Microbiome Taxon |
| eHOMD | expanded Human Oral Microbiome Database |
| 16S | 16S Ribosomal RNA |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
