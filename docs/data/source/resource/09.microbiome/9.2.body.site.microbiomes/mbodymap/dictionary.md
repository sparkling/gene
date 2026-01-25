# mBodyMap - Data Dictionary

## Overview

This data dictionary documents mBodyMap body site microbiome atlas records.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | mbodymap |
| **Name** | mBodyMap |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| sample_id | string | Yes | mBodyMap sample identifier |
| site_code | string | Yes | Hierarchical body site code |
| site_category | string | Yes | Major body region code |
| site_name | string | No | Body site common name |
| sequencing_method | string | No | 16S variant or WGS |
| health_status | string | No | healthy, disease, unknown |

---

## Site Code Format

```
{Region}_{Subregion}_{Specific_Site}

Examples:
GI_LI_COLON     - Gastrointestinal > Large Intestine > Colon
ORAL_BUCCAL     - Oral > Buccal mucosa
SKIN_ARM_FOSSA  - Skin > Arm > Antecubital fossa
URO_VAG         - Urogenital > Vagina
RESP_NASAL      - Respiratory > Nasal cavity
```

---

## Site Categories

| Code | Region | Subsites |
|------|--------|----------|
| GI | Gastrointestinal | Stomach, Small intestine, Large intestine |
| ORAL | Oral cavity | Tongue, Buccal, Gingiva, Palate |
| SKIN | Skin | Various body locations |
| URO | Urogenital | Vagina, Urethra |
| RESP | Respiratory | Nasal, Pharynx, Lung |
| OTHER | Other | Blood, Eye, Ear |

---

## Diversity Metrics

| Metric | Description |
|--------|-------------|
| shannon | Shannon diversity index |
| simpson | Simpson diversity index |
| chao1 | Chao1 richness estimator |
| observed_species | Count of observed species |
| evenness | Pielou evenness index |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| GI | Gastrointestinal |
| URO | Urogenital |
| RESP | Respiratory |
| WGS | Whole-Genome Shotgun |
| 16S | 16S Ribosomal RNA |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
