# Exposome-Explorer - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | exposome-explorer |
| **Name** | Exposome-Explorer |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| association_id | integer | 1:1 | Yes | Record identifier | `456` |
| biomarker_id | integer | 1:1 | Yes | Biomarker ID | `123` |
| exposure_id | integer | 1:1 | Yes | Exposure ID | `789` |
| biomarker_name | string | 1:1 | Yes | Biomarker name | `4-O-methylgallic acid` |
| exposure_name | string | 1:1 | Yes | Exposure name | `Tea consumption` |

### Cross-References

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| hmdb_id | string | 1:1 | No | HMDB identifier | `HMDB0002122` |
| pubchem_cid | integer | 1:1 | No | PubChem compound ID | `13159` |
| chebi_id | string | 1:1 | No | ChEBI identifier | `CHEBI:16243` |

### Statistical Association

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| effect_size | number | 1:1 | No | Association magnitude | `0.42` |
| effect_unit | string | 1:1 | No | Effect measure type | `Pearson correlation` |
| p_value | number | 1:1 | No | P-value | `0.00001` |
| sample_size | integer | 1:1 | No | Study size | `450` |

---

## Enumerations

### Exposure Categories

| Value | Examples |
|-------|----------|
| Dietary | Foods, beverages, nutrients |
| Environmental | Pollutants, metals, pesticides |
| Lifestyle | Smoking, alcohol, physical activity |
| Occupational | Industrial chemicals |

### Specimen Types

| Value | Description |
|-------|-------------|
| plasma | Blood plasma |
| serum | Blood serum |
| urine | Spot or 24h urine |
| whole_blood | Complete blood |
| erythrocytes | Red blood cells |
| adipose | Fat tissue |
| hair | Head hair |
| nails | Fingernails/toenails |

### Study Types

| Value | Description |
|-------|-------------|
| cross-sectional | Single time point |
| cohort | Longitudinal follow-up |
| intervention | Controlled feeding study |
| case-control | Cases vs controls |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| IARC | International Agency for Research on Cancer | Database host |
| WHO | World Health Organization | Parent organization |
| HMDB | Human Metabolome Database | Cross-reference |
| PAH | Polycyclic Aromatic Hydrocarbons | Pollutants |
| PCB | Polychlorinated Biphenyls | Pollutants |
| OR | Odds Ratio | Effect measure |
| CI | Confidence Interval | Statistical uncertainty |

---

## Data Quality Notes

1. **Manually curated**: All associations from peer-reviewed literature
2. **Systematic collection**: Comprehensive literature review
3. **Effect validation**: Multiple studies per association when available
4. **Open access**: CC BY 4.0 license

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
