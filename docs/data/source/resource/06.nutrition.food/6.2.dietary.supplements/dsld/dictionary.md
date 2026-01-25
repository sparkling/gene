# DSLD - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | dsld |
| **Name** | Dietary Supplement Label Database |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| dsld_id | integer | 1:1 | Yes | Internal product ID | `123456` |
| product_name | string | 1:1 | Yes | Full product name | `Vitamin D3 5000 IU` |
| brand_name | string | 1:1 | No | Brand name | `Nature's Way` |
| upc | string | 1:1 | No | Universal Product Code | `033674100233` |

### Product Details

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| net_contents | string | 1:1 | No | Package quantity | `120 Softgels` |
| serving_size | string | 1:1 | No | Serving description | `1 softgel` |
| product_type | string | 1:1 | No | Product form | `Softgel` |

### Ingredient Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ingredient_name | string | 1:N | Yes | Ingredient as labeled | `Vitamin D3` |
| amount | string | 1:1 | No | Amount per serving | `5000` |
| unit | string | 1:1 | No | Unit of measurement | `IU` |
| daily_value_percent | string | 1:1 | No | Percent Daily Value | `1250` |
| is_blend_component | boolean | 1:1 | No | Part of proprietary blend | `false` |

---

## Enumerations

### Product Types

| Value | Description |
|-------|-------------|
| Capsule | Hard or soft capsule |
| Tablet | Compressed tablet |
| Softgel | Soft gelatin capsule |
| Liquid | Liquid form |
| Powder | Powder form |
| Gummy | Gummy/chewable |
| Lozenge | Dissolving lozenge |
| Spray | Oral spray |

### Ingredient Groups

| Value | Examples |
|-------|----------|
| Vitamins | Vitamin A, B12, C, D, E, K |
| Minerals | Calcium, Iron, Magnesium |
| Herbs/Botanicals | Echinacea, Ginkgo, Turmeric |
| Amino Acids | L-Arginine, L-Glutamine |
| Probiotics | Lactobacillus, Bifidobacterium |
| Specialty | CoQ10, Glucosamine, Melatonin |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DSLD | Dietary Supplement Label Database | NIH database |
| ODS | Office of Dietary Supplements | NIH office |
| NIH | National Institutes of Health | US agency |
| NHP | Natural Health Product | Canadian identifier |
| UPC | Universal Product Code | Barcode |
| IU | International Unit | Vitamin measurement |
| DV | Daily Value | Nutrient reference |

---

## Data Quality Notes

1. **Public domain data**: US government resource
2. **Label accuracy**: Reflects labels, not verified content
3. **Proprietary blends**: Individual amounts not disclosed
4. **API access**: Free with registration

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
