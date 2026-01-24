# Category 06: Nutrition & Food - Data Dictionary

## Overview

This data dictionary documents all fields for Category 06: Nutrition & Food, which encompasses food composition databases, dietary supplement information, bioactive compound analysis, and metabolomics data.

| Attribute | Value |
|-----------|-------|
| Category ID | 06 |
| Category Name | Nutrition & Food |
| Total Data Sources | 8 |
| Extraction Date | 2026-01-24 |

## Subcategories

| ID | Name | Data Sources |
|----|------|--------------|
| 6.1 | Food Composition | Open Food Facts, FooDB |
| 6.2 | Dietary Supplements | ConsumerLab, DSLD, Natural Medicines |
| 6.3 | Bioactive Food Compounds | eBASIS |
| 6.4 | Metabolomics | Exposome-Explorer, HMDB |

---

## Cross-Category Common Fields

These fields appear across multiple subcategories and data sources, enabling cross-referencing and integration.

### Chemical Identifiers

| Field Name | Data Type | Description | Sources |
|------------|-----------|-------------|---------|
| cas_number | string | Chemical Abstracts Service registry number | FooDB, eBASIS, Natural Medicines, Exposome-Explorer, HMDB |
| pubchem_cid | integer | PubChem compound identifier | eBASIS, Exposome-Explorer, HMDB |
| inchikey | string | Hashed International Chemical Identifier | FooDB, eBASIS, HMDB |
| smiles | string | SMILES molecular structure notation | FooDB, HMDB |
| molecular_formula | string | Chemical elemental formula | FooDB, eBASIS, Exposome-Explorer, HMDB |
| molecular_weight | number | Molecular mass in g/mol | FooDB, eBASIS, HMDB |

### Food Identifiers

| Field Name | Data Type | Description | Sources |
|------------|-----------|-------------|---------|
| food_name | string | Common name of food item | Open Food Facts, FooDB, eBASIS |
| scientific_name | string | Taxonomic scientific name | FooDB, eBASIS, Natural Medicines |
| food_group | string | Category classification of food | Open Food Facts, FooDB, eBASIS, Exposome-Explorer |

### Concentration Fields

| Field Name | Data Type | Description | Sources |
|------------|-----------|-------------|---------|
| value/concentration | number | Measured concentration value | FooDB (orig_content), eBASIS (value), Exposome-Explorer (concentration_value), HMDB (concentration_value) |
| unit | string | Unit of measurement | FooDB (orig_unit), eBASIS (unit), DSLD (unit), Exposome-Explorer (concentration_units), HMDB (concentration_units) |
| min/max values | number | Range of measured values | FooDB, eBASIS, Exposome-Explorer, HMDB |
| n_samples/subject_count | integer | Number of samples/subjects analyzed | eBASIS, Exposome-Explorer, HMDB |

### Reference Fields

| Field Name | Data Type | Description | Sources |
|------------|-----------|-------------|---------|
| pubmed_id | integer | PubMed article identifier | eBASIS, Natural Medicines, Exposome-Explorer, HMDB |
| doi | string | Digital Object Identifier | eBASIS, Exposome-Explorer, HMDB |

### Supplement Identifiers

| Field Name | Data Type | Description | Sources |
|------------|-----------|-------------|---------|
| upc | string | Universal Product Code barcode | ConsumerLab, DSLD |
| brand_name | string | Commercial brand name | Open Food Facts, ConsumerLab, DSLD |
| product_name | string | Product name | Open Food Facts, ConsumerLab, DSLD |

---

## Key Cross-References

The following identifiers enable linking data across databases:

1. **HMDB ID** - Links metabolomics to food composition databases
2. **FooDB ID** - Provides cross-reference between HMDB and food compounds
3. **PubChem CID** - Enables linking across all chemical databases
4. **CAS Number** - Universal chemical identifier across all sources
5. **InChIKey** - Enables structural matching across databases

---

## Summary Statistics

| Subcategory | Unified Fields | Source-Specific Fields |
|-------------|---------------|------------------------|
| 6.1 Food Composition | 7 | 12 |
| 6.2 Dietary Supplements | 10 | 34 |
| 6.3 Bioactive Food Compounds | 10 | 18 |
| 6.4 Metabolomics | 13 | 40 |

---

## Related Documents

- [6.1 Food Composition Data Dictionary](./6.1.food.composition/_data-dictionary.md)
- [6.2 Dietary Supplements Data Dictionary](./6.2.dietary.supplements/_data-dictionary.md)
- [6.3 Bioactive Food Compounds Data Dictionary](./6.3.bioactive.food.compounds/_data-dictionary.md)
- [6.4 Metabolomics Data Dictionary](./6.4.metabolomics/_data-dictionary.md)
