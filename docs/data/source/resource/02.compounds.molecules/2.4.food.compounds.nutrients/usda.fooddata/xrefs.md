---
id: xrefs-usda-fooddata
title: "USDA FoodData Central Cross-References"
type: xrefs
parent: _index.md
last_updated: 2026-01-23
---

# USDA FoodData Central Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `02.compounds.molecules/2.4.food.compounds.nutrients/usda.fooddata` |
| Secondary | Symlink | `06.nutrition.food/6.1.food.composition/usda.fooddata` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| NDB Number | ndb_number | High (legacy) |
| FDC ID | fdc_id | High (current) |
| GTIN/UPC | gtin_upc | High (branded) |
| Nutrient ID | nutrient_id | High |
| Food Group | food_group_id | High |
| INFOODS | infoods_tagname | High |

## Integration Notes

USDA FoodData Central is the authoritative source for food composition data in the United States, bridging nutrient chemistry with food science applications.

**Primary use (Food Compounds):** Nutrient content per food, analytical chemistry data, food composition standards.

**Secondary use (Food Composition):** Dietary assessment, recipe analysis, food labeling, nutritional research.

FoodData Central integrates multiple data types: SR Legacy (research), Foundation Foods (analytical), FNDDS (survey), and Branded Foods (industry data).
