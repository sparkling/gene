---
id: schema-consumerlab
title: "ConsumerLab - Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [schema, supplements, testing, quality]
---

# ConsumerLab - Schema Documentation

**Document ID:** SCHEMA-CONSUMERLAB
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://www.consumerlab.com/

---

## TL;DR

ConsumerLab is an independent organization that tests health and nutrition products and provides quality evaluations of dietary supplements. The database contains test results, product reviews, ingredient comparisons, and evidence assessments. Data is proprietary and requires subscription access.

---

## Database Statistics

| Metric | Count |
|--------|-------|
| Products Tested | 6,500+ |
| Brands Covered | 800+ |
| Product Categories | 50+ |
| Ingredient Monographs | 300+ |
| Drug Interactions | 100+ |

---

## Data Schema

### Core Entities

#### 1. Product Reviews
```json
{
  "review_id": "CL-2024-VIT-D-001",
  "product_name": "Vitamin D3 5000 IU",
  "brand": "Nature Made",
  "category": "Vitamins & Minerals",
  "subcategory": "Vitamin D",
  "upc": "012345678901",
  "test_date": "2024-01-15",
  "approval_status": "Approved",
  "test_results": {
    "claimed_amount": "5000 IU",
    "actual_amount": "5120 IU",
    "percent_of_claim": 102.4,
    "contaminants": {
      "lead": "<0.5 mcg",
      "arsenic": "Not detected",
      "cadmium": "Not detected",
      "mercury": "Not detected"
    },
    "disintegration": "Pass",
    "bioavailability": "Adequate"
  }
}
```

#### 2. Ingredient Monographs
```json
{
  "ingredient_id": "ING-VIT-D",
  "name": "Vitamin D",
  "common_names": ["Vitamin D3", "Cholecalciferol", "Vitamin D2", "Ergocalciferol"],
  "category": "Vitamins",
  "evidence_summary": {
    "bone_health": "Strong evidence",
    "immune_function": "Moderate evidence",
    "cancer_prevention": "Limited evidence"
  },
  "recommended_intake": {
    "adults": "600-800 IU/day",
    "upper_limit": "4000 IU/day"
  },
  "forms": ["D3 (Cholecalciferol)", "D2 (Ergocalciferol)"],
  "drug_interactions": ["Corticosteroids", "Weight-loss drugs", "Seizure medications"]
}
```

#### 3. Test Categories
```json
{
  "category_id": "CAT-001",
  "name": "Vitamins & Minerals",
  "subcategories": [
    "Multivitamins",
    "Vitamin D",
    "Vitamin B12",
    "Calcium",
    "Magnesium",
    "Iron"
  ],
  "test_criteria": {
    "identity": "Confirm ingredient identity",
    "potency": "Verify claimed amounts",
    "purity": "Check for contaminants",
    "disintegration": "Dissolution testing"
  }
}
```

---

## Quality Assessment Framework

### Approval Criteria

| Criterion | Description | Threshold |
|-----------|-------------|-----------|
| Identity | Correct ingredient present | Required |
| Potency | Amount vs. label claim | 90-110% of claim |
| Purity | Heavy metal content | Below USP limits |
| Disintegration | Tablet/capsule breakdown | <30 minutes |
| Contamination | Microbial, pesticides | Not detected |

### Product Status

| Status | Description |
|--------|-------------|
| Approved | Meets all quality criteria |
| Not Approved | Failed one or more tests |
| Voluntarily Withdrawn | Removed by manufacturer |
| Not Tested | Awaiting evaluation |

---

## Test Result Schema

### Heavy Metal Testing
```json
{
  "heavy_metals": {
    "lead_mcg_daily": {
      "value": 0.3,
      "limit": 0.5,
      "status": "Pass"
    },
    "arsenic_mcg_daily": {
      "value": "ND",
      "limit": 10.0,
      "status": "Pass"
    },
    "cadmium_mcg_daily": {
      "value": 0.1,
      "limit": 4.1,
      "status": "Pass"
    },
    "mercury_mcg_daily": {
      "value": "ND",
      "limit": 2.0,
      "status": "Pass"
    }
  }
}
```

### Potency Testing
```json
{
  "potency_results": [
    {
      "ingredient": "Vitamin D3",
      "claimed": "5000 IU",
      "measured": "5120 IU",
      "percent_claim": 102.4,
      "acceptable_range": "90-110%",
      "status": "Pass"
    }
  ]
}
```

---

## Data Access

### Web Interface

| Access Level | Features |
|--------------|----------|
| Free | Limited product listings, headlines |
| Subscription | Full test results, comparisons, monographs |

### No API Available

ConsumerLab does not provide public API access. Data must be accessed through the web interface.

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Review ID | `CL-YYYY-CAT-NNN` | CL-2024-VIT-D-001 |
| Product UPC | Standard barcode | 012345678901 |
| Ingredient ID | `ING-XXX` | ING-VIT-D |

---

## Relationships

### Entity Relationships
```
products (N) ----< reviews (N)
    |
    +---- (N) ingredients (M)
    |
    +---- (1) brands (N)
    |
    +---- (N) categories (M)
```

---

## Use Cases

### 1. Find Approved Vitamin D Products
```
Navigate to: ConsumerLab.com > Reviews > Vitamins > Vitamin D
Filter by: Approved products
Sort by: Cost per dose
```

### 2. Check Drug Interactions
```
Navigate to: ConsumerLab.com > Encyclopedia > [Ingredient]
View: Drug Interactions section
Note: Severity ratings and mechanisms
```

### 3. Compare Products
```
Navigate to: ConsumerLab.com > Reviews > [Category]
Use: Comparison tool
Compare: Price, potency, purity across brands
```

---

## License

- **License Type:** Proprietary
- **Subscription Required:** Yes ($49.95/year individual)
- **Commercial Use:** Requires license agreement
- **Data Export:** Not permitted
- **Citation:** Allowed with attribution

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Approved | Product passed all quality tests | CL seal displayed |
| Not Approved | Product failed one or more tests | Detailed reasons provided |
| Potency | Amount of active ingredient present | 102% of label claim |
| Disintegration | Time for tablet to break down | <30 minutes |
| USP | United States Pharmacopeia standards | Heavy metal limits |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CL | ConsumerLab | Testing organization |
| USP | United States Pharmacopeia | Quality standards |
| UPC | Universal Product Code | Barcode identifier |
| IU | International Unit | Vitamin measurement |
| ND | Not Detected | Below detection limit |
| mcg | Microgram | Weight unit |

---

## Sample Data

### Example Product Review Record
```json
{
  "review_id": "CL-2024-VIT-D-001",
  "product_name": "Vitamin D3 5000 IU",
  "brand": "Nature Made",
  "upc": "012345678901",
  "category": "Vitamins & Minerals",
  "subcategory": "Vitamin D",
  "test_date": "2024-01-15",
  "approval_status": "Approved",
  "test_results": {
    "identity_test": {
      "status": "Pass",
      "notes": "Correct vitamin D3 form confirmed"
    },
    "potency_results": [
      {
        "ingredient": "Vitamin D3",
        "claimed": "5000 IU",
        "measured": "5120 IU",
        "percent_claim": 102.4,
        "acceptable_range": "90-110%",
        "status": "Pass"
      }
    ],
    "heavy_metals": {
      "lead_mcg_daily": {
        "value": 0.3,
        "limit": 0.5,
        "status": "Pass"
      },
      "arsenic_mcg_daily": {
        "value": "ND",
        "limit": 10.0,
        "status": "Pass"
      },
      "cadmium_mcg_daily": {
        "value": 0.1,
        "limit": 4.1,
        "status": "Pass"
      },
      "mercury_mcg_daily": {
        "value": "ND",
        "limit": 2.0,
        "status": "Pass"
      }
    },
    "disintegration": {
      "time_minutes": 18,
      "limit_minutes": 30,
      "status": "Pass"
    },
    "contaminant_screen": {
      "microbiological": "Pass",
      "pesticides": "Not detected",
      "solvents": "Not detected"
    }
  },
  "cost_analysis": {
    "price_paid": 12.99,
    "servings_per_container": 90,
    "cost_per_serving": 0.14,
    "cost_per_1000_iu": 0.03
  },
  "quality_seal": true,
  "date_approved": "2024-01-20",
  "last_updated": "2024-06-15"
}
```

### Example Ingredient Monograph Record
```json
{
  "ingredient_id": "ING-VIT-D",
  "name": "Vitamin D",
  "common_names": [
    "Vitamin D3",
    "Cholecalciferol",
    "Vitamin D2",
    "Ergocalciferol",
    "Calcifediol",
    "25-hydroxyvitamin D"
  ],
  "category": "Vitamins",
  "subcategory": "Fat-Soluble Vitamins",
  "overview": "Vitamin D is a fat-soluble vitamin essential for calcium absorption...",
  "evidence_summary": {
    "bone_health": {
      "rating": "Strong evidence",
      "summary": "Well-established role in calcium metabolism and bone mineralization"
    },
    "immune_function": {
      "rating": "Moderate evidence",
      "summary": "Growing evidence for immune modulation; optimal levels unclear"
    },
    "cancer_prevention": {
      "rating": "Limited evidence",
      "summary": "Observational data; clinical trials inconclusive"
    },
    "cardiovascular_health": {
      "rating": "Insufficient evidence",
      "summary": "Mixed results from supplementation trials"
    }
  },
  "recommended_intake": {
    "infants_0_12_mo": "400 IU/day",
    "children_1_18_yr": "600 IU/day",
    "adults_19_70_yr": "600 IU/day",
    "adults_over_70": "800 IU/day",
    "pregnancy_lactation": "600 IU/day",
    "upper_limit_adults": "4000 IU/day"
  },
  "forms": [
    {
      "name": "D3 (Cholecalciferol)",
      "source": "Animal/lichen",
      "notes": "More effective at raising blood levels"
    },
    {
      "name": "D2 (Ergocalciferol)",
      "source": "Plant/fungal",
      "notes": "Vegan option; may require higher doses"
    }
  ],
  "drug_interactions": [
    {
      "drug": "Corticosteroids",
      "severity": "Moderate",
      "effect": "May reduce vitamin D absorption"
    },
    {
      "drug": "Orlistat",
      "severity": "Moderate",
      "effect": "Reduces fat-soluble vitamin absorption"
    },
    {
      "drug": "Thiazide diuretics",
      "severity": "Moderate",
      "effect": "May increase calcium levels"
    }
  ],
  "safety_considerations": [
    "Toxicity possible at very high doses (>10,000 IU/day long-term)",
    "Monitor calcium levels with high-dose supplementation",
    "Caution with kidney disease, sarcoidosis"
  ],
  "testing_notes": "ConsumerLab tests for D3/D2 content and heavy metals"
}
```

**Note:** The sample data above represents the structure of data available through ConsumerLab's subscription interface. Actual data values are proprietary and for illustrative purposes only.

---

## Related Documents

- [Download Instructions](./download.md)
- [DSLD](../dsld/_index.md) - NIH supplement label database
- [Natural Medicines](../natural.medicines/_index.md) - Evidence database
