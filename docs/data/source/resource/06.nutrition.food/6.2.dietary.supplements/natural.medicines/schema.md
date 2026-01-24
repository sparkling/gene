---
id: schema-natural-medicines
title: "Natural Medicines - Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [schema, natural-products, evidence, efficacy, interactions]
---

# Natural Medicines Comprehensive Database - Schema Documentation

**Document ID:** SCHEMA-NATURAL-MEDICINES
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://naturalmedicines.therapeuticresearch.com/

---

## TL;DR

The Natural Medicines Comprehensive Database (formerly by TRC Healthcare) is a clinical decision support resource providing evidence-based information about natural products, dietary supplements, and complementary therapies. It includes detailed monographs with effectiveness ratings, safety grades, drug interactions, and dosing guidelines.

---

## Database Statistics

| Metric | Count |
|--------|-------|
| Monographs | 1,400+ |
| Drug Interactions | 175,000+ |
| Effectiveness Ratings | Systematic |
| Interaction Severity Levels | 5 |
| Safety Grades | 6 |

---

## Data Schema

### Core Entities

#### 1. Ingredient Monographs
```json
{
  "monograph_id": "NM-001234",
  "name": "Turmeric",
  "scientific_names": ["Curcuma longa", "Curcuma domestica"],
  "common_names": ["Curcumin", "Indian Saffron", "Haridra"],
  "family": "Zingiberaceae",
  "parts_used": ["Rhizome", "Root"],
  "overview": "Turmeric is a plant native to South Asia...",
  "safety_grade": "Likely Safe",
  "effectiveness_ratings": {
    "osteoarthritis": {
      "rating": "Possibly Effective",
      "evidence_level": "B",
      "summary": "Clinical evidence suggests..."
    },
    "depression": {
      "rating": "Insufficient Evidence",
      "evidence_level": "C",
      "summary": "Preliminary research..."
    }
  },
  "mechanism_of_action": "Curcumin inhibits...",
  "pharmacokinetics": {
    "absorption": "Poor oral bioavailability",
    "metabolism": "Hepatic glucuronidation",
    "half_life": "1-2 hours"
  }
}
```

#### 2. Drug Interactions
```json
{
  "interaction_id": "INT-56789",
  "natural_product": "St. John's Wort",
  "drug": "Warfarin",
  "severity": "Major",
  "documentation": "Well-documented",
  "mechanism": "CYP3A4 induction",
  "clinical_significance": "May reduce anticoagulant effect",
  "management": "Avoid combination or monitor INR closely",
  "references": ["PMID:12345678", "PMID:23456789"]
}
```

#### 3. Disease-Based Evidence
```json
{
  "condition_id": "COND-001",
  "condition_name": "Osteoarthritis",
  "natural_treatments": [
    {
      "ingredient": "Glucosamine",
      "effectiveness": "Possibly Effective",
      "typical_dose": "1500 mg daily",
      "evidence_summary": "Multiple RCTs..."
    },
    {
      "ingredient": "Chondroitin",
      "effectiveness": "Possibly Effective",
      "typical_dose": "800-1200 mg daily"
    }
  ]
}
```

---

## Effectiveness Rating Scale

| Rating | Description | Evidence Level |
|--------|-------------|----------------|
| Effective | Strong evidence from well-designed trials | A |
| Likely Effective | Good evidence from clinical trials | A-B |
| Possibly Effective | Some evidence, more research needed | B |
| Possibly Ineffective | Evidence suggests lack of benefit | B |
| Likely Ineffective | Good evidence of no benefit | A-B |
| Ineffective | Strong evidence of no benefit | A |
| Insufficient Evidence | Not enough data to evaluate | C-D |

---

## Safety Rating Scale

| Rating | Description |
|--------|-------------|
| Likely Safe | No significant safety concerns when used appropriately |
| Possibly Safe | Limited safety data, no significant concerns reported |
| Possibly Unsafe | Reports of adverse effects; use with caution |
| Likely Unsafe | Significant adverse effects documented |
| Unsafe | Do not use; serious adverse effects |
| Insufficient Information | Not enough safety data available |

---

## Drug Interaction Severity

| Severity | Description | Action |
|----------|-------------|--------|
| Contraindicated | Do not use together | Avoid combination |
| Major | Life-threatening or permanent effects possible | Use alternative or monitor closely |
| Moderate | May worsen condition or require adjustment | Monitor and consider alternatives |
| Minor | Limited clinical significance | Be aware; usually no action needed |
| Unknown | Insufficient data to assess | Use caution |

---

## Evidence Levels

| Level | Description |
|-------|-------------|
| A | High quality RCTs or meta-analyses |
| B | Moderate quality clinical studies |
| C | Low quality studies or case reports |
| D | Traditional use or expert opinion only |

---

## Data Access

### Web Interface (Subscription Required)

| Feature | Access Level |
|---------|--------------|
| Monograph search | Subscriber |
| Interaction checker | Subscriber |
| Patient handouts | Subscriber |
| Mobile app | Subscriber |

### API Access

- **Status:** Enterprise only
- **Contact:** TRC Healthcare sales
- **Integration:** EPIC, Cerner, other EHR systems

### No Bulk Downloads

Natural Medicines does not provide bulk data downloads. Data is proprietary and accessible only through subscription.

---

## Monograph Structure

### Standard Sections

| Section | Content |
|---------|---------|
| Overview | General description |
| Scientific Names | Taxonomic names |
| Common Names | Regional/trade names |
| Uses | Traditional and modern uses |
| Safety | Pregnancy, lactation, contraindications |
| Effectiveness | Evidence-based ratings by condition |
| Mechanism of Action | Pharmacological basis |
| Adverse Reactions | Side effects |
| Interactions | Drug, herb, food interactions |
| Dosing | Recommended doses by condition |
| References | Primary literature citations |

---

## Use Cases

### 1. Check Drug-Supplement Interaction
```
1. Navigate to Interaction Checker
2. Enter drug name (e.g., "Warfarin")
3. Enter supplement name (e.g., "Vitamin E")
4. Review interaction details and severity
```

### 2. Evaluate Evidence for Condition
```
1. Search by condition (e.g., "Depression")
2. View ranked treatments by effectiveness
3. Review evidence summaries
4. Download patient handout
```

### 3. Create Patient Education
```
1. Search supplement monograph
2. Select "Patient Handout"
3. Customize language/detail level
4. Print or email to patient
```

---

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Monograph ID | `NM-NNNNNN` | NM-001234 |
| Interaction ID | `INT-NNNNN` | INT-56789 |
| CAS Number | Standard | 458-37-7 |
| UNII | FDA code | IT942ZTH98 |

---

## Relationships

### Entity Relationships
```
ingredients (1) ----< interactions (N)
    |
    +---- (N) conditions (M) [effectiveness ratings]
    |
    +---- (N) drugs (M) [interactions]
    |
    +---- (N) references (M)
```

---

## License

- **License Type:** Proprietary
- **Subscription Required:** Yes (institutional/individual)
- **Commercial Use:** Licensed separately
- **Data Export:** Not permitted
- **EHR Integration:** Available for enterprise

---

## Citation Format

```
Natural Medicines. [Monograph Name].
Natural Medicines Comprehensive Database.
TRC Healthcare.
https://naturalmedicines.therapeuticresearch.com/
[Access date]
```

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Monograph | Comprehensive ingredient documentation | Turmeric monograph |
| Effectiveness Rating | Evidence-based efficacy assessment | Possibly Effective |
| Safety Grade | Overall safety classification | Likely Safe |
| Interaction Severity | Clinical significance of combination | Major |
| Evidence Level | Quality of supporting research | Level A |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| CYP Induction | Increased drug metabolism enzyme activity | Drug interactions |
| CYP Inhibition | Decreased drug metabolism enzyme activity | Drug interactions |
| Bioavailability | Fraction of dose reaching systemic circulation | Pharmacokinetics |
| RCT | Randomized Controlled Trial | Evidence level |
| Meta-analysis | Statistical combination of studies | Evidence level |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NM | Natural Medicines | Database name |
| TRC | Therapeutic Research Center | Publisher |
| CYP | Cytochrome P450 | Drug metabolism enzymes |
| RCT | Randomized Controlled Trial | Study design |
| UNII | Unique Ingredient Identifier | FDA identifier |
| EHR | Electronic Health Record | Integration target |
| INR | International Normalized Ratio | Warfarin monitoring |

---

## Sample Data

### Example Ingredient Monograph Record
```json
{
  "monograph_id": "NM-001234",
  "name": "Turmeric",
  "scientific_names": ["Curcuma longa", "Curcuma domestica"],
  "common_names": [
    "Curcumin",
    "Indian Saffron",
    "Haridra",
    "Haldi",
    "Yellow Root",
    "Curcuma"
  ],
  "family": "Zingiberaceae",
  "parts_used": ["Rhizome", "Root"],
  "overview": "Turmeric is a plant native to South Asia, widely used as a spice and traditional medicine. The primary active compounds are curcuminoids, particularly curcumin, which has antioxidant and anti-inflammatory properties.",
  "safety_grade": "Likely Safe",
  "pregnancy_lactation": {
    "pregnancy": "Possibly Unsafe in medicinal amounts",
    "lactation": "Insufficient Reliable Information"
  },
  "effectiveness_ratings": [
    {
      "condition": "Osteoarthritis",
      "rating": "Possibly Effective",
      "evidence_level": "B",
      "summary": "Clinical evidence suggests that taking curcumin 500 mg twice daily for 2-4 months can reduce pain and improve function in people with osteoarthritis.",
      "typical_dose": "500-2000 mg curcumin daily"
    },
    {
      "condition": "Hyperlipidemia",
      "rating": "Possibly Effective",
      "evidence_level": "B",
      "summary": "Some clinical research shows that taking turmeric extract can reduce total cholesterol and LDL cholesterol levels.",
      "typical_dose": "1000-2000 mg daily"
    },
    {
      "condition": "Depression",
      "rating": "Possibly Effective",
      "evidence_level": "B",
      "summary": "Curcumin may help reduce symptoms of major depressive disorder when used adjunctively with standard therapy.",
      "typical_dose": "500-1000 mg curcumin daily"
    },
    {
      "condition": "Inflammatory Bowel Disease",
      "rating": "Insufficient Evidence",
      "evidence_level": "C",
      "summary": "Preliminary research shows mixed results for ulcerative colitis. More research is needed.",
      "typical_dose": "Variable"
    }
  ],
  "mechanism_of_action": "Curcumin inhibits multiple inflammatory pathways including COX-2, LOX, NF-kB, and various cytokines. It also exhibits antioxidant activity through free radical scavenging.",
  "pharmacokinetics": {
    "absorption": "Poor oral bioavailability (1-2%)",
    "bioavailability_enhancement": "Improved with piperine, lipid formulations, or nanoparticle delivery",
    "metabolism": "Extensive hepatic glucuronidation and sulfation",
    "half_life": "1-2 hours for curcumin",
    "elimination": "Primarily fecal"
  },
  "adverse_effects": [
    {
      "effect": "GI upset",
      "frequency": "Common at high doses",
      "severity": "Mild"
    },
    {
      "effect": "Increased bleeding risk",
      "frequency": "Theoretical",
      "severity": "Potentially significant with anticoagulants"
    }
  ],
  "contraindications": [
    "Bile duct obstruction",
    "Gallstones (may stimulate bile flow)"
  ]
}
```

### Example Drug Interaction Record
```json
{
  "interaction_id": "INT-56789",
  "natural_product": {
    "id": "NM-002345",
    "name": "St. John's Wort",
    "scientific_name": "Hypericum perforatum"
  },
  "drug": {
    "generic_name": "Warfarin",
    "brand_names": ["Coumadin", "Jantoven"],
    "drug_class": "Anticoagulants"
  },
  "severity": "Major",
  "documentation": "Well-documented",
  "onset": "Delayed (1-2 weeks)",
  "mechanism": {
    "primary": "CYP3A4 induction",
    "secondary": "CYP2C9 induction",
    "description": "St. John's Wort induces hepatic cytochrome P450 enzymes, increasing warfarin metabolism and reducing anticoagulant effect."
  },
  "clinical_significance": "May significantly reduce anticoagulant effect, increasing risk of thromboembolism",
  "evidence": [
    {
      "type": "Clinical study",
      "summary": "INR decreased by average of 2 points in patients taking both",
      "reference": "PMID:12345678"
    },
    {
      "type": "Case report",
      "summary": "Multiple cases of subtherapeutic INR documented",
      "reference": "PMID:23456789"
    }
  ],
  "management": {
    "recommendation": "Avoid combination",
    "alternatives": "Consider other antidepressants that don't induce CYP enzymes",
    "monitoring": "If combination cannot be avoided, monitor INR closely for 2-4 weeks after starting or stopping St. John's Wort"
  }
}
```

### Example Condition-Based Record
```json
{
  "condition_id": "COND-001",
  "condition_name": "Osteoarthritis",
  "icd10_codes": ["M15", "M16", "M17", "M18", "M19"],
  "overview": "Osteoarthritis is a degenerative joint disease characterized by cartilage breakdown, pain, and functional impairment.",
  "natural_treatments": [
    {
      "ingredient": "Glucosamine Sulfate",
      "monograph_id": "NM-000567",
      "effectiveness": "Possibly Effective",
      "evidence_level": "B",
      "typical_dose": "1500 mg daily or 500 mg three times daily",
      "onset": "4-8 weeks for symptomatic benefit",
      "evidence_summary": "Multiple RCTs show modest benefit for pain and function. European formulations may be more effective than US products.",
      "notes": "Sulfate form preferred; glucosamine HCl shows less consistent benefit"
    },
    {
      "ingredient": "Chondroitin Sulfate",
      "monograph_id": "NM-000568",
      "effectiveness": "Possibly Effective",
      "evidence_level": "B",
      "typical_dose": "800-1200 mg daily",
      "onset": "2-4 months",
      "evidence_summary": "May reduce pain and slow structural progression. Often combined with glucosamine."
    },
    {
      "ingredient": "Turmeric/Curcumin",
      "monograph_id": "NM-001234",
      "effectiveness": "Possibly Effective",
      "evidence_level": "B",
      "typical_dose": "500-2000 mg curcumin extract daily",
      "onset": "4-8 weeks",
      "evidence_summary": "Anti-inflammatory effects may reduce pain comparable to NSAIDs in some studies."
    },
    {
      "ingredient": "SAMe",
      "monograph_id": "NM-000789",
      "effectiveness": "Possibly Effective",
      "evidence_level": "B",
      "typical_dose": "600-1200 mg daily",
      "onset": "1-2 weeks",
      "evidence_summary": "May be as effective as NSAIDs for pain relief with fewer GI side effects."
    },
    {
      "ingredient": "Avocado-Soybean Unsaponifiables",
      "monograph_id": "NM-000890",
      "effectiveness": "Possibly Effective",
      "evidence_level": "B",
      "typical_dose": "300-600 mg daily",
      "onset": "2-3 months",
      "evidence_summary": "May reduce pain and improve function; may have structure-modifying effects."
    }
  ]
}
```

**Note:** The sample data above represents the structure of data available through Natural Medicines subscription interface. Actual data values are proprietary and for illustrative purposes only. Data is based on publicly documented schema patterns.

---

## Related Documents

- [Download Instructions](./download.md)
- [DSLD](../dsld/_index.md) - NIH supplement label database
- [ConsumerLab](../consumerlab/_index.md) - Product testing
