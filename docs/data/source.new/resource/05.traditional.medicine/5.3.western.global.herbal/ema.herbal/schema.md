---
id: schema-ema.herbal
title: "EMA Herbal Medicines Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: final
tags: [schema, herbal-medicine, european, regulatory, monographs, pharmacovigilance]
---

# EMA Herbal Medicines Schema Documentation

**Document ID:** SCHEMA-EMA-HERBAL
**Version:** 1.0
**Source Version:** EMA HMPC (Current)

---

## TL;DR

EMA Herbal Medicines provides regulatory-grade monographs for herbal substances in the EU. Each monograph contains standardized sections covering botanical identity, pharmaceutical quality, clinical efficacy, and safety. Two regulatory pathways exist: "well-established use" (clinical evidence) and "traditional use" (long-standing use). Data available as PDF documents with structured sections.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Herbal Monographs | 200+ | HMPC assessments |
| Community List Entries | 80+ | Harmonized EU list |
| Assessment Reports | 200+ | Full evidence reviews |
| Herbal Substances | 250+ | Assessed materials |
| Traditional Use Registrations | 150+ | Historical use basis |
| Well-Established Use | 50+ | Clinical evidence basis |

---

## Document Structure Overview

```
                    +---------------------+
                    |   Herbal Substance  |
                    | (Latin binomial +   |
                    |    plant part)      |
                    +----------+----------+
                               |
              +----------------+----------------+
              |                                 |
              v                                 v
    +---------+---------+           +---------+---------+
    |   Assessment      |           |    Community      |
    |    Report         |           |      List         |
    |   (Full review)   |           |   (Harmonized)    |
    +---------+---------+           +-------------------+
              |
              v
    +---------+---------+
    |    Monograph      |
    +-------------------+
              |
    +---------+---------+
    |                   |
    v                   v
+---+---+         +-----+-----+
| Well- |         |Traditional|
| Est.  |         |   Use     |
| Use   |         +-----------+
+-------+
```

---

## Core Entities

### HerbalSubstance

**Description:** The primary herbal material under assessment

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| substance_id | String | Yes | EMA reference number |
| latin_name | String | Yes | Botanical binomial |
| plant_part | String | Yes | Standardized part (radix, folium, etc.) |
| common_names | Array[String] | No | Vernacular names by language |
| family | String | Yes | Botanical family |
| pharmacopoeia_reference | String | No | European Pharmacopoeia monograph |
| synonyms | Array[String] | No | Taxonomic synonyms |

### Monograph

**Description:** HMPC monograph document

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| monograph_id | String | Yes | EMA/HMPC/xxxxxx/xxxx |
| herbal_substance | HerbalSubstance | Yes | Subject of monograph |
| adoption_date | Date | Yes | HMPC adoption date |
| revision_date | Date | No | Latest revision |
| status | Enum | Yes | Final, Draft, Under revision |
| use_category | Enum | Yes | Well-established, Traditional |
| sections | Array[MonographSection] | Yes | Document sections |

### MonographSection

**Description:** Standardized sections within monograph

| Section ID | Title | Content |
|------------|-------|---------|
| 1 | NAME OF THE HERBAL SUBSTANCE | Latin name, plant part |
| 2 | COMPOSITION | Active components, markers |
| 3 | PHARMACEUTICAL FORM | Preparations, extracts |
| 4.1 | THERAPEUTIC INDICATIONS | Approved uses |
| 4.2 | POSOLOGY AND METHOD | Dosage, administration |
| 4.3 | CONTRAINDICATIONS | When not to use |
| 4.4 | SPECIAL WARNINGS | Precautions |
| 4.5 | INTERACTIONS | Drug interactions |
| 4.6 | PREGNANCY/LACTATION | Safety in pregnancy |
| 4.7 | EFFECTS ON DRIVING | Impairment warnings |
| 4.8 | UNDESIRABLE EFFECTS | Adverse reactions |
| 4.9 | OVERDOSE | Overdose management |
| 5 | PHARMACOLOGICAL PROPERTIES | Mechanism, pharmacokinetics |
| 6 | NON-CLINICAL DATA | Toxicology studies |

### AssessmentReport

**Description:** Full HMPC assessment with evidence review

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| report_id | String | Yes | EMA reference |
| herbal_substance | HerbalSubstance | Yes | Subject |
| adoption_date | Date | Yes | Publication date |
| total_pages | Integer | No | Document length |
| sections | Array[AssessmentSection] | Yes | Evidence sections |
| references | Array[Reference] | Yes | Literature cited |

### CommunityListEntry

**Description:** EU Community list harmonized entry

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| entry_id | String | Yes | Community list reference |
| herbal_substance | HerbalSubstance | Yes | Listed substance |
| herbal_preparations | Array[Preparation] | Yes | Specified preparations |
| indication | String | Yes | Accepted indication |
| route | String | Yes | Administration route |
| posology | String | Yes | Dosage specification |
| period_of_use | String | No | Maximum duration |
| restrictions | Array[String] | No | Use restrictions |

---

## Preparation Types

### Preparation

**Description:** Pharmaceutical preparation of herbal substance

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| preparation_type | Enum | Yes | Type of preparation |
| description | String | Yes | Detailed specification |
| extraction_solvent | String | No | Solvent used |
| drug_extract_ratio | String | No | DER specification |
| standardization | String | No | Marker compound content |

### Preparation Types Enum

| Type | Description | Example |
|------|-------------|---------|
| Herbal substance | Dried plant material | Valerian root, dried |
| Powdered | Ground plant material | Powdered ginger root |
| Dry extract | Concentrated extract | Extract DER 4-7:1 |
| Soft extract | Semi-solid extract | Soft extract, ethanol |
| Liquid extract | Fluid extract | Liquid extract 1:1 |
| Tincture | Alcoholic solution | Tincture 1:5, 70% ethanol |
| Expressed juice | Fresh plant juice | Fresh juice, stabilized |
| Essential oil | Steam distilled | Peppermint essential oil |

---

## Regulatory Categories

### Use Categories

| Category | Evidence Basis | Requirements |
|----------|----------------|--------------|
| Well-Established Use | Clinical trials, systematic reviews | Documented efficacy and safety |
| Traditional Use | Long-standing use (>=30 years, >=15 in EU) | Plausibility, safety acceptable |

### Well-Established Use Criteria

| Criterion | Requirement |
|-----------|-------------|
| Time | >=10 years EU/>=30 years global use |
| Evidence | Clinical trials or systematic reviews |
| Safety | Adequate safety documentation |
| Bibliographic | Comprehensive literature review |

### Traditional Use Criteria

| Criterion | Requirement |
|-----------|-------------|
| Time | >=30 years use (>=15 in EU) |
| Plausibility | Pharmacologically plausible |
| Safety | No evidence of harm |
| Indication | Minor, self-limiting conditions |

---

## Safety Data Structure

### SafetyProfile

**Description:** Comprehensive safety information

| Field | Type | Description |
|-------|------|-------------|
| contraindications | Array[Contraindication] | Absolute exclusions |
| warnings | Array[Warning] | Special precautions |
| interactions | Array[Interaction] | Drug interactions |
| pregnancy | PregnancyData | Pregnancy/lactation |
| adverse_effects | Array[AdverseEffect] | Known side effects |
| overdose | OverdoseInfo | Overdose guidance |

### Contraindication

| Field | Type | Description |
|-------|------|-------------|
| condition | String | Medical condition |
| severity | Enum | Absolute/Relative |
| evidence_level | Enum | Clinical/Preclinical/Theoretical |

### Interaction

| Field | Type | Description |
|-------|------|-------------|
| interacting_substance | String | Drug or substance |
| effect | String | Nature of interaction |
| severity | Enum | Major/Moderate/Minor |
| mechanism | String | Pharmacological basis |
| evidence | String | Supporting evidence |

### AdverseEffect

| Field | Type | Description |
|-------|------|-------------|
| effect | String | Adverse reaction |
| frequency | Enum | Very common/Common/Uncommon/Rare/Very rare |
| meddra_term | String | MedDRA preferred term |
| system_organ_class | String | SOC classification |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | PDF documents |
| Structured | Semi-structured text sections |
| Encoding | UTF-8 |
| Languages | English (primary), EU languages |

---

## Sample Records

### Monograph Summary (Valeriana officinalis)

```json
{
  "monograph_id": "EMA/HMPC/575871/2007",
  "substance": {
    "latin_name": "Valeriana officinalis L.",
    "plant_part": "radix",
    "family": "Caprifoliaceae",
    "common_name": "Valerian root"
  },
  "adoption_date": "2016-03-22",
  "status": "Final",
  "use_categories": ["Traditional use"],
  "indications": [
    {
      "category": "Traditional use",
      "indication": "Traditional herbal medicinal product for relief of mild nervous tension and to aid sleep",
      "route": "Oral",
      "duration": "2 weeks maximum"
    }
  ],
  "preparations": [
    {
      "type": "Herbal substance",
      "description": "Dried, whole or cut root",
      "form": "Herbal tea"
    },
    {
      "type": "Dry extract",
      "description": "DER 4-7:1, ethanol 40-70%",
      "standardization": "Valerenic acids 0.1-0.8%"
    }
  ],
  "posology": {
    "adults": "Single dose 0.3-1g as tea, up to 3 times daily",
    "children": "Not recommended under 12 years"
  },
  "safety": {
    "contraindications": ["Hypersensitivity"],
    "interactions": ["May enhance sedatives"],
    "pregnancy": "Not recommended (insufficient data)",
    "adverse_effects": ["GI disorders (rare)"]
  }
}
```

### Community List Entry

```json
{
  "entry_id": "COMM/2014/01",
  "substance": {
    "latin_name": "Matricaria recutita L.",
    "plant_part": "flos",
    "common_name": "Chamomile flower"
  },
  "use_category": "Traditional use",
  "indications": [
    {
      "indication": "Minor skin inflammations and minor wounds",
      "route": "Cutaneous",
      "preparation": "Liquid extract"
    },
    {
      "indication": "Minor inflammations of mouth and throat",
      "route": "Oromucosal",
      "preparation": "Herbal tea"
    }
  ],
  "restrictions": [
    "Not to be used in children under 12 years",
    "Maximum 2 weeks without medical advice"
  ]
}
```

---

## Quality Standards

### Pharmacopoeia References

| Standard | Description |
|----------|-------------|
| Ph. Eur. | European Pharmacopoeia |
| Herbal Drug | Whole or processed plant material |
| Herbal Preparation | Extract, tincture, powder |
| DER | Drug-Extract Ratio |
| Marker Compound | Standardization reference |

### Standardization Approaches

| Approach | Description | Example |
|----------|-------------|---------|
| Native extract | No standardization | Simple dry extract |
| Quantified | Defined marker range | Valerenic acids 0.1-0.8% |
| Standardized | Adjusted to marker | Hypericin 0.1-0.3% |

---

## Cross-References

| External Source | Reference Type |
|-----------------|----------------|
| European Pharmacopoeia | Quality monographs |
| WHO Monographs | Global herbal references |
| ESCOP | European scientific assessment |
| EFSA | Food safety opinions |
| MedDRA | Adverse effect terminology |
| ATC Code | Therapeutic classification |

---

## Glossary

| Term | Definition |
|------|------------|
| HMPC | Committee on Herbal Medicinal Products |
| Well-Established Use | Efficacy proven by clinical evidence |
| Traditional Use | Based on long-standing historical use |
| Community List | EU harmonized herbal substances |
| DER | Drug-Extract Ratio (quantity plant:extract) |
| Ph. Eur. | European Pharmacopoeia |
| Pharmacovigilance | Ongoing safety monitoring |
| Marker Compound | Chemical for standardization |

---

## References

1. EMA HMPC Guidance: https://www.ema.europa.eu/en/human-regulatory/herbal-medicinal-products
2. EU Community Herbal Monographs: EMA website
3. European Pharmacopoeia: EDQM
