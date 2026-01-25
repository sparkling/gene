---
id: schema-etcm
title: "ETCM Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: active
tags: [schema, tcm, traditional-chinese-medicine, pharmacology]
---

# ETCM Database Schema

**Document ID:** ETCM-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** ETCM

---

## TL;DR

ETCM (Encyclopedia of Traditional Chinese Medicine) provides comprehensive TCM data organized according to classical theory. Contains 403 herbs, 3,677 formulas, 7,274 compounds, 3,000+ targets, and 500+ diseases. Emphasizes TCM properties (nature, flavor, meridian tropism) alongside modern pharmacological data.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **TCM Herbs** | 403 |
| **TCM Formulas** | 3,677 |
| **Compounds** | 7,274 |
| **Target Proteins** | 3,000+ |
| **Diseases** | 500+ |
| **Herb-Compound Links** | 23,000+ |

---

## Data Model

### Entity Relationships

```
TCM Theory Properties
        |
    Herbs (403) -----> Formulas (3,677)
        |                    |
        v                    v
  Compounds (7,274) <--------+
        |
        v
   Targets (3,000+)
        |
        v
   Diseases (500+)
```

### Herb Entity

```json
{
  "herb_id": "string (ETCM internal, e.g., H001)",
  "herb_name_chinese": "string",
  "herb_name_pinyin": "string",
  "herb_name_english": "string",
  "herb_name_latin": "string (botanical name)",
  "tcm_properties": {
    "nature": "cold|cool|neutral|warm|hot",
    "flavor": ["sour", "bitter", "sweet", "pungent", "salty"],
    "meridian_tropism": ["lung", "heart", "liver", "spleen", "kidney", "stomach", "gallbladder", "bladder", "large_intestine", "small_intestine", "triple_burner", "pericardium"]
  },
  "actions": ["string array (therapeutic effects)"],
  "indications": ["string array (disease patterns)"],
  "formula_count": "integer",
  "compound_count": "integer"
}
```

### Formula Entity

```json
{
  "formula_id": "string (ETCM internal, e.g., F001)",
  "formula_name_chinese": "string",
  "formula_name_pinyin": "string",
  "formula_name_english": "string",
  "classical_source": "string (historical text reference)",
  "composition": [
    {
      "herb_id": "string",
      "herb_name": "string",
      "role": "jun|chen|zuo|shi (sovereign|minister|assistant|envoy)",
      "dosage": "string"
    }
  ],
  "therapeutic_category": "string",
  "actions": ["string array"],
  "indications": ["string array"]
}
```

### Compound Entity

```json
{
  "compound_id": "string (ETCM internal, e.g., C001)",
  "pubchem_cid": "integer",
  "compound_name": "string",
  "smiles": "string",
  "inchi": "string",
  "molecular_formula": "string",
  "molecular_weight": "float",
  "source_herbs": ["herb_id array"],
  "predicted_targets": ["target_id array"]
}
```

### Target Entity

```json
{
  "target_id": "string (UniProt accession)",
  "gene_symbol": "string",
  "gene_name": "string",
  "protein_class": "string",
  "source_compounds": ["compound_id array"],
  "disease_associations": ["disease_id array"]
}
```

### Disease Entity

```json
{
  "disease_id": "string (ETCM internal or external)",
  "disease_name": "string",
  "mesh_id": "string",
  "associated_targets": ["target_id array"],
  "associated_herbs": ["herb_id array"],
  "associated_formulas": ["formula_id array"]
}
```

---

## TCM Properties System

### Nature (Xing) Classification

| Nature | Description | Effect |
|--------|-------------|--------|
| Cold (Han) | Strongly cooling | Clears heat, purges fire |
| Cool (Liang) | Mildly cooling | Clears heat gently |
| Neutral (Ping) | Neither warm nor cold | Balanced effect |
| Warm (Wen) | Mildly warming | Warms interior |
| Hot (Re) | Strongly warming | Expels cold, warms yang |

### Flavor (Wei) Classification

| Flavor | TCM Function |
|--------|--------------|
| Sour (Suan) | Astringent, consolidates |
| Bitter (Ku) | Drains, dries dampness |
| Sweet (Gan) | Tonifies, harmonizes |
| Pungent (Xin) | Disperses, promotes circulation |
| Salty (Xian) | Softens, purges |

### Meridian Tropism (Gui Jing)

The 12 primary meridians plus pericardium:
- Lung, Heart, Liver, Spleen, Kidney
- Stomach, Gallbladder, Bladder
- Large Intestine, Small Intestine
- Triple Burner (San Jiao)
- Pericardium (Xin Bao)

---

## Formula Composition Roles

| Role | Chinese | Function |
|------|---------|----------|
| Sovereign (Jun) | Main therapeutic ingredient |
| Minister (Chen) | Enhances sovereign effect |
| Assistant (Zuo) | Reduces toxicity, treats secondary symptoms |
| Envoy (Shi) | Directs formula to target area |

---

## Cross-References

### Identifier Mappings

| ETCM Field | External Database | Notes |
|------------|-------------------|-------|
| compound_id | PubChem CID | Via structure |
| target_id | UniProt | Protein accession |
| disease_id | MESH | Disease concepts |
| herb_name_latin | Plant databases | Taxonomic reference |

---

## Sample Data

### Sample Herb Record

| Field | Value |
|-------|-------|
| herb_id | H001 |
| herb_name_pinyin | Huang Qi |
| nature | Warm |
| flavor | Sweet |
| meridian_tropism | Lung, Spleen |
| compound_count | 89 |

### Sample Formula Record

| Field | Value |
|-------|-------|
| formula_id | F001 |
| formula_name_pinyin | Bu Zhong Yi Qi Tang |
| herb_count | 8 |
| sovereign_herb | Huang Qi |

---

## Integration Recommendations

### Complementary Databases

| Database | Integration Purpose |
|----------|---------------------|
| BATMAN-TCM | Larger compound-target predictions |
| TCMBank | More compound structures |
| HERB | Gene expression validation |
| SymMap | Symptom mapping |

---

## License

| Aspect | Value |
|--------|-------|
| License | Academic use |
| Commercial Use | Contact maintainers |
| Attribution | Required for publications |

---

## Limitations

1. No public REST API
2. Manual download required for bulk data
3. Less comprehensive than BATMAN-TCM for predicted interactions
4. TCM theory concepts may require domain expertise

---

## Glossary

| Term | Definition |
|------|------------|
| Xing | TCM nature/thermal property |
| Wei | TCM flavor classification |
| Gui Jing | Meridian tropism (target organs) |
| Jun-Chen-Zuo-Shi | Classical formula composition roles |
| Zheng | TCM pattern/syndrome |

---

## References

1. Xu HY, et al. (2019). ETCM: an encyclopaedia of traditional Chinese medicine. Nucleic Acids Research.
2. Official Website: http://www.tcmip.cn/ETCM/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
