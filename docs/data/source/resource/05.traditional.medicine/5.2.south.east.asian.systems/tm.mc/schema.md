---
id: schema-tmmc
title: "TM-MC Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: active
tags: [schema, traditional-medicine, multi-system, molecular-targets]
---

# TM-MC Database Schema

**Document ID:** TMMC-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** TM-MC

---

## TL;DR

TM-MC (Traditional Medicine - Molecular Correlates) integrates traditional medicine knowledge from multiple Asian systems (Ayurveda, TCM, Unani, Siddha, Kampo) to modern molecular understanding. Contains 2,500+ medicinal plants, 5+ traditional systems, 8,000+ compounds, 1,500+ molecular targets, 500+ therapeutic indications, and 25,000+ plant-compound links. Enables cross-cultural comparative analysis and convergent validation.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **Medicinal Plants** | 2,500+ |
| **Traditional Systems** | 5+ |
| **Compounds** | 8,000+ |
| **Molecular Targets** | 1,500+ |
| **Therapeutic Indications** | 500+ |
| **Plant-Compound Links** | 25,000+ |

---

## Data Model

### Entity Relationships

```
Traditional Systems (Ayurveda, TCM, Unani, Siddha, etc.)
                    |
                    v
          Medicinal Plants (2,500+)
                    |
         +----------+----------+
         |                     |
         v                     v
  Compounds (8,000+)    Therapeutic Uses (500+)
         |
         v
  Molecular Targets (1,500+)
         |
         v
  Pathways & Mechanisms
```

### Plant Entity (Core)

```json
{
  "plant_id": "string (TMMC internal, e.g., TMMC_P_001234)",
  "botanical_name": "string (Latin binomial)",
  "family": "string",
  "synonyms": ["string array"],
  "common_names": {
    "english": "string",
    "regional": ["string array"]
  },
  "traditional_systems": [
    {
      "system_code": "AYU|TCM|UNA|SID|KAM|JAM",
      "system_name": "string",
      "local_name": "string",
      "therapeutic_category": "string",
      "uses": ["string array"]
    }
  ],
  "compounds_isolated": ["compound_id array"],
  "therapeutic_indications": ["indication_id array"]
}
```

### Traditional System Entity

```json
{
  "system_code": "string (AYU|TCM|UNA|SID|KAM|JAM)",
  "system_name": "string",
  "region": "string",
  "description": "string",
  "plant_count": "integer",
  "characteristics": {
    "diagnostic_approach": "string",
    "classification_system": "string",
    "preparation_methods": ["string array"]
  }
}
```

### Compound Entity

```json
{
  "compound_id": "string (TMMC internal, e.g., TMMC_C_001234)",
  "compound_name": "string",
  "identifiers": {
    "pubchem_cid": "integer",
    "chembl_id": "string",
    "cas_number": "string"
  },
  "structure": {
    "smiles": "string",
    "inchi_key": "string",
    "molecular_formula": "string",
    "molecular_weight": "float"
  },
  "source_plants": ["plant_id array"],
  "traditional_systems": ["system_code array"],
  "target_proteins": ["target_id array"]
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
  "source_systems": ["system_code array"],
  "pathway_associations": ["pathway_id array"]
}
```

### Therapeutic Indication Entity

```json
{
  "indication_id": "string",
  "indication_name": "string",
  "mesh_id": "string",
  "traditional_terms": [
    {
      "system": "string",
      "term": "string"
    }
  ],
  "associated_plants": ["plant_id array"],
  "associated_targets": ["target_id array"]
}
```

---

## Traditional Systems Covered

### System Profiles

| Code | System | Region | Coverage |
|------|--------|--------|----------|
| AYU | Ayurveda | India | Comprehensive |
| TCM | Traditional Chinese Medicine | China | Major herbs |
| UNA | Unani | Middle East/India | Selected |
| SID | Siddha | South India | Selected |
| KAM | Kampo | Japan | Selected |
| JAM | Jamu | Indonesia | Limited |

### System-Specific Fields

```json
{
  "ayurveda": {
    "rasa": ["sweet", "sour", "salty", "pungent", "bitter", "astringent"],
    "guna": ["heavy", "light", "hot", "cold"],
    "vipaka": "string",
    "virya": "string",
    "dosha_effect": ["vata", "pitta", "kapha"]
  },
  "tcm": {
    "nature": "cold|cool|neutral|warm|hot",
    "flavor": ["sweet", "sour", "bitter", "pungent", "salty"],
    "meridian_tropism": ["lung", "heart", "liver", "spleen", "kidney"]
  }
}
```

---

## Cross-System Analysis Features

### Shared Entity Detection

```json
{
  "cross_system_plant": {
    "plant_id": "string",
    "botanical_name": "string",
    "systems_present": ["AYU", "TCM", "KAM"],
    "convergent_uses": [
      {
        "indication": "string",
        "systems_agreeing": ["AYU", "TCM"],
        "confidence": "high|medium|low"
      }
    ]
  }
}
```

### Cross-System Queries

| Query Type | Description | Example |
|------------|-------------|---------|
| Plant overlap | Species in multiple systems | Curcuma longa in AYU, TCM |
| Compound sharing | Chemicals linking traditions | Curcumin across systems |
| Use convergence | Similar indications across cultures | Anti-inflammatory |
| Target overlap | Shared molecular mechanisms | COX-2 inhibition |

---

## Cross-References

### Identifier Mappings

| TM-MC Field | External Database | Notes |
|-------------|-------------------|-------|
| plant_id | MPNS | Medicinal Plant Names |
| compound_id | PubChem CID | Via structure |
| target_id | UniProt | Protein accession |
| indication_id | MESH | Medical subject headings |
| botanical_name | The Plant List | Taxonomic reference |

---

## Sample Data

### Sample Multi-System Plant

| Field | Value |
|-------|-------|
| plant_id | TMMC_P_000001 |
| botanical_name | Curcuma longa |
| family | Zingiberaceae |
| ayurveda_name | Haridra |
| tcm_name | Jiang Huang |
| systems_present | AYU, TCM, UNA |
| compounds | 58 |
| shared_uses | Anti-inflammatory, digestive |

### Sample Cross-System Compound

| Field | Value |
|-------|-------|
| compound_id | TMMC_C_000001 |
| compound_name | Curcumin |
| source_systems | AYU, TCM, KAM |
| source_plants | 3 |
| targets | 45 |

---

## Unique Value Proposition

### Cross-Cultural Validation

| Analysis Type | Insight |
|---------------|---------|
| Convergent uses | Multiple traditions agree on efficacy |
| Molecular basis | Shared targets explain traditional uses |
| Discovery opportunity | Novel uses from one system may apply to others |

### Comparative Analysis

```
System A (Ayurveda)         System B (TCM)
        |                        |
        v                        v
    Plant X                  Plant X
        |                        |
        v                        v
  Compound Y              Compound Y
        |                        |
        +------> Shared Target <------+
                      |
                      v
           Validated Mechanism
```

---

## Integration Recommendations

### Complementary Databases

| Database | Purpose |
|----------|---------|
| IMPPAT | Deep Ayurveda coverage |
| BATMAN-TCM | Deep TCM coverage |
| KampoDB | Deep Kampo coverage |
| HIT 2.0 | Experimental validation |

### Multi-System Workflow

```
1. Query TM-MC for cross-system plants
2. Identify convergent therapeutic uses
3. Extract shared compounds
4. Map to molecular targets
5. Validate with system-specific databases
6. Confirm with HIT 2.0 experimental data
```

---

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

---

## Limitations

1. Coverage varies by traditional system
2. Not all plants have molecular data
3. Therapeutic term standardization challenges
4. Updates depend on literature curation
5. Some system-specific nuances may be lost

---

## Glossary

| Term | Definition |
|------|------------|
| Convergent Use | Same therapeutic application across systems |
| Cross-system | Involving multiple traditional medicine systems |
| Molecular Correlate | Modern molecular explanation for traditional use |
| Traditional System | Codified medical tradition with diagnostic/therapeutic framework |
| Ethnobotanical | Related to traditional plant use knowledge |

---

## References

1. TM-MC Official Database
2. Official Website: http://tm-mc.org/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
