---
id: schema-tcmbank
title: "TCMBank Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: active
tags: [schema, tcm, compound-library, drug-discovery]
---

# TCMBank Database Schema

**Document ID:** TCMBANK-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** TCMBank

---

## TL;DR

TCMBank is a comprehensive TCM compound library designed for drug discovery. Contains 9,000+ herbs, 60,000+ compounds with structures, 75,000+ prescriptions, 15,000+ targets, and 8,000+ diseases. Emphasizes compound-centric data with drug-likeness assessment and ADMET predictions for identifying leads from natural sources.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **TCM Herbs** | 9,000+ |
| **Compounds** | 60,000+ |
| **Prescriptions/Formulas** | 75,000+ |
| **Target Proteins** | 15,000+ |
| **Diseases** | 8,000+ |
| **Compound-Target Pairs** | 250,000+ |

---

## Data Model

### Entity Relationships

```
Herbs (9,000+) -----> Formulas (75,000+)
     |
     v
Compounds (60,000+) -----> Drug-likeness Assessment
     |
     v
Targets (15,000+) -----> Diseases (8,000+)
```

### Compound Entity (Primary Focus)

```json
{
  "compound_id": "string (TCMBank internal, e.g., TCM001234)",
  "compound_name": "string",
  "compound_name_chinese": "string (if applicable)",
  "identifiers": {
    "pubchem_cid": "integer",
    "chembl_id": "string",
    "cas_number": "string"
  },
  "structure": {
    "smiles": "string",
    "inchi": "string",
    "inchi_key": "string (27 characters)"
  },
  "physicochemical_properties": {
    "molecular_formula": "string",
    "molecular_weight": "float (Da)",
    "logP": "float",
    "tpsa": "float (topological polar surface area)",
    "hydrogen_bond_donors": "integer",
    "hydrogen_bond_acceptors": "integer",
    "rotatable_bonds": "integer",
    "heavy_atom_count": "integer",
    "ring_count": "integer"
  },
  "drug_likeness": {
    "lipinski_violations": "integer (0-4)",
    "veber_compliant": "boolean",
    "pains_alerts": "integer",
    "lead_like": "boolean"
  },
  "admet_predictions": {
    "absorption": {
      "gi_absorption": "High|Low",
      "pgp_substrate": "Yes|No"
    },
    "distribution": {
      "bbb_permeant": "Yes|No"
    },
    "metabolism": {
      "cyp_inhibitors": ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
    },
    "toxicity": {
      "ames_mutagenicity": "Yes|No",
      "hepatotoxicity_risk": "High|Medium|Low"
    }
  },
  "source_herbs": ["herb_id array"],
  "predicted_targets": ["target_id array"]
}
```

### Herb Entity

```json
{
  "herb_id": "string (TCMBank internal, e.g., HERB001234)",
  "herb_name_chinese": "string",
  "herb_name_pinyin": "string",
  "herb_name_english": "string",
  "herb_name_latin": "string (botanical name)",
  "compound_count": "integer",
  "associated_formulas": ["formula_id array"]
}
```

### Formula Entity

```json
{
  "formula_id": "string (TCMBank internal)",
  "formula_name_chinese": "string",
  "formula_name_pinyin": "string",
  "classical_source": "string",
  "herb_composition": ["herb_id array"],
  "therapeutic_category": "string"
}
```

### Target Entity

```json
{
  "target_id": "string (UniProt or Gene Symbol)",
  "gene_symbol": "string",
  "gene_name": "string",
  "uniprot_id": "string",
  "protein_class": "string",
  "source_compounds": ["compound_id array"],
  "disease_associations": ["disease_id array"]
}
```

### Disease Entity

```json
{
  "disease_id": "string (internal or external)",
  "disease_name": "string",
  "mesh_id": "string",
  "associated_targets": ["target_id array"],
  "associated_herbs": ["herb_id array"]
}
```

---

## Drug-Likeness Filters

### Lipinski Rule of Five

| Rule | Criteria | Interpretation |
|------|----------|----------------|
| MW | <= 500 Da | Molecular weight |
| logP | <= 5 | Lipophilicity |
| HBD | <= 5 | Hydrogen bond donors |
| HBA | <= 10 | Hydrogen bond acceptors |

Violations: 0-4, with 0-1 typically acceptable for oral drugs.

### Veber Filter

| Rule | Criteria |
|------|----------|
| Rotatable Bonds | <= 10 |
| TPSA | <= 140 A^2 |

### PAINS (Pan-Assay Interference Compounds)

Compounds flagged for structural alerts that cause assay interference:
- Quinones
- Catechols
- Rhodanines
- Other reactive groups

### Lead-Like Properties

| Property | Criteria |
|----------|----------|
| MW | 250-350 Da |
| logP | <= 3.5 |
| HBD | <= 3 |
| HBA | <= 6 |
| Rotatable Bonds | <= 4 |

---

## Structure Files

| Format | Extension | Use Case |
|--------|-----------|----------|
| SDF | .sdf | Standard 2D/3D with properties |
| MOL | .mol | Single molecule format |
| SMILES | text | Linear notation |

---

## Cross-References

### Identifier Mappings

| TCMBank Field | External Database | Notes |
|---------------|-------------------|-------|
| compound_id | PubChem CID | Via structure |
| compound_id | ChEMBL ID | Via InChIKey |
| target_id | UniProt | Protein accession |
| disease_id | MESH | Disease concepts |

---

## Sample Data

### Sample Compound Record

| Field | Value |
|-------|-------|
| compound_id | TCM001234 |
| compound_name | Quercetin |
| pubchem_cid | 5280343 |
| molecular_weight | 302.24 |
| logP | 1.54 |
| lipinski_violations | 0 |
| pains_alerts | 0 |

---

## Integration Recommendations

### Virtual Screening Workflow

```
1. Download TCMBank compound library (SDF)
2. Filter by drug-likeness (Lipinski, PAINS)
3. Dock against target of interest
4. Identify hits from natural sources
5. Validate via HIT 2.0 experimental data
```

### Complementary Databases

| Database | Purpose |
|----------|---------|
| BATMAN-TCM | Target prediction validation |
| TCMSID | Additional ADME properties |
| HIT 2.0 | Experimental target validation |
| ChEMBL | Bioactivity comparison |

---

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

---

## Limitations

1. Focus on Chinese publications may limit coverage
2. ADMET predictions are computational
3. Regular updates needed for new literature
4. Some structures may require verification

---

## Glossary

| Term | Definition |
|------|------------|
| Drug-likeness | Properties suggesting oral drug potential |
| PAINS | Pan-Assay Interference Compounds |
| Lead-like | Properties suitable for lead optimization |
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity |
| TPSA | Topological Polar Surface Area |

---

## References

1. Wang R, et al. (2022). TCMBank: A comprehensive database for Traditional Chinese Medicine. Briefings in Bioinformatics.
2. Official Website: http://tcmbank.cn/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
