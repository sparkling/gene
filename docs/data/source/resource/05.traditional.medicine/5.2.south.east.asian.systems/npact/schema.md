---
id: schema-npact
title: "NPACT Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: active
tags: [schema, natural-products, anti-cancer, bioactivity, drug-discovery]
---

# NPACT Database Schema

**Document ID:** NPACT-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** NPACT

---

## TL;DR

NPACT (Naturally occurring Plant-based Anti-cancer Compound-activity-Target database) focuses on experimentally validated anti-cancer compounds from plant sources. Contains 1,574 anti-cancer compounds, 1,000+ plant sources, 150+ cancer cell lines, 200+ targets, 3,000+ literature references, and 5,000+ activity data points with IC50/EC50 values. Unlike prediction-based databases, NPACT prioritizes experimentally confirmed bioactivity data.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **Anti-cancer Compounds** | 1,574 |
| **Plant Sources** | 1,000+ |
| **Cancer Cell Lines** | 150+ |
| **Target Proteins** | 200+ |
| **Literature References** | 3,000+ |
| **Activity Data Points** | 5,000+ |

---

## Data Model

### Entity Relationships

```
Plant Sources (1,000+)
        |
        v
Anti-cancer Compounds (1,574) ---> Chemical Properties
        |                          - SMILES, MW, logP
        |                          - Drug-likeness
        v
Cancer Cell Lines (150+) <-----> Activity Data (IC50, etc.)
        |
        v
Target Proteins (200+)
```

### Compound Entity (Primary Focus)

```json
{
  "npact_id": "string (internal, e.g., NPACT001234)",
  "compound_name": "string",
  "identifiers": {
    "pubchem_cid": "integer",
    "cas_number": "string",
    "chembl_id": "string"
  },
  "structure": {
    "smiles": "string",
    "inchi": "string",
    "molecular_formula": "string",
    "molecular_weight": "float"
  },
  "physicochemical_properties": {
    "logP": "float",
    "tpsa": "float",
    "hydrogen_bond_donors": "integer",
    "hydrogen_bond_acceptors": "integer",
    "rotatable_bonds": "integer"
  },
  "drug_likeness": {
    "lipinski_violations": "integer",
    "bioavailability_score": "float"
  },
  "source_plants": ["plant_id array"],
  "activity_data": ["activity_record array"],
  "target_proteins": ["target_id array"]
}
```

### Activity Data Entity (Key Feature)

```json
{
  "activity_id": "string",
  "compound_id": "string",
  "activity_type": "string (IC50|EC50|GI50|LC50)",
  "activity_value": "float",
  "activity_unit": "string (nM|uM|mM|ug/mL)",
  "activity_relation": "string (=|<|>|~)",
  "cell_line": {
    "name": "string (e.g., MCF-7)",
    "cancer_type": "string (e.g., Breast)",
    "species": "string (Human|Mouse)",
    "tissue": "string"
  },
  "assay_details": {
    "assay_type": "string",
    "exposure_time": "string",
    "endpoint": "string"
  },
  "reference": {
    "pubmed_id": "integer",
    "doi": "string",
    "authors": "string",
    "year": "integer"
  }
}
```

### Plant Source Entity

```json
{
  "plant_id": "string",
  "botanical_name": "string (Latin binomial)",
  "family": "string",
  "common_names": ["string array"],
  "traditional_systems": ["Ayurveda", "TCM", "Western Herbal"],
  "plant_parts_used": ["root", "leaf", "bark", "seed"],
  "compounds_isolated": ["compound_id array"]
}
```

### Cancer Cell Line Entity

```json
{
  "cell_line_name": "string (e.g., MCF-7)",
  "cancer_type": "string",
  "tissue_origin": "string",
  "species": "string",
  "atcc_id": "string",
  "characteristics": ["string array"],
  "tested_compounds": "integer"
}
```

### Target Protein Entity

```json
{
  "target_id": "string",
  "gene_symbol": "string",
  "gene_name": "string",
  "uniprot_id": "string",
  "target_class": "string (Kinase|GPCR|Nuclear Receptor|etc.)",
  "cancer_relevance": "string",
  "compound_associations": ["compound_id array"]
}
```

---

## Activity Data Types

### Activity Measurements

| Type | Full Name | Description |
|------|-----------|-------------|
| IC50 | Half-maximal Inhibitory Concentration | Concentration for 50% inhibition |
| EC50 | Half-maximal Effective Concentration | Concentration for 50% effect |
| GI50 | Growth Inhibition 50 | Concentration for 50% growth inhibition |
| LC50 | Lethal Concentration 50 | Concentration killing 50% of cells |
| CC50 | Cytotoxic Concentration 50 | Concentration for 50% cytotoxicity |

### Unit Conversions

| Unit | Abbreviation | Conversion to nM |
|------|--------------|------------------|
| Nanomolar | nM | 1 |
| Micromolar | uM | 1,000 |
| Millimolar | mM | 1,000,000 |
| Microgram/mL | ug/mL | MW-dependent |

---

## Cancer Type Coverage

| Cancer Type | Compound Count | Cell Lines |
|-------------|----------------|------------|
| Breast | 400+ | MCF-7, MDA-MB-231, T47D |
| Lung | 300+ | A549, H460, H1299 |
| Colon | 250+ | HT-29, HCT-116, SW480 |
| Leukemia | 200+ | HL-60, K562, Jurkat |
| Liver | 150+ | HepG2, Hep3B, SMMC-7721 |
| Prostate | 150+ | PC-3, DU145, LNCaP |
| Cervical | 100+ | HeLa, SiHa, CaSki |
| Ovarian | 80+ | SKOV-3, A2780, OVCAR-3 |

---

## Cross-References

### Identifier Mappings

| NPACT Field | External Database | Notes |
|-------------|-------------------|-------|
| pubchem_cid | PubChem | Compound structure |
| chembl_id | ChEMBL | Bioactivity data |
| target_id | UniProt | Protein accession |
| pubmed_id | PubMed | Literature reference |
| cell_line | ATCC | Cell line repository |

---

## Sample Data

### Sample Compound Record

| Field | Value |
|-------|-------|
| npact_id | NPACT001234 |
| compound_name | Curcumin |
| pubchem_cid | 969516 |
| molecular_weight | 368.38 |
| source_plants | Curcuma longa |
| activity_count | 156 |

### Sample Activity Record

| Field | Value |
|-------|-------|
| compound | Curcumin |
| cell_line | MCF-7 |
| cancer_type | Breast |
| activity_type | IC50 |
| activity_value | 15.3 |
| activity_unit | uM |
| pubmed_id | 12345678 |

---

## Data Quality Features

### Experimental Validation

| Feature | Description |
|---------|-------------|
| All compounds | Tested activity (no predictions) |
| Quantitative data | IC50/EC50 values, not just +/- |
| Literature-linked | PubMed citations for validation |
| Cell line specified | Standardized cancer model names |

### Quality Indicators

- Primary literature sources
- Standard cell line nomenclature
- Reproducible assay conditions
- Activity value precision

---

## Integration Recommendations

### Complementary Databases

| Database | Purpose |
|----------|---------|
| IMPPAT | Broader Ayurvedic coverage |
| ChEMBL | Additional bioactivity data |
| NCI-60 | Standardized screening data |
| HIT 2.0 | Additional validated targets |

### Drug Discovery Workflow

```
1. Query NPACT for cancer type of interest
2. Filter by activity threshold (e.g., IC50 < 10 uM)
3. Identify common plant sources
4. Cross-reference with IMPPAT for additional compounds
5. Validate lead structures against ChEMBL
6. Perform ADMET filtering (via TCMSID)
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

1. Focus on anti-cancer limits scope
2. Activity data from diverse assays (heterogeneous)
3. No standardized target prediction method
4. Updates may lag behind literature
5. English literature bias

---

## Glossary

| Term | Definition |
|------|------------|
| IC50 | Half-maximal Inhibitory Concentration |
| EC50 | Half-maximal Effective Concentration |
| GI50 | Concentration for 50% growth inhibition |
| Cell Line | Immortalized cancer cells for testing |
| Bioactivity | Measurable biological effect |
| Lead Compound | Promising candidate for drug development |

---

## References

1. Mangal M, et al. (2013). NPACT: Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target database. Nucleic Acids Research.
2. Official Website: http://crdd.osdd.net/raghava/npact/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
