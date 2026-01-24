---
id: schema-phytohub
title: "PhytoHub Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, phytochemicals, metabolites, biomarkers, mass-spectrometry]
---

# PhytoHub - Dietary Phytochemical Metabolome Schema

**Document ID:** SCHEMA-PHYTOHUB
**Version:** 1.5
**Source Version:** 2024

---

## TL;DR

PhytoHub catalogs dietary phytochemicals and their human/microbial metabolites with mass spectrometry data. The schema organizes parent compounds, their biotransformation products, and analytical reference data (MS/MS spectra) for metabolomics identification of dietary biomarkers.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Total Compounds | 1,800+ | Compound database |
| Parent Phytochemicals | 600+ | Food-derived |
| Human Metabolites | 400+ | Phase I/II products |
| Microbial Metabolites | 600+ | Gut bacteria products |
| Dietary Sources | 350+ | Food associations |
| MS/MS Spectra | 1,200+ | Reference library |

---

## Entity Relationship Overview

```
Parent_Compounds (1) ←→ (many) Metabolites
        ↓                          ↓
    Food sources             Transformation type

Compounds (1) ←→ (many) MS_Spectra
                            ↓
                   m/z, intensity, conditions

Compounds (1) ←→ (many) Food_Sources
                            ↓
                    Food name, concentration

All Compounds ←→ External_IDs
                     ↓
              PubChem, HMDB, ChEBI
```

---

## Core Tables/Entities

### compounds

**Description:** All phytochemicals and metabolites with structural data.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| phytohub_id | string | Yes | PhytoHub ID (PHUB+digits) |
| name | string | Yes | Compound name |
| compound_type | string | Yes | parent, phase1, phase2, microbial |
| parent_id | string | No | Parent compound (for metabolites) |
| smiles | string | No | SMILES structure |
| inchi | string | No | InChI identifier |
| inchi_key | string | Yes | InChI Key |
| molecular_formula | string | No | Chemical formula |
| molecular_weight | decimal | No | Exact mass |
| chemical_class | string | No | Polyphenol, terpene, etc. |

### metabolites

**Description:** Biotransformation products of parent compounds.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| metabolite_id | string | Yes | Primary identifier |
| parent_id | string | Yes | Parent compound ID |
| metabolite_name | string | Yes | Metabolite name |
| transformation | string | Yes | Chemical transformation |
| enzyme | string | No | Metabolizing enzyme |
| location | string | No | liver, gut, etc. |
| biofluid | string | No | Detection matrix |

### transformations

**Description:** Types of metabolic transformations.

| Type | Description | Location |
|------|-------------|----------|
| Glucuronidation | UGT conjugation | Liver |
| Sulfation | SULT conjugation | Liver |
| Methylation | COMT conjugation | Liver |
| Oxidation | CYP450 | Liver |
| Reduction | Reductases | Liver/gut |
| Ring fission | Microbial | Gut |
| Dehydroxylation | Microbial | Gut |
| Demethylation | Microbial | Gut |

### ms_spectra

**Description:** Mass spectrometry reference spectra.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| spectrum_id | integer | Yes | Primary identifier |
| phytohub_id | string | Yes | Compound ID |
| precursor_mz | decimal | Yes | Precursor ion m/z |
| adduct | string | Yes | [M+H]+, [M-H]-, etc. |
| collision_energy | decimal | No | CE in eV |
| ionization_mode | string | Yes | positive, negative |
| instrument_type | string | No | Q-TOF, Orbitrap, etc. |
| peaks | array | Yes | m/z, intensity pairs |

### food_sources

**Description:** Links compounds to dietary sources.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| source_id | integer | Yes | Primary identifier |
| phytohub_id | string | Yes | Compound ID |
| food_name | string | Yes | Food name |
| food_group | string | No | Fruits, vegetables, etc. |
| typical_content | decimal | No | Amount in food |
| content_unit | string | No | mg/100g |

### external_ids

**Description:** Cross-references to other databases.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| phytohub_id | string | Yes | PhytoHub ID |
| database | string | Yes | PubChem, HMDB, ChEBI |
| external_id | string | Yes | External identifier |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | Web interface |
| MS Spectra | MSP, MGF formats |
| Structures | SMILES, InChI |
| Encoding | UTF-8 |

---

## Sample Record

```json
{
  "phytohub_id": "PHUB000001",
  "name": "Quercetin",
  "compound_type": "parent",
  "smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
  "inchi_key": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
  "molecular_weight": 302.236,
  "chemical_class": "Flavonoid",
  "metabolites": [
    {
      "metabolite_id": "PHUB000002",
      "name": "Quercetin-3-O-glucuronide",
      "transformation": "Glucuronidation",
      "location": "liver",
      "enzyme": "UGT1A9"
    },
    {
      "metabolite_id": "PHUB000010",
      "name": "3,4-dihydroxyphenylacetic acid",
      "transformation": "Ring fission",
      "location": "gut",
      "enzyme": "Microbial"
    }
  ],
  "food_sources": [
    {"food": "Onion", "content": 39.21},
    {"food": "Apple", "content": 4.42}
  ],
  "ms_spectrum": {
    "precursor_mz": 303.049,
    "adduct": "[M+H]+",
    "collision_energy": 30,
    "peaks": [[153.018, 100], [137.023, 45], [121.028, 30]]
  }
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| PhytoHub ID | Identifier format PHUB + 6 digits |
| Phase I | Oxidation/reduction reactions |
| Phase II | Conjugation reactions |
| Microbial metabolite | Gut bacteria product |
| Ring fission | Cleavage of aromatic ring |
| MS/MS | Tandem mass spectrometry |
| m/z | Mass-to-charge ratio |

---

## References

1. PhytoHub: http://phytohub.eu
2. Related: Phenol-Explorer (http://phenol-explorer.eu)
3. HMDB (https://hmdb.ca) for metabolite reference
