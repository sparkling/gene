# IMPPAT 2.0 Database Schema

**Document ID:** IMPPAT-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** IMPPAT 2.0 (June 2022)

---

## TL;DR

IMPPAT 2.0 (Indian Medicinal Plants, Phytochemistry And Therapeutics) is the largest digital database on phytochemicals of Indian medicinal plants with 4,010 plants, 17,967 phytochemicals, and 27,365 predicted protein-target interactions. Data is organized by plant parts with comprehensive ADMET properties, drug-likeness scores, and 1,875 molecular descriptors per compound. Licensed CC BY-NC 4.0.

---

## Database Statistics (Version 2.0)

| Entity | Count |
|--------|-------|
| **Indian Medicinal Plants** | 4,010 |
| **Phytochemicals** | 17,967 |
| **Therapeutic Uses** | 1,095 |
| **Plant-Part-Phytochemical Associations** | 189,386 |
| **Plant-Part-Therapeutic Use Associations** | 89,733 |
| **Predicted Human Target Proteins** | 5,042 |
| **Phytochemical-Target Interactions** | 27,365 |
| **2D/3D Molecular Descriptors per Compound** | 1,875 |
| **Drug-like Phytochemicals** | 1,335 |

---

## Data Access

### Primary URL
```
https://cb.imsc.res.in/imppat/
```

### API Availability
**No public REST API**. Data access via:
- Web interface queries
- TSV export from search results
- Structure downloads (SDF, MOL, MOL2, PDB, PDBQT)

### GitHub Repository
```
https://github.com/asamallab/IMPPAT2
```
Contains analysis scripts (not raw data).

---

## Data Model

### Core Entity Relationships

```
┌─────────────────────────────────────────────────────────────────┐
│                    INDIAN MEDICINAL PLANTS                       │
│  4,010 species from Ayurveda, Siddha, Unani, Homeopathy, etc.   │
│  Includes: taxonomy, conservation status, vernacular names       │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                        PLANT PARTS                               │
│  Granular associations at plant-part level                       │
│  Examples: roots, stems, leaves, bark, flowers, seeds, fruits    │
└──────────┬──────────────────────────────────┬───────────────────┘
           │                                  │
           ▼                                  ▼
┌──────────────────────────┐    ┌──────────────────────────────────┐
│     PHYTOCHEMICALS       │    │       THERAPEUTIC USES            │
│   17,967 compounds       │    │    1,095 traditional uses         │
│   189,386 associations   │    │    89,733 associations            │
└──────────┬───────────────┘    └──────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────┐
│               PREDICTED HUMAN TARGET PROTEINS                    │
│  5,042 proteins from STITCH database                            │
│  27,365 high-confidence interactions (score >= 700)             │
└─────────────────────────────────────────────────────────────────┘
```

---

## Entity Schemas

### Plant Entity

```json
{
  "imppat_plant_id": "IMPPAT_PLANT_XXXXX",
  "botanical_name": "string",
  "family": "string",
  "synonyms": ["string array"],
  "vernacular_names": {
    "hindi": "string",
    "tamil": "string",
    "telugu": "string",
    "kannada": "string",
    "malayalam": "string",
    "bengali": "string",
    "gujarati": "string",
    "marathi": "string",
    "punjabi": "string",
    "oriya": "string"
  },
  "taxonomic_classification": {
    "kingdom": "string",
    "phylum": "string",
    "class": "string",
    "order": "string",
    "family": "string",
    "genus": "string",
    "species": "string"
  },
  "traditional_medicine_systems": ["Ayurveda", "Siddha", "Unani", "Homeopathy", "Sowa-Rigpa"],
  "iucn_conservation_status": "string",
  "external_links": {
    "the_plant_list": "url",
    "tropicos": "url",
    "mpns": "url",
    "ipni": "url",
    "pow": "url",
    "wfo": "url"
  }
}
```

### Phytochemical Entity

Six data tabs per compound on dedicated pages:

#### Tab 1: Summary
```json
{
  "imppat_id": "IMPPAT_CHEM_XXXXX",
  "compound_name": "string",
  "iupac_name": "string",
  "smiles": "string",
  "inchi": "string",
  "inchi_key": "string",
  "deepsmiles": "string (for ML applications)",
  "molecular_formula": "string",
  "molecular_weight": "float",
  "source_plants": ["plant_id array"],
  "source_plant_parts": ["string array"]
}
```

#### Tab 2: Physicochemical Properties
```json
{
  "molecular_weight": "float (Da)",
  "logP": "float (partition coefficient)",
  "topological_polar_surface_area": "float (Å²)",
  "hydrogen_bond_donors": "integer",
  "hydrogen_bond_acceptors": "integer",
  "rotatable_bonds": "integer",
  "heavy_atom_count": "integer",
  "ring_count": "integer",
  "aromatic_ring_count": "integer"
}
```

#### Tab 3: Drug-Likeness
```json
{
  "lipinski_rule_of_five": {
    "passes": "boolean",
    "violations": "integer",
    "mw_ok": "boolean (<= 500)",
    "logp_ok": "boolean (<= 5)",
    "hbd_ok": "boolean (<= 5)",
    "hba_ok": "boolean (<= 10)"
  },
  "ghose_filter": {
    "passes": "boolean",
    "violations": "integer"
  },
  "veber_filter": {
    "passes": "boolean",
    "rotatable_bonds_ok": "boolean (<= 10)",
    "tpsa_ok": "boolean (<= 140)"
  },
  "egan_filter": {
    "passes": "boolean"
  },
  "pfizer_3_75_rule": {
    "passes": "boolean",
    "logp_ok": "boolean",
    "tpsa_ok": "boolean"
  },
  "gsk_4_400_rule": {
    "passes": "boolean",
    "mw_ok": "boolean",
    "logp_ok": "boolean"
  },
  "np_likeness_score": "float (natural product likeness)",
  "weighted_qed_score": "float (quantitative estimate of drug-likeness)"
}
```

#### Tab 4: ADMET Properties (from SwissADME)
```json
{
  "absorption": {
    "gi_absorption": "High/Low",
    "pgp_substrate": "Yes/No",
    "water_solubility": "string (ESOL class)",
    "caco2_permeability": "float"
  },
  "distribution": {
    "bbb_permeant": "Yes/No",
    "vdss": "float (volume of distribution)"
  },
  "metabolism": {
    "cyp1a2_inhibitor": "Yes/No",
    "cyp2c19_inhibitor": "Yes/No",
    "cyp2c9_inhibitor": "Yes/No",
    "cyp2d6_inhibitor": "Yes/No",
    "cyp3a4_inhibitor": "Yes/No"
  },
  "excretion": {
    "clearance": "float"
  },
  "toxicity": {
    "ames_mutagenicity": "Yes/No",
    "hepatotoxicity": "Yes/No",
    "herg_inhibition": "Yes/No",
    "skin_sensitization": "Yes/No"
  }
}
```

**Note:** ADMET predictions unavailable for 493 compounds due to SMILES length restrictions in SwissADME.

#### Tab 5: Molecular Descriptors
```json
{
  "descriptor_count": 1875,
  "descriptor_types": {
    "2d_descriptors": ["list of 2D molecular descriptors"],
    "3d_descriptors": ["list of 3D molecular descriptors"]
  },
  "note": "1,875 2D and 3D chemical descriptors computed per compound"
}
```

#### Tab 6: Predicted Human Target Proteins
```json
{
  "targets": [
    {
      "protein_name": "string",
      "gene_symbol": "string",
      "hgnc_id": "string (HUGO Gene Nomenclature Committee)",
      "combined_score": "integer (>= 700)",
      "source": "STITCH"
    }
  ],
  "interaction_count": "integer",
  "high_confidence_threshold": 700
}
```

---

## Structural Classification

### Molecular Scaffolds

Three levels of scaffold abstraction (Murcko decomposition):

| Level | Code | Description |
|-------|------|-------------|
| Generic/Node/Bond | G/N/B | Full Murcko scaffold |
| Generic/Node | G/N | Atom types removed |
| Graph | Graph | Only topology preserved |

### Chemical Classification (NPClassifier)

```json
{
  "biosynthetic_pathway": "string (natural product classification)",
  "superclass": "string (ClassyFire prediction)",
  "class": "string",
  "subclass": "string",
  "functional_groups": ["string array"]
}
```

---

## Download Formats

### Structure Files

| Format | Extension | Use Case |
|--------|-----------|----------|
| SDF | .sdf | Standard 2D/3D |
| MOL | .mol | Single molecule |
| MOL2 | .mol2 | With atom types |
| PDB | .pdb | For docking |
| PDBQT | .pdbqt | AutoDock ready |

### Data Export

| Format | Content |
|--------|---------|
| TSV | Plant-phytochemical associations |
| TSV | Plant-therapeutic use associations |
| Excel | Supplementary tables from publication |

---

## Cross-References

### Identifier Mappings

| IMPPAT Field | External Database | Link Pattern |
|--------------|-------------------|--------------|
| Plant Name | The Plant List | Direct URL |
| Plant Name | Tropicos | Direct URL |
| Plant Name | MPNS | Direct URL |
| Plant Name | IPNI | Direct URL |
| Plant Name | POW (Plants of the World) | Direct URL |
| Plant Name | WFO (World Flora Online) | Direct URL |
| InChIKey | PubChem | UniChem mapping |
| Target Protein | HGNC | Gene nomenclature |
| Target Interactions | STITCH | Score >= 700 |

### UniChem Integration

IMPPAT provides external links to standard chemical databases via UniChem service.

---

## Data Quality

### FAIR Compliance

- **Findable**: Unique IMPPAT identifiers for plants and compounds
- **Accessible**: Web interface with export capabilities
- **Interoperable**: Standard formats (InChI, SMILES, SDF)
- **Reusable**: CC BY-NC 4.0 license, documented schema

### Target Prediction Quality

| Metric | Value |
|--------|-------|
| Source | STITCH database |
| Confidence Threshold | Combined score >= 700 |
| Coverage | 1,294 phytochemicals have predictions |
| Protein Coverage | 5,042 unique human proteins |
| Total Interactions | 27,365 high-confidence pairs |

---

## Web Interface Queries

### Search Options

1. **Plant Search**
   - By botanical name
   - By vernacular name (10 Indian languages)
   - By family
   - By therapeutic use

2. **Phytochemical Search**
   - By compound name
   - By SMILES/InChI
   - By molecular formula
   - By physicochemical properties
   - By drug-likeness filters

3. **Advanced Search**
   - Substructure search
   - Similarity search
   - Property range filters
   - Target-based search

---

## Sample Data Records

### Sample Plant Record

**Ashwagandha (Withania somnifera)**

| Field | Value |
|-------|-------|
| IMPPAT ID | IMPPAT_PLANT_03985 |
| Botanical Name | Withania somnifera |
| Family | Solanaceae |
| Traditional Systems | Ayurveda, Unani |
| IUCN Status | Not Evaluated |
| Phytochemicals | 127 compounds |
| Therapeutic Uses | Adaptogen, immunomodulator, anti-stress |

### Sample Phytochemical Record

**Withaferin A**

| Field | Value |
|-------|-------|
| IMPPAT ID | IMPPAT_CHEM_15234 |
| Name | Withaferin A |
| Formula | C28H38O6 |
| Molecular Weight | 470.60 |
| logP | 3.75 |
| TPSA | 96.36 Å² |
| Lipinski | Passes (0 violations) |
| GI Absorption | High |
| BBB Permeant | Yes |
| Predicted Targets | 23 proteins |

---

## Integration Recommendations

### Priority Integration Steps

1. **Export plant-phytochemical associations** via web TSV export
2. **Download structure files** (SDF) for compound library
3. **Map to PubChem** via InChIKey for standardization
4. **Link targets** via HGNC to UniProt IDs

### Complementary Databases

| Database | Purpose |
|----------|---------|
| NPACT | Validated cancer targets |
| NPASS | Quantitative activity data |
| CMAUP | Pathway annotations |
| STRING | Extended target networks |
| ChEMBL | Bioactivity validation |

### Data Harmonization

```python
# Example: Map IMPPAT to UniProt via HGNC
import pandas as pd

def map_imppat_to_uniprot(hgnc_id, hgnc_mapping):
    """Map IMPPAT HGNC IDs to UniProt accessions."""
    return hgnc_mapping.get(hgnc_id, {}).get('uniprot_id')

# IMPPAT uses HGNC for target annotation
# UniProt provides HGNC -> UniProt mapping in ID mapping service
```

---

## License

**CC BY-NC 4.0** (Creative Commons Attribution-NonCommercial 4.0 International)

- **Permits**: Non-commercial use, sharing, adaptation
- **Requires**: Attribution to IMPPAT 2.0 and citation
- **Restricts**: Commercial use without separate agreement
- **Maintainer**: Institute of Mathematical Sciences, Chennai

---

## Limitations

1. **No REST API**: Manual export or web scraping required
2. **Commercial restriction**: CC BY-NC 4.0 limits commercial use
3. **ADMET gaps**: 493 compounds lack ADMET predictions
4. **Target predictions only**: No experimentally validated interactions
5. **Confidence threshold**: Only high-confidence (>=700) STITCH interactions included

---

## GitHub Repository Contents

The associated GitHub repository (asamallab/IMPPAT2) contains:

### Analysis Scripts

| Script | Purpose |
|--------|---------|
| ChemicalSimilarityNetwork.py | Tanimoto coefficient calculations |
| ChemicalStructureImages.py | SVG/PNG structure generation |
| DruglikenessProperties.py | Drug-likeness evaluation |
| MolecularProperties.py | Physicochemical computation |
| MolecularScaffolds.py | Scaffold analysis |
| MurckoScaffold.py | Modified RDKit scaffold code |

**Note**: Raw database files are NOT in the repository; access via web interface.

---

## References

1. Vivek-Ananth RP, et al. (2023). IMPPAT 2.0: An Enhanced and Expanded Phytochemical Atlas of Indian Medicinal Plants. ACS Omega, 8(9):8827-8845.

2. Mohanraj K, et al. (2018). IMPPAT: A curated database of Indian Medicinal Plants, Phytochemistry And Therapeutics. Scientific Reports, 8:4329.

3. GitHub Repository: https://github.com/asamallab/IMPPAT2

4. Official Website: https://cb.imsc.res.in/imppat/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation from publication and web interface analysis |
