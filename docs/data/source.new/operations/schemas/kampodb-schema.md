---
id: schemas-kampodb-schema
title: "KampoDB Schema"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, kampo, traditional-medicine, natural-products, docking, rest-api]
---

**Parent:** [Schema Documentation](./_index.md)

# KampoDB Schema

**Document ID:** KAMPODB-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** KampoDB (2018, Active)

---

## TL;DR

KampoDB provides a comprehensive REST API with 298 Kampo formulas, 180 crude drugs, 3,002 compounds, and 62,906 proteins across 3,063,505 docking simulations. The 4-layer hierarchy (Formula -> Crude Drug -> Compound -> Protein) is fully traversable via clean JSON endpoints. Licensed under CC BY-SA 4.0, making it commercial-friendly.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **Kampo Formulas** | 298 |
| **Crude Drugs** | 180 |
| **Natural Compounds** | 3,002 |
| **Proteins/Genes** | 62,906 |
| **Docking Simulations** | 3,063,505 |
| **Known Target Proteins** | 460 |
| **Predicted Target Proteins** | 1,369 |

---

## REST API Documentation

### Base URL

```
https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/
```

All endpoints return JSON responses.

---

## API Endpoint Categories

### 1. List Endpoints

Retrieve all IDs and names for each entity type.

| Endpoint | Description | Response Size |
|----------|-------------|---------------|
| `/api/formula/` | All 298 formula IDs and names | ~481 entries |
| `/api/crude/` | All 180 crude drug IDs and names | 180 entries |
| `/api/compound/` | All 3,002 compound IDs and names | 3,002 entries |
| `/api/protein/` | All 62,906 protein IDs and names | 62,906 entries |

**Sample Response** (`/api/formula/`):
```json
[
  {"id": "ACS", "name": "Anchusan"},
  {"id": "KT", "name": "Kakkonto"},
  {"id": "YKS", "name": "Yokukansankan"},
  {"id": "TSS", "name": "Tokishakuyakusan"},
  {"id": "ZMT", "name": "Zokumeito"}
]
```

---

### 2. Information Endpoints

Get detailed information for a specific entity.

| Endpoint Pattern | Example | Description |
|-----------------|---------|-------------|
| `/api/formula/{id}/info` | `/api/formula/KT/info` | Formula details |
| `/api/crude/{id}/info` | `/api/crude/40/info` | Crude drug details |
| `/api/compound/{id}/info` | `/api/compound/969516/info` | Compound details |
| `/api/protein/{id}/info` | `/api/protein/2475/info` | Protein details |

**Sample Response** (`/api/formula/KT/info`):
```json
{
  "id": "KT",
  "phonetic": "かっこんとう",
  "name": "Kakkonto",
  "name_jp": "葛根湯",
  "document": "kakkonto/",
  "synonyms": null
}
```

**Sample Response** (`/api/crude/40/info`):
```json
{
  "id": 40,
  "phonetic": "けいひ",
  "name": "Cinnamon Bark",
  "name_jp": "桂皮",
  "origin": "Cinnamomum cassia",
  "synonyms": null
}
```

**Sample Response** (`/api/compound/969516/info`):
```json
{
  "id": 969516,
  "name": "Curcumin",
  "formula": "C21H20O6"
}
```

**Sample Response** (`/api/protein/2475/info`):
```json
{
  "id": 2475,
  "name": "MTOR",
  "aliases": "FRAP, FRAP1, FRAP2, RAFT1, RAPT1, SKS",
  "description": "mechanistic target of rapamycin kinase"
}
```

---

### 3. Relational Endpoints

Navigate the 4-layer hierarchy by querying relationships.

#### Formula Relations

| Endpoint | Description | Example |
|----------|-------------|---------|
| `/api/formula/{id}/crude` | Crude drugs in formula | `/api/formula/KT/crude` |
| `/api/formula/{id}/compound` | Compounds in formula | `/api/formula/KT/compound` |
| `/api/formula/{id}/protein` | Target proteins of formula | `/api/formula/KT/protein` |

**Sample Response** (`/api/formula/KT/crude`):
```json
[
  {"id": 40, "name": "Cinnamon Bark"},
  {"id": 62, "name": "Ephedra Herb"},
  {"id": 76, "name": "Ginger"},
  {"id": 78, "name": "Glycyrrhiza"},
  {"id": 90, "name": "Jujube"},
  {"id": 126, "name": "Peony Root"},
  {"id": 147, "name": "Pueraria Root"}
]
```

**Sample Response** (`/api/formula/KT/compound`):
```
Returns 400+ compounds including:
- Ephedrine, pseudoephedrine (alkaloids)
- Quercetin, kaempferol (flavonoids)
- Curcumin, gingerols, shogaols
- Genistein, daidzein (isoflavones)
```

**Sample Response** (`/api/formula/KT/protein`):
```
Returns 1,295 target proteins including:
- CYP450 family (CYP3A4, CYP2D6)
- Kinases (JAK1, JAK2, MAPK1)
- Apoptosis regulators (BCL2, BAX, CASP3)
- ABC transporters (ABCB1, ABCC1)
```

#### Crude Drug Relations

| Endpoint | Description |
|----------|-------------|
| `/api/crude/{id}/formula` | Formulas containing this crude drug |
| `/api/crude/{id}/compound` | Compounds in this crude drug |
| `/api/crude/{id}/protein` | Target proteins via compounds |

**Sample Response** (`/api/crude/40/compound`):
```json
[
  {"id": 637511, "name": "Cinnamaldehyde"},
  {"id": 444539, "name": "Cinnamic acid"},
  {"id": 323, "name": "Coumarin"},
  {"id": 72, "name": "3,4-Dihydroxybenzoic acid"},
  {"id": 637520, "name": "Cinnamtannin A3"}
]
```
(31 compounds for Cinnamon Bark)

#### Compound Relations

| Endpoint | Description |
|----------|-------------|
| `/api/compound/{id}/formula` | Formulas containing this compound |
| `/api/compound/{id}/crude` | Crude drugs containing this compound |
| `/api/compound/{id}/protein` | Target proteins of this compound |

**Sample Response** (`/api/compound/969516/crude`):
```json
[
  {"id": 70, "name": "Fresh Ginger"},
  {"id": 76, "name": "Ginger"},
  {"id": 144, "name": "Processed Ginger"},
  {"id": 174, "name": "Turmeric"}
]
```

**Sample Response** (`/api/compound/969516/protein`):
```
Returns 273 target proteins including:
- Apoptosis: BAX, BAK1, CASP3, CASP8, BCL2, TP53
- Transporters: ABCB1, ABCC1, ABCG2
- Signal transduction: JAK2, STAT3, MAPK1, AKT1
- Inflammation: TNF, IL6, IL1B, ICAM1
```

#### Protein Relations

| Endpoint | Description |
|----------|-------------|
| `/api/protein/{id}/formula` | Formulas targeting this protein |
| `/api/protein/{id}/crude` | Crude drugs targeting this protein |
| `/api/protein/{id}/compound` | Compounds targeting this protein |

---

### 4. Enrichment & Annotation Endpoints

Functional annotations for proteins.

| Endpoint Pattern | Analysis Type | Example |
|-----------------|---------------|---------|
| `/api/{category}/{id}/process` | GO Biological Processes | `/api/protein/2475/process` |
| `/api/{category}/{id}/function` | GO Molecular Functions | `/api/protein/2475/function` |
| `/api/{category}/{id}/pathway` | KEGG Pathways | `/api/protein/2475/pathway` |
| `/api/{category}/{id}/disease` | KEGG Diseases | `/api/protein/2475/disease` |

**Sample Response** (`/api/protein/2475/pathway`):
```json
[
  {"id": "hsa04152", "name": "AMPK signaling pathway"},
  {"id": "hsa04150", "name": "mTOR signaling pathway"},
  {"id": "hsa04151", "name": "PI3K-Akt signaling pathway"},
  {"id": "hsa04140", "name": "Autophagy"},
  {"id": "hsa05210", "name": "Colorectal cancer"},
  {"id": "hsa05014", "name": "Amyotrophic lateral sclerosis"},
  {"id": "hsa05010", "name": "Alzheimer disease"}
]
```
(50 KEGG pathways for MTOR)

---

### 5. Prediction Endpoints

Target prediction methods.

| Endpoint | Method | Notes |
|----------|--------|-------|
| `/api/compound/{id}/lbp` | Ligand-based prediction | Uses chemical similarity |
| `/api/compound/{id}/sbp` | Structure-based prediction | Up to 1000 results |
| `/api/protein/{id}/lbp` | Ligand-based prediction | |
| `/api/protein/{id}/sbp` | Structure-based prediction | |

---

### 6. Docking Simulation Endpoint

Molecular docking data between compounds and proteins.

**Endpoint:** `/api/docking/compound/{compound_id}/protein/{protein_id}`

**Example:** `/api/docking/compound/969516/protein/2475`

**Sample Response:**
```json
{
  "compound_id": 969516,
  "protein_id": 2475,
  "domain_id": 14339,
  "target_name": "NP_004949.1_holo_2118-2483",
  "affinity_kcal_mol": -7.8,
  "ligand_atoms": 29,
  "protein_atoms": 1401
}
```

**Interpretation:**
- Binding affinity: -7.8 kcal/mol indicates moderate favorable binding
- Domain covers residues 2118-2483
- Contains over 3 million docking results across database

---

## Data Model

### 4-Layer Hierarchical Structure

The KampoDB data model implements a strict hierarchical organization enabling traversal from traditional formula to molecular targets:

```
Layer 1: Kampo Medicine (Formula)
    ├─ ID: String code (e.g., "KT")
    ├─ Name: Romanized (e.g., "Kakkonto")
    ├─ Name_jp: Japanese kanji (e.g., "葛根湯")
    └─ Count: 298 total formulas

    │ Formula-Crude Drug Relationship (many-to-many)
    ▼
Layer 2: Crude Drugs (Drug)
    ├─ ID: Integer (0-179)
    ├─ Name: English common name (e.g., "Ephedra Herb")
    ├─ Name_jp: Japanese name (e.g., "麻黄")
    ├─ Origin: Botanical source (e.g., "Ephedra sinica")
    └─ Count: 180 total crude drugs

    │ Crude Drug-Compound Relationship (many-to-many)
    ▼
Layer 3: Compounds (Compound)
    ├─ ID: Integer (PubChem CID)
    ├─ Name: Chemical name (e.g., "Curcumin")
    ├─ Formula: Molecular formula (e.g., "C21H20O6")
    ├─ Source: External IDs (KNApSAcK, ChEMBL)
    └─ Count: 3,002 total compounds

    │ Compound-Protein Relationship (many-to-many, via docking predictions)
    ▼
Layer 4: Target Proteins (Protein)
    ├─ ID: Integer (NCBI Gene ID)
    ├─ Name: Gene symbol (e.g., "MTOR")
    ├─ Aliases: Alternative names (e.g., "FRAP, RAFT1")
    ├─ Description: Functional annotation
    └─ Count: 62,906 total proteins

Key Relationships:
- Formula -> Crude Drug: 1 formula contains multiple crude drugs
- Crude Drug -> Compound: 1 crude drug contains multiple compounds
- Compound -> Protein: 1 compound targets multiple proteins (via docking)
- All relationships are many-to-many, fully traversable
```

### Hierarchy Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                       KAMPO FORMULAS (298)                       │
│  Traditional Japanese medicine preparations                      │
│  Example: Kakkonto (KT), Yokukansankan (YKS)                    │
└─────────────────────┬───────────────────────────────────────────┘
                      │ contains (many-to-many)
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                       CRUDE DRUGS (180)                          │
│  Raw medicinal materials (herbs, minerals, animal products)      │
│  Example: Cinnamon Bark (40), Ephedra Herb (62), Ginger (76)    │
└─────────────────────┬───────────────────────────────────────────┘
                      │ contains (many-to-many)
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                       COMPOUNDS (3,002)                          │
│  Natural chemical compounds with PubChem identifiers             │
│  Example: Curcumin (969516), Ephedrine (9294), Gingerol         │
└─────────────────────┬───────────────────────────────────────────┘
                      │ targets (many-to-many, predicted via docking)
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                       PROTEINS (62,906)                          │
│  Human proteins/genes with NCBI Gene identifiers                 │
│  Example: MTOR (2475), CYP3A4, BCL2, STAT3                      │
└─────────────────────────────────────────────────────────────────┘
```

### Implementation Details

#### Hierarchy Navigation

Each layer provides complete relational endpoints for traversal:

| From | To | Endpoint | Returns |
|------|-----|----------|---------|
| Formula | Crude Drug | `/api/formula/{id}/crude` | Array of crude drugs |
| Formula | Compound | `/api/formula/{id}/compound` | Array of compounds |
| Formula | Protein | `/api/formula/{id}/protein` | Array of target proteins |
| Crude Drug | Formula | `/api/crude/{id}/formula` | Parent formulas |
| Crude Drug | Compound | `/api/crude/{id}/compound` | Array of compounds |
| Crude Drug | Protein | `/api/crude/{id}/protein` | Target proteins via compounds |
| Compound | Crude Drug | `/api/compound/{id}/crude` | Source crude drugs |
| Compound | Protein | `/api/compound/{id}/protein` | Array of target proteins |
| Compound | Formula | `/api/compound/{id}/formula` | Parent formulas |
| Protein | Compound | `/api/protein/{id}/compound` | Ligand compounds |
| Protein | Crude Drug | `/api/protein/{id}/crude` | Source crude drugs |
| Protein | Formula | `/api/protein/{id}/formula` | Parent formulas |

#### Hierarchy Constraints

- **Directionality**: One-way from formula down; reverse queries available
- **All-to-All Connectivity**: Every layer can reach every other layer
- **No Intermediate Skipping**: Must traverse via proper relationships (no direct formula→protein shortcuts)
- **Many-to-Many Cardinality**: Single formula may contain same crude drug multiple times (rare), all relationships are N:N

#### Data Integration Points

1. **Formula Layer**: TradMPD, STORK cross-references
2. **Crude Drug Layer**: Botanical taxonomy via IPNI IDs
3. **Compound Layer**: KNApSAcK C_IDs, PubChem CIDs for external lookups
4. **Protein Layer**: KEGG GENES, UniProt, pathway/disease annotations

### Entity Schemas

#### Formula Entity
```json
{
  "id": "string (formula code, e.g., 'KT')",
  "phonetic": "string (Japanese kana)",
  "name": "string (English name)",
  "name_jp": "string (Japanese kanji)",
  "document": "string (documentation path)",
  "synonyms": "string or null"
}
```

#### Crude Drug Entity
```json
{
  "id": "integer (0-179)",
  "phonetic": "string (Japanese kana)",
  "name": "string (English name)",
  "name_jp": "string (Japanese kanji)",
  "origin": "string (botanical/scientific source)",
  "synonyms": "string or null"
}
```

#### Compound Entity
```json
{
  "id": "integer (PubChem CID)",
  "name": "string (compound name)",
  "formula": "string (molecular formula)"
}
```

#### Protein Entity
```json
{
  "id": "integer (NCBI Gene ID)",
  "name": "string (gene symbol)",
  "aliases": "string (comma-separated aliases)",
  "description": "string (gene/protein description)"
}
```

#### Docking Result Entity
```json
{
  "compound_id": "integer (PubChem CID)",
  "protein_id": "integer (NCBI Gene ID)",
  "domain_id": "integer (protein domain ID)",
  "target_name": "string (domain identifier with residue range)",
  "affinity_kcal_mol": "float (binding energy in kcal/mol)",
  "ligand_atoms": "integer (atom count)",
  "protein_atoms": "integer (atom count)"
}
```

---

## Cross-References

### Identifier Mappings

| KampoDB ID Type | External Database | Notes |
|-----------------|-------------------|-------|
| Compound ID | PubChem CID | Direct mapping |
| Compound ID | KNApSAcK C_ID | Integration partner |
| Protein ID | NCBI Gene ID | Standard gene identifiers |
| Pathway ID | KEGG pathway | hsa prefix for human |
| Disease ID | KEGG disease | Disease annotations |

### External Links

- **KNApSAcK**: http://www.knapsackfamily.com/knapsack_core/top.php
- **PubChem**: https://pubchem.ncbi.nlm.nih.gov/
- **KEGG**: https://www.kegg.jp/
- **ChEMBL**: https://www.ebi.ac.uk/chembl/
- **UniProt**: https://www.uniprot.org/

---

## Usage Examples

### Get All Compounds in a Formula
```bash
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/KT/compound"
```

### Get Target Proteins for a Compound
```bash
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/compound/969516/protein"
```

### Get Docking Score
```bash
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/docking/compound/969516/protein/2475"
```

### Python Integration Example
```python
import requests

BASE_URL = "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api"

def get_formula_compounds(formula_id):
    """Get all compounds in a Kampo formula."""
    response = requests.get(f"{BASE_URL}/formula/{formula_id}/compound")
    return response.json()

def get_compound_targets(compound_id):
    """Get target proteins for a compound."""
    response = requests.get(f"{BASE_URL}/compound/{compound_id}/protein")
    return response.json()

def get_docking_score(compound_id, protein_id):
    """Get molecular docking results."""
    response = requests.get(f"{BASE_URL}/docking/compound/{compound_id}/protein/{protein_id}")
    return response.json()

# Example: Analyze Kakkonto formula
compounds = get_formula_compounds("KT")
print(f"Kakkonto contains {len(compounds)} compounds")

# Get targets for curcumin
curcumin_targets = get_compound_targets(969516)
print(f"Curcumin targets {len(curcumin_targets)} proteins")
```

---

## License

**CC BY-SA 4.0** (Creative Commons Attribution-ShareAlike 4.0 International)

- **Permits:** Commercial use with attribution
- **Requires:** Attribution, ShareAlike for derivatives
- **Source:** Institute of Natural Medicine, University of Toyama

---

## Limitations

1. **No bulk download endpoint**: Must query individual entities
2. **Rate limiting**: Not documented, use reasonable delays
3. **Prediction confidence**: No confidence scores on compound-target relationships
4. **UniProt mapping**: Gene IDs used, manual UniProt mapping required

---

## Integration Recommendations

1. **Primary Use**: Formula-compound-target traversal for network pharmacology
2. **Supplement**: Add BATMAN-TCM 2.0 for TCM overlap and larger TTI predictions
3. **Cross-reference**: Map compound IDs to PubChem for structure data
4. **Validation**: Use HIT 2.0 for experimentally validated interactions

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `formula_id` | String identifier for Kampo medicine formula | `KT` |
| `crude_id` | Integer identifier for crude drug (0-179) | `40` |
| `compound_id` | PubChem CID for natural compound | `969516` |
| `protein_id` | NCBI Gene ID for target protein | `2475` |
| `affinity_kcal_mol` | Predicted binding energy from docking simulation | `-7.8` |
| `name_jp` | Japanese kanji name for formula or crude drug | `葛根湯` |
| `phonetic` | Japanese kana pronunciation guide | `かっこんとう` |
| `origin` | Botanical or biological source of crude drug | `Cinnamomum cassia` |
| `domain_id` | Identifier for protein structural domain used in docking | `14339` |
| `pathway_id` | KEGG pathway identifier for functional annotation | `hsa04150` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Kampo Medicine | Traditional Japanese herbal medicine system | formula_id |
| Crude Drug | Raw medicinal material (herb, mineral, animal product) | crude_id, origin |
| Molecular Docking | Computational prediction of ligand-protein binding | affinity_kcal_mol |
| 4-Layer Hierarchy | KampoDB organization: Formula -> Crude Drug -> Compound -> Protein | Data Model |
| Network Pharmacology | Systems approach analyzing multiple targets of traditional medicines | Compound-Protein relationships |
| Binding Affinity | Strength of molecular interaction (more negative = stronger) | affinity_kcal_mol |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| KampoDB | Kampo Database | Traditional Japanese medicine database |
| CID | Compound Identifier | PubChem compound ID |
| NCBI | National Center for Biotechnology Information | Gene ID source |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |
| GO | Gene Ontology | Functional annotation |
| PDB | Protein Data Bank | 3D structure repository |
| LBP | Ligand-Based Prediction | Chemical similarity method |
| SBP | Structure-Based Prediction | Docking-based method |
| TCM | Traditional Chinese Medicine | Related traditional medicine |
| CC BY-SA | Creative Commons Attribution-ShareAlike | License type |

---

## References

1. Sawada R, et al. (2018). KampoDB, database of predicted targets and functional annotations of natural medicines. Scientific Reports 8:11216.

2. Official Website: https://wakanmoview.inm.u-toyama.ac.jp/kampo/

3. API Documentation: https://wakanmoview.inm.u-toyama.ac.jp/kampo/api

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation with live API testing |
