---
id: schema-coconut
title: "COCONUT Natural Products Database Schema"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, natural-products, compounds, chemistry]
---

**Parent:** [Schema Documentation](./_index.md)

# COCONUT - COlleCtion of Open Natural prodUcTs Database Schema

## Overview
COCONUT is a comprehensive database of natural products (NPs) aggregated from various open sources. It provides structural, physicochemical, and biological information about natural compounds.

## Database Information
- **Type**: PostgreSQL
- **REST API**: Available
- **Web Interface**: https://coconut.naturalproducts.net
- **API Base URL**: https://coconut.naturalproducts.net/api

## Core Tables

### compounds
Main table containing natural product information.

```sql
CREATE TABLE compounds (
  id SERIAL PRIMARY KEY,
  coconut_id VARCHAR(20) UNIQUE NOT NULL,
  name TEXT,
  iupac_name TEXT,
  molecular_formula VARCHAR(255),
  molecular_weight DECIMAL(10,4),
  canonical_smiles TEXT,
  isomeric_smiles TEXT,
  inchi TEXT,
  inchi_key VARCHAR(27),
  sugar_free_smiles TEXT,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

#### Key Fields
- **coconut_id**: Unique COCONUT identifier (e.g., CNP0123456)
- **name**: Common or trivial name of the natural product
- **iupac_name**: IUPAC systematic name
- **molecular_formula**: Chemical formula (e.g., C21H30O2)
- **molecular_weight**: Molecular weight in Daltons
- **canonical_smiles**: Simplified molecular-input line-entry system
- **isomeric_smiles**: SMILES with stereochemistry
- **inchi**: International Chemical Identifier
- **inchi_key**: Hashed InChI for fast lookups
- **sugar_free_smiles**: SMILES without glycosidic moieties

### properties
Physicochemical properties of compounds.

```sql
CREATE TABLE properties (
  id SERIAL PRIMARY KEY,
  compound_id INTEGER REFERENCES compounds(id),
  total_atom_number INTEGER,
  heavy_atom_number INTEGER,
  aromatic_heavy_atom_number INTEGER,
  bond_number INTEGER,
  rotatable_bond_number INTEGER,
  hbond_acceptor INTEGER,
  hbond_donor INTEGER,
  murko_framework TEXT,
  number_of_carbons INTEGER,
  number_of_nitrogens INTEGER,
  number_of_oxygens INTEGER,
  number_of_sulphurs INTEGER,
  number_of_rings INTEGER,
  aromatic_rings_count INTEGER,
  qed_drug_likeliness DECIMAL(5,4),
  lipinski_rule_of_5 INTEGER,
  topological_polar_surface_area DECIMAL(10,4),
  alogp DECIMAL(10,4),
  xlogp DECIMAL(10,4),
  found_in_databases TEXT[],
  UNIQUE(compound_id)
);
```

#### Key Properties
- **Heavy atoms**: Non-hydrogen atoms
- **H-bond donors/acceptors**: Hydrogen bonding capability
- **Rotatable bonds**: Conformational flexibility
- **Murko framework**: Core scaffold structure
- **QED**: Quantitative Estimate of Drug-likeness (0-1)
- **Lipinski Rule of 5**: Drug-likeness criteria (0-4 violations)
- **TPSA**: Topological Polar Surface Area
- **aLogP/xLogP**: Partition coefficient (lipophilicity)

### sources
Data sources where compounds were found.

```sql
CREATE TABLE sources (
  id SERIAL PRIMARY KEY,
  name VARCHAR(255) UNIQUE NOT NULL,
  url TEXT,
  description TEXT,
  compound_count INTEGER DEFAULT 0
);
```

Common sources include:
- **ZINC**: ZINC Database
- **ChEMBL**: Bioactivity Database
- **PubChem**: NIH Chemical Database
- **LOTUS**: Natural Products Occurrence Database
- **NPASS**: Natural Products Activity and Species Source Database
- **NAPRALERT**: Natural Products Database
- **TCMDB**: Traditional Chinese Medicine Database
- **NORMAN**: NORMAN Suspect List Exchange
- **COCONUT**: Original COCONUT sources

### compound_sources
Many-to-many relationship between compounds and sources.

```sql
CREATE TABLE compound_sources (
  compound_id INTEGER REFERENCES compounds(id),
  source_id INTEGER REFERENCES sources(id),
  external_id VARCHAR(255),
  PRIMARY KEY (compound_id, source_id)
);
```

### organisms
Organisms from which natural products are derived.

```sql
CREATE TABLE organisms (
  id SERIAL PRIMARY KEY,
  name TEXT NOT NULL,
  taxonomic_rank VARCHAR(50),
  ncbi_taxonomy_id INTEGER,
  canonical_name TEXT,
  kingdom VARCHAR(50),
  phylum VARCHAR(100),
  class VARCHAR(100),
  order_name VARCHAR(100),
  family VARCHAR(100),
  genus VARCHAR(100),
  species VARCHAR(100)
);
```

### compound_organisms
Many-to-many relationship between compounds and organisms.

```sql
CREATE TABLE compound_organisms (
  compound_id INTEGER REFERENCES compounds(id),
  organism_id INTEGER REFERENCES organisms(id),
  PRIMARY KEY (compound_id, organism_id)
);
```

### biological_activities
Biological activity annotations.

```sql
CREATE TABLE biological_activities (
  id SERIAL PRIMARY KEY,
  compound_id INTEGER REFERENCES compounds(id),
  activity_type VARCHAR(255),
  target TEXT,
  organism TEXT,
  value DECIMAL(10,4),
  unit VARCHAR(50),
  reference TEXT
);
```

#### Activity Types
- **IC50**: Half maximal inhibitory concentration
- **EC50**: Half maximal effective concentration
- **Ki**: Inhibition constant
- **MIC**: Minimum inhibitory concentration
- **LD50**: Median lethal dose
- **ED50**: Median effective dose

### citations
Literature references for compounds.

```sql
CREATE TABLE citations (
  id SERIAL PRIMARY KEY,
  compound_id INTEGER REFERENCES compounds(id),
  doi VARCHAR(255),
  pmid INTEGER,
  title TEXT,
  authors TEXT,
  journal TEXT,
  year INTEGER,
  citation_text TEXT
);
```

### fragments
Molecular fragments for substructure searching.

```sql
CREATE TABLE fragments (
  id SERIAL PRIMARY KEY,
  compound_id INTEGER REFERENCES compounds(id),
  fragment_smiles TEXT,
  fragment_type VARCHAR(50)
);
```

### stereoisomers
Stereochemical variants of compounds.

```sql
CREATE TABLE stereoisomers (
  id SERIAL PRIMARY KEY,
  parent_compound_id INTEGER REFERENCES compounds(id),
  stereoisomer_compound_id INTEGER REFERENCES compounds(id),
  relationship VARCHAR(50)
);
```

## REST API Endpoints

### Search Endpoints

#### Text Search
```
GET /api/search?query={text}
```

#### Structure Search
```
POST /api/structure/search
Body: {
  "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "searchType": "exact|substructure|similarity",
  "threshold": 0.85
}
```

#### Property Search
```
GET /api/search/properties?mw_min=100&mw_max=500&alogp_max=5
```

### Retrieval Endpoints

#### Get Compound by ID
```
GET /api/compound/{coconut_id}
```

Response:
```json
{
  "coconut_id": "CNP0123456",
  "name": "Example Natural Product",
  "molecular_formula": "C21H30O2",
  "molecular_weight": 314.46,
  "smiles": "CC(C)CCCC(C)CCCC(C)CCCC(C)C",
  "inchi_key": "XXXXXXXXXXXXXXXXXXXXX-XXXXX",
  "properties": {
    "qed": 0.85,
    "lipinski_violations": 0,
    "tpsa": 40.46,
    "alogp": 6.2
  },
  "sources": ["ZINC", "PubChem"],
  "organisms": [
    {
      "name": "Artemisia annua",
      "taxonomy_id": 35608
    }
  ]
}
```

#### Get Compound Properties
```
GET /api/compound/{coconut_id}/properties
```

#### Get Compound Sources
```
GET /api/compound/{coconut_id}/sources
```

#### Get Compound Organisms
```
GET /api/compound/{coconut_id}/organisms
```

#### Get Biological Activities
```
GET /api/compound/{coconut_id}/activities
```

### Batch Operations

#### Batch Retrieval
```
POST /api/batch/compounds
Body: {
  "ids": ["CNP0123456", "CNP0234567", "CNP0345678"]
}
```

#### Batch Download
```
POST /api/download
Body: {
  "ids": [...],
  "format": "sdf|csv|json"
}
```

### Statistics Endpoints

#### Database Statistics
```
GET /api/stats
```

Response:
```json
{
  "total_compounds": 450000,
  "total_organisms": 35000,
  "total_sources": 52,
  "molecular_weight_distribution": {...},
  "most_common_organisms": [...],
  "kingdom_distribution": {...}
}
```

## Query Parameters for Property Search

### Molecular Properties
- `mw_min`, `mw_max`: Molecular weight range
- `alogp_min`, `alogp_max`: Lipophilicity range
- `hba_min`, `hba_max`: H-bond acceptor count
- `hbd_min`, `hbd_max`: H-bond donor count
- `rb_min`, `rb_max`: Rotatable bond count
- `tpsa_min`, `tpsa_max`: Topological polar surface area

### Structure Filters
- `rings_min`, `rings_max`: Number of rings
- `aromatic_rings_min`, `aromatic_rings_max`: Aromatic rings
- `carbons_min`, `carbons_max`: Carbon count
- `nitrogens_min`, `nitrogens_max`: Nitrogen count
- `oxygens_min`, `oxygens_max`: Oxygen count

### Drug-likeness
- `qed_min`, `qed_max`: QED score (0-1)
- `lipinski_violations`: 0-4

### Taxonomy
- `kingdom`: plant|bacteria|fungi|animal
- `organism`: Organism name or taxonomy ID

### Source
- `source`: Database name

## Example Queries

### Search for Artemisinin-like Compounds
```
GET /api/search?query=artemisinin
```

### Find Plant-derived Compounds
```
GET /api/search/properties?kingdom=plant&mw_max=500
```

### Substructure Search for Flavonoids
```
POST /api/structure/search
Body: {
  "smiles": "c1ccc(cc1)C2CC(=O)c3ccccc3O2",
  "searchType": "substructure"
}
```

### Drug-like Natural Products
```
GET /api/search/properties?lipinski_violations=0&qed_min=0.5
```

### Compounds from Specific Organism
```
GET /api/search?organism=Artemisia%20annua
```

## Data Relationships

```
compounds (1) ← (many) compound_sources (many) → (1) sources
compounds (1) ← (many) compound_organisms (many) → (1) organisms
compounds (1) ← (many) biological_activities
compounds (1) ← (many) citations
compounds (1) ← (1) properties
compounds (1) ← (many) fragments
compounds (1) ← (many) stereoisomers → (1) compounds
```

## Download Formats

### SDF (Structure-Data File)
- Complete 3D structure information
- All properties as data fields

### CSV
- Tabular format
- Selected properties

### JSON
- Complete compound data
- Nested relationships

### SMILES
- Text format
- One compound per line: SMILES + Name

## Integration with Other Databases

COCONUT cross-references:
- **ChEMBL**: Bioactivity data
- **PubChem**: Chemical information
- **ZINC**: Commercial availability
- **ChEBI**: Chemical ontology
- **HMDB**: Human metabolome
- **DrugBank**: Drug information
- **FooDB**: Food composition

## Citation
COCONUT online: Collection of Open Natural Products database.
Sorokina M, Merseburger P, Rajan K, Yirik MA, Steinbeck C.
Journal of Cheminformatics 13, 2 (2021).
DOI: 10.1186/s13321-020-00478-9

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `coconut_id` | Unique COCONUT natural product identifier | CNP0123456 |
| `canonical_smiles` | Standardized SMILES string without stereochemistry | CC(C)CCCC(C)C |
| `isomeric_smiles` | SMILES with stereochemistry information preserved | C[C@@H](O)c1ccccc1 |
| `inchi_key` | 27-character hash of InChI for fast structure lookup | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| `sugar_free_smiles` | SMILES representation with glycosidic moieties removed | Aglycone structure |
| `qed_drug_likeliness` | Quantitative Estimate of Drug-likeness score (0-1) | 0.85 |
| `lipinski_rule_of_5` | Number of Lipinski rule violations (0-4) | 0 (drug-like) |
| `murko_framework` | Core molecular scaffold with all side chains removed | Benzene ring |
| `alogp` | Calculated partition coefficient (lipophilicity) | 2.5 |
| `tpsa` | Topological Polar Surface Area in Angstrom squared | 40.46 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Natural Product | Secondary metabolite produced by living organisms | COCONUT scope |
| Molecular Framework | Core scaffold structure after removing substituents | murko_framework |
| Drug-likeness | Properties suggesting potential as oral drug | QED, Lipinski |
| Lipinski Rule of 5 | Criteria for oral bioavailability (MW, LogP, HBD, HBA) | lipinski_rule_of_5 |
| Glycoside | Compound with sugar moiety attached | sugar_free_smiles |
| Stereoisomer | Compounds with same connectivity but different 3D arrangement | isomeric_smiles |
| Substructure Search | Finding compounds containing a specific molecular fragment | API feature |
| Similarity Search | Finding compounds structurally similar to a query | Tanimoto threshold |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| COCONUT | COlleCtion of Open Natural prodUcTs | Database name |
| NP | Natural Product | Compound type |
| SMILES | Simplified Molecular Input Line Entry System | Structure notation |
| InChI | International Chemical Identifier | IUPAC identifier |
| QED | Quantitative Estimate of Drug-likeness | Drug-likeness score |
| TPSA | Topological Polar Surface Area | Absorption predictor |
| HBD | Hydrogen Bond Donor | Lipinski parameter |
| HBA | Hydrogen Bond Acceptor | Lipinski parameter |
| MW | Molecular Weight | Daltons |
| IC50 | Half Maximal Inhibitory Concentration | Activity measure |
| MIC | Minimum Inhibitory Concentration | Antimicrobial activity |
| SDF | Structure-Data File | Chemical file format |
| NCBI | National Center for Biotechnology Information | Taxonomy source |

---

## References
- [COCONUT Website](https://coconut.naturalproducts.net)
- [COCONUT GitHub](https://github.com/Steinbeck-Lab/coconut)
- [API Documentation](https://coconut.naturalproducts.net/api/documentation)
- [Publication](https://doi.org/10.1186/s13321-020-00478-9)