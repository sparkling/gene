---
id: schema-pubchem
title: "PubChem Compound Database Schema"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, compounds, chemistry, pubchem]
---

**Parent:** [Schema Documentation](./_index.md)

# PubChem Compound Database Schema

## Overview
PubChem is a comprehensive database containing information about chemical molecules and their biological activities. The compound database specifically contains information about small chemical molecules.

## API Access
- **REST API**: `https://pubchem.ncbi.nlm.nih.gov/rest/pug`
- **PUG View API**: `https://pubchem.ncbi.nlm.nih.gov/rest/pug_view`

## Core Data Structures

### Compound Record
```json
{
  "cid": "integer",
  "iupac_name": "string",
  "molecular_formula": "string",
  "molecular_weight": "float",
  "canonical_smiles": "string",
  "isomeric_smiles": "string",
  "inchi": "string",
  "inchi_key": "string",
  "synonyms": ["string"],
  "properties": {
    "XLogP": "float",
    "ExactMass": "float",
    "MonoisotopicMass": "float",
    "TPSA": "float",
    "Complexity": "float",
    "Charge": "integer",
    "HBondDonorCount": "integer",
    "HBondAcceptorCount": "integer",
    "RotatableBondCount": "integer",
    "HeavyAtomCount": "integer",
    "IsotopeAtomCount": "integer",
    "AtomStereoCount": "integer",
    "DefinedAtomStereoCount": "integer",
    "UndefinedAtomStereoCount": "integer",
    "BondStereoCount": "integer",
    "DefinedBondStereoCount": "integer",
    "UndefinedBondStereoCount": "integer",
    "CovalentUnitCount": "integer",
    "Volume3D": "float",
    "XStericQuadrupole3D": "float",
    "YStericQuadrupole3D": "float",
    "ZStericQuadrupole3D": "float",
    "FeatureCount3D": "integer",
    "FeatureAcceptorCount3D": "integer",
    "FeatureDonorCount3D": "integer",
    "FeatureAnionCount3D": "integer",
    "FeatureCationCount3D": "integer",
    "FeatureRingCount3D": "integer",
    "FeatureHydrophobeCount3D": "integer",
    "ConformerModelRMSD3D": "float",
    "EffectiveRotorCount3D": "integer",
    "ConformerCount3D": "integer",
    "Fingerprint2D": "string"
  }
}
```

### Bioactivity Data
```json
{
  "aid": "integer",
  "source_name": "string",
  "assay_type": "string",
  "target": {
    "name": "string",
    "gene_id": "integer",
    "accession": "string"
  },
  "activity": "string",
  "value": "float",
  "unit": "string"
}
```

### Patent Data
```json
{
  "patent_id": "string",
  "country": "string",
  "date": "date",
  "title": "string"
}
```

### Literature References
```json
{
  "pmid": "integer",
  "title": "string",
  "authors": ["string"],
  "journal": "string",
  "year": "integer",
  "doi": "string"
}
```

### Vendors and Suppliers
```json
{
  "vendor_name": "string",
  "catalog_number": "string",
  "availability": "string"
}
```

### Safety and Toxicity
```json
{
  "ghs_hazards": ["string"],
  "precautionary_statements": ["string"],
  "ld50": {
    "species": "string",
    "route": "string",
    "value": "float",
    "unit": "string"
  }
}
```

### Conformer Data (3D Structure)
```json
{
  "conformer_id": "string",
  "atoms": [
    {
      "element": "string",
      "x": "float",
      "y": "float",
      "z": "float"
    }
  ],
  "bonds": [
    {
      "atom1_id": "integer",
      "atom2_id": "integer",
      "order": "integer"
    }
  ],
  "energy": "float"
}
```

## API Endpoints

### Search Endpoints
- **Name Search**: `/compound/name/{name}/JSON`
- **SMILES Search**: `/compound/smiles/{smiles}/JSON`
- **InChI Search**: `/compound/inchi/JSON`
- **Formula Search**: `/compound/formula/{formula}/JSON`
- **Similarity Search**: `/compound/similarity/smiles/{smiles}/JSON`
- **Substructure Search**: `/compound/substructure/smiles/{smiles}/JSON`

### Retrieval Endpoints
- **By CID**: `/compound/cid/{cid}/JSON`
- **Properties**: `/compound/cid/{cid}/property/{properties}/JSON`
- **Synonyms**: `/compound/cid/{cid}/synonyms/JSON`
- **Classification**: `/compound/cid/{cid}/classification/JSON`
- **Assays**: `/compound/cid/{cid}/assaysummary/JSON`
- **Bioactivity**: `/compound/cid/{cid}/concise/JSON`
- **Patents**: `/compound/cid/{cid}/patents/JSON`
- **Literature**: `/compound/cid/{cid}/xrefs/PubMedID/JSON`

### Download Endpoints
- **SDF Download**: `/compound/cid/{cid}/SDF`
- **XML Download**: `/compound/cid/{cid}/XML`
- **PNG Image**: `/compound/cid/{cid}/PNG`

## Common Query Parameters

- **return_type**: JSON, XML, SDF, CSV, PNG
- **search_type**: exact, similarity, substructure
- **threshold**: 0-100 (for similarity searches)
- **max_records**: integer (limit results)

## Example Queries

### Get Compound by Name
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/JSON
```

### Get Molecular Properties
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,MolecularWeight,CanonicalSMILES/JSON
```

### Similarity Search
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/CC(=O)OC1=CC=CC=C1C(=O)O/JSON?Threshold=95
```

### Substructure Search
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/c1ccccc1/JSON
```

### Batch Retrieval
```
POST https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/JSON
Body: {"cid": [2244, 5352, 5362]}
```

## Data Relationships

```
Compound (CID)
├── Properties
│   ├── Physical Properties
│   ├── Chemical Properties
│   └── 3D Properties
├── Identifiers
│   ├── SMILES
│   ├── InChI/InChIKey
│   ├── Synonyms
│   └── Cross-references
├── Bioactivity
│   ├── Assays (AIDs)
│   ├── Targets
│   └── Activity Values
├── Literature
│   └── PubMed References
├── Patents
├── Vendors
├── Classification
│   ├── Chemical Hierarchy
│   └── Pharmacological Classes
└── Safety/Toxicity
    ├── GHS Hazards
    └── LD50 Values
```

## Rate Limits
- **Max 5 requests per second**
- **Max 400 requests per minute**
- For higher throughput, consider using FTP bulk downloads

## Bulk Downloads
Available via FTP: `ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/`

### Available Formats
- **SDF**: Structure-Data File
- **XML**: Full compound records
- **ASN.1**: Binary format
- **JSON**: Individual compound files
- **CSV**: Selected properties

## Common Use Cases

1. **Drug Discovery**: Search for similar compounds, retrieve bioactivity data
2. **Chemical Safety**: Get toxicity and hazard information
3. **Literature Mining**: Link compounds to research papers
4. **Structure Analysis**: Compare 3D conformations
5. **Commercial Sourcing**: Find vendors and suppliers

## Integration Points

- **DrugBank**: Cross-referenced drug identifiers
- **ChEMBL**: Bioactivity data
- **UniProt**: Target protein information
- **PubMed**: Literature references
- **Patent databases**: Patent information
- **Chemical vendors**: Commercial availability

## Notes

- CIDs are permanent identifiers
- Some compounds may have multiple conformers
- Stereoisomers get separate CIDs
- Regular updates with new compounds and data
- Historical versions available through versioning system

## Data Format

| Format | Description |
|--------|-------------|
| Primary | SDF (Structure-Data File) |
| Alternative | XML, JSON, CSV, ASN.1 |
| Compression | gzip (.gz) |
| Encoding | UTF-8 |
| API Response | JSON, XML, SDF, PNG |

---

## Download

### Bulk Data Downloads

| Format | Size (approx) | URL |
|--------|---------------|-----|
| SDF (Structure-Data File) | ~500 GB (uncompressed) | ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/ |
| XML | ~300 GB (uncompressed) | ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/ |
| JSON | ~400 GB (uncompressed) | ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/JSON/ |
| ASN.1 Binary | ~200 GB (compressed) | ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/ASN/ |
| CSV (selected properties) | ~50 GB (compressed) | ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/CSV/ |

### REST API

**Base URL:** https://pubchem.ncbi.nlm.nih.gov/rest/pug

**Rate Limits:**
- 5 requests/second maximum
- 400 requests/minute maximum
- Higher throughput: Use FTP bulk downloads

---

## Data Set Size

| Metric | Value |
|--------|-------|
| **Total Compounds** | 115,000,000+ |
| **Active Compounds** | 95,000,000+ with bioactivity data |
| **Unique Structures** | 95,000,000+ (CID) |
| **Substances Tracked** | 300,000,000+ (SID) |
| **Bioassays** | 2,000,000+ (AID) |
| **Bioactivity Records** | 500,000,000+ |
| **Patent References** | 50,000,000+ |
| **Literature References** | 100,000,000+ |
| **3D Conformers** | 50,000,000+ compounds |
| **API Response Rate** | ~1,000 queries/hour typical |
| **FTP Total Size** | ~1.5 TB (all formats, uncompressed) |
| **Update Frequency** | Daily incremental updates |
| **Last Full Rebuild** | Monthly |

---

## Sample Data

### Example Record
```json
{
  "cid": 2244,
  "compound_name": "Aspirin",
  "molecular_formula": "C9H8O4",
  "molecular_weight": 180.16,
  "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
}
```

### Sample Query Result
| cid | compound_name | molecular_formula | molecular_weight | canonical_smiles |
|-----|-------------|------------------|------------------|------------------|
| 2244 | Aspirin | C9H8O4 | 180.16 | CC(=O)Oc1ccccc1C(=O)O |
| 3672 | Ibuprofen | C13H18O2 | 206.28 | CC(C)Cc1ccc(cc1)C(C)C(=O)O |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "CID123456" |
| `name` | string | Entity name | "Aspirin" |
| `type` | string | Record type | "compound" / "substance" / "bioassay" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `has_property` | Property | N:M |
| `similar_to` | Compound | N:M |
| `tested_in` | Bioassay | N:M |

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| PubChem | CC0 (Public Domain) | Yes |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `CID` | Compound Identifier, unique numeric ID for each compound in PubChem | 2244 (Aspirin) |
| `AID` | Assay Identifier, unique ID for each bioassay experiment | 504327 |
| `canonical_smiles` | Standardized SMILES notation without stereochemistry information | CC(=O)OC1=CC=CC=C1C(=O)O |
| `isomeric_smiles` | SMILES notation preserving stereochemistry information | C[C@@H](O)c1ccccc1 |
| `inchi_key` | 27-character hash of InChI for fast compound lookup | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| `XLogP` | Computed partition coefficient (octanol/water), measures lipophilicity | 1.31 |
| `TPSA` | Topological Polar Surface Area, predictor of drug absorption | 63.6 A^2 |
| `HBondDonorCount` | Number of hydrogen bond donor atoms in the molecule | 1 |
| `HBondAcceptorCount` | Number of hydrogen bond acceptor atoms in the molecule | 4 |
| `Complexity` | Computed molecular complexity score | 212 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Conformer | 3D arrangement of atoms representing a low-energy molecular shape | Conformer3D fields |
| Fingerprint | Binary vector encoding structural features for similarity search | Fingerprint2D field |
| GHS Hazard | Globally Harmonized System classification for chemical hazards | Safety data |
| LD50 | Lethal Dose 50%, dose killing 50% of test animals | Toxicity data |
| Stereocenter | Atom with four different substituents creating chirality | AtomStereoCount field |
| Rotatable Bond | Single bond that can freely rotate (excludes rings, double bonds) | RotatableBondCount |
| Heavy Atom | Non-hydrogen atom in the molecular structure | HeavyAtomCount field |
| Bioassay | Experimental test measuring compound activity against a target | AID records |
| Pharmacophore | 3D arrangement of features required for biological activity | 3D features |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CID | Compound Identifier | PubChem unique compound ID |
| AID | Assay Identifier | PubChem unique assay ID |
| SID | Substance Identifier | PubChem substance record ID |
| SMILES | Simplified Molecular Input Line Entry System | Chemical notation |
| InChI | International Chemical Identifier | IUPAC standard identifier |
| SDF | Structure-Data File | Chemical file format |
| TPSA | Topological Polar Surface Area | Absorption predictor |
| GHS | Globally Harmonized System | Hazard classification |
| LD50 | Lethal Dose 50% | Toxicity measure |
| PUG | Power User Gateway | PubChem API system |
| MW | Molecular Weight | Daltons |
| XLogP | Computed LogP | Lipophilicity estimate |
| RMSD | Root Mean Square Deviation | Conformer comparison |
| HBA | Hydrogen Bond Acceptor | Lipinski parameter |
| HBD | Hydrogen Bond Donor | Lipinski parameter |
| FTP | File Transfer Protocol | Bulk download method |

---

## References
- [PubChem REST API Documentation](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
- [PUG View Documentation](https://pubchem.ncbi.nlm.nih.gov/docs/pug-view)
- [FTP Download Instructions](https://ftp.ncbi.nlm.nih.gov/pubchem/)