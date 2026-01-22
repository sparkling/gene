# PubChem Compound Database Schema

**Document ID:** PUBCHEM-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Data Source URL:** https://pubchem.ncbi.nlm.nih.gov/

---

## TL;DR

PubChem is the world's largest free chemistry database containing 118M+ compounds with computed properties, bioassay data, and extensive cross-references. Access via PUG REST API at 5 requests/second. Key endpoints: `/compound/cid/`, `/compound/name/`, `/compound/smiles/`. Output formats: JSON, XML, CSV, SDF, TXT. **License: Public Domain (no restrictions)**. Core data includes molecular properties (weight, formula, SMILES, InChI), 3D conformers, synonyms, and links to ChEBI, ChEMBL, and 800+ other sources.

---

## License Information

| Attribute | Value |
|-----------|-------|
| **License** | Public Domain |
| **Attribution Required** | No |
| **Commercial Use** | Allowed (unrestricted) |
| **Modifications** | Allowed (unrestricted) |
| **Redistribution** | Allowed (unrestricted) |
| **NCBI Policy** | https://www.ncbi.nlm.nih.gov/home/about/policies/ |

**Note:** PubChem data is in the public domain and may be freely used without restriction. However, users are requested to cite PubChem as the source of data.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **Compounds (CIDs)** | 118,000,000+ |
| **Substances (SIDs)** | 300,000,000+ |
| **BioAssays (AIDs)** | 1,500,000+ |
| **Protein Targets** | 15,000+ |
| **Data Sources** | 800+ |
| **Gene Targets** | 30,000+ |
| **Patents** | 40,000,000+ |

---

## API Overview

### Base URL

```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/
```

### URL Structure

```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/{domain}/{namespace}/{identifiers}/{operation}/{output}
```

### Rate Limits

| Constraint | Limit |
|------------|-------|
| **Requests per second** | 5 maximum |
| **Concurrent requests** | 5 maximum |
| **Request timeout** | 30 seconds |
| **Maximum CIDs per request** | 100 (properties), 10 (full records) |

**Policy:** Users exceeding rate limits may be temporarily blocked from accessing PubChem.

---

## Primary Endpoints

### Compound Retrieval Endpoints

| Endpoint Pattern | Description | Example |
|------------------|-------------|---------|
| `/compound/cid/{cid}` | By PubChem Compound ID | `/compound/cid/2244` |
| `/compound/name/{name}` | By compound name | `/compound/name/aspirin` |
| `/compound/smiles/{smiles}` | By SMILES string | `/compound/smiles/CC(=O)C` |
| `/compound/inchi/{inchi}` | By InChI string | `/compound/inchi/InChI=1S/...` |
| `/compound/inchikey/{key}` | By InChI Key | `/compound/inchikey/BSYNR...` |
| `/compound/formula/{formula}` | By molecular formula | `/compound/formula/C9H8O4` |

### Operation Types

| Operation | Description | Example |
|-----------|-------------|---------|
| `/property/{properties}` | Get specific properties | `/property/MolecularWeight,SMILES` |
| `/synonyms` | Get compound synonyms | `/synonyms/JSON` |
| `/cids` | Get matching CIDs | `/cids/JSON` |
| `/sids` | Get substance IDs | `/sids/JSON` |
| `/aids` | Get assay IDs | `/aids/JSON` |
| `/record` | Get full record | `/record/JSON` |
| `/conformers` | Get 3D conformer IDs | `/conformers/TXT` |
| `/xrefs` | Get cross-references | `/xrefs/SourceName/JSON` |

### Output Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| **JSON** | `.json` | Recommended for programmatic access |
| **XML** | `.xml` | Full structured data |
| **CSV** | `.csv` | Tabular property data |
| **TXT** | `.txt` | Single property only |
| **SDF** | `.sdf` | Structure-Data File (3D coordinates) |
| **ASN.1** | `.asn` | NCBI ASN.1 format |
| **PNG** | `.png` | 2D/3D structure images |

---

## Compound Property Dictionary

### Basic Molecular Properties

| Property | Type | Description |
|----------|------|-------------|
| `CID` | Integer | PubChem Compound Identifier (primary key) |
| `MolecularFormula` | String | Molecular formula (e.g., "C9H8O4") |
| `MolecularWeight` | Float | Molecular weight in g/mol |
| `ExactMass` | Float | Exact mass (most abundant isotope) |
| `MonoisotopicMass` | Float | Monoisotopic molecular mass |
| `IUPACName` | String | IUPAC systematic name |

### Structure Representations

| Property | Type | Description |
|----------|------|-------------|
| `CanonicalSMILES` | String | Canonical SMILES (no stereochemistry) |
| `IsomericSMILES` | String | SMILES with stereochemistry |
| `InChI` | String | IUPAC InChI identifier |
| `InChIKey` | String | 27-character InChI hash key |
| `Fingerprint2D` | String | Binary fingerprint (encoded) |

### Physicochemical Properties

| Property | Type | Description |
|----------|------|-------------|
| `XLogP` | Float | Computed octanol/water partition coefficient |
| `TPSA` | Float | Topological Polar Surface Area (A^2) |
| `Complexity` | Float | Bertz/Hendrickson/Ihlenfeldt complexity |
| `Charge` | Integer | Net formal charge |
| `HBondDonorCount` | Integer | Hydrogen bond donor count |
| `HBondAcceptorCount` | Integer | Hydrogen bond acceptor count |
| `RotatableBondCount` | Integer | Number of rotatable bonds |
| `HeavyAtomCount` | Integer | Non-hydrogen atom count |
| `IsotopeAtomCount` | Integer | Atoms with non-standard isotopes |
| `CovalentUnitCount` | Integer | Covalently bonded units |

### Stereochemistry Properties

| Property | Type | Description |
|----------|------|-------------|
| `AtomStereoCount` | Integer | Total stereocenters |
| `DefinedAtomStereoCount` | Integer | Defined atom stereocenters |
| `UndefinedAtomStereoCount` | Integer | Undefined atom stereocenters |
| `BondStereoCount` | Integer | Total bond stereocenters |
| `DefinedBondStereoCount` | Integer | Defined bond stereocenters (E/Z) |
| `UndefinedBondStereoCount` | Integer | Undefined bond stereocenters |

### 3D Conformer Properties

| Property | Type | Description |
|----------|------|-------------|
| `Volume3D` | Float | Molecular volume (A^3) |
| `XStericQuadrupole3D` | Float | X steric quadrupole |
| `YStericQuadrupole3D` | Float | Y steric quadrupole |
| `ZStericQuadrupole3D` | Float | Z steric quadrupole |
| `FeatureCount3D` | Integer | Total pharmacophore features |
| `FeatureAcceptorCount3D` | Integer | HB acceptor features |
| `FeatureDonorCount3D` | Integer | HB donor features |
| `FeatureAnionCount3D` | Integer | Anionic features |
| `FeatureCationCount3D` | Integer | Cationic features |
| `FeatureRingCount3D` | Integer | Ring features |
| `FeatureHydrophobeCount3D` | Integer | Hydrophobic features |
| `ConformerModelRMSD3D` | Float | Conformer model RMSD |
| `EffectiveRotorCount3D` | Float | Effective rotor count |
| `ConformerCount3D` | Integer | Number of conformers |

---

## Sample API Responses

### Property Request

**Request:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularWeight,MolecularFormula,CanonicalSMILES,InChI,IUPACName/JSON
```

**Response:**
```json
{
  "PropertyTable": {
    "Properties": [
      {
        "CID": 2244,
        "MolecularFormula": "C9H8O4",
        "MolecularWeight": 180.16,
        "CanonicalSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
        "IUPACName": "2-acetyloxybenzoic acid"
      }
    ]
  }
}
```

### Multiple Compounds Request

**Request:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244,5090,3672/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON
```

**Response:**
```json
{
  "PropertyTable": {
    "Properties": [
      {
        "CID": 2244,
        "MolecularWeight": 180.16,
        "XLogP": 1.2,
        "TPSA": 63.6,
        "HBondDonorCount": 1,
        "HBondAcceptorCount": 4
      },
      {
        "CID": 5090,
        "MolecularWeight": 194.19,
        "XLogP": 2.3,
        "TPSA": 37.3,
        "HBondDonorCount": 1,
        "HBondAcceptorCount": 3
      },
      {
        "CID": 3672,
        "MolecularWeight": 206.28,
        "XLogP": 3.5,
        "TPSA": 49.3,
        "HBondDonorCount": 1,
        "HBondAcceptorCount": 2
      }
    ]
  }
}
```

### CID Lookup by Name

**Request:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/cids/JSON
```

**Response:**
```json
{
  "IdentifierList": {
    "CID": [2244]
  }
}
```

### Synonyms Request

**Request:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/synonyms/JSON
```

**Response:**
```json
{
  "InformationList": {
    "Information": [
      {
        "CID": 2244,
        "Synonym": [
          "aspirin",
          "Acetylsalicylic acid",
          "50-78-2",
          "2-Acetoxybenzoic acid",
          "Acylpyrin",
          "Ecotrin",
          "Acetysal",
          "Acetosalic acid",
          "Bayer Aspirin",
          "O-Acetylsalicylic acid"
        ]
      }
    ]
  }
}
```

### Full Record Request (Abbreviated)

**Request:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/JSON
```

**Response (abbreviated):**
```json
{
  "PC_Compounds": [
    {
      "id": {
        "id": {
          "cid": 2244
        }
      },
      "atoms": {
        "aid": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        "element": [6, 6, 6, 6, 6, 6, 6, 8, 8, 6, 8, 8, 6]
      },
      "bonds": {
        "aid1": [1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 10, 10],
        "aid2": [2, 6, 3, 7, 4, 5, 6, 7, 8, 10, 11, 12],
        "order": [2, 1, 1, 1, 2, 1, 2, 1, 1, 1, 2, 1]
      },
      "coords": [...],
      "charge": 0,
      "props": [
        {
          "urn": {
            "label": "Molecular Formula"
          },
          "value": {
            "sval": "C9H8O4"
          }
        }
      ],
      "count": {
        "heavy_atom": 13,
        "atom_chiral": 0,
        "bond_chiral": 0,
        "isotope_atom": 0,
        "covalent_unit": 1,
        "tautomers": -1
      }
    }
  ]
}
```

---

## Cross-Reference Integration

### Linked Databases

PubChem provides cross-references to 800+ external databases. Key databases for knowledge integration:

| Database | Relationship | Access Pattern |
|----------|--------------|----------------|
| **ChEBI** | Bidirectional | Via SID or depositor synonyms |
| **ChEMBL** | Bidirectional via BioAssay | Via AID cross-references |
| **DrugBank** | Depositor source | Via SID mapping |
| **KEGG** | Depositor source | Via synonyms/SID |
| **UniProt** | Target links | Via BioAssay targets |
| **PDB** | Structure links | Via ligand mappings |
| **MeSH** | Vocabulary mapping | Via RDF/synonyms |
| **CAS Registry** | Identifier mapping | Via synonyms (not authoritative) |

### Cross-Reference API

**Get external references:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/xrefs/SourceName/JSON
```

**Get substances from a source:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceid/ChEBI/15365/sids/JSON
```

### ChEBI Mapping

ChEBI is both a depositor source and provides semantic annotations:

```
# Get ChEBI-deposited substances for a CID
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/sids/JSON?source=ChEBI

# Search for ChEBI ID in synonyms
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/CHEBI:15365/cids/JSON
```

### ChEMBL Mapping

ChEMBL integration is via BioAssay:

```
# Get BioAssays sourced from ChEMBL
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/source/ChEMBL/aids/JSON

# Get activities for a compound
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/assaysummary/JSON
```

---

## BioAssay Data Structure

### Assay Summary Response

**Request:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/assaysummary/JSON
```

**Response (abbreviated):**
```json
{
  "Table": {
    "Columns": {
      "Column": ["AID", "SourceName", "TargetGI", "TargetName", "Activity", "ActivityValue", "ActivityUnit"]
    },
    "Row": [
      {
        "Cell": [
          "1234",
          "ChEMBL",
          "166897622",
          "Cyclooxygenase-1",
          "Active",
          "0.3",
          "uM"
        ]
      }
    ]
  }
}
```

### Activity Outcome Codes

| Code | Meaning |
|------|---------|
| 1 | Inactive |
| 2 | Active |
| 3 | Inconclusive |
| 4 | Unspecified |
| 5 | Probe |

---

## Similarity and Substructure Search

### Similarity Search

```
# 2D fingerprint similarity (default threshold 90%)
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/CC(=O)OC1=CC=CC=C1C(=O)O/cids/JSON

# With custom threshold
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/CC(=O)OC1=CC=CC=C1C(=O)O/cids/JSON?Threshold=95
```

### Substructure Search

```
# Substructure search by SMILES
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/c1ccccc1C(=O)O/cids/JSON
```

### Async Operations

Large searches return a job key for polling:

```json
{
  "Waiting": {
    "ListKey": "624510668210559704"
  }
}
```

**Poll for results:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/624510668210559704/cids/JSON
```

---

## Common Use Case Examples

### 1. Get All Properties for Drug Discovery

```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount/JSON
```

### 2. Convert Name to Structure

```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/caffeine/property/CanonicalSMILES,InChI,InChIKey/JSON
```

### 3. Get 2D Structure Image

```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/PNG?image_size=large
```

### 4. Get 3D Conformer

```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/SDF?record_type=3d
```

### 5. Batch Property Retrieval (CSV)

```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244,5090,3672,5988,2519/property/MolecularFormula,MolecularWeight,CanonicalSMILES,XLogP,TPSA/CSV
```

---

## Data Quality Notes

### Computed vs. Deposited Data

| Data Type | Source | Reliability |
|-----------|--------|-------------|
| Structure (SMILES, InChI) | Depositors + standardization | High |
| Molecular properties | PubChem computed | High |
| Synonyms | Depositors | Variable (may include errors) |
| BioAssay data | Depositors (ChEMBL, etc.) | Depends on source |
| Cross-references | Depositors | Variable |

### Standardization

PubChem applies standardization to all deposited structures:
- Neutralization of charges where appropriate
- Standardization of tautomeric forms
- Generation of canonical SMILES
- InChI/InChIKey computation

---

## Integration Recommendations

### For Knowledge Graph Integration

1. **Use CID as primary identifier** - Stable, unique compound ID
2. **Use InChIKey for deduplication** - 27-character hash enables cross-database matching
3. **Map to ChEBI/ChEMBL** - Via depositor cross-references or InChIKey
4. **Cache synonyms** - High variability, useful for text mining
5. **Batch requests** - Up to 100 CIDs per property request

### Identifier Mapping Priority

```
1. InChIKey (universal structure hash)
2. InChI (full structure identifier)
3. PubChem CID (database-specific)
4. CanonicalSMILES (structure notation)
5. Name/Synonyms (ambiguous, use cautiously)
```

---

## References

1. Kim S, et al. (2023) "PubChem 2023 update." Nucleic Acids Res. 51(D1):D1373-D1380.

2. PUG REST Documentation: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

3. PubChem Source Information: https://pubchem.ncbi.nlm.nih.gov/sources/

4. IUPAC FAIR Chemistry Cookbook - PUG REST: https://iupac.github.io/WFChemCookbook/datasources/pubchem_pugrest1.html

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
