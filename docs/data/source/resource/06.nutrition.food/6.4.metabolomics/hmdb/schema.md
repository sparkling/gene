---
id: schema-hmdb
title: "HMDB - Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [schema, metabolome, metabolites, human, biomarkers]
---

# HMDB - Human Metabolome Database Schema Documentation

**Document ID:** SCHEMA-HMDB
**Status:** Final
**Last Updated:** January 2026
**Data Source URL:** https://hmdb.ca/

---

## TL;DR

The Human Metabolome Database (HMDB) is the most comprehensive database of human metabolites. It contains detailed information about small molecule metabolites found in the human body, including chemical data, clinical data, molecular biology data, and reference spectra for metabolite identification.

---

## Database Statistics

| Metric | Count |
|--------|-------|
| Metabolites | 220,000+ |
| With Concentration Data | 5,000+ |
| MS/MS Spectra | 650,000+ |
| NMR Spectra | 5,000+ |
| Disease Associations | 2,400+ |
| Pathways | 800+ |

---

## Database Schema

### Core Tables

#### 1. Metabolites Table
```sql
CREATE TABLE metabolites (
    id INTEGER PRIMARY KEY,
    accession VARCHAR(20) UNIQUE,     -- HMDB ID (e.g., HMDB0000001)
    version INTEGER,
    name VARCHAR(500),                 -- Primary name
    description TEXT,                  -- Detailed description
    chemical_formula VARCHAR(100),
    average_molecular_weight DECIMAL(15,6),
    monoisotopic_molecular_weight DECIMAL(15,6),
    iupac_name VARCHAR(1000),
    traditional_iupac VARCHAR(1000),
    cas_registry_number VARCHAR(20),
    smiles TEXT,
    inchi TEXT,
    inchikey VARCHAR(50),
    state VARCHAR(20),                 -- Solid, Liquid, Gas
    synthesis_reference TEXT,
    INDEX idx_accession (accession),
    INDEX idx_name (name),
    INDEX idx_inchikey (inchikey)
);
```

#### 2. Concentrations Table
```sql
CREATE TABLE concentrations (
    id INTEGER PRIMARY KEY,
    metabolite_id INTEGER REFERENCES metabolites(id),
    biofluid VARCHAR(50),             -- Blood, Urine, CSF, etc.
    biofluid_type VARCHAR(50),        -- Specific type
    concentration_value DECIMAL(15,6),
    concentration_units VARCHAR(30),   -- uM, ng/mL, etc.
    age VARCHAR(50),                   -- Age group
    sex VARCHAR(20),                   -- Male, Female, Both
    condition VARCHAR(100),            -- Normal, Disease name
    subject_count INTEGER,
    reference_id INTEGER,
    INDEX idx_metabolite (metabolite_id),
    INDEX idx_biofluid (biofluid)
);
```

#### 3. Diseases Table
```sql
CREATE TABLE diseases (
    id INTEGER PRIMARY KEY,
    metabolite_id INTEGER REFERENCES metabolites(id),
    disease_name VARCHAR(500),
    omim_id VARCHAR(20),
    pubmed_id INTEGER,
    concentration_change VARCHAR(50),  -- Elevated, Reduced
    reference_text TEXT,
    INDEX idx_metabolite (metabolite_id),
    INDEX idx_omim (omim_id)
);
```

#### 4. Pathways Table
```sql
CREATE TABLE pathways (
    id INTEGER PRIMARY KEY,
    metabolite_id INTEGER REFERENCES metabolites(id),
    pathway_name VARCHAR(500),
    kegg_id VARCHAR(20),
    smpdb_id VARCHAR(20),             -- SMPDB pathway ID
    INDEX idx_metabolite (metabolite_id)
);
```

#### 5. Spectra Tables
```sql
CREATE TABLE ms_spectra (
    id INTEGER PRIMARY KEY,
    metabolite_id INTEGER REFERENCES metabolites(id),
    spectrum_type VARCHAR(50),         -- MS, MS/MS
    ionization_mode VARCHAR(20),       -- Positive, Negative
    collision_energy VARCHAR(20),
    instrument VARCHAR(100),
    peaks TEXT,                        -- JSON array of peaks
    predicted BOOLEAN,
    INDEX idx_metabolite (metabolite_id)
);

CREATE TABLE nmr_spectra (
    id INTEGER PRIMARY KEY,
    metabolite_id INTEGER REFERENCES metabolites(id),
    nucleus VARCHAR(10),               -- 1H, 13C
    frequency VARCHAR(20),             -- 400 MHz, 500 MHz
    solvent VARCHAR(50),
    peaks TEXT,                        -- JSON array of chemical shifts
    INDEX idx_metabolite (metabolite_id)
);
```

#### 6. External Links Table
```sql
CREATE TABLE external_links (
    id INTEGER PRIMARY KEY,
    metabolite_id INTEGER REFERENCES metabolites(id),
    resource_name VARCHAR(100),
    external_id VARCHAR(100),
    INDEX idx_metabolite (metabolite_id),
    INDEX idx_resource (resource_name)
);
```

---

## Biofluid Types

| Biofluid | Description | Typical Metabolites |
|----------|-------------|---------------------|
| Blood | Whole blood | Amino acids, lipids |
| Serum | Blood without cells/clotting factors | Proteins, metabolites |
| Plasma | Blood without cells | Clinical chemistry |
| Urine | Urinary metabolites | Organic acids, drugs |
| CSF | Cerebrospinal fluid | Neurotransmitters |
| Saliva | Salivary metabolites | Hormones |
| Feces | Fecal metabolites | Bile acids, microbial |
| Breast milk | Maternal milk | Lipids, sugars |
| Sweat | Sweat metabolites | Electrolytes |

---

## JSON Schemas

### Metabolite Object
```json
{
  "accession": "HMDB0000001",
  "version": 5,
  "name": "1-Methylhistidine",
  "description": "1-Methylhistidine (1-MHis) is...",
  "chemical_formula": "C7H11N3O2",
  "average_molecular_weight": 169.181,
  "monoisotopic_molecular_weight": 169.085126611,
  "iupac_name": "(2S)-2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid",
  "cas_registry_number": "332-80-9",
  "smiles": "CN1C=NC(CC(N)C(O)=O)=C1",
  "inchi": "InChI=1S/C7H11N3O2/c1-10-4-9-5(3-10)2-6(8)7(11)12/h3-4,6H,2,8H2,1H3,(H,11,12)/t6-/m0/s1",
  "inchikey": "BRMWTNUJHUMWMS-LURJTMIESA-N",
  "state": "Solid",
  "biological_properties": {
    "cellular_locations": ["Cytoplasm"],
    "biofluid_locations": ["Blood", "Urine"],
    "tissue_locations": ["Muscle"]
  },
  "taxonomy": {
    "kingdom": "Organic compounds",
    "superclass": "Organic acids and derivatives",
    "class": "Carboxylic acids and derivatives"
  }
}
```

### Concentration Object
```json
{
  "biofluid": "Blood",
  "concentration_value": 5.2,
  "concentration_units": "uM",
  "age": "Adult (>18 years)",
  "sex": "Both",
  "condition": "Normal",
  "subject_count": 120,
  "reference": {
    "pubmed_id": 12345678,
    "citation": "Smith et al. (2020)"
  }
}
```

### Disease Association Object
```json
{
  "disease_name": "Histidinemia",
  "omim_id": "235800",
  "concentration_change": "Elevated",
  "reference": {
    "pubmed_id": 23456789,
    "description": "Elevated levels in affected individuals"
  }
}
```

### MS/MS Spectrum Object
```json
{
  "spectrum_type": "MS/MS",
  "ionization_mode": "Positive",
  "collision_energy": "20 eV",
  "instrument": "Q-TOF",
  "precursor_mz": 170.0924,
  "peaks": [
    {"mz": 56.0500, "intensity": 100},
    {"mz": 83.0609, "intensity": 85},
    {"mz": 95.0609, "intensity": 45},
    {"mz": 124.0869, "intensity": 70},
    {"mz": 170.0924, "intensity": 25}
  ],
  "predicted": false
}
```

---

## Data Access

### Web Interface

| Feature | URL |
|---------|-----|
| Browse metabolites | https://hmdb.ca/metabolites |
| Search | https://hmdb.ca/unearth |
| Spectra search | https://hmdb.ca/spectra |

### REST API

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/metabolites/{id}` | GET | Get metabolite by HMDB ID |
| `/metabolites/{id}.json` | GET | JSON format |
| `/metabolites/{id}.xml` | GET | XML format |

### Bulk Downloads

| Format | Size | Content |
|--------|------|---------|
| hmdb_metabolites.zip | ~2 GB | All metabolites (XML) |
| structures.zip | ~200 MB | SDF structures |
| hmdb_proteins.zip | ~500 MB | Associated proteins |
| spectra files | ~5 GB | MS/MS and NMR spectra |

**Download URL:** https://hmdb.ca/downloads

---

## Cross-References

### External Databases

| Database | ID Format | Example |
|----------|-----------|---------|
| PubChem | Numeric CID | 92105 |
| ChEBI | CHEBI:NNNNN | CHEBI:50599 |
| KEGG | CNNNNN | C01152 |
| DrugBank | DBNNNNN | DB00118 |
| Metlin | Numeric | 3741 |
| BioCyc | Alphanumeric | 1-METHYLHISTIDINE |
| FooDB | FDBNNNNN | FDB000001 |

---

## Relationships

### Entity Relationship Diagram
```
metabolites (1) ----< concentrations (N)
    |
    +----< diseases (N)
    |
    +----< pathways (N)
    |
    +----< ms_spectra (N)
    |
    +----< nmr_spectra (N)
    |
    +----< external_links (N)
```

---

## Use Cases

### 1. Get Metabolite Information
```bash
curl "https://hmdb.ca/metabolites/HMDB0000001.json" | jq
```

### 2. Find Disease-Associated Metabolites
```
Search: disease_name = "Diabetes"
Result: Glucose, HbA1c, etc.
```

### 3. MS/MS Spectrum Matching
```
Input: Unknown spectrum peaks
Match: Against predicted/experimental spectra
Result: Candidate metabolite IDs
```

### 4. Biofluid-Specific Query
```sql
SELECT m.name, c.concentration_value, c.concentration_units
FROM metabolites m
JOIN concentrations c ON m.id = c.metabolite_id
WHERE c.biofluid = 'Urine' AND c.condition = 'Normal';
```

---

## License

- **License Type:** Creative Commons Attribution-NonCommercial 4.0
- **Academic Use:** Free
- **Commercial Use:** Requires license
- **Attribution:** Required

---

## Citation

```
Wishart DS, Guo AC, Oler E, et al.
HMDB 5.0: the Human Metabolome Database for 2022.
Nucleic Acids Research. 2022;50(D1):D1285-D1292.
doi: 10.1093/nar/gkab1062

HMDB Version 5.0. https://hmdb.ca/
```

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| HMDB ID | Unique metabolite accession number | HMDB0000001 |
| InChIKey | Hashed chemical structure identifier | BRMWTNUJHUMWMS-... |
| Biofluid | Body fluid containing metabolite | Blood, Urine |
| MS/MS | Tandem mass spectrometry | Fragmentation pattern |
| NMR | Nuclear magnetic resonance spectroscopy | Chemical shifts |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Metabolite | Small molecule product of metabolism | Primary data entity |
| Metabolome | Complete set of metabolites in organism | Database scope |
| Reference spectrum | Known spectrum for identification | MS/MS, NMR |
| Biomarker | Metabolite indicating biological state | Clinical application |
| Endogenous | Produced within the body | Classification |
| Exogenous | From external sources (diet, drugs) | Classification |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HMDB | Human Metabolome Database | Database name |
| MS | Mass Spectrometry | Analytical method |
| MS/MS | Tandem Mass Spectrometry | Fragmentation |
| NMR | Nuclear Magnetic Resonance | Spectroscopy |
| CSF | Cerebrospinal Fluid | Biofluid |
| OMIM | Online Mendelian Inheritance in Man | Disease database |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |
| SMPDB | Small Molecule Pathway Database | Pathway database |

---

## Related Documents

- [Download Instructions](./download.md)
- [Exposome-Explorer](../exposome.explorer/_index.md) - Exposure biomarkers
- [FooDB](../../6.1.food.composition/foodb/_index.md) - Food compound database
