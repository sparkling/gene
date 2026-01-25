---
id: schema-chembl
title: "ChEMBL Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, bioactivity, drugs, compounds]
---

**Parent:** [Schema Documentation](./README.md)

# ChEMBL Database Schema

**Document ID:** CHEMBL-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** ChEMBL 36 (July 2025)

---

## TL;DR

ChEMBL 36 contains 24.3M bioactivity records for 2.9M compounds across 17,803 targets from 99,142 publications. The core entity chain is: Molecule -> Activity -> Assay -> Target. Key tables include MOLECULE_DICTIONARY, ACTIVITIES, ASSAYS, and TARGET_DICTIONARY with extensive cross-references to UniProt, PubChem, and DrugBank.

---

## Database Statistics (ChEMBL 36)

| Entity | Count |
|--------|-------|
| **Activities** | 24,267,312 |
| **Compound Records** | 3,774,137 |
| **Distinct Compounds** | 2,878,135 |
| **Targets** | 17,803 |
| **Publications** | 99,142 |
| **Assays** | ~1.6M (estimated) |
| **Natural Products (flagged)** | ~64,000 |

---

## Available Downloads

| File | Size | Description |
|------|------|-------------|
| chembl_36_sqlite.tar.gz | 5.2 GB | SQLite database |
| chembl_36_mysql.tar.gz | 1.9 GB | MySQL format |
| chembl_36_postgresql.tar.gz | 1.9 GB | PostgreSQL format |
| chembl_36.h5 | 309 MB | HDF5 format |
| chembl_36.sdf.gz | 893 MB | Structure data format |
| chembl_36_chemreps.txt.gz | 274 MB | Chemical representations |
| schema_documentation.html | 116 KB | HTML schema reference |
| schema_documentation.txt | 95 KB | Text schema reference |
| chembl_36_schema.png | 4.3 MB | Visual schema diagram |

**Download URL:** https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/

---

## Core Entity Relationship

```
MOLECULE_DICTIONARY
       |
       | (molregno)
       v
COMPOUND_RECORDS ----> DOCS (publications)
       |
       | (record_id)
       v
   ACTIVITIES
       |
       | (assay_id)
       v
    ASSAYS
       |
       | (tid)
       v
TARGET_DICTIONARY ----> COMPONENT_SEQUENCES (UniProt)
```

---

## Schema Tables by Domain

### Molecular Data Tables

#### MOLECULE_DICTIONARY

**Description:** Central repository for all chemical compounds

| Column | Type | Description |
|--------|------|-------------|
| molregno | INTEGER | Primary key, internal molecule ID |
| pref_name | VARCHAR(255) | Preferred compound name |
| chembl_id | VARCHAR(20) | ChEMBL identifier (e.g., CHEMBL25) |
| max_phase | NUMERIC | Highest clinical phase (0-4) |
| therapeutic_flag | SMALLINT | Is therapeutic (0/1) |
| dosed_ingredient | SMALLINT | Dosed ingredient flag |
| oral | SMALLINT | Oral bioavailability |
| parenteral | SMALLINT | Parenteral administration |
| topical | SMALLINT | Topical administration |
| black_box_warning | SMALLINT | Black box warning flag |
| natural_product | SMALLINT | **Natural product flag (key field)** |
| first_approval | SMALLINT | Year of first approval |
| indication_class | VARCHAR(1000) | Therapeutic indication |
| withdrawn_flag | SMALLINT | Withdrawn from market |
| chirality | SMALLINT | Chirality type |
| prodrug | SMALLINT | Is prodrug |
| inorganic_flag | SMALLINT | Inorganic compound |
| polymer_flag | SMALLINT | Polymeric compound |
| first_in_class | SMALLINT | First in class drug |
| orphan | SMALLINT | Orphan drug designation |

**Indexes:** PRIMARY KEY (molregno), UNIQUE (chembl_id)

---

#### COMPOUND_STRUCTURES

**Description:** Chemical structure representations

| Column | Type | Description |
|--------|------|-------------|
| molregno | INTEGER | Foreign key to MOLECULE_DICTIONARY |
| molfile | TEXT | MDL MOL file format |
| standard_inchi | TEXT | Standard InChI string |
| standard_inchi_key | VARCHAR(27) | InChI key (14+1+10+1+1) |
| canonical_smiles | TEXT | Canonical SMILES notation |

**Indexes:** PRIMARY KEY (molregno), INDEX (standard_inchi_key)

---

#### COMPOUND_PROPERTIES

**Description:** Physicochemical properties (computed)

| Column | Type | Description |
|--------|------|-------------|
| molregno | INTEGER | Foreign key to MOLECULE_DICTIONARY |
| mw_freebase | NUMERIC(9,2) | Molecular weight (freebase) |
| alogp | NUMERIC(9,2) | Calculated LogP (Wildman-Crippen) |
| hba | SMALLINT | Hydrogen bond acceptors |
| hbd | SMALLINT | Hydrogen bond donors |
| psa | NUMERIC(9,2) | Polar surface area |
| rtb | SMALLINT | Rotatable bonds |
| ro3_pass | VARCHAR(3) | Rule of 3 compliance |
| num_ro5_violations | SMALLINT | Lipinski violations |
| cx_most_apka | NUMERIC(9,2) | Most acidic pKa |
| cx_most_bpka | NUMERIC(9,2) | Most basic pKa |
| cx_logp | NUMERIC(9,2) | Calculated LogP (ChemAxon) |
| cx_logd | NUMERIC(9,2) | Calculated LogD |
| molecular_species | VARCHAR(50) | Species type |
| full_mwt | NUMERIC(9,2) | Full molecular weight |
| aromatic_rings | SMALLINT | Number of aromatic rings |
| heavy_atoms | SMALLINT | Heavy atom count |
| qed_weighted | NUMERIC(3,2) | QED drug-likeness score |
| mw_monoisotopic | NUMERIC(11,4) | Monoisotopic mass |
| full_molformula | VARCHAR(100) | Molecular formula |
| hba_lipinski | SMALLINT | HBA (Lipinski definition) |
| hbd_lipinski | SMALLINT | HBD (Lipinski definition) |
| np_likeness_score | NUMERIC(3,2) | **Natural product likeness score** |

---

#### MOLECULE_SYNONYMS

**Description:** Alternative names and identifiers

| Column | Type | Description |
|--------|------|-------------|
| molregno | INTEGER | Foreign key to MOLECULE_DICTIONARY |
| synonyms | VARCHAR(200) | Alternative name |
| syn_type | VARCHAR(50) | Type (TRADE_NAME, INN, USAN, etc.) |
| res_stem_id | INTEGER | Research stem ID |

---

### Bioactivity Tables

#### ACTIVITIES

**Description:** Core bioactivity measurements

| Column | Type | Description |
|--------|------|-------------|
| activity_id | INTEGER | Primary key |
| assay_id | INTEGER | Foreign key to ASSAYS |
| molregno | INTEGER | Foreign key to MOLECULE_DICTIONARY |
| record_id | INTEGER | Foreign key to COMPOUND_RECORDS |
| standard_type | VARCHAR(250) | Activity type (IC50, Ki, EC50, etc.) |
| standard_relation | VARCHAR(50) | Relation (=, <, >, ~) |
| standard_value | NUMERIC | Standardized value |
| standard_units | VARCHAR(100) | Units (nM, uM, etc.) |
| standard_flag | SMALLINT | Standardization applied |
| standard_upper_value | NUMERIC | Upper bound for ranges |
| pchembl_value | NUMERIC(4,2) | pChEMBL value (-log10) |
| data_validity_comment | VARCHAR(30) | Validity notes |
| potential_duplicate | SMALLINT | Possible duplicate flag |
| text_value | VARCHAR(1000) | Text result (for qualitative data) |
| src_id | SMALLINT | Source ID |
| bao_endpoint | VARCHAR(11) | BioAssay Ontology endpoint |
| activity_comment | VARCHAR(4000) | Activity notes |
| data_validity_description | VARCHAR(500) | Extended validity comment |
| uo_units | VARCHAR(10) | Unit Ontology units |
| qudt_units | VARCHAR(70) | QUDT units |
| doc_id | INTEGER | Document ID |

**Indexes:** PRIMARY KEY (activity_id), INDEX (assay_id), INDEX (molregno), INDEX (record_id)

---

#### ACTIVITY_PROPERTIES

**Description:** Additional activity parameters

| Column | Type | Description |
|--------|------|-------------|
| ap_id | INTEGER | Primary key |
| activity_id | INTEGER | Foreign key to ACTIVITIES |
| type | VARCHAR(250) | Property type |
| relation | VARCHAR(50) | Relation operator |
| value | NUMERIC | Property value |
| units | VARCHAR(100) | Units |
| text_value | VARCHAR(1000) | Text value |
| standard_type | VARCHAR(250) | Standardized type |
| standard_relation | VARCHAR(50) | Standardized relation |
| standard_value | NUMERIC | Standardized value |
| standard_units | VARCHAR(100) | Standardized units |

---

#### LIGAND_EFF

**Description:** Ligand efficiency metrics

| Column | Type | Description |
|--------|------|-------------|
| activity_id | INTEGER | Foreign key to ACTIVITIES |
| bei | NUMERIC(9,2) | Binding efficiency index |
| sei | NUMERIC(9,2) | Surface efficiency index |
| le | NUMERIC(9,2) | Ligand efficiency |
| lle | NUMERIC(9,2) | Lipophilic ligand efficiency |

---

### Assay Tables

#### ASSAYS

**Description:** Experimental assay definitions

| Column | Type | Description |
|--------|------|-------------|
| assay_id | INTEGER | Primary key |
| doc_id | INTEGER | Foreign key to DOCS |
| description | VARCHAR(4000) | Assay description |
| assay_type | VARCHAR(1) | Type: B=Binding, F=Functional, A=ADMET, T=Toxicity, P=Physicochemical, U=Unassigned |
| assay_test_type | VARCHAR(20) | In vivo / In vitro |
| assay_category | VARCHAR(20) | Screening / Confirmatory / etc. |
| assay_organism | VARCHAR(250) | Test organism |
| assay_tax_id | INTEGER | NCBI taxonomy ID |
| assay_strain | VARCHAR(200) | Organism strain |
| assay_tissue | VARCHAR(100) | Tissue type |
| assay_cell_type | VARCHAR(100) | Cell type |
| assay_subcellular_fraction | VARCHAR(100) | Subcellular fraction |
| tid | INTEGER | Foreign key to TARGET_DICTIONARY |
| relationship_type | VARCHAR(1) | D=Direct, H=Homologous, etc. |
| confidence_score | SMALLINT | Target confidence (0-9) |
| curated_by | VARCHAR(32) | Curation source |
| src_id | SMALLINT | Data source ID |
| src_assay_id | VARCHAR(50) | Source assay ID |
| chembl_id | VARCHAR(20) | ChEMBL assay ID |
| bao_format | VARCHAR(11) | BioAssay Ontology format |
| bao_label | VARCHAR(100) | BAO label |
| confidence_description | VARCHAR(250) | Confidence explanation |
| variant_id | INTEGER | Variant sequence ID |

**Assay Types:**
- B = Binding
- F = Functional
- A = ADMET
- T = Toxicity
- P = Physicochemical
- U = Unassigned

**Confidence Scores:**
- 9 = Direct single protein target
- 8 = Homologous single protein target
- 7 = Direct protein complex
- 6 = Homologous protein complex
- 5 = Direct selectivity assay
- 4 = Indirect protein complex
- 3 = Cell-based
- 2 = Tissue/whole organism
- 1 = Uncurated
- 0 = Target not assigned

---

#### ASSAY_PARAMETERS

**Description:** Additional assay conditions

| Column | Type | Description |
|--------|------|-------------|
| assay_param_id | INTEGER | Primary key |
| assay_id | INTEGER | Foreign key to ASSAYS |
| type | VARCHAR(250) | Parameter type |
| relation | VARCHAR(50) | Relation |
| value | NUMERIC | Parameter value |
| units | VARCHAR(100) | Units |
| text_value | VARCHAR(1000) | Text value |
| standard_type | VARCHAR(250) | Standardized type |
| standard_relation | VARCHAR(50) | Standardized relation |
| standard_value | NUMERIC | Standardized value |
| standard_units | VARCHAR(100) | Standardized units |
| standard_text_value | VARCHAR(1000) | Standardized text |
| comments | VARCHAR(4000) | Notes |

---

### Target Tables

#### TARGET_DICTIONARY

**Description:** Biological targets

| Column | Type | Description |
|--------|------|-------------|
| tid | INTEGER | Primary key |
| target_type | VARCHAR(30) | SINGLE PROTEIN, PROTEIN COMPLEX, etc. |
| pref_name | VARCHAR(200) | Preferred name |
| tax_id | INTEGER | NCBI taxonomy ID |
| organism | VARCHAR(150) | Organism name |
| chembl_id | VARCHAR(20) | ChEMBL target ID |
| species_group_flag | SMALLINT | Species group indicator |

**Target Types:**
- SINGLE PROTEIN
- PROTEIN COMPLEX
- PROTEIN FAMILY
- PROTEIN-PROTEIN INTERACTION
- CELL-LINE
- TISSUE
- ORGANISM
- SELECTIVITY GROUP
- NUCLEIC ACID
- SUBCELLULAR
- CHIMERIC PROTEIN
- NO TARGET
- UNKNOWN

---

#### TARGET_COMPONENTS

**Description:** Links targets to protein sequences

| Column | Type | Description |
|--------|------|-------------|
| tid | INTEGER | Foreign key to TARGET_DICTIONARY |
| component_id | INTEGER | Foreign key to COMPONENT_SEQUENCES |
| targcomp_id | INTEGER | Primary key |
| homologue | SMALLINT | Homologue indicator |

---

#### COMPONENT_SEQUENCES

**Description:** Protein and nucleic acid sequences

| Column | Type | Description |
|--------|------|-------------|
| component_id | INTEGER | Primary key |
| component_type | VARCHAR(50) | PROTEIN, DNA, RNA |
| accession | VARCHAR(25) | **UniProt accession** |
| sequence | TEXT | Amino acid sequence |
| sequence_md5sum | VARCHAR(32) | MD5 hash |
| description | VARCHAR(200) | Component description |
| tax_id | INTEGER | NCBI taxonomy ID |
| organism | VARCHAR(150) | Organism |
| db_source | VARCHAR(25) | Sequence source database |
| db_version | VARCHAR(10) | Source database version |

---

#### TARGET_RELATIONS

**Description:** Target-to-target relationships

| Column | Type | Description |
|--------|------|-------------|
| tid | INTEGER | Target ID |
| relationship | VARCHAR(20) | SUBSET, OVERLAP, etc. |
| related_tid | INTEGER | Related target ID |
| targrel_id | INTEGER | Primary key |

---

### Drug and Clinical Tables

#### DRUG_MECHANISM

**Description:** Mechanism of action

| Column | Type | Description |
|--------|------|-------------|
| mec_id | INTEGER | Primary key |
| molregno | INTEGER | Foreign key to MOLECULE_DICTIONARY |
| mechanism_of_action | VARCHAR(250) | MOA description |
| tid | INTEGER | Foreign key to TARGET_DICTIONARY |
| action_type | VARCHAR(50) | INHIBITOR, AGONIST, etc. |
| direct_interaction | SMALLINT | Direct target interaction |
| disease_efficacy | SMALLINT | Disease efficacy evidence |
| site_id | INTEGER | Binding site ID |
| molecular_mechanism | SMALLINT | Molecular mechanism flag |
| mechanism_comment | VARCHAR(2000) | Additional notes |
| selectivity_comment | VARCHAR(1000) | Selectivity notes |
| binding_site_comment | VARCHAR(1000) | Binding site notes |

---

#### DRUG_INDICATION

**Description:** Disease indications

| Column | Type | Description |
|--------|------|-------------|
| drugind_id | INTEGER | Primary key |
| molregno | INTEGER | Foreign key to MOLECULE_DICTIONARY |
| mesh_id | VARCHAR(7) | MeSH disease ID |
| mesh_heading | VARCHAR(200) | MeSH term |
| efo_id | VARCHAR(20) | EFO disease ID |
| efo_term | VARCHAR(200) | EFO term |
| max_phase_for_ind | SMALLINT | Max phase for this indication |

---

#### DRUG_WARNING

**Description:** Safety warnings and withdrawals

| Column | Type | Description |
|--------|------|-------------|
| warning_id | INTEGER | Primary key |
| molregno | INTEGER | Foreign key to MOLECULE_DICTIONARY |
| warning_type | VARCHAR(20) | BLACK_BOX, WITHDRAWN |
| warning_class | VARCHAR(10) | Warning class |
| warning_description | VARCHAR(1000) | Description |
| warning_country | VARCHAR(1000) | Affected countries |
| warning_year | SMALLINT | Year of warning |
| efo_term | VARCHAR(200) | EFO term |
| efo_id | VARCHAR(20) | EFO ID |
| efo_id_for_warning_class | VARCHAR(20) | Warning class EFO |

---

#### ATC_CLASSIFICATION

**Description:** WHO ATC drug classification

| Column | Type | Description |
|--------|------|-------------|
| level5 | VARCHAR(10) | ATC level 5 code |
| level4 | VARCHAR(10) | ATC level 4 code |
| level3 | VARCHAR(10) | ATC level 3 code |
| level2 | VARCHAR(10) | ATC level 2 code |
| level1 | VARCHAR(10) | ATC level 1 code |
| who_name | VARCHAR(200) | WHO drug name |
| level1_description | VARCHAR(150) | Level 1 description |
| level2_description | VARCHAR(150) | Level 2 description |
| level3_description | VARCHAR(200) | Level 3 description |
| level4_description | VARCHAR(200) | Level 4 description |

---

### Reference Tables

#### DOCS

**Description:** Source publications and datasets

| Column | Type | Description |
|--------|------|-------------|
| doc_id | INTEGER | Primary key |
| journal | VARCHAR(50) | Journal name |
| year | SMALLINT | Publication year |
| volume | VARCHAR(50) | Volume |
| issue | VARCHAR(50) | Issue |
| first_page | VARCHAR(50) | First page |
| last_page | VARCHAR(50) | Last page |
| pubmed_id | INTEGER | PubMed ID |
| doi | VARCHAR(100) | DOI |
| chembl_id | VARCHAR(20) | ChEMBL document ID |
| title | VARCHAR(500) | Document title |
| doc_type | VARCHAR(50) | PUBLICATION, PATENT, etc. |
| authors | VARCHAR(4000) | Author list |
| abstract | TEXT | Abstract text |
| patent_id | VARCHAR(20) | Patent ID |
| src_id | SMALLINT | Data source ID |

---

#### SOURCE

**Description:** Data origin

| Column | Type | Description |
|--------|------|-------------|
| src_id | SMALLINT | Primary key |
| src_description | VARCHAR(500) | Source description |
| src_short_name | VARCHAR(20) | Short name |

**Common Sources:**
- Scientific Literature
- PubChem BioAssays
- DrugBank
- BindingDB
- Patent data

---

## Cross-Reference Tables

#### COMPONENT_XREF

**Description:** External database cross-references for components

| Column | Type | Description |
|--------|------|-------------|
| comp_xref_id | INTEGER | Primary key |
| component_id | INTEGER | Foreign key to COMPONENT_SEQUENCES |
| xref_src_db | VARCHAR(50) | Source database |
| xref_id | VARCHAR(300) | External ID |
| xref_name | VARCHAR(200) | External name |

**Supported External Databases:**
- UniProt
- InterPro
- Pfam
- PDB
- GO (Gene Ontology)
- ChEBI
- Reactome

---

## REST API Structure

### Base URL
```
https://www.ebi.ac.uk/chembl/api/data/
```

### Key Endpoints

| Endpoint | Description |
|----------|-------------|
| `/molecule.json` | Molecule records |
| `/molecule/{chembl_id}.json` | Single molecule |
| `/assay.json` | Assay records |
| `/activity.json` | Activity records |
| `/target.json` | Target records |
| `/document.json` | Document records |
| `/status.json` | Database statistics |

### Molecule API Response Schema

```json
{
  "molecule_chembl_id": "CHEMBL25",
  "molecule_type": "Small molecule",
  "max_phase": 4.0,
  "first_approval": 1950,
  "oral": true,
  "parenteral": false,
  "topical": false,
  "black_box_warning": false,
  "natural_product": 0,
  "therapeutic_flag": true,
  "polymer_flag": false,
  "prodrug": false,
  "first_in_class": false,
  "inorganic_flag": false,
  "chirality": -1,
  "withdrawn_flag": false,

  "molecule_properties": {
    "alogp": 1.31,
    "aromatic_rings": 1,
    "cx_logp": 1.24,
    "full_molformula": "C9H8O4",
    "full_mwt": 180.16,
    "hba": 3,
    "hbd": 1,
    "heavy_atoms": 13,
    "mw_freebase": 180.16,
    "np_likeness_score": -0.86,
    "num_ro5_violations": 0,
    "psa": 63.60,
    "qed_weighted": 0.55,
    "ro3_pass": "Y",
    "rtb": 2
  },

  "molecule_structures": {
    "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "standard_inchi": "InChI=1S/C9H8O4/c...",
    "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
  },

  "molecule_synonyms": [
    {"synonym": "Aspirin", "syn_type": "TRADE_NAME"},
    {"synonym": "Acetylsalicylic acid", "syn_type": "INN"}
  ],

  "atc_classifications": [
    "B01AC06", "N02BA01", "A01AD05"
  ],

  "molecule_hierarchy": {
    "active_chembl_id": "CHEMBL25",
    "parent_chembl_id": "CHEMBL25"
  }
}
```

### Activity API Response Schema

```json
{
  "activity_id": 12345,
  "molecule_chembl_id": "CHEMBL25",
  "assay_chembl_id": "CHEMBL1234567",
  "target_chembl_id": "CHEMBL220",
  "document_chembl_id": "CHEMBL1123456",

  "assay_description": "Inhibition of cyclooxygenase-1",
  "assay_type": "B",

  "target_pref_name": "Cyclooxygenase-1",
  "target_organism": "Homo sapiens",
  "target_tax_id": 9606,

  "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",

  "standard_type": "IC50",
  "standard_relation": "=",
  "standard_value": 300.0,
  "standard_units": "nM",
  "pchembl_value": 6.52,

  "ligand_efficiency": {
    "bei": 15.31,
    "le": 0.51,
    "lle": 5.21,
    "sei": 8.27
  },

  "document_journal": "J Med Chem",
  "document_year": 1995,

  "data_validity_comment": null,
  "potential_duplicate": false
}
```

### Target API Response Schema

```json
{
  "target_chembl_id": "CHEMBL220",
  "pref_name": "Cyclooxygenase-1",
  "target_type": "SINGLE PROTEIN",
  "organism": "Homo sapiens",
  "tax_id": 9606,
  "species_group_flag": false,

  "target_components": [{
    "component_id": 220,
    "component_type": "PROTEIN",
    "accession": "P23219",
    "component_description": "Prostaglandin G/H synthase 1",
    "relationship": "SINGLE PROTEIN",

    "target_component_synonyms": [
      {"synonym": "COX-1", "syn_type": "GENE_SYMBOL"},
      {"synonym": "PTGS1", "syn_type": "GENE_SYMBOL"}
    ],

    "target_component_xrefs": [
      {"xref_id": "P23219", "xref_src_db": "UniProt"},
      {"xref_id": "IPR001128", "xref_src_db": "InterPro"},
      {"xref_id": "PF00067", "xref_src_db": "Pfam"}
    ]
  }],

  "cross_references": [
    {"xref_id": "P23219", "xref_src_db": "UniProt"}
  ]
}
```

---

## Key Foreign Key Relationships

| From Table | From Column | To Table | To Column |
|------------|-------------|----------|-----------|
| COMPOUND_STRUCTURES | molregno | MOLECULE_DICTIONARY | molregno |
| COMPOUND_PROPERTIES | molregno | MOLECULE_DICTIONARY | molregno |
| COMPOUND_RECORDS | molregno | MOLECULE_DICTIONARY | molregno |
| COMPOUND_RECORDS | doc_id | DOCS | doc_id |
| ACTIVITIES | molregno | MOLECULE_DICTIONARY | molregno |
| ACTIVITIES | record_id | COMPOUND_RECORDS | record_id |
| ACTIVITIES | assay_id | ASSAYS | assay_id |
| ACTIVITIES | doc_id | DOCS | doc_id |
| ASSAYS | tid | TARGET_DICTIONARY | tid |
| ASSAYS | doc_id | DOCS | doc_id |
| TARGET_COMPONENTS | tid | TARGET_DICTIONARY | tid |
| TARGET_COMPONENTS | component_id | COMPONENT_SEQUENCES | component_id |
| DRUG_MECHANISM | molregno | MOLECULE_DICTIONARY | molregno |
| DRUG_MECHANISM | tid | TARGET_DICTIONARY | tid |
| DRUG_INDICATION | molregno | MOLECULE_DICTIONARY | molregno |
| MOLECULE_SYNONYMS | molregno | MOLECULE_DICTIONARY | molregno |

---

## Natural Product Query Examples

### Find All Natural Products

```sql
SELECT md.chembl_id, md.pref_name, cs.canonical_smiles, cp.np_likeness_score
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
LEFT JOIN compound_properties cp ON md.molregno = cp.molregno
WHERE md.natural_product = 1;
```

### Natural Products with Bioactivity

```sql
SELECT DISTINCT md.chembl_id, md.pref_name,
       a.standard_type, a.standard_value, a.standard_units,
       td.pref_name AS target_name
FROM molecule_dictionary md
JOIN activities a ON md.molregno = a.molregno
JOIN assays ay ON a.assay_id = ay.assay_id
JOIN target_dictionary td ON ay.tid = td.tid
WHERE md.natural_product = 1
AND a.standard_type IN ('IC50', 'Ki', 'EC50')
AND a.standard_value IS NOT NULL;
```

### Find Compounds by NP-Likeness Score

```sql
SELECT md.chembl_id, md.pref_name, cp.np_likeness_score
FROM molecule_dictionary md
JOIN compound_properties cp ON md.molregno = cp.molregno
WHERE cp.np_likeness_score > 1.5
ORDER BY cp.np_likeness_score DESC
LIMIT 100;
```

---

## Integration Notes for Knowledge Base

### Key Identifiers for Cross-Referencing

| Identifier | Format | Example |
|------------|--------|---------|
| ChEMBL Compound | CHEMBL + digits | CHEMBL25 |
| ChEMBL Target | CHEMBL + digits | CHEMBL220 |
| ChEMBL Assay | CHEMBL + digits | CHEMBL1234567 |
| InChI Key | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| UniProt Accession | 6 or 10 char | P23219 |
| PubChem CID | Integer | 2244 |
| PubMed ID | Integer | 12345678 |

### Data Quality Indicators

- **pchembl_value**: Standardized -log10 activity, NULL if cannot be computed
- **data_validity_comment**: Flags suspicious data
- **confidence_score**: Target assignment confidence (0-9, higher = better)
- **potential_duplicate**: Indicates possible duplicate entry

### Licensing

**License:** CC BY-SA 3.0
**Attribution Required:** Yes
**Commercial Use:** Allowed with share-alike

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | SQLite, MySQL, PostgreSQL |
| Alternative | SDF (.sdf.gz), HDF5 (.h5), TSV |
| Compression | gzip (.gz), tar.gz |
| Encoding | UTF-8 |
| API Response | JSON, XML |

---

## Sample Data

### Example Record
```json
{
  "molregno": 2244,
  "chembl_id": "CHEMBL25",
  "pref_name": "Aspirin",
  "activity_id": 12345,
  "target_id": "CHEMBL220",
  "standard_type": "IC50",
  "standard_value": 300.0,
  "standard_units": "nM",
  "pchembl_value": 6.52
}
```

### Sample Query Result
| chembl_id | pref_name | standard_type | standard_value | target_name | pchembl_value |
|-----------|-----------|---------------|----------------|-------------|---------------|
| CHEMBL25 | Aspirin | IC50 | 300.0 | Cyclooxygenase-1 | 6.52 |
| CHEMBL210 | Ibuprofen | Ki | 150.0 | COX1 | 6.82 |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| ChEMBL | FTP | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/ |
| SQLite | FTP | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_36_sqlite.tar.gz |
| PostgreSQL | FTP | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_36_postgresql.tar.gz |
| SDF Structures | FTP | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_36.sdf.gz |

**Access Requirements:** Open access, no registration required

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| ChEMBL | CC BY-SA 3.0 | Yes (with share-alike) |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `molregno` | Internal ChEMBL molecule registration number, primary key for compound identification | 1234567 |
| `chembl_id` | Public ChEMBL identifier for molecules, targets, assays, or documents | CHEMBL25 |
| `pchembl_value` | Standardized activity value as -log10(molar IC50/Ki/EC50), enables comparison across assays | 6.52 (-log10 of 300nM) |
| `max_phase` | Highest clinical development phase reached (0=preclinical, 1-3=trials, 4=approved) | 4 |
| `standard_inchi_key` | 27-character hash of InChI for structure lookup and deduplication | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| `canonical_smiles` | Standardized SMILES notation for molecular structure representation | CC(=O)Oc1ccccc1C(=O)O |
| `tid` | Target identifier, internal primary key in TARGET_DICTIONARY | 100126 |
| `assay_id` | Assay identifier, internal primary key in ASSAYS table | 1234567 |
| `confidence_score` | Target assignment confidence (0-9), higher values indicate more reliable target identification | 9 (direct single protein) |
| `natural_product` | Binary flag (0/1) indicating whether compound is derived from natural sources | 1 |
| `np_likeness_score` | Computed score indicating structural similarity to known natural products | 1.5 (NP-like) |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Bioactivity | Measurable effect of a compound on a biological target or system | ACTIVITIES table |
| IC50 | Inhibitor concentration producing 50% inhibition of target activity | standard_type field |
| Ki | Inhibition constant, binding affinity of inhibitor for target | standard_type field |
| EC50 | Effective concentration producing 50% of maximum response | standard_type field |
| Ligand Efficiency (LE) | Binding energy per heavy atom, measures optimization efficiency | LIGAND_EFF table |
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity properties | assay_type = 'A' |
| Prodrug | Inactive compound converted to active drug in vivo | prodrug column |
| Lipinski Rule of 5 | Drug-likeness criteria: MW<500, LogP<5, HBD<5, HBA<10 | num_ro5_violations |
| QED | Quantitative Estimate of Drug-likeness, weighted score 0-1 | qed_weighted column |
| ATC Classification | WHO Anatomical Therapeutic Chemical drug classification system | ATC_CLASSIFICATION table |
| Mechanism of Action (MOA) | Description of how drug produces pharmacological effect | DRUG_MECHANISM table |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ChEMBL | Chemical database of bioactive molecules with drug-like properties | European Bioinformatics Institute |
| SMILES | Simplified Molecular Input Line Entry System | Chemical structure notation |
| InChI | International Chemical Identifier | IUPAC standard |
| HBA | Hydrogen Bond Acceptor | Lipinski parameter |
| HBD | Hydrogen Bond Donor | Lipinski parameter |
| PSA | Polar Surface Area | Topological descriptor (Angstrom^2) |
| RTB | Rotatable Bonds | Flexibility measure |
| MW | Molecular Weight | Daltons |
| LogP | Partition Coefficient | Lipophilicity measure |
| pKa | Acid Dissociation Constant | -log10 of Ka |
| BAO | BioAssay Ontology | Assay classification ontology |
| MOA | Mechanism of Action | Drug-target interaction type |
| INN | International Nonproprietary Name | WHO generic drug name |
| USAN | United States Adopted Name | US generic drug name |
| EFO | Experimental Factor Ontology | Disease classification |
| MeSH | Medical Subject Headings | NLM controlled vocabulary |

---

## References

1. Zdrazil B, et al. (2024) "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods." Nucleic Acids Res. 52(D1):D1180-D1192.

2. ChEMBL Web Services: https://www.ebi.ac.uk/chembl/api/data/docs

3. Schema Documentation: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 24,267,312 |
| Storage | 5.2 GB (compressed SQLite) |
| Last updated | January 2026 |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation from ChEMBL 36 |
