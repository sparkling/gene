---
id: schema-binding-affinity
title: "Binding Affinity Database Schema"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, binding-affinity, ttd, bindingdb, gtopdb, drug-discovery]
---

**Parent:** [Schema Documentation](./_index.md)

# Binding Affinity Database Schema

**Document ID:** BINDING-AFFINITY-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Scope:** TTD, BindingDB, GtoPdb

---

## TL;DR

This document covers three major binding affinity databases: TTD (Therapeutic Target Database) with 3,131 targets and 40,000+ drugs organized across 3 druggability perspectives; BindingDB with 2.9M binding measurements (IC50, Ki, EC50, Kd) for 1.3M compounds; and GtoPdb (Guide to Pharmacology) with 170,000+ curated interactions across 10 target families. All three provide REST APIs and bulk downloads with extensive cross-references to ChEMBL, PubChem, and UniProt.

---

## 1. TTD (Therapeutic Target Database)

**URL:** https://idrblab.net/ttd/ (or https://db.idrblab.net/ttd/)
**Schema URL:** https://db.idrblab.net/ttd/schema
**Downloads:** https://idrblab.net/ttd/full-data-download
**Reference:** TTD 2024 (NAR 2024) - DOI: 10.1093/nar/gkad751
**License:** Free for academic use

### 1.1 Database Statistics (2024)

| Entity | Count |
|--------|-------|
| **Successful Targets** | 426 |
| **Clinical Trial Targets** | 1,014 |
| **Preclinical/Patented Targets** | 212 |
| **Literature-reported Targets** | 1,479 |
| **Total Targets** | 3,131 |
| **Approved Drugs** | 2,895 |
| **Clinical Trial Drugs** | 11,796 |
| **Preclinical/Patented Drugs** | 5,041 |
| **Experimental Drugs** | 20,130 |
| **Total Drugs** | 39,862 |

### 1.2 Target Classification Schema

#### Three Druggability Perspectives

TTD organizes druggability data across three perspectives for comprehensive target assessment:

**Perspective 1: Molecular Interactions/Regulations**

| Data Category | Description | Data Points |
|---------------|-------------|-------------|
| Ligand binding pocket | 3D structure of drug binding site | PDB structures |
| Protein-protein interactions | Network properties | Interaction partners |
| Microbiota-drug regulation | Gut microbiome interactions | Regulation data |
| Post-translational modifications | Target modifications | PTM annotations |

**Perspective 2: Human System Profiles**

| Data Category | Description | Data Points |
|---------------|-------------|-------------|
| Protein similarity | Similarity to non-family proteins | % identity scores |
| Pathway involvement | Life-essential pathways | 241 pathways |
| Tissue distribution | Human tissue expression | 32 tissues |
| Disease associations | Target-disease links | ICD-11 codes |

**Perspective 3: Cell-Based Expression Variations**

| Data Category | Description | Data Points |
|---------------|-------------|-------------|
| Cell type expression | Single-cell expression data | 1,742 cell types |
| Exogenous stimuli | Environmental factor responses | 625 factors |
| Endogenous factors | Internal regulatory factors | 447 factors |
| Expression variation | Expression level changes | Fold change values |

#### Target Status Categories

| Status | Definition | Count (2024) |
|--------|------------|--------------|
| `Successful` | Targeted by at least one FDA-approved drug | 426 |
| `Clinical Trial` | Not approved, but in clinical trials | 1,014 |
| `Preclinical/Patented` | Not in trials, but preclinical/patented | 212 |
| `Literature-reported` | Research targets with experimental drugs only | 1,479 |

#### Target Validation Status Levels

| Level | Description | Evidence Required |
|-------|-------------|-------------------|
| `Validated` | Clinically validated target | Approved drug exists |
| `Advanced Clinical` | Phase 3 clinical trials | Phase 3 data |
| `Clinical` | Phase 1-2 clinical trials | Clinical trial data |
| `Preclinical` | Animal model validation | In vivo efficacy |
| `In vitro` | Cell-based validation | Cell assay data |
| `Predicted` | Computational prediction | Bioinformatics evidence |

#### Target Molecular Types

| Type | Description | Examples |
|------|-------------|----------|
| `Protein` | Most common target type | Kinases, GPCRs, enzymes |
| `Nucleic Acid` | DNA, mRNA, miRNA, lncRNA targets | Antisense targets |
| `Other Molecule` | Non-protein biomolecules | Uric acid, iron, ROS |

### 1.3 Drug Classification Schema

#### Drug Status Categories

| Status | Definition | Count (2024) |
|--------|------------|--------------|
| `Approved` | FDA/regulatory approved | 2,895 |
| `Clinical Trial` | In clinical development | 11,796 |
| `Preclinical/Patented` | Preclinical or patented | 5,041 |
| `Experimental` | Research compounds | 20,130 |

#### Seven Drug Types

| Type | Description | Examples |
|------|-------------|----------|
| `Small Molecule` | Traditional chemical drugs | Imatinib, aspirin |
| `Antibody` | mAbs, ADCs, bispecifics, IgG mixtures | Trastuzumab, adalimumab |
| `Peptide` | Therapeutic peptides | Insulin, GLP-1 agonists |
| `Nucleic Acid Drug` | ASOs, siRNAs, saRNAs, miRNAs, mRNAs | Patisiran, nusinersen |
| `Cell Therapy` | CAR-T, stem cells | Tisagenlecleucel |
| `Gene Therapy` | Gene delivery vectors | Zolgensma |
| `Vaccine` | Preventive and therapeutic vaccines | mRNA vaccines |

### 1.4 Data Model Structure

#### Target Data Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `TTD_TARGET_ID` | string | Internal target identifier | `TTDT00001` |
| `Target_Name` | string | Standard target name | `B-Raf proto-oncogene` |
| `Target_Type` | string | Molecular type | `Protein` |
| `Target_Status` | string | Development status | `Successful Target` |
| `UniProt_ID` | string | UniProt accession | `P15056` |
| `Gene_Name` | string | Gene symbol | `BRAF` |
| `PDB_ID` | string | PDB structure IDs (semicolon-delimited) | `4MNE;4MNF;5HI2` |
| `KEGG_Pathway` | string | KEGG pathway mappings | `hsa04010:MAPK` |
| `Biochemical_Class` | string | Protein family | `Kinase` |
| `Sequence` | text | Amino acid sequence | FASTA format |
| `EC_Number` | string | Enzyme Commission number | `2.7.11.1` |
| `Target_Function` | text | Functional description | Free text |
| `Disease_Association` | string | Associated diseases | ICD-11 codes |

#### Drug Data Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `TTD_DRUG_ID` | string | Internal drug identifier | `D0A9YA` |
| `Drug_Name` | string | Standard drug name | `Vemurafenib` |
| `Drug_Status` | string | Development status | `Approved` |
| `Drug_Type` | string | Drug modality | `Small Molecule` |
| `Therapeutic_Class` | string | ATC classification | `L01EC01` |
| `CAS_Number` | string | CAS registry number | `918504-65-1` |
| `PubChem_CID` | integer | PubChem compound ID | `42611257` |
| `DrugBank_ID` | string | DrugBank identifier | `DB08881` |
| `ChEMBL_ID` | string | ChEMBL identifier | `CHEMBL1667` |
| `SMILES` | string | Chemical structure | `CCCS...` |
| `InChIKey` | string | InChIKey (27 chars) | `JOBWBQQXX...` |
| `Molecular_Weight` | decimal | Molecular weight (Da) | `489.92` |
| `LogP` | decimal | Partition coefficient | `5.2` |

#### Target-Drug Relationship Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `TTD_TARGET_ID` | string | Target identifier | `TTDT00001` |
| `TTD_DRUG_ID` | string | Drug identifier | `D0A9YA` |
| `Activity_Type` | string | Mechanism type | `Inhibitor` |
| `MOA` | string | Mode of action | `Kinase inhibition` |
| `Indication` | string | Disease indication | `Melanoma` |
| `ICD_Code` | string | ICD-11 disease code | `2B90.0` |
| `Clinical_Status` | string | Trial phase | `Approved` |
| `Highest_Status` | string | Highest development status | `Phase 4` |
| `Reference` | string | PubMed IDs (semicolon-delimited) | `20823850;21639808` |

### 1.5 Download File Structure

TTD provides multiple download files in tab-delimited format:

| File | Content | Size (approx.) |
|------|---------|----------------|
| `P1-01-TTD_target_download.txt` | Target information | 15 MB |
| `P1-02-TTD_drug_download.txt` | Drug information | 25 MB |
| `P1-03-TTD_crossmatching.txt` | ID cross-references | 8 MB |
| `P1-04-Drug_disease.txt` | Drug-disease mappings | 12 MB |
| `P1-05-Target_disease.txt` | Target-disease mappings | 5 MB |
| `P1-06-Target_pathway.txt` | Pathway annotations | 3 MB |
| `P2-01-TTD_uniprot_all.txt` | UniProt mappings | 2 MB |
| `P2-02-TTD_pubchem_drug.txt` | PubChem mappings | 10 MB |

#### Sample Target Download Format (P1-01)

```
TTDTRTID	TTD Target ID	Field Name	Content
TTDT00001	T00001	TARGETNAME	B-Raf proto-oncogene serine/threonine kinase
TTDT00001	T00001	TARGETYPE	Successful target
TTDT00001	T00001	BIOCLASS	Kinase
TTDT00001	T00001	UNIPROID	P15056
TTDT00001	T00001	GENENAME	BRAF
TTDT00001	T00001	PDBIDS	4MNE;4MNF;5HI2;6UAN;7JVP
```

#### Sample Drug Download Format (P1-02)

```
TTDDRGID	TTD Drug ID	Field Name	Content
D0A9YA	D0A9YA	DRUGNAME	Vemurafenib
D0A9YA	D0A9YA	TRADNAME	Zelboraf
D0A9YA	D0A9YA	DRUGTYPE	Small molecule drug
D0A9YA	D0A9YA	DRUGSTAT	Approved
D0A9YA	D0A9YA	DRUGSMIL	CCCS(=O)(=O)Nc1ccc...
D0A9YA	D0A9YA	DRUGCOMP	C23H18ClF2N3O3S
```

### 1.6 ID Systems and Cross-References

| Entity | ID Format | Example |
|--------|-----------|---------|
| Target | TTD Target ID | `TTDT00001` |
| Drug | TTD Drug ID | `D0A9YA` |
| Disease | ICD-11 Code | `2B90.0` |
| UniProt | UniProt Accession | `P15056` |
| PubChem | PubChem CID | `42611257` |
| DrugBank | DrugBank ID | `DB08881` |
| ChEMBL | ChEMBL ID | `CHEMBL1667` |
| KEGG | KEGG Gene ID | `hsa:673` |
| PDB | PDB ID | `4MNE` |

---

## 2. BindingDB

**URL:** https://www.bindingdb.org
**API Documentation:** https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp
**Downloads:** https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp
**Reference:** BindingDB 2024 (NAR 2025) - DOI: 10.1093/nar/gkae1039
**License:** CC BY 3.0

### 2.1 Database Statistics (2024)

| Entity | Count |
|--------|-------|
| **Total Binding Measurements** | 2,900,000+ |
| **Distinct Compounds** | 1,300,000+ |
| **Protein Targets** | ~9,400 |
| **IC50 Measurements** | 1,800,000+ |
| **Ki Measurements** | 560,000+ |
| **EC50 Measurements** | 220,000+ |
| **Kd Measurements** | 100,000+ |
| **Publications** | 40,000+ |
| **Patents** | 35,000+ |

### 2.2 Four Affinity Measurement Types

| Type | Full Name | Count | Description | Typical Range |
|------|-----------|-------|-------------|---------------|
| `IC50` | Half-maximal inhibitory concentration | 1.8M | Concentration causing 50% inhibition | 1 nM - 100 uM |
| `Ki` | Inhibition constant | 560K | Inhibitor binding constant | 0.1 nM - 10 uM |
| `EC50` | Half-maximal effective concentration | 220K | Concentration for 50% effect | 1 nM - 100 uM |
| `Kd` | Dissociation constant | 100K | Binding affinity measure | 0.01 nM - 1 uM |

#### Additional Measurement Types

| Type | Description | Units |
|------|-------------|-------|
| `kon` | Association rate constant | M^-1 s^-1 |
| `koff` | Dissociation rate constant | s^-1 |
| `deltaG` | Free energy of binding | kcal/mol |
| `deltaH` | Enthalpy of binding | kcal/mol |
| `-TdeltaS` | Entropy term | kcal/mol |
| `pH` | Measurement pH | pH units |
| `Temperature` | Measurement temperature | Celsius |

### 2.3 REST API Endpoints

**Base URL:** `https://www.bindingdb.org/rest/` or `https://bindingdb.org/rest/`

#### getLigandsByPDBs

Retrieve binding data for PDB targets.

```
GET /getLigandsByPDBs?pdb={pdb_ids}&cutoff={affinity_nm}&identity={percent}&response={format}
```

| Parameter | Type | Description | Required | Example |
|-----------|------|-------------|----------|---------|
| `pdb` | string | PDB ID(s), comma-separated | Yes | `1Q0L,3ANM` |
| `cutoff` | integer | Affinity threshold (nM) | Yes | `100` |
| `identity` | integer | Sequence identity cutoff (%) | Yes | `92` |
| `response` | string | Response format | No | `application/json` |

**Example Request:**
```
https://bindingdb.org/rest/getLigandsByPDBs?pdb=1Q0L,3ANM&cutoff=100&identity=92&response=application/json
```

#### getLigandsByUniprots

Retrieve binding data for UniProt targets.

```
GET /getLigandsByUniprots?uniprot={uniprot_ids}&cutoff={affinity_nm}&code={filter}&response={format}
```

| Parameter | Type | Description | Required | Example |
|-----------|------|-------------|----------|---------|
| `uniprot` | string | UniProt ID(s), comma-separated | Yes | `P00176,P00183` |
| `cutoff` | integer | Affinity threshold (nM) | Yes | `10000` |
| `code` | integer | 0=all, 1=commercial, 2=FDA approved | No | `0` |
| `response` | string | Response format | No | `application/json` |

**Example Request:**
```
https://bindingdb.org/rest/getLigandsByUniprots?uniprot=P35355&cutoff=10000&response=application/json
```

#### getLigandsByUniprot (Single)

```
GET /getLigandsByUniprot?uniprot={uniprot_id};{cutoff}&response={format}
```

**Example Request:**
```
https://bindingdb.org/rest/getLigandsByUniprot?uniprot=P35355;100&response=application/json
```

#### getTargetByCompound

Find targets for a compound by structure similarity.

```
GET /getTargetByCompound?smiles={smiles}&cutoff={similarity}&response={format}
```

| Parameter | Type | Description | Required | Example |
|-----------|------|-------------|----------|---------|
| `smiles` | string | SMILES structure (URL-encoded) | Yes | `CC(C)Cc1ccc(cc1)...` |
| `cutoff` | decimal | Similarity threshold (0-1) | Yes | `0.85` |
| `response` | string | Response format | No | `application/json` |

### 2.4 API Response Schema (30+ Fields)

#### JSON Response Structure

```json
{
  "affinities": [
    {
      "monomerid": "50000001",
      "smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
      "inchi": "InChI=1S/C13H18O2/c1-9...",
      "inchikey": "HEFNNWSXXWATRW-JTQLQIEISA-N",
      "affinity_type": "IC50",
      "affinity_value": "100",
      "affinity_relation": "=",
      "affinity_unit": "nM",
      "target_name": "Cyclooxygenase-2",
      "target_source": "Homo sapiens",
      "uniprot_id": "P35354",
      "pdb_id": "1CX2",
      "pubmed_id": "12345678",
      "doi": "10.1021/jm000001",
      "patent": null,
      "ki_value": null,
      "kd_value": null,
      "ec50_value": null,
      "kon_value": null,
      "koff_value": null,
      "ph": "7.4",
      "temperature": "25",
      "assay_type": "Enzyme activity",
      "assay_description": "Inhibition of COX-2 enzyme activity",
      "ligand_name": "Ibuprofen",
      "cas_number": "15687-27-1",
      "chembl_id": "CHEMBL521",
      "pubchem_cid": "3672",
      "drugbank_id": "DB01050",
      "zinc_id": "ZINC000000001532",
      "institution": "University Example",
      "authors": "Smith J, et al.",
      "article_doi": "10.1021/jm000001",
      "entry_doi": "10.37126/bdb..."
    }
  ]
}
```

#### Complete Response Fields Dictionary

| Field | Type | Description | Always Present |
|-------|------|-------------|----------------|
| `monomerid` | string | BindingDB compound ID | Yes |
| `smiles` | string | SMILES structure | Yes |
| `inchi` | string | InChI string | Yes |
| `inchikey` | string | InChIKey (27 chars) | Yes |
| `affinity_type` | string | IC50, Ki, Kd, or EC50 | Yes |
| `affinity_value` | string | Numeric value | Yes |
| `affinity_relation` | string | =, <, >, ~ | Yes |
| `affinity_unit` | string | nM (nanomolar) | Yes |
| `target_name` | string | Protein target name | Yes |
| `target_source` | string | Species | Yes |
| `uniprot_id` | string | UniProt accession | No |
| `pdb_id` | string | PDB structure ID | No |
| `pubmed_id` | string | PubMed reference | No |
| `doi` | string | Article DOI | No |
| `patent` | string | Patent number | No |
| `ki_value` | string | Ki value (nM) | No |
| `kd_value` | string | Kd value (nM) | No |
| `ec50_value` | string | EC50 value (nM) | No |
| `ic50_value` | string | IC50 value (nM) | No |
| `kon_value` | string | Association rate | No |
| `koff_value` | string | Dissociation rate | No |
| `delta_g` | string | Free energy (kcal/mol) | No |
| `delta_h` | string | Enthalpy (kcal/mol) | No |
| `tds` | string | Entropy term (kcal/mol) | No |
| `ph` | string | Measurement pH | No |
| `temperature` | string | Temperature (C) | No |
| `assay_type` | string | Assay category | No |
| `assay_description` | string | Assay details | No |
| `ligand_name` | string | Compound name | No |
| `cas_number` | string | CAS registry number | No |
| `chembl_id` | string | ChEMBL ID | No |
| `pubchem_cid` | string | PubChem compound ID | No |
| `pubchem_sid` | string | PubChem substance ID | No |
| `drugbank_id` | string | DrugBank ID | No |
| `chebi_id` | string | ChEBI ID | No |
| `zinc_id` | string | ZINC ID | No |
| `institution` | string | Research institution | No |
| `authors` | string | Author list | No |
| `article_doi` | string | Publication DOI | No |
| `entry_doi` | string | BindingDB entry DOI | No |

### 2.5 TSV Download Format

Main TSV files contain one row per binding measurement:

| Column | Description | Example |
|--------|-------------|---------|
| `BindingDB MonomerID` | Compound identifier | `50000001` |
| `Ligand SMILES` | Chemical structure | `CC(C)Cc1ccc...` |
| `Ligand InChI` | InChI notation | `InChI=1S/...` |
| `Ligand InChI Key` | Hashed identifier | `HEFNNWSXX...` |
| `BindingDB Ligand Name` | Compound name | `Ibuprofen` |
| `Target Name Assigned by Curator` | Curated target name | `Cyclooxygenase-2` |
| `Target Source Organism` | Species | `Homo sapiens` |
| `Ki (nM)` | Ki value | `150` |
| `IC50 (nM)` | IC50 value | `320` |
| `Kd (nM)` | Kd value | `null` |
| `EC50 (nM)` | EC50 value | `null` |
| `kon (M-1 s-1)` | Association rate | `null` |
| `koff (s-1)` | Dissociation rate | `null` |
| `pH` | pH of measurement | `7.4` |
| `Temp (C)` | Temperature | `25` |
| `Curation/DataSource` | Data source | `ChEMBL` |
| `Article DOI` | Publication DOI | `10.1021/jm...` |
| `PMID` | PubMed ID | `12345678` |
| `Patent Number` | Patent reference | `US8012345` |
| `Authors` | Author list | `Smith J, et al.` |
| `Institution` | Research institution | `MIT` |
| `BindingDB Entry DOI` | Entry DOI | `10.37126/bdb...` |
| `Link to Ligand in BindingDB` | Ligand URL | `https://...` |
| `Link to Target in BindingDB` | Target URL | `https://...` |
| `Link to PDB` | PDB structure link | `https://...` |
| `UniProt (SwissProt) Entry Name` | UniProt name | `PGH2_HUMAN` |
| `UniProt (SwissProt) Primary ID` | UniProt accession | `P35354` |
| `UniProt (TrEMBL) Entry Name` | TrEMBL name | `null` |
| `UniProt (TrEMBL) Primary ID` | TrEMBL accession | `null` |
| `PubChem CID` | PubChem compound ID | `3672` |
| `PubChem SID` | PubChem substance ID | `123456` |
| `ChEBI ID` | ChEBI identifier | `CHEBI:5855` |
| `ChEMBL ID` | ChEMBL compound ID | `CHEMBL521` |
| `DrugBank ID` | DrugBank identifier | `DB01050` |
| `ZINC ID` | ZINC identifier | `ZINC000000001532` |
| `Number of Protein Chains in Target` | Chain count | `2` |
| `BindingDB Target Chain Sequence` | Protein sequence | `MSLKFLA...` |

#### Sample TSV Data

```tsv
BindingDB MonomerID	Ligand SMILES	Ligand InChI Key	Ki (nM)	IC50 (nM)	Target Name	UniProt ID	PMID
50000001	CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O	HEFNNWSXXWATRW	150	320	Cyclooxygenase-2	P35354	12345678
50000002	Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C	KTUFNOKKBVMGRW	1.2	3.5	ABL1 kinase	P00519	21639808
```

### 2.6 Download Files

| File | Description | Format |
|------|-------------|--------|
| `BindingDB_All.tsv` | Complete dataset | TSV |
| `BindingDB_All_2D.sdf` | 2D coordinates | SDF |
| `BindingDB_All_3D.sdf` | 3D coordinates (Vconf) | SDF |
| `BindingDB_All_terse.tsv` | One compound per row | TSV |
| `BindingDB_TargetSequences.fasta` | Target sequences | FASTA |
| `BindingDB_CID.txt` | PubChem CID mapping | TXT |
| `BindingDB_SID.txt` | PubChem SID mapping | TXT |
| `BindingDB_CHEBI_ID.txt` | ChEBI mapping | TXT |
| `BindingDB_DrugBankID.txt` | DrugBank mapping | TXT |
| `BindingDB_PubMed.txt` | PubMed ID collection | TXT |
| `BindingDB_UniProt.txt` | UniProt mapping | TXT |
| `BDB_Assays.tsv` | Assay descriptions | TSV |
| `BDB_rsid_eaids.txt` | Reaction Set mapping | TXT |

### 2.7 ID Systems

| Entity | ID Format | Example |
|--------|-----------|---------|
| Compound (Monomer) | Numeric MonomerID | `50000001` |
| Target (Polymer) | Numeric PolymerID | `3000` |
| Entry | DOI | `10.37126/bdb...` |
| PubChem CID | Numeric | `3672` |
| ChEMBL ID | CHEMBL prefix | `CHEMBL521` |
| DrugBank ID | DB prefix | `DB01050` |
| UniProt | Accession | `P35354` |

---

## 3. GtoPdb (Guide to Pharmacology)

**URL:** https://www.guidetopharmacology.org
**API Base:** https://www.guidetopharmacology.org/services/
**Documentation:** https://www.guidetopharmacology.org/webServices.jsp
**Downloads:** https://www.guidetopharmacology.org/download.jsp
**Reference:** GtoPdb 2024 (NAR 2024) - DOI: 10.1093/nar/gkad944
**License:** CC BY-SA 4.0

### 3.1 Database Statistics (2024)

| Entity | Count |
|--------|-------|
| **Total Targets** | ~3,000 |
| **Target Families** | 800+ |
| **Ligands** | ~13,000 |
| **Approved Drugs** | ~2,200 |
| **Interactions** | 170,000+ |
| **References** | 85,000+ |
| **Diseases** | 1,200+ |

### 3.2 Ten Target Family Types

| Type | Code | Description | Count |
|------|------|-------------|-------|
| G protein-coupled receptors | `GPCR` | 7TM receptors | 400+ |
| Nuclear hormone receptors | `NHR` | Ligand-activated TFs | 48 |
| Ligand-gated ion channels | `LGIC` | Cys-loop, glutamate, P2X | 80+ |
| Voltage-gated ion channels | `VGIC` | Na+, K+, Ca2+ channels | 140+ |
| Other ion channels | `OtherIC` | TRP, CLC, connexins | 60+ |
| Enzymes | `Enzyme` | Kinases, proteases, etc. | 600+ |
| Catalytic receptors | `CatalyticReceptor` | RTKs, RSKs, cytokine receptors | 200+ |
| Transporters | `Transporter` | SLC, ABC families | 400+ |
| Other proteins | `OtherProtein` | Various non-receptor proteins | 300+ |
| Accessory proteins | `AccessoryProtein` | Auxiliary subunits | 100+ |

#### Target Family Hierarchy

```
Target Families
├── G protein-coupled receptors
│   ├── Class A (Rhodopsin)
│   ├── Class B (Secretin)
│   ├── Class C (Glutamate)
│   ├── Class F (Frizzled)
│   └── Adhesion GPCRs
├── Ion channels
│   ├── Ligand-gated
│   │   ├── Cys-loop receptors
│   │   ├── Glutamate receptors
│   │   └── P2X receptors
│   ├── Voltage-gated
│   │   ├── Sodium channels
│   │   ├── Potassium channels
│   │   └── Calcium channels
│   └── Other channels
├── Nuclear hormone receptors
│   ├── Thyroid hormone receptors
│   ├── Retinoid receptors
│   ├── Steroid receptors
│   └── Orphan receptors
├── Catalytic receptors
│   ├── Receptor tyrosine kinases
│   ├── Receptor serine/threonine kinases
│   └── Cytokine receptors
├── Enzymes
│   ├── Kinases
│   ├── Phosphatases
│   ├── Proteases
│   └── Oxidoreductases
└── Transporters
    ├── SLC transporters
    └── ABC transporters
```

### 3.3 Seven Ligand Types

| Type | Description | Count |
|------|-------------|-------|
| `Synthetic organic` | Synthetic small molecules | ~8,000 |
| `Metabolite` | Endogenous metabolites | ~1,500 |
| `Natural product` | Plant/microbial compounds | ~800 |
| `Endogenous peptide` | Hormones, neuropeptides | ~600 |
| `Antibody` | Therapeutic antibodies | ~400 |
| `Inorganic` | Inorganic compounds | ~100 |
| `Labelled` | Radiolabeled/fluorescent | ~200 |

### 3.4 Six Interaction Affinity Types

GtoPdb uses logarithmic (pX) notation for affinity values:

| Type | Full Name | Description | Conversion |
|------|-----------|-------------|------------|
| `pKi` | Negative log of Ki | Inhibitor binding constant | pKi = -log10(Ki in M) |
| `pIC50` | Negative log of IC50 | 50% inhibition concentration | pIC50 = -log10(IC50 in M) |
| `pEC50` | Negative log of EC50 | 50% effective concentration | pEC50 = -log10(EC50 in M) |
| `pKd` | Negative log of Kd | Dissociation constant | pKd = -log10(Kd in M) |
| `pKB` | Negative log of KB | Antagonist equilibrium constant | pKB = -log10(KB in M) |
| `pA2` | Negative log of A2 | Schild analysis antagonist affinity | pA2 = -log10(A2 in M) |

#### Affinity Unit Conversions

| pX Value | nM | uM | pM |
|----------|----|----|-----|
| 9.0 | 1 | 0.001 | 1,000 |
| 8.0 | 10 | 0.01 | 10,000 |
| 7.0 | 100 | 0.1 | 100,000 |
| 6.0 | 1,000 | 1 | 1,000,000 |
| 5.0 | 10,000 | 10 | - |
| 4.0 | 100,000 | 100 | - |

**Conversion Formula:**
```
pX = -log10(X in M) = -log10(X in nM * 10^-9) = 9 - log10(X in nM)

Example: 10 nM = pX 8.0
         pX = 9 - log10(10) = 9 - 1 = 8
```

### 3.5 REST API Endpoints

**Base URL:** `https://www.guidetopharmacology.org/services/`
**Response Format:** JSON (default)

#### Target Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/targets` | GET | List all targets (filterable) |
| `/targets/{targetId}` | GET | Single target details |
| `/targets/families` | GET | List target families |
| `/targets/families/{familyId}` | GET | Single family details |
| `/targets/{targetId}/interactions` | GET | Interactions for target |
| `/targets/{targetId}/subunits` | GET | Component subunits |
| `/targets/{targetId}/geneProteinInformation` | GET | Gene/protein data |
| `/targets/{targetId}/databaseLinks` | GET | External database links |
| `/targets/{targetId}/diseases` | GET | Associated diseases |
| `/targets/{targetId}/variants` | GET | Genetic variants |
| `/targets/{targetId}/pdbStructure` | GET | PDB structures |

**Target Filter Parameters:**

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `type` | string | Target type filter | `GPCR` |
| `name` | string | Name search | `dopamine` |
| `geneSymbol` | string | Gene symbol | `DRD2` |
| `ecNumber` | string | EC number | `2.7.11.1` |
| `accession` | string | UniProt/RefSeq | `P14416` |

#### Ligand Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/ligands` | GET | List all ligands (filterable) |
| `/ligands/{ligandId}` | GET | Single ligand details |
| `/ligands/exact` | GET | Exact SMILES match |
| `/ligands/substructure` | GET | Substructure search |
| `/ligands/similarity` | GET | Similarity search |
| `/ligands/{ligandId}/interactions` | GET | Interactions for ligand |
| `/ligands/{ligandId}/molecularProperties` | GET | Physico-chemical properties |
| `/ligands/{ligandId}/synonyms` | GET | Name synonyms |
| `/ligands/{ligandId}/databaseLinks` | GET | External database links |
| `/ligands/{ligandId}/image` | GET | Structure image |

**Ligand Filter Parameters:**

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `type` | string | Ligand type | `Synthetic organic` |
| `name` | string | Name search | `aspirin` |
| `approved` | boolean | Approved drugs only | `true` |
| `immuno` | boolean | Immunopharmacology portal | `true` |
| `malaria` | boolean | Malaria portal | `true` |
| `antibacterial` | boolean | Antibacterial portal | `true` |

**Molecular Property Filters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `logP` | range | Partition coefficient |
| `molecularWeight` | range | Molecular weight (Da) |
| `hBondAcceptors` | integer | H-bond acceptors |
| `hBondDonors` | integer | H-bond donors |
| `rotatableBonds` | integer | Rotatable bonds |
| `polarSurfaceArea` | range | TPSA |
| `ruleOfFive` | integer | Lipinski violations |

#### Interaction Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/interactions` | GET | List all interactions (filterable) |
| `/interactions/{interactionId}` | GET | Single interaction details |

**Interaction Filter Parameters:**

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `targetId` | integer | Filter by target | `290` |
| `ligandId` | integer | Filter by ligand | `5085` |
| `type` | string | Interaction type | `Inhibitor` |
| `affinityType` | string | Affinity measurement | `pKi` |
| `species` | string | Species filter | `Human` |
| `primaryTarget` | boolean | Primary target only | `true` |

#### Disease Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/diseases` | GET | List diseases |
| `/diseases/{diseaseId}` | GET | Single disease |
| `/diseases/{diseaseId}/diseaseTargets` | GET | Associated targets |
| `/diseases/{diseaseId}/diseaseLigands` | GET | Associated ligands |

### 3.6 Response Data Models

#### Target Response

```json
{
  "targetId": 290,
  "name": "BRAF",
  "abbreviation": "BRAF",
  "systematicName": null,
  "type": "Enzyme",
  "familyIds": [839],
  "subunitIds": [],
  "complexIds": []
}
```

#### Target Fields Dictionary

| Field | Type | Description |
|-------|------|-------------|
| `targetId` | integer | GtoPdb target ID |
| `name` | string | Target name |
| `abbreviation` | string | Short name |
| `systematicName` | string | Systematic nomenclature |
| `type` | string | Target type |
| `familyIds` | array[int] | Parent family IDs |
| `subunitIds` | array[int] | Subunit target IDs |
| `complexIds` | array[int] | Complex target IDs |

#### Ligand Response

```json
{
  "ligandId": 5085,
  "name": "vemurafenib",
  "abbreviation": null,
  "inn": "vemurafenib",
  "type": "Synthetic organic",
  "species": null,
  "radioactive": false,
  "labelled": false,
  "approved": true,
  "withdrawn": false,
  "approvalSource": "FDA (2011)",
  "inchiKey": "DJPRUMZNCDNCBQ-UHFFFAOYSA-N",
  "smiles": "CCCS(=O)(=O)Nc1ccc(F)c(c1F)C(=O)c2cnc3ccc(Cl)cc3n2",
  "pubchemCid": "42611257",
  "chemblId": "CHEMBL1667"
}
```

#### Ligand Fields Dictionary

| Field | Type | Description |
|-------|------|-------------|
| `ligandId` | integer | GtoPdb ligand ID |
| `name` | string | Ligand name |
| `abbreviation` | string | Short name |
| `inn` | string | International nonproprietary name |
| `type` | string | Ligand type |
| `species` | string | Species-specific (if applicable) |
| `radioactive` | boolean | Radioactive ligand |
| `labelled` | boolean | Labelled compound |
| `approved` | boolean | Approved drug |
| `withdrawn` | boolean | Withdrawn drug |
| `approvalSource` | string | Approval details |
| `inchiKey` | string | InChI Key |
| `smiles` | string | SMILES structure |
| `pubchemCid` | string | PubChem CID |
| `chemblId` | string | ChEMBL ID |

#### Interaction Response

```json
{
  "interactionId": 12345,
  "targetId": 290,
  "ligandId": 5085,
  "type": "Inhibitor",
  "action": "Inhibition",
  "actionComment": null,
  "species": "Human",
  "endogenous": false,
  "selectivity": "Selective",
  "concentrationRange": "-",
  "affinity": "9.0",
  "affinityType": "pIC50",
  "originalAffinity": "1",
  "originalAffinityType": "IC50",
  "originalAffinityRelation": "=",
  "originalAffinityUnits": "nM",
  "assayDescription": "Inhibition of BRAF kinase activity",
  "assayConditions": null,
  "useDependent": false,
  "voltageDependent": false,
  "primaryTarget": true,
  "ligandContext": null,
  "references": [12345]
}
```

#### Interaction Fields Dictionary

| Field | Type | Description |
|-------|------|-------------|
| `interactionId` | integer | Interaction ID |
| `targetId` | integer | Target ID |
| `ligandId` | integer | Ligand ID |
| `type` | string | Interaction type |
| `action` | string | Action description |
| `actionComment` | string | Action notes |
| `species` | string | Species |
| `endogenous` | boolean | Endogenous ligand |
| `selectivity` | string | Selectivity class |
| `concentrationRange` | string | Concentration range |
| `affinity` | string | Affinity value (pX format) |
| `affinityType` | string | pKi, pIC50, pEC50, etc. |
| `originalAffinity` | string | Original value |
| `originalAffinityType` | string | Original unit type |
| `originalAffinityRelation` | string | =, <, > |
| `originalAffinityUnits` | string | Units (nM, uM) |
| `assayDescription` | string | Assay description |
| `assayConditions` | string | Assay conditions |
| `useDependent` | boolean | Use-dependent effect |
| `voltageDependent` | boolean | Voltage-dependent |
| `primaryTarget` | boolean | Primary target flag |
| `ligandContext` | string | Ligand context |
| `references` | array[int] | Reference IDs |

### 3.7 Download Formats

#### CSV Data Files

| File | Description | Records |
|------|-------------|---------|
| `targets_and_families.csv` | All targets with family | ~3,000 |
| `ligands.csv` | All ligands with properties | ~13,000 |
| `interactions.csv` | All interactions | ~170,000 |
| `approved_drugs.csv` | Approved drug ligands | ~2,200 |

#### Interaction Files by Target Type

| File | Target Type |
|------|-------------|
| `interactions_gpcr.csv` | GPCR interactions |
| `interactions_ic.csv` | Ion channel interactions |
| `interactions_cr.csv` | Catalytic receptor interactions |
| `interactions_enzyme.csv` | Enzyme interactions |
| `interactions_nhr.csv` | Nuclear hormone receptor |
| `interactions_transporter.csv` | Transporter interactions |

#### Sample CSV Data (interactions.csv)

```csv
target_id,target_name,ligand_id,ligand_name,type,action,species,affinity,affinity_type,original_affinity,original_units,pubmed_id
290,BRAF,5085,vemurafenib,Inhibitor,Inhibition,Human,9.0,pIC50,1,nM,20823850
100,Dopamine D2 receptor,384,haloperidol,Antagonist,Antagonism,Human,9.3,pKi,0.5,nM,2878367
```

#### Mapping Files

| File | Description |
|------|-------------|
| `hgnc_gene_ids.csv` | GtoPdb to HGNC gene ID |
| `uniprot_ids.csv` | GtoPdb to UniProt |

#### SQL Database

| File | Format | Description |
|------|--------|-------------|
| `gtopdb_pg.dmp` | PostgreSQL dump | Complete database |
| `gtopdb_schema.sql` | SQL | Schema only |

### 3.8 PostgreSQL Database Schema

**Version:** PostgreSQL 12.20

#### Core Tables

| Table | Description | Key Columns |
|-------|-------------|-------------|
| `target` | Target records | targetId, name, type |
| `target_family` | Family classification | familyId, name, type |
| `ligand` | Ligand records | ligandId, name, type |
| `interaction` | Target-ligand interactions | interactionId, targetId, ligandId |
| `reference` | Literature references | referenceId, pubmedId |
| `disease` | Disease annotations | diseaseId, name |
| `species` | Species lookup | speciesId, name |

#### Key Relationships

```
target ──┬── target_family (many-to-many)
         ├── interaction (one-to-many)
         ├── disease_target (many-to-many)
         └── database_link (one-to-many)

ligand ──┬── interaction (one-to-many)
         ├── disease_ligand (many-to-many)
         ├── ligand_synonym (one-to-many)
         └── database_link (one-to-many)

interaction ── reference (many-to-many)
```

### 3.9 ID Systems

| Entity | ID Format | Example |
|--------|-----------|---------|
| Target | Numeric targetId | `290` |
| Ligand | Numeric ligandId | `5085` |
| Family | Numeric familyId | `839` |
| Interaction | Numeric interactionId | `12345` |
| Reference | Numeric refId | `6789` |
| HGNC | HGNC ID | `HGNC:1097` |
| UniProt | Accession | `P15056` |
| PubChem | CID | `42611257` |
| ChEMBL | ChEMBL ID | `CHEMBL1667` |

### 3.10 Example API Calls

```bash
# Get all GPCR targets
curl "https://www.guidetopharmacology.org/services/targets?type=GPCR"

# Get interactions for a target with affinity filter
curl "https://www.guidetopharmacology.org/services/targets/290/interactions?affinityType=pKi&species=Human"

# Search ligands by approved status
curl "https://www.guidetopharmacology.org/services/ligands?approved=true"

# Get interactions by type
curl "https://www.guidetopharmacology.org/services/interactions?type=Inhibitor&species=Human"

# Similarity search
curl "https://www.guidetopharmacology.org/services/ligands/similarity?smiles=CC(C)Cc1ccc(cc1)C(C)C(=O)O&threshold=70"

# Get ligand with all cross-references
curl "https://www.guidetopharmacology.org/services/ligands/5085/databaseLinks"
```

---

## 4. Cross-Database Integration

### 4.1 Common Identifier Systems

| Database | Target/Gene ID | Compound/Drug ID | Protein ID |
|----------|----------------|------------------|------------|
| TTD | TTD Target ID | TTD Drug ID, PubChem, DrugBank, ChEMBL | UniProt |
| BindingDB | Polymer ID | Monomer ID, PubChem, ChEMBL, DrugBank | UniProt |
| GtoPdb | Target ID, HGNC | Ligand ID, PubChem, ChEMBL | UniProt |

### 4.2 Cross-Reference Mapping Strategy

**Recommended approach for integration:**

1. **Primary compound mapping:** Use InChIKey (27 characters) - exact structure match
2. **Secondary compound mapping:** Use ChEMBL ID or PubChem CID
3. **Target/protein mapping:** Use UniProt accession
4. **Gene mapping:** Use HGNC symbols

### 4.3 Affinity Unit Conversions

Converting between different affinity representations:

| From | To | Formula |
|------|-----|---------|
| nM | pX | pX = 9 - log10(nM) |
| uM | pX | pX = 6 - log10(uM) |
| pM | pX | pX = 12 - log10(pM) |
| pX | nM | nM = 10^(9-pX) |
| pX | uM | uM = 10^(6-pX) |
| Ki | pKi | pKi = -log10(Ki in M) |
| IC50 | pIC50 | pIC50 = -log10(IC50 in M) |

---

## 5. Integration Code Examples

### 5.1 Python: Query BindingDB REST API

```python
import requests
from typing import List, Dict, Optional

class BindingDBClient:
    """Client for BindingDB REST API."""

    BASE_URL = "https://bindingdb.org/rest"

    def get_ligands_by_uniprot(
        self,
        uniprot_id: str,
        cutoff_nm: int = 10000,
        response_format: str = "application/json"
    ) -> List[Dict]:
        """
        Retrieve binding data for a UniProt target.

        Args:
            uniprot_id: UniProt accession (e.g., P35354)
            cutoff_nm: Affinity cutoff in nM
            response_format: Response format

        Returns:
            List of binding affinity records
        """
        url = f"{self.BASE_URL}/getLigandsByUniprot"
        params = {
            "uniprot": f"{uniprot_id};{cutoff_nm}",
            "response": response_format
        }

        response = requests.get(url, params=params)
        response.raise_for_status()

        data = response.json()
        return data.get("affinities", [])

    def get_targets_by_smiles(
        self,
        smiles: str,
        similarity: float = 0.85
    ) -> List[Dict]:
        """
        Find targets for a compound by similarity.

        Args:
            smiles: SMILES structure
            similarity: Similarity threshold (0-1)

        Returns:
            List of matching targets
        """
        import urllib.parse

        url = f"{self.BASE_URL}/getTargetByCompound"
        params = {
            "smiles": urllib.parse.quote(smiles),
            "cutoff": similarity,
            "response": "application/json"
        }

        response = requests.get(url, params=params)
        response.raise_for_status()

        return response.json().get("targets", [])


# Example usage
client = BindingDBClient()

# Get all ligands for COX-2 with IC50 < 1 uM
cox2_ligands = client.get_ligands_by_uniprot("P35354", cutoff_nm=1000)

for ligand in cox2_ligands[:5]:
    print(f"Compound: {ligand.get('ligand_name', 'Unknown')}")
    print(f"  IC50: {ligand.get('ic50_value', 'N/A')} nM")
    print(f"  Ki: {ligand.get('ki_value', 'N/A')} nM")
    print(f"  ChEMBL: {ligand.get('chembl_id', 'N/A')}")
```

### 5.2 Python: Query GtoPdb REST API

```python
import requests
from typing import List, Dict, Optional
import math

class GtoPdbClient:
    """Client for GtoPdb REST API."""

    BASE_URL = "https://www.guidetopharmacology.org/services"

    def get_target(self, target_id: int) -> Dict:
        """Get target details by ID."""
        url = f"{self.BASE_URL}/targets/{target_id}"
        response = requests.get(url)
        response.raise_for_status()
        return response.json()

    def get_target_interactions(
        self,
        target_id: int,
        species: str = "Human",
        affinity_type: Optional[str] = None
    ) -> List[Dict]:
        """
        Get interactions for a target.

        Args:
            target_id: GtoPdb target ID
            species: Species filter (Human, Mouse, Rat)
            affinity_type: Filter by affinity type (pKi, pIC50, etc.)

        Returns:
            List of interaction records
        """
        url = f"{self.BASE_URL}/targets/{target_id}/interactions"
        params = {"species": species}
        if affinity_type:
            params["affinityType"] = affinity_type

        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json()

    def search_ligands(
        self,
        approved: bool = True,
        ligand_type: Optional[str] = None
    ) -> List[Dict]:
        """Search for ligands."""
        url = f"{self.BASE_URL}/ligands"
        params = {"approved": str(approved).lower()}
        if ligand_type:
            params["type"] = ligand_type

        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json()

    @staticmethod
    def convert_px_to_nm(px_value: float) -> float:
        """Convert pX value to nM."""
        return 10 ** (9 - px_value)

    @staticmethod
    def convert_nm_to_px(nm_value: float) -> float:
        """Convert nM value to pX."""
        return 9 - math.log10(nm_value)


# Example usage
client = GtoPdbClient()

# Get BRAF interactions
braf_interactions = client.get_target_interactions(
    target_id=290,
    species="Human",
    affinity_type="pIC50"
)

for interaction in braf_interactions[:5]:
    px = float(interaction.get("affinity", 0))
    nm = GtoPdbClient.convert_px_to_nm(px) if px else None

    print(f"Ligand ID: {interaction['ligandId']}")
    print(f"  Type: {interaction['type']}")
    print(f"  Affinity: {px} {interaction['affinityType']}")
    if nm:
        print(f"  Affinity: {nm:.2f} nM")
```

### 5.3 Python: Parse TTD Download Files

```python
from typing import Dict, List, Generator
from collections import defaultdict

def parse_ttd_target_file(filepath: str) -> Dict[str, Dict]:
    """
    Parse TTD target download file.

    Args:
        filepath: Path to P1-01-TTD_target_download.txt

    Returns:
        Dictionary of targets keyed by TTD_TARGET_ID
    """
    targets = defaultdict(dict)

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            ttd_id = parts[0]
            field_name = parts[2]
            content = parts[3]

            # Map TTD field names to standardized names
            field_mapping = {
                'TARGETNAME': 'name',
                'TARGETYPE': 'status',
                'BIOCLASS': 'biochemical_class',
                'UNIPROID': 'uniprot_id',
                'GENENAME': 'gene_name',
                'PDBIDS': 'pdb_ids',
                'ECNUMBER': 'ec_number',
                'FUNCTION': 'function',
                'KEGGPATH': 'kegg_pathway',
                'SEQUENCE': 'sequence'
            }

            if field_name in field_mapping:
                targets[ttd_id][field_mapping[field_name]] = content

    return dict(targets)


def parse_ttd_drug_file(filepath: str) -> Dict[str, Dict]:
    """
    Parse TTD drug download file.

    Args:
        filepath: Path to P1-02-TTD_drug_download.txt

    Returns:
        Dictionary of drugs keyed by TTD_DRUG_ID
    """
    drugs = defaultdict(dict)

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            ttd_id = parts[0]
            field_name = parts[2]
            content = parts[3]

            field_mapping = {
                'DRUGNAME': 'name',
                'TRADNAME': 'trade_name',
                'DRUGTYPE': 'drug_type',
                'DRUGSTAT': 'status',
                'DRUGSMIL': 'smiles',
                'DRUGINCH': 'inchi',
                'DRUGCOMP': 'formula',
                'PUBCHCID': 'pubchem_cid',
                'DRUGBKID': 'drugbank_id',
                'CHEMBDID': 'chembl_id',
                'CASNUMBE': 'cas_number'
            }

            if field_name in field_mapping:
                drugs[ttd_id][field_mapping[field_name]] = content

    return dict(drugs)


# Example usage
targets = parse_ttd_target_file('P1-01-TTD_target_download.txt')
drugs = parse_ttd_drug_file('P1-02-TTD_drug_download.txt')

# Find BRAF target
for ttd_id, target in targets.items():
    if target.get('gene_name') == 'BRAF':
        print(f"TTD ID: {ttd_id}")
        print(f"Name: {target.get('name')}")
        print(f"Status: {target.get('status')}")
        print(f"UniProt: {target.get('uniprot_id')}")
        print(f"PDB IDs: {target.get('pdb_ids')}")
        break
```

### 5.4 Python: Unified Affinity Data Aggregator

```python
import math
from typing import Dict, List, Optional
from dataclasses import dataclass
from enum import Enum

class AffinityType(Enum):
    IC50 = "IC50"
    KI = "Ki"
    KD = "Kd"
    EC50 = "EC50"
    PIC50 = "pIC50"
    PKI = "pKi"
    PKD = "pKd"
    PEC50 = "pEC50"
    PKB = "pKB"
    PA2 = "pA2"


@dataclass
class UnifiedAffinity:
    """Unified binding affinity record."""

    compound_id: str
    compound_name: str
    target_id: str
    target_name: str
    affinity_type: AffinityType
    affinity_value: float
    affinity_unit: str
    affinity_nm: Optional[float]  # Normalized to nM
    affinity_px: Optional[float]  # Normalized to pX
    source_database: str
    uniprot_id: Optional[str]
    chembl_id: Optional[str]
    pubchem_cid: Optional[str]
    inchikey: Optional[str]
    pubmed_id: Optional[str]

    @staticmethod
    def normalize_to_nm(value: float, unit: str) -> float:
        """Normalize affinity value to nanomolar."""
        unit_lower = unit.lower()

        if unit_lower == 'nm':
            return value
        elif unit_lower == 'um' or unit_lower == 'um':
            return value * 1000
        elif unit_lower == 'pm':
            return value / 1000
        elif unit_lower == 'm':
            return value * 1e9
        elif unit_lower.startswith('p'):
            # pX value: convert to nM
            return 10 ** (9 - value)
        else:
            return value

    @staticmethod
    def normalize_to_px(value: float, unit: str) -> float:
        """Normalize affinity value to pX (negative log10)."""
        unit_lower = unit.lower()

        if unit_lower.startswith('p'):
            return value
        else:
            # Convert to nM first, then to pX
            nm = UnifiedAffinity.normalize_to_nm(value, unit)
            if nm > 0:
                return 9 - math.log10(nm)
            return None


class AffinityAggregator:
    """Aggregate binding affinity data from multiple sources."""

    def __init__(self):
        self.records: List[UnifiedAffinity] = []

    def add_bindingdb_record(self, record: Dict) -> None:
        """Add a BindingDB record."""
        # Determine primary affinity type
        for aff_type in ['ic50_value', 'ki_value', 'kd_value', 'ec50_value']:
            value = record.get(aff_type)
            if value and value != 'null':
                type_name = aff_type.replace('_value', '').upper()

                unified = UnifiedAffinity(
                    compound_id=record.get('monomerid', ''),
                    compound_name=record.get('ligand_name', ''),
                    target_id=record.get('uniprot_id', ''),
                    target_name=record.get('target_name', ''),
                    affinity_type=AffinityType[type_name],
                    affinity_value=float(value),
                    affinity_unit='nM',
                    affinity_nm=float(value),
                    affinity_px=UnifiedAffinity.normalize_to_px(float(value), 'nM'),
                    source_database='BindingDB',
                    uniprot_id=record.get('uniprot_id'),
                    chembl_id=record.get('chembl_id'),
                    pubchem_cid=record.get('pubchem_cid'),
                    inchikey=record.get('inchikey'),
                    pubmed_id=record.get('pubmed_id')
                )
                self.records.append(unified)
                break

    def add_gtopdb_record(self, interaction: Dict, ligand: Dict) -> None:
        """Add a GtoPdb interaction record."""
        affinity = interaction.get('affinity')
        if not affinity:
            return

        affinity_type = interaction.get('affinityType', 'pIC50')

        unified = UnifiedAffinity(
            compound_id=str(interaction.get('ligandId', '')),
            compound_name=ligand.get('name', ''),
            target_id=str(interaction.get('targetId', '')),
            target_name='',  # Would need separate lookup
            affinity_type=AffinityType[affinity_type.upper().replace('P', 'P')],
            affinity_value=float(affinity),
            affinity_unit=affinity_type,
            affinity_nm=UnifiedAffinity.normalize_to_nm(float(affinity), affinity_type),
            affinity_px=float(affinity),
            source_database='GtoPdb',
            uniprot_id=None,  # Would need separate lookup
            chembl_id=ligand.get('chemblId'),
            pubchem_cid=ligand.get('pubchemCid'),
            inchikey=ligand.get('inchiKey'),
            pubmed_id=None
        )
        self.records.append(unified)

    def get_by_target(self, target_id: str) -> List[UnifiedAffinity]:
        """Get all records for a target."""
        return [r for r in self.records if r.target_id == target_id]

    def get_by_compound(self, compound_id: str) -> List[UnifiedAffinity]:
        """Get all records for a compound."""
        return [r for r in self.records if r.compound_id == compound_id]

    def get_by_inchikey(self, inchikey: str) -> List[UnifiedAffinity]:
        """Get all records by InChIKey (cross-database)."""
        return [r for r in self.records if r.inchikey == inchikey]


# Example usage
aggregator = AffinityAggregator()

# Add BindingDB records
bindingdb_record = {
    'monomerid': '50000001',
    'ligand_name': 'Ibuprofen',
    'target_name': 'Cyclooxygenase-2',
    'uniprot_id': 'P35354',
    'ic50_value': '320',
    'chembl_id': 'CHEMBL521',
    'pubchem_cid': '3672',
    'inchikey': 'HEFNNWSXXWATRW-JTQLQIEISA-N',
    'pubmed_id': '12345678'
}
aggregator.add_bindingdb_record(bindingdb_record)

# Query aggregated data
cox2_affinities = aggregator.get_by_target('P35354')
for aff in cox2_affinities:
    print(f"{aff.compound_name}: {aff.affinity_nm:.1f} nM ({aff.affinity_px:.2f} pX)")
```

---

## 6. Summary Comparison

| Feature | TTD | BindingDB | GtoPdb |
|---------|-----|-----------|--------|
| **Primary Focus** | Therapeutic targets | Binding affinities | Pharmacology |
| **API Type** | Web download | REST | REST |
| **Target Count** | 3,131 | ~9,400 | ~3,000 |
| **Compound Count** | 39,862 | 1,300,000+ | ~13,000 |
| **Interaction Count** | ~60,000 | 2,900,000+ | ~170,000 |
| **Affinity Data** | Limited | Yes (IC50, Ki, Kd, EC50) | Yes (pKi, pIC50, etc.) |
| **Affinity Units** | nM | nM | pX (negative log) |
| **Bulk Download** | TXT | TSV, SDF | CSV, PostgreSQL |
| **License** | Free academic | CC BY 3.0 | CC BY-SA 4.0 |
| **Cross-refs** | UniProt, PubChem, DrugBank, ChEMBL | UniProt, PubChem, ChEMBL, DrugBank | UniProt, PubChem, ChEMBL |
| **Update Frequency** | Annual | Weekly | Quarterly |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `IC50` | Half-maximal inhibitory concentration - concentration causing 50% inhibition | `100 nM` |
| `Ki` | Inhibition constant - binding affinity of inhibitor to target | `1.2 nM` |
| `Kd` | Dissociation constant - measure of binding affinity | `0.5 nM` |
| `EC50` | Half-maximal effective concentration - concentration for 50% effect | `50 nM` |
| `pKi` | Negative log10 of Ki in molar units | `9.0` (equals 1 nM) |
| `pIC50` | Negative log10 of IC50 in molar units | `8.0` (equals 10 nM) |
| `MonomerID` | BindingDB unique compound identifier | `50000001` |
| `TTD_TARGET_ID` | Therapeutic Target Database target identifier | `TTDT00001` |
| `TTD_DRUG_ID` | Therapeutic Target Database drug identifier | `D0A9YA` |
| `InChIKey` | 27-character hash of chemical structure for exact matching | `HEFNNWSXXWATRW-JTQLQIEISA-N` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Affinity | Strength of binding between drug and target | IC50, Ki, Kd |
| Druggability | Assessment of whether a target can be modulated by a drug | TTD Target Status |
| Successful Target | Target with at least one FDA-approved drug | TTD Target Status |
| Clinical Trial Target | Target in clinical development but not yet approved | TTD Target Status |
| pX Notation | Logarithmic scale where higher values indicate stronger binding | pKi, pIC50, pEC50 |
| Selectivity | Preference of a compound for one target over others | GtoPdb interactions |
| Primary Target | Main intended molecular target of a drug | GtoPdb primaryTarget |
| Ligand Type | Classification of chemical compound by origin or structure | Synthetic organic, Antibody |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| TTD | Therapeutic Target Database | Drug target database with 3,131 targets |
| GtoPdb | Guide to Pharmacology Database | IUPHAR/BPS pharmacology database |
| nM | Nanomolar | Concentration unit (10^-9 M) |
| uM | Micromolar | Concentration unit (10^-6 M) |
| pM | Picomolar | Concentration unit (10^-12 M) |
| GPCR | G Protein-Coupled Receptor | Major drug target family |
| NHR | Nuclear Hormone Receptor | Ligand-activated transcription factors |
| LGIC | Ligand-Gated Ion Channel | Ion channels activated by ligands |
| VGIC | Voltage-Gated Ion Channel | Ion channels activated by voltage |
| RTK | Receptor Tyrosine Kinase | Catalytic receptor family |
| SLC | Solute Carrier | Transporter superfamily |
| ABC | ATP-Binding Cassette | Transporter superfamily |
| INN | International Nonproprietary Name | Generic drug naming standard |
| MOA | Mechanism of Action | How drug achieves therapeutic effect |
| ATC | Anatomical Therapeutic Chemical | Drug classification system |

---

## 7. References

1. TTD 2024: Zhou Y, et al. Nucleic Acids Res. 2024;52(D1):D1465-D1477. https://doi.org/10.1093/nar/gkad751

2. BindingDB 2024: Liu T, et al. Nucleic Acids Res. 2025;53(D1):D1633-D1641. https://doi.org/10.1093/nar/gkae1039

3. GtoPdb 2024: Harding SD, et al. Nucleic Acids Res. 2024;52(D1):D1438-D1449. https://doi.org/10.1093/nar/gkad944

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation for TTD, BindingDB, GtoPdb |
