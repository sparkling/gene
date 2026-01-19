# Drug-Target Interaction Database Data Models

This document provides comprehensive technical documentation of the data models, schemas, and field descriptions for major drug-target interaction databases.

---

## Table of Contents

1. [DGIdb (Drug Gene Interaction Database)](#1-dgidb-drug-gene-interaction-database)
2. [Open Targets Platform](#2-open-targets-platform)
3. [TTD (Therapeutic Target Database)](#3-ttd-therapeutic-target-database)
4. [BindingDB](#4-bindingdb)
5. [GtoPdb (Guide to Pharmacology)](#5-gtopdb-guide-to-pharmacology)

---

## 1. DGIdb (Drug Gene Interaction Database)

**URL**: https://dgidb.org
**Current Version**: DGIdb 5.0
**Primary Reference**: [DGIdb 5.0: rebuilding the drug-gene interaction database (NAR 2024)](https://academic.oup.com/nar/article/52/D1/D1227/7416371)

### 1.1 Overview

DGIdb aggregates drug-gene interaction data from 44+ disparate sources to support hypothesis generation for clinicians and researchers. The database provides normalized drug and gene records with improved coverage (88.58% drug normalization, 94.91% gene normalization).

### 1.2 GraphQL API

**Endpoint**: `https://dgidb.org/api/graphql`
**Playground**: `https://dgidb.org/api/graphql` (interactive interface)

#### Core Entity Types

```graphql
# Primary entities in the DGIdb GraphQL schema
type Gene {
  name: String!
  conceptId: String       # Normalized gene identifier
  aliases: [String]
  geneCategories: [GeneCategory]
  interactions: [Interaction]
  attributes: [GeneAttribute]
}

type Drug {
  name: String!
  conceptId: String       # Normalized drug identifier
  aliases: [String]
  approvalStatus: String  # FDA approval status
  approvedIndications: [String]
  interactions: [Interaction]
  attributes: [DrugAttribute]
}

type Interaction {
  gene: Gene!
  drug: Drug!
  interactionTypes: [InteractionType]
  sources: [Source]
  publications: [Publication]
  score: Float            # Interaction score
  attributes: [InteractionAttribute]
}

type GeneCategory {
  name: String!
  sourceDbName: String
}

type Source {
  sourceDbName: String!
  sourceDbVersion: String
  citation: String
}
```

#### Example Queries

```graphql
# Query interactions for a specific gene
query {
  genes(names: ["BRAF"]) {
    nodes {
      name
      conceptId
      geneCategories {
        name
      }
      interactions {
        drug {
          name
          approvalStatus
        }
        interactionTypes {
          type
          directionality
        }
        sources {
          sourceDbName
        }
      }
    }
  }
}

# Query drugs by approval status
query {
  drugs(names: ["Vemurafenib"]) {
    nodes {
      name
      aliases
      approvalStatus
      interactions {
        gene {
          name
        }
        interactionTypes {
          type
        }
      }
    }
  }
}
```

### 1.3 Interaction Types

DGIdb defines 26 interaction types with associated directionality:

| Interaction Type | Directionality | Description |
|------------------|----------------|-------------|
| `activator` | Activating | Increases target activity |
| `agonist` | Activating | Binds and activates receptor |
| `antagonist` | Inhibiting | Blocks receptor activation |
| `inhibitor` | Inhibiting | Reduces target activity |
| `blocker` | Inhibiting | Prevents ion channel opening |
| `partial agonist` | Activating | Partial receptor activation |
| `inverse agonist` | Inhibiting | Reduces constitutive activity |
| `allosteric modulator` | Variable | Modifies target via allosteric site |
| `positive modulator` | Activating | Enhances target response |
| `negative modulator` | Inhibiting | Reduces target response |
| `inducer` | Activating | Increases target expression |
| `antibody` | Variable | Immunological binding |
| `antisense oligonucleotide` | Inhibiting | Gene expression knockdown |
| `binder` | N/A | Non-specific binding |
| `ligand` | Variable | Binds to target |
| `cofactor` | Activating | Required for target function |
| `chaperone` | Activating | Assists protein folding |
| `cleavage` | Inhibiting | Proteolytic processing |
| `adduct` | Variable | Covalent modification |
| `potentiator` | Activating | Enhances drug effect |
| `multitarget` | Variable | Multiple mechanism |
| `inhibitory allosteric modulator` | Inhibiting | Allosteric inhibition |
| `n/a` | Unknown | Not specified |
| `other/unknown` | Unknown | Unclassified |

### 1.4 Gene Druggability Categories

DGIdb assigns genes to druggability categories based on membership criteria:

| Category | Description | Sources |
|----------|-------------|---------|
| `CLINICALLY ACTIONABLE` | Genes informing clinical decisions | CIViC, OncoKB |
| `DRUGGABLE GENOME` | Union of druggable genome definitions | Hopkins & Groom, Russ & Lampel, dGene |
| `DRUG RESISTANCE` | Genes conferring resistance | GO:0042493, COSMIC |
| `KINASE` | Protein kinases | UniProt, Gene Ontology |
| `TUMOR SUPPRESSOR` | Tumor suppressor genes | CIViC, OncoKB |
| `TRANSCRIPTION FACTOR` | Transcriptional regulators | TF databases |
| `CELL SURFACE` | Surface-expressed proteins | Human Protein Atlas |
| `G PROTEIN COUPLED RECEPTOR` | GPCRs | Guide to Pharmacology |
| `ION CHANNEL` | Ion channels | Guide to Pharmacology |
| `NUCLEAR HORMONE RECEPTOR` | NHRs | Guide to Pharmacology |
| `TRANSPORTER` | Membrane transporters | SLC/ABC families |
| `PROTEASE` | Proteolytic enzymes | MEROPS |
| `PHOSPHATASE` | Protein phosphatases | DEPhOsphorylation |

### 1.5 Source Databases

DGIdb aggregates from 44+ sources (as of v5.0):

| Source Abbreviation | Full Name | Data Type |
|---------------------|-----------|-----------|
| CGI | Cancer Genome Interpreter | Interactions |
| ChEMBL | ChEMBL | Interactions |
| ChemIDplus | ChemIDplus | Drug claims |
| CIViC | Clinical Interpretation of Variants in Cancer | Interactions, Categories |
| COSMIC | Catalogue of Somatic Mutations in Cancer | Interactions, Categories |
| DTC | Drug Target Commons | Interactions |
| DrugBank | DrugBank | Interactions, Drug claims |
| Drugs@FDA | FDA Drug Database | Drug claims, Approval |
| GO | Gene Ontology | Gene categories |
| GToPdb | Guide to Pharmacology | Interactions |
| HGNC | HUGO Gene Nomenclature Committee | Gene claims |
| HemOnc | Hematology/Oncology | Drug claims |
| HPA | Human Protein Atlas | Gene categories |
| IDG | Illuminating the Druggable Genome | Gene categories |
| NCIt | NCI Thesaurus | Drug claims |
| OncoKB | Precision Oncology Knowledge Base | Interactions |
| PharmGKB | Pharmacogenomics Knowledge Base | Interactions |
| RxNorm | RxNorm | Drug claims |
| TALC | Targeted Agents in Lung Cancer | Interactions |
| Tempus | Tempus xT Panel | Gene categories |
| TTD | Therapeutic Target Database | Interactions |

### 1.6 Bulk Download TSV Format

**Download URL**: https://dgidb.org/downloads

#### Interactions TSV Columns

| Column | Description | Example |
|--------|-------------|---------|
| `gene_name` | HGNC gene symbol | `BRAF` |
| `gene_concept_id` | Normalized gene ID (NCBI) | `hgnc:1097` |
| `gene_aliases` | Alternative names (pipe-delimited) | `B-RAF1|BRAF1` |
| `drug_name` | Normalized drug name | `VEMURAFENIB` |
| `drug_concept_id` | Normalized drug ID | `rxcui:1147220` |
| `drug_aliases` | Alternative names | `PLX4032|RG7204` |
| `interaction_types` | Interaction types (pipe-delimited) | `inhibitor` |
| `interaction_score` | Computed interaction score | `12.5` |
| `source_db_names` | Contributing sources | `ChEMBL|DrugBank` |
| `pmids` | PubMed IDs | `20823850|21639808` |
| `approval_status` | FDA approval status | `Approved` |
| `approved_indications` | Approved uses | `Melanoma` |

#### Genes TSV Columns

| Column | Description |
|--------|-------------|
| `gene_name` | HGNC gene symbol |
| `gene_concept_id` | Normalized identifier |
| `gene_aliases` | Alternative names |
| `gene_categories` | Druggability categories |
| `interaction_count` | Number of interactions |

#### Drugs TSV Columns

| Column | Description |
|--------|-------------|
| `drug_name` | Normalized drug name |
| `drug_concept_id` | Normalized identifier |
| `drug_aliases` | Alternative names |
| `approval_status` | FDA approval status |
| `approved_indications` | Therapeutic indications |
| `active_andas_ndas` | Active ANDA/NDA applications |
| `interaction_count` | Number of interactions |

### 1.7 ID Systems

| Entity | ID Format | Example |
|--------|-----------|---------|
| Gene | HGNC ID | `hgnc:1097` |
| Drug (RxNorm) | RxCUI | `rxcui:1147220` |
| Drug (ChEMBL) | ChEMBL ID | `chembl:CHEMBL1336` |
| Drug (NCIt) | NCIt Code | `ncit:C64768` |
| Interaction | Internal UUID | Auto-generated |

---

## 2. Open Targets Platform

**URL**: https://platform.opentargets.org
**API Base**: https://api.platform.opentargets.org/api/v4/graphql
**Documentation**: https://platform-docs.opentargets.org
**Primary Reference**: [Open Targets Platform Documentation](https://platform-docs.opentargets.org/data-access/graphql-api)

### 2.1 Overview

The Open Targets Platform provides evidence-based target-disease associations integrating genetics, genomics, transcriptomics, drugs, and literature data to support drug target identification and prioritization.

### 2.2 GraphQL API Schema

**Endpoint**: `https://api.platform.opentargets.org/api/v4/graphql`
**Schema URL**: `https://api.platform.opentargets.org/api/v4/graphql/schema`
**Interactive Browser**: `https://api.platform.opentargets.org/api/v4/graphql/browser`

#### Core Entity Types

```graphql
type Target {
  id: String!                    # Ensembl Gene ID (e.g., ENSG00000157764)
  approvedSymbol: String!        # HGNC symbol
  approvedName: String!          # Full gene name
  biotype: String               # Gene biotype

  # Genomic location
  genomicLocation: GenomicLocation

  # Functional annotations
  functionDescriptions: [String]
  go: [GO]                      # Gene Ontology terms

  # Tractability assessments
  tractability: Tractability

  # Expression data
  expressions: [Expression]

  # Associated diseases
  associatedDiseases: AssociatedDiseases

  # Known drugs
  knownDrugs: KnownDrugs
}

type Disease {
  id: String!                   # EFO ID (e.g., EFO_0000311)
  name: String!                 # Disease name
  description: String

  # Ontology structure
  ancestors: [Disease]
  descendants: [Disease]
  parents: [Disease]
  children: [Disease]

  # Phenotypes
  phenotypes: [Phenotype]       # HPO annotations

  # Associated targets
  associatedTargets: AssociatedTargets

  # Known drugs
  knownDrugs: KnownDrugs
}

type Drug {
  id: String!                   # ChEMBL ID (e.g., CHEMBL1336)
  name: String!                 # Drug name

  # Clinical status
  maximumClinicalTrialPhase: Float
  isApproved: Boolean
  hasBeenWithdrawn: Boolean

  # Classification
  drugType: String              # Small molecule, Antibody, etc.

  # Mechanisms
  mechanismsOfAction: [MechanismOfAction]

  # Indications
  indications: [Indication]

  # Safety
  adverseEvents: [AdverseEvent]
  warnings: [DrugWarning]
}

type Evidence {
  id: String!
  targetId: String!
  diseaseId: String!
  score: Float!                 # Evidence score (0-1)

  # Source information
  datasourceId: String!
  datatypeId: String!

  # Evidence-specific fields (varies by source)
  variant: Variant
  drug: Drug
  study: Study
  literature: [String]          # PMIDs

  # Additional context
  clinicalSignificances: [String]
  confidence: String
}

type AssociatedDiseases {
  count: Int!
  rows: [AssociatedDisease]
}

type AssociatedDisease {
  disease: Disease!
  score: Float!                 # Overall association score

  # Score breakdown
  datatypeScores: [DatatypeScore]
  datasourceScores: [DatasourceScore]

  evidenceCount: Int!
}
```

### 2.3 Evidence Types and Data Sources

Evidence is categorized into data types, each containing multiple data sources:

#### Genetic Associations

| Data Source | Description | Scoring |
|-------------|-------------|---------|
| `ot_genetics_portal` | GWAS with L2G scores | L2G score > 0.05 |
| `gene_burden` | Rare variant burden | Scaled p-value (0.25-1.0) |
| `eva` | ClinVar variants | Clinical significance + review status |
| `eva_somatic` | ClinVar somatic | Same as eva |
| `gene2phenotype` | Gene-phenotype links | Confidence level |
| `genomics_england` | PanelApp genes | Panel confidence |
| `uniprot_variants` | UniProt variants | Fixed score |
| `uniprot_literature` | UniProt literature | Fixed score |
| `orphanet` | Rare diseases | Association type |
| `clingen` | ClinGen curations | Classification |

#### Somatic Mutations

| Data Source | Description |
|-------------|-------------|
| `cancer_gene_census` | COSMIC Cancer Gene Census |
| `intogen` | IntOGen driver genes |
| `cancer_biomarkers` | Cancer biomarker associations |

#### Drugs (Known Drug Target)

| Data Source | Description | Scoring |
|-------------|-------------|---------|
| `chembl` | ChEMBL clinical data | Clinical phase (0.05-1.0) |

#### Pathways & Systems Biology

| Data Source | Description |
|-------------|-------------|
| `crispr` | CRISPR screens |
| `slapenrich` | SLAPenrich |
| `progeny` | PROGENy pathway |
| `reactome` | Reactome pathways |
| `sysbio` | Systems biology |

#### RNA Expression

| Data Source | Description |
|-------------|-------------|
| `expression_atlas` | Differential expression |

#### Animal Models

| Data Source | Description |
|-------------|-------------|
| `impc` | Mouse phenotype data |

#### Literature

| Data Source | Description |
|-------------|-------------|
| `europepmc` | Europe PMC text mining |

### 2.4 Association Score Calculation

#### Score Formula

Association scores use a harmonic sum approach:

```
score = sum(evidence_score_i / i^2) / max_harmonic_sum

Where:
- evidence_score_i = score of i-th evidence (sorted descending)
- i = position in sorted list (1-indexed)
- max_harmonic_sum ≈ 1.644 (theoretical maximum)
```

#### Data Source Weights

| Data Source | Weight | Data Type |
|-------------|--------|-----------|
| Europe PMC | 0.2 | Literature |
| Expression Atlas | 0.2 | RNA Expression |
| IMPC | 0.2 | Animal Models |
| PROGENy | 0.5 | Pathways |
| SLAPenrich | 0.5 | Pathways |
| Cancer Biomarkers | 0.5 | Somatic |
| SysBio | 0.5 | Pathways |
| OTAR Projects | 0.5 | Various |
| All Others | 1.0 | Various |

#### Direct vs Indirect Associations

- **Direct**: Evidence directly links target to disease
- **Indirect**: Evidence applied through disease ontology hierarchy

### 2.5 Drug Annotation Format

```graphql
type Drug {
  id: String!                        # ChEMBL compound ID
  name: String!
  synonyms: [String]

  # Classification
  drugType: DrugType                 # SMALL_MOLECULE, ANTIBODY, PROTEIN, etc.

  # Clinical status
  maximumClinicalTrialPhase: Float   # 0, 0.5, 1, 2, 3, 4
  isApproved: Boolean
  hasBeenWithdrawn: Boolean
  withdrawnNotice: WithdrawnNotice

  # Chemistry
  tradeNames: [String]

  # Mechanisms of action
  mechanismsOfAction: [MechanismOfAction]

  # Linked diseases
  linkedDiseases: LinkedDiseases
  linkedTargets: LinkedTargets
}

type MechanismOfAction {
  mechanismOfAction: String          # e.g., "Kinase inhibitor"
  targetName: String
  targets: [Target]
  actionType: String                 # INHIBITOR, AGONIST, etc.
  references: [Reference]
}
```

### 2.6 Parquet File Schemas (Bulk Downloads)

**FTP Base**: `https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{version}/output/`

#### Targets Dataset Schema

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Ensembl gene ID |
| `approvedSymbol` | string | HGNC symbol |
| `approvedName` | string | Gene name |
| `biotype` | string | Gene biotype |
| `chromosome` | string | Chromosome |
| `start` | long | Start position |
| `end` | long | End position |
| `strand` | string | Strand (+/-) |
| `tpiClinical` | double | Clinical tractability |
| `tpiPreclinical` | double | Preclinical tractability |
| `functionDescriptions` | array[string] | Function descriptions |
| `synonyms` | array[string] | Alternative names |

#### Diseases Dataset Schema

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | EFO disease ID |
| `name` | string | Disease name |
| `description` | string | Description |
| `therapeuticAreas` | array[string] | Therapeutic area IDs |
| `synonyms` | array[string] | Alternative names |
| `ancestors` | array[string] | Parent term IDs |
| `descendants` | array[string] | Child term IDs |

#### Evidence Dataset Schema

| Field | Type | Description |
|-------|------|-------------|
| `targetId` | string | Ensembl gene ID |
| `diseaseId` | string | EFO disease ID |
| `datasourceId` | string | Data source name |
| `datatypeId` | string | Data type category |
| `score` | double | Evidence score (0-1) |
| `variantRsId` | string | rs ID (if applicable) |
| `studyId` | string | Study identifier |
| `literature` | array[string] | PubMed IDs |
| `clinicalSignificances` | array[string] | Clinical significance |
| `confidence` | string | Confidence level |

#### Associations Dataset Schema

| Field | Type | Description |
|-------|------|-------------|
| `targetId` | string | Ensembl gene ID |
| `diseaseId` | string | EFO disease ID |
| `score` | double | Overall association score |
| `evidenceCount` | long | Number of evidence items |
| `datatypeId` | string | Data type (if by datatype) |
| `datasourceId` | string | Data source (if by source) |

### 2.7 ID Systems

| Entity | ID Format | Example |
|--------|-----------|---------|
| Target/Gene | Ensembl Gene ID | `ENSG00000157764` |
| Disease | EFO ID | `EFO_0000311` |
| Drug | ChEMBL ID | `CHEMBL1336` |
| Variant | rsID or variant ID | `rs113488022` |
| Study | GCST ID | `GCST006979` |

### 2.8 Example Queries

```graphql
# Get target with associations
query TargetInfo($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    tractability {
      smallmolecule {
        topCategory
      }
      antibody {
        topCategory
      }
    }
    associatedDiseases(page: {index: 0, size: 10}) {
      count
      rows {
        disease {
          id
          name
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}

# Search for drugs
query DrugSearch($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    name
    drugType
    maximumClinicalTrialPhase
    mechanismsOfAction {
      mechanismOfAction
      actionType
      targets {
        id
        approvedSymbol
      }
    }
  }
}
```

---

## 3. TTD (Therapeutic Target Database)

**URL**: https://idrblab.net/ttd/ (or https://db.idrblab.net/ttd/)
**Schema URL**: https://db.idrblab.net/ttd/schema
**Downloads**: https://idrblab.net/ttd/full-data-download
**Primary Reference**: [TTD 2024 (NAR 2024)](https://academic.oup.com/nar/article/52/D1/D1465/7275004)

### 3.1 Overview

The Therapeutic Target Database (TTD) provides comprehensive druggability information for therapeutic targets, organized across three perspectives: molecular interactions, human system profiles, and cell-based expression variations.

### 3.2 Target Classification Schema

#### Target Status Categories

| Status | Definition | Count (2024) |
|--------|------------|--------------|
| `Successful` | Targeted by at least one FDA-approved drug | 426 |
| `Clinical Trial` | Not approved, but in clinical trials | 1,014 |
| `Preclinical/Patented` | Not in trials, but preclinical/patented | 212 |
| `Literature-reported` | Research targets with experimental drugs only | 1,479 |

#### Target Molecular Types

| Type | Description |
|------|-------------|
| `Protein` | Most common target type |
| `Nucleic Acid` | DNA, mRNA, miRNA, lncRNA targets |
| `Other Molecule` | Uric acid, iron, reactive oxygen species |

### 3.3 Drug Status Categories

| Status | Definition | Count (2024) |
|--------|------------|--------------|
| `Approved` | FDA/regulatory approved | 2,895 |
| `Clinical Trial` | In clinical development | 11,796 |
| `Preclinical/Patented` | Preclinical or patented | 5,041 |
| `Experimental` | Research compounds | 20,130 |

#### Drug Types

| Type | Description |
|------|-------------|
| `Small Molecule` | Traditional chemical drugs |
| `Antibody` | mAbs, ADCs, bispecifics, IgG mixtures |
| `Nucleic Acid Drug` | ASOs, siRNAs, saRNAs, miRNAs, mRNAs |
| `Cell Therapy` | CAR-T, stem cells |
| `Gene Therapy` | Gene delivery vectors |
| `Vaccine` | Preventive and therapeutic vaccines |

### 3.4 Data Model Structure

#### Target Data Fields

| Field | Description | Example |
|-------|-------------|---------|
| `TTD_TARGET_ID` | Internal target identifier | `TTDT00001` |
| `Target_Name` | Standard target name | `B-Raf proto-oncogene` |
| `Target_Type` | Molecular type | `Protein` |
| `Target_Status` | Development status | `Successful Target` |
| `UniProt_ID` | UniProt accession | `P15056` |
| `Gene_Name` | Gene symbol | `BRAF` |
| `PDB_ID` | PDB structure IDs | `4MNE;4MNF;5HI2` |
| `KEGG_Pathway` | KEGG pathway mappings | `hsa04010:MAPK` |
| `Biochemical_Class` | Protein family | `Kinase` |
| `Sequence` | Amino acid sequence | FASTA format |

#### Drug Data Fields

| Field | Description | Example |
|-------|-------------|---------|
| `TTD_DRUG_ID` | Internal drug identifier | `D0A9YA` |
| `Drug_Name` | Standard drug name | `Vemurafenib` |
| `Drug_Status` | Development status | `Approved` |
| `Drug_Type` | Drug modality | `Small Molecule` |
| `Therapeutic_Class` | ATC classification | `L01EC01` |
| `CAS_Number` | CAS registry number | `918504-65-1` |
| `PubChem_CID` | PubChem compound ID | `42611257` |
| `DrugBank_ID` | DrugBank identifier | `DB08881` |
| `ChEMBL_ID` | ChEMBL identifier | `CHEMBL1667` |
| `SMILES` | Chemical structure | `CCCS...` |
| `InChIKey` | InChIKey | `JOBWBQQXXXXXX-XXXX` |

#### Target-Drug Relationship Fields

| Field | Description |
|-------|-------------|
| `TTD_TARGET_ID` | Target identifier |
| `TTD_DRUG_ID` | Drug identifier |
| `Activity_Type` | Mechanism (inhibitor, agonist, etc.) |
| `MOA` | Mode of action |
| `Indication` | Disease indication |
| `ICD_Code` | ICD-11 disease code |
| `Clinical_Status` | Trial phase |
| `Reference` | PubMed IDs |

### 3.5 Druggability Perspectives

TTD organizes druggability data across three perspectives:

#### 1. Molecular Interactions/Regulations

| Data Category | Description |
|---------------|-------------|
| Ligand binding pocket | 3D structure of drug binding site |
| Protein-protein interactions | Network properties |
| Microbiota-drug regulation | Gut microbiome interactions |

#### 2. Human System Profiles

| Data Category | Description |
|---------------|-------------|
| Protein similarity | Similarity to non-family proteins |
| Pathway involvement | 241 life-essential pathways |
| Tissue distribution | 32 human tissues |

#### 3. Cell-Based Expression

| Data Category | Description |
|---------------|-------------|
| Cell type expression | 1,742 cell types |
| Exogenous stimuli | 625 environmental factors |
| Endogenous factors | 447 internal factors |

### 3.6 Download File Structure

TTD provides multiple download files:

| File | Content | Format |
|------|---------|--------|
| `P1-01-TTD_target_download.txt` | Target information | Tab-delimited |
| `P1-02-TTD_drug_download.txt` | Drug information | Tab-delimited |
| `P1-03-TTD_crossmatching.txt` | ID cross-references | Tab-delimited |
| `P1-04-Drug_disease.txt` | Drug-disease mappings | Tab-delimited |
| `P1-05-Target_disease.txt` | Target-disease mappings | Tab-delimited |
| `P1-06-Target_pathway.txt` | Pathway annotations | Tab-delimited |
| `P2-01-TTD_uniprot_all.txt` | UniProt mappings | Tab-delimited |
| `P2-02-TTD_pubchem_drug.txt` | PubChem mappings | Tab-delimited |

### 3.7 ID Systems

| Entity | ID Format | Example |
|--------|-----------|---------|
| Target | TTD Target ID | `TTDT00001` |
| Drug | TTD Drug ID | `D0A9YA` |
| Disease | ICD-11 Code | `2B90.0` |
| Cross-references | UniProt, PubChem, DrugBank, ChEMBL | Various |

---

## 4. BindingDB

**URL**: https://www.bindingdb.org
**API Documentation**: https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp
**Primary Reference**: [BindingDB 2024 (NAR 2025)](https://academic.oup.com/nar/article/53/D1/D1633/7906836)

### 4.1 Overview

BindingDB is a FAIR-compliant knowledgebase containing 2.9 million experimentally measured protein-small molecule binding affinities spanning 1.3 million compounds and thousands of protein targets.

### 4.2 REST API Endpoints

**Base URL**: `https://www.bindingdb.org/rest/` or `https://bindingdb.org/rest/`

#### getLigandsByPDBs

Retrieve binding data for PDB targets.

```
GET /getLigandsByPDBs?pdb={pdb_ids}&cutoff={affinity_nm}&identity={percent}&response={format}
```

| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| `pdb` | string | PDB ID(s), comma-separated | Yes |
| `cutoff` | integer | Affinity threshold (nM) | Yes |
| `identity` | integer | Sequence identity cutoff (%) | Yes |
| `response` | string | `application/json` or XML (default) | No |

**Example**:
```
https://bindingdb.org/rest/getLigandsByPDBs?pdb=1Q0L,3ANM&cutoff=100&identity=92&response=application/json
```

#### getLigandsByUniprots

Retrieve binding data for UniProt targets.

```
GET /getLigandsByUniprots?uniprot={uniprot_ids}&cutoff={affinity_nm}&code={filter}&response={format}
```

| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| `uniprot` | string | UniProt ID(s), comma-separated | Yes |
| `cutoff` | integer | Affinity threshold (nM) | Yes |
| `code` | integer | 0=all, 1=commercial, 2=FDA approved | No |
| `response` | string | Response format | No |

**Example**:
```
https://bindingdb.org/rest/getLigandsByUniprots?uniprot=P00176,P00183&cutoff=10000&response=application/json
```

#### getLigandsByUniprot (Single)

```
GET /getLigandsByUniprot?uniprot={uniprot_id};{cutoff}&response={format}
```

**Example**:
```
https://bindingdb.org/rest/getLigandsByUniprot?uniprot=P35355;100&response=application/json
```

#### getTargetByCompound

Find targets for a compound by structure similarity.

```
GET /getTargetByCompound?smiles={smiles}&cutoff={similarity}&response={format}
```

| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| `smiles` | string | SMILES structure (URL-encoded) | Yes |
| `cutoff` | decimal | Similarity threshold (0-1) | Yes |
| `response` | string | Response format | No |

### 4.3 Response Schema

#### JSON Response Structure

```json
{
  "affinities": [
    {
      "monomerid": "50000001",
      "smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
      "affinity_type": "IC50",
      "affinity_value": "100",
      "affinity_unit": "nM",
      "target_name": "Cyclooxygenase-2",
      "uniprot_id": "P35354",
      "organism": "Homo sapiens",
      "pubmed_id": "12345678",
      "patent": null
    }
  ]
}
```

#### Response Fields

| Field | Type | Description |
|-------|------|-------------|
| `monomerid` | string | BindingDB compound ID |
| `smiles` | string | SMILES structure |
| `inchi` | string | InChI string |
| `inchikey` | string | InChIKey |
| `affinity_type` | string | IC50, Ki, Kd, or EC50 |
| `affinity_value` | string | Numeric value |
| `affinity_unit` | string | nM (nanomolar) |
| `target_name` | string | Protein target name |
| `uniprot_id` | string | UniProt accession |
| `organism` | string | Species |
| `pubmed_id` | string | PubMed reference |
| `patent` | string | Patent number (if applicable) |

### 4.4 Binding Affinity Data Format

BindingDB stores four primary affinity measurement types:

| Type | Full Name | Count | Description |
|------|-----------|-------|-------------|
| `IC50` | Half-maximal inhibitory concentration | 1.8M | Concentration causing 50% inhibition |
| `Ki` | Inhibition constant | 560K | Inhibitor binding constant |
| `EC50` | Half-maximal effective concentration | 220K | Concentration for 50% effect |
| `Kd` | Dissociation constant | 100K | Binding affinity measure |

#### Additional Measurements

| Type | Description |
|------|-------------|
| `kon` | Association rate constant |
| `koff` | Dissociation rate constant |
| `ΔG` | Free energy of binding |
| `ΔH` | Enthalpy of binding |
| `-TΔS` | Entropy term |
| `pH` | Measurement pH |
| `Temperature` | Measurement temperature |

### 4.5 Download File Formats

**Download URL**: https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp

#### TSV Format

Main TSV files contain one row per binding measurement:

| Column | Description |
|--------|-------------|
| `BindingDB MonomerID` | Compound identifier |
| `Ligand SMILES` | Chemical structure |
| `Ligand InChI` | InChI notation |
| `Ligand InChI Key` | Hashed identifier |
| `BindingDB Ligand Name` | Compound name |
| `Target Name Assigned by Curator` | Curated target name |
| `Target Source Organism` | Species |
| `Ki (nM)` | Ki value |
| `IC50 (nM)` | IC50 value |
| `Kd (nM)` | Kd value |
| `EC50 (nM)` | EC50 value |
| `kon (M-1 s-1)` | Association rate |
| `koff (s-1)` | Dissociation rate |
| `pH` | pH of measurement |
| `Temp (C)` | Temperature |
| `Curation/DataSource` | Data source |
| `Article DOI` | Publication DOI |
| `PMID` | PubMed ID |
| `Patent Number` | Patent reference |
| `Authors` | Author list |
| `Institution` | Research institution |
| `BindingDB Entry DOI` | Entry DOI |
| `Link to Ligand in BindingDB` | Ligand URL |
| `Link to Target in BindingDB` | Target URL |
| `Link to PDB` | PDB structure link |
| `UniProt (SwissProt) Entry Name` | UniProt name |
| `UniProt (SwissProt) Primary ID` | UniProt accession |
| `UniProt (TrEMBL) Entry Name` | TrEMBL name |
| `UniProt (TrEMBL) Primary ID` | TrEMBL accession |
| `PubChem CID` | PubChem compound ID |
| `PubChem SID` | PubChem substance ID |
| `ChEBI ID` | ChEBI identifier |
| `ChEMBL ID` | ChEMBL compound ID |
| `DrugBank ID` | DrugBank identifier |
| `ZINC ID` | ZINC identifier |
| `Number of Protein Chains in Target` | Chain count |
| `BindingDB Target Chain Sequence` | Protein sequence |

#### SDF Format Variants

| Variant | Description |
|---------|-------------|
| `*_2D.sdf` | 2D coordinates |
| `*_3D.sdf` | 3D coordinates (Vconf computed) |
| `*_terse_*` | One compound per entry, multi-target in data blocks |

#### Supporting Files

| File | Description |
|------|-------------|
| `BindingDBTargetSequences.fasta` | Target protein sequences |
| `BindingDB_CID.txt` | BindingDB → PubChem CID mapping |
| `BindingDB_SID.txt` | BindingDB → PubChem SID mapping |
| `BindingDB_CHEBI_ID.txt` | BindingDB → ChEBI mapping |
| `BindingDB_DrugBankID.txt` | BindingDB → DrugBank mapping |
| `BindingDB_PubMed.txt` | PubMed ID collection |
| `BindingDB_UniProt.txt` | Polymer → UniProt mapping |
| `BDB_Assays.tsv` | Assay descriptions |
| `BDB_rsid_eaids.txt` | Reaction Set → Entry/Assay mapping |

### 4.6 ID Systems

| Entity | ID Format | Example |
|--------|-----------|---------|
| Compound (Monomer) | Numeric MonomerID | `50000001` |
| Target (Polymer) | Numeric PolymerID | `3000` |
| Entry | DOI | `10.37126/bdb...` |
| Cross-references | PubChem, ChEMBL, DrugBank, UniProt | Various |

---

## 5. GtoPdb (Guide to Pharmacology)

**URL**: https://www.guidetopharmacology.org
**API Base**: https://www.guidetopharmacology.org/services/
**Documentation**: https://www.guidetopharmacology.org/webServices.jsp
**Downloads**: https://www.guidetopharmacology.org/download.jsp
**Primary Reference**: [GtoPdb 2024 (NAR 2024)](https://academic.oup.com/nar/article/52/D1/D1438/7332061)

### 5.1 Overview

The IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb) is an expert-curated resource of ligand-activity-target relationships providing quantitative pharmacological data for drug targets.

### 5.2 REST API Endpoints

**Base URL**: `https://www.guidetopharmacology.org/services/`
**Response Format**: JSON (default)

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

**Target Filters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `type` | string | Target type filter |
| `name` | string | Name search |
| `geneSymbol` | string | Gene symbol |
| `ecNumber` | string | EC number |
| `accession` | string | UniProt/RefSeq |

**Target Types**:
- `GPCR` - G protein-coupled receptors
- `NHR` - Nuclear hormone receptors
- `LGIC` - Ligand-gated ion channels
- `VGIC` - Voltage-gated ion channels
- `OtherIC` - Other ion channels
- `Enzyme` - Enzymes
- `CatalyticReceptor` - Catalytic receptors
- `Transporter` - Transporters
- `OtherProtein` - Other proteins
- `AccessoryProtein` - Accessory proteins

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

**Ligand Filters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `type` | string | Ligand type |
| `name` | string | Name search |
| `approved` | boolean | Approved drugs only |
| `immuno` | boolean | Immunopharmacology |
| `malaria` | boolean | Malaria Portal |
| `antibacterial` | boolean | Antibacterial Portal |

**Ligand Types**:
- `Synthetic organic`
- `Metabolite`
- `Natural product`
- `Endogenous peptide`
- `Antibody`
- `Inorganic`
- `Labelled`

**Molecular Property Filters**:
- `logP` - Partition coefficient range
- `molecularWeight` - MW range
- `hBondAcceptors` - H-bond acceptor count
- `hBondDonors` - H-bond donor count
- `rotatableBonds` - Rotatable bond count
- `polarSurfaceArea` - TPSA range
- `ruleOfFive` - Lipinski violations

#### Interaction Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/interactions` | GET | List all interactions (filterable) |
| `/interactions/{interactionId}` | GET | Single interaction details |

**Interaction Filters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `targetId` | integer | Filter by target |
| `ligandId` | integer | Filter by ligand |
| `type` | string | Interaction type |
| `affinityType` | string | Affinity measurement type |
| `species` | string | Species (Human, Mouse, Rat) |
| `primaryTarget` | boolean | Primary target only |

**Interaction Types**:
- `Activator`
- `Agonist`
- `Allosteric modulator`
- `Antagonist`
- `Antibody`
- `Channel blocker`
- `Gating inhibitor`
- `Inhibitor`
- `Subunit-specific`

**Affinity Types**:
- `pA2` - Negative log of antagonist concentration
- `pEC50` - Negative log of EC50
- `pIC50` - Negative log of IC50
- `pKB` - Negative log of equilibrium dissociation constant
- `pKd` - Negative log of Kd
- `pKi` - Negative log of Ki

#### Disease Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/diseases` | GET | List diseases |
| `/diseases/{diseaseId}` | GET | Single disease |
| `/diseases/{diseaseId}/diseaseTargets` | GET | Associated targets |
| `/diseases/{diseaseId}/diseaseLigands` | GET | Associated ligands |

**Disease Databases**:
- OMIM
- Disease Ontology (DOID)
- Orphanet

#### Reference Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/refs` | GET | List references |
| `/refs/{referenceId}` | GET | Single reference |

### 5.3 Response Data Models

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

#### Target Fields

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
  "approvalSource": "FDA (2011)"
}
```

#### Ligand Fields

| Field | Type | Description |
|-------|------|-------------|
| `ligandId` | integer | GtoPdb ligand ID |
| `name` | string | Ligand name |
| `abbreviation` | string | Short name |
| `inn` | string | International nonproprietary name |
| `type` | string | Ligand type |
| `species` | string | Species-specific |
| `radioactive` | boolean | Radioactive ligand |
| `labelled` | boolean | Labelled compound |
| `approved` | boolean | Approved drug |
| `withdrawn` | boolean | Withdrawn drug |
| `approvalSource` | string | Approval details |

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

#### Interaction Fields

| Field | Type | Description |
|-------|------|-------------|
| `interactionId` | integer | Interaction ID |
| `targetId` | integer | Target ID |
| `ligandId` | integer | Ligand ID |
| `type` | string | Interaction type |
| `action` | string | Action description |
| `species` | string | Species |
| `endogenous` | boolean | Endogenous ligand |
| `selectivity` | string | Selectivity |
| `affinity` | string | Affinity value (pX format) |
| `affinityType` | string | pKi, pIC50, pEC50, etc. |
| `originalAffinity` | string | Original value |
| `originalAffinityType` | string | Original unit type |
| `originalAffinityUnits` | string | Units (nM, μM) |
| `assayDescription` | string | Assay description |
| `primaryTarget` | boolean | Primary target flag |
| `references` | array[int] | Reference IDs |

### 5.4 Target Family Classification

GtoPdb organizes targets into hierarchical families:

| Family Category | Subcategories |
|-----------------|---------------|
| **G protein-coupled receptors** | Class A (Rhodopsin), Class B (Secretin), Class C (Glutamate), Frizzled, Adhesion |
| **Ion channels** | Ligand-gated (Cys-loop, glutamate, P2X), Voltage-gated (sodium, potassium, calcium), Other |
| **Nuclear hormone receptors** | Thyroid, Retinoid, Steroid, Orphan |
| **Catalytic receptors** | Receptor tyrosine kinases, Receptor serine/threonine kinases, Cytokine receptors |
| **Enzymes** | Kinases, Proteases, Phosphatases, Oxidoreductases, Transferases |
| **Transporters** | SLC transporters, ABC transporters |

### 5.5 Download File Structure

**Download URL**: https://www.guidetopharmacology.org/download.jsp

#### CSV/TSV Data Files

| File | Description |
|------|-------------|
| `targets_and_families.csv` | All targets with family classification |
| `ligands.csv` | All ligands with properties |
| `interactions.csv` | Target-ligand interactions |
| `approved_drugs.csv` | Approved drug ligands |

#### Interaction Files by Target Type

| File | Target Type |
|------|-------------|
| `interactions_gpcr.csv` | GPCR interactions |
| `interactions_ic.csv` | Ion channel interactions |
| `interactions_cr.csv` | Catalytic receptor interactions |
| `interactions_enzyme.csv` | Enzyme interactions |
| `interactions_nhr.csv` | Nuclear hormone receptor interactions |
| `interactions_transporter.csv` | Transporter interactions |

#### Mapping Files

| File | Description |
|------|-------------|
| `hgnc_gene_ids.csv` | GtoPdb → HGNC gene ID mapping |
| `uniprot_ids.csv` | GtoPdb → UniProt mapping |

#### SDF Files

| File | Description |
|------|-------------|
| `ligand_structures.sdf` | All ligand structures with SMILES |

#### RDF/Linked Data

| File | Description |
|------|-------------|
| `gtopdb.n3` | Full database in Notation3 format |

### 5.6 PostgreSQL Database Schema

**Version**: PostgreSQL 12.20

#### Core Tables

| Table | Description |
|-------|-------------|
| `target` | Target records |
| `target_family` | Family classification |
| `ligand` | Ligand records |
| `interaction` | Target-ligand interactions |
| `reference` | Literature references |
| `disease` | Disease annotations |
| `species` | Species lookup |

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

### 5.7 ID Systems

| Entity | ID Format | Example |
|--------|-----------|---------|
| Target | Numeric targetId | `290` |
| Ligand | Numeric ligandId | `5085` |
| Family | Numeric familyId | `839` |
| Interaction | Numeric interactionId | `12345` |
| Reference | Numeric refId | `6789` |
| Cross-references | HGNC, UniProt, PubChem, ChEMBL | Various |

### 5.8 Example API Calls

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
```

---

## Cross-Database ID Mapping

### Common Identifier Systems

| Database | Gene ID | Drug/Compound ID | Target Protein ID |
|----------|---------|------------------|-------------------|
| DGIdb | HGNC | RxNorm, ChEMBL, NCIt | N/A |
| Open Targets | Ensembl | ChEMBL | Ensembl (gene-level) |
| TTD | Internal | Internal, PubChem, DrugBank | UniProt |
| BindingDB | N/A | Internal, PubChem, ChEMBL | UniProt |
| GtoPdb | HGNC | Internal, PubChem, ChEMBL | UniProt |

### Recommended Cross-Reference Strategy

1. **Gene/Target mapping**: Use HGNC symbols or Ensembl IDs as primary keys
2. **Drug/Compound mapping**: Use ChEMBL IDs or PubChem CIDs for cross-referencing
3. **Protein mapping**: Use UniProt accessions for protein-level data

---

## Summary Comparison

| Feature | DGIdb | Open Targets | TTD | BindingDB | GtoPdb |
|---------|-------|--------------|-----|-----------|--------|
| **API Type** | GraphQL | GraphQL | Web | REST | REST |
| **Focus** | Drug-gene interactions | Target-disease associations | Therapeutic targets | Binding affinities | Pharmacology |
| **Gene Coverage** | ~45,000 | ~60,000 | ~3,700 | Thousands | ~3,000 |
| **Drug Coverage** | ~10,000 | ~13,000 | ~40,000 | 1.3M compounds | ~13,000 |
| **Interaction Count** | ~100,000 | Millions of evidence | ~60,000 | 2.9M measurements | ~170,000 |
| **Affinity Data** | No | Indirect | Limited | Yes (IC50, Ki, Kd, EC50) | Yes (pKi, pIC50, etc.) |
| **Bulk Download** | TSV | Parquet | TXT | TSV, SDF | CSV, PostgreSQL |
| **License** | Open | Open | Open | Open | ODbL |

---

## References

1. DGIdb 5.0: Freshour SL, et al. Nucleic Acids Res. 2024;52(D1):D1227-D1235. https://doi.org/10.1093/nar/gkad1040
2. Open Targets Platform: Ochoa D, et al. Nucleic Acids Res. 2021;49(D1):D1302-D1310. https://doi.org/10.1093/nar/gkaa1027
3. TTD 2024: Zhou Y, et al. Nucleic Acids Res. 2024;52(D1):D1465-D1477. https://doi.org/10.1093/nar/gkad751
4. BindingDB 2024: Liu T, et al. Nucleic Acids Res. 2025;53(D1):D1633-D1641. https://doi.org/10.1093/nar/gkae1039
5. GtoPdb 2024: Harding SD, et al. Nucleic Acids Res. 2024;52(D1):D1438-D1449. https://doi.org/10.1093/nar/gkad944
