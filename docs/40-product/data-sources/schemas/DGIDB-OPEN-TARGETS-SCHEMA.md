# DGIdb and Open Targets Platform Schema Reference

## Overview

This document provides comprehensive technical documentation for two major drug-target interaction platforms: **DGIdb 5.0** (Drug Gene Interaction Database) and the **Open Targets Platform**. These complementary resources support drug discovery and target identification through different but synergistic approaches.

| Aspect | DGIdb 5.0 | Open Targets Platform |
|--------|-----------|----------------------|
| **Primary Focus** | Drug-gene interactions | Target-disease associations |
| **API Type** | GraphQL | GraphQL |
| **Gene Coverage** | ~45,000 genes | ~60,000 targets |
| **Drug Coverage** | ~10,000 drugs | ~13,000 drugs |
| **Interactions** | ~100,000 | Millions of evidence items |
| **Sources Aggregated** | 44 databases | 12+ data types |
| **License** | MIT | Apache 2.0 |

---

## Part 1: DGIdb 5.0 (Drug Gene Interaction Database)

### 1.1 Database Overview

| Property | Value |
|----------|-------|
| **Name** | Drug Gene Interaction Database |
| **Version** | 5.0 |
| **URL** | https://dgidb.org |
| **API Endpoint** | https://dgidb.org/api/graphql |
| **Interactive Playground** | https://dgidb.org/api/graphql |
| **Bulk Downloads** | https://dgidb.org/downloads |
| **Primary Reference** | [NAR 2024](https://academic.oup.com/nar/article/52/D1/D1227/7416371) |
| **License** | MIT |
| **Normalization Rate** | 88.58% drugs, 94.91% genes |

DGIdb aggregates drug-gene interaction data from 44+ disparate sources to support hypothesis generation for clinicians and researchers, particularly in precision oncology.

### 1.2 GraphQL API Schema

#### Core Entity Types

```graphql
# =============================================================================
# PRIMARY ENTITY: Gene
# =============================================================================
type Gene {
  # Core identifiers
  name: String!                    # HGNC gene symbol (e.g., "BRAF")
  conceptId: String                # Normalized gene ID (e.g., "hgnc:1097")

  # Alternative names
  aliases: [String]                # Alternative gene names/symbols

  # Classification
  geneCategories: [GeneCategory]   # Druggability categories

  # Relationships
  interactions: [Interaction]      # Drug interactions
  attributes: [GeneAttribute]      # Additional annotations

  # Gene claims from sources
  geneClaims: [GeneClaim]          # Source-specific gene records
}

# =============================================================================
# PRIMARY ENTITY: Drug
# =============================================================================
type Drug {
  # Core identifiers
  name: String!                    # Normalized drug name (e.g., "VEMURAFENIB")
  conceptId: String                # Normalized drug ID (e.g., "rxcui:1147220")

  # Alternative names
  aliases: [String]                # Trade names, synonyms

  # Regulatory status
  approvalStatus: String           # FDA approval status
  approvedIndications: [String]    # Approved therapeutic uses

  # Relationships
  interactions: [Interaction]      # Gene interactions
  attributes: [DrugAttribute]      # Additional annotations

  # Drug claims from sources
  drugClaims: [DrugClaim]          # Source-specific drug records
}

# =============================================================================
# PRIMARY ENTITY: Interaction
# =============================================================================
type Interaction {
  # Core relationships
  gene: Gene!                      # Target gene
  drug: Drug!                      # Interacting drug

  # Interaction characterization
  interactionTypes: [InteractionType]  # Mechanism types
  score: Float                     # Computed interaction score

  # Provenance
  sources: [Source]                # Contributing databases
  publications: [Publication]      # PubMed references
  attributes: [InteractionAttribute]   # Additional annotations

  # Interaction claims from sources
  interactionClaims: [InteractionClaim]
}

# =============================================================================
# SUPPORTING TYPES
# =============================================================================
type InteractionType {
  type: String!                    # Interaction type name
  directionality: String           # "Activating", "Inhibiting", "Variable", "Unknown"
  definition: String               # Type definition
}

type GeneCategory {
  name: String!                    # Category name (e.g., "KINASE")
  sourceDbName: String             # Source of category assignment
}

type Source {
  sourceDbName: String!            # Database name (e.g., "ChEMBL")
  sourceDbVersion: String          # Version identifier
  citation: String                 # Publication citation
  license: String                  # Data license
  licenseLink: String              # License URL
}

type Publication {
  pmid: String!                    # PubMed ID
  citation: String                 # Full citation
}

type GeneAttribute {
  name: String!                    # Attribute name
  value: String!                   # Attribute value
  sources: [Source]                # Source of attribute
}

type DrugAttribute {
  name: String!                    # Attribute name
  value: String!                   # Attribute value
  sources: [Source]                # Source of attribute
}

# =============================================================================
# CLAIMS (Source-specific records before normalization)
# =============================================================================
type GeneClaim {
  name: String!                    # Source-specific gene name
  sourceDbName: String!            # Source database
  aliases: [String]                # Source-specific aliases
}

type DrugClaim {
  name: String!                    # Source-specific drug name
  sourceDbName: String!            # Source database
  aliases: [String]                # Source-specific aliases
  primaryName: String              # Primary name from source
}

type InteractionClaim {
  gene: Gene
  drug: Drug
  interactionType: String
  sourceDbName: String!
  publications: [Publication]
}
```

#### Query Root Types

```graphql
type Query {
  # Gene queries
  genes(
    names: [String!]!              # Gene names to search
    sourceDbName: String           # Filter by source
  ): GeneConnection!

  gene(
    name: String!                  # Single gene lookup
  ): Gene

  # Drug queries
  drugs(
    names: [String!]!              # Drug names to search
    sourceDbName: String           # Filter by source
  ): DrugConnection!

  drug(
    name: String!                  # Single drug lookup
  ): Drug

  # Interaction queries
  interactions(
    genes: [String!]               # Filter by gene names
    drugs: [String!]               # Filter by drug names
    interactionTypes: [String!]    # Filter by interaction types
    sources: [String!]             # Filter by source databases
  ): InteractionConnection!

  # Category queries
  geneCategories: [GeneCategory!]!

  # Source queries
  sources: [Source!]!

  # Interaction type queries
  interactionTypes: [InteractionType!]!
}

type GeneConnection {
  nodes: [Gene!]!
  totalCount: Int!
  pageInfo: PageInfo!
}

type DrugConnection {
  nodes: [Drug!]!
  totalCount: Int!
  pageInfo: PageInfo!
}

type InteractionConnection {
  nodes: [Interaction!]!
  totalCount: Int!
  pageInfo: PageInfo!
}

type PageInfo {
  hasNextPage: Boolean!
  hasPreviousPage: Boolean!
  startCursor: String
  endCursor: String
}
```

### 1.3 Interaction Types (26 Types)

DGIdb defines 26 standardized interaction types with associated directionality:

| Interaction Type | Directionality | Description | Example |
|------------------|----------------|-------------|---------|
| `activator` | Activating | Increases target activity | Enzyme activators |
| `agonist` | Activating | Binds and activates receptor | GPCR agonists |
| `antagonist` | Inhibiting | Blocks receptor activation | Beta-blockers |
| `inhibitor` | Inhibiting | Reduces target activity | Kinase inhibitors |
| `blocker` | Inhibiting | Prevents ion channel opening | Calcium channel blockers |
| `partial agonist` | Activating | Partial receptor activation | Buprenorphine |
| `inverse agonist` | Inhibiting | Reduces constitutive activity | Beta-carboline |
| `allosteric modulator` | Variable | Modifies via allosteric site | Benzodiazepines |
| `positive modulator` | Activating | Enhances target response | PAMs |
| `negative modulator` | Inhibiting | Reduces target response | NAMs |
| `inducer` | Activating | Increases target expression | CYP inducers |
| `antibody` | Variable | Immunological binding | Monoclonal antibodies |
| `antisense oligonucleotide` | Inhibiting | Gene expression knockdown | ASOs |
| `binder` | N/A | Non-specific binding | General ligands |
| `ligand` | Variable | Binds to target | Receptor ligands |
| `cofactor` | Activating | Required for function | Vitamins |
| `chaperone` | Activating | Assists protein folding | HSP90 clients |
| `cleavage` | Inhibiting | Proteolytic processing | Protease substrates |
| `adduct` | Variable | Covalent modification | Covalent inhibitors |
| `potentiator` | Activating | Enhances drug effect | Drug combinations |
| `multitarget` | Variable | Multiple mechanisms | Polypharmacology |
| `inhibitory allosteric modulator` | Inhibiting | Allosteric inhibition | Specific NAMs |
| `suppressor` | Inhibiting | Suppresses target function | Gene suppressors |
| `modulator` | Variable | General modulation | Unspecified mechanism |
| `n/a` | Unknown | Not specified | Unknown mechanism |
| `other/unknown` | Unknown | Unclassified | Novel mechanisms |

### 1.4 Gene Druggability Categories (13 Categories)

| Category | Description | Primary Sources |
|----------|-------------|-----------------|
| `CLINICALLY ACTIONABLE` | Genes informing clinical decisions | CIViC, OncoKB |
| `DRUGGABLE GENOME` | Union of druggable genome definitions | Hopkins & Groom, Russ & Lampel, dGene |
| `DRUG RESISTANCE` | Genes conferring drug resistance | GO:0042493, COSMIC |
| `KINASE` | Protein kinases | UniProt, Gene Ontology |
| `TUMOR SUPPRESSOR` | Tumor suppressor genes | CIViC, OncoKB |
| `TRANSCRIPTION FACTOR` | Transcriptional regulators | TF databases |
| `CELL SURFACE` | Surface-expressed proteins | Human Protein Atlas |
| `G PROTEIN COUPLED RECEPTOR` | GPCRs | Guide to Pharmacology |
| `ION CHANNEL` | Ion channels | Guide to Pharmacology |
| `NUCLEAR HORMONE RECEPTOR` | Nuclear hormone receptors | Guide to Pharmacology |
| `TRANSPORTER` | Membrane transporters | SLC/ABC families |
| `PROTEASE` | Proteolytic enzymes | MEROPS |
| `PHOSPHATASE` | Protein phosphatases | DEPhOsphorylation |

### 1.5 Source Databases (44 Aggregated Sources)

| Source | Full Name | Data Type | URL |
|--------|-----------|-----------|-----|
| CGI | Cancer Genome Interpreter | Interactions | https://www.cancergenomeinterpreter.org |
| ChEMBL | ChEMBL | Interactions | https://www.ebi.ac.uk/chembl |
| ChemIDplus | ChemIDplus | Drug claims | https://chem.nlm.nih.gov |
| CIViC | Clinical Interpretation of Variants in Cancer | Interactions, Categories | https://civicdb.org |
| COSMIC | Catalogue of Somatic Mutations in Cancer | Interactions, Categories | https://cancer.sanger.ac.uk/cosmic |
| DTC | Drug Target Commons | Interactions | https://drugtargetcommons.fimm.fi |
| DrugBank | DrugBank | Interactions, Drug claims | https://go.drugbank.com |
| Drugs@FDA | FDA Drug Database | Drug claims, Approval | https://www.accessdata.fda.gov/scripts/cder/daf |
| GO | Gene Ontology | Gene categories | http://geneontology.org |
| GToPdb | Guide to Pharmacology | Interactions | https://www.guidetopharmacology.org |
| HGNC | HUGO Gene Nomenclature Committee | Gene claims | https://www.genenames.org |
| HemOnc | Hematology/Oncology | Drug claims | https://hemonc.org |
| HPA | Human Protein Atlas | Gene categories | https://www.proteinatlas.org |
| IDG | Illuminating the Druggable Genome | Gene categories | https://druggablegenome.net |
| NCIt | NCI Thesaurus | Drug claims | https://ncit.nci.nih.gov |
| OncoKB | Precision Oncology Knowledge Base | Interactions | https://www.oncokb.org |
| PharmGKB | Pharmacogenomics Knowledge Base | Interactions | https://www.pharmgkb.org |
| RxNorm | RxNorm | Drug claims | https://www.nlm.nih.gov/research/umls/rxnorm |
| TALC | Targeted Agents in Lung Cancer | Interactions | Clinical trial data |
| Tempus | Tempus xT Panel | Gene categories | https://www.tempus.com |
| TTD | Therapeutic Target Database | Interactions | https://db.idrblab.net/ttd |

### 1.6 Field Dictionaries

#### Gene Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `name` | String | HGNC gene symbol | `"BRAF"` |
| `conceptId` | String | Normalized HGNC ID | `"hgnc:1097"` |
| `aliases` | Array[String] | Alternative names | `["B-RAF1", "BRAF1", "RAFB1"]` |
| `geneCategories` | Array[GeneCategory] | Druggability categories | `[{name: "KINASE"}]` |
| `interactions` | Array[Interaction] | Drug interactions | See Interaction type |
| `attributes` | Array[GeneAttribute] | Additional annotations | Various |

#### Drug Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `name` | String | Normalized drug name | `"VEMURAFENIB"` |
| `conceptId` | String | Normalized ID (RxNorm preferred) | `"rxcui:1147220"` |
| `aliases` | Array[String] | Trade names, synonyms | `["Zelboraf", "PLX4032"]` |
| `approvalStatus` | String | FDA approval status | `"Approved"` |
| `approvedIndications` | Array[String] | Approved uses | `["Melanoma with BRAF V600E"]` |

#### Interaction Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `gene` | Gene | Target gene | See Gene type |
| `drug` | Drug | Interacting drug | See Drug type |
| `interactionTypes` | Array[InteractionType] | Mechanism types | `[{type: "inhibitor"}]` |
| `score` | Float | Interaction score | `12.5` |
| `sources` | Array[Source] | Contributing databases | `[{sourceDbName: "ChEMBL"}]` |
| `publications` | Array[Publication] | PubMed references | `[{pmid: "20823850"}]` |

### 1.7 Identifier Systems

| Entity | ID System | Format | Example |
|--------|-----------|--------|---------|
| Gene | HGNC ID | `hgnc:{id}` | `hgnc:1097` |
| Drug (RxNorm) | RxCUI | `rxcui:{id}` | `rxcui:1147220` |
| Drug (ChEMBL) | ChEMBL ID | `chembl:{id}` | `chembl:CHEMBL1336` |
| Drug (NCIt) | NCIt Code | `ncit:{code}` | `ncit:C64768` |
| Interaction | Internal UUID | Auto-generated | `interaction-uuid` |
| Publication | PubMed ID | Numeric | `20823850` |

### 1.8 GraphQL Query Examples

#### Query Gene Interactions

```graphql
query GetGeneInteractions($geneNames: [String!]!) {
  genes(names: $geneNames) {
    nodes {
      name
      conceptId
      geneCategories {
        name
        sourceDbName
      }
      interactions {
        drug {
          name
          conceptId
          approvalStatus
          approvedIndications
        }
        interactionTypes {
          type
          directionality
        }
        score
        sources {
          sourceDbName
          sourceDbVersion
        }
        publications {
          pmid
        }
      }
    }
    totalCount
  }
}

# Variables
{
  "geneNames": ["BRAF", "EGFR", "ALK"]
}
```

#### Query Drug Information

```graphql
query GetDrugDetails($drugNames: [String!]!) {
  drugs(names: $drugNames) {
    nodes {
      name
      conceptId
      aliases
      approvalStatus
      approvedIndications
      interactions {
        gene {
          name
          conceptId
          geneCategories {
            name
          }
        }
        interactionTypes {
          type
          directionality
        }
        score
      }
      attributes {
        name
        value
      }
    }
    totalCount
  }
}

# Variables
{
  "drugNames": ["Vemurafenib", "Imatinib", "Gefitinib"]
}
```

#### Query by Interaction Type

```graphql
query GetInhibitorInteractions {
  interactions(
    interactionTypes: ["inhibitor"]
    sources: ["ChEMBL", "DrugBank"]
  ) {
    nodes {
      gene {
        name
        geneCategories {
          name
        }
      }
      drug {
        name
        approvalStatus
      }
      interactionTypes {
        type
        directionality
      }
      score
      sources {
        sourceDbName
      }
    }
    totalCount
  }
}
```

#### Query Gene Categories

```graphql
query GetDruggableGenes {
  genes(names: ["*"]) {
    nodes {
      name
      conceptId
      geneCategories {
        name
        sourceDbName
      }
    }
  }
}
```

### 1.9 Sample JSON Response

```json
{
  "data": {
    "genes": {
      "nodes": [
        {
          "name": "BRAF",
          "conceptId": "hgnc:1097",
          "geneCategories": [
            {
              "name": "KINASE",
              "sourceDbName": "Gene Ontology"
            },
            {
              "name": "CLINICALLY ACTIONABLE",
              "sourceDbName": "CIViC"
            },
            {
              "name": "DRUGGABLE GENOME",
              "sourceDbName": "dGene"
            }
          ],
          "interactions": [
            {
              "drug": {
                "name": "VEMURAFENIB",
                "conceptId": "rxcui:1147220",
                "approvalStatus": "Approved",
                "approvedIndications": [
                  "Melanoma with BRAF V600E mutation"
                ]
              },
              "interactionTypes": [
                {
                  "type": "inhibitor",
                  "directionality": "Inhibiting"
                }
              ],
              "score": 15.2,
              "sources": [
                {
                  "sourceDbName": "ChEMBL",
                  "sourceDbVersion": "32"
                },
                {
                  "sourceDbName": "DrugBank",
                  "sourceDbVersion": "5.1.9"
                },
                {
                  "sourceDbName": "CIViC",
                  "sourceDbVersion": "2023-10"
                }
              ],
              "publications": [
                {
                  "pmid": "20823850"
                },
                {
                  "pmid": "21639808"
                }
              ]
            },
            {
              "drug": {
                "name": "DABRAFENIB",
                "conceptId": "rxcui:1424911",
                "approvalStatus": "Approved",
                "approvedIndications": [
                  "Melanoma with BRAF V600E/K mutation",
                  "NSCLC with BRAF V600E mutation"
                ]
              },
              "interactionTypes": [
                {
                  "type": "inhibitor",
                  "directionality": "Inhibiting"
                }
              ],
              "score": 14.8,
              "sources": [
                {
                  "sourceDbName": "ChEMBL",
                  "sourceDbVersion": "32"
                },
                {
                  "sourceDbName": "OncoKB",
                  "sourceDbVersion": "2023-11"
                }
              ],
              "publications": [
                {
                  "pmid": "22608338"
                }
              ]
            }
          ]
        }
      ],
      "totalCount": 1
    }
  }
}
```

### 1.10 Bulk Download TSV Format

#### Interactions TSV Schema

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `gene_name` | String | HGNC gene symbol | `BRAF` |
| `gene_concept_id` | String | Normalized gene ID | `hgnc:1097` |
| `gene_aliases` | String | Pipe-delimited aliases | `B-RAF1|BRAF1` |
| `drug_name` | String | Normalized drug name | `VEMURAFENIB` |
| `drug_concept_id` | String | Normalized drug ID | `rxcui:1147220` |
| `drug_aliases` | String | Pipe-delimited aliases | `PLX4032|RG7204|Zelboraf` |
| `interaction_types` | String | Pipe-delimited types | `inhibitor` |
| `interaction_score` | Float | Computed score | `15.2` |
| `source_db_names` | String | Pipe-delimited sources | `ChEMBL|DrugBank|CIViC` |
| `pmids` | String | Pipe-delimited PubMed IDs | `20823850|21639808` |
| `approval_status` | String | FDA approval status | `Approved` |
| `approved_indications` | String | Pipe-delimited indications | `Melanoma with BRAF V600E` |

#### Genes TSV Schema

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `gene_name` | String | HGNC gene symbol | `BRAF` |
| `gene_concept_id` | String | Normalized ID | `hgnc:1097` |
| `gene_aliases` | String | Pipe-delimited aliases | `B-RAF1|BRAF1` |
| `gene_categories` | String | Pipe-delimited categories | `KINASE|CLINICALLY ACTIONABLE` |
| `interaction_count` | Integer | Number of interactions | `45` |

#### Drugs TSV Schema

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `drug_name` | String | Normalized drug name | `VEMURAFENIB` |
| `drug_concept_id` | String | Normalized ID | `rxcui:1147220` |
| `drug_aliases` | String | Pipe-delimited aliases | `PLX4032|Zelboraf` |
| `approval_status` | String | FDA status | `Approved` |
| `approved_indications` | String | Pipe-delimited indications | `Melanoma` |
| `active_andas_ndas` | String | Active applications | `NDA203505` |
| `interaction_count` | Integer | Number of interactions | `12` |

---

## Part 2: Open Targets Platform

### 2.1 Platform Overview

| Property | Value |
|----------|-------|
| **Name** | Open Targets Platform |
| **URL** | https://platform.opentargets.org |
| **API Endpoint** | https://api.platform.opentargets.org/api/v4/graphql |
| **API Browser** | https://api.platform.opentargets.org/api/v4/graphql/browser |
| **Documentation** | https://platform-docs.opentargets.org |
| **Bulk Downloads** | https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/ |
| **Primary Reference** | [NAR 2021](https://doi.org/10.1093/nar/gkaa1027) |
| **License** | Apache 2.0 |

The Open Targets Platform integrates genetics, genomics, transcriptomics, drugs, and literature data to provide evidence-based target-disease associations for drug target identification and prioritization.

### 2.2 GraphQL API Schema

#### Core Entity Types

```graphql
# =============================================================================
# PRIMARY ENTITY: Target
# =============================================================================
type Target {
  # Core identifiers
  id: String!                         # Ensembl Gene ID (e.g., "ENSG00000157764")
  approvedSymbol: String!             # HGNC symbol (e.g., "BRAF")
  approvedName: String!               # Full gene name

  # Gene classification
  biotype: String                     # Gene biotype (e.g., "protein_coding")

  # Genomic location
  genomicLocation: GenomicLocation

  # Functional annotations
  functionDescriptions: [String]      # UniProt function descriptions
  go: [GO]                            # Gene Ontology terms

  # Druggability assessment
  tractability: Tractability

  # Expression data
  expressions: [Expression]

  # Cross-references
  proteinIds: [IdAndSource]           # UniProt IDs
  dbXrefs: [IdAndSource]              # External database references

  # Associations
  associatedDiseases(
    page: Pagination
    orderByScore: String
    BFilter: String
    facetFilters: [String]
  ): AssociatedDiseases

  # Known drugs targeting this gene
  knownDrugs(
    page: Pagination
    freeTextQuery: String
  ): KnownDrugs

  # Pathways
  pathways: [Pathway]

  # Protein interactions
  interactions: TargetInteractions

  # Mouse phenotypes
  mousePhenotypes: [MousePhenotype]

  # Genetic constraint scores
  geneticConstraint: [GeneticConstraint]

  # Homologues
  homologues: [Homologue]

  # Safety information
  safetyLiabilities: [SafetyLiability]

  # Subcellular location
  subcellularLocations: [SubcellularLocation]

  # Target enabling packages
  tep: TEP
}

# =============================================================================
# PRIMARY ENTITY: Disease
# =============================================================================
type Disease {
  # Core identifiers
  id: String!                         # EFO ID (e.g., "EFO_0000311")
  name: String!                       # Disease name
  description: String                 # Disease description

  # Ontology structure
  ancestors: [String]                 # Ancestor EFO IDs
  descendants: [String]               # Descendant EFO IDs
  parents: [Disease]                  # Parent disease objects
  children: [Disease]                 # Child disease objects

  # Classification
  therapeuticAreas: [TherapeuticArea]

  # Synonyms
  synonyms: [Synonym]

  # Cross-references
  dbXRefs: [String]

  # Phenotypes
  phenotypes(page: Pagination): DiseasePhenotypes

  # Associations
  associatedTargets(
    page: Pagination
    orderByScore: String
    BFilter: String
    facetFilters: [String]
  ): AssociatedTargets

  # Known drugs for this disease
  knownDrugs(
    page: Pagination
    freeTextQuery: String
  ): KnownDrugs

  # Ontology information
  isTherapeuticArea: Boolean
  ontology: DiseaseOntology
}

# =============================================================================
# PRIMARY ENTITY: Drug
# =============================================================================
type Drug {
  # Core identifiers
  id: String!                         # ChEMBL ID (e.g., "CHEMBL1336")
  name: String!                       # Drug name

  # Synonyms and trade names
  synonyms: [String]
  tradeNames: [String]

  # Classification
  drugType: String                    # Small molecule, Antibody, etc.

  # Clinical status
  maximumClinicalTrialPhase: Float    # 0, 0.5, 1, 2, 3, 4
  isApproved: Boolean
  hasBeenWithdrawn: Boolean
  withdrawnNotice: WithdrawnNotice

  # Pharmacology
  mechanismsOfAction: [MechanismOfAction]

  # Linked entities
  linkedDiseases: LinkedDiseases
  linkedTargets: LinkedTargets

  # Indications
  indications: Indications

  # Safety
  adverseEvents(page: Pagination): AdverseEvents
  warnings: [DrugWarning]

  # Chemistry
  blackBoxWarning: Boolean
  description: String
  parentId: String
  childChemblIds: [String]

  # Cross-references
  crossReferences: [DrugCrossReference]
}

# =============================================================================
# PRIMARY ENTITY: Evidence
# =============================================================================
type Evidence {
  id: String!                         # Evidence record ID
  targetId: String!                   # Ensembl gene ID
  diseaseId: String!                  # EFO disease ID

  # Scoring
  score: Float!                       # Evidence score (0-1)

  # Source information
  datasourceId: String!               # Data source identifier
  datatypeId: String!                 # Data type category

  # Evidence-specific fields (varies by source)

  # Genetic evidence
  variantId: String
  variantRsId: String
  studyId: String
  beta: Float
  betaConfidenceIntervalLower: Float
  betaConfidenceIntervalUpper: Float
  oddsRatio: Float
  oddsRatioConfidenceIntervalLower: Float
  oddsRatioConfidenceIntervalUpper: Float
  pValueMantissa: Float
  pValueExponent: Int

  # Drug evidence
  drug: Drug
  clinicalPhase: Float
  clinicalStatus: String
  clinicalUrls: [Url]

  # Literature
  literature: [String]                # PubMed IDs

  # Pathways
  pathwayId: String
  pathwayName: String

  # Clinical significance
  clinicalSignificances: [String]
  confidence: String

  # Allele information
  allelicRequirements: [String]

  # Target information
  targetFromSourceId: String

  # Disease information
  diseaseFromSourceId: String
  diseaseFromSourceMappedId: String

  # Resource information
  resourceScore: Float
  urls: [Url]
}

# =============================================================================
# ASSOCIATION TYPES
# =============================================================================
type AssociatedDiseases {
  count: Int!
  rows: [AssociatedDisease]
}

type AssociatedDisease {
  disease: Disease!
  score: Float!                       # Overall association score

  # Score breakdown by data type
  datatypeScores: [DatatypeScore]

  # Score breakdown by data source
  datasourceScores: [DatasourceScore]

  # Evidence count
  evidenceCount: Int!
}

type AssociatedTargets {
  count: Int!
  rows: [AssociatedTarget]
}

type AssociatedTarget {
  target: Target!
  score: Float!
  datatypeScores: [DatatypeScore]
  datasourceScores: [DatasourceScore]
  evidenceCount: Int!
}

type DatatypeScore {
  id: String!                         # Data type ID
  score: Float!                       # Aggregated score
  evidenceCount: Int
}

type DatasourceScore {
  id: String!                         # Data source ID
  score: Float!                       # Aggregated score
  evidenceCount: Int
}
```

#### Supporting Types

```graphql
# =============================================================================
# SUPPORTING TYPES
# =============================================================================
type GenomicLocation {
  chromosome: String!
  start: Long!
  end: Long!
  strand: Int!
}

type GO {
  id: String!                         # GO term ID
  term: String!                       # GO term name
  category: String                    # BP, MF, CC
}

type Tractability {
  smallmolecule: [TractabilityBucket]
  antibody: [TractabilityBucket]
  protac: [TractabilityBucket]
  otherModalities: [TractabilityBucket]
}

type TractabilityBucket {
  label: String!
  value: Boolean!
}

type Expression {
  tissue: Tissue!
  rna: ExpressionData
  protein: ExpressionData
}

type Tissue {
  id: String!
  label: String!
  anatomicalSystems: [String]
  organs: [String]
}

type ExpressionData {
  value: Float
  unit: String
  level: Int                          # 0-3 expression level
}

type MechanismOfAction {
  mechanismOfAction: String           # e.g., "Kinase inhibitor"
  targetName: String
  targets: [Target]
  actionType: String                  # INHIBITOR, AGONIST, etc.
  references: [Reference]
}

type Indication {
  disease: Disease!
  maxPhaseForIndication: Float
  references: [Reference]
}

type Pathway {
  id: String!
  name: String!
  pathway: String
  topLevelTerm: String
}

type Reference {
  ids: [String]
  source: String
  urls: [String]
}

type Pagination {
  index: Int!
  size: Int!
}

type KnownDrugs {
  count: Int!
  cursor: String
  rows: [KnownDrug]
}

type KnownDrug {
  approvedSymbol: String
  approvedName: String
  drugId: String
  drugFromSource: String
  drug: Drug
  phase: Float
  status: String
  targetClass: [String]
  mechanismOfAction: String
  disease: Disease
  ctIds: [String]
  urls: [Url]
}

type Url {
  url: String!
  name: String
}

type LinkedDiseases {
  count: Int!
  rows: [Disease]
}

type LinkedTargets {
  count: Int!
  rows: [Target]
}

type TherapeuticArea {
  id: String!
  name: String!
}
```

### 2.3 Evidence Data Types (12 Categories)

Evidence is organized into 12 data types, each containing multiple data sources:

#### 1. Genetic Associations

| Data Source | ID | Description | Scoring Method |
|-------------|-----|-------------|----------------|
| Open Targets Genetics | `ot_genetics_portal` | GWAS with L2G machine learning | L2G score (>0.05 threshold) |
| Gene Burden | `gene_burden` | Rare variant burden tests | Scaled p-value (0.25-1.0) |
| ClinVar (germline) | `eva` | ClinVar germline variants | Clinical significance + review |
| ClinVar (somatic) | `eva_somatic` | ClinVar somatic variants | Same as eva |
| Gene2Phenotype | `gene2phenotype` | Gene-phenotype associations | Confidence level mapping |
| Genomics England PanelApp | `genomics_england` | Panel gene annotations | Panel confidence |
| UniProt Variants | `uniprot_variants` | UniProt variant annotations | Fixed score |
| UniProt Literature | `uniprot_literature` | UniProt literature links | Fixed score |
| Orphanet | `orphanet` | Rare disease associations | Association type |
| ClinGen | `clingen` | ClinGen gene curation | Classification level |

#### 2. Somatic Mutations

| Data Source | ID | Description |
|-------------|-----|-------------|
| Cancer Gene Census | `cancer_gene_census` | COSMIC Cancer Gene Census |
| IntOGen | `intogen` | Cancer driver gene database |
| Cancer Biomarkers | `cancer_biomarkers` | Biomarker associations |

#### 3. Known Drugs

| Data Source | ID | Description | Scoring |
|-------------|-----|-------------|---------|
| ChEMBL | `chembl` | Clinical trial and drug data | Phase-based (0.05-1.0) |

#### 4. Pathways & Systems Biology

| Data Source | ID | Description |
|-------------|-----|-------------|
| CRISPR Screens | `crispr` | Genome-wide CRISPR data |
| SLAPenrich | `slapenrich` | Pathway enrichment |
| PROGENy | `progeny` | Pathway activity inference |
| Reactome | `reactome` | Pathway annotations |
| SysBio | `sysbio` | Systems biology data |

#### 5. RNA Expression

| Data Source | ID | Description |
|-------------|-----|-------------|
| Expression Atlas | `expression_atlas` | Differential expression |

#### 6. Animal Models

| Data Source | ID | Description |
|-------------|-----|-------------|
| IMPC | `impc` | Mouse phenotype data |

#### 7. Literature

| Data Source | ID | Description |
|-------------|-----|-------------|
| Europe PMC | `europepmc` | Text mining results |

### 2.4 Association Score Calculation

#### Harmonic Sum Formula

```
overall_score = sum(evidence_score[i] / i^2) / harmonic_sum_constant

Where:
- evidence_score[i] = score of i-th evidence (sorted descending)
- i = position in sorted list (1-indexed)
- harmonic_sum_constant ≈ 1.644 (sum of 1/n^2 for n=1 to infinity)
```

#### Data Source Weights

| Data Source | Weight | Rationale |
|-------------|--------|-----------|
| Europe PMC | 0.2 | Text mining, lower confidence |
| Expression Atlas | 0.2 | Expression correlation |
| IMPC | 0.2 | Animal model extrapolation |
| PROGENy | 0.5 | Pathway inference |
| SLAPenrich | 0.5 | Pathway enrichment |
| Cancer Biomarkers | 0.5 | Biomarker associations |
| SysBio | 0.5 | Systems biology |
| OTAR Projects | 0.5 | Open Targets research |
| **All Others** | **1.0** | Full weight (genetics, drugs, etc.) |

#### Evidence Scoring by Source

| Source | Score Range | Formula |
|--------|-------------|---------|
| ChEMBL | 0.05-1.0 | `phase / 4` (phase 0.5 = 0.05) |
| OT Genetics | 0.05-1.0 | L2G score (minimum 0.05) |
| Gene Burden | 0.25-1.0 | Scaled p-value |
| ClinVar | 0.0-1.0 | Significance + review status |
| Gene2Phenotype | 0.25-1.0 | Confidence mapping |
| IntOGen | Fixed | Based on evidence type |

### 2.5 GraphQL Query Examples

#### Query Target with Disease Associations

```graphql
query TargetAssociations($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    biotype

    genomicLocation {
      chromosome
      start
      end
      strand
    }

    tractability {
      smallmolecule {
        label
        value
      }
      antibody {
        label
        value
      }
    }

    associatedDiseases(page: {index: 0, size: 25}) {
      count
      rows {
        disease {
          id
          name
          therapeuticAreas {
            id
            name
          }
        }
        score
        datatypeScores {
          id
          score
        }
        evidenceCount
      }
    }

    knownDrugs(page: {index: 0, size: 10}) {
      count
      rows {
        drugId
        drug {
          name
          drugType
          maximumClinicalTrialPhase
          isApproved
        }
        phase
        mechanismOfAction
      }
    }
  }
}

# Variables
{
  "ensemblId": "ENSG00000157764"
}
```

#### Query Drug Details

```graphql
query DrugDetails($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    name
    synonyms
    tradeNames
    drugType

    maximumClinicalTrialPhase
    isApproved
    hasBeenWithdrawn

    mechanismsOfAction {
      mechanismOfAction
      actionType
      targets {
        id
        approvedSymbol
      }
    }

    linkedTargets {
      count
      rows {
        id
        approvedSymbol
        approvedName
      }
    }

    linkedDiseases {
      count
      rows {
        id
        name
      }
    }

    indications {
      count
      rows {
        disease {
          id
          name
        }
        maxPhaseForIndication
      }
    }

    adverseEvents(page: {index: 0, size: 10}) {
      count
      rows {
        name
        count
        logLR
      }
    }
  }
}

# Variables
{
  "chemblId": "CHEMBL1336"
}
```

#### Query Evidence for Target-Disease Pair

```graphql
query EvidenceQuery(
  $ensemblId: String!
  $efoId: String!
  $size: Int!
) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    evidences(
      efoIds: [$efoId]
      size: $size
    ) {
      count
      rows {
        id
        score
        datasourceId
        datatypeId

        # Genetic evidence fields
        variantId
        variantRsId
        studyId
        beta
        oddsRatio
        pValueMantissa
        pValueExponent

        # Drug evidence fields
        drug {
          id
          name
        }
        clinicalPhase
        clinicalStatus

        # Literature
        literature

        # Clinical significance
        clinicalSignificances
      }
    }
  }
}

# Variables
{
  "ensemblId": "ENSG00000157764",
  "efoId": "EFO_0000616",
  "size": 50
}
```

#### Search Query

```graphql
query SearchQuery($queryString: String!, $size: Int!) {
  search(
    queryString: $queryString
    entityNames: ["target", "disease", "drug"]
    page: {index: 0, size: $size}
  ) {
    total
    hits {
      id
      entity
      name
      description
      score
    }
  }
}

# Variables
{
  "queryString": "BRAF melanoma",
  "size": 20
}
```

### 2.6 Sample JSON Response

```json
{
  "data": {
    "target": {
      "id": "ENSG00000157764",
      "approvedSymbol": "BRAF",
      "approvedName": "B-Raf proto-oncogene, serine/threonine kinase",
      "biotype": "protein_coding",
      "genomicLocation": {
        "chromosome": "7",
        "start": 140719327,
        "end": 140924929,
        "strand": -1
      },
      "tractability": {
        "smallmolecule": [
          {
            "label": "Phase 4",
            "value": true
          },
          {
            "label": "Clinical Candidate",
            "value": true
          },
          {
            "label": "Discovery Precedence",
            "value": true
          }
        ],
        "antibody": [
          {
            "label": "Predicted Tractable - High confidence",
            "value": false
          }
        ]
      },
      "associatedDiseases": {
        "count": 245,
        "rows": [
          {
            "disease": {
              "id": "EFO_0000616",
              "name": "melanoma",
              "therapeuticAreas": [
                {
                  "id": "MONDO_0045024",
                  "name": "cancer or benign tumor"
                }
              ]
            },
            "score": 0.89,
            "datatypeScores": [
              {
                "id": "genetic_association",
                "score": 0.72
              },
              {
                "id": "somatic_mutation",
                "score": 0.95
              },
              {
                "id": "known_drug",
                "score": 1.0
              },
              {
                "id": "affected_pathway",
                "score": 0.45
              },
              {
                "id": "literature",
                "score": 0.68
              }
            ],
            "evidenceCount": 1247
          },
          {
            "disease": {
              "id": "EFO_0003060",
              "name": "non-small cell lung carcinoma",
              "therapeuticAreas": [
                {
                  "id": "MONDO_0045024",
                  "name": "cancer or benign tumor"
                }
              ]
            },
            "score": 0.76,
            "datatypeScores": [
              {
                "id": "genetic_association",
                "score": 0.55
              },
              {
                "id": "somatic_mutation",
                "score": 0.82
              },
              {
                "id": "known_drug",
                "score": 0.75
              }
            ],
            "evidenceCount": 523
          }
        ]
      },
      "knownDrugs": {
        "count": 12,
        "rows": [
          {
            "drugId": "CHEMBL1336",
            "drug": {
              "name": "VEMURAFENIB",
              "drugType": "Small molecule",
              "maximumClinicalTrialPhase": 4,
              "isApproved": true
            },
            "phase": 4,
            "mechanismOfAction": "BRAF protein kinase inhibitor"
          },
          {
            "drugId": "CHEMBL2028663",
            "drug": {
              "name": "DABRAFENIB",
              "drugType": "Small molecule",
              "maximumClinicalTrialPhase": 4,
              "isApproved": true
            },
            "phase": 4,
            "mechanismOfAction": "BRAF protein kinase inhibitor"
          },
          {
            "drugId": "CHEMBL3301610",
            "drug": {
              "name": "ENCORAFENIB",
              "drugType": "Small molecule",
              "maximumClinicalTrialPhase": 4,
              "isApproved": true
            },
            "phase": 4,
            "mechanismOfAction": "BRAF protein kinase inhibitor"
          }
        ]
      }
    }
  }
}
```

### 2.7 Parquet Bulk Download Schemas

**FTP Base URL**: `https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{version}/output/`

#### Targets Dataset (`targets/`)

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Ensembl gene ID |
| `approvedSymbol` | string | HGNC symbol |
| `approvedName` | string | Full gene name |
| `biotype` | string | Gene biotype |
| `chromosome` | string | Chromosome number |
| `start` | long | Genomic start position |
| `end` | long | Genomic end position |
| `strand` | int | Strand (+1/-1) |
| `tpiClinical` | double | Clinical tractability priority |
| `tpiPreclinical` | double | Preclinical tractability priority |
| `functionDescriptions` | array[string] | Function descriptions |
| `synonyms` | array[string] | Gene name synonyms |
| `proteinIds` | array[struct] | UniProt identifiers |
| `dbXrefs` | array[struct] | External database cross-references |

#### Diseases Dataset (`diseases/`)

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | EFO disease ID |
| `name` | string | Disease name |
| `description` | string | Disease description |
| `therapeuticAreas` | array[string] | Therapeutic area IDs |
| `synonyms` | array[struct] | Disease synonyms |
| `ancestors` | array[string] | Ancestor EFO IDs |
| `descendants` | array[string] | Descendant EFO IDs |
| `parents` | array[string] | Direct parent IDs |
| `children` | array[string] | Direct child IDs |
| `dbXRefs` | array[string] | External cross-references |

#### Drugs/Molecule Dataset (`molecule/`)

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | ChEMBL compound ID |
| `name` | string | Drug name |
| `synonyms` | array[string] | Name synonyms |
| `tradeNames` | array[string] | Trade/brand names |
| `drugType` | string | Drug modality |
| `maximumClinicalTrialPhase` | double | Max clinical phase |
| `isApproved` | boolean | Approved drug flag |
| `hasBeenWithdrawn` | boolean | Withdrawn flag |
| `blackBoxWarning` | boolean | Black box warning flag |
| `description` | string | Drug description |
| `parentId` | string | Parent molecule ChEMBL ID |
| `childChemblIds` | array[string] | Child molecule IDs |

#### Evidence Dataset (`evidence/`)

| Field | Type | Description |
|-------|------|-------------|
| `targetId` | string | Ensembl gene ID |
| `diseaseId` | string | EFO disease ID |
| `datasourceId` | string | Data source name |
| `datatypeId` | string | Data type category |
| `score` | double | Evidence score (0-1) |
| `literature` | array[string] | PubMed IDs |
| `variantId` | string | Variant ID |
| `variantRsId` | string | dbSNP rs ID |
| `studyId` | string | Study identifier |
| `beta` | double | Effect size (beta) |
| `oddsRatio` | double | Odds ratio |
| `pValueMantissa` | double | P-value mantissa |
| `pValueExponent` | int | P-value exponent |
| `clinicalPhase` | double | Clinical trial phase |
| `clinicalStatus` | string | Trial status |
| `clinicalSignificances` | array[string] | Clinical significance |
| `confidence` | string | Confidence level |

#### Associations Dataset (Various)

**Overall Associations** (`associationByOverallDirect/`):

| Field | Type | Description |
|-------|------|-------------|
| `targetId` | string | Ensembl gene ID |
| `diseaseId` | string | EFO disease ID |
| `score` | double | Overall association score |
| `evidenceCount` | long | Number of evidence items |

**By Data Type** (`associationByDatatypeDirect/`):

| Field | Type | Description |
|-------|------|-------------|
| `targetId` | string | Ensembl gene ID |
| `diseaseId` | string | EFO disease ID |
| `datatypeId` | string | Data type identifier |
| `score` | double | Data type score |
| `evidenceCount` | long | Evidence count for type |

**By Data Source** (`associationByDatasourceDirect/`):

| Field | Type | Description |
|-------|------|-------------|
| `targetId` | string | Ensembl gene ID |
| `diseaseId` | string | EFO disease ID |
| `datasourceId` | string | Data source identifier |
| `score` | double | Data source score |
| `evidenceCount` | long | Evidence count for source |

### 2.8 Identifier Systems

| Entity | ID System | Format | Example |
|--------|-----------|--------|---------|
| Target/Gene | Ensembl Gene ID | `ENSG{11 digits}` | `ENSG00000157764` |
| Disease | EFO ID | `EFO_{7 digits}` | `EFO_0000616` |
| Drug | ChEMBL ID | `CHEMBL{digits}` | `CHEMBL1336` |
| Variant | dbSNP rsID | `rs{digits}` | `rs113488022` |
| Study | GCST ID | `GCST{6 digits}` | `GCST006979` |
| Pathway | Reactome ID | `R-HSA-{digits}` | `R-HSA-162582` |
| GO Term | GO ID | `GO:{7 digits}` | `GO:0006468` |

---

## Part 3: Cross-References and Integration

### 3.1 Common Identifier Mapping

| Entity Type | DGIdb | Open Targets | Mapping Key |
|-------------|-------|--------------|-------------|
| Gene | HGNC ID | Ensembl Gene ID | UniProt, HGNC symbol |
| Drug | RxNorm, ChEMBL | ChEMBL | ChEMBL ID |
| Disease | - | EFO ID | - |
| Protein | - | UniProt | UniProt accession |

### 3.2 Cross-Reference Strategy

```
DGIdb Gene (hgnc:1097) ←→ UniProt (P15056) ←→ Open Targets (ENSG00000157764)
                         ↓
                   Ensembl mapping

DGIdb Drug (chembl:CHEMBL1336) ←→ Open Targets Drug (CHEMBL1336)
                                  ↓
                            Direct ChEMBL ID
```

### 3.3 Integration Code Examples

#### Python: Query Both APIs

```python
import requests
from typing import Dict, List, Optional

class DrugTargetIntegration:
    """Integration layer for DGIdb and Open Targets APIs."""

    DGIDB_ENDPOINT = "https://dgidb.org/api/graphql"
    OPENTARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

    def query_dgidb(self, gene_name: str) -> Dict:
        """Query DGIdb for drug-gene interactions."""
        query = """
        query GetGeneInteractions($geneName: [String!]!) {
            genes(names: $geneName) {
                nodes {
                    name
                    conceptId
                    geneCategories { name }
                    interactions {
                        drug {
                            name
                            conceptId
                            approvalStatus
                        }
                        interactionTypes {
                            type
                            directionality
                        }
                        score
                        sources { sourceDbName }
                    }
                }
            }
        }
        """

        response = requests.post(
            self.DGIDB_ENDPOINT,
            json={"query": query, "variables": {"geneName": [gene_name]}}
        )
        return response.json()

    def query_opentargets(self, ensembl_id: str) -> Dict:
        """Query Open Targets for target-disease associations."""
        query = """
        query TargetInfo($ensemblId: String!) {
            target(ensemblId: $ensemblId) {
                id
                approvedSymbol
                approvedName
                associatedDiseases(page: {index: 0, size: 25}) {
                    count
                    rows {
                        disease { id name }
                        score
                        datatypeScores { id score }
                    }
                }
                knownDrugs(page: {index: 0, size: 10}) {
                    count
                    rows {
                        drugId
                        drug { name isApproved }
                        mechanismOfAction
                        phase
                    }
                }
            }
        }
        """

        response = requests.post(
            self.OPENTARGETS_ENDPOINT,
            json={"query": query, "variables": {"ensemblId": ensembl_id}}
        )
        return response.json()

    def get_hgnc_to_ensembl_mapping(self, hgnc_symbol: str) -> Optional[str]:
        """Map HGNC symbol to Ensembl ID via Open Targets search."""
        query = """
        query SearchGene($symbol: String!) {
            search(queryString: $symbol, entityNames: ["target"], page: {index: 0, size: 5}) {
                hits {
                    id
                    name
                }
            }
        }
        """

        response = requests.post(
            self.OPENTARGETS_ENDPOINT,
            json={"query": query, "variables": {"symbol": hgnc_symbol}}
        )
        data = response.json()

        hits = data.get("data", {}).get("search", {}).get("hits", [])
        for hit in hits:
            if hit.get("name", "").upper() == hgnc_symbol.upper():
                return hit.get("id")
        return None

    def integrated_query(self, gene_name: str) -> Dict:
        """Get comprehensive drug-target data from both sources."""

        # Get DGIdb data
        dgidb_data = self.query_dgidb(gene_name)

        # Map to Ensembl ID and get Open Targets data
        ensembl_id = self.get_hgnc_to_ensembl_mapping(gene_name)
        opentargets_data = None
        if ensembl_id:
            opentargets_data = self.query_opentargets(ensembl_id)

        return {
            "gene_name": gene_name,
            "ensembl_id": ensembl_id,
            "dgidb": dgidb_data,
            "opentargets": opentargets_data
        }


# Example usage
if __name__ == "__main__":
    integration = DrugTargetIntegration()

    # Query BRAF
    result = integration.integrated_query("BRAF")

    # Extract DGIdb interactions
    dgidb_genes = result["dgidb"]["data"]["genes"]["nodes"]
    if dgidb_genes:
        print(f"DGIdb Interactions for {result['gene_name']}:")
        for interaction in dgidb_genes[0]["interactions"][:5]:
            drug = interaction["drug"]
            types = [t["type"] for t in interaction["interactionTypes"]]
            print(f"  - {drug['name']} ({drug['approvalStatus']}): {', '.join(types)}")

    # Extract Open Targets associations
    if result["opentargets"]:
        ot_target = result["opentargets"]["data"]["target"]
        print(f"\nOpen Targets Disease Associations for {ot_target['approvedSymbol']}:")
        for assoc in ot_target["associatedDiseases"]["rows"][:5]:
            disease = assoc["disease"]
            print(f"  - {disease['name']} (score: {assoc['score']:.3f})")
```

#### JavaScript/TypeScript: Combined Query

```typescript
interface DGIdbInteraction {
  drug: {
    name: string;
    conceptId: string;
    approvalStatus: string;
  };
  interactionTypes: Array<{
    type: string;
    directionality: string;
  }>;
  score: number;
}

interface OpenTargetsAssociation {
  disease: {
    id: string;
    name: string;
  };
  score: number;
  datatypeScores: Array<{
    id: string;
    score: number;
  }>;
}

class DrugTargetAPI {
  private readonly DGIDB_ENDPOINT = 'https://dgidb.org/api/graphql';
  private readonly OPENTARGETS_ENDPOINT = 'https://api.platform.opentargets.org/api/v4/graphql';

  async queryDGIdb(geneName: string): Promise<DGIdbInteraction[]> {
    const query = `
      query GetGeneInteractions($geneName: [String!]!) {
        genes(names: $geneName) {
          nodes {
            name
            interactions {
              drug {
                name
                conceptId
                approvalStatus
              }
              interactionTypes {
                type
                directionality
              }
              score
            }
          }
        }
      }
    `;

    const response = await fetch(this.DGIDB_ENDPOINT, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        query,
        variables: { geneName: [geneName] }
      })
    });

    const data = await response.json();
    return data.data?.genes?.nodes?.[0]?.interactions || [];
  }

  async queryOpenTargets(ensemblId: string): Promise<OpenTargetsAssociation[]> {
    const query = `
      query TargetAssociations($ensemblId: String!) {
        target(ensemblId: $ensemblId) {
          associatedDiseases(page: {index: 0, size: 50}) {
            rows {
              disease { id name }
              score
              datatypeScores { id score }
            }
          }
        }
      }
    `;

    const response = await fetch(this.OPENTARGETS_ENDPOINT, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        query,
        variables: { ensemblId }
      })
    });

    const data = await response.json();
    return data.data?.target?.associatedDiseases?.rows || [];
  }

  async getComprehensiveData(geneName: string, ensemblId: string) {
    const [interactions, associations] = await Promise.all([
      this.queryDGIdb(geneName),
      this.queryOpenTargets(ensemblId)
    ]);

    return {
      geneName,
      ensemblId,
      drugInteractions: interactions,
      diseaseAssociations: associations,
      summary: {
        totalDrugs: interactions.length,
        approvedDrugs: interactions.filter(i => i.drug.approvalStatus === 'Approved').length,
        totalDiseases: associations.length,
        topDisease: associations[0]?.disease.name || null
      }
    };
  }
}

// Usage
const api = new DrugTargetAPI();
api.getComprehensiveData('BRAF', 'ENSG00000157764')
  .then(data => {
    console.log(`Gene: ${data.geneName}`);
    console.log(`Total Drug Interactions: ${data.summary.totalDrugs}`);
    console.log(`Approved Drugs: ${data.summary.approvedDrugs}`);
    console.log(`Disease Associations: ${data.summary.totalDiseases}`);
    console.log(`Top Disease: ${data.summary.topDisease}`);
  });
```

#### R: Combined Analysis

```r
library(httr)
library(jsonlite)
library(dplyr)

# DGIdb Query Function
query_dgidb <- function(gene_name) {
  query <- '
    query GetGeneInteractions($geneName: [String!]!) {
      genes(names: $geneName) {
        nodes {
          name
          conceptId
          interactions {
            drug {
              name
              approvalStatus
            }
            interactionTypes {
              type
            }
            score
          }
        }
      }
    }
  '

  response <- POST(
    "https://dgidb.org/api/graphql",
    body = list(
      query = query,
      variables = list(geneName = list(gene_name))
    ),
    encode = "json"
  )

  content(response, as = "parsed")
}

# Open Targets Query Function
query_opentargets <- function(ensembl_id) {
  query <- '
    query TargetInfo($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        approvedSymbol
        associatedDiseases(page: {index: 0, size: 25}) {
          rows {
            disease { id name }
            score
          }
        }
        knownDrugs(page: {index: 0, size: 10}) {
          rows {
            drug { name isApproved }
            phase
          }
        }
      }
    }
  '

  response <- POST(
    "https://api.platform.opentargets.org/api/v4/graphql",
    body = list(
      query = query,
      variables = list(ensemblId = ensembl_id)
    ),
    encode = "json"
  )

  content(response, as = "parsed")
}

# Example: Analyze BRAF
dgidb_result <- query_dgidb("BRAF")
ot_result <- query_opentargets("ENSG00000157764")

# Extract and combine drug information
dgidb_drugs <- dgidb_result$data$genes$nodes[[1]]$interactions %>%
  lapply(function(x) {
    data.frame(
      source = "DGIdb",
      drug_name = x$drug$name,
      approval_status = x$drug$approvalStatus,
      interaction_score = x$score
    )
  }) %>%
  bind_rows()

ot_drugs <- ot_result$data$target$knownDrugs$rows %>%
  lapply(function(x) {
    data.frame(
      source = "OpenTargets",
      drug_name = x$drug$name,
      is_approved = x$drug$isApproved,
      max_phase = x$phase
    )
  }) %>%
  bind_rows()

print("DGIdb Drug Interactions:")
print(head(dgidb_drugs))

print("\nOpen Targets Known Drugs:")
print(head(ot_drugs))
```

---

## Part 4: Data Access Summary

### 4.1 API Endpoints

| Database | Endpoint | Type |
|----------|----------|------|
| DGIdb | `https://dgidb.org/api/graphql` | GraphQL |
| Open Targets | `https://api.platform.opentargets.org/api/v4/graphql` | GraphQL |

### 4.2 Bulk Downloads

| Database | URL | Formats |
|----------|-----|---------|
| DGIdb | https://dgidb.org/downloads | TSV |
| Open Targets | https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/ | Parquet, JSON |

### 4.3 Rate Limits and Best Practices

| Aspect | DGIdb | Open Targets |
|--------|-------|--------------|
| Rate Limit | No explicit limit | 100 requests/minute recommended |
| Batch Size | Up to 100 genes per query | Page size up to 10,000 |
| Caching | Recommended | Recommended |
| Bulk Data | Preferred for large analyses | Parquet files for full dataset |

### 4.4 License Information

| Database | License | Commercial Use |
|----------|---------|----------------|
| DGIdb 5.0 | MIT | Allowed |
| Open Targets | Apache 2.0 | Allowed |

---

## References

1. **DGIdb 5.0**: Freshour SL, et al. "DGIdb 5.0: rebuilding the drug-gene interaction database" Nucleic Acids Res. 2024;52(D1):D1227-D1235. https://doi.org/10.1093/nar/gkad1040

2. **Open Targets Platform**: Ochoa D, et al. "Open Targets Platform: supporting systematic drug-target identification and prioritisation" Nucleic Acids Res. 2021;49(D1):D1302-D1310. https://doi.org/10.1093/nar/gkaa1027

3. **Open Targets Genetics**: Ghoussaini M, et al. "Open Targets Genetics: systematic identification of trait-associated genes using large-scale genetics and functional genomics" Nucleic Acids Res. 2021;49(D1):D1311-D1320. https://doi.org/10.1093/nar/gkaa840

---

## Version History

| Date | Version | Changes |
|------|---------|---------|
| 2024-01 | 1.0 | Initial comprehensive schema documentation |
| 2024-01 | 1.1 | Added integration code examples and cross-references |
