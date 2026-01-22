---
id: schema-dgidb-open-targets
title: "DGIdb and Open Targets Platform Schema Reference"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, drug-gene-interactions, graphql, open-targets]
---

**Parent:** [Schema Documentation](./_index.md)

# DGIdb and Open Targets Platform Schema Reference

This document provides schema information for two major drug-gene interaction resources:
1. **DGIdb (Drug Gene Interaction Database)**
2. **Open Targets Platform**

## DGIdb (Drug Gene Interaction Database)

### Overview
DGIdb aggregates drug-gene interaction data from multiple sources, providing information about how drugs interact with genes and their products.

### GraphQL API
- **Endpoint**: `https://dgidb.org/api/graphql`
- **Version**: v5

### Core Types

#### Gene
```graphql
type Gene {
  id: ID!
  name: String!
  conceptId: String!
  longName: String
  description: String
  geneCategories: [GeneCategory!]!
  geneAliases: [GeneAlias!]!
  geneClaims: [GeneClaim!]!
  interactions: [Interaction!]!
}
```

#### Drug
```graphql
type Drug {
  id: ID!
  name: String!
  conceptId: String!
  approved: Boolean!
  immunotherapy: Boolean!
  antineoplastic: Boolean!
  drugAliases: [DrugAlias!]!
  drugAttributes: [DrugAttribute!]!
  drugApplications: [DrugApplication!]!
  drugApprovalRatings: [DrugApprovalRating!]!
  interactions: [Interaction!]!
}
```

#### Interaction
```graphql
type Interaction {
  id: ID!
  gene: Gene!
  drug: Drug!
  interactionScore: Float
  interactionTypes: [InteractionType!]!
  interactionClaims: [InteractionClaim!]!
  publications: [Publication!]!
  sources: [Source!]!
}
```

#### InteractionType
```graphql
type InteractionType {
  id: ID!
  type: String!
  directionality: String
}
```

#### InteractionClaim
```graphql
type InteractionClaim {
  id: ID!
  source: Source!
  interactionClaimTypes: [InteractionClaimType!]!
  drug: Drug!
  gene: Gene!
  publications: [Publication!]!
}
```

#### Source
```graphql
type Source {
  id: ID!
  fullName: String!
  sourceDbName: String!
  sourceDbVersion: String
  citationShort: String
  citationFull: String
  pmid: String
  pmcid: String
  doi: String
  sourceUrl: String
  license: String
  licenseLink: String
}
```

#### Publication
```graphql
type Publication {
  id: ID!
  pmid: Int
  citation: String
}
```

#### GeneCategory
```graphql
type GeneCategory {
  id: ID!
  name: String!
}
```

#### DrugApplication
```graphql
type DrugApplication {
  id: ID!
  appNo: String!
  appType: String
}
```

#### DrugApprovalRating
```graphql
type DrugApprovalRating {
  id: ID!
  rating: String!
  source: Source!
}
```

### Query Examples

#### Search for Gene Interactions
```graphql
query {
  genes(names: ["EGFR"]) {
    nodes {
      name
      longName
      interactions {
        drug {
          name
          approved
        }
        interactionTypes {
          type
          directionality
        }
        interactionScore
      }
    }
  }
}
```

#### Search for Drug Interactions
```graphql
query {
  drugs(names: ["Gefitinib"]) {
    nodes {
      name
      approved
      antineoplastic
      interactions {
        gene {
          name
          longName
        }
        interactionTypes {
          type
          directionality
        }
        sources {
          fullName
        }
      }
    }
  }
}
```

#### Search by Interaction Type
```graphql
query {
  interactions(
    interactionTypes: ["inhibitor"]
    geneName: "BRAF"
  ) {
    nodes {
      drug {
        name
        approved
      }
      gene {
        name
      }
      interactionTypes {
        type
      }
      publications {
        pmid
        citation
      }
    }
  }
}
```

### Interaction Types

Common interaction types in DGIdb:
- **inhibitor**: Drug inhibits gene/protein activity
- **activator**: Drug activates gene/protein activity
- **antagonist**: Drug blocks receptor activity
- **agonist**: Drug activates receptor
- **antibody**: Antibody-based therapeutic
- **modulator**: Drug modulates activity
- **binder**: Drug binds to target
- **inducer**: Drug induces expression
- **suppressor**: Drug suppresses expression
- **blocker**: Drug blocks function

### Gene Categories

Common gene categories:
- **KINASE**: Protein kinase genes
- **TUMOR SUPPRESSOR**: Tumor suppressor genes
- **TRANSCRIPTION FACTOR**: Transcription factor genes
- **CELL SURFACE**: Cell surface protein genes
- **DRUG RESISTANCE**: Drug resistance genes
- **DRUGGABLE GENOME**: Potentially druggable genes
- **CLINICALLY ACTIONABLE**: Clinically relevant genes

---

## Open Targets Platform

### Overview
The Open Targets Platform integrates evidence from multiple data sources to support systematic identification and prioritization of drug targets.

### GraphQL API
- **Endpoint**: `https://api.platform.opentargets.org/api/v4/graphql`
- **Version**: v4

### Core Types

#### Target
```graphql
type Target {
  id: ID!
  approvedSymbol: String!
  approvedName: String!
  biotype: String!
  geneticConstraint: GeneticConstraint
  genomicLocation: GenomicLocation!
  hallmarks: Hallmarks
  synonyms: [String!]
  symbolSynonyms: [String!]
  nameSynonyms: [String!]
  functionDescriptions: [String!]
  pathways: [Pathway!]!
  go: [GO!]!
  proteinAnnotations: ProteinAnnotations
  subcellularLocations: [SubcellularLocation!]
  tractability: Tractability
  safetyLiabilities: [SafetyLiability!]
  tep: TEP
  chemicalProbes: [ChemicalProbe!]
  knownDrugs: [KnownDrug!]
  interactions: [Interaction!]
}
```

#### Disease
```graphql
type Disease {
  id: ID!
  name: String!
  description: String
  synonyms: [String!]
  parents: [String!]
  ancestors: [String!]
  descendants: [String!]
  therapeuticAreas: [TherapeuticArea!]!
  ontology: Ontology!
  indirectLocations: [Tissue!]
  directLocations: [Tissue!]
}
```

#### Evidence
```graphql
type Evidence {
  id: ID!
  target: Target!
  disease: Disease!
  score: Float!
  datasourceId: String!
  datatypeId: String!
  publications: [Publication!]
  studyId: String
  cohortId: String
  clinicalPhase: Int
  clinicalStatus: String
  urls: [URL!]
}
```

#### Association
```graphql
type Association {
  id: ID!
  target: Target!
  disease: Disease!
  score: Float!
  datasourceScores: [DatasourceScore!]!
  datatypeScores: [DatatypeScore!]!
}
```

#### Drug
```graphql
type Drug {
  id: ID!
  name: String!
  synonyms: [String!]
  tradeNames: [String!]
  yearOfFirstApproval: Int
  drugType: String!
  maximumClinicalTrialPhase: Int
  hasBeenWithdrawn: Boolean!
  withdrawnNotice: WithdrawnNotice
  linkedTargets: [Target!]!
  linkedDiseases: [Disease!]!
  childChemblIds: [String!]
  knownDrugs: [KnownDrug!]
  mechanismsOfAction: [MechanismOfAction!]
  indications: [Indication!]
  adverseEvents: [AdverseEvent!]
}
```

#### KnownDrug
```graphql
type KnownDrug {
  approvedSymbol: String!
  approvedName: String!
  prefName: String!
  drugType: String!
  drugId: String!
  phase: Int!
  status: String
  urls: [URL!]
  mechanismOfAction: String
  targets: [Target!]!
  disease: Disease!
  references: [Reference!]
}
```

#### MechanismOfAction
```graphql
type MechanismOfAction {
  mechanismOfAction: String!
  actionType: String
  targetName: String
  targets: [Target!]!
  references: [Reference!]
}
```

#### GenomicLocation
```graphql
type GenomicLocation {
  chromosome: String!
  start: Long!
  end: Long!
  strand: Int!
}
```

#### Pathway
```graphql
type Pathway {
  pathway: String!
  pathwayId: String!
}
```

#### GO
```graphql
type GO {
  id: String!
  term: String!
  aspect: String!
}
```

#### Tractability
```graphql
type Tractability {
  smallmolecule: TractabilityAssessment
  antibody: TractabilityAssessment
  otherModalities: TractabilityAssessment
}

type TractabilityAssessment {
  topCategory: String!
  buckets: [Int!]!
}
```

#### SafetyLiability
```graphql
type SafetyLiability {
  event: String!
  effects: [SafetyEffect!]!
  biosample: String
  literature: String
}
```

### Query Examples

#### Search for Target-Disease Associations
```graphql
query {
  target(ensemblId: "ENSG00000146648") {
    approvedSymbol
    approvedName
    associatedDiseases(page: { size: 10 }) {
      count
      rows {
        disease {
          id
          name
        }
        score
        datatypeScores {
          componentId
          score
        }
      }
    }
  }
}
```

#### Get Known Drugs for a Target
```graphql
query {
  target(ensemblId: "ENSG00000146648") {
    approvedSymbol
    knownDrugs(size: 20) {
      uniqueDrugs
      uniqueDiseases
      count
      rows {
        drug {
          id
          name
          drugType
        }
        phase
        status
        mechanismOfAction
        disease {
          id
          name
        }
      }
    }
  }
}
```

#### Search for Evidence
```graphql
query {
  disease(efoId: "EFO_0000685") {
    name
    associatedTargets(page: { size: 10 }) {
      count
      rows {
        target {
          approvedSymbol
          approvedName
        }
        score
        datasourceScores {
          componentId
          score
        }
      }
    }
  }
}
```

#### Drug Mechanism of Action
```graphql
query {
  drug(chemblId: "CHEMBL941") {
    name
    mechanismsOfAction {
      mechanismOfAction
      actionType
      targets {
        approvedSymbol
        approvedName
      }
    }
  }
}
```

### Data Sources

#### Genetic Evidence
- **GWAS Catalog**: Genome-wide association studies
- **PheWAS**: Phenome-wide association studies
- **Europe PMC**: Text mining from literature
- **UniProt Literature**: Curated literature
- **Gene2Phenotype**: Clinical genetic databases
- **GeneBass**: Gene-based association studies
- **ClinGen**: Clinical genome resource

#### Somatic Evidence
- **Cancer Gene Census**: Curated cancer genes
- **IntOGen**: Cancer driver genes
- **EVA Somatic**: Somatic variants
- **Intogen Cancerdrivers**: Cancer driver mutations

#### Drug Evidence
- **ChEMBL**: Bioactivity database
- **DrugBank**: Drug information
- **DGIdb**: Drug-gene interactions

#### Pathways & Systems Biology
- **Reactome**: Pathway database
- **SLAPenrich**: Pathway enrichment
- **PROGENy**: Pathway responsive genes
- **CRISPR**: CRISPR screens

#### RNA Expression
- **Expression Atlas**: Gene expression patterns
- **HPA**: Human Protein Atlas

#### Animal Models
- **IMPC**: International Mouse Phenotyping Consortium
- **MGI**: Mouse Genome Informatics

### Evidence Scores

Evidence scores range from 0 to 1:
- **>0.7**: Strong evidence
- **0.5-0.7**: Moderate evidence
- **0.3-0.5**: Weak evidence
- **<0.3**: Very weak evidence

### Clinical Phases
- **0**: Preclinical
- **1**: Phase I
- **2**: Phase II
- **3**: Phase III
- **4**: Phase IV (approved)

## Integration Between DGIdb and Open Targets

Both platforms complement each other:
- **DGIdb**: Focused on drug-gene interactions, comprehensive interaction types
- **Open Targets**: Broader target-disease associations, extensive evidence integration

### Combined Use Case
```
1. Use Open Targets to find targets associated with a disease
2. Use DGIdb to find drugs that interact with those targets
3. Cross-reference known drugs from both platforms
4. Validate with clinical trial data from Open Targets
```

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `conceptId` | Normalized identifier for genes or drugs across sources | `HGNC:3236` |
| `interactionScore` | Numeric score indicating strength of drug-gene interaction | `0.85` |
| `interactionTypes` | Categories describing how drug affects gene/protein | `inhibitor` |
| `approved` | Boolean indicating FDA approval status of drug | `true` |
| `antineoplastic` | Boolean indicating anti-cancer activity | `true` |
| `phase` | Clinical trial phase (0-4) for drug-disease indication | `4` |
| `score` | Evidence score (0-1) for target-disease association | `0.7` |
| `datasourceId` | Identifier for evidence source contributing to association | `chembl` |
| `mechanismOfAction` | Description of how drug produces its therapeutic effect | `EGFR inhibitor` |
| `tractability` | Assessment of target druggability by modality | `smallmolecule` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Drug-Gene Interaction | Relationship where a drug affects gene/protein function | interactionTypes |
| Druggable Genome | Genes encoding proteins amenable to drug targeting | geneCategories |
| Evidence Score | Aggregated confidence score from multiple data sources | score, datasourceScores |
| Therapeutic Area | High-level disease classification (oncology, neurology, etc.) | therapeuticAreas |
| Safety Liability | Known adverse effects associated with target modulation | SafetyLiability |
| Chemical Probe | Small molecule tool compound for target validation | ChemicalProbe |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| DGIdb | Drug Gene Interaction Database | Drug-target interaction resource |
| EGFR | Epidermal Growth Factor Receptor | Common cancer target |
| BRAF | B-Raf Proto-Oncogene | Kinase cancer target |
| GWAS | Genome-Wide Association Study | Genetic association method |
| PheWAS | Phenome-Wide Association Study | Phenotype association method |
| ChEMBL | - | EBI bioactivity database |
| EFO | Experimental Factor Ontology | Disease ontology |
| GO | Gene Ontology | Functional annotation |
| IMPC | International Mouse Phenotyping Consortium | Mouse model resource |
| TEP | Target Enabling Package | Target validation data |

---

## References

### DGIdb
- [DGIdb Website](https://dgidb.org)
- [DGIdb API Documentation](https://dgidb.org/api)
- [DGIdb GitHub](https://github.com/dgidb/dgidb-v5)

### Open Targets
- [Open Targets Platform](https://platform.opentargets.org)
- [Open Targets API Documentation](https://platform-docs.opentargets.org)
- [Open Targets GitHub](https://github.com/opentargets)