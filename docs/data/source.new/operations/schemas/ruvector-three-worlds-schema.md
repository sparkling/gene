---
id: schema-ruvector-three-worlds
title: RuVector THREE WORLDS Schema Design
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, architecture, ruvector]
---

# RuVector THREE WORLDS Schema Design

**Document ID:** RUVECTOR-THREE-WORLDS-SCHEMA
**Status:** Architecture Decision
**Owner:** System Architecture
**Date:** January 21, 2026
**Version:** 1.0
**Parent:** [Schema Documentation](./_index.md)

---

## Executive Summary

This document defines the optimal schema for storing THREE WORLDS data (Genetics, Traditional Medicine, Nutritional Science) using RuVector's native capabilities. The design leverages RuVector's unique features: CodeGraph hyperedges, HNSW vector search, SONA continual learning, and distributed sharding via RuvectorCluster.

### Scale Requirements

| Domain | Entity Type | Estimated Count | Storage Strategy |
|--------|-------------|-----------------|------------------|
| **WORLD 1: Genetics** | SNPs/Variants | ~1.2B | Domain-sharded, quantized vectors |
| **WORLD 1: Genetics** | Genes | ~25K | Hot cache, full precision |
| **WORLD 2: Traditional Medicine** | Compounds (NP) | ~800K | HNSW-indexed, structure vectors |
| **WORLD 2: Traditional Medicine** | Formulas | ~60K | Hyperedge patterns |
| **WORLD 3: Nutrition** | Foods | ~380K | Multi-vector (composition + embedding) |
| **Cross-Domain** | Pathways | ~10K | Graph topology, hyperedge-native |
| **Cross-Domain** | Diseases | ~80K | Ontology hierarchy in graph |

---

## Part 1: RuVector Capability Mapping

### 1.1 Core RuVector Components Used

```
+------------------+     +------------------+     +------------------+
|    VectorDB      |     |    CodeGraph     |     | RuvectorCluster  |
|  (HNSW-indexed)  |     |  (Hyperedges)    |     |  (Distributed)   |
+------------------+     +------------------+     +------------------+
        |                        |                        |
        v                        v                        v
  Vector Search           N-ary Relations           Auto-Sharding
  Similarity Queries      Pathway Modeling          Raft Consensus
  Semantic Embeddings     Formula Composition       Cross-Node Queries
```

### 1.2 Component Selection Rationale

| RuVector Component | THREE WORLDS Use Case | Why This Component |
|--------------------|----------------------|-------------------|
| `VectorDB` | Gene/compound embeddings, semantic search | 150x faster HNSW, metadata filtering |
| `CodeGraph` | Pathways, gene networks, compound-target | Native hyperedges for N-ary relations |
| `RuvectorCluster` | Cross-node queries, domain sharding | Raft consensus, auto-rebalancing |
| `SONA Engine` | Learning from query patterns | Micro-LoRA adaptation (<0.1ms) |
| `NeuralSubstrate` | Semantic drift detection | Memory physics for cache management |
| `GNN Wrapper` | Graph neural network queries | Differentiable search, hierarchy traversal |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "GENE:12345" |
| `name` | string | Entity name | "TP53" |
| `type` | string | Record type | "gene" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## Part 2: Node Type Definitions

### 2.1 WORLD 1: Genetics Node Types

#### Gene Node
```typescript
interface GeneNode {
  // === Identity (stored in CodeGraph) ===
  id: string;                    // Format: "GENE:{ncbi_gene_id}"
  labels: ['Gene', 'WORLD1'];

  properties: {
    // Primary identifiers
    ncbi_gene_id: number;        // Entrez Gene ID (canonical)
    hgnc_id: string;             // HGNC:12345
    hgnc_symbol: string;         // Official symbol (e.g., "TP53")
    ensembl_id: string;          // ENSG00000000000
    uniprot_ids: string[];       // Multiple protein products

    // Genomic location
    chromosome: string;          // 1-22, X, Y, MT
    start_pos: number;           // GRCh38
    end_pos: number;
    strand: '+' | '-';
    cytoband: string;            // e.g., "17q21.31"

    // Annotations
    name: string;
    description: string;
    gene_type: string;           // protein_coding, lncRNA, etc.

    // Pharmacogenomics flags
    cpic_level: string;          // A, B, C, D
    is_druggable: boolean;
    has_haplotypes: boolean;
  };

  // === Embedding (stored in VectorDB) ===
  // Stored separately: GENE_EMB:{ncbi_gene_id}
  // Dimensions: 384 (MiniLM) or 768 (mpnet)
  // Content: Gene description + GO terms + pathway names
}
```

#### Variant Node (SNP/SNV)
```typescript
interface VariantNode {
  id: string;                    // Format: "VAR:{rsid}" or "VAR:{chr}-{pos}-{ref}-{alt}"
  labels: ['Variant', 'WORLD1'];

  properties: {
    // Primary identifiers
    rsid: string;                // rs123456789
    clinvar_vcv: string;         // VCV000123456
    gnomad_id: string;           // chr-pos-ref-alt

    // Location (SPDI format)
    assembly: 'GRCh38';
    chromosome: string;
    position: number;            // 1-based
    ref_allele: string;
    alt_allele: string;

    // HGVS nomenclature
    hgvs_g: string;              // Genomic
    hgvs_c: string;              // Coding
    hgvs_p: string;              // Protein

    // Clinical classification
    clinvar_sig: 'pathogenic' | 'likely_pathogenic' | 'VUS' |
                 'likely_benign' | 'benign' | null;
    clinvar_review: string;      // Review status

    // Functional predictions (quantized to save space)
    cadd_phred: number;          // 0-99
    revel: number;               // 0-1
    spliceai_max: number;        // 0-1

    // Consequence
    consequence: string;         // VEP consequence
    impact: 'HIGH' | 'MODERATE' | 'LOW' | 'MODIFIER';

    // Population frequency (quantized)
    gnomad_af: number;           // Global allele frequency
    gnomad_popmax: number;       // Maximum population AF
  };

  // === Population Frequencies (separate VectorDB entry) ===
  // Format: VAR_POP:{rsid}
  // Vector: [af_afr, af_amr, af_asj, af_eas, af_fin, af_mid, af_nfe, af_sas]
  // Enables: "Find variants with similar population distribution"
}
```

#### Haplotype Node
```typescript
interface HaplotypeNode {
  id: string;                    // Format: "HAP:{gene}:{star_allele}"
  labels: ['Haplotype', 'WORLD1'];

  properties: {
    gene_symbol: string;         // e.g., "CYP2D6"
    star_allele: string;         // e.g., "*4"
    activity_score: number;      // 0-2
    function: string;            // "No function", "Decreased function", etc.
    defining_variants: string[]; // rsIDs
    frequency_eur: number;
    frequency_afr: number;
    frequency_eas: number;
  };
}
```

### 2.2 WORLD 2: Traditional Medicine Node Types

#### Compound Node (Natural Product)
```typescript
interface CompoundNode {
  id: string;                    // Format: "COMP:{inchikey}"
  labels: ['Compound', 'WORLD2'];

  properties: {
    // Structure identifiers (InChIKey is canonical)
    inchikey: string;            // 27-char standard
    inchi: string;               // Full InChI
    smiles_canonical: string;
    smiles_isomeric: string;

    // Cross-references
    pubchem_cid: number;
    chembl_id: string;
    coconut_id: string;
    chebi_id: string;
    drugbank_id: string;

    // Physical properties
    molecular_weight: number;
    exact_mass: number;
    formula: string;

    // Lipinski
    alogp: number;
    hbd: number;                 // H-bond donors
    hba: number;                 // H-bond acceptors
    tpsa: number;                // Topological PSA
    rotatable_bonds: number;

    // Classification
    np_likeness: number;         // Natural product score
    superclass: string;          // ClassyFire
    class: string;
    subclass: string;

    // Names
    preferred_name: string;
    iupac_name: string;
    synonyms: string[];
  };

  // === Structure Embedding (VectorDB) ===
  // Format: COMP_FP:{inchikey}
  // Vector: Morgan fingerprint (2048-bit) or mol2vec (300-dim)
  // Enables: "Find structurally similar compounds"

  // === Semantic Embedding (VectorDB) ===
  // Format: COMP_SEM:{inchikey}
  // Vector: 384-dim from name + class + activities
  // Enables: "Find compounds with similar bioactivity profile"
}
```

#### Organism Node (Source Species)
```typescript
interface OrganismNode {
  id: string;                    // Format: "ORG:{ncbi_taxonomy_id}"
  labels: ['Organism', 'WORLD2'];

  properties: {
    ncbi_taxonomy_id: number;
    scientific_name: string;
    common_names: string[];
    kingdom: string;
    family: string;
    genus: string;
    species: string;

    // Traditional medicine context
    traditional_systems: string[];  // ['TCM', 'Kampo', 'Ayurveda']
    geographic_origin: string[];
    parts_used: string[];           // ['root', 'leaf', 'bark']
  };
}
```

#### Formula Node (Traditional Preparation)
```typescript
interface FormulaNode {
  id: string;                    // Format: "FORM:{db}:{id}"
  labels: ['Formula', 'WORLD2'];

  properties: {
    // Source database
    source_db: 'KampoDB' | 'BATMAN-TCM' | 'IMPPAT';
    source_id: string;

    // Names
    name_original: string;       // In original language/script
    name_pinyin: string;         // Romanized
    name_english: string;

    // Composition (stored as hyperedge - see Part 3)
    component_count: number;

    // Traditional use
    traditional_category: string;
    traditional_indications: string[];

    // Modern validation
    pubmed_citations: number;
    has_clinical_trial: boolean;
  };

  // === Formula Embedding (VectorDB) ===
  // Format: FORM_SEM:{id}
  // Vector: 384-dim from indications + category
  // Enables: "Find formulas for similar conditions"
}
```

### 2.3 WORLD 3: Nutritional Science Node Types

#### Food Node
```typescript
interface FoodNode {
  id: string;                    // Format: "FOOD:{db}:{id}"
  labels: ['Food', 'WORLD3'];

  properties: {
    // Source
    source_db: 'USDA' | 'FooDB' | 'OpenFoodFacts';
    source_id: string;

    // Identity
    name: string;
    description: string;
    food_group: string;
    food_category: string;

    // Nutrient summary (full nutrients in separate store)
    energy_kcal: number;
    protein_g: number;
    carbs_g: number;
    fat_g: number;
    fiber_g: number;

    // Flags
    is_whole_food: boolean;
    is_processed: boolean;
    allergens: string[];
  };

  // === Nutrient Profile Vector (VectorDB) ===
  // Format: FOOD_NUT:{id}
  // Vector: Normalized nutrient profile (50-dim: vitamins, minerals, etc.)
  // Enables: "Find foods with similar nutrient profile"

  // === Phytochemical Profile Vector (VectorDB) ===
  // Format: FOOD_PHYTO:{id}
  // Vector: Phytochemical fingerprint from FooDB
  // Enables: "Find foods with similar bioactive compounds"
}
```

#### Nutrient Node
```typescript
interface NutrientNode {
  id: string;                    // Format: "NUTR:{usda_id}"
  labels: ['Nutrient', 'WORLD3'];

  properties: {
    usda_nutrient_id: number;
    name: string;
    unit: string;

    // Classification
    category: 'Macronutrient' | 'Vitamin' | 'Mineral' |
              'Amino Acid' | 'Fatty Acid' | 'Other';
    subcategory: string;

    // Bioactivity (linked to compounds)
    compound_inchikey: string;   // If mappable to a compound

    // RDA
    rda_male: number;
    rda_female: number;
    upper_limit: number;
  };
}
```

### 2.4 Cross-Domain Node Types

#### Disease Node
```typescript
interface DiseaseNode {
  id: string;                    // Format: "DIS:{mondo_id}"
  labels: ['Disease', 'CrossDomain'];

  properties: {
    // Primary identifier (MONDO is hub)
    mondo_id: string;            // MONDO:0000000

    // Cross-references
    omim_id: string;
    orphanet_id: string;
    doid: string;
    mesh_id: string;
    umls_cui: string;
    icd10: string[];

    // Names
    name: string;
    synonyms: string[];

    // Classification
    category: string;            // Cancer, Cardiovascular, etc.
    is_rare: boolean;
    inheritance: string[];       // AD, AR, XL, etc.

    // Epidemiology
    prevalence: number;          // Per 100,000
  };

  // === Disease Embedding (VectorDB) ===
  // Format: DIS_SEM:{mondo_id}
  // Vector: 384-dim from name + synonyms + phenotypes
  // Enables: "Find similar diseases"
}
```

#### Pathway Node
```typescript
interface PathwayNode {
  id: string;                    // Format: "PATH:{reactome_id}"
  labels: ['Pathway', 'CrossDomain'];

  properties: {
    // Identifiers
    reactome_id: string;         // R-HSA-000000
    wikipathways_id: string;
    kegg_id: string;
    go_bp_id: string;

    // Names
    name: string;
    synonyms: string[];

    // Hierarchy
    top_level_pathway: string;
    parent_pathway: string;

    // Statistics
    gene_count: number;
    reaction_count: number;

    // Disease relevance
    disease_associated: boolean;
  };

  // === Pathway Embedding (VectorDB) ===
  // Derived from member gene embeddings (mean pooling)
}
```

#### Phenotype Node (HPO)
```typescript
interface PhenotypeNode {
  id: string;                    // Format: "HP:{hpo_id}"
  labels: ['Phenotype', 'CrossDomain'];

  properties: {
    hpo_id: string;              // HP:0000000
    name: string;
    definition: string;
    synonyms: string[];

    // Hierarchy
    parent_ids: string[];

    // Frequency
    frequency_category: string;  // HP:0040280 (frequency terms)
  };
}
```

---

## Part 3: Edge Type Definitions

### 3.1 Binary Edges (CodeGraph.createEdge)

#### Gene-Variant Relationships
```typescript
// Gene contains variant
{
  from: "GENE:7157",           // TP53
  to: "VAR:rs1042522",
  type: "HAS_VARIANT",
  properties: {
    consequence: "missense_variant",
    impact: "MODERATE",
    exon: 4,
    codon_change: "Arg>Pro"
  }
}
```

#### Variant-Disease Relationships
```typescript
{
  from: "VAR:rs1042522",
  to: "DIS:MONDO:0005072",     // Cancer
  type: "ASSOCIATED_WITH",
  properties: {
    clinvar_significance: "pathogenic",
    review_status: "reviewed by expert panel",
    inheritance: "AD",
    penetrance: "high",
    evidence_pmids: [12345678, 23456789]
  }
}
```

#### Compound-Target Relationships
```typescript
{
  from: "COMP:BQJCRHHNABKAKU-UHFFFAOYSA-N",  // Quercetin
  to: "GENE:5594",                            // MAPK1
  type: "INHIBITS",
  properties: {
    activity_type: "IC50",
    activity_value: 5.2,
    activity_unit: "uM",
    assay_type: "binding",
    source: "ChEMBL",
    confidence: 0.95
  }
}
```

#### Organism-Compound Relationships
```typescript
{
  from: "ORG:3760",            // Prunus persica (peach)
  to: "COMP:REFJWTPEDVJJIY-UHFFFAOYSA-N",  // Amygdalin
  type: "PRODUCES",
  properties: {
    plant_part: "seed",
    concentration_mg_g: 2.5,
    source: "LOTUS"
  }
}
```

#### Gene-Pathway Relationships
```typescript
{
  from: "GENE:7157",           // TP53
  to: "PATH:R-HSA-69473",      // G2/M DNA damage checkpoint
  type: "PARTICIPATES_IN",
  properties: {
    role: "regulator",
    evidence: "experimental"
  }
}
```

#### Disease-Phenotype Relationships
```typescript
{
  from: "DIS:MONDO:0007256",   // Breast cancer
  to: "HP:0100013",            // Neoplasm of breast
  type: "HAS_PHENOTYPE",
  properties: {
    frequency: "HP:0040281",   // Very frequent (99-80%)
    onset: "adult"
  }
}
```

#### Drug-Gene (Pharmacogenomics) Relationships
```typescript
{
  from: "COMP:CYQFCXCEBYINGO-UHFFFAOYSA-N",  // Warfarin
  to: "GENE:1565",                            // CYP2D6
  type: "METABOLIZED_BY",
  properties: {
    pharmgkb_level: "1A",
    cpic_guideline: true,
    dosing_recommendation: "Reduce dose in poor metabolizers"
  }
}
```

### 3.2 Full Edge Type Catalog

| Edge Type | From Node | To Node | Properties |
|-----------|-----------|---------|------------|
| `HAS_VARIANT` | Gene | Variant | consequence, impact, exon |
| `ASSOCIATED_WITH` | Variant | Disease | significance, inheritance, penetrance |
| `INHIBITS` | Compound | Gene/Protein | activity_type, activity_value, source |
| `ACTIVATES` | Compound | Gene/Protein | activity_type, activity_value, source |
| `BINDS_TO` | Compound | Gene/Protein | affinity, assay_type |
| `PRODUCES` | Organism | Compound | plant_part, concentration |
| `PARTICIPATES_IN` | Gene | Pathway | role, evidence |
| `HAS_PHENOTYPE` | Disease | Phenotype | frequency, onset |
| `METABOLIZED_BY` | Compound | Gene | pharmgkb_level, guideline |
| `TREATS` | Compound | Disease | evidence_level, mechanism |
| `CAUSES` | Compound | Disease | adverse_effect, frequency |
| `CONTAINS_NUTRIENT` | Food | Nutrient | amount, unit |
| `CONTAINS_COMPOUND` | Food | Compound | concentration |
| `IS_A` | Disease | Disease | ontology_source |
| `PART_OF` | Pathway | Pathway | hierarchy_level |
| `INTERACTS_WITH` | Gene | Gene | interaction_type, confidence |
| `CO_EXPRESSED` | Gene | Gene | correlation, tissue |
| `HAPLOTYPE_MEMBER` | Variant | Haplotype | allele_contribution |

---

## Part 4: Hyperedge Patterns (N-ary Relationships)

### 4.1 Why Hyperedges?

Traditional binary edges cannot capture relationships involving 3+ entities simultaneously. RuVector's `CodeGraph.createHyperedge()` enables native N-ary relationships critical for:

1. **Formula composition** - Formula + multiple ingredients + dosages
2. **Pathway reactions** - Substrates + enzymes + products + cofactors
3. **Pharmacogenomic interactions** - Drug + gene + variant + phenotype
4. **Food-nutrient matrices** - Food + multiple nutrients + amounts
5. **Polypharmacy interactions** - Drug1 + Drug2 + Gene + Effect

### 4.2 Hyperedge Type Definitions

#### Formula Composition Hyperedge
```typescript
// Traditional medicine formula with multiple ingredients
{
  id: "HE:FORM_COMP:liuwei_dihuang",
  nodes: [
    "FORM:KampoDB:K001",              // Liu Wei Di Huang Wan
    "COMP:LYGYBDGRHCYKKM-UHFFFAOYSA-N", // Ingredient 1
    "COMP:QHTVXCFVNLQJJT-UHFFFAOYSA-N", // Ingredient 2
    "COMP:XZWLCPWMRBFABZ-UHFFFAOYSA-N", // Ingredient 3
    // ... more ingredients
  ],
  type: "FORMULA_COMPOSITION",
  properties: {
    roles: {
      "FORM:KampoDB:K001": "formula",
      "COMP:LYGYBDGRHCYKKM-UHFFFAOYSA-N": "emperor",  // Jun
      "COMP:QHTVXCFVNLQJJT-UHFFFAOYSA-N": "minister", // Chen
      "COMP:XZWLCPWMRBFABZ-UHFFFAOYSA-N": "assistant" // Zuo
    },
    dosages: {
      "COMP:LYGYBDGRHCYKKM-UHFFFAOYSA-N": { amount: 24, unit: "g" },
      "COMP:QHTVXCFVNLQJJT-UHFFFAOYSA-N": { amount: 12, unit: "g" }
    },
    synergy_score: 0.85,
    traditional_indication: "Kidney Yin deficiency"
  }
}
```

#### Pathway Reaction Hyperedge
```typescript
// Metabolic reaction with multiple participants
{
  id: "HE:RXN:R-HSA-71493",
  nodes: [
    "COMP:WQZGKKKJIJFFOK-GASJEMHNSA-N", // Glucose (substrate)
    "COMP:BAWFJGJZGIEFAR-NNYOXOHSSA-O", // ATP (cofactor)
    "GENE:3098",                         // HK1 (enzyme)
    "COMP:NBSCHQHZLSJFNQ-GASJEMHNSA-L", // G6P (product)
    "COMP:XTWYTFMLZFPYCI-KQYNXXCUSA-K", // ADP (byproduct)
    "PATH:R-HSA-70171"                   // Glycolysis pathway
  ],
  type: "PATHWAY_REACTION",
  properties: {
    roles: {
      "COMP:WQZGKKKJIJFFOK-GASJEMHNSA-N": "substrate",
      "COMP:BAWFJGJZGIEFAR-NNYOXOHSSA-O": "cofactor",
      "GENE:3098": "catalyst",
      "COMP:NBSCHQHZLSJFNQ-GASJEMHNSA-L": "product",
      "COMP:XTWYTFMLZFPYCI-KQYNXXCUSA-K": "byproduct",
      "PATH:R-HSA-70171": "context"
    },
    reaction_type: "phosphorylation",
    reversible: false,
    stoichiometry: { substrate: 1, product: 1, cofactor: 1 }
  }
}
```

#### Pharmacogenomic Interaction Hyperedge
```typescript
// Drug-gene-variant-phenotype quadruple
{
  id: "HE:PGX:warfarin-cyp2c9-star3",
  nodes: [
    "COMP:CYQFCXCEBYINGO-UHFFFAOYSA-N", // Warfarin
    "GENE:1559",                         // CYP2C9
    "VAR:rs1057910",                     // *3 allele
    "HP:0001892"                         // Abnormal bleeding
  ],
  type: "PHARMACOGENOMIC_INTERACTION",
  properties: {
    roles: {
      "COMP:CYQFCXCEBYINGO-UHFFFAOYSA-N": "drug",
      "GENE:1559": "metabolizing_enzyme",
      "VAR:rs1057910": "variant",
      "HP:0001892": "phenotype"
    },
    effect: "increased_bleeding_risk",
    mechanism: "reduced_metabolism",
    clinical_action: "reduce_dose_30_50_percent",
    evidence_level: "1A",
    cpic_guideline: "CPIC-2017-warfarin"
  }
}
```

#### Food Nutrient Profile Hyperedge
```typescript
// Food with full nutrient composition
{
  id: "HE:FOOD_NUT:usda-09003",
  nodes: [
    "FOOD:USDA:09003",    // Apple, raw, with skin
    "NUTR:1003",          // Protein
    "NUTR:1004",          // Total fat
    "NUTR:1005",          // Carbohydrate
    "NUTR:1008",          // Energy
    "NUTR:1087",          // Calcium
    "NUTR:1089",          // Iron
    // ... more nutrients
  ],
  type: "NUTRIENT_PROFILE",
  properties: {
    serving_size_g: 100,
    amounts: {
      "NUTR:1003": 0.26,   // g protein
      "NUTR:1004": 0.17,   // g fat
      "NUTR:1005": 13.81,  // g carbs
      "NUTR:1008": 52,     // kcal
      "NUTR:1087": 6,      // mg calcium
      "NUTR:1089": 0.12    // mg iron
    },
    data_quality: "measured"
  }
}
```

#### Drug-Drug-Gene Interaction Hyperedge
```typescript
// Polypharmacy interaction
{
  id: "HE:DDI:omeprazole-clopidogrel-cyp2c19",
  nodes: [
    "COMP:SUBDBMMJDZJVKF-UHFFFAOYSA-N", // Omeprazole
    "COMP:GKTWGGQPFAXNFI-HNNXBMFYSA-N", // Clopidogrel
    "GENE:1557",                         // CYP2C19
    "HP:0001658"                         // Myocardial infarction risk
  ],
  type: "DRUG_DRUG_GENE_INTERACTION",
  properties: {
    roles: {
      "COMP:SUBDBMMJDZJVKF-UHFFFAOYSA-N": "perpetrator",
      "COMP:GKTWGGQPFAXNFI-HNNXBMFYSA-N": "victim",
      "GENE:1557": "shared_enzyme",
      "HP:0001658": "clinical_outcome"
    },
    interaction_type: "competitive_inhibition",
    clinical_significance: "major",
    recommendation: "avoid_combination_or_use_alternative_ppi"
  }
}
```

### 4.3 Hyperedge Querying Patterns

```typescript
// Find all hyperedges containing a specific compound
const hyperedges = codeGraph.getHyperedges(
  "COMP:CYQFCXCEBYINGO-UHFFFAOYSA-N",
  "PHARMACOGENOMIC_INTERACTION"
);

// Find formulas containing a specific ingredient
const formulas = codeGraph.getHyperedges(
  "COMP:QUERCETIN_INCHIKEY",
  "FORMULA_COMPOSITION"
).filter(he => he.properties.roles[compoundId] !== "formula");

// Cypher query for pathway reactions involving a gene
const reactions = codeGraph.cypher(`
  MATCH (n:Gene {ncbi_gene_id: 3098})-[:MEMBER_OF]->(he:Hyperedge)
  WHERE he.type = 'PATHWAY_REACTION' AND he.properties.roles[$nodeId] = 'catalyst'
  RETURN he
`, { nodeId: "GENE:3098" });
```

---

## Part 5: Vector Embedding Integration Points

### 5.1 Embedding Architecture Overview

```
+------------------+     +------------------+     +------------------+
|  Text Embeddings |     |Structure Embeddings|   | Profile Vectors  |
|    (384-dim)     |     |   (2048-bit FP)    |   |   (50-100 dim)   |
+------------------+     +------------------+     +------------------+
        |                        |                        |
        v                        v                        v
  Gene descriptions       Compound structures      Nutrient profiles
  Disease names           Morgan fingerprints      Population frequencies
  Pathway summaries       Pharmacophores           Activity profiles
  Phenotype terms
```

### 5.2 Embedding Storage Strategy

| Entity Type | VectorDB Key Format | Dimensions | Content Source | Use Case |
|-------------|---------------------|------------|----------------|----------|
| Gene semantic | `GENE_SEM:{ncbi_id}` | 384 | Description + GO + pathways | Functional similarity |
| Variant population | `VAR_POP:{rsid}` | 8 | Population frequencies | Population similarity |
| Compound structure | `COMP_FP:{inchikey}` | 2048 | Morgan fingerprint | Structure similarity |
| Compound semantic | `COMP_SEM:{inchikey}` | 384 | Name + class + activities | Bioactivity similarity |
| Food nutrient | `FOOD_NUT:{id}` | 50 | Nutrient values normalized | Nutrient similarity |
| Food phytochemical | `FOOD_PHYTO:{id}` | 300 | Compound presence | Phytochemical similarity |
| Disease semantic | `DIS_SEM:{mondo_id}` | 384 | Name + phenotypes | Disease similarity |
| Pathway semantic | `PATH_SEM:{id}` | 384 | Description + genes | Pathway similarity |
| Formula semantic | `FORM_SEM:{id}` | 384 | Indications + composition | Formula similarity |

### 5.3 Multi-Vector Queries

```typescript
// Find compounds structurally similar AND semantically similar to target
async function findSimilarCompounds(
  targetInchikey: string,
  structureWeight: number = 0.6,
  semanticWeight: number = 0.4,
  topK: number = 20
) {
  // Get both embeddings for target
  const structVec = await vectorDb.get(`COMP_FP:${targetInchikey}`);
  const semVec = await vectorDb.get(`COMP_SEM:${targetInchikey}`);

  // Search both spaces
  const structResults = await vectorDb.search({
    vector: structVec.vector,
    k: topK * 2,
    filter: { namespace: 'COMP_FP' }
  });

  const semResults = await vectorDb.search({
    vector: semVec.vector,
    k: topK * 2,
    filter: { namespace: 'COMP_SEM' }
  });

  // Fuse results with weighted scoring
  const fused = fuseResults(structResults, semResults, structureWeight, semanticWeight);
  return fused.slice(0, topK);
}
```

### 5.4 Embedding Generation Pipeline

```typescript
import { OnnxEmbedder, AdaptiveEmbedder } from 'ruvector';

// Initialize ONNX embedder (MiniLM-L6 for semantic)
const semanticEmbedder = new OnnxEmbedder({
  model: 'all-MiniLM-L6-v2',
  dimensions: 384
});

// Gene embedding generation
async function embedGene(gene: GeneNode): Promise<Float32Array> {
  const text = [
    gene.properties.hgnc_symbol,
    gene.properties.name,
    gene.properties.description,
    // Add GO terms
    ...(await getGOTerms(gene.id)).map(t => t.name),
    // Add pathway names
    ...(await getPathways(gene.id)).map(p => p.name)
  ].join(' ');

  return semanticEmbedder.embed(text);
}

// Compound structure embedding (Morgan fingerprint)
function embedCompoundStructure(smiles: string): Float32Array {
  // Using RDKit-js or mol2vec
  return computeMorganFingerprint(smiles, radius=2, nBits=2048);
}

// Population frequency vector (for variant similarity by population)
function embedVariantPopulation(variant: VariantNode): Float32Array {
  return new Float32Array([
    variant.properties.gnomad_populations.afr || 0,
    variant.properties.gnomad_populations.amr || 0,
    variant.properties.gnomad_populations.asj || 0,
    variant.properties.gnomad_populations.eas || 0,
    variant.properties.gnomad_populations.fin || 0,
    variant.properties.gnomad_populations.mid || 0,
    variant.properties.gnomad_populations.nfe || 0,
    variant.properties.gnomad_populations.sas || 0
  ]);
}
```

### 5.5 SONA Integration for Learning

```typescript
import { Sona } from 'ruvector';

// Initialize SONA for adaptive learning from queries
const sona = new Sona.Engine(384);

// Track query patterns for personalization
async function adaptiveQuery(queryText: string, userId: string) {
  const embedding = await semanticEmbedder.embed(queryText);

  // Begin trajectory for this query
  const trajId = sona.beginTrajectory(embedding);

  // Execute search
  const results = await vectorDb.search({
    vector: embedding,
    k: 20
  });

  // Apply micro-LoRA adaptation based on past success patterns
  const adaptedEmbedding = sona.applyMicroLora(embedding);

  // Re-rank with adapted embedding
  const rerankedResults = await vectorDb.search({
    vector: adaptedEmbedding,
    k: 20
  });

  // Record step with attention weights
  sona.addStep(trajId, embedding, new Float32Array(384), 0.5);

  return rerankedResults;
}

// When user clicks/selects a result, record feedback
function recordFeedback(trajId: number, selectedId: string, quality: number) {
  sona.endTrajectory(trajId, quality);
}
```

---

## Part 6: Sharding Strategy

### 6.1 Domain-Based Sharding

The THREE WORLDS data naturally partitions into domains with different access patterns and sizes:

```
+------------------------------------------------------------------+
|                      RuvectorCluster                              |
+------------------------------------------------------------------+
|                                                                   |
|  Shard 0-3: GENETICS_HOT         Shard 4-7: GENETICS_COLD        |
|  +-----------------------+       +-----------------------+        |
|  | ~10M variants         |       | ~1.19B variants       |        |
|  | Clinically relevant   |       | Rare/unknown signif.  |        |
|  | High query volume     |       | Low query volume      |        |
|  | Full precision        |       | Quantized vectors     |        |
|  +-----------------------+       +-----------------------+        |
|                                                                   |
|  Shard 8-9: GENES                Shard 10-11: COMPOUNDS           |
|  +-----------------------+       +-----------------------+        |
|  | ~25K genes            |       | ~800K compounds       |        |
|  | Hot cache             |       | Structure + semantic  |        |
|  | Full graph topology   |       | HNSW-indexed         |        |
|  +-----------------------+       +-----------------------+        |
|                                                                   |
|  Shard 12-13: FOODS              Shard 14-15: ONTOLOGIES          |
|  +-----------------------+       +-----------------------+        |
|  | ~380K foods           |       | Diseases, phenotypes  |        |
|  | Nutrient vectors      |       | Pathways              |        |
|  | Phytochemical maps    |       | Graph hierarchy       |        |
|  +-----------------------+       +-----------------------+        |
|                                                                   |
+------------------------------------------------------------------+
```

### 6.2 Shard Configuration

```typescript
import { RuvectorCluster, createCluster } from 'ruvector';

const cluster = createCluster({
  nodeId: 'node-1',
  address: 'localhost:8080',
  peers: ['localhost:8081', 'localhost:8082'],
  shards: 16,
  replicationFactor: 2
});

// Shard routing configuration
const SHARD_CONFIG = {
  // Genetics - Hot (clinical variants)
  'GENETICS_HOT': {
    shards: [0, 1, 2, 3],
    keyPattern: /^VAR:(rs\d+|VCV\d+)/,
    filter: (v: VariantNode) =>
      v.properties.clinvar_sig !== null ||
      v.properties.gnomad_af > 0.01
  },

  // Genetics - Cold (rare/VUS)
  'GENETICS_COLD': {
    shards: [4, 5, 6, 7],
    keyPattern: /^VAR:/,
    filter: (v: VariantNode) =>
      v.properties.clinvar_sig === null &&
      v.properties.gnomad_af <= 0.01
  },

  // Genes
  'GENES': {
    shards: [8, 9],
    keyPattern: /^GENE:/
  },

  // Compounds
  'COMPOUNDS': {
    shards: [10, 11],
    keyPattern: /^COMP:/
  },

  // Foods
  'FOODS': {
    shards: [12, 13],
    keyPattern: /^FOOD:/
  },

  // Ontologies (diseases, phenotypes, pathways)
  'ONTOLOGIES': {
    shards: [14, 15],
    keyPattern: /^(DIS|HP|PATH|FORM):/
  }
};

// Route key to shard
function routeToShard(key: string): number {
  for (const [domain, config] of Object.entries(SHARD_CONFIG)) {
    if (config.keyPattern.test(key)) {
      // Hash within domain's shard range
      const hash = cyrb53(key);
      const shardIndex = hash % config.shards.length;
      return config.shards[shardIndex];
    }
  }
  // Default: hash across all shards
  return cyrb53(key) % 16;
}
```

### 6.3 Hot/Cold Data Strategy

```typescript
// Variant access patterns determine hot/cold placement
interface VariantAccessProfile {
  queryFrequency: number;      // Queries per day
  lastAccess: number;          // Timestamp
  clinicalRelevance: number;   // 0-1 score
}

class HotColdManager {
  private accessProfiles: Map<string, VariantAccessProfile> = new Map();

  // Track access for promotion/demotion
  recordAccess(variantId: string, isClinical: boolean) {
    const profile = this.accessProfiles.get(variantId) || {
      queryFrequency: 0,
      lastAccess: 0,
      clinicalRelevance: isClinical ? 1 : 0
    };

    profile.queryFrequency++;
    profile.lastAccess = Date.now();
    this.accessProfiles.set(variantId, profile);

    // Check for promotion to hot shard
    if (this.shouldPromote(profile)) {
      this.promoteToHot(variantId);
    }
  }

  shouldPromote(profile: VariantAccessProfile): boolean {
    return profile.queryFrequency > 100 || profile.clinicalRelevance > 0.8;
  }

  shouldDemote(profile: VariantAccessProfile): boolean {
    const daysSinceAccess = (Date.now() - profile.lastAccess) / (1000 * 60 * 60 * 24);
    return daysSinceAccess > 90 && profile.clinicalRelevance < 0.5;
  }

  async promoteToHot(variantId: string) {
    const coldShard = cluster.getShardForKey(variantId);
    if (SHARD_CONFIG.GENETICS_COLD.shards.includes(coldShard.id)) {
      // Move to hot shard with full precision
      const data = await cluster.get(variantId);
      const hotShardId = SHARD_CONFIG.GENETICS_HOT.shards[
        cyrb53(variantId) % SHARD_CONFIG.GENETICS_HOT.shards.length
      ];
      await cluster.put(variantId, data); // Will route to hot
    }
  }
}
```

### 6.4 Vector Quantization for Cold Storage

```typescript
import { TensorCompress } from 'ruvector';

const compressor = new TensorCompress();

// Store cold variants with compressed vectors
async function storeVariantCold(variant: VariantNode, embedding: Float32Array) {
  // Quantize embedding for cold storage (8-bit)
  const accessFreq = 0.1; // Low access frequency
  const compressed = compressor.compress(embedding, accessFreq);

  // Store in cold shard
  await cluster.put(`VAR:${variant.properties.rsid}`, {
    node: variant,
    embedding_compressed: compressed,
    compression_level: 'int8'
  });
}

// Decompress on access
async function getVariantEmbedding(variantId: string): Promise<Float32Array> {
  const stored = await cluster.get(variantId);

  if (stored.embedding_compressed) {
    return new Float32Array(compressor.decompress(stored.embedding_compressed));
  }
  return stored.embedding;
}
```

### 6.5 Cross-Domain Query Patterns

```typescript
// Query spanning multiple domains (shards)
async function findDrugTargetsForDisease(diseaseId: string) {
  // 1. Get disease from ONTOLOGIES shard
  const disease = await cluster.get(`DIS:${diseaseId}`);

  // 2. Find associated genes from GENES shard
  const geneEdges = await codeGraph.getIncomingEdges(diseaseId, 'ASSOCIATED_WITH');
  const geneIds = geneEdges.map(e => e.from);

  // 3. Find compounds targeting these genes from COMPOUNDS shard
  const compounds = [];
  for (const geneId of geneIds) {
    const targetEdges = await codeGraph.getIncomingEdges(geneId, 'INHIBITS');
    compounds.push(...targetEdges.map(e => e.from));
  }

  // 4. Fetch compound details
  const compoundDetails = await Promise.all(
    compounds.map(id => cluster.get(id))
  );

  return compoundDetails;
}

// Parallel cross-shard query
async function parallelCrossShardQuery(geneIds: string[]) {
  // Group by shard
  const byShardMap = new Map<number, string[]>();
  for (const geneId of geneIds) {
    const shardId = routeToShard(geneId);
    if (!byShardMap.has(shardId)) {
      byShardMap.set(shardId, []);
    }
    byShardMap.get(shardId)!.push(geneId);
  }

  // Query each shard in parallel
  const results = await Promise.all(
    Array.from(byShardMap.entries()).map(async ([shardId, ids]) => {
      return Promise.all(ids.map(id => cluster.get(id)));
    })
  );

  return results.flat();
}
```

---

## Part 7: Implementation Roadmap

### 7.1 Phase 1: Core Infrastructure (Weeks 1-4)

| Task | RuVector Component | Deliverable |
|------|-------------------|-------------|
| Initialize CodeGraph | `CodeGraph` | Graph database with hyperedge support |
| Setup VectorDB instances | `VectorDB` | Per-domain vector stores |
| Configure cluster | `RuvectorCluster` | 16-shard distributed setup |
| Implement shard routing | Custom | Key-to-shard routing logic |

### 7.2 Phase 2: Node Types (Weeks 5-8)

| Task | Entities | Storage |
|------|----------|---------|
| Gene nodes | 25K genes | Graph + Vector (semantic) |
| Variant nodes (hot) | 10M variants | Graph + Vector (population) |
| Compound nodes | 800K compounds | Graph + Vector (structure + semantic) |
| Disease nodes | 80K diseases | Graph + Vector (semantic) |

### 7.3 Phase 3: Relationships (Weeks 9-12)

| Task | Edge Types | Count Estimate |
|------|------------|----------------|
| Gene-Variant edges | `HAS_VARIANT` | 1M+ |
| Variant-Disease edges | `ASSOCIATED_WITH` | 2M+ |
| Compound-Target edges | `INHIBITS`, `BINDS_TO` | 20M+ |
| Pathway membership | `PARTICIPATES_IN` | 50K+ |

### 7.4 Phase 4: Hyperedges (Weeks 13-16)

| Task | Hyperedge Type | Count Estimate |
|------|----------------|----------------|
| Formula compositions | `FORMULA_COMPOSITION` | 60K |
| Pathway reactions | `PATHWAY_REACTION` | 100K+ |
| Pharmacogenomic | `PHARMACOGENOMIC_INTERACTION` | 10K+ |
| Food nutrients | `NUTRIENT_PROFILE` | 380K |

### 7.5 Phase 5: Cold Storage & Optimization (Weeks 17-20)

| Task | Description |
|------|-------------|
| Variant cold migration | Move 1.19B rare variants to cold shards |
| Vector quantization | Compress cold vectors to int8 |
| Query optimization | Implement cross-shard parallel queries |
| SONA integration | Enable adaptive learning from query patterns |

---

## Part 8: Query Patterns & Examples

### 8.1 Semantic Similarity Queries

```typescript
// Find genes similar to TP53 by function
const similarGenes = await vectorDb.search({
  vector: await getEmbedding('GENE_SEM:7157'),
  k: 20,
  filter: { labels: ['Gene'] }
});

// Find compounds similar to quercetin by structure
const similarCompounds = await vectorDb.search({
  vector: await getEmbedding('COMP_FP:REFJWTPEDVJJIY-UHFFFAOYSA-N'),
  k: 50,
  filter: { labels: ['Compound'] }
});

// Find foods with similar nutrient profile to spinach
const similarFoods = await vectorDb.search({
  vector: await getEmbedding('FOOD_NUT:USDA:11457'),
  k: 30,
  filter: { labels: ['Food'] }
});
```

### 8.2 Graph Traversal Queries

```typescript
// Find all pathways affected by a drug
const drugPathways = codeGraph.cypher(`
  MATCH (c:Compound {inchikey: $inchikey})-[:INHIBITS|:ACTIVATES]->(g:Gene)
        -[:PARTICIPATES_IN]->(p:Pathway)
  RETURN DISTINCT p.name, p.reactome_id, count(g) as affected_genes
  ORDER BY affected_genes DESC
`, { inchikey: 'CYQFCXCEBYINGO-UHFFFAOYSA-N' });

// Find pharmacogenomic interactions for a patient's variants
const pgxInteractions = codeGraph.cypher(`
  MATCH (v:Variant)-[:MEMBER_OF]->(he:Hyperedge)
  WHERE v.rsid IN $rsids AND he.type = 'PHARMACOGENOMIC_INTERACTION'
  RETURN he.properties as interaction
`, { rsids: ['rs1042522', 'rs1057910', 'rs4244285'] });
```

### 8.3 Hybrid Queries (Vector + Graph)

```typescript
// Find natural products that might treat a disease
async function findNaturalProductsForDisease(diseaseId: string) {
  // 1. Get disease embedding
  const diseaseEmb = await vectorDb.get(`DIS_SEM:${diseaseId}`);

  // 2. Find semantically similar compounds (by bioactivity description)
  const similarCompounds = await vectorDb.search({
    vector: diseaseEmb.vector,
    k: 100,
    filter: { labels: ['Compound'], is_natural_product: true }
  });

  // 3. For each compound, trace to source organisms
  const results = [];
  for (const comp of similarCompounds) {
    const sources = codeGraph.getIncomingEdges(comp.id, 'PRODUCES');
    const targets = codeGraph.getOutgoingEdges(comp.id, 'INHIBITS');

    results.push({
      compound: comp,
      sources: sources.map(s => s.from),
      targets: targets.map(t => t.to),
      traditional_use: await getTraditionalUse(comp.id)
    });
  }

  return results;
}
```

### 8.4 Hyperedge Queries

```typescript
// Find all formulas containing a specific ingredient with its role
const formulasWithIngredient = codeGraph.cypher(`
  MATCH (c:Compound {inchikey: $inchikey})-[:MEMBER_OF]->(he:Hyperedge)
  WHERE he.type = 'FORMULA_COMPOSITION'
  RETURN he.nodes, he.properties.roles[$compId] as role,
         he.properties.dosages[$compId] as dosage
`, { inchikey: 'QUERCETIN_INCHIKEY', compId: 'COMP:QUERCETIN_INCHIKEY' });

// Find all reactions where a gene is the catalyst
const enzymeReactions = codeGraph.cypher(`
  MATCH (g:Gene {ncbi_gene_id: $geneId})-[:MEMBER_OF]->(he:Hyperedge)
  WHERE he.type = 'PATHWAY_REACTION' AND he.properties.roles[$nodeId] = 'catalyst'
  RETURN he.nodes, he.properties.reaction_type
`, { geneId: 3098, nodeId: 'GENE:3098' });
```

---

## Part 9: Performance Considerations

### 9.1 Expected Query Performance

| Query Type | Expected Latency | Notes |
|------------|------------------|-------|
| Single node lookup | <5ms | Direct key lookup via cluster |
| Vector similarity (k=20) | <50ms | HNSW search with ef_search=100 |
| Graph traversal (2 hops) | <100ms | Cached adjacency |
| Hyperedge lookup | <20ms | Direct hyperedge retrieval |
| Cross-shard query (10 shards) | <500ms | Parallel execution |
| Full-text + vector hybrid | <200ms | Combined index access |

### 9.2 Memory Footprint Estimates

| Component | Memory Estimate | Notes |
|-----------|-----------------|-------|
| Gene embeddings (25K x 384) | ~37 MB | Hot, full precision |
| Hot variant embeddings (10M x 8) | ~80 MB | Population vectors only |
| Compound structure FP (800K x 2048) | ~1.6 GB | Binary fingerprints |
| Compound semantic (800K x 384) | ~1.2 GB | Dense vectors |
| Food nutrient (380K x 50) | ~76 MB | Normalized profiles |
| HNSW indices | ~2x vector size | Graph overhead |
| Graph topology | ~500 MB | Edges + hyperedges |
| **Total Hot Tier** | **~6-8 GB** | Per node |

### 9.3 Storage Estimates

| Data Tier | Storage Estimate | Compression |
|-----------|------------------|-------------|
| Hot vectors | ~20 GB | None (float32) |
| Cold vectors (1.2B variants) | ~150 GB | int8 quantized |
| Graph storage | ~50 GB | LZ4 compressed |
| Hyperedge storage | ~20 GB | JSON + compression |
| Indices | ~100 GB | HNSW + B-tree |
| **Total per node** | **~340 GB** | With replication |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [schemas-index.md](./schemas-index.md) | Parent navigation |
| [unified-schema-analysis.md](./unified-schema-analysis.md) | Entity definitions source |
| [world1-schema-research.md](../research/world1-schema-research.md) | Genetics data sources |

---

## Appendix A: RuVector API Quick Reference

```typescript
// VectorDB
const db = new VectorDB({ dimensions: 384, storagePath: './data' });
await db.insert({ id: 'key', vector: [...], metadata: {...} });
const results = await db.search({ vector: [...], k: 20, filter: {...} });

// CodeGraph
const graph = new CodeGraph({ storagePath: './graph' });
graph.createNode('id', ['Label1', 'Label2'], { prop: 'value' });
graph.createEdge('from', 'to', 'TYPE', { prop: 'value' });
graph.createHyperedge(['n1', 'n2', 'n3'], 'TYPE', { roles: {...} });
const path = graph.shortestPath('from', 'to', maxDepth);
const result = graph.cypher('MATCH (n) RETURN n');

// RuvectorCluster
const cluster = createCluster({ nodeId: 'n1', shards: 16 });
await cluster.put('key', value);
const data = await cluster.get('key');
const shard = cluster.getShardForKey('key');

// SONA
const sona = new Sona.Engine(384);
const trajId = sona.beginTrajectory(embedding);
sona.addStep(trajId, activations, attention, reward);
sona.endTrajectory(trajId, quality);
const adapted = sona.applyMicroLora(input);
```

---

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | Variable |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | JSON |
| Alternative | Hyperbolic embeddings |
| Encoding | UTF-8 |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| RuVector | Repository | See main database |
| Three Worlds Model | Research | Research publication |

**Access Requirements:** Open access, research use

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| RuVector | Open Source | Yes (see repository) |
| Three Worlds | Research | See publication terms |

---

## Sample Data

### Example Record
```json
{
  "id": "GENE:7157",
  "labels": ["Gene", "WORLD1"],
  "ncbi_gene_id": 7157,
  "hgnc_symbol": "TP53",
  "chromosome": "17",
  "start_pos": 7571720
}
```

### Sample Query Result
| id | labels | hgnc_symbol | ncbi_gene_id | chromosome |
|----|--------|-------------|------|-----------|
| GENE:7157 | Gene, WORLD1 | TP53 | 7157 | 17 |
| GENE:3156 | Gene, WORLD1 | HMGCR | 3156 | 5 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `node` | Graph database entity representing a concept | GeneNode, CompoundNode |
| `edge` | Binary relationship between two nodes | HAS_VARIANT, INHIBITS |
| `hyperedge` | N-ary relationship connecting 3+ nodes simultaneously | FORMULA_COMPOSITION |
| `shard` | Partition of distributed database for scalability | GENETICS_HOT shard |
| `embedding` | Vector representation of an entity for similarity search | 384-dim gene embedding |
| `hot_tier` | Frequently accessed data stored at full precision | Clinical variants |
| `cold_tier` | Rarely accessed data stored with compression | Rare variants |
| `quantization` | Compression technique reducing vector precision | int8 quantized vectors |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| VectorDB | RuVector component for HNSW-indexed vector storage | Similarity search |
| CodeGraph | RuVector component for graph with hyperedge support | Relationships |
| RuvectorCluster | RuVector component for distributed sharding | Scalability |
| SONA | Self-Optimizing Neural Architecture for adaptive learning | Query personalization |
| HNSW | Hierarchical Navigable Small World graph for fast similarity | Vector index |
| NeuralSubstrate | RuVector memory physics for cache management | Semantic drift |
| GNN_Wrapper | Graph Neural Network integration for differentiable search | ML queries |
| Morgan_fingerprint | Circular molecular fingerprint for structure similarity | 2048-bit vectors |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HNSW | Hierarchical Navigable Small World | Vector index |
| SONA | Self-Optimizing Neural Architecture | RuVector learning |
| GNN | Graph Neural Network | ML on graphs |
| FP | Fingerprint | Molecular structure |
| LoRA | Low-Rank Adaptation | Efficient fine-tuning |
| DDI | Drug-Drug Interaction | Polypharmacy |
| PGx | Pharmacogenomics | Drug-gene variants |
| SPDI | Sequence-Position-Deletion-Insertion | Variant notation |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 21, 2026 | System Architecture | Initial schema design |
