# Unified Schema Analysis

**Document ID:** UNIFIED-SCHEMA-ANALYSIS
**Status:** Final
**Owner:** Architecture Team
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../index.md](../index.md)

---

## Executive Summary

This document synthesizes schema patterns across 250+ databases documented in the data sources collection to propose a unified data model for the knowledge base. The analysis identifies five core entity types (Gene, Variant, Compound, Disease, Pathway), maps 15+ identifier systems, analyzes relationship patterns, and proposes a normalized schema with provenance tracking.

### Key Findings

1. **Hub Identifiers**: UniProt (proteins), NCBI Gene ID (genes), InChIKey (compounds), MONDO (diseases), and Reactome (pathways) serve as optimal canonical identifiers
2. **Common Patterns**: All databases follow a similar entity-relationship model with Gene -> Variant -> Disease and Compound -> Target -> Pathway chains
3. **Evidence Complexity**: Evidence quality varies dramatically (experimental vs. predicted vs. text-mined) requiring robust provenance tracking
4. **Integration Challenges**: Identifier fragmentation, naming inconsistencies, and schema heterogeneity require careful normalization

---

## Part 1: Common Entity Type Analysis

### 1.1 Gene Entity Comparison

| Database | Gene ID Type | Gene Count | Key Attributes | Strengths |
|----------|--------------|------------|----------------|-----------|
| **dbSNP** | NCBI Gene ID | ~30K human | Symbol, location, variants | Variant linkage |
| **PharmGKB** | PA ID + HGNC | 1,761 | Haplotypes, PGx annotations | Clinical annotations |
| **Reactome** | UniProt ID | 11,020 | Pathway membership, reactions | Pathway context |
| **DisGeNET** | NCBI Gene ID | 21,671 | Gene-disease associations | Disease linkage |
| **Gene Ontology** | UniProt/NCBI | 20K+ human | Functions, processes, locations | Functional annotations |

**Proposed Gene Schema:**

```yaml
Gene:
  identifiers:
    internal_id: string          # UUID for internal reference
    ncbi_gene_id: integer        # Primary: NCBI Gene (Entrez)
    hgnc_id: string              # HGNC:12345
    hgnc_symbol: string          # Official symbol
    ensembl_gene_id: string      # ENSG00000000000
    uniprot_ids: array[string]   # Multiple protein products
    pharmgkb_id: string          # PA123

  core_attributes:
    symbol: string               # Primary gene symbol
    name: string                 # Full gene name
    aliases: array[string]       # Alternative symbols
    description: string          # Function summary

  genomic_location:
    chromosome: string           # 1-22, X, Y, MT
    start_position: integer      # GRCh38
    end_position: integer
    strand: enum[+, -]
    cytogenetic_band: string     # e.g., "17q21.31"

  functional_annotations:
    go_terms: array[GOAnnotation]    # Molecular function, process, component
    pathways: array[PathwayRef]      # Reactome, KEGG, WikiPathways
    protein_domains: array[string]   # Pfam, InterPro

  pharmacogenomics:
    haplotypes: array[Haplotype]     # Star alleles (*1, *2, etc.)
    cpic_level: string               # CPIC gene level
    druggable: boolean

  metadata:
    source_databases: array[string]
    last_updated: datetime
    data_quality_score: float
```

### 1.2 Variant Entity Comparison

| Database | Variant ID | Variant Count | Key Attributes | Strengths |
|----------|------------|---------------|----------------|-----------|
| **dbSNP** | rsID | 1B+ | SPDI format, frequencies | Comprehensive coverage |
| **ClinVar** | VCV/RCV/SCV | 2M+ | Clinical significance, conditions | Clinical interpretation |
| **gnomAD** | chrom-pos-ref-alt | 807M+ | Population frequencies (8 groups) | Population data |
| **GWAS Catalog** | rsID + trait | 400K+ | P-values, effect sizes | GWAS associations |
| **PharmGKB** | rsID + haplotype | 10K+ | Drug response annotations | PGx relevance |

**Proposed Variant Schema:**

```yaml
Variant:
  identifiers:
    internal_id: string          # UUID
    rsid: string                 # rs123456789 (primary for SNVs)
    clinvar_vcv: string          # VCV000123456
    gnomad_id: string            # 1-12345678-A-G

  location:
    assembly: string             # GRCh38 (canonical)
    chromosome: string
    position: integer            # 1-based
    reference_allele: string
    alternate_allele: string

  nomenclature:
    hgvs_genomic: string         # NC_000001.11:g.12345A>G
    hgvs_coding: string          # NM_000123.4:c.456A>G
    hgvs_protein: string         # NP_000123.1:p.Lys152Glu
    spdi: string                 # Sequence:Position:Deletion:Insertion

  classification:
    variant_type: enum[SNV, insertion, deletion, indel, CNV, SV]
    consequence: string          # VEP/SnpEff consequence
    impact: enum[HIGH, MODERATE, LOW, MODIFIER]
    biotype: string              # protein_coding, regulatory, etc.

  population_frequencies:
    gnomad_af_global: float
    gnomad_af_popmax: float
    gnomad_populations: map[string, PopulationFrequency]
      # afr, amr, asj, eas, fin, mid, nfe, sas
    dbsnp_caf: float             # 1000 Genomes

  clinical_annotations:
    clinvar_significance: enum[pathogenic, likely_pathogenic, VUS, likely_benign, benign]
    clinvar_review_status: string
    clinvar_conditions: array[DiseaseRef]
    snpedia_magnitude: float     # 0-10 scale
    snpedia_repute: enum[good, bad, neutral]

  functional_predictions:
    cadd_phred: float            # 1-99
    revel_score: float           # 0-1
    spliceai_max: float          # 0-1
    phylop: float                # Conservation
    sift_prediction: string
    polyphen_prediction: string

  pharmacogenomics:
    pharmgkb_level: enum[1A, 1B, 2A, 2B, 3, 4]
    affected_drugs: array[DrugAnnotation]
    haplotype_membership: array[HaplotypeRef]

  gene_context:
    gene_id: string              # Reference to Gene
    gene_symbol: string
    transcript_id: string
    exon_number: integer

  metadata:
    sources: array[SourceEvidence]
    first_observed: datetime
    last_updated: datetime
```

### 1.3 Compound Entity Comparison

| Database | Compound ID | Count | Key Attributes | Strengths |
|----------|-------------|-------|----------------|-----------|
| **ChEMBL** | CHEMBL ID | 2.4M | Bioactivity, assays | Validated activity data |
| **COCONUT** | COCONUT ID | 695K | NP-likeness, source DBs | Natural products focus |
| **FooDB** | FooDB ID | 26K+ | Food sources, concentrations | Dietary compounds |
| **Dr. Duke's** | - | 50K+ | Ethnobotanical uses, LD50 | Traditional use data |
| **PubChem** | CID | 116M+ | Structures, properties | Comprehensive coverage |
| **DrugBank** | DB ID | 15K | Drug properties, interactions | Drug focus |

**Proposed Compound Schema:**

```yaml
Compound:
  identifiers:
    internal_id: string          # UUID
    inchikey: string             # PRIMARY: 27-char standard InChIKey
    pubchem_cid: integer
    chembl_id: string            # CHEMBL123456
    drugbank_id: string          # DB00000
    coconut_id: string
    foodb_id: string
    chebi_id: string

  structure:
    smiles_canonical: string     # Canonical SMILES
    smiles_isomeric: string      # With stereochemistry
    inchi: string                # Full InChI
    molecular_formula: string

  properties:
    molecular_weight: float
    exact_mass: float
    alogp: float                 # Lipophilicity
    hbd: integer                 # H-bond donors
    hba: integer                 # H-bond acceptors
    psa: float                   # Polar surface area
    rotatable_bonds: integer
    aromatic_rings: integer

  classification:
    compound_class: string       # Natural product, synthetic, semi-synthetic
    superclass: string           # ClassyFire superclass
    class: string                # ClassyFire class
    subclass: string
    is_natural_product: boolean
    is_drug: boolean
    is_food_compound: boolean

  drug_properties:
    np_likeness_score: float     # Natural product likeness
    drug_likeness_score: float   # QED or similar
    synthetic_accessibility: float
    lipinski_violations: integer

  admet_predictions:
    absorption:
      intestinal_absorption: float
      caco2_permeability: float
      pgp_substrate: boolean
    distribution:
      vdss: float
      bbb_permeability: boolean
      cns_permeability: boolean
    metabolism:
      cyp_inhibition: map[string, boolean]  # CYP1A2, CYP2C9, etc.
      cyp_substrate: map[string, boolean]
    excretion:
      total_clearance: float
    toxicity:
      ames_mutagenicity: boolean
      herg_inhibition: boolean
      hepatotoxicity: boolean
      ld50_oral: float
      toxicity_class: integer    # 1-6

  source_organisms:
    organisms: array[OrganismSource]
    plant_parts: array[string]   # leaf, root, bark, etc.
    geographic_regions: array[string]

  names:
    preferred_name: string
    iupac_name: string
    common_names: array[string]
    trade_names: array[string]   # For drugs

  metadata:
    sources: array[string]
    data_quality: float
    last_updated: datetime
```

### 1.4 Disease Entity Comparison

| Database | Disease ID | Count | Key Attributes | Strengths |
|----------|------------|-------|----------------|-----------|
| **MONDO** | MONDO ID | 24K+ | Unified ontology, cross-refs | Best integration hub |
| **HPO** | HP ID | 18K | Phenotypes, inheritance | Phenotype detail |
| **Orphanet** | ORPHA ID | 10K+ | Rare diseases, prevalence | Rare disease expertise |
| **DisGeNET** | UMLS CUI | 30K+ | Gene-disease scores | Association strength |
| **OMIM** | MIM ID | 27K+ | Mendelian genetics | Genetic basis |
| **MedGen** | CUI | 100K+ | NCBI integration | NCBI ecosystem |

**Proposed Disease Schema:**

```yaml
Disease:
  identifiers:
    internal_id: string          # UUID
    mondo_id: string             # PRIMARY: MONDO:0000000
    omim_id: string              # OMIM:123456
    orphanet_id: string          # ORPHA:123
    mesh_id: string              # D000000
    doid: string                 # DOID:0000
    efo_id: string               # EFO:0000000
    umls_cui: string             # C0000000
    icd10: array[string]         # Multiple ICD-10 codes
    icd11: array[string]
    medgen_id: string

  names:
    preferred_name: string       # From MONDO
    synonyms: array[string]
    abbreviations: array[string]

  classification:
    disease_category: string     # Cancer, cardiovascular, etc.
    inheritance_pattern: array[enum[AD, AR, XL, XLD, XLR, MT, sporadic]]
    onset_age: string            # Pediatric, adult, etc.

  phenotypes:
    hpo_terms: array[HPOAnnotation]
    clinical_features: array[string]

  epidemiology:
    prevalence: float            # Per 100,000
    prevalence_class: string     # >1/1000, 1-5/10000, etc.
    incidence: float
    geographic_distribution: array[string]

  genetic_basis:
    causative_genes: array[GeneAssociation]
    associated_variants: array[VariantAssociation]

  therapeutic_options:
    approved_treatments: array[DrugRef]
    clinical_trials_count: integer

  metadata:
    sources: array[string]
    curation_status: string
    last_updated: datetime
```

### 1.5 Pathway Entity Comparison

| Database | Pathway ID | Count | Key Attributes | Strengths |
|----------|------------|-------|----------------|-----------|
| **Reactome** | R-HSA-ID | 2,700+ | Curated reactions, events | Quality curation |
| **WikiPathways** | WP ID | 3,000+ | Community pathways | Community coverage |
| **KEGG** | hsa ID | 350+ | Metabolic focus | Metabolic completeness |
| **GO (BP)** | GO:BP | 30K+ | Biological processes | Standardization |

**Proposed Pathway Schema:**

```yaml
Pathway:
  identifiers:
    internal_id: string          # UUID
    reactome_id: string          # PRIMARY: R-HSA-000000
    wikipathways_id: string      # WP000
    kegg_id: string              # hsa00000
    go_bp_id: string             # GO:0000000 (biological process)

  names:
    name: string
    synonyms: array[string]

  classification:
    category: string             # Metabolism, signaling, disease, etc.
    subcategory: string

  components:
    genes: array[GeneRef]
    proteins: array[ProteinRef]
    compounds: array[CompoundRef]  # Small molecule participants
    reactions: array[ReactionRef]

  hierarchy:
    parent_pathway: string       # Reference to parent
    child_pathways: array[string]

  disease_associations:
    diseases: array[DiseaseRef]

  metadata:
    species: string              # Human, mouse, etc.
    evidence_type: string        # Curated, inferred
    source: string
    last_updated: datetime
```

---

## Part 2: Identifier Cross-Reference Matrix

### 2.1 Gene Identifier Mapping

| ID Type | Coverage | Primary Source | Cross-Referenced In |
|---------|----------|----------------|---------------------|
| **NCBI Gene ID** | Universal | NCBI Gene | dbSNP, ClinVar, DisGeNET, KEGG |
| **HGNC ID** | Human genes | HGNC | PharmGKB, ClinVar, Reactome |
| **Ensembl Gene** | Multi-species | Ensembl | gnomAD, UniProt |
| **UniProt ID** | Proteins | UniProt | Reactome, GO, STRING |
| **PharmGKB PA** | PGx genes | PharmGKB | - |

**Hub Recommendation: NCBI Gene ID (Entrez) with UniProt for protein-level**

```
NCBI Gene ID (7157)
  |-> HGNC ID (HGNC:11998)
  |-> Ensembl (ENSG00000141510)
  |-> UniProt (P04637)
  |-> PharmGKB (PA36679)
```

### 2.2 Variant Identifier Mapping

| ID Type | Coverage | Primary Source | Cross-Referenced In |
|---------|----------|----------------|---------------------|
| **rsID** | SNVs | dbSNP | ClinVar, gnomAD, PharmGKB, GWAS |
| **ClinVar VCV** | Clinical | ClinVar | dbSNP |
| **gnomAD variant_id** | Population | gnomAD | - |
| **SPDI** | Any variant | NCBI | dbSNP |
| **HGVS** | Standard | Various | Universal |

**Hub Recommendation: rsID for SNVs, SPDI for normalization**

### 2.3 Compound Identifier Mapping

| ID Type | Coverage | Primary Source | Cross-Referenced In |
|---------|----------|----------------|---------------------|
| **InChIKey** | Universal | IUPAC | PubChem, ChEMBL, COCONUT, all NP DBs |
| **PubChem CID** | 116M+ | PubChem | ChEMBL, DrugBank, BindingDB |
| **ChEMBL ID** | 2.4M | ChEMBL | UniProt, PubChem |
| **DrugBank ID** | 15K drugs | DrugBank | PharmGKB |
| **ChEBI ID** | 195K | ChEBI | Reactome, GO |
| **CAS Number** | Commercial | CAS | Dr. Duke's, many |
| **COCONUT ID** | 695K NPs | COCONUT | - |

**Hub Recommendation: InChIKey (structure-based, universal)**

### 2.4 Disease Identifier Mapping

| ID Type | Coverage | Primary Source | Cross-Referenced In |
|---------|----------|----------------|---------------------|
| **MONDO** | 24K unified | MONDO | Wide adoption |
| **OMIM** | 27K genetic | OMIM | ClinVar, DisGeNET |
| **ORPHA** | 10K rare | Orphanet | MONDO |
| **MeSH** | Medical | NLM | PubMed |
| **ICD-10/11** | Clinical | WHO | EHR systems |
| **UMLS CUI** | Unified | UMLS | DisGeNET |
| **HPO** | Phenotypes | HPO | Orphanet, ClinVar |
| **EFO** | Experimental | EBI | GWAS Catalog |

**Hub Recommendation: MONDO (best cross-reference coverage)**

### 2.5 Pathway Identifier Mapping

| ID Type | Coverage | Primary Source | Cross-Referenced In |
|---------|----------|----------------|---------------------|
| **Reactome** | 2,700 | Reactome | WikiPathways |
| **WikiPathways** | 3,000+ | WikiPathways | - |
| **KEGG** | 350 metabolic | KEGG | Many |
| **GO:BP** | 30K processes | GO | Universal |

**Hub Recommendation: Reactome ID (best curation quality)**

### 2.6 Master Cross-Reference Matrix

```
                    UniProt    NCBI    ChEMBL   PubChem   MONDO   Reactome
                   ID Mapping  ELink   API      PUG       OLS     API
                   =========  ======  =======  ========  ======  =========
Gene IDs             286 DBs    38      Yes      No        No      Yes
Variant IDs          Partial   Yes      No       No        No      No
Compound IDs         ChEBI     No      Yes      116M+     No      ChEBI
Disease IDs          Partial   OMIM    No       No        All     Disease
Pathway IDs          Reactome  No      No       No        GO      All
Literature           PubMed    Yes     Yes      Yes       No      Yes
```

---

## Part 3: Relationship Pattern Analysis

### 3.1 Core Relationship Chains

**Chain 1: Genetic Causation**
```
Gene --[has_variant]--> Variant --[associated_with]--> Disease
  |                       |                              |
  |                       |--[affects]--> Phenotype (HPO)
  |                       |
  |                       |--[found_in]--> Population
  |
  |--[participates_in]--> Pathway
```

**Chain 2: Pharmacological Action**
```
Compound --[binds_to/inhibits]--> Target (Protein/Gene)
    |                                |
    |                                |--[participates_in]--> Pathway
    |                                |
    |                                |--[associated_with]--> Disease
    |
    |--[treats]--> Disease
    |--[causes]--> Adverse Effect
```

**Chain 3: Natural Product Discovery**
```
Organism --[produces]--> Compound --[predicted_target]--> Target
    |                       |                               |
    |                       |--[has_activity]--> Bioassay  |
    |                                                      |
    |--[traditional_use]--> Condition <--[modern_mapping]--+
```

### 3.2 Evidence Type Hierarchy

```yaml
Evidence_Types:
  experimental:
    - binding_assay: "IC50, Ki, Kd measurements"
    - functional_assay: "Cellular/enzymatic activity"
    - clinical_trial: "Human intervention study"
    - genetic_study: "GWAS, family studies"

  computational:
    - predicted_interaction: "Target prediction (SwissTarget, SEA)"
    - structural_similarity: "2D/3D molecular similarity"
    - network_inference: "Pathway/network analysis"
    - machine_learning: "ML-based predictions"

  literature:
    - text_mined: "NLP extraction from papers"
    - curated: "Expert manual curation"
    - meta_analysis: "Systematic review"

  traditional:
    - ethnobotanical: "Traditional use records"
    - pharmacopoeia: "Official monograph"
```

### 3.3 Relationship Cardinalities

| Relationship | Cardinality | Example |
|--------------|-------------|---------|
| Gene -> Variant | 1:N | TP53 has ~30K known variants |
| Variant -> Disease | N:M | One variant, multiple diseases |
| Compound -> Target | N:M | Multi-target drugs |
| Compound -> Organism | N:M | Same compound in multiple species |
| Gene -> Pathway | N:M | Genes in multiple pathways |
| Disease -> Phenotype | 1:N | Disease has multiple phenotypes |
| Drug -> Gene (PGx) | N:M | Drug response genes |

---

## Part 4: Proposed Unified Schema

### 4.1 Entity-Relationship Diagram

```
+------------+       +------------+       +------------+
|    Gene    |-------|  Variant   |-------|  Disease   |
+------------+       +------------+       +------------+
      |                    |                    |
      |                    |                    |
      v                    v                    v
+------------+       +------------+       +------------+
|  Pathway   |       |  Evidence  |       | Phenotype  |
+------------+       +------------+       +------------+
      |                    ^                    ^
      |                    |                    |
      v                    |                    |
+------------+       +------------+       +------------+
|   Target   |<------|  Compound  |-------|  Organism  |
+------------+       +------------+       +------------+
```

### 4.2 Unified Evidence/Provenance Schema

```yaml
Evidence:
  id: string

  source:
    database: string           # e.g., "ChEMBL", "ClinVar"
    database_version: string   # Release version
    record_id: string          # Original record ID
    access_date: datetime

  publication:
    pmid: integer
    pmcid: string
    doi: string
    title: string
    authors: array[string]
    journal: string
    year: integer

  evidence_type:
    category: enum[experimental, computational, literature, traditional]
    method: string             # Specific method used

  quality_metrics:
    confidence_score: float    # 0-1 normalized
    evidence_level: string     # Database-specific (e.g., PharmGKB 1A-4)
    review_status: string      # ClinVar-style review status
    validation_status: enum[validated, predicted, inferred]

  quantitative_data:
    value: float
    unit: string
    value_type: string         # IC50, Ki, p-value, OR, etc.
    error: float
    n_samples: integer

  context:
    species: string
    tissue: string
    cell_line: string
    experimental_conditions: map[string, string]
```

### 4.3 Association Tables

**Gene-Variant Association:**
```yaml
GeneVariantAssociation:
  gene_id: string
  variant_id: string
  consequence: string
  transcript_id: string
  protein_position: integer
  amino_acid_change: string
  sources: array[Evidence]
```

**Variant-Disease Association:**
```yaml
VariantDiseaseAssociation:
  variant_id: string
  disease_id: string
  clinical_significance: string
  inheritance_mode: string
  penetrance: string
  condition_specifics: string
  sources: array[Evidence]
```

**Compound-Target Association:**
```yaml
CompoundTargetAssociation:
  compound_id: string
  target_id: string            # Usually protein/gene
  action_type: enum[inhibitor, agonist, antagonist, modulator, substrate, binder]
  activity_type: string        # IC50, Ki, EC50, etc.
  activity_value: float
  activity_unit: string
  sources: array[Evidence]
```

**Target-Pathway Association:**
```yaml
TargetPathwayAssociation:
  target_id: string
  pathway_id: string
  role: string                 # Enzyme, receptor, transporter, etc.
  sources: array[Evidence]
```

**Compound-Organism Association:**
```yaml
CompoundOrganismAssociation:
  compound_id: string
  organism_id: string          # NCBI Taxonomy ID
  organism_name: string
  plant_part: string
  concentration: float
  concentration_unit: string
  sources: array[Evidence]
```

### 4.4 Canonical ID Resolution Service

```yaml
IdentifierResolution:
  description: "Service to resolve and normalize identifiers"

  supported_types:
    gene: [ncbi_gene, hgnc, ensembl, uniprot, symbol]
    variant: [rsid, clinvar_vcv, hgvs, spdi, gnomad]
    compound: [inchikey, pubchem, chembl, drugbank, chebi]
    disease: [mondo, omim, orphanet, mesh, icd10, umls]
    pathway: [reactome, wikipathways, kegg, go_bp]

  resolution_priority:
    gene: [ncbi_gene_id, hgnc_id, ensembl_gene_id]
    variant: [rsid, clinvar_vcv, hgvs_genomic]
    compound: [inchikey, pubchem_cid, chembl_id]
    disease: [mondo_id, omim_id, orphanet_id]
    pathway: [reactome_id, wikipathways_id, kegg_id]

  cross_reference_sources:
    - UniProt ID Mapping (286 databases)
    - NCBI ELink (38 Entrez databases)
    - PubChem ID Exchange
    - OLS (Ontology Lookup Service)
    - Wikidata SPARQL
```

---

## Part 5: Integration Architecture

### 5.1 Data Flow Architecture

```
+------------------+     +------------------+     +------------------+
|   Source Layer   |     |  Transform Layer |     |   Target Layer   |
+------------------+     +------------------+     +------------------+

  dbSNP (JSON)    --+                          +--> Gene Collection
  ClinVar (TSV)   --+--> ID Resolution   --+   |
  gnomAD (VCF)    --+    Normalization     +---+--> Variant Collection
                   |     Deduplication     |   |
  ChEMBL (SQL)    --+                      |   +--> Compound Collection
  COCONUT (SDF)   --+--> Entity Mapping  --+   |
  PharmGKB (TSV)  --+    Evidence Merge    |   +--> Disease Collection
                   |                       |   |
  Reactome (API)  --+                      +---+--> Pathway Collection
  DisGeNET (TSV)  --+--> Quality Scoring       |
  MONDO (OWL)     --+                          +--> Evidence Collection
                                               |
                                               +--> Association Tables
```

### 5.2 Storage Recommendations

| Collection | Estimated Records | Storage (raw) | Storage (indexed) |
|------------|-------------------|---------------|-------------------|
| Gene | ~25K (human focus) | 500 MB | 1 GB |
| Variant | ~10M (clinically relevant) | 50 GB | 100 GB |
| Compound | ~1M (NP + drugs) | 10 GB | 20 GB |
| Disease | ~50K | 500 MB | 1 GB |
| Pathway | ~10K | 200 MB | 500 MB |
| Evidence | ~100M | 200 GB | 400 GB |
| Associations | ~500M | 100 GB | 200 GB |
| **Total** | - | **~360 GB** | **~720 GB** |

### 5.3 Update Strategy

| Data Type | Update Frequency | Strategy |
|-----------|------------------|----------|
| Variants (dbSNP) | Monthly | Incremental merge |
| Clinical (ClinVar) | Weekly | Full refresh |
| Frequencies (gnomAD) | Major releases | Full refresh |
| Compounds (ChEMBL) | Quarterly | Incremental merge |
| Natural Products | As available | Additive only |
| Pathways | Monthly | Full refresh |
| Evidence | Continuous | Event-driven |

---

## Part 6: Implementation Priorities

### 6.1 Phase 1: Core Infrastructure (Weeks 1-4)

| Task | Source | Priority | Effort |
|------|--------|----------|--------|
| Gene entity model | NCBI Gene, HGNC | Critical | 1 week |
| Variant entity model | dbSNP, ClinVar | Critical | 2 weeks |
| ID resolution service | UniProt, NCBI ELink | Critical | 1 week |

### 6.2 Phase 2: Clinical Layer (Weeks 5-8)

| Task | Source | Priority | Effort |
|------|--------|----------|--------|
| Disease entity model | MONDO, HPO | High | 1 week |
| Variant-Disease associations | ClinVar, DisGeNET | High | 2 weeks |
| Population frequencies | gnomAD | High | 1 week |

### 6.3 Phase 3: Pharmacological Layer (Weeks 9-12)

| Task | Source | Priority | Effort |
|------|--------|----------|--------|
| Compound entity model | ChEMBL, PubChem | High | 2 weeks |
| Compound-Target associations | ChEMBL, BindingDB | High | 2 weeks |

### 6.4 Phase 4: Natural Products Layer (Weeks 13-16)

| Task | Source | Priority | Effort |
|------|--------|----------|--------|
| Natural product integration | COCONUT, LOTUS | Medium | 2 weeks |
| Traditional medicine | TCM DBs, Dr. Duke's | Medium | 2 weeks |

### 6.5 Phase 5: Pathway & Network Layer (Weeks 17-20)

| Task | Source | Priority | Effort |
|------|--------|----------|--------|
| Pathway entity model | Reactome, WikiPathways | Medium | 1 week |
| Gene-Pathway associations | Reactome, GO | Medium | 1 week |
| Network analysis capabilities | STRING, IntAct | Lower | 2 weeks |

---

## Part 7: Quality Assurance Framework

### 7.1 Data Quality Dimensions

| Dimension | Metric | Target |
|-----------|--------|--------|
| Completeness | Required fields populated | >95% |
| Accuracy | Cross-database concordance | >98% |
| Timeliness | Days since source update | <30 days |
| Consistency | ID resolution success rate | >99% |
| Provenance | Evidence traceability | 100% |

### 7.2 Validation Rules

```yaml
Validation_Rules:
  gene:
    - ncbi_gene_id must be valid integer
    - symbol must match HGNC approved symbol
    - chromosome must be valid (1-22, X, Y, MT)

  variant:
    - rsid must match pattern rs\d+
    - position must be positive integer
    - reference/alternate must be valid nucleotides
    - hgvs must be parseable

  compound:
    - inchikey must be 27 characters
    - smiles must be chemically valid
    - molecular_weight must match structure

  disease:
    - mondo_id must exist in MONDO ontology
    - icd10 codes must be valid format

  evidence:
    - pmid must exist in PubMed (if provided)
    - confidence_score must be 0-1
```

### 7.3 Conflict Resolution Strategy

| Conflict Type | Resolution Strategy |
|---------------|---------------------|
| ID mapping conflicts | Use primary source priority list |
| Clinical significance disagreement | Aggregate with ClinVar review status |
| Activity value discrepancies | Report range and mean |
| Name variations | Canonical name from primary source + synonyms |

---

## Appendix A: Database-Specific Field Mappings

### A.1 Gene Field Mappings

| Unified Field | dbSNP | PharmGKB | Reactome | DisGeNET | GO |
|---------------|-------|----------|----------|----------|-----|
| ncbi_gene_id | gene_id | geneId | - | geneId | gene_id |
| hgnc_symbol | symbol | symbol | name | geneSymbol | symbol |
| ensembl_gene_id | - | - | stableId | - | - |
| uniprot_id | - | - | uniprotId | uniprotId | db_object_id |

### A.2 Variant Field Mappings

| Unified Field | dbSNP | ClinVar | gnomAD | PharmGKB |
|---------------|-------|---------|--------|----------|
| rsid | refsnp_id | RS# | - | variantId |
| chromosome | seq_id | Chromosome | chrom | location |
| position | position | Start | pos | - |
| reference | deleted | ReferenceAllele | ref | - |
| alternate | inserted | AlternateAllele | alt | - |
| clinical_significance | clinical | ClinicalSignificance | - | - |

### A.3 Compound Field Mappings

| Unified Field | ChEMBL | COCONUT | FooDB | DrugBank |
|---------------|--------|---------|-------|----------|
| inchikey | standard_inchi_key | inchikey | inchikey | inchi_key |
| smiles | canonical_smiles | smiles | smiles | smiles |
| molecular_weight | full_mwt | molecular_weight | moldb_mono_mass | average_mass |
| pubchem_cid | - | pubchem_cid | pubchem_compound_id | - |

---

## Appendix B: API Integration Summary

| Database | API Type | Rate Limit | Auth | Format |
|----------|----------|------------|------|--------|
| dbSNP | REST | 1/sec | None | JSON |
| ClinVar | E-utilities | 3/sec (10 w/key) | API key | XML/JSON |
| gnomAD | GraphQL | Fair use | None | JSON |
| PharmGKB | REST | 2/sec | None | JSON |
| ChEMBL | REST | None stated | None | JSON |
| UniProt | REST | 1/sec (3 w/email) | None | JSON/TSV |
| Reactome | REST | None stated | None | JSON |
| COCONUT | REST | None stated | None | JSON |
| PubChem | PUG REST | 5/sec | None | JSON/XML |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [index.md](../index.md) | Parent navigation |
| [primary.md](../genetics/primary.md) | Gene/Variant source details |
| [primary.md](../pathways/primary.md) | Pathway source details |
| [pharmaceuticals.md](../compounds/pharmaceuticals.md) | Drug/compound source details |
| [natural-products.md](../compounds/natural-products.md) | Natural product source details |
| [xrefs.md](../integration/xrefs.md) | Cross-reference details |
| [world1-schema-research.md](../research/world1-schema-research.md) | Genetics schema research |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Architecture Team | Initial unified schema analysis |
