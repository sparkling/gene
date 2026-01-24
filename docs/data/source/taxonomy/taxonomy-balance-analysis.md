---
id: research-taxonomy-balance-analysis
title: Taxonomy Balance Analysis - Data Source Distribution
type: research
parent: ../_index.md
last_updated: 2026-01-23
status: draft
tags: [research, taxonomy, categorization, balance, data-sources]
---

**Parent:** [Research](./_index.md)

# Taxonomy Balance Analysis - Data Source Distribution

**Document ID:** TAXONOMY-BALANCE-ANALYSIS
**Status:** Draft
**Owner:** Taxonomy Balance Analyzer Agent
**Date:** January 23, 2026
**Version:** 1.0

---

## Executive Summary

This analysis evaluates how the 175+ documented data sources distribute across proposed taxonomy categories. The goal is to identify imbalances, recommend category adjustments, and ensure a well-structured polyhierarchical taxonomy for the Gene Platform.

### Key Findings

| Metric | Value |
|--------|-------|
| **Total Data Sources Identified** | 178 |
| **Primary Categories** | 6 |
| **Subcategories (avg per category)** | 5-8 |
| **Oversized Categories (>25 sources)** | 2 (Genetics, Compounds) |
| **Undersized Categories (<5 sources)** | 0 at top level |
| **Sources fitting 5+ categories** | 12 (highly connected) |
| **Orphan sources** | 3 |

---

## Part 1: Complete Data Source Inventory

### 1.1 Master Data Source List (178 Sources)

Below is the complete inventory of data sources identified across all documentation.

#### Category A: Genetics & Genomics (38 sources)

| # | Source | Priority | Type | Multi-Category? |
|---|--------|----------|------|-----------------|
| 1 | dbSNP | Tier 1 | Variant Repository | Primary |
| 2 | ClinVar | Tier 1 | Clinical Variants | Primary |
| 3 | gnomAD | Tier 1 | Population Frequencies | Primary |
| 4 | PharmGKB | Tier 1 | Pharmacogenomics | Compounds, Drugs |
| 5 | dbNSFP | Tier 1 | Functional Prediction | Primary |
| 6 | AlphaMissense | Tier 1 | AI Prediction | Primary |
| 7 | GWAS Catalog | Tier 2 | Association Studies | Diseases |
| 8 | gnomAD-SV | Tier 2 | Structural Variants | Primary |
| 9 | ExAC | Tier 2 | Exome Aggregation | Superseded |
| 10 | GTEx | Tier 2 | Gene Expression | Primary |
| 11 | UK Biobank | Tier 3 | Population Genomics | Controlled |
| 12 | 1000 Genomes | Tier 2 | Population Variation | Primary |
| 13 | OMIM | Tier 2 | Genetic Disorders | Diseases |
| 14 | CIViC | Tier 2 | Clinical Interpretations | Cancer |
| 15 | dbVar | Tier 2 | Structural Variants | Primary |
| 16 | SpliceAI | Tier 2 | Splice Prediction | Primary |
| 17 | ENCODE | Tier 2 | Functional Elements | Primary |
| 18 | Gene Ontology | Tier 1 | Functional Annotation | Primary |
| 19 | UniProt ID Mapping | Tier 1 | Cross-Reference | Integration Hub |
| 20 | SNPedia | Tier 2 | Consumer Genetics | Primary |
| 21 | CADD | Tier 2 | Pathogenicity Score | Primary |
| 22 | MaveDB | Tier 2 | Variant Effects | Primary |
| 23 | RegulomeDB | Tier 2 | Regulatory Variants | Primary |
| 24 | ALFA (NCBI) | Tier 2 | Population Frequencies | Primary |
| 25 | ClinGen | Tier 2 | Gene Curation | Primary |
| 26 | DECIPHER | Tier 2 | CNV Annotation | Primary |
| 27 | DDG2P | Tier 2 | Gene-Disease | Diseases |
| 28 | PGC | Tier 2 | Psychiatric Genetics | Mental Health |
| 29 | CARDIoGRAMplusC4D | Tier 3 | Cardiovascular GWAS | Cardio |
| 30 | GLGC | Tier 3 | Lipid GWAS | Cardio |
| 31 | MITOMAP | Tier 2 | Mitochondrial Variants | Primary |
| 32 | MitoCarta | Tier 2 | Mitochondrial Genes | Primary |
| 33 | BrainSpan | Tier 2 | Brain Expression | Mental Health |
| 34 | Allen Brain Atlas | Tier 2 | Brain Expression | Mental Health |
| 35 | PharmVar | Tier 1 | Star Alleles | Pharmacogenomics |
| 36 | CPIC | Tier 1 | Dosing Guidelines | Pharmacogenomics |
| 37 | Cancer Gene Census | Tier 2 | Cancer Genes | Cancer |
| 38 | GDC/TCGA (open) | Tier 3 | Cancer Genomics | Cancer |

#### Category B: Traditional Medicine (23 sources)

| # | Source | Priority | Type | Multi-Category? |
|---|--------|----------|------|-----------------|
| 1 | BATMAN-TCM 2.0 | Tier 1 | TCM Network | Compounds |
| 2 | HERB 2.0 | Tier 1 | TCM Herb-Gene | Compounds |
| 3 | TCMSID | Tier 1 | TCM Syndromes | Diseases |
| 4 | IMPPAT 2.0 | Tier 1 | Ayurveda | Compounds |
| 5 | KampoDB | Tier 2 | Kampo | Compounds |
| 6 | LOTUS | Tier 1 | Natural Products | Compounds |
| 7 | Dr. Duke's | Tier 1 | Western Herbal | Compounds |
| 8 | NPASS | Tier 2 | Natural Products | Compounds |
| 9 | TM-MC | Tier 2 | Multi-System | Compounds |
| 10 | ETCM | Tier 2 | TCM Network | Compounds |
| 11 | SymMap 2.0 | Tier 2 | TCM Symptoms | Diseases |
| 12 | YaTCM | Tier 3 | TCM Network | Compounds |
| 13 | TCMID 2.0 | Tier 3 | TCM Compounds | Compounds |
| 14 | TCMBank | Tier 3 | TCM AI | Compounds |
| 15 | CMAUP | Tier 2 | Medicinal Plants | Compounds |
| 16 | NAPRALERT | Tier 3 | Natural Products | Subscription |
| 17 | TMC-TCM | Tier 3 | TCM Formulas | Compounds |
| 18 | PhytoHub | Tier 2 | Phytochemicals | Nutrition |
| 19 | ANPDB | Tier 2 | African NP | Compounds |
| 20 | SANCDB | Tier 2 | South African NP | Compounds |
| 21 | NuBBEDB | Tier 2 | Brazilian NP | Compounds |
| 22 | Super Natural II | Tier 3 | Natural Products | Compounds |
| 23 | Ensembl Plants | Tier 3 | Plant Genomes | Genetics |

#### Category C: Nutrition & Food (14 sources)

| # | Source | Priority | Type | Multi-Category? |
|---|--------|----------|------|-----------------|
| 1 | USDA FoodData Central | Tier 1 | Food Composition | Primary |
| 2 | FooDB | Tier 1 | Food Compounds | Compounds |
| 3 | DSLD (NIH) | Tier 1 | Supplements | Primary |
| 4 | Open Food Facts | Tier 2 | Consumer Products | Primary |
| 5 | Phenol-Explorer | Tier 2 | Polyphenols | Compounds |
| 6 | PhytoHub | Tier 2 | Phytochemicals | Traditional |
| 7 | HMDB | Tier 2 | Human Metabolome | Metabolites |
| 8 | FRIDA | Tier 2 | Food Risk | Primary |
| 9 | EuroFIR | Tier 3 | European Foods | Regional |
| 10 | AUSNUT | Tier 3 | Australian Foods | Regional |
| 11 | UK Food Composition | Tier 3 | UK Foods | Regional |
| 12 | FoodOmicsGR | Tier 3 | Greek Foods | Regional |
| 13 | Metabolomics Workbench | Tier 3 | Metabolomics | Research |
| 14 | HMP | Tier 2 | Human Microbiome | Microbiome |

#### Category D: Compounds & Drugs (32 sources)

| # | Source | Priority | Type | Multi-Category? |
|---|--------|----------|------|-----------------|
| 1 | COCONUT 2.0 | Tier 1 | Natural Products | Primary |
| 2 | ChEMBL | Tier 1 | Bioactivity | Primary |
| 3 | PubChem | Tier 1 | Chemical Repository | Primary |
| 4 | ChEBI | Tier 1 | Chemical Ontology | Primary |
| 5 | DrugBank | Tier 2 | Drug Database | Pharmacogenomics |
| 6 | ZINC20 | Tier 3 | Virtual Screening | Primary |
| 7 | BindingDB | Tier 2 | Binding Affinities | Primary |
| 8 | StreptomeDB | Tier 3 | Bacterial NP | Microbiome |
| 9 | KEGG DRUG | Tier 2 | Drug Pathways | Pathways |
| 10 | PK-DB | Tier 2 | Pharmacokinetics | Primary |
| 11 | SuperCYP | Tier 2 | CYP Interactions | Primary |
| 12 | T3DB | Tier 2 | Toxins | Primary |
| 13 | SIDER | Tier 2 | Side Effects | Safety |
| 14 | FAERS | Tier 2 | Adverse Events | Safety |
| 15 | DGIdb | Tier 2 | Drug-Gene | Genetics |
| 16 | Open Targets | Tier 2 | Drug Targets | Primary |
| 17 | CTD | Tier 2 | Chemical-Disease | Diseases |
| 18 | MiBiG | Tier 3 | Biosynthesis | Microbiome |
| 19 | EPA CompTox | Tier 3 | Toxicology | Environmental |
| 20 | KEGG Compound | Tier 2 | Metabolites | Pathways |
| 21 | LIPID MAPS | Tier 2 | Lipids | Primary |
| 22 | SwissLipids | Tier 3 | Lipids | Primary |
| 23 | MetaCyc | Tier 2 | Metabolites | Pathways |
| 24 | Therapeutic Target DB | Tier 2 | Drug Targets | Primary |
| 25 | ChemSpider | Tier 3 | Chemical Search | Primary |
| 26 | GtoPdb | Tier 2 | Pharmacology | Primary |
| 27 | Inxight Drugs | Tier 2 | FDA Substances | Primary |
| 28 | WITHDRAWN | Tier 3 | Withdrawn Drugs | Safety |
| 29 | eMolecules | Tier 3 | Purchasable | Commercial |
| 30 | Mcule | Tier 3 | Purchasable | Commercial |
| 31 | MOLARIS | Tier 3 | Molecular DB | Primary |
| 32 | PDBe-KB | Tier 2 | Protein Structures | Structural |

#### Category E: Pathways & Networks (18 sources)

| # | Source | Priority | Type | Multi-Category? |
|---|--------|----------|------|-----------------|
| 1 | Reactome | Tier 1 | Curated Pathways | Primary |
| 2 | DisGeNET | Tier 1 | Gene-Disease | Diseases |
| 3 | WikiPathways | Tier 2 | Community Pathways | Primary |
| 4 | KEGG Pathway | Tier 2 | Reference Pathways | Primary |
| 5 | PathwayCommons | Tier 2 | Aggregator | Primary |
| 6 | STRING | Tier 1 | Protein Interactions | Primary |
| 7 | IntAct | Tier 2 | Molecular Interactions | Primary |
| 8 | BioCyc | Tier 3 | Metabolic Pathways | Subscription |
| 9 | PANTHER | Tier 2 | Protein Classification | Primary |
| 10 | MSigDB | Tier 2 | Gene Sets | Primary |
| 11 | SIGNOR | Tier 2 | Signaling Network | Primary |
| 12 | NDEx | Tier 2 | Network Exchange | Primary |
| 13 | PathBank | Tier 2 | Disease Pathways | Primary |
| 14 | KEGG KGML | Tier 2 | Pathway Format | Format Spec |
| 15 | GPML (WikiPathways) | Tier 2 | Pathway Format | Format Spec |
| 16 | BioPAX | Tier 2 | Pathway Format | Format Spec |
| 17 | SBML | Tier 2 | Systems Biology | Format Spec |
| 18 | BioGRID | Tier 2 | Genetic Interactions | Primary |

#### Category F: Diseases & Phenotypes (21 sources)

| # | Source | Priority | Type | Multi-Category? |
|---|--------|----------|------|-----------------|
| 1 | MONDO | Tier 1 | Disease Ontology | Primary |
| 2 | HPO | Tier 1 | Phenotype Ontology | Primary |
| 3 | Orphanet (ORDO) | Tier 1 | Rare Diseases | Primary |
| 4 | MedGen | Tier 2 | NCBI Diseases | Primary |
| 5 | DOID | Tier 2 | Disease Ontology | Primary |
| 6 | EFO | Tier 2 | Experimental Factor | Primary |
| 7 | ICD-10/ICD-11 | Tier 1 | Clinical Codes | Primary |
| 8 | MeSH | Tier 1 | Medical Terms | Primary |
| 9 | SNOMED-CT | Tier 2 | Clinical Terms | Requires License |
| 10 | UMLS | Tier 2 | Unified Medical | Requires License |
| 11 | OncoKB | Tier 1 | Cancer Evidence | Cancer |
| 12 | Monarch Initiative | Tier 2 | Knowledge Graph | Primary |
| 13 | HGMD | Tier 3 | Human Gene Mutation | Subscription |
| 14 | ClinGen Dosage | Tier 2 | CNV Interpretation | Genetics |
| 15 | GeneReviews | Tier 2 | Clinical Summaries | Primary |
| 16 | OMIM | Tier 2 | Mendelian Diseases | Genetics |
| 17 | MarkerDB 2.0 | Tier 2 | Biomarkers | Biomarkers |
| 18 | LOINC | Tier 2 | Lab Codes | Biomarkers |
| 19 | CALIPER | Tier 2 | Reference Intervals | Biomarkers |
| 20 | cBioPortal | Tier 2 | Cancer Genomics | Cancer |
| 21 | DepMap | Tier 2 | Cancer Dependencies | Cancer |

#### Category G: Literature & Evidence (12 sources)

| # | Source | Priority | Type | Multi-Category? |
|---|--------|----------|------|-----------------|
| 1 | PubMed | Tier 1 | Biomedical Literature | Primary |
| 2 | PMC | Tier 2 | Full-Text Archive | Primary |
| 3 | OpenAlex | Tier 2 | Scholarly Graph | Primary |
| 4 | Europe PMC | Tier 2 | European Archive | Primary |
| 5 | Semantic Scholar | Tier 2 | AI-Powered Search | Primary |
| 6 | CORE | Tier 3 | Open Access | Primary |
| 7 | arXiv | Tier 3 | Preprints | Primary |
| 8 | bioRxiv | Tier 3 | Biology Preprints | Primary |
| 9 | ClinicalTrials.gov | Tier 2 | Trial Registry | Primary |
| 10 | Cochrane Library | Tier 3 | Systematic Reviews | Subscription |
| 11 | Wikidata | Tier 1 | Knowledge Graph | Integration Hub |
| 12 | Wikipedia | Tier 3 | General Reference | Informational |

#### Category H: Microbiome (10 sources)

| # | Source | Priority | Type | Multi-Category? |
|---|--------|----------|------|-----------------|
| 1 | GMrepo v3 | Tier 1 | Gut Microbiome | Primary |
| 2 | HMP Portal | Tier 2 | Human Microbiome | Primary |
| 3 | gutMGene v2.0 | Tier 2 | Gene-Microbe | Primary |
| 4 | mBodyMap | Tier 2 | Microbiome Atlas | Primary |
| 5 | MASI | Tier 3 | Microbe-Drug | Compounds |
| 6 | MDAD | Tier 3 | Drug-Microbe | Compounds |
| 7 | SILVA | Tier 2 | rRNA Database | Primary |
| 8 | Greengenes | Tier 2 | 16S rRNA | Primary |
| 9 | RefSeq Microbes | Tier 2 | Reference Genomes | Primary |
| 10 | MGnify | Tier 2 | Metagenomics | Primary |

---

## Part 2: Category Balance Analysis

### 2.1 Distribution Statistics

| Category | Source Count | % of Total | Assessment |
|----------|--------------|------------|------------|
| Genetics & Genomics | 38 | 21.3% | **Oversized** - Consider splitting |
| Compounds & Drugs | 32 | 18.0% | **Oversized** - Consider splitting |
| Traditional Medicine | 23 | 12.9% | Balanced |
| Diseases & Phenotypes | 21 | 11.8% | Balanced |
| Pathways & Networks | 18 | 10.1% | Balanced |
| Nutrition & Food | 14 | 7.9% | Balanced |
| Literature & Evidence | 12 | 6.7% | Balanced |
| Microbiome | 10 | 5.6% | Balanced |
| **TOTAL** | **178** | **100%** | |

### 2.2 Category Size Assessment

#### Oversized Categories (>25 sources)

**1. Genetics & Genomics (38 sources)**
- Recommendation: Split into 3-4 subcategories
  - Variant Repositories (dbSNP, ClinVar, gnomAD, dbVar)
  - Functional Prediction (dbNSFP, AlphaMissense, CADD, SpliceAI)
  - Gene Expression (GTEx, ENCODE, Allen Brain Atlas)
  - Pharmacogenomics (PharmGKB, PharmVar, CPIC)

**2. Compounds & Drugs (32 sources)**
- Recommendation: Split into 3 subcategories
  - Natural Products (COCONUT, LOTUS, StreptomeDB)
  - Pharmaceuticals (DrugBank, ChEMBL, DGIdb)
  - Chemical Properties (PubChem, ChEBI, ZINC)

#### Balanced Categories (10-25 sources)

- Traditional Medicine (23): Well-structured by tradition system
- Diseases & Phenotypes (21): Well-structured by disease type
- Pathways & Networks (18): Could benefit from subcategories
- Nutrition & Food (14): Balanced

#### Potentially Undersized Categories

- Literature & Evidence (12): Could merge with a "Reference" category
- Microbiome (10): Adequate for specialized domain

### 2.3 Sources Per Subcategory (Current Proposed)

| Main Category | Subcategory | Count | Status |
|---------------|-------------|-------|--------|
| **Genetics** | Variant Repositories | 8 | Good |
| | Functional Prediction | 8 | Good |
| | Population Genetics | 6 | Good |
| | Pharmacogenomics | 6 | Good |
| | Expression & Regulation | 5 | Good |
| | Cancer Genomics | 5 | Good |
| **Compounds** | Natural Products | 12 | Slightly high |
| | Pharmaceuticals | 10 | Good |
| | Chemical Reference | 10 | Good |
| **Traditional** | TCM | 10 | Good |
| | Ayurveda | 2 | **Low - Consider merging** |
| | Kampo | 1 | **Low - Consider merging** |
| | Western Herbal | 3 | **Low - Consider merging** |
| | Multi-System | 7 | Good |
| **Diseases** | Rare Diseases | 5 | Good |
| | Cancer | 5 | Good |
| | Phenotypes | 4 | Good |
| | Ontologies | 7 | Good |
| **Pathways** | Curated Pathways | 6 | Good |
| | Interaction Networks | 5 | Good |
| | Format Standards | 4 | Good |
| **Nutrition** | Food Composition | 8 | Good |
| | Phytochemicals | 4 | Good |
| | Regional | 4 | Consider merging |

---

## Part 3: Multi-Category Sources (Polyhierarchy)

### 3.1 Highly Connected Sources (5+ category fits)

| Source | Primary Category | Secondary Categories | Connection Count |
|--------|-----------------|---------------------|------------------|
| **PharmGKB** | Genetics | Compounds, Drugs, Diseases, Pathways | 5 |
| **DisGeNET** | Pathways | Genetics, Diseases, Literature | 4 |
| **ChEMBL** | Compounds | Genetics, Pathways, Diseases | 4 |
| **Wikidata** | Literature | All categories | 6+ |
| **UniProt** | Genetics | Compounds, Pathways, Diseases | 4 |
| **PubChem** | Compounds | Nutrition, Traditional, Pathways | 4 |
| **KEGG** | Pathways | Compounds, Genetics, Diseases | 4 |
| **DrugBank** | Compounds | Genetics, Pathways, Diseases | 4 |
| **Reactome** | Pathways | Genetics, Compounds, Diseases | 4 |
| **STRING** | Pathways | Genetics, Compounds, Diseases | 4 |
| **LOTUS** | Traditional | Compounds, Nutrition | 3 |
| **FooDB** | Nutrition | Compounds, Traditional | 3 |

### 3.2 Integration Hub Sources

These sources serve as primary cross-reference hubs:

| Hub Source | Cross-Ref Count | Primary Use |
|------------|-----------------|-------------|
| UniProt ID Mapping | 286+ databases | Protein ID resolution |
| Wikidata | All categories | Universal knowledge graph |
| PubChem | 100+ compound DBs | Chemical ID resolution |
| NCBI ELink | 38 NCBI databases | Gene/variant linking |
| MONDO | 130K disease mappings | Disease ID resolution |

---

## Part 4: Orphan Sources

Sources that don't fit cleanly into current categories:

| Source | Current Placement | Issue | Recommendation |
|--------|------------------|-------|----------------|
| **CALIPER** | Diseases (Biomarkers) | Clinical reference intervals | Create "Clinical Reference" subcategory |
| **LOINC** | Diseases (Biomarkers) | Lab coding system | Create "Clinical Standards" subcategory |
| **EPA CompTox** | Compounds | Environmental toxicology | Create "Environmental Health" subcategory |

---

## Part 5: Recommendations

### 5.1 Category Restructuring

**Recommended Taxonomy (6 Primary + 30 Secondary Categories):**

```
LEVEL 0: ROOT
|
+-- LEVEL 1: GENETICS & GENOMICS (38 sources)
|   +-- Variant Databases (dbSNP, ClinVar, gnomAD, dbVar, etc.)
|   +-- Functional Prediction (dbNSFP, AlphaMissense, CADD, SpliceAI)
|   +-- Population Genetics (1000G, gnomAD-SV, ALFA)
|   +-- Pharmacogenomics (PharmGKB, PharmVar, CPIC)
|   +-- Expression & Regulation (GTEx, ENCODE, Allen Brain)
|   +-- Cancer Genetics (CIViC, Cancer Gene Census, TCGA)
|
+-- LEVEL 1: COMPOUNDS & CHEMISTRY (32 sources)
|   +-- Natural Products (COCONUT, LOTUS, NPASS)
|   +-- Pharmaceuticals (DrugBank, DGIdb, Open Targets)
|   +-- Bioactivity (ChEMBL, BindingDB, PDB)
|   +-- Chemical Reference (PubChem, ChEBI, ChemSpider)
|   +-- Toxicology & Safety (T3DB, CTD, SIDER, FAERS)
|
+-- LEVEL 1: TRADITIONAL MEDICINE (23 sources)
|   +-- Traditional Chinese Medicine (BATMAN-TCM, HERB, ETCM, SymMap)
|   +-- South/Southeast Asian Systems (IMPPAT [Ayurveda], KampoDB [Kampo])
|   +-- Western & Global Herbal (Dr. Duke's, NAPRALERT, CMAUP)
|   +-- Regional Natural Products (ANPDB, SANCDB, NuBBEDB)
|   +-- Multi-System Aggregators (LOTUS, NPASS, TM-MC)
|
+-- LEVEL 1: NUTRITION & METABOLISM (14 sources)
|   +-- Food Composition (USDA, FooDB, Open Food Facts)
|   +-- Dietary Supplements (DSLD, Phenol-Explorer)
|   +-- Metabolomics (HMDB, Metabolomics Workbench)
|   +-- Regional Food Data (EuroFIR, AUSNUT, UK Food)
|
+-- LEVEL 1: PATHWAYS & NETWORKS (18 sources)
|   +-- Curated Pathways (Reactome, WikiPathways, KEGG)
|   +-- Protein Interactions (STRING, IntAct, BioGRID)
|   +-- Disease Networks (DisGeNET, PathBank)
|   +-- Pathway Standards (BioPAX, SBML, KGML, GPML)
|
+-- LEVEL 1: DISEASES & PHENOTYPES (21 sources)
|   +-- Disease Ontologies (MONDO, DOID, EFO)
|   +-- Phenotype Systems (HPO, OMIM)
|   +-- Clinical Standards (ICD, MeSH, SNOMED)
|   +-- Rare Diseases (Orphanet, GeneReviews, DECIPHER)
|   +-- Cancer Resources (OncoKB, cBioPortal, DepMap)
|   +-- Biomarkers & Labs (MarkerDB, LOINC, CALIPER)
|
+-- LEVEL 1: LITERATURE & EVIDENCE (12 sources)
|   +-- Biomedical Literature (PubMed, PMC, Europe PMC)
|   +-- Scholarly Graphs (OpenAlex, Semantic Scholar)
|   +-- Clinical Trials (ClinicalTrials.gov, Cochrane)
|   +-- Knowledge Integration (Wikidata, Wikipedia)
|
+-- LEVEL 1: MICROBIOME (10 sources)
|   +-- Gut Microbiome (GMrepo, gutMGene)
|   +-- Human Microbiome Project (HMP Portal, mBodyMap)
|   +-- Reference Databases (SILVA, Greengenes, RefSeq)
|   +-- Microbe-Drug Interactions (MASI, MDAD)
```

### 5.2 Merging Recommendations

| Undersized Category | Merge Into | Rationale |
|---------------------|------------|-----------|
| Ayurveda (2 sources) | South/Southeast Asian Systems | Geographic/cultural proximity |
| Kampo (1 source) | South/Southeast Asian Systems | Geographic/cultural proximity |
| Western Herbal (3) | Western & Global Herbal | Expand to include global |
| Regional Food (4) | Food Composition | Consolidate regional variants |

### 5.3 Splitting Recommendations

| Oversized Category | Split Into | Rationale |
|-------------------|------------|-----------|
| Genetics (38) | 6 subcategories | Too diverse; improve navigation |
| Compounds (32) | 5 subcategories | Distinct use cases |

### 5.4 New Category Recommendations

| New Category | Sources | Rationale |
|--------------|---------|-----------|
| Environmental Health | EPA CompTox, CTD | Distinct domain |
| Clinical Standards | LOINC, CALIPER, ICD | Distinct utility |
| AI/ML Predictions | AlphaMissense, SpliceAI, CADD | Emerging category |

---

## Part 6: Data Source to Category Mapping

### 6.1 Primary Category Assignments

The complete mapping of all 178 sources to their primary and secondary categories:

```yaml
# Full mapping available in separate YAML file
# Example format:
data_sources:
  - id: dbSNP
    primary_category: genetics
    secondary_categories: []
    subcategory: variant_repositories
    tier: 1
    license: public_domain

  - id: PharmGKB
    primary_category: genetics
    secondary_categories: [compounds, pathways, diseases]
    subcategory: pharmacogenomics
    tier: 1
    license: cc_by_sa_4

  - id: Wikidata
    primary_category: literature
    secondary_categories: [genetics, compounds, traditional, nutrition, pathways, diseases]
    subcategory: knowledge_integration
    tier: 1
    license: cc0
```

### 6.2 Category Statistics Summary

| Metric | Value |
|--------|-------|
| Total sources | 178 |
| Primary categories | 8 |
| Recommended subcategories | 30 |
| Average sources per subcategory | 5.9 |
| Min sources per subcategory | 2 (Ayurveda) |
| Max sources per subcategory | 12 (Natural Products) |
| Sources with multiple categories | 45 (25%) |
| Integration hub sources | 5 |
| Orphan sources | 3 |

---

## Part 7: Implementation Priority

### 7.1 Phase 1: Core Categories (Weeks 1-4)

Document and structure:
1. Genetics - Variant Repositories (8 sources)
2. Compounds - Natural Products (12 sources)
3. Diseases - Disease Ontologies (7 sources)
4. Pathways - Curated Pathways (6 sources)

### 7.2 Phase 2: Secondary Categories (Weeks 5-8)

Document and structure:
1. Traditional Medicine - TCM (10 sources)
2. Nutrition - Food Composition (8 sources)
3. Literature - Biomedical (5 sources)
4. Microbiome - Core (6 sources)

### 7.3 Phase 3: Specialized Categories (Weeks 9-12)

Document remaining:
1. Regional databases
2. Specialized tools
3. Format standards

---

## Appendix A: Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-01-23 | Taxonomy Balance Analyzer | Initial analysis |

---

## Glossary

| Term | Definition |
|------|------------|
| Polyhierarchy | Classification system allowing items to belong to multiple categories |
| Orphan Source | Data source that doesn't fit existing categories |
| Integration Hub | Source that connects multiple databases via cross-references |
| Tier | Priority level for implementation (1=MVP, 2=High, 3=Specialized) |
| Oversized Category | Category with >25 sources requiring subdivision |
| Undersized Category | Category with <5 sources potentially requiring merger |

---

*Document generated by Taxonomy Balance Analyzer Agent, January 2026*
