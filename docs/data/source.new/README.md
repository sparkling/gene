# Data Source Catalog

**Version:** 2.0.0
**Last Updated:** January 2026
**Status:** Active
**Total Files:** 116 markdown documents across 14 directories

---

## Overview

This catalog organizes all data sources for the Gene Platform by **data category**—a unified approach integrating genetics, traditional medicine, and nutrition data into a coherent knowledge system.

The documentation follows a dual categorization approach:
1. **Physical Organization** - Folder structure organizing files by data type and function
2. **Content Categories** - Data categories organizing knowledge by domain (genetics, traditional, nutrition, shared)

---

## Data Categories

The Gene Platform bridges complementary knowledge systems:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                           DATA CATEGORIES                               │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│   ┌─────────────┐    ┌─────────────┐    ┌─────────────┐                │
│   │  GENETICS   │    │ TRADITIONAL │    │  NUTRITION  │                │
│   │             │    │  MEDICINE   │    │             │                │
│   │             │    │             │    │             │                │
│   ├─────────────┤    ├─────────────┤    ├─────────────┤                │
│   │ • Variants  │    │ • TCM       │    │ • Foods     │                │
│   │ • Genes     │    │ • Ayurveda  │    │ • Nutrients │                │
│   │ • Proteins  │    │ • Kampo     │    │ • Metabolites│               │
│   │ • Pathways  │    │ • Western   │    │ • Bioactives│                │
│   │             │    │   Herbal    │    │             │                │
│   │             │    │ • Global    │    │             │                │
│   └──────┬──────┘    └──────┬──────┘    └──────┬──────┘                │
│          │                  │                  │                        │
│          └──────────────────┼──────────────────┘                        │
│                             │                                           │
│                    ┌────────▼────────┐                                  │
│                    │     SHARED      │                                  │
│                    │   RESOURCES     │                                  │
│                    ├─────────────────┤                                  │
│                    │ • Pathways      │                                  │
│                    │ • Literature    │                                  │
│                    │ • Compounds     │                                  │
│                    └─────────────────┘                                  │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### Genetics

Modern genomic and molecular biology data including human genome sequences, variants, gene expression, protein structures, and clinical genomic associations.

| Tier | Databases | Coverage |
|------|-----------|----------|
| **Tier 1 (MVP)** | dbSNP, ClinVar, gnomAD v4.1, dbNSFP v4.9, AlphaMissense | Core variant annotation |
| **Tier 2** | TOPMed, UK Biobank, CADD, SpliceAI, ENCODE 4 | Population & functional |
| **Tier 3** | gnomAD-SV, DGV, dbVar, DECIPHER | Structural variants |

**Location:** `databases/genetics/`
**Frontmatter:** `category: genetics`

### Traditional Medicine

Time-tested healing knowledge from diverse cultures, organized into five subcategories:

| Subcategory | Key Databases | Content |
|-------------|---------------|---------|
| **TCM** | TCMSP, TCMID, SymMap, ETCM, HERB, TCMBank | 8,000+ compounds, 50+ herbs, formulas |
| **Ayurveda** | IMPPAT 2.0, AyurMedBase, InDiaMed | 9,500+ phytochemicals, 1,700+ plants |
| **Kampo** | KCONSORT, Kampo Medicine Database | 200+ standardized formulas |
| **Western Herbal** | NAPRALERT, Dr. Duke's, HerbMedPro | 200,000+ articles on natural products |
| **Global** | ANPDB, SANCDB, NuBBEDB, BIOFACQUIM | African & Latin American ethnobotany |

**Location:** `databases/traditional/`
**Frontmatter:** `category: traditional`, `subcategory: tcm|ayurveda|kampo|western-herbal|global`

### Nutrition

Evidence-based nutritional science including nutrient composition, dietary patterns, nutrigenomics, and metabolic effects.

| Tier | Databases | Content |
|------|-----------|---------|
| **Tier 1 (MVP)** | FooDB, USDA FoodData Central | 800+ foods, 350K+ food items |
| **Tier 2** | Phenol-Explorer, PhytoHub, HMDB | 114K+ metabolites |
| **Tier 3** | EuroFIR, AUSNUT, FoodOmicsGR | Regional composition data |

**Location:** `databases/nutrition/`
**Frontmatter:** `category: nutrition`

### Shared Resources

Databases that bridge categories through shared identifiers and relationships:

| Category | Databases | Purpose |
|----------|-----------|---------|
| **Pathways** | Reactome, DisGeNET, KEGG, WikiPathways, STRING | Gene → Pathway → Disease mapping |
| **Literature** | PubMed (36M+ citations), PMC, OpenAlex | Evidence synthesis across categories |
| **Compounds** | COCONUT (400K+), LOTUS (750K+), ChEMBL, DrugBank | Chemical structures & bioactivity |

**Location:** `databases/pathways/`, `databases/literature/`, `databases/compounds/`
**Frontmatter:** `category: shared`

---

## Folder Structure

```
source.new/
├── _index.md                    # Main catalog index
├── README.md                    # This file
│
├── databases/                   # PRIMARY DATA SOURCE COLLECTIONS
│   ├── _index.md               # Database overview (26 files total)
│   │
│   ├── genetics/               # Genetics: Genetic variation
│   │   ├── _index.md          # Genetics overview
│   │   ├── primary.md         # Core variant databases (dbSNP, ClinVar, gnomAD)
│   │   └── population.md      # Population genetics (TOPMed, UK Biobank)
│   │
│   ├── traditional/            # Traditional: Traditional medicine systems
│   │   ├── _index.md          # Traditional medicine overview
│   │   ├── tcm.md             # Traditional Chinese Medicine
│   │   ├── ayurveda.md        # Ayurvedic medicine
│   │   ├── kampo.md           # Japanese Kampo
│   │   ├── western-herbal.md  # Western herbal medicine
│   │   └── global.md          # Global/indigenous systems
│   │
│   ├── nutrition/              # Nutrition: Nutritional data
│   │   └── _index.md          # Nutrition overview (schemas in operations/)
│   │
│   ├── pathways/               # Shared: Biological pathways
│   │   ├── _index.md
│   │   ├── primary.md         # Reactome, KEGG
│   │   ├── disease.md         # DisGeNET, MONDO
│   │   └── processes.md       # GO, WikiPathways
│   │
│   ├── literature/             # Shared: Scientific publications
│   │   ├── _index.md
│   │   ├── sources.md         # PubMed, PMC, OpenAlex
│   │   ├── public-sources.md  # Open access sources
│   │   ├── pipeline-design.md # Literature processing pipeline
│   │   ├── data-structures.md # Document schemas
│   │   ├── coverage-analysis.md
│   │   └── abstracts-vs-fulltext.md
│   │
│   └── compounds/              # Shared: Chemical compounds
│       ├── _index.md
│       ├── natural-products.md # COCONUT, LOTUS
│       ├── pharmaceuticals.md  # ChEMBL, DrugBank
│       └── drug-metabolism.md  # ADMET, metabolism
│
├── domains/                     # HEALTH APPLICATION VIEWS
│   ├── _index.md               # Domain overview (13 domain files)
│   │
│   │  # Primary Health Domains
│   ├── mental-cognitive.md     # Psychiatric, neurodevelopmental (35+ DBs)
│   ├── cardio-metabolic.md     # Heart, diabetes, metabolic (~88GB)
│   ├── cancer-oncology.md      # COSMIC 38M+ mutations, TCGA 2.5PB
│   ├── autoimmune.md           # HLA 43K+ alleles, inflammation
│   ├── rare.md                 # Orphanet 6,500+ diseases, HPO, MONDO
│   ├── womens-pediatric.md     # Reproductive, developmental
│   ├── microbiome.md           # Gut health, microbiota
│   ├── allergy-pain.md         # Allergies, pain syndromes
│   ├── sleep-longevity-nutri.md # Sleep, aging, nutrigenomics
│   ├── oral-skin-sensory.md    # Dental, skin, vision, hearing
│   │
│   │  # Additional Domains
│   ├── clinical-biomarkers-labs.md    # LOINC, MarkerDB (~22GB)
│   ├── clinical-environmental-mito.md # Environmental, mitochondrial
│   └── community-patient-networks.md  # Patient networks
│
├── operations/                  # DATA MANAGEMENT & GOVERNANCE
│   ├── _index.md               # Operations overview (70 files total)
│   │
│   ├── downloads/              # Data acquisition procedures (5 files)
│   │   ├── _index.md
│   │   ├── processing-pipeline.md  # Master ETL pipeline
│   │   ├── traditional-medicine.md # Traditional medicine downloads
│   │   ├── pathways-targets.md     # Pathway data acquisition
│   │   ├── pharmaceuticals.md      # Drug data downloads
│   │   └── wikidata-bulk.md        # Wikidata/Wikipedia/DBpedia
│   │
│   ├── governance/             # Quality, compliance, stewardship (2 files)
│   │   ├── _index.md
│   │   ├── curation-framework.md   # 240 DBs evaluated, 126 retained
│   │   └── data-access-legal.md    # Legal risk assessment
│   │
│   ├── integration/            # ETL and harmonization (11 files)
│   │   ├── _index.md
│   │   ├── integration-guide.md    # Cross-reference strategy
│   │   ├── xrefs.md                # ID mapping (286 DBs via UniProt)
│   │   ├── compound-pathway-linking.md
│   │   ├── pathway-target-mapping.md
│   │   ├── size-estimates.md
│   │   ├── alt-sources.md
│   │   ├── wikidata-master-reference.md
│   │   ├── wikidata-pharma.md
│   │   ├── wikidata-supplements.md
│   │   ├── wikidata-traditional.md
│   │   └── wikipedia-wikidata.md
│   │
│   └── schemas/                # Data models and validation (47 files)
│       ├── _index.md
│       ├── ruvector-schema.md              # Master architecture
│       ├── unified-schema-analysis.md       # 5 core entity types
│       ├── schemas-index.md
│       ├── schemas-navigation.md
│       │
│       │  # Individual Database Schemas (43 files)
│       ├── dbsnp-schema.md
│       ├── clinvar-schema.md
│       ├── gnomad-schema.md
│       ├── pharmgkb-schema.md
│       ├── batman-tcm-schema.md
│       ├── imppat-schema.md
│       ├── kampodb-schema.md
│       ├── foodb.md
│       ├── usda-fooddata-central.md
│       ├── reactome-schema.md
│       ├── disgenet-schema.md
│       ├── coconut-schema.md
│       ├── lotus-schema.md
│       ├── chembl-schema.md
│       └── ... (28 more schema files)
│
└── research/                    # RESEARCH PRIORITIES (4 files)
    ├── _index.md
    ├── interventions-priority.md   # Intervention research priorities
    ├── literature-priority.md      # Literature review priorities
    ├── schema-gaps.md              # Database schema gap analysis
    └── genetics-schema-research.md # Genetics schema research
```

---

## Domain Views vs Primary Data

An important distinction in this catalog:

| Aspect | Primary Data (`databases/`) | Domain Views (`domains/`) |
|--------|----------------------------|---------------------------|
| **Purpose** | Store raw data source documentation | Aggregate by health application |
| **Organization** | By data type (genetics, nutrition, etc.) | By medical condition/specialty |
| **Content** | Individual database details | Cross-references to multiple databases |
| **Updates** | When source databases change | When clinical applications evolve |
| **Example** | `databases/genetics/primary.md` | `domains/cardio-metabolic.md` |

**Key Insight:** Domain files are "view" documents that reference databases from multiple categories by health application area. For example, `domains/mental-cognitive.md` references:
- Genetics: PGC psychiatric GWAS, neuroimaging genetics
- Traditional: Adaptogenic herbs, nootropic compounds
- Nutrition: Omega-3s, B vitamins, gut-brain nutrients

---

## Cross-Domain Considerations

Many biological processes span domains:

| Process | Affected Domains |
|---------|------------------|
| **Inflammation** | Cardiovascular, Autoimmune, Metabolic |
| **Methylation** | Mental Health, Cancer, Aging |
| **Oxidative Stress** | Cardiovascular, Neurodegenerative, Cancer |
| **Gut Microbiome** | Mental Health, Autoimmune, Metabolic |

Users should consult multiple domain views when researching complex conditions with multi-system involvement.

---

## Documentation Standards

Every content file in this catalog contains 7 standardized sections:

### Required Sections

| Section | Purpose | Example Content |
|---------|---------|-----------------|
| **Download** | Data acquisition instructions | URLs, API commands, FTP paths |
| **Schema** | Data model definition | Core fields, relationships, cardinality |
| **Glossary** | Term definitions | Domain terms, acronyms |
| **Sample Data** | Example records | JSON examples, query results |
| **License** | Usage rights | License type, commercial use, attribution |
| **Data Set Size** | Storage estimates | Record counts, file sizes |
| **Data Format** | File specifications | Primary format, compression, encoding |

### Frontmatter Schema

All files use YAML frontmatter:

```yaml
---
id: unique-identifier
title: "Human-Readable Title"
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [category, subcategory, keywords]
tier: 1 | 2 | 3             # Implementation priority
category: genetics | traditional | nutrition | shared
subcategory: tcm | ayurveda  # Traditional only
---
```

---

## Data Integration Strategy

### Hub Identifiers

Cross-category integration uses standardized hub identifiers:

| Entity Type | Hub Identifier | Coverage |
|-------------|----------------|----------|
| Proteins | UniProt ID | 286 databases |
| Genes | NCBI Gene ID | Standard reference |
| Compounds | InChIKey | Chemical uniqueness |
| Diseases | MONDO ID | Unified disease ontology |
| Pathways | Reactome ID | Pathway reference |

### Integration Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                    INTEGRATION PIPELINE                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  GENETICS              TRADITIONAL           NUTRITION          │
│                        MEDICINE                                  │
│     │                     │                     │               │
│     │                     │                     │               │
│     ▼                     ▼                     ▼               │
│  ┌─────┐              ┌─────┐              ┌─────┐             │
│  │Gene │              │Herb │              │Food │             │
│  │ID   │              │     │              │     │             │
│  └──┬──┘              └──┬──┘              └──┬──┘             │
│     │                    │                    │                 │
│     │    ┌───────────────┼────────────────────┘                │
│     │    │               │                                      │
│     ▼    ▼               ▼                                      │
│  ┌────────────────────────────┐                                 │
│  │       COMPOUND HUB         │                                 │
│  │      (InChIKey/PubChem)    │                                 │
│  └─────────────┬──────────────┘                                 │
│                │                                                 │
│                ▼                                                 │
│  ┌────────────────────────────┐                                 │
│  │       TARGET GENES         │                                 │
│  │      (UniProt/NCBI)        │                                 │
│  └─────────────┬──────────────┘                                 │
│                │                                                 │
│                ▼                                                 │
│  ┌────────────────────────────┐                                 │
│  │        PATHWAYS            │                                 │
│  │    (Reactome/KEGG)         │                                 │
│  └─────────────┬──────────────┘                                 │
│                │                                                 │
│                ▼                                                 │
│  ┌────────────────────────────┐                                 │
│  │        DISEASES            │                                 │
│  │      (MONDO/DisGeNET)      │                                 │
│  └────────────────────────────┘                                 │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Storage Requirements

### By Scope

| Scope | Download | Working | Processed | Total |
|-------|----------|---------|-----------|-------|
| **MVP** | 20 GB | 40 GB | 3 GB | ~65 GB |
| **Standard** | 50 GB | 100 GB | 12 GB | ~165 GB |
| **Comprehensive** | 200 GB | 400 GB | 30 GB | ~650 GB |
| **Full (with gnomAD)** | 700 GB | 1.4 TB | 120 GB | ~2.2 TB |

### By Category

| Category | Key Sources | Estimated Size |
|----------|-------------|----------------|
| Genetics | gnomAD, dbSNP, ClinVar | 500 GB - 1.5 TB |
| Traditional | BATMAN-TCM, IMPPAT, HERB | 5 - 10 GB |
| Nutrition | FooDB, USDA, HMDB | 10 - 20 GB |
| Shared | PubMed, Reactome, COCONUT | 100 - 500 GB |

---

## Governance Summary

### Curation Results

| Category | Count |
|----------|-------|
| Databases Evaluated | 240 |
| Databases Retained | 126 |
| Databases Pruned | 114 |

### Legal Risk Levels

| Risk | Approach | Examples |
|------|----------|----------|
| **LOW** | Direct API/partnerships | ClinicalTrials.gov, AACT, PubMed |
| **MEDIUM** | Restricted access with consent | Wearable APIs, Health Connect |
| **HIGH** | Avoid or partnership-only | Reddit, ConsumerLab, HealthUnlocked |

---

## Quick Reference

### File Counts by Directory

| Directory | Files | Purpose |
|-----------|-------|---------|
| `databases/` | 26 | Primary data sources |
| `domains/` | 14 | Health application views |
| `operations/downloads/` | 6 | Acquisition guides |
| `operations/governance/` | 3 | Compliance docs |
| `operations/integration/` | 12 | ETL processes |
| `operations/schemas/` | 48 | Data models |
| `research/` | 5 | Research priorities |
| **TOTAL** | **116** | |

### Category Assignment Quick Guide

```
category: genetics    → databases/genetics/
category: traditional → databases/traditional/
category: nutrition   → databases/nutrition/
category: shared      → databases/pathways/, literature/, compounds/
```

### Implementation Tiers

| Tier | Priority | Timeline |
|------|----------|----------|
| **Tier 1** | MVP Core | Phase 1 |
| **Tier 2** | Post-MVP | Phase 2-3 |
| **Tier 3** | Future | As needed |

---

## Navigation

- **Main Index:** [_index.md](./_index.md)
- **Databases:** [databases/_index.md](./databases/_index.md)
- **Domains:** [domains/_index.md](./domains/_index.md)
- **Operations:** [operations/_index.md](./operations/_index.md)
- **Research:** [research/_index.md](./research/_index.md)

---

## Maintenance

This README and the catalog documentation are updated when:
- New databases are added to primary storage
- Research reveals new domain-specific applications
- Clinical guidelines change for domain-specific interventions
- Cross-domain relationships are identified
- Storage requirements or access methods change

**Last Comprehensive Review:** January 2026

---

*Generated from swarm research findings. Data Team Contact: data-team@geneplatform.org*
