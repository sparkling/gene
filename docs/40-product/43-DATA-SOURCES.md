# Data Sources — THREE WORLDS Database Inventory

**Document ID:** 43-DATA-SOURCES
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 2.0

---

## TL;DR

The platform integrates data from THREE WORLDS: Modern Genetics (1.2B+ SNPs), Traditional Medicine (100K+ compounds across TCM, Ayurveda, Kampo, Western Herbal), and Nutritional Science (380K+ foods). **250+ databases cataloged** across 26 specialized documents covering genetics, traditional medicine, interventions, health domains, pathways, literature, and technical integration. This document serves as the master index; see detailed files in [`data-sources/`](./data-sources/) for complete database catalogs.

---

## Quick Navigation — Detailed Documentation

### Genetics (World 1)

| Doc | Title | Databases | Size |
|-----|-------|-----------|------|
| [43-11](./data-sources/43-11-GENETICS-PRIMARY.md) | **Genetics Primary** | Population, functional annotation, epigenetics, structural variants, 30 total | 28KB |

### Traditional Medicine (World 2)

| Doc | Title | Databases | Size |
|-----|-------|-----------|------|
| [43-21](./data-sources/43-21-TCM.md) | **TCM** | BATMAN-TCM, HERB, SymMap, TCMBank, 20 total | 27KB |
| [43-22](./data-sources/43-22-AYURVEDA.md) | **Ayurveda** | IMPPAT, OSADHI, GRAYU, TKDL, 16 total | 22KB |
| [43-23](./data-sources/43-23-KAMPO.md) | **Kampo** | KampoDB, STORK, TM-MC, 15+ total | 22KB |
| [43-24](./data-sources/43-24-WESTERN-HERBAL.md) | **Western Herbal** | DSLD, Dr. Duke's, EMA, 27 total | 27KB |
| [43-25](./data-sources/43-25-GLOBAL-TRADITIONAL.md) | **Global Traditional** | African, Latin American medicine, 15 total | 18KB |

### Interventions (Cross-World)

| Doc | Title | Databases | Size |
|-----|-------|-----------|------|
| [43-51](./data-sources/43-51-PHARMACEUTICALS.md) | **Pharmaceuticals** | PharmGKB, CPIC, DrugBank, ChEMBL, 22 total | 27KB |
| [43-52](./data-sources/43-52-NATURAL-PRODUCTS.md) | **Natural Products** | COCONUT, LOTUS, NPASS, NPAtlas | 28KB |
| [43-53](./data-sources/43-53-DRUG-METABOLISM.md) | **Drug Metabolism** | CYP450, drug interactions, transporters, 21 total | 24KB |

### Health Domains

| Doc | Title | Focus | Size |
|-----|-------|-------|------|
| [43-71](./data-sources/43-71-MENTAL-COGNITIVE.md) | **Mental-Cognitive** | Psychiatry, neurology, cognition | 25KB |
| [43-72](./data-sources/43-72-CARDIO-METABOLIC.md) | **Cardio-Metabolic** | Heart, diabetes, metabolic | 17KB |
| [43-73](./data-sources/43-73-CANCER-ONCOLOGY.md) | **Cancer-Oncology** | COSMIC, cBioPortal, OncoKB | 15KB |
| [43-74](./data-sources/43-74-AUTOIMMUNE.md) | **Autoimmune** | Autoimmune, hormones | 21KB |
| [43-75](./data-sources/43-75-RARE-DISEASES.md) | **Rare Diseases** | Orphanet, OMIM, HPO | 18KB |
| [43-76](./data-sources/43-76-WOMENS-PEDIATRIC.md) | **Women's-Pediatric** | Reproductive, pediatric | 24KB |
| [43-77](./data-sources/43-77-MICROBIOME.md) | **Microbiome** | Gut, oral, skin microbiome | 15KB |
| [43-78](./data-sources/43-78-ALLERGY-PAIN.md) | **Allergy-Pain** | Allergy, histamine, mast cell, pain | 25KB |
| [43-79](./data-sources/43-79-SLEEP-LONGEVITY-NUTRI.md) | **Sleep-Longevity-Nutri** | Circadian, aging, nutrigenomics, 30 total | 22KB |

### Pathways & Literature

| Doc | Title | Content | Size |
|-----|-------|---------|------|
| [43-41](./data-sources/43-41-PATHWAYS-PRIMARY.md) | **Pathways Primary** | Reactome, WikiPathways, KEGG | 23KB |
| [43-43](./data-sources/43-43-PATHWAYS-DISEASE.md) | **Disease Pathways** | DisGeNET, OMIM, HPO | 16KB |
| [43-61](./data-sources/43-61-LITERATURE-SOURCES.md) | **Literature** | PubMed, PMC, OpenAlex | 19KB |

### Integration & Technical

| Doc | Title | Content | Size |
|-----|-------|---------|------|
| [43-81](./data-sources/43-81-BIOMARKERS-LABS.md) | **Biomarkers-Labs** | Lab standards, biomarker databases, 21 total | 22KB |
| [43-83](./data-sources/43-83-ENVIRONMENTAL-MITO.md) | **Environmental-Mito** | Toxicogenomics, mitochondrial, 14 total | 18KB |
| [43-84](./data-sources/43-84-WIKIPEDIA-WIKIDATA.md) | **Wikipedia-Wikidata** | Semantic web, Wikidata extraction | 20KB |
| [43-85](./data-sources/43-85-ALT-DATA-SOURCES.md) | **Alt Data Sources** | Alternative access methods, 15 sources | 16KB |

**Full Index:** [data-sources/43-00-INDEX.md](./data-sources/43-00-INDEX.md)

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary nutrition source | USDA FoodData Central | CC0 license, 380K foods, research-grade | Jan 2026 |
| Genetics priority | dbSNP + ClinVar + PharmGKB | Foundational, well-maintained, free | Jan 2026 |
| Traditional medicine approach | TCMSP + IMPPAT + KampoDB | Only comprehensive open sources available | Jan 2026 |
| No homeopathy | Excluded | Lacks evidence basis | Jan 2026 |

---

## THREE WORLDS Overview

| World | Description | Primary Sources | Record Count |
|-------|-------------|-----------------|--------------|
| **1. Modern Genetics** | Science of DNA, SNPs, pathways | dbSNP, ClinVar, PharmGKB, SNPedia | 1.2B+ SNPs |
| **2. Traditional Medicine** | Ancient healing wisdom | TCMSP, IMPPAT, KampoDB, HERB | 100K+ compounds |
| **3. Nutritional Science** | Food, nutrients, effects | USDA FoodData Central, FooDB | 380K+ foods |

---

## World 1: Modern Genetics

### Primary Sources

#### dbSNP (NCBI)
| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/snp/ |
| **Content** | Reference SNP database |
| **Records** | 1+ billion SNPs |
| **License** | Public Domain |
| **API** | E-utilities REST API |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |

#### ClinVar (NCBI)
| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | Clinical variant interpretations |
| **Records** | 2M+ submissions |
| **License** | Public Domain |
| **API** | E-utilities REST API |
| **Update Frequency** | Weekly |
| **Priority** | Tier 1 (MVP) |

#### PharmGKB
| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | Pharmacogenomics relationships |
| **Records** | 700+ drugs, 1000+ genes |
| **License** | CC BY-SA 4.0 |
| **API** | REST API (free registration) |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |

#### SNPedia
| Field | Value |
|-------|-------|
| **URL** | https://www.snpedia.com/ |
| **Content** | SNP interpretations (wiki-style) |
| **Records** | 100K+ SNPs documented |
| **License** | CC BY-NC-SA 3.0 |
| **API** | MediaWiki API |
| **Update Frequency** | Community-driven |
| **Priority** | Tier 1 (MVP) |

### Secondary Genetics Sources

| Source | Content | Records | License | Priority |
|--------|---------|---------|---------|----------|
| **Reactome** | Pathway database | 2,600+ pathways | CC0 | Tier 1 |
| **WikiPathways** | Curated pathways | 3,000+ pathways | CC0 | Tier 1 |
| **KEGG** | Metabolic pathways | 500+ pathways | Academic free | Tier 2 |
| **GeneCards** | Gene information | 20K+ genes | Non-commercial | Tier 2 |
| **OMIM** | Genetic disorders | 25K+ entries | Subscription | Tier 3 |

---

## World 2: Traditional Medicine

### TCM (Traditional Chinese Medicine)

#### TCMSP (TCM Systems Pharmacology)
| Field | Value |
|-------|-------|
| **URL** | https://old.tcmsp-e.com/tcmsp.php |
| **Content** | Chinese herbs, compounds, targets |
| **Records** | 500+ herbs, 30K+ compounds |
| **License** | Academic/Research use |
| **API** | Web scraping required |
| **Update Frequency** | Static (2014) |
| **Priority** | Tier 2 |

#### HERB Database
| Field | Value |
|-------|-------|
| **URL** | http://herb.ac.cn/ |
| **Content** | High-throughput herb-disease associations |
| **Records** | 12K+ targets, 7K+ diseases |
| **License** | Academic use |
| **Priority** | Tier 2 |

### Ayurveda

#### IMPPAT (Indian Medicinal Plants)
| Field | Value |
|-------|-------|
| **URL** | https://cb.imsc.res.in/imppat/ |
| **Content** | Indian medicinal plants, phytochemicals |
| **Records** | 1,700+ plants, 10K+ phytochemicals |
| **License** | Academic use |
| **API** | Download available |
| **Priority** | Tier 2 |

### Kampo (Japanese Herbal Medicine)

#### KampoDB
| Field | Value |
|-------|-------|
| **URL** | http://wakanmoview.inm.u-toyama.ac.jp/kampo/ |
| **Content** | Kampo formulas and ingredients |
| **Records** | 50+ major formulas |
| **License** | Academic use |
| **Priority** | Tier 2 |

### Western Herbal

| Source | Content | Records | License | Priority |
|--------|---------|---------|---------|----------|
| **Dr. Duke's Phytochemical** | Plant compounds + activities | 15K+ compounds | Public | Tier 2 |
| **NAPRALERT** | Natural products | 200K+ records | Subscription | Tier 3 |

---

## World 3: Nutritional Science

### Primary Source: USDA FoodData Central

| Field | Value |
|-------|-------|
| **URL** | https://fdc.nal.usda.gov/ |
| **Content** | Comprehensive food nutrient profiles |
| **Records** | 380,000+ foods |
| **License** | CC0 1.0 (Public Domain) |
| **API** | REST API (free, 1,000 requests/hour) |
| **Documentation** | https://fdc.nal.usda.gov/api-guide/ |
| **Download** | https://fdc.nal.usda.gov/download-datasets/ |
| **Update Frequency** | Branded Foods: monthly; Full: every 6 months |
| **Priority** | Tier 1 (MVP) |

**Data Types Available:**
- Foundation Foods (research-grade reference data)
- SR Legacy (historical USDA data)
- FNDDS (food consumption surveys)
- Branded Foods (commercial products)
- Experimental Foods

**Why Primary:**
- Government-verified data from official USDA testing
- Free with no usage limits
- Research-grade accuracy including detailed micronutrients
- Covers US foods comprehensively

### Secondary Source: FooDB

| Field | Value |
|-------|-------|
| **URL** | https://foodb.ca/ |
| **Content** | Chemical composition of foods |
| **Records** | Detailed compound data |
| **License** | Free non-commercial; commercial requires permission |
| **API** | Beta API (https://foodb.ca/api_doc) |
| **Download** | https://foodb.ca/downloads |
| **Priority** | Tier 1 (MVP) |

**Unique Value:**
- Deeper chemical composition than USDA
- Links compounds to biological targets
- Good for mechanistic understanding (how foods affect pathways)

### Secondary Source: Open Food Facts

| Field | Value |
|-------|-------|
| **URL** | https://world.openfoodfacts.org/ |
| **Content** | Crowdsourced product database |
| **Records** | 2.8M+ products from 150+ countries |
| **License** | Open Database License (ODbL) |
| **API** | REST API (free, open source) |
| **Download** | Daily data dumps |
| **Priority** | Tier 2 |

**Unique Value:**
- Best for international/branded product coverage
- Community-curated, constantly updated
- Includes environmental impact scores

### Tertiary Nutrient Sources

| Source | Content | Use Case | Priority |
|--------|---------|----------|----------|
| **DSLD (NIH)** | 140K+ supplement labels | Supplement ingredients | Tier 2 |
| **NIH ODS** | Vitamin/mineral fact sheets | Evidence summaries | Tier 2 |
| **Phenol-Explorer** | Polyphenol content | Antioxidant analysis | Tier 3 |
| **USDA Flavonoid Database** | Flavonoid content | Specific nutrient queries | Tier 3 |

---

## Knowledge Store Domains

The platform integrates knowledge across 10 domains:

| Domain | Content | Sources |
|--------|---------|---------|
| **Biochemical Pathways** | Methylation, detox, neurotransmitters, energy, inflammation, BH4 | Reactome, WikiPathways, KEGG |
| **Genetics/SNPs** | Gene functions, variant impacts, research | dbSNP, ClinVar, SNPedia |
| **Pharmaceuticals** | Drug mechanisms, interactions, pharmacogenomics | PharmGKB, DrugBank |
| **Nutritional Biochemistry** | Vitamins, minerals, amino acids, cofactors | USDA, FooDB |
| **TCM** | Patterns, herbs, formulas, meridians, diagnostics | TCMSP, HERB |
| **Kampo** | Japanese herbal formulas | KampoDB |
| **Ayurveda** | Doshas, herbs, lifestyle recommendations | IMPPAT |
| **Western Herbal** | European/American botanical medicine | Dr. Duke's |
| **Folk Medicine** | Regional/indigenous practices | Literature review |
| **Research Literature** | PubMed, clinical trials across all modalities | PubMed API |

*Note: No homeopathy (lacks evidence basis)*

---

## User Data Integration

### Inputs the Platform Accepts

| Category | Examples | Phase |
|----------|----------|-------|
| **Genetic Data** | 23andMe, AncestryDNA, VCF files | MVP |
| **Lab Results** | Blood, hormones, stool, urine | MVP |
| **Symptoms** | Physical and mental health complaints | MVP |
| **Diet** | Food logs, patterns, restrictions | MVP |
| **Lifestyle** | Sleep, exercise, stress, environment | MVP |
| **Mental Health** | History, cognitive patterns, mood | MVP |
| **Diagnoses** | Conditions, family history | MVP |
| **Treatments Tried** | What's worked, what hasn't | MVP |
| **Microbiome** | Gut, oral, skin genetics | Phase 2+ |
| **Epigenetics** | Gene expression data | Phase 2+ |
| **Wearables** | Continuous health data | Phase 2+ |

---

## Mental Health to Biochemistry Connections

| Pathway/Gene | Mental Health Link |
|--------------|-------------------|
| Methylation | Schizophrenia, bipolar, depression |
| Dopamine/norepinephrine | ADHD, motivation, focus |
| Serotonin synthesis | Anxiety, depression, OCD |
| GABA/glutamate balance | Anxiety, excitotoxicity |
| MAOA/MAOB/COMT | Aggression, stress response, "warrior gene" |
| Inflammation | Depression, brain fog, fatigue |
| BH4 cycle | Neurotransmitter synthesis |
| Histamine | Anxiety, insomnia, brain inflammation |
| Glutathione | Detox, oxidative stress, brain health |

---

## Integration Priority Matrix

### Tier 1: MVP (Months 1-4)

| Source | World | Records | License | Status |
|--------|-------|---------|---------|--------|
| USDA FoodData Central | Nutrition | 380K foods | CC0 | Ready |
| FooDB | Nutrition | Compounds | Non-commercial | Ready |
| dbSNP | Genetics | 1B+ SNPs | Public Domain | Ready |
| ClinVar | Genetics | 2M+ variants | Public Domain | Ready |
| PharmGKB | Genetics | 700+ drugs | CC BY-SA 4.0 | Ready |
| SNPedia | Genetics | 100K+ SNPs | CC BY-NC-SA | Ready |
| Reactome | Pathways | 2,600+ pathways | CC0 | Ready |
| WikiPathways | Pathways | 3,000+ pathways | CC0 | Ready |

### Tier 2: Post-MVP (Months 5-8)

| Source | World | Records | License | Status |
|--------|-------|---------|---------|--------|
| TCMSP | TCM | 500+ herbs | Academic | Needs scraping |
| IMPPAT | Ayurveda | 1,700 plants | Academic | Ready |
| KampoDB | Kampo | 50+ formulas | Academic | Ready |
| HERB | TCM | 12K+ targets | Academic | Ready |
| Open Food Facts | Nutrition | 2.8M products | ODbL | Ready |
| DSLD | Supplements | 140K labels | Public | Ready |

### Tier 3: Future (Months 9+)

| Source | World | Notes |
|--------|-------|-------|
| KEGG | Pathways | Academic license review needed |
| DrugBank | Pharma | Commercial license needed |
| OMIM | Genetics | Subscription required |
| NAPRALERT | Herbal | Subscription required |

---

## Data Quality Assessment

| Source | Accuracy | Completeness | Currency | Reliability |
|--------|----------|--------------|----------|-------------|
| USDA FoodData | High | High | Updated regularly | Government-verified |
| dbSNP | High | High | Continuous updates | Gold standard |
| ClinVar | Medium-High | Growing | Weekly updates | Expert-curated |
| PharmGKB | High | Medium | Monthly | Expert-curated |
| TCMSP | Medium | Medium | Static (2014) | Research-derived |
| IMPPAT | Medium | Medium | Updated 2021 | Research-derived |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [44-ARCHITECTURE](./44-ARCHITECTURE.md) | Informs database design |
| [45-DATA-MODEL](./45-DATA-MODEL.md) | Informs entity structure |
| [42-ROADMAP](./42-ROADMAP.md) | Informed by integration priority |
| [data-sources/43-00-INDEX.md](./data-sources/43-00-INDEX.md) | Detailed navigation index |

---

## Open Questions

- [ ] Initial SNP panel — which 100-200 to start with for MVP?
- [ ] TCM/Kampo/Ayurveda data quality — needs validation
- [ ] KEGG license — is academic use sufficient?
- [ ] Data update pipeline — automated vs manual refresh?

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Complete THREE WORLDS inventory |
| 2.0 | January 2026 | Data Engineering | Added navigation to 17 detailed data-sources/ documents (150+ databases migrated from research.old/) |
