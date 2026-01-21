# Data Sources Index

**Document ID:** 43-00-INDEX
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-DATA-SOURCES.md](../43-DATA-SOURCES.md)

---

## TL;DR

Navigation index for detailed data source documentation. **27 specialized documents** (20 Final, 7 Pending) organized by THREE WORLDS (genetics, traditional medicine, nutrition) plus cross-cutting concerns (pathways, interventions, literature, health domains, integration). Covers **256+ databases** total. Each document follows the 43-XX numbering convention where the second digit indicates category.

---

## Numbering Convention

| Pattern | Category | Example |
|---------|----------|---------|
| `43-1x` | World 1: Modern Genetics | `43-11-GENETICS-PRIMARY.md` |
| `43-2x` | World 2: Traditional Medicine | `43-21-TCM.md` |
| `43-3x` | World 3: Nutritional Science | `43-31-FOODS.md` |
| `43-4x` | Pathways & Mechanisms | `43-41-PATHWAYS-PRIMARY.md` |
| `43-5x` | Interventions | `43-51-PHARMACEUTICALS.md` |
| `43-6x` | Literature & Evidence | `43-61-SOURCES.md` |
| `43-7x` | Health Domains | `43-71-MENTAL-COGNITIVE.md` |
| `43-8x` | Integration & Technical | `43-81-APIS.md` |

---

## World 1: Modern Genetics (43-1x)

| Doc | Title | Content | Status |
|-----|-------|---------|--------|
| 43-10 | [GENETICS-OVERVIEW](./43-10-GENETICS-OVERVIEW.md) | Section summary | Pending |
| 43-11 | [GENETICS-PRIMARY](./43-11-GENETICS-PRIMARY.md) | Population genetics, functional annotation, epigenetics, structural variants (30 databases) | Final |
| 43-12 | [GENETICS-POPULATION](./43-12-GENETICS-POPULATION.md) | gnomAD, TOPMed, All of Us, ALFA | Pending |
| 43-13 | [GENETICS-FUNCTIONAL](./43-13-GENETICS-FUNCTIONAL.md) | AlphaMissense, dbNSFP, SpliceAI, CADD | Pending |
| 43-14 | [GENETICS-EPIGENETICS](./43-14-GENETICS-EPIGENETICS.md) | ENCODE, IHEC, MethBank, FANTOM5 | Pending |
| 43-15 | [GENETICS-STRUCTURAL](./43-15-GENETICS-STRUCTURAL.md) | gnomAD-SV, dbVar, DECIPHER | Pending |

---

## World 2: Traditional Medicine (43-2x)

| Doc | Title | Content | Status |
|-----|-------|---------|--------|
| 43-20 | [TRADITIONAL-OVERVIEW](./43-20-TRADITIONAL-OVERVIEW.md) | Section summary | Pending |
| 43-21 | [TCM](./43-21-TCM.md) | BATMAN-TCM, HERB, ETCM, TCMID, SymMap (20 databases) | Final |
| 43-22 | [AYURVEDA](./43-22-AYURVEDA.md) | IMPPAT, OSADHI, GRAYU, TKDL (16 databases) | Final |
| 43-23 | [KAMPO](./43-23-KAMPO.md) | KampoDB, TM-MC, STORK, WAKANYAKU (15 databases) | Final |
| 43-24 | [WESTERN-HERBAL](./43-24-WESTERN-HERBAL.md) | Dr. Duke's, DSLD, ODS, EMA Monographs (27 databases) | Final |
| 43-25 | [GLOBAL-TRADITIONAL](./43-25-GLOBAL-TRADITIONAL.md) | African, Latin American medicine (15 databases) | Final |

---

## World 3: Nutritional Science (43-3x)

| Doc | Title | Content | Status |
|-----|-------|---------|--------|
| 43-30 | [NUTRITION-OVERVIEW](./43-30-NUTRITION-OVERVIEW.md) | Section summary | Pending |
| 43-31 | [FOODS](./43-31-FOODS.md) | USDA FoodData, FooDB, Open Food Facts | Pending |
| 43-32 | [SUPPLEMENTS](./43-32-SUPPLEMENTS.md) | DSLD, NIH ODS, Natural Medicines | Pending |
| 43-33 | [BIOACTIVE-COMPOUNDS](./43-33-BIOACTIVE-COMPOUNDS.md) | Phenol-Explorer, Flavonoid DB | Pending |

---

## Cross-Cutting: Pathways & Mechanisms (43-4x)

| Doc | Title | Content | Status |
|-----|-------|---------|--------|
| 43-40 | [PATHWAYS-OVERVIEW](./43-40-PATHWAYS-OVERVIEW.md) | Section summary | Pending |
| 43-41 | [PATHWAYS-PRIMARY](./43-41-PATHWAYS-PRIMARY.md) | Reactome, WikiPathways, KEGG | Final |
| 43-42 | [PATHWAYS-METABOLISM](./43-42-PATHWAYS-METABOLISM.md) | KEGG, MetaCyc, SMPDB | Pending |
| 43-43 | [PATHWAYS-DISEASE](./43-43-PATHWAYS-DISEASE.md) | DisGeNET, OMIM, HPO (12 databases) | Final |

---

## Cross-Cutting: Interventions (43-5x)

| Doc | Title | Content | Status |
|-----|-------|---------|--------|
| 43-50 | [INTERVENTIONS-OVERVIEW](./43-50-INTERVENTIONS-OVERVIEW.md) | Section summary, storage estimates | Pending |
| 43-51 | [PHARMACEUTICALS](./43-51-PHARMACEUTICALS.md) | PharmGKB, CPIC, DrugBank, ChEMBL (22 databases) | Final |
| 43-52 | [NATURAL-PRODUCTS](./43-52-NATURAL-PRODUCTS.md) | COCONUT, LOTUS, NuBBE, NPASS | Final |
| 43-53 | [DRUG-METABOLISM](./43-53-DRUG-METABOLISM.md) | CYP450, drug interactions, transporters (21 databases) | Final |

---

## Cross-Cutting: Literature & Evidence (43-6x)

| Doc | Title | Content | Status |
|-----|-------|---------|--------|
| 43-60 | [LITERATURE-OVERVIEW](./43-60-LITERATURE-OVERVIEW.md) | Section summary | Pending |
| 43-61 | [LITERATURE-SOURCES](./43-61-LITERATURE-SOURCES.md) | PubMed, PMC, OpenAlex, Europe PMC | Final |
| 43-62 | [LITERATURE-PIPELINE](./43-62-LITERATURE-PIPELINE.md) | Processing architecture and strategy | Pending |

---

## Health Domain Applications (43-7x)

| Doc | Title | Content | Status |
|-----|-------|---------|--------|
| 43-70 | [DOMAINS-OVERVIEW](./43-70-DOMAINS-OVERVIEW.md) | Section summary | Pending |
| 43-71 | [MENTAL-COGNITIVE](./43-71-MENTAL-COGNITIVE.md) | Depression, anxiety, ADHD, cognition | Final |
| 43-72 | [CARDIO-METABOLIC](./43-72-CARDIO-METABOLIC.md) | Heart, diabetes, metabolic syndrome (10 databases) | Final |
| 43-73 | [CANCER-ONCOLOGY](./43-73-CANCER-ONCOLOGY.md) | COSMIC, cBioPortal (8 databases) | Final |
| 43-74 | [AUTOIMMUNE](./43-74-AUTOIMMUNE.md) | Autoimmune, hormones, thyroid (23 databases) | Final |
| 43-75 | [RARE-DISEASES](./43-75-RARE-DISEASES.md) | Orphanet, OMIM (9 databases) | Final |
| 43-76 | [WOMENS-PEDIATRIC](./43-76-WOMENS-PEDIATRIC.md) | Reproductive, pediatric health | Final |
| 43-77 | [MICROBIOME](./43-77-MICROBIOME.md) | HMP, GMrepo (12 databases) | Final |
| 43-78 | [ALLERGY-PAIN](./43-78-ALLERGY-PAIN.md) | Allergy, histamine, mast cell, pain (25 databases) | Final |
| 43-79 | [SLEEP-LONGEVITY-NUTRI](./43-79-SLEEP-LONGEVITY-NUTRI.md) | Circadian, aging, nutrigenomics (30 databases) | Final |

---

## Integration & Technical (43-8x)

| Doc | Title | Content | Status |
|-----|-------|---------|--------|
| 43-80 | [INTEGRATION-OVERVIEW](./43-80-INTEGRATION-OVERVIEW.md) | Section summary | Pending |
| 43-81 | [BIOMARKERS-LABS](./43-81-BIOMARKERS-LABS.md) | Lab standards, biomarker databases (21 databases) | Final |
| 43-82 | [STORAGE](./43-82-STORAGE.md) | Size estimates, formats, optimization | Pending |
| 43-83 | [ENVIRONMENTAL-MITO](./43-83-ENVIRONMENTAL-MITO.md) | Toxicogenomics, mitochondrial (14 databases) | Final |
| 43-84 | [WIKIPEDIA-WIKIDATA](./43-84-WIKIPEDIA-WIKIDATA.md) | Semantic web, Wikidata extraction (12+ resources) | Final |
| 43-85 | [ALT-DATA-SOURCES](./43-85-ALT-DATA-SOURCES.md) | Alternative access methods (15 sources) | Final |
| 43-86 | [INTEGRATION-XREFS](./43-86-INTEGRATION-XREFS.md) | UniProt ID mapping, Gene Ontology, cross-reference infrastructure (6 databases) | Final |

---

## MECE Compliance

| Category | Scope | Does NOT Include |
|----------|-------|------------------|
| Genetics (43-1x) | DNA variants, frequencies, annotations, epigenetics | Pathway mechanisms |
| Traditional Medicine (43-2x) | Historical healing systems by tradition | Pharmaceutical drugs |
| Nutritional Science (43-3x) | Foods, supplements, bioactive compounds | Traditional formulas |
| Pathways (43-4x) | Biological mechanisms, reactions, diseases | Raw database catalogs |
| Interventions (43-5x) | Therapeutic agents, targets across modalities | Disease databases |
| Literature (43-6x) | Research papers, evidence, processing | Database sources |
| Domains (43-7x) | Condition-specific applications | General databases |
| Integration (43-8x) | Technical implementation, APIs, storage | Scientific content |

---

## Migration Source Mapping

| research.old File | Target Document |
|-------------------|-----------------|
| data-sources-genetics-expanded.md | 43-11 through 43-15 |
| data-sources-traditional-medicine-expanded.md | 43-21 through 43-25 |
| interventions-tcm.md | 43-21-TCM.md |
| interventions-ayurveda.md | 43-22-AYURVEDA.md |
| interventions-kampo.md | 43-23-KAMPO.md |
| interventions-western-herbal.md | 43-24-WESTERN-HERBAL.md |
| interventions-pharma.md | 43-51-PHARMACEUTICALS.md |
| interventions-natural-products.md | 43-52-NATURAL-PRODUCTS.md |
| pathways-databases.md | 43-41, 43-42 |
| biological-processes-databases.md | 43-40, 43-43 |
| papers-*.md | 43-61, 43-62 |
| data-sources-mental-cognitive.md | 43-71-MENTAL-COGNITIVE.md |
| data-sources-cardio-metabolic.md | 43-72-CARDIO-METABOLIC.md |
| data-sources-cancer-oncology.md | 43-73-CANCER-ONCOLOGY.md |
| data-sources-hormones-autoimmune.md | 43-74-AUTOIMMUNE.md |
| data-sources-rare-diseases.md | 43-75-RARE-DISEASES.md |
| data-sources-womens-pediatric.md | 43-76-WOMENS-PEDIATRIC.md |
| data-sources-microbiome.md | 43-77-MICROBIOME.md |
| data-size-estimates.md | 43-82-STORAGE.md |
| interventions-recommendation.md | 43-50-INTERVENTIONS-OVERVIEW.md |
| compound-pathway-linking.md | 43-53-TARGETS-LINKING.md |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial index with migration mapping |
