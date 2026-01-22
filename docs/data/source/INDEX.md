# Data Sources Index

**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 2.0

---

## TL;DR

Navigation index for detailed data source documentation. **27 specialized documents** (20 Final, 7 Pending) organized by domain: genetics, traditional medicine, nutrition, pathways, compounds, literature, diseases, and integration. Covers **256+ databases** total.

---

## Directory Structure

```
docs/data/source/
├── INDEX.md                 # This file
├── genetics/                # World 1: Modern Genetics
├── traditional/             # World 2: Traditional Medicine
├── pathways/                # Biological pathways & mechanisms
├── compounds/               # Pharmaceuticals & natural products
├── literature/              # Research papers & evidence
├── diseases/                # Health domain applications
├── clinical/                # Biomarkers & labs
├── integration/             # Cross-references & APIs
├── schemas/                 # Database schemas
└── research/                # Deep-dive research
```

---

## World 1: Modern Genetics

| Document | Content | Status |
|----------|---------|--------|
| [genetics/PRIMARY.md](./genetics/PRIMARY.md) | Population genetics, functional annotation, epigenetics, structural variants (30 databases) | Final |

---

## World 2: Traditional Medicine

| Document | Content | Status |
|----------|---------|--------|
| [traditional/TCM.md](./traditional/TCM.md) | BATMAN-TCM, HERB, ETCM, TCMID, SymMap (20 databases) | Final |
| [traditional/AYURVEDA.md](./traditional/AYURVEDA.md) | IMPPAT, OSADHI, GRAYU, TKDL (16 databases) | Final |
| [traditional/KAMPO.md](./traditional/KAMPO.md) | KampoDB, TM-MC, STORK, WAKANYAKU (15 databases) | Final |
| [traditional/WESTERN-HERBAL.md](./traditional/WESTERN-HERBAL.md) | Dr. Duke's, DSLD, ODS, EMA Monographs (27 databases) | Final |
| [traditional/GLOBAL.md](./traditional/GLOBAL.md) | African, Latin American medicine (15 databases) | Final |

---

## Pathways & Mechanisms

| Document | Content | Status |
|----------|---------|--------|
| [pathways/PRIMARY.md](./pathways/PRIMARY.md) | Reactome, WikiPathways, KEGG | Final |
| [pathways/DISEASE.md](./pathways/DISEASE.md) | DisGeNET, OMIM, HPO (12 databases) | Final |

---

## Compounds & Interventions

| Document | Content | Status |
|----------|---------|--------|
| [compounds/PHARMACEUTICALS.md](./compounds/PHARMACEUTICALS.md) | PharmGKB, CPIC, DrugBank, ChEMBL (22 databases) | Final |
| [compounds/NATURAL-PRODUCTS.md](./compounds/NATURAL-PRODUCTS.md) | COCONUT, LOTUS, NuBBE, NPASS | Final |
| [compounds/DRUG-METABOLISM.md](./compounds/DRUG-METABOLISM.md) | CYP450, drug interactions, transporters (21 databases) | Final |

---

## Literature & Evidence

| Document | Content | Status |
|----------|---------|--------|
| [literature/SOURCES.md](./literature/SOURCES.md) | PubMed, PMC, OpenAlex, Europe PMC | Final |

---

## Health Domains (Diseases)

| Document | Content | Status |
|----------|---------|--------|
| [diseases/MENTAL-COGNITIVE.md](./diseases/MENTAL-COGNITIVE.md) | Depression, anxiety, ADHD, cognition | Final |
| [diseases/CARDIO-METABOLIC.md](./diseases/CARDIO-METABOLIC.md) | Heart, diabetes, metabolic syndrome (10 databases) | Final |
| [diseases/CANCER-ONCOLOGY.md](./diseases/CANCER-ONCOLOGY.md) | COSMIC, cBioPortal (8 databases) | Final |
| [diseases/AUTOIMMUNE.md](./diseases/AUTOIMMUNE.md) | Autoimmune, hormones, thyroid (23 databases) | Final |
| [diseases/RARE.md](./diseases/RARE.md) | Orphanet, OMIM (9 databases) | Final |
| [diseases/WOMENS-PEDIATRIC.md](./diseases/WOMENS-PEDIATRIC.md) | Reproductive, pediatric health | Final |
| [diseases/MICROBIOME.md](./diseases/MICROBIOME.md) | HMP, GMrepo (12 databases) | Final |
| [diseases/ALLERGY-PAIN.md](./diseases/ALLERGY-PAIN.md) | Allergy, histamine, mast cell, pain (25 databases) | Final |
| [diseases/SLEEP-LONGEVITY-NUTRI.md](./diseases/SLEEP-LONGEVITY-NUTRI.md) | Circadian, aging, nutrigenomics (30 databases) | Final |

---

## Clinical & Biomarkers

| Document | Content | Status |
|----------|---------|--------|
| [clinical/BIOMARKERS-LABS.md](./clinical/BIOMARKERS-LABS.md) | Lab standards, biomarker databases (21 databases) | Final |
| [clinical/ENVIRONMENTAL-MITO.md](./clinical/ENVIRONMENTAL-MITO.md) | Toxicogenomics, mitochondrial (14 databases) | Final |

---

## Integration & Technical

| Document | Content | Status |
|----------|---------|--------|
| [integration/WIKIPEDIA-WIKIDATA.md](./integration/WIKIPEDIA-WIKIDATA.md) | Semantic web, Wikidata extraction (12+ resources) | Final |
| [integration/ALT-SOURCES.md](./integration/ALT-SOURCES.md) | Alternative access methods (15 sources) | Final |
| [integration/XREFS.md](./integration/XREFS.md) | UniProt ID mapping, Gene Ontology, cross-references (6 databases) | Final |

---

## Schemas

Database schema documentation: [schemas/INDEX.md](./schemas/INDEX.md)

45 detailed schema files covering CLINVAR, DBSNP, CHEBI, UNIPROT, Reactome, and more.

---

## Research

Deep-dive analysis: [research/](./research/)

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 2.0 | January 2026 | Data Engineering | Reorganized into domain subdirectories, removed 43-XX prefixes |
| 1.0 | January 2026 | Data Engineering | Initial index with migration mapping |
