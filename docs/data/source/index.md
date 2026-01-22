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
├── index.md                 # This file
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
| [genetics/primary.md](./genetics/primary.md) | Population genetics, functional annotation, epigenetics, structural variants (30 databases) | Final |

---

## World 2: Traditional Medicine

| Document | Content | Status |
|----------|---------|--------|
| [traditional/tcm.md](./traditional/tcm.md) | BATMAN-TCM, HERB, ETCM, TCMID, SymMap (20 databases) | Final |
| [traditional/ayurveda.md](./traditional/ayurveda.md) | IMPPAT, OSADHI, GRAYU, TKDL (16 databases) | Final |
| [traditional/kampo.md](./traditional/kampo.md) | KampoDB, TM-MC, STORK, WAKANYAKU (15 databases) | Final |
| [traditional/western-herbal.md](./traditional/western-herbal.md) | Dr. Duke's, DSLD, ODS, EMA Monographs (27 databases) | Final |
| [traditional/global.md](./traditional/global.md) | African, Latin American medicine (15 databases) | Final |

---

## Pathways & Mechanisms

| Document | Content | Status |
|----------|---------|--------|
| [pathways/primary.md](./pathways/primary.md) | Reactome, WikiPathways, KEGG | Final |
| [pathways/disease.md](./pathways/disease.md) | DisGeNET, OMIM, HPO (12 databases) | Final |

---

## Compounds & Interventions

| Document | Content | Status |
|----------|---------|--------|
| [compounds/pharmaceuticals.md](./compounds/pharmaceuticals.md) | PharmGKB, CPIC, DrugBank, ChEMBL (22 databases) | Final |
| [compounds/natural-products.md](./compounds/natural-products.md) | COCONUT, LOTUS, NuBBE, NPASS | Final |
| [compounds/drug-metabolism.md](./compounds/drug-metabolism.md) | CYP450, drug interactions, transporters (21 databases) | Final |

---

## Literature & Evidence

| Document | Content | Status |
|----------|---------|--------|
| [literature/sources.md](./literature/sources.md) | PubMed, PMC, OpenAlex, Europe PMC | Final |

---

## Health Domains (Diseases)

| Document | Content | Status |
|----------|---------|--------|
| [diseases/mental-cognitive.md](./diseases/mental-cognitive.md) | Depression, anxiety, ADHD, cognition | Final |
| [diseases/cardio-metabolic.md](./diseases/cardio-metabolic.md) | Heart, diabetes, metabolic syndrome (10 databases) | Final |
| [diseases/cancer-oncology.md](./diseases/cancer-oncology.md) | COSMIC, cBioPortal (8 databases) | Final |
| [diseases/autoimmune.md](./diseases/autoimmune.md) | Autoimmune, hormones, thyroid (23 databases) | Final |
| [diseases/rare.md](./diseases/rare.md) | Orphanet, OMIM (9 databases) | Final |
| [diseases/womens-pediatric.md](./diseases/womens-pediatric.md) | Reproductive, pediatric health | Final |
| [diseases/microbiome.md](./diseases/microbiome.md) | HMP, GMrepo (12 databases) | Final |
| [diseases/allergy-pain.md](./diseases/allergy-pain.md) | Allergy, histamine, mast cell, pain (25 databases) | Final |
| [diseases/sleep-longevity-nutri.md](./diseases/sleep-longevity-nutri.md) | Circadian, aging, nutrigenomics (30 databases) | Final |

---

## Clinical & Biomarkers

| Document | Content | Status |
|----------|---------|--------|
| [clinical/biomarkers-labs.md](./clinical/biomarkers-labs.md) | Lab standards, biomarker databases (21 databases) | Final |
| [clinical/environmental-mito.md](./clinical/environmental-mito.md) | Toxicogenomics, mitochondrial (14 databases) | Final |

---

## Integration & Technical

| Document | Content | Status |
|----------|---------|--------|
| [integration/wikipedia-wikidata.md](./integration/wikipedia-wikidata.md) | Semantic web, Wikidata extraction (12+ resources) | Final |
| [integration/alt-sources.md](./integration/alt-sources.md) | Alternative access methods (15 sources) | Final |
| [integration/xrefs.md](./integration/xrefs.md) | UniProt ID mapping, Gene Ontology, cross-references (6 databases) | Final |

---

## Schemas

Database schema documentation: [schemas/index.md](./schemas/index.md)

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
