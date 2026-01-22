# Data Sources Index

**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 2.1

---

## TL;DR

Navigation index for detailed data source documentation. **45 specialized documents** organized by domain: genetics, traditional medicine, nutrition, pathways, compounds, literature, diseases, integration, downloads, and research. Covers **256+ databases** total.

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
├── research/                # Deep-dive research
└── downloads/               # Bulk download guides
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
| [pathways/processes.md](./pathways/processes.md) | Biological process databases catalog | Final |

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
| [literature/public-sources.md](./literature/public-sources.md) | Public paper sources & access | Final |
| [literature/coverage-analysis.md](./literature/coverage-analysis.md) | Literature coverage analysis | Final |
| [literature/abstracts-vs-fulltext.md](./literature/abstracts-vs-fulltext.md) | Abstract vs full-text comparison | Final |
| [literature/pipeline-design.md](./literature/pipeline-design.md) | Paper processing pipeline | Final |
| [literature/data-structures.md](./literature/data-structures.md) | Literature data structures | Final |

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
| [integration/size-estimates.md](./integration/size-estimates.md) | Data size estimates per source | Final |
| [integration/wikidata-pharma.md](./integration/wikidata-pharma.md) | Wikidata pharmaceutical integration | Final |
| [integration/integration-guide.md](./integration/integration-guide.md) | Data integration guide | Final |
| [integration/compound-pathway-linking.md](./integration/compound-pathway-linking.md) | Compound-pathway linking | Final |

---

## Schemas

Database schema documentation: [schemas/index.md](./schemas/index.md)

45 detailed schema files covering CLINVAR, DBSNP, CHEBI, UNIPROT, Reactome, and more.

---

## Research

| Document | Content | Status |
|----------|---------|--------|
| [research/interventions-priority.md](./research/interventions-priority.md) | Intervention source priorities | Final |
| [research/literature-priority.md](./research/literature-priority.md) | Literature source priorities | Final |

---

## Downloads

Bulk download guides for major data sources.

| Document | Content | Status |
|----------|---------|--------|
| [downloads/processing-pipeline.md](./downloads/processing-pipeline.md) | Data processing pipeline | Final |
| [downloads/traditional-medicine.md](./downloads/traditional-medicine.md) | Traditional medicine downloads | Final |
| [downloads/pharmaceuticals.md](./downloads/pharmaceuticals.md) | Pharmaceutical data downloads | Final |
| [downloads/pathways-targets.md](./downloads/pathways-targets.md) | Pathway & target downloads | Final |
| [downloads/wikidata-bulk.md](./downloads/wikidata-bulk.md) | Wikidata bulk downloads | Final |

---

## Gap Analysis

| Document | Content | Status |
|----------|---------|--------|
| [schema-gaps.md](./schema-gaps.md) | Missing data source details | Final |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 2.1 | January 2026 | Data Engineering | Added 18 files from research.old migration (downloads, literature, integration, research) |
| 2.0 | January 2026 | Data Engineering | Reorganized into domain subdirectories, removed 43-XX prefixes |
| 1.0 | January 2026 | Data Engineering | Initial index with migration mapping |
