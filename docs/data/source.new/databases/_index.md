---
title: "Database Sources"
parent: ../_index.md
last_updated: 2026-01-22
status: draft
---

# Database Sources

This directory contains comprehensive documentation for all database sources organized by domain.

## Directory Structure

| Directory | Description | Category | Primary Focus |
|-----------|-------------|----------|---------------|
| [genetics/](./genetics/_index.md) | Genetic variation, clinical annotations, pharmacogenomics | genetics | SNPs, variants, drug-gene interactions |
| [traditional/](./traditional/_index.md) | Traditional medicine systems and herbal knowledge | traditional | TCM, Ayurveda, Kampo, Western Herbal |
| [nutrition/](./nutrition/_index.md) | Nutritional composition and food chemistry | nutrition | Food compounds, nutritional data |
| [pathways/](./pathways/_index.md) | Biological pathways and disease mechanisms | shared | Metabolic pathways, disease networks |
| [literature/](./literature/_index.md) | Scientific publications and research articles | shared | PubMed, PMC, OpenAlex |
| [compounds/](./compounds/_index.md) | Natural products and pharmaceutical compounds | shared | Chemical structures, bioactivity |

## Implementation Tiers

### Tier 1 (MVP)
- **Genetics**: dbSNP, ClinVar, gnomAD, PharmGKB
- **Traditional**: BATMAN-TCM, HERB, TCMSID, IMPPAT
- **Nutrition**: FooDB, USDA
- **Pathways**: Reactome, DisGeNET
- **Literature**: PubMed
- **Compounds**: COCONUT, LOTUS

### Tier 2 (Post-MVP)
- **Genetics**: dbNSFP, ExAC, GTEx
- **Traditional**: NPASS, KampoDB, Regional databases
- **Nutrition**: Phenol-Explorer, PhytoHub
- **Pathways**: KEGG, WikiPathways
- **Literature**: PMC, OpenAlex
- **Compounds**: DrugBank, ChEMBL

### Tier 3 (Future)
- Specialized genetics databases
- Regional traditional medicine databases
- Advanced nutritional databases
- Proprietary pathway databases
- Commercial compound databases

## Data Integration Strategy

### Genetics Integration
- Primary: dbSNP → ClinVar → PharmGKB
- Annotation: gnomAD → dbNSFP
- Expression: GTEx → UK Biobank

### Traditional Medicine Integration
- TCM: BATMAN-TCM → HERB → TCMSID
- Ayurveda: IMPPAT
- Cross-system: NPASS (natural product activity)

### Nutrition Integration
- Composition: FooDB → USDA
- Phytochemicals: Phenol-Explorer → PhytoHub

### Cross-Category Integration
- Compounds: COCONUT/LOTUS → DrugBank/ChEMBL
- Pathways: Reactome → KEGG → WikiPathways
- Literature: PubMed → PMC → OpenAlex
- Disease: DisGeNET bridges all categories

## Access Methods

| Method | Databases | Notes |
|--------|-----------|-------|
| REST API | dbSNP, ClinVar, PubMed, Reactome, DrugBank | Rate limits apply |
| Bulk Download | gnomAD, dbNSFP, FooDB, USDA, COCONUT | Large files (GB-TB) |
| SPARQL | WikiPathways, DisGeNET | RDF/semantic queries |
| MySQL/PostgreSQL | KEGG (subscription), local mirrors | Requires local setup |
| Web Scraping | BATMAN-TCM, HERB | Rate limiting required |

## Update Frequencies

- **Real-time**: PubMed, PMC
- **Monthly**: dbSNP, ClinVar, PharmGKB
- **Quarterly**: gnomAD, FooDB, USDA
- **Biannual**: KEGG, Reactome, DrugBank
- **Annual**: dbNSFP, GTEx, traditional medicine databases

## Navigation

- **Parent**: [Data Sources](../_index.md)
- **Subdirectories**: See table above
