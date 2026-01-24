---
id: natural.products
title: "Natural Products"
type: subcategory
parent: ../_index.md
last_updated: 2026-01-23
status: active
tags: [natural-products, plants, microbial, marine, secondary-metabolites]
---

# Natural Products

**Parent:** [Compounds & Molecules](../_index.md)

## Overview

Natural product databases catalog secondary metabolites from plants, microorganisms, and marine organisms. These resources are essential for drug discovery from nature, dereplication of known compounds, and understanding the chemical diversity of natural sources.

Key resources include LOTUS (comprehensive open natural products), COCONUT (aggregated natural products), Dr. Duke's (phytochemical database), NPASS (natural products), and NPAtlas (microbial natural products). Together they cover millions of natural compounds.

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [LOTUS](./lotus/_index.md) | 1 | Open natural products database |
| [COCONUT](./coconut/_index.md) | 1 | Collection of Open Natural Products |
| [Dr. Duke's](./dr.dukes/_index.md) | 2 | Phytochemical and ethnobotanical database |
| [NPASS](./npass/_index.md) | 2 | Natural Product Activity & Species Source |
| [NPAtlas](./npatlas/_index.md) | 2 | Atlas of microbial natural products |

## Integration Notes

LOTUS provides the most comprehensive open dataset with taxonomic source information. COCONUT aggregates multiple sources. Cross-reference with PubChem/ChEMBL for additional annotations. Use NPClassifier for chemical classification. Link to traditional medicine databases for ethnobotanical context.

---

## Database Selection Guide

### TL;DR

Natural products databases provide structural, bioactivity, and source organism data for 750K+ unique compounds from plants, microbes, and marine organisms. Priority sources are COCONUT (695K structures, CC0), LOTUS (750K structure-organism pairs, CC0), and NPASS (204K compounds with quantitative activity data). Target prediction tools (SwissTargetPrediction, PharmMapper) and ADMET predictors (ProTox 3.0, pkCSM) complete the compound-to-mechanism pipeline.

### Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Primary structure database | COCONUT 2.0 | Largest open (CC0), 695K structures, REST API |
| Structure-organism pairs | LOTUS/Wikidata | CC0 license, community-curated, SPARQL access |
| Quantitative bioactivity | NPASS | 1M+ activity records with IC50/Ki values |
| Microbial NPs | NPAtlas | Comprehensive microbial coverage, CC BY 4.0 |
| Target prediction | SwissTargetPrediction + PharmMapper | Complementary methods (2D/3D similarity + pharmacophore) |
| ADMET prediction | ProTox 3.0 + pkCSM | Free, comprehensive endpoints, no registration |
| Commercial databases | DNP excluded | Subscription-only, not suitable for open KB |

### Comparative Statistics

| Database | Compounds | Source Type | Key Feature | License |
|----------|-----------|-------------|-------------|---------|
| COCONUT 2.0 | 695K | Multi-source | Largest open NP collection | CC0 |
| LOTUS | 750K pairs | Taxonomic | Structure-organism linking | CC0 |
| NPASS | 204K | Bioactivity | Quantitative IC50/Ki data | Academic |
| NPAtlas | 35K+ | Microbial | Microbial NPs | CC BY 4.0 |
| SuperNatural 3.0 | 400K+ | Multi-source | Drug-likeness filtering | Academic |
| GNPS | Spectra | MS/MS | Spectral identification | Open |

### Tier Priorities

| Tier | Databases | Priority Reason |
|------|-----------|-----------------|
| **1 (Critical)** | COCONUT 2.0, LOTUS, NPASS | Comprehensive coverage; open licenses; bioactivity data |
| **2 (High Value)** | NPAtlas, GNPS, MIBiG | Microbial NPs; spectral data; biosynthesis |
| **3 (Enrichment)** | FooDB, Phenol-Explorer, KNApSAcK | Food compounds; polyphenols; specialized collections |

### Licensing Summary

| License Type | Databases | Commercial Use |
|--------------|-----------|----------------|
| CC0 (Public Domain) | COCONUT, LOTUS | Yes |
| CC BY 4.0 | NPAtlas, MIBiG | Yes (with attribution) |
| Open Access | FooDB, Phenol-Explorer, GNPS | Yes |
| Academic Only | NPASS, SuperNatural, UNPD, KNApSAcK | Research only |
