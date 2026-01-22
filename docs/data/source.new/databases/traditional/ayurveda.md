---
id: traditional-ayurveda
title: Ayurveda Data Sources
category: traditional
tier: 2
subcategory: ayurveda
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [traditional, ayurveda, herbs, plants]
---

# Ayurveda Data Sources

**Document ID:** 43-22-AYURVEDA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

Ayurveda databases provide 4,000+ medicinal plants, 170,000+ phytochemicals, and 450,000+ traditional formulations. IMPPAT 2.0 is the best overall resource with CC BY 4.0 licensing and predicted protein targets. TKDL contains the largest formulation repository but has restricted access. Integration priority focuses on open-access sources (IMPPAT, OSADHI, GRAYU, NPACT) with supplementary data from specialized databases.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary data source | IMPPAT 2.0 | Best combination of coverage, targets, CC BY 4.0 license | Jan 2026 |
| Formulation data source | GRAYU | Only source with formulation-disease mappings | Jan 2026 |
| Target prediction source | IMPPAT + NPACT | Predicted (IMPPAT) + validated (NPACT) targets | Jan 2026 |
| Geographic coverage | OSADHI | Best state-level distribution data | Jan 2026 |
| Restricted sources | TKDL deferred | Patent office access only, future subscription | Jan 2026 |
| Literature source | DHARA | 90% of Ayurveda literature not in PubMed | Jan 2026 |

---

## Database Catalog

### 1. IMPPAT 2.0 (Indian Medicinal Plants, Phytochemistry And Therapeutics)

| Field | Value |
|-------|-------|
| **URL** | https://cb.imsc.res.in/imppat/ |
| **Content** | Comprehensive phytochemical database with predicted protein targets |
| **Records** | 4,010 plants; 17,967 phytochemicals; 1,095 therapeutic uses; 27,365 predicted interactions; 5,042 human target proteins |
| **License** | CC BY 4.0 (permits use, sharing, adaptation with attribution) |
| **API** | No public API; web export to TSV; structure downloads (SDF, MOL, MOL2, PDB, PDBQT) |
| **Update Frequency** | Version 2.0 released June 2022 |
| **Priority** | Tier 1 (Primary) |
| **Storage Estimate** | ~500 MB (compounds + structures + associations) |

**Data Structure:**
- Plant -> Plant Part -> Phytochemical associations
- Plant -> Plant Part -> Therapeutic Use associations
- Six data tabs per compound: Summary, Physicochemical, Drug-likeness, ADMET, Descriptors, Predicted Targets

**Key Metrics:**
- 5-fold increase in plant-phytochemical associations over v1.0
- 1,335 drug-like phytochemicals (filtered subset)
- High-confidence target predictions (combined score >= 700)

**GitHub:** https://github.com/asamallab/IMPPAT2

**Publication:** ACS Omega, 8(9):8827-8845, 2023

---

### 2. OSADHI (Online Structural and Analytics Database for Herbs of India)

| Field | Value |
|-------|-------|
| **URL** | https://neist.res.in/osadhi/ |
| **Content** | 4Ds approach: Documentation, Digitization, Deposition, Data Science |
| **Records** | 6,959 unique plants; 27,440 phytochemicals; 348 plant families; 2,477 therapeutic uses; 28 states + 8 UTs coverage |
| **License** | Open access (terms not explicitly stated) |
| **API** | No documented API |
| **Update Frequency** | Active (2022 publication) |
| **Priority** | Tier 1 (Geographic coverage) |
| **Storage Estimate** | ~800 MB (includes geographic data) |

**Unique Features:**
- Largest unique plant coverage (6,959 species)
- Geographic distribution across all Indian states/UTs
- Chemoinformatics: physicochemical properties, ADMET, NPClassifier
- Antiviral potency predictions (Random Forest & XGBoost)
- 2D and 3D structure downloads

**Data Fields:**
- Traditional knowledge: vernacular names, plant parts, taxonomy
- Phytochemicals: SMILES, InChIKey, IUPAC, structures
- Geographic classification by state

**Publication:** Computational Biology and Chemistry (2022)

---

### 3. GRAYU (Graph-based Database integrating Ayurvedic formulations)

| Field | Value |
|-------|-------|
| **URL** | https://caps.ncbs.res.in/GRAYU/ |
| **Content** | Graph-based Ayurvedic formulation database with disease ontology mapping |
| **Records** | 1,039 formulations; 12,743 plants; 129,542 phytochemicals; 13,480 diseases; 1,382,362 plant-compound associations; 4,087 formulation-disease associations |
| **License** | Research & Education Only (no clinical use; citation required) |
| **API** | No documented API; Interactive Knowledge Graph available |
| **Update Frequency** | Publication 2025 |
| **Priority** | Tier 1 (Formulations) |
| **Storage Estimate** | ~2 GB (large association network) |

**Unique Value:**
- Only database linking formulations to diseases via ontology
- Maps Sanskrit nosology to MeSH and DOID
- Graph-based architecture for network pharmacology
- Aggregates from: CMAUP, FooDB, HMDB, IMPPAT, NPASS, PubChem

**Association Counts:**
- Plant-Phytochemical: 1,382,362
- Plant-Disease: 116,824
- Plant-Formulation: 2,405
- Formulation-Disease: 4,087

**Publication:** Frontiers in Pharmacology (2025)

---

### 4. TKDL (Traditional Knowledge Digital Library)

| Field | Value |
|-------|-------|
| **URL** | https://www.tkdl.res.in/ |
| **Content** | Largest repository of traditional Indian medicine formulations |
| **Records** | 454,885+ total (119,269 Ayurveda; 236,399 Unani; 54,689 Siddha; 4,151 Yoga; 4,377 Sowa Rigpa) |
| **License** | Highly Restricted (Patent Offices only via TKDL Access Agreement) |
| **API** | No public API |
| **Update Frequency** | As of March 2022 |
| **Priority** | Tier 3 (Future - paid subscription model planned) |
| **Storage Estimate** | N/A (restricted access) |

**Access Details:**
- Current: 16 Patent Offices (EPO, USPTO, JPO, UKIPO, CIPO, DPMA, IP Australia, Indian Patent Office, Chile, Malaysia, etc.)
- Future: Paid subscription for businesses, research institutions, educational institutions, practitioners

**Classification System (TKRC):**
- ~5,000 subgroups for Ayurveda, Unani, Siddha, Yoga
- Correlates traditional to modern terminology
- Recognized by WIPO, incorporated into IPC under A61K 36/00

**Languages:** English, Japanese, French, German, Spanish

**Search Features:** Full text, Boolean, proximity, phrase, IPC/TKRC codes

**Contact:** headtkdl@csir.res.in

---

### 5. DHARA (Digital Helpline for Ayurveda Research Articles)

| Field | Value |
|-------|-------|
| **URL** | https://dharaonline.org/ |
| **Content** | Ayurveda-specific literature index |
| **Records** | 10,583+ articles (4,064 full text; 24,129 unique authors); 4,100+ journals indexed |
| **License** | Free public service for academic research |
| **API** | No API |
| **Update Frequency** | Active (maintained by AVP Research Foundation) |
| **Priority** | Tier 2 (Literature discovery) |
| **Storage Estimate** | ~50 MB (bibliographic index) |

**Coverage:**
- 17 Ayurveda-specific journals
- 18 CAM-related journals
- 846 mainstream medical journals

**Critical Value:** <10% of DHARA-indexed articles retrievable via PubMed "Ayurveda" search

**Search Features:** Keyword, advanced Boolean, field tags

---

### 6. AYUSH Research Portal

| Field | Value |
|-------|-------|
| **URL** | https://ayushportal.nic.in/ or https://arp.ayush.gov.in/ |
| **Content** | Government research repository for AYUSH systems |
| **Records** | 2.1M+ article views; 292K+ downloads |
| **License** | Academic use only |
| **API** | No documented API |
| **Update Frequency** | Active (Government maintained) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB (article metadata) |

**Systems Covered:** Ayurveda, Yoga, Naturopathy, Unani, Siddha, Sowa-Rigpa, Homeopathy

**Research Categories:** Clinical, Pre-clinical, Drug, Fundamental

---

### 7. NPACT (Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target)

| Field | Value |
|-------|-------|
| **URL** | http://crdd.osdd.net/raghava/npact/ |
| **Content** | Anti-cancer compounds with experimentally validated targets |
| **Records** | 1,574 compounds; 353 cell lines; 27 cancer types; ~5,214 compound-cell interactions; ~1,980 validated compound-target interactions |
| **License** | Open access |
| **API** | No API; MOL structure downloads; bulk via ZINC |
| **Update Frequency** | Static (2013 publication) |
| **Priority** | Tier 1 (Validated targets) |
| **Storage Estimate** | ~100 MB |

**Data Fields:**
- Compound: NPACT ID, IUPAC, SMILES, InChI, CAS
- Activity: IC50, ED50, EC50, GI50 values
- Targets: Molecular targets (experimentally validated)
- Suppliers: 50 commercial sources

**Bulk Download:** ZINC Catalog NPACT

**Publication:** Nucleic Acids Research, 41, D1124-D1129, 2013

---

### 8. CMAUP (Collective Molecular Activities of Useful Plants)

| Field | Value |
|-------|-------|
| **URL** | https://bidd.group/CMAUP/ |
| **Content** | Global medicinal plant database including Indian species |
| **Records** | 5,645 plants; ~48,000 phytochemicals; 646 targets; 234 KEGG pathways; 2,473 Gene Ontologies; 656 diseases |
| **License** | Freely accessible; all data downloadable |
| **API** | No API; free downloads from entry pages |
| **Update Frequency** | 2024 update (clinical trials, DNA barcodes, oral bioavailability) |
| **Priority** | Tier 1 (Target + pathway data) |
| **Storage Estimate** | ~400 MB |

**Plant Categories:** 2,567 medicinal; 170 food; 1,567 edible; 3 agricultural; 119 garden

**Geographic Coverage:** 153 countries with interactive world map

**Search Options:** Keywords, usage classes, families, targets, KEGG, GO, ICD diseases, locations

**Publication:** Nucleic Acids Research, 47(D1), D1118-D1127, 2019

---

### 9. NPASS (Natural Product Activity and Species Source)

| Field | Value |
|-------|-------|
| **URL** | http://bidd.group/NPASS |
| **Content** | Natural products with quantitative activity data |
| **Records** | 35,032+ NPs; 25,041+ species; ~7,700 targets; 222,092+ NP-target pairs; ~95,000 composition records |
| **License** | Freely accessible |
| **API** | No API; free downloads |
| **Update Frequency** | 2023 (3rd version) |
| **Priority** | Tier 1 (Activity data) |
| **Storage Estimate** | ~600 MB |

**Species Breakdown:** Plants 67.8%; Animals 9.4%; Fungi 7.9%; Bacteria 6.7%

**Activity Values:** IC50, Ki, EC50, GI50, MIC (quantitative)

**Additional Data:** Drug-likeness, ADMET, Chemical Checker activity profiles

**Publication:** Nucleic Acids Research (2023)

---

### 10. AromaDb (CIMAP Aroma Molecules Database)

| Field | Value |
|-------|-------|
| **URL** | https://aromadb.cimapbioinfo.in/ or https://bioinfo.cimap.res.in/aromadb/ |
| **Content** | Aroma molecules from Indian medicinal and aromatic plants |
| **Records** | 1,321 structures; 357 fragrance types; 166 commercial plants; 148 high-yielding varieties |
| **License** | Open access |
| **API** | No API; 2D/3D structure downloads |
| **Update Frequency** | Active |
| **Priority** | Tier 2 (Aromatherapy/essential oils) |
| **Storage Estimate** | ~50 MB |

**Data Fields:**
- IUPAC, CAS, chemical classification, functional groups
- 2D/3D structures (volatile <300 MW, medium <500 MW)
- ADMET properties, bioactivities, pathways

**Search Types:** Simple, advanced, similarity-based, property-based

**Publication:** Frontiers in Plant Science, 2018

---

### 11. sCentInDB (Essential Oil Chemical Profiles Database)

| Field | Value |
|-------|-------|
| **URL** | https://cb.imsc.res.in/scentindb/ |
| **Content** | Essential oil chemical profiles of Indian medicinal plants |
| **Records** | 554 plants; 2,170 EO profiles; 3,420 chemicals; 471 therapeutic use associations; 778 source articles |
| **License** | Academic research use |
| **API** | No API; CSV downloads |
| **Update Frequency** | 2025 publication |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~30 MB |

**Unique Features:**
- Plant-part level data compilation
- Sensory attributes (odor, color)
- FAIR-compliant design
- Pharmacokinetic data

**Publication:** Molecular Diversity, 2025

---

### 12. Indian Medicinal Plants Database (FRLHT/NMPB)

| Field | Value |
|-------|-------|
| **URL** | https://www.medicinalplants.in/ |
| **Content** | Botanical and vernacular names with medicinal uses |
| **Records** | 7,263 botanical names; 150,000+ vernacular names (10 Indian languages); 5,000+ plant images |
| **License** | Research and Development Only (commercial use prohibited) |
| **API** | No API; registration required |
| **Update Frequency** | Active |
| **Priority** | Tier 2 (Vernacular name mapping) |
| **Storage Estimate** | ~200 MB (includes images) |

**Systems Covered:** Ayurveda, Unani, Siddha, Homeopathy, + 2 others

---

### 13. IMPDB (Indian Medicinal Phytochemical Database)

| Field | Value |
|-------|-------|
| **URL** | https://impdb.org.in/ |
| **Content** | Phytochemicals from "Database on Medicinal Plants Used in Ayurveda" |
| **Records** | Based on CCRAS publication |
| **License** | Not explicitly stated |
| **API** | No API |
| **Update Frequency** | 2022 publication |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~20 MB |

**Key Finding:** 99.96-99.99% of phytochemicals are novel vs 2,068 FDA-approved drugs

**Technical Stack:** PHP front-end, MySQL back-end

**Publication:** Journal of Computational Biology and Chemistry, 2022

---

### 14. ABIM (Annotated Bibliography of Indian Medicine)

| Field | Value |
|-------|-------|
| **URL** | https://indianmedicine.eldoc.ub.rug.nl/ |
| **Content** | Literature citations for Indian medicine |
| **Records** | 50,000+ citations |
| **License** | Academic/research use |
| **API** | No API |
| **Update Frequency** | Static |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~30 MB |

**Features:** Simple search, diacritical marks for Sanskrit terms

**Platform:** EPrints 3, University of Southampton

---

### 15. e-Charak Portal

| Field | Value |
|-------|-------|
| **URL** | https://echarak.ayush.gov.in/ |
| **Content** | Classical and modern literature references for medicinal plants |
| **Records** | 22,000+ references for 16 selected plants |
| **License** | Research and development only (commercial use prohibited) |
| **API** | No API; registration required |
| **Update Frequency** | Active (Government maintained) |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~50 MB |

**Categories:** Botany, Chemistry, Pharmacology, Pharmacy, Miscellaneous, Reprints

---

### 16. Dr. Duke's Phytochemical and Ethnobotanical Databases

| Field | Value |
|-------|-------|
| **URL** | https://phytochem.nal.usda.gov/ |
| **Content** | Global ethnobotanical data including Ayurvedic applications |
| **Records** | Comprehensive species, phytochemicals, bioactivities |
| **License** | CC0 Public Domain (unrestricted use) |
| **API** | No API; bulk CSV downloads available |
| **Update Frequency** | Static (historical) |
| **Priority** | Tier 2 (CC0 license advantage) |
| **Storage Estimate** | ~150 MB |

**Download Sources:**
- Ag Data Commons: https://data.nal.usda.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases
- Figshare: https://agdatacommons.nal.usda.gov/articles/dataset/Dr_Duke_s_Phytochemical_and_Ethnobotanical_Databases/24660351
- Data.gov: https://catalog.data.gov/dataset/dr-dukes-phytochemical-and-ethnobotanical-databases-0849e

**Data Dictionary:** DrDukesDatabaseDataDictionary-prelim.csv

---

## Siddha Medicine Databases

### TKDL Siddha

| Field | Value |
|-------|-------|
| **Records** | 54,689 formulations |
| **License** | Restricted (via TKDL) |
| **Priority** | Tier 3 |

### Expert System for Siddha (eSS)

| Field | Value |
|-------|-------|
| **Status** | Pilot project |
| **Content** | 110 preparations from 2 Siddhars (Agathiyar, Theran) |
| **Purpose** | Pattern extraction from Siddha prescriptions |
| **Challenge** | Entire literature in Tamil (poems/padal format) |

---

## Unani Medicine Databases

### TKDL Unani

| Field | Value |
|-------|-------|
| **Records** | 236,399 formulations (largest in TKDL) |
| **License** | Restricted (via TKDL) |
| **Priority** | Tier 3 |

### National Formulary of Unani Medicine (NFUM)

| Field | Value |
|-------|-------|
| **Content** | 6 volumes; 1,228 standardized formulations |

### Unani Pharmacopoeia of India (UPI)

| Field | Value |
|-------|-------|
| **Content** | 298 single drugs (6 volumes, Part I); 100 compound drugs (2 volumes, Part II) |

### CCRUM

| Field | Value |
|-------|-------|
| **URL** | https://ccrum.res.in/ |
| **Network** | 31 institutes/units across India |

---

## Sowa Rigpa (Tibetan Medicine) Databases

### TKDL Sowa Rigpa

| Field | Value |
|-------|-------|
| **Records** | 4,377 formulations |
| **License** | Restricted (via TKDL) |

### Dataset of Materia Medica in Sowa Rigpa

| Field | Value |
|-------|-------|
| **Content** | ~700 commonly used materia medica |
| **Format** | Unicode Tibetan font, Wylie transliteration, phonetic transcription |
| **Publication** | ScienceDirect, Data in Brief (2020) |

**Note:** Government of India officially recognized Sowa Rigpa in 2010 as a "science of healing."

---

## Comparison Matrix

| Database | Plants | Compounds | Formulations | Targets | Open Access | Bulk Download | API |
|----------|--------|-----------|--------------|---------|-------------|---------------|-----|
| **IMPPAT 2.0** | 4,010 | 17,967 | - | 27,365 predicted | CC BY 4.0 | TSV export | No |
| **OSADHI** | 6,959 | 27,440 | - | No | Yes | Partial | No |
| **GRAYU** | 12,743 | 129,542 | 1,039 | Indirect | Research only | No | No |
| **TKDL** | - | - | 454,885+ | No | Restricted | No | No |
| **DHARA** | - | - | - | No | Yes | No | No |
| **NPACT** | - | 1,574 | - | 1,980 validated | Yes | MOL files | No |
| **CMAUP** | 5,645 | 48,000 | - | 646 | Yes | Yes | No |
| **NPASS** | - | 35,000+ | - | 222,092 pairs | Yes | Yes | No |
| **Dr. Duke's** | Many | Many | - | Limited | CC0 | CSV | No |
| **AromaDb** | 166 | 1,321 | - | Limited | Yes | 2D/3D | No |
| **sCentInDB** | 554 | 3,420 | - | No | Yes | CSV | No |
| **FRLHT/NMPB** | 7,263 | - | - | No | Registration | No | No |

---

## Schema Overview

### Core Data Model

| Entity | Description | Key Fields |
|--------|-------------|------------|
| Plant | Medicinal plant record | plant_id, botanical_name, common_names, family, vernacular_names |
| Compound | Phytochemical compound | compound_id, iupac_name, smiles, inchi_key, cas_number |
| Formulation | Traditional preparation (Yoga) | formulation_id, sanskrit_name, components, dosage_form |
| Target | Biological/protein target | target_id, protein_name, uniprot_id, gene_symbol |
| Therapeutic Use | Traditional indication | use_id, sanskrit_term, english_term, body_system |
| Plant Part | Botanical organ used | part_id, part_name (root, leaf, seed, etc.), preparation_method |

### Key Tables (from Primary Databases)

#### IMPPAT 2.0 Schema

| Table | Description | Key Fields |
|-------|-------------|------------|
| `plants` | Indian medicinal plants | plant_id, botanical_name, family, common_names |
| `phytochemicals` | Chemical compounds | compound_id, name, smiles, inchi_key, pubchem_cid |
| `plant_parts` | Plant organ associations | plant_id, part_id, part_name |
| `plant_phytochemical` | Plant-compound links | plant_id, part_id, compound_id |
| `therapeutic_uses` | Traditional indications | use_id, use_term, plant_id, part_id |
| `predicted_targets` | Target predictions | compound_id, target_id, combined_score |

#### GRAYU Schema (Graph-Based)

| Node Type | Description | Key Attributes |
|-----------|-------------|----------------|
| `Formulation` | Ayurvedic yoga/preparation | formulation_id, name, source_text |
| `Plant` | Medicinal plant | plant_id, botanical_name |
| `Phytochemical` | Chemical compound | compound_id, smiles, inchi_key |
| `Disease` | Disease entity | disease_id, mesh_id, doid, sanskrit_name |

| Edge Type | Description | Properties |
|-----------|-------------|------------|
| `CONTAINS` | Formulation contains plant | quantity, unit |
| `PRODUCES` | Plant produces compound | plant_part, concentration |
| `TREATS` | Formulation/plant treats disease | evidence_type, source |
| `TARGETS` | Compound targets protein | binding_score, mechanism |

#### OSADHI Schema

| Table | Description | Key Fields |
|-------|-------------|------------|
| `plants` | Plant records with geographic data | plant_id, botanical_name, family, state_distribution |
| `compounds` | Phytochemicals | compound_id, smiles, inchi_key, molecular_weight |
| `admet_properties` | ADME-Tox predictions | compound_id, oral_bioavailability, bbb_permeability |
| `therapeutic_uses` | Traditional uses | plant_id, use_term, part_used |

### Identifier Cross-References

| Database | Plant ID | Compound ID | Target ID |
|----------|----------|-------------|-----------|
| IMPPAT | IMPPAT_P#### | IMPPAT_C#### | UniProt |
| OSADHI | Internal | PubChem CID | N/A |
| GRAYU | Internal | PubChem CID | UniProt |
| NPACT | N/A | NPACT_ID | UniProt |
| CMAUP | CMAUP_P#### | CMAUP_C#### | UniProt |

### Common Data Formats

| Format | Use Case | Source Databases |
|--------|----------|------------------|
| TSV | Tabular exports | IMPPAT, OSADHI |
| SDF/MOL | Chemical structures | IMPPAT, CMAUP, NPACT |
| CSV | Bulk downloads | Dr. Duke's, NPASS |
| JSON | API responses | PubChem integration |

---

## Integration Recommendations

### Priority 1: Primary Data Sources (Highest Value)

| Database | Rationale | Integration Effort |
|----------|-----------|-------------------|
| **IMPPAT 2.0** | Best compound-target coverage, CC BY 4.0 license | Medium (web export) |
| **NPASS** | Quantitative activity data, free bulk download | Medium |
| **CMAUP** | Global coverage, target + pathway data | Medium |
| **NPACT** | Validated anti-cancer targets | Low |

### Priority 2: Supplementary Sources

| Database | Rationale | Integration Effort |
|----------|-----------|-------------------|
| **OSADHI** | Geographic + traditional knowledge context | High (no bulk) |
| **GRAYU** | Formulation-level data (unique value) | High (no bulk) |
| **Dr. Duke's** | CC0 license, ethnobotanical context | Low (CSV) |
| **AromaDb** | Essential oil specifics | Low |

### Priority 3: Reference Sources

| Database | Rationale | Integration Effort |
|----------|-----------|-------------------|
| **DHARA** | Literature discovery (90% not in PubMed) | Medium |
| **AYUSH Portal** | Government research articles | Medium |
| **sCentInDB** | Essential oil chemical profiles | Low |
| **FRLHT/NMPB** | Vernacular name mapping | Medium |

### Data Integration Strategy

1. **Core Compound-Target Data**
   - Start with IMPPAT 2.0 predicted targets
   - Supplement with NPACT validated cancer targets
   - Add NPASS quantitative activity data
   - Cross-reference with CMAUP pathway data

2. **Traditional Knowledge Layer**
   - Plant-compound associations from IMPPAT + OSADHI
   - Formulation data from GRAYU
   - Therapeutic use mappings

3. **Disease Ontology Mapping**
   - Leverage GRAYU's MeSH/DOID mappings
   - Connect traditional indications to modern disease terms

4. **Literature Support**
   - Index via DHARA and AYUSH Portal
   - Link compounds to supporting research

---

## Technical Considerations

| Consideration | Details |
|---------------|---------|
| **No APIs available** | All databases require web scraping or manual export |
| **Data harmonization** | Different identifier systems across databases |
| **PubChem integration** | Most databases reference PubChem for compound standardization |
| **STRING/ChEMBL linkage** | Required for target prediction validation |
| **Language barriers** | Siddha (Tamil), Sowa Rigpa (Tibetan) require special handling |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| IMPPAT 2.0 | 500 MB (Compounds + structures + associations) |
| OSADHI | 800 MB (Includes geographic data) |
| GRAYU | 2 GB (Large association network) |
| NPACT | 100 MB (Compounds + targets) |
| CMAUP | 400 MB (Global plant data) |
| NPASS | 600 MB (Activity data) |
| AromaDb | 50 MB (Aroma molecules) |
| sCentInDB | 30 MB (EO profiles) |
| FRLHT/NMPB | 200 MB (Includes images) |
| Dr. Duke's | 150 MB (Ethnobotanical) |
| Literature DBs | 200 MB (Metadata only) |
| Total storage estimate | ~5 GB (Core Ayurveda data) |
| Last updated | January 2026 |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [index.md](./../index.md) | Parent index |
| [natural-products.md](./../compounds/natural-products.md) | Cross-references COCONUT, LOTUS |

---

## References

1. Vivek-Ananth RP, et al. (2023). IMPPAT 2.0: An Enhanced and Expanded Phytochemical Atlas of Indian Medicinal Plants. ACS Omega, 8(9):8827-8845.
2. Mangal M, et al. (2013). NPACT: Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target database. Nucleic Acids Research, 41, D1124-D1129.
3. Zeng X, et al. (2019). CMAUP: a database of collective molecular activities of useful plants. Nucleic Acids Research, 47(D1), D1118-D1127.
4. Kumar Y, et al. (2018). AromaDb: A Database of Medicinal and Aromatic Plant's Aroma Molecules. Frontiers in Plant Science.
5. Baskaran SP, et al. (2025). sCentInDB: a database of essential oil chemical profiles of Indian medicinal plants. Molecular Diversity.
6. TKDL Official Website: https://www.tkdl.res.in/
7. OSADHI: Computational Biology and Chemistry (2022)
8. GRAYU: Frontiers in Pharmacology (2025)

---

## License

This document catalogs multiple databases with varying license terms:

| Database | License | Commercial Use | Attribution | Access |
|----------|---------|----------------|-------------|--------|
| IMPPAT 2.0 | CC BY 4.0 | Yes | Required | Open |
| OSADHI | Open access | Yes | Citation | Open |
| GRAYU | Research & Education Only | No | Required | Open |
| TKDL | Highly Restricted | No | N/A | Patent Offices only |
| DHARA | Free public service | Research only | Citation | Open |
| AYUSH Research Portal | Academic use | Research only | Required | Open |
| NPACT | Open access | Yes | Citation | Open |
| CMAUP | Freely accessible | Yes | Citation | Open (downloads) |
| NPASS | Freely accessible | Yes | Citation | Open (downloads) |
| AromaDb | Open access | Yes | Citation | Open |
| sCentInDB | Academic research use | Research only | Required | Open |
| FRLHT/NMPB | Research & Development Only | No | Required | Registration required |
| IMPDB | Not explicitly stated | Unknown | Citation | Open |
| ABIM | Academic/research use | Research only | Citation | Open |
| e-Charak Portal | Research & Development Only | No | Required | Registration required |
| Dr. Duke's | CC0 (Public Domain) | Yes | None required | Open (bulk CSV) |

**Key Considerations:**
- **Fully Open (Commercial OK):** IMPPAT 2.0 (CC BY 4.0), Dr. Duke's (CC0), NPACT, CMAUP, NPASS, OSADHI, AromaDb
- **Research/Academic Only:** GRAYU, DHARA, AYUSH Portal, sCentInDB, FRLHT/NMPB, ABIM, e-Charak
- **Highly Restricted:** TKDL (Patent Office access only)
- **CC0 Public Domain:** Dr. Duke's allows unrestricted commercial use

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **IMPPAT 2.0** | Web | `https://cb.imsc.res.in/imppat/` |
| **Dr. Duke's** | Bulk CSV | `https://phytochem.nal.usda.gov/phytochem/search` |
| **AYUSH Research Portal** | Web | `https://ayushportal.nic.in/` (registration required) |
| **NPACT** | Web | `http://crdd.osdd.net/raghava/npact/` |
| **CMAUP** | Download | `https://www.cmaup.cn/` |
| **NPASS** | Download | `https://bidd.group/NPASS/` |
| **TKDL** | Restricted | Patent office access only |

**Access Requirements:** Most databases are freely accessible; TKDL requires patent office authorization; AYUSH Portal and e-Charak require registration.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | SDF, MOL2, CSV |
| Alternative | TSV, JSON, XML |
| Chemical structures | SMILES, InChI, InChIKey |
| Encoding | UTF-8 |

## Sample Data

### Example Compound Record (IMPPAT 2.0)
```json
{
  "compound_id": "IMPPAT001234",
  "name": "Withaferin A",
  "smiles": "CC1(C)CCC2C(CC(=O)C3(C)C2CCC2C4CC(O)C(C(C)C(=O)O4)C(O)C23)C1",
  "plant_source": "Withania somnifera",
  "traditional_use": "Rasayana (rejuvenation)",
  "pharmacological_activity": ["anti-inflammatory", "immunomodulatory"]
}
```

### Sample Query Result
| compound_id | name | plant_source | dosha_effect |
|-------------|------|--------------|--------------|
| IMPPAT001234 | Withaferin A | Withania somnifera | Vata-balancing |
| IMPPAT002567 | Curcumin | Curcuma longa | Pitta-pacifying |

---

## Glossary

### Core Ayurvedic Terms

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Dosha` | The three fundamental bodily humors that govern physiological and psychological functions | Vata, Pitta, Kapha |
| `Vata` | Dosha governing movement, nervous system, and elimination; composed of air and space elements | Associated with dry, cold, light qualities |
| `Pitta` | Dosha governing metabolism, digestion, and transformation; composed of fire and water elements | Associated with hot, sharp, oily qualities |
| `Kapha` | Dosha governing structure, stability, and lubrication; composed of water and earth elements | Associated with heavy, slow, cool qualities |
| `Prakriti` | Constitutional type; an individual's inherent doshic balance determined at conception | Vata-Pitta prakriti |
| `Vikriti` | Current state of doshic imbalance; deviation from one's prakriti | Used in diagnosis |

### Pharmacological Properties

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Rasa` | Taste; one of six tastes that influence doshic balance | Madhura (sweet), Amla (sour), Lavana (salty), Katu (pungent), Tikta (bitter), Kashaya (astringent) |
| `Guna` | Quality or attribute of a substance | Guru (heavy), Laghu (light), Snigdha (oily), Ruksha (dry) |
| `Virya` | Potency; the heating or cooling effect of a substance | Ushna (hot), Shita (cold) |
| `Vipaka` | Post-digestive effect; the taste that emerges after digestion | Madhura, Amla, or Katu vipaka |
| `Prabhava` | Special potency; unique therapeutic action not explainable by rasa, guna, virya, or vipaka | Unexplained specific actions |

### Therapeutic Categories

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Rasayana` | Rejuvenation therapy; substances that promote longevity and vitality | Ashwagandha, Amalaki, Brahmi |
| `Vajikarana` | Aphrodisiac therapy; substances that enhance reproductive health | Shatavari, Kapikacchu |
| `Panchakarma` | Five purification therapies for detoxification | Vamana, Virechana, Basti, Nasya, Raktamokshana |
| `Shodhana` | Purification or detoxification therapy | Part of Panchakarma |
| `Shamana` | Palliative therapy; symptom relief without purification | Herbal decoctions for symptom management |

### Preparation Types (Dosage Forms)

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Churna` | Fine powder preparation | Triphala churna, Ashwagandha churna |
| `Kwatha` (Kashaya) | Decoction; water extract prepared by boiling | Dashamoola kwatha |
| `Asava` | Self-generated fermented preparation (no external heat) | Ashwagandharishta |
| `Arishta` | Fermented preparation with decoction as base | Draksharishta |
| `Taila` | Medicated oil preparation | Brahmi taila, Bhringaraj taila |
| `Ghrita` | Medicated ghee (clarified butter) preparation | Brahmi ghrita |
| `Vati/Gutika` | Tablet or pill form | Arogyavardhini vati |
| `Lepa` | Paste for external application | Face masks, poultices |
| `Avaleha/Lehya` | Confection or jam-like preparation | Chyawanprash |
| `Bhasma` | Calcined mineral/metal ash preparation | Abhrak bhasma, Swarna bhasma |

### Anatomical and Physiological Terms

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Agni` | Digestive fire; metabolic capacity | Jatharagni (digestive fire), Dhatvagni (tissue metabolism) |
| `Ama` | Toxins; undigested metabolic waste | Accumulates when agni is weak |
| `Dhatu` | Seven bodily tissues | Rasa, Rakta, Mamsa, Meda, Asthi, Majja, Shukra |
| `Srotas` | Channels or pathways in the body | Pranavaha srotas (respiratory), Annavaha srotas (digestive) |
| `Mala` | Waste products | Purisha (feces), Mutra (urine), Sveda (sweat) |
| `Ojas` | Vital essence; immunity and vitality | Product of proper digestion and metabolism |

### Database and Technical Terms

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Phytochemical` | Chemical compound produced by plants | Alkaloids, flavonoids, terpenoids |
| `ADMET` | Absorption, Distribution, Metabolism, Excretion, Toxicity | Drug-likeness properties |
| `TTI` | Target-Target Interaction or Target-Ingredient Interaction | Compound-protein binding data |
| `InChIKey` | International Chemical Identifier hash key | Unique 27-character compound identifier |
| `SMILES` | Simplified Molecular Input Line Entry System | Linear notation for chemical structures |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| IMPPAT | Indian Medicinal Plants, Phytochemistry And Therapeutics | Primary Ayurveda database |
| OSADHI | Online Structural and Analytics Database for Herbs of India | Geographic coverage focus |
| GRAYU | Graph-based Ayurveda | Formulation-disease network |
| TKDL | Traditional Knowledge Digital Library | Restricted patent office access |
| DHARA | Digital Helpline for Ayurveda Research Articles | Literature index |
| AYUSH | Ayurveda, Yoga, Unani, Siddha, Homeopathy | Indian Ministry designation |
| NPACT | Naturally Occurring Plant-based Anti-cancer Compound-Activity-Target | Cancer targets database |
| CMAUP | Collective Molecular Activities of Useful Plants | Global medicinal plant database |
| NPASS | Natural Product Activity and Species Source | Activity data source |
| CC BY 4.0 | Creative Commons Attribution 4.0 | Open license with attribution |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial comprehensive Ayurveda database catalog migrated from research.old |
