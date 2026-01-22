---
id: traditional-kampo
title: Kampo Medicine Data Sources
world: 2
category: traditional
tier: 2
subcategory: kampo
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [traditional, kampo, japanese, herbs]
---

# Kampo Medicine Data Sources

**Document ID:** 43-23-KAMPO
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

Kampo is Japanese traditional medicine derived from TCM, with 148 officially approved formulas in Japan. Primary sources include KampoDB (298 formulas, 62,906 protein targets, CC BY-SA 4.0), STORK (official reference for all 148 approved formulas), and TM-MC 2.0 (192 Japanese Pharmacopoeia herbs). Cross-reference with TCM databases (BATMAN-TCM, HERB) extends compound-target coverage significantly.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary Kampo source | KampoDB | Only comprehensive open database with target predictions, CC BY-SA 4.0 license | Jan 2026 |
| Official formula reference | STORK | Authoritative source for all 148 approved ethical Kampo drugs | Jan 2026 |
| Japanese herb coverage | TM-MC 2.0 | 192 JP Pharmacopoeia herbs with compound/target data, downloadable | Jan 2026 |
| Compound-species relationships | KNApSAcK KAMPO | 336 formulas, 278 plants, integrates with KampoDB IDs | Jan 2026 |
| Target prediction supplementation | BATMAN-TCM 2.0 | API access, 2.3M TTIs, compensates for Kampo-specific data gaps | Jan 2026 |
| Clinical evidence | EKAT | 500+ RCT abstracts, Cochrane-linked specialized register | Jan 2026 |

---

## Database Catalog

### Primary Kampo Databases

#### KampoDB

| Field | Value |
|-------|-------|
| **URL** | https://wakanmoview.inm.u-toyama.ac.jp/kampo/ |
| **Content** | Kampo formula-compound-target predictions with docking simulations |
| **Records** | 298 formulas, 180 crude drugs, 3,002 compounds, 62,906 proteins, 3,063,505 docking results |
| **License** | CC BY-SA 4.0 |
| **API** | No (web interface only) |
| **Update Frequency** | Active (current version) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 MB (formulas + compounds + target mappings) |
| **Maintainer** | Institute of Natural Medicine, University of Toyama |

**Key Features:**
- 460 known target proteins (experimental evidence)
- 1,369 predicted target proteins (docking/ML)
- Hierarchical structure: Formula -> Crude drugs -> Compounds -> Targets
- Compound IDs correspond to KNApSAcK IDs
- Chemical structures linked to KNApSAcK and PubChem

**Data Access Notes:**
- Web-based queries return structured data
- Bulk download not explicitly documented; contact maintainers
- Publication: Sawada R, et al. (2018) Scientific Reports 8:11216

---

#### STORK (Standards of Reporting Kampo Products)

| Field | Value |
|-------|-------|
| **URL** | http://mpdb.nibiohn.go.jp/stork/ |
| **Content** | Official reference for all 148 approved ethical Kampo drugs in Japan |
| **Records** | 148 prescription Kampo drugs, 241 approved crude drugs, 5 crude drug preparations |
| **License** | Open reference (free citation in research) |
| **API** | No |
| **Update Frequency** | As regulations change |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~5 MB (formula compositions + package insert data) |
| **Maintainer** | NIBIOHN + NIHS |

**Key Features:**
- Corresponds to CONSORT Statement 2010 Item 5 (Intervention)
- Formula composition with crude drug ratios
- Manufacturer information
- Package insert equivalent data (English)
- Glycyrrhiza (licorice) dose information

**Data Access Notes:**
- Reference use only; no bulk download
- Web scraping required for systematic extraction

---

#### TradMPD (Traditional Medical & Pharmaceutical Database)

| Field | Value |
|-------|-------|
| **URL** | https://dentomed.toyama-wakan.net/index_en |
| **Content** | Experimental data for Kampo natural medicines (first such repository) |
| **Records** | ~80 compounds (Wakanyaku Library), 80 crude drug extracts, 80 Kampo extracts |
| **License** | Academic/research use |
| **API** | No |
| **Update Frequency** | Static (2008-2012 grants) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~20 MB (profiles + analytical data) |
| **Maintainer** | Section of Pharmacognosy, Institute of Natural Medicine, University of Toyama |

**Key Features:**
- LC-MS profiling data
- Biological activity data for:
  - Neurodegenerative diseases
  - Allergy & inflammation
  - Cancer
  - Lifestyle diseases
- Genetic information on medicinal plants
- Literature searches on Kampo formulas

---

### Umbrella Portal

#### WAKANYAKU Database Portal

| Field | Value |
|-------|-------|
| **URL** | https://www.inm.u-toyama.ac.jp/en/database/ |
| **Content** | Umbrella portal hosting multiple traditional medicine databases |
| **Maintainer** | Institute of Natural Medicine, University of Toyama |

**Component Databases:**

| Database | URL | Content |
|----------|-----|---------|
| TradMPD | https://dentomed.toyama-wakan.net/index_en | Genetic, chemical, biological data on crude drugs & Kampo |
| WAKANYAKU Wiki | https://www.inm.u-toyama.ac.jp/wiki/ | Natural medicine reviews, LC-MS profiles, Edo-period texts |
| ETHMEDmmm | https://ethmed.toyama-wakan.net/SearchEn/ | Ethnomedicine database, Museum specimens |
| KampoDB | https://wakanmoview.inm.u-toyama.ac.jp/kampo/ | Compound-protein interactions |
| Shouruihonzou DB | https://ethmed.toyama-wakan.net/honzou | Historical Chinese materia medica text |

---

### Compound-Species Databases

#### KNApSAcK KAMPO

| Field | Value |
|-------|-------|
| **URL** | http://www.knapsackfamily.com/knapsack_core/top.php |
| **Content** | Metabolite-species relationships for Kampo formulas |
| **Records** | 336 Kampo formulae, 278 medicinal plants |
| **License** | Non-commercial (contact skanaya@gtc.naist.jp for commercial use) |
| **API** | URL pattern access: `info.php?sname=[type]&word=[query]` |
| **Update Frequency** | Active |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~30 MB (KAMPO subset) |
| **Maintainer** | Nara Institute of Science and Technology (NAIST) |

**Related KNApSAcK Databases:**

| Database | Records | Content |
|----------|---------|---------|
| KNApSAcK Core | 63,715 metabolites, 159,095 metabolite-species pairs, 24,749 species | Main metabolite database |
| KNApSAcK WorldMap | 41,548 GZ-plant pairs, 222 geographic zones, 15,240 plants | Geographic distribution |
| Biological Activity | 2,418 activities, 33,706 plant-activity relationships | Activity annotations |

**Search Capabilities:**
- Organism names
- Metabolite identifiers
- Molecular formulas
- C_ID, CAS_ID, INCHI-KEY, INCHI-CODE, SMILES

---

### Regional Integration Databases

#### TM-MC 2.0

| Field | Value |
|-------|-------|
| **URL** | https://tm-mc.kr |
| **Content** | Northeast Asian traditional medicine (Korean, Chinese, Japanese pharmacopoeias) |
| **Records** | 635 medicinal materials, 34,107 compounds (21,306 de-duplicated), 13,992 targets, 27,997 diseases, 5,075 prescriptions |
| **License** | Open access |
| **API** | No |
| **Update Frequency** | Active |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 MB (full database download) |
| **Maintainer** | Korean research community |

**Japanese Pharmacopoeia Coverage:**
- **192 medicinal materials** from Japanese Pharmacopoeia included
- Cross-referenced with Korean (454) and Chinese (556) pharmacopoeias
- Compounds deposited in PubChem (22,643)
- Identifiers and pharmacokinetic properties provided

**Download:** Downloadable files available directly from website

---

### Clinical Evidence Database

#### EKAT (Evidence Reports of Kampo Treatment)

| Field | Value |
|-------|-------|
| **URL** | http://www.jsom.or.jp/medical/ebm/ere/index.html |
| **Content** | Clinical evidence database for Kampo medicines |
| **Records** | 500+ RCT structured abstracts (345 in EKAT 2010 version) |
| **License** | Free access for physicians and researchers |
| **API** | No |
| **Update Frequency** | Annual updates |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~10 MB |
| **Maintainer** | Japan Society for Oriental Medicine (JSOM) |

**Key Features:**
- Structured abstracts with 12 components each
- Available in Japanese, English, and Korean
- Linked to Cochrane Library (CENTRAL) as Specialized Register

---

### Ethnomedicine Database

#### ETHMEDmmm

| Field | Value |
|-------|-------|
| **URL** | https://ethmed.toyama-wakan.net/SearchEn/ |
| **Content** | Herbal medicine specimens and ethnomedicine data |
| **Records** | ~30,000 herbal medicine specimens |
| **License** | Restricted (authorization required for imagery data) |
| **API** | No |
| **Update Frequency** | Active |
| **Priority** | Tier 3 |
| **Storage Estimate** | Reference only |
| **Maintainer** | Museum of Materia Medica, Institute of Natural Medicine, University of Toyama |

**Access Notes:**
- Free consultation for reference
- Contents and photographs cannot be diverted without authorization
- Must follow prescribed procedures for imagery data use

---

### Reference Database

#### NIBIOHN Medicinal Plant Database (MPDB)

| Field | Value |
|-------|-------|
| **URL** | http://mpdb.nibiohn.go.jp/ |
| **Content** | Common names, crude drug information, plant species, compounds, samples |
| **Records** | Variable |
| **License** | Open reference |
| **API** | No |
| **Update Frequency** | Active |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~15 MB |
| **Maintainer** | Research Center for Medicinal Plant Resources, NIBIOHN |

**Note:** Server reliability may be inconsistent.

---

### Wiki-Based Repository

#### Metabolomics.jp Wiki - Kampo Section

| Field | Value |
|-------|-------|
| **URL** | http://metabolomics.jp/wiki/ |
| **Content** | Wiki-based repository for crude drugs and Kampo medicine |
| **Records** | 158 crude drugs, 348 Kampo prescriptions |
| **License** | Open access (wiki-based) |
| **API** | Partial (MediaWiki API for page contents) |
| **Update Frequency** | Community-driven |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~10 MB |
| **Maintainer** | Research community |

**Key Features:**
- Taxonomic information
- Chemical information
- MediaWiki structure with SQL backend
- Japanese Pharmacopoeia entries on "Index:JPR" page

---

## Cross-Reference TCM Databases

Since Kampo derives from TCM, these databases provide essential supplementary compound-target data.

### BATMAN-TCM 2.0

| Field | Value |
|-------|-------|
| **URL** | http://bionet.ncpsb.org/batman-tcm |
| **Content** | Most comprehensive TCM database with API access |
| **Records** | 54,832 formulae, 8,404 herbs, 39,171 ingredients, 17,068 known TTIs, 2,319,272 predicted TTIs |
| **License** | Academic/research use |
| **API** | **Yes** (JSON format, URL-based parameters, ~3 sec/herb query) |
| **Update Frequency** | Active |
| **Priority** | Tier 1 (supplementary) |
| **Storage Estimate** | ~500 MB (full TTI data) |
| **Maintainer** | National Center for Protein Sciences Beijing |

**API Features:**
- Build URL with selected parameters
- Request type, output format, input item
- Returns JSON (computer-readable) or hypertext
- KEGG/GO/disease enrichment analysis
- Network visualization

**Publication:** Liu Z, et al. (2024) Nucleic Acids Research 52(D1):D1110-D1116

---

### HERB 2.0

| Field | Value |
|-------|-------|
| **URL** | http://herb.ac.cn/ |
| **Content** | High-throughput experiment-guided TCM database |
| **Records** | 7,263 herbs, 49,258 ingredients, 12,933 targets, 28,212 diseases, 2,231 experiments, 8,558 clinical trials |
| **License** | Open access |
| **API** | No (custom scripts downloadable: probe2gene) |
| **Update Frequency** | Active |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~200 MB |
| **Maintainer** | Research community |

**Publication:** Fang S, et al. (2021) Nucleic Acids Research 49(D1):D1197-D1206

---

### TCMBank

| Field | Value |
|-------|-------|
| **URL** | https://TCMBank.cn/ |
| **Content** | Largest downloadable non-commercial TCM database |
| **Records** | 9,192 herbs, 61,966 ingredients, 15,179 targets, 32,529 diseases |
| **License** | Non-commercial (free for research) |
| **API** | Partial (`/api` endpoint exists) |
| **Update Frequency** | Active |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~300 MB |
| **Maintainer** | AI Medical Center, Sun Yat-sen University |

**Data Format:**
- 3D structures in mol2 format
- Cross-references: CAS, DrugBank, PubChem, MeSH, OMIM, DO, ETCM, HERB

---

### SymMap 2.0

| Field | Value |
|-------|-------|
| **URL** | http://www.symmap.org/ |
| **Content** | TCM database with symptom mapping |
| **Records** | 1,717 TCM symptoms, 499 herbs, 961 modern symptoms, 5,235 diseases, 19,595 ingredients, 4,302 target genes |
| **License** | Open access |
| **API** | No (download button available) |
| **Update Frequency** | Active |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~80 MB |
| **Maintainer** | Research community |

---

### HIT 2.0 (Herbal Ingredients' Targets)

| Field | Value |
|-------|-------|
| **URL** | http://hit2.badd-cao.net |
| **Content** | Herbal ingredient-target relationships |
| **Records** | 10,031 compound-target activity pairs, 2,208 targets, 1,237 ingredients, 1,250+ herbs |
| **License** | Open access |
| **API** | No ("My-target" feature for personal curation) |
| **Update Frequency** | Active |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~40 MB |
| **Maintainer** | Research community |

**Crosslinks:** TTD, DrugBank, KEGG, PDB, UniProt, Pfam, NCBI, TCM-ID

---

### Additional TCM References

#### TCMID 2.0

| Field | Value |
|-------|-------|
| **URL** | http://www.megabionet.org/tcmid/ |
| **Download** | http://47.100.169.139/tcmid/download/ |
| **Records** | ~47,000 prescriptions, 8,159 herbs, 25,210+ compounds, 17,521+ targets |
| **License** | Free for research (citation requested) |
| **Storage Estimate** | ~150 MB |

#### YaTCM

| Field | Value |
|-------|-------|
| **URL** | http://cadd.pharmacy.nankai.edu.cn/yatcm/home |
| **Content** | Prescriptions, herbs, ingredients, diseases, targets, pathways |
| **Features** | MV-SEA target prediction, similarity search, substructure search, network analysis |
| **Maintainer** | Nankai University |

#### SuperTCM

| Field | Value |
|-------|-------|
| **URL** | http://tcm.charite.de/supertcm |
| **Records** | 6,516 herbs, 5,372 species, 55,772 ingredients, 543 targets, 254 KEGG pathways, 8,634 diseases |
| **Features** | KEGG Global Maps visualization, iPath3.0 pathway projection |
| **License** | Open access |
| **Maintainer** | Charite - Universitatsmedizin Berlin |

---

## General Reference Databases

### COCONUT 2.0

| Field | Value |
|-------|-------|
| **URL** | https://coconut.naturalproducts.net |
| **Content** | Largest open natural products database (63+ sources integrated) |
| **License** | CC0 (Public domain) |
| **API** | **Yes** (REST API, OpenAPI compliant) |
| **Download** | CSV, SDF, SQL dump; monthly Zenodo releases |
| **Priority** | Tier 1 (compound structures) |
| **Storage Estimate** | ~2 GB (full database) |

### PubChem

| Field | Value |
|-------|-------|
| **URL** | https://pubchem.ncbi.nlm.nih.gov/ |
| **Content** | World's largest chemical information collection |
| **License** | Public domain |
| **API** | **Yes** (programmatic access) |
| **Download** | FTP bulk download (millions of structures) |
| **Priority** | Tier 1 (ID mapping) |

### ChEMBL

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/chembl/ |
| **Content** | 1.7M compounds, 14M activities, 24K+ natural products |
| **License** | Open access (academic) |
| **API** | **Yes** (REST API with Python client) |
| **Download** | FTP (RDF, flatfiles) |
| **Priority** | Tier 1 (bioactivity data) |

---

## Comparison Matrix

### Coverage Comparison

| Database | Formulas | Herbs | Compounds | Targets | Kampo-Specific |
|----------|----------|-------|-----------|---------|----------------|
| KampoDB | 298 | 180 | 3,002 | 62,906 | **Yes** |
| STORK | 148 | - | - | - | **Yes** (official) |
| KNApSAcK KAMPO | 336 | 278 | - | - | **Yes** |
| TradMPD | - | 80 | 80 | - | **Yes** |
| TM-MC 2.0 | 5,075 | 635 (192 JP) | 34,107 | 13,992 | Partial |
| BATMAN-TCM 2.0 | 54,832 | 8,404 | 39,171 | 2.3M TTIs | No (TCM) |
| TCMBank | - | 9,192 | 61,966 | 15,179 | No (TCM) |
| HERB 2.0 | - | 7,263 | 49,258 | 12,933 | No (TCM) |

### Access Method Comparison

| Database | Web | Download | API | Commercial OK |
|----------|-----|----------|-----|---------------|
| KampoDB | Yes | Limited | No | **Yes** (CC BY-SA) |
| STORK | Yes | No | No | **Yes** |
| KNApSAcK | Yes | Partial | URL patterns | No |
| TradMPD | Yes | No | No | Academic only |
| TM-MC 2.0 | Yes | **Yes** | No | Open |
| BATMAN-TCM 2.0 | Yes | Yes | **Yes** | Academic |
| TCMBank | Yes | **Yes** | Partial | Non-commercial |
| COCONUT | Yes | **Yes** | **Yes** | **Yes** (CC0) |
| ChEMBL | Yes | **Yes** | **Yes** | Academic |
| PubChem | Yes | **Yes** | **Yes** | **Yes** |

---

## Schema Overview

### Core Data Model

| Entity | Description | Key Fields |
|--------|-------------|------------|
| Formula | Kampo prescription | formula_id, japanese_name, romaji, composition |
| Crude Drug | Raw medicinal material (Shoyaku) | drug_id, japanese_name, latin_name, plant_source |
| Compound | Chemical constituent | compound_id, name, smiles, inchi_key, knapsack_id |
| Target | Protein target | target_id, protein_name, uniprot_id |
| Sho | Constitutional pattern | sho_id, japanese_term, description, indications |

### Key Tables (from Primary Databases)

#### KampoDB Schema

| Table | Description | Key Fields |
|-------|-------------|------------|
| `formulas` | Kampo prescriptions | formula_id, name, romaji, crude_drug_composition |
| `crude_drugs` | Raw materials | drug_id, japanese_name, latin_name, family |
| `compounds` | Phytochemicals | compound_id, name, smiles, knapsack_id |
| `targets` | Protein targets | target_id, uniprot_id, protein_name |
| `docking_results` | Compound-target predictions | compound_id, target_id, docking_score |
| `formula_crude_drug` | Formula composition | formula_id, drug_id, amount, unit |
| `drug_compound` | Crude drug constituents | drug_id, compound_id |

#### STORK Schema (Reference)

| Table | Description | Key Fields |
|-------|-------------|------------|
| `ethical_formulas` | 148 approved Kampo drugs | formula_id, product_name, manufacturer |
| `crude_drugs` | 241 approved crude drugs | drug_id, jp_pharmacopoeia_name, specifications |
| `formula_composition` | Ingredient ratios | formula_id, drug_id, ratio, glycyrrhiza_dose |
| `package_insert` | Regulatory information | formula_id, indications, warnings, dosage |

#### TM-MC 2.0 Schema (Japanese Subset)

| Table | Description | Key Fields |
|-------|-------------|------------|
| `medicinal_materials` | JP Pharmacopoeia materials (192) | material_id, jp_name, latin_name |
| `compounds` | Chemical constituents | compound_id, smiles, pubchem_cid |
| `targets` | Predicted targets | target_id, uniprot_id |
| `diseases` | Disease associations | disease_id, name, icd_code |
| `prescriptions` | Traditional formulas | prescription_id, name, source |

#### KNApSAcK KAMPO Schema

| Table | Description | Key Fields |
|-------|-------------|------------|
| `kampo_formulas` | 336 Kampo prescriptions | formula_id, name, composition |
| `medicinal_plants` | 278 source plants | plant_id, species_name, family |
| `metabolites` | Chemical compounds | c_id, name, molecular_formula, smiles |
| `formula_plant` | Formula-plant relationships | formula_id, plant_id |
| `plant_metabolite` | Plant-compound relationships | plant_id, c_id |

### Identifier Cross-References

| Database | Formula ID | Crude Drug ID | Compound ID | Target ID |
|----------|------------|---------------|-------------|-----------|
| KampoDB | Internal | Internal | KNApSAcK C_ID | UniProt |
| STORK | NPN | JP Pharmacopoeia | N/A | N/A |
| TM-MC 2.0 | Internal | Internal | PubChem CID | UniProt |
| KNApSAcK | Internal | Species name | C_ID | N/A |
| BATMAN-TCM | TCM equivalent | TCM equivalent | PubChem CID | UniProt |

### Common Data Formats

| Format | Use Case | Source Databases |
|--------|----------|------------------|
| Web tables | Manual extraction | KampoDB, STORK |
| JSON | API queries | BATMAN-TCM (supplementary) |
| SDF | Structure downloads | TM-MC 2.0, PubChem |
| TSV | Bulk downloads | TM-MC 2.0 |
| MediaWiki | Community data | Metabolomics.jp Wiki |

### Kampo-Specific Attributes

| Attribute | Description | Values/Format |
|-----------|-------------|---------------|
| `sho` | Constitutional pattern | Example: Kyo-sho (deficiency), Jitsu-sho (excess) |
| `glycyrrhiza_dose` | Licorice content (safety) | mg per daily dose |
| `ethical_status` | Regulatory approval | Approved/Not approved for prescription use |
| `source_text` | Classical reference | Shang Han Lun, Jin Gui Yao Lue, etc. |

---

## Integration Recommendations

### Data Integration Strategy

```
Layer 1: Kampo Formulas (STORK + KampoDB)
    |
Layer 2: Crude Drugs/Herbs (TM-MC 2.0 + KNApSAcK)
    |
Layer 3: Compounds (TCMBank + COCONUT + PubChem)
    |
Layer 4: Targets (BATMAN-TCM 2.0 + KampoDB + ChEMBL)
    |
Layer 5: Pathways/Diseases (KEGG + DisGeNET)
```

### Priority Integration Order

**Phase 1 (MVP):**
1. **STORK** - Reference for all 148 approved Kampo formulas (web scraping)
2. **KampoDB** - Formula-compound-target relationships (contact maintainers for bulk access)
3. **TM-MC 2.0** - 192 JP herbs with compound/target data (direct download)

**Phase 2 (Post-MVP):**
4. **BATMAN-TCM 2.0** - Large-scale ingredient-target interactions (REST API)
5. **KNApSAcK KAMPO** - Metabolite-species relationships (URL pattern extraction)
6. **EKAT** - Clinical evidence abstracts

**Phase 3 (Enhancement):**
7. **TCMBank** - Comprehensive compound/target mappings
8. **COCONUT** - CC0 chemical structures
9. **TradMPD** - Experimental activity data

### ID Mapping Strategy

| Source | Primary ID | Cross-references |
|--------|------------|------------------|
| KampoDB | Internal | KNApSAcK, PubChem |
| KNApSAcK | C_ID | CAS, InChI, SMILES |
| PubChem | CID | CAS, InChI, SMILES |
| ChEMBL | ChEMBL_ID | PubChem, InChI |
| UniProt | UniProt ID | NCBI Gene, Ensembl |

### Recommended Data Formats

- **Compounds**: SDF, SMILES, InChI (from PubChem/COCONUT)
- **Targets**: UniProt IDs (for gene integration)
- **Interactions**: JSON (from BATMAN-TCM API)
- **Pathways**: KEGG IDs

---

## Storage Estimates Summary

| Category | Database(s) | Estimated Size |
|----------|-------------|----------------|
| Core Kampo | KampoDB, STORK, TradMPD | ~75 MB |
| Metabolite-Species | KNApSAcK KAMPO | ~30 MB |
| Regional Integration | TM-MC 2.0 | ~100 MB |
| TCM Cross-Reference | BATMAN-TCM, HERB, TCMBank | ~1 GB |
| Chemical Structures | COCONUT subset | ~200 MB |
| Clinical Evidence | EKAT | ~10 MB |
| **Total Kampo Subsystem** | - | **~1.5 GB** |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [tcm.md](./tcm.md) | Cross-references for TCM-derived compounds |
| [natural-products.md](./../compounds/natural-products.md) | COCONUT, chemical structure data |

---

## Technical Notes

### Data Access Limitations

1. **KampoDB**: Bulk download not documented; web scraping or maintainer contact required
2. **STORK**: Reference-only; requires web scraping for systematic extraction
3. **KNApSAcK**: Non-commercial restriction; contact required for commercial use
4. **ETHMEDmmm**: Authorization required for imagery data
5. **TradMPD**: No bulk download; academic use only

### API Availability

| Database | API Type | Rate Limit | Notes |
|----------|----------|------------|-------|
| BATMAN-TCM 2.0 | REST/JSON | ~3 sec/query | Best programmatic access for targets |
| KNApSAcK | URL patterns | Unknown | Semi-structured URL-based queries |
| COCONUT | REST/OpenAPI | Unknown | Best for chemical structures |
| PubChem | E-utilities | 3 req/sec | Standard NCBI limits |
| ChEMBL | REST | 1000/page | Pagination with caching |

---

## References

### Key Publications

1. Sawada R, et al. (2018) "KampoDB, database of predicted targets and functional annotations of natural medicines." *Scientific Reports* 8:11216.

2. Liu Z, et al. (2024) "BATMAN-TCM 2.0: an enhanced integrative database for known and predicted interactions between traditional Chinese medicine ingredients and target proteins." *Nucleic Acids Research* 52(D1):D1110-D1116.

3. Kim SK, et al. (2024) "TM-MC 2.0: an enhanced chemical database of medicinal materials in Northeast Asian traditional medicine." *BMC Complementary Medicine and Therapies* 24:37.

4. Fang S, et al. (2021) "HERB: a high-throughput experiment- and reference-guided database of traditional Chinese medicine." *Nucleic Acids Research* 49(D1):D1197-D1206.

### URL Summary

| Database | URL |
|----------|-----|
| KampoDB | https://wakanmoview.inm.u-toyama.ac.jp/kampo/ |
| WAKANYAKU Portal | https://www.inm.u-toyama.ac.jp/en/database/ |
| TradMPD | https://dentomed.toyama-wakan.net/index_en |
| STORK | http://mpdb.nibiohn.go.jp/stork/ |
| KNApSAcK | http://www.knapsackfamily.com/knapsack_core/top.php |
| TM-MC 2.0 | https://tm-mc.kr |
| EKAT | http://www.jsom.or.jp/medical/ebm/ere/index.html |
| BATMAN-TCM 2.0 | http://bionet.ncpsb.org/batman-tcm |
| TCMBank | https://TCMBank.cn/ |
| HERB | http://herb.ac.cn/ |
| SymMap | http://www.symmap.org/ |
| COCONUT | https://coconut.naturalproducts.net |
| PubChem | https://pubchem.ncbi.nlm.nih.gov/ |
| ChEMBL | https://www.ebi.ac.uk/chembl/ |

---

## Glossary

### Core Kampo Concepts

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Kampo` | Japanese traditional medicine derived from ancient Chinese medicine, systematized in Japan | Literally "Chinese method" (Kan = Chinese, Po = method) |
| `Sho` | Constitutional pattern or symptom-sign complex guiding formula selection | Similar to TCM Zheng but uniquely Japanese interpretation |
| `Kyo-sho` | Deficiency pattern; weak constitution requiring tonification | Treated with warming, nourishing formulas |
| `Jitsu-sho` | Excess pattern; robust constitution with pathogenic factors | Treated with draining, clearing formulas |
| `Chu-kan-sho` | Intermediate pattern; between deficiency and excess | Requires balanced approach |

### Formula Classification

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Shoyaku` | Crude drug; raw medicinal material | Individual herb components |
| `Yakumi` | Medicinal taste; pharmacological property | Similar to TCM five flavors |
| `Hozai` | Tonic formula | Strengthening formulas |
| `Shiazai` | Purging/draining formula | Clearing excess |
| `Wazai` | Harmonizing formula | Balancing approach |

### Diagnostic Terms

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Fukushin` | Abdominal diagnosis; unique to Kampo | Palpation patterns guide formula selection |
| `Myakushin` | Pulse diagnosis | Assess Qi, Blood, constitution |
| `Zetsushin` | Tongue diagnosis | Coating, color, shape analysis |
| `Monshin` | Inquiry; patient questioning | Symptoms, history, preferences |
| `Boshin` | Inspection; visual observation | Complexion, posture, demeanor |

### Constitutional Types

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Taishitsu` | Constitution; inherent body type | Determines formula suitability |
| `Kiketsu-sui` | Qi-Blood-Water framework | Three vital substances in Kampo theory |
| `Ki` | Vital energy (equivalent to Qi) | Ki deficiency, Ki stagnation |
| `Ketsu` | Blood in Kampo context | Ketsu stasis, Ketsu deficiency |
| `Sui` | Water/fluid metabolism | Sui stagnation (edema, dampness) |

### Common Formula Names (Romanized)

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Kakkonto` | Pueraria Combination | Cold/flu, stiff shoulders |
| `Shosaikoto` | Minor Bupleurum Combination | Liver conditions, alternating fever/chills |
| `Hochuekkito` | Tonify the Middle and Augment Qi | Fatigue, poor appetite, prolapse |
| `Juzentaihoto` | Ten Perfect Great Supplement | Severe fatigue, post-illness recovery |
| `Rikkunshito` | Six Gentlemen Combination | Digestive weakness, GERD |
| `Daikenchuto` | Major Construct the Middle | Abdominal cold, intestinal obstruction |

### Preparation Types

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Ekisu-zai` | Extract preparation; modern granule form | Concentrated powder from decoction |
| `Tou-zai` | Decoction; traditional boiled preparation | Classical preparation method |
| `Maru-zai` | Pill form | Traditional pill preparation |
| `San-zai` | Powder form | Fine powder for ingestion |

### Regulatory Terms

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Ethical Kampo` | Prescription Kampo drugs covered by national health insurance | 148 approved formulas |
| `OTC Kampo` | Over-the-counter Kampo products | Available without prescription |
| `Japanese Pharmacopoeia` (JP) | Official compendium of drug standards | Quality specifications |
| `Glycyrrhiza content` | Licorice dose tracked for safety | Monitored due to pseudoaldosteronism risk |

### Database and Technical Terms

| Term | Definition | Example/Context |
|------|------------|-----------------|
| `Docking simulation` | Computational prediction of compound-protein binding | KampoDB methodology |
| `KNApSAcK` | Species-metabolite relationship database | Compound identification system |
| `C_ID` | KNApSAcK compound identifier | Cross-reference ID |
| `InChIKey` | International Chemical Identifier hash | Unique compound fingerprint |
| `SMILES` | Simplified Molecular Input Line Entry System | Chemical structure notation |
| `UniProt` | Universal Protein Resource | Target protein database |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| KampoDB | Kampo Database | Primary Kampo compound-target database |
| STORK | Standards of Reporting Kampo Products | 148 approved formula reference |
| TradMPD | Traditional Medical & Pharmaceutical Database | Experimental Kampo data |
| TM-MC | Traditional Medicine Molecular Chemistry | Northeast Asian traditional medicine |
| EKAT | Evidence Reports of Kampo Treatment | Clinical evidence database |
| KNApSAcK | Species-metabolite relationship database | Metabolomics resource |
| NIBIOHN | National Institutes of Biomedical Innovation, Health and Nutrition | Japanese research institute |
| WAKANYAKU | Japanese/Chinese medicine | Portal umbrella term |
| JP | Japanese Pharmacopoeia | Quality standards |
| NPN | Natural Product Number | Health Canada license ID |
| BATMAN-TCM | Bioinformatics Analysis Tool for Molecular mechANism of TCM | Supplementary target data |
| CC BY-SA 4.0 | Creative Commons Attribution ShareAlike 4.0 | Open license |
| CC0 | Creative Commons Zero | Public domain dedication |
| REST API | Representational State Transfer API | Web service interface |
| CADD | Computer-Aided Drug Design | Computational drug discovery |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial document migrated from research.old/interventions-kampo.md |
