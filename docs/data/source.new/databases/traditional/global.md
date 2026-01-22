---
id: traditional-global
title: Global Traditional Medicine Databases (African & Latin American)
world: 2
category: traditional
tier: 3
subcategory: global
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [traditional, global, african, latin-american]
---

# Global Traditional Medicine Databases (African & Latin American)

**Document ID:** 43-25-GLOBAL-TRADITIONAL
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

Comprehensive inventory of 15 databases covering African Traditional Medicine (6 databases) and Latin American Traditional Medicine (5 databases), plus 4 global multi-regional natural product databases. Primary sources include ANPDB (11,063 African compounds with REST API), COCONUT 2.0 (695K compounds, CC0 license), and NuBBEDB (Brazilian biodiversity). SANCDB provides unique commercial analog links for African compounds. Estimated total storage: ~25 GB with CC0/open licenses for major resources.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary African source | ANPDB | Largest African NP repository (11K molecules), REST API available | Jan 2026 |
| African physical library | p-ANAPL | 500+ physical samples available for screening | Jan 2026 |
| Primary Brazilian source | NuBBEDB | First and most comprehensive Brazilian NP database | Jan 2026 |
| Open license priority | COCONUT 2.0 | 695K compounds, CC0 license, includes UNPD data | Jan 2026 |
| Linked data source | LOTUS | Wikidata integration, community curation, SPARQL endpoint | Jan 2026 |
| Activity-based source | NPASS | Quantitative IC50/Ki/EC50/MIC values (not just binary) | Jan 2026 |
| South African focus | SANCDB | REST API + commercial analog links (Mcule, MolPort) | Jan 2026 |

---

## Database Catalog

### Category 1: African Traditional Medicine

#### ANPDB - African Natural Products Database (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://african-compounds.org/anpdb/ or https://phabidb.vm.uni-freiburg.de/anpdb/ |
| **Content** | Natural products from African medicinal plants |
| **Records** | 11,063 molecules, 6 kingdoms, 264 families, 1,757 species, 2,148 references (1961-2024) |
| **License** | Free for academic use |
| **API** | RESTful API available |
| **Data Format** | Web interface, SDF downloads, API access |
| **Update Frequency** | 2024 (published in Nucleic Acids Research) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB |

**Key Features:**
- Largest curated African compound repository
- Integrates experimental and ML-predicted NMR/MS data
- RESTful API for programmatic access
- SDF structure downloads

**Unique Value:** Only African NP database with REST API and predicted spectroscopic data for dereplication.

---

#### SANCDB - South African Natural Compounds Database

| Field | Value |
|-------|-------|
| **URL** | https://sancdb.rubi.ru.ac.za/ |
| **Content** | Natural compounds from South African plants and marine life |
| **Records** | 1,012 compounds |
| **License** | Free for academic use |
| **API** | REST API documented |
| **Data Format** | Web interface, API, SDF downloads |
| **Update Frequency** | 2021 (added 412 compounds, commercial analog links) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~200 MB |

**Key Features:**
- REST API documented and functional
- Links to Mcule and MolPort for commercial analogs
- South African biodiversity focus
- Integration with commercial sourcing

**Unique Value:** Only web-based NP database in Africa with commercial analog links for procurement.

---

#### p-ANAPL - Pan-African Natural Products Library

| Field | Value |
|-------|-------|
| **URL** | University of Botswana - Cesriki |
| **Content** | Physical compound samples from African medicinal plants |
| **Records** | 500+ physical compound samples, ~30 compound classes |
| **License** | Collaborative access |
| **API** | Virtual library available for computational screening |
| **Data Format** | Physical samples, virtual library for docking |
| **Update Frequency** | Established 2009 (consortium from 6 African countries) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB (virtual library) |

**Key Features:**
- Physical library available for screening (unique among African databases)
- Consortium from 6 African countries
- Virtual library for computational screening
- ~30 compound classes represented

**Unique Value:** Largest physical collection of African medicinal plant compounds available for HTS screening.

---

#### ETM-DB - Ethiopian Traditional Medicine Database

| Field | Value |
|-------|-------|
| **URL** | Regional database |
| **Content** | Ethiopian medicinal herbs, compounds, and targets |
| **Records** | 1,054 Ethiopian medicinal herbs, 16,426 herb-compound relations, 162,632 compound-target relations |
| **License** | Academic use |
| **API** | None documented |
| **Data Format** | Web interface |
| **Update Frequency** | 2019 publication |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |

**Key Features:**
- Ethiopian medicinal plant focus
- Herb-compound-target relationship mapping
- Regional ethnobotanical coverage

**Unique Value:** Comprehensive Ethiopian traditional medicine coverage with target predictions.

---

#### AfroDb - Potent African Natural Products

| Field | Value |
|-------|-------|
| **URL** | Research database |
| **Content** | Curated potent compounds from African sources |
| **Records** | Select highly potent natural products |
| **License** | Academic use |
| **API** | None documented |
| **Data Format** | Research data |
| **Update Frequency** | 2013 publication |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- Focus on bioactive compounds with demonstrated potency
- Quality-filtered African NPs
- Lead compound discovery focus

**Unique Value:** Pre-filtered African compounds with demonstrated bioactivity for lead discovery.

---

#### NANPDB - Northern African Natural Products Database (MERGED)

| Field | Value |
|-------|-------|
| **URL** | Now merged into ANPDB |
| **Content** | Natural products from Northern African sources |
| **Status** | Merged with EANPDB to form ANPDB |
| **Priority** | N/A - Use ANPDB |
| **Publication** | Ntie-Kang F, et al. J Nat Prod. 2017;80(7):2067-2076 |

**Note:** Historical database now integrated into ANPDB. Use ANPDB for current data access.

---

### Category 2: Latin American Traditional Medicine

#### NuBBEDB - Brazilian Natural Products Database (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | http://nubbe.iq.unesp.br/nubbedb.html |
| **Content** | Compounds from Brazilian biodiversity |
| **Records** | 2,000+ compounds |
| **License** | Free for academic research |
| **API** | Not documented |
| **Data Format** | Web interface, SDF downloads |
| **Update Frequency** | Established 2013, active |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |

**Key Features:**
- First Brazilian biodiversity NP library
- ADMET property predictions included
- Pharmacological and toxicological data
- SDF structure downloads available

**Unique Value:** First and most comprehensive Brazilian NP database with drug-likeness assessments.

---

#### NuBBE[KG] - NuBBE Knowledge Graph

| Field | Value |
|-------|-------|
| **URL** | https://nubbekg.aksw.org/ |
| **Content** | Biochemical knowledge graph based on NuBBE database |
| **Records** | Knowledge graph representation of NuBBEDB |
| **License** | Open access |
| **API** | SPARQL endpoint |
| **Data Format** | RDF, SPARQL queries |
| **Update Frequency** | Active |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~200 MB |

**Key Features:**
- Semantic web format
- Linked data integration
- SPARQL endpoint for queries
- RDF data model

**Unique Value:** Only Brazilian NP database with semantic web/linked data capabilities.

---

#### BIOFACQUIM - Mexican Natural Products Database

| Field | Value |
|-------|-------|
| **URL** | https://www.biofacquim.unam.mx/ |
| **Content** | Mexican natural products |
| **Records** | Mexican biodiversity compounds |
| **License** | Academic use |
| **API** | Not documented |
| **Data Format** | Web interface |
| **Update Frequency** | Active (2019 publication) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~300 MB |

**Key Features:**
- Mexican biodiversity focus
- UNAM institutional support
- Chemical diversity analysis

**Unique Value:** Primary resource for Mexican traditional medicine compounds.

---

#### MAMPDB - Mesoamerican Medicinal Plant Database

| Field | Value |
|-------|-------|
| **URL** | Regional database |
| **Content** | Ethnobotanical data from Mexico and Central America |
| **Records** | Ethnopharmacological studies 1975-2016 |
| **License** | Academic use |
| **API** | None documented |
| **Data Format** | Literature compilation |
| **Update Frequency** | Static (study compilation) |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~100 MB |

**Key Features:**
- Integration focus for national health systems
- 40+ years of ethnopharmacological studies
- Mesoamerican coverage (Mexico + Central America)

**Unique Value:** Systematic compilation of Mesoamerican ethnopharmacological research.

---

#### Oaxaca Medicinal Plants Inventory

| Field | Value |
|-------|-------|
| **URL** | Publication-based |
| **Content** | Medicinal plants from Oaxaca, Mexico |
| **Records** | 1,032 medicinal plants from 164 families in 84 municipalities |
| **License** | Open access (Ethnobiomed 2020) |
| **API** | None |
| **Data Format** | Publication data |
| **Update Frequency** | 2020 publication |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~20 MB |

**Key Features:**
- 8 ethnic groups documented
- Top families: Asteraceae, Fabaceae, Rubiaceae
- Municipal-level distribution data
- Ethnobotanical documentation

**Unique Value:** Detailed regional Mexican ethnobotany with multi-ethnic documentation.

---

### Category 3: Global Multi-Regional Natural Product Databases

#### COCONUT 2.0 - Collection of Open Natural Products (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://coconut.naturalproducts.net/ |
| **GitHub** | https://github.com/Steinbeck-Lab/coconut |
| **Content** | Aggregated natural products from 63 data collections |
| **Records** | 695,133 unique structures |
| **License** | CC0 (public domain) |
| **API** | REST API (OpenAPI compliant) - https://coconut.naturalproducts.net/api-documentation |
| **Data Format** | Web, API, SDF, CSV, database dump |
| **Update Frequency** | COCONUT 2.0 released September 2024 |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~15 GB |

**Key Features:**
- Largest open NP database (695K unique structures)
- CC0 license allows unrestricted commercial use
- REST API (OpenAPI compliant)
- Integrates 63 source databases including UNPD data
- FAIR principles compliant

**Unique Value:** Most permissive license (CC0) + largest compound collection + excellent API.

---

#### LOTUS - Linking Open Natural Products Datasets

| Field | Value |
|-------|-------|
| **URL** | https://lotus.naturalproducts.net/ and Wikidata |
| **Content** | Referenced structure-organism pairs |
| **Records** | 750,000+ structure-organism pairs |
| **License** | CC0 (via Wikidata) |
| **API** | SPARQL endpoint via Wikidata |
| **Data Format** | Wikidata, SPARQL, web interface |
| **Update Frequency** | Active (community curated) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~3 GB |

**Key Features:**
- Full Wikidata integration
- Community curation model
- SPARQL endpoint for semantic queries
- Referenced structure-organism pairs
- Open science principles

**Unique Value:** Only NP database fully integrated with Wikidata; enables linked data applications.

---

#### NPASS - Natural Product Activity and Species Source

| Field | Value |
|-------|-------|
| **URL** | http://bidd.group/NPASS |
| **Content** | Natural products with quantitative activity values |
| **Records** | 94,400 NPs, 7,700 targets, 31,600 species sources |
| **License** | Free for academic use |
| **API** | Web interface |
| **Data Format** | Web, downloadable |
| **Update Frequency** | 2023 update |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 GB |

**Key Features:**
- Quantitative activity values (IC50, Ki, EC50, MIC)
- Not just binary activity - actual potency data
- ADMET predictions included
- Species source tracking

**Unique Value:** Focus on quantitative activity values rather than binary yes/no activity data.

---

#### SuperNatural 3.0

| Field | Value |
|-------|-------|
| **URL** | http://bioinf-applied.charite.de/supernatural_3/ |
| **Content** | Natural product database with mechanism annotations |
| **Records** | 449,058 unique compounds (790,096 including isomers) |
| **License** | Free access (no registration required) |
| **API** | Not documented |
| **Data Format** | Web interface |
| **Update Frequency** | SuperNatural 3.0 released 2022 |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 GB |

**Key Features:**
- Mechanism of action annotations
- Pathway associations
- Toxicity predictions
- Taste prediction (unique feature)
- Disease-specific searches

**Unique Value:** Broad annotations including mechanism, pathway, toxicity, and taste data.

---

## Comparative Analysis

### Regional Coverage Comparison

| Database | Region | Compounds | Species/Herbs | Targets | API |
|----------|--------|-----------|---------------|---------|-----|
| **African** |
| ANPDB | Pan-Africa | 11,063 | 1,757 | - | Yes (REST) |
| SANCDB | South Africa | 1,012 | - | - | Yes (REST) |
| p-ANAPL | Pan-Africa | 500+ | Physical | - | Virtual only |
| ETM-DB | Ethiopia | - | 1,054 | 162,632 | No |
| **Latin American** |
| NuBBEDB | Brazil | 2,000+ | - | - | No |
| BIOFACQUIM | Mexico | - | - | - | No |
| NuBBE[KG] | Brazil | KG | - | - | SPARQL |
| **Global** |
| COCONUT 2.0 | Global | 695,133 | - | - | Yes (REST) |
| LOTUS | Global | 750,000+ | - | - | SPARQL |
| NPASS | Global | 94,400 | 31,600 | 7,700 | No |
| SuperNatural 3.0 | Global | 449,058 | - | - | No |

### License Summary

| License Type | Databases |
|--------------|-----------|
| **CC0 (Public Domain)** | COCONUT 2.0, LOTUS |
| **Academic Free** | ANPDB, SANCDB, NuBBEDB, BIOFACQUIM, NPASS, ETM-DB |
| **Free (No Registration)** | SuperNatural 3.0 |
| **Collaborative** | p-ANAPL |

### API Availability

| Database | API Type | Documentation |
|----------|----------|---------------|
| COCONUT 2.0 | REST (OpenAPI) | https://coconut.naturalproducts.net/api-documentation |
| ANPDB | REST | Available |
| SANCDB | REST | Documented |
| LOTUS | SPARQL | Via Wikidata |
| NuBBE[KG] | SPARQL | Available |

---

## Storage Summary

| Category | Sources | Est. Size |
|----------|---------|-----------|
| African Traditional Medicine | 5 | ~3 GB |
| Latin American Traditional Medicine | 5 | ~1 GB |
| Global Multi-Regional | 4 | ~25 GB |
| **Total** | **14** | **~29 GB** |

*Note: COCONUT 2.0 accounts for majority of storage; deduplicated estimate with overlap: ~25 GB*

---

## Integration Recommendations

### Tier 1 - Essential (High Value, Open/API Access)

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | COCONUT 2.0 | Largest open NP database, CC0 license, REST API |
| 2 | ANPDB | Best African coverage, REST API, predicted spectroscopy |
| 3 | SANCDB | REST API, commercial analog links for procurement |
| 4 | LOTUS | Wikidata integration, linked data, CC0 license |
| 5 | NuBBEDB | Brazilian coverage, ADMET predictions |

### Tier 2 - Important

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | NPASS | Quantitative activity data (IC50, Ki values) |
| 2 | NuBBE[KG] | Semantic web integration, SPARQL queries |
| 3 | p-ANAPL | Physical screening library access |
| 4 | ETM-DB | Ethiopian coverage with target predictions |
| 5 | BIOFACQUIM | Mexican biodiversity coverage |

### Tier 3 - Supplementary

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | SuperNatural 3.0 | Mechanism/toxicity annotations |
| 2 | AfroDb | Pre-filtered potent African compounds |
| 3 | MAMPDB | Mesoamerican ethnopharmacology compilation |
| 4 | Oaxaca Inventory | Regional Mexican ethnobotany |

### Integration Architecture

1. **API-First Approach**: Start with COCONUT 2.0, ANPDB, SANCDB for programmatic access
2. **Bulk Download Strategy**: COCONUT (SDF/CSV), NuBBEDB for local database
3. **Linked Data Integration**: LOTUS/Wikidata SPARQL for semantic queries
4. **Commercial Sourcing**: Use SANCDB Mcule/MolPort links for compound procurement

---

## Data Harmonization

### Identifier Mapping

| Source Identifier | Target Standard |
|-------------------|-----------------|
| Compound names | InChI / InChIKey |
| Structures | SMILES (canonical) |
| Species | NCBI Taxonomy ID |
| Geographic origin | ISO 3166 country codes |

### Regional Coverage Gaps

| Region | Current Coverage | Gap Assessment |
|--------|------------------|----------------|
| Sub-Saharan Africa | Good (ANPDB, SANCDB, ETM-DB) | East Africa underrepresented |
| North Africa | Moderate (merged into ANPDB) | Adequate |
| South America | Moderate (NuBBEDB for Brazil) | Argentina, Chile, Peru gaps |
| Central America | Limited (MAMPDB, Oaxaca) | Guatemala, Honduras gaps |
| Caribbean | Limited | Significant gap |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [index.md](./../index.md) | Parent index |
| [tcm.md](./tcm.md) | TCM databases (complementary) |
| [ayurveda.md](./ayurveda.md) | Ayurveda databases (complementary) |
| [kampo.md](./kampo.md) | Kampo databases (complementary) |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `natural product` | Chemical compound produced by living organisms | Artemisinin from Artemisia annua |
| `SMILES` | Simplified Molecular Input Line Entry System - text chemical notation | CC(=O)OC1=CC=CC=C1C(=O)O |
| `InChI` | International Chemical Identifier - standardized structure representation | InChI=1S/C9H8O4/... |
| `InChIKey` | Fixed-length hash of InChI for database lookups | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| `ethnobotany` | Study of traditional plant use by indigenous cultures | Oaxaca medicinal plants |
| `SDF` | Structure Data File - chemical format containing 3D coordinates | Molecular structure file |
| `HTS` | High-Throughput Screening - automated large-scale compound testing | p-ANAPL physical library |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `ANPDB` | African Natural Products Database - 11K+ African compounds with REST API | African TM |
| `SANCDB` | South African Natural Compounds Database with commercial analog links | South African TM |
| `COCONUT` | Collection of Open Natural Products - 695K compounds with CC0 license | Global NPs |
| `LOTUS` | Linked Open Unified Taxonomic/occurrence Score via Wikidata | Linked data |
| `NPASS` | Natural Product Activity and Species Source with quantitative activity | Bioactivity data |
| `NuBBEDB` | Brazilian Natural Products Database with ADMET predictions | Brazilian TM |
| `p-ANAPL` | Pan-African Natural Products Library - physical compound samples | Screening library |
| `SPARQL` | Query language for RDF/linked data (Wikidata, LOTUS) | Semantic queries |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NP | Natural Product | Compound from organism |
| TM | Traditional Medicine | Ethnopharmacology |
| HTS | High-Throughput Screening | Automated testing |
| SDF | Structure Data File | Chemical format |
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity | Drug-likeness |
| API | Application Programming Interface | Programmatic access |
| REST | Representational State Transfer | API architecture |
| SPARQL | SPARQL Protocol and RDF Query Language | Linked data queries |
| CC0 | Creative Commons Zero | Public domain |
| FAIR | Findable, Accessible, Interoperable, Reusable | Data principles |
| IC50 | Half-maximal Inhibitory Concentration | Potency measure |
| Ki | Inhibition constant | Binding affinity |

---

## Open Questions

- [ ] p-ANAPL access agreement - what are collaboration requirements?
- [ ] Caribbean traditional medicine coverage - identify potential sources
- [ ] South American coverage expansion - Argentina, Chile, Peru databases?
- [ ] ANPDB commercial licensing terms - contact for clarification
- [ ] NuBBEDB API development - any plans for programmatic access?

---

## References

1. ANPDB (2024). African Natural Products Database. Nucleic Acids Res. (advance article)

2. Hatherley R, et al. (2021). SANCDB update. J Cheminform. 2021.

3. Ntie-Kang F, et al. (2014). p-ANAPL: Pan-African Natural Products Library. PLoS ONE. 9(3):e90655.

4. Valli M, et al. (2013). NuBBEDB. J Nat Prod. 76(3):439-444.

5. Chandrasekhar V, et al. (2024). COCONUT 2.0. Nucleic Acids Res. 53(D1):D634.

6. Rutz A, et al. (2022). LOTUS. eLife. 11:e70780.

7. Zeng X, et al. (2023). NPASS 2023. Nucleic Acids Res. 51(D1):D621-D628.

8. Gallo K, et al. (2023). SuperNatural 3.0. Nucleic Acids Res. 51(D1):D654-D659.

9. Pilon-Jimenez BA, et al. (2019). BIOFACQUIM. Biomolecules. 9(1):31.

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial creation from research.old/data-sources-traditional-medicine-expanded.md |
