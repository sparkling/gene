# Traditional Medicine and Herbal Databases - Expanded Research

**Research Date:** January 2026
**Purpose:** Gene Platform data source expansion beyond IMPPAT, KampoDB

---

## Table of Contents

1. [Additional TCM Databases](#1-additional-tcm-databases)
2. [Ayurveda Databases Beyond IMPPAT](#2-ayurveda-databases-beyond-imppat)
3. [Western Herbal Medicine Databases](#3-western-herbal-medicine-databases)
4. [African Traditional Medicine](#4-african-traditional-medicine)
5. [Latin American Traditional Medicine](#5-latin-american-traditional-medicine)
6. [Herb-Drug Interaction Databases](#6-herb-drug-interaction-databases)
7. [Global/Multi-Regional Natural Product Databases](#7-globalmulti-regional-natural-product-databases)
8. [Summary Comparison Table](#8-summary-comparison-table)
9. [Integration Recommendations](#9-integration-recommendations)

---

## 1. Additional TCM Databases

### 1.1 SymMap

| Attribute | Details |
|-----------|---------|
| **Full Name** | SymMap: Integrative Database of TCM Enhanced by Symptom Mapping |
| **URL** | http://www.symmap.org/ or https://www.bioinfo.org/symmap |
| **Content** | 1,717 TCM symptoms, 499 herbs, 961 modern medicine symptoms, 5,235 diseases, 19,595 herbal ingredients, 4,302 target genes |
| **Key Feature** | Maps TCM symptoms to modern medicine symptoms - unique bridging approach |
| **API** | No public REST API documented |
| **License** | Free for academic use |
| **Data Format** | Web interface, downloadable tables |
| **Latest Update** | SymMap 2.0 (2020 edition based on Chinese Pharmacopoeia) |
| **Publication** | Wu Y, et al. Nucleic Acids Res. 2019;47(D1):D1110-D1117 |
| **Use Case** | Symptom-based drug discovery, TCM-modern medicine bridging |

**Unique Value:** SymMap is the only database that systematically maps TCM symptoms to modern medicine disease classifications, curated by 17 leading TCM experts.

---

### 1.2 ETCM (Encyclopedia of Traditional Chinese Medicine)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Encyclopedia of Traditional Chinese Medicine |
| **URL** | http://www.tcmip.cn/ETCM2/front/#/ (v2.0) |
| **Content** | 48,442 TCM formulas, 9,872 Chinese patent drugs, 2,079 medicinal materials, 38,298 ingredients |
| **Key Feature** | Comprehensive formula database with predicted targets and JavaScript network visualization |
| **API** | No public API; web-based systematic analysis tools |
| **License** | Free for academic research |
| **Data Format** | Web interface, network downloads |
| **Latest Update** | ETCM v2.0 (2023) |
| **Publication** | Xu HY, et al. Nucleic Acids Res. 2019;47(D1):D976-D982 |
| **Use Case** | Formula composition analysis, quality marker identification, drug repurposing |

**Unique Value:** Largest formula database with 48,442 formulas from ancient Chinese medical books plus 9,872 modern Chinese patent drugs.

---

### 1.3 TCMID (Traditional Chinese Medicine Integrated Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Traditional Chinese Medicine Integrated Database |
| **URL** | http://www.megabionet.org/tcmid/ |
| **Content** | ~47,000 prescriptions, 8,159 herbs, 25,210 compounds, 6,828 drugs, 3,791 diseases, 17,521 targets |
| **Key Feature** | Largest integrated dataset linking to DrugBank, OMIM, PubChem |
| **API** | No documented public API |
| **License** | Free for academic use |
| **Data Format** | Web interface, network visualization tool |
| **Latest Update** | TCMID 2.0 (2017) - added mass spectrometry spectra |
| **Publication** | Xue R, et al. Nucleic Acids Res. 2013;41(D1):D1089-D1095 |
| **Use Case** | Herb-disease relationship mining, combination therapy research |

**Unique Value:** Contains 778 herbal mass spectrometry spectra related to 170 herbs showing quality variation data.

---

### 1.4 TCMBank

| Attribute | Details |
|-----------|---------|
| **Full Name** | TCMBank - Largest TCM Database |
| **URL** | https://tcmbank.cn/ |
| **Content** | 9,192 herbs, 61,966 ingredients, 15,179 targets, 32,529 diseases |
| **Key Feature** | Deep learning-based Chinese-Western medicine exclusion prediction |
| **API** | Not documented |
| **License** | Non-commercial use; downloadable |
| **Data Format** | Web interface, bulk downloads available |
| **Latest Update** | 2024 (integrated into TCMM) |
| **Publication** | Huang L, et al. Signal Transduct Target Ther. 2023;8:127 |
| **Use Case** | Drug-drug interaction prediction, adverse effect prediction |

**Unique Value:** Uses graph neural networks for functional group extraction and causal learning for Chinese-Western medicine interaction prediction.

---

### 1.5 HERB (High-throughput Experiment Reference-guided Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | HERB - BenCaoZuJian |
| **URL** | http://herb.ac.cn/ |
| **Content** | 7,263 herbs, 49,258 ingredients, 12,933 targets, 28,212 diseases, 6,164 gene expression profiles |
| **Key Feature** | Pharmacotranscriptomics data with CMap drug connections |
| **API** | Not documented |
| **License** | Free for academic research |
| **Data Format** | Web interface with transcriptomic data downloads |
| **Latest Update** | HERB 2.0 (2024) - 8,558 clinical trials, 8,032 meta-analyses |
| **Publication** | Fang S, et al. Nucleic Acids Res. 2021;49(D1):D1197-D1206 |
| **Use Case** | Drug repositioning, mechanism validation via transcriptomics |

**Unique Value:** Only TCM database with systematic high-throughput transcriptomic data (2,231 experiments) and CMap drug connectivity maps.

---

### 1.6 BATMAN-TCM

| Attribute | Details |
|-----------|---------|
| **Full Name** | Bioinformatics Analysis Tool for Molecular mechANism of TCM |
| **URL** | http://bionet.ncpsb.org/batman-tcm/ |
| **Content** | Online analysis tool with ingredient-target prediction |
| **Key Feature** | Similarity-based target prediction + pathway enrichment |
| **API** | Web service for batch analysis |
| **License** | Free for academic use |
| **Data Format** | Web interface, result downloads |
| **Latest Update** | BATMAN-TCM 2.0 (2023) - bidirectional query |
| **Publication** | Liu Z, et al. Sci Rep. 2016;6:21146 |
| **Use Case** | Network pharmacology analysis, target prediction |

**Unique Value:** Analysis tool rather than static database - provides GO, KEGG, disease enrichment and network visualization.

---

### 1.7 CMAUP (Collective Molecular Activities of Useful Plants)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Collective Molecular Activities of Useful Plants |
| **URL** | https://bidd.group/CMAUP/ |
| **Content** | 60,222 ingredients, 7,865 plants, 758 targets, 1,399 diseases, 238 KEGG pathways |
| **Key Feature** | Global medicinal plant coverage (153 countries), transcriptomic disease associations |
| **API** | Download available, interactive world map |
| **License** | Free download on individual pages |
| **Data Format** | Web interface, CSV downloads |
| **Latest Update** | 2024 - added clinical trial data, DNA barcodes, oral bioavailability |
| **Publication** | Zeng X, et al. Nucleic Acids Res. 2024;52(D1):D1508-D1517 |
| **Use Case** | Global ethnobotany research, drug-producing plant identification |

**Unique Value:** Geographic coverage of 153 countries with world map navigation; includes 691 clinical trials for 185 plants.

---

## 2. Ayurveda Databases Beyond IMPPAT

### 2.1 TKDL (Traditional Knowledge Digital Library)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Traditional Knowledge Digital Library |
| **URL** | https://www.tkdl.res.in/ |
| **Content** | 418,885+ formulations (119,269 Ayurveda, 236,399 Unani, 54,689 Siddha, 4,151 Yoga, 4,377 Sowa Rigpa) |
| **Key Feature** | Largest Indian traditional medicine formulation database |
| **API** | No public API - restricted access via agreements |
| **License** | Restricted to Patent Offices via TKDL Access Agreement |
| **Data Format** | Proprietary search interface |
| **Access** | EPO, USPTO, JPO, UKPO, CIPO, DPMA, IP Australia, Indian Patent Office, Chile Patent Office |
| **Publication** | Government of India initiative |
| **Use Case** | Prior art search, biopiracy prevention |

**Limitation:** Not publicly accessible; requires formal agreement with CSIR. Primary use is patent examination.

---

### 2.2 IMPLAD (Indian Medicinal Plants Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Indian Medicinal Plants Database |
| **URL** | Part of FRLHT/TDU infrastructure |
| **Content** | 6,500 medicinal plants, 200,000 vernacular names in 32 Indian languages, GIS maps |
| **Key Feature** | Botanical names, distribution, threat status, state inventories |
| **API** | Not documented |
| **License** | Research access |
| **Data Format** | Web interface |
| **Coverage** | Flora, fauna, metals, minerals from 1500 BC to 1900 AD texts |
| **Use Case** | Ethnobotanical research, conservation assessment |

**Unique Value:** Contains vernacular names in 32 Indian languages; maps threatened medicinal plant species.

---

### 2.3 BSI Medicinal Plant Database

| Attribute | Details |
|-----------|---------|
| **Full Name** | Botanical Survey of India Medicinal Plant Database |
| **URL** | https://bsi.gov.in/page/en/medicinal-plant-database |
| **Content** | 1,915 species (Phase 1), ~1,000 more planned |
| **Key Feature** | Official government database with herbarium specimen images |
| **API** | Not available |
| **License** | Free public access |
| **Data Format** | Web interface |
| **Content Fields** | Scientific name, family, vernacular name, medicinal uses, location, herbarium images |
| **Use Case** | Species identification, medicinal use verification |

---

### 2.4 NeMedPlant

| Attribute | Details |
|-----------|---------|
| **Full Name** | Medicinal Plants of North-East India |
| **URL** | Regional database |
| **Content** | Therapeutic applications and chemical constituents from NE India |
| **Key Feature** | Regional focus on biodiversity hotspot |
| **License** | Academic use |
| **Use Case** | NE India ethnobotany research |

---

### 2.5 Ayurvedic Pharmacopoeia of India

| Attribute | Details |
|-----------|---------|
| **Full Name** | Ayurvedic Pharmacopoeia of India |
| **URL** | Available via Government of India archives |
| **Content** | Official monographs for Ayurvedic drugs |
| **Key Feature** | Quality standards, identification, therapeutic use |
| **License** | Government publication |
| **Data Format** | PDF volumes |
| **Use Case** | Quality control, authentication |

---

## 3. Western Herbal Medicine Databases

### 3.1 ESCOP Monographs

| Attribute | Details |
|-----------|---------|
| **Full Name** | European Scientific Cooperative on Phytotherapy Monographs |
| **URL** | https://www.escop.com/ |
| **Content** | 120+ evidence-based monographs |
| **Key Feature** | EU regulatory standard for phytomedicines |
| **API** | No API - monographs available for purchase |
| **License** | Commercial (monograph purchase required) |
| **Data Format** | PDF monographs |
| **Content Fields** | Composition, dosages, warnings, interactions, clinical data, toxicity |
| **Use Case** | EU herbal medicine registration, evidence assessment |

**Unique Value:** Official European standard; monographs accepted by EMA for herbal medicine registration.

---

### 3.2 PhytoNet (ESCOP Adverse Reactions Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | PhytoNet Herbal Adverse Drug Reactions Database |
| **URL** | Associated with ESCOP |
| **Content** | Adverse event reports for herbal medicines |
| **Key Feature** | Pharmacovigilance for phytomedicines |
| **Data Format** | Report submission forms, database queries |
| **Use Case** | Safety surveillance, interaction tracking |

---

### 3.3 EMA HMPC Monographs

| Attribute | Details |
|-----------|---------|
| **Full Name** | European Medicines Agency - Committee on Herbal Medicinal Products |
| **URL** | https://www.ema.europa.eu/en/human-regulatory-overview/herbal-medicinal-products |
| **Content** | Community herbal monographs and assessment reports |
| **Key Feature** | Official EU regulatory documents |
| **API** | EMA has APIs for some data |
| **License** | Free public access |
| **Data Format** | PDF documents |
| **Use Case** | EU market authorization, traditional use registration |

**Unique Value:** Official EU regulatory standard; required for herbal medicine authorization in EU.

---

### 3.4 AMED (Allied and Complementary Medicine Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Allied and Complementary Medicine Database |
| **URL** | British Library Health Care Information Service |
| **Content** | Literature database for complementary medicine |
| **Key Feature** | Bibliographic database covering herbal medicine literature |
| **License** | Subscription required |
| **Data Format** | Bibliographic records |
| **Use Case** | Literature search, evidence review |

---

## 4. African Traditional Medicine

### 4.1 ANPDB (African Natural Products Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | African Natural Products Database |
| **URL** | https://african-compounds.org/anpdb/ or https://phabidb.vm.uni-freiburg.de/anpdb/ |
| **Content** | 11,063 molecules, 6 kingdoms, 264 families, 1,757 species, 2,148 references (1961-2024) |
| **Key Feature** | Largest African NP repository with predicted NMR/MS data |
| **API** | RESTful API available |
| **License** | Free for academic use |
| **Data Format** | Web interface, SDF downloads, API access |
| **Latest Update** | 2024 (published in Nucleic Acids Research) |
| **Publication** | Nucleic Acids Res. 2024 (advance article) |
| **Use Case** | African NP drug discovery, dereplication |

**Unique Value:** Integrates experimental and ML-predicted NMR/MS data; largest curated African compound repository.

---

### 4.2 SANCDB (South African Natural Compounds Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | South African Natural Compounds Database |
| **URL** | https://sancdb.rubi.ru.ac.za/ |
| **Content** | 1,012 compounds from South African plants and marine life |
| **Key Feature** | REST API, links to Mcule/MolPort commercial analogs |
| **API** | REST API documented |
| **License** | Free for academic use |
| **Data Format** | Web interface, API, SDF downloads |
| **Latest Update** | 2021 (added 412 compounds, commercial analog links) |
| **Publication** | Hatherley R, et al. J Cheminform. 2015;7:29 |
| **Use Case** | South African biodiversity screening, analog sourcing |

**Unique Value:** Only web-based NP database in Africa with commercial analog links (Mcule, MolPort).

---

### 4.3 p-ANAPL (Pan-African Natural Products Library)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Pan-African Natural Products Library |
| **URL** | University of Botswana - Cesriki |
| **Content** | 500+ physical compound samples, ~30 compound classes |
| **Key Feature** | Physical library available for screening + virtual library |
| **API** | Virtual library available for computational screening |
| **License** | Collaborative access |
| **Data Format** | Physical samples, virtual library for docking |
| **Established** | 2009 (consortium from 6 African countries) |
| **Publication** | Ntie-Kang F, et al. PLoS ONE. 2014;9(3):e90655 |
| **Use Case** | Physical screening, virtual screening campaigns |

**Unique Value:** Largest physical collection of African medicinal plant compounds available for screening.

---

### 4.4 ETM-DB (Ethiopian Traditional Medicine Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Ethiopian Traditional Herbal Medicine and Phytochemicals Database |
| **URL** | Regional database |
| **Content** | 1,054 Ethiopian medicinal herbs, 16,426 herb-compound relations, 162,632 compound-target relations |
| **Key Feature** | Ethiopian medicinal plant focus |
| **License** | Academic use |
| **Publication** | Yigezu Y, et al. BMC Complement Med Ther. 2019 |
| **Use Case** | Ethiopian traditional medicine research |

---

### 4.5 AfroDb

| Attribute | Details |
|-----------|---------|
| **Full Name** | Select Highly Potent Natural Product Library from African Medicinal Plants |
| **URL** | Research database |
| **Content** | Curated potent compounds from African sources |
| **Key Feature** | Focus on bioactive compounds with demonstrated potency |
| **Publication** | Ntie-Kang F, et al. PLoS ONE. 2013;8(10):e78085 |
| **Use Case** | Lead compound discovery |

---

### 4.6 NANPDB (Northern African Natural Products Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Northern African Natural Products Database |
| **URL** | Now merged into ANPDB |
| **Content** | NPs from Northern African sources |
| **Status** | Merged with EANPDB to form ANPDB |
| **Publication** | Ntie-Kang F, et al. J Nat Prod. 2017;80(7):2067-2076 |

---

## 5. Latin American Traditional Medicine

### 5.1 NuBBEDB (Brazilian Natural Products Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Nuclei of Bioassays, Ecophysiology and Biosynthesis of Natural Products Database |
| **URL** | http://nubbe.iq.unesp.br/nubbedb.html |
| **Content** | 2,000+ compounds from Brazilian biodiversity |
| **Key Feature** | First Brazilian biodiversity NP library; ADMET properties |
| **API** | Not documented |
| **License** | Free for academic research |
| **Data Format** | Web interface, SDF downloads |
| **Established** | 2013 |
| **Publication** | Valli M, et al. J Nat Prod. 2013;76(3):439-444 |
| **Use Case** | Brazilian biodiversity drug discovery, dereplication |

**Unique Value:** First and most comprehensive Brazilian NP database; includes pharmacological and toxicological data.

---

### 5.2 NuBBE[KG] (Knowledge Graph)

| Attribute | Details |
|-----------|---------|
| **Full Name** | NuBBE Knowledge Graph |
| **URL** | https://nubbekg.aksw.org/ |
| **Content** | Biochemical knowledge graph based on NuBBE database |
| **Key Feature** | Semantic web format, linked data |
| **Data Format** | RDF, SPARQL endpoint |
| **Use Case** | Semantic queries, data integration |

---

### 5.3 BIOFACQUIM

| Attribute | Details |
|-----------|---------|
| **Full Name** | Mexican Compound Database of Natural Products |
| **URL** | https://www.biofacquim.unam.mx/ |
| **Content** | Mexican natural products |
| **Key Feature** | Focus on Mexican biodiversity |
| **License** | Academic use |
| **Publication** | Pilón-Jiménez BA, et al. Biomolecules. 2019;9(1):31 |
| **Use Case** | Mexican NP research |

---

### 5.4 MAMPDB (Mesoamerican Medicinal Plant Database)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Mesoamerican Medicinal Plant Database |
| **Content** | Ethnobotanical data from Mexico and Central America |
| **Key Feature** | Integration focus for national health systems |
| **Publication** | Ethnopharmacological studies 1975-2016 |
| **Use Case** | Mesoamerican traditional medicine integration |

---

### 5.5 Oaxaca Medicinal Plants Inventory

| Attribute | Details |
|-----------|---------|
| **Content** | 1,032 medicinal plants from 164 families in 84 municipalities |
| **Key Feature** | 8 ethnic groups documented |
| **Top Families** | Asteraceae, Fabaceae, Rubiaceae |
| **Publication** | Ethnobiomed 2020 |

---

## 6. Herb-Drug Interaction Databases

### 6.1 Natural Medicines Database (NatMed Pro)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Natural Medicines Comprehensive Database |
| **URL** | https://naturalmedicines.therapeuticresearch.com/ |
| **Content** | Evidence-based monographs, interaction checker |
| **Key Feature** | Multi-ingredient interaction checking, severity ratings |
| **API** | Commercial API available |
| **License** | Subscription required (commercial) |
| **Data Format** | Web interface, API |
| **Use Case** | Clinical decision support, interaction checking |

**Unique Value:** Leading clinical tool; allows checking multiple supplements against multiple drugs simultaneously.

---

### 6.2 Memorial Sloan Kettering "About Herbs"

| Attribute | Details |
|-----------|---------|
| **Full Name** | MSK About Herbs, Botanicals & Other Products |
| **URL** | https://www.mskcc.org/cancer-care/diagnosis-treatment/symptom-management/integrative-medicine/herbs |
| **Content** | 290 entries on herbs, botanicals, supplements |
| **Key Feature** | Oncology-focused; professional and consumer formats |
| **API** | No API; mobile app available |
| **License** | Free public access |
| **Data Format** | Web interface, iOS/Android app |
| **Use Case** | Oncology herb-drug interactions, patient education |

**Unique Value:** Oncology specialist focus; dual formatting for healthcare professionals and patients.

---

### 6.3 DrugBank

| Attribute | Details |
|-----------|---------|
| **Full Name** | DrugBank Knowledgebase |
| **URL** | https://go.drugbank.com/ |
| **Content** | 4,563 FDA approved drugs, 1,413,413 drug-drug interactions, 2,475 drug-food interactions |
| **Key Feature** | Clinical API, comprehensive drug data |
| **API** | Yes - Clinical API and Discovery API documented |
| **License** | Free academic access; commercial licensing available |
| **Data Format** | Web, API (JSON), bulk downloads |
| **Latest Update** | DrugBank 6.0 (2024) |
| **Publication** | Knox C, et al. Nucleic Acids Res. 2024;52(D1):D1265-D1275 |
| **Use Case** | Drug interaction checking, drug discovery |

**Note:** Limited herb-specific data; primarily drug-focused but includes drug-food interactions.

---

### 6.4 PharmMapper

| Attribute | Details |
|-----------|---------|
| **Full Name** | PharmMapper Server |
| **URL** | https://www.lilab-ecust.cn/pharmmapper/ |
| **Content** | 23,236 proteins (16,159 druggable, 51,431 ligandable pharmacophore models) |
| **Key Feature** | Reverse pharmacophore mapping for target prediction |
| **API** | Web service for batch submission |
| **License** | Free for academic use |
| **Data Format** | Web interface, downloadable results |
| **Latest Update** | 2017 (6x database expansion) |
| **Publication** | Wang X, et al. Nucleic Acids Res. 2017;45(W1):W356-W360 |
| **Use Case** | Natural product target identification, mechanism elucidation |

---

## 7. Global/Multi-Regional Natural Product Databases

### 7.1 COCONUT (Collection of Open Natural Products)

| Attribute | Details |
|-----------|---------|
| **Full Name** | COlleCtion of Open Natural prodUcTs |
| **URL** | https://coconut.naturalproducts.net/ |
| **Content** | 695,133 unique structures from 63 data collections |
| **Key Feature** | Open data (CC0), REST API, FAIR principles |
| **API** | REST API (OpenAPI compliant) - https://coconut.naturalproducts.net/api-documentation |
| **License** | CC0 (public domain) |
| **Data Format** | Web, API, SDF, CSV, database dump |
| **Latest Update** | COCONUT 2.0 (September 2024) |
| **Publication** | Chandrasekhar V, et al. Nucleic Acids Res. 2024 |
| **GitHub** | https://github.com/Steinbeck-Lab/coconut |
| **Use Case** | Universal NP screening, data integration |

**Unique Value:** Largest open NP database; CC0 license allows unrestricted use; integrates 63 source databases including UNPD data.

---

### 7.2 LOTUS (Wikidata Integration)

| Attribute | Details |
|-----------|---------|
| **Full Name** | Linking Open natural prodUcTs dataSets |
| **URL** | https://lotus.naturalproducts.net/ and Wikidata |
| **Content** | 750,000+ referenced structure-organism pairs |
| **Key Feature** | Wikidata integration, community curation |
| **API** | SPARQL endpoint via Wikidata |
| **License** | CC0 (via Wikidata) |
| **Data Format** | Wikidata, SPARQL, web interface |
| **Publication** | Rutz A, et al. eLife. 2022;11:e70780 |
| **Use Case** | Open science, linked data applications |

**Unique Value:** Only NP database fully integrated with Wikidata; enables community curation and semantic queries.

---

### 7.3 NPASS

| Attribute | Details |
|-----------|---------|
| **Full Name** | Natural Product Activity and Species Source |
| **URL** | http://bidd.group/NPASS |
| **Content** | 94,400 NPs, 7,700 targets, 31,600 species sources |
| **Key Feature** | Quantitative activity values (IC50, Ki, EC50, MIC) |
| **API** | Web interface |
| **License** | Free for academic use |
| **Data Format** | Web, downloadable |
| **Latest Update** | 2023 |
| **Publication** | Zeng X, et al. Nucleic Acids Res. 2023;51(D1):D621-D628 |
| **Use Case** | Activity-based screening, SAR studies |

**Unique Value:** Focus on quantitative activity values rather than binary activity; includes ADMET predictions.

---

### 7.4 SuperNatural 3.0

| Attribute | Details |
|-----------|---------|
| **Full Name** | SuperNatural Database |
| **URL** | http://bioinf-applied.charite.de/supernatural_3/ |
| **Content** | 449,058 unique compounds (790,096 including isomers) |
| **Key Feature** | Mechanism of action, pathways, toxicity, taste prediction |
| **API** | Not documented |
| **License** | Free access (no registration) |
| **Data Format** | Web interface |
| **Latest Update** | SuperNatural 3.0 (2022) |
| **Publication** | Gallo K, et al. Nucleic Acids Res. 2023;51(D1):D654-D659 |
| **Use Case** | NP screening, disease-specific searches |

---

### 7.5 NAPRALERT

| Attribute | Details |
|-----------|---------|
| **Full Name** | NAtural PRoducts ALERT |
| **URL** | https://pharmacognosy.pharmacy.uic.edu/napralert/ |
| **Content** | 200,000+ published studies on natural products |
| **Key Feature** | Comprehensive literature from 1975+, ethnomedicine records |
| **API** | Limited free queries; fee-based extensive searches |
| **License** | Freemium (limited free, paid for extensive) |
| **Data Format** | Web interface |
| **Established** | 1975 (University of Illinois Chicago) |
| **Use Case** | Literature search, ethnomedicine research |

**Unique Value:** World's largest herb/medicinal plant research database with 50+ years of literature coverage.

---

### 7.6 KampoDB (Japanese Kampo)

| Attribute | Details |
|-----------|---------|
| **Full Name** | KampoDB - Predicted Targets and Annotations |
| **URL** | http://wakanmoview.inm.u-toyama.ac.jp/kampo/ |
| **Content** | 158 crude drugs, 348 prescriptions, predicted targets |
| **Key Feature** | Machine learning target prediction, Japanese Pharmacopoeia alignment |
| **API** | Web interface |
| **License** | Free for academic use |
| **Publication** | Sawada R, et al. Sci Rep. 2018;8:11216 |
| **Use Case** | Kampo mechanism research, target prediction |

---

## 8. Summary Comparison Table

| Database | Region | Compounds | Plants/Herbs | Targets | API | License | Updated |
|----------|--------|-----------|--------------|---------|-----|---------|---------|
| **TCM Databases** |
| SymMap | China | 19,595 | 499 | 4,302 | No | Academic | 2020 |
| ETCM v2.0 | China | 38,298 | 2,079 | Predicted | No | Academic | 2023 |
| TCMID 2.0 | China | 25,210 | 8,159 | 17,521 | No | Academic | 2017 |
| TCMBank | China | 61,966 | 9,192 | 15,179 | No | Non-commercial | 2024 |
| HERB 2.0 | China | 49,258 | 7,263 | 12,933 | No | Academic | 2024 |
| CMAUP | Global | 60,222 | 7,865 | 758 | Download | Free | 2024 |
| **Ayurveda Databases** |
| TKDL | India | - | 418,885 formulas | - | No | Restricted | 2022 |
| IMPLAD | India | - | 6,500 | - | No | Academic | Active |
| BSI | India | - | 1,915 | - | No | Free | Active |
| **Western Herbal** |
| ESCOP | Europe | - | 120+ monographs | - | No | Commercial | Active |
| EMA HMPC | Europe | - | Official monographs | - | Some | Free | Active |
| NatMed Pro | Global | - | Comprehensive | - | Commercial | Subscription | Active |
| **African** |
| ANPDB | Africa | 11,063 | 1,757 species | - | Yes | Academic | 2024 |
| SANCDB | S. Africa | 1,012 | - | - | Yes | Academic | 2021 |
| p-ANAPL | Africa | 500+ | Physical library | - | Virtual | Collaborative | Active |
| **Latin American** |
| NuBBEDB | Brazil | 2,000+ | - | - | No | Academic | Active |
| BIOFACQUIM | Mexico | - | - | - | No | Academic | Active |
| **Global NP** |
| COCONUT 2.0 | Global | 695,133 | - | - | REST | CC0 | 2024 |
| LOTUS | Global | 750,000+ | - | - | SPARQL | CC0 | Active |
| NPASS | Global | 94,400 | 31,600 spp | 7,700 | No | Academic | 2023 |
| SuperNatural 3.0 | Global | 449,058 | - | - | No | Free | 2022 |
| **Interaction DBs** |
| DrugBank 6.0 | Global | 4,563 drugs | - | - | Yes | Freemium | 2024 |
| PharmMapper | Global | - | - | 23,236 | Web | Academic | 2017 |

---

## 9. Integration Recommendations

### Priority 1: High-Value, API-Accessible Databases

| Database | Rationale | Integration Effort |
|----------|-----------|-------------------|
| **COCONUT 2.0** | Largest open NP collection, REST API, CC0 license | Medium |
| **ANPDB** | African coverage gap, API available | Medium |
| **SANCDB** | REST API, commercial analog links | Low |

### Priority 2: Unique Content, Web Scraping Required

| Database | Rationale | Integration Effort |
|----------|-----------|-------------------|
| **HERB 2.0** | Unique transcriptomic data, clinical trials | High |
| **ETCM v2.0** | Largest formula database | High |
| **TCMBank** | DL-based interaction prediction | High |
| **SymMap** | Unique symptom mapping | Medium |
| **CMAUP** | Global coverage, clinical trials | Medium |

### Priority 3: Complementary Resources

| Database | Rationale | Integration Effort |
|----------|-----------|-------------------|
| **NPASS** | Quantitative activity data | Medium |
| **LOTUS/Wikidata** | Semantic web integration | Medium |
| **NuBBEDB** | Latin American coverage | Medium |
| **PharmMapper** | Target prediction tool | Low (tool use) |

### Integration Architecture Recommendations

1. **API-First Approach**: Start with COCONUT, STRING, SANCDB, ANPDB for programmatic access
2. **Bulk Download Strategy**: COCONUT (SDF/CSV), CMAUP, NPASS for local database
3. **Web Service Integration**: PharmMapper, BATMAN-TCM for analysis workflows
4. **Linked Data**: LOTUS/Wikidata SPARQL for semantic queries

### Data Licensing Summary

| License Type | Databases |
|--------------|-----------|
| CC0 (Public Domain) | COCONUT, LOTUS |
| Creative Commons | STRING |
| Academic Free | SymMap, ETCM, TCMID, HERB, CMAUP, ANPDB, SANCDB, NuBBEDB, NPASS |
| Restricted | TKDL (Patent Offices only) |
| Commercial | ESCOP, NatMed Pro, DrugBank (some features) |

---

## Sources

### TCM Databases
- [SymMap - Nucleic Acids Research](https://academic.oup.com/nar/article/47/D1/D1110/5150228)
- [ETCM v2.0 - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2211383523001016)
- [TCMID - Database Commons](https://ngdc.cncb.ac.cn/databasecommons/database/id/437)
- [TCMBank - Nature Signal Transduction](https://www.nature.com/articles/s41392-023-01339-1)
- [HERB - Nucleic Acids Research](https://academic.oup.com/nar/article/49/D1/D1197/6017358)
- [BATMAN-TCM - Scientific Reports](https://www.nature.com/articles/srep21146)
- [CMAUP 2024 - Nucleic Acids Research](https://academic.oup.com/nar/article/52/D1/D1508/7332076)

### Ayurveda Databases
- [TKDL Official Website](https://www.tkdl.res.in/)
- [IMPPAT - Nature Scientific Reports](https://www.nature.com/articles/s41598-018-22631-z)
- [BSI Medicinal Plant Database](https://bsi.gov.in/page/en/medicinal-plant-database)

### Western Herbal
- [ESCOP Official Website](https://www.escop.com/)
- [Natural Medicines Database](https://naturalmedicines.therapeuticresearch.com/)
- [MSK About Herbs](https://www.mskcc.org/cancer-care/diagnosis-treatment/symptom-management/integrative-medicine/herbs)

### African Databases
- [ANPDB - Nucleic Acids Research 2024](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaf1186/8324951)
- [SANCDB - Journal of Cheminformatics](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00514-2)
- [p-ANAPL - PLoS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090655)

### Latin American Databases
- [NuBBEDB - Nature Scientific Reports](https://www.nature.com/articles/s41598-017-07451-x)
- [BIOFACQUIM - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6358837/)

### Global Databases
- [COCONUT 2.0 - Nucleic Acids Research](https://academic.oup.com/nar/article/53/D1/D634/7908792)
- [LOTUS - eLife](https://elifesciences.org/articles/70780)
- [NPASS 2023 - Nucleic Acids Research](https://academic.oup.com/nar/article/51/D1/D621/6851105)
- [SuperNatural 3.0 - Nucleic Acids Research](https://academic.oup.com/nar/article/51/D1/D654/6833249)
- [NAPRALERT - UIC](https://pharmacognosy.pharmacy.uic.edu/napralert/)

### Interaction Databases
- [DrugBank 6.0 - Nucleic Acids Research](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767804/)
- [STRING - Nucleic Acids Research](https://academic.oup.com/nar/article/49/D1/D483/5983621)
- [PharmMapper 2017 - Nucleic Acids Research](https://academic.oup.com/nar/article/45/W1/W356/3791213)

---

*Document generated for Gene Platform data source expansion*
