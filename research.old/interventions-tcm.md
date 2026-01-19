# Traditional Chinese Medicine (TCM) Databases - Comprehensive Research

**Research Date:** January 2024
**Purpose:** Data sources for genetics/health knowledge base platform
**Status:** Comprehensive review of all major TCM databases

---

## Executive Summary

Traditional Chinese Medicine has the most extensive collection of specialized databases among all intervention categories. This document catalogs **30+ databases** covering herbs, formulas, compounds, targets, diseases, symptoms, and gene expression data. The TCM database landscape has seen significant growth, with major updates in 2023-2024 including BATMAN-TCM 2.0, HERB 2.0, ETCM v2.0, TCMM, and CMAUP 2024.

### Key Findings

1. **Largest Databases by Coverage:**
   - TCM-Mesh: 383,840 compounds
   - TCMBank: 61,966 ingredients, 9,192 herbs
   - HERB: 49,258 ingredients with gene expression data
   - BATMAN-TCM 2.0: 2.3M predicted target interactions

2. **Best for Molecular Targets:**
   - BATMAN-TCM 2.0 (API available, bulk download)
   - TCMSID (ADME properties, CC BY 4.0)
   - HIT 2.0 (curated target data)

3. **Best for Gene Expression:**
   - HERB 2.0 (2,231 high-throughput experiments)
   - ITCM (pharmacotranscriptomic data)

4. **Best for Clinical Evidence:**
   - CMAUP 2024 (691 clinical trials)
   - HERB 2.0 (8,558 clinical trials, 8,032 meta-analyses)

5. **Most Accessible (API/Bulk Download):**
   - BATMAN-TCM 2.0 (REST API, tab-delimited downloads)
   - TCMBank (downloadable, CC BY 4.0)
   - TCMNP (R package, open source)

---

## Database Catalog

### 1. TCMBank

| Attribute | Details |
|-----------|---------|
| **URL** | https://tcmbank.cn/ |
| **Maintainer** | China Medical University, Taiwan (extended from TCM Database@Taiwan) |
| **Content** | Herbs, ingredients, targets, diseases, deep learning predictions |
| **Coverage** | 9,192 herbs, 61,966 ingredients, 15,179 targets, 32,529 diseases |
| **Access Method** | Web interface; Bulk download available |
| **Data Format** | Downloadable datasets |
| **Schema** | Herb-ingredient-target-disease mapping |
| **Licensing** | CC BY 4.0 - allows commercial use with attribution |
| **Full Content** | Yes - comprehensive mapping data |
| **Gene/Protein Targets** | Yes - 15,179 targets |
| **Key Features** | Deep learning for Chinese-Western medicine exclusion prediction; largest herb-ingredient-target-disease mapping |

**Notes:** Extends TCM Database@Taiwan. Non-commercial designation may be outdated - verify current terms.

---

### 3. TCMID 2.0 (Traditional Chinese Medicine Integrated Database)

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.megabionet.org/tcmid/ |
| **Maintainer** | Megabionet, China |
| **Content** | Prescriptions, herbs, ingredients, targets, diseases, drugs |
| **Coverage** | 49,000 prescriptions, 8,159 herbs, 25,210 ingredients, 17,521 targets, 3,791 diseases, 6,828 drugs |
| **Access Method** | Web interface |
| **Data Format** | Web tables; linked to DrugBank, OMIM, PubChem |
| **Schema** | Six interconnected parts: prescriptions, herbs, ingredients, diseases, drugs, targets |
| **Licensing** | Academic use; contact for commercial |
| **Full Content** | Yes - includes MS spectra for 170 herbs |
| **Gene/Protein Targets** | Yes - 17,521 targets |
| **Key Features** | 778 herbal mass spectrometry spectra; text-mining integrated data |

---

### 4. HERB 2.0 (High-throughput Experiment and Reference-guided Database)

| Attribute | Details |
|-----------|---------|
| **URL** | http://herb.ac.cn/ |
| **Maintainer** | Tsinghua University, China |
| **Content** | Herbs, ingredients, targets, diseases, gene expression, clinical trials, meta-analyses |
| **Coverage** | 7,263 herbs, 49,258 ingredients, 12,933 targets, 28,212 diseases, 2,231 high-throughput experiments, 8,558 clinical trials, 8,032 meta-analyses |
| **Access Method** | Web interface; probe2gene script downloadable |
| **Data Format** | Gene expression matrices; web tables |
| **Schema** | Herbs/Ingredients -> DEGs -> Pathways -> Diseases |
| **Licensing** | Academic use |
| **Full Content** | Yes - includes RNA-seq and microarray data (2001-2024) |
| **Gene/Protein Targets** | Yes - differentially expressed genes from experiments |
| **Key Features** | First TCM database with high-throughput gene expression data from GEO; clinical evidence integration |

**Notes:** Best source for transcriptomic evidence in TCM research.

---

### 5. YaTCM (Yet another Traditional Chinese Medicine)

| Attribute | Details |
|-----------|---------|
| **URL** | http://cadd.pharmacy.nankai.edu.cn/yatcm/home |
| **Maintainer** | Nankai University, China |
| **Content** | Prescriptions, herbs, ingredients, targets, pathways, diseases |
| **Coverage** | 47,696 compounds, 6,220 herbs, 18,697 targets, 390 pathways |
| **Access Method** | Web interface with analysis tools |
| **Data Format** | Web tables; network visualization |
| **Schema** | Prescription -> Herbs -> Ingredients -> Targets -> Pathways -> Diseases |
| **Licensing** | Free for academic use |
| **Full Content** | Yes - 50 ADME properties per compound |
| **Gene/Protein Targets** | Yes - multi-voting chemical similarity ensemble approach for target prediction |
| **Key Features** | Similarity search, substructure search, pathway analysis, network pharmacology |

---

### 6. SymMap 2.0

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.symmap.org/ and https://www.bioinfo.org/symmap |
| **Maintainer** | Peking University, China |
| **Content** | TCM symptoms, herbs, ingredients, targets, modern symptoms, diseases |
| **Coverage** | 1,717 TCM symptoms, 499 herbs, 961 modern symptoms, 5,235 diseases, 19,595 ingredients, 4,302 targets |
| **Access Method** | Web interface; bulk download with P-value filtering |
| **Data Format** | Tab-delimited files; full/filtered sets available |
| **Schema** | TCM Symptoms <-> Herbs <-> Ingredients <-> Targets <-> Modern Symptoms <-> Diseases |
| **Licensing** | Academic use |
| **Full Content** | Yes - heterogeneous network with 32,281 nodes, 403,318 edges |
| **Gene/Protein Targets** | Yes - 4,302 target genes |
| **Key Features** | Expert-curated symptom mapping (17 TCM experts); bridges TCM and modern medicine phenotypes |

**Notes:** Unique for symptom-based drug discovery approach.

---

### 7. ETCM v2.0 (Encyclopedia of Traditional Chinese Medicine)

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.tcmip.cn/ETCM2/front/#/ |
| **Maintainer** | Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences |
| **Content** | Formulas, Chinese patent drugs, herbs, ingredients, targets, diseases |
| **Coverage** | 48,442 TCM formulas, 9,872 Chinese patent drugs, 2,079 medicinal materials, 38,298 ingredients, 1,040 targets, 8,045 diseases |
| **Access Method** | Web interface; downloadable knowledge graphs |
| **Data Format** | Network visualization; downloadable graphs |
| **Schema** | Formulas/Patents -> Herbs -> Ingredients -> Targets -> Diseases |
| **Licensing** | Academic use |
| **Full Content** | Yes - includes binding activity data for targets |
| **Gene/Protein Targets** | Yes - 2D ligand similarity-based target identification |
| **Key Features** | Ancient Chinese medical book coverage; JavaScript network visualization tool |

---

### 8. BATMAN-TCM 2.0

| Attribute | Details |
|-----------|---------|
| **URL** | http://bionet.ncpsb.org.cn/batman-tcm/ |
| **Maintainer** | National Center for Protein Sciences Beijing |
| **Content** | Formulas, herbs, ingredients, known/predicted target interactions |
| **Coverage** | 54,832 formulas, 8,404 herbs, 39,171 ingredients, 17,068 known TTIs, 2,319,272 predicted TTIs |
| **Access Method** | Web interface; REST API; Bulk download |
| **Data Format** | Tab-delimited text files; JSON via API |
| **Schema** | Ingredients -> Targets (with confidence scores) |
| **Licensing** | CC BY-NC 4.0 - non-commercial; contact for commercial |
| **Full Content** | Yes - all known and predicted interactions downloadable |
| **Gene/Protein Targets** | Yes - largest predicted interaction dataset (ROC AUC = 0.9663) |
| **Key Features** | API access; bioinformatics analysis workflow (KEGG, GO, OMIM/TTD); disease signature query |

**Notes:** Best programmatic access among TCM databases. Excellent for computational analysis.

---

### 9. HIT 2.0 (Herbal Ingredients' Targets)

| Attribute | Details |
|-----------|---------|
| **URL** | http://hit.cnu.edu.cn/ |
| **Maintainer** | Capital Normal University, Beijing |
| **Content** | Herbs, ingredients, experimentally validated targets |
| **Coverage** | 1,250+ herbs, 1,237 ingredients, 2,208 targets, 10,031 compound-target activity pairs |
| **Access Method** | Web interface |
| **Data Format** | Web tables with quality indicators |
| **Schema** | Herbs -> Ingredients -> Targets (with activity type: direct/indirect, activated/inhibited) |
| **Licensing** | Academic use |
| **Full Content** | Yes - includes activity types and experimental evidence |
| **Gene/Protein Targets** | Yes - manually curated from 2000-2020 literature |
| **Key Features** | Quality indicators for each interaction; experimental validation focus |

---

### 10. TCM-ID (Traditional Chinese Medicine Information Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://bidd.group/TCMID/ |
| **Maintainer** | BIDD Group, National University of Singapore |
| **Content** | Prescriptions, herbs, ingredients, 3D structures |
| **Coverage** | 1,588 prescriptions, 1,313 herbs, 5,669 ingredients, 3,725 3D structures |
| **Access Method** | Web interface |
| **Data Format** | 3D structures available |
| **Schema** | Prescriptions -> Herbs -> Ingredients -> Targets |
| **Licensing** | Academic use |
| **Full Content** | Yes - includes INVDOCK computational targets |
| **Gene/Protein Targets** | Yes - computational and experimental |
| **Key Features** | One of earliest TCM databases (2006); 3D structure availability |

---

### 11. CMAUP 2024 (Collective Molecular Activities of Useful Plants)

| Attribute | Details |
|-----------|---------|
| **URL** | https://bidd.group/CMAUP/ |
| **Maintainer** | BIDD Group, National University of Singapore |
| **Content** | Plants, ingredients, targets, diseases, clinical trials, drug development, transcriptomics |
| **Coverage** | 7,865 plants, 60,222 ingredients, 758 targets, 1,399 diseases, 691 clinical trials, 14,516 ingredient clinical records |
| **Access Method** | Web interface |
| **Data Format** | Web tables; phylogenetic tree visualization |
| **Schema** | Plants -> Ingredients -> Targets -> Diseases; Clinical trial associations |
| **Licensing** | Academic use |
| **Full Content** | Yes - includes transcriptomic reversion analysis from 20,027 patient samples |
| **Gene/Protein Targets** | Yes - targets with disease associations |
| **Key Features** | Clinical trial integration; DNA barcodes for 3,949 plants; SwissADME/HobPre bioavailability predictions |

**Notes:** Excellent for clinical translation evidence.

---

### 12. TCMM (Traditional Chinese Medicine Modernization)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.tcmm.net.cn/ |
| **Maintainer** | Peng Cheng Laboratory, Shenzhen |
| **Content** | 20 types of TCM concepts with 46 biological relations |
| **Coverage** | 3,447,023 records; integrates TCMBank, CPMCP, SymMap, PharMeBINet, PrimeKG |
| **Access Method** | Web interface; contact authors for complete data |
| **Data Format** | Knowledge graph; web interface |
| **Schema** | Unified schema across prescriptions, symptoms, herbs, ingredients, targets, diseases |
| **Licensing** | CC BY (open access article) |
| **Full Content** | Partial via web; complete data on request |
| **Gene/Protein Targets** | Yes - integrated from multiple sources |
| **Key Features** | Largest unified TCM modernization database (2024); Rx Gen module for prescription generation; GPT-based symptom descriptions |

**Notes:** Newest comprehensive database (April 2024). Bridges TCM theory with modern medicine.

---

### 13. TCM Database@Taiwan

| Attribute | Details |
|-----------|---------|
| **URL** | http://tcm.cmu.edu.tw/ |
| **Maintainer** | China Medical University, Taiwan |
| **Content** | Herbs, compounds with 2D/3D structures |
| **Coverage** | 453 herbs, ~20,000+ compounds |
| **Access Method** | Web interface; structure downloads |
| **Data Format** | CDX (2D) and Tripos mol2 (3D) formats |
| **Schema** | Herbs -> Compounds (with structures) |
| **Licensing** | CC BY - allows commercial use |
| **Full Content** | Yes - virtual screening ready |
| **Gene/Protein Targets** | Limited |
| **Key Features** | First large-scale TCM database (2011); included in ZINC database; CADD-ready structures |

---

### 14. TCM-Mesh

| Attribute | Details |
|-----------|---------|
| **URL** | http://mesh.tcm.microbioinformatics.org |
| **Maintainer** | Microbioinformatics Lab |
| **Content** | Herbs, compounds, genes, diseases, side effects, toxicity |
| **Coverage** | 6,235 herbs, 383,840 compounds, 14,298 genes, 6,204 diseases, 163,221 side effects, 71 toxicity records |
| **Access Method** | Web interface |
| **Data Format** | Network visualization |
| **Schema** | Integrates TCMID, STRING, OMIM, GAD, TOXNET, SIDER |
| **Licensing** | Academic use |
| **Full Content** | Yes - largest compound collection |
| **Gene/Protein Targets** | Yes - 14,298 genes from STRING/OMIM |
| **Key Features** | Toxicity and side effects data unique among TCM databases |

---

### 15. TCMSID (Traditional Chinese Medicine Simplified Integrated Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://tcm.scbdd.com/ |
| **Maintainer** | Sichuan University, China |
| **Content** | Herbs, ingredients, targets, ADME/T properties |
| **Coverage** | 499 herbs, 20,015 ingredients, 3,270 targets |
| **Access Method** | Web interface; simplification tools |
| **Data Format** | Web tables |
| **Schema** | Herbs -> Key Ingredients -> Targets |
| **Licensing** | CC BY 4.0 - allows commercial use |
| **Full Content** | Yes - 14 ADME/T properties; structural reliability scores |
| **Gene/Protein Targets** | Yes - multi-tool target prediction |
| **Key Features** | Unique simplification function to extract key pharmacologically active ingredients |

---

### 16. DCABM-TCM (Constituents Absorbed into Blood and Metabolites)

| Attribute | Details |
|-----------|---------|
| **URL** | http://bionet.ncpsb.org.cn/dcabm-tcm/ |
| **Maintainer** | National Center for Protein Sciences Beijing |
| **Content** | Blood-absorbed constituents (prototypes and metabolites) |
| **Coverage** | 192 prescriptions, 194 herbs, 1,816 blood constituents with structures, 970 prototypes, 385 metabolites |
| **Access Method** | Web interface |
| **Data Format** | Web tables; linked to PubChem |
| **Schema** | Prescriptions/Herbs -> Blood Constituents -> Targets -> Pathways -> Diseases |
| **Licensing** | No derivatives for commercial |
| **Full Content** | Yes - includes detection conditions |
| **Gene/Protein Targets** | Yes - via ADMET associations |
| **Key Features** | First database of actually bioavailable TCM compounds; experimental blood detection data |

**Notes:** Critical for understanding which compounds actually reach systemic circulation.

---

### 17. LTM-TCM (Linking Traditional Chinese Medicine with Modern Medicine)

| Attribute | Details |
|-----------|---------|
| **URL** | http://cloud.tasly.com/#/tcm/home |
| **Maintainer** | Tasly Pharmaceutical Group |
| **Content** | Symptoms, prescriptions, herbs, ingredients, targets |
| **Coverage** | 1,928 symptoms, 48,126 prescriptions, 9,122 plants, 34,967 ingredients, 13,109 targets, 1,170,133 interactions |
| **Access Method** | Web interface (Chinese/English) |
| **Data Format** | Web tables; network visualization |
| **Schema** | Symptoms <-> Prescriptions <-> Herbs <-> Ingredients <-> Targets |
| **Licensing** | Academic use |
| **Full Content** | Yes - includes 41,025 clinical treatment records, 213 ancient books |
| **Gene/Protein Targets** | Yes - 13,109 targets |
| **Key Features** | BioNLP-based interaction correction from 30M articles; virtual screening pipelines |

---

### 18. CPMCP (Chinese Patent Medicine and Compound Prescription)

| Attribute | Details |
|-----------|---------|
| **URL** | http://cpmcp.top |
| **Maintainer** | Academic consortium |
| **Content** | Chinese patent medicines, ancient prescriptions, symptoms |
| **Coverage** | 2,279 prescriptions, 71,414 prescription-symptom associations |
| **Access Method** | Web interface |
| **Data Format** | MariaDB backend; React frontend |
| **Schema** | Prescriptions -> Herbs -> TCM Symptoms -> Modern Symptoms |
| **Licensing** | Academic use |
| **Full Content** | Yes - expert-mapped symptom vocabularies |
| **Gene/Protein Targets** | Via linked databases |
| **Key Features** | Drug combination principles analysis; prescription compatibility |

---

### 19. TCMIO (Traditional Chinese Medicine on Immuno-Oncology)

| Attribute | Details |
|-----------|---------|
| **URL** | http://tcmio.xielab.net |
| **Maintainer** | Xie Lab |
| **Content** | Immuno-oncology focused: prescriptions, herbs, ingredients, IO targets |
| **Coverage** | 1,493 prescriptions, 618 TCM medicines, 16,437 ingredients, 400 IO targets, 120,000+ small molecules |
| **Access Method** | Web interface with cheminformatics tools |
| **Data Format** | Web tables; network visualization |
| **Schema** | Herbs -> Ingredients -> IO Targets -> Cancer pathways |
| **Licensing** | Academic use |
| **Full Content** | Yes - immuno-oncology specific |
| **Gene/Protein Targets** | Yes - 400 immuno-oncology targets |
| **Key Features** | Bingo cheminformatics integration; pathway enrichment; cancer immune microenvironment focus |

---

### 20. TM-MC 2.0 (Traditional Medicine Molecular Chemistry)

| Attribute | Details |
|-----------|---------|
| **URL** | https://tm-mc.kr |
| **Maintainer** | Korea Institute of Oriental Medicine |
| **Content** | Medicinal materials from Korean/Chinese/Japanese pharmacopoeias |
| **Coverage** | 635 materials, 21,306 unique compounds, 13,992 targets, 27,997 diseases, 5,075 prescriptions |
| **Access Method** | Web interface; SDF/JPG downloads |
| **Data Format** | SDF (structures), JPG images |
| **Schema** | Materials -> Compounds (with PubMed evidence) |
| **Licensing** | Free academic use |
| **Full Content** | Yes - all compounds with literature evidence |
| **Gene/Protein Targets** | Yes - 13,992 targets |
| **Key Features** | Northeast Asian traditional medicine focus; ChemDraw-verified structures; pharmacopoeia cross-reference |

---

### 21. NPASS 2023/2026 (Natural Product Activity and Species Source)

| Attribute | Details |
|-----------|---------|
| **URL** | https://bidd.group/NPASS/ |
| **Maintainer** | BIDD Group, National University of Singapore |
| **Content** | Natural products with quantitative activity values |
| **Coverage** | 94,400 NPs, 7,700 targets, 31,600 species sources, 446,552 activity records |
| **Access Method** | Web interface; integrated with PubChem |
| **Data Format** | Quantitative data (IC50, Ki, EC50, GI50, MIC) |
| **Schema** | Species -> Natural Products -> Targets (with activity values) |
| **Licensing** | Academic use |
| **Full Content** | Yes - quantitative bioactivity data |
| **Gene/Protein Targets** | Yes - with experimental activity values |
| **Key Features** | Quantitative activity focus; Chemical Checker integration; ADMET predictions |

**Notes:** Not TCM-specific but includes TCM plants. Excellent for quantitative pharmacology.

---

### 22. NPACT (Natural Products Activity Targets - Anti-cancer)

| Attribute | Details |
|-----------|---------|
| **URL** | http://crdd.osdd.net/raghava/npact/ |
| **Maintainer** | CSIR-IMTECH, India |
| **Content** | Plant-derived anti-cancer compounds |
| **Coverage** | 1,574 compounds, 353 cancer cell lines, 5,214 compound-cell line interactions, 1,980 compound-target interactions |
| **Access Method** | Web interface |
| **Data Format** | Web tables; linked to SuperNatural, HIT, CTD |
| **Schema** | Compounds -> Cancer Cell Lines -> Targets |
| **Licensing** | Academic use |
| **Full Content** | Yes - inhibitory values (IC50, ED50, EC50, GI50) |
| **Gene/Protein Targets** | Yes - experimentally validated |
| **Key Features** | Cancer cell line activity data; anti-cancer focus |

---

### 23. ITCM (Integrating Traditional Chinese Medicine)

| Attribute | Details |
|-----------|---------|
| **URL** | http://itcm.biotcm.net |
| **Maintainer** | Shanghai University of TCM / Second Military Medical University |
| **Content** | Active ingredients with pharmacotranscriptomic profiles |
| **Coverage** | 496 representative active ingredients (5x larger than previous) |
| **Access Method** | Web interface; TCM-Query screening |
| **Data Format** | Gene expression profiles |
| **Schema** | Ingredients -> Transcriptomic Signatures -> Diseases |
| **Licensing** | Academic use |
| **Full Content** | Yes - unified high-throughput experiments |
| **Gene/Protein Targets** | Yes - transcriptomic signatures |
| **Key Features** | Pharmacotranscriptomic screening; six signature search methods; COVID-19 and prostate cancer applications |

---

### 24. SuperTCM

| Attribute | Details |
|-----------|---------|
| **URL** | http://tcm.charite.de/supertcm |
| **Maintainer** | Charite - Universitatsmedizin Berlin |
| **Content** | TCM drugs, botanical species, ingredients, targets, pathways, diseases |
| **Coverage** | 6,516 TCM drugs, 5,372 botanical species, 55,772 ingredients, 543 targets, 254 KEGG pathways, 8,634 diseases |
| **Access Method** | Web interface |
| **Data Format** | Web tables |
| **Schema** | TCM Drugs -> Species -> Ingredients -> Targets -> Pathways -> Diseases |
| **Licensing** | Academic use |
| **Full Content** | Yes - historical linguistic data included |
| **Gene/Protein Targets** | Yes - 543 targets in KEGG pathways |
| **Key Features** | Biocultural approach combining biological pathways with historical/linguistic data |

---

### 25. TCMPG 2.0 (Traditional Chinese Medicine Plant Genomes)

| Attribute | Details |
|-----------|---------|
| **URL** | http://cbcb.cdutcm.edu.cn/TCMPG2/ |
| **Maintainer** | Chengdu University of TCM |
| **Content** | Medicinal plant genomes, herbs, ingredients |
| **Coverage** | 274 plants, 324 genomes, 414 herbs, 13,868 ingredients |
| **Access Method** | Web interface with analysis tools |
| **Data Format** | Genomic data; analysis outputs |
| **Schema** | Plants -> Genomes -> Herbs -> Ingredients |
| **Licensing** | CC BY 4.0 |
| **Full Content** | Yes - complete genomes |
| **Gene/Protein Targets** | Genomic level |
| **Key Features** | Heatmap, Primer3, PlantiSMASH, CRISPRCasFinder tools |

---

### 26. IGTCM (Integrative Genome Database of TCM Plants)

| Attribute | Details |
|-----------|---------|
| **URL** | http://yeyn.group:96/ |
| **Maintainer** | Academic consortium |
| **Content** | TCM plant genomes, genes, proteins, components |
| **Coverage** | 83 annotated genomes, 3,610,350 genes, 3,534,314 proteins, 4,032,242 RNAs, 1,033 components |
| **Access Method** | Web interface; sequence similarity search |
| **Data Format** | Genomic sequences; eggNOG annotations |
| **Schema** | Genomes -> Genes/Proteins -> KEGG Pathways |
| **Licensing** | Academic use |
| **Full Content** | Yes - integrated from GenBank/RefSeq |
| **Gene/Protein Targets** | Yes - full genomic annotation |
| **Key Features** | Bridges genomics to pharmacology; KEGG pathway mapping |

---

### 27. TCMNP (TCM Network Pharmacology)

| Attribute | Details |
|-----------|---------|
| **URL** | https://tcmlab.com.cn/tcmnp |
| **Maintainer** | Academic consortium |
| **Content** | Integrated network pharmacology data with R package |
| **Coverage** | 571 herbs, 17,118 ingredients, 10,013 diseases, 15,956 targets |
| **Access Method** | Web interface; R package (GitHub: tcmlab/TCMNP) |
| **Data Format** | R data structures; visualizations |
| **Schema** | Herbs -> Ingredients -> Targets -> Diseases/Pathways |
| **Licensing** | Open source R package |
| **Full Content** | Yes - integrates DrugBank, ETCM, DisGeNET, OMIM |
| **Gene/Protein Targets** | Yes - 15,956 targets |
| **Key Features** | R package for programmatic analysis; visualization with ggplot2/circlize; PPI and TF analysis |

**Notes:** Excellent for computational researchers preferring R environment.

---

### 28. TCM-Suite

| Attribute | Details |
|-----------|---------|
| **URL** | Not publicly listed |
| **Maintainer** | Huazhong University of Science and Technology |
| **Content** | Biological ingredient identification (Holmes-Suite) + Network pharmacology (Watson-Suite) |
| **Coverage** | 1,251,548 marker gene sequences, 235,470 biological ingredients |
| **Access Method** | Web interface |
| **Data Format** | DNA marker sequences; network data |
| **Schema** | DNA Markers -> Species Identification -> Ingredients -> Targets |
| **Licensing** | Academic use |
| **Full Content** | Yes - six marker gene types |
| **Gene/Protein Targets** | Yes - via Watson-Suite |
| **Key Features** | DNA-based TCM authentication; integrated identification + pharmacology pipeline |

---

### 29. TCMNPAS (TCM Network Pharmacology Analysis System)

| Attribute | Details |
|-----------|---------|
| **URL** | Contact authors |
| **Maintainer** | Academic consortium |
| **Content** | Network formulaology + network pharmacology |
| **Coverage** | Comprehensive integration |
| **Access Method** | Web interface |
| **Data Format** | Analysis outputs |
| **Schema** | Formulas -> Herbs -> Ingredients -> Targets |
| **Licensing** | Academic use |
| **Full Content** | Yes |
| **Gene/Protein Targets** | Yes |
| **Key Features** | Integrates network formulaology concept |

---

### 30. SysNatMed (2025)

| Attribute | Details |
|-----------|---------|
| **URL** | Research framework (no standalone database) |
| **Maintainer** | Multiple institutions |
| **Content** | Systems genetics framework for natural medicine discovery |
| **Coverage** | 8,399 herbs analyzed for Alzheimer's disease case study |
| **Access Method** | Methodology/framework |
| **Data Format** | Analysis pipelines |
| **Schema** | Gene-Disease associations -> NP-Target relationships -> Herb prioritization |
| **Licensing** | Open research |
| **Full Content** | Yes - methodology |
| **Gene/Protein Targets** | Yes - SMR-based gene-disease associations |
| **Key Features** | Barabasi systems medicine framework; GSEA-based herb prioritization |

---

## Comparative Analysis

### Coverage Comparison (Top Databases)

| Database | Herbs | Ingredients | Targets | Diseases | Formulas |
|----------|-------|-------------|---------|----------|----------|
| TCM-Mesh | 6,235 | 383,840 | 14,298 | 6,204 | - |
| TCMBank | 9,192 | 61,966 | 15,179 | 32,529 | - |
| HERB 2.0 | 7,263 | 49,258 | 12,933 | 28,212 | 9 |
| TCMID 2.0 | 8,159 | 25,210 | 17,521 | 3,791 | 49,000 |
| YaTCM | 6,220 | 47,696 | 18,697 | - | - |
| ETCM v2.0 | 2,079 | 38,298 | 1,040 | 8,045 | 48,442 |
| BATMAN-TCM 2.0 | 8,404 | 39,171 | 2.3M pairs | - | 54,832 |
| TCMSID | 499 | 20,015 | 3,270 | - | - |

### Access Methods Comparison

| Database | Web UI | API | Bulk Download | Data Format |
|----------|--------|-----|---------------|-------------|
| BATMAN-TCM 2.0 | Yes | Yes (REST) | Yes | TSV, JSON |
| TCMBank | Yes | No | Yes | Various |
| HERB 2.0 | Yes | No | Partial | Gene matrices |
| SymMap 2.0 | Yes | No | Yes | TSV with filtering |
| TCMSID | Yes | No | Via web | Web tables |
| ETCM v2.0 | Yes | No | Graphs only | Network files |
| TCMNP | Yes | R package | Yes | R objects |
| TM-MC 2.0 | Yes | No | SDF/JPG | Structures |

### Licensing Summary

| License Type | Databases |
|--------------|-----------|
| **CC BY 4.0** (Commercial OK) | TCMBank, TCMPG 2.0, TCMSID, TCM Database@Taiwan |
| **CC BY-NC** (Non-commercial) | BATMAN-TCM 2.0, DCABM-TCM |
| **Academic Use** | HERB, YaTCM, SymMap, ETCM, TCMID, HIT, CMAUP, most others |

### Unique Features by Database

| Feature | Best Database(s) |
|---------|-----------------|
| **ADME Properties** | TCMSID (14 ADME/T), YaTCM (50 properties) |
| **Gene Expression** | HERB 2.0 (2,231 experiments), ITCM (pharmacotranscriptomics) |
| **Target Prediction** | BATMAN-TCM 2.0 (2.3M predictions, AUC=0.97) |
| **Clinical Trials** | CMAUP (691 trials), HERB 2.0 (8,558 trials) |
| **Symptom Mapping** | SymMap (1,717 TCM symptoms mapped) |
| **Blood Constituents** | DCABM-TCM (actual bioavailable compounds) |
| **Formulas** | ETCM v2.0 (48,442), BATMAN-TCM 2.0 (54,832), TCMID (49,000) |
| **Genomics** | TCMPG 2.0, IGTCM |
| **Toxicity** | TCM-Mesh (163,221 side effect records) |
| **API Access** | BATMAN-TCM 2.0, TCMNP (R package) |

---

## Recommendations for Integration

### Priority 1: Core Databases (Recommended for Initial Integration)

1. **BATMAN-TCM 2.0** - Best API access, largest target predictions
2. **TCMSID** - ADME properties, CC BY 4.0 license
3. **HERB 2.0** - Gene expression and clinical trial data
4. **TCMBank** - Comprehensive coverage, downloadable

### Priority 2: Specialized Databases

1. **SymMap 2.0** - Symptom mapping (unique)
2. **DCABM-TCM** - Bioavailable compounds (unique)
3. **CMAUP 2024** - Clinical trial integration
4. **TM-MC 2.0** - Northeast Asian medicines, good structures

### Priority 3: Supporting Databases

1. **ETCM v2.0** - Ancient formulas
2. **TCMID 2.0** - Large prescription collection
3. **TCMSID** - Simplified/key compounds
4. **NPASS** - Quantitative activity data

### Data Harmonization Challenges

1. **Herb naming**: Different databases use different naming conventions (Chinese, Latin, common names)
2. **Compound identification**: Varying levels of structural verification
3. **Target annotation**: Mix of UniProt, Gene Symbol, and proprietary IDs
4. **Evidence quality**: Mix of predicted vs. experimentally validated interactions

### Suggested Approach

1. Use BATMAN-TCM 2.0 API as primary data source for programmatic access
2. Supplement with TCMSID for ADME properties
3. Add HERB 2.0 for transcriptomic evidence
4. Map identifiers using PubChem, UniProt, and ChEMBL as intermediaries
5. Validate critical interactions against HIT 2.0 (experimentally validated)

---

## References

1. Wang Y, et al. (2024). A critical assessment of Traditional Chinese Medicine databases as a source for drug discovery. Front. Pharmacol. 15:1303693. https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2024.1303693/full

2. Kong X, et al. (2024). BATMAN-TCM 2.0: an enhanced integrative database. Nucleic Acids Res. 52(D1):D1110-D1117. https://academic.oup.com/nar/article/52/D1/D1110/7334089

3. Fang S, et al. (2025). HERB 2.0: an updated database integrating clinical and experimental evidence. Nucleic Acids Res. 53(D1):D1404. https://academic.oup.com/nar/article/53/D1/D1404/7903361

4. Ren Z, et al. (2024). TCMM: A unified database for traditional Chinese medicine modernization. Comput Struct Biotechnol J. 23:1619-1630. https://www.sciencedirect.com/science/article/pii/S2001037024001053

5. Zeng X, et al. (2024). CMAUP database update 2024. Nucleic Acids Res. 52(D1):D1508. https://academic.oup.com/nar/article/52/D1/D1508/7332076

---

*Last updated: January 2024*
