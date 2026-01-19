# Environmental and Mitochondrial Databases for Gene Platform

Research compiled: January 2026

This document catalogs databases relevant to environmental exposure, toxicogenomics, mitochondrial variants, detoxification pathways, and heavy metal interactions for integration into the Gene Platform.

---

## Table of Contents

1. [Toxicogenomics Databases](#1-toxicogenomics-databases)
2. [Environmental Exposure Databases](#2-environmental-exposure-databases)
3. [Mitochondrial Variant Databases](#3-mitochondrial-variant-databases)
4. [Detoxification Pathway Databases](#4-detoxification-pathway-databases)
5. [Heavy Metal Interaction Databases](#5-heavy-metal-interaction-databases)
6. [Integration Summary](#6-integration-summary)

---

## 1. Toxicogenomics Databases

### 1.1 Comparative Toxicogenomics Database (CTD)

| Attribute | Details |
|-----------|---------|
| **URL** | https://ctdbase.org |
| **Description** | Premier resource harmonizing cross-species heterogeneous data for chemical exposures and biological repercussions. Manually curates and interrelates chemical, gene, phenotype, anatomy, disease, taxa, and exposure content from published literature. |
| **Content** | As of August 2024: 3.8 million direct interactions from 149,000+ scientific articles; 17,700+ chemicals; 55,400 genes; 6,700 phenotypes; 7,200 diseases; 214,000 exposures; 980 anatomical terms; 630+ species. Total: 50+ million toxicogenomic relationships. |
| **Key Features** | CTD Tetramers (four-unit information blocks connecting chemical-gene-phenotype-disease); inferred relationships; batch query tools; Venn analysis |
| **API Access** | Batch Query Tool available; Data downloads at http://ctdbase.org/downloads/; Web services available (may require CAPTCHA verification for automated access) |
| **Data Formats** | CSV, TSV, XML downloads |
| **License** | Free for academic use; commercial licensing may apply |
| **Database Size** | 50+ million relationships; monthly updates |
| **Integration Priority** | HIGH - Core resource for chemical-gene-disease relationships |

**Use Cases for Gene Platform:**
- Chemical exposure impact on gene expression
- Disease mechanisms from environmental exposures
- Drug-chemical interactions
- Air pollution/metal exposure linked to diseases (e.g., Alzheimer's)

### 1.2 T3DB - Toxin and Toxin Target Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.t3db.ca |
| **Description** | Comprehensive database combining detailed toxin data with toxin target information. Also known as the Toxic Exposome Database. |
| **Content** | 3,678 toxins; 41,602 synonyms; 2,073 toxin target records; 42,374 toxin-target associations. Each ToxCard contains 90+ data fields including chemical properties, toxicity values, molecular interactions, medical information. |
| **Categories** | Pollutants, pesticides, drugs, food toxins, household toxins, industrial/workplace toxins, cigarette toxins, uremic toxins, heavy metals and metal salts |
| **API Access** | Web interface with search capabilities; downloadable library files |
| **Data Formats** | Web interface; downloadable datasets |
| **License** | Non-proprietary; freely accessible; EULA required for downloads |
| **Database Size** | 3,678+ toxins with comprehensive annotations |
| **Integration Priority** | HIGH - Essential for toxin-target relationships |

**Use Cases for Gene Platform:**
- Toxin metabolism prediction
- Toxin/drug interaction prediction
- Environmental exposure studies
- Household/workplace toxin awareness

---

## 2. Environmental Exposure Databases

### 2.1 EPA CompTox Chemicals Dashboard

| Attribute | Details |
|-----------|---------|
| **URL** | https://comptox.epa.gov/dashboard/ |
| **Description** | Integrated chemistry, toxicity, and exposure information for environmental chemicals. Combines physicochemical properties, environmental fate/transport, exposure data, usage information, in vivo toxicity, and in vitro bioassay data. |
| **Content** | 1+ million chemicals; 300+ chemical lists by structure/category; physicochemical properties; exposure predictions; bioactivity data |
| **Key Features** | Chemical search; batch queries; exposure predictions; hazard assessments; integration with ToxCast bioassay data |
| **API Access** | **CTX APIs** (Computational Toxicology and Exposure APIs) - FREE; requires API key (contact ccte_api@epa.gov) |
| **API Endpoints** | Chemical API, Exposure API (via CPDat), Hazard API, Bioactivity API |
| **Data Formats** | JSON via REST API; CSV/Excel downloads |
| **License** | **Public Domain** - Free of all copyright restrictions; fully available for commercial and non-commercial use |
| **Database Size** | 1+ million chemicals |
| **Integration Priority** | HIGH - Government-backed, freely licensed exposure data |

**API Details:**
- Base URL: APIs hosted on cloud.gov
- Authentication: Individual API key required (free)
- Rate Limits: Contact EPA for specifics
- Documentation: https://www.epa.gov/comptox-tools/computational-toxicology-and-exposure-apis

**Use Cases for Gene Platform:**
- Environmental chemical exposure assessment
- Chemical property lookups
- Exposure prediction models
- Hazard characterization

### 2.2 HERCULES Exposome Research Center (Emory)

| Attribute | Details |
|-----------|---------|
| **URL** | https://emoryhercules.com |
| **Description** | NIEHS-funded center (since 2013) providing infrastructure for exposome research - the sum of all environmental exposures across the lifespan. |
| **Content** | Untargeted metabolomics covering 20,000+ chemical signals; environmental chemicals; food metabolites; microbiome-derived metabolites |
| **Key Features** | Exposome workbench (cloud-based database); mass spectrometry platforms; data sciences core |
| **API Access** | No public API documented; collaborative access through facility cores |
| **Data Formats** | Research collaboration-based |
| **License** | Academic research collaboration |
| **Database Size** | 20,000+ chemical signals per sample |
| **Integration Priority** | MEDIUM - Primarily for research collaborations |

### 2.3 NEXUS (NIH Exposome Coordinating Center)

| Attribute | Details |
|-----------|---------|
| **URL** | Via NIEHS (https://www.niehs.nih.gov) |
| **Description** | $7.7 million center (September 2024) at Columbia, Harvard, and USC to transform environmental health research. Building NIH Real World Data Platform. |
| **Content** | Environmental monitoring data; wearable sensor data; biospecimens; geospatial tools; social determinants data |
| **Key Features** | Data harmonization; AI/ML infrastructure; exposome-wide association studies (ExWAS) |
| **API Access** | Under development as part of NIH Real World Data Platform |
| **Data Formats** | Standardized formats in development |
| **License** | Expected to be publicly accessible |
| **Integration Priority** | MEDIUM - Emerging resource (monitor for availability) |

---

## 3. Mitochondrial Variant Databases

### 3.1 MITOMAP

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.mitomap.org |
| **Description** | Compendium of polymorphisms and mutations in human mitochondrial DNA. Catalogues published research on mtDNA variation. Operated by Center for Mitochondrial & Epigenomic Medicine, Children's Hospital of Philadelphia. |
| **Content** | 62,556 full-length GenBank sequences; 81,778 control region sequences; 19,892 single nucleotide variants (SNVs); disease mutation reports; mtDNA deletions, inversions, insertions, complex rearrangements; nuclear genes in mitochondrial disease |
| **Key Features** | Allele Search tool; MITOMASTER sequence analysis (variant identification, haplogroup prediction); haplogroup distribution data |
| **Haplogroup Distribution** | N: 68%; L: 10%; M: 20% (of 62,556 FL sequences) |
| **API Access** | No formal REST API; web-based tools; GenBank sequence downloads; **Alternative: MSeqDR mvTool API** (https://mseqdr.org/mvtool.php) |
| **Data Formats** | HTML tables, GenBank format, FASTA |
| **License** | Content is property of contributing authors; no explicit open license stated |
| **Database Size** | 62,556+ full-length sequences; updated every 4-6 months |
| **Integration Priority** | HIGH - Primary mtDNA variant resource |

**MSeqDR mvTool API (Alternative Access):**
- URL: https://mseqdr.org/mvtool.php
- Input formats: VCF, HGVS, classical mtDNA nomenclature
- Output formats: HTML, JSON, Excel, VCF
- Integrates: MITOMAP, HmtDB, GeneDx, 1000 Genomes, MSeqDR datasets
- FREE access

### 3.2 HmtDB (Human Mitochondrial Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.hmtdb.uniba.it (or via Database Commons) |
| **Description** | Genomic resource for mitochondrion-based human variability studies. Supports population genetics and biomedical research with healthy and disease phenotype sequences. |
| **Content** | Human mitochondrial genome sequences with phenotype annotations; site-specific nucleotide variability; continental subsets; haplogroup classifications (Phylotree-based) |
| **Key Features** | Query interface; MToolBox pipeline integration; haplogroup prediction; multi-aligned genome downloads |
| **API Access** | Web interface; downloadable alignments and variability data; MToolBox command-line tool |
| **Data Formats** | FASTA alignments; variability tables; relational DB (IBM DB2) |
| **License** | Academic research use |
| **Database Size** | Tripled since original release (1,255 genomes initially) |
| **Integration Priority** | MEDIUM - Complementary to MITOMAP |

**Use Cases for Gene Platform:**
- mtDNA variant pathogenicity assessment
- Population-specific haplogroup analysis
- mtDNA-disease associations
- Mitochondrial health scoring

---

## 4. Detoxification Pathway Databases

### 4.1 KEGG (Kyoto Encyclopedia of Genes and Genomes)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.kegg.jp |
| **API URL** | https://rest.kegg.jp |
| **Description** | Comprehensive pathway database including drug metabolism, xenobiotic biotransformation, and detoxification pathways. |
| **Content** | Pathway maps; compound database; enzyme database; reaction database; drug metabolism pathways (CYP450, Phase I/II) |
| **Key Pathways** | Drug metabolism - cytochrome P450; Drug metabolism - other enzymes; Glutathione metabolism; Xenobiotics biodegradation |
| **API Access** | **REST API** at https://rest.kegg.jp; Operations: info, list, find, get, conv, link |
| **Data Formats** | KGML (pathway markup); flat files; images |
| **License** | **Academic use only** for API; Commercial use requires license from Pathway Solutions (https://www.pathway.jp) |
| **Database Size** | Comprehensive (exact size varies by database) |
| **Integration Priority** | HIGH - Essential pathway resource (with license considerations) |

**API Examples:**
```
https://rest.kegg.jp/list/pathway/hsa  # List human pathways
https://rest.kegg.jp/get/hsa00980      # Drug metabolism - CYP450
https://rest.kegg.jp/find/compound/aspirin  # Find compound
```

**Libraries:**
- Python: Biopython Bio.KEGG.REST, kegg_pull
- R: KEGGREST (Bioconductor)

### 4.2 Reactome

| Attribute | Details |
|-----------|---------|
| **URL** | https://reactome.org |
| **Downloads** | https://reactome.org/download-data |
| **Description** | Open-source, manually curated, peer-reviewed pathway database. Strong coverage of biological oxidations and detoxification. |
| **Content** | Human pathways, reactions, processes; Phase I/II metabolism; Cytochrome P450 pathways; Xenobiotic metabolism; Detoxification of reactive oxygen species |
| **Key Pathways** | Biological oxidations (R-HSA-211859); Xenobiotics (R-HSA-211981); Metabolism (R-HSA-1430728); Detoxification of ROS |
| **API Access** | **Content Service API** (pathway data access); **Analysis Service API** (pathway analysis) |
| **Data Formats** | BioPAX level 3; SBML Level 3 v1; SBGN; PSI-MITAB; Neo4j GraphDB; MySQL dumps; SVG/PNG diagrams |
| **Mapping Files** | UniProt, ChEBI, ENSEMBL, miRBase, NCBI, GtoP identifiers |
| **License** | **Creative Commons** - Open source, open data |
| **Database Size** | Quarterly releases; Version 89+ (June 2024 onwards) on Zenodo |
| **Integration Priority** | HIGH - Open license, comprehensive pathways |

**API Documentation:** https://reactome.org/ContentService/
**GitHub:** https://github.com/reactome

### 4.3 DrugBank

| Attribute | Details |
|-----------|---------|
| **URL** | https://go.drugbank.com |
| **Description** | Gold standard for drug, drug-target, and pharmaceutical information. Extensive drug metabolism data including enzymes, transporters, and metabolic pathways. |
| **Content** | DrugBank 6.0 (2024): 4,563 FDA-approved drugs; 6,231 investigational drugs; 3,037 drug metabolites with structures; 1.4+ million drug-drug interactions; 2,475 drug-food interactions; 2,550 new drug-enzyme entries; 1,560 drug-transporter entries |
| **Metabolism Data** | Substrates, products, involved enzymes, reaction types, metabolite activity; manually curated from product monographs and literature |
| **API Access** | **Clinical API**; **Discovery API**; plugins for drug-drug interactions |
| **Data Formats** | XML, JSON, CSV downloads; protein identifiers linked to UniProt/PDB |
| **License** | **CC BY-NC 4.0** for datasets; **CC0 (Public Domain)** for Open Data; Academic: free license for certain uses; Commercial: paid license required |
| **Database Size** | 4,900+ distinct bioentities; 50,000+ drug entries |
| **Integration Priority** | HIGH - Drug metabolism reference (license considerations) |

**API Documentation:** https://docs.drugbank.com/

### 4.4 SuperCYP

| Attribute | Details |
|-----------|---------|
| **URL** | http://bioinformatics.charite.de/supercyp |
| **Prediction Tool** | http://insilico-cyp.charite.de/SuperCYPsPred/ |
| **Description** | Comprehensive database focused on Cytochrome P450 enzymes and drug metabolism interactions. |
| **Content** | 1,170 drugs; 3,800+ CYP interactions; 2,000+ SNPs and mutations with expression/activity effects; 3D homology models of all 48 human CYPs |
| **Key Features** | CYP-drug interaction analysis; drug cocktail compatibility checking; alternative combination finder; SNP effect viewer; CYP alignment tools |
| **CYP Coverage** | 1A2, 2C9, 2C19, 2D6, 3A4 (responsible for 90%+ of clinical drug metabolism) |
| **API Access** | No formal API; web interface; downloadable structure files |
| **Data Formats** | PDB (protein structures); MOL (drug files) |
| **License** | **CC BY-NC-SA 3.0** |
| **Integration Priority** | MEDIUM - Specialized CYP resource |

**SuperCYPsPred Features:**
- Input: 2D chemical structure
- Output: CYP inhibition profile, confidence scores, similar compounds, DDI tables, radar charts
- Free, no registration required

### 4.5 PharmGKB/ClinPGx

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pharmgkb.org (redirects to https://www.clinpgx.org) |
| **Description** | Pharmacogenomics knowledge resource with drug-gene relationships, clinical annotations, and dosing guidelines. |
| **Content** | Drug-gene associations; clinical annotations; dosing guidelines; pharmacokinetic pathways; variant annotations |
| **API Access** | Data access page available (requires JavaScript) |
| **License** | Academic use typically permitted; check current terms |
| **Integration Priority** | MEDIUM - Pharmacogenomics focus |

---

## 5. Heavy Metal Interaction Databases

### 5.1 CTD (Heavy Metals Focus)

The Comparative Toxicogenomics Database (Section 1.1) provides extensive heavy metal-gene interaction data.

| Heavy Metal | CTD Coverage |
|-------------|--------------|
| Lead (Pb) | Gene interactions, disease associations, exposure data |
| Cadmium (Cd) | Toxicogenomic relationships, phenotype effects |
| Mercury (Hg) | Methylmercury and inorganic forms |
| Arsenic (As) | Pathway mechanisms to diseases |
| Aluminum (Al) | Neurological associations |
| Chromium (Cr) | Carcinogenicity data |

**CTD Tetramer Example:** Air pollution/metal exposure -> intermediate genes/phenotypes -> Alzheimer disease

### 5.2 T3DB (Heavy Metals Focus)

The T3DB (Section 1.2) includes comprehensive heavy metal toxin data.

| Attribute | Details |
|-----------|---------|
| **Metal Coverage** | Heavy metals and metal salts as toxins |
| **Data Fields** | Chemical properties, toxicity values (LD50, etc.), molecular targets, medical information |
| **Target Information** | Protein targets affected by heavy metals |

### 5.3 EPA CompTox (Heavy Metals)

| Attribute | Details |
|-----------|---------|
| **Coverage** | Heavy metals in environmental exposure context |
| **Data Types** | Exposure estimates, hazard data, bioactivity |
| **API Access** | Via CTX APIs (free with API key) |

### 5.4 Genome-Wide Association Studies (GWAS) Data

| Attribute | Details |
|-----------|---------|
| **Key Findings** | Genetic variants associated with heavy metal blood levels |
| **Metals Studied** | Aluminum, cadmium, cobalt, copper, chromium, mercury, manganese, molybdenum, nickel, lead, zinc |
| **Notable Gene** | SLC39A8 (rs13107325) - manganese and zinc transport |
| **Data Sources** | GWAS Catalog (https://www.ebi.ac.uk/gwas/); individual study datasets |

---

## 6. Integration Summary

### Priority Matrix

| Database | Priority | License | API | Key Use Case |
|----------|----------|---------|-----|--------------|
| CTD | HIGH | Academic free | Batch/Download | Chemical-gene-disease |
| EPA CompTox | HIGH | **Public Domain** | REST API (free) | Environmental exposure |
| MITOMAP | HIGH | Author-owned | MSeqDR API | mtDNA variants |
| Reactome | HIGH | **CC Open** | REST API | Detox pathways |
| T3DB | HIGH | Free (EULA) | Web/Download | Toxin targets |
| DrugBank | HIGH | CC BY-NC 4.0 | REST API | Drug metabolism |
| KEGG | HIGH | Academic only | REST API | Pathways (license needed) |
| HmtDB | MEDIUM | Academic | Web/Download | mtDNA population |
| SuperCYP | MEDIUM | CC BY-NC-SA | Web/Download | CYP450 focus |
| HERCULES | MEDIUM | Collaborative | None | Exposome research |

### Recommended Integration Order

1. **Phase 1 - Core (Open/Free License)**
   - EPA CompTox (public domain, full API)
   - Reactome (CC, full API)
   - T3DB (free, comprehensive toxins)
   - MITOMAP + MSeqDR (mtDNA variants)

2. **Phase 2 - Academic License**
   - CTD (academic free, essential relationships)
   - KEGG (academic free, pathways)
   - HmtDB (mtDNA complement)

3. **Phase 3 - Commercial Considerations**
   - DrugBank (academic free / commercial paid)
   - SuperCYP (specialized CYP data)
   - PharmGKB (pharmacogenomics)

### Data Volume Estimates

| Database | Estimated Size | Update Frequency |
|----------|----------------|------------------|
| CTD | 50M+ relationships | Monthly |
| EPA CompTox | 1M+ chemicals | Quarterly |
| MITOMAP | 62K sequences, 20K SNVs | 4-6 months |
| Reactome | ~15,000 pathways | Quarterly |
| T3DB | 42K associations | Periodic |
| DrugBank | 50K+ entries | Periodic |
| KEGG | Varies by DB | Ongoing |

### API Integration Notes

**Free API Keys Required:**
- EPA CompTox: Contact ccte_api@epa.gov

**No API Key Needed:**
- Reactome Content/Analysis Services
- KEGG REST (academic use)
- MSeqDR mvTool

**Download-Based Access:**
- CTD (batch downloads)
- T3DB (library files)
- DrugBank (account required)

---

## References

1. CTD 2025 Update: https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkae883/7816860
2. DrugBank 6.0: https://academic.oup.com/nar/article/52/D1/D1265/7416367
3. MITOMAP: https://www.mitomap.org/MITOMAP
4. Reactome: https://reactome.org/
5. EPA CompTox APIs: https://www.epa.gov/comptox-tools/computational-toxicology-and-exposure-apis
6. T3DB: https://pmc.ncbi.nlm.nih.gov/articles/PMC2808899/
7. KEGG API: https://www.kegg.jp/kegg/rest/keggapi.html
8. SuperCYP: https://academic.oup.com/nar/article/38/suppl_1/D237/3112315
9. NIEHS Exposome: https://www.niehs.nih.gov/research/supported/exposure/bio
10. HmtDB 2016: https://academic.oup.com/nar/article/45/D1/D698/2605728
