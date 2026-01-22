# Genetics Primary Databases

**Document ID:** 43-11-GENETICS-PRIMARY
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../INDEX.md](./../INDEX.md)

---

## TL;DR

Primary genetics databases provide population allele frequencies, functional variant annotations, epigenetic marks, and structural variant catalogs essential for personalized health intelligence. dbNSFP v4.9 (35 prediction scores) and AlphaMissense (71M missense variants, CC BY 4.0) are recommended as Tier 1 priorities for functional annotation. gnomAD-SV v4.1 provides population-scale structural variant frequencies, while TOPMed BRAVO and ALFA R4 offer diverse population allele frequencies. ENCODE 4 and Roadmap Epigenomics provide regulatory element annotations critical for non-coding variant interpretation.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary functional annotation | dbNSFP v4.9 | 35 prediction scores in single download; includes AlphaMissense, SpliceAI, CADD, REVEL | Jan 2026 |
| Missense pathogenicity | AlphaMissense | 71M variants scored; CC BY 4.0 allows commercial use; state-of-art accuracy | Jan 2026 |
| Population frequencies (diverse) | TOPMed BRAVO | 868M+ variants from 150K genomes; ~60% non-European ancestry | Jan 2026 |
| Population frequencies (NCBI) | ALFA R4 | Public domain; 200K+ subjects; excellent ClinVar coverage (960K+ RSIDs) | Jan 2026 |
| Structural variants | gnomAD-SV v4.1 | 1.2M SVs from 527K individuals; GraphQL API; open access | Jan 2026 |
| Regulatory annotation | ENCODE 4 | 926K human cCREs; comprehensive ChIP-seq, DNase-seq, ATAC-seq | Jan 2026 |
| Splice prediction | SpliceAI | Deep learning splice predictions; pre-computed scores available | Jan 2026 |
| African diversity | H3Africa/AGVP | 3M+ novel variants; addresses critical gap in reference panels | Jan 2026 |
| Asian diversity | GenomeAsia 100K | 66M+ variants from 219 Asian population groups | Jan 2026 |
| Primary compound identifier | ChEBI | Cross-database compatibility with pathway databases | Jan 2026 |

---

## Database Catalog

### 1. Population Genetics Databases

#### 1.1 TOPMed BRAVO

| Field | Value |
|-------|-------|
| **URL** | https://bravo.sph.umich.edu/ |
| **Content** | Allele frequencies from deep whole-genome sequencing |
| **Records** | 868M+ variants from 150,000+ whole genomes |
| **Populations** | ~60% non-European ancestry |
| **Version** | Freeze 10 (current) |
| **License** | Public (summary data), controlled access (individual) |
| **API** | REST API via bravo_api (https://github.com/statgen/bravo_api) |
| **Download** | Web browser, BioData Catalyst for controlled data |
| **Data Formats** | VCF, JSON via API |
| **Update Frequency** | Periodic freezes |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~10 GB (summary data) |

**Key Features:**
- Best resource for rare variants in diverse populations
- Trans-Omics for Precision Medicine (TOPMed) program
- Deep sequencing (30x average coverage)
- Comprehensive variant calling pipeline

**Access Notes:**
- Summary statistics publicly available via web browser
- Individual-level data requires dbGaP application
- BioData Catalyst platform for cloud-based analysis

---

#### 1.2 All of Us Researcher Workbench

| Field | Value |
|-------|-------|
| **URL** | https://www.researchallofus.org/ |
| **Content** | Genetic variants from diverse US population cohort |
| **Records** | 1B+ genetic variants from 245K+ individuals |
| **Novel Variants** | 275M+ previously unreported variants |
| **Diversity** | ~50% non-European ancestry |
| **License** | Data Use Agreement required |
| **API** | Researcher Workbench APIs (controlled tier) |
| **Download** | Researcher Workbench (free registration, $300 credits) |
| **Data Formats** | VDS, Hail MT, VCF, BGEN, PLINK |
| **Update Frequency** | Continuous enrollment; periodic data releases |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 TB (full), summary data available publicly |

**Key Features:**
- NIH precision medicine initiative
- Emphasis on historically underrepresented populations
- Linked to EHR, survey, and physical measurement data
- Fall 2025: Verily Workbench integration

**Access Notes:**
- Free registration for researchers
- Tiered access: registered, controlled, authorized
- Cloud-based analysis environment

---

#### 1.3 ALFA (Allele Frequency Aggregator) - NCBI

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/ |
| **Content** | Aggregated allele frequencies from dbGaP studies |
| **Records** | Allele frequencies from 200K+ subjects |
| **Version** | Release 4 (June 2025) - doubled cohort size from R3 |
| **Populations** | 12 populations (European, African, Asian, Latin American, etc.) |
| **ClinVar Coverage** | 960K+ ClinVar RS IDs (74% increase from R3) |
| **License** | Public domain |
| **API** | E-utilities API, Variation Services API |
| **Download** | https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/ |
| **Data Formats** | VCF, JSON via API |
| **Update Frequency** | Periodic releases |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~5 GB |

**Key Features:**
- NCBI-integrated, easy dbSNP cross-referencing
- Public domain license ideal for commercial use
- Excellent ClinVar overlap for clinical variant interpretation
- 12 population breakdowns

**API Example:**
```
GET https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid}
```

---

#### 1.4 UK Biobank Allele Frequency Browser (AFB)

| Field | Value |
|-------|-------|
| **URL** | https://afb.ukbiobank.ac.uk/ |
| **Content** | SNP and indel frequencies from large-scale WGS |
| **Records** | Frequencies from 150,119 WGS individuals (AFB); 490,640 WGS (full dataset, Sept 2025) |
| **Browser** | Global Biobank Engine (https://biobankengine.stanford.edu/) |
| **License** | Public (AFB), controlled (full data) |
| **API** | Limited public access; full via Research Analysis Platform |
| **Download** | AFB publicly available; full data requires application |
| **Data Formats** | VCF, summary statistics |
| **Update Frequency** | Major releases |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~30 TB (full), AFB summary ~10 GB |

**Key Features:**
- Largest single-population WGS cohort
- Linked phenotype and health record data
- PheWAS associations available
- European ancestry focus (limitation for diversity)

---

#### 1.5 GenomeAsia 100K

| Field | Value |
|-------|-------|
| **URL** | https://browser.genomeasia100k.org/ |
| **Content** | Whole-genome sequences from Asian populations |
| **Records** | 1,739 WGS from 219 Asian population groups |
| **Countries** | 64 countries across Asia |
| **Variants** | 66M+ catalogued variants |
| **License** | Open for research (data access agreement) |
| **API** | Data access via application |
| **Download** | Application to dataaccess@genomeasia100k.org |
| **Data Formats** | VCF (individual-level) |
| **Update Frequency** | Periodic updates |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~2 TB (full VCFs) |

**Key Features:**
- Best resource for Asian population-specific variants
- Covers diverse Asian ethnic groups
- Addresses critical gap in reference databases
- South Asian, East Asian, Southeast Asian representation

---

#### 1.6 H3Africa / African Genome Variation Project

| Field | Value |
|-------|-------|
| **URL** | https://h3africa.org/ |
| **Content** | WGS and dense genotypes from African populations |
| **Records** | 426+ WGS individuals, 50 ethnolinguistic groups; AGVP: 1,481 dense genotypes + 320 WGS |
| **Novel Variants** | 3M+ previously undescribed variants |
| **Repository** | EGA (EGAS00001002976) |
| **License** | Controlled access (bona fide researchers) |
| **API** | H3ABioNet tools (https://github.com/h3abionet/h3agwas) |
| **Download** | Application to dbac@h3africa.org |
| **Data Formats** | VCF, H3Africa 2.3M SNP array |
| **Update Frequency** | Project-based releases |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~500 GB |

**Key Features:**
- Critical resource for African genetic diversity
- Addresses major gap in global reference panels
- 50 ethnolinguistic groups across sub-Saharan Africa
- H3Africa SNP array (2.3M SNPs) optimized for African populations

**Related Resources:**
- AGenDA (Assessing Genetic Diversity in Africa): WGS from 9 African countries
- H3ABioNet: Analysis tools and pipelines

---

### 2. Functional Annotation Databases

#### 2.1 dbNSFP v4.9

| Field | Value |
|-------|-------|
| **URL** | https://www.dbnsfp.org/home |
| **Content** | Comprehensive functional annotations for all possible coding variants |
| **Records** | 83M non-synonymous SNVs + 2.4M splice-site SNVs |
| **Scores Included** | 35 prediction scores: SIFT, PolyPhen2, CADD, REVEL, AlphaMissense, SpliceAI, EVE, PrimateAI, etc. |
| **Conservation** | PhyloP, phastCons, GERP++, bStatistic |
| **Version** | v4.9a (academic), v4.9c (commercial) |
| **License** | Academic (v4.9a free), Commercial (v4.9c - excludes some scores) |
| **API** | Bulk download (no REST API) |
| **Download** | Amazon S3, Box |
| **Data Formats** | TSV, VCF-ready |
| **Update Frequency** | Periodic (October 2025: GENCODE 49, Ensembl 115) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~35 GB (compressed) |

**Included Prediction Scores:**
| Score | Type | Range | Threshold |
|-------|------|-------|-----------|
| SIFT | Conservation | 0-1 | <0.05 deleterious |
| PolyPhen2 | Structure/Conservation | 0-1 | >0.85 damaging |
| CADD | Ensemble | 0-99 (PHRED) | >20 pathogenic |
| REVEL | Ensemble | 0-1 | >0.5 pathogenic |
| AlphaMissense | Deep learning | 0-1 | >0.564 pathogenic |
| SpliceAI | Deep learning | 0-1 | >0.5 splice-altering |
| EVE | Evolutionary | 0-1 | >0.5 pathogenic |
| PrimateAI | Deep learning | 0-1 | >0.8 pathogenic |

**Key Features:**
- Single download captures 35 scores
- Pre-computed for all possible coding variants
- VEP plugin compatible
- Regular updates with new scores

---

#### 2.2 AlphaMissense (Google DeepMind)

| Field | Value |
|-------|-------|
| **URL** | https://alphamissense.hegelab.org/ |
| **Content** | Deep learning pathogenicity scores for missense variants |
| **Records** | 71M missense variant predictions |
| **Classification** | 57% likely benign, 32% likely pathogenic, 11% uncertain |
| **License** | CC BY 4.0 (academic and commercial use) |
| **API** | No REST API; bulk download |
| **Download** | Google Cloud Public Dataset |
| **R Package** | AlphaMissenseR (Bioconductor) |
| **Data Formats** | TSV, VEP plugin compatible |
| **Update Frequency** | Model updates |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~10 GB (all predictions) |

**Key Features:**
- State-of-art missense pathogenicity prediction
- Trained on protein structure (AlphaFold) and evolutionary data
- CC BY 4.0 license enables commercial use
- Ensembl VEP plugin available
- GitHub: https://github.com/google-deepmind/alphamissense

**Score Thresholds:**
| Category | Score Range | Proportion |
|----------|-------------|------------|
| Likely pathogenic | >0.564 | 32% |
| Uncertain | 0.34-0.564 | 11% |
| Likely benign | <0.34 | 57% |

---

#### 2.3 CADD (Combined Annotation Dependent Depletion)

| Field | Value |
|-------|-------|
| **URL** | https://cadd.gs.washington.edu/ |
| **Content** | Deleteriousness scores integrating multiple annotations |
| **Records** | All possible SNVs and indels |
| **Version** | v1.7 (current) |
| **Scores** | Raw C-scores and PHRED-scaled scores |
| **License** | Academic use free, commercial license available |
| **API** | Web API, offline scoring available |
| **Download** | Pre-scored files (coding and non-coding) |
| **Data Formats** | TSV, VCF |
| **Update Frequency** | Major version releases |
| **Priority** | Tier 1 (via dbNSFP) |
| **Storage Estimate** | ~300 GB (genome-wide), ~15 GB (coding only) |

**Key Features:**
- Integrates 60+ genomic features
- PHRED-like scoring (>20 top 1% most deleterious)
- Includes non-coding variants
- VEP plugin, included in dbNSFP

**Score Interpretation:**
| PHRED Score | Interpretation | Percentile |
|-------------|----------------|------------|
| 10 | Top 10% most deleterious | 90th |
| 20 | Top 1% most deleterious | 99th |
| 30 | Top 0.1% most deleterious | 99.9th |

---

#### 2.4 SpliceAI

| Field | Value |
|-------|-------|
| **URL** | https://spliceailookup.broadinstitute.org/ |
| **Content** | Deep learning splice-altering predictions |
| **Records** | All SNVs and indels |
| **Scores** | DS_AG (acceptor gain), DS_AL (acceptor loss), DS_DG (donor gain), DS_DL (donor loss) |
| **License** | Free for academic, commercial license from Illumina |
| **API** | Web lookup |
| **Download** | Illumina BaseSpace (https://basespace.illumina.com/s/otSPW8hnhaZR) |
| **PyPI** | `pip install spliceai` |
| **Data Formats** | TSV, VCF |
| **Update Frequency** | Model updates |
| **Priority** | Tier 1 (via dbNSFP) |
| **Storage Estimate** | ~50 GB (pre-computed scores) |

**Key Features:**
- Deep neural network trained on splice junctions
- Four separate scores for different splice effects
- 2025 update: SpliceAI-splint addresses precomputed limitations
- Included in dbNSFP

**Score Thresholds:**
| Threshold | Precision | Recall |
|-----------|-----------|--------|
| 0.2 | Lower | High |
| 0.5 | Balanced | Balanced |
| 0.8 | High | Lower |

---

#### 2.5 REVEL

| Field | Value |
|-------|-------|
| **URL** | https://sites.google.com/site/revelgenomics/ |
| **Content** | Ensemble missense pathogenicity scores |
| **Components** | MutPred, FATHMM, VEST, PolyPhen, SIFT, PROVEAN, etc. |
| **Range** | 0-1 (higher = more pathogenic) |
| **License** | Academic use |
| **API** | No REST API |
| **Download** | Direct download available |
| **Data Formats** | TSV |
| **Update Frequency** | Periodic |
| **Priority** | Tier 1 (via dbNSFP) |
| **Storage Estimate** | ~5 GB |

**Key Features:**
- Combines 13 individual prediction tools
- Trained to distinguish pathogenic from rare neutral variants
- Higher discriminative power than individual tools
- Included in dbNSFP, VEP plugin available

---

#### 2.6 EVE (Evolutionary model of Variant Effect)

| Field | Value |
|-------|-------|
| **URL** | https://evemodel.org/ |
| **Content** | Pathogenicity predictions using deep evolutionary models |
| **Records** | 36M+ variants across 3,219 genes |
| **VUS Classified** | 256K+ variants of unknown significance |
| **Method** | Bayesian VAE trained on evolutionary data |
| **License** | Academic use |
| **API** | No REST API |
| **Download** | evemodel.org, Supplementary Data |
| **GitHub** | https://github.com/OATML-Markslab/EVE |
| **Data Formats** | TSV |
| **Update Frequency** | Model updates |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 GB |

**Key Features:**
- Unsupervised learning from evolutionary sequences
- No reliance on known pathogenic variants for training
- Good for genes lacking clinical data
- Included in dbNSFP

---

#### 2.7 PrimateAI-3D

| Field | Value |
|-------|-------|
| **URL** | https://primateai3d.basespace.illumina.com/download |
| **Content** | Pathogenicity scores using primate variation and protein structures |
| **Method** | Semi-supervised 3D-CNN trained on 4.5M benign variants from primates |
| **License** | Free academic, commercial license from Illumina |
| **API** | No REST API |
| **Download** | Illumina BaseSpace |
| **GitHub** | https://github.com/Illumina/PrimateAI |
| **Data Formats** | TSV |
| **Update Frequency** | Model updates |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 GB |

**Key Features:**
- Uses primate common variants as benign training set
- Incorporates AlphaFold protein structures
- Structural context for variant interpretation

**Score Thresholds:**
| Score | Interpretation |
|-------|----------------|
| >0.8 | Likely pathogenic |
| <0.6 | Likely benign |
| 0.6-0.8 | Uncertain |

---

#### 2.8 MaveDB

| Field | Value |
|-------|-------|
| **URL** | https://www.mavedb.org/ |
| **Content** | Experimental variant effect measurements from MAVE assays |
| **Records** | 7M+ variant effect measurements |
| **Datasets** | 1,884 datasets (Nov 2024) |
| **Data Types** | Deep mutational scanning (DMS), MPRA, saturation genome editing |
| **License** | Open access |
| **API** | REST API (JSON format) |
| **API Documentation** | https://github.com/VariantEffect/mavedb-api |
| **Download** | CSV per score set |
| **Data Formats** | CSV, JSON |
| **Update Frequency** | Continuous community submissions |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 GB |

**Key Features:**
- Experimental rather than computational predictions
- Gold standard for functional validation
- Growing rapidly with new assay technologies
- Publication: Genome Biology Jan 2025

**API Example:**
```
GET https://api.mavedb.org/api/v1/score-sets/
```

---

#### 2.9 RegulomeDB v2

| Field | Value |
|-------|-------|
| **URL** | http://regulomedb.org |
| **Content** | Regulatory potential scores for non-coding variants |
| **Data Sources** | ENCODE, Roadmap Epigenomics, GTEx |
| **Scoring** | 1 (strongest) to 6 (weakest) regulatory evidence |
| **Records** | 650M intervals (hg19), 1.5B intervals (GRCh38) |
| **Builds** | hg19 and GRCh38 |
| **License** | Open access |
| **API** | Web interface |
| **Download** | Bulk download available |
| **Data Formats** | BED, TSV |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~10 GB |

**Key Features:**
- Integrates multiple regulatory evidence types
- SURF algorithm scores included
- Essential for non-coding variant interpretation
- Links to eQTL data

**Score Categories:**
| Score | Evidence Level |
|-------|----------------|
| 1a-1f | Strong evidence (TF binding + motif + expression) |
| 2a-2c | Good evidence (TF binding + motif OR expression) |
| 3a-3b | Moderate evidence (TF binding OR motif) |
| 4-6 | Minimal evidence |

---

### 3. Epigenetics Databases

#### 3.1 ENCODE 4

| Field | Value |
|-------|-------|
| **URL** | https://www.encodeproject.org/ |
| **Content** | Encyclopedia of DNA elements - regulatory annotations |
| **Records** | 926,535 human cCREs, 339,815 mouse cCREs |
| **Data Types** | ChIP-seq (H3K4me1/3, H3K27ac/me3, H3K36me3, H3K9me3, CTCF), DNase-seq, ATAC-seq, RNA-seq |
| **Biosamples** | 234 human biosamples, 1,794 experiments |
| **License** | Open access |
| **API** | REST API |
| **API Documentation** | https://www.encodeproject.org/help/rest-api/ |
| **Download** | FTP, rsync (recommended) |
| **Browser** | SCREEN web interface |
| **Data Formats** | BED, bigWig, BAM |
| **Update Frequency** | Continuous data production |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 TB (full), ~50 GB (processed cCREs) |

**Key Features:**
- Gold standard for regulatory element annotation
- Candidate cis-regulatory elements (cCREs) catalog
- SCREEN browser for visualization
- Phase 4 focus on single-cell data

**cCRE Categories:**
| Type | Description | Count (human) |
|------|-------------|---------------|
| PLS | Promoter-like signature | ~35K |
| pELS | Proximal enhancer-like | ~150K |
| dELS | Distal enhancer-like | ~700K |
| CTCF | CTCF-bound | ~40K |

---

#### 3.2 Roadmap Epigenomics

| Field | Value |
|-------|-------|
| **URL** | https://egg2.wustl.edu/roadmap/web_portal/ |
| **Content** | Reference human epigenomes across tissues and cell types |
| **Records** | 111 reference human epigenomes |
| **Data Types** | Histone ChIP-seq (H3K4me1/3, H3K27me3, H3K36me3, H3K9me3, H3K27ac, H3K9ac), DNA methylation, RNA-seq |
| **Methylation** | 104 datasets (WGBS, RRBS, MeDIP-seq, MRE-seq) |
| **License** | Open access |
| **API** | FTP download |
| **Download** | http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics |
| **Release** | Human Epigenome Atlas Release 9 |
| **Data Formats** | BED, bigWig, WIG |
| **Update Frequency** | Completed project; static |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 TB |

**Key Features:**
- 111 primary human cell/tissue types
- Chromatin state annotations (ChromHMM)
- DNA methylation at single-base resolution
- Foundational resource for tissue-specific regulation

---

#### 3.3 IHEC Data Portal

| Field | Value |
|-------|-------|
| **URL** | https://epigenomesportal.ca/ihec/ |
| **Content** | International Human Epigenome Consortium data |
| **Records** | 7,500+ reference epigenomic datasets |
| **Tissues** | 600+ tissue types |
| **Contributors** | ENCODE, Roadmap, CEEHRC, Blueprint, DEEP, AMED-CREST, KNIH |
| **License** | Open access |
| **API** | REST API |
| **Download** | Bulk download available |
| **Browser** | UCSC mirror with preloaded tracks |
| **Data Formats** | BED, bigWig |
| **Update Frequency** | Continuous contributions |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~10 TB (full) |

**Key Features:**
- Aggregates international epigenomics efforts
- Harmonized data processing
- Quality-controlled reference epigenomes
- Cross-species comparisons

---

#### 3.4 MethBank 4.0

| Field | Value |
|-------|-------|
| **URL** | https://ngdc.cncb.ac.cn/methbank/ |
| **Content** | Single-base resolution DNA methylation data |
| **Records** | 1,449 methylomes across 23 species |
| **Human Data** | 53,680 age-specific DMCs, 1,716 age-specific DMRs |
| **License** | Open access |
| **API** | Web interface |
| **Download** | Web interface, bulk download |
| **Tools** | DMR Toolkit for analysis |
| **Data Formats** | BED, TSV |
| **Update Frequency** | Aug 2025: 30 new human projects; Cancer Module Jan 2025 |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~500 GB |

**Key Features:**
- Multi-species methylation resource
- Age-specific methylation changes
- Cancer methylation module (2025)
- Chinese-centric but internationally relevant

---

#### 3.5 4D Nucleome Data Portal

| Field | Value |
|-------|-------|
| **URL** | https://data.4dnucleome.org/ |
| **Content** | 3D genome organization and nuclear architecture |
| **Records** | 1,800+ experiment sets, 36,000+ files |
| **Data Types** | Hi-C, ChIA-PET, Capture Hi-C, PLAC-Seq |
| **High-Resolution** | 40+ datasets with >1B read pairs |
| **License** | Open access |
| **API** | REST API |
| **Download** | Bulk download, filtering interface |
| **Visualization** | HiGlass, 3D Genome Browser, WashU Browser |
| **Data Formats** | .hic, .cool, .mcool, pairs |
| **Update Frequency** | Continuous (consortium active) |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~5 TB |

**Key Features:**
- 3D chromatin organization data
- Topologically associated domains (TADs)
- Chromatin loops and interactions
- Publication: Nature 2025 (final consortium paper)

---

#### 3.6 FANTOM5

| Field | Value |
|-------|-------|
| **URL** | https://fantom.gsc.riken.jp/5/ |
| **Content** | Promoter and enhancer atlas from CAGE data |
| **Records** | 210K human promoters, 63K human enhancers |
| **Samples** | 1,800 human CAGE profiles, 1,000 mouse profiles |
| **Data Types** | CAGE (cap analysis gene expression) |
| **License** | CC BY 4.0 |
| **API** | BioMart interface |
| **Download** | https://dbarchive.biosciencedbc.jp/data/fantom5/ |
| **Tools** | UCSC DataHub, Enhancer Selector |
| **Data Formats** | BED, GTF, TSV |
| **Update Frequency** | NAR 2025 update for lncRNA interfaces |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~200 GB |

**Key Features:**
- Most comprehensive promoter atlas
- Tissue-specific enhancer activity
- Transcription start site (TSS) resolution
- lncRNA expression atlas

---

### 4. Structural Variant Databases

#### 4.1 gnomAD-SV v4.1

| Field | Value |
|-------|-------|
| **URL** | https://gnomad.broadinstitute.org/ |
| **Content** | Structural variants from population-scale sequencing |
| **Records** | 1.2M SVs from 464,297 exomes + 63,046 genomes |
| **SV Types** | Deletions, duplications, insertions, inversions, mCNVs, translocations, complex |
| **Per-Person** | ~11,844 SVs average |
| **Detection** | GATK-SV (genomes), GATK-gCNV (exomes) |
| **License** | Open access |
| **API** | GraphQL API |
| **Download** | BED files, Hail Tables (AWS) |
| **Data Formats** | VCF, BED, Hail |
| **Update Frequency** | Major version releases |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~5 GB (summary), ~500 GB (full) |

**Key Features:**
- Gold standard for population SV frequencies
- Ancestry-stratified frequencies
- Quality-filtered call set
- Essential for CNV/SV interpretation

**SV Type Distribution:**
| Type | Approximate Count |
|------|-------------------|
| Deletions | ~500K |
| Duplications | ~300K |
| Insertions | ~200K |
| Inversions | ~50K |
| Complex | ~100K |

---

#### 4.2 Database of Genomic Variants (DGV)

| Field | Value |
|-------|-------|
| **URL** | http://dgv.tcag.ca/ |
| **Content** | Structural variation in healthy control populations |
| **Records** | 2.5M+ SV entries from 55 studies |
| **SV Types** | CNVs >50 bp, inversions |
| **Population** | Control/healthy individuals |
| **License** | Open access |
| **API** | Web interface |
| **Download** | FTP, web interface |
| **Data Formats** | GFF, BED, VCF |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |

**Key Features:**
- Benchmark for benign CNV catalog
- Multiple platforms and studies integrated
- Essential for filtering common SVs
- Long-standing reference resource

---

#### 4.3 dbVar (NCBI)

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/dbvar/ |
| **Content** | Structural variants from submitted studies |
| **Records** | 6M+ structural variants from 185+ studies |
| **Sources** | 1000 Genomes, gnomAD, ClinVar, ClinGen |
| **Common SVs** | nstd186 (curated common SVs) |
| **License** | Public domain |
| **API** | E-utilities |
| **Download** | FTP |
| **Data Formats** | VCF, GVF, tab-delimited |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~10 GB |

**Key Features:**
- NCBI official SV archive
- Links to ClinVar for clinical SVs
- Integration with NCBI variant resources
- ClinGen dosage sensitivity annotations

---

#### 4.4 DECIPHER

| Field | Value |
|-------|-------|
| **URL** | https://www.deciphergenomics.org/ |
| **Content** | Clinical structural variants with phenotypes |
| **Records** | 51,894 patients, 51K+ variants, 172K+ phenotype terms |
| **SV Types** | CNVs, aneuploidy, UPD, inversions, insertions, STRs |
| **Phenotypes** | HPO terms with quantitative data |
| **Projects** | 250+ projects from 40 countries |
| **License** | Open (consented data), Data Display Agreement for third parties |
| **API** | Data export |
| **Download** | Bulk download (data access agreement) |
| **Data Formats** | TSV, VCF |
| **Update Frequency** | Continuous clinical submissions |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~1 GB |

**Key Features:**
- Clinical CNV-phenotype associations
- Sanger Institute curation
- Patient consent for data sharing
- Integration with Ensembl and UCSC

---

#### 4.5 SV4GD (Structural Variation for Genetic Diseases)

| Field | Value |
|-------|-------|
| **URL** | https://sv4gd.com/ |
| **Content** | Curated disease-associated structural variants |
| **Records** | 10,305 germline SV records |
| **Diseases** | 58 neoplastic, 232 non-neoplastic genetic diseases |
| **Disease-Related** | 2,695 disease-related SVs |
| **Common Types** | Deletions (1,462), duplications (743) = 81.8% |
| **Curation** | Manual from ~300 publications |
| **License** | Open access |
| **API** | Web interface |
| **Download** | Web download |
| **Data Formats** | TSV |
| **Update Frequency** | Periodic |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- Disease-focused SV catalog
- Literature-curated associations
- Publication: NAR Jan 2025
- Useful for rare disease diagnostics

---

#### 4.6 HGSVC (Human Genome Structural Variation Consortium)

| Field | Value |
|-------|-------|
| **URL** | https://www.hgsvc.org/ |
| **Content** | Complete haplotype-resolved structural variants |
| **Records** | Complete haplotype sequences from 65 diverse individuals |
| **HGSVC3** | Near T2T completeness, fully resolved centromeres |
| **Data Types** | PacBio, Strand-Seq, Bionano maps |
| **License** | Open access |
| **API** | Data portal |
| **Download** | https://www.internationalgenome.org/data-portal/data-collection/hgsvc2 |
| **GitHub** | https://github.com/hgsvc |
| **Data Formats** | VCF, FASTA |
| **Update Frequency** | Major releases |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~1 TB |

**Key Features:**
- Gold standard SV reference
- Haplotype-resolved variants
- Long-read sequencing based
- Publications: Nature Aug 2025, Nature Jul 2025

---

#### 4.7 ClinGen Dosage Sensitivity

| Field | Value |
|-------|-------|
| **URL** | https://dosage.clinicalgenome.org/ |
| **Content** | Gene-level haploinsufficiency and triplosensitivity scores |
| **Records** | Scores for thousands of genes |
| **Scoring** | Haploinsufficiency (HI) and Triplosensitivity (TS) on 0-3 scale |
| **License** | Open access |
| **API** | REST API |
| **Download** | https://search.clinicalgenome.org/kb/downloads |
| **Data Formats** | CSV, BED, TSV |
| **Update Frequency** | Daily |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~10 MB |

**Key Features:**
- Clinical interpretation of CNVs
- Expert-curated dosage sensitivity
- Essential for CNV pathogenicity assessment
- ACMG/ClinGen framework aligned

**Scoring System:**
| Score | Evidence Level |
|-------|----------------|
| 3 | Sufficient evidence |
| 2 | Emerging evidence |
| 1 | Little evidence |
| 0 | No evidence |
| 40 | Autosomal recessive |

---

## Comparative Analysis

### Coverage Comparison

| Database | Variants | Populations | Primary Use | License |
|----------|----------|-------------|-------------|---------|
| dbNSFP | 83M nsSNVs | N/A | Functional annotation | Academic |
| AlphaMissense | 71M missense | N/A | Pathogenicity | CC BY 4.0 |
| gnomAD-SV | 1.2M SVs | 8 ancestry groups | SV frequencies | Open |
| TOPMed BRAVO | 868M | ~60% non-European | Rare variants | Public/Controlled |
| ALFA | 200K+ subjects | 12 populations | Allele frequencies | Public domain |
| ENCODE 4 | 926K cCREs | N/A | Regulatory elements | Open |

### Licensing Summary

| License Type | Databases | Commercial Use |
|--------------|-----------|----------------|
| **CC BY 4.0** | AlphaMissense, FANTOM5, ENCODE | Yes (with attribution) |
| **Public Domain** | ALFA | Yes |
| **Open Access** | gnomAD, gnomAD-SV, ENCODE, MaveDB, RegulomeDB | Yes |
| **Academic** | dbNSFP (v4.9a), CADD, SpliceAI, EVE, REVEL | No |
| **Controlled** | TOPMed (individual), All of Us, H3Africa, GenomeAsia 100K | Research only |

### API Availability

| Database | REST API | GraphQL | Bulk DL | Rate Limit |
|----------|----------|---------|---------|------------|
| gnomAD | No | Yes | Yes | None stated |
| TOPMed BRAVO | Yes | No | Limited | Unknown |
| ALFA | Yes | No | Yes | 3-10/sec |
| ENCODE | Yes | No | Yes | Reasonable |
| MaveDB | Yes | No | Yes | None stated |
| dbNSFP | No | No | Yes | N/A |

---

## Integration Priority

### Tier 1: Critical (Weeks 1-2)

| Database | Priority Reason | Size | License |
|----------|-----------------|------|---------|
| **dbNSFP v4.9** | 35 scores in one download; comprehensive | ~35 GB | Academic |
| **AlphaMissense** | State-of-art missense predictions; CC BY 4.0 | ~10 GB | CC BY 4.0 |
| **gnomAD-SV v4.1** | Population SV frequencies; open access | ~5 GB | Open |
| **ALFA R4** | NCBI-integrated; public domain | ~5 GB | Public domain |

### Tier 2: High Value (Weeks 3-4)

| Database | Priority Reason | Size | License |
|----------|-----------------|------|---------|
| **TOPMed BRAVO** | Diverse populations; rare variants | ~10 GB | Public |
| **ENCODE 4** | Regulatory element annotations | ~50 GB | Open |
| **MaveDB** | Experimental variant effects | ~5 GB | Open |
| **ClinGen Dosage** | CNV interpretation | ~10 MB | Open |

### Tier 3: Enrichment (Weeks 5-8)

| Database | Priority Reason | Size | License |
|----------|-----------------|------|---------|
| **GenomeAsia 100K** | Asian-specific variants | ~2 TB | Controlled |
| **H3Africa/AGVP** | African diversity | ~500 GB | Controlled |
| **4D Nucleome** | 3D genome context | ~5 TB | Open |
| **FANTOM5** | Enhancer/promoter annotations | ~200 GB | CC BY 4.0 |

---

## Integration Tools & Annotation Pipelines

### Ensembl VEP (Variant Effect Predictor)

| Field | Value |
|-------|-------|
| **URL** | https://www.ensembl.org/vep |
| **Version** | Release 115 (Sept 2025) |
| **Plugins** | AlphaMissense, SpliceAI, CADD, MaveDB, UTRAnnotator, NMD, OpenTargets |
| **License** | Apache 2.0 |
| **API** | REST API (15 req/sec) |
| **Installation** | Conda, Docker, web interface |
| **Updates 2025** | gnomAD-SV allele frequencies, ClinVar somatic classifications, All of Us frequencies |

### ANNOVAR

| Field | Value |
|-------|-------|
| **URL** | https://annovar.openbioinformatics.org/ |
| **Content** | Gene-based, region-based, filter-based annotation |
| **Databases** | gnomAD, dbSNP, ClinVar, dbNSFP |
| **License** | Free for academic use |
| **Updates** | March 2024 version |

---

## Size Estimates for Full Integration

| Category | Compressed Size | Uncompressed | Notes |
|----------|-----------------|--------------|-------|
| Population genetics | ~100 GB | ~500 GB | Summary-level only |
| Functional annotation | ~60 GB | ~300 GB | dbNSFP largest |
| Epigenetics | ~20 GB | ~100 GB | ENCODE cCREs focus |
| Structural variants | ~15 GB | ~75 GB | gnomAD-SV primary |
| **Total** | **~195 GB** | **~975 GB** | Summary/processed data |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-00-INDEX](./../INDEX.md) | Parent index |
| [43-10-GENETICS-INDEX](./43-10-GENETICS-INDEX.md) | Category index |
| [43-21-CLINICAL-VARIANTS](./43-21-CLINICAL-VARIANTS.md) | Clinical variant interpretation (ClinVar, etc.) |
| [43-31-PHENOTYPES](./43-31-PHENOTYPES.md) | Phenotype ontologies (HPO, OMIM) |
| [45-DATA-MODEL](../45-DATA-MODEL.md) | Entity structure for variant tables |
| [44-ARCHITECTURE](../44-ARCHITECTURE.md) | Storage and API architecture decisions |

---

## References

1. Taliun D, et al. (2021). Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program. Nature. 590:290-299. https://doi.org/10.1038/s41586-021-03205-y

2. All of Us Research Program Investigators. (2023). The "All of Us" Research Program. N Engl J Med. 381:668-676. https://doi.org/10.1056/NEJMsr1809937

3. Phan L, et al. (2020). ALFA: Allele Frequency Aggregator. NCBI. https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/

4. Cheng J, et al. (2023). AlphaMissense predicts the pathogenicity of missense variants. Science. 381:eadg7492. https://doi.org/10.1126/science.adg7492

5. Liu X, et al. (2020). dbNSFP v4: a comprehensive database of transcript-specific functional predictions. Genome Med. 12:103. https://doi.org/10.1186/s13073-020-00803-9

6. Rentzsch P, et al. (2019). CADD: predicting the deleteriousness of variants. Nucleic Acids Res. 47(D1):D886-D894. https://doi.org/10.1093/nar/gky1016

7. Jaganathan K, et al. (2019). Predicting Splicing from Primary Sequence with Deep Learning. Cell. 176:535-548. https://doi.org/10.1016/j.cell.2018.12.015

8. Collins RL, et al. (2020). A structural variation reference for medical and population genetics. Nature. 581:444-451. https://doi.org/10.1038/s41586-020-2287-8

9. ENCODE Project Consortium. (2020). Expanded encyclopaedias of DNA elements. Nature. 583:699-710. https://doi.org/10.1038/s41586-020-2493-4

10. Kundaje A, et al. (2015). Integrative analysis of 111 reference human epigenomes. Nature. 518:317-330. https://doi.org/10.1038/nature14248

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial catalog: 30 genetics databases migrated from research.old/data-sources-genetics-expanded.md |
