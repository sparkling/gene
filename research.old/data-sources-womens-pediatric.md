# Women's Health and Pediatric Genomics Databases

Research compiled for Gene Platform integration. This document catalogs databases relevant to reproductive genetics, pregnancy/prenatal testing, hormone-related conditions, pediatric pharmacogenomics, growth/development, and fertility genetics.

---

## Table of Contents

1. [Reproductive Genetics Databases](#1-reproductive-genetics-databases)
2. [Pregnancy/Prenatal Databases](#2-pregnancyprenatal-databases)
3. [Hormone-Related Databases](#3-hormone-related-databases)
4. [Pediatric Pharmacogenomics](#4-pediatric-pharmacogenomics)
5. [Growth and Development Databases](#5-growth-and-development-databases)
6. [Fertility Genetics](#6-fertility-genetics)
7. [Cross-Cutting Resources](#7-cross-cutting-resources)
8. [Integration Recommendations](#8-integration-recommendations)

---

## 1. Reproductive Genetics Databases

### 1.1 ClinGen (Clinical Genome Resource)

| Attribute | Details |
|-----------|---------|
| **URL** | https://clinicalgenome.org/ |
| **Description** | NIH-funded resource for defining clinical relevance of genes and variants. FDA-recognized for human variant curation. |
| **Content** | Gene-disease validity curations, variant pathogenicity, dosage sensitivity, actionability assessments |
| **API** | REST APIs via ClinGen Allele Registry and Linked Data Hub. JSON-LD format. 650+ million registered variants. |
| **Download** | File downloads available at https://clinicalgenome.org/tools/ |
| **License** | Public domain, freely available |
| **Size** | 650+ million registered variants including gnomAD, ExAC, dbSNP, ClinVar |
| **Relevance** | Variant interpretation for carrier screening, prenatal diagnosis |

**API Endpoints:**
- Allele Registry: https://reg.clinicalgenome.org/
- Linked Data Hub (LDH): REST APIs for gene and variant linking

### 1.2 GTR (NIH Genetic Testing Registry)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/gtr/ |
| **Description** | Centralized registry of genetic tests including carrier screening, prenatal testing, pharmacogenomics |
| **Content** | Test purpose, methodology, validity, lab credentials, gene-disease associations |
| **API** | NCBI E-utilities (esearch, efetch). XML data format. |
| **Download** | FTP: https://ftp.ncbi.nlm.nih.gov/pub/GTR/data/ |
| **License** | Public domain (attribution requested) |
| **Size** | Thousands of registered tests |
| **Relevance** | Comprehensive registry of reproductive and prenatal genetic tests |

**Key Data Files:**
- `test_condition_gene.txt` - Test-condition-gene relationships

### 1.3 GeneReviews

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/books/NBK1116/ |
| **Description** | Expert-authored, peer-reviewed disease descriptions covering diagnosis, management, genetic counseling |
| **Content** | Standardized clinical summaries for heritable conditions, prenatal testing availability |
| **API** | NCBI Bookshelf API; UCSC Genome Browser track |
| **Download** | Via NCBI Bookshelf |
| **License** | Public (NCBI) |
| **Size** | ~800+ disease chapters |
| **Relevance** | Authoritative clinical guidance for reproductive genetics counseling |

### 1.4 OMIM (Online Mendelian Inheritance in Man)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.omim.org/ |
| **Description** | Comprehensive catalog of human genes and genetic disorders |
| **Content** | Gene-phenotype relationships, inheritance patterns, clinical features |
| **API** | REST API (requires registration and annual key renewal). GET requests only. 20 entries/request limit. |
| **Download** | https://www.omim.org/downloads/ (requires agreement) |
| **License** | Free for non-commercial research; commercial license required from JHU |
| **Size** | 15,000+ genes, all known Mendelian disorders |
| **Relevance** | Core reference for genetic disorder information |

**API Requirements:**
- Registration required
- Annual key renewal
- Weekly data refresh required if downloaded

---

## 2. Pregnancy/Prenatal Databases

### 2.1 ClinVar

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Description** | Public archive of variant-phenotype relationships with clinical significance classifications |
| **Content** | Variant pathogenicity (Pathogenic, Likely pathogenic, VUS, Benign), gene-disease associations |
| **API** | NCBI E-utilities; Clinical Tables API: https://clinicaltables.nlm.nih.gov/api/variants/v4/search |
| **Download** | FTP: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/ (VCF, XML formats) |
| **License** | Public domain |
| **Size** | Millions of variant submissions |
| **Relevance** | Essential for prenatal variant interpretation |

**Key Files:**
- `clinvar.vcf.gz` (GRCh38)
- `variant_summary.txt.gz`
- `ClinVarVariationRelease_00-latest.xml.gz`

**Simple ClinVar Tool:** http://simple-clinvar.broadinstitute.org/ - Web interface for filtering

### 2.2 gnomAD (Genome Aggregation Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://gnomad.broadinstitute.org/ |
| **Description** | Population variant frequencies from exome and genome sequencing |
| **Content** | Allele frequencies across diverse populations (European, Latino, African, South Asian, East Asian, Ashkenazi Jewish) |
| **API** | GraphQL API (syntax updated, check current docs); gnomadR package for R |
| **Download** | AWS Open Data: https://registry.opendata.aws/broad-gnomad/ |
| **License** | Open access |
| **Size** | v4.1: 730,947 exomes + 76,215 genomes; ~120GB/chromosome raw; SQLite version ~19-91GB |
| **Relevance** | Critical for interpreting variant rarity in prenatal findings |

**Data Statistics (v4.1):**
- 241 million small variants
- 335,470 structural variants
- Diverse ancestral populations

### 2.3 DECIPHER

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.deciphergenomics.org/ |
| **Description** | Database of chromosomal imbalance and phenotype for developmental disorders |
| **Content** | Copy number variants, sequence variants, phenotypes from ~40,000+ rare disease patients |
| **API** | Web interface; data sharing via EGA for DDD study |
| **Download** | Registered access; patient data via EGA (EGAS00001000775) |
| **License** | Academic/research use; registration required |
| **Size** | ~40,000 openly shared patient records; 63,000+ with limited sharing |
| **Relevance** | Prenatal diagnosis of developmental disorders, CNV interpretation |

### 2.4 dbGaP (Database of Genotypes and Phenotypes)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/gap/ |
| **Description** | Archive of genotype-phenotype studies including prenatal and pediatric cohorts |
| **Content** | Individual-level genotype data, phenotypes, sequencing data |
| **API** | SRA Toolkit for download; dbgap2x R package |
| **Download** | Authorized access required via dbGaP system |
| **License** | Controlled access (IRB approval required) |
| **Size** | 1000+ studies |
| **Relevance** | Research cohorts for prenatal/pediatric genomics studies |

**Access Tools:**
- SRA Toolkit (v2.10.2+)
- dbgap2x R package: https://github.com/gversmee/dbgap2x

---

## 3. Hormone-Related Databases

### 3.1 ERGDB (Estrogen Responsive Genes Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://pmc.ncbi.nlm.nih.gov/articles/PMC308817/ (historical) |
| **Description** | Database of genes whose expression is regulated by estrogen |
| **Content** | Up/down-regulated genes, experimental conditions, tissue types, estrogen response elements |
| **API** | Not available (static resource) |
| **Download** | Via publication supplementary materials |
| **License** | Academic use |
| **Size** | Hundreds of curated genes |
| **Relevance** | Understanding estrogen-related reproductive conditions |

### 3.2 ERTargetDB

| Attribute | Details |
|-----------|---------|
| **URL** | https://jme.bioscientifica.com/view/journals/jme/35/2/0350225.xml |
| **Description** | Estrogen receptor target gene database with regulatory mechanisms |
| **Content** | 40 genes with verified EREs, 381 microarray-identified genes, 2948 computationally predicted genes |
| **API** | Not available |
| **Download** | Via publication |
| **License** | Academic use |
| **Size** | ~3,400 genes total |
| **Relevance** | Hormone receptor biology for reproductive health |

### 3.3 EstroGene2.0

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.nature.com/articles/s41523-024-00709-4 |
| **Description** | Database for endocrine therapy response and resistance in breast cancer |
| **Content** | Transcriptomic and cistromic data for ER biology |
| **API** | Web interface |
| **Download** | Via publication/database |
| **License** | Academic use |
| **Size** | Aggregated public datasets |
| **Relevance** | Hormone therapy pharmacogenomics |

### 3.4 PCOS Genetics Resources

| Attribute | Details |
|-----------|---------|
| **URL** | Multiple GWAS consortia; EGA: phs000368 |
| **Description** | Genetic studies of Polycystic Ovary Syndrome |
| **Content** | 19 GWAS loci (THADA, FSHR, INS-VNTR, DENND1A), gene-hormone interactions |
| **API** | Via GWAS Catalog and dbGaP |
| **Download** | EGA/dbGaP controlled access |
| **License** | Controlled access |
| **Size** | Multiple cohorts; ~70% heritability estimated |
| **Relevance** | Fertility and hormonal disorder genetics |

**Key Genes:** THADA, FSHR, INS-VNTR, DENND1A, PA2G4, ERBB3

### 3.5 Endometriosis Genetics

| Attribute | Details |
|-----------|---------|
| **URL** | GWAS Catalog; https://www.nature.com/articles/s41588-023-01323-z |
| **Description** | Genetic basis of endometriosis and comorbidities |
| **Content** | 42 genome-wide significant loci (49 signals); gene-disease associations |
| **API** | Via GWAS Catalog |
| **Download** | Summary statistics available |
| **License** | Open access (summary stats) |
| **Size** | 60,674 cases; 701,926 controls |
| **Relevance** | Understanding hormonal disorder genetics |

**Key Genes:** ESR1, CYP19A1, HSD17B1, WNT4, VEZT

### 3.6 Endomet Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.nature.com/articles/s41597-020-00623-x |
| **Description** | Gene expression database for endometrium and endometriosis lesions |
| **Content** | Differentially expressed genes, microarray data |
| **API** | Web interface |
| **Download** | Via publication/GEO |
| **License** | Academic use |
| **Size** | Most extensive collection of lesion expression data |
| **Relevance** | Molecular basis of endometriosis |

---

## 4. Pediatric Pharmacogenomics

### 4.1 PharmGKB (Pharmacogenomics Knowledge Base)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pharmgkb.org/ |
| **Description** | Curated resource for genetic variation and drug response |
| **Content** | Drug-gene associations, clinical annotations, pathways, VIP gene summaries |
| **API** | Web services for bulk download (Support Protocol 4) |
| **Download** | Registration required; zipped spreadsheets |
| **License** | Free for research; no redistribution without permission |
| **Size** | 1700+ genes, 700+ drugs, 9000+ literature annotations, 153 pathways |
| **Relevance** | Core pediatric pharmacogenomics reference |

**Pediatric-Relevant Content:**
- Age-dependent gene expression considerations
- Ontogeny-genetics interactions
- Pediatric dosing considerations

### 4.2 CPIC (Clinical Pharmacogenetics Implementation Consortium)

| Attribute | Details |
|-----------|---------|
| **URL** | https://cpicpgx.org/ |
| **Description** | Evidence-based gene/drug clinical practice guidelines |
| **Content** | Prescribing recommendations by genotype, pediatric assessments |
| **API** | JSON download via PharmGKB integration |
| **Download** | https://cpicpgx.org/guidelines/ (PDF, supplement files) |
| **License** | Creative Commons public domain |
| **Size** | 25+ guidelines covering 100+ drugs |
| **Relevance** | Actionable pediatric prescribing guidance |

**Key Pediatric Guidelines:**
- TPMT/NUDT15 for thiopurines (childhood ALL)
- CYP2C9/VKORC1 for warfarin (pediatric dosing)
- HLA-B*57:01 for abacavir (HIV-infected children)
- CYP2D6 for SSRIs (citalopram, escitalopram, sertraline)

### 4.3 FDA Pharmacogenomic Biomarkers

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling |
| **Description** | FDA-recognized pharmacogenomic associations in drug labels |
| **Content** | Biomarker-drug pairs with labeling recommendations |
| **API** | Not available (PDF table) |
| **Download** | PDF from FDA website |
| **License** | Public domain |
| **Size** | 261 drugs; 136 (52%) approved for pediatric use |
| **Relevance** | Regulatory-approved PGx for children |

### 4.4 DPWG (Dutch Pharmacogenetics Working Group)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.knmp.nl/ (via PharmGKB) |
| **Description** | Dutch pharmacogenetics guidelines complementing CPIC |
| **Content** | Therapeutic recommendations, dose adjustments |
| **API** | Via PharmGKB |
| **Download** | Via PharmGKB |
| **License** | Open access |
| **Size** | 100+ gene-drug pairs |
| **Relevance** | Additional pediatric PGx guidance |

---

## 5. Growth and Development Databases

### 5.1 DDG2P / Gene2Phenotype (G2P)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/gene2phenotype |
| **Description** | Curated gene-disease associations for developmental disorders |
| **Content** | 1000+ genes with inheritance patterns, mutation mechanisms, phenotypes |
| **API** | REST API via EBI |
| **Download** | Static releases; PanelApp integration |
| **License** | Open access |
| **Size** | 1000+ genes |
| **Relevance** | Core resource for pediatric developmental disorder diagnosis |

**Key Features:**
- Autosomal dominant/recessive/X-linked annotations
- Loss of function vs. activating mutation mechanisms
- Integrated with DECIPHER

### 5.2 HPO (Human Phenotype Ontology)

| Attribute | Details |
|-----------|---------|
| **URL** | https://hpo.jax.org/ |
| **Description** | Standardized vocabulary of phenotypic abnormalities |
| **Content** | Hierarchical phenotype terms linked to 7000+ diseases |
| **API** | https://clinicaltables.nlm.nih.gov/api/hpo/v3/search; FHIR ValueSet |
| **Download** | OBO Foundry, BioPortal, Monarch Initiative |
| **License** | Open source |
| **Size** | Thousands of phenotype terms |
| **Relevance** | Essential for pediatric phenotype coding |

**Major Adopters:**
- UK 100,000 Genomes Project
- NIH Undiagnosed Diseases Program
- RD-Connect GPAP

### 5.3 Orphanet / ORDO

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.orpha.net/ |
| **Data Portal** | https://www.orphadata.com/ |
| **Description** | Portal for rare diseases and orphan drugs |
| **Content** | 6,100+ rare diseases, 5,400+ genes, clinical/diagnostic information |
| **API** | ORDO ontology via BioPortal |
| **Download** | Orphadata datasets in multiple formats |
| **License** | Free access |
| **Size** | 6,100+ diseases |
| **Relevance** | Pediatric rare disease reference |

**Integrations:**
- MeSH, SNOMED CT, UMLS, MedDRA
- OMIM, UniProtKB, HGNC, Ensembl

### 5.4 GA4K (Genomics Answers for Kids)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.childrensmercy.org/childrens-mercy-research-institute/research-areas/genomic-medicine-center/ |
| **Description** | Children's Mercy pediatric genomics initiative |
| **Content** | Pediatric genome data, HiFi sequencing |
| **API** | Not publicly available |
| **Download** | Research collaboration required |
| **License** | Restricted |
| **Size** | Target: 30,000 children (100,000 genomes); 1000+ HiFi genomes produced |
| **Relevance** | Pediatric-focused genomic resource |

### 5.5 PGDP (Pediatric Genomics Discovery Program) - Yale

| Attribute | Details |
|-----------|---------|
| **URL** | https://medicine.yale.edu/trial/pediatric-genomics-discovery-program-pgdp/ |
| **Description** | Gene discovery for childhood diseases |
| **Content** | Novel gene-disease associations |
| **API** | Not publicly available |
| **Download** | Research collaboration required |
| **License** | Research collaboration |
| **Size** | Active research program |
| **Relevance** | Novel pediatric disease gene discovery |

---

## 6. Fertility Genetics

### 6.1 Preimplantation Genetic Testing Resources

**PGT-A/M/SR Testing Information:**

| Attribute | Details |
|-----------|---------|
| **Guidelines** | SART/ASRM 2024 recommendations |
| **Testing Labs** | Registered in GTR (https://www.ncbi.nlm.nih.gov/gtr/) |
| **Technology** | Next-Generation Sequencing (NGS) for all 23 chromosome pairs |
| **Content** | Aneuploidy detection, monogenic disorder testing, structural rearrangements |
| **Relevance** | IVF embryo selection |

**Key Statistics:**
- Women <35: up to 37% aneuploidy rate
- Women >40: 75% aneuploidy rate
- Biopsy: 5-10 cells at blastocyst stage
- ~5% embryo loss from handling/biopsy/freezing

### 6.2 NIPT Data Resources

| Attribute | Details |
|-----------|---------|
| **Studies** | PEGASUS (Canada), TRIDENT-2 (Netherlands) |
| **Content** | Cell-free DNA screening data, validation studies |
| **Technology** | Ultra-low-pass sequencing (0.06-0.5x depth) |
| **API** | Research collaboration required |
| **Download** | Via dbGaP/EGA for specific studies |
| **License** | Controlled access |
| **Relevance** | Prenatal screening validation |

**Key Findings:**
- 89.99% maternal DNA, 10.01% fetal DNA
- TRIDENT-2: 149,318 pregnancies; 1/275 additional findings beyond common aneuploidies

### 6.3 ISPD (International Society for Prenatal Diagnosis)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ispdhome.org/ |
| **Description** | Professional society with guidelines and resources |
| **Content** | Clinical guidelines, best practices, educational materials |
| **API** | Not applicable |
| **Download** | Member resources |
| **License** | Professional membership |
| **Relevance** | Clinical practice standards |

---

## 7. Cross-Cutting Resources

### 7.1 Newborn Screening Resources

#### NBSTRN (Newborn Screening Translational Research Network)

| Attribute | Details |
|-----------|---------|
| **URL** | https://nbstrn.org/ |
| **Description** | Tools and resources for NBS research |
| **Content** | LPDR (Longitudinal Pediatric Data Resource), NBS-CR, NBS-VR, ELSI resources |
| **API** | REDCap integration |
| **Download** | Registered access at https://nbstrn.org/login |
| **License** | Research registration required |
| **Relevance** | Newborn screening research infrastructure |

#### BeginNGS

| Attribute | Details |
|-----------|---------|
| **URL** | https://radygenomics.org/begin-ngs-newborn-sequencing/ |
| **Description** | Genome-based newborn screening platform |
| **Content** | 53,575 variants for 412 genetic diseases with 1,603 therapies |
| **Technology** | Query federation via TileDB |
| **License** | Research/clinical collaboration |
| **Relevance** | Next-generation newborn screening |

#### GTRx (Genome to Treatment)

| Attribute | Details |
|-----------|---------|
| **URL** | https://gtrx.rbsapp.net/ |
| **Code** | https://github.com/rao-madhavrao-rcigm/gtrx |
| **Description** | Genome-to-treatment decision support |
| **License** | Open source |
| **Relevance** | Clinical decision support for NBS |

### 7.2 Deciphering Developmental Disorders (DDD) Study

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.sanger.ac.uk/collaboration/deciphering-developmental-disorders-ddd/ |
| **Data Access** | EGA: EGAS00001000775 |
| **Description** | UK/Ireland study of undiagnosed developmental disorders |
| **Content** | ~14,000 children with trio exome/genome data |
| **API** | Via EGA |
| **Download** | Controlled access via EGA |
| **License** | Research application required |
| **Relevance** | Large pediatric developmental disorder cohort |

---

## 8. Integration Recommendations

### 8.1 Priority Databases for Gene Platform

| Priority | Database | Rationale |
|----------|----------|-----------|
| **High** | ClinVar | Essential variant pathogenicity data |
| **High** | gnomAD | Population frequencies critical for interpretation |
| **High** | PharmGKB/CPIC | Pediatric drug dosing |
| **High** | HPO | Phenotype standardization |
| **High** | ClinGen | Variant curation standards |
| **Medium** | OMIM | Gene-disease relationships |
| **Medium** | GTR | Test availability |
| **Medium** | Orphanet | Rare disease information |
| **Medium** | DDG2P/G2P | Developmental disorder genes |
| **Lower** | DECIPHER | Specialized CNV interpretation |

### 8.2 API Integration Approach

```
Tier 1 (Direct API):
- ClinVar (E-utilities)
- ClinGen Allele Registry (REST)
- HPO (REST)
- gnomAD (GraphQL)

Tier 2 (Bulk Download + Local):
- gnomAD (SQLite conversion)
- OMIM (weekly refresh)
- PharmGKB (registered download)

Tier 3 (Controlled Access):
- dbGaP studies
- DDD/DECIPHER patient data
```

### 8.3 Data Size Planning

| Database | Approximate Size | Storage Notes |
|----------|------------------|---------------|
| gnomAD v4.1 | 19-91 GB (SQLite) | Consider cloud storage |
| ClinVar | ~2-5 GB | VCF + metadata |
| PharmGKB | ~500 MB | Structured spreadsheets |
| HPO | ~100 MB | Ontology files |
| OMIM | ~500 MB | Text + relationships |

### 8.4 License Compliance Matrix

| Database | Commercial Use | Attribution | Redistribution |
|----------|----------------|-------------|----------------|
| ClinVar | Yes | Yes | Yes |
| gnomAD | Yes | Yes | Yes |
| ClinGen | Yes | Yes | Yes |
| OMIM | License required | Yes | No |
| PharmGKB | Research only | Yes | No |
| HPO | Yes | Yes | Yes |
| Orphanet | Yes | Yes | Check terms |

---

## Sources

### Prenatal/Reproductive Genetics
- [SMFM Prenatal Genetic Screening Guidelines](https://www.smfm.org/news/a-brief-guide-to-smfms-updated-prenatal-genetic-screening-recommendations)
- [ISPD Home](https://www.ispdhome.org/)
- [NIPT Overview - MedlinePlus](https://medlineplus.gov/genetics/understanding/testing/nipt/)
- [Prenatal Genetic Screening - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK557702/)

### Databases
- [ClinGen](https://clinicalgenome.org/)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [gnomAD](https://gnomad.broadinstitute.org/)
- [PharmGKB](https://www.pharmgkb.org/)
- [CPIC](https://cpicpgx.org/)
- [OMIM](https://www.omim.org/)
- [GTR](https://www.ncbi.nlm.nih.gov/gtr/)
- [HPO](https://hpo.jax.org/)
- [Orphanet](https://www.orpha.net/)
- [DECIPHER](https://www.deciphergenomics.org/)
- [Gene2Phenotype](https://www.ebi.ac.uk/gene2phenotype)
- [dbGaP](https://www.ncbi.nlm.nih.gov/gap/)

### Research Publications
- [Pediatric Pharmacogenomics - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8259917/)
- [PCOS Genetics - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6935309/)
- [Endometriosis GWAS - Nature Genetics](https://www.nature.com/articles/s41588-023-01323-z)
- [NIPT Research Data - Cell Genomics](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00288-X)
- [DDD Study - Sanger Institute](https://www.sanger.ac.uk/collaboration/deciphering-developmental-disorders-ddd/)

---

*Document compiled: January 2026*
*For Gene Platform genomics database integration*
