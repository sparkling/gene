---
id: domain-womens-pediatric
title: "Women's Health and Pediatric Genomics Databases"
type: health-domain
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [disease, health-domain, databases]
---

# Women's Health and Pediatric Genomics Databases

**Document ID:** 43-76-WOMENS-PEDIATRIC
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [_index.md](./_index.md)

---

## TL;DR

Comprehensive catalog of 30+ databases for women's health and pediatric genomics, covering reproductive genetics, pregnancy/prenatal testing, hormone-related conditions, pediatric pharmacogenomics, growth/development, and fertility genetics. Priority integration targets include ClinVar, gnomAD, PharmGKB/CPIC, HPO, and ClinGen for variant interpretation and clinical guidance. Total storage estimate: ~100-150GB for high-priority sources.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Variant interpretation priority | ClinVar + ClinGen | FDA-recognized, public domain, 650M+ variants | Jan 2026 |
| Population frequency source | gnomAD v4.1 | Diverse populations, open access, critical for prenatal interpretation | Jan 2026 |
| Pediatric PGx standard | PharmGKB + CPIC | Evidence-based guidelines with pediatric-specific content | Jan 2026 |
| Phenotype standardization | HPO | Universal adoption (UK 100K Genomes, NIH UDP) | Jan 2026 |
| Controlled access strategy | EGA/dbGaP pathway | IRB approval for DDD, DECIPHER patient data | Jan 2026 |
| Hormone database approach | Multiple sources | No single comprehensive resource; aggregate GWAS + expression data | Jan 2026 |

---

## 1. Reproductive Genetics Databases

### 1.1 ClinGen (Clinical Genome Resource)

| Field | Value |
|-------|-------|
| **URL** | https://clinicalgenome.org/ |
| **Content** | Gene-disease validity curations, variant pathogenicity, dosage sensitivity, actionability assessments |
| **Records** | 650+ million registered variants (including gnomAD, ExAC, dbSNP, ClinVar) |
| **License** | Public domain, freely available |
| **API** | REST APIs via ClinGen Allele Registry and Linked Data Hub (JSON-LD format) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~5-10 GB |

**API Endpoints:**
- Allele Registry: https://reg.clinicalgenome.org/
- Linked Data Hub (LDH): REST APIs for gene and variant linking

**Relevance:** Variant interpretation for carrier screening, prenatal diagnosis. FDA-recognized for human variant curation.

### 1.2 GTR (NIH Genetic Testing Registry)

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/gtr/ |
| **Content** | Test purpose, methodology, validity, lab credentials, gene-disease associations |
| **Records** | Thousands of registered tests |
| **License** | Public domain (attribution requested) |
| **API** | NCBI E-utilities (esearch, efetch) - XML format |
| **Download** | FTP: https://ftp.ncbi.nlm.nih.gov/pub/GTR/data/ |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |

**Key Data Files:**
- `test_condition_gene.txt` - Test-condition-gene relationships

**Relevance:** Comprehensive registry of reproductive and prenatal genetic tests.

### 1.3 GeneReviews

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/books/NBK1116/ |
| **Content** | Expert-authored disease descriptions covering diagnosis, management, genetic counseling |
| **Records** | ~800+ disease chapters |
| **License** | Public (NCBI) |
| **API** | NCBI Bookshelf API; UCSC Genome Browser track |
| **Download** | Via NCBI Bookshelf |
| **Update Frequency** | Ongoing revisions |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~200 MB |

**Relevance:** Authoritative clinical guidance for reproductive genetics counseling, prenatal testing availability information.

### 1.4 OMIM (Online Mendelian Inheritance in Man)

| Field | Value |
|-------|-------|
| **URL** | https://www.omim.org/ |
| **Content** | Gene-phenotype relationships, inheritance patterns, clinical features |
| **Records** | 15,000+ genes, all known Mendelian disorders |
| **License** | Free for non-commercial research; commercial license required from JHU |
| **API** | REST API (registration required, annual key renewal, GET only, 20 entries/request limit) |
| **Download** | https://www.omim.org/downloads/ (requires agreement) |
| **Update Frequency** | Weekly refresh required if downloaded |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |

**Relevance:** Core reference for genetic disorder information in reproductive counseling.

---

## 2. Pregnancy/Prenatal Databases

### 2.1 ClinVar

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | Variant pathogenicity classifications (Pathogenic, Likely pathogenic, VUS, Benign), gene-disease associations |
| **Records** | Millions of variant submissions |
| **License** | Public domain |
| **API** | NCBI E-utilities; Clinical Tables API: https://clinicaltables.nlm.nih.gov/api/variants/v4/search |
| **Download** | FTP: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/ (VCF, XML formats) |
| **Update Frequency** | Weekly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2-5 GB |

**Key Files:**
- `clinvar.vcf.gz` (GRCh38)
- `variant_summary.txt.gz`
- `ClinVarVariationRelease_00-latest.xml.gz`

**Tools:**
- Simple ClinVar: http://simple-clinvar.broadinstitute.org/ - Web interface for filtering

**Relevance:** Essential for prenatal variant interpretation.

### 2.2 gnomAD (Genome Aggregation Database)

| Field | Value |
|-------|-------|
| **URL** | https://gnomad.broadinstitute.org/ |
| **Content** | Allele frequencies across diverse populations (European, Latino, African, South Asian, East Asian, Ashkenazi Jewish) |
| **Records** | v4.1: 730,947 exomes + 76,215 genomes; 241M small variants; 335,470 structural variants |
| **License** | Open access |
| **API** | GraphQL API; gnomadR package for R |
| **Download** | AWS Open Data: https://registry.opendata.aws/broad-gnomad/ |
| **Update Frequency** | Major releases |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | 19-91 GB (SQLite conversion) |

**Relevance:** Critical for interpreting variant rarity in prenatal findings. Population-specific frequencies essential for carrier screening interpretation.

### 2.3 DECIPHER

| Field | Value |
|-------|-------|
| **URL** | https://www.deciphergenomics.org/ |
| **Content** | Copy number variants, sequence variants, phenotypes from rare disease patients |
| **Records** | ~40,000 openly shared patient records; 63,000+ with limited sharing |
| **License** | Academic/research use; registration required |
| **API** | Web interface; data sharing via EGA for DDD study |
| **Download** | Registered access; patient data via EGA (EGAS00001000775) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 3 |
| **Storage Estimate** | Varies by access level |

**Relevance:** Prenatal diagnosis of developmental disorders, CNV interpretation.

### 2.4 dbGaP (Database of Genotypes and Phenotypes)

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/gap/ |
| **Content** | Individual-level genotype data, phenotypes, sequencing data for prenatal/pediatric cohorts |
| **Records** | 1000+ studies |
| **License** | Controlled access (IRB approval required) |
| **API** | SRA Toolkit for download; dbgap2x R package |
| **Download** | Authorized access required via dbGaP system |
| **Update Frequency** | Continuous |
| **Priority** | Tier 3 |
| **Storage Estimate** | Study-dependent |

**Access Tools:**
- SRA Toolkit (v2.10.2+)
- dbgap2x R package: https://github.com/gversmee/dbgap2x

**Relevance:** Research cohorts for prenatal/pediatric genomics studies.

---

## 3. Hormone-Related Databases

### 3.1 ERGDB (Estrogen Responsive Genes Database)

| Field | Value |
|-------|-------|
| **URL** | https://pmc.ncbi.nlm.nih.gov/articles/PMC308817/ (historical) |
| **Content** | Up/down-regulated genes, experimental conditions, tissue types, estrogen response elements |
| **Records** | Hundreds of curated genes |
| **License** | Academic use |
| **API** | Not available (static resource) |
| **Download** | Via publication supplementary materials |
| **Update Frequency** | Static |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~10 MB |

**Relevance:** Understanding estrogen-related reproductive conditions.

### 3.2 ERTargetDB

| Field | Value |
|-------|-------|
| **URL** | https://jme.bioscientifica.com/view/journals/jme/35/2/0350225.xml |
| **Content** | 40 genes with verified EREs, 381 microarray-identified genes, 2948 computationally predicted genes |
| **Records** | ~3,400 genes total |
| **License** | Academic use |
| **API** | Not available |
| **Download** | Via publication |
| **Update Frequency** | Static |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~20 MB |

**Relevance:** Hormone receptor biology for reproductive health.

### 3.3 EstroGene2.0

| Field | Value |
|-------|-------|
| **URL** | https://www.nature.com/articles/s41523-024-00709-4 |
| **Content** | Transcriptomic and cistromic data for ER biology, endocrine therapy response and resistance |
| **Records** | Aggregated public datasets |
| **License** | Academic use |
| **API** | Web interface |
| **Download** | Via publication/database |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~50 MB |

**Relevance:** Hormone therapy pharmacogenomics.

### 3.4 PCOS Genetics Resources

| Field | Value |
|-------|-------|
| **URL** | Multiple GWAS consortia; EGA: phs000368 |
| **Content** | 19 GWAS loci, gene-hormone interactions |
| **Records** | Multiple cohorts; ~70% heritability estimated |
| **License** | Controlled access |
| **API** | Via GWAS Catalog and dbGaP |
| **Download** | EGA/dbGaP controlled access |
| **Update Frequency** | As studies publish |
| **Priority** | Tier 3 |
| **Storage Estimate** | Study-dependent |

**Key Genes:** THADA, FSHR, INS-VNTR, DENND1A, PA2G4, ERBB3

**Relevance:** Fertility and hormonal disorder genetics.

### 3.5 Endometriosis Genetics

| Field | Value |
|-------|-------|
| **URL** | GWAS Catalog; https://www.nature.com/articles/s41588-023-01323-z |
| **Content** | 42 genome-wide significant loci (49 signals); gene-disease associations |
| **Records** | 60,674 cases; 701,926 controls |
| **License** | Open access (summary stats) |
| **API** | Via GWAS Catalog |
| **Download** | Summary statistics available |
| **Update Frequency** | As studies publish |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~100 MB |

**Key Genes:** ESR1, CYP19A1, HSD17B1, WNT4, VEZT

**Relevance:** Understanding hormonal disorder genetics.

### 3.6 Endomet Database

| Field | Value |
|-------|-------|
| **URL** | https://www.nature.com/articles/s41597-020-00623-x |
| **Content** | Differentially expressed genes, microarray data for endometrium and endometriosis lesions |
| **Records** | Most extensive collection of lesion expression data |
| **License** | Academic use |
| **API** | Web interface |
| **Download** | Via publication/GEO |
| **Update Frequency** | Static |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~50 MB |

**Relevance:** Molecular basis of endometriosis.

---

## 4. Pediatric Pharmacogenomics

### 4.1 PharmGKB (Pharmacogenomics Knowledge Base)

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | Drug-gene associations, clinical annotations, pathways, VIP gene summaries |
| **Records** | 1700+ genes, 700+ drugs, 9000+ literature annotations, 153 pathways |
| **License** | Free for research; no redistribution without permission |
| **API** | Web services for bulk download |
| **Download** | Registration required; zipped spreadsheets |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |

**Pediatric-Relevant Content:**
- Age-dependent gene expression considerations
- Ontogeny-genetics interactions
- Pediatric dosing considerations

**Relevance:** Core pediatric pharmacogenomics reference.

### 4.2 CPIC (Clinical Pharmacogenetics Implementation Consortium)

| Field | Value |
|-------|-------|
| **URL** | https://cpicpgx.org/ |
| **Content** | Evidence-based gene/drug clinical practice guidelines, pediatric assessments |
| **Records** | 25+ guidelines covering 100+ drugs |
| **License** | Creative Commons public domain |
| **API** | JSON download via PharmGKB integration |
| **Download** | https://cpicpgx.org/guidelines/ (PDF, supplement files) |
| **Update Frequency** | As guidelines update |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 MB |

**Key Pediatric Guidelines:**
- TPMT/NUDT15 for thiopurines (childhood ALL)
- CYP2C9/VKORC1 for warfarin (pediatric dosing)
- HLA-B*57:01 for abacavir (HIV-infected children)
- CYP2D6 for SSRIs (citalopram, escitalopram, sertraline)

**Relevance:** Actionable pediatric prescribing guidance.

### 4.3 FDA Pharmacogenomic Biomarkers

| Field | Value |
|-------|-------|
| **URL** | https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling |
| **Content** | Biomarker-drug pairs with labeling recommendations |
| **Records** | 261 drugs; 136 (52%) approved for pediatric use |
| **License** | Public domain |
| **API** | Not available (PDF table) |
| **Download** | PDF from FDA website |
| **Update Frequency** | As labels update |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~10 MB |

**Relevance:** Regulatory-approved pharmacogenomics for children.

### 4.4 DPWG (Dutch Pharmacogenetics Working Group)

| Field | Value |
|-------|-------|
| **URL** | https://www.knmp.nl/ (via PharmGKB) |
| **Content** | Therapeutic recommendations, dose adjustments complementing CPIC |
| **Records** | 100+ gene-drug pairs |
| **License** | Open access |
| **API** | Via PharmGKB |
| **Download** | Via PharmGKB |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

**Relevance:** Additional pediatric pharmacogenomics guidance.

---

## 5. Growth and Development Databases

### 5.1 DDG2P / Gene2Phenotype (G2P)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gene2phenotype |
| **Content** | Curated gene-disease associations with inheritance patterns, mutation mechanisms, phenotypes |
| **Records** | 1000+ genes |
| **License** | Open access |
| **API** | REST API via EBI |
| **Download** | Static releases; PanelApp integration |
| **Update Frequency** | Ongoing curation |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- Autosomal dominant/recessive/X-linked annotations
- Loss of function vs. activating mutation mechanisms
- Integrated with DECIPHER

**Relevance:** Core resource for pediatric developmental disorder diagnosis.

### 5.2 HPO (Human Phenotype Ontology)

| Field | Value |
|-------|-------|
| **URL** | https://hpo.jax.org/ |
| **Content** | Hierarchical phenotype terms linked to 7000+ diseases |
| **Records** | Thousands of phenotype terms |
| **License** | Open source |
| **API** | https://clinicaltables.nlm.nih.gov/api/hpo/v3/search; FHIR ValueSet |
| **Download** | OBO Foundry, BioPortal, Monarch Initiative |
| **Update Frequency** | Regular releases |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 MB |

**Major Adopters:**
- UK 100,000 Genomes Project
- NIH Undiagnosed Diseases Program
- RD-Connect GPAP

**Relevance:** Essential for pediatric phenotype coding and standardization.

### 5.3 Orphanet / ORDO

| Field | Value |
|-------|-------|
| **URL** | https://www.orpha.net/ |
| **Data Portal** | https://www.orphadata.com/ |
| **Content** | 6,100+ rare diseases, 5,400+ genes, clinical/diagnostic information |
| **Records** | 6,100+ diseases |
| **License** | Free access |
| **API** | ORDO ontology via BioPortal |
| **Download** | Orphadata datasets in multiple formats |
| **Update Frequency** | Regular |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~200 MB |

**Integrations:**
- MeSH, SNOMED CT, UMLS, MedDRA
- OMIM, UniProtKB, HGNC, Ensembl

**Relevance:** Pediatric rare disease reference.

### 5.4 GA4K (Genomics Answers for Kids)

| Field | Value |
|-------|-------|
| **URL** | https://www.childrensmercy.org/childrens-mercy-research-institute/research-areas/genomic-medicine-center/ |
| **Content** | Pediatric genome data, HiFi sequencing |
| **Records** | Target: 30,000 children (100,000 genomes); 1000+ HiFi genomes produced |
| **License** | Restricted |
| **API** | Not publicly available |
| **Download** | Research collaboration required |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 3 |
| **Storage Estimate** | Collaboration-dependent |

**Relevance:** Pediatric-focused genomic resource.

### 5.5 PGDP (Pediatric Genomics Discovery Program) - Yale

| Field | Value |
|-------|-------|
| **URL** | https://medicine.yale.edu/trial/pediatric-genomics-discovery-program-pgdp/ |
| **Content** | Novel gene-disease associations for childhood diseases |
| **Records** | Active research program |
| **License** | Research collaboration |
| **API** | Not publicly available |
| **Download** | Research collaboration required |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 3 |
| **Storage Estimate** | Collaboration-dependent |

**Relevance:** Novel pediatric disease gene discovery.

---

## 6. Fertility Genetics

### 6.1 Preimplantation Genetic Testing Resources

| Field | Value |
|-------|-------|
| **Guidelines** | SART/ASRM 2024 recommendations |
| **Testing Labs** | Registered in GTR (https://www.ncbi.nlm.nih.gov/gtr/) |
| **Technology** | Next-Generation Sequencing (NGS) for all 23 chromosome pairs |
| **Content** | Aneuploidy detection (PGT-A), monogenic disorder testing (PGT-M), structural rearrangements (PGT-SR) |
| **Priority** | Tier 2 |

**Key Statistics:**
- Women <35: up to 37% aneuploidy rate
- Women >40: 75% aneuploidy rate
- Biopsy: 5-10 cells at blastocyst stage
- ~5% embryo loss from handling/biopsy/freezing

**Relevance:** IVF embryo selection.

### 6.2 NIPT Data Resources

| Field | Value |
|-------|-------|
| **Studies** | PEGASUS (Canada), TRIDENT-2 (Netherlands) |
| **Content** | Cell-free DNA screening data, validation studies |
| **Technology** | Ultra-low-pass sequencing (0.06-0.5x depth) |
| **License** | Controlled access |
| **API** | Research collaboration required |
| **Download** | Via dbGaP/EGA for specific studies |
| **Priority** | Tier 3 |

**Key Findings:**
- 89.99% maternal DNA, 10.01% fetal DNA
- TRIDENT-2: 149,318 pregnancies; 1/275 additional findings beyond common aneuploidies

**Relevance:** Prenatal screening validation.

### 6.3 ISPD (International Society for Prenatal Diagnosis)

| Field | Value |
|-------|-------|
| **URL** | https://www.ispdhome.org/ |
| **Content** | Clinical guidelines, best practices, educational materials |
| **Records** | Professional resources |
| **License** | Professional membership |
| **API** | Not applicable |
| **Download** | Member resources |
| **Update Frequency** | Ongoing |
| **Priority** | Reference only |

**Relevance:** Clinical practice standards for prenatal diagnosis.

---

## 7. Cross-Cutting Resources

### 7.1 Newborn Screening Resources

#### NBSTRN (Newborn Screening Translational Research Network)

| Field | Value |
|-------|-------|
| **URL** | https://nbstrn.org/ |
| **Content** | LPDR (Longitudinal Pediatric Data Resource), NBS-CR, NBS-VR, ELSI resources |
| **License** | Research registration required |
| **API** | REDCap integration |
| **Download** | Registered access at https://nbstrn.org/login |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB |

**Relevance:** Newborn screening research infrastructure.

#### BeginNGS

| Field | Value |
|-------|-------|
| **URL** | https://radygenomics.org/begin-ngs-newborn-sequencing/ |
| **Content** | 53,575 variants for 412 genetic diseases with 1,603 therapies |
| **Technology** | Query federation via TileDB |
| **License** | Research/clinical collaboration |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |

**Relevance:** Next-generation newborn screening.

#### GTRx (Genome to Treatment)

| Field | Value |
|-------|-------|
| **URL** | https://gtrx.rbsapp.net/ |
| **Code** | https://github.com/rao-madhavrao-rcigm/gtrx |
| **License** | Open source |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

**Relevance:** Clinical decision support for newborn screening.

### 7.2 Deciphering Developmental Disorders (DDD) Study

| Field | Value |
|-------|-------|
| **URL** | https://www.sanger.ac.uk/collaboration/deciphering-developmental-disorders-ddd/ |
| **Data Access** | EGA: EGAS00001000775 |
| **Content** | ~14,000 children with trio exome/genome data from UK/Ireland |
| **License** | Research application required |
| **API** | Via EGA |
| **Download** | Controlled access via EGA |
| **Update Frequency** | Study-dependent |
| **Priority** | Tier 3 |
| **Storage Estimate** | Multi-TB (full data) |

**Relevance:** Large pediatric developmental disorder cohort for research.

---

## Integration Recommendations

### Priority Tier Matrix

| Tier | Databases | Rationale |
|------|-----------|-----------|
| **Tier 1 (MVP)** | ClinVar, gnomAD, PharmGKB, CPIC, HPO, ClinGen, FDA PGx | Essential variant interpretation + pediatric prescribing |
| **Tier 2** | GTR, GeneReviews, OMIM, DDG2P, Orphanet, DPWG, NBSTRN, BeginNGS | Clinical guidance + rare disease coverage |
| **Tier 3** | DECIPHER, dbGaP, hormone databases, GA4K, PGDP, DDD | Specialized research + controlled access |

### API Integration Approach

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
- Orphanet (Orphadata)

Tier 3 (Controlled Access):
- dbGaP studies
- DDD/DECIPHER patient data
- EGA submissions
```

### Data Size Planning

| Database | Approximate Size | Storage Notes |
|----------|------------------|---------------|
| gnomAD v4.1 | 19-91 GB | SQLite conversion recommended |
| ClinVar | 2-5 GB | VCF + metadata |
| PharmGKB | ~500 MB | Structured spreadsheets |
| HPO | ~100 MB | Ontology files |
| OMIM | ~500 MB | Text + relationships |
| ClinGen | 5-10 GB | Variant registry |
| Orphanet | ~200 MB | Multiple formats |
| **Total Tier 1** | ~30-110 GB | Varies by gnomAD format |

### License Compliance Matrix

| Database | Commercial Use | Attribution | Redistribution |
|----------|----------------|-------------|----------------|
| ClinVar | Yes | Yes | Yes |
| gnomAD | Yes | Yes | Yes |
| ClinGen | Yes | Yes | Yes |
| OMIM | License required | Yes | No |
| PharmGKB | Research only | Yes | No |
| HPO | Yes | Yes | Yes |
| Orphanet | Yes | Yes | Check terms |
| CPIC | Yes | Yes | Yes (CC0) |

---

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `Pathogenic` | Classification indicating a variant causes disease | ClinVar pathogenic variants |
| `VUS` | Variant of Uncertain Significance - variant with insufficient evidence for classification | ~50% of ClinVar submissions |
| `CNV` | Copy Number Variation - gain or loss of DNA segments | DECIPHER CNV database |
| `Carrier screening` | Testing to identify carriers of recessive disease alleles | Reproductive genetics panels |
| `Heritability` | Proportion of phenotypic variation attributable to genetic factors | PCOS ~70% heritability |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `PGT-A` | Preimplantation Genetic Testing for Aneuploidy - IVF embryo chromosome screening | 37-75% aneuploidy detection |
| `PGT-M` | Preimplantation Genetic Testing for Monogenic disorders - single gene disease testing | Carrier screening follow-up |
| `NIPT` | Non-Invasive Prenatal Testing - cell-free DNA screening during pregnancy | TRIDENT-2 study |
| `cfDNA` | Cell-free DNA - fragmented DNA circulating in blood | ~10% fetal origin in pregnancy |
| `Ontogeny` | Development of an organism from embryo to adult | Pediatric PGx consideration |
| `Dosage sensitivity` | Genes where copy number changes cause disease | ClinGen curation |
| `ERE` | Estrogen Response Element - DNA sequence where estrogen receptors bind | ERGDB, ERTargetDB |
| `Endometriosis` | Condition where endometrial tissue grows outside uterus | 42 GWAS loci identified |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ClinGen | Clinical Genome Resource | 650M+ registered variants |
| ClinVar | Clinical Variation database | NCBI variant pathogenicity |
| gnomAD | Genome Aggregation Database | 730K exomes, 76K genomes |
| CPIC | Clinical Pharmacogenetics Implementation Consortium | 25+ guidelines |
| DPWG | Dutch Pharmacogenetics Working Group | Complements CPIC |
| HPO | Human Phenotype Ontology | Pediatric phenotype coding |
| GTR | NIH Genetic Testing Registry | Test registry |
| OMIM | Online Mendelian Inheritance in Man | Genetic disorder reference |
| DECIPHER | Database of Chromosomal Imbalance and Phenotype | CNV interpretation |
| EGA | European Genome-phenome Archive | Controlled data access |
| dbGaP | Database of Genotypes and Phenotypes | NIH controlled access |
| PCOS | Polycystic Ovary Syndrome | 19 GWAS loci |
| TPMT | Thiopurine S-methyltransferase | Pediatric PGx gene |
| HLA | Human Leukocyte Antigen | Immune/drug response genes |
| ACMG | American College of Medical Genetics | Variant classification guidelines |
| DDD | Deciphering Developmental Disorders | 14K children study |
| NBSTRN | Newborn Screening Translational Research Network | NBS research infrastructure |
| CC0 | Creative Commons Zero (public domain) | CPIC license |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [_index.md](./_index.md) | Parent navigation index |
| [primary.md](./../genetics/primary.md) | Shared variant databases (ClinVar, dbSNP) |
| [pharmaceuticals.md](./../compounds/pharmaceuticals.md) | Pharmacogenomics overlap (PharmGKB, CPIC) |
| [rare.md](./rare.md) | Rare disease overlap (Orphanet, DECIPHER) |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial migration from research.old/data-sources-womens-pediatric.md |
