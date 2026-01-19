# Mental Health, Nootropics, and Cognitive Databases for Gene Platform

Research compiled: 2026-01-19

This document catalogs databases relevant to psychiatric genetics, cognitive enhancement compounds, neurotransmitter pathways, brain gene expression, cognitive trait genetics, and sleep/circadian genetics.

---

## Table of Contents

1. [Psychiatric Genetics Databases](#1-psychiatric-genetics-databases)
2. [Nootropics Compound Databases](#2-nootropics-compound-databases)
3. [Neurotransmitter Pathway Databases](#3-neurotransmitter-pathway-databases)
4. [Brain Expression Databases](#4-brain-expression-databases)
5. [Cognitive Trait Genetics](#5-cognitive-trait-genetics)
6. [Sleep Genetics Databases](#6-sleep-genetics-databases)
7. [Summary Comparison Table](#7-summary-comparison-table)

---

## 1. Psychiatric Genetics Databases

### 1.1 Psychiatric Genomics Consortium (PGC)

The largest biological investigation in psychiatric genetics history, operating since 2007.

| Attribute | Details |
|-----------|---------|
| **URL** | https://pgc.unc.edu/ |
| **Download Page** | https://pgc.unc.edu/for-researchers/download-results/ |
| **Data Access Portal** | https://pgc.unc.edu/for-researchers/data-access-committee/data-access-portal/ |
| **GitHub** | https://github.com/psychiatric-genomics-consortium |
| **dbGaP Collection** | phs001254.v1.p1 |
| **Content** | GWAS meta-analysis data for major depression, bipolar disorder, ADHD, schizophrenia, Tourette syndrome, OCD, autism, anorexia nervosa, PTSD, substance use disorders |
| **API** | No dedicated API; data access through web portal and controlled-access requests |
| **License** | Academic research use; commercial use requires PGC DAC permission |
| **Size** | Summary statistics: varies by disorder (tens of MB to GB); Individual-level data: multiple TB across all cohorts |

**Data Types Available:**
- **Open Access**: GWAS summary statistics (downloadable without restriction)
- **Controlled Access**: Individual-level genotype data (requires approval, accessed via LISA/GCC cluster in Netherlands)

**Key Studies:**
- Schizophrenia: 150+ loci identified (N > 300,000)
- Major Depression: 100+ significant loci (N > 1 million)
- Bipolar Disorder: 64 genome-wide significant loci

**Access Requirements for Individual Data:**
1. LISA/GCC cluster account
2. Signed analyst memo
3. Signed WTCCC memo
4. Approved research proposal
5. Workgroup chair approval email

---

### 1.2 dbGaP Psychiatric Collections

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/gap/ |
| **PGC Collection** | https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/collection.cgi?study_id=phs001254.v1.p1 |
| **Content** | Individual-level genotype and phenotype data for psychiatric disorders |
| **API** | dbGaP API available for metadata; data access requires DAR approval |
| **License** | NIH Data Use Certification required |
| **Size** | Varies by study; total psychiatric collections span hundreds of GB |

---

### 1.3 NIMH Genetics Repository (NRGR)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.nimhgenetics.org/ |
| **PGC Bipolar Data** | https://www.nimhgenetics.org/resources/genomic-data/pgc-bp |
| **Content** | DNA samples, cell lines, and genomic data for psychiatric research |
| **API** | No public API; web-based ordering and data request system |
| **License** | Academic research; requires NIMH approval |
| **Size** | 200,000+ samples across multiple psychiatric phenotypes |

---

## 2. Nootropics Compound Databases

### 2.1 ChEMBL

The primary resource for bioactive compound data including nootropics and cognitive enhancers.

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/chembl/ |
| **API Documentation** | https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services |
| **API Live Docs** | https://www.ebi.ac.uk/chembl/api/data/docs |
| **Content** | 2.4M+ compounds, 1.5M+ assays, 20M+ activities, drug-target interactions |
| **API** | REST API with JSON/XML/CSV output; supports filtering, substructure search, similarity search |
| **License** | Creative Commons Attribution-Share Alike 3.0 Unported |
| **Size** | ~15 GB (MySQL dump); ~3 GB (SDF file) |

**Nootropics-Relevant Queries:**
```
# Approved cognitive drugs
https://www.ebi.ac.uk/chembl/api/data/molecule?max_phase=4&molecule_properties__mw_freebase__lte=500

# Compounds targeting acetylcholinesterase
https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id=CHEMBL220
```

**Key Targets for Nootropics:**
- Acetylcholinesterase (CHEMBL220)
- NMDA receptors (CHEMBL1907602)
- Dopamine receptors (D1-D5)
- Serotonin receptors (5-HT subtypes)
- Alpha-7 nicotinic receptor (CHEMBL3193)

---

### 2.2 DrugBank

| Attribute | Details |
|-----------|---------|
| **URL** | https://go.drugbank.com/ |
| **API Documentation** | https://docs.drugbank.com/ |
| **Content** | 15,000+ drug entries with detailed pharmacology, targets, interactions |
| **API** | Clinical API (commercial); Academic license available |
| **License** | Free for academic non-commercial use; commercial license required for API |
| **Size** | ~2 GB (full XML download) |

**Nootropic-Related Drugs in DrugBank:**
- Piracetam and racetam derivatives
- Modafinil/Armodafinil
- Methylphenidate
- Donepezil, Rivastigmine (AChE inhibitors)
- Memantine (NMDA antagonist)

---

### 2.3 PubChem

| Attribute | Details |
|-----------|---------|
| **URL** | https://pubchem.ncbi.nlm.nih.gov/ |
| **Downloads** | https://pubchem.ncbi.nlm.nih.gov/docs/downloads |
| **Content** | 118M+ compounds, 313M+ substances, biological assay data |
| **API** | PUG REST API, PUG-View, E-utilities |
| **License** | Public domain (US government work) |
| **Size** | Full database: ~100 GB+ (compound, substance, bioassay) |

**Nootropics Search Strategy:**
- Search by MeSH term "Nootropic Agents"
- Search by pharmacological class
- Cross-reference with cognitive/memory assays

---

### 2.4 ChEBI (Chemical Entities of Biological Interest)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/chebi/ |
| **Nootropic Ontology** | https://www.ebi.ac.uk/chebi/searchId.do?chebiId=66980 |
| **Content** | Chemical ontology with hierarchical classification |
| **API** | Web services available |
| **License** | Creative Commons CC0 1.0 |
| **Size** | ~500 MB |

**Nootropic Classification (CHEBI:66980):**
- Definition: Compound improving mental functions (cognition, memory, intelligence, motivation, attention, concentration)
- Includes natural and synthetic cognitive enhancers

---

### 2.5 Natural Products Databases (Nootropic Herbs)

#### NAPRALERT
| Attribute | Details |
|-----------|---------|
| **URL** | https://napralert.org/ |
| **Content** | 200,000+ compounds from natural sources including cognitive-enhancing plants |
| **License** | Subscription-based |

**Key Nootropic Plants:**
- Bacopa monnieri (Brahmi)
- Ginkgo biloba
- Panax ginseng
- Rhodiola rosea
- Ashwagandha (Withania somnifera)
- Huperzine A source (Huperzia serrata)

---

## 3. Neurotransmitter Pathway Databases

### 3.1 KEGG Pathway Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Dopaminergic Synapse** | https://www.kegg.jp/pathway/hsa04728 |
| **Serotonergic Synapse** | https://www.kegg.jp/pathway/hsa04726 |
| **Content** | Manually curated molecular interaction pathways |
| **API** | REST API available; bulk FTP requires license |
| **License** | Free for academic use; commercial license required for FTP |
| **Size** | ~5 GB (full pathway database) |

**Key Neurotransmitter Pathways:**
| Pathway ID | Name | Genes |
|------------|------|-------|
| hsa04728 | Dopaminergic synapse | DRD1-5, DAT1, TH, COMT, MAOA/B |
| hsa04726 | Serotonergic synapse | HTR1-7 subtypes, SLC6A4, TPH1/2 |
| hsa04724 | Glutamatergic synapse | GRIN1-3, GRM1-8, SLC1A1-3 |
| hsa04727 | GABAergic synapse | GABRA1-6, GABBR1-2, GAD1/2 |
| hsa04725 | Cholinergic synapse | CHRM1-5, CHRNA/B, ACHE |
| hsa00350 | Tyrosine metabolism | Dopamine synthesis enzymes |
| hsa00380 | Tryptophan metabolism | Serotonin synthesis enzymes |

---

### 3.2 Reactome

| Attribute | Details |
|-----------|---------|
| **URL** | https://reactome.org/ |
| **Download** | https://reactome.org/download-data |
| **Content Service** | https://reactome.org/dev/content-service |
| **Analysis Service** | https://reactome.org/dev/analysis |
| **Content** | Peer-reviewed pathway maps with molecular interactions |
| **API** | REST API for Content Service and Analysis Service |
| **License** | Creative Commons CC0 |
| **Size** | ~2 GB (full database) |

**Neurotransmitter Pathways in Reactome:**
- Neurotransmitter Release Cycle (R-HSA-112310)
- Neurotransmitter Receptors and Postsynaptic Signal Transmission
- Neurotransmitter Clearance
- G-protein coupled receptor signaling
- Ion channel activity

**API Example:**
```
# Get pathway details
https://reactome.org/ContentService/data/query/R-HSA-112310

# Pathway analysis
https://reactome.org/AnalysisService/identifiers/?pageSize=50&page=1
```

---

### 3.3 Gene Ontology (GO)

| Attribute | Details |
|-----------|---------|
| **URL** | https://geneontology.org/ |
| **Content** | Standardized gene function annotations |
| **API** | REST API, BioMart |
| **License** | CC BY 4.0 |
| **Size** | ~3 GB (full annotations) |

**Relevant GO Terms:**
- GO:0007212 - dopamine receptor signaling pathway
- GO:0007210 - serotonin receptor signaling pathway
- GO:0042136 - neurotransmitter biosynthetic process
- GO:0001505 - regulation of neurotransmitter levels
- GO:0007268 - chemical synaptic transmission

---

### 3.4 SynGO (Synaptic Gene Ontology)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.syngoportal.org/ |
| **Download** | Bulk ZIP files for all releases |
| **Content** | Expert-curated synaptic gene annotations |
| **API** | JSON/Excel bulk downloads |
| **License** | CC BY 4.0 |
| **Size** | ~100 MB |

**Features:**
- 1,100+ annotated synaptic genes
- Synaptic location ontology
- Synaptic function annotations
- Brain-specific background gene sets

---

## 4. Brain Expression Databases

### 4.1 Allen Brain Atlas (Human)

| Attribute | Details |
|-----------|---------|
| **URL** | https://human.brain-map.org/ |
| **API Documentation** | https://help.brain-map.org/display/humanbrain/API |
| **Download Page** | https://human.brain-map.org/static/download |
| **Gene Expression Portal** | https://portal.brain-map.org/gene-expression |
| **Content** | Genome-wide gene expression across 6 adult human brains, ~1000 sampling sites |
| **API** | REST API with JSON/XML/CSV output |
| **License** | Allen Institute Terms of Use (free for research) |
| **Size** | ~50 GB (microarray data for all donors) |

**API Structure:**
- RMA (RESTful Model Access) queries
- High-resolution image downloads
- Grid-level expression statistics
- Reference atlas integration

**Key Features:**
- 400-1000 sampling sites per brain
- MNI space registration
- Cytoarchitectural delineation
- Multi-histological stain data

---

### 4.2 Allen Brain Atlas (Mouse)

| Attribute | Details |
|-----------|---------|
| **URL** | https://mouse.brain-map.org/ |
| **API Documentation** | https://brain-map.org/support/documentation/api-for-mouse-brain-atlas |
| **Content** | In situ hybridization data for ~20,000 genes in adult mouse brain |
| **API** | REST API |
| **License** | Allen Institute Terms of Use |
| **Size** | ~2 TB (full ISH image data); ~500 GB (expression energy volumes) |

**Downloadable Files:**
- atlasVolume: 8-bit grayscale Nissl (25 um resolution)
- annotation: 32-bit structural annotation (25 um)
- gridAnnotation: 200 um resolution for gene expression

---

### 4.3 Allen Brain Cell Atlas (Single Cell)

| Attribute | Details |
|-----------|---------|
| **URL** | https://alleninstitute.github.io/abc_atlas_access/ |
| **AWS Bucket** | Public dataset on AWS S3 |
| **Content** | Single-cell RNA-seq across brain regions |
| **API** | AWS S3 access, anndata h5ad format |
| **License** | Allen Institute Terms of Use |
| **Size** | ~100+ GB (single-cell matrices) |

---

### 4.4 BrainSpan (Developmental Transcriptome)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.brainspan.org/ |
| **Download Page** | https://www.brainspan.org/static/download.html |
| **API Documentation** | https://help.brain-map.org/display/devhumanbrain/API |
| **Content** | RNA-seq and exon microarray across human brain development (prenatal to adult) |
| **API** | REST API with expression data access |
| **License** | Allen Institute Terms of Use |
| **Size** | ~10 GB (developmental transcriptome data) |

**Data Types:**
- RNA-Seq RPKM (gene-level): `rna_seq_genes`
- RNA-Seq RPKM (exon-level): `rna_seq_exons`
- Exon microarray (gene and probeset levels)

**Search Types:**
1. Gene Search - expression for specific genes
2. Differential Search - compare brain structures/stages
3. Correlative Search - find co-expressed genes

---

### 4.5 GTEx (Genotype-Tissue Expression)

| Attribute | Details |
|-----------|---------|
| **URL** | https://gtexportal.org/ |
| **Datasets** | https://www.gtexportal.org/home/datasets |
| **API Documentation** | https://gtexportal.org/home/apiPage |
| **API V2 Docs** | https://gtexportal.org/api/v2/redoc |
| **Content** | Gene expression across 54 tissues including 13 brain regions |
| **API** | GTEx API V2, GA4GH RNAget API |
| **License** | Open access (portal); protected access (raw data via dbGaP) |
| **Size** | Open access: ~10 GB; eQTL data: ~260 GB |

**Brain Tissues Available (13 regions):**
- Amygdala, Anterior cingulate cortex
- Caudate, Nucleus accumbens, Putamen
- Cerebellar Hemisphere, Cerebellum
- Cortex, Frontal Cortex (BA9)
- Hippocampus, Hypothalamus
- Spinal cord (cervical c-1)
- Substantia nigra

**R Packages:**
- `gtexr` - Query GTEx Portal API V2
- `gtexRNA` - Retrieve tissue-specific expression

---

### 4.6 Brain.GMT (Curated Brain Gene Sets)

| Attribute | Details |
|-----------|---------|
| **Reference** | PMC11030436 |
| **Content** | 918 curated gene sets for nervous system function, tissue, cell types |
| **Species** | Human, Mouse, Rat |
| **Format** | GMT format for GSEA/limma/edgeR |
| **License** | Academic use |
| **Size** | ~10 MB |

---

## 5. Cognitive Trait Genetics

### 5.1 COGENT (Cognitive Genomics Consortium)

| Attribute | Details |
|-----------|---------|
| **URL** | No dedicated portal; data shared via publications |
| **Key Publications** | Nature Communications (2018), Molecular Psychiatry (2017) |
| **Content** | GWAS meta-analysis of general cognitive function |
| **API** | No API; summary statistics via publication supplements |
| **License** | Academic research use |
| **Size** | Summary statistics: ~500 MB per study |

**Key Studies:**

**COGENT 2017 (Molecular Psychiatry):**
- N = 35,298 individuals
- 24 cohorts, European ancestry
- 2 genome-wide significant SNP loci
- 3 gene-based significant loci
- SNP heritability: 21.5%

**COGENT + CHARGE + UK Biobank 2018:**
- N = 300,486 individuals
- 148 genome-wide significant loci
- 709 genes associated
- Up to 4.3% variance explained by polygenic scores

**Latest 2025 Study:**
- Genomic-SEM analysis
- 3,842 genome-wide significant loci
- 275 novel loci
- 13 high-confidence causal genes

**Genetic Correlations Found:**
- Educational attainment (positive)
- Schizophrenia (negative with cognitive function)
- Bipolar disorder, depression
- Birth length/weight
- Smoking behavior
- Openness personality trait

---

### 5.2 UK Biobank Cognitive Data

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ukbiobank.ac.uk/ |
| **Application** | https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/ |
| **Content** | Cognitive assessments + genetic data for 500,000 participants |
| **API** | Research Analysis Platform access |
| **License** | UK Biobank Access Agreement |
| **Size** | Cognitive phenotypes: ~1 GB; Genetic data: multiple TB |

**Cognitive Measures:**
- Fluid intelligence
- Reaction time
- Numeric memory
- Prospective memory
- Pairs matching
- Trail making
- Symbol digit substitution

---

### 5.3 GWAS Catalog - Cognitive Traits

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Search** | https://www.ebi.ac.uk/gwas/search?query=cognitive |
| **Content** | Curated GWAS associations for cognitive phenotypes |
| **API** | REST API, bulk downloads |
| **License** | Open access |
| **Size** | Full catalog: ~2 GB |

**Cognitive Trait Studies in GWAS Catalog:**
- Intelligence
- Educational attainment
- Cognitive performance
- Memory
- Processing speed
- Executive function

---

## 6. Sleep Genetics Databases

### 6.1 Sleep Disorders Knowledge Portal

| Attribute | Details |
|-----------|---------|
| **URL** | http://sleepdisordergenetics.org/ |
| **Data Download** | http://sleepdisordergenetics.org/informational/data |
| **Content** | GWAS summary statistics for sleep phenotypes |
| **API** | Web-based queries; bulk downloads |
| **License** | Academic research use |
| **Size** | ~5 GB (summary statistics) |

**Available Phenotypes:**
- Insomnia
- Sleep duration
- Chronotype (morningness/eveningness)
- Daytime sleepiness
- Sleep apnea
- Restless legs syndrome
- Narcolepsy

---

### 6.2 GWAS Catalog - Sleep/Circadian

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/gwas/search?query=sleep |
| **Insomnia** | https://www.ebi.ac.uk/gwas/search?query=insomnia |
| **Content** | 4,236 genes with SNPs associated with circadian traits |
| **API** | REST API |
| **License** | Open access |

**Key Sleep GWAS Studies:**
| Study | N | Loci | Phenotype |
|-------|---|------|-----------|
| Jansen et al. 2019 | 1,331,010 | 202 | Insomnia |
| Watanabe et al. 2022 | 2,365,010 | 554 | Insomnia |
| Jones et al. 2019 | 697,828 | 351 | Chronotype |
| Dashti et al. 2019 | 446,118 | 78 | Sleep duration |

---

### 6.3 CNCR Summary Statistics

| Attribute | Details |
|-----------|---------|
| **URL** | https://cncr.nl/research/summary_statistics/ |
| **Content** | GWAS summary statistics from VU Amsterdam |
| **License** | Academic use; 23andMe data requires separate DTA |
| **Size** | ~2 GB per study |

**Available Studies:**
- Insomnia (excluding 23andMe)
- Chronotype
- Sleep duration
- Mental health phenotypes

---

### 6.4 UK Biobank Sleep Data

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ukbiobank.ac.uk/ |
| **Sleep Category** | Category 100057 |
| **Content** | 160+ sleep questions from ~180,000 respondents |
| **API** | Research Analysis Platform |
| **License** | UK Biobank Access Agreement |
| **Size** | ~500 MB (sleep phenotypes) |

**Sleep Phenotypes:**
- Sleep duration
- Chronotype (morning/evening person)
- Insomnia/sleeplessness
- Snoring
- Daytime napping/dozing
- Sleep apnea indicators
- Accelerometer-derived sleep measures

---

### 6.5 OpenGWAS

| Attribute | Details |
|-----------|---------|
| **URL** | https://opengwas.io/ |
| **Content** | 50,000+ GWAS summary datasets including sleep traits |
| **API** | REST API with account access |
| **License** | Varies by dataset |
| **Size** | Full database: multiple TB |

**Sleep-Related Datasets:**
- Search "insomnia", "sleep", "chronotype", "circadian"
- Programmatic batch queries available
- PheWAS capabilities

---

### 6.6 Circadian Gene Resources

| Resource | URL | Content |
|----------|-----|---------|
| KEGG Circadian Rhythm | hsa04710 | Core clock genes pathway |
| CircaDB | http://circadb.hogeneschlab.org/ | Circadian expression profiles |
| CGDB | Circadian Gene Database | Time-series expression |

**Core Circadian Genes:**
- CLOCK, BMAL1 (ARNTL)
- PER1, PER2, PER3
- CRY1, CRY2
- REV-ERBa/b (NR1D1/2)
- RORa/b/c
- CKI epsilon/delta (CSNK1E/D)

---

## 7. Summary Comparison Table

| Database | Category | API | License | Size | Primary Use |
|----------|----------|-----|---------|------|-------------|
| **PGC** | Psychiatric GWAS | No | Academic | TB (ind.) / GB (sumstats) | Psychiatric disorder genetics |
| **ChEMBL** | Compounds | REST | CC BY-SA 3.0 | 15 GB | Nootropic compound data |
| **DrugBank** | Drugs | Commercial | Academic free | 2 GB | Drug pharmacology |
| **PubChem** | Compounds | REST | Public domain | 100+ GB | Chemical structures |
| **KEGG** | Pathways | REST/FTP | Academic free | 5 GB | Neurotransmitter pathways |
| **Reactome** | Pathways | REST | CC0 | 2 GB | Molecular interactions |
| **Allen Brain** | Expression | REST | Free research | 50+ GB | Brain gene expression |
| **BrainSpan** | Expression | REST | Free research | 10 GB | Developmental expression |
| **GTEx** | Expression | REST | Open/dbGaP | 260+ GB | Multi-tissue expression |
| **COGENT** | Cognition GWAS | No | Academic | 500 MB | Cognitive genetics |
| **UK Biobank** | Multi-phenotype | Platform | Agreement | TB+ | Comprehensive phenotypes |
| **Sleep Portal** | Sleep GWAS | Web | Academic | 5 GB | Sleep genetics |
| **GWAS Catalog** | All traits | REST | Open | 2 GB | Curated associations |
| **OpenGWAS** | All traits | REST | Varies | TB+ | GWAS aggregation |

---

## Recommended Integration Strategy

### Priority 1: Core Databases
1. **PGC** - Psychiatric GWAS summary statistics
2. **ChEMBL** - Nootropic compounds and targets
3. **GTEx** - Brain expression with eQTLs
4. **GWAS Catalog** - Curated cognitive/sleep associations

### Priority 2: Specialized Resources
5. **Allen Brain Atlas** - Detailed brain expression maps
6. **Reactome** - Neurotransmitter pathways
7. **Sleep Disorders Knowledge Portal** - Sleep phenotypes
8. **COGENT summary statistics** - Cognitive trait genetics

### Priority 3: Supporting Databases
9. **BrainSpan** - Developmental context
10. **SynGO** - Synaptic gene annotations
11. **DrugBank** - Drug interactions
12. **UK Biobank** (if access obtained) - Deep phenotyping

### API Integration Notes

**Well-documented REST APIs:**
- ChEMBL, Reactome, GTEx, Allen Brain Atlas, GWAS Catalog

**Bulk download only:**
- PGC summary statistics, COGENT, BrainSpan

**Requires special access:**
- UK Biobank (application), dbGaP (DAR), PGC individual data (DAC approval)

---

## References

1. Psychiatric Genomics Consortium. https://pgc.unc.edu/
2. ChEMBL Database. https://www.ebi.ac.uk/chembl/
3. KEGG Pathway Database. https://www.genome.jp/kegg/
4. Reactome Pathway Database. https://reactome.org/
5. Allen Brain Atlas. https://brain-map.org/
6. GTEx Portal. https://gtexportal.org/
7. GWAS Catalog. https://www.ebi.ac.uk/gwas/
8. Trampush et al. (2017). GWAS meta-analysis reveals novel loci for cognitive function. Mol Psychiatry.
9. Davies et al. (2018). Study of 300,486 individuals identifies 148 genetic loci. Nat Commun.
10. Jansen et al. (2019). Genome-wide analysis of insomnia. Nat Genet.
11. Jones et al. (2019). Genome-wide association analyses of chronotype. Nat Commun.
