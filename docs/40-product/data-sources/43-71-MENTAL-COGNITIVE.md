# Mental Health & Cognitive Data Sources

**Document ID:** 43-71-MENTAL-COGNITIVE
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-00-INDEX.md](./43-00-INDEX.md)

---

## TL;DR

Comprehensive catalog of 35+ databases covering psychiatric genetics (PGC, dbGaP), nootropic compounds (ChEMBL, DrugBank, PubChem), neurotransmitter pathways (KEGG, Reactome, SynGO), brain expression atlases (Allen Brain, GTEx, BrainSpan), cognitive trait genetics (COGENT, UK Biobank), and sleep/circadian genetics (Sleep Disorders Knowledge Portal, OpenGWAS). Combined storage estimate: ~500 GB for summary statistics, 3+ TB for individual-level data.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary psychiatric genetics | PGC summary statistics | Largest meta-analyses, open access summaries | Jan 2026 |
| Nootropics compound source | ChEMBL + PubChem | CC licenses, comprehensive APIs | Jan 2026 |
| Brain expression priority | GTEx + Allen Brain Atlas | REST APIs, research-friendly licenses | Jan 2026 |
| Sleep genetics source | GWAS Catalog + OpenGWAS | Open access, programmatic queries | Jan 2026 |
| Individual-level data | Defer to Tier 3 | DAC approval complexity | Jan 2026 |

---

## Database Catalog

### 1. Psychiatric Genetics Databases

#### 1.1 Psychiatric Genomics Consortium (PGC)

| Field | Value |
|-------|-------|
| **URL** | https://pgc.unc.edu/ |
| **Download** | https://pgc.unc.edu/for-researchers/download-results/ |
| **Data Portal** | https://pgc.unc.edu/for-researchers/data-access-committee/data-access-portal/ |
| **GitHub** | https://github.com/psychiatric-genomics-consortium |
| **Content** | GWAS meta-analysis for major depression, bipolar, ADHD, schizophrenia, Tourette, OCD, autism, anorexia, PTSD, substance use |
| **Records** | Summary stats: varies by disorder; Individual-level: multiple TB |
| **License** | Academic research; commercial requires PGC DAC permission |
| **API** | No dedicated API; web portal + controlled access requests |
| **Update Frequency** | New studies published periodically |
| **Priority** | Tier 1 (summary stats) / Tier 3 (individual data) |
| **Storage Estimate** | Summary: 5-10 GB; Individual: 2+ TB |

**Key Studies:**
- Schizophrenia: 150+ loci (N > 300,000)
- Major Depression: 100+ loci (N > 1 million)
- Bipolar Disorder: 64 genome-wide significant loci

**Access Requirements (Individual Data):**
1. LISA/GCC cluster account
2. Signed analyst memo
3. Signed WTCCC memo
4. Approved research proposal
5. Workgroup chair approval

---

#### 1.2 dbGaP Psychiatric Collections

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/gap/ |
| **PGC Collection** | https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/collection.cgi?study_id=phs001254.v1.p1 |
| **Content** | Individual-level genotype and phenotype data for psychiatric disorders |
| **Records** | Hundreds of GB across psychiatric collections |
| **License** | NIH Data Use Certification required |
| **API** | dbGaP API for metadata; data requires DAR approval |
| **Update Frequency** | Continuous submissions |
| **Priority** | Tier 3 |
| **Storage Estimate** | 500+ GB |

---

#### 1.3 NIMH Genetics Repository (NRGR)

| Field | Value |
|-------|-------|
| **URL** | https://www.nimhgenetics.org/ |
| **Bipolar Data** | https://www.nimhgenetics.org/resources/genomic-data/pgc-bp |
| **Content** | DNA samples, cell lines, genomic data for psychiatric research |
| **Records** | 200,000+ samples across psychiatric phenotypes |
| **License** | Academic research; requires NIMH approval |
| **API** | No public API; web-based ordering |
| **Update Frequency** | Continuous collection |
| **Priority** | Tier 3 |
| **Storage Estimate** | N/A (sample repository) |

---

### 2. Nootropics & Cognitive Enhancement Compounds

#### 2.1 ChEMBL

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/chembl/ |
| **API Documentation** | https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services |
| **API Live Docs** | https://www.ebi.ac.uk/chembl/api/data/docs |
| **Content** | 2.4M+ compounds, 1.5M+ assays, 20M+ activities, drug-target interactions |
| **Records** | 2.4M+ compounds |
| **License** | CC BY-SA 3.0 |
| **API** | REST API with JSON/XML/CSV; substructure & similarity search |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 |
| **Storage Estimate** | 15 GB (MySQL); 3 GB (SDF) |

**Key Nootropic Targets:**
- Acetylcholinesterase (CHEMBL220)
- NMDA receptors (CHEMBL1907602)
- Dopamine receptors D1-D5
- Serotonin receptors (5-HT subtypes)
- Alpha-7 nicotinic receptor (CHEMBL3193)

**Example Queries:**
```
# Approved cognitive drugs
https://www.ebi.ac.uk/chembl/api/data/molecule?max_phase=4

# Compounds targeting acetylcholinesterase
https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id=CHEMBL220
```

---

#### 2.2 DrugBank

| Field | Value |
|-------|-------|
| **URL** | https://go.drugbank.com/ |
| **API Documentation** | https://docs.drugbank.com/ |
| **Content** | 15,000+ drugs with pharmacology, targets, interactions |
| **Records** | 15,000+ drug entries |
| **License** | Free academic non-commercial; commercial license for API |
| **API** | Clinical API (commercial); Academic download available |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 2 |
| **Storage Estimate** | 2 GB (XML download) |

**Nootropic-Related Entries:**
- Piracetam and racetam derivatives
- Modafinil/Armodafinil
- Methylphenidate
- Donepezil, Rivastigmine (AChE inhibitors)
- Memantine (NMDA antagonist)

---

#### 2.3 PubChem

| Field | Value |
|-------|-------|
| **URL** | https://pubchem.ncbi.nlm.nih.gov/ |
| **Downloads** | https://pubchem.ncbi.nlm.nih.gov/docs/downloads |
| **Content** | 118M+ compounds, 313M+ substances, biological assay data |
| **Records** | 118M+ compounds |
| **License** | Public Domain (US government) |
| **API** | PUG REST API, PUG-View, E-utilities |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 |
| **Storage Estimate** | 100+ GB (full database) |

**Nootropics Search Strategy:**
- Search MeSH term "Nootropic Agents"
- Filter by pharmacological class
- Cross-reference cognitive/memory assays

---

#### 2.4 ChEBI (Chemical Entities of Biological Interest)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/chebi/ |
| **Nootropic Ontology** | https://www.ebi.ac.uk/chebi/searchId.do?chebiId=66980 |
| **Content** | Chemical ontology with hierarchical classification |
| **Records** | 60,000+ entities |
| **License** | CC0 1.0 |
| **API** | Web services available |
| **Update Frequency** | Monthly |
| **Priority** | Tier 2 |
| **Storage Estimate** | 500 MB |

**CHEBI:66980 (Nootropic):**
- Definition: Compounds improving mental functions (cognition, memory, intelligence, motivation, attention)
- Includes natural and synthetic cognitive enhancers

---

#### 2.5 NAPRALERT (Natural Products)

| Field | Value |
|-------|-------|
| **URL** | https://napralert.org/ |
| **Content** | 200,000+ compounds from natural sources including cognitive-enhancing plants |
| **Records** | 200,000+ compounds |
| **License** | Subscription-based |
| **API** | No public API |
| **Update Frequency** | Regular updates |
| **Priority** | Tier 3 |
| **Storage Estimate** | N/A (subscription access) |

**Key Nootropic Plants:**
- Bacopa monnieri (Brahmi)
- Ginkgo biloba
- Panax ginseng
- Rhodiola rosea
- Ashwagandha (Withania somnifera)
- Huperzine A source (Huperzia serrata)

---

### 3. Neurotransmitter Pathway Databases

#### 3.1 KEGG Pathway Database

| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Dopaminergic Synapse** | https://www.kegg.jp/pathway/hsa04728 |
| **Serotonergic Synapse** | https://www.kegg.jp/pathway/hsa04726 |
| **Content** | Manually curated molecular interaction pathways |
| **Records** | 500+ human pathways |
| **License** | Free academic use; commercial license for FTP |
| **API** | REST API; bulk FTP requires license |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 |
| **Storage Estimate** | 5 GB (full database) |

**Key Neurotransmitter Pathways:**

| Pathway ID | Name | Key Genes |
|------------|------|-----------|
| hsa04728 | Dopaminergic synapse | DRD1-5, DAT1, TH, COMT, MAOA/B |
| hsa04726 | Serotonergic synapse | HTR1-7, SLC6A4, TPH1/2 |
| hsa04724 | Glutamatergic synapse | GRIN1-3, GRM1-8, SLC1A1-3 |
| hsa04727 | GABAergic synapse | GABRA1-6, GABBR1-2, GAD1/2 |
| hsa04725 | Cholinergic synapse | CHRM1-5, CHRNA/B, ACHE |
| hsa00350 | Tyrosine metabolism | Dopamine synthesis enzymes |
| hsa00380 | Tryptophan metabolism | Serotonin synthesis enzymes |

---

#### 3.2 Reactome

| Field | Value |
|-------|-------|
| **URL** | https://reactome.org/ |
| **Download** | https://reactome.org/download-data |
| **Content Service** | https://reactome.org/dev/content-service |
| **Analysis Service** | https://reactome.org/dev/analysis |
| **Content** | Peer-reviewed pathway maps with molecular interactions |
| **Records** | 2,600+ human pathways |
| **License** | CC0 |
| **API** | REST API for Content and Analysis services |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 |
| **Storage Estimate** | 2 GB |

**Neurotransmitter Pathways:**
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

#### 3.3 Gene Ontology (GO)

| Field | Value |
|-------|-------|
| **URL** | https://geneontology.org/ |
| **Content** | Standardized gene function annotations |
| **Records** | 45,000+ GO terms |
| **License** | CC BY 4.0 |
| **API** | REST API, BioMart |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 |
| **Storage Estimate** | 3 GB (full annotations) |

**Relevant GO Terms:**
- GO:0007212 - dopamine receptor signaling pathway
- GO:0007210 - serotonin receptor signaling pathway
- GO:0042136 - neurotransmitter biosynthetic process
- GO:0001505 - regulation of neurotransmitter levels
- GO:0007268 - chemical synaptic transmission

---

#### 3.4 SynGO (Synaptic Gene Ontology)

| Field | Value |
|-------|-------|
| **URL** | https://www.syngoportal.org/ |
| **Download** | Bulk ZIP files for all releases |
| **Content** | Expert-curated synaptic gene annotations |
| **Records** | 1,100+ annotated synaptic genes |
| **License** | CC BY 4.0 |
| **API** | JSON/Excel bulk downloads |
| **Update Frequency** | Annual |
| **Priority** | Tier 2 |
| **Storage Estimate** | 100 MB |

**Features:**
- Synaptic location ontology
- Synaptic function annotations
- Brain-specific background gene sets
- Integration with GO

---

### 4. Brain Expression Databases

#### 4.1 Allen Brain Atlas (Human)

| Field | Value |
|-------|-------|
| **URL** | https://human.brain-map.org/ |
| **API Documentation** | https://help.brain-map.org/display/humanbrain/API |
| **Download** | https://human.brain-map.org/static/download |
| **Gene Expression Portal** | https://portal.brain-map.org/gene-expression |
| **Content** | Genome-wide expression across 6 adult human brains, ~1000 sampling sites |
| **Records** | 6 donor brains, 400-1000 sites each |
| **License** | Allen Institute Terms of Use (free research) |
| **API** | REST API (RMA queries) with JSON/XML/CSV |
| **Update Frequency** | Static reference dataset |
| **Priority** | Tier 1 |
| **Storage Estimate** | 50 GB (microarray data) |

**Key Features:**
- MNI space registration
- Cytoarchitectural delineation
- Multi-histological stain data
- High-resolution image downloads

---

#### 4.2 Allen Brain Atlas (Mouse)

| Field | Value |
|-------|-------|
| **URL** | https://mouse.brain-map.org/ |
| **API Documentation** | https://brain-map.org/support/documentation/api-for-mouse-brain-atlas |
| **Content** | In situ hybridization for ~20,000 genes in adult mouse brain |
| **Records** | 20,000+ genes |
| **License** | Allen Institute Terms of Use |
| **API** | REST API |
| **Update Frequency** | Static reference |
| **Priority** | Tier 2 |
| **Storage Estimate** | 2 TB (full ISH); 500 GB (expression volumes) |

**Downloadable Files:**
- atlasVolume: 8-bit grayscale Nissl (25 um)
- annotation: 32-bit structural (25 um)
- gridAnnotation: 200 um for gene expression

---

#### 4.3 Allen Brain Cell Atlas (Single Cell)

| Field | Value |
|-------|-------|
| **URL** | https://alleninstitute.github.io/abc_atlas_access/ |
| **AWS Bucket** | Public S3 dataset |
| **Content** | Single-cell RNA-seq across brain regions |
| **Records** | Millions of cells |
| **License** | Allen Institute Terms of Use |
| **API** | AWS S3 access, anndata h5ad format |
| **Update Frequency** | Periodic releases |
| **Priority** | Tier 2 |
| **Storage Estimate** | 100+ GB |

---

#### 4.4 BrainSpan (Developmental Transcriptome)

| Field | Value |
|-------|-------|
| **URL** | https://www.brainspan.org/ |
| **Download** | https://www.brainspan.org/static/download.html |
| **API Documentation** | https://help.brain-map.org/display/devhumanbrain/API |
| **Content** | RNA-seq and exon microarray across human brain development (prenatal to adult) |
| **Records** | 16 brain structures, 31 developmental stages |
| **License** | Allen Institute Terms of Use |
| **API** | REST API |
| **Update Frequency** | Static reference |
| **Priority** | Tier 2 |
| **Storage Estimate** | 10 GB |

**Data Types:**
- RNA-Seq RPKM (gene-level): `rna_seq_genes`
- RNA-Seq RPKM (exon-level): `rna_seq_exons`
- Exon microarray (gene and probeset)

**Search Types:**
1. Gene Search - expression for specific genes
2. Differential Search - compare structures/stages
3. Correlative Search - find co-expressed genes

---

#### 4.5 GTEx (Genotype-Tissue Expression)

| Field | Value |
|-------|-------|
| **URL** | https://gtexportal.org/ |
| **Datasets** | https://www.gtexportal.org/home/datasets |
| **API Documentation** | https://gtexportal.org/home/apiPage |
| **API V2 Docs** | https://gtexportal.org/api/v2/redoc |
| **Content** | Gene expression across 54 tissues including 13 brain regions |
| **Records** | 17,382 samples, 948 donors |
| **License** | Open access (portal); protected (raw via dbGaP) |
| **API** | GTEx API V2, GA4GH RNAget API |
| **Update Frequency** | Major releases annually |
| **Priority** | Tier 1 |
| **Storage Estimate** | Open: 10 GB; eQTL: 260 GB |

**Brain Tissues (13 regions):**
- Amygdala
- Anterior cingulate cortex
- Caudate
- Cerebellar Hemisphere
- Cerebellum
- Cortex
- Frontal Cortex (BA9)
- Hippocampus
- Hypothalamus
- Nucleus accumbens
- Putamen
- Spinal cord (cervical c-1)
- Substantia nigra

**R Packages:**
- `gtexr` - Query GTEx Portal API V2
- `gtexRNA` - Retrieve tissue-specific expression

---

#### 4.6 Brain.GMT (Curated Brain Gene Sets)

| Field | Value |
|-------|-------|
| **Reference** | PMC11030436 |
| **Content** | 918 curated gene sets for nervous system function |
| **Records** | 918 gene sets |
| **Species** | Human, Mouse, Rat |
| **License** | Academic use |
| **API** | GMT format downloads |
| **Update Frequency** | Static |
| **Priority** | Tier 2 |
| **Storage Estimate** | 10 MB |

---

### 5. Cognitive Trait Genetics

#### 5.1 COGENT (Cognitive Genomics Consortium)

| Field | Value |
|-------|-------|
| **URL** | No dedicated portal; via publications |
| **Publications** | Nature Communications (2018), Molecular Psychiatry (2017) |
| **Content** | GWAS meta-analysis of general cognitive function |
| **Records** | Summary statistics for 300K+ individuals |
| **License** | Academic research use |
| **API** | No API; supplementary data downloads |
| **Update Frequency** | Publication-driven |
| **Priority** | Tier 1 |
| **Storage Estimate** | 500 MB per study |

**Key Studies:**

| Study | N | Loci | Key Finding |
|-------|---|------|-------------|
| COGENT 2017 | 35,298 | 2 GWS SNP, 3 gene-based | SNP heritability 21.5% |
| COGENT+UK Biobank 2018 | 300,486 | 148 | 709 genes, 4.3% variance PRS |
| Latest 2025 | N/A | 3,842 GWS | 275 novel, 13 causal genes |

**Genetic Correlations:**
- Educational attainment (positive)
- Schizophrenia (negative)
- Bipolar disorder, depression
- Birth length/weight
- Smoking behavior
- Openness personality trait

---

#### 5.2 UK Biobank Cognitive Data

| Field | Value |
|-------|-------|
| **URL** | https://www.ukbiobank.ac.uk/ |
| **Application** | https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/ |
| **Content** | Cognitive assessments + genetic data for 500,000 participants |
| **Records** | 500,000 participants |
| **License** | UK Biobank Access Agreement |
| **API** | Research Analysis Platform |
| **Update Frequency** | Ongoing collection |
| **Priority** | Tier 2 (requires application) |
| **Storage Estimate** | Cognitive: 1 GB; Genetic: multiple TB |

**Cognitive Measures:**
- Fluid intelligence
- Reaction time
- Numeric memory
- Prospective memory
- Pairs matching
- Trail making
- Symbol digit substitution

---

#### 5.3 GWAS Catalog - Cognitive Traits

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Search** | https://www.ebi.ac.uk/gwas/search?query=cognitive |
| **Content** | Curated GWAS associations for cognitive phenotypes |
| **Records** | 500+ cognitive trait studies |
| **License** | Open access |
| **API** | REST API, bulk downloads |
| **Update Frequency** | Continuous curation |
| **Priority** | Tier 1 |
| **Storage Estimate** | 2 GB (full catalog) |

**Cognitive Traits Indexed:**
- Intelligence
- Educational attainment
- Cognitive performance
- Memory
- Processing speed
- Executive function

---

### 6. Sleep & Circadian Genetics

#### 6.1 Sleep Disorders Knowledge Portal

| Field | Value |
|-------|-------|
| **URL** | http://sleepdisordergenetics.org/ |
| **Data Download** | http://sleepdisordergenetics.org/informational/data |
| **Content** | GWAS summary statistics for sleep phenotypes |
| **Records** | Multiple large-scale GWAS |
| **License** | Academic research use |
| **API** | Web-based queries; bulk downloads |
| **Update Frequency** | Publication-driven |
| **Priority** | Tier 1 |
| **Storage Estimate** | 5 GB |

**Available Phenotypes:**
- Insomnia
- Sleep duration
- Chronotype (morningness/eveningness)
- Daytime sleepiness
- Sleep apnea
- Restless legs syndrome
- Narcolepsy

---

#### 6.2 GWAS Catalog - Sleep/Circadian

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/search?query=sleep |
| **Insomnia** | https://www.ebi.ac.uk/gwas/search?query=insomnia |
| **Content** | 4,236 genes with SNPs associated with circadian traits |
| **Records** | 4,236+ associated genes |
| **License** | Open access |
| **API** | REST API |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 |
| **Storage Estimate** | Included in GWAS Catalog |

**Key Sleep GWAS Studies:**

| Study | N | Loci | Phenotype |
|-------|---|------|-----------|
| Jansen et al. 2019 | 1,331,010 | 202 | Insomnia |
| Watanabe et al. 2022 | 2,365,010 | 554 | Insomnia |
| Jones et al. 2019 | 697,828 | 351 | Chronotype |
| Dashti et al. 2019 | 446,118 | 78 | Sleep duration |

---

#### 6.3 CNCR Summary Statistics

| Field | Value |
|-------|-------|
| **URL** | https://cncr.nl/research/summary_statistics/ |
| **Content** | GWAS summary statistics from VU Amsterdam |
| **Records** | Multiple phenotypes |
| **License** | Academic use; 23andMe data requires DTA |
| **API** | Download only |
| **Update Frequency** | Publication-driven |
| **Priority** | Tier 2 |
| **Storage Estimate** | 2 GB per study |

**Available Studies:**
- Insomnia (excluding 23andMe)
- Chronotype
- Sleep duration
- Mental health phenotypes

---

#### 6.4 UK Biobank Sleep Data

| Field | Value |
|-------|-------|
| **URL** | https://www.ukbiobank.ac.uk/ |
| **Sleep Category** | Category 100057 |
| **Content** | 160+ sleep questions from ~180,000 respondents |
| **Records** | 180,000+ respondents |
| **License** | UK Biobank Access Agreement |
| **API** | Research Analysis Platform |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 2 |
| **Storage Estimate** | 500 MB (sleep phenotypes) |

**Sleep Phenotypes:**
- Sleep duration
- Chronotype (morning/evening)
- Insomnia/sleeplessness
- Snoring
- Daytime napping/dozing
- Sleep apnea indicators
- Accelerometer-derived measures

---

#### 6.5 OpenGWAS

| Field | Value |
|-------|-------|
| **URL** | https://opengwas.io/ |
| **Content** | 50,000+ GWAS summary datasets including sleep traits |
| **Records** | 50,000+ datasets |
| **License** | Varies by dataset |
| **API** | REST API with account |
| **Update Frequency** | Continuous aggregation |
| **Priority** | Tier 1 |
| **Storage Estimate** | Multiple TB (full) |

**Sleep-Related Searches:**
- insomnia
- sleep
- chronotype
- circadian
- Programmatic batch queries
- PheWAS capabilities

---

#### 6.6 Circadian Gene Resources

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

## Summary Comparison Table

| Database | Category | API | License | Size | Primary Use |
|----------|----------|-----|---------|------|-------------|
| **PGC** | Psychiatric GWAS | No | Academic | GB-TB | Psychiatric genetics |
| **ChEMBL** | Compounds | REST | CC BY-SA 3.0 | 15 GB | Nootropic compounds |
| **DrugBank** | Drugs | Commercial | Academic free | 2 GB | Drug pharmacology |
| **PubChem** | Compounds | REST | Public domain | 100+ GB | Chemical structures |
| **KEGG** | Pathways | REST/FTP | Academic free | 5 GB | Neurotransmitter pathways |
| **Reactome** | Pathways | REST | CC0 | 2 GB | Molecular interactions |
| **SynGO** | Synaptic | Download | CC BY 4.0 | 100 MB | Synaptic genes |
| **Allen Brain** | Expression | REST | Free research | 50+ GB | Brain expression |
| **BrainSpan** | Expression | REST | Free research | 10 GB | Developmental expression |
| **GTEx** | Expression | REST | Open/dbGaP | 260+ GB | Multi-tissue eQTL |
| **COGENT** | Cognition GWAS | No | Academic | 500 MB | Cognitive genetics |
| **UK Biobank** | Multi-phenotype | Platform | Agreement | TB+ | Comprehensive phenotypes |
| **Sleep Portal** | Sleep GWAS | Web | Academic | 5 GB | Sleep genetics |
| **GWAS Catalog** | All traits | REST | Open | 2 GB | Curated associations |
| **OpenGWAS** | All traits | REST | Varies | TB+ | GWAS aggregation |

---

## Integration Recommendations

### Priority 1: MVP (Tier 1)

| Source | Category | Rationale |
|--------|----------|-----------|
| PGC Summary Stats | Psychiatric GWAS | Largest psychiatric genetics resource |
| ChEMBL | Nootropic Compounds | CC license, comprehensive API |
| PubChem | Compounds | Public domain, complete coverage |
| GTEx | Brain Expression | 13 brain regions with eQTL |
| Reactome | Pathways | CC0, excellent API |
| KEGG | Pathways | Essential neurotransmitter maps |
| GWAS Catalog | Cognitive/Sleep | Open access, curated |
| Sleep Portal | Sleep GWAS | Specialized sleep genetics |
| OpenGWAS | Aggregated GWAS | Programmatic batch access |

### Priority 2: Post-MVP (Tier 2)

| Source | Category | Rationale |
|--------|----------|-----------|
| DrugBank | Drugs | Detailed pharmacology |
| Allen Brain Atlas | Expression | Spatial brain expression |
| BrainSpan | Development | Developmental trajectories |
| SynGO | Synaptic | Specialized annotations |
| UK Biobank | Deep phenotyping | Application required |
| CNCR | Sleep GWAS | VU Amsterdam studies |
| ChEBI | Ontology | Chemical classification |

### Priority 3: Future (Tier 3)

| Source | Category | Notes |
|--------|----------|-------|
| dbGaP Individual Data | Psychiatric | DAR approval required |
| NIMH Repository | Samples | Research application |
| NAPRALERT | Natural Products | Subscription |
| Allen Mouse Brain | Mouse Expression | 2 TB storage |

---

## API Integration Summary

**Well-documented REST APIs:**
- ChEMBL, PubChem, Reactome, GTEx, Allen Brain, GWAS Catalog, OpenGWAS

**Bulk download only:**
- PGC summary statistics, COGENT, BrainSpan, SynGO, CNCR

**Requires special access:**
- UK Biobank (application)
- dbGaP (DAR approval)
- PGC individual data (DAC)
- NAPRALERT (subscription)
- DrugBank API (commercial)

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-00-INDEX](./43-00-INDEX.md) | Parent navigation |
| [43-DATA-SOURCES](../43-DATA-SOURCES.md) | Overview context |
| [43-41-PATHWAYS-PRIMARY](./43-41-PATHWAYS-PRIMARY.md) | Shared pathway sources |
| [43-51-PHARMACEUTICALS](./43-51-PHARMACEUTICALS.md) | Drug/compound overlap |
| [44-ARCHITECTURE](../44-ARCHITECTURE.md) | Database design |
| [45-DATA-MODEL](../45-DATA-MODEL.md) | Entity relationships |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial migration from research.old/data-sources-mental-cognitive.md |
