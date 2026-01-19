# Sleep, Longevity, and Nutrigenomics Data Sources

Research compiled for the Gene Platform - specialized genetics databases covering circadian rhythms, aging/longevity, nutrigenomics, vitamin metabolism, and food-gene interactions.

---

## 1. Circadian Rhythm Genetics Databases

### 1.1 GWAS Catalog (Chronotype/Sleep Data)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Content** | Comprehensive catalog of GWAS studies including 351+ loci for chronotype ("morning person" trait), sleep duration, insomnia. Contains data from 697,828+ UK Biobank and 23andMe participants. |
| **API** | REST API: https://www.ebi.ac.uk/gwas/rest/docs/api - Returns JSON format. R package (gwasrapidd) and Python package (pandasGWAS) available for programmatic access. |
| **License** | Open access (CC0 for catalog, study-specific for summary statistics) |
| **Size** | 300,000+ variant-trait associations; chronotype-specific: 351 loci |
| **Key Genes** | CLOCK, BMAL1 (ARNTL), PER1-3, CRY1-2, NPAS2, CSNK1D, CSNK1E |
| **Download** | https://www.ebi.ac.uk/gwas/docs/file-downloads |

### 1.2 CircaDB (Circadian Gene Expression)

| Attribute | Details |
|-----------|---------|
| **URL** | http://circadb.hogeneschlab.org/ |
| **Content** | Database of circadian transcriptional profiles from time-course expression experiments in mice and humans. Evaluates transcripts using JTK_Cycle, Lomb Scargle, and DeLichtenberg algorithms. |
| **API** | No REST API. Web interface with Google Charts visualization. |
| **License** | Open source (GitHub: https://github.com/itmat/circadb) |
| **Size** | ~3,000 cycling genes with expression profiles |
| **Download** | Raw data downloadable; processed results via web interface |
| **Notes** | Last major update 2014; integrates with BioGPS |

### 1.3 CirGRDB (Circadian Genes and Regulators)

| Attribute | Details |
|-----------|---------|
| **URL** | http://cirgrdb.biols.ac.cn |
| **Content** | Genome-wide deciphering of circadian genes and regulators. Temporal expression patterns across 37 human/mouse tissues. Includes associations with sleep disorders, aging, and tumors. Eight types of regulatory data: transcription factors, histone modifications, chromatin accessibility, enhancer RNAs, miRNAs, RNA-binding proteins, RNA editing, RNA methylation. |
| **API** | No dedicated API |
| **License** | Academic use |
| **Size** | 4,936+ genome-wide assays integrated |
| **Download** | Expression and regulatory data as plain text files |

### 1.4 Sleep Disorders Knowledge Portal

| Attribute | Details |
|-----------|---------|
| **URL** | https://sleep.hugeamp.org/ |
| **Content** | GWAS summary statistics for sleep phenotypes including insomnia, sleep duration, chronotype, sleep apnea. Contains SNVs across genomic positions with p-values reaching genome-wide significance. |
| **API** | Check portal for API documentation |
| **License** | Academic research |
| **Size** | Hundreds to thousands of SNVs per gene; 57+ loci for insomnia |
| **Download** | GWAS summary statistics available |

### 1.5 UK Biobank Sleep Traits GWAS

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.kp4cd.org/node/235 |
| **Content** | Self-reported sleep traits GWAS from UK Biobank participants |
| **API** | Via Knowledge Portal Network |
| **License** | UK Biobank data access agreement required |
| **Size** | 453,379 participants |

---

## 2. Longevity/Aging Genetics Databases

### 2.1 Human Ageing Genomic Resources (HAGR)

| Attribute | Details |
|-----------|---------|
| **URL** | https://genomics.senescence.info/ |
| **Content** | Collection of databases: GenAge, AnAge, DrugAge, LongevityMap, CellAge, GenDR |
| **API** | No REST API; downloadable datasets |
| **License** | Freely available under HAGR license terms |
| **Format** | Tab-delimited ASCII files |

#### 2.1.1 GenAge (Ageing Gene Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://genomics.senescence.info/genes/ |
| **Content** | 307 genes linked to human ageing; 2,205 genes from model organisms (yeast, worms, flies, mice) |
| **Download** | Zipped tab-delimited ASCII datasets |
| **Key Data** | Gene name, organism, lifespan effect, observation, references |

#### 2.1.2 LongevityMap

| Attribute | Details |
|-----------|---------|
| **URL** | https://genomics.senescence.info/longevity/ |
| **Content** | 3,144 genetic variants across 884 genes associated with human longevity. Includes both positive and negative results. Details on study design, population ethnicity, sample size, age of probands. |
| **Size** | 550 entries (Build 3, 2017) |
| **Download** | Zipped tab-delimited ASCII dataset |
| **Key Genes** | FOXO3, APOE/TOMM40, CDKN2B, SH2B3/ATXN2, ABO, CETP |

#### 2.1.3 DrugAge

| Attribute | Details |
|-----------|---------|
| **URL** | https://genomics.senescence.info/drugs/ |
| **Content** | 1,316 entries featuring 418 compounds that extend lifespan in model organisms (worms, flies, yeast, mice). From 324 research articles. |
| **Download** | Available via Downloads section |
| **Key Data** | Compound name, species, lifespan effect, dosage, references |

#### 2.1.4 GenDR (Dietary Restriction Gene Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://genomics.senescence.info/diet/ |
| **Content** | 214 genes associated with dietary restriction-induced life extension |
| **Download** | Tab-delimited ASCII and Excel files |

#### 2.1.5 CellAge

| Attribute | Details |
|-----------|---------|
| **URL** | https://genomics.senescence.info/cells/ |
| **Content** | 866 genes associated with cellular senescence |
| **Download** | Tab-delimited format |

### 2.2 Open Genes Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://open-genes.com/ |
| **Content** | 2,402 genes associated with aging. Comprehensive data on lifespan-extending interventions, age-related changes, longevity associations, gene evolution, disease associations, hallmarks of aging. |
| **API** | REST API (FastAPI/Python 3): https://github.com/open-genes/open-genes-api - Returns JSON format |
| **License** | Mozilla Public License 2.0 |
| **Size** | ~80 database tables in MySQL 8 |
| **Download** | https://open-genes.com/download - TSV files and MySQL dump |
| **GitHub** | https://github.com/open-genes |

---

## 3. Nutrigenomics Databases

### 3.1 NutriGenomeDB

| Attribute | Details |
|-----------|---------|
| **URL** | https://nutrigenomedb.org/ |
| **Content** | First platform for nutrigenomics data exploration. Gene expression data from 231 manually curated nutrigenomics experiments showing effects of nutrients and bioactive food compounds on gene expression. Linked to GEO database. |
| **API** | Web services to query Panther database for molecular functions |
| **License** | Academic research |
| **Download** | Excel and PDF formats |
| **Key Features** | Gene expression browser, pathway enrichment analysis |

### 3.2 PharmGKB (Pharmacogenomics Knowledge Base)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | While focused on pharmacogenomics, contains relevant drug-gene and metabolite-gene interactions. 715 drugs, 1,761 genes, 227 diseases, 165 clinical guidelines. Some overlap with nutrient metabolism (e.g., folate pathways). |
| **API** | REST API available with registration |
| **License** | Free for research; registration required; no redistribution |
| **Download** | Zipped spreadsheets with annotations, pathways, clinical guidelines |
| **Key Data** | Drug-gene associations, variant annotations, metabolic pathways |

### 3.3 SNPedia

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.snpedia.com/ |
| **Content** | Wiki-based database of 109,729+ SNPs with medical condition associations (537+). Links to personal genomics services. Includes nutrigenomics-relevant variants (MTHFR, FTO, etc.). |
| **API** | MediaWiki API (semantic wiki) |
| **License** | Creative Commons; acquired by MyHeritage (2019) |
| **Size** | 109,729 SNPs, 537 medical conditions |
| **Tools** | Promethease for personal genome analysis |

### 3.4 FooDB (Food Compound Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://foodb.ca/ |
| **Content** | Comprehensive database of food constituents, chemistry, and biology. Food-specific compounds with nutritional information. |
| **API** | Contact required for API access |
| **License** | Free for academic use; commercial requires permission |
| **Download** | https://foodb.ca/downloads |
| **Key Data** | Food constituents, chemical structures, nutritional content |

### 3.5 Human Metabolome Database (HMDB)

| Attribute | Details |
|-----------|---------|
| **URL** | https://hmdb.ca/ |
| **Content** | World's largest organism-specific metabolomic database. 220,945 metabolite entries, 8,610 protein sequences (enzymes/transporters). Includes nutrient metabolites, enzymatic reactions, metabolic pathways. |
| **API** | Contact eponine@ualberta.ca or samackay@ualberta.ca for API access |
| **License** | Free for academic; commercial requires permission |
| **Download** | https://hmdb.ca/downloads |
| **Key Data** | Metabolite structures, concentrations, tissue locations, disease associations, genetic associations, MS/NMR spectra |

---

## 4. Vitamin Metabolism Genetics

### 4.1 Key Genes and Variants

| Vitamin | Key Genes | Notable SNPs | Database Resources |
|---------|-----------|--------------|-------------------|
| **Folate (B9)** | MTHFR, MTR, MTRR, SHMT1, FOLR3 | rs1801133 (C677T), rs1801131 (A1298C) | dbSNP, ClinVar, OMIM |
| **B12** | TCN1, TCN2, CUBN, CD320, FUT2, FUT6 | rs601338 (FUT2) | GWAS Catalog, dbSNP |
| **Vitamin D** | VDR, CYP2R1, CYP27B1, CYP24A1, GC, DHCR7 | rs11234027, rs1790349, rs12785878 | GWAS Catalog, VDR database |
| **Vitamin A** | BCMO1, RBP4 | Various | dbSNP |
| **Vitamin E** | TTPA, CYP4F2 | Various | dbSNP |

### 4.2 MTHFR Variant Resources

| Attribute | Details |
|-----------|---------|
| **OMIM Entry** | https://omim.org/entry/607093 |
| **Content** | 34 rare mutations, 9 common polymorphisms documented |
| **Key Variants** | C677T (rs1801133): 35% (het) to 70% (hom) reduced enzyme activity; A1298C (rs1801131) |
| **Population Frequency** | 25% Hispanics, 10-15% North American whites (homozygous C677T) |

### 4.3 dbNSFP (Functional Predictions Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.dbnsfp.org/ or https://sites.google.com/site/jpopgen/dbNSFP |
| **Content** | 83+ million nsSNVs and 2.4+ million splice-site SNVs. Includes metabolic pathway annotations via ConsensusPathDB and KEGG. Gene-disease relationships from GenCC, OMIM, Orphanet. |
| **Download** | Amazon, Box, Google Drive links available |
| **License** | Free for academic research |
| **Size** | Several GB |
| **Key Data** | 46 prediction algorithm scores, 9 conservation scores, pathway data |

### 4.4 Genes & Nutrition Journal Database

| Attribute | Details |
|-----------|---------|
| **Source** | https://genesandnutrition.biomedcentral.com/ |
| **Content** | Peer-reviewed research on vitamin-gene interactions |
| **Key Studies** | B12-related gene polymorphisms, folate metabolism variants |

---

## 5. Food-Gene Interaction Databases

### 5.1 GeneCards Suite

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.genecards.org/ |
| **Content** | Integrates 190+ data sources for gene-centric information. Includes metabolic pathways, compound associations, gene functions. |
| **API** | Limited batch queries via GeneALaCart (100 genes/day); full API requires academic agreement or commercial license |
| **License** | Contact LifeMap Sciences; scraping prohibited |
| **Tools** | GeneAnalytics for compound-gene relationships (83,000+ compounds) |

### 5.2 PhenoScanner

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.phenoscanner.medschl.cam.ac.uk/ |
| **Content** | 350+ million genotype-phenotype associations, 10+ million genetic variants. Includes dietary/metabolic phenotypes. Catalogs: GWAS, eQTL, pQTL, mQTL, methQTL. |
| **API** | R package: https://github.com/phenoscanner/phenoscanner; Python: https://github.com/phenoscanner/phenoscannerpy |
| **License** | Academic research |
| **Limits** | 10,000 SNPs (batches of 100), 1,000 genes (batches of 10) per hour |
| **Key Features** | Proxy SNP lookup via 1000 Genomes/HapMap |

### 5.3 KEGG (Kyoto Encyclopedia of Genes and Genomes)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.genome.jp/kegg/ |
| **Content** | Most widely used pathway database. Metabolite, reaction, enzyme, and gene information across species. Extensive nutrient metabolism pathways. |
| **API** | KEGG API: https://rest.kegg.jp/ |
| **License** | Free for academic; commercial license required for business use |
| **Key Pathways** | Vitamin digestion/absorption, one-carbon metabolism, amino acid metabolism |

### 5.4 USDA FoodData Central

| Attribute | Details |
|-----------|---------|
| **URL** | https://fdc.nal.usda.gov/ |
| **Content** | Comprehensive food and nutrient database. While not genetic, useful for linking food composition to gene-nutrient interactions. |
| **API** | REST API available for developers |
| **License** | Public domain (US government) |
| **Key Data** | Nutrient content, food composition, serving sizes |

---

## 6. Comprehensive Gene Annotation Resources

### 6.1 ClinVar

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | Aggregates genetic variants and their relationships to human health. Includes vitamin metabolism variants, drug response variants. |
| **API** | E-utilities API |
| **License** | Public domain |
| **Download** | FTP downloads available |

### 6.2 dbSNP

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/snp/ |
| **Content** | Central repository for SNPs. Contains all referenced vitamin metabolism SNPs, nutrigenomics variants. |
| **API** | E-utilities API |
| **License** | Public domain |
| **Download** | VCF, JSON formats |

### 6.3 Ensembl

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ensembl.org/ |
| **Content** | Genome browser with variant annotations, regulatory features, pathway links |
| **API** | REST API: https://rest.ensembl.org/ |
| **License** | Open access |
| **Download** | BioMart bulk download |

---

## 7. Summary Table: Database Comparison

| Database | Category | API | License | Format | Size |
|----------|----------|-----|---------|--------|------|
| GWAS Catalog | Circadian/Sleep | REST | CC0/Open | JSON/TSV | 300K+ associations |
| CircaDB | Circadian | None | Open Source | Web/Raw | ~3K genes |
| CirGRDB | Circadian | None | Academic | Text | 4.9K+ assays |
| Sleep Portal | Sleep | TBD | Academic | GWAS stats | 57+ loci |
| GenAge | Longevity | None | HAGR License | TSV | 2.5K genes |
| LongevityMap | Longevity | None | HAGR License | TSV | 3.1K variants |
| Open Genes | Longevity | REST (FastAPI) | MPL 2.0 | JSON/TSV/MySQL | 2.4K genes |
| DrugAge | Longevity | None | HAGR License | TSV | 1.3K entries |
| NutriGenomeDB | Nutrigenomics | Web Services | Academic | Excel/PDF | 231 experiments |
| PharmGKB | Nutrigenomics | REST | Research Only | Spreadsheets | 715 drugs |
| SNPedia | Nutrigenomics | MediaWiki | CC | Wiki | 109K SNPs |
| FooDB | Food-Gene | Contact | Academic | Download | Large |
| HMDB | Metabolism | Contact | Academic | XML/SDF | 220K metabolites |
| PhenoScanner | Food-Gene | R/Python | Academic | JSON/TSV | 350M associations |
| dbNSFP | Annotation | None | Academic | TSV | 83M+ variants |
| GeneCards | Annotation | Limited | License | HTML/API | 190+ sources |

---

## 8. Integration Recommendations for Gene Platform

### Priority 1: Core Databases
1. **GWAS Catalog** - Chronotype, sleep, metabolic traits (REST API)
2. **Open Genes** - Longevity genes (REST API, MPL 2.0)
3. **PhenoScanner** - Phenotype associations (R/Python API)
4. **HMDB** - Metabolite-gene links (API available)

### Priority 2: Specialized Databases
1. **LongevityMap/GenAge** - Aging variants (downloadable)
2. **CirGRDB** - Circadian regulators (downloadable)
3. **NutriGenomeDB** - Gene expression (web services)
4. **DrugAge** - Longevity compounds (downloadable)

### Priority 3: Supporting Resources
1. **dbNSFP** - Variant annotations
2. **ClinVar/dbSNP** - Clinical significance
3. **FooDB** - Food compounds
4. **KEGG/Reactome** - Metabolic pathways

### Key SNPs to Prioritize

| Category | Gene | SNP | Clinical Relevance |
|----------|------|-----|-------------------|
| Circadian | CLOCK | rs1801260 | Sleep timing |
| Circadian | PER3 | rs228697 | Chronotype |
| Longevity | FOXO3 | rs2802292 | Longevity |
| Longevity | APOE | rs429358, rs7412 | Alzheimer's/longevity |
| Folate | MTHFR | rs1801133 | Folate metabolism |
| B12 | TCN2 | rs1801198 | B12 transport |
| Vitamin D | VDR | rs2228570 | Vitamin D response |
| Vitamin D | GC | rs4588, rs7041 | Vitamin D binding |
| Diet | FTO | rs9939609 | Obesity/diet response |
| Diet | PPARG | rs1801282 | Fat metabolism |
| Lactose | LCT/MCM6 | rs4988235 | Lactose tolerance |
| Caffeine | CYP1A2 | rs762551 | Caffeine metabolism |

---

## References

1. Jones SE, et al. (2019). Genome-wide association analyses of chronotype in 697,828 individuals. Nat Commun. 10:343.
2. Budovsky A, et al. (2013). LongevityMap: a database of human genetic variants associated with longevity. Trends Genet.
3. Barardo D, et al. (2017). The DrugAge database of aging-related drugs. Aging Cell. 16(3):594-597.
4. Pizarro A, et al. (2013). CircaDB: a database of mammalian circadian gene expression profiles. Nucleic Acids Res.
5. Li X, et al. (2018). CirGRDB: a database for genome-wide deciphering circadian genes and regulators. Nucleic Acids Res.
6. Kamat MA, et al. (2019). PhenoScanner V2. Bioinformatics. 35(22):4851-4853.
7. Wishart DS, et al. (2022). HMDB 5.0: the Human Metabolome Database for 2022. Nucleic Acids Res.
8. Open Genes Consortium (2023). Open Genes - a comprehensive database of human genes associated with aging and longevity. Nucleic Acids Res.

---

*Document generated: 2026-01-19*
*For Gene Platform data integration planning*
