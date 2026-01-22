# Sleep, Longevity & Nutrigenomics Data Sources

**Document ID:** 43-79-SLEEP-LONGEVITY-NUTRI
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../index.md](./../index.md)

---

## TL;DR

Comprehensive catalog of 30+ databases covering circadian rhythm genetics (GWAS Catalog chronotype, CircaDB, CirGRDB, Sleep Disorders Knowledge Portal), longevity/aging genetics (HAGR suite including GenAge, LongevityMap, DrugAge, CellAge, GenDR; Open Genes), nutrigenomics (NutriGenomeDB, PharmGKB, FooDB, HMDB), vitamin metabolism (MTHFR/VDR variants via ClinVar, dbNSFP), and food-gene interactions (PhenoScanner, KEGG, GeneCards, USDA FoodData Central). Combined storage estimate: ~150 GB for curated data and summary statistics.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary circadian GWAS source | GWAS Catalog + Sleep Disorders Portal | 351+ chronotype loci, REST API | Jan 2026 |
| Circadian expression data | CircaDB + CirGRDB | Comprehensive time-series profiles | Jan 2026 |
| Primary longevity database | Open Genes + HAGR suite | REST API (MPL 2.0), 2,400+ genes | Jan 2026 |
| Nutrigenomics platform | NutriGenomeDB + PharmGKB | Curated experiments, pathway context | Jan 2026 |
| Metabolite data source | HMDB | 220K+ metabolites, comprehensive | Jan 2026 |
| Vitamin variant source | ClinVar + dbNSFP | Clinical significance, functional predictions | Jan 2026 |
| Food-gene interaction | PhenoScanner + KEGG | R/Python APIs, pathway context | Jan 2026 |

---

## Database Catalog

### 1. Circadian Rhythm Genetics Databases

#### 1.1 GWAS Catalog (Chronotype/Sleep Data)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **API Documentation** | https://www.ebi.ac.uk/gwas/rest/docs/api |
| **Content** | 351+ loci for chronotype ("morning person" trait), sleep duration, insomnia; data from 697,828+ UK Biobank and 23andMe participants |
| **Records** | 300,000+ variant-trait associations; chronotype-specific: 351 loci |
| **License** | Open access (CC0 for catalog, study-specific for summary statistics) |
| **API** | REST API (JSON); R package (gwasrapidd), Python package (pandasGWAS) |
| **Download** | https://www.ebi.ac.uk/gwas/docs/file-downloads |
| **Update Frequency** | Continuous curation |
| **Priority** | Tier 1 |
| **Storage Estimate** | 2 GB (catalog); 10+ GB (sleep summary stats) |

**Key Circadian Genes:**
CLOCK, BMAL1 (ARNTL), PER1-3, CRY1-2, NPAS2, CSNK1D, CSNK1E

---

#### 1.2 CircaDB (Circadian Gene Expression)

| Field | Value |
|-------|-------|
| **URL** | http://circadb.hogeneschlab.org/ |
| **GitHub** | https://github.com/itmat/circadb |
| **Content** | Database of circadian transcriptional profiles from time-course expression experiments in mice and humans |
| **Algorithms** | JTK_Cycle, Lomb Scargle, DeLichtenberg |
| **Records** | ~3,000 cycling genes with expression profiles |
| **License** | Open source |
| **API** | No REST API; web interface with Google Charts visualization |
| **Update Frequency** | Static (last major update 2014) |
| **Priority** | Tier 2 |
| **Storage Estimate** | 500 MB |

**Features:**
- Raw data downloadable
- Processed results via web interface
- Integrates with BioGPS

---

#### 1.3 CirGRDB (Circadian Genes and Regulators)

| Field | Value |
|-------|-------|
| **URL** | http://cirgrdb.biols.ac.cn |
| **Content** | Genome-wide deciphering of circadian genes and regulators; temporal expression patterns across 37 human/mouse tissues |
| **Regulatory Data** | Transcription factors, histone modifications, chromatin accessibility, enhancer RNAs, miRNAs, RNA-binding proteins, RNA editing, RNA methylation |
| **Records** | 4,936+ genome-wide assays integrated |
| **License** | Academic use |
| **API** | No dedicated API |
| **Download** | Expression and regulatory data as plain text files |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | 2 GB |

**Associations:**
- Sleep disorders
- Aging
- Tumors

---

#### 1.4 Sleep Disorders Knowledge Portal

| Field | Value |
|-------|-------|
| **URL** | https://sleep.hugeamp.org/ |
| **Content** | GWAS summary statistics for sleep phenotypes including insomnia, sleep duration, chronotype, sleep apnea |
| **Records** | Hundreds to thousands of SNVs per gene; 57+ loci for insomnia |
| **License** | Academic research |
| **API** | Check portal for API documentation |
| **Download** | GWAS summary statistics available |
| **Update Frequency** | Publication-driven |
| **Priority** | Tier 1 |
| **Storage Estimate** | 5 GB |

---

#### 1.5 UK Biobank Sleep Traits GWAS

| Field | Value |
|-------|-------|
| **URL** | https://www.kp4cd.org/node/235 |
| **Content** | Self-reported sleep traits GWAS from UK Biobank participants |
| **Records** | 453,379 participants |
| **License** | UK Biobank data access agreement required |
| **API** | Via Knowledge Portal Network |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 2 (requires application) |
| **Storage Estimate** | 2 GB |

---

### 2. Longevity/Aging Genetics Databases

#### 2.1 Human Ageing Genomic Resources (HAGR)

| Field | Value |
|-------|-------|
| **URL** | https://genomics.senescence.info/ |
| **Content** | Collection of databases: GenAge, AnAge, DrugAge, LongevityMap, CellAge, GenDR |
| **License** | Freely available under HAGR license terms |
| **API** | No REST API; downloadable datasets |
| **Format** | Tab-delimited ASCII files |
| **Priority** | Tier 1 |
| **Storage Estimate** | 500 MB (combined) |

---

##### 2.1.1 GenAge (Ageing Gene Database)

| Field | Value |
|-------|-------|
| **URL** | https://genomics.senescence.info/genes/ |
| **Content** | 307 genes linked to human ageing; 2,205 genes from model organisms (yeast, worms, flies, mice) |
| **Records** | 2,512 total genes |
| **Download** | Zipped tab-delimited ASCII datasets |
| **Key Data** | Gene name, organism, lifespan effect, observation, references |
| **Priority** | Tier 1 |

---

##### 2.1.2 LongevityMap

| Field | Value |
|-------|-------|
| **URL** | https://genomics.senescence.info/longevity/ |
| **Content** | 3,144 genetic variants across 884 genes associated with human longevity; includes positive and negative results |
| **Records** | 550 entries (Build 3, 2017) |
| **Details** | Study design, population ethnicity, sample size, age of probands |
| **Download** | Zipped tab-delimited ASCII dataset |
| **Priority** | Tier 1 |

**Key Longevity Genes:**
FOXO3, APOE/TOMM40, CDKN2B, SH2B3/ATXN2, ABO, CETP

---

##### 2.1.3 DrugAge

| Field | Value |
|-------|-------|
| **URL** | https://genomics.senescence.info/drugs/ |
| **Content** | 1,316 entries featuring 418 compounds that extend lifespan in model organisms (worms, flies, yeast, mice) |
| **Sources** | 324 research articles |
| **Download** | Available via Downloads section |
| **Key Data** | Compound name, species, lifespan effect, dosage, references |
| **Priority** | Tier 2 |

---

##### 2.1.4 GenDR (Dietary Restriction Gene Database)

| Field | Value |
|-------|-------|
| **URL** | https://genomics.senescence.info/diet/ |
| **Content** | 214 genes associated with dietary restriction-induced life extension |
| **Records** | 214 genes |
| **Download** | Tab-delimited ASCII and Excel files |
| **Priority** | Tier 2 |

---

##### 2.1.5 CellAge

| Field | Value |
|-------|-------|
| **URL** | https://genomics.senescence.info/cells/ |
| **Content** | 866 genes associated with cellular senescence |
| **Records** | 866 genes |
| **Download** | Tab-delimited format |
| **Priority** | Tier 2 |

---

#### 2.2 Open Genes Database

| Field | Value |
|-------|-------|
| **URL** | https://open-genes.com/ |
| **API** | REST API (FastAPI/Python 3): https://github.com/open-genes/open-genes-api |
| **GitHub** | https://github.com/open-genes |
| **Content** | 2,402 genes associated with aging; lifespan-extending interventions, age-related changes, longevity associations, gene evolution, disease associations, hallmarks of aging |
| **Records** | ~80 database tables in MySQL 8; 2,402 genes |
| **License** | Mozilla Public License 2.0 |
| **Download** | https://open-genes.com/download - TSV files and MySQL dump |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 |
| **Storage Estimate** | 2 GB |

---

### 3. Nutrigenomics Databases

#### 3.1 NutriGenomeDB

| Field | Value |
|-------|-------|
| **URL** | https://nutrigenomedb.org/ |
| **Content** | First platform for nutrigenomics data exploration; gene expression data from 231 manually curated nutrigenomics experiments |
| **Records** | 231 experiments linked to GEO database |
| **License** | Academic research |
| **API** | Web services to query Panther database for molecular functions |
| **Download** | Excel and PDF formats |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 |
| **Storage Estimate** | 500 MB |

**Key Features:**
- Gene expression browser
- Pathway enrichment analysis
- Effects of nutrients and bioactive food compounds on gene expression

---

#### 3.2 PharmGKB (Pharmacogenomics Knowledge Base)

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | Drug-gene and metabolite-gene interactions; 715 drugs, 1,761 genes, 227 diseases, 165 clinical guidelines |
| **Records** | 715 drugs, 1,761 genes |
| **License** | Free for research; registration required; no redistribution |
| **API** | REST API available with registration |
| **Download** | Zipped spreadsheets with annotations, pathways, clinical guidelines |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 |
| **Storage Estimate** | 500 MB |

**Nutrigenomics Overlap:**
- Folate pathways
- Vitamin metabolism variants
- Drug response annotations

---

#### 3.3 SNPedia

| Field | Value |
|-------|-------|
| **URL** | https://www.snpedia.com/ |
| **Content** | Wiki-based database of 109,729+ SNPs with medical condition associations (537+) |
| **Nutrigenomics Variants** | MTHFR, FTO, etc. |
| **Records** | 109,729 SNPs, 537 medical conditions |
| **License** | Creative Commons; acquired by MyHeritage (2019) |
| **API** | MediaWiki API (semantic wiki) |
| **Tools** | Promethease for personal genome analysis |
| **Update Frequency** | Continuous community edits |
| **Priority** | Tier 2 |
| **Storage Estimate** | 500 MB |

---

#### 3.4 FooDB (Food Compound Database)

| Field | Value |
|-------|-------|
| **URL** | https://foodb.ca/ |
| **Download** | https://foodb.ca/downloads |
| **Content** | Comprehensive database of food constituents, chemistry, and biology |
| **Records** | Food-specific compounds with nutritional information |
| **License** | Free for academic use; commercial requires permission |
| **API** | Contact required for API access |
| **Update Frequency** | Regular |
| **Priority** | Tier 2 |
| **Storage Estimate** | 5 GB |

---

#### 3.5 Human Metabolome Database (HMDB)

| Field | Value |
|-------|-------|
| **URL** | https://hmdb.ca/ |
| **Download** | https://hmdb.ca/downloads |
| **Content** | World's largest organism-specific metabolomic database; nutrient metabolites, enzymatic reactions, metabolic pathways |
| **Records** | 220,945 metabolite entries, 8,610 protein sequences (enzymes/transporters) |
| **License** | Free for academic; commercial requires permission |
| **API** | Contact eponine@ualberta.ca or samackay@ualberta.ca |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 |
| **Storage Estimate** | 20 GB |

**Key Data:**
- Metabolite structures and concentrations
- Tissue locations
- Disease associations
- Genetic associations
- MS/NMR spectra

---

### 4. Vitamin Metabolism Genetics

#### 4.1 Key Genes and Variants

| Vitamin | Key Genes | Notable SNPs | Database Resources |
|---------|-----------|--------------|-------------------|
| **Folate (B9)** | MTHFR, MTR, MTRR, SHMT1, FOLR3 | rs1801133 (C677T), rs1801131 (A1298C) | dbSNP, ClinVar, OMIM |
| **B12** | TCN1, TCN2, CUBN, CD320, FUT2, FUT6 | rs601338 (FUT2) | GWAS Catalog, dbSNP |
| **Vitamin D** | VDR, CYP2R1, CYP27B1, CYP24A1, GC, DHCR7 | rs11234027, rs1790349, rs12785878 | GWAS Catalog, VDR database |
| **Vitamin A** | BCMO1, RBP4 | Various | dbSNP |
| **Vitamin E** | TTPA, CYP4F2 | Various | dbSNP |

---

#### 4.2 MTHFR Variant Resources

| Field | Value |
|-------|-------|
| **OMIM Entry** | https://omim.org/entry/607093 |
| **Content** | 34 rare mutations, 9 common polymorphisms documented |
| **Key Variants** | C677T (rs1801133): 35% (het) to 70% (hom) reduced enzyme activity; A1298C (rs1801131) |
| **Population Frequency** | 25% Hispanics, 10-15% North American whites (homozygous C677T) |
| **Priority** | Tier 1 |

---

#### 4.3 dbNSFP (Functional Predictions Database)

| Field | Value |
|-------|-------|
| **URL** | https://www.dbnsfp.org/ or https://sites.google.com/site/jpopgen/dbNSFP |
| **Content** | 83+ million nsSNVs and 2.4+ million splice-site SNVs; metabolic pathway annotations via ConsensusPathDB and KEGG |
| **Records** | 83+ million variants |
| **License** | Free for academic research |
| **Download** | Amazon, Box, Google Drive links available |
| **Update Frequency** | Major releases |
| **Priority** | Tier 2 |
| **Storage Estimate** | 50 GB |

**Key Data:**
- 46 prediction algorithm scores
- 9 conservation scores
- Gene-disease relationships from GenCC, OMIM, Orphanet

---

#### 4.4 Genes & Nutrition Journal Database

| Field | Value |
|-------|-------|
| **Source** | https://genesandnutrition.biomedcentral.com/ |
| **Content** | Peer-reviewed research on vitamin-gene interactions |
| **Key Studies** | B12-related gene polymorphisms, folate metabolism variants |
| **Priority** | Tier 3 (literature) |

---

### 5. Food-Gene Interaction Databases

#### 5.1 GeneCards Suite

| Field | Value |
|-------|-------|
| **URL** | https://www.genecards.org/ |
| **Content** | Integrates 190+ data sources for gene-centric information; metabolic pathways, compound associations |
| **Records** | 190+ integrated sources |
| **License** | Contact LifeMap Sciences; scraping prohibited |
| **API** | Limited batch queries via GeneALaCart (100 genes/day); full API requires agreement |
| **Tools** | GeneAnalytics for compound-gene relationships (83,000+ compounds) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | N/A (web resource) |

---

#### 5.2 PhenoScanner

| Field | Value |
|-------|-------|
| **URL** | http://www.phenoscanner.medschl.cam.ac.uk/ |
| **Content** | 350+ million genotype-phenotype associations, 10+ million genetic variants |
| **Catalogs** | GWAS, eQTL, pQTL, mQTL, methQTL |
| **Records** | 350M associations, 10M variants |
| **License** | Academic research |
| **API** | R: https://github.com/phenoscanner/phenoscanner; Python: https://github.com/phenoscanner/phenoscannerpy |
| **Rate Limits** | 10,000 SNPs (batches of 100), 1,000 genes (batches of 10) per hour |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 |
| **Storage Estimate** | 50 GB |

**Features:**
- Proxy SNP lookup via 1000 Genomes/HapMap
- Dietary/metabolic phenotypes included

---

#### 5.3 KEGG (Kyoto Encyclopedia of Genes and Genomes)

| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/ |
| **API** | KEGG API: https://rest.kegg.jp/ |
| **Content** | Most widely used pathway database; metabolite, reaction, enzyme, and gene information |
| **Records** | 500+ human pathways |
| **License** | Free for academic; commercial license required |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 |
| **Storage Estimate** | 5 GB |

**Key Nutrigenomics Pathways:**
- Vitamin digestion/absorption
- One-carbon metabolism
- Amino acid metabolism

---

#### 5.4 USDA FoodData Central

| Field | Value |
|-------|-------|
| **URL** | https://fdc.nal.usda.gov/ |
| **Content** | Comprehensive food and nutrient database |
| **Records** | Extensive food composition data |
| **License** | Public domain (US government) |
| **API** | REST API available for developers |
| **Update Frequency** | Regular |
| **Priority** | Tier 2 |
| **Storage Estimate** | 2 GB |

**Key Data:**
- Nutrient content
- Food composition
- Serving sizes
- Useful for linking food composition to gene-nutrient interactions

---

### 6. Comprehensive Gene Annotation Resources

#### 6.1 ClinVar

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | Aggregates genetic variants and their relationships to human health; vitamin metabolism variants, drug response variants |
| **License** | Public domain |
| **API** | E-utilities API |
| **Download** | FTP downloads available |
| **Update Frequency** | Weekly |
| **Priority** | Tier 1 |
| **Storage Estimate** | 5 GB |

---

#### 6.2 dbSNP

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/snp/ |
| **Content** | Central repository for SNPs; all referenced vitamin metabolism SNPs, nutrigenomics variants |
| **License** | Public domain |
| **API** | E-utilities API |
| **Download** | VCF, JSON formats |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 |
| **Storage Estimate** | 50+ GB |

---

#### 6.3 Ensembl

| Field | Value |
|-------|-------|
| **URL** | https://www.ensembl.org/ |
| **API** | REST API: https://rest.ensembl.org/ |
| **Content** | Genome browser with variant annotations, regulatory features, pathway links |
| **License** | Open access |
| **Download** | BioMart bulk download |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 |
| **Storage Estimate** | 10 GB (human annotations) |

---

## Key SNPs to Prioritize

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

## Summary Comparison Table

| Database | Category | API | License | Size | Primary Use |
|----------|----------|-----|---------|------|-------------|
| **GWAS Catalog** | Circadian/Sleep | REST | CC0/Open | 300K+ assoc. | Sleep/chronotype GWAS |
| **CircaDB** | Circadian | None | Open Source | ~3K genes | Circadian expression |
| **CirGRDB** | Circadian | None | Academic | 4.9K+ assays | Circadian regulators |
| **Sleep Portal** | Sleep | TBD | Academic | 57+ loci | Sleep GWAS |
| **GenAge** | Longevity | None | HAGR License | 2.5K genes | Aging genes |
| **LongevityMap** | Longevity | None | HAGR License | 3.1K variants | Longevity variants |
| **Open Genes** | Longevity | REST (FastAPI) | MPL 2.0 | 2.4K genes | Comprehensive aging |
| **DrugAge** | Longevity | None | HAGR License | 1.3K entries | Lifespan compounds |
| **CellAge** | Longevity | None | HAGR License | 866 genes | Senescence genes |
| **GenDR** | Longevity | None | HAGR License | 214 genes | Dietary restriction |
| **NutriGenomeDB** | Nutrigenomics | Web Services | Academic | 231 experiments | Gene expression |
| **PharmGKB** | Nutrigenomics | REST | Research Only | 715 drugs | Drug-gene |
| **SNPedia** | Nutrigenomics | MediaWiki | CC | 109K SNPs | Annotated SNPs |
| **FooDB** | Food-Gene | Contact | Academic | Large | Food compounds |
| **HMDB** | Metabolism | Contact | Academic | 220K metabolites | Metabolome |
| **PhenoScanner** | Food-Gene | R/Python | Academic | 350M assoc. | Phenotype associations |
| **dbNSFP** | Annotation | None | Academic | 83M+ variants | Functional predictions |
| **GeneCards** | Annotation | Limited | License | 190+ sources | Gene information |
| **KEGG** | Pathways | REST | Academic free | 500+ pathways | Nutrient pathways |
| **ClinVar** | Clinical | E-utilities | Public domain | Weekly updates | Clinical significance |
| **dbSNP** | Reference | E-utilities | Public domain | Billions | Reference SNPs |
| **Ensembl** | Annotation | REST | Open | Comprehensive | Genome browser |

---

## Integration Priority

### Priority 1: MVP (Tier 1)

| Source | Category | Rationale |
|--------|----------|-----------|
| GWAS Catalog | Chronotype/Sleep | 351+ loci, REST API |
| Open Genes | Longevity | REST API, MPL 2.0, 2,402 genes |
| PhenoScanner | Phenotype Associations | R/Python API, 350M associations |
| HMDB | Metabolite-Gene | 220K+ metabolites, comprehensive |
| GenAge/LongevityMap | Aging Variants | Curated longevity data |
| NutriGenomeDB | Gene Expression | 231 curated experiments |
| Sleep Disorders Portal | Sleep GWAS | Specialized sleep genetics |
| PharmGKB | Drug-Gene | Clinical guidelines included |
| KEGG | Nutrient Pathways | Essential metabolic maps |
| ClinVar/dbSNP | Clinical Significance | Public domain, authoritative |
| Ensembl | Genome Annotation | REST API, open access |

### Priority 2: Post-MVP (Tier 2)

| Source | Category | Rationale |
|--------|----------|-----------|
| CircaDB | Circadian Expression | Time-series profiles |
| CirGRDB | Circadian Regulators | Multi-tissue data |
| DrugAge | Longevity Compounds | Lifespan-extending compounds |
| CellAge | Senescence Genes | Cellular aging |
| GenDR | Dietary Restriction | Diet-longevity genes |
| FooDB | Food Compounds | Food chemistry |
| SNPedia | Annotated SNPs | User-friendly annotations |
| dbNSFP | Functional Predictions | 46 prediction algorithms |
| GeneCards | Gene Information | 190+ sources integrated |
| USDA FoodData | Food Composition | Nutrient content |
| UK Biobank Sleep | Sleep Deep Phenotyping | Requires application |

### Priority 3: Future (Tier 3)

| Source | Category | Notes |
|--------|----------|-------|
| Genes & Nutrition Journal | Literature | Peer-reviewed studies |
| UK Biobank Full | Multi-phenotype | Application required |
| GeneCards Full API | Gene Annotation | Commercial agreement |

---

## API Integration Summary

**Well-documented REST APIs:**
- GWAS Catalog, Open Genes, PhenoScanner (R/Python), KEGG, ClinVar, dbSNP, Ensembl, PharmGKB

**Bulk download only:**
- GenAge, LongevityMap, DrugAge, CellAge, GenDR, CircaDB, CirGRDB, NutriGenomeDB, dbNSFP

**Requires special access:**
- HMDB (contact for API)
- FooDB (contact for API)
- GeneCards (academic agreement)
- UK Biobank (application)

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [index.md](./../index.md) | Parent navigation |
| [primary.md](./../pathways/primary.md) | Shared pathway sources |
| [pharmaceuticals.md](./../compounds/pharmaceuticals.md) | Drug/compound overlap |
| [mental-cognitive.md](./mental-cognitive.md) | Sleep/circadian overlap |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial migration from research.old/data-sources-sleep-longevity-nutri.md |
