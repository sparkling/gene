# Allergy, Histamine, and Pain Genetics Data Sources

**Research Date:** January 19, 2026
**Purpose:** Comprehensive catalog of databases for Gene Platform integration

---

## Table of Contents

1. [Allergy Genetics Databases](#1-allergy-genetics-databases)
2. [Histamine Intolerance Resources](#2-histamine-intolerance-resources)
3. [Mast Cell Activation Databases](#3-mast-cell-activation-databases)
4. [Pain Genetics Databases](#4-pain-genetics-databases)
5. [Inflammation Pathway Databases](#5-inflammation-pathway-databases)
6. [Autoimmune Genetics Databases](#6-autoimmune-genetics-databases)
7. [Multi-Purpose Genomics Platforms](#7-multi-purpose-genomics-platforms)
8. [Summary Table](#8-summary-table)

---

## 1. Allergy Genetics Databases

### 1.1 NHGRI-EBI GWAS Catalog

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Content** | Comprehensive repository of published GWAS data; >45,000 published GWAS across >5,000 human traits; >40,000 full P-value summary statistics datasets |
| **Allergy Data** | ~40 asthma GWAS, 3 atopy GWAS, 3 atopic dermatitis GWAS; 10 loci for allergic sensitization (TLR6, C11orf30, STAT6, SLC25A46, HLA-DQB1, IL1RL1, LPP, MYC, IL2, HLA-B) |
| **API** | REST API at https://www.ebi.ac.uk/gwas/rest/docs/api; Summary statistics API at https://www.ebi.ac.uk/gwas/summary-statistics/docs/ |
| **Programmatic Access** | R package: gwasrapidd (https://github.com/ramiromagno/gwasrapidd); Python package: pandasGWAS |
| **License** | Open access, FAIR compliant |
| **Data Format** | JSON, TSV downloads, REST API |
| **Size** | 1,924+ publications, 13,403+ SNPs cataloged |

**Key Allergy Findings:**
- 5 food allergy loci at genome-wide significance: SERPINB cluster (18q21.3), cytokine cluster (5q31.1), filaggrin gene, C11orf30/LRRC32, HLA region
- 76 genetic variants identified for allergic disease using age-of-onset information
- 2,038 genome-wide significant SNP loci including 31 novel loci (2025 genomic SEM study)

---

### 1.2 Allergome

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.allergome.org/ |
| **Content** | Comprehensive repository of IgE-binding compounds; WHO/IUIS-approved allergens plus non-recognized allergens; biochemical and clinical data |
| **API** | No public REST API; web interface search |
| **Programmatic Access** | Basic and advanced search via web; RefArray module for references; ReTiME for real-time IgE sensitization data |
| **License** | Free, open resource (maintained by Allergen Data Laboratories, Italy) |
| **Data Format** | Web-based, links to external databases |
| **Size** | Most comprehensive collection of allergen source and molecule data |

**Key Features:**
- Links to PDB, UniProt, NCBI databases
- Extensive literature references grouped by topic
- Support modules: RefArray, ReTiME

---

### 1.3 WHO/IUIS Allergen Nomenclature Database

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.allergen.org/ |
| **Content** | Official nomenclature for allergen proteins; systematic naming conventions |
| **API** | Limited programmatic access |
| **License** | Open access |
| **Data Format** | Structured nomenclature data |
| **Size** | All officially recognized allergens |

---

### 1.4 COMPARE (Comprehensive Protein Allergen Resource)

| Attribute | Details |
|-----------|---------|
| **URL** | https://comparedatabase.org/ |
| **Content** | Allergen sequences for food safety assessment |
| **API** | Downloadable sequence lists |
| **Programmatic Access** | Annually updated freely downloadable allergen sequences |
| **License** | Free (maintained by HESI) |
| **Data Format** | FASTA sequences |
| **Size** | Comprehensive allergen sequence database |

---

## 2. Histamine Intolerance Resources

### 2.1 Key Genes and Variants

**DAO (Diamine Oxidase) - AOC1 Gene:**

| Variant | rs ID | Effect |
|---------|-------|--------|
| C-4107T | rs2052129 | Reduced DAO activity |
| T-692G | rs2268999 | Reduced DAO activity |
| His663Asn | rs10156191 | Reduced DAO activity |
| Ser332Phe | rs1049742 | Reduced DAO activity |
| rs2071514 | rs2071514 | Moderate protective effect |
| rs1049748 | rs1049748 | Moderate protective effect |
| rs2071517 | rs2071517 | Moderate protective effect |

**HNMT (Histamine N-Methyltransferase) Gene:**

| Variant | rs ID | Effect |
|---------|-------|--------|
| Thr105Ile | rs11558538 | 30-50% reduction in HNMT activity; found in ~10% Europeans |
| C314T | - | Associated with aspirin intolerance |

### 2.2 Database Resources for Histamine Genetics

| Database | URL | Content |
|----------|-----|---------|
| **ClinVar** | https://www.ncbi.nlm.nih.gov/clinvar/ | Variant-phenotype relationships for DAO/HNMT |
| **dbSNP** | https://www.ncbi.nlm.nih.gov/snp/ | SNP records for histamine-related genes |
| **SNPedia** | https://www.snpedia.com/ | Wiki-based SNP annotations including histamine variants |
| **PharmGKB** | https://www.pharmgkb.org/ | Pharmacogenomics of antihistamines |

### 2.3 ClinVar (for Histamine Variants)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | Variant-phenotype relationships; clinical significance classifications |
| **API** | E-utilities (esearch, esummary, elink, efetch); Clinical Table Search API |
| **Programmatic Access** | Simple ClinVar (http://simple-clinvar.broadinstitute.org/) for easy queries |
| **License** | Public domain (NIH) |
| **Data Format** | XML, JSON, TSV downloads |
| **Size** | Updated weekly |
| **Query Example** | Search "histamine intolerance" or "diamine oxidase" |

---

## 3. Mast Cell Activation Databases

### 3.1 Key Genetic Associations

**Primary Genes Implicated in MCAS:**

| Gene | Function | Variants |
|------|----------|----------|
| **KIT** | Mast cell growth regulator | D816V (asp816val) - most common in mastocytosis |
| **TPSAB1** | Alpha-tryptase | Duplications cause Hereditary alpha-Tryptasemia (HaT) - 4-6% prevalence |
| **MT-CYB** | Mitochondrial cytochrome b | Associated with hEDS+MCAS |
| **HTT** | Huntingtin | Associated with hEDS+MCAS |
| **MUC3A** | Mucin | Associated with hEDS+MCAS |
| **HLA-B** | MHC Class I | Associated with hEDS+MCAS |
| **HLA-DRB1** | MHC Class II | Associated with hEDS+MCAS |

**Additional Mast Cell Regulatory Genes with Variants:**
ASXL1, CBL, DNMT3B, ETV6, EZH2, HNMT, IDH1, IL13, JAK2, KMT2A, KRAS, MS4A2, NLRP3, RASGRP4, SETBP1, SF3B1, TBXA2R, TET2, TP53

### 3.2 GARD (Genetic and Rare Diseases Information Center)

| Attribute | Details |
|-----------|---------|
| **URL** | https://rarediseases.info.nih.gov/diseases/12981/mast-cell-activation-syndrome |
| **Content** | Clinical information, genetic information, prevalence data |
| **API** | Limited; primarily web-based resource |
| **License** | Public (NIH) |
| **Data Format** | Structured clinical information |

### 3.3 OMIM (for Mastocytosis/MCAS)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.omim.org/ |
| **Content** | >27,000 gene/phenotype entries; detailed mastocytosis entries |
| **API** | OMIM API available (requires registration) |
| **Programmatic Access** | FTP download and API batch queries |
| **License** | Free for academic use; commercial license required |
| **Data Format** | JSON API output, structured text |
| **Size** | 7,400+ disorders, 4,800+ genes |

---

## 4. Pain Genetics Databases

### 4.1 Human Pain Genes Database (HPGDB)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.humanpaingenetics.ca/resources/ |
| **Content** | Comprehensive compilation of SNPs associated with pain; manually-curated literature review |
| **Structure** | Six columns: genetic locus, SNP rsID, allele/haplotype, direction of effect, associated phenotype, citation |
| **API** | Web-based search with hover details and NCBI/PubMed links |
| **License** | Open access |
| **Data Format** | Tabular format, regularly updated |
| **Size** | Comprehensive SNP-phenotype associations for pain |

### 4.2 Key Pain-Associated Genetic Variants

| Gene | Variant | rs ID | Association |
|------|---------|-------|-------------|
| **SCN9A** | Arg1150Trp | rs6746030 | Enhanced DRG excitation; increased pain in OA, sciatica, phantom limb |
| **P2RX7** | Arg270His | rs7958311 | Impaired pore formation; pain modulation |
| **CACNA2D3** | Intronic | rs6777055 | Reduced thermal pain; diminished chronic pain post-discectomy |

### 4.3 UK Biobank Pain GWAS

| Attribute | Details |
|-----------|---------|
| **Study** | Genome-wide association study of multisite chronic pain |
| **Sample Size** | ~380,000 participants |
| **Findings** | SNP heritability 10.2%; 76 independent lead SNPs at 39 risk loci |
| **Key Associations** | Brain function/development genes; mental health correlation (depression, PTSD); autoimmune traits (asthma) |
| **Data Access** | UK Biobank application required |
| **URL** | https://www.ukbiobank.ac.uk/ |

### 4.4 Pain Genetics Genomic SEM Study

| Attribute | Details |
|-----------|---------|
| **Content** | Shared genetic architecture across 24 chronic pain conditions |
| **Findings** | 184 novel targets for cross-condition chronic pain; general factor explaining genetic variance |
| **Data** | Published in peer-reviewed literature |

### 4.5 PharmGKB (for Pain Pharmacogenomics)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | Pharmacogenomics of pain medications (opioids, NSAIDs, etc.); variant-drug response annotations |
| **API** | Web services for bulk download |
| **Programmatic Access** | PharmCAT tool for clinical annotation |
| **License** | Requires account; license agreement for download |
| **Data Format** | TSV, BioPax pathways |
| **Size** | Comprehensive drug-gene-variant annotations |

---

## 5. Inflammation Pathway Databases

### 5.1 Reactome Pathway Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://reactome.org/ |
| **Content** | Curated, peer-reviewed pathways; cytokine signaling (R-HSA-1280215); inflammasomes (R-HSA-622312); TNF signaling (R-HSA-75893) |
| **API** | Graph Database, Analysis Service, Content Service API |
| **Programmatic Access** | ReactomeFIViz Cytoscape plugin; direct API access |
| **License** | Creative Commons (CC BY 4.0) |
| **Data Format** | SBML, BioPAX (Level 2/3), PDF, SVG, PNG |
| **Size** | 306 proteins in Cytokine Signaling pathway alone |
| **Download** | https://reactome.org/download-data |

**Key Inflammation Pathways:**
- Cytokine Signaling in Immune System (R-HSA-1280215)
- Inflammasomes (R-HSA-622312)
- TNF Signaling (R-HSA-75893)
- Cell Recruitment (Pro-inflammatory Response) (R-HSA-9664424)

### 5.2 KEGG Pathway Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Content** | Molecular interaction, reaction, and relation networks |
| **Inflammation Pathways** | MAPK signaling, NF-kB signaling, cytokine-cytokine receptor interaction, JAK-STAT signaling |
| **API** | KEGG REST API |
| **License** | Academic use free; commercial license required |
| **Data Format** | KGML, various image formats |

### 5.3 InnateDB

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.innatedb.com/ |
| **Content** | Innate immunity interactions and pathways; p38, ERK, JNK signaling; NF-kB transcriptional module |
| **API** | Web services available |
| **License** | Open access |
| **Data Format** | Various download formats |
| **Size** | 761 genes including TNF-alpha and IFN-alpha pathways |

### 5.4 WFINFLAM Panel

| Attribute | Details |
|-----------|---------|
| **Content** | Assembly of inflammation-related genes for pathway-focused analysis |
| **Coverage** | 1,027 candidate genes with primary subpathway assignments |
| **Pathways Included** | MAPK, NF-kB, PI3K/Akt, GPCR, cytokine, leukocyte signaling |
| **Source** | Published in PLOS ONE |

---

## 6. Autoimmune Genetics Databases

### 6.1 IPD-IMGT/HLA Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/ipd/imgt/hla/ |
| **Content** | Official sequences of human MHC; WHO Nomenclature Committee named alleles |
| **API** | Alignment tool; allele query tool |
| **License** | Open access (EBI) |
| **Data Format** | Sequence alignments, FASTA |
| **Size** | All officially named HLA alleles |

**HLA-Disease Associations:**
- Associated with 100+ diseases
- Accounts for ~50% of genetic susceptibility to Type 1 diabetes
- Key in: rheumatoid arthritis, psoriasis, asthma, autoimmune disorders

### 6.2 ImmunoBase

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.immunobase.org/ |
| **Content** | Curated GWAS data for immunologically-related diseases; ImmunoChip consortium data for 12 diseases |
| **API** | Summary statistics download |
| **License** | Open access |
| **Data Format** | Summary statistics TSV |
| **Size** | Example: 133,352 SNPs for celiac disease |

**Diseases Covered:**
- Celiac disease
- Rheumatoid arthritis
- Systemic lupus erythematosus
- Multiple sclerosis
- Primary biliary cirrhosis
- Ulcerative colitis
- Crohn's disease

### 6.3 GAAD (Gene and Autoimmune Disease Association Database)

| Attribute | Details |
|-----------|---------|
| **URL** | Published in Genomics, Proteomics & Bioinformatics |
| **Content** | Gene-autoimmune disease associations from public databases and MEDLINE |
| **Size** | 44,762 associations between 49 autoimmune diseases and 4,249 genes |
| **License** | Open access |

### 6.4 ADEx (Autoimmune Diseases Explorer)

| Attribute | Details |
|-----------|---------|
| **URL** | https://adex.genyo.es |
| **Content** | Integrated transcriptomics and methylation studies for autoimmune diseases |
| **GitHub** | https://github.com/GENyO-BioInformatics/ADEx_public |
| **License** | Open access |
| **Size** | 82 curated studies, 5,609 samples |

---

## 7. Multi-Purpose Genomics Platforms

### 7.1 Open Targets Genetics

| Attribute | Details |
|-----------|---------|
| **URL** | https://genetics.opentargets.org/ (now integrated into Platform) |
| **Content** | GWAS and functional genomics integration; gene expression, protein abundance, chromatin data |
| **API** | GraphQL API; BigQuery access |
| **Programmatic Access** | R package: otargen; FTP and Google Cloud downloads |
| **License** | Open access |
| **Data Format** | JSON, Parquet, TSV |
| **Documentation** | https://genetics-docs.opentargets.org/ |

**Relevant for:**
- Autoimmune and inflammatory diseases
- IBD, ulcerative colitis, Crohn's disease
- Cross-disease colocalization analysis

### 7.2 DisGeNET

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.disgenet.org/ |
| **Content** | Disease-gene associations; integrates multiple sources including literature |
| **API** | REST API for programmatic access |
| **Programmatic Access** | disgenet2r R package; Cytoscape App |
| **License** | Free for academic (Academic license); commercial license required for for-profit |
| **Data Format** | TSV, API responses |
| **Size** | 24,000+ diseases/traits, 17,000 genes, 117,000 genomic variants |

### 7.3 dbSNP

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/snp/ |
| **Content** | SNP and short variation database |
| **API** | E-utilities; NIH Clinical Table Search API |
| **Programmatic Access** | R packages: rsnps, biomaRt; Ensembl REST API |
| **License** | Public domain (NCBI) |
| **Data Format** | JSON, XML, VCF |
| **FTP** | ftp://ftp.ncbi.nih.gov/snp/ |

### 7.4 SNPedia

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.snpedia.com/ |
| **Content** | Wiki-based SNP database with phenotype annotations |
| **API** | MediaWiki API (500 hits/query max) |
| **Programmatic Access** | SNPediaR (Bioconductor); Python wikitools |
| **License** | Creative Commons |
| **Size** | 109,729 SNPs, 537 medical conditions |

### 7.5 European Variation Archive (EVA)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/eva/ |
| **Content** | All types of genetic variation from all species |
| **API** | REST API, Variant Browser |
| **License** | Open access (EBI) |
| **Data Format** | VCF, JSON |
| **Note** | Continuation of dbSNP identifiers (SS/RS) |

---

## 8. Summary Table

| Database | Focus Area | API | License | Data Size |
|----------|-----------|-----|---------|-----------|
| **GWAS Catalog** | Allergy GWAS | REST API | Open | 45,000+ studies |
| **Allergome** | Allergen molecules | Web only | Open | Comprehensive |
| **COMPARE** | Allergen sequences | Download | Free | Updated annually |
| **ClinVar** | Clinical variants | E-utilities | Public domain | Weekly updates |
| **HPGDB** | Pain genetics | Web | Open | Comprehensive |
| **PharmGKB** | Pharmacogenomics | Web services | License req. | Comprehensive |
| **Reactome** | Inflammation pathways | GraphQL | CC BY 4.0 | 2,600+ pathways |
| **KEGG** | Metabolic pathways | REST API | Academic free | 500+ pathways |
| **InnateDB** | Innate immunity | Web services | Open | 761+ genes |
| **IPD-IMGT/HLA** | HLA alleles | Query tools | Open | All HLA alleles |
| **ImmunoBase** | Autoimmune GWAS | Download | Open | 12 diseases |
| **GAAD** | Autoimmune genes | - | Open | 44,762 associations |
| **ADEx** | Autoimmune omics | GitHub | Open | 5,609 samples |
| **Open Targets** | Drug targets | GraphQL | Open | Comprehensive |
| **DisGeNET** | Gene-disease | REST API | Academic free | 24,000+ diseases |
| **dbSNP** | All SNPs | E-utilities | Public domain | Billions of variants |
| **SNPedia** | Annotated SNPs | MediaWiki | CC | 109,729 SNPs |
| **OMIM** | Mendelian diseases | OMIM API | Academic free | 27,000+ entries |
| **EVA** | Variation archive | REST API | Open | All species |

---

## Integration Recommendations for Gene Platform

### Priority 1 (Core Integration)
1. **GWAS Catalog** - Primary source for allergy/inflammation GWAS SNPs
2. **ClinVar** - Clinical significance for all variants
3. **dbSNP** - Reference SNP data
4. **Open Targets** - Target validation and drug discovery

### Priority 2 (Specialized Data)
1. **HPGDB** - Pain-specific genetic associations
2. **PharmGKB** - Drug response variants
3. **Reactome** - Pathway context
4. **IPD-IMGT/HLA** - HLA typing for autoimmune risk

### Priority 3 (Enhanced Coverage)
1. **DisGeNET** - Broader disease-gene associations
2. **ImmunoBase** - Autoimmune-specific GWAS
3. **Allergome** - Allergen-specific data
4. **SNPedia** - User-friendly annotations

### API Integration Notes

- **REST APIs**: GWAS Catalog, ClinVar, dbSNP, Open Targets, DisGeNET, EVA
- **GraphQL**: Open Targets Platform
- **Download-only**: ImmunoBase, COMPARE, HPGDB
- **License Required**: PharmGKB (academic), DisGeNET (commercial), OMIM (API access)

---

## References

1. GWAS Catalog: https://www.ebi.ac.uk/gwas/
2. Allergome: https://www.allergome.org/
3. ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
4. Human Pain Genes Database: https://www.humanpaingenetics.ca/
5. PharmGKB: https://www.pharmgkb.org/
6. Reactome: https://reactome.org/
7. KEGG: https://www.genome.jp/kegg/
8. InnateDB: https://www.innatedb.com/
9. IPD-IMGT/HLA: https://www.ebi.ac.uk/ipd/imgt/hla/
10. ImmunoBase: https://www.immunobase.org/
11. Open Targets: https://platform.opentargets.org/
12. DisGeNET: https://www.disgenet.org/
13. dbSNP: https://www.ncbi.nlm.nih.gov/snp/
14. SNPedia: https://www.snpedia.com/
15. OMIM: https://www.omim.org/
16. EVA: https://www.ebi.ac.uk/eva/
