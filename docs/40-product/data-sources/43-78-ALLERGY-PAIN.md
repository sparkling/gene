# Allergy, Histamine & Pain Genetics Databases

**Document ID:** 43-78-ALLERGY-PAIN
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-DATA-SOURCES.md](../43-DATA-SOURCES.md)

---

## TL;DR

Comprehensive inventory of 25+ databases covering allergy genetics, histamine intolerance, mast cell activation, pain genetics, and inflammation pathways. Primary sources include GWAS Catalog (45K+ studies), HPGDB (pain SNP associations), Reactome (inflammation pathways), and IPD-IMGT/HLA (43K+ alleles). Key allergy findings: 5 food allergy loci at genome-wide significance, 76 variants identified using age-of-onset. Pain: 76 independent lead SNPs at 39 risk loci from UK Biobank (380K participants). Estimated total storage: ~120 GB with majority from GWAS summary statistics.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary allergy GWAS source | NHGRI-EBI GWAS Catalog | 45K+ studies, open access, REST API, R/Python packages | Jan 2026 |
| Pain genetics primary | HPGDB + UK Biobank | Curated SNP-phenotype associations, largest cohort (380K) | Jan 2026 |
| Inflammation pathways | Reactome | CC BY 4.0, GraphQL API, peer-reviewed curation | Jan 2026 |
| Allergen nomenclature | WHO/IUIS + Allergome | Official nomenclature + comprehensive molecule data | Jan 2026 |
| Histamine variant source | ClinVar + dbSNP | Clinical significance + reference SNP data | Jan 2026 |
| HLA associations | IPD-IMGT/HLA | WHO official, 43K+ alleles, CC BY-ND | Jan 2026 |

---

## Database Catalog

### Category 1: Allergy Genetics

#### NHGRI-EBI GWAS Catalog (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Content** | Comprehensive repository of published GWAS data |
| **Records** | 45,000+ studies, 5,000+ traits, 40,000+ summary statistics datasets |
| **Allergy Coverage** | ~40 asthma GWAS, 3 atopy GWAS, 3 atopic dermatitis GWAS |
| **Key Loci** | TLR6, C11orf30, STAT6, SLC25A46, HLA-DQB1, IL1RL1, LPP, MYC, IL2, HLA-B |
| **License** | Open access, FAIR compliant |
| **API** | REST API (https://www.ebi.ac.uk/gwas/rest/docs/api); Summary statistics API |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 GB (with summary statistics) |
| **R Package** | gwasrapidd (https://github.com/ramiromagno/gwasrapidd) |
| **Python Package** | pandasGWAS |

**Key Allergy Findings:**
- 5 food allergy loci at genome-wide significance: SERPINB (18q21.3), cytokine cluster (5q31.1), filaggrin gene, C11orf30/LRRC32, HLA region
- 76 genetic variants identified for allergic disease using age-of-onset information
- 2,038 genome-wide significant SNP loci including 31 novel loci (2025 genomic SEM study)

---

#### Allergome

| Field | Value |
|-------|-------|
| **URL** | https://www.allergome.org/ |
| **Content** | Comprehensive repository of IgE-binding compounds; WHO/IUIS-approved allergens plus non-recognized allergens |
| **Records** | Most comprehensive collection of allergen source and molecule data |
| **Features** | Links to PDB, UniProt, NCBI; RefArray module for references; ReTiME for real-time IgE sensitization |
| **License** | Free, open resource (Allergen Data Laboratories, Italy) |
| **API** | None (web interface search) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |

---

#### WHO/IUIS Allergen Nomenclature Database

| Field | Value |
|-------|-------|
| **URL** | http://www.allergen.org/ |
| **Content** | Official nomenclature for allergen proteins; systematic naming conventions |
| **Records** | All officially recognized allergens |
| **License** | Open access |
| **API** | Limited programmatic access |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

---

#### COMPARE (Comprehensive Protein Allergen Resource)

| Field | Value |
|-------|-------|
| **URL** | https://comparedatabase.org/ |
| **Content** | Allergen sequences for food safety assessment |
| **Records** | Comprehensive allergen sequence database |
| **License** | Free (maintained by HESI) |
| **API** | Downloadable sequence lists |
| **Update Frequency** | Annual |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB |
| **Format** | FASTA sequences |

---

### Category 2: Histamine & Mast Cell Genetics

#### Key Histamine Intolerance Variants

**DAO (Diamine Oxidase) - AOC1 Gene:**

| Variant | rs ID | Effect | Population |
|---------|-------|--------|------------|
| C-4107T | rs2052129 | Reduced DAO activity | European |
| T-692G | rs2268999 | Reduced DAO activity | European |
| His663Asn | rs10156191 | Reduced DAO activity | European |
| Ser332Phe | rs1049742 | Reduced DAO activity | European |
| - | rs2071514 | Moderate protective effect | European |
| - | rs1049748 | Moderate protective effect | European |
| - | rs2071517 | Moderate protective effect | European |

**HNMT (Histamine N-Methyltransferase) Gene:**

| Variant | rs ID | Effect | Prevalence |
|---------|-------|--------|------------|
| Thr105Ile | rs11558538 | 30-50% reduction in HNMT activity | ~10% Europeans |
| C314T | - | Associated with aspirin intolerance | - |

---

#### ClinVar (for Histamine Variants)

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | Variant-phenotype relationships; clinical significance classifications |
| **Records** | 2M+ submissions, weekly updates |
| **Query** | Search "histamine intolerance" or "diamine oxidase" |
| **License** | Public domain (NIH) |
| **API** | E-utilities (esearch, esummary, elink, efetch); Clinical Table Search API |
| **Update Frequency** | Weekly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~5 GB |
| **Alt Access** | Simple ClinVar (http://simple-clinvar.broadinstitute.org/) |

---

#### Mast Cell Activation - Key Genetic Associations

**Primary MCAS Genes:**

| Gene | Function | Key Variants |
|------|----------|--------------|
| **KIT** | Mast cell growth regulator | D816V (asp816val) - most common in mastocytosis |
| **TPSAB1** | Alpha-tryptase | Duplications cause Hereditary alpha-Tryptasemia (HaT); 4-6% prevalence |
| **MT-CYB** | Mitochondrial cytochrome b | Associated with hEDS+MCAS |
| **HTT** | Huntingtin | Associated with hEDS+MCAS |
| **MUC3A** | Mucin | Associated with hEDS+MCAS |
| **HLA-B** | MHC Class I | Associated with hEDS+MCAS |
| **HLA-DRB1** | MHC Class II | Associated with hEDS+MCAS |

**Additional Mast Cell Regulatory Genes:**
ASXL1, CBL, DNMT3B, ETV6, EZH2, HNMT, IDH1, IL13, JAK2, KMT2A, KRAS, MS4A2, NLRP3, RASGRP4, SETBP1, SF3B1, TBXA2R, TET2, TP53

---

#### GARD (Genetic and Rare Diseases Information Center)

| Field | Value |
|-------|-------|
| **URL** | https://rarediseases.info.nih.gov/diseases/12981/mast-cell-activation-syndrome |
| **Content** | Clinical information, genetic information, prevalence data for MCAS |
| **Records** | Comprehensive rare disease entries |
| **License** | Public (NIH) |
| **API** | Limited; primarily web-based resource |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB |

---

#### OMIM (for Mastocytosis/MCAS)

| Field | Value |
|-------|-------|
| **URL** | https://www.omim.org/ |
| **Content** | 27,000+ gene/phenotype entries; detailed mastocytosis entries |
| **Records** | 7,400+ disorders, 4,800+ genes |
| **License** | Free for academic use; commercial license required |
| **API** | OMIM API (requires registration) |
| **Update Frequency** | Daily |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |
| **Format** | JSON API output, structured text |

---

### Category 3: Pain Genetics

#### Human Pain Genes Database (HPGDB) (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://www.humanpaingenetics.ca/resources/ |
| **Content** | Comprehensive compilation of SNPs associated with pain; manually-curated literature review |
| **Records** | Comprehensive SNP-phenotype associations for pain |
| **Structure** | Six columns: genetic locus, SNP rsID, allele/haplotype, direction of effect, associated phenotype, citation |
| **License** | Open access |
| **API** | Web-based search with hover details and NCBI/PubMed links |
| **Update Frequency** | Regularly updated |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 MB |

---

#### UK Biobank Pain GWAS (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://www.ukbiobank.ac.uk/ |
| **Content** | Genome-wide association study of multisite chronic pain |
| **Records** | ~380,000 participants |
| **Findings** | SNP heritability 10.2%; 76 independent lead SNPs at 39 risk loci |
| **Key Associations** | Brain function/development genes; mental health correlation (depression, PTSD); autoimmune traits (asthma) |
| **License** | Application required |
| **API** | Data access via approved applications |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~10 GB (summary statistics) |

---

#### Key Pain-Associated Genetic Variants

| Gene | Variant | rs ID | Association |
|------|---------|-------|-------------|
| **SCN9A** | Arg1150Trp | rs6746030 | Enhanced DRG excitation; increased pain in OA, sciatica, phantom limb |
| **P2RX7** | Arg270His | rs7958311 | Impaired pore formation; pain modulation |
| **CACNA2D3** | Intronic | rs6777055 | Reduced thermal pain; diminished chronic pain post-discectomy |

---

#### Pain Genetics Genomic SEM Study

| Field | Value |
|-------|-------|
| **Content** | Shared genetic architecture across 24 chronic pain conditions |
| **Findings** | 184 novel targets for cross-condition chronic pain; general factor explaining genetic variance |
| **Data** | Published in peer-reviewed literature |
| **Priority** | Tier 2 |

---

#### PharmGKB (for Pain Pharmacogenomics)

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | Pharmacogenomics of pain medications (opioids, NSAIDs, etc.); variant-drug response annotations |
| **Records** | 700+ drugs, 1000+ genes |
| **License** | Requires account; license agreement for download |
| **API** | Web services for bulk download |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~1 GB |
| **Tool** | PharmCAT for clinical annotation |
| **Format** | TSV, BioPax pathways |

---

### Category 4: Inflammation Pathways

#### Reactome Pathway Database (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://reactome.org/ |
| **Content** | Curated, peer-reviewed inflammation pathways |
| **Key Pathways** | Cytokine Signaling (R-HSA-1280215), Inflammasomes (R-HSA-622312), TNF Signaling (R-HSA-75893), Cell Recruitment (R-HSA-9664424) |
| **Records** | 306 proteins in Cytokine Signaling pathway alone |
| **License** | Creative Commons (CC BY 4.0) |
| **API** | Graph Database, Analysis Service, Content Service API |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB |
| **Download** | https://reactome.org/download-data |
| **Formats** | SBML, BioPAX (Level 2/3), PDF, SVG, PNG |
| **Plugin** | ReactomeFIViz Cytoscape plugin |

---

#### KEGG Pathway Database

| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Content** | Molecular interaction, reaction, and relation networks |
| **Inflammation Pathways** | MAPK signaling, NF-kB signaling, cytokine-cytokine receptor interaction, JAK-STAT signaling |
| **License** | Academic use free; commercial license required |
| **API** | KEGG REST API |
| **Update Frequency** | Monthly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |
| **Format** | KGML, various image formats |

---

#### InnateDB

| Field | Value |
|-------|-------|
| **URL** | https://www.innatedb.com/ |
| **Content** | Innate immunity interactions and pathways; p38, ERK, JNK signaling; NF-kB transcriptional module |
| **Records** | 761 genes including TNF-alpha and IFN-alpha pathways |
| **License** | Open access |
| **API** | Web services available |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~200 MB |

---

#### WFINFLAM Panel

| Field | Value |
|-------|-------|
| **Content** | Assembly of inflammation-related genes for pathway-focused analysis |
| **Records** | 1,027 candidate genes with primary subpathway assignments |
| **Pathways Included** | MAPK, NF-kB, PI3K/Akt, GPCR, cytokine, leukocyte signaling |
| **Source** | Published in PLOS ONE |
| **License** | Open access |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~50 MB |

---

### Category 5: Autoimmune & HLA Associations

#### IPD-IMGT/HLA Database (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/ipd/imgt/hla/ |
| **Content** | Official sequences of human MHC; WHO Nomenclature Committee named alleles |
| **Records** | 43,000+ unique alleles from 47 genes |
| **License** | Creative Commons Attribution-NoDerivs (CC BY-ND) |
| **API** | Alignment tool; allele query tool |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB |
| **Formats** | FASTA, MSF, PIR, XML |

**HLA-Disease Associations:**
- Associated with 100+ diseases
- Accounts for ~50% of genetic susceptibility to Type 1 diabetes
- Key in: rheumatoid arthritis, psoriasis, asthma, autoimmune disorders

---

#### ImmunoBase

| Field | Value |
|-------|-------|
| **URL** | https://www.immunobase.org/ |
| **Content** | Curated GWAS data for immunologically-related diseases; ImmunoChip consortium data |
| **Records** | 12 diseases; Example: 133,352 SNPs for celiac disease |
| **Diseases** | Celiac disease, RA, SLE, MS, primary biliary cirrhosis, UC, Crohn's disease |
| **License** | Open access |
| **API** | Summary statistics download |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 GB |

---

#### GAAD (Gene and Autoimmune Disease Association Database)

| Field | Value |
|-------|-------|
| **URL** | Published in Genomics, Proteomics & Bioinformatics |
| **Content** | Gene-autoimmune disease associations from public databases and MEDLINE |
| **Records** | 44,762 associations between 49 autoimmune diseases and 4,249 genes |
| **License** | Open access |
| **API** | None |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

---

#### ADEx (Autoimmune Diseases Explorer)

| Field | Value |
|-------|-------|
| **URL** | https://adex.genyo.es |
| **Content** | Integrated transcriptomics and methylation studies for autoimmune diseases |
| **Records** | 82 curated studies, 5,609 samples |
| **License** | Open access |
| **API** | None |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~5 GB |
| **GitHub** | https://github.com/GENyO-BioInformatics/ADEx_public |

---

### Category 6: Multi-Purpose Genomics Platforms

#### Open Targets Genetics

| Field | Value |
|-------|-------|
| **URL** | https://genetics.opentargets.org/ (integrated into Platform) |
| **Content** | GWAS and functional genomics integration; gene expression, protein abundance, chromatin data |
| **Use Cases** | Autoimmune and inflammatory diseases, IBD, UC, Crohn's, cross-disease colocalization |
| **License** | Open access |
| **API** | GraphQL API; BigQuery access |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~20 GB |
| **R Package** | otargen |
| **Formats** | JSON, Parquet, TSV |

---

#### DisGeNET

| Field | Value |
|-------|-------|
| **URL** | https://www.disgenet.org/ |
| **Content** | Disease-gene associations; integrates multiple sources including literature |
| **Records** | 24,000+ diseases/traits, 17,000 genes, 117,000 genomic variants |
| **License** | Free for academic; commercial license required |
| **API** | REST API for programmatic access |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 GB |
| **R Package** | disgenet2r |
| **Plugin** | Cytoscape App |

---

#### dbSNP

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/snp/ |
| **Content** | SNP and short variation database |
| **Records** | 1+ billion SNPs |
| **License** | Public domain (NCBI) |
| **API** | E-utilities; NIH Clinical Table Search API |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 GB |
| **R Packages** | rsnps, biomaRt |
| **FTP** | ftp://ftp.ncbi.nih.gov/snp/ |
| **Formats** | JSON, XML, VCF |

---

#### SNPedia

| Field | Value |
|-------|-------|
| **URL** | https://www.snpedia.com/ |
| **Content** | Wiki-based SNP database with phenotype annotations |
| **Records** | 109,729 SNPs, 537 medical conditions |
| **License** | Creative Commons |
| **API** | MediaWiki API (500 hits/query max) |
| **Update Frequency** | Community-driven |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |
| **R Package** | SNPediaR (Bioconductor) |
| **Python** | wikitools |

---

#### European Variation Archive (EVA)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/eva/ |
| **Content** | All types of genetic variation from all species |
| **Records** | Continuation of dbSNP identifiers (SS/RS) |
| **License** | Open access (EBI) |
| **API** | REST API, Variant Browser |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 GB |
| **Formats** | VCF, JSON |

---

## Storage Summary

| Category | Sources | Est. Size |
|----------|---------|-----------|
| Allergy Genetics | 4 | ~55 GB |
| Histamine & Mast Cell | 4 | ~6 GB |
| Pain Genetics | 5 | ~15 GB |
| Inflammation Pathways | 4 | ~3 GB |
| Autoimmune & HLA | 4 | ~12 GB |
| Multi-Purpose Platforms | 5 | ~175 GB |
| **Total** | **26** | **~120 GB** |

*Note: dbSNP full data is ~100 GB; using summary/subset data reduces to ~20 GB for allergy/pain-specific variants*

---

## Summary Table

| Database | URL | Content | Records | License | API | Update | Priority | Storage |
|----------|-----|---------|---------|---------|-----|--------|----------|---------|
| **GWAS Catalog** | ebi.ac.uk/gwas | Allergy GWAS | 45K+ studies | Open | REST | Continuous | Tier 1 | 50 GB |
| **Allergome** | allergome.org | Allergen molecules | Comprehensive | Open | Web only | Continuous | Tier 2 | 500 MB |
| **WHO/IUIS** | allergen.org | Allergen nomenclature | All official | Open | Limited | Periodic | Tier 2 | 50 MB |
| **COMPARE** | comparedatabase.org | Allergen sequences | Comprehensive | Free | Download | Annual | Tier 2 | 100 MB |
| **ClinVar** | ncbi.nlm.nih.gov/clinvar | Clinical variants | 2M+ | Public domain | E-utilities | Weekly | Tier 1 | 5 GB |
| **GARD** | rarediseases.info.nih.gov | MCAS clinical info | Comprehensive | Public | Limited | Continuous | Tier 2 | 100 MB |
| **OMIM** | omim.org | Mendelian diseases | 27K+ entries | Academic free | OMIM API | Daily | Tier 2 | 500 MB |
| **HPGDB** | humanpaingenetics.ca | Pain genetics | Comprehensive | Open | Web | Regular | Tier 1 | 100 MB |
| **UK Biobank** | ukbiobank.ac.uk | Pain GWAS | 380K participants | Application | Limited | Continuous | Tier 1 | 10 GB |
| **PharmGKB** | pharmgkb.org | Pharmacogenomics | 700+ drugs | License req. | Web services | Monthly | Tier 1 | 1 GB |
| **Reactome** | reactome.org | Inflammation pathways | 2,600+ pathways | CC BY 4.0 | GraphQL | Quarterly | Tier 1 | 2 GB |
| **KEGG** | genome.jp/kegg | Metabolic pathways | 500+ pathways | Academic free | REST | Monthly | Tier 2 | 500 MB |
| **InnateDB** | innatedb.com | Innate immunity | 761+ genes | Open | Web services | Periodic | Tier 2 | 200 MB |
| **IPD-IMGT/HLA** | ebi.ac.uk/ipd/imgt/hla | HLA alleles | 43K+ alleles | CC BY-ND | Query tools | Quarterly | Tier 1 | 2 GB |
| **ImmunoBase** | immunobase.org | Autoimmune GWAS | 12 diseases | Open | Download | Periodic | Tier 2 | 5 GB |
| **GAAD** | - | Autoimmune genes | 44,762 assoc. | Open | None | Static | Tier 2 | 50 MB |
| **ADEx** | adex.genyo.es | Autoimmune omics | 5,609 samples | Open | None | Periodic | Tier 3 | 5 GB |
| **Open Targets** | opentargets.org | Drug targets | Comprehensive | Open | GraphQL | Quarterly | Tier 1 | 20 GB |
| **DisGeNET** | disgenet.org | Gene-disease | 24K+ diseases | Academic free | REST | Periodic | Tier 2 | 2 GB |
| **dbSNP** | ncbi.nlm.nih.gov/snp | All SNPs | 1B+ variants | Public domain | E-utilities | Continuous | Tier 1 | 100 GB |
| **SNPedia** | snpedia.com | Annotated SNPs | 109K SNPs | CC | MediaWiki | Community | Tier 2 | 500 MB |
| **EVA** | ebi.ac.uk/eva | Variation archive | All species | Open | REST | Continuous | Tier 2 | 50 GB |

---

## Integration Recommendations

### Tier 1 - Essential (High Value, Open Access)

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | GWAS Catalog | Primary source for allergy/inflammation GWAS SNPs, open access, REST API |
| 2 | ClinVar | Clinical significance for all variants, public domain |
| 3 | dbSNP | Reference SNP data, public domain |
| 4 | Open Targets | Target validation and drug discovery, GraphQL API |
| 5 | HPGDB | Pain-specific genetic associations, curated |
| 6 | Reactome | Pathway context for inflammation, CC BY 4.0 |
| 7 | IPD-IMGT/HLA | HLA typing for autoimmune risk, WHO official |

### Tier 2 - Important

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | PharmGKB | Drug response variants for pain medications |
| 2 | DisGeNET | Broader disease-gene associations |
| 3 | ImmunoBase | Autoimmune-specific GWAS |
| 4 | Allergome | Allergen-specific data |
| 5 | SNPedia | User-friendly annotations |
| 6 | OMIM | Mendelian disease context |

### Tier 3 - Supplementary

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | KEGG | Pathway maps (academic use) |
| 2 | InnateDB | Innate immunity detail |
| 3 | ADEx | Autoimmune expression data |
| 4 | GAAD | Gene-autoimmune associations |
| 5 | WFINFLAM | Inflammation gene panel |

---

## API Integration Notes

### REST APIs Available

| Database | API Endpoint |
|----------|--------------|
| GWAS Catalog | https://www.ebi.ac.uk/gwas/rest/docs/api |
| ClinVar | E-utilities (https://www.ncbi.nlm.nih.gov/books/NBK25501/) |
| dbSNP | E-utilities |
| Open Targets | GraphQL (https://api.platform.opentargets.org/api/v4/graphql) |
| DisGeNET | https://www.disgenet.org/api/ |
| EVA | https://www.ebi.ac.uk/eva/webservices/rest/ |
| Reactome | https://reactome.org/ContentService/ |

### R Packages

| Package | Database |
|---------|----------|
| gwasrapidd | GWAS Catalog |
| clinvarR | ClinVar |
| rsnps | dbSNP |
| otargen | Open Targets |
| disgenet2r | DisGeNET |
| SNPediaR | SNPedia |
| ReactomePA | Reactome pathway analysis |

### Python Packages

| Package | Database |
|---------|----------|
| pandasGWAS | GWAS Catalog |
| pyensembl | Ensembl/EVA |
| opentargets | Open Targets Platform |

### Download-Only Resources

| Database | Access Method |
|----------|---------------|
| ImmunoBase | Summary statistics download |
| COMPARE | Annual FASTA download |
| HPGDB | Web interface with links |
| UK Biobank | Application required |
| WFINFLAM | Publication supplementary |

---

## License Summary

### Fully Open Access

| Database | License | Download Method |
|----------|---------|-----------------|
| GWAS Catalog | Open | FTP, API |
| ClinVar | Public domain | FTP, API |
| dbSNP | Public domain | FTP, API |
| Reactome | CC BY 4.0 | FTP, API |
| Open Targets | Open | FTP, API, BigQuery |
| IPD-IMGT/HLA | CC BY-ND | FTP, GitHub |
| InnateDB | Open | Web services |
| EVA | Open | API, FTP |
| HPGDB | Open | Web |
| Allergome | Open | Web |
| WHO/IUIS | Open | Web |
| COMPARE | Free | Download |
| ImmunoBase | Open | Download |
| GAAD | Open | Web |
| ADEx | Open | GitHub |
| SNPedia | CC | MediaWiki API |

### Academic Use Only

| Database | License | Notes |
|----------|---------|-------|
| KEGG | Academic free | Commercial license required for business use |
| OMIM | Academic free | API key required, yearly renewal |
| PharmGKB | License required | Account needed, agreement for download |
| DisGeNET | Academic free | Commercial license for for-profit use |
| UK Biobank | Application | Research project approval required |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-00-INDEX.md](./43-00-INDEX.md) | Parent index |
| [43-74-AUTOIMMUNE.md](./43-74-AUTOIMMUNE.md) | Overlaps with HLA, autoimmune GWAS |
| [43-41-PATHWAYS-PRIMARY.md](./43-41-PATHWAYS-PRIMARY.md) | Overlaps with Reactome |
| [43-51-PHARMACEUTICALS.md](./43-51-PHARMACEUTICALS.md) | Overlaps with PharmGKB |
| [44-ARCHITECTURE.md](../44-ARCHITECTURE.md) | Informs database design |
| [45-DATA-MODEL.md](../45-DATA-MODEL.md) | Informs entity structure |

---

## Open Questions

- [ ] Histamine variant panel - which SNPs to include in MVP (DAO, HNMT)?
- [ ] UK Biobank application - timeline for pain GWAS access?
- [ ] MCAS genetic panel - standardize variant list across sources?
- [ ] Allergen cross-reactivity mapping - best approach for integration?
- [ ] Pain pharmacogenomics - opioid response variants for initial panel?

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial migration from research.old/data-sources-allergy-pain.md |
