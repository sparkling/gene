# Hormone and Autoimmune Disease Databases

**Last Updated:** January 2026
**Purpose:** Comprehensive inventory of hormone, endocrine, and autoimmune disease genetics databases for the Gene Platform.

---

## Summary

| Category | Sources | Key Resources | Est. Size |
|----------|---------|---------------|-----------|
| Thyroid Genetics | 3 | ThyroidOmics, GWAS Catalog, TCGA | ~5 GB |
| Hormone Pathways | 4 | KEGG, Reactome, Hmrbase2, HormoneBase | ~10 GB |
| Autoimmune Genetics | 5 | GWAS Catalog, OpenGWAS, GAAD, ADEx, Open Targets | ~50 GB |
| HLA/Immune System | 3 | IPD-IMGT/HLA, IPD-KIR, gnomAD | ~15 GB |
| Endocrine Disorders | 3 | EndoGene, OMIM, MalaCards | ~5 GB |
| Steroid Metabolism | 3 | Reactome, KEGG, PharmVar | ~3 GB |

---

## 1. Thyroid Genetics Databases

### 1.1 ThyroidOmics Consortium - PRIMARY

| Field | Value |
|-------|-------|
| **URL** | https://transfer.sysepi.medizin.uni-greifswald.de/thyroidomics/datasets/ |
| **Content** | GWAS summary statistics for thyroid function traits |
| **Traits** | TSH, FT4, FT3, TT3, FT3/FT4 ratio, TPOAb, goiter, thyroid volume, Hashimoto's thyroiditis |
| **Sample Size** | Up to 271,040 euthyroid individuals from 46 cohorts |
| **API** | None (direct download) |
| **Formats** | Summary statistics files (TSV/TXT) |
| **License** | Academic use (check individual datasets) |
| **Alt Access** | CHARGE dbGaP (phs000930) |
| **Size** | ~2 GB |
| **Note** | Largest thyroid-specific GWAS consortium |

### 1.2 GWAS Catalog - Thyroid Traits

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/efotraits/EFO_0006812 |
| **Content** | Curated thyroid disease associations |
| **Traits** | Autoimmune thyroid disease, hypothyroidism, hyperthyroidism, thyroid cancer |
| **Associations** | 260+ genetic associations for TSH alone |
| **API** | REST API (https://www.ebi.ac.uk/gwas/rest/docs/api) |
| **R Package** | gwasrapidd |
| **Formats** | TSV, JSON |
| **License** | Open access |
| **Size** | ~100 MB (thyroid subset) |

### 1.3 TCGA Thyroid Cancer (THCA)

| Field | Value |
|-------|-------|
| **URL** | https://portal.gdc.cancer.gov/ |
| **Content** | Comprehensive genomic characterization of thyroid cancer |
| **Samples** | 507 papillary thyroid carcinoma cases |
| **Key Findings** | BRAFV600E, RAS mutations, RET/NTRK/ALK fusions |
| **Mutation Rate** | 97% of papillary thyroid cancers have driver mutations |
| **API** | GDC API |
| **Formats** | MAF, VCF, TSV |
| **License** | Open access (controlled access for patient data) |
| **Size** | ~500 GB (full THCA project) |

---

## 2. Hormone Pathway Databases

### 2.1 Reactome - Steroid/Hormone Pathways - PRIMARY

| Field | Value |
|-------|-------|
| **URL** | https://reactome.org/content/detail/R-HSA-196071 |
| **Content** | Curated hormone metabolism pathways |
| **Pathways** | Steroid hormone metabolism, thyroid hormone synthesis, insulin signaling |
| **Components** | Pregnenolone, glucocorticoid, mineralocorticoid, androgen, estrogen biosynthesis |
| **API** | RESTful API, Content Service endpoints |
| **Formats** | SBML, BioPAX (L2/L3), PDF, SVG, PNG, PPTX, SBGN |
| **License** | **Creative Commons Attribution 4.0 International (CC BY 4.0)** |
| **Size** | ~2 GB (full database) |
| **Note** | Excellent for pathway visualization |

### 2.2 KEGG - Endocrine System Pathways

| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Content** | Manually curated pathway maps |
| **Endocrine Pathways** | Thyroid hormone synthesis (04918), Thyroid signaling (04919), Insulin secretion/signaling, GnRH, Estrogen, Prolactin, Oxytocin, Growth hormone, Parathyroid, Renin-angiotensin |
| **Steroid Pathway** | Steroid hormone biosynthesis (map00140) |
| **API** | REST API (academic use only) |
| **License** | **Academic: Free (FTP subscription available); Commercial: License required** |
| **Commercial Contact** | https://www.pathway.jp/en/licensing.html |
| **Size** | ~500 MB (pathway data) |
| **Note** | Non-academic use requires commercial license |

### 2.3 Hmrbase2 - Hormone Receptor Database

| Field | Value |
|-------|-------|
| **URL** | https://webs.iiitd.edu.in/raghava/hmrbase2/ |
| **Legacy URL** | http://crdd.osdd.net/raghava/hmrbase/ |
| **Content** | Hormones, receptors, and hormone-receptor pairs |
| **Entries** | 12,056 total (7,406 peptide hormones, 753 non-peptide hormones, 3,897 receptors) |
| **Pairs** | 5,662 hormone-receptor pairs |
| **Organisms** | 803 species |
| **API** | None (web interface) |
| **Search** | Keyword, organism, MW, amino acids, PDB/KEGG/DrugBank IDs |
| **License** | Free for academic use |
| **Size** | ~100 MB |
| **Backend** | MySQL, PHP, Apache |

### 2.4 HormoneBase - Steroid Levels Database

| Field | Value |
|-------|-------|
| **URL** | https://www.nature.com/articles/sdata201897 |
| **Content** | Population-level steroid hormone data across vertebrates |
| **Entries** | 6,580+ entries from 476 species |
| **Hormones** | Glucocorticoids (baseline/stress-induced), androgens |
| **Sources** | 648 publications (1967-2015) + unpublished data |
| **Data Type** | Mean, variation, range in free-living adult vertebrates |
| **API** | None |
| **Formats** | CSV/Excel (supplementary data) |
| **License** | Open access (Scientific Data) |
| **Size** | ~10 MB |

---

## 3. Autoimmune Disease Genetics Databases

### 3.1 NHGRI-EBI GWAS Catalog - PRIMARY

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Content** | Curated GWAS associations and summary statistics |
| **Statistics** | 7,400+ publications, 1,040,000+ SNP-trait associations |
| **Autoimmune Traits** | RA, SLE, MS, T1D, IBD, celiac disease, psoriasis, etc. |
| **API** | REST API (https://www.ebi.ac.uk/gwas/rest/docs/api) |
| **R Package** | gwasrapidd |
| **Formats** | TSV, JSON, Summary statistics |
| **License** | Open access |
| **FTP** | ftp://ftp.ebi.ac.uk/pub/databases/gwas/ |
| **Size** | ~50 GB (with summary statistics) |

### 3.2 OpenGWAS (MRC IEU) - PRIMARY

| Field | Value |
|-------|-------|
| **URL** | https://gwas.mrcieu.ac.uk/ |
| **API URL** | https://api.opengwas.io/api/ |
| **Content** | Complete GWAS summary datasets |
| **Statistics** | 126 billion genetic associations from 14,582 datasets |
| **API** | REST API with authentication |
| **R Package** | ieugwasr |
| **Python Package** | ieugwaspy |
| **Formats** | VCF, TSV |
| **License** | Open access |
| **Size** | ~100 GB+ |
| **Note** | Query limits apply; register for increased access |

### 3.3 GAAD - Gene and Autoimmune Disease Association Database

| Field | Value |
|-------|-------|
| **URL** | http://gaad.medgenius.info |
| **Genes List** | http://gaad.medgenius.info/genes |
| **Content** | Gene-autoimmune disease associations |
| **Associations** | 44,762 associations between 49 ADs and 4,249 genes |
| **Sources** | 19,299 MEDLINE documents + NCBI Gene + GeneCards |
| **Features** | Co-occurring gene pairs, polymorphism data, expression changes |
| **API** | None (web interface) |
| **License** | Free access |
| **Size** | ~50 MB |
| **Published** | Genomics, Proteomics & Bioinformatics (2018) |

### 3.4 ADEx - Autoimmune Diseases Explorer

| Field | Value |
|-------|-------|
| **URL** | https://adex.genyo.es |
| **Content** | Integrated transcriptomics and methylation data |
| **Studies** | 82 curated studies, 5,609 samples |
| **Diseases** | SLE, Rheumatoid Arthritis, Sjogren's Syndrome, Systemic Sclerosis, Type 1 Diabetes |
| **Analysis** | Differential expression, pathway analysis, meta-analysis (Rank Product) |
| **API** | None |
| **Download** | Expression/methylation data as text files |
| **License** | Academic use (contact for details) |
| **Size** | ~5 GB |

### 3.5 Open Targets Platform

| Field | Value |
|-------|-------|
| **URL** | https://platform.opentargets.org/ |
| **Genetics URL** | https://genetics.opentargets.org/ (redirects to Platform) |
| **Content** | Drug target identification with genetic evidence |
| **Features** | Target-disease associations, GWAS loci to gene mapping |
| **Autoimmune Data** | IBD, RA, SLE, MS, and other immune diseases integrated |
| **API** | GraphQL API |
| **Download** | EMBL-EBI FTP, Google BigQuery, Google Cloud Storage |
| **License** | Open access |
| **Size** | ~20 GB |

### 3.6 ImmunoBase (Historical)

| Field | Value |
|-------|-------|
| **URL** | https://www.immunobase.org/ (may be deprecated) |
| **Content** | ImmunoChip consortium GWAS data for 12 immune diseases |
| **Diseases** | RA, SLE, MS, T1D, UC, celiac disease, primary biliary cirrhosis |
| **Status** | **Data migrated to GWAS Catalog** |
| **Note** | Previously supported by Eli Lilly; use GWAS Catalog for current access |

---

## 4. HLA/Immune System Databases

### 4.1 IPD-IMGT/HLA Database - PRIMARY

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/ipd/imgt/hla/ |
| **Content** | Official WHO nomenclature HLA sequences |
| **Alleles** | 43,000+ unique alleles from 47 genes |
| **Sequences** | 77,000+ component entries |
| **Release** | 3.63 (January 2026) |
| **Growth** | 8,148 new alleles since 2023 |
| **API** | Allele Query Tool, Alignment Tool |
| **Download** | FTP site, GitHub (ANHIG/IMGTHLA) |
| **Formats** | FASTA, MSF, PIR, XML |
| **License** | **Creative Commons Attribution-NoDerivs** |
| **Size** | ~2 GB |
| **Clinical Use** | Transplant matching, disease association |

### 4.2 IPD-KIR Database

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/ipd/kir/ |
| **GitHub** | https://github.com/ANHIG/IPDKIR |
| **Content** | Killer Immunoglobulin-like Receptor sequences |
| **Features** | Allele alignments, haplotypes, ethnicity data, primer search |
| **Files** | KIR.dat (EMBL format), KIR_gen.fasta, KIR_nuc.fasta, KIR_prot.fasta |
| **API** | Download via FTP/GitHub |
| **License** | **Creative Commons Attribution-NoDerivs** |
| **Size** | ~100 MB |

### 4.3 gnomAD - Immune Variant Frequencies

| Field | Value |
|-------|-------|
| **URL** | https://gnomad.broadinstitute.org/ |
| **Content** | Allele frequencies including HLA/immune regions |
| **Version** | v4.1 (734,947 exomes + 76,215 genomes) |
| **Features** | Filtering allele frequency (FAF), local ancestry inference |
| **API** | GraphQL API |
| **Download** | AWS Open Data (s3://gnomad-public-us-east-1/) |
| **License** | Open access |
| **Size** | ~30 TB (full), ~50 GB (summary) |
| **Note** | Use AWS to avoid egress fees |

---

## 5. Endocrine Disorder Databases

### 5.1 EndoGene Database

| Field | Value |
|-------|-------|
| **URL** | https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2025.1472754/full |
| **Content** | Genetic variants in endocrine disorder patients |
| **Patients** | 5,926 Russian patients |
| **Diseases** | 450 endocrine and concomitant diseases |
| **Variants** | 2,711 clinically relevant variants |
| **Panels** | 4 panels (220-382 genes) + WES (31,969 genes) |
| **Period** | November 2017 - January 2024 |
| **API** | None |
| **License** | Open access (Frontiers publication) |
| **Size** | ~50 MB |
| **Published** | February 2025 |

### 5.2 OMIM - Endocrine Entries

| Field | Value |
|-------|-------|
| **URL** | https://www.omim.org/ |
| **Content** | Mendelian disorders including endocrine conditions |
| **Entries** | 15,000+ genes, extensive endocrine disease entries |
| **Endocrine Examples** | Diabetes, thyroid disorders, adrenal insufficiency, pheochromocytoma |
| **API** | OMIM API (requires registration, yearly renewal) |
| **Download** | FTP (requires agreement) |
| **License** | Free for academic; commercial license required |
| **Updates** | Daily |
| **Size** | ~500 MB |

### 5.3 MalaCards - Endocrine Diseases

| Field | Value |
|-------|-------|
| **URL** | https://www.malacards.org/ |
| **Content** | Integrated human disease database |
| **Total Disorders** | 22,560 disorders |
| **Endocrine Diseases** | 1,665 endocrine diseases |
| **Genetic Diseases** | 8,439 genetic diseases |
| **Sources** | 78 integrated web sources |
| **API** | Limited (web interface primary) |
| **License** | Free for academic use |
| **Size** | ~1 GB |

### 5.4 NIH Genetic Testing Registry (GTR) - Endocrine

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/gtr/conditions/C0014130/ |
| **Content** | Genetic tests for endocrine disorders |
| **Key Genes** | CYP11B1, CYP17A1, CYP21A2, HSD3B2, STAR |
| **Disorders** | Neoplastic and non-neoplastic endocrine conditions |
| **API** | E-utilities API |
| **License** | Public domain |

---

## 6. Steroid Metabolism Databases

### 6.1 Reactome - Steroid Metabolism - PRIMARY

| Field | Value |
|-------|-------|
| **URL** | https://reactome.org/content/detail/R-HSA-8957322 |
| **Pathway** | R-HSA-196071 (Metabolism of steroid hormones) |
| **Content** | Complete steroid biosynthesis and metabolism |
| **Enzymes** | CYP11A1, CYP17A1, CYP19A1, and related enzymes |
| **Steps** | Cholesterol to pregnenolone, glucocorticoids, mineralocorticoids, androgens, estrogens |
| **API** | RESTful API |
| **Formats** | SBML, BioPAX, SBGN |
| **License** | **CC BY 4.0** |
| **Reference** | Payne & Hales (2004) Endocrine Reviews |

### 6.2 KEGG - Steroid Hormone Biosynthesis

| Field | Value |
|-------|-------|
| **URL** | https://www.kegg.jp/pathway/map00140 |
| **Human Pathway** | https://www.kegg.jp/pathway/hsa00140 |
| **Content** | Steroid hormone biosynthesis pathway map |
| **Compounds** | C27 cholesterol to C21, C19, C18 steroids |
| **Key Enzymes** | CYP11A1 (cholesterol side-chain cleavage), CYP17, CYP19 |
| **API** | KEGG REST API (academic use) |
| **License** | Academic free; commercial license required |

### 6.3 PharmVar - Pharmacogene Variation

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmvar.org/ |
| **Content** | Pharmacogene allele nomenclature |
| **Genes** | 29+ CYP enzymes + POR (steroid-related: CYP17A1, CYP19A1, CYP21A2) |
| **Steroid CYPs** | CYP11A1, CYP11B1, CYP11B2, CYP17A1, CYP19A1, CYP21A2 |
| **API** | Data files available for download |
| **Formats** | VCF, FASTA, standardized definitions |
| **License** | Open access |
| **Partners** | PharmGKB, CPIC |
| **Size** | ~200 MB |
| **Note** | Essential for steroid pharmacogenomics |

### 6.4 Harmonizome - Steroid Gene Sets

| Field | Value |
|-------|-------|
| **URL** | https://maayanlab.cloud/Harmonizome/gene_set/c21+steroid+hormone+metabolism/KEGG+Pathways |
| **Content** | Gene sets for steroid metabolism pathways |
| **C21 Pathway** | 11 proteins in c21 steroid hormone metabolism |
| **API** | REST API |
| **License** | Open access |

---

## 7. Cross-Disease & Integration Resources

### 7.1 Endo-ERN (European Reference Network)

| Field | Value |
|-------|-------|
| **URL** | https://endo-ern.eu/ |
| **Content** | Rare endocrine conditions guidelines and resources |
| **Focus** | NGS implementation, genetic testing standards |
| **Value** | Clinical interpretation guidelines for endocrine variants |

### 7.2 Gene Ontology - Steroid/Hormone Terms

| Field | Value |
|-------|-------|
| **URL** | https://www.informatics.jax.org/vocab/gene_ontology/GO:0008202 |
| **Term** | GO:0008202 (steroid metabolic process) |
| **Content** | Standardized gene function annotations |
| **API** | AmiGO API |
| **License** | CC BY 4.0 |

---

## Data Access Summary

### Fully Open Access (No Restrictions)

| Database | License | Download |
|----------|---------|----------|
| Reactome | CC BY 4.0 | FTP, API |
| GWAS Catalog | Open | FTP, API |
| OpenGWAS | Open | API, Files |
| IPD-IMGT/HLA | CC BY-ND | FTP, GitHub |
| gnomAD | Open | AWS S3 |
| GAAD | Free | Web |
| PharmVar | Open | Files |

### Academic Use Only

| Database | License | Notes |
|----------|---------|-------|
| KEGG | Academic free | Commercial license required for business |
| Hmrbase2 | Academic | Web interface |
| OMIM | Academic free | API key required |
| MalaCards | Academic | Web interface |
| GeneCards | Academic | JSON dumps available |

### Deprecated/Limited

| Database | Status | Alternative |
|----------|--------|-------------|
| ImmunoBase | Deprecated | GWAS Catalog |
| T1DBase | Deprecated | GWAS Catalog, Open Targets |

---

## Recommended Integration Priority

### Tier 1 - Essential (High Value, Open Access)

1. **IPD-IMGT/HLA** - Gold standard for HLA nomenclature
2. **GWAS Catalog** - Comprehensive autoimmune GWAS data
3. **Reactome** - Hormone pathway data with excellent API
4. **ThyroidOmics** - Best thyroid-specific GWAS resource
5. **PharmVar** - Steroid pharmacogenomics

### Tier 2 - Important

1. **OpenGWAS** - Bulk GWAS data access
2. **GAAD** - Autoimmune gene associations
3. **Hmrbase2** - Hormone-receptor relationships
4. **Open Targets** - Drug target integration

### Tier 3 - Supplementary

1. **ADEx** - Autoimmune expression data
2. **EndoGene** - Endocrine variant collection
3. **KEGG** - Pathway maps (academic use)
4. **MalaCards** - Disease integration

---

## API Integration Notes

### REST APIs Available

```
GWAS Catalog: https://www.ebi.ac.uk/gwas/rest/docs/api
OpenGWAS: https://api.opengwas.io/api/
Reactome: https://reactome.org/ContentService/
gnomAD: GraphQL at https://gnomad.broadinstitute.org/api
Open Targets: GraphQL at https://api.platform.opentargets.org/api/v4/graphql
```

### R Packages

- `gwasrapidd` - GWAS Catalog
- `ieugwasr` - OpenGWAS
- `ReactomePA` - Reactome pathway analysis

### Python Packages

- `ieugwaspy` - OpenGWAS
- `opentargets` - Open Targets Platform

---

## Sources

- [ThyroidOmics Consortium](https://transfer.sysepi.medizin.uni-greifswald.de/thyroidomics/datasets/)
- [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/)
- [OpenGWAS](https://gwas.mrcieu.ac.uk/)
- [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/)
- [Reactome Pathway Database](https://reactome.org/)
- [KEGG Pathway Database](https://www.kegg.jp/kegg/pathway.html)
- [Hmrbase2](https://webs.iiitd.edu.in/raghava/hmrbase2/)
- [GAAD Database](http://gaad.medgenius.info)
- [ADEx Database](https://adex.genyo.es)
- [Open Targets Platform](https://platform.opentargets.org/)
- [PharmVar Consortium](https://www.pharmvar.org/)
- [EndoGene Database](https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2025.1472754/full)
- [OMIM](https://www.omim.org/)
- [MalaCards](https://www.malacards.org/)
- [gnomAD](https://gnomad.broadinstitute.org/)
