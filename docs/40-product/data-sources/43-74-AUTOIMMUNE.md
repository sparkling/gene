# Autoimmune & Hormone Disease Databases

**Document ID:** 43-74-AUTOIMMUNE
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-00-INDEX.md](./43-00-INDEX.md)

---

## TL;DR

Comprehensive inventory of 23 databases across 7 categories covering autoimmune diseases, hormone pathways, thyroid genetics, HLA/immune system, endocrine disorders, and steroid metabolism. Primary sources include IPD-IMGT/HLA (43K+ alleles), GWAS Catalog (1M+ associations), Reactome (hormone pathways), ThyroidOmics (271K individuals), and PharmVar (steroid CYP enzymes). Estimated total storage: ~88 GB with majority from GWAS summary statistics. Key APIs: GWAS Catalog REST, OpenGWAS REST, Reactome REST, gnomAD GraphQL, Open Targets GraphQL.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary HLA source | IPD-IMGT/HLA | WHO official nomenclature, 43K+ alleles, CC BY-ND license | Jan 2026 |
| Autoimmune GWAS primary | GWAS Catalog + OpenGWAS | Open access, comprehensive coverage, REST APIs | Jan 2026 |
| Hormone pathway source | Reactome | CC BY 4.0 license, excellent API, visualization | Jan 2026 |
| Thyroid genetics primary | ThyroidOmics | Largest consortium (271K individuals), 46 cohorts | Jan 2026 |
| Steroid pharmacogenomics | PharmVar | Open access, CYP enzyme nomenclature, VCF format | Jan 2026 |
| ImmunoBase status | Deprecated | Data migrated to GWAS Catalog | Jan 2026 |

---

## Database Catalog

### Category 1: Thyroid Genetics

#### ThyroidOmics Consortium (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://transfer.sysepi.medizin.uni-greifswald.de/thyroidomics/datasets/ |
| **Content** | GWAS summary statistics for thyroid function traits |
| **Records** | 271,040 euthyroid individuals from 46 cohorts |
| **Traits** | TSH, FT4, FT3, TT3, FT3/FT4 ratio, TPOAb, goiter, thyroid volume, Hashimoto's thyroiditis |
| **License** | Academic use (check individual datasets) |
| **API** | None (direct download) |
| **Update Frequency** | Periodic consortium releases |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB |
| **Alt Access** | CHARGE dbGaP (phs000930) |

#### GWAS Catalog - Thyroid Traits

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/efotraits/EFO_0006812 |
| **Content** | Curated thyroid disease associations |
| **Records** | 260+ genetic associations for TSH; autoimmune thyroid, hypothyroidism, hyperthyroidism, thyroid cancer |
| **License** | Open access |
| **API** | REST API (https://www.ebi.ac.uk/gwas/rest/docs/api) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 MB (thyroid subset) |
| **R Package** | gwasrapidd |

#### TCGA Thyroid Cancer (THCA)

| Field | Value |
|-------|-------|
| **URL** | https://portal.gdc.cancer.gov/ |
| **Content** | Comprehensive genomic characterization of thyroid cancer |
| **Records** | 507 papillary thyroid carcinoma cases |
| **Key Findings** | BRAFV600E, RAS mutations, RET/NTRK/ALK fusions; 97% driver mutation rate |
| **License** | Open access (controlled access for patient data) |
| **API** | GDC API |
| **Update Frequency** | Project-based |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 GB (full THCA project); ~5 GB (summary data) |

---

### Category 2: Hormone Pathways

#### Reactome - Steroid/Hormone Pathways (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://reactome.org/content/detail/R-HSA-196071 |
| **Content** | Curated hormone metabolism pathways |
| **Records** | Steroid hormone metabolism, thyroid hormone synthesis, insulin signaling pathways |
| **Components** | Pregnenolone, glucocorticoid, mineralocorticoid, androgen, estrogen biosynthesis |
| **License** | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| **API** | RESTful API, Content Service endpoints |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB (full database) |
| **Formats** | SBML, BioPAX (L2/L3), PDF, SVG, PNG, PPTX, SBGN |

#### KEGG - Endocrine System Pathways

| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Content** | Manually curated pathway maps |
| **Records** | Thyroid hormone synthesis (04918), Thyroid signaling (04919), Insulin secretion/signaling, GnRH, Estrogen, Prolactin, Oxytocin, Growth hormone, Parathyroid, Renin-angiotensin, Steroid biosynthesis (map00140) |
| **License** | Academic: Free (FTP subscription available); Commercial: License required |
| **API** | REST API (academic use only) |
| **Update Frequency** | Monthly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB (pathway data) |
| **Commercial Contact** | https://www.pathway.jp/en/licensing.html |

#### Hmrbase2 - Hormone Receptor Database

| Field | Value |
|-------|-------|
| **URL** | https://webs.iiitd.edu.in/raghava/hmrbase2/ |
| **Content** | Hormones, receptors, and hormone-receptor pairs |
| **Records** | 12,056 total (7,406 peptide hormones, 753 non-peptide hormones, 3,897 receptors, 5,662 pairs) |
| **Organisms** | 803 species |
| **License** | Free for academic use |
| **API** | None (web interface) |
| **Update Frequency** | Static |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB |
| **Search** | Keyword, organism, MW, amino acids, PDB/KEGG/DrugBank IDs |

#### HormoneBase - Steroid Levels Database

| Field | Value |
|-------|-------|
| **URL** | https://www.nature.com/articles/sdata201897 |
| **Content** | Population-level steroid hormone data across vertebrates |
| **Records** | 6,580+ entries from 476 species; 648 publications (1967-2015) |
| **Hormones** | Glucocorticoids (baseline/stress-induced), androgens |
| **License** | Open access (Scientific Data) |
| **API** | None |
| **Update Frequency** | Static (2018 publication) |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~10 MB |
| **Formats** | CSV/Excel (supplementary data) |

---

### Category 3: Autoimmune Disease Genetics

#### NHGRI-EBI GWAS Catalog (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Content** | Curated GWAS associations and summary statistics |
| **Records** | 7,400+ publications, 1,040,000+ SNP-trait associations |
| **Autoimmune Coverage** | RA, SLE, MS, T1D, IBD, celiac disease, psoriasis, and more |
| **License** | Open access |
| **API** | REST API (https://www.ebi.ac.uk/gwas/rest/docs/api) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 GB (with summary statistics) |
| **FTP** | ftp://ftp.ebi.ac.uk/pub/databases/gwas/ |
| **R Package** | gwasrapidd |

#### OpenGWAS (MRC IEU) (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://gwas.mrcieu.ac.uk/ |
| **API URL** | https://api.opengwas.io/api/ |
| **Content** | Complete GWAS summary datasets |
| **Records** | 126 billion genetic associations from 14,582 datasets |
| **License** | Open access |
| **API** | REST API with authentication |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 GB+ |
| **R Package** | ieugwasr |
| **Python Package** | ieugwaspy |
| **Note** | Query limits apply; register for increased access |

#### GAAD - Gene and Autoimmune Disease Association Database

| Field | Value |
|-------|-------|
| **URL** | http://gaad.medgenius.info |
| **Genes List** | http://gaad.medgenius.info/genes |
| **Content** | Gene-autoimmune disease associations |
| **Records** | 44,762 associations between 49 autoimmune diseases and 4,249 genes |
| **Sources** | 19,299 MEDLINE documents + NCBI Gene + GeneCards |
| **License** | Free access |
| **API** | None (web interface) |
| **Update Frequency** | Static (2018) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |
| **Features** | Co-occurring gene pairs, polymorphism data, expression changes |

#### ADEx - Autoimmune Diseases Explorer

| Field | Value |
|-------|-------|
| **URL** | https://adex.genyo.es |
| **Content** | Integrated transcriptomics and methylation data |
| **Records** | 82 curated studies, 5,609 samples |
| **Diseases** | SLE, Rheumatoid Arthritis, Sjogren's Syndrome, Systemic Sclerosis, Type 1 Diabetes |
| **License** | Academic use (contact for details) |
| **API** | None |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 GB |
| **Analysis** | Differential expression, pathway analysis, meta-analysis (Rank Product) |

#### Open Targets Platform

| Field | Value |
|-------|-------|
| **URL** | https://platform.opentargets.org/ |
| **Content** | Drug target identification with genetic evidence |
| **Records** | IBD, RA, SLE, MS, and other immune diseases integrated |
| **Features** | Target-disease associations, GWAS loci to gene mapping |
| **License** | Open access |
| **API** | GraphQL API |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~20 GB |
| **Download** | EMBL-EBI FTP, Google BigQuery, Google Cloud Storage |

#### ImmunoBase (Historical - DEPRECATED)

| Field | Value |
|-------|-------|
| **URL** | https://www.immunobase.org/ (deprecated) |
| **Content** | ImmunoChip consortium GWAS data for 12 immune diseases |
| **Diseases** | RA, SLE, MS, T1D, UC, celiac disease, primary biliary cirrhosis |
| **Status** | Data migrated to GWAS Catalog |
| **Priority** | N/A - Use GWAS Catalog |
| **Note** | Previously supported by Eli Lilly |

---

### Category 4: HLA/Immune System

#### IPD-IMGT/HLA Database (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/ipd/imgt/hla/ |
| **Content** | Official WHO nomenclature HLA sequences |
| **Records** | 43,000+ unique alleles from 47 genes; 77,000+ component entries |
| **Release** | 3.63 (January 2026); 8,148 new alleles since 2023 |
| **License** | Creative Commons Attribution-NoDerivs (CC BY-ND) |
| **API** | Allele Query Tool, Alignment Tool |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB |
| **Download** | FTP site, GitHub (ANHIG/IMGTHLA) |
| **Formats** | FASTA, MSF, PIR, XML |
| **Clinical Use** | Transplant matching, disease association |

#### IPD-KIR Database

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/ipd/kir/ |
| **GitHub** | https://github.com/ANHIG/IPDKIR |
| **Content** | Killer Immunoglobulin-like Receptor sequences |
| **Records** | Allele alignments, haplotypes, ethnicity data |
| **License** | Creative Commons Attribution-NoDerivs (CC BY-ND) |
| **API** | Download via FTP/GitHub |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB |
| **Files** | KIR.dat (EMBL format), KIR_gen.fasta, KIR_nuc.fasta, KIR_prot.fasta |

#### gnomAD - Immune Variant Frequencies

| Field | Value |
|-------|-------|
| **URL** | https://gnomad.broadinstitute.org/ |
| **Content** | Allele frequencies including HLA/immune regions |
| **Records** | v4.1: 734,947 exomes + 76,215 genomes |
| **Features** | Filtering allele frequency (FAF), local ancestry inference |
| **License** | Open access |
| **API** | GraphQL API |
| **Update Frequency** | Major version releases |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~30 TB (full); ~50 GB (summary) |
| **Download** | AWS Open Data (s3://gnomad-public-us-east-1/) |
| **Note** | Use AWS to avoid egress fees |

---

### Category 5: Endocrine Disorders

#### EndoGene Database

| Field | Value |
|-------|-------|
| **URL** | https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2025.1472754/full |
| **Content** | Genetic variants in endocrine disorder patients |
| **Records** | 5,926 Russian patients, 450 diseases, 2,711 clinically relevant variants |
| **Panels** | 4 panels (220-382 genes) + WES (31,969 genes) |
| **License** | Open access (Frontiers publication) |
| **API** | None |
| **Update Frequency** | Static (February 2025) |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~50 MB |
| **Period** | November 2017 - January 2024 |

#### OMIM - Endocrine Entries

| Field | Value |
|-------|-------|
| **URL** | https://www.omim.org/ |
| **Content** | Mendelian disorders including endocrine conditions |
| **Records** | 15,000+ genes; extensive endocrine disease entries |
| **Examples** | Diabetes, thyroid disorders, adrenal insufficiency, pheochromocytoma |
| **License** | Free for academic; commercial license required |
| **API** | OMIM API (requires registration, yearly renewal) |
| **Update Frequency** | Daily |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |
| **Download** | FTP (requires agreement) |

#### MalaCards - Endocrine Diseases

| Field | Value |
|-------|-------|
| **URL** | https://www.malacards.org/ |
| **Content** | Integrated human disease database |
| **Records** | 22,560 total disorders; 1,665 endocrine diseases; 8,439 genetic diseases |
| **Sources** | 78 integrated web sources |
| **License** | Free for academic use |
| **API** | Limited (web interface primary) |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~1 GB |

#### NIH Genetic Testing Registry (GTR) - Endocrine

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/gtr/conditions/C0014130/ |
| **Content** | Genetic tests for endocrine disorders |
| **Key Genes** | CYP11B1, CYP17A1, CYP21A2, HSD3B2, STAR |
| **License** | Public domain |
| **API** | E-utilities API |
| **Update Frequency** | Continuous |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~100 MB |

---

### Category 6: Steroid Metabolism

#### Reactome - Steroid Metabolism (PRIMARY)

| Field | Value |
|-------|-------|
| **URL** | https://reactome.org/content/detail/R-HSA-8957322 |
| **Pathway** | R-HSA-196071 (Metabolism of steroid hormones) |
| **Content** | Complete steroid biosynthesis and metabolism |
| **Enzymes** | CYP11A1, CYP17A1, CYP19A1, and related enzymes |
| **Steps** | Cholesterol to pregnenolone, glucocorticoids, mineralocorticoids, androgens, estrogens |
| **License** | CC BY 4.0 |
| **API** | RESTful API |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | Included in Reactome total (~2 GB) |
| **Formats** | SBML, BioPAX, SBGN |

#### KEGG - Steroid Hormone Biosynthesis

| Field | Value |
|-------|-------|
| **URL** | https://www.kegg.jp/pathway/map00140 |
| **Human Pathway** | https://www.kegg.jp/pathway/hsa00140 |
| **Content** | Steroid hormone biosynthesis pathway map |
| **Compounds** | C27 cholesterol to C21, C19, C18 steroids |
| **Key Enzymes** | CYP11A1 (cholesterol side-chain cleavage), CYP17, CYP19 |
| **License** | Academic free; commercial license required |
| **API** | KEGG REST API (academic use) |
| **Update Frequency** | Monthly |
| **Priority** | Tier 2 |
| **Storage Estimate** | Included in KEGG total (~500 MB) |

#### PharmVar - Pharmacogene Variation

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmvar.org/ |
| **Content** | Pharmacogene allele nomenclature |
| **Records** | 29+ CYP enzymes + POR |
| **Steroid CYPs** | CYP11A1, CYP11B1, CYP11B2, CYP17A1, CYP19A1, CYP21A2 |
| **License** | Open access |
| **API** | Data files available for download |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~200 MB |
| **Formats** | VCF, FASTA, standardized definitions |
| **Partners** | PharmGKB, CPIC |

#### Harmonizome - Steroid Gene Sets

| Field | Value |
|-------|-------|
| **URL** | https://maayanlab.cloud/Harmonizome/gene_set/c21+steroid+hormone+metabolism/KEGG+Pathways |
| **Content** | Gene sets for steroid metabolism pathways |
| **Records** | 11 proteins in c21 steroid hormone metabolism |
| **License** | Open access |
| **API** | REST API |
| **Update Frequency** | Periodic |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~50 MB |

---

### Category 7: Cross-Disease & Integration Resources

#### Endo-ERN (European Reference Network)

| Field | Value |
|-------|-------|
| **URL** | https://endo-ern.eu/ |
| **Content** | Rare endocrine conditions guidelines and resources |
| **Focus** | NGS implementation, genetic testing standards |
| **Value** | Clinical interpretation guidelines for endocrine variants |
| **Priority** | Reference only |

#### Gene Ontology - Steroid/Hormone Terms

| Field | Value |
|-------|-------|
| **URL** | https://www.informatics.jax.org/vocab/gene_ontology/GO:0008202 |
| **Term** | GO:0008202 (steroid metabolic process) |
| **Content** | Standardized gene function annotations |
| **License** | CC BY 4.0 |
| **API** | AmiGO API |
| **Priority** | Tier 2 |

---

## Storage Summary

| Category | Sources | Est. Size |
|----------|---------|-----------|
| Thyroid Genetics | 3 | ~5 GB |
| Hormone Pathways | 4 | ~10 GB |
| Autoimmune Genetics | 5 | ~50 GB |
| HLA/Immune System | 3 | ~15 GB |
| Endocrine Disorders | 4 | ~5 GB |
| Steroid Metabolism | 4 | ~3 GB |
| **Total** | **23** | **~88 GB** |

*Note: gnomAD full data is ~30 TB but only summary data (~50 GB) recommended for initial integration*

---

## Integration Recommendations

### Tier 1 - Essential (High Value, Open Access)

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | IPD-IMGT/HLA | Gold standard for HLA nomenclature, WHO official, CC BY-ND |
| 2 | GWAS Catalog | Comprehensive autoimmune GWAS data, open access, REST API |
| 3 | Reactome | Hormone pathway data, CC BY 4.0, excellent API |
| 4 | ThyroidOmics | Largest thyroid-specific GWAS resource (271K individuals) |
| 5 | PharmVar | Steroid pharmacogenomics, open access, VCF format |
| 6 | gnomAD | Population frequencies, HLA region coverage |

### Tier 2 - Important

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | OpenGWAS | Bulk GWAS data access, 126B associations |
| 2 | GAAD | Autoimmune gene associations, 44K associations |
| 3 | Hmrbase2 | Hormone-receptor relationships, 12K entries |
| 4 | Open Targets | Drug target integration, GraphQL API |
| 5 | IPD-KIR | KIR sequences, complements HLA data |
| 6 | OMIM | Endocrine Mendelian disorders |

### Tier 3 - Supplementary

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | ADEx | Autoimmune expression data, 82 studies |
| 2 | EndoGene | Endocrine variant collection |
| 3 | KEGG | Pathway maps (academic use) |
| 4 | MalaCards | Disease integration |
| 5 | HormoneBase | Comparative steroid data |
| 6 | Harmonizome | Gene set enrichment |

---

## API Integration Notes

### REST APIs Available

| Database | API Endpoint |
|----------|--------------|
| GWAS Catalog | https://www.ebi.ac.uk/gwas/rest/docs/api |
| OpenGWAS | https://api.opengwas.io/api/ |
| Reactome | https://reactome.org/ContentService/ |
| gnomAD | GraphQL at https://gnomad.broadinstitute.org/api |
| Open Targets | GraphQL at https://api.platform.opentargets.org/api/v4/graphql |

### R Packages

| Package | Database |
|---------|----------|
| gwasrapidd | GWAS Catalog |
| ieugwasr | OpenGWAS |
| ReactomePA | Reactome pathway analysis |

### Python Packages

| Package | Database |
|---------|----------|
| ieugwaspy | OpenGWAS |
| opentargets | Open Targets Platform |

---

## License Summary

### Fully Open Access

| Database | License | Download Method |
|----------|---------|-----------------|
| Reactome | CC BY 4.0 | FTP, API |
| GWAS Catalog | Open | FTP, API |
| OpenGWAS | Open | API, Files |
| IPD-IMGT/HLA | CC BY-ND | FTP, GitHub |
| gnomAD | Open | AWS S3 |
| GAAD | Free | Web |
| PharmVar | Open | Files |
| Gene Ontology | CC BY 4.0 | FTP, API |

### Academic Use Only

| Database | License | Notes |
|----------|---------|-------|
| KEGG | Academic free | Commercial license required for business use |
| Hmrbase2 | Academic | Web interface only |
| OMIM | Academic free | API key required, yearly renewal |
| MalaCards | Academic | Web interface primary |
| ThyroidOmics | Academic | Check individual datasets |
| ADEx | Academic | Contact for commercial use |

### Deprecated/Limited

| Database | Status | Alternative |
|----------|--------|-------------|
| ImmunoBase | Deprecated | GWAS Catalog |
| T1DBase | Deprecated | GWAS Catalog, Open Targets |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-00-INDEX.md](./43-00-INDEX.md) | Parent index |
| [43-11-GENETICS-PRIMARY.md](./43-11-GENETICS-PRIMARY.md) | Overlaps with dbSNP, ClinVar |
| [43-12-GENETICS-POPULATION.md](./43-12-GENETICS-POPULATION.md) | Overlaps with gnomAD |
| [43-41-PATHWAYS-PRIMARY.md](./43-41-PATHWAYS-PRIMARY.md) | Overlaps with Reactome |
| [43-42-PATHWAYS-METABOLISM.md](./43-42-PATHWAYS-METABOLISM.md) | Overlaps with KEGG |
| [43-51-PHARMACEUTICALS.md](./43-51-PHARMACEUTICALS.md) | Overlaps with PharmVar |
| [44-ARCHITECTURE.md](../44-ARCHITECTURE.md) | Informs database design |
| [45-DATA-MODEL.md](../45-DATA-MODEL.md) | Informs entity structure |

---

## Open Questions

- [ ] HLA typing depth - 2-field vs 4-field resolution for disease associations?
- [ ] gnomAD subset - which immune regions to prioritize for initial load?
- [ ] ThyroidOmics data sharing agreements - review per-cohort requirements
- [ ] KEGG license - confirm academic use sufficient for platform purposes
- [ ] TCGA controlled access - NIH dbGaP application timeline?

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial migration from research.old/data-sources-hormones-autoimmune.md |
