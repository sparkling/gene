---
id: domain-cancer-oncology
title: "Cancer & Oncology Data Sources"
type: health-domain
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [disease, health-domain, databases]
---

# Cancer & Oncology Data Sources

**Document ID:** 43-73-CANCER-ONCOLOGY
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [_index.md](./_index.md)

---

## TL;DR

Cancer and oncology databases provide somatic mutation data (COSMIC: 38M+ mutations), clinical interpretation (OncoKB, CIViC), multi-omics research data (GDC/TCGA: 2.5PB), and hereditary cancer variant classifications (ClinVar, BRCA Exchange). CIViC offers fully open CC0 data ideal for research; OncoKB and COSMIC require licensing for clinical/commercial use. Priority integration: ClinVar (hereditary) + CIViC (open clinical) + GDC (research backbone).

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary somatic mutation reference | COSMIC | Largest resource (38M+ mutations), comprehensive coverage | Jan 2026 |
| Open clinical interpretation source | CIViC | CC0 license, community-curated, GraphQL API | Jan 2026 |
| Clinical actionability source | OncoKB | FDA-recognized, MSK-curated, evidence levels | Jan 2026 |
| Research data backbone | GDC/TCGA | 2.5PB multi-omics, harmonized pipelines | Jan 2026 |
| Hereditary cancer variants | ClinVar + BRCA Exchange | Public domain, expert classifications | Jan 2026 |
| Lynch syndrome focus | InSiGHT Database | Specialized MMR gene expertise | Jan 2026 |

---

## Database Catalog

### 1. COSMIC (Catalogue Of Somatic Mutations In Cancer)

| Field | Value |
|-------|-------|
| **URL** | https://cancer.sanger.ac.uk/cosmic |
| **Content** | World's largest somatic mutation database: coding mutations, gene fusions, mutational signatures, drug resistance, CNVs, non-coding mutations |
| **Records** | 38M+ mutations, 1.4M+ samples, 25M+ genomic variants across 6,800+ cancer types |
| **License** | Academic research: Free (registration required); Patient testing/Commercial: License fee via QIAGEN |
| **API** | GA4GH Beacon; NIH Clinical Tables API (v89, GRCh37 only, limited); Download files (TSV/VCF) |
| **Update Frequency** | Quarterly (current: v103, November 2025) |
| **Priority** | Tier 2 (licensing considerations) |
| **Storage Estimate** | ~1 GB RAM to load; multi-GB download files |

**Data Products:**
- COSMIC Core: Complete somatic mutation dataset
- Cell Lines Project: Cancer cell line mutations
- Cancer Mutation Census (CMC): Clinically actionable mutations
- Actionability: Drug-variant associations
- Cancer Gene Census (CGC): Curated cancer driver genes

**Download Formats:** TSV, VCF (GRCh37/GRCh38), 3 previous releases available

**Integration Notes:**
- Python access via `gget cosmic` package
- Google Cloud Life Sciences: Public dataset available
- Commercial licensing may limit redistribution
- Contact: cosmic@sanger.ac.uk

---

### 2. GDC / TCGA (Genomic Data Commons)

| Field | Value |
|-------|-------|
| **URL** | https://portal.gdc.cancer.gov/ |
| **Content** | NCI harmonized genomic and clinical data: TCGA (33 cancer types), TARGET (pediatric), CPTAC (proteomics), CGCI, HCMI |
| **Records** | 2.5 petabytes; TCGA: 11,000+ patients across 33 cancer types |
| **License** | Open Access: No authentication; Controlled Access: dbGaP authorization + eRA Commons |
| **API** | Full REST API (JSON), GraphQL; R: `GenomicDataCommons`; Python: `gdclient` |
| **Update Frequency** | Continuous (current: Data Release 43.0, May 2025) |
| **Priority** | Tier 1 (research backbone) |
| **Storage Estimate** | Individual files: 50 MB to 200-300 GB; Total: 2.5 PB |

**Data Types:**
| Type | Description |
|------|-------------|
| WGS/WXS | Whole genome/exome sequencing |
| RNA-Seq | Gene expression data |
| miRNA-Seq | MicroRNA expression |
| DNA Methylation | Epigenetic profiles (SeSAMe pipeline) |
| SNP Arrays | Copy number data |
| RPPA | Protein expression |
| Clinical | Patient demographics, outcomes |

**Download Methods:**
1. GDC Data Portal: Web-based browsing and cart system
2. GDC Data Transfer Tool (DTT) v2.3.0: CLI for large files (Linux, macOS, Windows)
3. API Downloads: Programmatic access
4. TCGADownloadHelper: Simplified extraction (GitHub: alex-baumann-ur/TCGADownloadHelper)

**Integration Notes:**
- Standard formats: BAM, VCF, MAF, TSV
- GDC harmonization pipelines ensure consistency
- Metadata-rich with standardized ontologies

---

### 3. Cancer Gene Census (CGC)

| Field | Value |
|-------|-------|
| **URL** | https://cancer.sanger.ac.uk/census |
| **Content** | High-confidence, expertly curated catalog of genes causally implicated in cancer |
| **Records** | ~700 genes (Tier 1: strong evidence; Tier 2: mutational patterns) |
| **License** | Same as COSMIC (Academic free with registration; Commercial via QIAGEN) |
| **API** | Part of COSMIC data downloads; CSV/TSV export from web interface |
| **Update Frequency** | Regular updates by COSMIC curation team |
| **Priority** | Tier 2 (gene-level filtering) |
| **Storage Estimate** | <10 MB (curated gene list) |

**Data Fields:**
| Field | Description |
|-------|-------------|
| Gene Symbol | HGNC symbol |
| Tier | 1 (high confidence) or 2 |
| Hallmark | Cancer hallmarks affected |
| Mutation Types | Point, amplification, translocation, etc. |
| Role in Cancer | Oncogene, TSG, or both |
| Cancer Types | Associated malignancies |

**Acceptance Criteria:**
- Literature-based evidence required
- Independent confirmation from 2+ groups
- Clear mutation patterns in specified diseases
- Statistical interpretations excluded (conservative approach)

**Integration Notes:**
- Often used as filter/annotation for variant analysis
- Compare with OncoKB gene list (1,216 genes as of Dec 2025)
- Foundation for many cancer gene panels

---

### 4. OncoKB

| Field | Value |
|-------|-------|
| **URL** | https://www.oncokb.org/ |
| **Content** | MSK Precision Oncology Knowledge Base; FDA-recognized human genetic variant database for clinical interpretation |
| **Records** | 978 genes, 10,623 alterations, 144 cancer types, 159 drugs, 1,216 cancer gene list |
| **License** | Academic research: Free (registration); Commercial/Clinical: Fee-based; AI/ML training: Prohibited |
| **API** | Full REST API (token required); Demo: https://demo.oncokb.org (BRAF, TP53, ROS1 only) |
| **Update Frequency** | Frequent (17 novel biomarkers added in 2024); Germline variants expected 2025 |
| **Priority** | Tier 2 (licensing for clinical use) |
| **Storage Estimate** | ~100 MB (focused curation) |

**Evidence Levels:**
| Level | Description |
|-------|-------------|
| Level 1 | FDA-recognized biomarker |
| Level 2 | Standard care biomarker |
| Level 3A | Compelling clinical evidence |
| Level 3B | Standard care in another cancer type |
| Level 4 | Compelling biological evidence |
| Level R1 | Standard care resistance |
| Level R2 | Compelling resistance evidence |

**API Endpoints:**
| Type | Description |
|------|-------------|
| Protein Change | Gene + alteration + tumor type |
| Copy Number | Amplification, deletion, gain, loss |
| Structural Variants | Fusions, rearrangements |
| Genomic Coordinates | TCGA MAF or HGVS format |

**Public Downloads (no registration):**
- Cancer Gene List
- All Curated Genes List
- Biomarker-Drug Association List

**Integration Notes:**
- Strongly recommends API over downloads (always current)
- Downloads can become outdated quickly
- Umbrella terms use logic that requires API
- Contact: contact@oncokb.org

---

### 5. CIViC (Clinical Interpretation of Variants in Cancer)

| Field | Value |
|-------|-------|
| **URL** | https://civicdb.org/ |
| **Content** | Open-access, community-driven knowledgebase for clinical variant interpretations |
| **Records** | 3,200+ variants, 470+ genes, 3,100+ publications, 300+ contributors |
| **License** | Data: CC0 1.0 Universal (Public Domain); Source Code: MIT License |
| **API** | GraphQL API (primary), REST API v2, GraphiQL interactive, CIViCpy Python SDK |
| **Update Frequency** | Nightly TSV/CSV dumps; Monthly full releases |
| **Priority** | Tier 1 (fully open) |
| **Storage Estimate** | ~50 MB (focused expert curation) |

**Evidence Types:**
| Type | Description |
|------|-------------|
| Predictive | Drug response predictions |
| Diagnostic | Supports diagnosis |
| Prognostic | Outcome predictions |
| Predisposing | Germline risk factors |
| Oncogenic | Functional oncogenicity |

**Download Options:**
- Nightly TSV/CSV dumps on AWS Open Data
- Monthly releases: Full data exports
- Formats: JSON, TSV, VCF

**Integration Notes:**
- Completely open with no restrictions
- MCP integration available (CIViC MCP, 2025)
- Excellent for research and clinical decision support
- Community contribution model
- AWS Registry: https://registry.opendata.aws/civic/

---

### 6. ClinVar (Hereditary Cancer Variants)

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | NCBI public database of human genetic variants and disease relationships; primary resource for hereditary cancer variant classification |
| **Records** | 3M+ variants, 2,800+ submitting organizations; includes oncogenicity classifications (new 2025) |
| **License** | Public domain (US Government work) |
| **API** | E-utilities, Clinical Tables API |
| **Update Frequency** | Weekly/monthly releases |
| **Priority** | Tier 1 (foundational) |
| **Storage Estimate** | Multi-GB download files |

**Hereditary Cancer Genes Covered:**
- BRCA1, BRCA2 (breast/ovarian)
- MLH1, MSH2, MSH6, PMS2, EPCAM (Lynch syndrome)
- APC, MUTYH (polyposis syndromes)
- TP53 (Li-Fraumeni)
- PTEN (Cowden syndrome)
- CDH1 (hereditary diffuse gastric cancer)
- PALB2, CHEK2, ATM (moderate penetrance)

**Download Formats:**
- VCF (GRCh37, GRCh38)
- XML (VCV, RCV formats)
- Tab-delimited summaries

**Integration Notes:**
- Simple ClinVar available: https://simple-clinvar.broadinstitute.org/
- Foundation for all hereditary variant classification
- Free to use, download, redistribute

---

### 7. InSiGHT Database

| Field | Value |
|-------|-------|
| **URL** | https://www.insight-group.org/mutations |
| **Content** | International Society for Gastrointestinal Hereditary Tumours database; Lynch syndrome and related hereditary GI cancer syndromes |
| **Records** | ~3,000 unique germline variants; MLH1 (40%), MSH2 (34%), MSH6 (18%), PMS2 (8%); also APC, MUTYH |
| **License** | Academic access (registration may be required); data sharing with ClinVar |
| **API** | Web-based query interface; export capabilities |
| **Update Frequency** | Regular updates |
| **Priority** | Tier 3 (specialized) |
| **Storage Estimate** | <50 MB |

**Classification System:**
- 5-tiered IARC classification
- Linked to clinical recommendations
- Based on variant characteristics, family data, functional assay results

**Related Resources:**
- PLSD (Prospective Lynch Syndrome Database): Cancer risk estimates by age, gender, gene

---

### 8. BRCA Exchange

| Field | Value |
|-------|-------|
| **URL** | https://brcaexchange.org/ |
| **Content** | Global resource for BRCA1/BRCA2 variant data aggregation and expert classification (GA4GH BRCA Challenge) |
| **Records** | 20,000+ unique BRCA1/2 variants, 6,100+ expert-classified, ~3,700 pathogenic |
| **License** | Open access; freely queryable, downloadable, redistributable |
| **API** | REST API; PostgreSQL backend; VR specification support |
| **Update Frequency** | Monthly updates with change tracking |
| **Priority** | Tier 2 (BRCA-specific expertise) |
| **Storage Estimate** | ~100 MB |

**Data Sources Aggregated:**
| Source | Type |
|--------|------|
| ClinVar | Clinical submissions |
| LOVD | Locus-specific database |
| BIC | Breast Cancer Information Core |
| ExAC/gnomAD | Population frequencies |
| 1000 Genomes | Population data |
| ENIGMA | Expert classifications |

**Download Options:**
- Full variant data via API or web portal
- Historical data releases (tarballs)
- Pipeline intermediate outputs included

---

## Summary Comparison

| Database | Focus | Records | License | API | Best For |
|----------|-------|---------|---------|-----|----------|
| COSMIC | Somatic mutations | 38M+ | Academic free, Commercial paid | Limited | Somatic variant annotation |
| GDC/TCGA | Multi-omics cancer | 2.5 PB | Open/Controlled | Full REST | Research data access |
| Cancer Gene Census | Cancer driver genes | ~700 genes | Same as COSMIC | Via COSMIC | Gene filtering |
| OncoKB | Clinical actionability | 10K+ alterations | Academic free, Clinical paid | Full REST | Treatment decisions |
| CIViC | Clinical interpretation | 3K+ variants | CC0 (open) | GraphQL | Research, clinical support |
| ClinVar | Germline classification | 3M+ | Public domain | E-utilities | Variant classification |
| InSiGHT | Lynch syndrome | 3K+ variants | Academic | Web | Lynch syndrome |
| BRCA Exchange | BRCA1/2 | 20K+ variants | Open | REST | BRCA interpretation |

---

## Integration Recommendations

### Tier 1: Essential (MVP)

| Database | Rationale |
|----------|-----------|
| **ClinVar** | Foundation for all variant classification; public domain; hereditary cancer focus |
| **CIViC** | Open CC0 license; clinical interpretations; excellent API |
| **GDC/TCGA** | Research data backbone; harmonized multi-omics |

### Tier 2: High Value (Post-MVP)

| Database | Rationale |
|----------|-----------|
| **OncoKB** | Clinical actionability; FDA-recognized; requires licensing for clinical use |
| **COSMIC** | Somatic mutation reference; licensing considerations for redistribution |
| **BRCA Exchange** | BRCA-specific expertise; aggregated expert classifications |

### Tier 3: Specialized

| Database | Rationale |
|----------|-----------|
| **Cancer Gene Census** | Gene-level filtering; part of COSMIC ecosystem |
| **InSiGHT** | Lynch syndrome focus; specialized MMR gene expertise |

---

## Technical Integration Considerations

### Data Harmonization Challenges

- **Variant nomenclature**: HGVS, protein, genomic coordinate differences
- **Reference genome**: GRCh37 vs GRCh38 liftover requirements
- **Classification systems**: ACMG, AMP/ASCO/CAP oncology-specific guidelines
- **Update frequency**: Varies from daily (CIViC) to quarterly (COSMIC)

### Recommended Annotation Tools

| Tool | Description |
|------|-------------|
| VEP | Variant Effect Predictor (Ensembl annotation) |
| ANNOVAR | Multi-database annotation |
| OpenCRAVAT | Modular annotation system |
| ClinGen Allele Registry | Variant normalization |

### API Rate Limits

| Database | Limit Notes |
|----------|-------------|
| COSMIC | Registration required |
| GDC | Generally unrestricted |
| OncoKB | Token-based, contact for limits |
| CIViC | Open, reasonable use |
| ClinVar | E-utilities guidelines apply |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `Somatic mutation` | A mutation acquired during an organism's lifetime, not inherited | BRAF V600E in melanoma |
| `Germline` | Inherited genetic variations present in egg/sperm cells | BRCA1/2 mutations |
| `Pathogenic` | Disease-causing variant classification under ACMG guidelines | ClinVar P/LP variants |
| `VUS` | Variant of Uncertain Significance - insufficient evidence to classify | ~40% of ClinVar variants |
| `CNV` | Copy Number Variation - duplication or deletion of DNA segments | Gene amplifications |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `Oncogene` | Gene that promotes cancer when mutated or overexpressed | KRAS, BRAF, MYC |
| `TSG` | Tumor Suppressor Gene - gene that normally prevents cancer | TP53, BRCA1, APC |
| `Driver mutation` | Mutation that directly contributes to cancer development | Identified in 97% of THCA |
| `Actionability` | Whether a variant has therapeutic implications | OncoKB evidence levels |
| `Lynch syndrome` | Hereditary cancer syndrome from MMR gene mutations | MLH1, MSH2, MSH6, PMS2 |
| `MMR genes` | Mismatch Repair genes that correct DNA replication errors | InSiGHT database focus |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| COSMIC | Catalogue Of Somatic Mutations In Cancer | 38M+ mutations, Sanger Institute |
| TCGA | The Cancer Genome Atlas | NCI/NHGRI cancer genomics project |
| GDC | Genomic Data Commons | NCI harmonized data portal |
| CIViC | Clinical Interpretation of Variants in Cancer | CC0 open-access database |
| OncoKB | Oncology Knowledge Base | MSK precision oncology, FDA-recognized |
| CGC | Cancer Gene Census | ~700 curated cancer driver genes |
| InSiGHT | International Society for Gastrointestinal Hereditary Tumours | Lynch syndrome database |
| ACMG | American College of Medical Genetics | Variant classification guidelines |
| AMP | Association for Molecular Pathology | Oncology variant guidelines |
| CAP | College of American Pathologists | Clinical laboratory standards |
| MAF | Mutation Annotation Format | Standard somatic mutation file |
| WGS | Whole Genome Sequencing | Complete genome analysis |
| WXS | Whole Exome Sequencing | Protein-coding region analysis |
| HGVS | Human Genome Variation Society | Variant nomenclature standard |
| GA4GH | Global Alliance for Genomics and Health | Data sharing standards |
| VEP | Variant Effect Predictor | Ensembl annotation tool |
| CC0 | Creative Commons Zero | Public domain dedication |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [primary.md](./../genetics/primary.md) | ClinVar also covered there (germline focus) |
| [disease.md](./../pathways/disease.md) | Disease ontology linking |
| [pharmaceuticals.md](./../compounds/pharmaceuticals.md) | Drug-variant associations |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial migration from research.old/data-sources-cancer-oncology.md |
