---
id: domain-rare
title: "Rare Diseases & Biobank Data Sources"
type: health-domain
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [disease, health-domain, databases]
---

# Rare Diseases & Biobank Data Sources

**Document ID:** 43-75-RARE-DISEASES
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [_index.md](./_index.md)

---

## TL;DR

Rare disease data integration centers on Orphanet (6,500+ diseases, CC-BY-4.0), ClinGen (11,000+ variant curations), and DECIPHER (51,000+ patients) for clinical genomics. Phenotype standardization uses HPO (13,000+ terms) and MONDO (26,000+ diseases). Biobank data from UK Biobank (500K participants) and BioBank Japan (270K patients) provide population-scale GWAS statistics, with individual-level data requiring formal access agreements.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary rare disease ontology | Orphanet (ORDO) | Most comprehensive, CC-BY-4.0, ELIXIR Core Resource | Jan 2026 |
| Phenotype standard | HPO + MONDO | Industry standard, multi-language support, precise equivalence | Jan 2026 |
| Variant interpretation source | ClinGen | NIH-funded, ACMG/AMP guidelines, daily updates | Jan 2026 |
| Biobank access | GWAS summary statistics only | Individual-level requires formal agreements | Jan 2026 |
| Cross-species integration | Monarch Initiative | Enables phenotype-to-genotype mapping across species | Jan 2026 |

---

## Database Catalog

### 1. Orphanet (Rare Disease Portal)

| Field | Value |
|-------|-------|
| **URL** | https://www.orpha.net/ |
| **Content** | Rare disease classifications, gene associations, epidemiology, expert centres |
| **Records** | 6,528 rare diseases, 4,512 linked genes, 36,595 diagnostic tests |
| **License** | CC-BY-4.0 |
| **API** | REST API (https://api.orphadata.com/) |
| **Update Frequency** | Bi-annual (July and December) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB (ORDO ontology + datasets) |

**Data Types:**
- Disease classifications (hierarchical)
- Gene-disease associations
- Epidemiological data (prevalence, incidence)
- Expert centres directory (8,722 centres)
- Orphan drug information
- Diagnostic test catalog

**Download Options:**
| Dataset | Format | Size | Description |
|---------|--------|------|-------------|
| ORDO (Ontology) | OWL | ~50 MB | Structured vocabulary for rare diseases |
| Disease Classifications | XML, JSON | Variable | Hierarchical disease organization |
| Gene-Disease Associations | XML | Variable | Gene linkages to rare diseases |
| HPO-ORDO Module (HOOM) | OWL | Variable | Phenotype mappings |
| Orphapackets | JSON | Variable | Structured disease packets |

**Integration Notes:**
- Integrates with: OMIM, MONDO, UniProtKB, HGNC, Ensembl, Reactome, IUPHAR
- Recognized by: ELIXIR Core Data Resource, Global Core Biodata Resource, IRDiRC
- Current ORDO version: 4.8 (December 2025)
- Languages: English, Czech, Dutch, French, German, Italian, Spanish, Polish

---

### 2. DECIPHER (Clinical Genomics)

| Field | Value |
|-------|-------|
| **URL** | https://www.deciphergenomics.org/ |
| **Content** | Chromosomal imbalances, phenotype-variant mappings |
| **Records** | 51,894 open-access patients, 250+ contributing projects |
| **License** | Academic use with Data Access Agreement |
| **API** | Deposition API (client keys available) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 GB (bulk download) |

**Data Types:**
- Copy Number Variants (CNVs)
- Sequence variants (SNVs, insertions, deletions, indels)
- Phenotype data (linked to HPO terms)
- Gene information with GeneReviews annotations
- CNV syndromes
- De-identified patient cases

**Access Tiers:**
| Tier | Visibility | Requirements |
|------|------------|--------------|
| Center-Only | Single center | Default upload |
| Open Access | Public/World | Explicit patient consent |
| Data Display Agreement | External genome browsers | Signed agreement |
| Bulk Research Access | Academic researchers | Data Access Agreement |

**Search Capabilities:**
- Genomic coordinates (GRCh38 and GRCh37)
- Gene symbols and transcripts
- HPO phenotype terms and identifiers
- Variant nomenclature (HGVSg, HGVSc, HGVSp)
- Patient IDs and cytogenetic bands

**Integration Notes:**
- Displays in: Ensembl Genome Browser, UCSC Genome Browser
- Located at: EMBL-EBI, Wellcome Genome Campus, UK
- Moved to EMBL-EBI in July 2023

---

### 3. ClinGen (Clinical Genome Resource)

| Field | Value |
|-------|-------|
| **URL** | https://clinicalgenome.org/ |
| **Content** | Gene-disease validity, variant pathogenicity, clinical actionability |
| **Records** | 3,289 gene-disease curations, 11,133 variant curations, 2,719 genes |
| **License** | Open access (citation required) |
| **API** | Multiple REST APIs (Evidence Repository, Allele Registry, etc.) |
| **Update Frequency** | Daily |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |

**Data Types:**
- Gene-disease validity assessments
- Variant pathogenicity (ACMG/AMP classifications)
- Dosage sensitivity (1,557 genes)
- Clinical actionability (447 gene-condition pairs)
- Somatic cancer variants
- ACMG secondary findings curations

**Download Formats:**
| Format | Use Case |
|--------|----------|
| CSV | Spreadsheet analysis |
| TSV | Tab-delimited processing |
| BED | Genomic coordinate visualization |
| AED | ChAS-compatible array analysis |
| JSON (nested/flat) | Programmatic access |

**API Endpoints:**
| API | Description |
|-----|-------------|
| Evidence Repository API | Variant pathogenicity data |
| Allele Registry API | Variant data retrieval |
| Criteria Specification Registry API | Specification data |
| Linked Data Hub API | Aggregated variant annotations |
| Clinical Actionability APIs | Context-specific JSON |

**Integration Notes:**
- Works with: ClinVar, GenomeConnect
- Expert panels use ACMG/AMP guidelines with gene-specific specifications
- 2,800+ active contributors from 74+ countries
- FTP access retains historical downloads for 60 days

---

### 4. UK Biobank

| Field | Value |
|-------|-------|
| **URL** | https://www.ukbiobank.ac.uk/ |
| **Content** | Large-scale biomedical database with genetic, imaging, health data |
| **Records** | 500,000 participants, 10,000+ data variables, 30+ petabytes |
| **License** | Restricted (controlled access only) |
| **API** | None for individual-level data; summary statistics may be public |
| **Update Frequency** | Continuous additions |
| **Priority** | Tier 3 (summary stats only) |
| **Storage Estimate** | N/A (cloud-based access via UKB-RAP) |

**Data Types:**
- Genetic data (whole genome sequencing, genotyping arrays)
- Metabolomic data (~250 blood metabolites, all 500K participants)
- Imaging data (brain, cardiac, body - targeting 100,000 by end of 2025)
- Biological samples (blood, urine, saliva)
- Questionnaires (health, lifestyle, diet, environment)
- Physical measurements (body composition, bone density)
- Activity tracking (accelerometer data)
- Follow-up data (health outcomes, hospital records, death registry)

**Access Model:**
- Primary: UK Biobank Research Analysis Platform (UKB-RAP)
- Secure cloud-based platform (default access since 2024)
- Requires: researcher approval, institutional affiliation, access committee approval, data usage agreement

**Public Data Availability:**
| Type | Availability |
|------|--------------|
| Summary Statistics | Publicly accessible |
| Individual-Level Data | Controlled access only |
| Published Results | Available via publications |

**Recent Updates (2025):**
- Metabolomic data expanded to all 500,000 participants
- 55,000 3D heart models added
- Additional 30,000 bone density measures
- Brain structure longitudinal data
- Imaging project approaching 100,000 participants

**Integration Notes:**
- No public API for individual-level data
- 22,000 researchers actively using data worldwide
- 30,000+ registered researchers from 90+ countries
- 9,000+ peer-reviewed publications
- Charity registered in England/Wales (1101332) and Scotland (SC039230)

---

### 5. BioBank Japan

| Field | Value |
|-------|-------|
| **URL** | https://biobankjp.org/en/ |
| **Content** | Disease-oriented prospective biobank for Japanese population |
| **Records** | ~270,000 patients, 51 diseases, 800K DNA tubes, 1.7M serum samples |
| **License** | GWAS summary statistics public; individual data controlled |
| **API** | None for individual data; GWAS via PheWeb.jp |
| **Update Frequency** | Periodic releases |
| **Priority** | Tier 2 (GWAS stats) |
| **Storage Estimate** | ~50 GB (GWAS summary statistics) |

**Data Types:**
- DNA samples (800,000 tubes)
- Serum samples (1,700,000 tubes from 200,000 patients)
- Whole genome sequencing (16,000 patients)
- SNP genotyping data (270,000 patients)
- Clinical information linked to samples
- GWAS summary statistics (publicly available)

**Public GWAS Data (PheWeb.jp):**
| Metric | Value |
|--------|-------|
| GWAS Results | 220+ |
| Study Population | ~260,000 participants |
| Ancestry | Japanese |
| Available Since | March 2021 |

**Data Categories:**
- Metabolic disorders (Type 2 diabetes)
- Autoimmune conditions (rheumatoid arthritis, psoriasis)
- Cancer types (breast cancer)
- Infectious disease responses (COVID-19)
- Neurological conditions
- Dietary phenotypes

**Key Publications:**
| Study | Sample Size | Findings |
|-------|-------------|----------|
| Quantitative Traits GWAS | 162,255 | 1,407 trait-associated loci (679 novel) |
| 42 Diseases GWAS | 212,453 | 320 independent signals, 25 novel loci |

**Certifications:**
- ISO 9001:2015 (Quality Management)
- ISO/IEC 27001:2013 (Information Security Management)

**Integration Notes:**
- NBDC Human Database hosts controlled data
- April 2024: URL changes for NBDC-hosted data
- Primarily serves Japanese population research
- 12 contributing medical centers

---

### 6. Human Phenotype Ontology (HPO)

| Field | Value |
|-------|-------|
| **URL** | https://hpo.jax.org/ |
| **Content** | Standardized vocabulary of phenotypic abnormalities |
| **Records** | 13,000+ terms, 156,000+ disease annotations |
| **License** | Open access (https://hpo.jax.org/app/license) |
| **API** | GitHub releases, OBO Foundry |
| **Update Frequency** | ~Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 MB |

**Download Formats:**
| Format | URL |
|--------|-----|
| OWL | http://purl.obolibrary.org/obo/hp.owl |
| OBO | http://purl.obolibrary.org/obo/hp.obo |
| JSON | Available via GitHub releases |

**Statistics (2024):**
| Metric | Value |
|--------|-------|
| Total Terms | 13,000+ |
| Annotations to Diseases | 156,000+ |
| New Terms Added (recent) | 2,239 |
| New Annotations Added | 49,235 |
| GitHub Commits | 8,308 |
| Releases | 61 |
| Contributors | 24+ |

**Languages Supported:**
English, German, Spanish, French, Italian, Dutch, Portuguese, Turkish, Japanese, Chinese + more

**Recent Developments:**
- Medical Action Ontology (MAxO) for treatment modeling
- GA4GH Phenopacket Schema integration
- EHR integration efforts
- Internationalization with 10+ language translations

**Integration Notes:**
- Maintained by Monarch Initiative and Jackson Laboratory
- Available via EBI OLS: https://www.ebi.ac.uk/ols/ontologies/hp
- GitHub: https://github.com/obophenotype/human-phenotype-ontology

---

### 7. Monarch Initiative

| Field | Value |
|-------|-------|
| **URL** | https://monarchinitiative.org/ |
| **Content** | Integrative platform connecting phenotypes to genotypes across species |
| **Records** | 33 integrated data sources |
| **License** | Open access (individual datasets may vary) |
| **API** | REST API v3 (https://api-v3.monarchinitiative.org/v3/docs) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 GB (knowledge graph) |

**Data Sources Integrated:**
- Model organisms: MGD, ZFIN, WormBase, FlyBase, Xenbase, SGD, PomBase, DictyBase, BGeeDB
- Human disease: OMIM, Orphanet
- Pathways: Reactome
- Interactions: STRING
- Ontologies: HPO, Gene Ontology, PHENIO

**Available Datasets:**
| Dataset | Description |
|---------|-------------|
| Monarch KG | Main knowledge graph |
| HPO | Human Phenotype Ontology |
| Exomiser | Variant prioritization resources |
| Phenologs | Phenotype ortholog data |
| Semantic Similarity | Phenotype similarity metrics |
| Mappings | Cross-reference mappings |
| UPHENO2 | Unified Phenotype Ontology |

**API Features:**
- Version: API v3 (launched October 2023)
- Documentation: OpenAPI/Swagger
- Concise abstracted endpoints for all entity combinations

**Integration Notes:**
- ChatGPT plugin available for phenotypic data queries
- Previous API (Biolink, OWLSim) discontinued March 2024
- R package available: monarchr
- Downloads: https://data.monarchinitiative.org/

---

### 8. MONDO Disease Ontology

| Field | Value |
|-------|-------|
| **URL** | https://mondo.monarchinitiative.org/ |
| **Content** | Unified disease ontology with precise equivalence axioms |
| **Records** | 25,938 disease concepts (22,977 human, 2,960 non-human) |
| **License** | CC-BY-4.0 |
| **API** | GitHub releases |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~200 MB |

**Download Formats:**
| Format | File | Description |
|--------|------|-------------|
| OWL | mondo-with-equivalents.owl | With equivalence axioms |
| OBO | mondo.obo | Using xrefs for linking |
| JSON | mondo-with-equivalents.json | Equivalent to OWL |

**Key Feature:**
Provides "precise 1:1 equivalence axioms" validated by OWL reasoning, distinguishing it from loose cross-references in other resources.

**Integration Notes:**
- GitHub: https://github.com/monarch-initiative/mondo
- Latest releases: https://github.com/monarch-initiative/mondo/releases/latest

---

### 9. Related Resources (Reference)

#### ClinVar
| Attribute | Value |
|-----------|-------|
| URL | https://www.ncbi.nlm.nih.gov/clinvar/ |
| FTP | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/ |
| Formats | VCF (GRCh37/38), XML, Tab-delimited |
| License | Public domain (NCBI) |
| Updates | Daily |

#### OMIM
| Attribute | Value |
|-----------|-------|
| URL | https://www.omim.org/ |
| API | Registration required |
| License | Johns Hopkins University trademark |
| Focus | Mendelian disorders and genes |

#### GWAS Catalog
| Attribute | Value |
|-----------|-------|
| URL | https://www.ebi.ac.uk/gwas/ |
| API | https://www.ebi.ac.uk/gwas/rest/api |
| Summary Stats API | https://www.ebi.ac.uk/gwas/summary-statistics/api |
| FTP | http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/ |
| Formats | TSV, OWL/RDF |
| License | Open access |

#### DisGeNET
| Attribute | Value |
|-----------|-------|
| URL | https://disgenet.com/ |
| Content | Gene-disease associations |
| Access | Registration for full access |

#### HGNC (Gene Nomenclature)
| Attribute | Value |
|-----------|-------|
| URL | https://www.genenames.org/ |
| Downloads | Google Storage Bucket |
| Formats | TSV, JSON, OWL |
| Content | Approved human gene symbols |

---

## Integration Recommendations

### Priority Order for Gene Platform

| Priority | Source | Rationale |
|----------|--------|-----------|
| 1 | HPO | Foundation for phenotype annotation |
| 2 | Orphanet/ORDO | Comprehensive rare disease classification |
| 3 | ClinGen | Clinical variant interpretation |
| 4 | MONDO | Unified disease ontology |
| 5 | Monarch KG | Cross-species integration |
| 6 | DECIPHER | Clinical CNV/variant data |
| 7 | BioBank Japan GWAS | Population-specific summary statistics |

### Recommended Integration Strategy

**1. Ontology Layer**
- Load HPO, ORDO, MONDO as base ontologies
- Map terms across ontologies using equivalence axioms
- Support multiple languages via HPO translations

**2. Gene-Disease Associations**
- Primary: ClinGen gene-disease validity
- Secondary: Orphanet gene associations
- Cross-reference: OMIM, MONDO

**3. Variant Data**
- Clinical: ClinGen variant pathogenicity, ClinVar
- Research: DECIPHER (with appropriate agreements)
- Population: GWAS Catalog summary statistics

**4. Biobank Data**
- UK Biobank: Summary statistics (publicly available)
- BioBank Japan: GWAS summary statistics via PheWeb
- Note: Individual-level data requires formal applications

**5. Knowledge Graph**
- Leverage Monarch KG for cross-species phenotype mapping
- Use for variant prioritization (Exomiser)
- Enable semantic similarity calculations

### License Compatibility Matrix

| Source | License | Commercial Use | Attribution |
|--------|---------|----------------|-------------|
| HPO | Custom (open) | Check terms | Required |
| ORDO | CC-BY-4.0 | Yes | Required |
| MONDO | CC-BY-4.0 | Yes | Required |
| ClinGen | Open | Yes | Required |
| ClinVar | Public Domain | Yes | Recommended |
| DECIPHER | Agreement | Contact | Required |
| UK Biobank | Restricted | No | N/A |
| BBJ (GWAS) | Open | Check terms | Required |
| GWAS Catalog | Open | Yes | Recommended |

### Data Refresh Recommendations

| Source | Update Frequency | Method |
|--------|------------------|--------|
| HPO | ~Monthly | GitHub releases |
| ORDO | Bi-annual | Orphadata downloads |
| MONDO | Monthly | GitHub releases |
| ClinGen | Daily | API/Downloads |
| ClinVar | Daily | FTP |
| GWAS Catalog | Weekly | API/FTP |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `Rare disease` | Disease affecting fewer than 200,000 people in the US or 1 in 2,000 in EU | Orphanet 6,528 rare diseases |
| `Ontology` | Structured vocabulary with hierarchical relationships for data standardization | HPO, MONDO, ORDO |
| `OWL` | Web Ontology Language - semantic web language for representing ontologies | HPO hp.owl file |
| `OBO` | Open Biomedical Ontologies format - simpler ontology format | HPO hp.obo file |
| `Equivalence axiom` | Formal statement that two ontology terms represent the same concept | MONDO precise equivalences |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `HPO` | Human Phenotype Ontology - standardized vocabulary of phenotypic abnormalities | 13,000+ terms |
| `MONDO` | Mondo Disease Ontology - unified disease classification | 25,938 disease concepts |
| `ORDO` | Orphanet Rare Disease Ontology - rare disease classification | Orphanet vocabulary |
| `CNV` | Copy Number Variant - large-scale genomic duplication or deletion | DECIPHER CNV data |
| `Pathogenicity` | Classification of whether a variant causes disease | ACMG 5-tier classification |
| `Biobank` | Repository of biological samples linked to phenotypic and genetic data | UK Biobank, BioBank Japan |
| `Phenotype` | Observable characteristics resulting from genotype and environment | HPO-standardized terms |
| `Dosage sensitivity` | Gene's tolerance to increased or decreased copy number | ClinGen haploinsufficiency |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HPO | Human Phenotype Ontology | 13,000+ phenotype terms |
| MONDO | Mondo Disease Ontology | Unified disease ontology |
| ORDO | Orphanet Rare Disease Ontology | Part of Orphanet |
| DECIPHER | DatabasE of genomiC varIation and Phenotype in Humans using Ensembl Resources | 51K+ patients |
| ClinGen | Clinical Genome Resource | NIH-funded variant curation |
| ACMG | American College of Medical Genetics | Variant classification guidelines |
| AMP | Association for Molecular Pathology | Co-authors ACMG guidelines |
| GWAS | Genome-Wide Association Study | UK Biobank, BBJ summary stats |
| BBJ | BioBank Japan | 270K patients, GWAS stats |
| IRDiRC | International Rare Diseases Research Consortium | Rare disease coordination |
| ELIXIR | European life-science Infrastructure for biological Information | Orphanet is Core Resource |
| GRCh37/38 | Genome Reference Consortium Human Build 37/38 | Reference genome versions |
| SNV | Single Nucleotide Variant | Point mutation |
| HGVS | Human Genome Variation Society | Variant nomenclature standard |
| CC BY | Creative Commons Attribution | Open license |
| CC BY-4.0 | Creative Commons Attribution 4.0 | Orphanet, MONDO license |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [_index.md](./_index.md) | Parent index |
| [primary.md](./../genetics/primary.md) | Overlaps ClinVar, dbSNP |
| [disease.md](./../pathways/disease.md) | Disease ontology overlap |

---

## Storage Estimate Summary

| Source | Estimated Size | Notes |
|--------|----------------|-------|
| Orphanet/ORDO | ~500 MB | Ontology + datasets |
| DECIPHER | ~2 GB | Bulk download |
| ClinGen | ~500 MB | All formats |
| HPO | ~100 MB | All formats |
| MONDO | ~200 MB | All formats |
| Monarch KG | ~5 GB | Full knowledge graph |
| BioBank Japan GWAS | ~50 GB | Summary statistics |
| **Total (Core)** | **~8.3 GB** | Excluding UK Biobank |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial document migrated from research.old/data-sources-rare-diseases.md |
