# Rare Disease, Phenotype, and Biobank Database Research

> Research compiled for Gene Platform data integration
> Date: 2026-01-19

## Table of Contents

1. [Orphanet (Rare Diseases)](#1-orphanet-rare-diseases)
2. [DECIPHER Database](#2-decipher-database)
3. [ClinGen (Clinical Genome Resource)](#3-clingen-clinical-genome-resource)
4. [UK Biobank](#4-uk-biobank)
5. [BioBank Japan](#5-biobank-japan)
6. [Phenotype Databases (HPO, Monarch)](#6-phenotype-databases)
7. [Related Resources](#7-related-resources)
8. [Integration Recommendations](#8-integration-recommendations)

---

## 1. Orphanet (Rare Diseases)

### Overview

Orphanet is the reference portal for information on rare diseases and orphan drugs, serving as an ELIXIR Core Data Resource and Global Core Biodata Resource.

### URLs

| Resource | URL |
|----------|-----|
| Main Portal | https://www.orpha.net/ |
| Orphadata (Downloads) | https://sciences.orphadata.com/ |
| ORDO Ontology | https://sciences.orphadata.com/ordo/ |
| API Documentation | https://api.orphadata.com/ |

### Content Statistics

| Metric | Count |
|--------|-------|
| Rare Diseases Catalogued | 6,528 |
| Unique Rare Diseases | ~6,172 |
| Genes Linked | 4,512 |
| Expert Centres | 8,722 |
| Diagnostic Tests | 36,595 |
| Healthcare Professionals | 30,456 |

### Data Types Available

1. **Disease Classifications**: Hierarchical organization of rare diseases
2. **Gene-Disease Associations**: Genes linked to rare conditions
3. **Epidemiological Data**: Prevalence and incidence figures
4. **Expert Centres**: Directory of specialized healthcare centers
5. **Orphan Drugs**: Information on rare disease treatments
6. **Diagnostic Tests**: Available testing resources

### Download Options

**Orphadata Science Datasets:**

| Dataset | Format | Size | Description |
|---------|--------|------|-------------|
| ORDO (Orphanet Rare Disease Ontology) | OWL | ~44-50 MB | Structured vocabulary for rare diseases |
| Disease Classifications | XML, JSON | Variable | Hierarchical disease organization |
| Gene-Disease Associations | XML | Variable | Gene linkages to rare diseases |
| HPO-ORDO Module (HOOM) | OWL | Variable | Phenotype mappings |
| Orphapackets | JSON | Variable | Structured disease packets |

**ORDO Details:**
- Current Version: 4.8 (December 2025)
- Languages: English, Czech, Dutch, French, German, Italian, Spanish, Polish
- Release Schedule: Bi-annual (July and December)

### API Access

- **Swagger Documentation**: https://api.orphadata.com/
- **Interactive Documentation**: Available via OpenAPI specification
- **Try-It-Out**: Disabled in documentation interface
- **Authentication**: Check official documentation for requirements

### License

**CC-BY-4.0 (Creative Commons Attribution 4.0)**
- Free to use with proper attribution
- Commercial and non-commercial use permitted
- Modifications allowed with attribution

### Integration Notes

- ORDO integrates with: OMIM, MONDO, UniProtKB, HGNC, Ensembl, Reactome, IUPHAR
- Recognized by: ELIXIR, Global Biodata Coalition, IRDiRC
- Disease classifications follow decreasing extension: groups > disorders > sub-types

---

## 2. DECIPHER Database

### Overview

DECIPHER (DatabasE of genomiC varIation and Phenotype in Humans using Ensembl Resources) is a clinical genomics resource mapping chromosomal imbalances and phenotypes, maintained at EMBL-EBI.

### URLs

| Resource | URL |
|----------|-----|
| Main Portal | https://www.deciphergenomics.org/ |
| Downloads | https://www.deciphergenomics.org/about/downloads |
| Data Files | https://www.deciphergenomics.org/about/downloads/data |
| Contact | contact@deciphergenomics.org |

### Content Statistics

| Metric | Count |
|--------|-------|
| Open-Access Patients | 51,894 |
| Contributing Centers | 250+ projects |
| Total Cases Uploaded | 19,000+ |
| Consented Patients | ~9,000+ |

### Data Types Available

1. **Copy Number Variants (CNVs)**: Microdeletions and duplications
2. **Sequence Variants**: SNVs, insertions, deletions, indels
3. **Phenotype Data**: Linked to HPO terms
4. **Gene Information**: Including GeneReviews annotations
5. **CNV Syndromes**: Characterized structural variations
6. **Patient Cases**: De-identified clinical genomic data

### Search Capabilities

- Genomic coordinates (GRCh38 and GRCh37)
- Gene symbols and transcripts
- HPO phenotype terms and identifiers
- Variant nomenclature (HGVSg, HGVSc, HGVSp)
- Patient IDs and cytogenetic bands

### Download Options

**Bulk Data Access:**
- Encrypted downloads for academic researchers
- Requires formal Data Access Agreement
- Updated frequently
- Available for: analytical method development, polymorphism pattern studies, critical interval mapping

**API Access:**
- Deposition API for programmatic data upload
- API Client Keys available
- Synchronization capabilities with external systems (e.g., EHRs)

### Access Tiers

| Tier | Visibility | Requirements |
|------|------------|--------------|
| Center-Only | Single center | Default upload |
| Open Access | Public/World | Explicit patient consent |
| Data Display Agreement | External genome browsers | Signed agreement |
| Bulk Research Access | Academic researchers | Data Access Agreement |

### License & Terms

- **Academic Use**: Free with Data Access Agreement
- **Display**: Requires Data Display Agreement
- **Commercial Use**: Contact DECIPHER team
- **Ethical Framework**: Available as downloadable PDF

### Integration Notes

- Displays in: Ensembl Genome Browser, UCSC Genome Browser
- Located at: EMBL-EBI, Wellcome Genome Campus, UK
- Moved to EMBL-EBI in July 2023

---

## 3. ClinGen (Clinical Genome Resource)

### Overview

ClinGen is an NIH-funded resource building a central resource defining the clinical relevance of genes and variants for precision medicine and research.

### URLs

| Resource | URL |
|----------|-----|
| Main Portal | https://clinicalgenome.org/ |
| Knowledge Base | https://search.clinicalgenome.org/ |
| Downloads & APIs | https://search.clinicalgenome.org/kb/downloads |
| Statistics | https://search.clinicalgenome.org/kb/reports/stats |

### Content Statistics (September 2025)

| Metric | Count |
|--------|-------|
| Gene-Disease Validity Curations | 3,289 |
| Genes Curated | 2,719 |
| Variant Pathogenicity Curations | 11,133 |
| Active Contributors | 2,800+ |
| Countries Represented | 74+ |
| Active GCEPs | 56 |

### Data Types Available

1. **Gene-Disease Validity**: Strength of evidence for gene-disease relationships
2. **Variant Pathogenicity**: ACMG/AMP variant classifications
3. **Dosage Sensitivity**: Gene copy number effects (1,557 genes)
4. **Clinical Actionability**: Clinical intervention evaluations (447 gene-condition pairs)
5. **Somatic Cancer Variants**: Cancer-specific annotations
6. **ACMG SF Curations**: Secondary findings recommendations

### Download Formats

| Format | Use Case |
|--------|----------|
| CSV | Spreadsheet analysis |
| TSV | Tab-delimited processing |
| BED | Genomic coordinate visualization |
| AED | ChAS-compatible array analysis |
| JSON (nested/flat) | Programmatic access |

### API Endpoints

| API | Description |
|-----|-------------|
| Evidence Repository API | Variant pathogenicity data |
| Allele Registry API | Variant data retrieval |
| Criteria Specification Registry API | Specification data |
| Linked Data Hub API | Aggregated variant annotations |
| Clinical Actionability APIs | Context-specific JSON |

### FTP Access

- Historical downloads retained for 60 days
- Dosage sensitivity files for GRCh37 and GRCh38
- Multiple genome build support

### License

- Open access for research and clinical use
- Citation required (see "Citing ClinGen & Terms of Use")
- Disclaimer: Not intended for direct diagnostic use without professional review

### Integration Notes

- Works with: ClinVar, GenomeConnect
- Expert panels use ACMG/AMP guidelines with gene-specific specifications
- Statistics updated daily

---

## 4. UK Biobank

### Overview

UK Biobank is a large-scale biomedical database containing genetic, imaging, and health data from 500,000 UK participants, recruited between 2006-2010.

### URLs

| Resource | URL |
|----------|-----|
| Main Portal | https://www.ukbiobank.ac.uk/ |
| Research Analysis Platform | UKB-RAP (cloud-based) |
| Data Releases | https://www.ukbiobank.ac.uk/about-our-data/data-releases/ |
| Access Information | https://www.ukbiobank.ac.uk/about-us/how-we-work/access-to-uk-biobank-data/ |

### Content Statistics

| Metric | Value |
|--------|-------|
| Participants | 500,000 |
| Original Age Range | 40-69 years |
| Recruitment Period | 2006-2010 |
| Data Variables | 10,000+ |
| Biological Samples | 15+ million |
| Database Size | 30+ petabytes |
| Registered Researchers | 30,000+ |
| Countries Represented | 90+ |
| Peer-Reviewed Publications | 9,000+ |

### Data Types Available

1. **Genetic Data**: Whole genome sequencing, genotyping arrays
2. **Metabolomic Data**: ~250 blood metabolites (all 500K participants as of 2025)
3. **Imaging Data**: Brain, cardiac, body (targeting 100,000 by end of 2025)
4. **Biological Samples**: Blood, urine, saliva
5. **Questionnaires**: Health, lifestyle, diet, environment
6. **Physical Measurements**: Body composition, bone density
7. **Activity Tracking**: Accelerometer data
8. **Follow-up Data**: Health outcomes, hospital records, death registry

### Access Model

**Primary Access**: UK Biobank Research Analysis Platform (UKB-RAP)
- Secure cloud-based platform
- Default access method since 2024
- Democratizes access for early-career researchers and lower-income countries

**Access Requirements:**
1. Researcher approval process
2. Institutional affiliation verification
3. Access committee approval
4. Data usage agreement

### Public Data Availability

| Type | Availability |
|------|--------------|
| Summary Statistics | Publicly accessible |
| Individual-Level Data | Controlled access only |
| Published Results | Available via publications |

### Recent Updates (2025)

- Metabolomic data expanded to all 500,000 participants
- 55,000 3D heart models added
- Additional 30,000 bone density measures
- Brain structure longitudinal data
- Imaging project approaching 100,000 participants

### License & Terms

- **Research Use Only**: Applications from insurance companies no longer accepted (Jan 2025)
- Charity registered in England/Wales (1101332) and Scotland (SC039230)
- Requires formal application and approval

### Integration Notes

- No public API for individual-level data
- Summary statistics may be available publicly
- 22,000 researchers actively using data worldwide

---

## 5. BioBank Japan

### Overview

BioBank Japan (BBJ) is a disease-oriented prospective biobank managed by the Institute of Medical Science, University of Tokyo, containing samples from approximately 270,000 Japanese patients across 51 diseases.

### URLs

| Resource | URL |
|----------|-----|
| Main Portal | https://biobankjp.org/en/ |
| PheWeb (GWAS Results) | https://pheweb.jp/ |
| JENGER Portal | http://jenger.riken.jp/en/ |
| Sample Search | BBJ Online Biological Samples Search System |

### Content Statistics

| Metric | Value |
|--------|-------|
| Total Patients | ~270,000 |
| Target Diseases | 51 |
| DNA Samples | 800,000 tubes |
| Serum Samples | 1,700,000 tubes |
| Contributing Institutions | 12 medical centers |
| WGS Samples | 16,000 patients |
| SNP Data | 270,000 patients |

### Data Types Available

1. **DNA Samples**: Stored biological material
2. **Serum Samples**: 200,000 patients
3. **Whole Genome Sequencing**: 16,000 patients
4. **SNP Genotyping Data**: 270,000 patients
5. **Clinical Information**: Linked to samples
6. **GWAS Summary Statistics**: Publicly available

### Public GWAS Data Access (PheWeb.jp)

**Available Since**: March 2021

| Metric | Value |
|--------|-------|
| Initial GWAS Results | 220 |
| Study Population | ~260,000 participants |
| Ancestry | Japanese |

**Data Categories:**
- Metabolic disorders (Type 2 diabetes)
- Autoimmune conditions (rheumatoid arthritis, psoriasis)
- Cancer types (breast cancer)
- Infectious disease responses (COVID-19)
- Neurological conditions
- Dietary phenotypes

**Download Options:**
- Full GWAS summary statistics (public)
- Fine-mapping results
- Raw eQTL data (via NBDC website)

### Controlled Data Access

**Requirements:**
- Formal application to BBJ Sample and Data Access Committee
- Rigorous screening process
- Compliance with Japanese regulations
- Security guidelines and contractual agreements

### Certifications

- ISO 9001:2015 (Quality Management)
- ISO/IEC 27001:2013 (Information Security Management)

### Key Publications

| Study | Sample Size | Findings |
|-------|-------------|----------|
| Quantitative Traits GWAS | 162,255 | 1,407 trait-associated loci (679 novel) |
| 42 Diseases GWAS | 212,453 | 320 independent signals, 25 novel loci |

### License & Terms

- **GWAS Summary Statistics**: Publicly available
- **Individual-Level Data**: Controlled access
- **Sample Requests**: Application required

### Integration Notes

- NBDC Human Database hosts controlled data
- April 2024: URL changes for NBDC-hosted data
- Primarily serves Japanese population research

---

## 6. Phenotype Databases

### 6.1 Human Phenotype Ontology (HPO)

#### Overview

HPO is a standardized vocabulary of phenotypic abnormalities observed in human disease, maintained by the Monarch Initiative and Jackson Laboratory.

#### URLs

| Resource | URL |
|----------|-----|
| Main Portal | https://hpo.jax.org/ |
| GitHub Repository | https://github.com/obophenotype/human-phenotype-ontology |
| OBO Foundry | http://purl.obolibrary.org/obo/hp.owl |
| EBI OLS | https://www.ebi.ac.uk/ols/ontologies/hp |

#### Content Statistics (2024)

| Metric | Value |
|--------|-------|
| Total Terms | 13,000+ |
| Annotations to Diseases | 156,000+ |
| New Terms Added (recent) | 2,239 |
| New Annotations Added | 49,235 |
| GitHub Commits | 8,308 |
| Releases | 61 |
| Contributors | 24+ |

#### Download Formats

| Format | URL |
|--------|-----|
| OWL | http://purl.obolibrary.org/obo/hp.owl |
| OBO | http://purl.obolibrary.org/obo/hp.obo |
| JSON | Available via GitHub releases |

#### Languages Supported

English, German, Spanish, French, Italian, Dutch, Portuguese, Turkish, Japanese, Chinese + more

#### License

- Open access with attribution requirements
- See: https://hpo.jax.org/app/license

#### Recent Developments

- Medical Action Ontology (MAxO) for treatment modeling
- GA4GH Phenopacket Schema integration
- EHR integration efforts
- Internationalization with 10+ language translations

---

### 6.2 Monarch Initiative

#### Overview

Monarch Initiative is an integrative data and analytic platform connecting phenotypes to genotypes across species.

#### URLs

| Resource | URL |
|----------|-----|
| Main Portal | https://monarchinitiative.org/ |
| Data Downloads | https://data.monarchinitiative.org/ |
| API (v3) | https://api-v3.monarchinitiative.org/v3/docs |
| GitHub | https://github.com/monarch-initiative |

#### Knowledge Graph Statistics

**Data Sources Integrated**: 33 heterogeneous sources including:
- MGD, ZFIN, WormBase, FlyBase, Xenbase, SGD, PomBase, DictyBase, BGeeDB
- OMIM, Orphanet (human disease catalogs)
- Reactome (pathways)
- STRING (protein-protein interactions)
- HPO, Gene Ontology, PHENIO

#### Primary Data Types

| Type | Description |
|------|-------------|
| Genes | Cross-species gene information |
| Diseases | Integrated disease ontologies |
| Phenotypes | Cross-species phenotype mappings |

#### Available Datasets

| Dataset | Description |
|---------|-------------|
| Monarch KG | Main knowledge graph |
| HPO | Human Phenotype Ontology |
| Exomiser | Variant prioritization resources |
| Phenologs | Phenotype ortholog data |
| Semantic Similarity | Phenotype similarity metrics |
| Mappings | Cross-reference mappings |
| UPHENO2 | Unified Phenotype Ontology |

#### API Access

- **Version**: API v3 (launched October 2023)
- **Documentation**: OpenAPI/Swagger
- **Endpoint**: https://api-v3.monarchinitiative.org/
- **Features**: Concise abstracted endpoints for all entity combinations

#### License

- Open access
- Individual datasets may have specific licenses

#### Integration Notes

- ChatGPT plugin available for phenotypic data queries
- Previous API (Biolink, OWLSim) discontinued March 2024
- R package available: monarchr

---

### 6.3 MONDO Disease Ontology

#### Overview

MONDO is a unified disease ontology providing precise equivalence axioms across disease resources.

#### URLs

| Resource | URL |
|----------|-----|
| Main Portal | https://mondo.monarchinitiative.org/ |
| GitHub | https://github.com/monarch-initiative/mondo |
| Latest Release | https://github.com/monarch-initiative/mondo/releases/latest |

#### Content Statistics

| Metric | Value |
|--------|-------|
| Total Disease Concepts | 25,938 |
| Human Diseases | 22,977 |
| Non-Human Diseases | 2,960 |

#### Download Formats

| Format | File | Description |
|--------|------|-------------|
| OWL | mondo-with-equivalents.owl | With equivalence axioms |
| OBO | mondo.obo | Using xrefs for linking |
| JSON | mondo-with-equivalents.json | Equivalent to OWL |

#### Key Feature

Provides "precise 1:1 equivalence axioms" validated by OWL reasoning, distinguishing it from loose cross-references in other resources.

#### License

**CC-BY-4.0** - Creative Commons Attribution 4.0

---

## 7. Related Resources

### 7.1 ClinVar

| Attribute | Value |
|-----------|-------|
| URL | https://www.ncbi.nlm.nih.gov/clinvar/ |
| FTP | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/ |
| Formats | VCF (GRCh37/38), XML, Tab-delimited |
| License | Public domain (NCBI) |
| Updates | Daily |

### 7.2 OMIM

| Attribute | Value |
|-----------|-------|
| URL | https://www.omim.org/ |
| API | Registration required |
| License | Johns Hopkins University trademark |
| Focus | Mendelian disorders and genes |

### 7.3 GWAS Catalog

| Attribute | Value |
|-----------|-------|
| URL | https://www.ebi.ac.uk/gwas/ |
| API | https://www.ebi.ac.uk/gwas/rest/api |
| Summary Stats API | https://www.ebi.ac.uk/gwas/summary-statistics/api |
| FTP | http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/ |
| Formats | TSV, OWL/RDF |
| License | Open access |

### 7.4 DisGeNET

| Attribute | Value |
|-----------|-------|
| URL | https://disgenet.com/ |
| Content | Gene-disease associations |
| Access | Registration for full access |

### 7.5 HGNC (Gene Nomenclature)

| Attribute | Value |
|-----------|-------|
| URL | https://www.genenames.org/ |
| Downloads | Google Storage Bucket |
| Formats | TSV, JSON, OWL |
| Content | Approved human gene symbols |

---

## 8. Integration Recommendations

### Priority Data Sources for Gene Platform

| Priority | Source | Rationale |
|----------|--------|-----------|
| 1 | HPO | Foundation for phenotype annotation |
| 2 | Orphanet/ORDO | Comprehensive rare disease classification |
| 3 | ClinGen | Clinical variant interpretation |
| 4 | MONDO | Unified disease ontology |
| 5 | Monarch KG | Cross-species integration |
| 6 | DECIPHER | Clinical CNV/variant data |

### Recommended Integration Strategy

1. **Ontology Layer**
   - Load HPO, ORDO, MONDO as base ontologies
   - Map terms across ontologies using equivalence axioms
   - Support multiple languages via HPO translations

2. **Gene-Disease Associations**
   - Primary: ClinGen gene-disease validity
   - Secondary: Orphanet gene associations
   - Cross-reference: OMIM, MONDO

3. **Variant Data**
   - Clinical: ClinGen variant pathogenicity, ClinVar
   - Research: DECIPHER (with appropriate agreements)
   - Population: GWAS Catalog summary statistics

4. **Biobank Data**
   - UK Biobank: Summary statistics (publicly available)
   - BioBank Japan: GWAS summary statistics via PheWeb
   - Note: Individual-level data requires formal applications

5. **Knowledge Graph**
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

## References

### Primary Sources

1. Orphanet - https://www.orpha.net/
2. DECIPHER - https://www.deciphergenomics.org/
3. ClinGen - https://clinicalgenome.org/
4. UK Biobank - https://www.ukbiobank.ac.uk/
5. BioBank Japan - https://biobankjp.org/en/
6. HPO - https://hpo.jax.org/
7. Monarch Initiative - https://monarchinitiative.org/

### Key Publications

1. "The Human Phenotype Ontology in 2024: phenotypes around the world" - Nucleic Acids Research (2024)
2. "The Monarch Initiative in 2024: an analytic platform integrating phenotypes, genes and diseases across species" - Nucleic Acids Research (2024)
3. "The Clinical Genome Resource (ClinGen): Advancing Genomic Knowledge through Global Curation" - Genetics in Medicine (2024)
4. "Estimating cumulative point prevalence of rare diseases: analysis of the Orphanet database" - European Journal of Human Genetics (2020)
5. "DECIPHER: Supporting the interpretation and sharing of rare disease phenotype-linked variant data" - Human Mutation (2022)

---

*Document compiled for Gene Platform integration planning*
*Last updated: 2026-01-19*
