# Pharmaceutical and Pharmacogenomics Data Sources Research

**Research Date:** January 2026
**Purpose:** Comprehensive evaluation of pharmaceutical and pharmacogenomics databases for genetics/health knowledge base integration

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Primary Pharmacogenomics Databases](#primary-pharmacogenomics-databases)
   - [PharmGKB / ClinPGx](#1-pharmgkb--clinpgx)
   - [CPIC Guidelines](#2-cpic-clinical-pharmacogenetics-implementation-consortium)
   - [DPWG Guidelines](#3-dpwg-dutch-pharmacogenetics-working-group)
   - [PharmVar](#4-pharmvar-pharmacogene-variation-consortium)
3. [Comprehensive Drug Databases](#comprehensive-drug-databases)
   - [DrugBank](#5-drugbank)
   - [ChEMBL](#6-chembl)
   - [PubChem](#7-pubchem)
   - [DrugCentral](#8-drugcentral)
4. [Drug-Target Interaction Databases](#drug-target-interaction-databases)
   - [Open Targets](#9-open-targets-platform)
   - [DGIdb](#10-dgidb-drug-gene-interaction-database)
   - [TTD](#11-ttd-therapeutic-target-database)
   - [BindingDB](#12-bindingdb)
   - [Guide to Pharmacology](#13-iupharbps-guide-to-pharmacology-gtopdb)
5. [Regulatory and Clinical Sources](#regulatory-and-clinical-sources)
   - [FDA DailyMed](#14-fda-dailymed)
   - [OpenFDA](#15-openfda)
   - [FDA Pharmacogenomic Biomarkers Table](#16-fda-table-of-pharmacogenomic-biomarkers)
   - [RxNorm](#17-rxnorm)
6. [Specialized Databases](#specialized-databases)
   - [CIViC](#18-civic-clinical-interpretation-of-variants-in-cancer)
   - [PGxDB](#19-pgxdb)
   - [SuperDRUG2](#20-superdrug2)
   - [KEGG Drug](#21-kegg-drug)
   - [UniProt](#22-uniprot)
7. [Comparison Matrix](#comparison-matrix)
8. [Recommendations](#recommendations)

---

## Executive Summary

This document provides a comprehensive analysis of pharmaceutical and pharmacogenomics databases suitable for integration into a genetics/health knowledge base platform. The research covers **23+ databases** spanning pharmacogenomics guidelines, drug information, drug-target interactions, and regulatory data.

### Key Findings

| Category | Top Recommendation | Rationale |
|----------|-------------------|-----------|
| **Pharmacogenomics** | PharmGKB/ClinPGx + CPIC | Gold standard for SNP-drug relationships and clinical dosing |
| **Drug Information** | DrugBank + ChEMBL | Most comprehensive drug data with target information |
| **Drug-Gene Interactions** | DGIdb + Open Targets | Aggregated, open-access interaction data |
| **Regulatory** | OpenFDA + DailyMed | Authoritative FDA drug labeling and adverse events |
| **Binding Affinity** | BindingDB + ChEMBL | Quantitative binding data for drug discovery |

---

## Primary Pharmacogenomics Databases

### 1. PharmGKB / ClinPGx

**URL:** https://www.pharmgkb.org / https://www.clinpgx.org
**Maintainer:** Stanford University (NIH-funded)

#### Content
- Clinical annotations linking variants to drug responses
- Variant annotations from published literature
- Drug-centered pharmacokinetic/pharmacodynamic pathways
- Dosing guidelines (CPIC, DPWG)
- Annotated FDA drug labels
- Gene summaries for pharmacogenes

#### Coverage
- **Genes:** 700+ pharmacogenes
- **Drugs:** 1,000+ with pharmacogenomic relevance
- **Clinical Annotations:** 3,000+
- **Variant Annotations:** 20,000+
- **Pathways:** 150+ drug-centered pathways

#### Access Method
- **Web Interface:** No login required for viewing
- **API:** RESTful API serving JSON (https://api.pharmgkb.org/)
- **Bulk Download:** Requires account and license agreement
- **Swagger Documentation:** https://api.pharmgkb.org/swagger/

#### Data Format
- TSV/CSV spreadsheets (downloads)
- JSON (API)
- BioPax (pathways)

#### Schema/Structure
Key entities and relationships:
```
Gene <---> Variant <---> Clinical Annotation <---> Drug
                              |
                              v
                    Level of Evidence (1A-4)
                              |
                              v
                    Phenotype Categories
```

**Clinical Annotation Fields:**
- Variant (rsID, star allele)
- Drug name
- Phenotype categories
- Level of evidence (1A, 1B, 2A, 2B, 3, 4)
- Genotype-phenotype associations
- Supporting publications

**Levels of Evidence:**
- **1A:** CPIC/medical society guideline or implemented clinically
- **1B:** Preponderance of evidence shows association
- **2A:** Level 2B + Very Important Pharmacogene (VIP)
- **2B:** Moderate evidence, replicated but some non-significant studies
- **3:** Single significant study or lacking clear evidence
- **4:** Case report, non-significant study, or in vitro evidence only

#### Licensing
- **Creative Commons Attribution-ShareAlike 4.0** (CC BY-SA 4.0)
- Free for academic and commercial use with attribution
- Download requires account creation and license agreement

#### Pharmacogenomics Focus
- **SNP-Drug Relationships:** Primary focus with clinical annotations
- **Dosing Recommendations:** Via CPIC/DPWG guidelines
- **Genotype-to-Phenotype:** Standardized metabolizer status assignments

#### Full Content vs Abstract
- **FULL CONTENT:** Complete clinical summaries, variant details, mechanism explanations
- Includes original publication references and curated evidence

---

### 2. CPIC (Clinical Pharmacogenetics Implementation Consortium)

**URL:** https://cpicpgx.org
**Maintainer:** PharmGKB/PGRN consortium

#### Content
- Evidence-based gene/drug clinical practice guidelines
- Allele functionality tables
- Genotype-to-phenotype translation tables
- Dosing recommendations by genotype
- Drug-gene pairs with clinical actionability

#### Coverage
- **Genes:** 34 genes
- **Drugs:** 164 drugs covered
- **Guidelines:** 27 published guidelines
- **Healthcare Institutions Using CPIC:** 128+
- **Commercial Laboratories:** 40+

#### Access Method
- **Web Interface:** Free access, no login
- **API:** Over 80,000 monthly queries
- **Bulk Download:** JSON format available on PharmGKB
- **Integration:** Epic's foundational genomics module

#### Data Format
- PDF guidelines
- Excel/TSV (allele definition tables, dosing tables)
- JSON (downloadable via PharmGKB)

#### Schema/Structure
```
Gene --> Allele Definitions --> Diplotype --> Phenotype --> Dosing Recommendation
                                                   |
                                                   v
                                        Strength of Recommendation (Strong/Moderate/Optional)
```

**Key Tables:**
1. Allele Definition Table (variant positions, defining variants)
2. Allele Functionality Table (function assignments)
3. Diplotype-Phenotype Table
4. Dosing Recommendation Table

#### Licensing
- **Creative Commons Public Domain (CC0)**
- Completely free for any use

#### Pharmacogenomics Focus
- **SNP-Based Dosing:** PRIMARY PURPOSE - detailed dosing by genotype
- **Clinical Implementation:** Ready for EHR integration
- **Actionable Recommendations:** Each guideline has specific dosing adjustments

#### Full Content vs Abstract
- **FULL CONTENT:** Complete guidelines with evidence summaries, dosing tables, implementation guidance

---

### 3. DPWG (Dutch Pharmacogenetics Working Group)

**URL:** https://www.knmp.nl (Royal Dutch Pharmacists Association)
**Also available via:** PharmGKB, ClinPGx
**Maintainer:** Royal Dutch Pharmacists Association (KNMP)

#### Content
- Gene-drug interaction guidelines
- Therapeutic recommendations based on genotype
- Risk assessments for gene-drug combinations

#### Coverage
- Multiple guidelines for major pharmacogenes including:
  - CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP1A2
  - DPYD, TPMT, NUDT15
  - UGT1A1
  - HLA-A, HLA-B

#### Recent Guidelines (2024)
- CYP2D6, CYP3A4, CYP1A2 and antipsychotics (March 2024)
- CYP2D6, CYP2C19 and non-SSRI/non-TCA antidepressants (November 2024)
- CYP2C9, HLA-A, HLA-B with anti-epileptic drugs (August 2024)

#### Access Method
- Published in European Journal of Human Genetics
- Available through PharmGKB
- Implemented in Dutch pharmacy/physician systems

#### Data Format
- Published guidelines (PDF)
- Integrated into PharmGKB annotations

#### Licensing
- Published openly in peer-reviewed journals
- Free for clinical implementation

#### Pharmacogenomics Focus
- **SNP-Based Dosing:** European perspective on pharmacogenomics implementation
- **Clinical Implementation:** Integrated into Dutch healthcare IT systems

---

### 4. PharmVar (Pharmacogene Variation Consortium)

**URL:** https://www.pharmvar.org
**Maintainer:** PharmVar Consortium

#### Content
- Official pharmacogene allele nomenclature
- Star allele definitions for pharmacogenes
- Variant-to-allele mappings
- Haplotype definitions

#### Coverage
- All major CYP450 enzymes (CYP2D6, CYP2C19, CYP2C9, CYP3A4, etc.)
- Other pharmacogenes (UGT1A1, DPYD, TPMT, etc.)
- Monthly releases with version numbering

#### Access Method
- **Web Interface:** www.PharmVar.org
- **Downloads:** Gene-specific allele definition files
- **Submissions:** Via submission form

#### Data Format
- Allele definition tables
- Read Me and Change Log documents
- Downloadable files per gene

#### Schema/Structure
```
Gene --> Star Allele --> Defining Variants (rsIDs) --> Function Assignment
                              |
                              v
                    Reference Sequence Positions
```

#### Licensing
- Open access for research and clinical use
- Works closely with PharmGKB and CPIC

#### Pharmacogenomics Focus
- **Allele Nomenclature:** Standardized naming system for pharmacogene variants
- **Essential for:** Interpreting pharmacogenomic test results

#### Recent Updates (2024-2025)
- CYP1A2 completed transition from legacy nomenclature (December 2024)
- Monthly database releases (currently at version 6.2+)

---

## Comprehensive Drug Databases

### 5. DrugBank

**URL:** https://go.drugbank.com
**Maintainer:** DrugBank (OMx Personal Health Analytics)

#### Content
- Comprehensive drug information (approved, experimental, investigational)
- Drug targets, enzymes, carriers, transporters
- Drug-drug interactions
- Chemical structures and properties
- Pharmacokinetic data
- Drug metabolism pathways

#### Coverage (Version 6.0, 2024)
- **Approved Drugs:** 2,700+
- **Experimental Drugs:** 6,700+
- **Total Drug Entries:** 15,000+
- **Drug Targets:** 5,000+
- **Drug-Drug Interactions:** 2.5 million+
- **SNP-Drug Associations:** Included

#### Access Method
- **Web Interface:** Free browsing
- **API:** REST API (subscription required for full access)
- **Bulk Download:** Multiple licensing tiers

#### Data Format
- XML (comprehensive database export)
- CSV (specific datasets)
- SDF (chemical structures)
- FASTA (protein sequences)

#### Schema/Structure
```
Drug --> Targets (proteins) --> Gene
  |          |
  |          v
  |     UniProt ID, Actions, Known Action
  |
  +--> Drug Interactions
  |
  +--> Metabolizing Enzymes
  |
  +--> Pharmacokinetics (ADMET)
```

**Key Fields:**
- DrugBank ID, Name, Type
- Description, Indication, Pharmacodynamics
- Mechanism of Action
- Targets (with actions: inhibitor, agonist, antagonist, etc.)
- Enzymes, Carriers, Transporters
- Drug Interactions
- Food Interactions
- SNP Effects (pharmacogenomics)

#### Licensing

| License Type | Access | Cost |
|--------------|--------|------|
| **Open Data (CC0)** | Vocabulary, Structures | Free |
| **Academic (CC BY-NC 4.0)** | Full XML, most data | Free for non-commercial |
| **Commercial** | Full access + support | Contact for pricing |

#### Pharmacogenomics Focus
- **SNP Effects:** Documented for many drugs
- **Metabolizing Enzymes:** CYP450 and other enzymes linked
- **Drug-Drug Interactions:** Metabolic pathway-based

#### Full Content vs Abstract
- **FULL CONTENT:** Complete drug monographs, target details, interaction mechanisms

---

### 6. ChEMBL

**URL:** https://www.ebi.ac.uk/chembl
**Maintainer:** EMBL-EBI

#### Content
- Bioactivity data from medicinal chemistry literature
- Drug-target binding data
- Assay results
- Compound structures
- ADMET properties
- Drug mechanism of action

#### Coverage (ChEMBL 36, July 2025)
- **Compounds:** 2.8 million distinct
- **Bioactivities:** 21+ million
- **Assays:** 1.6+ million
- **Targets:** 17,803
- **Documents:** 95,000+ publications

#### Access Method
- **Web Interface:** Free, no login
- **REST API:** Full programmatic access
- **Bulk Download:** FTP/HTTPS

#### Data Format
- SQLite, MySQL, PostgreSQL (database dumps)
- SDF (structures)
- FASTA (target sequences)
- RDF/Turtle (linked data)
- JSON, XML, YAML (API)

#### Schema/Structure
```
Document --> Assay --> Activity --> Compound
               |                        |
               v                        v
            Target                  CHEMBL ID
               |                   Structure (SMILES, InChI)
               v
          UniProt ID
```

**Key Tables:**
- molecule_dictionary (compounds)
- assays (experimental setup)
- activities (measured values: IC50, Ki, EC50, etc.)
- target_dictionary (proteins, organisms)
- mechanism (drug mechanisms of action)

#### Licensing
- **Open Access** - Free for all uses
- Apache 2.0 license for software
- Data freely available

#### Pharmacogenomics Focus
- **Drug-Target Binding:** Quantitative affinity data
- **Mechanism of Action:** Detailed target information
- **Drug Discovery:** Pre-clinical and clinical compound data

#### Full Content vs Abstract
- **FULL CONTENT:** Complete bioactivity measurements, assay details, compound properties

---

### 7. PubChem

**URL:** https://pubchem.ncbi.nlm.nih.gov
**Maintainer:** NIH/NCBI

#### Content
- Chemical structures
- Biological activities from bioassays
- Patent data
- Literature references
- Physical/chemical properties
- Safety data

#### Coverage
- **Compounds:** 118+ million
- **Substances:** 326+ million
- **BioAssays:** 1.5+ million
- **Patent Compounds:** 19+ million

#### Access Method
- **Web Interface:** Free, no login
- **PUG REST API:** Programmatic access
- **Bulk Download:** FTP site

#### Data Format
- SDF (structures)
- JSON, XML (API responses)
- CSV (tabular data)
- ASN.1 (legacy)

#### Schema/Structure
```
Compound (CID) <--> Substance (SID) <--> Source
       |
       v
   BioAssay (AID) --> Activity Data
       |
       v
   Target (Gene/Protein)
```

#### Licensing
- **Public Domain** - No restrictions
- Free for any use

#### Pharmacogenomics Focus
- **Limited direct PGx:** More focused on chemical/bioassay data
- **Cross-references:** Links to PharmGKB, DrugBank

#### Full Content vs Abstract
- **FULL CONTENT:** Complete chemical data, bioassay results

---

### 8. DrugCentral

**URL:** https://drugcentral.org
**Maintainer:** University of New Mexico, Translational Informatics Division

#### Content
- FDA, EMA, PMDA approved drugs
- Drug structures and properties
- Mechanisms of action
- Drug indications
- Pharmacological actions
- Adverse events (FAERS analysis)
- Veterinary drugs (new in 2023)

#### Coverage
- Active pharmaceutical ingredients with approval status
- Drug-target interactions with bioactivity
- Age-stratified adverse event data (pediatric, geriatric)

#### Access Method
- **Web Interface:** Free, open access
- **API:** PostgreSQL database access
- **Bulk Download:** PostgreSQL dump, Docker container

**Database Connection (public):**
- Host: unmtid-dbs.net
- Port: 5433
- Database: drugcentral
- User: drugman

#### Data Format
- PostgreSQL database
- API (JSON)
- Docker container

#### Schema/Structure
```
Drug --> Target --> Gene
  |          |
  |          v
  |     Action Type, Bioactivity
  |
  +--> Indications (diseases)
  |
  +--> Adverse Events (FAERS)
```

#### Licensing
- **Open Access** - Free for all uses

#### Pharmacogenomics Focus
- **Drug-Target Interactions:** With bioactivity data
- **FAERS Integration:** Adverse event analysis

#### Full Content vs Abstract
- **FULL CONTENT:** Complete drug profiles with targets and adverse events

---

## Drug-Target Interaction Databases

### 9. Open Targets Platform

**URL:** https://platform.opentargets.org
**Maintainer:** Open Targets consortium (Wellcome Sanger, EMBL-EBI, GSK, etc.)

#### Content
- Target-disease associations
- Drug-target relationships
- Genetic associations
- Pathway information
- Literature mining results

#### Coverage (2025)
- Evidence from multiple data sources aggregated
- Genetic associations from GWAS, functional genomics
- Drug data from ChEMBL

#### Access Method
- **Web Interface:** Free, interactive visualizations
- **GraphQL API:** Flexible queries
- **Bulk Download:** FTP, Google Cloud, AWS, Azure

#### Data Format
- Parquet files (bulk downloads)
- JSON (API)
- TSV exports from web interface

#### Schema/Structure
```
Target (gene) <--> Disease <--> Drug
       |               |
       v               v
  Evidence Score   Evidence Sources
```

#### Licensing
- **Open Access** (CC0 for data)
- Open source tools

#### Pharmacogenomics Focus
- **Genetic Evidence:** GWAS associations for drug targets
- **Drug Tractability:** Assessment of target druggability

---

### 10. DGIdb (Drug Gene Interaction Database)

**URL:** https://dgidb.org
**Maintainer:** Washington University School of Medicine

#### Content
- Drug-gene interactions from 40+ sources
- Druggable gene categories
- Drug groupings

#### Coverage (v5.0, 2024)
- **Genes:** 10,000+
- **Drugs:** 20,000+
- **Interactions:** 70,000+

#### Access Method
- **Web Interface:** Search and browse
- **GraphQL API:** Customizable queries
- **Bulk Download:** TSV files

#### Data Format
- TSV (downloads)
- JSON (API)
- PostgreSQL backend

#### Schema/Structure
```
Gene <--> Interaction <--> Drug
              |
              v
        Interaction Type
        (inhibitor, agonist, etc.)
              |
              v
         Source Database
```

#### Licensing
- **Open Access** - Free for all uses

#### Pharmacogenomics Focus
- **Drug-Gene Interactions:** Aggregated from multiple sources
- **Druggable Genome:** Categorization of genes by druggability

#### Full Content vs Abstract
- **AGGREGATED CONTENT:** Combines data from DrugBank, PharmGKB, ChEMBL, etc.

---

### 11. TTD (Therapeutic Target Database)

**URL:** https://idrblab.net/ttd
**Maintainer:** Zhejiang University

#### Content
- Therapeutic protein and nucleic acid targets
- Target druggability information
- Drugs/ligands for targets
- Pathway information
- Disease associations

#### Coverage (2024)
- **Successful Targets:** 426
- **Clinical Trial Targets:** 1,014
- **Preclinical Targets:** 212
- **Literature Targets:** 1,479

#### Access Method
- **Web Interface:** Free, no login
- **Bulk Download:** https://idrblab.net/ttd/full-data-download

#### Data Format
- Downloadable data files
- MySQL backend

#### Schema/Structure
```
Target --> Drugs/Ligands
   |
   v
Disease --> Pathway
   |
   v
Druggability Metrics
```

#### Licensing
- **Free Access** - No login required

#### Pharmacogenomics Focus
- **Target Druggability:** Molecular interactions, system features, expression variations

---

### 12. BindingDB

**URL:** https://www.bindingdb.org
**Maintainer:** UC San Diego

#### Content
- Experimentally measured binding affinities
- Protein-ligand interactions
- Quantitative data (IC50, Ki, Kd, EC50)

#### Coverage (2024)
- **Binding Measurements:** 2.9 million
- **Compounds:** 1.3 million
- **Targets:** Thousands of proteins

#### Access Method
- **Web Interface:** Search and browse
- **REST API:** JSON/XML responses
- **Bulk Download:** SDF, TSV files
- **KNIME Workflows:** For data retrieval

#### Data Format
- SDF (structures)
- TSV (data tables)
- FASTA (protein sequences)
- JSON/XML (API)

#### API Examples
```
# By UniProt IDs:
http://bindingdb.org/rest/getLigandsByUniprots?uniprot=P00176,P00183&cutoff=10000&response=application/json

# By SMILES structure:
http://bindingdb.org/rest/getLigandsBySmiles?smiles=...&similarity=0.8
```

#### Schema/Structure
```
Compound (SMILES) <--> Binding Data <--> Target (UniProt)
                           |
                           v
                    Affinity Value (IC50, Ki, Kd)
                    Assay Type
                    Publication
```

#### Licensing
- **Open Access** - FAIR principles compliant

#### Pharmacogenomics Focus
- **Quantitative Binding:** Precise affinity measurements
- **Drug Discovery:** Supports computational drug design

---

### 13. IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb)

**URL:** https://www.guidetopharmacology.org
**Maintainer:** IUPHAR (International Union of Basic and Clinical Pharmacology) / BPS

#### Content
- Expert-curated pharmacological targets
- Ligand information
- Quantitative interaction data
- Target families (GPCRs, ion channels, enzymes, etc.)

#### Coverage (Version 2025.4)
- **Human Targets:** 3,097
- **Targets with Ligand Interactions:** 1,758
- **Ligands:** 13,131
- **Ligands with Target Interactions:** 9,485

#### Access Method
- **Web Interface:** Free access
- **REST API:** JSON format
- **Bulk Download:** CSV/TSV, PostgreSQL dump, RDF

#### Data Format
- CSV/TSV files
- PostgreSQL database
- RDF/N3 (linked data)
- JSON (API)

#### Licensing
- **Open Access** - Free for all uses

#### Pharmacogenomics Focus
- **Receptor Pharmacology:** Detailed target characterization
- **Ligand-Target Interactions:** Quantitative affinity data

---

## Regulatory and Clinical Sources

### 14. FDA DailyMed

**URL:** https://dailymed.nlm.nih.gov
**Maintainer:** NIH/NLM

#### Content
- FDA-approved drug labels (package inserts)
- Structured Product Labeling (SPL) files
- Over-the-counter drug labels
- Animal drug labels

#### Coverage
- **Drug Labels:** 155,000+
- Updated daily with new submissions

#### Access Method
- **Web Interface:** Search and browse
- **Bulk Download:** Daily, weekly, monthly updates
- **FTP:** All SPL files

#### Data Format
- XML (SPL format)
- ZIP archives

#### Licensing
- **Public Domain** - No restrictions

#### Pharmacogenomics Focus
- **Drug Labels:** Include pharmacogenomics sections when applicable
- **Warnings:** Genetic testing requirements

---

### 15. OpenFDA

**URL:** https://open.fda.gov
**Maintainer:** FDA

#### Content
- Drug adverse events (FAERS)
- Drug labels
- Drug enforcement reports (recalls)
- NDC directory
- Drug shortages

#### Coverage
- **API Calls:** 200+ million total
- **Adverse Events:** Millions of reports
- **Drug Labels:** From DailyMed

#### Access Method
- **REST API:** Free, rate-limited
- **Bulk Download:** Zipped JSON

#### Data Format
- JSON (API and downloads)

#### API Examples
```
# Adverse events for a drug:
https://api.fda.gov/drug/event.json?search=patient.drug.medicinalproduct:"aspirin"

# Drug labels:
https://api.fda.gov/drug/label.json?search=openfda.brand_name:"Prozac"
```

#### Licensing
- **Public Domain** - No restrictions

#### Pharmacogenomics Focus
- **Adverse Events:** Pharmacovigilance data
- **Drug Labels:** Include PGx information

---

### 16. FDA Table of Pharmacogenomic Biomarkers

**URL:** https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling
**Maintainer:** FDA Division of Translational and Precision Medicine

#### Content
- List of drugs with pharmacogenomic biomarkers in labeling
- Gene-drug associations requiring testing or providing information
- Labeling sections affected

#### Coverage (January 2025)
- **Drugs:** 541 with PGx biomarkers
- Updated periodically (last: June 2024)

#### Access Method
- **Web Page:** Browsable table
- **PDF Download:** Detailed version with labeling text

#### Data Format
- PDF
- Web table (HTML)

#### Licensing
- **Public Domain**

#### Pharmacogenomics Focus
- **AUTHORITATIVE SOURCE:** FDA-recognized biomarkers
- **Clinical Relevance:** Required/recommended testing

---

### 17. RxNorm

**URL:** https://www.nlm.nih.gov/research/umls/rxnorm
**Maintainer:** NIH/NLM

#### Content
- Normalized drug names
- Links between drug vocabularies
- Drug-drug interactions (via DrugBank)
- Prescribable content

#### Coverage
- Links to: First Databank, Micromedex, Gold Standard, Multum, DrugBank

#### Access Method
- **Web Browser (RxNav):** https://lhncbc.nlm.nih.gov/RxNav/
- **REST API:** RxNorm API, RxClass API
- **Bulk Download:** RRF files (Rich Release Format)
- **RxNav-in-a-Box:** Local installation

#### Data Format
- RRF (pipe-delimited text)
- UTF-8 encoding
- MySQL/Oracle load scripts included

#### Licensing
- **Free for prescribable subset** (SAB=RXNORM, SAB=MTHSPL)
- Some source data may require additional licensing

#### Pharmacogenomics Focus
- **Drug Terminology:** Standardization for EHR integration
- **Interoperability:** Links clinical systems

---

## Specialized Databases

### 18. CIViC (Clinical Interpretation of Variants in Cancer)

**URL:** https://civicdb.org
**Maintainer:** Washington University, Griffith Lab

#### Content
- Clinical interpretation of cancer variants
- Therapeutic, prognostic, diagnostic relevance
- Drug-variant relationships
- Evidence items from literature

#### Coverage
- **Variants:** 3,200+
- **Genes:** 470+
- **Publications:** 3,100+
- **Contributors:** 300+

#### Access Method
- **Web Interface:** No login required
- **REST API:** Full programmatic access
- **Bulk Downloads:** Nightly updates, monthly stable releases
- **AWS Open Data:** Registry listing

#### Data Format
- JSON (API)
- TSV (downloads)

#### Licensing
- **Public Domain (CC0)**
- MIT license for software

#### Pharmacogenomics Focus
- **Cancer PGx:** Variant-drug sensitivity/resistance
- **Actionable Variants:** Clinical significance

---

### 19. PGxDB

**URL:** https://pgx-db.org
**Maintainer:** Academic consortium (launched August 2024)

#### Content
- Integrated pharmacogenomics data
- Drug-gene-variant associations
- Adverse drug reactions
- Gene-based association statistics

#### Data Sources Integrated
- DrugBank, PharmGKB, ChEMBL
- UniProt, Ensembl
- SIDER, DisGeNet
- ClinicalTrials.gov

#### Access Method
- **Web Interface:** Interactive platform
- **API:** Programmatic access
- **GitHub:** Source code and data

#### Data Format
- Web tools
- Downloadable data files

#### Licensing
- **Open Access** - FAIR compliant

#### Pharmacogenomics Focus
- **Integrated PGx Platform:** Combines multiple data types
- **Research-Focused:** Supports translational studies

---

### 20. SuperDRUG2

**URL:** http://cheminfo.charite.de/superdrug2
**Maintainer:** Charite Berlin

#### Content
- Approved/marketed drugs
- Chemical structures (2D/3D)
- Drug targets
- Pharmacokinetic data
- Drug-drug interactions
- Side effects

#### Coverage
- **Active Pharmaceutical Ingredients:** 4,587

#### Access Method
- **Web Interface:** Search and browse
- **Download:** Customized download links

#### Data Format
- MySQL database
- Downloadable files

#### Licensing
- **Open Access**

#### Pharmacogenomics Focus
- **Drug Properties:** Comprehensive drug profiles
- **PK Simulation:** Plasma concentration curves

---

### 21. KEGG Drug

**URL:** https://www.genome.jp/kegg/drug/
**Maintainer:** Kanehisa Laboratories

#### Content
- Approved drugs (Japan, USA, Europe)
- Chemical structures
- Target molecules
- Metabolizing enzymes
- Drug interactions

#### Access Method
- **KEGG REST API:** For academic use
- **Web Interface:** Free browsing
- **Bulk Download:** Paid subscription for most data

#### Data Format
- KEGG flat file format
- API responses

#### Licensing
- **API:** Academic use only
- **Medicus directory:** Free (includes KEGG Drug)
- **Full FTP:** Subscription required

#### Pharmacogenomics Focus
- **Drug-Enzyme Relationships:** Metabolism pathways
- **Drug-Target Networks:** Molecular interactions

---

### 22. UniProt

**URL:** https://www.uniprot.org
**Maintainer:** UniProt Consortium (EBI, SIB, PIR)

#### Content
- Protein sequences and functions
- Post-translational modifications
- Disease associations
- Drug binding sites (via cross-references)

#### Coverage (2024)
- **Sequences:** 246 million in UniProtKB

#### Access Method
- **Web Interface:** Free, powerful search
- **REST API:** Full programmatic access
- **Bulk Download:** FTP site

#### Data Format
- FASTA, XML, JSON, RDF
- TSV exports

#### Licensing
- **Open Access (CC BY 4.0)**

#### Pharmacogenomics Focus
- **Drug Targets:** Protein information for drug targets
- **Variants:** Disease and functional variant annotations
- **Cross-References:** Links to DrugBank, ChEMBL, etc.

---

## Comparison Matrix

### Data Coverage Comparison

| Database | Drugs | Targets | Interactions | SNP-Drug | Dosing |
|----------|-------|---------|--------------|----------|--------|
| PharmGKB | 1,000+ | 700+ genes | 20,000+ annotations | YES | YES |
| CPIC | 164 | 34 genes | - | YES | YES |
| DrugBank | 15,000+ | 5,000+ | 2.5M DDI | YES | Partial |
| ChEMBL | 2.8M compounds | 17,803 | 21M bioactivities | No | No |
| DGIdb | 20,000+ | 10,000+ | 70,000+ | No | No |
| Open Targets | Integrated | Genome-wide | Evidence scores | Genetic | No |
| BindingDB | 1.3M | Thousands | 2.9M affinity | No | No |
| GtoPdb | 13,131 | 3,097 | Quantitative | No | No |
| TTD | Ligands | 3,131 | Target-drug | No | No |

### Access and Licensing Comparison

| Database | API | Bulk Download | Academic | Commercial | Format |
|----------|-----|---------------|----------|------------|--------|
| PharmGKB | REST/JSON | Yes (license) | Free | Free (CC BY-SA) | TSV, JSON |
| CPIC | Via PharmGKB | Yes | Free (CC0) | Free (CC0) | PDF, JSON |
| DrugBank | REST | Yes | Free | Paid | XML, CSV |
| ChEMBL | REST | Yes | Free | Free | SQL, SDF |
| PubChem | PUG REST | Yes | Free | Free | SDF, JSON |
| DGIdb | GraphQL | Yes | Free | Free | TSV, JSON |
| Open Targets | GraphQL | Yes | Free | Free | Parquet |
| OpenFDA | REST | Yes | Free | Free | JSON |
| BindingDB | REST | Yes | Free | Free | SDF, TSV |
| GtoPdb | REST | Yes | Free | Free | CSV, PostgreSQL |

### Pharmacogenomics Focus Comparison

| Database | SNP-Drug Pairs | Dosing Guidelines | Metabolizer Status | Clinical Implementation |
|----------|----------------|-------------------|--------------------|-----------------------|
| PharmGKB | PRIMARY | Yes (CPIC/DPWG) | Yes | Yes |
| CPIC | Yes | PRIMARY | Yes | Yes (EHR ready) |
| DPWG | Yes | Yes | Yes | Yes (Dutch systems) |
| PharmVar | Allele definitions | No | Function tables | Nomenclature |
| DrugBank | Partial | No | Enzyme info | Partial |
| FDA PGx Table | Biomarkers | In labels | No | FDA guidance |
| CIViC | Cancer variants | Therapeutic | No | Cancer treatment |

---

## Recommendations

### For SNP-Drug Relationships and Clinical Dosing

**PRIMARY:** PharmGKB + CPIC + DPWG
- Most comprehensive pharmacogenomics data
- Clinically validated dosing guidelines
- Ready for clinical implementation

### For Comprehensive Drug Information

**PRIMARY:** DrugBank (academic license) + ChEMBL
- Complete drug profiles with targets
- Quantitative bioactivity data
- Well-structured schemas

### For Drug-Gene Interactions

**PRIMARY:** DGIdb + Open Targets
- Aggregated from multiple sources
- Free, open access
- Good API support

### For Drug Targets and Binding Data

**PRIMARY:** BindingDB + ChEMBL + GtoPdb
- Quantitative affinity measurements
- Expert-curated target information
- Research-grade data

### For Regulatory Compliance

**PRIMARY:** OpenFDA + DailyMed + FDA PGx Biomarkers Table
- Official FDA data
- Drug labels with PGx information
- Adverse event data

### Integration Priority for Health Knowledge Base

1. **Tier 1 (Essential):**
   - PharmGKB/ClinPGx (pharmacogenomics gold standard)
   - CPIC guidelines (clinical dosing)
   - DrugBank (drug information)
   - OpenFDA (adverse events, labels)

2. **Tier 2 (Important):**
   - ChEMBL (bioactivity data)
   - DGIdb (drug-gene interactions)
   - Open Targets (target-disease associations)
   - PharmVar (allele nomenclature)

3. **Tier 3 (Supplementary):**
   - BindingDB (binding affinity)
   - GtoPdb (receptor pharmacology)
   - TTD (target druggability)
   - RxNorm (terminology standardization)
   - CIViC (cancer variants)

### Data Pipeline Recommendations

```
User Query (Gene/Drug/Variant)
         |
         v
+------------------+
| Normalize IDs    | <-- RxNorm, PharmVar, dbSNP
+------------------+
         |
         v
+------------------+
| PharmGKB/CPIC    | --> SNP-Drug associations, Dosing
+------------------+
         |
         v
+------------------+
| DrugBank/ChEMBL  | --> Drug info, Targets, Mechanisms
+------------------+
         |
         v
+------------------+
| DGIdb/Open Targets| --> Drug-Gene interactions
+------------------+
         |
         v
+------------------+
| OpenFDA/DailyMed | --> Labels, Adverse events
+------------------+
         |
         v
    Integrated Response
```

---

## References

1. PharmGKB: https://www.pharmgkb.org
2. CPIC: https://cpicpgx.org
3. DrugBank: https://go.drugbank.com
4. ChEMBL: https://www.ebi.ac.uk/chembl
5. PubChem: https://pubchem.ncbi.nlm.nih.gov
6. Open Targets: https://platform.opentargets.org
7. DGIdb: https://dgidb.org
8. TTD: https://idrblab.net/ttd
9. BindingDB: https://www.bindingdb.org
10. GtoPdb: https://www.guidetopharmacology.org
11. OpenFDA: https://open.fda.gov
12. DailyMed: https://dailymed.nlm.nih.gov
13. RxNorm: https://www.nlm.nih.gov/research/umls/rxnorm
14. PharmVar: https://www.pharmvar.org
15. CIViC: https://civicdb.org
16. DrugCentral: https://drugcentral.org
17. UniProt: https://www.uniprot.org
18. KEGG Drug: https://www.genome.jp/kegg/drug
19. PGxDB: https://pgx-db.org
20. SuperDRUG2: http://cheminfo.charite.de/superdrug2

---

*Document generated: January 2026*
*Last updated: January 17, 2026*
