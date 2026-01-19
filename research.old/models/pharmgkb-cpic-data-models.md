# Pharmacogenomics Database Data Models

This document provides comprehensive data model documentation for major pharmacogenomics databases used in clinical implementation and research.

## Table of Contents

1. [PharmGKB (ClinPGx)](#1-pharmgkb-clinpgx)
2. [CPIC Guidelines Database](#2-cpic-guidelines-database)
3. [PharmVar](#3-pharmvar)
4. [FDA Pharmacogenomic Biomarkers Table](#4-fda-pharmacogenomic-biomarkers-table)
5. [Cross-Database Relationships](#5-cross-database-relationships)

---

## 1. PharmGKB (ClinPGx)

PharmGKB (Pharmacogenomics Knowledge Base) is an NIH-funded resource that provides information about how human genetic variation affects response to medications. The database has been rebranded as ClinPGx with API access at `https://api.pharmgkb.org/`.

### 1.1 API Overview

**Base URL:** `https://api.pharmgkb.org/v1/`

**Documentation:** OpenAPI 3.0 specification available at `https://api.pharmgkb.org/swagger/`

**Rate Limit:** Maximum 2 requests per second

**License:** Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)

### 1.2 Core Entity Types

#### 1.2.1 Genes

| Field | Type | Description |
|-------|------|-------------|
| `accessionId` | String | PharmGKB accession ID (format: `PA####`) |
| `symbol` | String | HGNC official gene symbol |
| `name` | String | Full gene name |
| `alternateNames` | Array[String] | Alternative names and aliases |
| `alternateSymbols` | Array[String] | Alternative gene symbols |
| `chromosome` | String | Chromosome location (e.g., "chr1") |
| `chromosomeStart` | Integer | GRCh38 start position |
| `chromosomeEnd` | Integer | GRCh38 end position |
| `strand` | String | Strand orientation (+/-) |
| `ensemblId` | String | Ensembl gene ID |
| `ncbiGeneId` | String | NCBI Gene ID |
| `hgncId` | String | HGNC ID |
| `hasVip` | Boolean | Is Very Important Pharmacogene |
| `hasHaplotypes` | Boolean | Has defined star alleles |
| `crossReferences` | Array[Object] | External database references |

**Example:**
```json
{
  "accessionId": "PA128",
  "symbol": "CYP2D6",
  "name": "cytochrome P450 family 2 subfamily D member 6",
  "chromosome": "chr22",
  "chromosomeStart": 42126499,
  "chromosomeEnd": 42130810,
  "strand": "-",
  "ensemblId": "ENSG00000100197",
  "ncbiGeneId": "1565",
  "hgncId": "HGNC:2625",
  "hasVip": true,
  "hasHaplotypes": true
}
```

#### 1.2.2 Drugs/Chemicals

| Field | Type | Description |
|-------|------|-------------|
| `accessionId` | String | PharmGKB accession ID (format: `PA#####`) |
| `name` | String | Generic drug name |
| `genericNames` | Array[String] | All generic names |
| `tradeNames` | Array[String] | Brand/trade names |
| `type` | String | Chemical type (Drug, Metabolite, etc.) |
| `rxnormId` | String | RxNorm concept ID (RxCUI) |
| `atcCodes` | Array[String] | ATC classification codes |
| `drugBankId` | String | DrugBank accession |
| `pubchemCid` | String | PubChem Compound ID |
| `meshId` | String | MeSH descriptor ID |
| `smiles` | String | SMILES chemical structure |
| `inchi` | String | InChI identifier |
| `molecularWeight` | Float | Molecular weight |
| `indication` | String | Therapeutic indication |
| `mechanismOfAction` | String | Drug mechanism |
| `hasRx` | Boolean | Is prescription drug |
| `crossReferences` | Array[Object] | External database links |

**Example:**
```json
{
  "accessionId": "PA449015",
  "name": "codeine",
  "genericNames": ["codeine"],
  "tradeNames": ["Tylenol with Codeine"],
  "type": "Drug",
  "rxnormId": "2670",
  "atcCodes": ["R05DA04", "N02AA08"],
  "drugBankId": "DB00318",
  "pubchemCid": "5284371"
}
```

#### 1.2.3 Variants

| Field | Type | Description |
|-------|------|-------------|
| `accessionId` | String | PharmGKB accession ID |
| `name` | String | Variant name |
| `rsid` | String | dbSNP rsID (e.g., "rs1045642") |
| `location` | Object | Genomic location data |
| `location.chromosome` | String | Chromosome |
| `location.position` | Integer | GRCh38 position |
| `location.assembly` | String | Genome assembly version |
| `referenceAllele` | String | Reference allele |
| `alternateAlleles` | Array[String] | Alternate alleles |
| `variantType` | String | SNP, insertion, deletion, etc. |
| `hgvs` | Object | HGVS nomenclature |
| `hgvs.genomic` | String | Genomic HGVS (g.) |
| `hgvs.coding` | String | Coding HGVS (c.) |
| `hgvs.protein` | String | Protein HGVS (p.) |
| `genes` | Array[Object] | Associated genes |
| `functionalConsequence` | String | Variant effect |
| `frequencies` | Array[Object] | Population allele frequencies |
| `clinicalAnnotations` | Array[Object] | Clinical annotation links |

**Example:**
```json
{
  "accessionId": "PA166153735",
  "name": "rs4149056",
  "rsid": "rs4149056",
  "location": {
    "chromosome": "chr12",
    "position": 21176804,
    "assembly": "GRCh38"
  },
  "referenceAllele": "T",
  "alternateAlleles": ["C"],
  "variantType": "SNP",
  "hgvs": {
    "genomic": "NC_000012.12:g.21176804T>C",
    "coding": "NM_006446.5:c.521T>C",
    "protein": "NP_006437.3:p.Val174Ala"
  }
}
```

#### 1.2.4 Clinical Annotations

| Field | Type | Description |
|-------|------|-------------|
| `id` | String | Clinical annotation ID |
| `type` | String | "variant-level" or "gene-level" |
| `variant` | Object | Associated variant/allele |
| `drug` | Object | Associated drug |
| `gene` | Object | Associated gene |
| `levelOfEvidence` | String | Evidence level (1A, 1B, 2A, 2B, 3, 4) |
| `score` | Float | Calculated evidence score |
| `phenotypeCategories` | Array[String] | Phenotype categories |
| `genotypes` | Array[Object] | Genotype-specific annotations |
| `genotypes[].genotype` | String | Genotype (e.g., "CC", "CT", "TT") |
| `genotypes[].phenotype` | String | Associated phenotype |
| `genotypes[].description` | String | Clinical description |
| `variantAnnotations` | Array[Object] | Supporting variant annotations |
| `literature` | Array[Object] | Literature citations |
| `lastModified` | DateTime | Last update timestamp |

**Genotype Annotation Structure:**
```json
{
  "genotype": "CC",
  "phenotype": "Increased myopathy risk",
  "description": "Patients with the CC genotype at rs4149056...",
  "strength": "Strong",
  "direction": "Increased risk"
}
```

### 1.3 Evidence Level Codes

PharmGKB uses a 6-tier evidence classification system:

| Level | Name | Criteria |
|-------|------|----------|
| **1A** | Guideline-based | Variant-drug in CPIC guideline, medical society guideline, or implemented at PGRN site/major health system |
| **1B** | High Evidence | Replicated in multiple cohorts with significant p-values and preferably strong effect size |
| **2A** | Moderate (VIP) | Level 2B criteria + variant in Very Important Pharmacogene |
| **2B** | Moderate Evidence | Replicated association but may have inconsistent significance or small effect size |
| **3** | Low Evidence | Single study, conflicting results, case reports, in vitro studies, or unreplicated GWAS |
| **4** | Unsupported | Preponderance of evidence does not support association (negative score) |

**Scoring Formula:**
```
Variant Score = (Step1 + Step2 + Step3 + Step4) x (Step5A x Step5B)
```

| Step | Factor | Points |
|------|--------|--------|
| Step 1 | Phenotype Category | 0-1 (clinical outcomes = 1, PK = 0.5) |
| Step 2 | P-value | 0-1 (p<0.01 = 1, p<=0.05 = 0.5) |
| Step 3 | Cohort Size | 0-1 (>500 = 1, <=50 = 0) |
| Step 4 | Effect Size | 0-0.5 (OR/HR >2 or <0.5 = 0.5) |
| Step 5A | Study Type | Multiplier 0-1 (in vitro = 0) |
| Step 5B | Association | Multiplier -1 to 1 (negative = -1) |

**Score Thresholds:**

| Level | Standard Range | Rare Variant Range |
|-------|----------------|-------------------|
| 1A | >=80 (+ guideline) | >=80 (+ guideline) |
| 1B | 25-79.9375 | 10-79.9375 |
| 2A/2B | 8-24.9375 | 3-9.9375 |
| 3 | 0-7.9375 | 0-2.9375 |
| 4 | <0 | <0 |

### 1.4 Phenotype Categories

| Category | Description |
|----------|-------------|
| `Toxicity` | Adverse drug reactions |
| `Efficacy` | Drug effectiveness |
| `Dosage` | Dose requirements |
| `Metabolism/PK` | Pharmacokinetic parameters |
| `PD` | Pharmacodynamic effects |
| `Other` | Other phenotypes |

### 1.5 Variant Annotations

| Field | Type | Description |
|-------|------|-------------|
| `id` | String | Annotation ID |
| `variant` | Object | Variant reference |
| `drug` | Object | Drug reference |
| `literature` | Object | Publication reference (PMID) |
| `phenotypeCategory` | String | Phenotype type |
| `significance` | String | Statistical significance |
| `pValue` | Float | P-value from study |
| `cohortSize` | Integer | Study cohort size |
| `ethnicity` | String | Study population |
| `alleles` | String | Alleles studied |
| `effectSize` | Object | OR/HR/RR with confidence interval |
| `sentence` | String | Key finding sentence |
| `studyParameters` | Object | Study metadata |

### 1.6 Pathways

| Field | Type | Description |
|-------|------|-------------|
| `accessionId` | String | Pathway accession ID |
| `name` | String | Pathway name |
| `type` | String | "PK" (pharmacokinetics) or "PD" (pharmacodynamics) |
| `drug` | Object | Central drug |
| `genes` | Array[Object] | Pathway genes |
| `tissues` | Array[String] | Relevant tissues |
| `description` | String | Pathway description |
| `nodes` | Array[Object] | Pathway diagram nodes |
| `edges` | Array[Object] | Pathway connections |
| `literature` | Array[Object] | Supporting publications |

**Download Formats:**
- TSV (tabular data)
- BioPAX Level 3 (OWL/RDF format)
- GPML (GenMAPP Pathway Markup Language)
- PDF (diagram)
- Adobe Illustrator (.ai)

### 1.7 Bulk Download TSV Files

#### genes.tsv

| Column | Description |
|--------|-------------|
| PharmGKB Accession Id | PA#### format |
| NCBI Gene ID | Entrez Gene ID |
| HGNC ID | HGNC identifier |
| Ensembl Id | Ensembl gene ID |
| Name | Gene name |
| Symbol | HGNC symbol |
| Alternate Names | Pipe-delimited aliases |
| Alternate Symbols | Pipe-delimited symbols |
| Is VIP | true/false |
| Has Variant Annotation | true/false |
| Cross-references | External IDs |
| Has CPIC Dosing Guideline | true/false |
| Chromosome | Chromosome name |
| Chromosomal Start - GRCh37 | Start position |
| Chromosomal Stop - GRCh37 | End position |
| Chromosomal Start - GRCh38 | Start position |
| Chromosomal Stop - GRCh38 | End position |

#### drugs.tsv

| Column | Description |
|--------|-------------|
| PharmGKB Accession Id | PA##### format |
| Name | Drug name |
| Generic Names | Pipe-delimited |
| Trade Names | Pipe-delimited |
| Brand Mixtures | Combination products |
| Type | Drug, Metabolite, etc. |
| Cross-references | External IDs |
| SMILES | Chemical structure |
| InChI | InChI identifier |
| Dosing Guideline | Has guideline |
| External Vocabulary | Term mappings |
| Clinical Annotation Count | Number of annotations |
| Variant Annotation Count | Number of annotations |
| Pathway Count | Number of pathways |
| VIP Count | Associated VIP genes |
| Dosing Guideline Sources | CPIC, DPWG, etc. |
| Top Clinical Annotation Level | Highest evidence level |
| Top FDA Label Testing Level | FDA recommendation |

#### clinical_annotations.tsv

| Column | Description |
|--------|-------------|
| Clinical Annotation ID | Unique identifier |
| Variant/Haplotypes | rsID or star alleles |
| Gene | Gene symbol |
| Level of Evidence | 1A, 1B, 2A, 2B, 3, 4 |
| Level Override | Manual override reason |
| Level Modifiers | Additional context |
| Score | Calculated score |
| Phenotype Category | Toxicity, Efficacy, etc. |
| PMID Count | Literature count |
| Evidence Count | Supporting evidence |
| Drug(s) | Associated drugs |
| Phenotype(s) | Clinical phenotypes |
| Latest History Date | Last update |
| URL | Web page link |
| Specialty Population | If applicable |

#### relationships.tsv

| Column | Description |
|--------|-------------|
| Entity1_id | First entity accession |
| Entity1_name | First entity name |
| Entity1_type | Gene, Drug, Disease, Variant |
| Entity2_id | Second entity accession |
| Entity2_name | Second entity name |
| Entity2_type | Gene, Drug, Disease, Variant |
| Evidence | Supporting evidence |
| Association | Relationship type |
| PK | Pharmacokinetics flag |
| PD | Pharmacodynamics flag |
| PMIDs | Literature references |

### 1.8 Relationship Types

| Relationship | Description |
|--------------|-------------|
| `associated` | General association |
| `affects` | Entity affects another |
| `is affected by` | Passive relationship |
| `metabolizes` | Gene metabolizes drug |
| `is metabolized by` | Drug is metabolized |
| `transports` | Gene transports drug |
| `is transported by` | Drug is transported |
| `targets` | Drug targets gene |
| `is targeted by` | Gene is targeted |

### 1.9 PGx Levels (Drug Labels)

| Level | Description |
|-------|-------------|
| `Testing required` | Genetic testing must be conducted before using drug |
| `Testing recommended` | Genetic testing recommended but not required |
| `Actionable PGx` | Describes variant impact without mandating testing |
| `Informative PGx` | States variants have no clinically significant effect |

---

## 2. CPIC Guidelines Database

The Clinical Pharmacogenetics Implementation Consortium (CPIC) provides standardized clinical guidelines for pharmacogenetic implementation.

### 2.1 API Overview

**Base URL:** `https://api.cpicpgx.org/v1/`

**Technology:** PostgreSQL database with PostgREST API

**Documentation:** `https://github.com/cpicpgx/cpic-data/wiki`

### 2.2 Core Tables

#### 2.2.1 Gene Table

| Field | Type | Description |
|-------|------|-------------|
| `symbol` | String (PK) | HGNC official gene symbol |
| `chr` | String | Chromosome |
| `genesequenceid` | String | NCBI reference sequence |
| `proteinsequenceid` | String | Protein sequence ID |
| `chromosequenceid` | String | Chromosome sequence ID |
| `mrnastart` | Integer | mRNA start position |
| `mrnaend` | Integer | mRNA end position |
| `lookupmethod` | Enum | "PHENOTYPE", "ACTIVITY_SCORE", or "ALLELE_STATUS" |
| `url` | String | PharmGKB gene page URL |
| `hgncid` | String | HGNC ID |
| `ncbiid` | String | NCBI Gene ID |
| `ensemblid` | String | Ensembl ID |

**Example API Query:**
```bash
curl 'https://api.cpicpgx.org/v1/gene?symbol=eq.CYP2D6'
```

#### 2.2.2 Drug Table

| Field | Type | Description |
|-------|------|-------------|
| `drugid` | String (PK) | Format: "source:id" (e.g., "RxNorm:2670") |
| `name` | String | Drug name |
| `rxnormid` | String | RxNorm CUI |
| `atcid` | Array[String] | ATC codes |
| `drugbankid` | String | DrugBank ID |
| `flowchart` | String | CDS flowchart URL |
| `guidelineid` | UUID | Parent guideline |

#### 2.2.3 Allele Definition Table

The allele definition model is split into two tables to handle copy number variations:

**allele_definition Table:**

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Unique definition ID |
| `genesymbol` | String (FK) | Gene symbol |
| `name` | String | Allele name (e.g., "*2") |
| `pharmvarid` | String | PharmVar ID |
| `function` | String | Allele function category |
| `activityvalue` | Float | Activity score value |

**allele Table:**

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Unique allele ID |
| `definitionid` | Integer (FK) | Reference to allele_definition |
| `name` | String | Full allele name (e.g., "*2x2") |
| `copynumber` | Integer | Gene copy number |
| `frequency` | Object | Population frequencies |
| `isreference` | Boolean | Is reference allele |

**allele_location_value Table:**

| Field | Type | Description |
|-------|------|-------------|
| `definitionid` | Integer (FK) | Allele definition reference |
| `locationid` | Integer (FK) | Sequence location reference |
| `variantallele` | String | Variant allele at position |

**sequence_location Table:**

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Location ID |
| `genesymbol` | String (FK) | Gene symbol |
| `chromosomelocation` | String | GRCh38 coordinate |
| `genelocation` | String | Gene coordinate |
| `proteinlocation` | String | Protein position |
| `dbsnpid` | String | rsID |
| `referenceallele` | String | Reference allele |

### 2.3 Allele Definition Table Schema (Excel Format)

CPIC distributes allele definitions as Excel files with a matrix structure:

**Columns:**
- Column A: Allele name (e.g., *1, *2, *3)
- Columns B onwards: Variant positions (rsIDs or genomic coordinates)

**Rows:**
- Row 1: Header with rsID/position identifiers
- Row 2: Reference allele (typically *1)
- Rows 3+: Named alleles with variant calls

**Cell Values:**
- Reference base at position (blank = same as reference)
- Alternate allele (e.g., "A", "G", "del", "ins")

**Example Structure:**

| Allele | rs12248560 | rs28399504 | rs41291556 | rs4244285 | Score |
|--------|------------|------------|------------|-----------|-------|
| *1 | C | A | G | G | 4 |
| *2 | T | | | | 1 |
| *3 | | G | | | 1 |
| *4 | | | A | | 1 |
| *17 | T | | | | 1 |

### 2.4 Allele Functionality Table

| Field | Type | Description |
|-------|------|-------------|
| `genesymbol` | String | Gene symbol |
| `allele` | String | Star allele name |
| `function` | String | Functional category |
| `activityvalue` | Float | Numeric activity score |
| `clinicalfunction` | String | Clinical interpretation |
| `citations` | Array[String] | Supporting PMIDs |

**Functional Categories:**

| Function | Description | Typical Activity Value |
|----------|-------------|----------------------|
| `Normal function` | Full enzyme activity | 1.0 |
| `Decreased function` | Reduced activity | 0.5 |
| `No function` | No enzyme activity | 0 |
| `Increased function` | Enhanced activity | >1.0 |
| `Uncertain function` | Unknown effect | null |
| `Unknown function` | Not characterized | null |

### 2.5 Gene Result Tables (3-Level Hierarchy)

#### gene_result Table

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Result ID |
| `genesymbol` | String (FK) | Gene symbol |
| `result` | String | Phenotype name (e.g., "Poor Metabolizer") |
| `ehrpriority` | Integer | EHR display priority |
| `consultationtext` | String | Example consultation text |
| `activityscore` | String | Activity score range |

#### gene_result_lookup Table

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Lookup ID |
| `generesultid` | Integer (FK) | Reference to gene_result |
| `function1` | String | First allele function |
| `function2` | String | Second allele function |
| `activityscore` | Float | Combined activity score |
| `description` | String | Combination description |

#### gene_result_diplotype Table

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Diplotype ID |
| `generesultid` | Integer (FK) | Reference to gene_result |
| `diplotype` | String | Diplotype notation (e.g., "*1/*2") |
| `ehrpriority` | Integer | EHR display priority |

### 2.6 Diplotype View

Pre-joined view for convenient querying:

| Field | Description |
|-------|-------------|
| `genesymbol` | Gene symbol |
| `diplotype` | Diplotype notation |
| `phenotype` | Assigned phenotype |
| `ehrpriority` | Display priority |
| `consultationtext` | Clinical text |
| `activityscore` | Activity score |
| `function1` | Allele 1 function |
| `function2` | Allele 2 function |

**Example Query:**
```bash
curl 'https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.CYP2D6&diplotype=eq.*1/*4'
```

**Example Response:**
```json
{
  "genesymbol": "CYP2D6",
  "diplotype": "*1/*4",
  "phenotype": "Intermediate Metabolizer",
  "activityscore": "1.0",
  "function1": "Normal function",
  "function2": "No function"
}
```

### 2.7 Diplotype-Phenotype Table Schema

| Field | Type | Description |
|-------|------|-------------|
| `genesymbol` | String | Gene symbol |
| `diplotype` | String | Diplotype (e.g., "*1/*1") |
| `ehr_priority` | Integer | EHR display priority |
| `activity_score` | String | Activity score value/range |
| `phenotype` | String | Metabolizer status |
| `description` | String | Clinical description |

**Standard Phenotypes:**

| Gene Type | Phenotypes |
|-----------|------------|
| CYP Enzymes | Ultrarapid, Rapid, Normal, Intermediate, Poor |
| Transporters | Increased, Normal, Decreased, Poor function |
| Receptors | Normal, Reduced, Absent function |

### 2.8 Recommendation Table

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Recommendation ID |
| `guidelineid` | UUID (FK) | Parent guideline |
| `drugid` | String (FK) | Drug reference |
| `lookupkey` | JSONB | Combined gene status keys |
| `population` | String | Target population (null = general) |
| `classification` | String | Recommendation strength |
| `implication` | String | Clinical implication text |
| `recommendation` | String | Dosing recommendation |
| `comments` | String | Additional comments |
| `phenotypes` | JSONB | Gene phenotype requirements |
| `activityscore` | String | Required activity score |

**Lookup Methods:**

| Method | Description | Example Key |
|--------|-------------|-------------|
| `PHENOTYPE` | Phenotype-based | `{"CYP2D6": "Poor Metabolizer"}` |
| `ACTIVITY_SCORE` | Activity score-based | `{"CYP2D6": "0"}` |
| `ALLELE_STATUS` | Specific allele | `{"HLA-B": "*57:01 positive"}` |

**Example Query:**
```bash
curl 'https://api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:2670&select=*'
```

### 2.9 Test Alert Table

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Alert ID |
| `drugid` | String (FK) | Drug reference |
| `genesymbol` | String (FK) | Gene symbol |
| `phenotype` | String | Triggering phenotype |
| `activityscore` | String | Triggering activity score |
| `population` | String | Target population |
| `alerttext` | String | Alert message text |
| `alerttype` | String | Pre-test or post-test |

### 2.10 Pair View (Gene-Drug Pairs)

| Field | Type | Description |
|-------|------|-------------|
| `genesymbol` | String | Gene symbol |
| `drugid` | String | Drug ID |
| `drugname` | String | Drug name |
| `cpiclevel` | String | CPIC evidence level |
| `cpicstatus` | String | Guideline status |
| `pgkblevel` | String | PharmGKB annotation level |
| `pgxlevel` | String | PGx testing level |
| `fdalevel` | String | FDA label level |

**CPIC Levels:**

| Level | Description |
|-------|-------------|
| A | Prescribing action recommended |
| A/B | Mixed evidence |
| B | Prescribing action may be recommended |
| C | No prescribing action recommended |
| D | Insufficient evidence |

### 2.11 Frequency Data

| Field | Type | Description |
|-------|------|-------------|
| `genesymbol` | String | Gene symbol |
| `allelename` | String | Star allele |
| `population` | String | Population/ethnicity |
| `frequency` | Float | Allele frequency (0-1) |
| `frequencymethod` | String | Determination method |
| `samplesize` | Integer | Study sample size |
| `pmid` | String | Literature reference |

### 2.12 JSON Download Format

CPIC provides JSON exports with nested structure:

```json
{
  "guideline": {
    "id": "uuid",
    "name": "CPIC Guideline for codeine and CYP2D6",
    "version": "3.0",
    "genes": ["CYP2D6"],
    "drugs": ["codeine"],
    "url": "https://cpicpgx.org/guidelines/..."
  },
  "citations": [
    {"pmid": "12345678", "citation": "..."}
  ],
  "recommendations": [
    {
      "drug": {
        "name": "codeine",
        "rxnormId": "2670"
      },
      "genes": {
        "CYP2D6": {
          "phenotype": "Poor Metabolizer",
          "activityScore": "0"
        }
      },
      "recommendation": "Avoid codeine use...",
      "implication": "Greatly reduced morphine formation...",
      "classification": "Strong"
    }
  ],
  "alleleDefinitions": {...},
  "alleleFunctionality": {...},
  "diplotypePhenotypes": {...}
}
```

---

## 3. PharmVar

The Pharmacogene Variation (PharmVar) Consortium maintains the official repository for pharmacogene allele nomenclature.

### 3.1 Overview

**Website:** `https://www.pharmvar.org/`

**Purpose:** Standardized star allele nomenclature and definitions

**Data Updates:** Monthly releases with versioned database

### 3.2 Star Allele Nomenclature

**Format:** `GENE*Allele.Suballele`

| Component | Description | Example |
|-----------|-------------|---------|
| Gene symbol | HGNC gene symbol | CYP2D6 |
| Asterisk | Separator | * |
| Major allele | Primary allele number | 4 |
| Suballele | Variant-specific | .001, .002 |

**Examples:**
- `CYP2D6*1` - Reference allele
- `CYP2D6*4` - Core allele (all *4 suballeles)
- `CYP2D6*4.001` - Specific suballele
- `CYP2D6*4x2` - Gene duplication with *4

### 3.3 PharmVar ID Format

Each allele receives a unique, immutable PharmVar ID:

| Format | Description |
|--------|-------------|
| `PV#####` | Unique PharmVar identifier |

**Key Properties:**
- Immutable: ID never changes
- Version-tracked: Definition changes trigger new ID
- Retired IDs: Marked as retired, not deleted

### 3.4 Allele Definition Structure

#### Core Allele Definition

| Field | Description |
|-------|-------------|
| `pharmvarId` | Unique PV##### identifier |
| `gene` | Gene symbol |
| `allele` | Star allele name |
| `isCore` | Is core allele definition |
| `function` | Functional classification |
| `clinicalFunction` | Clinical interpretation |
| `activity` | Activity value |
| `definitionSource` | Definition source (CPIC, literature) |
| `variants` | Array of defining variants |

#### Variant Definition

| Field | Type | Description |
|-------|------|-------------|
| `rsId` | String | dbSNP rsID |
| `hgvsGenomic` | String | Genomic HGVS (GRCh38) |
| `hgvsCoding` | String | Coding HGVS |
| `hgvsProtein` | String | Protein HGVS |
| `chromosome` | String | Chromosome |
| `position` | Integer | GRCh38 position |
| `referenceAllele` | String | Reference base |
| `variantAllele` | String | Variant base |
| `variantType` | String | SNV, insertion, deletion |
| `impact` | String | Variant impact |
| `isDefining` | Boolean | Core defining variant |

### 3.5 Page Formats

PharmVar provides two display formats:

**Star Format:**
- Organized by star allele
- Shows all variants defining each allele
- Includes suballele hierarchy

**rsID Format:**
- Organized by rsID
- Shows all alleles containing each variant
- Useful for VCF-based lookups

### 3.6 Haplotype Definition Format

| Field | Description |
|-------|-------------|
| `gene` | Gene symbol |
| `allele` | Star allele name |
| `coreAllele` | Parent core allele (for suballeles) |
| `definingVariants` | Array of required variants |
| `additionalVariants` | Non-defining variants (suballele-specific) |
| `copyNumber` | Gene copy number (1, 2, >2) |
| `structuralVariant` | CNV/hybrid information |

### 3.7 Variant-to-Allele Mapping

| rsID | Gene | Alleles | Function Effect |
|------|------|---------|-----------------|
| rs3892097 | CYP2D6 | *4, *4.001, *4.013 | Splice defect |
| rs1065852 | CYP2D6 | *4, *10, *36 | Decreased activity |
| rs5030655 | CYP2D6 | *6 | Frameshift deletion |

### 3.8 Reference Sequence Standards

| Standard | Description |
|----------|-------------|
| RefSeqGene | Gene-specific reference |
| MANE Select | Matched Annotation from NCBI and EBI |
| GRCh38 | Genome assembly for coordinates |
| GRCh37 | Legacy coordinates (also provided) |

### 3.9 HGVS Nomenclature

PharmVar provides HGVS annotations at multiple levels:

| Level | Format | Example |
|-------|--------|---------|
| Genomic | NC_######.#:g.pos | NC_000022.11:g.42128945G>A |
| Coding | NM_######.#:c.pos | NM_000106.6:c.100C>T |
| Protein | NP_######.#:p.pos | NP_000097.3:p.Pro34Ser |

### 3.10 Version Control

| Entity | Version Format | Description |
|--------|---------------|-------------|
| Database | X.Y.Z | Major.Minor.Patch releases |
| Gene | Gene-version | Gene-specific updates |
| Allele | Allele-version | Allele definition version |

**Change Tracking:**
- New allele additions
- Allele mergers (legacy names preserved)
- Allele retirements
- Definition updates (new PharmVar ID issued)

---

## 4. FDA Pharmacogenomic Biomarkers Table

The FDA maintains a table of drugs with pharmacogenomic information in their labeling.

### 4.1 Table Structure

| Column | Data Type | Description |
|--------|-----------|-------------|
| Drug | String | Drug name (hyperlinked to Drugs@FDA) |
| Therapeutic Area | String | Clinical category |
| Biomarker | String | Genetic biomarker name |
| Labeling Sections | String | Where PGx info appears in label |

### 4.2 Therapeutic Areas

| Area | Example Drugs |
|------|---------------|
| Oncology | Pembrolizumab, Olaparib |
| Psychiatry | Atomoxetine, Citalopram |
| Neurology | Carbamazepine, Phenytoin |
| Cardiology | Warfarin, Clopidogrel |
| Infectious Diseases | Abacavir, Efavirenz |
| Gastroenterology | Pantoprazole, Omeprazole |
| Pulmonary | Ivacaftor |
| Rheumatology | Azathioprine |
| Anesthesiology | Codeine, Tramadol |
| Dermatology | Allopurinol |

### 4.3 Biomarker Types

| Type | Description | Examples |
|------|-------------|----------|
| Gene variants | Germline or somatic mutations | EGFR, BRAF, KRAS |
| HLA alleles | HLA typing | HLA-B*57:01, HLA-B*15:02 |
| Enzyme deficiencies | Functional deficiencies | G6PD deficiency |
| Chromosomal abnormalities | Cytogenetic markers | Chromosome 17p deletion |
| Protein markers | Protein expression | PD-L1, HER2 |
| Gene signatures | Expression patterns | Oncotype DX |

### 4.4 Labeling Sections

| Section | When Used |
|---------|-----------|
| Boxed Warning | Serious safety concern requiring test |
| Indications and Usage | Patient selection based on biomarker |
| Dosage and Administration | Dose adjustment based on genotype |
| Warnings and Precautions | Safety information |
| Use in Specific Populations | Population-specific guidance |
| Clinical Pharmacology | PK/PD impact of genetics |
| Clinical Studies | Efficacy by genetic subgroup |

### 4.5 PGx Information Categories

| Category | Description | Testing Requirement |
|----------|-------------|---------------------|
| Required | Must test before prescribing | Mandatory |
| Recommended | Should test before prescribing | Strong suggestion |
| Actionable | Info on dose/selection | Optional testing |
| Informative | Background information | Not required |

### 4.6 Example Data Rows

| Drug | Therapeutic Area | Biomarker | Labeling Sections |
|------|------------------|-----------|-------------------|
| Abacavir | Infectious Diseases | HLA-B*57:01 | Boxed Warning, Contraindications, Warnings and Precautions |
| Carbamazepine | Neurology | HLA-B*15:02, HLA-A*31:01 | Boxed Warning, Warnings and Precautions |
| Clopidogrel | Cardiology | CYP2C19 | Boxed Warning, Dosage and Administration, Clinical Pharmacology |
| Warfarin | Cardiology | CYP2C9, VKORC1, CYP4F2 | Dosage and Administration, Clinical Pharmacology |
| Codeine | Anesthesiology | CYP2D6 | Boxed Warning, Contraindications, Warnings and Precautions |
| Azathioprine | Rheumatology | TPMT, NUDT15 | Dosage and Administration, Warnings and Precautions |
| Irinotecan | Oncology | UGT1A1 | Dosage and Administration, Warnings and Precautions |

### 4.7 Biomarker Nomenclature

| Convention | Usage | Example |
|------------|-------|---------|
| HUGO symbol | Standard gene names | CYP2D6, CYP2C19 |
| HLA nomenclature | HLA alleles | HLA-B*57:01 |
| Descriptive | Non-genetic markers | G6PD deficiency |
| Protein names | Protein biomarkers | CD274 (PD-L1) |

### 4.8 Data Access

**Downloadable Formats:**
- PDF with detailed labeling text
- Online HTML table (updated periodically)

**Update Frequency:**
- Quarterly updates
- Blue text indicates recent changes

---

## 5. Cross-Database Relationships

### 5.1 Identifier Cross-References

| Database | Gene ID | Drug ID | Variant ID |
|----------|---------|---------|------------|
| PharmGKB | PA#### | PA##### | PA########## |
| CPIC | HGNC symbol | RxNorm:#### | rsID |
| PharmVar | HGNC symbol | N/A | PV##### |
| FDA Table | HUGO symbol | Drug name | Biomarker name |

### 5.2 External ID Mappings

| System | Gene | Drug | Variant |
|--------|------|------|---------|
| NCBI | Entrez Gene ID | - | dbSNP rsID |
| Ensembl | ENSG ID | - | - |
| HGNC | HGNC:#### | - | - |
| RxNorm | - | RxCUI | - |
| ATC | - | ATC code | - |
| DrugBank | - | DB##### | - |
| ClinVar | - | - | VCV ID |

### 5.3 Data Flow Between Systems

```
PharmVar (Allele Definitions)
    |
    v
CPIC (Clinical Guidelines)
    |
    v
PharmGKB (Annotations & Relationships)
    |
    v
FDA Labels (Regulatory Information)
```

### 5.4 Integration Points

| Source | Target | Relationship |
|--------|--------|--------------|
| PharmVar | CPIC | Allele nomenclature |
| PharmVar | PharmGKB | Star allele definitions |
| CPIC | PharmGKB | Guideline content |
| PharmGKB | FDA | Label annotations |
| dbSNP | All | Variant identifiers |
| RxNorm | All | Drug identifiers |

### 5.5 Common Data Exchange Formats

| Format | Use Case |
|--------|----------|
| JSON | API responses, structured data |
| TSV/CSV | Bulk downloads, spreadsheets |
| VCF | Variant calls |
| HGVS | Variant nomenclature |
| BioPAX | Pathway data |
| FHIR | Clinical implementation |

---

## References

### Primary Sources

1. PharmGKB: https://www.pharmgkb.org/ (https://www.clinpgx.org/)
2. CPIC: https://cpicpgx.org/
3. PharmVar: https://www.pharmvar.org/
4. FDA Pharmacogenomic Biomarkers: https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling

### Key Publications

1. Whirl-Carrillo M, et al. "An Evidence-Based Framework for Evaluating Pharmacogenomics Knowledge for Personalized Medicine." Clinical Pharmacology & Therapeutics. 2021. PMC8457105.

2. Whirl-Carrillo M, et al. "PharmGKB, an Integrated Resource of Pharmacogenomic Knowledge." Clinical Pharmacology & Therapeutics. 2021. PMC8650697.

3. Gaedigk A, et al. "PharmVar: A Global Resource and Repository for Pharmacogene Variation." Clinical Pharmacology & Therapeutics. 2021. PMC8725060.

4. Relling MV, Klein TE. "CPIC: Clinical Pharmacogenetics Implementation Consortium of the Pharmacogenomics Research Network." Clinical Pharmacology & Therapeutics. 2011.

### API Documentation

1. PharmGKB API: https://api.pharmgkb.org/swagger/
2. CPIC API Wiki: https://github.com/cpicpgx/cpic-data/wiki
3. PharmGKB Postman: https://www.postman.com/pharmgkb/pharmgkb-api/

### Data Downloads

1. PharmGKB Downloads: https://www.clinpgx.org/downloads
2. CPIC Data Releases: https://github.com/cpicpgx/cpic-data/releases
3. FDA Biomarkers PDF: https://www.fda.gov/media/124784/download
