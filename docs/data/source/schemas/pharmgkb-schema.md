# PharmGKB Schema Documentation

**Document ID:** SCHEMA-PHARMGKB
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** PharmGKB/ClinPGx API, CPIC API, S3 Downloads

---

## Overview

PharmGKB (Pharmacogenomics Knowledge Base) provides curated pharmacogenomic data including clinical annotations, drug labels, dosing guidelines, and variant-drug associations. The database has rebranded to ClinPGx but maintains the PharmGKB data model.

**API Base:** `https://api.pharmgkb.org/v1/` (redirects to ClinPGx)
**S3 Downloads:** `https://s3.pgkb.org/`
**Rate Limit:** 2 requests/second
**License:** CC BY-SA 4.0

---

## Database Statistics

Based on CPIC API and available documentation:

| Metric | Value | Source |
|--------|-------|--------|
| **CPIC Guidelines** | 33 guidelines | CPIC API |
| **Unique Genes** | ~25 pharmacogenes | CPIC guidelines |
| **Multi-gene Guidelines** | 8 guidelines | 2-6 genes each |
| **S3 File Categories** | 1000+ files | Various formats |

---

## API Endpoints

### Core Data Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/data/clinicalAnnotation` | GET | Clinical annotations |
| `/data/gene` | GET | Gene information |
| `/data/chemical` | GET | Drug/chemical data |
| `/data/variant` | GET | Genetic variants |
| `/data/label` | GET | Drug label annotations |
| `/data/guideline` | GET | Dosing guidelines |
| `/data/pathway` | GET | Pharmacokinetic pathways |
| `/data/automaticAnnotation` | GET | Automated literature annotations |

### Reference Data

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/data/haplotype` | GET | Star allele definitions |
| `/data/diplotype` | GET | Diplotype-phenotype mappings |
| `/data/phenotype` | GET | Phenotype definitions |

---

## Core Data Schemas

### Clinical Annotation

Clinical annotations link genetic variants to drug responses with evidence levels.

```json
{
  "id": "1449309937",
  "objCls": "ClinicalAnnotation",
  "name": "Annotation of CPIC Guideline for warfarin and CYP2C9, VKORC1, CYP4F2",
  "location": {
    "id": "1449309937",
    "_url": "/clinicalAnnotation/1449309937"
  },
  "relatedGenes": [
    {
      "id": "PA126",
      "name": "CYP2C9",
      "symbol": "CYP2C9"
    }
  ],
  "relatedChemicals": [
    {
      "id": "PA451906",
      "name": "warfarin"
    }
  ],
  "relatedVariants": [
    {
      "id": "PA166155091",
      "name": "rs1799853"
    }
  ],
  "evidenceLevel": "1A",
  "phenotypeCategories": [
    {
      "id": "PA166182818",
      "name": "Dosage"
    }
  ],
  "literatureCitations": [
    {
      "id": "15258624",
      "type": "PubMed"
    }
  ]
}
```

#### Clinical Annotation Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | Yes | PharmGKB accession (PA format) |
| `objCls` | string | Yes | Object class ("ClinicalAnnotation") |
| `name` | string | Yes | Annotation title |
| `location` | object | Yes | Self-reference URL |
| `relatedGenes` | array | Yes | Associated genes |
| `relatedChemicals` | array | Yes | Associated drugs/chemicals |
| `relatedVariants` | array | No | Associated variants |
| `evidenceLevel` | string | Yes | CPIC evidence level (1A-4) |
| `phenotypeCategories` | array | Yes | Effect categories |
| `literatureCitations` | array | No | Supporting references |

#### Evidence Levels

| Level | Description | Strength |
|-------|-------------|----------|
| `1A` | CPIC Level A - Strong evidence | Highest |
| `1B` | CPIC Level B - Moderate evidence | High |
| `2A` | Guideline annotation | Moderate |
| `2B` | PharmGKB annotation | Moderate |
| `3` | FDA label annotation | Supporting |
| `4` | Literature annotation | Lowest |

### Gene

```json
{
  "id": "PA126",
  "objCls": "Gene",
  "symbol": "CYP2C9",
  "name": "cytochrome P450 family 2 subfamily C member 9",
  "alternateName": ["CYPIIC9", "CYP2C10"],
  "chromosome": "10",
  "chromosomeSequenceId": "NC_000010.11",
  "location": {
    "_url": "/gene/PA126"
  },
  "hasNonStandardHaplotypes": false,
  "hasVariantAnnotation": true,
  "version": 75,
  "crossReferences": [
    {
      "resource": "HGNC",
      "resourceId": "HGNC:2623"
    },
    {
      "resource": "NCBI Gene",
      "resourceId": "1559"
    },
    {
      "resource": "Ensembl",
      "resourceId": "ENSG00000138109"
    }
  ]
}
```

#### Gene Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | Yes | PharmGKB ID (PA format) |
| `objCls` | string | Yes | "Gene" |
| `symbol` | string | Yes | HGNC symbol |
| `name` | string | Yes | Full gene name |
| `alternateName` | array | No | Alternative symbols |
| `chromosome` | string | Yes | Chromosome location |
| `chromosomeSequenceId` | string | No | RefSeq chromosome ID |
| `hasNonStandardHaplotypes` | boolean | Yes | Complex haplotype flag |
| `hasVariantAnnotation` | boolean | Yes | Variant data available |
| `version` | integer | Yes | Record version |
| `crossReferences` | array | Yes | External database links |

### Variant

```json
{
  "id": "PA166155091",
  "objCls": "Variant",
  "name": "rs1799853",
  "type": "snp",
  "gene": {
    "id": "PA126",
    "symbol": "CYP2C9"
  },
  "location": {
    "chromosome": "10",
    "position": 94942290,
    "assembly": "GRCh38"
  },
  "hgvs": {
    "genomic": "NC_000010.11:g.94942290C>T",
    "coding": "NM_000771.4:c.430C>T",
    "protein": "NP_000762.2:p.Arg144Cys"
  },
  "alleles": {
    "reference": "C",
    "alternate": ["T"]
  },
  "frequency": {
    "global": 0.12,
    "populations": [
      {"name": "European", "frequency": 0.15},
      {"name": "African", "frequency": 0.03}
    ]
  },
  "functionalStatus": "decreased function"
}
```

#### Variant Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | Yes | PharmGKB variant ID |
| `objCls` | string | Yes | "Variant" |
| `name` | string | Yes | rsID or other identifier |
| `type` | string | Yes | snp, indel, repeat, etc. |
| `gene` | object | No | Associated gene |
| `location` | object | Yes | Genomic coordinates |
| `hgvs` | object | No | HGVS nomenclature |
| `alleles` | object | Yes | Reference and alternate alleles |
| `frequency` | object | No | Population frequencies |
| `functionalStatus` | string | No | Functional classification |

### Chemical (Drug)

```json
{
  "id": "PA451906",
  "objCls": "Chemical",
  "name": "warfarin",
  "genericNames": ["warfarin sodium", "coumadin"],
  "tradeName": ["Coumadin", "Jantoven"],
  "drugClasses": [
    {
      "id": "PA164712363",
      "name": "Vitamin K antagonists"
    }
  ],
  "atcCodes": ["B01AA03"],
  "rxNormId": "11289",
  "inChIKey": "PJVWKTKQMONHTI-UHFFFAOYSA-N",
  "smiles": "CC(=O)CC(C1=CC=CC=C1)C2=C(O)C3=CC=CC=C3OC2=O",
  "types": ["Drug"],
  "hasRxAnnotation": true,
  "hasCpicGuideline": true,
  "hasFdaLabel": true
}
```

#### Chemical Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | Yes | PharmGKB chemical ID |
| `objCls` | string | Yes | "Chemical" |
| `name` | string | Yes | Primary name |
| `genericNames` | array | No | Generic names |
| `tradeName` | array | No | Brand names |
| `drugClasses` | array | No | Drug classifications |
| `atcCodes` | array | No | ATC classification codes |
| `rxNormId` | string | No | RxNorm identifier |
| `inChIKey` | string | No | InChI key |
| `smiles` | string | No | SMILES structure |
| `types` | array | Yes | "Drug", "Metabolite", etc. |
| `hasRxAnnotation` | boolean | Yes | Has prescribing info |
| `hasCpicGuideline` | boolean | Yes | Has CPIC guideline |
| `hasFdaLabel` | boolean | Yes | Has FDA label |

### Haplotype (Star Allele)

```json
{
  "id": "PA165816543",
  "objCls": "Haplotype",
  "name": "*2",
  "gene": {
    "id": "PA126",
    "symbol": "CYP2C9"
  },
  "function": "Decreased function",
  "activityScore": 0.5,
  "variants": [
    {
      "id": "PA166155091",
      "name": "rs1799853",
      "allele": "T"
    }
  ],
  "frequency": {
    "European": 0.13,
    "African": 0.03,
    "Asian": 0.02
  }
}
```

#### Haplotype Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | Yes | PharmGKB haplotype ID |
| `objCls` | string | Yes | "Haplotype" |
| `name` | string | Yes | Star allele name |
| `gene` | object | Yes | Associated gene |
| `function` | string | Yes | Functional classification |
| `activityScore` | number | No | CPIC activity score |
| `variants` | array | Yes | Defining variants |
| `frequency` | object | No | Population frequencies |

### Functional Classifications

| Classification | Activity Score | Description |
|----------------|---------------|-------------|
| `Normal function` | 1.0 | Wild-type activity |
| `Decreased function` | 0.5 | Reduced activity |
| `No function` | 0.0 | Complete loss |
| `Increased function` | 1.5+ | Enhanced activity |
| `Unknown function` | null | Not characterized |

---

## CPIC Guideline Schema

From CPIC API (`https://api.cpicpgx.org/v1/guideline`):

```json
{
  "id": "1",
  "version": 62,
  "name": "CYP2D6 and Codeine",
  "url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/",
  "pharmgkbid": ["PA166104996"],
  "genes": ["CYP2D6"],
  "clinpgxid": "CA123456",
  "notesonusage": null
}
```

### CPIC Guideline Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | Yes | CPIC guideline ID |
| `version` | integer | Yes | Version number |
| `name` | string | Yes | Gene-drug pair name |
| `url` | string | Yes | CPIC website URL |
| `pharmgkbid` | array | No | PharmGKB guideline IDs |
| `genes` | array | Yes | Associated genes (1-6) |
| `clinpgxid` | string | No | ClinPGx identifier |
| `notesonusage` | string | No | Special usage notes |

### Current CPIC Guidelines (33 Total)

| Category | Guidelines | Examples |
|----------|------------|----------|
| **Single-gene** | 25 | CYP2D6+codeine, CYP2C19+clopidogrel |
| **Multi-gene (2)** | 5 | CYP2C9+VKORC1+warfarin |
| **Multi-gene (3+)** | 3 | HLA-A+HLA-B+carbamazepine |
| **Multi-gene (6)** | 1 | Beta-blockers (6 genes) |

---

## Downloadable Data Files

### S3 Bucket Structure (s3.pgkb.org)

| Category | File Types | Examples |
|----------|------------|----------|
| **Clinical Evidence** | PDF, XLSX | Drug labels, regulatory docs |
| **Allele Definitions** | XLSX, TSV | Star allele tables |
| **Gene Resources** | XLSX, PNG | VIP summaries, pathways |
| **Regulatory Documents** | PDF | FDA, EMA, PMDA labels |
| **Molecular Data** | PNG, SVG | Structures, metabolite diagrams |

### Common Download Files

| File | Format | Description |
|------|--------|-------------|
| `clinicalAnnotations.zip` | ZIP/TSV | Clinical annotation data |
| `genes.zip` | ZIP/TSV | Gene information |
| `chemicals.zip` | ZIP/TSV | Drug data |
| `variants.zip` | ZIP/TSV | Variant information |
| `relationships.zip` | ZIP/TSV | Entity relationships |
| `pathways.zip` | ZIP/TSV | Pathway data |
| `dosingGuidelines.zip` | ZIP/TSV | CPIC dosing tables |
| `drugLabels.zip` | ZIP/TSV | FDA/EMA label data |

### TSV File Schema (clinicalAnnotations.tsv)

| Column | Type | Description |
|--------|------|-------------|
| `Clinical Annotation ID` | string | PharmGKB annotation ID |
| `Variant/Haplotypes` | string | Variant identifiers |
| `Gene` | string | Gene symbol |
| `Level of Evidence` | string | 1A, 1B, 2A, 2B, 3, 4 |
| `Level Override` | string | Manual evidence override |
| `Level Modifiers` | string | Evidence modifiers |
| `Score` | string | Annotation score |
| `Phenotype Category` | string | Dosage, Efficacy, Toxicity |
| `Drug(s)` | string | Associated drugs |
| `Phenotype(s)` | string | Associated phenotypes |
| `Variant Annotations` | string | Variant-specific annotations |
| `Variant Annotations (Continued)` | string | Overflow text |
| `PMIDs` | string | PubMed IDs |
| `Evidence Count` | integer | Supporting evidence count |
| `Related Chemicals` | string | Related drug IDs |
| `Related Diseases` | string | Disease associations |
| `Biogeographical Groups` | string | Population groups |
| `Chromosome` | string | Chromosome location |

### TSV File Schema (variants.tsv)

| Column | Type | Description |
|--------|------|-------------|
| `Variant ID` | string | PharmGKB variant ID |
| `Variant Name` | string | Common name (rsID) |
| `Gene IDs` | string | Associated gene IDs |
| `Gene Symbols` | string | Gene symbols |
| `Location` | string | Genomic coordinates |
| `Variant Annotation count` | integer | Annotation count |
| `Clinical Annotation count` | integer | Clinical annotation count |
| `Level 1/2 Clinical Annotation count` | integer | High-evidence count |
| `Guideline Annotation count` | integer | Guideline mentions |
| `Label Annotation count` | integer | Label mentions |
| `Synonyms` | string | Alternative identifiers |

---

## Cross-Reference Mappings

### External Database Links

| Database | ID Pattern | Example |
|----------|------------|---------|
| **dbSNP** | rs[0-9]+ | rs1799853 |
| **ClinVar** | RCV[0-9]+ | RCV000018137 |
| **HGNC** | HGNC:[0-9]+ | HGNC:2623 |
| **NCBI Gene** | [0-9]+ | 1559 |
| **Ensembl** | ENSG[0-9]+ | ENSG00000138109 |
| **UniProt** | [A-Z0-9]+ | P11712 |
| **RxNorm** | [0-9]+ | 11289 |
| **ATC** | [A-Z][0-9]{2}[A-Z]{2}[0-9]{2} | B01AA03 |
| **PubChem** | [0-9]+ | 54678486 |
| **ChEBI** | CHEBI:[0-9]+ | CHEBI:10033 |
| **DrugBank** | DB[0-9]+ | DB00682 |

### PharmGKB ID Patterns

| Entity | Pattern | Example |
|--------|---------|---------|
| Gene | PA[0-9]+ | PA126 |
| Chemical | PA[0-9]+ | PA451906 |
| Variant | PA[0-9]+ | PA166155091 |
| Haplotype | PA[0-9]+ | PA165816543 |
| Clinical Annotation | [0-9]+ | 1449309937 |
| Guideline | PA[0-9]+ | PA166104996 |

---

## API Response Format

### Success Response

```json
{
  "status": "success",
  "data": {
    // Entity data
  },
  "metadata": {
    "version": "1.0",
    "timestamp": "2026-01-21T12:00:00Z"
  }
}
```

### Error Response

```json
{
  "status": "error",
  "code": 400,
  "message": "Invalid parameter",
  "details": "Gene symbol not found"
}
```

### Pagination

```json
{
  "status": "success",
  "data": [],
  "pagination": {
    "offset": 0,
    "limit": 25,
    "total": 150,
    "hasMore": true
  }
}
```

---

## Integration Examples

### Python: Fetch Clinical Annotations

```python
import requests
import time

BASE_URL = "https://api.pharmgkb.org/v1"
RATE_LIMIT = 0.5  # 2 requests/second

def get_clinical_annotations(gene_symbol: str) -> list:
    """Fetch clinical annotations for a gene."""
    url = f"{BASE_URL}/data/clinicalAnnotation"
    params = {"gene": gene_symbol}

    response = requests.get(url, params=params)
    response.raise_for_status()

    time.sleep(RATE_LIMIT)  # Rate limiting
    return response.json().get("data", [])

# Example
annotations = get_clinical_annotations("CYP2D6")
for ann in annotations:
    print(f"{ann['name']}: Level {ann['evidenceLevel']}")
```

### Python: Parse TSV Downloads

```python
import csv
import zipfile
from io import TextIOWrapper

def parse_clinical_annotations_zip(filepath: str):
    """Parse PharmGKB clinical annotations ZIP file."""
    with zipfile.ZipFile(filepath, 'r') as zf:
        with zf.open('clinical_annotations.tsv') as f:
            reader = csv.DictReader(
                TextIOWrapper(f, 'utf-8'),
                delimiter='\t'
            )
            for row in reader:
                yield {
                    'id': row['Clinical Annotation ID'],
                    'gene': row['Gene'],
                    'drug': row['Drug(s)'],
                    'evidence_level': row['Level of Evidence'],
                    'phenotype': row['Phenotype Category'],
                    'pmids': row['PMIDs'].split(';') if row['PMIDs'] else []
                }
```

### Python: CPIC API Integration

```python
import requests

def get_cpic_guidelines() -> list:
    """Fetch all CPIC guidelines."""
    url = "https://api.cpicpgx.org/v1/guideline"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

def get_guidelines_for_gene(gene: str) -> list:
    """Filter guidelines by gene."""
    all_guidelines = get_cpic_guidelines()
    return [g for g in all_guidelines if gene in g['genes']]

# Example
cyp2d6_guidelines = get_guidelines_for_gene("CYP2D6")
print(f"Found {len(cyp2d6_guidelines)} CYP2D6 guidelines")
```

---

## PharmCAT Integration

PharmCAT (Pharmacogenomics Clinical Annotation Tool) processes VCF files using PharmGKB data.

### PharmCAT Data Files

| File | Description |
|------|-------------|
| `pharmcat_positions.vcf` | Curated PGx variant positions |
| `pharmcat_regions.bed` | Gene regions for analysis |

### PharmCAT Output Schema

```json
{
  "metaData": {
    "pharmcatVersion": "2.x.x",
    "timestamp": "2026-01-21T12:00:00Z"
  },
  "genotypeCalls": [
    {
      "gene": "CYP2D6",
      "diplotype": "*1/*4",
      "phenotype": "Intermediate Metabolizer",
      "activityScore": 1.0
    }
  ],
  "recommendations": [
    {
      "drug": "codeine",
      "recommendation": "Use label-recommended dosage",
      "source": "CPIC",
      "strength": "Strong"
    }
  ]
}
```

---

## Implementation Details

### Star Allele Naming Conventions (PharmVar Standard)

**Format:** `GENE*Allele.Suballele`

| Component | Description | Example | Notes |
|-----------|-------------|---------|-------|
| Gene symbol | HGNC official symbol | CYP2D6, CYP2C9 | Required, immutable |
| Asterisk | Separator | * | Standard notation |
| Major allele | Primary allele number | 1, 2, 4, 17 | Ordered by discovery |
| Suballele | Variant-specific | .001, .002 | Optional, precise definition |
| Copy number | Gene duplication | x2, x3, xN | For CNV alleles |

**Examples:**
- `CYP2D6*1` - Reference allele (normal function)
- `CYP2D6*4` - Core allele with decreased function
- `CYP2D6*4.001` - Specific suballele variant
- `CYP2D6*4x2` - Duplication of *4 allele
- `CYP2C9*2` - Single base substitution allele

**PharmVar ID Format:** `PV#####` (immutable, version-tracked)

### CPIC Guideline Structure

#### Guideline Organization

| Element | Type | Description | Notes |
|---------|------|-------------|-------|
| **ID** | UUID | Unique guideline identifier | From CPIC API |
| **Name** | String | Gene-drug pair name | e.g., "CYP2D6 and Codeine" |
| **Version** | Integer | Release version | Updated quarterly |
| **Genes** | Array | Associated genes | 1-6 genes per guideline |
| **Drugs** | Array | Associated drugs | CPIC or FDA contraindicated drugs |
| **URL** | String | CPIC website reference | Published guideline page |
| **Status** | String | A, A/B, B, C, or D | Prescribing action strength |

#### Guideline Categories (33 Total)

| Category | Count | Examples |
|----------|-------|----------|
| **Single-gene** | 25 | CYP2D6+codeine, CYP2C19+clopidogrel, TPMT+azathioprine |
| **Multi-gene (2)** | 5 | CYP2C9+VKORC1+warfarin, CYP3A4+CYP3A5+tacrolimus |
| **Multi-gene (3+)** | 3 | HLA-A+HLA-B+carbamazepine |
| **Multi-gene (6)** | 1 | Beta-blockers (ADRB1, CYP2B6, CYP2C19, CYP2D6, CYP3A4, CYP3A5) |

#### Recommendation Lookup Methods

| Method | Description | Lookup Key Format | Use Case |
|--------|-------------|-------------------|----------|
| **PHENOTYPE** | Phenotype-based recommendations | `{"CYP2D6": "Poor Metabolizer"}` | Most common for CYP enzymes |
| **ACTIVITY_SCORE** | Activity score thresholds | `{"CYP2D6": "0"}` | Precision dosing calculations |
| **ALLELE_STATUS** | Specific allele presence/absence | `{"HLA-B": "*57:01 positive"}` | HLA-based contraindications |

### Phenotype Mapping Examples

#### Metabolizer Phenotypes (CYP Enzymes)

| Phenotype | Activity Score | Diplotype Example | Functional Impact | Clinical Example |
|-----------|----------------|--------------------|-------------------|------------------|
| **Ultrarapid** | >1.5 | *1x2/*1 | Enhanced enzyme activity | Increased drug clearance, subtherapeutic levels |
| **Rapid** | 1.25-1.5 | *1/*1x2 | Moderately enhanced | Normal to elevated clearance |
| **Normal** | 1.0 | *1/*1 | Full enzyme activity | Standard dosing appropriate |
| **Intermediate** | 0.5-0.75 | *1/*2 | Reduced activity | Moderate dose reduction may be needed |
| **Poor** | 0 | *3/*4 | Complete loss | Avoid drug or significantly reduce dose |

#### Activity Score Calculation

**Formula:** Combine allele 1 function + allele 2 function

| Allele Function | Activity Value |
|-----------------|-----------------|
| Normal function | 1.0 |
| Decreased function | 0.5 |
| No function | 0.0 |
| Increased function | >1.0 |
| Unknown/Uncertain | null |

**Diplotype Scoring Examples:**
- `*1/*1` = 1.0 + 1.0 = 2.0 (Normal/Rapid) → Normal Metabolizer
- `*1/*2` = 1.0 + 0.5 = 1.5 (Normal/Decreased) → Intermediate Metabolizer
- `*2/*4` = 0.5 + 0.0 = 0.5 (Decreased/No) → Poor Metabolizer

#### HLA Allele Phenotypes (Transporters/Receptors)

| Phenotype | Typical Alleles | Activity | Example Drug |
|-----------|-----------------|----------|--------------|
| **Positive** | HLA-B*57:01 | N/A | Abacavir contraindicated |
| **Negative** | Any other HLA-B | N/A | Abacavir acceptable |
| **Normal** | Standard alleles | Expected | Default phenotype |
| **Reduced** | Non-canonical alleles | Decreased | Special consideration |

### Allele Definition Structure

#### CPIC Allele Definition Table Format

| Column | Value Type | Description | Example |
|--------|-----------|-------------|---------|
| **Allele Name** | String | Star allele identifier | *1, *2, *4 |
| **Variants** | Array[String] | Defining rsIDs or positions | rs1065852, rs5030655 |
| **Function** | String | Functional category | "Normal function", "No function" |
| **Activity Score** | Float | Numeric activity value | 1.0, 0.5, 0 |
| **PharmVar ID** | String | Immutable variant identifier | PV00123 |
| **Copy Number** | Integer | Gene copy count | 1, 2, 3+ |

#### Variant Specification Example

| Allele | rs1065852 | rs5030655 | rs3892097 | Function | Score |
|--------|-----------|-----------|-----------|----------|-------|
| *1 | C | A | G | Normal | 1.0 |
| *2 | T | A | G | Decreased | 0.5 |
| *4 | C | A | A | No | 0.0 |
| *6 | C | del | G | No | 0.0 |
| *10 | T | A | G | Decreased | 0.5 |

### CPIC Evidence Levels and Scoring

| Level | Name | Score Range | Strength | Action |
|-------|------|-------------|----------|--------|
| **A** | Strong | ≥80 | Highest evidence | Prescribing action recommended |
| **B** | Moderate | 25-79 | Moderate evidence | Prescribing action may be recommended |
| **C** | Supportive | 8-24 | Supporting evidence | No prescribing action recommended |
| **D** | Insufficient | <8 | Low/no evidence | Insufficient evidence |

#### Evidence Scoring Factors

| Factor | Category | Points |
|--------|----------|--------|
| **Phenotype Type** | Clinical outcomes | 1.0 (highest weight) |
|  | Pharmacokinetic | 0.5 |
| **P-value** | p < 0.01 | 1.0 |
|  | p ≤ 0.05 | 0.5 |
| **Cohort Size** | > 500 | 1.0 |
|  | ≤ 50 | 0.0 |
| **Effect Size** | OR/HR > 2 or < 0.5 | 0.5 |
| **Study Type** | In vivo studies | 1.0 |
|  | In vitro studies | 0.0 |

---

## References

1. PharmGKB: https://www.pharmgkb.org/
2. CPIC: https://cpicpgx.org/
3. PharmCAT: https://pharmcat.org/
4. PharmVar: https://www.pharmvar.org/
5. Whirl-Carrillo M, et al. (2021). An Evidence-Based Framework for Evaluating Pharmacogenomics Knowledge for Personalized Medicine. Clin Pharmacol Ther. 110(3):563-572.

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema from API and CPIC data |
