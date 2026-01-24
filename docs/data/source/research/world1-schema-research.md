---
id: research-world1-schema
title: WORLD 1 Modern Genetics - Schema & Sample Data Research
type: research
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [research, genetics, schema, dbsnp, clinvar, gnomad, pharmgkb, snpedia]
---

**Parent:** [Research](./_index.md)

# WORLD 1: Modern Genetics - Schema & Sample Data Research

**Document ID:** WORLD1-SCHEMA-RESEARCH
**Status:** Research Complete
**Owner:** Research Agent
**Date:** January 21, 2026
**Version:** 1.0

---

## Executive Summary

This research report documents schema availability, sample data access, and integration value for five priority databases in WORLD 1 (Modern Genetics). All five databases have PUBLIC schemas or sample data available, making them suitable for unified schema development.

### Quick Reference: Database Rankings

| Rank | Database | Schema Availability | Sample Data | API Access | Integration Value |
|------|----------|---------------------|-------------|------------|-------------------|
| **1** | **dbSNP (NCBI)** | Excellent | Yes | Yes | Critical |
| **2** | **ClinVar** | Excellent | Yes | Yes | Critical |
| **3** | **gnomAD** | Good | Yes | Yes | High |
| **4** | **PharmGKB** | Good | Yes | Yes | High |
| **5** | **SNPedia** | Moderate | Limited | Yes | Medium |

---

## 1. dbSNP (NCBI)

### Overview
- **URL:** https://www.ncbi.nlm.nih.gov/snp/
- **Content:** 1B+ single nucleotide polymorphisms
- **License:** Public Domain
- **Priority:** Tier 1 (Critical for MVP)

### Schema Documentation

**Status: EXCELLENT - Full OpenAPI specification available**

The dbSNP schema is formally documented through the NCBI Variation Services API with a complete OpenAPI specification.

**Schema Location:**
- OpenAPI YAML: https://api.ncbi.nlm.nih.gov/variation/v0/var_service.yaml
- Interactive Docs: https://api.ncbi.nlm.nih.gov/variation/v0/

### Core Data Model

```yaml
RefSNP_Snapshot:
  description: "Each JSON object represents a single RefSNP"
  properties:
    rsid: string  # RS identifier (e.g., rs123456)
    create_date: datetime
    last_update: datetime
    citations: array[PMID]

SPDI_Representation:
  description: "Sequence, Position, Deletion, Insertion format"
  properties:
    seq_id: string        # Sequence identifier (NC_000001.11)
    position: integer     # 0-based position
    deleted_sequence: string
    inserted_sequence: string

Allele_Annotation:
  properties:
    clinical_significance: enum
    frequency_data: object
    gene_annotations: array
    protein_consequences: array

Population_Frequency:
  properties:
    bioproject_id: string
    population_name: string
    allele_count: integer
    genotype_counts: object
    hardy_weinberg_pvalue: float
```

### Key Entity Types

1. **RefSNP** - Primary variant identifier
2. **SPDI** - Normalized sequence representation
3. **Allele** - Alternative sequences at a position
4. **Population** - Frequency data by population
5. **Clinical Annotation** - Significance and citations

### Sample Data Access

| Resource | Location | Format | Size |
|----------|----------|--------|------|
| Sample JSON | https://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-sample.json.bz2 | JSON | ~1 MB |
| Demo Script | https://ftp.ncbi.nih.gov/snp/latest_release/JSON/rsjson_demo.py | Python | 10 KB |
| Chr1 Data | https://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chr1.json.bz2 | JSON | ~5 GB |
| VCF (GRCh38) | https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz | VCF | 28 GB |

### VCF INFO Fields

```
CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO

INFO Tags:
- NSM: Missense variant
- SYN: Synonymous variant
- NSN: Nonsense variant
- CAF: Allele frequencies (1000 Genomes Phase III)
- G5: MAF >5% in general populations
- VLD: Validated variants
- COMMON: MAF >= 0.01
- OM: Has OMIM information
```

### API Response Format

```json
{
  "refsnp_id": "rs123456",
  "primary_snapshot_data": {
    "placements_with_allele": [{
      "seq_id": "NC_000001.11",
      "position": 12345678,
      "deleted_sequence": "A",
      "inserted_sequence": "G"
    }],
    "allele_annotations": [{
      "frequency": [{
        "study_name": "1000Genomes",
        "allele_count": 500,
        "total_count": 5000
      }],
      "clinical": {
        "clinical_significance": ["pathogenic"]
      }
    }]
  }
}
```

### Integration Value

- **Cross-references:** Links to ClinVar, OMIM, PubMed
- **Population data:** 12 populations via ALFA integration
- **Clinical overlap:** 960K+ ClinVar RS IDs
- **Rate Limit:** 1 request/second recommended

---

## 2. ClinVar

### Overview
- **URL:** https://www.ncbi.nlm.nih.gov/clinvar/
- **Content:** 2M+ clinical variant interpretations
- **License:** Public Domain
- **Priority:** Tier 1 (Critical for MVP)

### Schema Documentation

**Status: EXCELLENT - XSD schemas + comprehensive documentation**

ClinVar provides formal XML Schema Definition (XSD) files and detailed field documentation.

**Schema Locations:**
- Public XSD: https://ftp.ncbi.nih.gov/pub/clinvar/xsd_public/
- Submission XSD: https://ftp.ncbi.nih.gov/pub/clinvar/clinvar_submission.xsd
- README: https://ftp.ncbi.nih.gov/pub/clinvar/README.txt

### Core Data Model

```yaml
# Three-tier accession model
Accession_Types:
  SCV: "Submitted Clinical Variant - individual submission"
  RCV: "Reference ClinVar - variant + condition aggregation"
  VCV: "Variation in ClinVar - variant-level aggregation"

Variation:
  properties:
    variation_id: integer
    name: string  # HGVS nomenclature
    variation_type: enum[SNV, deletion, insertion, indel, duplication]

ClinicalSignificance:
  # 2025 update: Split into three classification types
  germline_classification: enum[pathogenic, likely_pathogenic, uncertain, likely_benign, benign]
  somatic_clinical_impact: enum  # For cancer variants
  oncogenicity_classification: enum  # For oncogenes

ReviewStatus:
  values:
    - "practice guideline"
    - "reviewed by expert panel"
    - "criteria provided, multiple submitters, no conflicts"
    - "criteria provided, conflicting interpretations"
    - "criteria provided, single submitter"
    - "no assertion criteria provided"

SubmittedRecord:
  properties:
    scv_accession: string
    submitter: string
    classification: ClinicalSignificance
    condition: string  # Disease/phenotype
    evidence: array[Citation]
    last_evaluated: date
```

### Key Entity Types & Relationships

```
Variant (VCV) ──1:N──> Variant-Condition Pair (RCV) ──1:N──> Submission (SCV)
     │                         │                              │
     └── Genomic Location      └── Disease (HPO/OMIM)        └── Evidence
     └── Gene                  └── Clinical Significance     └── Submitter
     └── HGVS                  └── Review Status             └── Method
```

### Sample Data Access

| Resource | Location | Format | Size |
|----------|----------|--------|------|
| variant_summary.txt | https://ftp.ncbi.nih.gov/pub/clinvar/tab_delimited/ | TSV | 421 MB |
| submission_summary.txt | https://ftp.ncbi.nih.gov/pub/clinvar/tab_delimited/ | TSV | 359 MB |
| VCV XML (latest) | https://ftp.ncbi.nih.gov/pub/clinvar/xml/ClinVarVCVRelease_00-latest.xml.gz | XML | 5.4 GB |
| VCF (GRCh38) | https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/ | VCF | ~500 MB |
| Weekly updates | https://ftp.ncbi.nih.gov/pub/clinvar/xml/weekly_release/ | XML | Variable |

### Tab-Delimited Fields (variant_summary.txt)

```
AlleleID, Type, Name, GeneID, GeneSymbol, HGNC_ID, ClinicalSignificance,
ClinSigSimple, LastEvaluated, RS# (dbSNP), nsv/esv (dbVar), RCVaccession,
PhenotypeIDS, PhenotypeList, Origin, OriginSimple, Assembly, ChromosomeAccession,
Chromosome, Start, Stop, ReferenceAllele, AlternateAllele, Cytogenetic,
ReviewStatus, NumberSubmitters, Guidelines, TestedInGTR, OtherIDs,
SubmitterCategories, VariationID, PositionVCF, ReferenceAlleleVCF,
AlternateAlleleVCF
```

### API Access

```bash
# E-utilities access
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=BRCA1

# Variation Services API
GET https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid}
```

### Integration Value

- **Phenotype links:** HPO, OMIM, MedGen, MeSH
- **Gene links:** NCBI Gene, HGNC
- **Variant links:** dbSNP (RS#), dbVar (SV)
- **Evidence:** PubMed citations, submitter methods
- **Update frequency:** Weekly (Mondays)

---

## 3. gnomAD (Genome Aggregation Database)

### Overview
- **URL:** https://gnomad.broadinstitute.org/
- **Content:** Variants from 807,162 individuals (v4.0)
- **License:** Open Access
- **Priority:** Tier 1 (Critical for population frequencies)

### Schema Documentation

**Status: GOOD - VCF headers + code documentation**

gnomAD schema is documented through VCF header specifications, Hail table structures, and the gnomad_methods Python library.

**Schema Locations:**
- VCF Utils Docs: https://broadinstitute.github.io/gnomad_methods/api_reference/utils/vcf.html
- Hail Tables: Available on AWS/GCP/Azure

### Core Data Model

```yaml
Variant:
  properties:
    chrom: string
    pos: integer
    ref: string
    alt: string
    variant_id: string  # chrom-pos-ref-alt format

AlleleFrequency:
  properties:
    AC: integer        # Allele count
    AN: integer        # Allele number (total chromosomes)
    AF: float          # Allele frequency
    nhomalt: integer   # Number of homozygotes

PopulationStratification:
  genetic_ancestry_groups:
    - afr: "African/African-American"
    - amr: "Admixed American"
    - asj: "Ashkenazi Jewish"
    - eas: "East Asian"
    - fin: "Finnish"
    - mid: "Middle Eastern"  # New in v4
    - nfe: "Non-Finnish European"
    - oth: "Other"
    - sas: "South Asian"

QualityAnnotations:
  allele_specific:
    - AS_FS: "Fisher strand bias"
    - AS_MQ: "Mapping quality"
    - AS_QD: "Quality by depth"
    - AS_SOR: "Strand odds ratio"
  site_level:
    - FS, MQ, QD, SOR, VarDP

InSilicoAnnotations:
  properties:
    cadd_phred: float      # CADD Phred-scaled score (1-99)
    cadd_raw_score: float
    revel_score: float
    spliceai_ds_max: float
    pangolin_largest_ds: float
    phylop: float          # Conservation score
```

### VCF INFO Fields

```
Standard fields:
AC, AN, AF, nhomalt, AC_XX, AC_XY

Population-specific (per genetic ancestry):
AC_afr, AN_afr, AF_afr, AC_nfe, AN_nfe, AF_nfe, etc.

Filtering allele frequency (FAF):
faf95, faf99 (95th/99th percentile)

Quality flags:
lcr (low complexity region)
segdup (segmental duplication)
nonpar (non-pseudoautosomal region)

In silico predictions:
cadd_phred, revel, spliceai, pangolin, phylop
```

### Sample Data Access

| Resource | Location | Format | Size |
|----------|----------|--------|------|
| VCF (exomes v4) | gs://gcp-public-data--gnomad/release/4.0/vcf/exomes/ | VCF | ~200 GB |
| VCF (genomes v4) | gs://gcp-public-data--gnomad/release/4.0/vcf/genomes/ | VCF | ~1 TB |
| Hail Tables | AWS/GCP/Azure (see downloads page) | Hail MT | Variable |
| Azure Parquet | https://datasetgnomad.blob.core.windows.net/ | Parquet | ~30 TB |

### GraphQL API

```graphql
query {
  variant(variantId: "1-55516888-G-A", dataset: GNOMAD_R4) {
    variantId
    chrom
    pos
    ref
    alt
    exome {
      ac
      an
      af
      populations {
        id
        ac
        an
        af
      }
    }
    genome {
      ac
      an
      af
    }
  }
}
```

### Integration Value

- **Population diversity:** 8 genetic ancestry groups + sex stratification
- **Quality metrics:** Comprehensive QC annotations
- **In silico scores:** CADD, REVEL, SpliceAI, Pangolin included
- **VRS annotations:** GA4GH Variant Representation Spec compliant
- **Cloud-native:** Available on AWS, GCP, Azure

---

## 4. PharmGKB (Pharmacogenomics Knowledge Base)

### Overview
- **URL:** https://www.pharmgkb.org/ (redirects to https://www.clinpgx.org/)
- **Content:** 715 drugs, 1,761 genes, 227 diseases
- **License:** CC BY-SA 4.0
- **Priority:** Tier 1 (Critical for pharmacogenomics)

### Schema Documentation

**Status: GOOD - Swagger API + entity documentation**

PharmGKB provides RESTful API documentation through Swagger, with entity relationships documented in publications.

**Schema Locations:**
- Swagger Docs: https://api.clinpgx.org/swagger/
- Data Description: API v1 documentation
- Rate Limit: 2 requests/second max

### Core Data Model

```yaml
Gene:
  properties:
    id: string           # PharmGKB ID (PA123)
    symbol: string       # HGNC symbol
    name: string
    chromosome: string
    haplotypes: array[Haplotype]
    variants: array[Variant]

Drug:
  properties:
    id: string           # PharmGKB ID (PA448XXX)
    name: string
    generic_names: array[string]
    trade_names: array[string]
    drugbank_id: string
    atc_codes: array[string]

Variant:
  properties:
    id: string           # rsID or PharmGKB ID
    gene: Gene
    location: string     # HGVS
    type: enum[SNP, haplotype, CNV, indel]

ClinicalAnnotation:
  description: "Aggregated variant-drug-phenotype associations"
  properties:
    id: string
    variant: Variant
    drug: Drug
    phenotype: Phenotype
    level_of_evidence: enum[1A, 1B, 2A, 2B, 3, 4]
    clinical_significance: string
    summary: text

VariantAnnotation:
  description: "Individual literature findings"
  properties:
    publication: PMID
    variant: Variant
    drug: Drug
    phenotype: Phenotype
    study_parameters:
      population: string
      study_size: integer
      p_value: float
      odds_ratio: float

DosingGuideline:
  properties:
    gene: Gene
    drug: Drug
    recommendation: text
    source: enum[CPIC, DPWG, FDA]
    strength: enum[strong, moderate, optional]
```

### Key Entity Relationships

```
Gene ──1:N──> Variant ──N:M──> Drug
  │              │              │
  └── Haplotype  │              └── Drug Label
                 │
           VariantAnnotation ──> ClinicalAnnotation
                 │                      │
                 └── Publication        └── Level of Evidence
                                        └── DosingGuideline
```

### Levels of Evidence

| Level | Description |
|-------|-------------|
| **1A** | Variant-drug combinations with CPIC or DPWG dosing guideline |
| **1B** | Variant-drug in FDA-required PGx testing |
| **2A** | Variant-drug with strong evidence + annotation in drug label |
| **2B** | Moderate evidence with replicated studies |
| **3** | Single significant study OR multiple non-significant studies |
| **4** | Case reports, in vitro, or functional assay evidence only |

### Sample Data Access

| Resource | Description | Format | Access |
|----------|-------------|--------|--------|
| Clinical Annotations | Curated gene-drug associations | TSV | Download page |
| Drug Labels | FDA/EMA label annotations | TSV | Download page |
| Dosing Guidelines | CPIC/DPWG recommendations | TSV | Download page |
| Pathways | PK/PD pathway diagrams | PDF/SVG | Web interface |
| Haplotypes | Star allele definitions | TSV | Gene pages |

### API Endpoints

```
GET /v1/data/gene/{id}
GET /v1/data/drug/{id}
GET /v1/data/variant/{id}
GET /v1/data/clinicalAnnotation
GET /v1/data/dosingGuideline
GET /v1/data/drugLabel
GET /v1/data/pathway
```

### Integration Value

- **Clinical actionability:** CPIC guidelines implementation-ready
- **Drug labels:** FDA, EMA, PMDA, HCSC, SwissMedic annotations
- **Cross-references:** DrugBank, ATC, RxNorm, ChEBI
- **Pathways:** Pharmacokinetic and pharmacodynamic visualization
- **Evidence hierarchy:** Clear LOE system for prioritization

---

## 5. SNPedia

### Overview
- **URL:** https://www.snpedia.com/
- **Content:** SNP interpretations (wiki-based)
- **License:** CC BY-NC-SA 3.0
- **Priority:** Tier 2 (Valuable for consumer interpretations)
- **Owner:** MyHeritage (acquired September 2019)

### Schema Documentation

**Status: MODERATE - MediaWiki API + semantic templates**

SNPedia uses MediaWiki with Semantic MediaWiki extensions. Schema is implicit in page templates.

**Schema Locations:**
- MediaWiki API: http://bots.snpedia.com/api.php
- R Package: https://bioconductor.org/packages/SNPediaR
- Python Tools: wikitools library

### Core Data Model

```yaml
SNP_Page:
  template_fields:
    rsid: string           # e.g., rs123456
    chromosome: string
    position: integer
    orientation: enum[plus, minus]
    reference: string
    gene: string
    summary: text

Genotype_Page:
  template_fields:
    rsid: string
    genotype: string       # e.g., (A;A), (A;G), (G;G)
    magnitude: float       # 0-10 scale of importance
    repute: enum[good, bad, neutral]
    summary: text

GeneInfo:
  properties:
    symbol: string
    name: string
    description: text
    associated_snps: array[rsid]

Genoset:
  description: "Combinations of multiple SNPs"
  properties:
    name: string
    criteria: array[Genotype]
    magnitude: float
    summary: text
```

### Key Entity Types

1. **SNP Pages** - Individual rs# documentation
2. **Genotype Pages** - Per-genotype interpretations with magnitude/repute
3. **Gene Pages** - Gene-level summaries
4. **Genosets** - Multi-SNP combinations
5. **Conditions** - Disease/trait associations

### Magnitude Scale

| Magnitude | Interpretation |
|-----------|----------------|
| 0 | Common genotype, nothing special |
| 1-2 | Slightly noteworthy |
| 3-4 | Interesting |
| 5-6 | Moderately significant |
| 7-8 | Significant |
| 9-10 | Very significant (rare) |

### Sample Data Access

| Resource | Description | Format | Access |
|----------|-------------|--------|--------|
| API Access | MediaWiki Action API | JSON/XML | bots.snpedia.com/api.php |
| R Package | SNPediaR (Bioconductor) | R | bioconductor.org |
| GFF Dump | Semi-regular exports | GFF | Web interface |
| Bulk Export | All page titles/content | Wiki XML | Special:Export |

### API Access Examples

```python
# Python with wikitools
import wikitools

site = wikitools.Wiki("http://bots.snpedia.com/api.php")
page = wikitools.Page(site, "Rs1234567")
text = page.getWikiText()
```

```r
# R with SNPediaR
library(SNPediaR)
pages <- getPages(c("rs123456", "rs789012"))
```

```bash
# MediaWiki API query
curl "http://bots.snpedia.com/api.php?action=query&titles=Rs1234567&prop=revisions&rvprop=content&format=json"
```

### OAuth Access (2016+)

Token creation: http://bots.snpedia.com/index.php/Special:OAuthConsumerRegistration

### Integration Value

- **Consumer-friendly:** Plain language interpretations
- **Magnitude scoring:** Quick importance assessment
- **Repute system:** Good/bad/neutral classification
- **Genosets:** Multi-SNP trait combinations
- **Wiki format:** Community-maintained, evolving

### Limitations

- Non-commercial license (CC BY-NC-SA 3.0)
- Wiki format requires parsing
- No formal API rate limits documented
- Owned by MyHeritage (potential access changes)

---

## Unified Schema Recommendations

### Common Entity Model

Based on the five databases analyzed, the following unified entity model is recommended:

```yaml
UnifiedVariant:
  identifiers:
    rsid: string              # dbSNP RS#
    clinvar_id: integer       # ClinVar VariationID
    gnomad_id: string         # chrom-pos-ref-alt
    pharmgkb_id: string       # PA ID
    snpedia_page: string      # Wiki page name

  location:
    chromosome: string
    position: integer         # GRCh38
    reference: string
    alternate: string
    hgvs_genomic: string
    hgvs_coding: string
    hgvs_protein: string

  frequency:
    gnomad_af: float
    gnomad_popmax_af: float
    population_frequencies: map[population, float]

  clinical:
    clinvar_significance: string
    clinvar_review_status: string
    clinvar_conditions: array[string]
    snpedia_magnitude: float
    snpedia_repute: string

  pharmacogenomics:
    pharmgkb_level: string    # 1A, 1B, 2A, etc.
    affected_drugs: array[Drug]
    dosing_recommendations: array[Guideline]

  functional:
    gene_symbol: string
    consequence: string       # VEP consequence
    cadd_score: float
    revel_score: float
    spliceai_score: float
```

### Cross-Reference Mapping

| Field | dbSNP | ClinVar | gnomAD | PharmGKB | SNPedia |
|-------|-------|---------|--------|----------|---------|
| Variant ID | RS# | VCV/RCV | variant_id | PA ID | rs# |
| Gene | gene_id | GeneID | gene_symbol | gene.id | gene |
| Position | pos | Start/Stop | pos | location | position |
| Reference | REF | ReferenceAllele | ref | - | reference |
| Alternate | ALT | AlternateAllele | alt | - | - |
| Frequency | CAF | - | AF | - | - |
| Clinical | clinical | ClinicalSignificance | - | LOE | magnitude |
| Condition | - | PhenotypeList | - | phenotype | - |

### API Integration Strategy

```yaml
Priority_1_APIs:
  - NCBI_Variation_Services:
      covers: [dbSNP, ClinVar]
      rate_limit: 1/sec
      format: JSON

  - gnomAD_GraphQL:
      covers: [gnomAD]
      rate_limit: "reasonable use"
      format: GraphQL/JSON

Priority_2_APIs:
  - PharmGKB_REST:
      covers: [PharmGKB]
      rate_limit: 2/sec
      format: JSON

  - SNPedia_MediaWiki:
      covers: [SNPedia]
      rate_limit: "polite crawling"
      format: JSON/wikitext

Bulk_Downloads_Recommended:
  - dbSNP: JSON by chromosome
  - ClinVar: variant_summary.txt (421 MB)
  - gnomAD: VCF per chromosome (large)
  - PharmGKB: All TSV downloads
  - SNPedia: GFF dump + API crawl
```

---

## Data Size Estimates

| Database | Bulk Download | Processed/Indexed | Update Frequency |
|----------|---------------|-------------------|------------------|
| dbSNP | ~60 GB (JSON) | ~20 GB | Monthly |
| ClinVar | ~6 GB (XML) | ~2 GB | Weekly |
| gnomAD | ~1 TB (VCF) | ~50 GB (summary) | Major releases |
| PharmGKB | ~500 MB | ~200 MB | Continuous |
| SNPedia | ~1 GB (estimated) | ~500 MB | Continuous |
| **Total** | **~1.1 TB** | **~73 GB** | Variable |

---

## Implementation Priority

### Week 1-2: Core Infrastructure
1. **dbSNP** - Load JSON schema, sample data, build RS# lookup
2. **ClinVar** - Parse variant_summary.txt, establish clinical significance mapping

### Week 3-4: Population Data
3. **gnomAD** - Integrate VCF parsing, population frequency extraction
4. Cross-reference dbSNP RS# to gnomAD variant_id

### Week 5-6: Pharmacogenomics
5. **PharmGKB** - Load clinical annotations, dosing guidelines
6. Link PharmGKB variants to dbSNP RS#

### Week 7-8: Consumer Layer
7. **SNPedia** - Crawl pages, extract magnitude/repute
8. Build unified API combining all sources

---

## Sources

- [dbSNP VCF Documentation](https://www.ncbi.nlm.nih.gov/snp/docs/products/vcf/)
- [ClinVar FTP Primer](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/)
- [ClinVar XML Updates 2025](https://academic.oup.com/nar/article/53/D1/D1313/7907366)
- [gnomAD v4.0 Release](https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/)
- [gnomAD VCF Utilities](https://broadinstitute.github.io/gnomad_methods/api_reference/utils/vcf.html)
- [PharmGKB Integrated Resource](https://pmc.ncbi.nlm.nih.gov/articles/PMC8650697/)
- [ClinPGx API](https://api.clinpgx.org/)
- [SNPedia Bulk Download](https://snpedia.com/index.php/Bulk)
- [SNPedia BioStars Discussion](https://www.biostars.org/p/369924/)
- [NCBI Variation Services API](https://api.ncbi.nlm.nih.gov/variation/v0/)
- [gnomAD Azure Dataset](https://learn.microsoft.com/en-us/azure/open-datasets/dataset-gnomad)

---

## Download

| Resource | Method | URL |
|----------|--------|-----|
| dbSNP | FTP/JSON | https://ftp.ncbi.nih.gov/snp/ |
| ClinVar | FTP | https://ftp.ncbi.nih.gov/pub/clinvar/ |
| gnomAD | Cloud/VCF | https://gnomad.broadinstitute.org/downloads |
| PharmGKB | REST API | https://api.clinpgx.org/ |
| SNPedia | MediaWiki API | http://bots.snpedia.com/api.php |

**Access Requirements:** All resources publicly available; no authentication required.

## Data Format

| Format | Description |
|--------|-------------|
| JSON | dbSNP and PharmGKB APIs |
| XML | ClinVar records |
| VCF | Standard genomics format |
| MediaWiki | SNPedia wiki pages |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `rsid` | string | dbSNP reference ID | "rs1801133" |
| `chromosome` | string | Chromosomal location | "1" |
| `position` | integer | 0-based genomic position | 11796321 |
| `allele_frequency` | float | Population allele frequency | 0.35 |
| `clinical_significance` | string | Disease association | "pathogenic" |

## Sample Data

### Example Unified Variant Record
```json
{
  "identifiers": {
    "rsid": "rs1801133",
    "clinvar_id": 17387,
    "gnomad_id": "1-11796321-C-T"
  },
  "location": {
    "chromosome": "1",
    "position": 11796321,
    "reference": "C",
    "alternate": "T",
    "hgvs_genomic": "NC_000001.11:g.11796321C>T"
  },
  "frequency": {
    "gnomad_af": 0.358,
    "populations": {
      "nfe": 0.32,
      "asj": 0.42,
      "eas": 0.28
    }
  },
  "clinical": {
    "clinvar_significance": "benign",
    "conditions": ["Homocystinuria"]
  },
  "functional": {
    "gene": "MTHFR",
    "consequence": "missense_variant",
    "cadd_score": 25.3
  }
}
```

## License

| Resource | License | Notes |
|----------|---------|-------|
| dbSNP | Public Domain | Freely usable |
| ClinVar | Public Domain | Freely usable |
| gnomAD | Open Access | CC0 |
| PharmGKB | CC BY-SA 4.0 | Attribution required |
| SNPedia | CC BY-NC-SA 3.0 | Non-commercial |

## Data Set Size

| Metric | Value |
|--------|-------|
| dbSNP variants | ~1 Billion+ |
| ClinVar records | ~2.5 Million |
| gnomAD variants | ~820 Million |
| PharmGKB drugs | ~715 |
| SNPedia pages | ~50,000+ |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Schema | Formal structure defining data organization and relationships | JSON Schema, OpenAPI |
| SNP | Single Nucleotide Polymorphism - single base pair variation in DNA | rs1800497 |
| Clinical Significance | Classification of a variant's health impact | Pathogenic, Benign |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| dbSNP | NCBI's database of genetic variants | 1B+ SNPs |
| ClinVar | Database of clinically relevant variants | Clinical significance |
| gnomAD | Genome Aggregation Database | Population frequencies |
| PharmGKB | Pharmacogenomics knowledge base | Drug-gene relationships |
| SNPedia | Wiki for SNP annotations | Consumer genetics |
| RefSNP | Reference SNP identifier in dbSNP | rs numbers |
| SPDI | Sequence-Position-Deletion-Insertion notation | Variant representation |
| VCF | Variant Call Format | Variant data standard |
| OpenAPI | Specification for describing REST APIs | Schema definition |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | Data access |
| FTP | File Transfer Protocol | Bulk downloads |
| JSON | JavaScript Object Notation | Data format |
| NCBI | National Center for Biotechnology Information | US agency |
| PMID | PubMed Identifier | Citation reference |
| SNP | Single Nucleotide Polymorphism | Genetic variant |
| SPDI | Sequence Position Deletion Insertion | Variant notation |
| VCF | Variant Call Format | Genomics standard |
| XML | Extensible Markup Language | Data format |
| YAML | YAML Ain't Markup Language | Config format |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 21, 2026 | Research Agent | Initial research report |
