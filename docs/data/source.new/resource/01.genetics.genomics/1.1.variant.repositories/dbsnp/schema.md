---
id: schema-dbsnp
title: "dbSNP Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./_index.md)

# dbSNP Schema Documentation

**Document ID:** SCHEMA-DBSNP
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** NCBI Variation Services API v0.1.10

---

## TL;DR

dbSNP (Database of Single Nucleotide Polymorphisms) is NCBI's public archive for genetic variation data. This document details the actual API schemas and data structures retrieved from the NCBI Variation Services API.

**Base URL:** `https://api.ncbi.nlm.nih.gov/variation/v0`
**Rate Limit:** 1 request/second recommended

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| **ALFA R4 Variants** | 900M+ variants | May 2025 release |
| **ALFA Subjects** | 400K+ subjects | May 2025 release |
| **ALFA Populations** | 12 populations | R4 specification |
| **Target Cohort** | 1M+ dbGaP subjects | Ongoing expansion |
| **Build** | GRCh38.p14 | Current reference |

---

## API Endpoint Categories

### 1. SPDI Endpoints (Sequence-Position-Deletion-Insertion)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/spdi/{spdi}/` | GET | Retrieve validated SPDI with associated resources |
| `/spdi/{spdi}/contextual` | GET | Get contextual allele (precision-corrected) |
| `/spdi/{spdi}/vcf_fields` | GET | Convert to VCF format (CHROM, POS, REF, ALT) |
| `/spdi/{spdi}/canonical_representative` | GET | Retrieve canonical allele |
| `/spdi/{spdi}/all_equivalent_contextual` | GET | Find equivalent alleles across sequences |
| `/spdi/{spdi}/rsids` | GET | Lookup associated RefSNP IDs |
| `/spdi/{spdi}/hgvs` | GET | Convert to HGVS notation |
| `/spdi/{spdi}/ga4gh-vr` | GET | Convert to GA4GH VR format (normalized) |
| `/spdi/{spdi}/ga4gh-vr-unnormalized` | GET | Convert to GA4GH VR (unnormalized) |
| `/spdi/` | POST | Convert GA4GH VR to SPDI |

### 2. HGVS Endpoints

| Endpoint | Method | Limit | Description |
|----------|--------|-------|-------------|
| `/hgvs/{hgvs}/contextuals` | GET | N/A | Retrieve contextual SPDI from HGVS |
| `/hgvs/batch/contextuals` | POST | 50,000 | Batch conversion |

### 3. VCF Endpoints

| Endpoint | Method | Limit | Description |
|----------|--------|-------|-------------|
| `/vcf/{chrom}/{pos}/{ref}/{alts}/contextuals` | GET | N/A | Convert VCF to contextual SPDI |
| `/vcf/file/set_rsids` | POST | 50,000 rows | Annotate VCF with RefSNP IDs |

### 4. RefSNP Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/refsnp/{rsid}` | GET | Retrieve RefSNP snapshot object |

### 5. ALFA Frequency Endpoints

| Endpoint | Method | Limit | Description |
|----------|--------|-------|-------------|
| `/metadata/frequency` | GET | N/A | Retrieve frequency study metadata |
| `/interval/{seq_id}:{position}:{length}/overlapping_frequency_records` | GET | 250 | Query variants by interval |
| `/refsnp/{rsid}/frequency` | GET | N/A | Lookup ALFA frequencies for RSID |

---

## Core Data Schemas

### SPDI Object

The canonical variant representation format.

```json
{
  "seq_id": "NC_000011.10",
  "position": 5227001,
  "deleted_sequence": "T",
  "inserted_sequence": "A",
  "warnings": {}
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `seq_id` | string | Yes | Accession.version (RefSeq identifier) |
| `position` | integer | Yes | 0-based position on sequence |
| `deleted_sequence` | string | Yes | IUPAC nucleotide sequence deleted |
| `inserted_sequence` | string | Yes | IUPAC nucleotide sequence inserted |
| `warnings` | object | No | Optional validation warnings |

### VCF Fields Object

Standard VCF representation.

```json
{
  "chrom": "11",
  "pos": 5227002,
  "ref": "T",
  "alt": "A",
  "warnings": {}
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `chrom` | string | Yes | Chromosome identifier |
| `pos` | integer | Yes | 1-based position (VCF standard) |
| `ref` | string | Yes | Reference allele |
| `alt` | string | Yes | Alternate allele |
| `warnings` | object | No | Optional validation warnings |

### HGVS Object

Human Genome Variation Society nomenclature.

```json
{
  "hgvs": "NC_000011.10:g.5227002T>A",
  "warnings": {}
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `hgvs` | string | Yes | Complete HGVS expression |
| `warnings` | object | No | Optional validation warnings |

### GA4GH VR Allele Object

Global Alliance for Genomics and Health Variant Representation.

```json
{
  "type": "Allele",
  "location": {
    "type": "SequenceLocation",
    "sequence_id": "refseq:NC_000011.10",
    "interval": {
      "type": "SequenceInterval",
      "start": {"type": "Number", "value": 5227001},
      "end": {"type": "Number", "value": 5227002}
    }
  },
  "state": {
    "type": "LiteralSequenceExpression",
    "sequence": "A"
  },
  "_id": "ga4gh:VA.xxx"
}
```

---

## RefSNP Snapshot Schema (Full Record)

Sample response from `/refsnp/{rsid}` (rs334 - sickle cell variant):

### Top-Level Structure

```json
{
  "refsnp_id": "334",
  "create_date": "2000-09-19T17:02:00Z",
  "last_update_date": "2025-01-15T10:30:00Z",
  "last_update_build_id": "156",
  "dbsnp1_merges": [],
  "citations": [12345, 67890],
  "lost_obs_movements": [],
  "present_obs_movements": [],
  "primary_snapshot_data": {},
  "support": [],
  "anchor": "NC_000011.10",
  "variant_type": "snv",
  "ga4gh": {},
  "mane_select_ids": ["NM_000518.5"]
}
```

### Key Fields

| Field | Type | Description |
|-------|------|-------------|
| `refsnp_id` | string | Unique RS identifier (without "rs" prefix) |
| `create_date` | ISO8601 | Initial submission date |
| `last_update_date` | ISO8601 | Most recent modification |
| `last_update_build_id` | string | dbSNP build number |
| `dbsnp1_merges` | array | Historical RS merges |
| `citations` | array[int] | PubMed IDs |
| `present_obs_movements` | array | Current observation mappings |
| `primary_snapshot_data` | object | Main variant data |
| `support` | array | Submission records |
| `anchor` | string | Primary reference sequence |
| `variant_type` | enum | snv, mnv, ins, del, delins, identity |
| `mane_select_ids` | array | MANE Select transcript IDs |

### Primary Snapshot Data Structure

```json
{
  "placements_with_allele": [
    {
      "seq_id": "NC_000011.10",
      "is_ptlp": true,
      "placement_annot": {
        "seq_type": "refseq_chromosome",
        "seq_id_traits_by_assembly": []
      },
      "alleles": [
        {
          "allele": {
            "spdi": {
              "seq_id": "NC_000011.10",
              "position": 5227001,
              "deleted_sequence": "T",
              "inserted_sequence": "T"
            }
          },
          "hgvs": "NC_000011.10:g.5227002="
        }
      ]
    }
  ],
  "allele_annotations": []
}
```

### Allele Annotations Structure

Includes population frequencies and clinical data:

```json
{
  "frequency": [
    {
      "study_name": "1000Genomes",
      "study_version": "Phase3",
      "observation": {
        "seq_id": "NC_000011.10",
        "position": 5227001,
        "deleted_sequence": "T",
        "inserted_sequence": "A"
      },
      "allele_count": 54,
      "total_count": 5008
    }
  ],
  "clinical": [
    {
      "accession_version": "RCV000018137.6",
      "disease_names": ["Sickle cell anemia"],
      "clinical_significances": ["Pathogenic"],
      "review_status": "criteria provided, single submitter"
    }
  ],
  "sequence_ontology": [
    {
      "accession": "SO:0001583",
      "name": "missense_variant"
    }
  ]
}
```

---

## ALFA Frequency Data Schema

### Frequency Response Structure

```json
{
  "build_id": "20250407153717",
  "project": "PRJNA507278",
  "reference_allele": "C",
  "results": [
    {
      "sample_id": "SAMN10492705",
      "allele_counts": {
        "T": 22382,
        "C": 258364
      }
    }
  ]
}
```

### ALFA Population Codes

| Sample ID | Population | Description |
|-----------|------------|-------------|
| SAMN10492695 | AFR | African |
| SAMN10492696 | AMR | American (Latin) |
| SAMN10492697 | ASN | Asian |
| SAMN10492698 | EAS | East Asian |
| SAMN10492699 | EUR | European |
| SAMN10492700 | SAS | South Asian |
| SAMN10492701 | OTH | Other |
| SAMN10492702 | Total | All populations |
| SAMN10492703 | AFO | African Other |
| SAMN10492704 | AFA | African American |
| SAMN10492705 | LAC | Latin American |
| SAMN11605645 | LEN | Latino/Hispanic |

---

## VCF File Format (dbSNP 2.0)

### Fixed Columns

| Column | Position | Description |
|--------|----------|-------------|
| CHROM | 1 | RefSeq chromosome identifier |
| POS | 2 | 1-based position (sorted) |
| ID | 3 | RS identifier (with suffix for multi-mapping) |
| REF | 4 | Reference base(s): A, C, G, T, N |
| ALT | 5 | Comma-separated alternate alleles |
| QUAL | 6 | Always '.' (not used) |
| FILTER | 7 | Always '.' (not used) |
| INFO | 8 | Semicolon-separated key-value pairs |

### INFO Field Tags

#### Clinical Attributes
| Tag | Description |
|-----|-------------|
| CLNACC | ClinVar accession(s) |
| CLNDN | ClinVar disease name(s) |
| CLNSIG | Clinical significance |
| CLNREVSTAT | ClinVar review status |

#### Sequence Ontology Location
| Tag | SO ID | Description |
|-----|-------|-------------|
| ASS | SO:0001574 | Splice acceptor variant |
| DSS | SO:0001575 | Splice donor variant |
| INT | SO:0001627 | Intron variant |
| R3 | SO:0001631 | 3' UTR variant |
| R5 | SO:0001623 | 5' UTR variant |
| U3 | SO:0001632 | Downstream gene variant |
| U5 | SO:0001636 | Upstream gene variant |

#### Functional Tags
| Tag | SO ID | Description |
|-----|-------|-------------|
| NSF | SO:0001587 | Frameshift variant |
| NSM | SO:0001583 | Missense variant |
| NSN | SO:0001587 | Nonsense variant |
| SYN | SO:0001819 | Synonymous variant |

#### Frequency Data
| Tag | Description |
|-----|-------------|
| CAF | Comma-delimited allele frequencies (1000Genomes) |

### Variant Types
| Type | SO ID | Description |
|------|-------|-------------|
| SNV | SO:0001483 | Single nucleotide variant |
| INS | SO:0000667 | Insertion |
| DEL | SO:0000159 | Deletion |
| INDEL | SO:1000032 | Insertion-deletion |
| MNV | SO:0002007 | Multiple nucleotide variant |
| NOALT | N/A | No alternate allele |

---

## Response Formats

### Success Responses
- **200 OK**: JSON payload with `data` field
- **206 Partial Content**: Frequency endpoint (limited to 250 results)

### Error Responses

```json
{
  "error": {
    "code": 400,
    "message": "Invalid SPDI format",
    "errors": [
      {
        "reason": "invalidParameter",
        "message": "Sequence ID not found"
      }
    ]
  }
}
```

| Code | Description |
|------|-------------|
| 400 | Bad Request (syntax errors) |
| 404 | Not Found (unsupported sequences/RSIDs) |
| 405 | Method Not Allowed |

---

## Content Negotiation

| Endpoint Type | Content-Type |
|---------------|--------------|
| SPDI resources | `application/hal+json` |
| GA4GH VR input | `application/vnd.ga4gh.vr+json` |
| VCF file operations | `text/plain; charset=utf-8` |
| Default responses | `application/json` |

---

## Batch Request Limits

| Endpoint | Maximum |
|----------|---------|
| HGVS batch | 50,000 expressions |
| VCF batch | 50,000 rows |
| Frequency interval | 250 results |

---

## Data Files (FTP Downloads)

### ALFA Frequency Files

| File | Size | Description |
|------|------|-------------|
| `freq.vcf.gz` | 16 GB | Complete population frequencies |
| `freq.vcf.gz.tbi` | 2.8 MB | Tabix index |

**FTP Location:** `ftp://ftp.ncbi.nlm.nih.gov/snp/population_frequency/latest_release/`

---

## Integration Examples

### Python: Fetch RefSNP

```python
import requests

def get_refsnp(rsid: str) -> dict:
    """Fetch RefSNP data from NCBI Variation Services API."""
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid}"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

# Example: rs334 (sickle cell)
data = get_refsnp("334")
print(f"Variant type: {data['variant_type']}")
print(f"Last updated: {data['last_update_date']}")
```

### Python: Batch HGVS Conversion

```python
import requests

def batch_hgvs_to_spdi(hgvs_list: list[str]) -> dict:
    """Convert multiple HGVS expressions to SPDI."""
    url = "https://api.ncbi.nlm.nih.gov/variation/v0/hgvs/batch/contextuals"
    response = requests.post(url, json={"hgvs": hgvs_list})
    response.raise_for_status()
    return response.json()

# Example
expressions = [
    "NM_000518.5:c.20A>T",
    "NM_005228.5:c.2573T>G"
]
results = batch_hgvs_to_spdi(expressions)
```

---

## Sample Data

### Example Record
```json
{
  "rs_id": "rs334",
  "refsnp_id": 334,
  "chromosome": "11",
  "position": 5227002,
  "reference_allele": "A",
  "alternate_allele": "T",
  "allele_frequency": 0.0108,
  "clinical_significance": "Pathogenic"
}
```

### Sample Query Result
| rs_id | chromosome | position | reference_allele | alternate_allele | allele_frequency | clinical_significance |
|-------|-----------|----------|------------------|------------------|-----------------|----------------------|
| rs334 | 11 | 5227002 | A | T | 0.0108 | Pathogenic |
| rs1805377 | 1 | 11796321 | G | A | 0.0045 | Likely pathogenic |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "rs12345" |
| `name` | string | Entity name | "dbSNP Entry" |
| `type` | string | Record type | "polymorphism" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| dbSNP | API | https://api.ncbi.nlm.nih.gov/variation/v0 |
| ALFA Frequencies | FTP | ftp://ftp.ncbi.nlm.nih.gov/snp/population_frequency/latest_release/ |
| VCF Files | FTP | ftp://ftp.ncbi.nlm.nih.gov/snp/latest_release/ |

**Access Requirements:** Open access via NCBI, no registration required

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| dbSNP | Public Domain (NCBI) | Yes |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 1,000,000,000+ |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | VCF (gzip compressed) |
| Alternative | JSON (API) |
| Encoding | UTF-8 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `refsnp_id` | Unique RefSNP identifier (without "rs" prefix) | 334 |
| `seq_id` | RefSeq sequence accession with version | NC_000011.10 |
| `position` | 0-based genomic coordinate in SPDI format | 5227001 |
| `deleted_sequence` | Reference nucleotide(s) at variant position | T |
| `inserted_sequence` | Alternate allele nucleotide(s) | A |
| `variant_type` | Classification of variant (snv, ins, del, etc.) | snv |
| `allele_count` | Number of observations of specific allele | 54 |
| `total_count` | Total number of chromosomes sampled | 5008 |
| `CLNSIG` | ClinVar clinical significance classification | Pathogenic |
| `mane_select_ids` | MANE Select transcript identifiers | NM_000518.5 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| SPDI | Sequence-Position-Deletion-Insertion notation | Canonical variant format |
| RefSNP | Reference SNP cluster in dbSNP | Primary variant record |
| Contextual Allele | Precision-corrected allele representation | SPDI endpoint |
| Canonical Allele | Standardized representative allele form | Normalization |
| ALFA | Allele Frequency Aggregator | Population frequency data |
| Clinical Significance | Pathogenicity classification from ClinVar | Clinical annotation |
| Sequence Ontology | Standardized vocabulary for sequence features | SO:XXXXXXX |
| Observation Movement | Tracking of variant positions across builds | Build updates |
| Population Frequency | Allele frequency in specific population | ALFA data |
| MANE Select | Matched Annotation from NCBI and EMBL | Transcript standard |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| dbSNP | Database of Single Nucleotide Polymorphisms | NCBI variant archive |
| SPDI | Sequence-Position-Deletion-Insertion | Variant notation |
| HGVS | Human Genome Variation Society | Nomenclature standard |
| VCF | Variant Call Format | File format |
| ALFA | Allele Frequency Aggregator | Frequency database |
| GA4GH | Global Alliance for Genomics and Health | Standards organization |
| VR | Variant Representation | GA4GH specification |
| SNV | Single Nucleotide Variant | Variant type |
| MNV | Multiple Nucleotide Variant | Variant type |
| INDEL | Insertion-Deletion | Variant type |
| SO | Sequence Ontology | Feature vocabulary |
| RefSeq | Reference Sequence | NCBI sequence database |
| GRCh38 | Genome Reference Consortium Human Build 38 | Reference genome |
| NCBI | National Center for Biotechnology Information | Database host |
| HAL | Hypertext Application Language | API response format |
| MANE | Matched Annotation from NCBI and EMBL | Transcript standard |

---

## References

1. NCBI Variation Services API: https://api.ncbi.nlm.nih.gov/variation/v0/
2. dbSNP Documentation: https://www.ncbi.nlm.nih.gov/snp/docs/
3. ALFA Documentation: https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/
4. Sherry ST, et al. (2001). dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 29(1):308-11.

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation from API spec |
