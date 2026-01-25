---
id: download-civic
title: "CIViC Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# CIViC Download Instructions

## Quick Start

```bash
# Download nightly release
wget https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv
```

## Prerequisites

- wget or curl
- ~50 MB disk space
- GraphQL client (for API)
- jq for JSON processing

## Download Methods

### Primary: Nightly TSV Downloads

```bash
# Clinical Evidence Summaries (most commonly used)
wget https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv

# Variant Summaries
wget https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv

# Gene Summaries
wget https://civicdb.org/downloads/nightly/nightly-GeneSummaries.tsv

# Assertion Summaries
wget https://civicdb.org/downloads/nightly/nightly-AssertionSummaries.tsv

# All files (full data dump)
wget https://civicdb.org/downloads/nightly/nightly-civic_all.tsv
```

### Alternative: GraphQL API

```bash
# Query via GraphQL
curl -X POST https://civicdb.org/api/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ genes(first: 10) { nodes { id name } } }"}' > genes.json

# Get specific variant
curl -X POST https://civicdb.org/api/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ variant(id: 12) { id name gene { name } evidenceItems { nodes { id evidenceType } } } }"}' > braf_v600e.json

# Get evidence for gene
curl -X POST https://civicdb.org/api/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ gene(id: 5) { name variants { nodes { name evidenceItems { nodes { id clinicalSignificance evidenceType } } } } } }"}' > braf_evidence.json
```

### Alternative: JSON Downloads

```bash
# JSON format nightly dumps
wget https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.json
wget https://civicdb.org/downloads/nightly/nightly-VariantSummaries.json
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| nightly-ClinicalEvidenceSummaries.tsv | ~10 MB | All evidence items |
| nightly-VariantSummaries.tsv | ~2 MB | Variant information |
| nightly-GeneSummaries.tsv | ~500 KB | Gene information |
| nightly-AssertionSummaries.tsv | ~200 KB | Expert assertions |
| nightly-civic_all.tsv | ~20 MB | Complete data dump |

## Post-Download Processing

```bash
# View column headers
head -1 nightly-ClinicalEvidenceSummaries.tsv | tr '\t' '\n' | nl

# Filter for specific gene
grep "BRAF" nightly-ClinicalEvidenceSummaries.tsv > braf_evidence.tsv

# Filter for Predictive evidence only
awk -F'\t' '$10 == "Predictive"' nightly-ClinicalEvidenceSummaries.tsv > predictive_only.tsv

# Extract Level A evidence
awk -F'\t' '$11 == "A"' nightly-ClinicalEvidenceSummaries.tsv > level_a_evidence.tsv

# Get unique drugs
cut -f15 nightly-ClinicalEvidenceSummaries.tsv | sort -u > drugs_list.txt
```

## Verification

```bash
# Count evidence items
wc -l nightly-ClinicalEvidenceSummaries.tsv
# Expected: ~8,000+

# Count by evidence type
cut -f10 nightly-ClinicalEvidenceSummaries.tsv | sort | uniq -c

# Count by evidence level
cut -f11 nightly-ClinicalEvidenceSummaries.tsv | sort | uniq -c
```

## Dataset Versions

### Current Release: Nightly (2026-01-25)

| Property | Value |
|----------|-------|
| Version | Nightly build |
| Release Date | Daily |
| Total Size | ~35 MB |
| Evidence Items | 8,000+ |
| Variants | 3,500+ |
| Genes | 500+ |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| nightly-ClinicalEvidenceSummaries.tsv | ~10 MB | 8,000+ | All evidence items |
| nightly-VariantSummaries.tsv | ~2 MB | 3,500+ | Variant information |
| nightly-GeneSummaries.tsv | ~500 KB | 500+ | Gene information |
| nightly-AssertionSummaries.tsv | ~200 KB | 200+ | Expert assertions |
| nightly-civic_all.tsv | ~20 MB | varies | Complete data dump |

### Evidence Type Distribution

| Type | Count | Description |
|------|-------|-------------|
| Predictive | 4,000+ | Drug response |
| Diagnostic | 1,000+ | Disease diagnosis |
| Prognostic | 1,500+ | Outcome prediction |
| Predisposing | 500+ | Cancer susceptibility |
| Oncogenic | 800+ | Functional evidence |

### Evidence Levels

| Level | Count | Description |
|-------|-------|-------------|
| A | 500+ | Validated association |
| B | 2,000+ | Clinical evidence |
| C | 1,500+ | Case study |
| D | 2,000+ | Preclinical data |
| E | 2,000+ | Inferential association |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://civicdb.org/api/graphql |
| Rate Limit | 10 req/sec |
| Auth Required | No |
| Response Format | GraphQL (JSON) |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Nightly dumps | Daily |
| Community curation | Continuous |
| Expert assertions | As reviewed |

## Integration Notes

- CC0 Public Domain license - no restrictions
- GraphQL preferred over deprecated REST API
- VICC harmonization with OncoKB, CGI, etc.
- Use DOID for disease mapping
- Coordinates are GRCh37 (check for GRCh38 liftover needs)
- Can contribute interpretations via web interface
