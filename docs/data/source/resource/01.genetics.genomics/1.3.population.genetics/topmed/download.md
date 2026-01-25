---
id: download-topmed
title: "TOPMed Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# TOPMed Download Instructions

## Quick Start

```bash
# Query Bravo API for variant frequency (public access)
curl "https://bravo.sph.umich.edu/freeze8/hg38/api/v1/variant?variant_id=1-10177-A-AC"
```

## Prerequisites

- curl or browser (Bravo public access)
- dbGaP authorization (individual-level data)
- ~50 TB disk space (full individual data)

## Download Methods

### Primary: Bravo Browser (Public Frequencies)

```bash
# Query single variant
curl "https://bravo.sph.umich.edu/freeze8/hg38/api/v1/variant?variant_id=chr-pos-ref-alt"

# Example: Query BRAF V600E region
curl "https://bravo.sph.umich.edu/freeze8/hg38/api/v1/variant?variant_id=7-140753336-A-T"

# Web interface for browsing
# https://bravo.sph.umich.edu/
```

### Alternative: TOPMed Imputation Server

```bash
# Submit VCF for imputation using TOPMed reference panel
# https://imputation.biodatacatalyst.nhlbi.nih.gov/

# No download required - cloud-based service
# Outputs imputed VCF with TOPMed-based frequencies
```

### Individual-Level Data: dbGaP (Authorized Access)

```bash
# 1. Apply for access at dbGaP
# https://www.ncbi.nlm.nih.gov/gap/

# 2. Once approved, use SRA toolkit
prefetch SRR_accession
fastq-dump --split-files SRR_accession

# 3. Or access via BioData Catalyst
# https://biodatacatalyst.nhlbi.nih.gov/
```

## File Inventory

| File | Size | Access | Description |
|------|------|--------|-------------|
| Bravo API | N/A | Public | Variant frequencies |
| TOPMed VCFs | ~50 TB | dbGaP | Individual genotypes |
| Imputation panel | ~100 GB | Server | Reference panel |
| Phenotype files | Varies | dbGaP | Study-specific |

## Post-Download Processing

```bash
# Parse Bravo API response (public)
curl -s "https://bravo.sph.umich.edu/freeze8/hg38/api/v1/variant?variant_id=1-10177-A-AC" | jq '.'

# For dbGaP data, follow study-specific protocols
# Typically involves:
# 1. Decryption using dbGaP tools
# 2. QC filtering
# 3. Ancestry assignment
```

## Verification

```bash
# Test Bravo API
curl -s "https://bravo.sph.umich.edu/freeze8/hg38/api/v1/variant?variant_id=1-10177-A-AC" | jq '.allele_freq'
# Should return frequency value
```

## Dataset Versions

### Current Release: Freeze 10b

| Property | Value |
|----------|-------|
| Version | Freeze 10b (GRCh38) |
| Release Date | 2024-Q4 |
| Total Size | ~50 TB (full WGS) |
| Samples | ~206,000 individuals |
| SNVs | 781 million |
| Indels | 62 million |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| Bravo API | N/A | 781M | Variant frequencies (public) |
| TOPMed VCFs | ~50 TB | 843M | Individual genotypes (dbGaP) |
| Imputation panel | ~100 GB | - | Reference panel (server) |

### Freeze History

| Version | Release | Samples | Status |
|---------|---------|---------|--------|
| Freeze 10b | 2024 | 206K | Current |
| Freeze 9b | 2023 | 180K | Available |
| Freeze 8 | 2021 | 150K | Available |
| Freeze 5b | 2018 | 62K | Archived |

### Ancestry Breakdown (Bravo)

| Ancestry | Proportion |
|----------|------------|
| European | 52% |
| African American | 25% |
| Hispanic/Latino | 12% |
| Asian | 8% |
| Other | 3% |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://bravo.sph.umich.edu/freeze8/hg38/api/v1/ |
| Rate Limit | 10 req/sec |
| Auth Required | No (public frequencies) |
| Response Format | JSON |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Freeze versions | Every 1-2 years |
| Bravo updates | With each freeze |
| Current | Freeze 10b (2024) |

## Integration Notes

- Bravo provides lookup only, no bulk download
- For bulk frequencies, consider gnomAD (includes TOPMed)
- Imputation server is free for academic use
- Individual data requires IRB approval and dbGaP authorization
