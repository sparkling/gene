---
id: download-topmed
title: "TOPMed Download Instructions"
type: download
parent: _index.md
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

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Freeze versions | Every 1-2 years |
| Bravo updates | With each freeze |
| Current | Freeze 8 (2021) |

## Integration Notes

- Bravo provides lookup only, no bulk download
- For bulk frequencies, consider gnomAD (includes TOPMed)
- Imputation server is free for academic use
- Individual data requires IRB approval and dbGaP authorization
