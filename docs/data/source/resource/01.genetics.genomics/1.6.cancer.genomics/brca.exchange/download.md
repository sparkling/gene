---
id: download-brca-exchange
title: "BRCA Exchange Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# BRCA Exchange Download Instructions

## Quick Start

```bash
# Download full variant list via API
curl "https://brcaexchange.org/backend/downloads/releases/current_release.tar.gz" -o brca_exchange.tar.gz
```

## Prerequisites

- curl or web browser
- ~100 MB disk space
- tar/gzip for extraction

## Download Methods

### Primary: BRCA Exchange Downloads

```bash
# Current release (all variants)
curl "https://brcaexchange.org/backend/downloads/releases/current_release.tar.gz" -o brca_exchange.tar.gz

# Extract
tar -xzf brca_exchange.tar.gz

# Files include:
# - built_with_change_types.tsv (full data)
# - built.tsv (simplified)
# - release_notes.txt
```

### Alternative: API Query

```bash
# Search for specific variant
curl "https://brcaexchange.org/backend/data/?format=json&search=c.68_69delAG"

# Get variant by ID
curl "https://brcaexchange.org/backend/data/?format=json&variant_id=12345"

# Filter by gene
curl "https://brcaexchange.org/backend/data/?format=json&Gene_Symbol=BRCA1&page_size=100"

# Filter by classification
curl "https://brcaexchange.org/backend/data/?format=json&Pathogenicity_expert=Pathogenic&page_size=100"
```

### Alternative: Web Interface

```bash
# Browse and download via web
# https://brcaexchange.org/variants

# Advanced search allows filtering by:
# - Gene (BRCA1/BRCA2)
# - Classification
# - Variant type
# - Source database
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| current_release.tar.gz | ~50 MB | Compressed release |
| built_with_change_types.tsv | ~200 MB | Full variant data |
| built.tsv | ~100 MB | Simplified version |
| release_notes.txt | <1 KB | Release information |

## Post-Download Processing

```bash
# Extract archive
tar -xzf brca_exchange.tar.gz
cd output/

# View columns
head -1 built_with_change_types.tsv | tr '\t' '\n' | nl

# Filter pathogenic variants
awk -F'\t' '$5 == "Pathogenic"' built_with_change_types.tsv > pathogenic_only.tsv

# Extract BRCA1 only
awk -F'\t' '$1 == "BRCA1"' built_with_change_types.tsv > brca1_variants.tsv

# Get ENIGMA-classified variants
awk -F'\t' '$5 != ""' built_with_change_types.tsv > enigma_classified.tsv
```

## Verification

```bash
# Count total variants
wc -l built_with_change_types.tsv

# Count by gene
cut -f1 built_with_change_types.tsv | sort | uniq -c

# Count by classification
cut -f5 built_with_change_types.tsv | sort | uniq -c
```

## Dataset Versions

### Current Release

| Property | Value |
|----------|-------|
| Version | Release 15 |
| Release Date | 2025-12 |
| Total Size | ~200 MB |
| Variants | 45,000+ |
| Sources | 15 databases |
| ENIGMA Classified | 8,000+ |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| current_release.tar.gz | ~50 MB | 45K | Compressed release |
| built_with_change_types.tsv | ~200 MB | 45K | Full variant data |
| built.tsv | ~100 MB | 45K | Simplified version |
| release_notes.txt | <1 KB | - | Release information |

### Variant Classification

| Classification | Count | Description |
|----------------|-------|-------------|
| Pathogenic | 8,000+ | Disease-causing |
| Likely Pathogenic | 2,000+ | Probably pathogenic |
| VUS | 20,000+ | Uncertain significance |
| Likely Benign | 3,000+ | Probably benign |
| Benign | 10,000+ | Not pathogenic |

### Gene Distribution

| Gene | Variants | Description |
|------|----------|-------------|
| BRCA1 | 20,000+ | Chromosome 17 |
| BRCA2 | 25,000+ | Chromosome 13 |

### Data Sources Aggregated

| Source | Contribution | Description |
|--------|--------------|-------------|
| ClinVar | Major | Clinical submissions |
| LOVD | Major | Locus-specific database |
| BIC | Historical | Breast cancer database |
| ENIGMA | Expert | Expert classifications |
| GnomAD | Frequencies | Population data |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://brcaexchange.org/backend/ |
| Rate Limit | 10 req/sec |
| Auth Required | No |
| Response Format | JSON |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Data updates | Monthly |
| ENIGMA reviews | Ongoing |
| Major releases | Quarterly |

## Integration Notes

- ENIGMA classifications are gold standard
- Cross-reference with ClinVar for additional evidence
- Use HGVS nomenclature for variant matching
- API supports pagination for large queries
- Free for academic and clinical use
