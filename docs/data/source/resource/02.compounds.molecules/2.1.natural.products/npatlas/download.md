---
id: download-npatlas
title: "NPAtlas Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# NPAtlas Download Instructions

## Quick Start

```bash
# Download latest NPAtlas data
wget https://www.npatlas.org/static/downloads/npatlas_v2.json.gz
gunzip npatlas_v2.json.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for decompression
- ~200 MB disk space
- JSON parser (jq recommended) or SDF viewer

## No Registration Required

Data is CC BY 4.0 licensed and freely downloadable.

## Download Methods

### Method 1: JSON Export (Recommended)

```bash
# Complete database in JSON format
wget https://www.npatlas.org/static/downloads/npatlas_v2.json.gz
gunzip npatlas_v2.json.gz

# Preview structure
head -100 npatlas_v2.json | jq '.[0]'
```

### Method 2: SDF Chemical Structures

```bash
# Download SDF with all structures
wget https://www.npatlas.org/static/downloads/npatlas_v2.sdf.gz
gunzip npatlas_v2.sdf.gz

# Count compounds
grep -c '^\$\$\$\$' npatlas_v2.sdf
```

### Method 3: REST API

```bash
# Get single compound by ID
curl "https://www.npatlas.org/api/v1/compound/NPA012345" | jq

# Search compounds
curl -X POST "https://www.npatlas.org/api/v1/search" \
  -H "Content-Type: application/json" \
  -d '{"query": "streptomycin"}'

# Get compounds by organism
curl "https://www.npatlas.org/api/v1/organism/1911/compounds" | jq
```

### Method 4: Web Interface Export

```bash
# 1. Navigate to https://www.npatlas.org
# 2. Use search to find compounds of interest
# 3. Export results in JSON or SDF format
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| npatlas_v2.json.gz | ~50 MB | Complete database (JSON) |
| npatlas_v2.sdf.gz | ~80 MB | Chemical structures (SDF) |
| npatlas_bacteria.json | ~30 MB | Bacterial compounds only |
| npatlas_fungi.json | ~25 MB | Fungal compounds only |

## Data Content

| Category | Records | Description |
|----------|---------|-------------|
| Bacterial | 15,000+ | Actinomycetes, cyanobacteria |
| Fungal | 18,000+ | Ascomycetes, basidiomycetes |
| Marine | 5,000+ | Marine-derived compounds |
| Total references | 25,000+ | Literature citations |

## Post-Download Processing

```bash
# Parse JSON with jq
jq 'length' npatlas_v2.json
# Returns total compound count

# Extract compound names and SMILES
jq -r '.[] | [.npa_id, .name, .smiles] | @tsv' npatlas_v2.json > compounds.tsv

# Filter bacterial compounds
jq '[.[] | select(.origin_type == "bacterial")]' npatlas_v2.json > bacterial.json

# Extract compounds by genus
jq '[.[] | select(.organisms[].genus == "Streptomyces")]' npatlas_v2.json > streptomyces.json

# Convert SDF to other formats (requires Open Babel)
obabel npatlas_v2.sdf -O npatlas_v2.smi
obabel npatlas_v2.sdf -O npatlas_v2.inchi

# Load into database
sqlite3 npatlas.db << 'EOF'
CREATE TABLE compounds (
  npa_id TEXT PRIMARY KEY,
  name TEXT,
  smiles TEXT,
  inchi_key TEXT,
  molecular_weight REAL,
  origin_type TEXT
);
EOF

jq -r '.[] | [.npa_id, .name, .smiles, .inchi_key, .molecular_weight, .origin_type] | @csv' \
  npatlas_v2.json | sqlite3 npatlas.db '.mode csv' '.import /dev/stdin compounds'
```

## Verification

```bash
# Check JSON integrity
jq 'length' npatlas_v2.json

# Verify compound count in SDF
grep -c '^\$\$\$\$' npatlas_v2.sdf

# Validate sample entries
jq '.[0:3] | .[].npa_id' npatlas_v2.json

# Check organism coverage
jq '[.[].organisms[].genus] | unique | length' npatlas_v2.json
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| NPAtlas 3.0 | 2025-01-06 | ~130 MB | Current |
| NPAtlas 2.x | 2022-2024 | ~100 MB | Archived |

### Version Notes

NPAtlas 3.0 contains 36,545 compounds from microbially-derived natural products:
- 1,347 newly curated papers
- 590 structural corrections and revisions
- Covers bacterial, fungal, and cyanobacterial compounds

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://www.npatlas.org/api/v1` |
| Rate Limit | 60-120 req/min |
| Auth Required | No |
| Documentation | https://www.npatlas.org/api/v1/docs |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major updates | Semi-annual |
| Data additions | Monthly |
| Last update | 2025-01 |

## API Rate Limits

| Endpoint | Limit |
|----------|-------|
| Search | 60 requests/minute |
| Compound lookup | 120 requests/minute |
| Bulk operations | Contact maintainers |

## Common Issues

- **Large JSON files**: Use streaming parsers for memory efficiency
- **SDF variants**: Some records have extended metadata fields
- **API pagination**: Large result sets require pagination handling
- **Structure validation**: Verify SMILES validity before use

## JSON Structure

```json
{
  "npa_id": "NPA012345",
  "name": "Compound Name",
  "smiles": "CC(=O)...",
  "inchi_key": "XXXXX-XXXXX-X",
  "molecular_formula": "C21H30O2",
  "molecular_weight": 314.46,
  "origin_type": "bacterial",
  "compound_class": "Polyketide",
  "organisms": [
    {
      "name": "Streptomyces sp.",
      "ncbi_taxon_id": 1234,
      "phylum": "Actinobacteria"
    }
  ],
  "references": [
    {
      "doi": "10.1038/xxx",
      "pmid": 12345678
    }
  ]
}
```

## License

CC BY 4.0 - Free for any use with attribution
