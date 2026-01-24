---
id: download-npass
title: "NPASS Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# NPASS Download Instructions

## Quick Start

```bash
# Download from NPASS website
# Navigate to http://bidd.group/NPASS/download.html
# Select desired data files and download
```

## Prerequisites

- Web browser for download page access
- **wget** or **curl** (if direct links available)
- ~500 MB disk space for full dataset
- SDF viewer (e.g., Open Babel) for structure files

## Registration

Academic use is free. Contact the BIDD group for commercial access.

## Download Methods

### Method 1: Web Interface (Recommended)

```bash
# 1. Navigate to NPASS download page
# http://bidd.group/NPASS/download.html

# 2. Available downloads include:
#    - Compound data (structures, properties)
#    - Activity data (IC50, Ki values)
#    - Species data (organism sources)
#    - Full database export

# 3. Click desired files to download
```

### Method 2: Individual Data Files

```bash
# Compound structures (SDF format)
# Contains SMILES, InChI, molecular properties

# Activity data (TSV format)
# Contains compound-target-activity triplets

# Species associations (TSV format)
# Contains compound-species links
```

### Method 3: API Access (Limited)

```bash
# Contact BIDD group for API access
# http://bidd.group/

# Basic compound lookup (if available)
curl "http://bidd.group/NPASS/api/compound/NPC12345"
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| npass_compounds.sdf | ~200 MB | Chemical structures |
| npass_activities.tsv | ~150 MB | Bioactivity measurements |
| npass_species.tsv | ~50 MB | Species source data |
| npass_targets.tsv | ~10 MB | Target information |
| npass_references.tsv | ~20 MB | Literature citations |

## Data Categories

| Category | Records | Description |
|----------|---------|-------------|
| Compounds | 35,000+ | Natural product structures |
| Activities | 460,000+ | Quantitative bioactivities |
| Species | 25,000+ | Source organisms |
| Targets | 5,000+ | Protein targets |

## Post-Download Processing

```bash
# Convert SDF to SMILES (requires Open Babel)
obabel npass_compounds.sdf -O npass_compounds.smi

# Preview TSV data
head -20 npass_activities.tsv

# Count activity records
wc -l npass_activities.tsv

# Filter for specific activity type
awk -F'\t' '$3=="IC50"' npass_activities.tsv > ic50_data.tsv

# Load into database
sqlite3 npass.db << 'EOF'
.mode tabs
.import npass_activities.tsv activities
CREATE INDEX idx_compound ON activities(npc_id);
CREATE INDEX idx_target ON activities(target_id);
SELECT COUNT(*) FROM activities;
EOF

# Extract compounds with IC50 < 100 nM
awk -F'\t' '$3=="IC50" && $4<100' npass_activities.tsv > potent_compounds.tsv
```

## Verification

```bash
# Check SDF file structure
grep -c '^\$\$\$\$' npass_compounds.sdf
# Should match compound count

# Validate TSV headers
head -1 npass_activities.tsv

# Count unique compounds with activity data
cut -f1 npass_activities.tsv | sort | uniq | wc -l

# Verify species coverage
cut -f1 npass_species.tsv | sort | uniq | wc -l
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major updates | Annual |
| Data additions | Periodic |
| Last major update | 2023 |

## Common Issues

- **Large files**: Use download manager for large SDF files
- **TSV parsing**: Tab-separated; handle quoted fields properly
- **SDF variants**: Some records may have extended properties
- **Activity units**: Standardize units before analysis (nM vs uM)
- **Missing values**: Some activity records lack complete metadata

## Activity Data Structure

```
npc_id    target_id    activity_type    value    unit    reference
NPC12345  1234         IC50             4.8      nM      PMID:12345678
NPC12345  1235         Ki               12.5     nM      PMID:12345679
```

## Integration Notes

NPASS cross-references with:
- **PubChem**: Via InChI Key matching
- **ChEMBL**: Via compound structure
- **UniProt**: Via target identifiers
- **NCBI Taxonomy**: Via species taxon IDs

## License

Free for academic use. Contact for commercial licensing.
