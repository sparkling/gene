---
id: download-stitch
title: "STITCH Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# STITCH Download Instructions

## Quick Start

```bash
# Download human chemical-protein interactions (high confidence)
wget http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz
gunzip 9606.protein_chemical.links.v5.0.tsv.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for decompression
- 10-50GB disk space (varies by species/data type)
- For API: No authentication required

## No Registration Required

All STITCH data is freely available under CC BY 4.0 license without registration.

## Download Methods

### Method 1: Species-Specific Downloads (Recommended)

```bash
# Create download directory
mkdir -p stitch && cd stitch

# Download human (9606) data
# Chemical-protein interactions
wget http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz

# Detailed scores (includes transfer evidence)
wget http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0/9606.protein_chemical.links.detailed.v5.0.tsv.gz

# Action annotations (activation, inhibition, etc.)
wget http://stitch.embl.de/download/actions.v5.0/9606.actions.v5.0.tsv.gz

# Chemical aliases for ID mapping
wget http://stitch.embl.de/download/chemical.aliases.v5.0/9606.chemical.aliases.v5.0.tsv.gz
```

### Method 2: Download All Species

```bash
# Full database - all species
# WARNING: Very large files (10-50GB compressed)

# All chemical-protein interactions
wget http://stitch.embl.de/download/protein_chemical.links.v5.0.tsv.gz

# All detailed interactions
wget http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0.tsv.gz

# All chemical-chemical links
wget http://stitch.embl.de/download/chemical_chemical.links.v5.0.tsv.gz

# All actions
wget http://stitch.embl.de/download/actions.v5.0.tsv.gz

# Chemical information
wget http://stitch.embl.de/download/chemicals.v5.0.tsv.gz
wget http://stitch.embl.de/download/chemical.sources.v5.0.tsv.gz
```

### Method 3: Download Multiple Species

```bash
# Common model organisms
SPECIES="9606 10090 10116 7955 7227 6239 4932"

for taxid in $SPECIES; do
  echo "Downloading species $taxid..."
  wget -q "http://stitch.embl.de/download/protein_chemical.links.v5.0/${taxid}.protein_chemical.links.v5.0.tsv.gz"
  wget -q "http://stitch.embl.de/download/actions.v5.0/${taxid}.actions.v5.0.tsv.gz"
done
```

### Method 4: API Bulk Query

```bash
# Query interactions for specific chemicals
# Aspirin interactions
curl "http://stitch.embl.de/api/tsv/interactions?identifier=CIDm00002244&species=9606&limit=100" \
  > aspirin_interactions.tsv

# Query by chemical name
curl "http://stitch.embl.de/api/tsv/interactions?identifier=aspirin&species=9606" \
  > aspirin_by_name.tsv

# Multiple chemicals
curl "http://stitch.embl.de/api/tsv/interactions?identifiers=aspirin%0dibuprofen%0dmetformin&species=9606" \
  > multiple_drugs.tsv

# Get JSON format
curl "http://stitch.embl.de/api/json/interactions?identifier=aspirin&species=9606&limit=50" \
  > aspirin_interactions.json

# Network query (chemicals + proteins)
curl "http://stitch.embl.de/api/tsv/network?identifiers=aspirin%0dTP53%0dPTGS2&species=9606" \
  > mixed_network.tsv
```

### Method 5: FTP Download

```bash
# Connect to STITCH FTP
ftp stitch.embl.de
# Username: anonymous
# Password: your@email.com

cd download
ls

# Download specific files
get protein_chemical.links.v5.0.tsv.gz
get chemicals.v5.0.tsv.gz
quit

# Or via command line
wget -r -np ftp://stitch.embl.de/download/protein_chemical.links.v5.0/
```

## File Inventory

### Interaction Files

| File | Size | Description |
|------|------|-------------|
| protein_chemical.links.v5.0.tsv.gz | ~10 GB | All interactions, basic scores |
| protein_chemical.links.detailed.v5.0.tsv.gz | ~15 GB | All interactions, detailed scores |
| 9606.protein_chemical.links.v5.0.tsv.gz | ~200 MB | Human interactions only |
| 9606.protein_chemical.links.detailed.v5.0.tsv.gz | ~300 MB | Human detailed |

### Chemical Files

| File | Size | Description |
|------|------|-------------|
| chemicals.v5.0.tsv.gz | ~50 MB | Chemical names and properties |
| chemical.aliases.v5.0.tsv.gz | ~500 MB | ID mappings (PubChem, DrugBank, etc.) |
| chemical.sources.v5.0.tsv.gz | ~2 GB | Evidence source annotations |
| chemical_chemical.links.v5.0.tsv.gz | ~5 GB | Chemical similarity network |

### Action Files

| File | Size | Description |
|------|------|-------------|
| actions.v5.0.tsv.gz | ~3 GB | All action annotations |
| 9606.actions.v5.0.tsv.gz | ~50 MB | Human actions only |

### Species-Specific Files

| Taxon | Links File | Detailed File | Actions |
|-------|------------|---------------|---------|
| 9606 (Human) | ~200 MB | ~300 MB | ~50 MB |
| 10090 (Mouse) | ~150 MB | ~250 MB | ~40 MB |
| 10116 (Rat) | ~100 MB | ~180 MB | ~30 MB |
| 7955 (Zebrafish) | ~80 MB | ~140 MB | ~20 MB |

## Post-Download Processing

### Filter by Confidence Score

```bash
# Extract high-confidence interactions (score >= 700)
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  awk -F'\t' 'NR==1 || $NF >= 700' > human_high_confidence.tsv

# Extract highest confidence only (score >= 900)
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  awk -F'\t' 'NR==1 || $NF >= 900' > human_highest_confidence.tsv
```

### Map Chemical IDs to PubChem

```bash
# Extract PubChem CID from STITCH ID
# CIDm00002244 -> 2244
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  awk -F'\t' 'NR>1 {
    cid = $1
    gsub(/^CID[ms]0*/, "", cid)
    print cid, $0
  }' > interactions_with_pubchem.tsv
```

### Map Proteins to Gene Symbols

```bash
# Download STRING aliases for gene symbol mapping
wget http://string-db.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz

# Create mapping file
gunzip -c 9606.protein.aliases.v11.5.txt.gz | \
  awk -F'\t' '$3 ~ /BioMart_HUGO/' | \
  cut -f1,2 > string_to_symbol.map

# Apply mapping to STITCH data
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  awk -F'\t' 'NR==FNR {map[$1]=$2; next}
    {print $0, (map[$2] ? map[$2] : "NA")}' \
    string_to_symbol.map - > interactions_with_symbols.tsv
```

### Extract Specific Interaction Types

```bash
# Get only inhibition actions
gunzip -c 9606.actions.v5.0.tsv.gz | \
  awk -F'\t' '$3 == "inhibition"' > inhibition_actions.tsv

# Get all drug-like compound interactions
# (requires chemical.aliases with DrugBank IDs)
gunzip -c chemical.aliases.v5.0.tsv.gz | \
  grep "DrugBank" | cut -f1 | sort -u > drugbank_chemicals.txt

gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  grep -Ff drugbank_chemicals.txt > drug_interactions.tsv
```

### Build Network File

```bash
# Create Cytoscape-compatible network
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  awk -F'\t' 'NR>1 && $NF >= 700 {print $1, "interacts_with", $2, $NF}' OFS='\t' \
  > stitch_network.sif

# Add edge attributes
gunzip -c 9606.actions.v5.0.tsv.gz | \
  awk -F'\t' '{print $1, $2, $3, $NF}' OFS='\t' \
  > stitch_edge_attributes.txt
```

### Python Processing

```python
import pandas as pd

# Load interactions
interactions = pd.read_csv(
    "9606.protein_chemical.links.v5.0.tsv.gz",
    sep='\t', compression='gzip'
)

# Filter high confidence
high_conf = interactions[interactions['combined_score'] >= 700]
print(f"High confidence interactions: {len(high_conf)}")

# Convert chemical IDs
high_conf['pubchem_cid'] = high_conf['chemical'].str.replace(
    r'^CID[ms]0*', '', regex=True
).astype(int)

# Extract gene from protein ID
high_conf['taxid'] = high_conf['protein'].str.split('.').str[0]
high_conf['ensp'] = high_conf['protein'].str.split('.').str[1]

# Save processed data
high_conf.to_csv('stitch_human_processed.tsv', sep='\t', index=False)
```

## Verification

```bash
# Check file integrity
gunzip -t 9606.protein_chemical.links.v5.0.tsv.gz && echo "File OK"

# Count total interactions
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | wc -l
# Expected: ~16 million for human

# Verify column structure
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | head -1
# Should show: chemical protein experimental prediction database textmining combined_score

# Check score distribution
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  awk -F'\t' 'NR>1 {
    if ($NF >= 900) high++
    else if ($NF >= 700) medium++
    else if ($NF >= 400) low++
    else vlow++
  }
  END {
    print "Highest (900+):", high
    print "High (700-899):", medium
    print "Medium (400-699):", low
    print "Low (<400):", vlow
  }'

# Verify species ID
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  awk -F'\t' 'NR>1 {print $2}' | head -100 | grep -c "^9606\."
# Should be 100 (all human proteins)

# Check unique chemicals count
gunzip -c 9606.protein_chemical.links.v5.0.tsv.gz | \
  awk -F'\t' 'NR>1 {print $1}' | sort -u | wc -l
# Expected: ~300,000-500,000
```

---

## Dataset Versions

### Current Release: STITCH v5.0

| Property | Value |
|----------|-------|
| Version | 5.0 |
| Release Date | 2016-01-01 |
| Total Size | ~100 GB (all species) |
| Chemicals | ~500,000 |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| 9606.protein_chemical.links.v5.0.tsv.gz | ~200 MB | ~16M | Human interactions |
| 9606.actions.v5.0.tsv.gz | ~50 MB | ~2M | Human actions |
| chemicals.v5.0.tsv.gz | ~50 MB | ~500K | Chemical info |
| chemical.aliases.v5.0.tsv.gz | ~500 MB | ~10M | ID mappings |

### Previous Versions

| Version | Release | Status |
|---------|---------|--------|
| 4.0 | 2014-01-01 | Archived |
| 3.1 | 2012-06-01 | Archived |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| Base URL | `http://stitch.embl.de/api` |
| Authentication | None required |
| Rate Limit | 1 request/second recommended |
| Response Format | TSV, JSON, PNG, SVG |

### API Endpoints

| Operation | Endpoint | Example |
|-----------|----------|---------|
| Interactions | `/tsv/interactions` | `?identifier=aspirin&species=9606` |
| Network | `/tsv/network` | `?identifiers=aspirin%0dTP53&species=9606` |
| Resolve | `/json/resolve` | `?identifier=metformin` |
| Enrichment | `/json/enrichment` | `?identifiers=CIDm00004091&species=9606` |
| Image | `/image/network` | `?identifier=aspirin&species=9606` |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major version | Every 2-3 years |
| Data refresh | Annual |
| STRING sync | With STRING releases |

## Common Issues

- **Large files**: Use `-c` flag with wget to resume interrupted downloads
- **Memory limits**: Process line-by-line instead of loading full file
- **ID mapping**: Use chemical.aliases for PubChem/DrugBank conversion
- **Missing proteins**: Some older Ensembl IDs may be deprecated
- **Score 0 channels**: Zero means no evidence (not missing data)
- **Transfer scores**: May inflate confidence for understudied organisms

## API Rate Limits

```bash
# Respect rate limits: ~1 request/second
# For bulk queries, use downloads instead of API

# Batch API calls with delay
for chem in aspirin ibuprofen metformin; do
  curl -s "http://stitch.embl.de/api/json/interactions?identifier=$chem&species=9606" \
    >> all_interactions.json
  sleep 1
done
```

## Integration with STRING

```bash
# Download STRING protein interactions
wget https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz

# Combine STITCH + STRING for full network
# Note: STRING protein IDs are compatible with STITCH

# Create unified network
cat stitch_network.sif string_network.sif > combined_network.sif
```

## API Query Examples

```bash
# Resolve chemical name to STITCH ID
curl "http://stitch.embl.de/api/json/resolve?identifier=metformin"

# Get functional enrichment for chemical targets
curl "http://stitch.embl.de/api/json/enrichment?identifiers=CIDm00004091&species=9606"

# Get network image
curl "http://stitch.embl.de/api/image/network?identifier=aspirin&species=9606" \
  -o aspirin_network.png

# Get SVG format
curl "http://stitch.embl.de/api/svg/network?identifier=aspirin&species=9606" \
  -o aspirin_network.svg
```

## Data Format Reference

### protein_chemical.links format

```
chemical	protein	experimental	prediction	database	textmining	combined_score
CIDm00002244	9606.ENSP00000269305	900	0	700	650	976
```

### actions format

```
item_id_a	item_id_b	mode	action	is_directional	a_is_acting	score
CIDm00002244	9606.ENSP00000269305	inhibition	inhibitor	t	t	900
```

### chemical.aliases format

```
chemical	alias	source
CIDm00002244	aspirin	DrugBank
CIDm00002244	acetylsalicylic acid	ChEBI
CIDm00002244	DB00945	DrugBank_ID
```

## License

CC BY 4.0 - Free for any use with attribution to STITCH.
