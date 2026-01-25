---
id: download-decipher
title: "DECIPHER Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# DECIPHER Download Instructions

## Quick Start

```bash
# Access DECIPHER gene browser (public data)
curl "https://www.deciphergenomics.org/browser/dl/ddg2p" -o ddg2p.csv

# Download haploinsufficiency scores
curl "https://www.deciphergenomics.org/browser/dl/hi" -o hi_scores.tsv
```

## Prerequisites

- **Web browser** for portal access
- **curl** or **wget** for API downloads
- **DECIPHER account** for detailed variant data
- **Data Sharing Agreement** for patient-level data

## Registration Requirements

| Data Type | Access Level | Requirements |
|-----------|--------------|--------------|
| Gene scores (HI, pLI) | Open | None |
| DDG2P gene list | Open | None |
| Aggregate variant data | Open | None |
| Individual variants | Registered | DECIPHER account |
| Patient-level data | Controlled | Data Sharing Agreement |

## Download Methods

### Method 1: Public Gene Data Downloads

```bash
# DDG2P (Developmental Disorders Gene2Phenotype)
curl "https://www.deciphergenomics.org/browser/dl/ddg2p" \
  -o ddg2p.csv

# Haploinsufficiency scores
curl "https://www.deciphergenomics.org/browser/dl/hi" \
  -o hi_scores.tsv

# Probability of LoF intolerance
curl "https://www.deciphergenomics.org/browser/dl/pLI" \
  -o pli_scores.tsv

# LOEUF scores (constraint)
curl "https://www.deciphergenomics.org/browser/dl/oe_lof" \
  -o loeuf_scores.tsv
```

### Method 2: DECIPHER API (Registered Users)

```bash
# Search for variants by gene
# Requires API token from DECIPHER account
TOKEN="your_api_token"

# Search CNVs overlapping a gene
curl "https://www.deciphergenomics.org/api/v1/patients?gene=MECP2" \
  -H "Authorization: Bearer ${TOKEN}" \
  -o mecp2_patients.json

# Get variant details
curl "https://www.deciphergenomics.org/api/v1/patients/<patient_id>" \
  -H "Authorization: Bearer ${TOKEN}" \
  -o patient_details.json

# Search by genomic region
curl "https://www.deciphergenomics.org/api/v1/patients?chr=X&start=153000000&end=154000000" \
  -H "Authorization: Bearer ${TOKEN}" \
  -o region_patients.json
```

### Method 3: Genome Browser Export

```bash
# Export from DECIPHER browser view
# 1. Navigate to https://www.deciphergenomics.org/browser
# 2. Search for gene or region
# 3. Use "Export" button for visible data

# BED format export for coordinates
curl "https://www.deciphergenomics.org/api/v1/regions/bed?gene=MECP2" \
  -H "Authorization: Bearer ${TOKEN}" \
  -o mecp2_cnvs.bed
```

### Method 4: DDD Study Data

```bash
# DDD (Deciphering Developmental Disorders) data
# Available via European Genome-phenome Archive (EGA)

# Summary statistics (open)
# https://www.ebi.ac.uk/ega/studies/EGAS00001000775

# Individual data requires EGA application
```

### Method 5: Gene Panels

```bash
# DDG2P panels by category
curl "https://www.deciphergenomics.org/browser/dl/ddg2p/organ/brain" \
  -o ddg2p_brain.csv

curl "https://www.deciphergenomics.org/browser/dl/ddg2p/organ/heart" \
  -o ddg2p_heart.csv

# All categories
for organ in brain eye heart skeletal skin; do
  curl "https://www.deciphergenomics.org/browser/dl/ddg2p/organ/${organ}" \
    -o "ddg2p_${organ}.csv"
done
```

## File Inventory

### Public Downloads

| File | Size | Description |
|------|------|-------------|
| ddg2p.csv | ~2 MB | DD gene-phenotype associations |
| hi_scores.tsv | ~5 MB | Haploinsufficiency predictions |
| pli_scores.tsv | ~5 MB | pLI constraint scores |
| loeuf_scores.tsv | ~5 MB | LoF observed/expected |

### DDG2P Columns

| Column | Description |
|--------|-------------|
| gene symbol | HGNC gene symbol |
| gene mim | OMIM gene ID |
| disease name | Associated disorder |
| disease mim | OMIM disease ID |
| DDD category | Diagnostic confidence |
| allelic requirement | Inheritance mode |
| mutation consequence | Effect type |
| organ specificity | Affected systems |
| pmids | Supporting publications |

## Post-Download Processing

```bash
# Parse DDG2P data
python3 << 'EOF'
import pandas as pd

# Read DDG2P
ddg2p = pd.read_csv('ddg2p.csv')

print(f"Total gene-disease associations: {len(ddg2p)}")
print(f"Unique genes: {ddg2p['gene symbol'].nunique()}")

# Filter by confidence
confirmed = ddg2p[ddg2p['DDD category'] == 'confirmed']
print(f"Confirmed associations: {len(confirmed)}")

# Group by inheritance
print("\nInheritance patterns:")
print(ddg2p['allelic requirement'].value_counts())
EOF

# Merge constraint scores
python3 << 'EOF'
import pandas as pd

# Load scores
hi = pd.read_csv('hi_scores.tsv', sep='\t')
pli = pd.read_csv('pli_scores.tsv', sep='\t')

# Merge on gene
merged = hi.merge(pli, on='gene', how='outer')
merged.to_csv('gene_constraint_scores.tsv', sep='\t', index=False)

# Find high-confidence constrained genes
constrained = merged[(merged['pLI'] > 0.9) & (merged['HI'] < 10)]
print(f"Highly constrained genes: {len(constrained)}")
EOF

# Create BED file from CNV data
python3 << 'EOF'
import json

with open('mecp2_patients.json') as f:
    data = json.load(f)

with open('cnvs.bed', 'w') as f:
    for patient in data.get('patients', []):
        for var in patient.get('variants', []):
            chrom = var.get('chromosome', '')
            start = var.get('start', '')
            end = var.get('end', '')
            vtype = var.get('type', 'CNV')
            f.write(f"chr{chrom}\t{start}\t{end}\t{vtype}\n")
EOF
```

## Verification

```bash
# Check DDG2P format
head -5 ddg2p.csv

# Count genes per category
cut -d',' -f5 ddg2p.csv | sort | uniq -c

# Verify constraint scores
head hi_scores.tsv
wc -l hi_scores.tsv

# Check gene overlap
cut -f1 hi_scores.tsv | sort -u > hi_genes.txt
cut -d',' -f1 ddg2p.csv | sort -u > ddg2p_genes.txt
comm -12 hi_genes.txt ddg2p_genes.txt | wc -l
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| DECIPHER 2025 | 2025-01 | ~50 MB public | Current |
| DDG2P v4.0 | 2024-12 | ~2 MB | Current |
| Continuous updates | Real-time | Varies | Active |

### Version Notes

DECIPHER current statistics:
- 50,000+ patient submissions worldwide
- 2,500+ DD genes in DDG2P
- Haploinsufficiency scores for 19,000+ genes
- pLI constraint scores for protein-coding genes
- Active DDD (Deciphering Developmental Disorders) cohort

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://www.deciphergenomics.org/api/v1` |
| Rate Limit | Account-based |
| Auth Required | Yes (for patient data) |
| Documentation | https://www.deciphergenomics.org/about/data-sharing |

## Update Schedule

| Data Type | Frequency |
|-----------|-----------|
| Patient submissions | Continuous |
| DDG2P | Quarterly |
| Constraint scores | Annual |

## Common Issues

- **Data access tiers**: Different access levels for different data
- **HPO terms**: Use HPO IDs for standardized phenotypes
- **Variant coordinates**: Uses GRCh38 coordinates
- **CNV interpretation**: Combine with population frequency data
- **Patient consent**: Follow data sharing agreement terms

## DDD Categories

| Category | Meaning |
|----------|---------|
| confirmed | Strong gene-disease evidence |
| probable | Good evidence, needs more |
| possible | Emerging evidence |
| both RD and IF | Rare disease and cancer predisposition |

## Constraint Score Interpretation

| Score | Range | Interpretation |
|-------|-------|----------------|
| pLI | 0-1 | >0.9 = LoF intolerant |
| HI score | 0-100 | <10 = likely haploinsufficient |
| LOEUF | 0-2 | <0.35 = constrained |

## Integration Examples

```bash
# Map DECIPHER genes to HPO
python3 << 'EOF'
import pandas as pd

ddg2p = pd.read_csv('ddg2p.csv')

# Extract unique HPO terms
hpo_terms = set()
for terms in ddg2p['hpo terms'].dropna():
    for term in terms.split(';'):
        hpo_terms.add(term.strip())

print(f"Unique HPO terms in DDG2P: {len(hpo_terms)}")
EOF

# Cross-reference with ClinVar
python3 << 'EOF'
import pandas as pd

# Load DDG2P genes
ddg2p = pd.read_csv('ddg2p.csv')
ddg2p_genes = set(ddg2p['gene symbol'].dropna())

# Load ClinVar gene summary (if available)
# clinvar = pd.read_csv('clinvar_gene_specific.txt', sep='\t')
# overlap = ddg2p_genes.intersection(set(clinvar['Symbol']))
# print(f"Genes in both: {len(overlap)}")
EOF
```

## Related Resources

- [Orphanet](../orphanet/) - Rare disease portal
- [OMIM](../../3.2.phenotype.databases/omim/) - Mendelian genetics
- [ClinVar](../../../../01.genetics.genomics/1.1.variant.repositories/clinvar/) - Clinical variants
- [gnomAD](../../../../01.genetics.genomics/1.3.population.genetics/gnomad/) - Population constraint
