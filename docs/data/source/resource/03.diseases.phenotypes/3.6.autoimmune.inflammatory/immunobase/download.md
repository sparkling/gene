---
id: download-immunobase
title: "ImmunoBase Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# ImmunoBase Download Instructions

## Quick Start

```bash
# Download association data
curl "https://www.immunobase.org/downloads/regions-meta.tsv" -o immunobase_regions.tsv

# Download disease summary
curl "https://www.immunobase.org/downloads/diseases.tsv" -o immunobase_diseases.tsv
```

## Prerequisites

- **curl** or **wget** for downloads
- **TSV/CSV parser** (pandas recommended)
- Approximately 100MB disk space

## No Registration Required

ImmunoBase data is openly accessible.

## Download Methods

### Method 1: Direct Downloads

```bash
# Association regions
curl "https://www.immunobase.org/downloads/regions-meta.tsv" \
  -o immunobase_regions.tsv

# Disease information
curl "https://www.immunobase.org/downloads/diseases.tsv" \
  -o immunobase_diseases.tsv

# SNP associations
curl "https://www.immunobase.org/downloads/snps.tsv" \
  -o immunobase_snps.tsv

# Gene annotations
curl "https://www.immunobase.org/downloads/genes.tsv" \
  -o immunobase_genes.tsv

# Credible sets (fine-mapping)
curl "https://www.immunobase.org/downloads/credible-sets.tsv" \
  -o immunobase_credible_sets.tsv
```

### Method 2: Web Interface Export

```bash
# 1. Navigate to https://www.immunobase.org/
# 2. Use Region Browser or Disease views
# 3. Export visible data via "Download" buttons
# 4. Select format (TSV, BED, etc.)
```

### Method 3: Region-Specific Downloads

```bash
# Download data for specific chromosomal region
REGION="1p13"
curl "https://www.immunobase.org/region/${REGION}/download" \
  -o "region_${REGION}.tsv"

# Download by gene
GENE="PTPN22"
curl "https://www.immunobase.org/gene/${GENE}/download" \
  -o "gene_${GENE}.tsv"
```

### Method 4: Disease-Specific Downloads

```bash
# Available diseases
# T1D (Type 1 Diabetes), RA (Rheumatoid Arthritis), CEL (Celiac)
# MS (Multiple Sclerosis), CD (Crohn's), UC (Ulcerative Colitis)
# AS (Ankylosing Spondylitis), PSO (Psoriasis), etc.

for disease in T1D RA CEL MS CD UC AS PSO; do
    curl "https://www.immunobase.org/disease/${disease}/associations" \
      -o "immunobase_${disease}_associations.tsv"
done
```

### Method 5: GWAS Catalog Integration

```bash
# ImmunoBase studies are also in GWAS Catalog
# Download related studies

# Type 1 Diabetes studies
curl "https://www.ebi.ac.uk/gwas/api/search/downloads?q=efoTraits.shortForm:EFO_0001359" \
  -o gwas_catalog_t1d.tsv

# Rheumatoid Arthritis studies
curl "https://www.ebi.ac.uk/gwas/api/search/downloads?q=efoTraits.shortForm:EFO_0000685" \
  -o gwas_catalog_ra.tsv
```

## File Inventory

### Core Files

| File | Size | Description |
|------|------|-------------|
| regions-meta.tsv | ~5 MB | Associated genomic regions |
| diseases.tsv | ~100 KB | Disease metadata |
| snps.tsv | ~10 MB | SNP associations |
| genes.tsv | ~2 MB | Gene annotations |
| credible-sets.tsv | ~5 MB | Fine-mapping results |

### Per-Disease Files

| Disease | Abbreviation | Approximate Loci |
|---------|--------------|------------------|
| Type 1 Diabetes | T1D | 60+ |
| Rheumatoid Arthritis | RA | 100+ |
| Celiac Disease | CEL | 40+ |
| Multiple Sclerosis | MS | 200+ |
| Crohn's Disease | CD | 140+ |
| Ulcerative Colitis | UC | 100+ |
| Ankylosing Spondylitis | AS | 30+ |
| Psoriasis | PSO | 60+ |

## Post-Download Processing

```bash
# Parse association regions
python3 << 'EOF'
import pandas as pd

regions = pd.read_csv('immunobase_regions.tsv', sep='\t')

print(f"Total regions: {len(regions)}")
print(f"\nColumns: {list(regions.columns)}")

# Count regions per disease
if 'disease' in regions.columns:
    print(f"\nRegions per disease:")
    print(regions['disease'].value_counts())

# Find shared regions (multiple diseases)
if 'diseases' in regions.columns:
    shared = regions[regions['diseases'].str.contains(',', na=False)]
    print(f"\nShared regions: {len(shared)}")
EOF

# Extract lead SNPs
python3 << 'EOF'
import pandas as pd

snps = pd.read_csv('immunobase_snps.tsv', sep='\t')

# Filter for lead/index SNPs
if 'is_lead' in snps.columns:
    lead_snps = snps[snps['is_lead'] == True]
    print(f"Lead SNPs: {len(lead_snps)}")
    lead_snps.to_csv('immunobase_lead_snps.tsv', sep='\t', index=False)
EOF

# Analyze credible sets
python3 << 'EOF'
import pandas as pd

credible = pd.read_csv('immunobase_credible_sets.tsv', sep='\t')

# Group by region
if 'region_id' in credible.columns:
    region_stats = credible.groupby('region_id').agg({
        'snp': 'count',
        'posterior_probability': 'sum'
    })
    print(f"Credible sets: {len(region_stats)}")
    print(f"\nVariants per credible set (median): {region_stats['snp'].median()}")
EOF

# Create cross-disease comparison
python3 << 'EOF'
import pandas as pd
from collections import defaultdict

regions = pd.read_csv('immunobase_regions.tsv', sep='\t')

# Build region-disease matrix
disease_regions = defaultdict(set)
for _, row in regions.iterrows():
    diseases = row.get('diseases', row.get('disease', '')).split(',')
    region = row.get('region_id', row.get('region', ''))
    for d in diseases:
        disease_regions[d.strip()].add(region)

# Find overlaps
diseases = list(disease_regions.keys())
print("Cross-disease region sharing:")
for i, d1 in enumerate(diseases):
    for d2 in diseases[i+1:]:
        overlap = len(disease_regions[d1] & disease_regions[d2])
        if overlap > 5:
            print(f"  {d1} & {d2}: {overlap} shared regions")
EOF
```

## Verification

```bash
# Check file format
head -5 immunobase_regions.tsv

# Count records
wc -l immunobase_*.tsv

# Check for specific SNP
grep "rs2476601" immunobase_snps.tsv

# Verify disease coverage
cut -f2 immunobase_regions.tsv | sort | uniq -c | sort -rn
```

## Update Schedule

| Data Type | Frequency |
|-----------|-----------|
| New associations | As published |
| Fine-mapping updates | Periodic |
| Gene annotations | Annual |

## Common Issues

- **Genome build**: Data may be in GRCh37 or GRCh38; check metadata
- **SNP IDs**: Use dbSNP rsIDs for consistency
- **P-value formats**: May be log-transformed or scientific notation
- **Multi-disease regions**: Some regions associated with multiple diseases
- **Credible sets**: Posterior probabilities should sum to ~1

## Data Fields

### Association Data

| Field | Description |
|-------|-------------|
| rsid | dbSNP identifier |
| chr | Chromosome |
| position | Genomic position |
| p_value | Association p-value |
| odds_ratio | Effect size |
| disease | Disease abbreviation |
| study | Source study |

### Credible Set Data

| Field | Description |
|-------|-------------|
| snp | Variant identifier |
| posterior_probability | Fine-mapping probability |
| region_id | Associated region |
| is_lead | Lead/index variant flag |

## Integration Examples

```bash
# Map ImmunoBase to GWAS Catalog
python3 << 'EOF'
import pandas as pd

immunobase = pd.read_csv('immunobase_snps.tsv', sep='\t')
gwas = pd.read_csv('gwas_catalog_associations.tsv', sep='\t')

# Find overlapping SNPs
ib_snps = set(immunobase['rsid'].dropna())
gwas_snps = set(gwas['SNPS'].dropna())

overlap = ib_snps.intersection(gwas_snps)
print(f"ImmunoBase SNPs: {len(ib_snps)}")
print(f"GWAS Catalog SNPs: {len(gwas_snps)}")
print(f"Overlap: {len(overlap)}")
EOF

# Annotate with functional data
python3 << 'EOF'
import pandas as pd

snps = pd.read_csv('immunobase_snps.tsv', sep='\t')

# Add functional annotations from Ensembl VEP or similar
# This is a placeholder for functional annotation workflow
print(f"SNPs to annotate: {len(snps)}")
EOF

# Create BED file for credible sets
python3 << 'EOF'
import pandas as pd

credible = pd.read_csv('immunobase_credible_sets.tsv', sep='\t')

with open('credible_sets.bed', 'w') as f:
    for _, row in credible.iterrows():
        chrom = row.get('chr', row.get('chromosome', ''))
        pos = row.get('position', row.get('pos', 0))
        snp = row.get('snp', row.get('rsid', ''))
        pp = row.get('posterior_probability', 0)
        f.write(f"chr{chrom}\t{pos}\t{pos+1}\t{snp}\t{pp}\n")
EOF
```

## Related Resources

- [GWAS Catalog](../../../../01.genetics.genomics/1.5.expression.regulation/gwas.catalog/) - Full GWAS database
- [IPD-IMGT/HLA](../ipd.imgt.hla/) - HLA immunogenetics
- [Open Targets](../../3.3.disease.gene.associations/open.targets/) - Disease-gene associations
