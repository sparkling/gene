---
id: download-disgenet
title: "DisGeNET Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# DisGeNET Download Instructions

## Quick Start

```bash
# Download curated gene-disease associations (requires API key)
curl -X GET "https://www.disgenet.org/api/gda/gene/1017" \
  -H "Authorization: Bearer YOUR_API_KEY"
```

## Prerequisites

- **Free DisGeNET account** (required for downloads and API)
- **curl** for API access
- Approximately 1-5GB disk space

## Registration Process

### Step 1: Create Account

1. Navigate to https://www.disgenet.org/signup
2. Provide email and institution details
3. Select account type (Academic/Commercial)
4. Verify email address

### Step 2: Get API Key

1. Log in to DisGeNET
2. Go to Profile > API Key
3. Generate and copy your API key

## Access Tiers

| Tier | Access |
|------|--------|
| Academic (free) | Full database access, API, downloads |
| Commercial | Requires license agreement |

## Download Methods

### Method 1: Web Download (After Login)

1. Log in at https://www.disgenet.org
2. Navigate to Downloads: https://www.disgenet.org/downloads
3. Select file type and format
4. Download directly

### Method 2: REST API

```bash
# Set API key
DISGENET_KEY="your_api_key_here"

# Get gene-disease associations for a gene
curl -X GET "https://www.disgenet.org/api/gda/gene/1017" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -o gene_1017_gda.json

# Get disease-gene associations for a disease
curl -X GET "https://www.disgenet.org/api/gda/disease/C0024115" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -o disease_gda.json

# Get variant-disease associations
curl -X GET "https://www.disgenet.org/api/vda/variant/rs7412" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -o variant_vda.json

# Search diseases by name
curl -X GET "https://www.disgenet.org/api/disease/search/diabetes" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -o diabetes_diseases.json
```

### Method 3: Bulk Data Files (via API)

```bash
# Download all gene-disease associations
curl -X GET "https://www.disgenet.org/api/gda/source/ALL" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -H "Accept: text/tab-separated-values" \
  -o all_gda.tsv

# Download curated associations only
curl -X GET "https://www.disgenet.org/api/gda/source/CURATED" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -H "Accept: text/tab-separated-values" \
  -o curated_gda.tsv

# Download variant-disease associations
curl -X GET "https://www.disgenet.org/api/vda/source/ALL" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -H "Accept: text/tab-separated-values" \
  -o all_vda.tsv
```

### Method 4: Disease-Specific Downloads

```bash
# Get all genes associated with a disease
curl -X GET "https://www.disgenet.org/api/gda/disease/C0024115?limit=500" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -H "Accept: text/tab-separated-values" \
  -o breast_cancer_genes.tsv

# Get disease classes
curl -X GET "https://www.disgenet.org/api/disease/C0024115" \
  -H "Authorization: Bearer ${DISGENET_KEY}" \
  -o disease_info.json
```

### Method 5: R Package

```r
# Install disgenet2r
install.packages("disgenet2r")
library(disgenet2r)

# Set API key
disgenet_api_key(api_key = "your_api_key")

# Get gene-disease associations
gda <- disease2gene(disease = "C0024115", database = "ALL")

# Download and save
write.csv(gda, "breast_cancer_gda.csv")
```

### Method 6: Python Package

```python
# Install
# pip install disgenet-py

from disgenet import DisGeNET

client = DisGeNET(api_key="your_api_key")

# Get gene-disease associations
gda = client.get_gda_by_gene(gene_id=1017)

# Get disease associations
dga = client.get_gda_by_disease(disease_id="C0024115")

# Save to file
import pandas as pd
pd.DataFrame(gda).to_csv("gene_associations.csv")
```

## File Inventory

### Gene-Disease Associations

| File Type | Size | Description |
|-----------|------|-------------|
| All GDA (TSV) | ~500 MB | Complete associations |
| Curated GDA | ~50 MB | Expert curated only |
| Animal models | ~100 MB | Model organism data |

### Variant-Disease Associations

| File Type | Size | Description |
|-----------|------|-------------|
| All VDA (TSV) | ~300 MB | Complete VDA |
| GWAS Catalog VDA | ~100 MB | GWAS associations |
| ClinVar VDA | ~50 MB | ClinVar variants |

### Supplementary Data

| File Type | Size | Description |
|-----------|------|-------------|
| Disease mappings | ~10 MB | UMLS/DO mappings |
| Gene annotations | ~20 MB | Gene metadata |
| Source statistics | ~1 MB | Data source stats |

## Post-Download Processing

```bash
# Parse JSON response
cat gene_1017_gda.json | jq '.[] | {gene: .gene_symbol, disease: .disease_name, score: .score}'

# Filter by score (GDA score > 0.3)
awk -F'\t' 'NR==1 || $5 > 0.3' all_gda.tsv > high_confidence_gda.tsv

# Extract specific disease class
grep "Neoplasms" all_gda.tsv > cancer_associations.tsv

# Create gene list for enrichment
cut -f2 high_confidence_gda.tsv | sort -u > disease_genes.txt

# Merge with gene annotations
python3 << 'EOF'
import pandas as pd

gda = pd.read_csv('all_gda.tsv', sep='\t')

# Group by disease
disease_summary = gda.groupby('diseaseId').agg({
    'geneId': 'count',
    'score': 'mean'
}).rename(columns={'geneId': 'gene_count', 'score': 'avg_score'})

disease_summary.to_csv('disease_summary.csv')
EOF
```

## Verification

```bash
# Check TSV structure
head -5 all_gda.tsv

# Count associations
wc -l all_gda.tsv

# Count unique genes
cut -f2 all_gda.tsv | sort -u | wc -l

# Check score distribution
awk -F'\t' 'NR>1 {print int($5*10)/10}' all_gda.tsv | sort -n | uniq -c
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| DisGeNET v7.0 | 2020-05 | ~1.5 GB | Current |
| DisGeNET v8.0 | 2024-Q4 (planned) | TBD | In Development |
| v6.0 | 2019-05 | ~1 GB | Archived |

### Version Notes

DisGeNET v7.0 current statistics:
- 1.1M+ gene-disease associations
- 30,000+ genes
- 24,000+ diseases (UMLS concepts)
- 117,000+ variant-disease associations
- 14 source databases integrated

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://www.disgenet.org/api` |
| Rate Limit | Account-based |
| Auth Required | Yes (API key required) |
| Documentation | https://www.disgenet.org/api |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major releases | Annually |
| Data updates | Quarterly |
| Bug fixes | As needed |

## Common Issues

- **API rate limits**: Respect rate limits; batch requests appropriately
- **Large queries**: Use pagination for large result sets (limit parameter)
- **UMLS CUI format**: Disease IDs use UMLS CUI format (C0000000)
- **Score interpretation**: Higher GDA score = stronger evidence
- **Source filtering**: Filter by source for specific evidence types

## GDA Score Interpretation

| Score Range | Interpretation |
|-------------|----------------|
| 0.0 - 0.1 | Low evidence |
| 0.1 - 0.3 | Moderate evidence |
| 0.3 - 0.5 | Good evidence |
| 0.5 - 0.7 | Strong evidence |
| > 0.7 | Very strong evidence |

## Data Sources

| Source | Type |
|--------|------|
| CTD | Curated |
| UniProt | Curated |
| ClinGen | Curated |
| GWAS Catalog | GWAS |
| ClinVar | Clinical |
| BEFREE | Text mining |
| PsyGeNET | Psychiatric |

## API Endpoints Reference

| Endpoint | Description |
|----------|-------------|
| /gda/gene/{id} | GDA by gene |
| /gda/disease/{id} | GDA by disease |
| /vda/variant/{id} | VDA by variant |
| /gene/search/{term} | Search genes |
| /disease/search/{term} | Search diseases |

## Related Resources

- [Open Targets](../open.targets/) - Drug target validation
- [ClinVar](../../01.genetics.genomics/1.1.variant.repositories/clinvar/) - Clinical variants
- [OMIM](../../3.2.phenotype.databases/omim/) - Mendelian diseases
