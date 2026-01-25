---
id: download-dgidb
title: "DGIdb Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# DGIdb Download Instructions

## Quick Start

```bash
# Download interactions TSV
wget https://dgidb.org/data/monthly_tsvs/interactions.tsv

# Or use GraphQL API
curl -X POST https://dgidb.org/api/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ genes(names: [\"EGFR\"]) { name interactions { drug { name } } } }"}'
```

## Prerequisites

- **wget** or **curl** for downloads
- ~500 MB disk space for full dataset
- GraphQL client or curl for API access
- R/Bioconductor for rDGIdb package

## No Registration Required

Data is CC BY 4.0 licensed and freely downloadable.

## Download Methods

### Method 1: TSV Downloads (Recommended)

```bash
# Download all interactions
wget https://dgidb.org/data/monthly_tsvs/interactions.tsv

# Download genes with categories
wget https://dgidb.org/data/monthly_tsvs/genes.tsv

# Download drugs
wget https://dgidb.org/data/monthly_tsvs/drugs.tsv

# Download gene categories
wget https://dgidb.org/data/monthly_tsvs/categories.tsv

# Download sources
wget https://dgidb.org/data/monthly_tsvs/sources.tsv
```

### Method 2: GraphQL API

```bash
# Interactive GraphQL explorer
# https://dgidb.org/api/graphiql

# Query genes and interactions
curl -X POST https://dgidb.org/api/graphql \
  -H "Content-Type: application/json" \
  -d '{
    "query": "{ genes(names: [\"EGFR\", \"BRAF\"]) { name longName geneCategories { name } interactions { drug { name approved } interactionTypes { type } sources { sourceDbName } } } }"
  }' | jq

# Query drugs
curl -X POST https://dgidb.org/api/graphql \
  -H "Content-Type: application/json" \
  -d '{
    "query": "{ drugs(names: [\"gefitinib\"]) { name approved chemblId interactions { gene { name } interactionTypes { type } } } }"
  }' | jq

# Query by gene category
curl -X POST https://dgidb.org/api/graphql \
  -H "Content-Type: application/json" \
  -d '{
    "query": "{ genes(geneCategories: [\"KINASE\"]) { name longName interactionCount } }"
  }' | jq
```

### Method 3: R/Bioconductor Package

```r
# Install rDGIdb
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rDGIdb")

# Load and query
library(rDGIdb)

# Query gene interactions
result <- queryDGIdb(genes = c("EGFR", "BRAF"))
interactions <- resultSummary(result)
print(interactions)

# Get detailed results
details <- detailedResults(result)

# Query by category
kinases <- queryDGIdb(geneCategories = "KINASE")
```

### Method 4: Web Interface Export

```bash
# 1. Go to: https://dgidb.org
# 2. Enter gene names in search box
# 3. View results and click "Download" for TSV export
# 4. Use filters for specific interaction types or sources
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| interactions.tsv | ~50 MB | All drug-gene interactions |
| genes.tsv | ~5 MB | Gene information |
| drugs.tsv | ~3 MB | Drug information |
| categories.tsv | ~1 MB | Gene category assignments |
| sources.tsv | ~50 KB | Source database info |

## TSV Column Reference

### interactions.tsv

| Column | Description |
|--------|-------------|
| gene_name | Gene symbol |
| gene_claim_name | Source gene name |
| entrez_id | NCBI Entrez ID |
| interaction_claim_source | Source database |
| interaction_types | Interaction type(s) |
| drug_claim_name | Source drug name |
| drug_name | Normalized drug name |
| drug_chembl_id | ChEMBL ID |
| PMIDs | PubMed references |

## Post-Download Processing

```bash
# Preview data
head -20 interactions.tsv

# Count total interactions
wc -l interactions.tsv

# Count unique genes
cut -f1 interactions.tsv | sort | uniq | wc -l

# Count unique drugs
cut -f7 interactions.tsv | sort | uniq | wc -l

# Filter for specific interaction type
head -1 interactions.tsv > inhibitors.tsv
grep -i "inhibitor" interactions.tsv >> inhibitors.tsv

# Filter for specific source
head -1 interactions.tsv > chembl_data.tsv
grep "ChEMBL" interactions.tsv >> chembl_data.tsv

# Load into SQLite
sqlite3 dgidb.db << 'EOF'
CREATE TABLE interactions (
  gene_name TEXT,
  gene_claim_name TEXT,
  entrez_id TEXT,
  source TEXT,
  interaction_types TEXT,
  drug_claim_name TEXT,
  drug_name TEXT,
  chembl_id TEXT,
  pmids TEXT
);
.mode tabs
.import interactions.tsv interactions
DELETE FROM interactions WHERE gene_name='gene_name';
CREATE INDEX idx_gene ON interactions(gene_name);
CREATE INDEX idx_drug ON interactions(drug_name);
CREATE INDEX idx_source ON interactions(source);
EOF

# Query for gene
sqlite3 dgidb.db "SELECT DISTINCT drug_name, interaction_types FROM interactions WHERE gene_name='EGFR';"
```

## GraphQL Schema

```graphql
type Gene {
  name: String!
  longName: String
  conceptId: String
  geneCategories: [GeneCategory]
  interactions: [Interaction]
  interactionCount: Int
}

type Drug {
  name: String!
  conceptId: String
  approved: Boolean
  chemblId: String
  drugbankId: String
  interactions: [Interaction]
}

type Interaction {
  id: String!
  gene: Gene
  drug: Drug
  interactionTypes: [InteractionType]
  sources: [Source]
  pmids: [Int]
  score: Float
}
```

## Verification

```bash
# Check file integrity
file interactions.tsv
# Should be: ASCII text

# Verify column count
head -1 interactions.tsv | tr '\t' '\n' | wc -l

# Check for expected genes
grep -c "EGFR" interactions.tsv
grep -c "BRAF" interactions.tsv
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| DGIdb 5.0 | 2023-10 | ~100 MB | Current |
| DGIdb 4.x | 2021-2023 | ~80 MB | Archived |

### Version Notes

DGIdb 5.0 contains comprehensive drug-gene interaction data:
- 100,000+ drug-gene interactions
- 40,000+ genes with druggability info
- 10,000+ drugs
- 40+ source databases integrated
- GraphQL API for flexible queries

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://dgidb.org/api/graphql` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://dgidb.org/api/graphiql |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database updates | Monthly |
| TSV downloads | Monthly |
| API data | Real-time |

## API Usage Notes

```bash
# Batch gene queries (up to 100 genes)
curl -X POST https://dgidb.org/api/graphql \
  -H "Content-Type: application/json" \
  -d '{
    "query": "{ genes(names: [\"EGFR\", \"BRAF\", \"KRAS\", \"TP53\", \"PIK3CA\"]) { name interactionCount interactions { drug { name } } } }"
  }'

# Include pagination for large results
# Use first/after parameters in GraphQL
```

## Common Issues

- **Gene symbols**: Use HGNC symbols for best matching
- **Drug names**: Normalized names preferred over brand names
- **Multiple sources**: Same interaction may appear from multiple sources
- **Interaction types**: May be multiple types per interaction

## License

CC BY 4.0 - Free for any use with attribution
