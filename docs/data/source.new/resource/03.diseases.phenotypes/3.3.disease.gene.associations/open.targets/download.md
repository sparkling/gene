---
id: download-open-targets
title: "Open Targets Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# Open Targets Download Instructions

## Quick Start

```bash
# Download association scores (parquet format)
wget https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/associationByOverallDirect/

# Or use BigQuery for cloud analysis
```

## Prerequisites

- **wget** or **gsutil** for downloads
- **Apache Spark** or **DuckDB** for Parquet processing
- **Google Cloud SDK** for BigQuery access (optional)
- 50-500GB disk space depending on data scope

## No Registration Required

Open Targets data is freely available under CC0 (public domain).

## Download Methods

### Method 1: FTP Parquet Files (Recommended)

```bash
# Create download directory
mkdir -p open_targets && cd open_targets

# Download association data
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/associationByOverallDirect/

# Download target data
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/targets/

# Download disease data
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/diseases/

# Download evidence data
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/evidence/
```

### Method 2: JSON Lines Format

```bash
# Download JSON format (easier to parse, larger files)
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/json/associationByOverallDirect/

# Download targets as JSON
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/json/targets/
```

### Method 3: Google BigQuery (Cloud Analysis)

```sql
-- Access Open Targets in BigQuery (no download needed)
-- Project: open-targets-prod
-- Dataset: platform

-- Query associations
SELECT
  targetId,
  diseaseId,
  score
FROM `open-targets-prod.platform.associationByOverallDirect`
WHERE score > 0.5
LIMIT 1000;

-- Query targets
SELECT
  id,
  approvedSymbol,
  approvedName
FROM `open-targets-prod.platform.targets`
WHERE approvedSymbol = 'TP53';
```

### Method 4: GraphQL API (For Specific Queries)

```bash
# Query target information
curl -X POST https://api.platform.opentargets.org/api/v4/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ target(ensemblId: \"ENSG00000141510\") { id approvedSymbol approvedName } }"}' \
  -o tp53_info.json

# Query associations
curl -X POST https://api.platform.opentargets.org/api/v4/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ target(ensemblId: \"ENSG00000141510\") { associatedDiseases { rows { disease { id name } score } } } }"}' \
  -o tp53_diseases.json

# Query disease
curl -X POST https://api.platform.opentargets.org/api/v4/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ disease(efoId: \"EFO_0000311\") { id name associatedTargets { rows { target { id approvedSymbol } score } } } }"}' \
  -o cancer_targets.json
```

### Method 5: Google Cloud Storage

```bash
# Using gsutil (faster for large downloads)
gsutil -m cp -r \
  gs://open-targets-data-releases/latest/output/etl/parquet/associationByOverallDirect \
  ./associations/

# List available versions
gsutil ls gs://open-targets-data-releases/
```

### Method 6: Specific Data Types

```bash
# Evidence by data source
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/evidence/sourceId=chembl/

# Drug data
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/molecule/

# Known drugs
wget -r -np -nH --cut-dirs=7 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/knownDrugsAggregated/
```

## File Inventory

### Core Data Files

| Dataset | Size | Description |
|---------|------|-------------|
| associationByOverallDirect | ~5 GB | Target-disease associations |
| associationByDatasourceDirect | ~20 GB | Associations by source |
| evidence | ~100 GB | Raw evidence records |
| targets | ~2 GB | Target annotations |
| diseases | ~500 MB | Disease annotations |

### Drug Data

| Dataset | Size | Description |
|---------|------|-------------|
| molecule | ~500 MB | Drug/compound info |
| knownDrugsAggregated | ~200 MB | Approved drugs |
| mechanismOfAction | ~100 MB | Drug mechanisms |

### Supplementary

| Dataset | Size | Description |
|---------|------|-------------|
| interaction | ~1 GB | Protein interactions |
| reactome | ~500 MB | Pathway data |
| go | ~200 MB | GO annotations |

## Post-Download Processing

```bash
# Read Parquet with DuckDB
duckdb << 'EOF'
CREATE TABLE associations AS
SELECT * FROM parquet_scan('associationByOverallDirect/*.parquet');

-- Query high-confidence associations
SELECT targetId, diseaseId, score
FROM associations
WHERE score > 0.7
ORDER BY score DESC
LIMIT 100;

-- Export to CSV
COPY (SELECT * FROM associations WHERE score > 0.5)
TO 'high_conf_associations.csv' (HEADER, DELIMITER ',');
EOF

# Read with PySpark
python3 << 'EOF'
from pyspark.sql import SparkSession

spark = SparkSession.builder.appName("OpenTargets").getOrCreate()

# Load associations
associations = spark.read.parquet("associationByOverallDirect/")

# Filter and save
associations.filter("score > 0.5") \
  .select("targetId", "diseaseId", "score") \
  .write.csv("filtered_associations", header=True)
EOF

# Read with pandas (for smaller files)
python3 << 'EOF'
import pandas as pd

# Read single parquet partition
targets = pd.read_parquet("targets/part-00000.parquet")
print(targets.head())

# Filter and save
targets[['id', 'approvedSymbol', 'approvedName']].to_csv('targets.csv', index=False)
EOF
```

## Verification

```bash
# Check parquet files
parquet-tools schema associationByOverallDirect/part-00000.parquet

# Count records with DuckDB
duckdb -c "SELECT COUNT(*) FROM parquet_scan('associationByOverallDirect/*.parquet')"

# Check data integrity
python3 << 'EOF'
import pandas as pd
df = pd.read_parquet("associationByOverallDirect/")
print(f"Total associations: {len(df)}")
print(f"Unique targets: {df['targetId'].nunique()}")
print(f"Unique diseases: {df['diseaseId'].nunique()}")
EOF
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major releases | Quarterly |
| Data updates | Monthly |
| Bug fixes | As needed |

## Common Issues

- **Large downloads**: Use gsutil or rsync for reliability
- **Parquet reading**: Requires pyarrow, duckdb, or spark
- **ID formats**: Targets use Ensembl IDs; diseases use EFO IDs
- **Score interpretation**: 0-1 scale; higher = stronger evidence
- **Memory issues**: Process parquet files partition by partition

## Data Source IDs

| Source | Evidence Type |
|--------|---------------|
| chembl | Drug mechanisms |
| eva | ClinVar variants |
| gwas_catalog | GWAS associations |
| europepmc | Literature mining |
| expression_atlas | Expression data |
| reactome | Pathway evidence |
| sysbio | Systems biology |

## GraphQL Schema Exploration

```bash
# Get schema
curl -X POST https://api.platform.opentargets.org/api/v4/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ __schema { types { name } } }"}' | jq '.data.__schema.types[].name'

# Introspect target type
curl -X POST https://api.platform.opentargets.org/api/v4/graphql \
  -H "Content-Type: application/json" \
  -d '{"query": "{ __type(name: \"Target\") { fields { name type { name } } } }"}' | jq
```

## Related Resources

- [DisGeNET](../disgenet/) - Additional disease-gene data
- [ChEMBL](../../02.compounds.molecules/2.2.pharmaceuticals/chembl/) - Drug bioactivity
- [GWAS Catalog](../../01.genetics.genomics/1.5.expression.regulation/gwas.catalog/) - GWAS data
