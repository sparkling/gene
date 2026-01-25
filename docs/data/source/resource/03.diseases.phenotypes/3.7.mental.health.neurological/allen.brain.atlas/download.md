---
id: download-allen-brain-atlas
title: "Allen Brain Atlas Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# Allen Brain Atlas Download Instructions

## Quick Start

```bash
# Human Brain Atlas - Gene expression matrix
curl "https://api.brain-map.org/api/v2/well_known_file_download/178238387" \
  -o human_brain_expression.zip

# Get all available data via API
curl "https://api.brain-map.org/api/v2/data/query.json?criteria=model::Product" \
  -o allen_products.json
```

## Prerequisites

- **curl** or **wget** for downloads
- **Python** with pandas, numpy for data processing
- **Allen SDK** (optional but recommended)
- Approximately 50GB for full human brain data

## No Registration Required

Allen Brain Atlas data is openly accessible.

## Download Methods

### Method 1: Allen Brain Atlas API

```bash
# List available products
curl "https://api.brain-map.org/api/v2/data/query.json?criteria=model::Product" \
  -H "Accept: application/json" \
  -o allen_products.json

# Get human brain microarray data
curl "https://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,[organism_id\$eq2]" \
  -o human_donors.json

# Get gene information
curl "https://api.brain-map.org/api/v2/data/query.json?criteria=model::Gene,rma::criteria,[acronym\$eq'BDNF']" \
  -o bdnf_gene.json

# Get expression data for a gene
GENE_ID=18339  # BDNF
curl "https://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[gene_ids\$eq${GENE_ID}]" \
  -o bdnf_expression.json

# Get structure ontology
curl "https://api.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,[ontology_id\$eq1]" \
  -o brain_structures.json
```

### Method 2: Well-Known File Downloads

```bash
# Human Brain Atlas expression matrix
curl "https://api.brain-map.org/api/v2/well_known_file_download/178238387" \
  -o human_microarray_expression.zip

# Structure annotations
curl "https://api.brain-map.org/api/v2/well_known_file_download/178238316" \
  -o human_brain_annotations.zip

# Probe information
curl "https://api.brain-map.org/api/v2/well_known_file_download/178238359" \
  -o human_probe_info.zip
```

### Method 3: Allen SDK (Python)

```python
# Install: pip install allensdk

from allensdk.api.queries.rma_api import RmaApi
from allensdk.api.queries.human_brain_atlas_api import HumanBrainAtlasApi

# Initialize API
hba_api = HumanBrainAtlasApi()

# Get donors
donors = hba_api.get_donors()
print(f"Human brain donors: {len(donors)}")

# Get genes
genes = hba_api.get_genes()
print(f"Genes profiled: {len(genes)}")

# Get expression for specific gene
gene_symbol = 'BDNF'
expression = hba_api.get_microarray_expression([gene_symbol])
print(f"Expression data points: {len(expression)}")

# Save to file
import pandas as pd
pd.DataFrame(expression).to_csv('bdnf_expression.csv', index=False)
```

### Method 4: BrainSpan (Developmental)

```bash
# BrainSpan API
BASE_URL="https://api.brain-map.org/api/v2"

# Get BrainSpan samples
curl "${BASE_URL}/data/query.json?criteria=model::Specimen,rma::criteria,[data_set_id\$eq3]" \
  -o brainspan_samples.json

# Get developmental expression
# Download from: https://www.brainspan.org/static/download.html
wget "https://www.brainspan.org/api/v2/well_known_file_download/267666525" \
  -O brainspan_expression.zip
```

### Method 5: Cell Types Database

```bash
# Single-cell RNA-seq data
# Mouse brain cell types
curl "https://api.brain-map.org/api/v2/data/query.json?criteria=model::CellType" \
  -o cell_types.json

# Human cell types (MTG)
curl "https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq" \
  # Download via web interface

# Human M1 cell types
# Available at: https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
```

### Method 6: Reference Atlases

```bash
# Allen Human Brain Atlas reference volumes (NIfTI)
# Download from: https://human.brain-map.org/

# Mouse Brain Atlas
curl "https://api.brain-map.org/api/v2/well_known_file_download/197642854" \
  -o mouse_brain_atlas.nrrd

# Structure graph (ontology)
curl "https://api.brain-map.org/api/v2/structure_graph_download/1.json" \
  -o structure_graph.json
```

## File Inventory

### Human Brain Atlas

| File | Size | Description |
|------|------|-------------|
| MicroarrayExpression.csv | ~2 GB | Full expression matrix |
| SampleAnnot.csv | ~1 MB | Sample metadata |
| Probes.csv | ~5 MB | Probe information |
| Ontology.csv | ~100 KB | Structure ontology |

### BrainSpan

| File | Size | Description |
|------|------|-------------|
| expression_matrix.csv | ~500 MB | Developmental expression |
| columns_metadata.csv | ~5 MB | Sample information |
| rows_metadata.csv | ~5 MB | Gene information |

### Cell Types

| Dataset | Description |
|---------|-------------|
| Human MTG | Temporal cortex single-cell |
| Human M1 | Motor cortex single-cell |
| Mouse VISp | Visual cortex single-cell |

## Post-Download Processing

```bash
# Parse expression data
python3 << 'EOF'
import pandas as pd
import zipfile
import os

# Extract downloaded files
with zipfile.ZipFile('human_microarray_expression.zip', 'r') as z:
    z.extractall('human_brain_data')

# Load expression matrix
expression = pd.read_csv('human_brain_data/MicroarrayExpression.csv',
                         header=None, index_col=0)
samples = pd.read_csv('human_brain_data/SampleAnnot.csv')
probes = pd.read_csv('human_brain_data/Probes.csv')

print(f"Expression matrix: {expression.shape}")
print(f"Samples: {len(samples)}")
print(f"Probes: {len(probes)}")
EOF

# Get expression by brain region
python3 << 'EOF'
import pandas as pd

samples = pd.read_csv('human_brain_data/SampleAnnot.csv')
expression = pd.read_csv('human_brain_data/MicroarrayExpression.csv',
                         header=None, index_col=0)

# Map sample IDs to structures
region_means = samples.groupby('structure_name').apply(
    lambda x: expression[x.index].mean(axis=1)
)

print(f"Brain regions: {len(region_means.columns)}")
region_means.to_csv('expression_by_region.csv')
EOF

# Extract gene-specific expression
python3 << 'EOF'
import pandas as pd

probes = pd.read_csv('human_brain_data/Probes.csv')
expression = pd.read_csv('human_brain_data/MicroarrayExpression.csv',
                         header=None, index_col=0)

# Find probes for specific gene
gene = 'BDNF'
gene_probes = probes[probes['gene_symbol'] == gene]['probe_id'].values

if len(gene_probes) > 0:
    gene_expr = expression.loc[gene_probes]
    gene_expr.to_csv(f'{gene}_expression.csv')
    print(f"{gene} probes: {len(gene_probes)}")
EOF

# Analyze regional expression patterns
python3 << 'EOF'
from allensdk.api.queries.human_brain_atlas_api import HumanBrainAtlasApi
import pandas as pd

hba_api = HumanBrainAtlasApi()

# Get top expressed genes per region
structures = hba_api.get_structures()
structure_df = pd.DataFrame(structures)
print(f"Brain structures: {len(structure_df)}")
print(structure_df[['id', 'name', 'acronym']].head(20))
EOF
```

## Verification

```bash
# Check downloaded files
ls -la human_brain_data/

# Verify expression matrix dimensions
head -1 human_brain_data/MicroarrayExpression.csv | tr ',' '\n' | wc -l

# Check sample count
wc -l human_brain_data/SampleAnnot.csv

# Verify probe mapping
head human_brain_data/Probes.csv
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| Allen Brain Atlas 2024 | 2024 | ~100 GB | Current |
| BrainSpan | 2013 (updated) | ~5 GB | Stable |
| Human Cell Types | 2023 | ~50 GB | Current |
| Mouse Brain Atlas | Continuous | ~200 GB | Active |

### Version Notes

Allen Brain Atlas current resources:
- 6 human brain donors with microarray data
- 50,000+ gene expression profiles per brain
- 3,700+ sampling sites across brain regions
- BrainSpan: prenatal to adult development
- Single-cell datasets: 1M+ cells profiled

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://api.brain-map.org/api/v2` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://help.brain-map.org/display/api |

## Update Schedule

| Resource | Frequency |
|----------|-----------|
| Human Brain Atlas | Stable (archived) |
| BrainSpan | Stable |
| Cell Types | Ongoing additions |
| Mouse Brain | Periodic updates |

## Common Issues

- **Large files**: Expression matrices can be several GB
- **Probe redundancy**: Multiple probes per gene
- **Sample variability**: Different donors, regions
- **Coordinate systems**: Multiple reference spaces
- **Missing data**: Not all genes in all samples

## Structure Ontology

| Level | Example |
|-------|---------|
| Division | Telencephalon |
| Region | Cerebral cortex |
| Area | Frontal lobe |
| Structure | Hippocampus |
| Substructure | CA1 field |

## Integration Examples

```bash
# Map to disease genes
python3 << 'EOF'
import pandas as pd

# Load expression data
probes = pd.read_csv('human_brain_data/Probes.csv')
expression = pd.read_csv('human_brain_data/MicroarrayExpression.csv',
                         header=None, index_col=0)
samples = pd.read_csv('human_brain_data/SampleAnnot.csv')

# Disease genes (example: psychiatric)
psych_genes = ['BDNF', 'DISC1', 'NRG1', 'COMT', 'DRD2', 'SLC6A4']

# Get expression for disease genes
gene_probes = probes[probes['gene_symbol'].isin(psych_genes)]
disease_expr = expression.loc[gene_probes['probe_id']]

# Merge with probe info
result = gene_probes[['probe_id', 'gene_symbol']].merge(
    disease_expr.reset_index(),
    left_on='probe_id',
    right_on='index'
)
result.to_csv('psychiatric_genes_expression.csv', index=False)
print(f"Disease gene probes: {len(result)}")
EOF

# Cross-reference with GTEx
python3 << 'EOF'
# Compare Allen Brain expression with GTEx brain tissues
# Requires GTEx brain expression data

# This is a placeholder for cross-database comparison
print("Cross-reference analysis placeholder")
print("Combine Allen Brain (spatial resolution) with GTEx (population variation)")
EOF
```

## Related Resources

- [GTEx](../../../../01.genetics.genomics/1.5.expression.regulation/gtex/) - Multi-tissue expression
- [PGC](../pgc/) - Psychiatric genetics
- [SynGO](../syngo/) - Synaptic gene ontology
