---
id: download-syngo
title: "SynGO Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# SynGO Download Instructions

## Quick Start

```bash
# Download SynGO gene annotations
wget https://syngoportal.org/data/SynGO_bulk_download_release_20231201.zip

# Download ontology
wget https://syngoportal.org/data/syngo.obo
```

## Prerequisites

- **wget** or **curl** for downloads
- **Python** with pandas for data processing
- Approximately 50MB disk space

## No Registration Required

SynGO is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: Bulk Downloads

```bash
# Complete bulk download (all data)
wget https://syngoportal.org/data/SynGO_bulk_download_release_20231201.zip

# SynGO ontology (OBO format)
wget https://syngoportal.org/data/syngo.obo

# Gene annotations (GAF format)
wget https://syngoportal.org/data/syngo_annotations.gaf

# Gene list (TSV)
wget https://syngoportal.org/data/syngo_genes.tsv
```

### Method 2: GitHub Repository

```bash
# Clone SynGO data repository
git clone https://github.com/syngoportal/SynGO-data.git
cd SynGO-data

# List available files
ls -la

# Ontology files
ls ontology/

# Annotation files
ls annotations/
```

### Method 3: Web Portal Export

```bash
# 1. Navigate to https://syngoportal.org/
# 2. Use Browse or Search functions
# 3. Select genes or terms of interest
# 4. Export via "Download" buttons

# Gene-specific export
# Search for gene, then download annotations
```

### Method 4: API Access (if available)

```bash
# Query specific genes
curl "https://syngoportal.org/api/genes/SYN1" \
  -H "Accept: application/json" \
  -o syn1_annotations.json

# Query ontology terms
curl "https://syngoportal.org/api/terms/SYNGO:synapse" \
  -H "Accept: application/json" \
  -o synapse_term.json
```

### Method 5: Enrichment Analysis Files

```bash
# Background gene sets for enrichment analysis
wget https://syngoportal.org/data/syngo_background_genes.txt

# Pre-computed enrichment datasets
wget https://syngoportal.org/data/syngo_enrichment_sets.zip
```

## File Inventory

### Core Files

| File | Size | Description |
|------|------|-------------|
| syngo.obo | ~500 KB | SynGO ontology |
| syngo_annotations.gaf | ~2 MB | Gene annotations (GAF) |
| syngo_genes.tsv | ~500 KB | Annotated gene list |
| SynGO_bulk_download.zip | ~20 MB | Complete download |

### Ontology Structure

| Branch | Description | Terms |
|--------|-------------|-------|
| Cellular Component | Synaptic localization | ~90 |
| Biological Process | Synaptic function | ~90 |

### Evidence Types

| Code | Description |
|------|-------------|
| ECO:0005593 | Biological assay |
| ECO:0005589 | Immunolocalization |
| ECO:0005644 | Electron microscopy |
| ECO:0006063 | Biochemical assay |

## Post-Download Processing

```bash
# Extract bulk download
unzip SynGO_bulk_download_release_20231201.zip -d syngo_data

# Parse gene annotations
python3 << 'EOF'
import pandas as pd

# Read gene list
genes = pd.read_csv('syngo_data/syngo_genes.tsv', sep='\t')

print(f"Total annotated genes: {len(genes)}")
print(f"\nColumns: {list(genes.columns)}")
print(f"\nSample genes:")
print(genes.head(10))
EOF

# Parse GAF annotations
python3 << 'EOF'
import pandas as pd

# GAF format columns
gaf_cols = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID',
            'DB_Reference', 'Evidence_Code', 'With_From', 'Aspect',
            'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type',
            'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']

# Read GAF (skip comment lines)
gaf = pd.read_csv('syngo_data/syngo_annotations.gaf', sep='\t',
                  comment='!', names=gaf_cols, low_memory=False)

print(f"Total annotations: {len(gaf)}")
print(f"Unique genes: {gaf['DB_Object_Symbol'].nunique()}")
print(f"Unique SynGO terms: {gaf['GO_ID'].nunique()}")

# Count by cellular component vs biological process
print(f"\nAspect distribution:")
print(gaf['Aspect'].value_counts())
EOF

# Parse ontology
python3 << 'EOF'
import pronto

# Load SynGO ontology
syngo = pronto.Ontology("syngo_data/syngo.obo")

print(f"Total terms: {len(list(syngo.terms()))}")

# List top-level terms
for term in syngo.terms():
    if len(list(term.superclasses(with_self=False, distance=1))) == 0:
        print(f"Root term: {term.id} - {term.name}")

# Find presynaptic terms
presynaptic = [t for t in syngo.terms() if 'presynaptic' in str(t.name).lower()]
print(f"\nPresynaptic terms: {len(presynaptic)}")
for t in presynaptic[:5]:
    print(f"  {t.id}: {t.name}")
EOF

# Build gene-term matrix
python3 << 'EOF'
import pandas as pd
from collections import defaultdict

# Read annotations
gaf = pd.read_csv('syngo_data/syngo_annotations.gaf', sep='\t',
                  comment='!', header=None, low_memory=False)
gaf.columns = ['DB', 'ID', 'Symbol', 'Qual', 'Term', 'Ref', 'Evid',
               'With', 'Aspect', 'Name', 'Syn', 'Type', 'Tax', 'Date',
               'By', 'Ext', 'Form']

# Build gene to terms mapping
gene_terms = defaultdict(set)
for _, row in gaf.iterrows():
    gene_terms[row['Symbol']].add(row['Term'])

# Convert to dataframe
gene_term_df = pd.DataFrame([
    {'gene': gene, 'term': term}
    for gene, terms in gene_terms.items()
    for term in terms
])

gene_term_df.to_csv('syngo_gene_terms.tsv', sep='\t', index=False)
print(f"Gene-term pairs: {len(gene_term_df)}")
EOF

# Extract synaptic compartment genes
python3 << 'EOF'
import pandas as pd

gaf = pd.read_csv('syngo_data/syngo_annotations.gaf', sep='\t',
                  comment='!', header=None, low_memory=False)
gaf.columns = ['DB', 'ID', 'Symbol', 'Qual', 'Term', 'Ref', 'Evid',
               'With', 'Aspect', 'Name', 'Syn', 'Type', 'Tax', 'Date',
               'By', 'Ext', 'Form']

# Filter for cellular component (C aspect)
cc = gaf[gaf['Aspect'] == 'C']

# Identify presynaptic vs postsynaptic
# This requires mapping to ontology structure
presynaptic_genes = cc[cc['Term'].str.contains('presynap', case=False, na=False)]['Symbol'].unique()
postsynaptic_genes = cc[cc['Term'].str.contains('postsynap', case=False, na=False)]['Symbol'].unique()

print(f"Presynaptic genes: {len(presynaptic_genes)}")
print(f"Postsynaptic genes: {len(postsynaptic_genes)}")

# Save gene lists
with open('presynaptic_genes.txt', 'w') as f:
    f.write('\n'.join(presynaptic_genes))
with open('postsynaptic_genes.txt', 'w') as f:
    f.write('\n'.join(postsynaptic_genes))
EOF
```

## Verification

```bash
# Check ontology format
head -50 syngo_data/syngo.obo

# Count ontology terms
grep "^\[Term\]" syngo_data/syngo.obo | wc -l

# Check GAF format
head -20 syngo_data/syngo_annotations.gaf

# Count unique genes
grep -v "^!" syngo_data/syngo_annotations.gaf | cut -f3 | sort -u | wc -l

# Verify gene list
wc -l syngo_data/syngo_genes.tsv
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major releases | Annual |
| New annotations | Continuous |
| Ontology updates | As needed |

## Common Issues

- **GO vs SynGO terms**: SynGO extends Gene Ontology for synapse
- **Evidence codes**: Different codes for different assay types
- **Multiple annotations**: Genes can have multiple terms
- **Human vs model**: Most evidence from mouse/rat models
- **Nomenclature**: Use HGNC symbols for human genes

## SynGO Hierarchy

### Cellular Component (CC)

| Level | Example Terms |
|-------|---------------|
| Synapse | synapse |
| Presynaptic | presynaptic membrane |
| Postsynaptic | postsynaptic density |
| Vesicle | synaptic vesicle |

### Biological Process (BP)

| Level | Example Terms |
|-------|---------------|
| Synaptic signaling | synaptic transmission |
| Vesicle cycle | synaptic vesicle cycle |
| Plasticity | synaptic plasticity |
| Organization | synapse organization |

## Enrichment Analysis

```bash
# Gene set enrichment with SynGO
python3 << 'EOF'
import pandas as pd
from scipy import stats

# Load your gene list (e.g., disease genes)
disease_genes = ['SHANK3', 'NRXN1', 'NLGN3', 'SYN1', 'GRIN2A']

# Load SynGO annotations
gaf = pd.read_csv('syngo_data/syngo_annotations.gaf', sep='\t',
                  comment='!', header=None)
syngo_genes = set(gaf[2])

# Overlap
overlap = set(disease_genes) & syngo_genes
print(f"Disease genes in SynGO: {len(overlap)}/{len(disease_genes)}")
print(f"Overlapping genes: {overlap}")

# For proper enrichment, use SynGO portal:
# https://syngoportal.org/enrichment.html
EOF
```

## Integration Examples

```bash
# Map SynGO to psychiatric GWAS
python3 << 'EOF'
import pandas as pd

# Load SynGO genes
syngo_genes = pd.read_csv('syngo_data/syngo_genes.tsv', sep='\t')
syngo_set = set(syngo_genes['gene_symbol'])

# Load psychiatric GWAS genes (from MAGMA or similar)
# gwas_genes = pd.read_csv('scz_magma_genes.tsv', sep='\t')
# gwas_set = set(gwas_genes['GENE'])

# overlap = syngo_set & gwas_set
# print(f"GWAS genes in SynGO: {len(overlap)}")
print("Placeholder for GWAS integration")
EOF

# Cross-reference with Allen Brain Atlas
python3 << 'EOF'
import pandas as pd

syngo_genes = pd.read_csv('syngo_data/syngo_genes.tsv', sep='\t')

# SynGO genes expressed in brain regions
# Requires Allen Brain Atlas expression data

print("Cross-reference with brain expression data")
print(f"SynGO genes for brain expression lookup: {len(syngo_genes)}")
EOF
```

## Related Resources

- [Allen Brain Atlas](../allen.brain.atlas/) - Brain expression data
- [PGC](../pgc/) - Psychiatric genetics
- [Gene Ontology](../../../../05.standards.ontologies/5.1.core.ontologies/gene.ontology/) - Parent ontology
- [UniProt](../../../../02.proteins.structures/2.1.protein.sequences/uniprot/) - Protein data
