---
id: download-symmap
title: "SymMap 2.0 Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# SymMap 2.0 Download Instructions

## Quick Start

```bash
# Visit the web interface for bulk downloads
# URL: http://www.symmap.org/
```

## Prerequisites

- Web browser for data access
- Python 3.8+ for data processing (optional)
- Cytoscape for network visualization (optional)

## Download Methods

### Primary: Web Interface

Access SymMap through the main portal:
- **URL**: http://www.symmap.org/
- **Data Access**: Bulk TSV files available

### Available Downloads

SymMap provides comprehensive symptom-phenotype mapping data:

| File Type | Format | Description |
|-----------|--------|-------------|
| Symptom Data | TSV | TCM symptoms with HPO mappings |
| Herb Data | TSV | Herbs with associations |
| Compound Data | TSV | Compounds with targets |
| Target Data | TSV | Gene targets |
| Disease Data | TSV | MESH/OMIM associations |
| Association Files | TSV | All relationship tables |

## File Inventory

| Category | Records | Description |
|----------|---------|-------------|
| TCM Symptoms | 1,717 | Traditional symptom patterns |
| Modern Phenotypes | 5,000+ | HPO mappings |
| TCM Herbs | 961 | Herb metadata |
| Compounds | 19,595 | Chemical compounds |
| Target Genes | 4,302 | Gene targets |
| Diseases | 5,235 | Disease associations |
| Symptom-Herb Links | 28,212 | Associations |
| Symptom-Gene Links | 322,073 | Associations |

## Data Structure

### Symptom Entity Fields
- Symptom ID (SMSY internal)
- TCM symptom name (Chinese)
- TCM symptom name (English)
- HPO mappings
- Associated herbs
- Associated genes
- Associated diseases

### Association Fields
- Source ID
- Target ID
- Association type
- Evidence level

## Download Steps

1. Navigate to http://www.symmap.org/
2. Access the "Download" section
3. Select data tables of interest
4. Download TSV files
5. Import into analysis tools

## Python Processing Example

```python
import pandas as pd

# Load symptom data
symptoms = pd.read_csv('symmap_symptoms.tsv', sep='\t')
print(f"TCM Symptoms: {len(symptoms)}")

# Load symptom-gene associations
symptom_genes = pd.read_csv('symmap_symptom_gene.tsv', sep='\t')
print(f"Symptom-Gene Links: {len(symptom_genes)}")

# Filter high-confidence associations
# (if confidence scores available)
```

## Network Analysis

SymMap provides Cytoscape-compatible network files:

```python
# Export for Cytoscape visualization
import networkx as nx

# Create symptom-target network
G = nx.Graph()
for _, row in symptom_genes.iterrows():
    G.add_edge(row['symptom_id'], row['gene_symbol'])

# Export to various formats
nx.write_graphml(G, 'symmap_network.graphml')
```

## Verification

Check against published statistics:

| Entity | Expected Count |
|--------|----------------|
| TCM Symptoms | 1,717 |
| HPO Phenotypes | 5,000+ |
| Herbs | 961 |
| Compounds | 19,595 |
| Targets | 4,302 |
| Diseases | 5,235 |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Current Version | SymMap 2.0 |
| Update Frequency | Periodic |
| Notification | Check website |

---

## Dataset Versions

### Current Release: SymMap 2.0

| Property | Value |
|----------|-------|
| Version | 2.0 |
| Release Date | 2020-01-01 |
| Total Size | ~200 MB |
| Focus | Symptom-Phenotype Mapping |

### Version Contents

| Component | Records | Description |
|-----------|---------|-------------|
| TCM Symptoms | 1,717 | Traditional symptom patterns |
| Modern Phenotypes | 5,000+ | HPO mappings |
| TCM Herbs | 961 | Herb metadata |
| Compounds | 19,595 | Chemical compounds |
| Target Genes | 4,302 | Gene targets |
| Diseases | 5,235 | Disease associations |
| Symptom-Gene Links | 322,073 | Associations |

### Previous Versions

| Version | Release | Status |
|---------|---------|--------|
| SymMap 1.0 | 2018-01-01 | Archived |

---

## Notes

- Unique symptom-centric perspective on TCM
- HPO mappings bridge traditional and modern medicine
- Network visualization supported
- Complements compound-centric databases
- Manual curation limits expansion speed
