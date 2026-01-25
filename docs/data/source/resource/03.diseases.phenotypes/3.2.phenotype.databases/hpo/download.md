---
id: download-hpo
title: "Human Phenotype Ontology Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# Human Phenotype Ontology (HPO) Download Instructions

## Quick Start

```bash
# Download HPO ontology (OBO format)
wget https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo

# Download gene-phenotype annotations
wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa
```

## Prerequisites

- **wget** or **curl** for downloads
- **OWL/OBO parser** (Pronto, OWLAPI, or similar)
- Approximately 500MB disk space

## No Registration Required

HPO data is freely available under public domain (CC0).

## Download Methods

### Method 1: Ontology Files

```bash
# OBO format (recommended for most use cases)
wget https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo

# OWL format (full semantic representation)
wget https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.owl

# JSON format
wget https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.json

# OBO with international translations
wget https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp-international.obo
```

### Method 2: Annotations (HPOA)

```bash
# Gene-to-phenotype annotations
wget http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt

# Phenotype-to-gene annotations
wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt

# Disease-to-phenotype annotations
wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa

# Medical action annotations
wget http://purl.obolibrary.org/obo/hp/hpoa/hpo_medical_actions.hpoa
```

### Method 3: GitHub Releases

```bash
# List available releases
curl -s https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases | \
  jq '.[].tag_name'

# Download specific version
VERSION="v2024-01-01"
wget "https://github.com/obophenotype/human-phenotype-ontology/releases/download/${VERSION}/hp.obo"
wget "https://github.com/obophenotype/human-phenotype-ontology/releases/download/${VERSION}/hp.owl"
```

### Method 4: Purl Stable URIs

```bash
# Official stable URIs (always point to latest)
wget http://purl.obolibrary.org/obo/hp.obo
wget http://purl.obolibrary.org/obo/hp.owl
wget http://purl.obolibrary.org/obo/hp.json
```

### Method 5: HPO API

```bash
# Search for terms
curl "https://hpo.jax.org/api/hpo/search/?q=seizure" -o seizure_terms.json

# Get term details
curl "https://hpo.jax.org/api/hpo/term/HP:0001250" -o hp_0001250.json

# Get genes for phenotype
curl "https://hpo.jax.org/api/hpo/term/HP:0001250/genes" -o seizure_genes.json

# Get diseases for phenotype
curl "https://hpo.jax.org/api/hpo/term/HP:0001250/diseases" -o seizure_diseases.json
```

## File Inventory

### Ontology Files

| File | Size | Description |
|------|------|-------------|
| hp.obo | ~10 MB | OBO format ontology |
| hp.owl | ~50 MB | OWL format ontology |
| hp.json | ~30 MB | JSON format ontology |
| hp-international.obo | ~15 MB | With translations |

### Annotation Files

| File | Size | Description |
|------|------|-------------|
| phenotype.hpoa | ~20 MB | Disease-phenotype annotations |
| genes_to_phenotype.txt | ~5 MB | Gene-phenotype mapping |
| phenotype_to_genes.txt | ~5 MB | Phenotype-gene mapping |

## Post-Download Processing

```bash
# Parse OBO with Python (pronto)
python3 << 'EOF'
import pronto

# Load ontology
hpo = pronto.Ontology("hp.obo")

# Get term information
term = hpo["HP:0001250"]  # Seizure
print(f"Name: {term.name}")
print(f"Definition: {term.definition}")
print(f"Parents: {[p.name for p in term.superclasses(distance=1)]}")
print(f"Children: {[c.name for c in term.subclasses(distance=1)]}")
EOF

# Extract all term names
grep "^name:" hp.obo | cut -d: -f2 | sort > hpo_term_names.txt

# Parse gene-phenotype annotations
python3 << 'EOF'
import pandas as pd

# Load gene-phenotype mapping
gp = pd.read_csv('genes_to_phenotype.txt', sep='\t', comment='#')
print(gp.head())

# Get genes for specific phenotype
seizure_genes = gp[gp['hpo_id'] == 'HP:0001250']['gene_symbol'].unique()
print(f"Genes associated with seizure: {len(seizure_genes)}")
EOF

# Create phenotype hierarchy visualization
python3 << 'EOF'
import pronto
import networkx as nx

hpo = pronto.Ontology("hp.obo")
G = nx.DiGraph()

# Build graph from root term
root = hpo["HP:0000118"]  # Phenotypic abnormality
for term in hpo.terms():
    for parent in term.superclasses(distance=1):
        if parent.id != term.id:
            G.add_edge(parent.id, term.id)

print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
nx.write_gexf(G, "hpo_hierarchy.gexf")
EOF
```

## Verification

```bash
# Check OBO format
head -50 hp.obo

# Count terms
grep "^\[Term\]" hp.obo | wc -l

# Verify annotations
wc -l genes_to_phenotype.txt

# Check annotation columns
head -1 genes_to_phenotype.txt
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| HPO 2025-01-13 | 2025-01-13 | ~100 MB | Current |
| HPO 2024-12-12 | 2024-12-12 | ~95 MB | Archived |
| Monthly releases | Second week | Varies | Rolling |

### Version Notes

HPO current release features:
- 18,000+ phenotypic terms
- 260,000+ phenotype annotations (diseases-HPO links)
- 5,400+ genes with HPO annotations
- Multiple languages supported
- CC0 public domain license

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://hpo.jax.org/api/hpo` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://hpo.jax.org/app/help/annotations |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Ontology updates | Monthly |
| Annotation updates | Quarterly |
| Major releases | Annually |

## Common Issues

- **ID format**: HPO IDs use format HP:NNNNNNN (7 digits)
- **Obsolete terms**: Check "is_obsolete: true" in OBO
- **Alternative IDs**: Some terms have been merged; check alt_id
- **Namespace confusion**: HP: prefix vs full URI
- **Annotation evidence**: Check evidence codes for reliability

## HPO ID Structure

| Level | Example | Description |
|-------|---------|-------------|
| Root | HP:0000001 | All |
| Top Level | HP:0000118 | Phenotypic abnormality |
| System | HP:0000707 | Nervous system |
| Category | HP:0001250 | Seizure |
| Specific | HP:0001263 | Global developmental delay |

## Annotation Evidence Codes

| Code | Meaning |
|------|---------|
| PCS | Published clinical study |
| IEA | Inferred from electronic annotation |
| TAS | Traceable author statement |
| ICE | Individual clinical experience |

## Integration with Other Resources

```bash
# Map HPO to OMIM
grep "OMIM:" phenotype.hpoa | cut -f1 | sort -u > hpo_with_omim.txt

# Cross-reference with Orphanet
grep "ORPHA:" phenotype.hpoa | cut -f1 | sort -u > hpo_with_orphanet.txt

# Link to UniProt via gene symbols
python3 << 'EOF'
import pandas as pd

gp = pd.read_csv('genes_to_phenotype.txt', sep='\t', comment='#')
# Gene symbols can be mapped to UniProt using UniProt ID mapping service
genes = gp['gene_symbol'].unique()
pd.Series(genes).to_csv('hpo_genes.txt', index=False, header=False)
EOF
```

## Related Resources

- [OMIM](../omim/) - Mendelian disease catalog
- [Orphanet](../../3.5.rare.orphan.diseases/orphanet/) - Rare disease data
- [MONDO](../../3.1.disease.ontologies/mondo/) - Disease ontology
