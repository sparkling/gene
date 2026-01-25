---
id: download-orphanet
title: "Orphanet Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# Orphanet Download Instructions

## Quick Start

```bash
# Download disease nomenclature (Product 1)
wget https://www.orphadata.com/data/xml/en_product1.xml

# Download gene-disease associations (Product 6)
wget https://www.orphadata.com/data/xml/en_product6.xml
```

## Prerequisites

- **wget** or **curl** for downloads
- **XML parser** (lxml, ElementTree)
- Approximately 100MB disk space

## No Registration Required

Orphanet data is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: Orphadata XML Products

```bash
# Product 1: Disease nomenclature
wget https://www.orphadata.com/data/xml/en_product1.xml -O orphanet_diseases.xml

# Product 1 with cross-references
wget https://www.orphadata.com/data/xml/en_product1_cross_ref.xml -O orphanet_xrefs.xml

# Product 3: Disease classifications
wget https://www.orphadata.com/data/xml/en_product3.xml -O orphanet_classifications.xml

# Product 4: Phenotype annotations (HPO)
wget https://www.orphadata.com/data/xml/en_product4.xml -O orphanet_phenotypes.xml

# Product 6: Gene-disease associations
wget https://www.orphadata.com/data/xml/en_product6.xml -O orphanet_genes.xml

# Product 9: Epidemiology data
wget https://www.orphadata.com/data/xml/en_product9_prev.xml -O orphanet_prevalence.xml
wget https://www.orphadata.com/data/xml/en_product9_ages.xml -O orphanet_ages.xml

# All products listing
curl https://www.orphadata.com/data/xml/ -o orphadata_listing.html
```

### Method 2: ORDO Ontology (OWL)

```bash
# ORDO from OLS
wget "https://www.ebi.ac.uk/ols/ontologies/ordo/download" \
  -O ordo.owl

# ORDO from BioPortal
wget "https://data.bioontology.org/ontologies/ORDO/download?apikey=YOUR_API_KEY&download_format=csv" \
  -O ordo_bioportal.csv

# ORDO from OBO Foundry
wget http://purl.obolibrary.org/obo/ordo.owl -O ordo_purl.owl
```

### Method 3: REST API

```bash
# API base
BASE_URL="https://api.orphanet.org/rd-api"

# Get API documentation
curl "${BASE_URL}" -o orphanet_api_info.json

# Search diseases
curl "${BASE_URL}/diseases?name=marfan" -o search_marfan.json

# Get disease by ORPHA code
curl "${BASE_URL}/diseases/558" -o disease_558.json

# Get genes for a disease
curl "${BASE_URL}/diseases/558/genes" -o disease_558_genes.json

# Get phenotypes for a disease
curl "${BASE_URL}/diseases/558/phenotypes" -o disease_558_phenotypes.json

# List all gene associations
curl "${BASE_URL}/genes" -o all_genes.json
```

### Method 4: Multiple Languages

```bash
# Available: en, de, es, fr, it, nl, pl, pt
for lang in en de es fr; do
    wget "https://www.orphadata.com/data/xml/${lang}_product1.xml" \
      -O "orphanet_diseases_${lang}.xml"
done
```

### Method 5: Classifications by Type

```bash
# Download specific classification hierarchies
# Product 3 contains all classifications

# Parse and extract specific classification
python3 << 'EOF'
import xml.etree.ElementTree as ET

tree = ET.parse('orphanet_classifications.xml')
root = tree.getroot()

# Find classification types
for classification in root.findall('.//ClassificationList/Classification'):
    name = classification.find('Name').text if classification.find('Name') is not None else 'Unknown'
    print(f"Classification: {name}")
EOF
```

## File Inventory

### Orphadata Products

| Product | File | Description |
|---------|------|-------------|
| Product 1 | en_product1.xml | Disease nomenclature (~20MB) |
| Product 3 | en_product3.xml | Disease classifications (~30MB) |
| Product 4 | en_product4.xml | HPO phenotype annotations (~15MB) |
| Product 6 | en_product6.xml | Gene-disease associations (~10MB) |
| Product 9 | en_product9_*.xml | Epidemiology (~5MB) |

### ORDO Files

| File | Size | Description |
|------|------|-------------|
| ordo.owl | ~50MB | Full ORDO ontology |
| ordo.obo | ~30MB | OBO format |

## Post-Download Processing

```bash
# Parse disease nomenclature
python3 << 'EOF'
import xml.etree.ElementTree as ET
import pandas as pd

tree = ET.parse('orphanet_diseases.xml')
root = tree.getroot()

diseases = []
for disorder in root.findall('.//Disorder'):
    orpha_code = disorder.find('OrphaCode').text
    name = disorder.find('Name').text

    # Get synonyms
    synonyms = []
    for syn in disorder.findall('.//Synonym'):
        synonyms.append(syn.text)

    diseases.append({
        'orpha_code': orpha_code,
        'name': name,
        'synonyms': ';'.join(synonyms)
    })

df = pd.DataFrame(diseases)
df.to_csv('orphanet_diseases.tsv', sep='\t', index=False)
print(f"Total diseases: {len(df)}")
EOF

# Parse gene-disease associations
python3 << 'EOF'
import xml.etree.ElementTree as ET
import pandas as pd

tree = ET.parse('orphanet_genes.xml')
root = tree.getroot()

associations = []
for disorder in root.findall('.//Disorder'):
    orpha_code = disorder.find('OrphaCode').text
    disease_name = disorder.find('Name').text

    for gene_assoc in disorder.findall('.//DisorderGeneAssociation'):
        gene = gene_assoc.find('.//Gene')
        if gene is not None:
            symbol = gene.find('Symbol').text if gene.find('Symbol') is not None else ''
            gene_name = gene.find('Name').text if gene.find('Name') is not None else ''

            # Get association type
            assoc_type = gene_assoc.find('.//DisorderGeneAssociationType/Name')
            assoc_type_text = assoc_type.text if assoc_type is not None else ''

            associations.append({
                'orpha_code': orpha_code,
                'disease_name': disease_name,
                'gene_symbol': symbol,
                'gene_name': gene_name,
                'association_type': assoc_type_text
            })

df = pd.DataFrame(associations)
df.to_csv('orphanet_gene_disease.tsv', sep='\t', index=False)
print(f"Total associations: {len(df)}")
print(f"Unique genes: {df['gene_symbol'].nunique()}")
EOF

# Parse cross-references
python3 << 'EOF'
import xml.etree.ElementTree as ET
import pandas as pd

tree = ET.parse('orphanet_xrefs.xml')
root = tree.getroot()

xrefs = []
for disorder in root.findall('.//Disorder'):
    orpha_code = disorder.find('OrphaCode').text

    for xref in disorder.findall('.//ExternalReference'):
        source = xref.find('Source').text if xref.find('Source') is not None else ''
        ref_id = xref.find('Reference').text if xref.find('Reference') is not None else ''

        xrefs.append({
            'orpha_code': orpha_code,
            'source': source,
            'reference': ref_id
        })

df = pd.DataFrame(xrefs)
df.to_csv('orphanet_xrefs.tsv', sep='\t', index=False)

# Count by source
print(df['source'].value_counts())
EOF

# Build disease hierarchy from classifications
python3 << 'EOF'
import xml.etree.ElementTree as ET

tree = ET.parse('orphanet_classifications.xml')
root = tree.getroot()

def extract_hierarchy(node, parent_id='', depth=0):
    results = []
    for disorder in node.findall('./Disorder'):
        orpha_code = disorder.find('OrphaCode').text
        name = disorder.find('Name').text

        results.append({
            'orpha_code': orpha_code,
            'name': name,
            'parent': parent_id,
            'depth': depth
        })

        # Recurse into children
        for child_list in disorder.findall('./ClassificationNodeChildList'):
            results.extend(extract_hierarchy(child_list, orpha_code, depth + 1))

    return results

# Find first classification
class_list = root.find('.//ClassificationList')
if class_list:
    hierarchy = extract_hierarchy(class_list)
    print(f"Nodes in hierarchy: {len(hierarchy)}")
EOF
```

## Verification

```bash
# Check XML structure
head -50 orphanet_diseases.xml

# Count diseases
grep "<Disorder " orphanet_diseases.xml | wc -l

# Check gene associations
grep "<Gene " orphanet_genes.xml | wc -l

# Verify cross-references
grep "<ExternalReference>" orphanet_xrefs.xml | wc -l

# Check specific disease
grep -A 20 "OrphaCode>558<" orphanet_diseases.xml
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| Orphanet Jan 2026 | 2026-01 | ~100 MB | Current |
| Monthly releases | First week | Varies | Rolling |
| ORDO v4.5 | 2025-12 | ~50 MB | Current |

### Version Notes

Orphanet current database statistics:
- 6,300+ rare diseases catalogued
- 4,200+ genes with disease associations
- 12,000+ disease-gene associations
- 10,000+ disease classifications
- Multi-language support (8 languages)
- Expert-curated by medical professionals

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://api.orphanet.org/rd-api` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://www.orphadata.com |

## Update Schedule

| Product | Frequency |
|---------|-----------|
| Disease nomenclature | Monthly |
| Gene associations | Monthly |
| Phenotypes | Monthly |
| Classifications | Quarterly |
| ORDO ontology | Quarterly |

## Common Issues

- **XML namespaces**: Some products may have XML namespaces
- **Unicode**: Disease names may contain special characters
- **Obsolete codes**: Some ORPHA codes are deprecated
- **Version tracking**: Use product metadata for versioning
- **Language codes**: Ensure correct language prefix

## Cross-Reference Sources

| Source | Description |
|--------|-------------|
| OMIM | Mendelian diseases |
| ICD-10 | Clinical codes |
| UMLS | Unified Medical Language |
| MeSH | Medical Subject Headings |
| MedDRA | Regulatory terminology |
| SNOMED CT | Clinical terminology |

## Association Types

| Type | Meaning |
|------|---------|
| Disease-causing germline mutation(s) in | Causal gene |
| Disease-causing somatic mutation(s) in | Somatic causal |
| Modifying germline mutation in | Modifier gene |
| Major susceptibility factor in | Risk factor |
| Role in the phenotype of | Phenotype modifier |
| Candidate gene tested in | Research candidate |
| Biomarker tested in | Biomarker |

## Integration Examples

```bash
# Map Orphanet to OMIM
python3 << 'EOF'
import pandas as pd

xrefs = pd.read_csv('orphanet_xrefs.tsv', sep='\t')
omim_map = xrefs[xrefs['source'] == 'OMIM']
omim_map.to_csv('orphanet_omim_mapping.tsv', sep='\t', index=False)
print(f"OMIM mappings: {len(omim_map)}")
EOF

# Create disease-gene network
python3 << 'EOF'
import pandas as pd
import networkx as nx

assoc = pd.read_csv('orphanet_gene_disease.tsv', sep='\t')

# Create bipartite graph
G = nx.Graph()
for _, row in assoc.iterrows():
    G.add_node(row['orpha_code'], bipartite=0, type='disease')
    G.add_node(row['gene_symbol'], bipartite=1, type='gene')
    G.add_edge(row['orpha_code'], row['gene_symbol'])

print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")

# Find highly connected genes
gene_degrees = [(n, d) for n, d in G.degree() if G.nodes[n].get('type') == 'gene']
gene_degrees.sort(key=lambda x: x[1], reverse=True)
print("Top connected genes:")
for gene, degree in gene_degrees[:10]:
    print(f"  {gene}: {degree} diseases")
EOF
```

## Related Resources

- [Orphanet Phenotypes](../../3.2.phenotype.databases/orphanet.phenotypes/) - HPO annotations
- [OMIM](../../3.2.phenotype.databases/omim/) - Mendelian inheritance
- [HPO](../../3.2.phenotype.databases/hpo/) - Phenotype ontology
- [MONDO](../../3.1.disease.ontologies/mondo/) - Disease ontology with Orphanet mappings
