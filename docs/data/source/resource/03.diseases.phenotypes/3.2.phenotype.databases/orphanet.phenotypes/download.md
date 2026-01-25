---
id: download-orphanet-phenotypes
title: "Orphanet Phenotype Annotations Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# Orphanet Phenotype Annotations Download Instructions

## Quick Start

```bash
# Download phenotype annotations (Product 4)
wget https://www.orphadata.com/data/xml/en_product4.xml

# Download with HOOM format
wget https://www.orphadata.com/data/HOOM/en_product4_HPO.xml
```

## Prerequisites

- **wget** or **curl** for downloads
- **XML parser** (lxml, ElementTree)
- Approximately 50MB disk space

## No Registration Required

Orphanet data is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: Orphadata Direct Downloads

```bash
# Phenotype annotations (Product 4)
wget https://www.orphadata.com/data/xml/en_product4.xml -O orphanet_phenotypes.xml

# Alternative: Gzipped version
wget https://www.orphadata.com/data/xml/en_product4.xml.gz

# HOOM format (HPO-Orphanet Ontological Module)
wget https://www.orphadata.com/data/HOOM/en_product4_HPO.xml -O orphanet_hoom.xml

# All Orphadata products listing
curl https://www.orphadata.com/data/ -o orphadata_listing.html
```

### Method 2: Related Orphanet Products

```bash
# Disease list (Product 1)
wget https://www.orphadata.com/data/xml/en_product1.xml -O orphanet_diseases.xml

# Gene-disease associations (Product 6)
wget https://www.orphadata.com/data/xml/en_product6.xml -O orphanet_genes.xml

# Epidemiology (Product 9)
wget https://www.orphadata.com/data/xml/en_product9_prev.xml -O orphanet_epidemiology.xml

# Classifications (Product 3)
wget https://www.orphadata.com/data/xml/en_product3.xml -O orphanet_classifications.xml

# Cross-references (Product 1 alternate)
wget https://www.orphadata.com/data/xml/en_product1_cross_ref.xml -O orphanet_xrefs.xml
```

### Method 3: Orphanet REST API

```bash
# API base URL
BASE_URL="https://api.orphanet.org/rd-api"

# List available endpoints
curl "${BASE_URL}" -o orphanet_api_endpoints.json

# Get disease information
ORPHA_CODE="558"  # Marfan syndrome
curl "${BASE_URL}/diseases/${ORPHA_CODE}" -o orphanet_disease_558.json

# Get phenotypes for a disease
curl "${BASE_URL}/diseases/${ORPHA_CODE}/phenotypes" -o orphanet_phenotypes_558.json

# Search diseases
curl "${BASE_URL}/diseases?name=marfan" -o orphanet_search_marfan.json
```

### Method 4: SPARQL/RDF Access

```bash
# ORDO (Orphanet Rare Disease Ontology) via OLS
wget https://www.ebi.ac.uk/ols/ontologies/ordo/download -O ordo.owl

# Query OLS SPARQL
curl "https://www.ebi.ac.uk/ols4/api/ontologies/ordo/terms" \
  -H "Accept: application/json" \
  -o ordo_terms.json
```

### Method 5: Multiple Languages

```bash
# Available languages: en, de, es, fr, it, nl, pl, pt
for lang in en de es fr it nl pl pt; do
    wget "https://www.orphadata.com/data/xml/${lang}_product4.xml" -O "orphanet_phenotypes_${lang}.xml"
done
```

## File Inventory

### Phenotype Files (Product 4)

| File | Size | Description |
|------|------|-------------|
| en_product4.xml | ~15 MB | English phenotype annotations |
| en_product4_HPO.xml | ~20 MB | HOOM format with HPO terms |

### Related Files

| Product | File | Description |
|---------|------|-------------|
| Product 1 | en_product1.xml | Disease nomenclature |
| Product 3 | en_product3.xml | Disease classifications |
| Product 6 | en_product6.xml | Gene-disease associations |
| Product 9 | en_product9_prev.xml | Epidemiology data |

## Post-Download Processing

```bash
# Parse phenotype annotations
python3 << 'EOF'
import xml.etree.ElementTree as ET
import pandas as pd

tree = ET.parse('orphanet_phenotypes.xml')
root = tree.getroot()

# Extract disease-phenotype associations
annotations = []
for disorder in root.findall('.//Disorder'):
    orpha_code = disorder.find('OrphaCode').text
    name = disorder.find('Name').text

    for hpo_assoc in disorder.findall('.//HPODisorderAssociation'):
        hpo_id = hpo_assoc.find('.//HPOId').text
        hpo_term = hpo_assoc.find('.//HPOTerm').text
        frequency = hpo_assoc.find('.//HPOFrequency/Name')
        freq_text = frequency.text if frequency is not None else 'Unknown'

        annotations.append({
            'orpha_code': orpha_code,
            'disease_name': name,
            'hpo_id': hpo_id,
            'hpo_term': hpo_term,
            'frequency': freq_text
        })

df = pd.DataFrame(annotations)
df.to_csv('orphanet_phenotypes.tsv', sep='\t', index=False)
print(f"Total annotations: {len(df)}")
print(f"Unique diseases: {df['orpha_code'].nunique()}")
print(f"Unique HPO terms: {df['hpo_id'].nunique()}")
EOF

# Extract frequency statistics
python3 << 'EOF'
import pandas as pd

df = pd.read_csv('orphanet_phenotypes.tsv', sep='\t')

# Frequency distribution
print("\nFrequency distribution:")
print(df['frequency'].value_counts())

# Map to percentage ranges
frequency_map = {
    'Obligate (100%)': 100,
    'Very frequent (99-80%)': 89.5,
    'Frequent (79-30%)': 54.5,
    'Occasional (29-5%)': 17,
    'Very rare (<4-1%)': 2.5,
    'Excluded (0%)': 0
}

df['freq_pct'] = df['frequency'].map(lambda x: next((v for k, v in frequency_map.items() if k in str(x)), None))
EOF

# Create disease phenotype profile
python3 << 'EOF'
import pandas as pd
from collections import defaultdict

df = pd.read_csv('orphanet_phenotypes.tsv', sep='\t')

# Get top phenotypes per disease
profiles = defaultdict(list)
for _, row in df.iterrows():
    if 'Obligate' in str(row['frequency']) or 'Very frequent' in str(row['frequency']):
        profiles[row['orpha_code']].append(row['hpo_id'])

# Save profiles
with open('disease_hpo_profiles.tsv', 'w') as f:
    f.write("orpha_code\ttop_hpo_terms\n")
    for orpha, hpos in profiles.items():
        f.write(f"{orpha}\t{';'.join(hpos)}\n")

print(f"Diseases with profiles: {len(profiles)}")
EOF

# Parse HOOM format for ontology relationships
python3 << 'EOF'
import xml.etree.ElementTree as ET

tree = ET.parse('orphanet_hoom.xml')
root = tree.getroot()

# HOOM includes ontological relationships
# Extract direct associations with evidence
for assoc in root.findall('.//Annotation')[:10]:
    print(ET.tostring(assoc, encoding='unicode')[:500])
EOF
```

## Verification

```bash
# Check XML structure
head -50 orphanet_phenotypes.xml

# Count diseases
grep "<Disorder " orphanet_phenotypes.xml | wc -l

# Count phenotype associations
grep "<HPODisorderAssociation>" orphanet_phenotypes.xml | wc -l

# Check specific disease (Marfan syndrome)
grep -A 50 "OrphaCode>558<" orphanet_phenotypes.xml | head -60

# Verify HPO term format
grep "HPOId" orphanet_phenotypes.xml | head -10
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| Orphanet Jan 2026 | 2026-01 | ~50 MB | Current |
| Monthly releases | First week | Varies | Rolling |

### Version Notes

Orphanet phenotype annotations current statistics:
- 6,700+ rare diseases with phenotype annotations
- 150,000+ disease-phenotype associations
- HPO terms used for standardization
- Frequency annotations (Obligate to Excluded)
- Multi-language support (8 languages)

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://api.orphanet.org/rd-api` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://www.orphadata.com |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Data exports | Monthly |
| API data | Real-time |
| ORDO ontology | Quarterly |

## Common Issues

- **XML namespace**: Product files may have XML namespaces
- **Language codes**: Use correct language prefix (en_, de_, etc.)
- **HPO version**: Phenotype terms align with HPO version
- **Frequency terms**: Use HPO frequency terms for standardization
- **API rate limits**: Be respectful with API queries

## Frequency Categories Mapping

| Orphanet Category | HPO Term | Range |
|-------------------|----------|-------|
| Obligate | HP:0040280 | 100% |
| Very frequent | HP:0040281 | 80-99% |
| Frequent | HP:0040282 | 30-79% |
| Occasional | HP:0040283 | 5-29% |
| Very rare | HP:0040284 | 1-4% |
| Excluded | HP:0040285 | 0% |

## Integration with HPO

```bash
# Download HPO for term details
wget http://purl.obolibrary.org/obo/hp.obo

# Map Orphanet phenotypes to HPO hierarchy
python3 << 'EOF'
import pronto
import pandas as pd

# Load HPO
hpo = pronto.Ontology("hp.obo")

# Load Orphanet phenotypes
df = pd.read_csv('orphanet_phenotypes.tsv', sep='\t')

# Enrich with HPO parents
def get_ancestors(hpo_id, ontology, max_depth=3):
    term = ontology.get(hpo_id)
    if term is None:
        return []
    ancestors = []
    for parent in term.superclasses(with_self=False, distance=max_depth):
        ancestors.append(parent.id)
    return ancestors

# Sample enrichment
sample = df.head(100).copy()
sample['hpo_ancestors'] = sample['hpo_id'].apply(lambda x: ';'.join(get_ancestors(x, hpo)))
sample.to_csv('orphanet_phenotypes_enriched.tsv', sep='\t', index=False)
EOF
```

## Related Resources

- [HPO](../hpo/) - Human Phenotype Ontology
- [Orphanet](../../3.5.rare.orphan.diseases/orphanet/) - Main Orphanet database
- [OMIM](../omim/) - Mendelian disease phenotypes
