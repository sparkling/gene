---
id: download-mondo
title: "MONDO Disease Ontology Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# MONDO Disease Ontology Download Instructions

## Quick Start

```bash
# Download MONDO ontology (OBO format)
wget http://purl.obolibrary.org/obo/mondo.obo

# Download with cross-references
wget http://purl.obolibrary.org/obo/mondo.owl
```

## Prerequisites

- **wget** or **curl** for downloads
- **OWL/OBO parser** (Pronto, OWLAPI, or ROBOT)
- Approximately 200MB disk space

## No Registration Required

MONDO is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: Official PURL Downloads

```bash
# OBO format (compact, widely supported)
wget http://purl.obolibrary.org/obo/mondo.obo

# OWL format (full OWL semantics)
wget http://purl.obolibrary.org/obo/mondo.owl

# JSON format
wget http://purl.obolibrary.org/obo/mondo.json

# OBO Graph JSON
wget http://purl.obolibrary.org/obo/mondo.owl -O mondo-base.owl
```

### Method 2: GitHub Releases

```bash
# List available releases
curl -s https://api.github.com/repos/monarch-initiative/mondo/releases | \
  jq '.[].tag_name'

# Download specific version
VERSION="v2024-01-01"
wget "https://github.com/monarch-initiative/mondo/releases/download/${VERSION}/mondo.obo"
wget "https://github.com/monarch-initiative/mondo/releases/download/${VERSION}/mondo.owl"

# Download with imports merged
wget "https://github.com/monarch-initiative/mondo/releases/download/${VERSION}/mondo-with-equivalents.owl"
```

### Method 3: Supplementary Files

```bash
# SSSOM mappings (to other disease resources)
wget https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.sssom.tsv

# Mondo to OMIM mapping
wget https://github.com/monarch-initiative/mondo/releases/latest/download/mondo_omim_mapping.sssom.tsv

# Mondo to DOID mapping
wget https://github.com/monarch-initiative/mondo/releases/latest/download/mondo_doid_mapping.sssom.tsv

# Cross-reference reports
wget https://github.com/monarch-initiative/mondo/releases/latest/download/mondo-xrefs.tsv
```

### Method 4: OLS API

```bash
# Search for terms
curl "https://www.ebi.ac.uk/ols/api/search?q=breast+cancer&ontology=mondo" \
  -o breast_cancer_search.json

# Get term details
curl "https://www.ebi.ac.uk/ols/api/ontologies/mondo/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FMONDO_0007254" \
  -o mondo_term.json

# Get children
curl "https://www.ebi.ac.uk/ols/api/ontologies/mondo/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FMONDO_0007254/children" \
  -o mondo_children.json
```

### Method 5: BioPortal API

```bash
# Search terms (requires API key)
BIOPORTAL_KEY="your_api_key"
curl "https://data.bioontology.org/search?q=diabetes&ontologies=MONDO&apikey=${BIOPORTAL_KEY}" \
  -o diabetes_search.json
```

## File Inventory

### Core Ontology Files

| File | Size | Description |
|------|------|-------------|
| mondo.obo | ~30 MB | OBO format |
| mondo.owl | ~100 MB | OWL format |
| mondo.json | ~60 MB | JSON format |
| mondo-with-equivalents.owl | ~150 MB | With merged imports |

### Mapping Files

| File | Size | Description |
|------|------|-------------|
| mondo.sssom.tsv | ~5 MB | All mappings (SSSOM) |
| mondo_omim_mapping.sssom.tsv | ~2 MB | OMIM mappings |
| mondo_doid_mapping.sssom.tsv | ~1 MB | DOID mappings |
| mondo-xrefs.tsv | ~10 MB | Cross-reference report |

## Post-Download Processing

```bash
# Parse OBO with Python
python3 << 'EOF'
import pronto

# Load ontology
mondo = pronto.Ontology("mondo.obo")

# Get term information
term = mondo["MONDO:0007254"]  # Breast cancer
print(f"Name: {term.name}")
print(f"Definition: {term.definition}")

# Get cross-references
for xref in term.xrefs:
    print(f"XRef: {xref.id}")

# Get synonyms
for syn in term.synonyms:
    print(f"Synonym: {syn.description}")
EOF

# Extract all disease names
grep "^name:" mondo.obo | cut -d: -f2- | sort > mondo_disease_names.txt

# Extract cross-references
python3 << 'EOF'
import pronto

mondo = pronto.Ontology("mondo.obo")

with open('mondo_xrefs.tsv', 'w') as f:
    f.write("mondo_id\tmondo_name\txref\n")
    for term in mondo.terms():
        for xref in term.xrefs:
            f.write(f"{term.id}\t{term.name}\t{xref.id}\n")
EOF

# Parse SSSOM mappings
python3 << 'EOF'
import pandas as pd

sssom = pd.read_csv('mondo.sssom.tsv', sep='\t', comment='#')
print(sssom.head())

# Filter for OMIM mappings
omim_maps = sssom[sssom['object_id'].str.startswith('OMIM:')]
print(f"OMIM mappings: {len(omim_maps)}")

# Filter for exact matches
exact_maps = sssom[sssom['predicate_id'] == 'skos:exactMatch']
print(f"Exact matches: {len(exact_maps)}")
EOF
```

## Verification

```bash
# Check OBO format
head -100 mondo.obo

# Count terms
grep "^\[Term\]" mondo.obo | wc -l

# Check for specific term
grep -A 10 "^id: MONDO:0007254" mondo.obo

# Verify mappings
head -5 mondo.sssom.tsv
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| v2025-01-06 | 2025-01-06 | ~150 MB | Current |
| v2024-12-03 | 2024-12-03 | ~145 MB | Archived |
| Monthly releases | First Monday | Varies | Rolling |

### Version Notes

MONDO current release features:
- 25,000+ disease concepts
- Comprehensive mappings to OMIM, Orphanet, DOID, EFO
- SSSOM standard for cross-ontology mappings
- Automated classification using Elk reasoner

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://www.ebi.ac.uk/ols/api/ontologies/mondo` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://oboacademy.github.io/obook |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Ontology updates | Monthly |
| Major revisions | Quarterly |
| Mapping updates | Monthly |

## Common Issues

- **ID format**: MONDO IDs use format MONDO:NNNNNNN (7 digits)
- **Obsolete terms**: Check is_obsolete tag; use replaced_by
- **Mapping predicates**: exactMatch vs closeMatch vs relatedMatch
- **Equivalence classes**: Use mondo-with-equivalents.owl for reasoning
- **Name spaces**: MONDO vs other ontology namespaces (OMIM, DOID)

## MONDO Structure

| Level | Example |
|-------|---------|
| Root | MONDO:0000001 (disease) |
| Category | MONDO:0045024 (cancer or benign tumor) |
| Type | MONDO:0005070 (neoplasm) |
| Specific | MONDO:0007254 (breast cancer) |
| Subtype | MONDO:0004989 (triple-negative breast cancer) |

## Cross-Reference Sources

| Source | Prefix | Description |
|--------|--------|-------------|
| OMIM | OMIM: | Mendelian diseases |
| Orphanet | Orphanet: | Rare diseases |
| DOID | DOID: | Disease Ontology |
| MESH | MESH: | Medical Subject Headings |
| ICD-10 | ICD10CM: | Clinical codes |
| UMLS | UMLS: | Unified Medical Language |

## SSSOM Mapping Predicates

| Predicate | Meaning |
|-----------|---------|
| skos:exactMatch | Equivalent concepts |
| skos:closeMatch | Very similar |
| skos:broadMatch | MONDO is more specific |
| skos:narrowMatch | MONDO is more general |
| skos:relatedMatch | Related but not equivalent |

## Integration Examples

```bash
# Map MONDO to OMIM
grep "OMIM:" mondo.obo | grep "xref:" > mondo_omim_xrefs.txt

# Create HPO-MONDO linkage
python3 << 'EOF'
import pronto

mondo = pronto.Ontology("mondo.obo")

# Find terms with HPO annotations
for term in mondo.terms():
    for xref in term.xrefs:
        if xref.id.startswith("HP:"):
            print(f"{term.id}\t{term.name}\t{xref.id}")
EOF
```

## Related Resources

- [HPO](../../3.2.phenotype.databases/hpo/) - Human Phenotype Ontology
- [OMIM](../../3.2.phenotype.databases/omim/) - Online Mendelian Inheritance
- [Orphanet](../../3.5.rare.orphan.diseases/orphanet/) - Rare disease database
