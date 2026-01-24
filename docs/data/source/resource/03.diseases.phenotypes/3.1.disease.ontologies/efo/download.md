---
id: download-efo
title: "Experimental Factor Ontology (EFO) Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# Experimental Factor Ontology (EFO) Download Instructions

## Quick Start

```bash
# Download EFO ontology (OWL format)
wget http://www.ebi.ac.uk/efo/efo.owl

# Download OBO format
wget http://www.ebi.ac.uk/efo/efo.obo
```

## Prerequisites

- **wget** or **curl** for downloads
- **OWL/OBO parser** (Pronto, OWLAPI, or ROBOT)
- Approximately 250MB disk space

## No Registration Required

EFO is freely available under Apache 2.0 license.

## Download Methods

### Method 1: Official EBI Downloads (PURL)

```bash
# OWL format (full semantics)
wget http://www.ebi.ac.uk/efo/efo.owl -O efo.owl

# OBO format (compact)
wget http://www.ebi.ac.uk/efo/efo.obo -O efo.obo

# JSON format
wget http://www.ebi.ac.uk/efo/efo.json -O efo.json

# EFO-lite (reduced version for GWAS Catalog)
wget http://www.ebi.ac.uk/efo/efo_lite.owl -O efo_lite.owl
```

### Method 2: GitHub Releases

```bash
# Clone repository for full history
git clone https://github.com/EBISPOT/efo.git
cd efo

# List releases
git tag -l | tail -20

# Checkout specific version
VERSION="v3.65.0"
git checkout ${VERSION}

# Download release assets
wget "https://github.com/EBISPOT/efo/releases/download/${VERSION}/efo.owl"
```

### Method 3: OBO Foundry

```bash
# Official OBO Foundry mirror
wget http://purl.obolibrary.org/obo/efo.owl -O efo_obofoundry.owl

# OBO format via PURL
wget http://purl.obolibrary.org/obo/efo.obo -O efo_obofoundry.obo
```

### Method 4: OLS API (Ontology Lookup Service)

```bash
# Search for terms
curl "https://www.ebi.ac.uk/ols4/api/search?q=diabetes&ontology=efo" \
  -o efo_diabetes_search.json

# Get term details
TERM_ID="EFO_0000685"
curl "https://www.ebi.ac.uk/ols4/api/ontologies/efo/terms?iri=http://www.ebi.ac.uk/efo/${TERM_ID}" \
  -o efo_term.json

# Get children of a term
curl "https://www.ebi.ac.uk/ols4/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000685/children" \
  -o efo_children.json

# Get ancestors
curl "https://www.ebi.ac.uk/ols4/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000685/ancestors" \
  -o efo_ancestors.json
```

### Method 5: SPARQL Endpoint

```bash
# Query EFO via SPARQL
curl -X POST "https://www.ebi.ac.uk/ols4/api/ontologies/efo/sparql" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "query=SELECT ?term ?label WHERE { ?term rdfs:label ?label } LIMIT 100" \
  -o efo_sparql_result.json
```

## File Inventory

### Core Ontology Files

| File | Size | Description |
|------|------|-------------|
| efo.owl | ~150 MB | Full OWL format with imports |
| efo.obo | ~50 MB | OBO format (no imports) |
| efo.json | ~80 MB | JSON format |
| efo_lite.owl | ~30 MB | Reduced version for GWAS Catalog |

### Supplementary Files

| File | Size | Description |
|------|------|-------------|
| efo_inferred.owl | ~180 MB | With inferred axioms |
| efo-edit.owl | ~100 MB | Editor version |
| mappings/ | Variable | Cross-reference mappings |

## Post-Download Processing

```bash
# Parse EFO with Python
python3 << 'EOF'
import pronto

# Load ontology
efo = pronto.Ontology("efo.obo")

# Get term information
term = efo["EFO:0000685"]  # Depression
print(f"Name: {term.name}")
print(f"Definition: {term.definition}")

# Get cross-references
for xref in term.xrefs:
    print(f"XRef: {xref.id}")

# Find disease terms
diseases = [t for t in efo.terms() if "disease" in str(t.name).lower()]
print(f"Disease terms: {len(diseases)}")
EOF

# Extract GWAS-relevant traits
python3 << 'EOF'
import pronto

efo = pronto.Ontology("efo.obo")

with open('efo_traits.tsv', 'w') as f:
    f.write("efo_id\tname\txrefs\n")
    for term in efo.terms():
        if term.id.startswith("EFO:"):
            xrefs = ";".join([x.id for x in term.xrefs])
            f.write(f"{term.id}\t{term.name}\t{xrefs}\n")
EOF

# Extract cross-references by source
grep "xref:" efo.obo | sort | uniq -c | sort -rn | head -20
```

## Verification

```bash
# Check file integrity
head -100 efo.owl

# Count classes (OBO format)
grep "^\[Term\]" efo.obo | wc -l

# Check for specific term
grep -A 15 "^id: EFO:0000685" efo.obo

# Verify cross-references
grep "xref: Orphanet:" efo.obo | wc -l
grep "xref: DOID:" efo.obo | wc -l
grep "xref: HP:" efo.obo | wc -l
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major releases | Monthly |
| Minor updates | Weekly |
| GWAS Catalog sync | Continuous |

## Common Issues

- **Large file size**: Use efo_lite.owl for lighter applications
- **Import errors**: OWL version includes imports; OBO does not
- **Obsolete terms**: Check `is_obsolete` flag; use `replaced_by`
- **ID format**: EFO IDs use format EFO_NNNNNNN (7 digits)
- **URL encoding**: OLS API requires double URL encoding for IRIs

## Cross-Reference Sources

| Source | Prefix | Description |
|--------|--------|-------------|
| Orphanet | Orphanet: | Rare diseases |
| DOID | DOID: | Disease Ontology |
| HP | HP: | Human Phenotype Ontology |
| MESH | MESH: | Medical Subject Headings |
| OMIM | OMIM: | Mendelian diseases |
| NCIt | NCIT: | NCI Thesaurus |

## Integration Examples

```bash
# Map EFO to GWAS Catalog traits
python3 << 'EOF'
import pronto

efo = pronto.Ontology("efo.obo")

# Find all terms with GWAS-related annotations
for term in efo.terms():
    if hasattr(term, 'annotations'):
        for ann in term.annotations:
            if 'gwas' in str(ann).lower():
                print(f"{term.id}\t{term.name}")
EOF

# Extract disease-phenotype hierarchy
python3 << 'EOF'
import pronto

efo = pronto.Ontology("efo.obo")

def get_ancestors(term, depth=0):
    indent = "  " * depth
    print(f"{indent}{term.id}: {term.name}")
    for parent in term.superclasses(with_self=False, distance=1):
        if depth < 3:
            get_ancestors(parent, depth + 1)

# Get hierarchy for a disease term
term = efo.get("EFO:0000685")
if term:
    get_ancestors(term)
EOF
```

## Related Resources

- [GWAS Catalog](../../../../01.genetics.genomics/1.5.expression.regulation/gwas.catalog/) - Primary EFO user
- [MONDO](../mondo/) - Disease ontology with EFO mappings
- [HPO](../../3.2.phenotype.databases/hpo/) - Phenotype ontology
