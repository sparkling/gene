---
id: download-chebi
title: "ChEBI Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# ChEBI (Chemical Entities of Biological Interest) Download Instructions

## Quick Start

```bash
# Download complete ChEBI (SDF format)
wget https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete.sdf.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **RDKit** or **Open Babel** for structure processing
- **OWL/OBO parser** for ontology files
- Approximately 5-10GB disk space

## No Registration Required

ChEBI data is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: Structure Files

```bash
# Complete ChEBI (all compounds)
wget https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete.sdf.gz

# 3-star compounds only (high quality)
wget https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete_3star.sdf.gz

# SMILES format
wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz

# InChI format
wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv
```

### Method 2: Ontology Files

```bash
# OBO format
wget https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz

# OWL format (full)
wget https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl.gz

# OWL lite (smaller)
wget https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.owl.gz

# OBOXML format
wget https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.xml.gz
```

### Method 3: Flat Files

```bash
# Chemical data
wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chemical_data.tsv.gz

# Compound origins (natural sources)
wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compound_origins.tsv.gz

# Names and synonyms
wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz

# Database links
wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv.gz

# Relations (ontology relations)
wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/relation.tsv.gz

# Vertice (compound attributes)
wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/vertice.tsv.gz
```

### Method 4: Web Services (LibChEBI)

```bash
# REST API - Get entity
curl "https://www.ebi.ac.uk/chebi/webServices/rest/ChEBI/getEntity/CHEBI:15377" \
  -o water.xml

# REST API - Search
curl "https://www.ebi.ac.uk/chebi/webServices/rest/ChEBI/getLiteEntity?search=aspirin&searchCategory=ALL&maximumResults=10" \
  -o aspirin_search.xml

# REST API - Get compound by InChI
curl "https://www.ebi.ac.uk/chebi/webServices/rest/ChEBI/getCompleteEntity?chebiId=CHEBI:15377" \
  -o water_complete.xml
```

### Method 5: Monthly Archives

```bash
# Download specific monthly release
wget https://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel224/SDF/ChEBI_complete.sdf.gz
```

## File Inventory

### Structure Files

| File | Size | Description |
|------|------|-------------|
| ChEBI_complete.sdf.gz | ~500 MB | All structures |
| ChEBI_complete_3star.sdf.gz | ~200 MB | High quality only |

### Ontology Files

| File | Size | Description |
|------|------|-------------|
| chebi.obo.gz | ~50 MB | OBO format |
| chebi.owl.gz | ~300 MB | Full OWL |
| chebi_lite.owl.gz | ~100 MB | Lite OWL |

### Flat Files

| File | Size | Description |
|------|------|-------------|
| compounds.tsv.gz | ~10 MB | Compound info |
| chemical_data.tsv.gz | ~20 MB | Chemical properties |
| names.tsv.gz | ~30 MB | Names/synonyms |
| database_accession.tsv.gz | ~10 MB | External links |
| relation.tsv.gz | ~5 MB | Ontology relations |

## Post-Download Processing

```bash
# Decompress
gunzip ChEBI_complete.sdf.gz

# Parse with RDKit
python3 << 'EOF'
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

suppl = Chem.SDMolSupplier('ChEBI_complete.sdf')

data = []
for mol in suppl:
    if mol is None:
        continue
    try:
        chebi_id = mol.GetProp('ChEBI ID')
        name = mol.GetProp('ChEBI Name')
        data.append({
            'chebi_id': chebi_id,
            'name': name,
            'mw': Descriptors.MolWt(mol),
            'smiles': Chem.MolToSmiles(mol)
        })
    except:
        continue

df = pd.DataFrame(data)
df.to_csv('chebi_compounds.tsv', sep='\t', index=False)
print(f"Processed {len(df)} compounds")
EOF

# Parse ontology with pronto
python3 << 'EOF'
import pronto

chebi = pronto.Ontology("chebi.obo.gz")

# Get term information
term = chebi["CHEBI:15377"]  # water
print(f"Name: {term.name}")
print(f"Definition: {term.definition}")

# Get all role relationships
for parent in term.superclasses(distance=1):
    print(f"is_a: {parent.name}")
EOF

# Extract flat file data
zcat compounds.tsv.gz | head -1000 > sample_compounds.tsv

# Parse relations
python3 << 'EOF'
import pandas as pd

relations = pd.read_csv('relation.tsv.gz', sep='\t', compression='gzip')
print(relations.columns)

# Filter for "is_a" relationships
is_a = relations[relations['TYPE'] == 'is_a']
print(f"is_a relations: {len(is_a)}")

# Filter for "has_role" relationships
has_role = relations[relations['TYPE'] == 'has_role']
print(f"has_role relations: {len(has_role)}")
EOF
```

## Verification

```bash
# Check SDF format
zcat ChEBI_complete.sdf.gz | head -100

# Count compounds
grep -c '^\$\$\$\$' ChEBI_complete.sdf

# Check ontology terms
grep "^\[Term\]" chebi.obo | wc -l

# Verify flat file structure
zcat compounds.tsv.gz | head -5
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| ChEBI Jan 2026 | 2026-01-06 | ~1 GB total | Current |
| Monthly releases | First week | Varies | Rolling |

### Version Notes

ChEBI current release features:
- Ontology provided in FULL, CORE, and LITE variants
- Available in OWL, OBO, and JSON formats (9 files total)
- Monthly release cycle with nightly ontology updates
- TSV, SDF, OBO, and OWL export formats

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://www.ebi.ac.uk/chebi/webServices/rest/ChEBI` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://www.ebi.ac.uk/chebi/webServices.do |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Monthly releases | Monthly |
| Nightly ontology | Daily |
| Bug fixes | As needed |

## Common Issues

- **Large ontology**: Use lite version for basic classification
- **Missing structures**: Not all ChEBI entries have structures
- **Star rating**: 3-star compounds have highest quality curation
- **ID format**: ChEBI IDs have format CHEBI:NNNNNNN
- **Ontology cycles**: Handle carefully in graph traversals

## ChEBI Ontology Structure

| Relationship | Description |
|--------------|-------------|
| is_a | Subsumption (parent class) |
| has_role | Biological/chemical role |
| has_part | Part-whole relationship |
| is_conjugate_base_of | Acid-base relationship |
| is_conjugate_acid_of | Acid-base relationship |
| is_tautomer_of | Tautomerism |
| has_functional_parent | Chemical derivation |

## Star Ratings

| Rating | Meaning |
|--------|---------|
| 3 stars | Fully manually annotated |
| 2 stars | Partially annotated |
| 1 star | Auto-submitted |

## Integration Examples

```bash
# Map ChEBI to PubChem
python3 << 'EOF'
import pandas as pd

xrefs = pd.read_csv('database_accession.tsv.gz', sep='\t', compression='gzip')
pubchem_map = xrefs[xrefs['TYPE'] == 'PUBCHEM'][['COMPOUND_ID', 'ACCESSION_NUMBER']]
pubchem_map.to_csv('chebi_pubchem_map.tsv', sep='\t', index=False)
EOF

# Extract drug roles
python3 << 'EOF'
import pronto

chebi = pronto.Ontology("chebi.obo.gz")
drug_role = chebi["CHEBI:23888"]  # drug

drugs = []
for term in chebi.terms():
    if drug_role in term.superclasses():
        drugs.append((term.id, term.name))

print(f"Found {len(drugs)} drug compounds")
EOF
```

## Related Resources

- [PubChem](../pubchem/) - Chemical database
- [DrugBank](../../2.2.pharmaceuticals/drugbank/) - Drug information
- [ChEMBL](../../2.2.pharmaceuticals/chembl/) - Bioactivity data
