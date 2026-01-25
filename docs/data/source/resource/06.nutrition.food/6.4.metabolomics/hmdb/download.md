---
id: download-hmdb
title: "HMDB Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# HMDB (Human Metabolome Database) Download Instructions

## Quick Start

```bash
# Download metabolite XML database
wget https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip
unzip hmdb_metabolites.zip
```

## Prerequisites

- **wget** or **curl** for downloads
- **unzip** for extraction
- **XML parser** for processing
- Approximately 10-20GB disk space

## No Registration Required

HMDB data is freely available for academic and non-commercial use under a custom license.

## Download Methods

### Method 1: Metabolite Database

```bash
# Complete metabolite database (XML)
wget https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip
unzip hmdb_metabolites.zip

# SDF structures
wget https://hmdb.ca/system/downloads/current/structures.zip
unzip structures.zip
```

### Method 2: Protein Database

```bash
# Complete protein database
wget https://hmdb.ca/system/downloads/current/hmdb_proteins.zip
unzip hmdb_proteins.zip
```

### Method 3: Specific Data Files

```bash
# FASTA sequences
wget https://hmdb.ca/system/downloads/current/protein_sequences.zip

# MS/MS spectra
wget https://hmdb.ca/system/downloads/current/hmdb_predicted_msms_spectra.zip

# NMR spectra
wget https://hmdb.ca/system/downloads/current/hmdb_nmr_spectra.zip

# Urine metabolites
wget https://hmdb.ca/system/downloads/current/urine_metabolites.zip

# Serum metabolites
wget https://hmdb.ca/system/downloads/current/serum_metabolites.zip

# CSF metabolites
wget https://hmdb.ca/system/downloads/current/csf_metabolites.zip
```

### Method 4: Tabular Exports

```bash
# Download CSV/TSV exports from web interface
# Navigate to: https://hmdb.ca/downloads
# Select "Metabolite Data Table" and export format
```

### Method 5: API Access

```bash
# Get metabolite by HMDB ID
curl "https://hmdb.ca/metabolites/HMDB0000001.xml" -o HMDB0000001.xml

# Get metabolite JSON
curl "https://hmdb.ca/metabolites/HMDB0000001.json" -o HMDB0000001.json

# Search metabolites
curl "https://hmdb.ca/unearth/q?query=glucose&searcher=metabolites&button=" \
  -o glucose_search.html
```

### Method 6: Biofluid-Specific Downloads

```bash
# Download by biofluid
for biofluid in urine serum saliva csf feces sweat breast_milk; do
  wget "https://hmdb.ca/system/downloads/current/${biofluid}_metabolites.zip"
done
```

## File Inventory

### Core Database Files

| File | Size | Description |
|------|------|-------------|
| hmdb_metabolites.zip | ~2 GB | All metabolites (XML) |
| hmdb_proteins.zip | ~500 MB | All proteins (XML) |
| structures.zip | ~200 MB | SDF structures |

### Spectra Files

| File | Size | Description |
|------|------|-------------|
| hmdb_predicted_msms_spectra.zip | ~5 GB | Predicted MS/MS |
| hmdb_nmr_spectra.zip | ~1 GB | NMR spectra |
| hmdb_experimental_msms_spectra.zip | ~500 MB | Experimental MS/MS |

### Biofluid-Specific

| File | Size | Description |
|------|------|-------------|
| serum_metabolites.zip | ~1 GB | Serum metabolites |
| urine_metabolites.zip | ~1.5 GB | Urine metabolites |
| csf_metabolites.zip | ~500 MB | CSF metabolites |

## Post-Download Processing

```bash
# Parse XML with Python
python3 << 'EOF'
import xml.etree.ElementTree as ET
import pandas as pd
import os

data = []

# Parse main XML file
tree = ET.parse('hmdb_metabolites.xml')
root = tree.getroot()

ns = {'hmdb': 'http://www.hmdb.ca'}

for metabolite in root.findall('.//hmdb:metabolite', ns):
    accession = metabolite.find('hmdb:accession', ns)
    name = metabolite.find('hmdb:name', ns)
    formula = metabolite.find('hmdb:chemical_formula', ns)
    smiles = metabolite.find('hmdb:smiles', ns)

    data.append({
        'hmdb_id': accession.text if accession is not None else None,
        'name': name.text if name is not None else None,
        'formula': formula.text if formula is not None else None,
        'smiles': smiles.text if smiles is not None else None
    })

df = pd.DataFrame(data)
df.to_csv('hmdb_metabolites.tsv', sep='\t', index=False)
print(f"Processed {len(df)} metabolites")
EOF

# Extract specific fields using xmlstarlet
xmlstarlet sel -N h="http://www.hmdb.ca" \
  -t -m "//h:metabolite" \
  -v "h:accession" -o $'\t' \
  -v "h:name" -o $'\t' \
  -v "h:chemical_formula" -n \
  hmdb_metabolites.xml > metabolites_basic.tsv

# Process SDF structures
python3 << 'EOF'
from rdkit import Chem
from rdkit.Chem import Descriptors

suppl = Chem.SDMolSupplier('structures.sdf')

with open('hmdb_structures.tsv', 'w') as f:
    f.write("hmdb_id\tsmiles\tmw\n")
    for mol in suppl:
        if mol is None:
            continue
        try:
            hmdb_id = mol.GetProp('DATABASE_ID')
            smiles = Chem.MolToSmiles(mol)
            mw = Descriptors.MolWt(mol)
            f.write(f"{hmdb_id}\t{smiles}\t{mw:.2f}\n")
        except:
            continue
EOF
```

## Verification

```bash
# Check XML structure
head -100 hmdb_metabolites.xml

# Count metabolites
grep -c "<accession>" hmdb_metabolites.xml

# Check SDF
grep -c '^\$\$\$\$' structures.sdf

# Verify specific metabolite
grep -A 20 "HMDB0000001" hmdb_metabolites.xml
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major updates | Annually |
| New metabolites | Continuous |
| Spectra updates | Quarterly |

## Common Issues

- **Large XML files**: Use streaming parser (iterparse) for memory efficiency
- **Namespace handling**: HMDB XML uses namespace; include in queries
- **Missing data**: Not all fields populated for every metabolite
- **ID format**: Current format is HMDB00XXXXX (7 digits after HMDB)
- **Encoding**: UTF-8 with some special characters

## HMDB ID Format

| Format | Example | Notes |
|--------|---------|-------|
| Current | HMDB0000001 | 7 digits |
| Legacy | HMDB00001 | 5 digits |
| Secondary | HMDB0001234 | Cross-reference |

## Data Categories

| Category | Count (approx) |
|----------|----------------|
| Metabolites | 220,000+ |
| Proteins | 5,700+ |
| Diseases | 2,400+ |
| Pathways | 800+ |
| Biofluid locations | 16 |

## XML Element Reference

```xml
<metabolite>
  <accession>HMDB0000001</accession>
  <name>1-Methylhistidine</name>
  <chemical_formula>C7H11N3O2</chemical_formula>
  <smiles>CN1C=NC(CC(N)C(O)=O)=C1</smiles>
  <inchi>InChI=1S/C7H11N3O2/c...</inchi>
  <inchikey>BRMWTNUJHUMWMS-UHFFFAOYSA-N</inchikey>
  <state>Solid</state>
  <biofluid_locations>
    <biofluid>Blood</biofluid>
    <biofluid>Urine</biofluid>
  </biofluid_locations>
  <diseases>...</diseases>
  <pathways>...</pathways>
</metabolite>
```

## Integration Examples

```bash
# Cross-reference with PubChem
python3 << 'EOF'
import xml.etree.ElementTree as ET

ns = {'hmdb': 'http://www.hmdb.ca'}
tree = ET.parse('hmdb_metabolites.xml')
root = tree.getroot()

with open('hmdb_pubchem_map.tsv', 'w') as f:
    f.write("hmdb_id\tpubchem_cid\n")
    for metabolite in root.findall('.//hmdb:metabolite', ns):
        accession = metabolite.find('hmdb:accession', ns)
        pubchem = metabolite.find('.//hmdb:pubchem_compound_id', ns)
        if accession is not None and pubchem is not None and pubchem.text:
            f.write(f"{accession.text}\t{pubchem.text}\n")
EOF

# Extract disease associations
python3 << 'EOF'
import xml.etree.ElementTree as ET

ns = {'hmdb': 'http://www.hmdb.ca'}
tree = ET.parse('hmdb_metabolites.xml')
root = tree.getroot()

with open('hmdb_diseases.tsv', 'w') as f:
    f.write("hmdb_id\tmetabolite_name\tdisease\n")
    for metabolite in root.findall('.//hmdb:metabolite', ns):
        acc = metabolite.find('hmdb:accession', ns).text
        name = metabolite.find('hmdb:name', ns).text
        for disease in metabolite.findall('.//hmdb:disease/hmdb:name', ns):
            if disease.text:
                f.write(f"{acc}\t{name}\t{disease.text}\n")
EOF
```

## Related Resources

- [FooDB](../../6.1.food.composition/foodb/) - Food compounds
- [PubChem](../../02.compounds.molecules/2.6.chemical.ontology.classification/pubchem/) - Chemical database
- [KEGG](../../04.pathways.networks/4.1.metabolic.pathways/kegg/) - Metabolic pathways
