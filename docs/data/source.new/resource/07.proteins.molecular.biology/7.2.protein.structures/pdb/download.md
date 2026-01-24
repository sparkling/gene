---
id: download-pdb
title: "PDB Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# PDB (Protein Data Bank) Download Instructions

## Quick Start

```bash
# Download single structure (mmCIF format, recommended)
wget https://files.rcsb.org/download/1tup.cif.gz

# Download single structure (PDB format, legacy)
wget https://files.rcsb.org/download/1tup.pdb.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **rsync** for bulk mirroring
- **PyMOL**, **ChimeraX**, or **VMD** for visualization
- 100GB-2TB storage for bulk downloads

## No Registration Required

PDB data is freely available under CC0 1.0 (public domain).

## Download Methods

### Method 1: Single Structure Download

```bash
# mmCIF format (recommended)
wget https://files.rcsb.org/download/1tup.cif.gz
gunzip 1tup.cif.gz

# PDB format (legacy, limited atom count)
wget https://files.rcsb.org/download/1tup.pdb.gz

# PDBML (XML format)
wget https://files.rcsb.org/download/1tup.xml.gz

# Biological assembly
wget https://files.rcsb.org/download/1tup-assembly1.cif.gz
```

### Method 2: Batch Download via File List

```bash
# Create list of PDB IDs
echo -e "1tup\n2hbs\n1a0a" > pdb_list.txt

# Download all structures
while read pdb; do
  wget "https://files.rcsb.org/download/${pdb}.cif.gz"
done < pdb_list.txt

# Using xargs for parallel downloads
cat pdb_list.txt | xargs -P 4 -I {} wget "https://files.rcsb.org/download/{}.cif.gz"
```

### Method 3: rsync Mirror (Recommended for Large Collections)

```bash
# Mirror entire PDB (mmCIF, divided format)
rsync -rlpt -v -z --delete \
  rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./pdb_mmcif/

# Mirror PDB format
rsync -rlpt -v -z --delete \
  rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./pdb_legacy/

# Mirror biological assemblies
rsync -rlpt -v -z --delete \
  rsync.rcsb.org::ftp_data/biounit/coordinates/divided/ ./pdb_biounit/
```

### Method 4: FTP Download

```bash
# Access FTP
ftp ftp.wwpdb.org
# cd /pub/pdb/data/structures/divided/mmCIF/
# cd tu/
# get 1tup.cif.gz

# Or via HTTP
wget https://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/tu/1tup.cif.gz
```

### Method 5: RCSB Search API + Download

```bash
# Search for human kinases
curl -X POST "https://search.rcsb.org/rcsbsearch/v2/query" \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "type": "group",
      "logical_operator": "and",
      "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entity_source_organism.taxonomy_lineage.name", "operator": "exact_match", "value": "Homo sapiens"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "struct_keywords.pdbx_keywords", "operator": "contains_words", "value": "kinase"}}
      ]
    },
    "return_type": "entry"
  }' | jq -r '.result_set[].identifier' > kinase_pdbs.txt

# Download results
while read pdb; do
  wget "https://files.rcsb.org/download/${pdb}.cif.gz"
done < kinase_pdbs.txt
```

### Method 6: Specific Data Types

```bash
# Structure factors
wget https://files.rcsb.org/download/1tup-sf.cif.gz

# Electron density maps (2Fo-Fc)
wget https://edmaps.rcsb.org/maps/1tup_2fofc.ccp4

# Validation reports
wget https://files.rcsb.org/pub/pdb/validation_reports/tu/1tup/1tup_validation.pdf

# Chemical component dictionary
wget https://files.rcsb.org/ligands/download/components.cif.gz
```

## File Inventory

### Coordinate Files

| Format | Size (typical) | Description |
|--------|----------------|-------------|
| .cif.gz | 100 KB - 10 MB | mmCIF (recommended) |
| .pdb.gz | 50 KB - 5 MB | Legacy PDB format |
| .xml.gz | 200 KB - 20 MB | PDBML/XML |
| .bcif.gz | 20 KB - 2 MB | Binary CIF |

### Full Database

| Directory | Size | Description |
|-----------|------|-------------|
| divided/mmCIF | ~100 GB | All mmCIF files |
| divided/pdb | ~80 GB | All PDB files |
| biounit | ~200 GB | Biological assemblies |
| obsolete | ~10 GB | Obsolete entries |

### Supplementary Data

| File Type | Size | Description |
|-----------|------|-------------|
| Structure factors | ~50 GB total | Experimental data |
| Validation reports | ~30 GB total | Quality metrics |
| Chemical components | ~500 MB | Ligand dictionary |

## Post-Download Processing

```bash
# Decompress
gunzip 1tup.cif.gz

# Extract FASTA sequence
python3 << 'EOF'
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder

parser = MMCIFParser()
structure = parser.get_structure("1TUP", "1tup.cif")
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
EOF

# Convert mmCIF to PDB
gemmi convert 1tup.cif 1tup.pdb

# Extract ligands
gemmi cif2json 1tup.cif | jq '.["_chem_comp"]'

# Extract coordinates as CSV
python3 << 'EOF'
from Bio.PDB import MMCIFParser
parser = MMCIFParser()
structure = parser.get_structure("1TUP", "1tup.cif")
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(f"{chain.id},{residue.resname},{residue.id[1]},{atom.name},{atom.coord[0]:.3f},{atom.coord[1]:.3f},{atom.coord[2]:.3f}")
EOF
```

## Verification

```bash
# Check file integrity
gzip -t 1tup.cif.gz

# Validate mmCIF
gemmi validate 1tup.cif

# Check resolution and method
grep "_refine.ls_d_res_high" 1tup.cif
grep "_exptl.method" 1tup.cif
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| New structures | Weekly (Wednesday) |
| Updates/corrections | Weekly |
| Obsolete structures | As needed |

## Common Issues

- **Large structures**: Use mmCIF; PDB format limited to 99,999 atoms
- **Multiple models**: NMR structures have multiple models; select appropriately
- **Biological assembly**: Download biounit for functional complex
- **Missing atoms**: Some structures have disordered regions
- **Ligand names**: Use Chemical Component Dictionary for definitions

## Quality Filters

```bash
# Filter by resolution (X-ray, < 2.5 A)
curl -X POST "https://search.rcsb.org/rcsbsearch/v2/query" \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "type": "terminal",
      "service": "text",
      "parameters": {
        "attribute": "rcsb_entry_info.resolution_combined",
        "operator": "less",
        "value": 2.5
      }
    }
  }'

# Filter by experimental method
# X-ray: "X-RAY DIFFRACTION"
# Cryo-EM: "ELECTRON MICROSCOPY"
# NMR: "SOLUTION NMR"
```

## REST API Access

```bash
# Get structure metadata
curl "https://data.rcsb.org/rest/v1/core/entry/1tup"

# Get polymer entities
curl "https://data.rcsb.org/rest/v1/core/polymer_entity/1tup/1"

# Get ligand information
curl "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/1tup/1"

# GraphQL query
curl -X POST "https://data.rcsb.org/graphql" \
  -H "Content-Type: application/json" \
  -d '{"query": "{ entry(entry_id: \"1TUP\") { rcsb_entry_info { resolution_combined } } }"}'
```

## Related Resources

- [AlphaFold DB](../alphafold.db/) - Predicted structures
- [SWISS-MODEL](../swiss.model/) - Homology models
- [UniProt](../../7.1.protein.sequences.annotations/uniprot/) - Sequence data
