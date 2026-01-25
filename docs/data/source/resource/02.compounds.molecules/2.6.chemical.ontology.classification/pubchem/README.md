---
id: pubchem
title: "PubChem - Chemical Information Repository"
type: source
parent: ../README.md
tier: 1
status: active
category: compounds.molecules
subcategory: chemical.ontology.classification
tags:
  - chemical-database
  - compound-repository
  - bioassay
  - structure
  - ncbi
---

# PubChem - Chemical Information Repository

## Overview

PubChem is the world's largest collection of freely accessible chemical information, maintained by the National Center for Biotechnology Information (NCBI). It provides information on chemical structures, identifiers, chemical and physical properties, biological activities, patents, health, safety, toxicity data, and links to external databases.

PubChem organizes data into three linked databases: Substance (depositor-provided chemical samples), Compound (unique chemical structures), and BioAssay (biological test results). This structure enables tracking of a compound's testing history across assays and identification of active compounds.

With over 100 million compounds and billions of bioactivity data points, PubChem serves as a critical hub for chemical and biological information, integrating data from hundreds of sources worldwide.

## Key Statistics

| Metric | Value |
|--------|-------|
| Compounds | 115+ million |
| Substances | 300+ million |
| BioAssays | 1.5+ million |
| Bioactivity Data Points | 290+ million |
| Data Sources | 850+ |
| Gene Targets | 70,000+ |
| Patents | 35+ million |

## Primary Use Cases

1. **Compound Information** - Look up chemical properties and identifiers
2. **Bioactivity Research** - Find biological test results for compounds
3. **Cross-Reference Hub** - Link chemical IDs across databases
4. **Structure Search** - Find similar or substructure-matching compounds
5. **Literature Mining** - Compound-publication associations

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| CID (Compound) | Integer | 2244 (Aspirin) |
| SID (Substance) | Integer | Depositor submission |
| AID (BioAssay) | Integer | Assay identifier |
| InChI Key | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| SMILES | Variable | CC(=O)OC1=CC=CC=C1C(=O)O |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://pubchem.ncbi.nlm.nih.gov | Search and browse |
| PUG-REST | https://pubchem.ncbi.nlm.nih.gov/rest/pug | REST API |
| PUG-View | https://pubchem.ncbi.nlm.nih.gov/rest/pug_view | Formatted data |
| FTP Download | https://ftp.ncbi.nlm.nih.gov/pubchem/ | Bulk data |
| Entrez | E-utilities | NCBI integration |

## Data Categories

| Category | Description |
|----------|-------------|
| Identifiers | Names, synonyms, database IDs |
| Properties | Physicochemical, computed, experimental |
| Bioactivity | Assay results, targets, pathways |
| Safety | Hazards, GHS classification |
| Literature | Patents, publications |

## Limitations

- Depositor-submitted data quality varies
- Bioassay results not standardized across depositors
- Large-scale downloads require FTP; API has rate limits
- Some compound records may have incomplete data

## Related Resources

- [ChEBI](../chebi/_index.md) - Chemical ontology
- [ChEMBL](../../2.2.pharmaceuticals/chembl/_index.md) - Curated bioactivity
- [DrugBank](../../2.2.pharmaceuticals/drugbank/_index.md) - Drug information

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## References

1. Kim S, et al. (2023) "PubChem 2023 update." Nucleic Acids Res. 51(D1):D1373-D1380.
