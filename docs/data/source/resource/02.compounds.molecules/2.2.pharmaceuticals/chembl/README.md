---
id: chembl
title: "ChEMBL - Bioactivity Database"
type: source
parent: ../README.md
tier: 1
status: active
category: compounds.molecules
subcategory: pharmaceuticals
tags:
  - bioactivity
  - drug-discovery
  - pharmacology
  - targets
  - assays
---

# ChEMBL - Bioactivity Database

## Overview

ChEMBL is a manually curated database of bioactive molecules with drug-like properties, maintained by the European Bioinformatics Institute (EMBL-EBI). It brings together chemical, bioactivity, and genomic data to aid the translation of genomic information into effective new drugs. ChEMBL is one of the most comprehensive resources for drug discovery and medicinal chemistry research.

The database contains bioactivity data extracted from the scientific literature, covering binding constants, pharmacology, and ADMET properties. Each compound is linked to its targets (proteins, cell lines, organisms) through standardized assay descriptions, enabling systematic analysis of structure-activity relationships across therapeutic areas.

ChEMBL provides extensive data on approved drugs, clinical candidates, and research compounds, making it essential for lead identification, target validation, and computational drug discovery workflows.

## Key Statistics

| Metric | Value |
|--------|-------|
| Compounds | 2.4+ million |
| Activities | 20+ million |
| Assays | 1.6+ million |
| Targets | 15,000+ |
| Documents | 90,000+ |
| Approved Drugs | 4,000+ |
| Last Update | ChEMBL 34 (2024) |

## Primary Use Cases

1. **Drug Discovery** - Identify active compounds against targets of interest
2. **Target Validation** - Assess druggability with existing compound data
3. **SAR Analysis** - Explore structure-activity relationships systematically
4. **ADMET Prediction** - Train models on curated property data
5. **Competitive Intelligence** - Track drug development landscape

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| ChEMBL ID (Compound) | CHEMBL + digits | CHEMBL25 |
| ChEMBL ID (Target) | CHEMBL + digits | CHEMBL2842 |
| ChEMBL ID (Assay) | CHEMBL + digits | CHEMBL1217643 |
| InChI Key | 27 characters | Standard format |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.ebi.ac.uk/chembl/ | Search and browse |
| REST API | https://www.ebi.ac.uk/chembl/api/data | Programmatic access |
| FTP Download | https://ftp.ebi.ac.uk/pub/databases/chembl/ | Bulk data |
| PostgreSQL | Available | Database dumps |
| ChEMBL-Python | pip install chembl_webresource_client | Python client |

## Data Categories

| Category | Description |
|----------|-------------|
| Compound Properties | Structure, physicochemical, drug-likeness |
| Bioactivities | IC50, Ki, EC50, potency measurements |
| Assays | Experimental protocols and conditions |
| Targets | Proteins, cell lines, organisms |
| Drug Mechanisms | Approved drug target interactions |

## Limitations

- Activity values extracted from literature may have measurement variability
- Assay conditions vary across sources
- Some historical data may use outdated target nomenclature
- Non-human target data may not translate to human biology

## Related Resources

- [DrugBank](../drugbank/README.md) - Comprehensive drug data
- [BindingDB](../../2.7.compound.target.interactions/bindingdb/README.md) - Binding affinities
- [PubChem](../../2.6.chemical.ontology.classification/pubchem/README.md) - Chemical repository

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## References

1. Zdrazil B, et al. (2024) "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods." Nucleic Acids Res. 52(D1):D1180-D1192. DOI: 10.1093/nar/gkad1004
