---
id: intact
title: "IntAct Molecular Interaction Database"
type: source
parent: ../README.md
tier: 1
status: active
category: pathways.networks
subcategory: protein.protein.interactions
tags:
  - ppi
  - interactions
  - imex
  - psi-mi
  - embl-ebi
  - open-access
---

# IntAct Molecular Interaction Database

**Category:** [Pathways & Networks](../../README.md) > [Protein-Protein Interactions](../README.md)

## Overview

IntAct is a freely available, open-source molecular interaction database maintained by EMBL-EBI. It is a founding member of the IMEx (International Molecular Exchange) consortium, which provides standardized curation of interaction data across multiple databases. IntAct focuses on high-quality manual curation with deep annotation using PSI-MI controlled vocabularies.

The database captures molecular interactions from literature with detailed experimental evidence, including interaction detection methods, participant identification methods, and confidence scoring. IntAct's MI-score provides a standardized confidence measure based on publication count, method diversity, and annotation depth.

As part of the IMEx consortium, IntAct data follows the MIMIx (Minimum Information about Molecular Interaction experiments) guidelines and is available in standard PSI-MI formats.

## Key Statistics

| Metric | Value |
|--------|-------|
| Binary Interactions | 1,233,546 |
| Interactor Proteins | 147,892 |
| Publications | 122,341 |
| Organisms | 623 |
| Human Interactions | ~450,000 |
| Human Proteins | ~19,500 |

## Primary Use Cases

1. **High-confidence PPI networks** - Build networks with experimental evidence
2. **IMEx federated queries** - Access multiple databases via PSICQUIC
3. **Interaction detection methods** - Filter by experimental technique
4. **Complex annotation** - Study protein complex composition
5. **Standards-compliant integration** - Use PSI-MI standard formats

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| IntAct Interaction | `EBI-{#######}` | EBI-77734 |
| IntAct Interactor | `EBI-{#######}` | EBI-366083 |
| IMEx ID | `IM-{#####}` | IM-12345 |
| UniProt | `[A-Z][0-9]{5}` | P04637 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.ebi.ac.uk/intact/ | Browse/search |
| PSICQUIC | https://www.ebi.ac.uk/Tools/webservices/psicquic | REST API |
| FTP | ftp://ftp.ebi.ac.uk/pub/databases/intact/current/ | Bulk downloads |
| Complex Portal | https://www.ebi.ac.uk/complexportal/ | Curated complexes |

### PSICQUIC API Examples

```bash
# Search by protein identifier
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/query/identifier:P04637?format=tab25"

# Search by gene name
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/query/alias:TP53?format=tab27"

# Human-human interactions
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/query/taxidA:9606%20AND%20taxidB:9606?format=tab25"

# By detection method (two-hybrid)
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/query/detmethod:%22MI:0018%22?format=tab25"

# By publication
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/query/pubid:1535557?format=tab25"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| PSI-MITAB 2.5 | 15 columns | Basic interaction data |
| PSI-MITAB 2.7 | 42 columns | Full annotation |
| PSI-MI XML 2.5/3.0 | Rich XML | Complete detail |
| JSON | Structured | Programmatic access |

### PSI-MITAB 2.5 Columns

| Column | Description | Example |
|--------|-------------|---------|
| 1-2 | Interactor IDs | uniprotkb:P04637 |
| 3-4 | Alternative IDs | intact:EBI-366083 |
| 5-6 | Aliases | uniprotkb:TP53(gene name) |
| 7 | Detection method | psi-mi:"MI:0018"(two hybrid) |
| 8 | First author | Momand J (1992) |
| 9 | Publication ID | pubmed:1535557 |
| 10-11 | Taxonomy IDs | taxid:9606(human) |
| 12 | Interaction type | psi-mi:"MI:0915"(physical association) |
| 13 | Source database | psi-mi:"MI:0469"(IntAct) |
| 14 | Interaction ID | intact:EBI-77734 |
| 15 | Confidence score | intact-miscore:0.73 |

## MI-Score Confidence

| Score Range | Confidence | Description |
|-------------|------------|-------------|
| 0.75 - 1.00 | High | Multiple evidence types |
| 0.50 - 0.74 | Medium | Good evidence |
| 0.25 - 0.49 | Low | Limited evidence |
| 0.00 - 0.24 | Very Low | Sparse annotation |

## PSI-MI Controlled Vocabulary

### Detection Methods (Selected)

| MI ID | Term |
|-------|------|
| MI:0018 | Two hybrid |
| MI:0019 | Coimmunoprecipitation |
| MI:0096 | Pull down |
| MI:0114 | X-ray crystallography |
| MI:0676 | Tandem affinity purification |

### Interaction Types

| MI ID | Term |
|-------|------|
| MI:0407 | Direct interaction |
| MI:0915 | Physical association |
| MI:0403 | Colocalization |
| MI:0195 | Covalent binding |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution | Required |
| Software License | Apache 2.0 |

## Limitations

- Literature curation may lag behind recent publications
- PSI-MI format requires specialized parsing
- MI-scores depend on available annotation depth
- Coverage biased toward well-studied proteins

## IMEx Consortium Partners

| Database | MI ID |
|----------|-------|
| IntAct | MI:0469 |
| MINT | MI:0471 |
| DIP | MI:0465 |
| MatrixDB | MI:0917 |
| InnateDB | MI:0974 |

## Cross-References

| Database | Relationship |
|----------|--------------|
| UniProt | Primary protein IDs |
| Ensembl | Gene/protein IDs |
| RefSeq | Sequence IDs |
| ChEBI | Small molecule IDs |
| GO | Functional annotations |
| PubMed | Literature citations |
| Reactome | Pathway context |

## See Also

- [Schema Documentation](./schema.md)
- [STRING](../string/README.md) - Functional associations
- [BioGRID](../biogrid/README.md) - Genetic/physical interactions
