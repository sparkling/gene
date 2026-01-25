---
id: biogrid
title: "BioGRID - Biological General Repository for Interaction Datasets"
type: source
parent: ../README.md
tier: 1
status: active
category: pathways.networks
subcategory: protein.protein.interactions
tags:
  - ppi
  - interactions
  - genetic
  - physical
  - open-access
---

# BioGRID - Biological General Repository for Interaction Datasets

**Category:** [Pathways & Networks](../../_index.md) > [Protein-Protein Interactions](../_index.md)

## Overview

BioGRID (Biological General Repository for Interaction Datasets) is a public database that archives and disseminates genetic and protein interaction data from model organisms and humans. It contains over 2 million curated interactions extracted from both high-throughput datasets and individual focused studies from the scientific literature.

Unlike some PPI databases that focus solely on physical interactions, BioGRID captures both physical protein-protein interactions AND genetic interactions, making it valuable for understanding functional relationships. The database is manually curated from primary literature by expert biocurators who extract interaction evidence with detailed experimental annotations.

BioGRID is freely available under an MIT license, making it one of the most permissive PPI resources for both academic and commercial applications.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Interactions | 2,200,000+ |
| Physical Interactions | 1,600,000+ |
| Genetic Interactions | 600,000+ |
| Unique Proteins | 85,000+ |
| Organisms | 70+ |
| Publications | 85,000+ |
| Human Interactions | ~700,000 |

## Primary Use Cases

1. **Protein interaction network construction** - Build PPI networks for analysis
2. **Functional annotation** - Infer gene function via interaction partners
3. **Drug target discovery** - Identify protein complexes and interaction hubs
4. **Genetic interaction analysis** - Study synthetic lethality and epistasis
5. **Cross-species comparison** - Compare interaction networks across organisms

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| BioGRID ID | Numeric | 112315 |
| Interaction ID | Numeric | 103 |
| Entrez Gene ID | Numeric | 7157 |
| UniProt | `[A-Z][0-9]{5}` | P04637 |
| Systematic Name | Organism-specific | YDL225W |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://thebiogrid.org | Browse/search |
| REST API | https://webservice.thebiogrid.org | JSON/Tab queries |
| Downloads | https://downloads.thebiogrid.org | Bulk files |
| Cytoscape App | BioGRID plugin | Network visualization |

### API Examples

```bash
# Search interactions for gene (requires access key)
curl "https://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=TP53&taxId=9606&accesskey=YOUR_KEY&format=json"

# Get interactions by BioGRID ID
curl "https://webservice.thebiogrid.org/interactions/?interactorList=112315&accesskey=YOUR_KEY&format=tab2"

# Search by PubMed ID
curl "https://webservice.thebiogrid.org/interactions/?pubmedList=1535557&accesskey=YOUR_KEY"

# Multi-gene query
curl "https://webservice.thebiogrid.org/interactions/?geneList=TP53|MDM2|BRCA1&taxId=9606&accesskey=YOUR_KEY"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| TAB 2.0 | Tab-delimited (26 columns) | Standard exchange |
| TAB 3.0 | Extended (35 columns) | Full annotation |
| PSI-MI TAB 2.5 | PSI-MI standard | Integration |
| JSON | Structured | Programmatic access |
| MITAB | Minimal tab | Simple networks |

### TAB 2.0 Key Columns

| Column | Description |
|--------|-------------|
| BioGRID Interaction ID | Unique interaction identifier |
| Entrez Gene ID A/B | Gene identifiers |
| BioGRID ID A/B | BioGRID gene identifiers |
| Systematic Name A/B | Systematic gene names |
| Official Symbol A/B | HGNC symbols |
| Experimental System | Detection method |
| Experimental System Type | physical/genetic |
| Throughput | high/low throughput |
| Score | Confidence score (if available) |
| Modification | PTM information |
| Phenotypes | Phenotype annotations |
| Qualifications | Additional details |
| Publication Source | PubMed ID |

## Experimental Systems

### Physical Interactions

| System | Description |
|--------|-------------|
| Two-hybrid | Yeast two-hybrid |
| Affinity Capture-MS | Mass spec pull-down |
| Affinity Capture-Western | Western blot pull-down |
| Co-fractionation | Co-purification |
| Co-crystal Structure | Crystal structure |
| Reconstituted Complex | In vitro reconstitution |
| PCA | Protein complementation |
| FRET | Fluorescence resonance |
| Proximity Label-MS | BioID, APEX |

### Genetic Interactions

| System | Description |
|--------|-------------|
| Synthetic Lethality | Lethal combination |
| Synthetic Growth Defect | Growth impairment |
| Synthetic Rescue | Suppressor interaction |
| Dosage Lethality | Gene dosage effect |
| Dosage Rescue | Dosage suppression |
| Phenotypic Enhancement | Phenotype worsening |
| Phenotypic Suppression | Phenotype improvement |

## License

| Aspect | Value |
|--------|-------|
| License | MIT License |
| Commercial Use | Yes |
| Attribution | Required |
| Redistribution | Allowed |
| Modification | Allowed |

## Cross-References

| Database | Relationship |
|----------|--------------|
| Entrez Gene | Primary gene IDs |
| UniProt | Protein accessions |
| RefSeq | Sequence IDs |
| Ensembl | Gene/protein IDs |
| PubMed | Literature references |
| SGD | Yeast systematic names |
| FlyBase | Drosophila IDs |
| WormBase | C. elegans IDs |

## Limitations

- API access requires registration and access key
- High-throughput studies may include false positives
- Genetic interactions may not reflect physical binding
- Coverage varies significantly across organisms

## Organisms (Selected)

| Organism | Tax ID | Interactions |
|----------|--------|--------------|
| Homo sapiens | 9606 | ~700,000 |
| Saccharomyces cerevisiae | 4932 | ~500,000 |
| Drosophila melanogaster | 7227 | ~100,000 |
| Mus musculus | 10090 | ~50,000 |
| Caenorhabditis elegans | 6239 | ~30,000 |
| Arabidopsis thaliana | 3702 | ~20,000 |
| Schizosaccharomyces pombe | 4896 | ~50,000 |

## See Also

- [IntAct](../intact/_index.md) - IMEx consortium member
- [STRING](../string/_index.md) - Functional associations
- [Pathway Commons](../../4.2.signaling.pathways/{pathwaycommons}/_index.md) - Pathway aggregator
