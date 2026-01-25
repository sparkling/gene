---
id: string
title: "STRING - Search Tool for Retrieval of Interacting Genes/Proteins"
type: source
parent: ../README.md
tier: 1
status: active
category: pathways.networks
subcategory: protein.protein.interactions
tags:
  - ppi
  - interactions
  - functional-associations
  - network
  - open-access
---

# STRING - Search Tool for Retrieval of Interacting Genes/Proteins

**Category:** [Pathways & Networks](../../_index.md) > [Protein-Protein Interactions](../_index.md)

## Overview

STRING (Search Tool for the Retrieval of Interacting Genes/Proteins) is a comprehensive database of known and predicted protein-protein interactions. Unlike curated databases that focus solely on experimentally verified physical interactions, STRING integrates multiple evidence channels including experimental data, text mining, genomic context (gene neighborhood, fusion, co-occurrence), and co-expression to compute functional association scores.

STRING covers 14,094 organisms with 59.3 million proteins and 12.5 billion interactions. Each interaction is assigned a combined confidence score (0-1000) computed from individual evidence channel scores using a probabilistic framework. The database provides network visualization, functional enrichment analysis, and homology mapping tools.

STRING is widely used for network-based functional analysis and is freely available under CC BY 4.0.

## Key Statistics

| Metric | Value |
|--------|-------|
| Organisms | 14,094 |
| Proteins | 59,309,604 |
| Total Interactions | 12,535,845,684 |
| Human Proteins | ~20,000 |
| Human Interactions | ~11.8 million |
| High-confidence (score >= 700) | 500,000+ |

## Primary Use Cases

1. **Functional association networks** - Explore protein functional relationships
2. **Pathway enrichment analysis** - Identify enriched GO/KEGG/Reactome terms
3. **Network visualization** - Generate publication-ready network images
4. **Homology-based transfer** - Map interactions across species
5. **Protein function prediction** - Infer function via network context

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| STRING ID | `{taxid}.{ENSP}` | 9606.ENSP00000269305 |
| Ensembl Protein | `ENSP{###########}` | ENSP00000269305 |
| Gene Symbol | Text | TP53 |
| UniProt | `[A-Z][0-9]{5}` | P04637 |

## Score Channels

STRING combines evidence from multiple sources:

| Score | Channel | Description |
|-------|---------|-------------|
| nscore | Neighborhood | Gene neighborhood in genome |
| fscore | Fusion | Gene fusion events |
| pscore | Phylogenetic | Co-occurrence across genomes |
| ascore | Co-expression | Expression correlation |
| escore | Experimental | Laboratory evidence |
| dscore | Database | Curated pathway databases |
| tscore | Text-mining | Literature co-mentions |
| score | Combined | Probabilistic integration |

### Score Confidence Levels

| Combined Score | Confidence |
|---------------|------------|
| 0.900 - 1.000 | Highest |
| 0.700 - 0.899 | High |
| 0.400 - 0.699 | Medium |
| 0.150 - 0.399 | Low |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://string-db.org | Interactive networks |
| REST API | https://string-db.org/api/ | JSON/TSV/XML |
| Downloads | https://string-db.org/cgi/download | Bulk data |
| Cytoscape App | stringApp | Network visualization |

### API Examples

```bash
# Get interaction network
curl "https://string-db.org/api/json/network?identifiers=TP53%0dBRCA1&species=9606&required_score=400"

# Get interaction partners
curl "https://string-db.org/api/json/interaction_partners?identifiers=TP53&species=9606&limit=10"

# Functional enrichment
curl "https://string-db.org/api/json/enrichment?identifiers=TP53%0dBRCA1%0dBRCA2%0dATM&species=9606"

# Map identifiers to STRING IDs
curl "https://string-db.org/api/json/get_string_ids?identifiers=p53&species=9606&limit=1"

# Get network image
curl "https://string-db.org/api/image/network?identifiers=TP53%0dMDM2%0dCDKN1A&species=9606" > network.png

# Homology mapping
curl "https://string-db.org/api/json/homology?identifiers=TP53&species=9606&species_b=10090"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| JSON | Structured | API responses |
| TSV | Tab-delimited | Bulk processing |
| TSV-no-header | Tab without header | Streaming |
| XML | Structured | Legacy integration |
| PSI-MI | Standard | PPI exchange |
| PNG/SVG | Images | Visualization |

## API Response Fields

### Network Response

| Field | Type | Description |
|-------|------|-------------|
| stringId_A | String | STRING ID protein A |
| stringId_B | String | STRING ID protein B |
| preferredName_A | String | Gene symbol A |
| preferredName_B | String | Gene symbol B |
| ncbiTaxonId | Integer | Taxonomy ID |
| score | Float | Combined score |
| nscore-tscore | Float | Individual channel scores |

### Enrichment Response

| Field | Type | Description |
|-------|------|-------------|
| category | String | GO, KEGG, Reactome, etc. |
| term | String | Term identifier |
| description | String | Term description |
| p_value | Float | Hypergeometric p-value |
| fdr | Float | Benjamini-Hochberg FDR |
| preferredNames | String | Genes matching term |

## Enrichment Categories

| Category | Source |
|----------|--------|
| Process | GO Biological Process |
| Function | GO Molecular Function |
| Component | GO Cellular Component |
| KEGG | KEGG Pathways |
| Reactome | Reactome pathways |
| WikiPathways | WikiPathways |
| Pfam | Protein families |
| InterPro | Protein domains |
| HPO | Human Phenotype Ontology |
| DisGeNET | Disease associations |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution | Required |

## Cross-References

| Database | Relationship |
|----------|--------------|
| Ensembl | Primary protein IDs |
| UniProt | Protein accessions |
| NCBI Gene | Gene identifiers |
| RefSeq | Sequence IDs |
| GO | Enrichment annotations |
| KEGG | Pathway enrichment |
| Reactome | Pathway enrichment |

## Network Types

| Type | Description |
|------|-------------|
| functional | All association types (default) |
| physical | Only physical binding interactions |

## Limitations

- Combined scores may overestimate confidence for text-mined interactions
- Physical interactions not distinguished from functional associations by default
- Very large bulk files require significant storage and processing
- Text-mining may introduce false positives from co-mentions

## Bulk Download Files

| File | Description |
|------|-------------|
| protein.links.v12.0.txt.gz | All interactions |
| protein.links.detailed.v12.0.txt.gz | With channel scores |
| protein.info.v12.0.txt.gz | Protein annotations |
| protein.aliases.v12.0.txt.gz | Name aliases |
| protein.sequences.v12.0.fa.gz | Protein sequences |

## See Also

- [Schema Documentation](./schema.md)
- [IntAct](../intact/_index.md) - Experimental interactions
- [BioGRID](../biogrid/_index.md) - Curated interactions
- [STITCH](../../4.4.drug.target.interactions/{stitch}/_index.md) - Protein-chemical interactions
