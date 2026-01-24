---
id: pathwaycommons
title: "Pathway Commons"
type: data-source
category: pathways
subcategory: signaling-pathways
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [pathways, signaling, aggregation, biopax, sif, open-access]
---

# Pathway Commons

**Category:** [Pathways & Networks](../../_index.md) > [Signaling Pathways](../_index.md)

## Overview

Pathway Commons is a network biology resource and web application that aggregates biological pathway and molecular interaction data from multiple public databases. It provides a unified, integrated view of pathway information using the BioPAX (Biological Pathway Exchange) standard format, enabling researchers to access data from diverse sources through a single interface.

The resource integrates data from 22+ pathway and interaction databases including Reactome, KEGG, WikiPathways, PANTHER, HumanCyc, and multiple protein interaction databases. Pathway Commons provides both rich BioPAX data and simplified interaction networks (SIF format) suitable for network analysis tools.

Pathway Commons is maintained by the Computational Biology Center at Memorial Sloan Kettering Cancer Center and is freely available for academic and commercial use.

## Key Statistics

| Metric | Value |
|--------|-------|
| Integrated Databases | 22+ |
| Pathways | 5,000+ |
| Interactions | 2.3+ million |
| Physical Entities | 1.4+ million |
| Species | 20+ |
| Publications | 41,000+ |

## Data Sources

| Source | Type | Content |
|--------|------|---------|
| Reactome | Pathway | Curated human pathways |
| KEGG | Pathway | Metabolic/signaling |
| WikiPathways | Pathway | Community-curated |
| PANTHER | Pathway | Gene family pathways |
| HumanCyc | Pathway | Human metabolic |
| PhosphoSitePlus | PTM | Phosphorylation sites |
| BioGRID | Interaction | PPI data |
| IntAct | Interaction | Molecular interactions |
| HPRD | Interaction | Human PPI |
| DIP | Interaction | PPI data |

## Primary Use Cases

1. **Unified pathway access** - Query multiple databases through one interface
2. **Network analysis** - Build interaction networks from integrated data
3. **Pathway enrichment** - Perform enrichment using comprehensive pathway sets
4. **Data integration** - Combine pathway data from diverse sources
5. **Systems biology** - Model biological systems with rich pathway context

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| URI | `http://pathwaycommons.org/pc12/{type}/{id}` | Complex123 |
| BioPAX ID | Source-specific | REACT:R-HSA-109581 |
| HGNC Symbol | Gene symbols | TP53 |
| UniProt | Protein accession | P04637 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.pathwaycommons.org | Browse/search |
| REST API | https://www.pathwaycommons.org/pc2/ | BioPAX/SIF queries |
| Downloads | https://www.pathwaycommons.org/archives/ | Bulk BioPAX/SIF |
| cPath2 Client | Java library | Programmatic access |

### API Examples

```bash
# Search pathways
curl "https://www.pathwaycommons.org/pc2/search?q=BRCA1&type=pathway"

# Get pathway BioPAX
curl "https://www.pathwaycommons.org/pc2/get?uri=http://identifiers.org/reactome/R-HSA-109581"

# Get neighborhood network
curl "https://www.pathwaycommons.org/pc2/graph?source=TP53&kind=neighborhood&format=SIF"

# Find paths between genes
curl "https://www.pathwaycommons.org/pc2/graph?source=TP53&target=MDM2&kind=pathsbetween&format=SIF"

# Top pathways query
curl "https://www.pathwaycommons.org/pc2/top_pathways?q=DNA+repair&datasource=reactome"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| BioPAX Level 3 | Rich pathway representation | Full pathway detail |
| SIF | Simple Interaction Format | Network analysis |
| GSEA GMT | Gene set format | Enrichment analysis |
| SBGN-ML | Graphical notation | Visualization |
| JSON-LD | Linked data | Semantic web |

### SIF Format

Simple tab-delimited format for network analysis:

```
TP53	INTERACTS_WITH	MDM2
BRCA1	IN_COMPLEX_WITH	BARD1
ATM	CONTROLS-PHOSPHORYLATION-OF	TP53
```

### SIF Interaction Types

| Type | Description |
|------|-------------|
| INTERACTS_WITH | Physical interaction |
| IN_COMPLEX_WITH | Complex membership |
| CONTROLS-STATE-CHANGE-OF | State control |
| CONTROLS-TRANSPORT-OF | Transport control |
| CONTROLS-PHOSPHORYLATION-OF | Phosphorylation |
| CONTROLS-EXPRESSION-OF | Expression control |
| CATALYSIS-PRECEDES | Enzyme cascade |
| NEIGHBOR_OF | Neighborhood relation |

## License

| Aspect | Value |
|--------|-------|
| License | Free for all uses |
| Commercial Use | Yes |
| Attribution | Encouraged (cite sources) |
| Data Sources | Inherit source licenses |

## Cross-References

| Database | Relationship |
|----------|--------------|
| UniProt | Primary protein IDs |
| HGNC | Gene symbols |
| ChEBI | Small molecules |
| GO | Biological process |
| PubMed | Literature references |
| Reactome | Pathway source |
| KEGG | Pathway source |

## Graph Query Types

| Query | Description |
|-------|-------------|
| neighborhood | All interactions for source genes |
| pathsbetween | Paths connecting two gene sets |
| pathsfromto | Directed paths from sources to targets |
| commonstream | Common regulators/targets |

## See Also

- [Schema Documentation](./schema.md) - Technical schema and data formats
- [Download Instructions](./download.md) - Bulk data acquisition guide
- [Reactome](../../4.1.metabolic.pathways/reactome/_index.md) - Primary data source
- [BioGRID](../../4.3.protein.protein.interactions/biogrid/_index.md) - Interaction source
- [IntAct](../../4.3.protein.protein.interactions/intact/_index.md) - Interaction source
