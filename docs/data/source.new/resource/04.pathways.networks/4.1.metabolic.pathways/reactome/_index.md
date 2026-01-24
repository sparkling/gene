---
id: reactome
title: "Reactome Pathway Database"
type: data-source
category: pathways
subcategory: metabolic-pathways
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [pathways, metabolic, reactome, neo4j, biopax, open-access]
---

# Reactome Pathway Database

**Category:** [Pathways & Networks](../../_index.md) > [Metabolic Pathways](../_index.md)

## Overview

Reactome is an open-source, open-access, manually curated and peer-reviewed pathway database. It provides intuitive bioinformatics tools for visualization, interpretation, and analysis of pathway knowledge. Reactome represents biological pathways as an interconnected graph where molecular reactions convert input physical entities to output entities.

The database is maintained by EMBL-EBI, OICR, and CSHL, with quarterly releases containing expert-curated human pathways and computationally inferred ortholog pathways for 24 species. Reactome uses a Neo4j graph database backend and provides data in multiple formats including BioPAX, SBML, and custom JSON.

Reactome is one of the most comprehensive open-access pathway resources, making it an excellent alternative to licensed databases like KEGG for both academic and commercial use.

## Key Statistics

| Metric | Value |
|--------|-------|
| Human Pathways | 2,712 |
| Human Reactions | 13,872 |
| Human Proteins | 11,196 |
| Small Molecules | 1,925 |
| Species | 24 |
| Literature References | 35,000+ |

## Primary Use Cases

1. **Pathway enrichment analysis** - Identify over-represented pathways in gene/protein lists
2. **Reaction network visualization** - Explore molecular mechanisms graphically
3. **Cross-species pathway comparison** - Compare pathways across 24 organisms
4. **Drug target contextualization** - Place drug targets in biological context
5. **Systems biology modeling** - Export to SBML for mathematical modeling

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Stable ID | `R-{species}-{#######}` | R-HSA-1430728 |
| Database ID (dbId) | Integer | 1430728 |
| Versioned ID | `R-{species}-{#######}.{version}` | R-HSA-1430728.16 |
| DOI | `10.3180/R-{species}-{#######}.{#}` | 10.3180/R-HSA-1430728.15 |

### Species Codes

| Species | Code | Tax ID |
|---------|------|--------|
| Homo sapiens | HSA | 9606 |
| Mus musculus | MMU | 10090 |
| Rattus norvegicus | RNO | 10116 |
| Danio rerio | DRE | 7955 |
| Drosophila melanogaster | DME | 7227 |
| Saccharomyces cerevisiae | SCE | 4932 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| REST API | https://reactome.org/ContentService | JSON responses |
| Neo4j Browser | https://reactome.org/dev/graph-database | Graph queries |
| Web Interface | https://reactome.org | Browse/search |
| FTP Downloads | https://reactome.org/download/current/ | Bulk data |

### API Examples

```bash
# Get pathway by ID
curl "https://reactome.org/ContentService/data/query/R-HSA-109581" \
  -H "Accept: application/json"

# Get pathway events
curl "https://reactome.org/ContentService/data/pathway/R-HSA-109581/containedEvents"

# Map UniProt to pathways
curl "https://reactome.org/ContentService/data/mapping/UniProt/P04637/pathways"

# Pathway enrichment analysis
curl -X POST "https://reactome.org/AnalysisService/identifiers/" \
  -H "Content-Type: text/plain" \
  -d "P04637,P53350,Q00987"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| BioPAX Level 3 | OWL-based pathway exchange | Pathway integration |
| SBML Level 3 | Systems biology models | Mathematical simulation |
| SBGN-ML | Graphical notation | Visualization |
| PSI-MITAB | Protein interactions | Network analysis |
| JSON | API responses | Programmatic access |
| Neo4j | Graph database | Complex queries |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution | Required |
| Redistribution | Allowed with attribution |

## Cross-References

| Database | Relationship |
|----------|--------------|
| UniProt | Protein identifiers (primary) |
| ChEBI | Small molecule references |
| GO | Biological process/compartment |
| Ensembl | Gene identifiers |
| KEGG | Pathway mappings |
| PubMed | Literature citations |

## Neo4j Query Examples

```cypher
-- Get pathway with sub-events
MATCH (p:Pathway {stId: 'R-HSA-109581'})-[:hasEvent*]->(e:Event)
RETURN p.displayName AS Pathway, collect(DISTINCT e.displayName) AS Events

-- Get all proteins in a pathway
MATCH (p:Pathway {stId: 'R-HSA-109581'})-[:hasEvent*]->(r:ReactionLikeEvent)
MATCH (r)-[:input|output]->(pe:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
WHERE re.databaseName = 'UniProt'
RETURN DISTINCT re.identifier AS UniProtID, re.geneName AS GeneName
```

## See Also

- [Schema Documentation](./schema.md)
- [WikiPathways](../wikipathways/_index.md) - Community-curated pathways
- [KEGG](../kegg/_index.md) - Alternative (licensed)
