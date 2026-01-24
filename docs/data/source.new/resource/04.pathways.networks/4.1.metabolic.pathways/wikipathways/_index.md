---
id: wikipathways
title: "WikiPathways"
type: data-source
category: pathways
subcategory: metabolic-pathways
parent: ../_index.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [pathways, metabolic, wikipathways, gpml, community, open-access, cc0]
---

# WikiPathways

**Category:** [Pathways & Networks](../../_index.md) > [Metabolic Pathways](../_index.md)

## Overview

WikiPathways is an open, collaborative platform for capturing and disseminating biological pathway knowledge. Built on MediaWiki technology, it allows researchers worldwide to contribute and curate pathway information. The database uses GPML (Graphical Pathway Markup Language) as its native format, storing both biological semantics and visual layout information.

WikiPathways is distinguished by its community-driven curation model and complete public domain dedication (CC0). Pathways cover metabolic, signaling, and disease-related processes across 48 organisms. The platform integrates with analysis tools like Cytoscape, PathVisio, and enrichment analysis services.

As a truly open resource with no licensing restrictions, WikiPathways is ideal for both academic and commercial applications where data freedom is essential.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Pathways | 3,100+ |
| Human Pathways | 955+ |
| Supported Organisms | 48 |
| Genes/Proteins Referenced | 200,000+ |
| Metabolites Referenced | 50,000+ |
| Total Interactions | 1,000,000+ |
| Active Contributors | 1,500+ |

## Primary Use Cases

1. **Community pathway curation** - Contribute and edit pathway knowledge
2. **Pathway enrichment analysis** - Identify enriched pathways in gene sets
3. **Disease pathway exploration** - Browse disease-specific pathways
4. **Education and teaching** - Use pathways for learning molecular biology
5. **Data integration** - Combine with other open resources

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| WikiPathways ID | `WP{####}` | WP254 |
| Revision ID | `WP{####}_r{#####}` | WP254_r129662 |
| GraphId | Alphanumeric | abc123 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.wikipathways.org | Browse/edit pathways |
| REST API | https://webservice.wikipathways.org | Programmatic access |
| Bulk Downloads | https://www.wikipathways.org/download | GPML, BioPAX, JSON |
| SPARQL Endpoint | https://sparql.wikipathways.org | RDF queries |

### API Examples

```bash
# Get pathway GPML
curl "https://webservice.wikipathways.org/getPathway?pwId=WP254&format=json"

# Search pathways by text
curl "https://webservice.wikipathways.org/findPathwaysByText?query=apoptosis&species=Homo%20sapiens&format=json"

# List pathways by organism
curl "https://webservice.wikipathways.org/listPathways?organism=Homo%20sapiens&format=json"

# Get pathway info
curl "https://webservice.wikipathways.org/getPathwayInfo?pwId=WP254&format=json"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| GPML | Graphical Pathway Markup Language | Native format with layout |
| BioPAX Level 3 | OWL pathway exchange | Integration/reasoning |
| JSON | Structured data | Programmatic access |
| SVG/PNG | Images | Visualization |
| RDF | Semantic web | SPARQL queries |

## License

| Aspect | Value |
|--------|-------|
| License | CC0 1.0 Universal (Public Domain) |
| Commercial Use | Yes, unrestricted |
| Attribution | Not legally required (appreciated) |
| Redistribution | Unlimited |
| Derivative Works | Unlimited |

## Cross-References

| Database | Relationship |
|----------|--------------|
| Ensembl | Gene identifiers |
| Entrez Gene | Gene identifiers |
| UniProt | Protein identifiers |
| ChEBI | Metabolite identifiers |
| HMDB | Metabolite identifiers |
| KEGG | Compound/pathway mappings |
| Reactome | Pathway cross-references |

## GPML DataNode Types

| Type | Description | Typical Databases |
|------|-------------|-------------------|
| GeneProduct | Genes and proteins | Ensembl, Entrez, UniProt |
| Metabolite | Small molecules | ChEBI, HMDB, KEGG |
| Protein | Specific proteins | UniProt |
| Rna | RNA molecules | Ensembl, miRBase |
| Complex | Protein complexes | Complex Portal |
| Pathway | Linked pathways | WikiPathways, Reactome |

## Interaction Arrow Types (MIM Notation)

| ArrowHead | Meaning |
|-----------|---------|
| mim-stimulation | Stimulation/Activation |
| mim-inhibition | Inhibition |
| mim-catalysis | Catalysis |
| mim-conversion | Conversion |
| mim-binding | Non-covalent binding |
| mim-modification | Covalent modification |

## See Also

- [Schema Documentation](./schema.md)
- [Reactome](../reactome/_index.md) - Expert-curated pathways
- [KEGG](../kegg/_index.md) - Reference pathways (licensed)
