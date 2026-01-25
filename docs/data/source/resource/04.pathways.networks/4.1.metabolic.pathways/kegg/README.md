---
id: kegg
title: "KEGG - Kyoto Encyclopedia of Genes and Genomes"
type: source
parent: ../README.md
tier: 1
status: active
category: pathways.networks
subcategory: metabolic.pathways
tags:
  - pathways
  - metabolic
  - kegg
  - kgml
  - licensed
  - reference-genome
---

# KEGG - Kyoto Encyclopedia of Genes and Genomes

**Category:** [Pathways & Networks](../../_index.md) > [Metabolic Pathways](../_index.md)

## Overview

KEGG (Kyoto Encyclopedia of Genes and Genomes) is a comprehensive database resource for understanding high-level functions and utilities of biological systems from molecular-level information. It integrates genomic, chemical, and systemic functional information to provide pathway maps representing molecular interaction networks.

KEGG pathways are manually curated and cover metabolic pathways, genetic information processing, environmental information processing, cellular processes, organismal systems, human diseases, and drug development. The database uses its own KGML (KEGG Markup Language) XML format for pathway representation.

**Important:** KEGG requires licensing for commercial use. Academic users have limited free API access. Consider open alternatives like Reactome or WikiPathways for unrestricted use.

## Key Statistics

| Metric | Value |
|--------|-------|
| Human Pathways | ~400 |
| Total Pathways (all organisms) | 10,000+ |
| Organisms | 8,000+ |
| Genes | 30+ million |
| Compounds | 18,000+ |
| Reactions | 11,000+ |
| KO (Orthology) Groups | 25,000+ |

## Primary Use Cases

1. **Metabolic pathway mapping** - Map genes/proteins to metabolic pathways
2. **Pathway enrichment analysis** - Identify enriched pathways in gene sets
3. **Drug target identification** - Link compounds to pathway targets
4. **Cross-species comparison** - Compare pathways across organisms via KO

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Pathway ID | `path:{org}{#####}` | path:hsa00010 |
| Gene ID | `{org}:{id}` | hsa:7157 |
| Compound ID | `C{#####}` | C00001 |
| KO Number | `K{#####}` | K00001 |
| Reaction ID | `R{#####}` | R00001 |
| EC Number | `{#}.{#}.{#}.{#}` | 1.1.1.1 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| REST API | https://rest.kegg.jp/ | Limited free access |
| Web Interface | https://www.kegg.jp/kegg/pathway.html | Browse pathways |
| FTP | ftp://ftp.genome.jp/pub/kegg/ | License required |

### API Examples

```bash
# List human pathways
curl "https://rest.kegg.jp/list/pathway/hsa"

# Get pathway KGML
curl "https://rest.kegg.jp/get/hsa04115/kgml"

# Get pathway image
curl "https://rest.kegg.jp/get/hsa00010/image"

# Find pathways by gene
curl "https://rest.kegg.jp/link/pathway/hsa:7157"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| KGML | KEGG Markup Language (XML) | Programmatic pathway access |
| Flat files | Tab-delimited | Bulk data processing |
| PNG/SVG | Pathway images | Visualization |

## License

| Aspect | Value |
|--------|-------|
| License | Academic/Commercial (Requires License) |
| Academic Use | Limited free API access |
| Commercial Use | Paid license from Kanehisa Laboratories |
| Attribution | Required |

## Limitations

- Requires commercial license for commercial use and bulk downloads
- API rate limits may restrict high-throughput queries
- Pathway updates may lag behind primary literature
- Organism-specific identifiers require cross-referencing

## Open Alternatives

Due to licensing restrictions, consider these open alternatives:

| Database | License | Coverage |
|----------|---------|----------|
| Reactome | CC BY 4.0 | 2,712 human pathways |
| WikiPathways | CC0 | 3,100+ pathways |
| Pathway Commons | Open | Aggregated pathways |

## Cross-References

| Database | Relationship |
|----------|--------------|
| UniProt | Protein identifiers |
| ChEBI | Compound cross-references |
| GO | Biological process mapping |
| Reactome | Pathway equivalents |
| EC | Enzyme classification |

## See Also

- [Schema Documentation](./schema.md)
- [Reactome](../reactome/_index.md) - Open alternative
- [WikiPathways](../wikipathways/_index.md) - Community-curated alternative
