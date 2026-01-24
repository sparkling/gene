---
id: guides-integration
title: "Integration Guides"
type: directory
parent: ../_index.md
description: "Cross-database integration patterns and API usage guides"
last_updated: 2026-01-23
status: active
children:
  - compound-pathway-linking.md
  - pathway-target-mapping.md
  - pathway-downloads.md
tags: [integration, pathways, targets, cross-references, apis, guides]
---

**Parent:** [Guides](../_index.md)

# Integration Guides

Practical guides for linking compounds to genes, pathways, and targets across multiple databases.

## Overview

These guides provide:
- Identifier mapping strategies between databases
- API usage examples for major pathway databases
- Code samples for data integration pipelines
- Cross-reference resolution patterns

## Guides in This Section

| Guide | Description | Databases Covered |
|-------|-------------|-------------------|
| [compound-pathway-linking.md](./compound-pathway-linking.md) | Linking compounds to genes and biochemical pathways | PubChem, ChEMBL, DrugBank, KEGG, Reactome |
| [pathway-target-mapping.md](./pathway-target-mapping.md) | Pathway and target database data models | Reactome, KEGG, WikiPathways, UniProt, STRING |
| [pathway-downloads.md](./pathway-downloads.md) | Bulk download methods for pathway databases | Reactome, KEGG, WikiPathways, UniProt, STRING, GO |

## Identifier Mapping Strategy

### Compound Hub (PubChem CID)

PubChem CID serves as the primary compound identifier, linking to:
- ChEMBL (~2M compounds)
- DrugBank (~15K drugs)
- KEGG COMPOUND (~18K compounds)
- ChEBI (~150K entities)
- HMDB (~220K metabolites)

### Protein Hub (UniProt)

UniProt Accession links to:
- Ensembl (ENSP proteins, ENSG genes)
- Entrez Gene
- RefSeq
- PDB structures
- Reactome pathways

### Gene Hub (HGNC/Ensembl/Entrez)

| Hub | Best For | Cross-references |
|-----|----------|------------------|
| HGNC ID | Human gene nomenclature | All major databases |
| Ensembl Gene ID | Genomic context, variants | UniProt, RefSeq, Entrez |
| Entrez Gene ID | NCBI ecosystem, PubMed | UniProt, Ensembl, KEGG |

## Integration Flow

```
Compound (PubChem CID)
    |
    +-> ChEMBL ID -> Bioactivities -> Targets (UniProt)
    |
    +-> KEGG Compound -> Pathways -> Genes
    |
    +-> DrugBank ID -> Drug-Target Interactions

Target (UniProt)
    |
    +-> Ensembl Gene -> Genomic Data
    |
    +-> Reactome -> Pathway Membership
    |
    +-> STRING -> Protein Interactions
```

## Key APIs

| Database | Base URL | Format |
|----------|----------|--------|
| PubChem | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/` | JSON |
| ChEMBL | `https://www.ebi.ac.uk/chembl/api/data/` | JSON |
| UniProt | `https://rest.uniprot.org/` | JSON/TSV |
| Reactome | `https://reactome.org/ContentService/` | JSON |
| KEGG | `https://rest.kegg.jp/` | Text |
| STRING | `https://string-db.org/api/` | JSON/TSV |

## Related Resources

- [Wikidata Guides](../wikidata/_index.md) - SPARQL queries for knowledge graph integration
- [Resource Documentation](../../resource/_index.md) - Individual database schemas

---

*Last Updated: January 2026*
