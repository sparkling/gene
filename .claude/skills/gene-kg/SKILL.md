---
name: gene-kg
description: |
  Query biomedical data sources naturally. Find databases for genetics, drugs, diseases, proteins,
  pathways, traditional medicine (TCM, Ayurveda), nutrition, microbiome, literature. Get APIs,
  licenses, formats, dataset sizes. Compare sources. Discover connections. 145 curated sources.

  Examples: "10 largest databases" "Find TCM sources" "production-ready genetics databases"
  "compare DrugBank and ChEMBL" "what connects to ClinVar" "databases with REST APIs"
allowed-tools: Read, Bash
---

# Gene Knowledge Graph

Query 145 biomedical data sources across 9 categories.

## How to Use

Run the query script with natural language:
```bash
python3 scripts/kg-query.py query "<natural language question>"
```

## Response Guidelines

When presenting results to the user:
- Use **numbered lists** for multiple sources
- **Bold** the source name, add maintainer in parentheses
- Explain tier levels: 1=Production-ready, 2=Beta, 3=Experimental
- Include API URL and license if available
- **Never expose**: URIs, SPARQL queries, internal IDs, prefixes
- Add helpful context (e.g., "Tier 1 means the API is stable")

## Example Questions

| Question | What it finds |
|----------|---------------|
| "10 largest databases" | Sources by dataset size |
| "databases for TCM" | Traditional Chinese Medicine sources |
| "production-ready genetics databases" | Tier 1 genetics sources |
| "compare DrugBank and ChEMBL" | Side-by-side comparison |
| "what connects to ClinVar" | Related databases |
| "databases with REST APIs" | Sources with programmatic access |

## Categories (9 total, 145 sources)

| Category | Sources | Examples |
|----------|---------|----------|
| Genetics & Genomics | 32 | ClinVar, gnomAD, dbSNP |
| Compounds & Molecules | 31 | DrugBank, ChEMBL, PubChem |
| Diseases & Phenotypes | 32 | OMIM, HPO, DisGeNET |
| Pathways & Networks | 19 | KEGG, Reactome, WikiPathways |
| Traditional Medicine | 18 | TCMSP, HerbMed, CMAUP |
| Nutrition & Food | 13 | FoodDB, USDA, Phenol-Explorer |
| Proteins | 9 | UniProt, PDB, InterPro |
| Literature | 17 | PubMed, Europe PMC |
| Microbiome | 14 | HMDB, gutMGene |

## Tier System

| Tier | Meaning | Count |
|------|---------|-------|
| 1 | Production-ready (stable API, maintained) | ~40 |
| 2 | Beta (usable but may change) | ~60 |
| 3 | Experimental (limited access) | ~45 |

## Behind the Scenes

This skill queries a knowledge graph with 12,098 triples using:
- Semantic search (RuVector embeddings)
- Structured queries (SPARQL on Fuseki)

You don't need to know SPARQL - just ask questions naturally.
