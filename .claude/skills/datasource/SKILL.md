---
name: datasource
description: "Data source, datasource, database queries - find TCM, genomics, compounds databases. Use for: 'what data sources', 'find databases', 'TCM sources', 'list sources', 'compare databases', 'which databases have'. Natural language queries over biomedical data catalog."
---

## EXECUTION INSTRUCTIONS (MANDATORY)

**STOP. Before doing anything else, follow these steps exactly:**

### Step 1: Generate SPARQL Query
1. Read `resources/ontology-context.md` (in this skill's directory) for prefixes and schema
2. Read `resources/query-guide.md` for query patterns matching user's intent
3. Generate a SPARQL query for the user's question

### Step 2: Execute Query Against Knowledge Graph
Execute SPARQL via the Fuseki endpoint:
```bash
curl -s "http://localhost:3030/gene/sparql" \
  --data-urlencode "query=[YOUR SPARQL]" \
  -H "Accept: application/sparql-results+json"
```

### Step 3: Format and Return Results
- Present clean, human-readable output
- Never show SPARQL queries, URIs, or prefixes to user
- Group by category, sort by tier

### FORBIDDEN ACTIONS
- **NEVER** read files in `docs/data/source/resource/` - that's raw documentation, not the KG
- **NEVER** manually browse markdown files to answer queries
- **ALWAYS** query the knowledge graph via qlever/SPARQL

---

# Data Source Knowledge Graph

## What This Skill Does

Enables natural language queries over the Gene biomedical data source catalog. Ask questions in plain English and get structured answers without needing to know SPARQL, RDF, or the underlying ontology.

**92 data sources** across **9 categories** and **44 subcategories** are available.

## Quick Start

Just ask your question naturally:

- "What TCM data sources are available?"
- "Find the 10 largest genomics databases"
- "Which sources have REST API access?"
- "Compare DrugBank and ChEMBL"
- "Show me open-license sources for cancer research"

---

## What You Can Ask

### Browse by Category

```
Show me all Traditional Chinese Medicine data sources
List genetics databases
What's in the Pathways category?
```

### Search and Filter

```
Find sources with SPARQL endpoints
Which databases have open licenses?
Show tier 1 primary sources
Find sources updated monthly
```

### Compare Sources

```
Compare UniProt and PDB
What's the difference between ClinVar and COSMIC?
```

### Get Details

```
Tell me about DrugBank
What formats does PubChem support?
How do I access the KEGG database?
```

### Statistics

```
How many data sources are there?
Count sources by category
Show tier distribution
```

---

## Categories Available

| # | Category | Description |
|---|----------|-------------|
| 01 | **Genetics and Genomics** | Variants, population genetics, pharmacogenomics, expression |
| 02 | **Compounds and Molecules** | Natural products, pharmaceuticals, drug metabolism |
| 03 | **Diseases and Phenotypes** | Disease ontologies, phenotypes, gene-disease associations |
| 04 | **Pathways and Networks** | Metabolic/signaling pathways, protein interactions |
| 05 | **Traditional Medicine** | TCM, Ayurveda, herbal medicine databases |
| 06 | **Nutrition and Food** | Food composition, supplements, metabolomics |
| 07 | **Proteins and Molecular Biology** | Protein sequences, structures, interactions |
| 08 | **Literature and Knowledge** | Scientific literature, knowledge bases, identifiers |
| 09 | **Microbiome** | Gut, body site, and host-microbe interactions |

---

## How It Works

This skill uses a **hybrid query architecture**:

1. **Semantic Understanding**: Your natural language query is analyzed to understand intent
2. **Knowledge Graph Query**: The LLM generates SPARQL queries using ontology definitions
3. **Vector Search**: Semantic search enriches results with related concepts
4. **Result Synthesis**: Results are formatted into clean, readable output

The underlying ontology (OWL), taxonomy (SKOS), and validation shapes (SHACL) provide the schema context for intelligent query generation.

---

## Output Examples

### Category Browse

```
Traditional Chinese Medicine Data Sources (5 sources)

Tier 1 (Primary):
  - TCMSP: TCM Systems Pharmacology Database
    Access: REST API, Direct Download | License: Open Access

  - TCMID: TCM Integrated Database
    Access: Web Interface | License: Academic

Tier 2 (Secondary):
  - HIT 2.0: Herb Ingredients' Targets Database
    Access: REST API | License: CC BY
```

### Source Details

```
DrugBank

Category: Compounds and Molecules > Pharmaceuticals
Tier: 1 (Primary)
Status: Active

Description: Comprehensive drug and drug target database combining
detailed drug data with comprehensive drug target information.

Access Methods:
  - REST API (authentication required)
  - Direct Download (XML, CSV)

License: CC BY-NC 4.0
  - Commercial use: Requires license
  - Attribution: Required

Website: https://go.drugbank.com/
```

---

## Advanced Usage

### Combining Filters

```
Find tier 1 sources in genomics with REST API and open license
```

### Cross-References

```
What sources link to UniProt?
Show cross-references for ClinVar
```

### Related Sources

```
Find sources related to DrugBank
What complements KEGG Pathway?
```

---

## Tips for Best Results

1. **Be specific**: "TCM herb databases" works better than "Chinese medicine"
2. **Use category names**: Reference the 9 main categories when browsing
3. **Combine criteria**: "tier 1 + REST API + open license" narrows results
4. **Ask for comparisons**: Great for choosing between similar sources

---

## Troubleshooting

### "No results found"

- Try broader search terms
- Check category spelling
- Use semantic search: "find sources similar to..."

### Too many results

- Add filters: tier, license, access method
- Specify subcategory
- Ask for "top 5" or "most relevant"

---

## Technical Notes

This skill orchestrates:
- **qlever** for high-performance SPARQL execution
- **sparql** skill for query optimization
- **ruvector-search** for semantic enrichment
- **owl/skos/shacl** skills for ontology-aware query generation

The LLM dynamically generates SPARQL queries using the ontology definitions as context, enabling flexible natural language understanding without hardcoded templates.

---

## Internal Resources (For LLM Query Generation)

When generating SPARQL queries, load these resources as context:

- **Ontology Context**: [resources/ontology-context.md](resources/ontology-context.md) - Prefixes, classes, properties, taxonomy concepts
- **Query Guide**: [resources/query-guide.md](resources/query-guide.md) - Intent classification and few-shot examples
- **Implementation**: [docs/IMPLEMENTATION.md](docs/IMPLEMENTATION.md) - Architecture and orchestration details

## Project Resources

- Ontology files: `docs/data/source/ontology/` (OWL, SKOS, SHACL)
- Taxonomy docs: `docs/data/source/taxonomy/`
- Source documentation: `docs/data/source/resource/`

## Related Skills

- [ruvector-search](../ruvector-search/) - Semantic vector search
- [sparql](link) - SPARQL query optimization
- [qlever](link) - High-performance SPARQL engine
- [owl](link) - OWL ontology understanding
- [skos](link) - SKOS taxonomy navigation
- [shacl](link) - SHACL constraint validation
