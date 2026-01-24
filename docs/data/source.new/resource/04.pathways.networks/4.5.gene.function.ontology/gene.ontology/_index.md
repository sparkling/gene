---
id: gene-ontology
title: "Gene Ontology (GO)"
type: data-source
category: pathways
subcategory: gene-function-ontology
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [ontology, gene-function, go-terms, annotation, open-access]
---

# Gene Ontology (GO)

**Category:** [Pathways & Networks](../../_index.md) > [Gene Function & Ontology](../_index.md)

## Overview

The Gene Ontology (GO) is the most widely used structured vocabulary for describing gene product functions across all species. It provides a computational representation of current scientific knowledge about gene and protein function, enabling consistent annotation and cross-species comparison.

GO comprises three ontologies: Molecular Function (what a gene product does), Biological Process (the larger biological goal), and Cellular Component (where in the cell it acts). The GO Consortium maintains both the ontology (~47,000 terms) and standardized gene-to-term annotations for thousands of species. GO-CAM (Causal Activity Models) extends traditional annotations by capturing how gene products work together in pathways.

GO is freely available under CC BY 4.0 and is the foundation for functional enrichment analysis in genomics research.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total GO Terms | ~47,000 |
| Molecular Function Terms | ~12,000 |
| Biological Process Terms | ~30,000 |
| Cellular Component Terms | ~4,500 |
| Total Annotations | 500,000,000+ |
| Human Annotations | ~800,000 |
| Species with Annotations | 4,000+ |
| GO-CAM Models | 1,500+ |

## Three Ontologies

| Ontology | Aspect Code | Root Term | Description |
|----------|-------------|-----------|-------------|
| Molecular Function | F | GO:0003674 | Molecular-level activities |
| Biological Process | P | GO:0008150 | Multi-step biological objectives |
| Cellular Component | C | GO:0005575 | Cellular locations |

## Primary Use Cases

1. **Functional enrichment analysis** - Identify over-represented functions in gene lists
2. **Gene annotation** - Describe gene/protein function systematically
3. **Cross-species comparison** - Compare functions across organisms
4. **Knowledge representation** - Formal biological knowledge modeling
5. **Data integration** - Common vocabulary for diverse databases

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| GO Term | `GO:{#######}` | GO:0008152 |
| GO-CAM Model | `gomodel:{hex}` | gomodel:5fce9b7300002411 |
| Evidence Code | 2-3 letters | IDA, IMP, IEA |
| RO Relation | `RO:{#######}` | RO:0002413 |

## Evidence Codes

### Experimental Evidence (High Weight)

| Code | Name | Description |
|------|------|-------------|
| EXP | Experiment | General experimental |
| IDA | Direct Assay | Direct biochemical assay |
| IPI | Physical Interaction | Physical interaction |
| IMP | Mutant Phenotype | Mutant analysis |
| IGI | Genetic Interaction | Genetic interaction |
| IEP | Expression Pattern | Expression evidence |

### Computational Evidence (Medium Weight)

| Code | Name | Description |
|------|------|-------------|
| ISS | Sequence Similarity | Sequence-based |
| ISO | Sequence Orthology | Ortholog evidence |
| IBA | Biological Ancestor | Phylogenetic inference |

### Electronic Evidence (Low Weight)

| Code | Name | Description |
|------|------|-------------|
| IEA | Electronic Annotation | Automated annotation |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| AmiGO Browser | https://amigo.geneontology.org | Browse/search |
| QuickGO | https://www.ebi.ac.uk/QuickGO | EBI interface |
| REST API | https://api.geneontology.org/api | JSON responses |
| SPARQL | https://sparql.geneontology.org/sparql | RDF queries |
| Downloads | https://current.geneontology.org/ | Bulk files |

### API Examples

```bash
# Get GO term details
curl "https://api.geneontology.org/api/ontology/term/GO:0008152"

# Search terms
curl "https://api.geneontology.org/api/search/entity?q=apoptosis&category=ontology_class"

# Get annotations for gene
curl "https://api.geneontology.org/api/bioentity/gene/HGNC:11998/function"

# AmiGO term page
# https://amigo.geneontology.org/amigo/term/GO:0006915
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| OBO | Text-based ontology | Standard use |
| OWL | Web Ontology Language | Reasoning/semantic web |
| GAF 2.2 | Gene Association File | Annotations |
| GO-CAM JSON | Causal activity models | Pathway-level |
| GPAD | Gene Product Association Data | Compact annotations |

### GAF 2.2 Columns

| Column | Field | Description |
|--------|-------|-------------|
| 1 | DB | Source database |
| 2 | DB_Object_ID | Identifier |
| 3 | DB_Object_Symbol | Gene symbol |
| 5 | GO_ID | GO term |
| 7 | Evidence_Code | Evidence type |
| 9 | Aspect | F, P, or C |
| 13 | Taxon | Species |
| 14 | Date | Annotation date |
| 15 | Assigned_By | Annotating group |

## GO Slims

Pre-defined subsets for high-level analysis:

| Slim | Terms | Use Case |
|------|-------|----------|
| goslim_generic | ~150 | General-purpose |
| goslim_chembl | ~100 | Drug discovery |
| goslim_metagenomics | ~150 | Metagenomics |
| goslim_plant | ~200 | Plant biology |
| goslim_yeast | ~180 | Yeast research |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution | Required |

## Cross-References

| Database | Relationship |
|----------|--------------|
| UniProt | Primary protein annotations |
| Ensembl | Gene annotations |
| NCBI Gene | Gene annotations |
| Reactome | Pathway cross-references |
| ChEBI | Chemical entity refs |
| EC | Enzyme classification |
| InterPro | Domain annotations |

## Download Files

| File | Format | Description |
|------|--------|-------------|
| go-basic.obo | OBO (35 MB) | Standard ontology |
| go.obo | OBO (45 MB) | Full with relationships |
| go.owl | OWL (150 MB) | Semantic web format |
| goa_human.gaf.gz | GAF | Human annotations |
| goa_uniprot_all.gaf.gz | GAF | All species |

## SPARQL Query Example

```sparql
PREFIX GO: <http://purl.obolibrary.org/obo/GO_>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT ?term ?label
WHERE {
  ?term rdfs:subClassOf* GO:0008150 .
  ?term rdfs:label ?label .
}
LIMIT 100
```

## See Also

- [Schema Documentation](./schema.md)
- [MSigDB](../msigdb/_index.md) - Gene set collections using GO
- [Reactome](../../4.1.metabolic.pathways/reactome/_index.md) - GO-annotated pathways
