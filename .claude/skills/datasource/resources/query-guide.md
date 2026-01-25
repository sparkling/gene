# SPARQL Query Generation Guide

This guide helps generate SPARQL queries from natural language questions about data sources.

## Query Intent Classification

First, classify the user's intent:

| Intent | Keywords | Query Pattern |
|--------|----------|---------------|
| **browse** | "list", "show", "what's in" | SELECT with category filter |
| **search** | "find", "search", "which" | SELECT with text/property filter |
| **details** | "tell me about", "describe" | SELECT all properties for one source |
| **compare** | "compare", "difference", "vs" | SELECT properties for multiple sources |
| **stats** | "count", "how many", "statistics" | SELECT with COUNT/GROUP BY |
| **related** | "related", "similar", "connected" | SELECT via relationship properties |

## Natural Language → SPARQL Examples

### Example 1: Category Browse

**User**: "Show me all TCM data sources"

**Intent**: browse (Traditional Chinese Medicine category)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX tax: <https://gene.ai/taxonomy/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?id ?title ?tierNum WHERE {
  ?source a ds:DataSource ;
          ds:id ?id ;
          ds:title ?title ;
          ds:belongsToSubcategory ?subcat ;
          ds:tier ?tier .
  ?subcat ds:belongsToCategory tax:TraditionalMedicine .
  ?tier skos:notation ?tierNum .
}
ORDER BY ?tierNum ?title
```

### Example 2: Multi-Filter Search

**User**: "Find tier 1 sources with REST API"

**Intent**: search (tier filter + access method filter)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX tax: <https://gene.ai/taxonomy/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?id ?title ?catLabel WHERE {
  ?source a ds:DataSource ;
          ds:id ?id ;
          ds:title ?title ;
          ds:tier tax:Tier1 ;
          ds:hasAccessMethod ?access ;
          ds:belongsToSubcategory/ds:belongsToCategory ?cat .
  ?access ds:accessMethodType tax:RestApi .
  ?cat skos:prefLabel ?catLabel .
}
ORDER BY ?catLabel ?title
```

### Example 3: Source Details

**User**: "Tell me about DrugBank"

**Intent**: details (single source)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?title ?description ?website ?tierLabel ?statusLabel
       ?catLabel ?subcatLabel ?licenseName ?licenseType WHERE {
  ?source ds:id "drugbank" ;
          ds:title ?title ;
          ds:tier ?tier ;
          ds:status ?status ;
          ds:belongsToSubcategory ?subcat .
  ?subcat ds:belongsToCategory ?cat ;
          skos:prefLabel ?subcatLabel .
  ?cat skos:prefLabel ?catLabel .
  ?tier skos:prefLabel ?tierLabel .
  ?status skos:prefLabel ?statusLabel .
  OPTIONAL { ?source ds:description ?description }
  OPTIONAL { ?source ds:website ?website }
  OPTIONAL {
    ?source ds:hasLicense ?lic .
    ?lic ds:licenseName ?licenseName ;
         ds:licenseType/skos:prefLabel ?licenseType .
  }
}
```

### Example 4: Access Methods for Source

**User**: "How do I access UniProt?"

**Intent**: details (access methods)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?methodLabel ?endpoint ?authRequired ?formats WHERE {
  ?source ds:id "uniprot" ;
          ds:hasAccessMethod ?method .
  ?method ds:accessMethodType ?type .
  ?type skos:prefLabel ?methodLabel .
  OPTIONAL { ?method ds:apiEndpoint ?endpoint }
  OPTIONAL { ?method ds:authenticationRequired ?authRequired }
  OPTIONAL {
    SELECT ?method (GROUP_CONCAT(?fmtLabel; separator=", ") AS ?formats) WHERE {
      ?method ds:hasFormat/ds:formatType/skos:prefLabel ?fmtLabel .
    }
    GROUP BY ?method
  }
}
```

### Example 5: Compare Sources

**User**: "Compare ClinVar and COSMIC"

**Intent**: compare (two sources)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?id ?title ?category ?tier ?license ?accessMethods WHERE {
  VALUES ?id { "clinvar" "cosmic" }
  ?source ds:id ?id ;
          ds:title ?title ;
          ds:tier/skos:prefLabel ?tier ;
          ds:belongsToSubcategory/ds:belongsToCategory/skos:prefLabel ?category .
  OPTIONAL { ?source ds:hasLicense/ds:licenseName ?license }
  OPTIONAL {
    SELECT ?source (GROUP_CONCAT(?methLabel; separator=", ") AS ?accessMethods) WHERE {
      ?source ds:hasAccessMethod/ds:accessMethodType/skos:prefLabel ?methLabel .
    }
    GROUP BY ?source
  }
}
ORDER BY ?id
```

### Example 6: License Filter

**User**: "Which sources have open licenses?"

**Intent**: search (license filter)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX tax: <https://gene.ai/taxonomy/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?id ?title ?category ?licenseName WHERE {
  ?source a ds:DataSource ;
          ds:id ?id ;
          ds:title ?title ;
          ds:belongsToSubcategory/ds:belongsToCategory/skos:prefLabel ?category ;
          ds:hasLicense ?lic .
  ?lic ds:licenseName ?licenseName ;
       ds:commercialUseAllowed true .
}
ORDER BY ?category ?title
```

### Example 7: Statistics

**User**: "How many sources are in each category?"

**Intent**: stats (count by category)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?category (COUNT(DISTINCT ?source) AS ?count) WHERE {
  ?source a ds:DataSource ;
          ds:belongsToSubcategory/ds:belongsToCategory ?cat .
  ?cat skos:prefLabel ?category .
}
GROUP BY ?category
ORDER BY DESC(?count)
```

### Example 8: Find Related

**User**: "What sources are related to KEGG?"

**Intent**: related (relationships)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?relatedId ?relatedTitle ?relationship ?category WHERE {
  ?source ds:id "kegg" .
  {
    ?source ds:relatedTo ?related .
    BIND("related to" AS ?relationship)
  } UNION {
    ?source ds:complementedBy ?related .
    BIND("complemented by" AS ?relationship)
  } UNION {
    ?related ds:complementedBy ?source .
    BIND("complements" AS ?relationship)
  }
  ?related ds:id ?relatedId ;
           ds:title ?relatedTitle ;
           ds:belongsToSubcategory/ds:belongsToCategory/skos:prefLabel ?category .
}
```

### Example 9: Keyword Search

**User**: "Find sources about cancer genomics"

**Intent**: search (keyword/text)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

SELECT ?id ?title ?category ?tier WHERE {
  ?source a ds:DataSource ;
          ds:id ?id ;
          ds:title ?title ;
          ds:tier/skos:notation ?tier ;
          ds:belongsToSubcategory/ds:belongsToCategory/skos:prefLabel ?category .
  {
    ?source ds:tag ?tag .
    FILTER(CONTAINS(LCASE(?tag), "cancer"))
  } UNION {
    ?source ds:description ?desc .
    FILTER(CONTAINS(LCASE(?desc), "cancer"))
  } UNION {
    ?source ds:belongsToSubcategory/skos:prefLabel ?subcatLabel .
    FILTER(CONTAINS(LCASE(?subcatLabel), "cancer"))
  }
}
ORDER BY ?tier ?title
```

### Example 10: Top N by Property

**User**: "Show the 10 largest databases"

**Intent**: search (sorted, limited)

**Generated SPARQL**:
```sparql
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

SELECT ?id ?title ?category ?recordCount WHERE {
  ?source a ds:DataSource ;
          ds:id ?id ;
          ds:title ?title ;
          ds:belongsToSubcategory/ds:belongsToCategory/skos:prefLabel ?category ;
          ds:currentVersion/ds:recordCount ?recordCount .
}
ORDER BY DESC(?recordCount)
LIMIT 10
```

## Category Name Mapping

When users mention categories, map to taxonomy URIs:

| User Says | Maps To |
|-----------|---------|
| "TCM", "Traditional Chinese Medicine", "Chinese medicine" | `tax:TraditionalMedicine` or `tax:TraditionalChineseMedicine` |
| "genomics", "genetics", "variants" | `tax:GeneticsGenomics` |
| "drugs", "pharmaceuticals", "compounds" | `tax:CompoundsMolecules` |
| "diseases", "phenotypes" | `tax:DiseasesPhenotypes` |
| "pathways", "networks", "interactions" | `tax:PathwaysNetworks` |
| "nutrition", "food" | `tax:NutritionFood` |
| "proteins", "molecular biology" | `tax:ProteinsMolecularBiology` |
| "literature", "publications", "knowledge" | `tax:LiteratureKnowledge` |
| "microbiome", "gut bacteria" | `tax:Microbiome` |

## Tier Mapping

| User Says | Maps To |
|-----------|---------|
| "primary", "tier 1", "authoritative" | `tax:Tier1` |
| "secondary", "tier 2", "supplementary" | `tax:Tier2` |
| "specialized", "tier 3", "niche" | `tax:Tier3` |

## Access Method Mapping

| User Says | Maps To |
|-----------|---------|
| "API", "REST", "web service" | `tax:RestApi` |
| "SPARQL", "RDF endpoint" | `tax:SparqlEndpoint` |
| "download", "bulk" | `tax:DirectDownload` or `tax:FtpDownload` |
| "FTP" | `tax:FtpDownload` |

## Result Formatting Rules

After executing the SPARQL query:

1. **Never show raw URIs** - Convert to human-readable labels
2. **Never show SPARQL** - Hide implementation details
3. **Group by category** when showing multiple sources
4. **Sort by tier** (1 before 2 before 3) within categories
5. **Use markdown tables** for comparisons
6. **Use bullet lists** for enumeration
7. **Show key properties** only (id, title, tier, category) in lists
8. **Show full details** only when asked about a specific source
