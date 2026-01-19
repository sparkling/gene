# Extracting Pathway and Biological Process Data from Wikidata

## Overview

Wikidata serves as a comprehensive knowledge graph that contains structured data about biological pathways, processes, and their relationships to genes and compounds. This document provides detailed guidance on extracting pathway-related information from Wikidata, including SPARQL queries, bulk extraction methods, and integration strategies.

---

## 1. Pathway Items in Wikidata

### Core Pathway Classes

Wikidata organizes pathway-related concepts using a class hierarchy. The following Q-identifiers represent the primary pathway types:

| Q-ID | Label | Description |
|------|-------|-------------|
| Q4915012 | biological pathway | A series of actions among molecules in a cell that leads to a product or change |
| Q2996394 | biological process | A process specifically pertinent to the functioning of integrated living units |
| Q4915059 | metabolic pathway | A series of chemical reactions occurring within a cell |
| Q189093 | signal transduction pathway | Process by which a cell converts an extracellular signal to a response |

### Class Hierarchy

```
Q2996394 (biological process)
├── Q4915012 (biological pathway)
│   ├── Q4915059 (metabolic pathway)
│   │   ├── Q27141728 (catabolic pathway)
│   │   ├── Q27141731 (anabolic pathway)
│   │   └── Q27141734 (amphibolic pathway)
│   ├── Q189093 (signal transduction pathway)
│   │   ├── Q14860489 (MAPK signaling pathway)
│   │   ├── Q420949 (Wnt signaling pathway)
│   │   └── Q14860445 (Notch signaling pathway)
│   └── Q14599311 (regulatory pathway)
```

### Counting Items Per Type

Use this SPARQL query to count items for each pathway type:

```sparql
SELECT ?type ?typeLabel (COUNT(?item) AS ?count)
WHERE {
  VALUES ?type { wd:Q4915012 wd:Q2996394 wd:Q4915059 wd:Q189093 }
  ?item wdt:P31 ?type .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
GROUP BY ?type ?typeLabel
ORDER BY DESC(?count)
```

#### Expected Results (Approximate Counts)

| Type | Approximate Count |
|------|-------------------|
| Q2996394 (biological process) | ~30,000+ |
| Q4915012 (biological pathway) | ~3,000+ |
| Q4915059 (metabolic pathway) | ~1,500+ |
| Q189093 (signal transduction pathway) | ~500+ |

### Query Including Subclasses

To get comprehensive counts including subclasses:

```sparql
SELECT ?type ?typeLabel (COUNT(DISTINCT ?item) AS ?count)
WHERE {
  VALUES ?type { wd:Q4915012 wd:Q2996394 wd:Q4915059 wd:Q189093 }
  ?item wdt:P31/wdt:P279* ?type .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
GROUP BY ?type ?typeLabel
ORDER BY DESC(?count)
```

---

## 2. Pathway Properties

### External Identifier Properties

| Property | Label | Description | Format |
|----------|-------|-------------|--------|
| P3937 | Reactome pathway ID | Identifier in Reactome database | R-HSA-nnnnnn |
| P2888 | exact match | Used for KEGG pathway IDs | External URI |
| P2410 | WikiPathways ID | Identifier in WikiPathways | WPnnn |
| P686 | Gene Ontology ID | GO term identifier | GO:nnnnnnn |

### Structural Properties

| Property | Label | Description | Usage |
|----------|-------|-------------|-------|
| P361 | part of | Indicates parent pathway | Pathway hierarchy |
| P527 | has part | Lists pathway components | Genes, compounds, sub-pathways |
| P1542 | has effect | Downstream effects | Regulatory relationships |
| P1479 | has contributing factor | Upstream regulators | Input signals |

### Query: Pathways with External Identifiers

```sparql
SELECT ?pathway ?pathwayLabel ?reactome ?wikipathways ?go
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .
  OPTIONAL { ?pathway wdt:P3937 ?reactome . }
  OPTIONAL { ?pathway wdt:P2410 ?wikipathways . }
  OPTIONAL { ?pathway wdt:P686 ?go . }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 1000
```

### Query: Pathway Hierarchy Using P361/P527

```sparql
SELECT ?parent ?parentLabel ?child ?childLabel
WHERE {
  ?child wdt:P31/wdt:P279* wd:Q4915012 ;
         wdt:P361 ?parent .
  ?parent wdt:P31/wdt:P279* wd:Q4915012 .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 500
```

---

## 3. Gene-Pathway Links in Wikidata

### Gene Ontology Annotation Properties

| Property | Label | GO Aspect | Description |
|----------|-------|-----------|-------------|
| P682 | biological process | BP | Links gene to biological process involvement |
| P680 | molecular function | MF | Links gene to its molecular activities |
| P681 | cellular component | CC | Links gene to cellular locations |

### How Genes Link to Pathways

Genes connect to pathways through multiple mechanisms:

1. **Direct GO Annotations**: Genes have P682 (biological process) statements
2. **Pathway Components**: Pathways list genes via P527 (has part)
3. **Protein Involvement**: Proteins encoded by genes participate in pathways

### Query: Genes in a Specific Pathway via GO

```sparql
SELECT DISTINCT ?gene ?geneLabel ?geneSymbol ?process ?processLabel
WHERE {
  # Get human genes
  ?gene wdt:P31 wd:Q7187 ;      # instance of: gene
        wdt:P703 wd:Q15978631 ; # found in taxon: Homo sapiens
        wdt:P682 ?process .     # biological process

  # Filter for specific pathway/process
  ?process wdt:P31/wdt:P279* wd:Q4915012 .  # is a biological pathway

  OPTIONAL { ?gene wdt:P353 ?geneSymbol . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 500
```

### Query: All Genes with Pathway Annotations

```sparql
SELECT ?gene ?geneLabel ?symbol ?pathway ?pathwayLabel
WHERE {
  ?gene wdt:P31 wd:Q7187 ;
        wdt:P703 wd:Q15978631 ;
        wdt:P353 ?symbol ;
        wdt:P682 ?pathway .

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 10000
```

### Query: Genes as Pathway Components (P527)

```sparql
SELECT ?pathway ?pathwayLabel ?gene ?geneLabel ?geneSymbol
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 ;
           wdt:P527 ?gene .
  ?gene wdt:P31 wd:Q7187 ;
        wdt:P703 wd:Q15978631 .
  OPTIONAL { ?gene wdt:P353 ?geneSymbol . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 1000
```

### Query: Proteins in Pathways

```sparql
SELECT ?pathway ?pathwayLabel ?protein ?proteinLabel ?gene ?geneLabel
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 ;
           wdt:P527 ?protein .
  ?protein wdt:P31 wd:Q8054 .  # instance of: protein

  OPTIONAL {
    ?protein wdt:P702 ?gene .  # encoded by
    ?gene wdt:P703 wd:Q15978631 .
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 1000
```

---

## 4. Compound-Pathway Links

### Relevant Properties for Compounds

| Property | Label | Description |
|----------|-------|-------------|
| P703 | found in taxon | Species where compound is found |
| P2868 | subject has role | Functional role (e.g., metabolite, drug) |
| P3780 | active ingredient in | Drug products containing compound |
| P129 | physically interacts with | Molecular interactions |
| P361 | part of | Membership in pathway/process |

### Metabolite Roles (Q-IDs)

| Q-ID | Label | Description |
|------|-------|-------------|
| Q11173 | chemical compound | General compound |
| Q113145171 | metabolite | Compound involved in metabolism |
| Q407595 | cofactor | Essential for enzyme activity |
| Q213449 | intermediate | Metabolic intermediate |

### Query: Metabolites in Pathways

```sparql
SELECT ?pathway ?pathwayLabel ?compound ?compoundLabel ?role ?roleLabel
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915059 ;  # metabolic pathway
           wdt:P527 ?compound .

  ?compound wdt:P31/wdt:P279* wd:Q11173 .  # chemical compound

  OPTIONAL { ?compound wdt:P2868 ?role . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 1000
```

### Query: Compounds with Specific Roles in Pathways

```sparql
SELECT ?pathway ?pathwayLabel ?compound ?compoundLabel ?inchi ?smiles
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915059 ;
           wdt:P527 ?compound .

  ?compound wdt:P31/wdt:P279* wd:Q11173 ;
            wdt:P2868 wd:Q113145171 .  # has role: metabolite

  OPTIONAL { ?compound wdt:P234 ?inchi . }
  OPTIONAL { ?compound wdt:P233 ?smiles . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 500
```

### Query: Compounds Found in Human (P703)

```sparql
SELECT ?compound ?compoundLabel ?pathway ?pathwayLabel
WHERE {
  ?compound wdt:P703 wd:Q15978631 ;  # found in: Homo sapiens
            wdt:P361 ?pathway .      # part of pathway

  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 500
```

---

## 5. SPARQL Queries

### 5.1 All Human Metabolic Pathways

```sparql
SELECT DISTINCT ?pathway ?pathwayLabel ?reactome ?wikipathways ?description
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915059 .  # metabolic pathway

  # Human-specific (has human genes or explicitly linked to human)
  {
    ?pathway wdt:P703 wd:Q15978631 .
  }
  UNION
  {
    ?pathway wdt:P527 ?component .
    ?component wdt:P703 wd:Q15978631 .
  }
  UNION
  {
    ?pathway wdt:P3937 ?reactome .
    FILTER(STRSTARTS(?reactome, "R-HSA-"))  # Human Reactome IDs start with R-HSA
  }

  OPTIONAL { ?pathway wdt:P3937 ?reactome . }
  OPTIONAL { ?pathway wdt:P2410 ?wikipathways . }
  OPTIONAL { ?pathway schema:description ?description . FILTER(LANG(?description) = "en") }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel
```

### 5.2 All Signaling Pathways with Genes

```sparql
SELECT ?pathway ?pathwayLabel ?gene ?geneLabel ?symbol
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q189093 .  # signal transduction pathway

  {
    ?pathway wdt:P527 ?gene .
    ?gene wdt:P31 wd:Q7187 ;
          wdt:P703 wd:Q15978631 .
  }
  UNION
  {
    ?gene wdt:P682 ?pathway ;
          wdt:P31 wd:Q7187 ;
          wdt:P703 wd:Q15978631 .
  }

  OPTIONAL { ?gene wdt:P353 ?symbol . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel ?symbol
LIMIT 5000
```

### 5.3 Compounds Involved in Specific Pathway

```sparql
# Example: Glycolysis pathway (adjust Q-ID as needed)
SELECT ?compound ?compoundLabel ?inchi ?chebi ?pubchem ?role ?roleLabel
WHERE {
  # Glycolysis - use appropriate Q-ID
  BIND(wd:Q172788 AS ?pathway)

  {
    ?pathway wdt:P527 ?compound .
    ?compound wdt:P31/wdt:P279* wd:Q11173 .
  }
  UNION
  {
    ?compound wdt:P361 ?pathway ;
              wdt:P31/wdt:P279* wd:Q11173 .
  }

  OPTIONAL { ?compound wdt:P234 ?inchi . }
  OPTIONAL { ?compound wdt:P683 ?chebi . }
  OPTIONAL { ?compound wdt:P662 ?pubchem . }
  OPTIONAL { ?compound wdt:P2868 ?role . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
```

### 5.4 Pathway Hierarchy Extraction

```sparql
SELECT ?pathway ?pathwayLabel ?parent ?parentLabel ?level
WHERE {
  # Get pathways with their hierarchy
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .

  # Optional parent pathway
  OPTIONAL {
    ?pathway wdt:P361 ?parent .
    ?parent wdt:P31/wdt:P279* wd:Q4915012 .
  }

  # Calculate hierarchy level (simplified)
  BIND(IF(BOUND(?parent), 2, 1) AS ?level)

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?parent ?pathway
LIMIT 2000
```

### Full Recursive Hierarchy

```sparql
SELECT ?pathway ?pathwayLabel ?parent ?parentLabel ?depth
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .

  # Get hierarchy depth using property path
  ?pathway wdt:P361* ?ancestor .
  ?ancestor wdt:P31/wdt:P279* wd:Q4915012 .

  OPTIONAL { ?pathway wdt:P361 ?parent . }

  # Count depth
  {
    SELECT ?pathway (COUNT(?mid) AS ?depth)
    WHERE {
      ?pathway wdt:P361* ?mid .
      ?mid wdt:P31/wdt:P279* wd:Q4915012 .
    }
    GROUP BY ?pathway
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?depth ?parent ?pathway
LIMIT 1000
```

### 5.5 Cross-Reference to Reactome/KEGG/WikiPathways

```sparql
SELECT ?pathway ?pathwayLabel ?reactome ?kegg ?wikipathways ?go
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .

  # Must have at least one external ID
  {
    { ?pathway wdt:P3937 ?reactome . }
    UNION
    { ?pathway wdt:P2410 ?wikipathways . }
    UNION
    { ?pathway wdt:P686 ?go . }
  }

  OPTIONAL { ?pathway wdt:P3937 ?reactome . }
  OPTIONAL { ?pathway wdt:P2410 ?wikipathways . }
  OPTIONAL { ?pathway wdt:P686 ?go . }

  # KEGG via exact match (P2888)
  OPTIONAL {
    ?pathway wdt:P2888 ?keggUri .
    FILTER(CONTAINS(STR(?keggUri), "kegg.jp/pathway"))
    BIND(REPLACE(STR(?keggUri), ".*pathway/", "") AS ?kegg)
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel
LIMIT 5000
```

### 5.6 Comprehensive Pathway Export Query

```sparql
SELECT ?pathway ?pathwayLabel ?type ?typeLabel
       ?reactome ?wikipathways ?go ?kegg
       (GROUP_CONCAT(DISTINCT ?geneSymbol; separator="|") AS ?genes)
       (GROUP_CONCAT(DISTINCT ?compoundLabel; separator="|") AS ?compounds)
WHERE {
  ?pathway wdt:P31 ?type .
  ?type wdt:P279* wd:Q4915012 .

  OPTIONAL { ?pathway wdt:P3937 ?reactome . }
  OPTIONAL { ?pathway wdt:P2410 ?wikipathways . }
  OPTIONAL { ?pathway wdt:P686 ?go . }

  OPTIONAL {
    ?pathway wdt:P2888 ?keggUri .
    FILTER(CONTAINS(STR(?keggUri), "kegg"))
    BIND(STR(?keggUri) AS ?kegg)
  }

  OPTIONAL {
    ?pathway wdt:P527 ?gene .
    ?gene wdt:P31 wd:Q7187 ;
          wdt:P353 ?geneSymbol .
  }

  OPTIONAL {
    ?pathway wdt:P527 ?compound .
    ?compound wdt:P31/wdt:P279* wd:Q11173 .
    ?compound rdfs:label ?compoundLabel .
    FILTER(LANG(?compoundLabel) = "en")
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
GROUP BY ?pathway ?pathwayLabel ?type ?typeLabel ?reactome ?wikipathways ?go ?kegg
LIMIT 2000
```

---

## 6. Federated Queries

### 6.1 Wikidata + WikiPathways SPARQL

WikiPathways has its own SPARQL endpoint at `https://sparql.wikipathways.org/sparql`. You can create federated queries to combine data from both sources.

```sparql
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX wpTypes: <http://vocabularies.wikipathways.org/wpTypes/>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT ?pathway ?pathwayLabel ?wpId ?wpTitle ?wpOrganism ?geneCount
WHERE {
  # Get WikiPathways ID from Wikidata
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 ;
           wdt:P2410 ?wpId .

  # Federated query to WikiPathways
  SERVICE <https://sparql.wikipathways.org/sparql> {
    ?wpPathway dcterms:identifier ?wpId ;
               dcterms:title ?wpTitle ;
               wp:organism ?wpOrganism .

    # Count genes in pathway
    {
      SELECT ?wpPathway (COUNT(DISTINCT ?gene) AS ?geneCount)
      WHERE {
        ?wpPathway a wp:Pathway .
        ?gene a wp:GeneProduct ;
              dcterms:isPartOf ?wpPathway .
      }
      GROUP BY ?wpPathway
    }
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 100
```

### 6.2 Get Pathway Details from WikiPathways via Wikidata Link

```sparql
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX gpml: <http://vocabularies.wikipathways.org/gpml#>

SELECT ?wikidataPathway ?wikidataLabel ?wpId ?wpTitle ?revision
       ?geneId ?geneLabel ?entrezId
WHERE {
  # Start with Wikidata pathway
  ?wikidataPathway wdt:P31/wdt:P279* wd:Q4915012 ;
                   wdt:P2410 ?wpId ;
                   wdt:P703 wd:Q15978631 .  # Human pathways

  # Get details from WikiPathways
  SERVICE <https://sparql.wikipathways.org/sparql> {
    ?wpPathway dcterms:identifier ?wpId ;
               dcterms:title ?wpTitle ;
               wp:organism <http://purl.obolibrary.org/obo/NCBITaxon_9606> ;  # Human
               dcterms:hasVersion ?revision .

    # Get genes
    ?geneProduct a wp:GeneProduct ;
                 dcterms:isPartOf ?wpPathway ;
                 rdfs:label ?geneLabel .

    OPTIONAL { ?geneProduct wp:bdbEntrezGene ?entrezId . }
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 500
```

### 6.3 Cross-Reference Validation Query

```sparql
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT ?wikidataId ?wikidataLabel ?wpId ?wpTitle ?match
WHERE {
  ?wikidataId wdt:P2410 ?wpId ;
              rdfs:label ?wikidataLabel .
  FILTER(LANG(?wikidataLabel) = "en")

  SERVICE <https://sparql.wikipathways.org/sparql> {
    ?wpPathway dcterms:identifier ?wpId ;
               dcterms:title ?wpTitle .
  }

  # Check if names match
  BIND(IF(LCASE(?wikidataLabel) = LCASE(?wpTitle), "exact",
       IF(CONTAINS(LCASE(?wikidataLabel), LCASE(?wpTitle)) ||
          CONTAINS(LCASE(?wpTitle), LCASE(?wikidataLabel)), "partial", "mismatch")) AS ?match)
}
ORDER BY ?match
LIMIT 200
```

### 6.4 Federated Query with Reactome (via Content Negotiation)

Note: Reactome does not have a public SPARQL endpoint, but you can use their REST API or RDF downloads. For SPARQL-based integration, use Wikidata's cross-references:

```sparql
SELECT ?pathway ?pathwayLabel ?reactomeId
       (CONCAT("https://reactome.org/content/detail/", ?reactomeId) AS ?reactomeUrl)
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 ;
           wdt:P3937 ?reactomeId .

  # Filter for human pathways
  FILTER(STRSTARTS(?reactomeId, "R-HSA-"))

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel
```

---

## 7. Bulk Extraction from Dump

### 7.1 Wikidata Dump Sources

Wikidata provides regular dumps at:
- **JSON dumps**: `https://dumps.wikimedia.org/wikidatawiki/entities/`
- **Compressed files**: `latest-all.json.bz2` (~100GB compressed)
- **Truthy dumps**: `latest-truthy.nt.bz2` (statements only)

### 7.2 Filter for Pathway Q-IDs

```python
#!/usr/bin/env python3
"""
Extract pathway-related items from Wikidata JSON dump.

This script filters the Wikidata dump for items that are:
1. Instances of pathway-related classes
2. Have pathway-related properties (P3937, P2410, P686)
3. Participate in biological processes
"""

import bz2
import json
import sys
from pathlib import Path
from typing import Iterator, Dict, Any, Set
from collections import defaultdict

# Pathway-related Q-IDs (classes)
PATHWAY_CLASSES = {
    'Q4915012',   # biological pathway
    'Q2996394',   # biological process
    'Q4915059',   # metabolic pathway
    'Q189093',    # signal transduction pathway
    'Q27141728',  # catabolic pathway
    'Q27141731',  # anabolic pathway
    'Q14860489',  # MAPK signaling pathway
    'Q420949',    # Wnt signaling pathway
}

# Properties indicating pathway involvement
PATHWAY_PROPERTIES = {
    'P3937',  # Reactome pathway ID
    'P2410',  # WikiPathways ID
    'P686',   # Gene Ontology ID
    'P682',   # biological process
    'P680',   # molecular function
    'P681',   # cellular component
}

# Properties for extracting relationships
RELATIONSHIP_PROPERTIES = {
    'P361',  # part of
    'P527',  # has part
    'P703',  # found in taxon
    'P702',  # encoded by
    'P688',  # encodes
}


def stream_wikidata_dump(dump_path: str) -> Iterator[Dict[str, Any]]:
    """Stream items from a bz2-compressed Wikidata JSON dump."""
    with bz2.open(dump_path, 'rt', encoding='utf-8') as f:
        # Skip opening bracket
        f.readline()

        for line in f:
            line = line.strip()
            if line in ['[', ']', '']:
                continue

            # Remove trailing comma if present
            if line.endswith(','):
                line = line[:-1]

            try:
                item = json.loads(line)
                yield item
            except json.JSONDecodeError:
                continue


def is_pathway_item(item: Dict[str, Any]) -> bool:
    """Check if an item is pathway-related."""
    if 'claims' not in item:
        return False

    claims = item['claims']

    # Check if instance of (P31) includes pathway classes
    if 'P31' in claims:
        for claim in claims['P31']:
            try:
                qid = claim['mainsnak']['datavalue']['value']['id']
                if qid in PATHWAY_CLASSES:
                    return True
            except (KeyError, TypeError):
                continue

    # Check if subclass of (P279) pathway classes
    if 'P279' in claims:
        for claim in claims['P279']:
            try:
                qid = claim['mainsnak']['datavalue']['value']['id']
                if qid in PATHWAY_CLASSES:
                    return True
            except (KeyError, TypeError):
                continue

    # Check for pathway-specific properties
    for prop in ['P3937', 'P2410']:  # Reactome, WikiPathways
        if prop in claims:
            return True

    return False


def has_pathway_annotation(item: Dict[str, Any]) -> bool:
    """Check if item has pathway/GO annotations (for genes)."""
    if 'claims' not in item:
        return False

    claims = item['claims']

    # Check for GO annotations
    for prop in ['P682', 'P680', 'P681', 'P686']:
        if prop in claims:
            return True

    return False


def extract_pathway_data(item: Dict[str, Any]) -> Dict[str, Any]:
    """Extract relevant pathway data from an item."""
    qid = item.get('id', '')

    result = {
        'qid': qid,
        'type': item.get('type', ''),
        'labels': {},
        'descriptions': {},
        'aliases': {},
        'claims': defaultdict(list),
    }

    # Extract labels
    labels = item.get('labels', {})
    for lang in ['en', 'de', 'fr', 'es', 'ja', 'zh']:
        if lang in labels:
            result['labels'][lang] = labels[lang].get('value', '')

    # Extract descriptions
    descriptions = item.get('descriptions', {})
    for lang in ['en']:
        if lang in descriptions:
            result['descriptions'][lang] = descriptions[lang].get('value', '')

    # Extract relevant claims
    claims = item.get('claims', {})

    properties_to_extract = (
        PATHWAY_PROPERTIES |
        RELATIONSHIP_PROPERTIES |
        {'P31', 'P279', 'P353', 'P352', 'P234', 'P233', 'P683', 'P662'}
    )

    for prop in properties_to_extract:
        if prop in claims:
            for claim in claims[prop]:
                try:
                    snak = claim['mainsnak']
                    if snak['snaktype'] == 'value':
                        datavalue = snak['datavalue']

                        if datavalue['type'] == 'wikibase-entityid':
                            value = datavalue['value']['id']
                        elif datavalue['type'] == 'string':
                            value = datavalue['value']
                        elif datavalue['type'] == 'monolingualtext':
                            value = datavalue['value']['text']
                        else:
                            value = str(datavalue['value'])

                        result['claims'][prop].append(value)
                except (KeyError, TypeError):
                    continue

    return result


def extract_pathways_from_dump(
    dump_path: str,
    output_dir: str,
    include_genes: bool = True
) -> Dict[str, int]:
    """
    Extract all pathway-related data from Wikidata dump.

    Args:
        dump_path: Path to Wikidata JSON dump (bz2 compressed)
        output_dir: Directory for output files
        include_genes: Also extract genes with pathway annotations

    Returns:
        Statistics about extracted data
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    stats = defaultdict(int)

    # Output files
    pathways_file = open(output_path / 'pathways.jsonl', 'w')
    genes_file = open(output_path / 'genes_with_pathways.jsonl', 'w') if include_genes else None
    relationships_file = open(output_path / 'pathway_relationships.jsonl', 'w')

    print(f"Processing dump: {dump_path}")

    for i, item in enumerate(stream_wikidata_dump(dump_path)):
        if i % 1000000 == 0:
            print(f"Processed {i:,} items, found {stats['pathways']} pathways, "
                  f"{stats['genes']} genes")

        # Check if pathway item
        if is_pathway_item(item):
            data = extract_pathway_data(item)
            pathways_file.write(json.dumps(data) + '\n')
            stats['pathways'] += 1

            # Extract relationships
            claims = item.get('claims', {})
            qid = item['id']

            for prop in ['P361', 'P527']:
                if prop in claims:
                    for claim in claims[prop]:
                        try:
                            target = claim['mainsnak']['datavalue']['value']['id']
                            rel = {
                                'source': qid,
                                'property': prop,
                                'target': target
                            }
                            relationships_file.write(json.dumps(rel) + '\n')
                            stats['relationships'] += 1
                        except (KeyError, TypeError):
                            continue

        # Check if gene with pathway annotation
        elif include_genes and has_pathway_annotation(item):
            claims = item.get('claims', {})

            # Verify it's a gene
            is_gene = False
            if 'P31' in claims:
                for claim in claims['P31']:
                    try:
                        qid = claim['mainsnak']['datavalue']['value']['id']
                        if qid == 'Q7187':  # gene
                            is_gene = True
                            break
                    except (KeyError, TypeError):
                        continue

            if is_gene:
                data = extract_pathway_data(item)
                genes_file.write(json.dumps(data) + '\n')
                stats['genes'] += 1

    # Close files
    pathways_file.close()
    if genes_file:
        genes_file.close()
    relationships_file.close()

    # Write statistics
    with open(output_path / 'extraction_stats.json', 'w') as f:
        json.dump(dict(stats), f, indent=2)

    return dict(stats)


def build_pathway_gene_matrix(
    pathways_file: str,
    genes_file: str,
    output_file: str
):
    """Build a pathway-gene relationship matrix."""

    # Load pathways
    pathways = {}
    with open(pathways_file) as f:
        for line in f:
            item = json.loads(line)
            qid = item['qid']
            label = item['labels'].get('en', qid)
            pathways[qid] = {
                'label': label,
                'genes': set(),
                'reactome': item['claims'].get('P3937', []),
                'wikipathways': item['claims'].get('P2410', []),
            }

    # Load genes and their pathway annotations
    genes = {}
    with open(genes_file) as f:
        for line in f:
            item = json.loads(line)
            qid = item['qid']
            symbol = item['claims'].get('P353', [''])[0]

            genes[qid] = {
                'symbol': symbol,
                'label': item['labels'].get('en', qid),
                'pathways': set()
            }

            # Link to pathways via P682 (biological process)
            for pathway_qid in item['claims'].get('P682', []):
                if pathway_qid in pathways:
                    pathways[pathway_qid]['genes'].add(qid)
                    genes[qid]['pathways'].add(pathway_qid)

    # Build matrix output
    result = {
        'pathways': {},
        'genes': {},
        'matrix': []
    }

    for qid, data in pathways.items():
        result['pathways'][qid] = {
            'label': data['label'],
            'reactome': data['reactome'],
            'wikipathways': data['wikipathways'],
            'gene_count': len(data['genes'])
        }

        for gene_qid in data['genes']:
            result['matrix'].append({
                'pathway': qid,
                'gene': gene_qid,
                'gene_symbol': genes[gene_qid]['symbol']
            })

    for qid, data in genes.items():
        if data['pathways']:
            result['genes'][qid] = {
                'symbol': data['symbol'],
                'label': data['label'],
                'pathway_count': len(data['pathways'])
            }

    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)

    return len(result['matrix'])


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python extract_pathways.py <dump_path> <output_dir>")
        sys.exit(1)

    dump_path = sys.argv[1]
    output_dir = sys.argv[2]

    stats = extract_pathways_from_dump(dump_path, output_dir)

    print("\nExtraction complete!")
    print(f"Pathways: {stats.get('pathways', 0)}")
    print(f"Genes with annotations: {stats.get('genes', 0)}")
    print(f"Relationships: {stats.get('relationships', 0)}")
```

### 7.3 Extract Pathway-Gene Relationships

```python
#!/usr/bin/env python3
"""
Extract pathway-gene relationships from pre-extracted Wikidata files.
"""

import json
from pathlib import Path
from collections import defaultdict
import csv


def load_jsonl(filepath: str):
    """Load a JSONL file."""
    items = []
    with open(filepath) as f:
        for line in f:
            items.append(json.loads(line))
    return items


def extract_relationships(data_dir: str, output_dir: str):
    """Extract and organize pathway-gene relationships."""

    data_path = Path(data_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading pathways...")
    pathways = {item['qid']: item for item in load_jsonl(data_path / 'pathways.jsonl')}

    print("Loading genes...")
    genes = {item['qid']: item for item in load_jsonl(data_path / 'genes_with_pathways.jsonl')}

    print("Loading relationships...")
    relationships = load_jsonl(data_path / 'pathway_relationships.jsonl')

    # Build relationship maps
    pathway_genes = defaultdict(set)  # pathway -> genes
    pathway_compounds = defaultdict(set)  # pathway -> compounds
    pathway_hierarchy = defaultdict(set)  # pathway -> parent pathways

    # Process P527 (has part) and P361 (part of) relationships
    for rel in relationships:
        source = rel['source']
        target = rel['target']
        prop = rel['property']

        if source in pathways:
            if prop == 'P527':  # has part
                if target in genes:
                    pathway_genes[source].add(target)
            elif prop == 'P361':  # part of
                if target in pathways:
                    pathway_hierarchy[source].add(target)

    # Process gene P682 (biological process) annotations
    for gene_qid, gene_data in genes.items():
        for process_qid in gene_data['claims'].get('P682', []):
            if process_qid in pathways:
                pathway_genes[process_qid].add(gene_qid)

    # Output: Pathway-Gene TSV
    print("Writing pathway-gene relationships...")
    with open(output_path / 'pathway_genes.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['pathway_qid', 'pathway_label', 'gene_qid', 'gene_symbol',
                        'reactome_id', 'wikipathways_id'])

        for pathway_qid, gene_set in sorted(pathway_genes.items()):
            pathway = pathways.get(pathway_qid, {})
            pathway_label = pathway.get('labels', {}).get('en', '')
            reactome = ','.join(pathway.get('claims', {}).get('P3937', []))
            wikipathways = ','.join(pathway.get('claims', {}).get('P2410', []))

            for gene_qid in sorted(gene_set):
                gene = genes.get(gene_qid, {})
                gene_symbol = gene.get('claims', {}).get('P353', [''])[0]

                writer.writerow([
                    pathway_qid, pathway_label, gene_qid, gene_symbol,
                    reactome, wikipathways
                ])

    # Output: Pathway hierarchy TSV
    print("Writing pathway hierarchy...")
    with open(output_path / 'pathway_hierarchy.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['child_qid', 'child_label', 'parent_qid', 'parent_label'])

        for child_qid, parent_set in sorted(pathway_hierarchy.items()):
            child = pathways.get(child_qid, {})
            child_label = child.get('labels', {}).get('en', '')

            for parent_qid in sorted(parent_set):
                parent = pathways.get(parent_qid, {})
                parent_label = parent.get('labels', {}).get('en', '')

                writer.writerow([child_qid, child_label, parent_qid, parent_label])

    # Output: Summary statistics
    stats = {
        'total_pathways': len(pathways),
        'pathways_with_genes': len(pathway_genes),
        'total_genes': len(genes),
        'total_pathway_gene_links': sum(len(g) for g in pathway_genes.values()),
        'pathways_with_parents': len(pathway_hierarchy),
    }

    print("Writing statistics...")
    with open(output_path / 'relationship_stats.json', 'w') as f:
        json.dump(stats, f, indent=2)

    print("\nExtraction complete!")
    for key, value in stats.items():
        print(f"  {key}: {value:,}")


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 3:
        print("Usage: python extract_relationships.py <data_dir> <output_dir>")
        sys.exit(1)

    extract_relationships(sys.argv[1], sys.argv[2])
```

### 7.4 Incremental Updates Using Recent Changes API

```python
#!/usr/bin/env python3
"""
Monitor Wikidata for pathway-related changes using the Recent Changes API.
"""

import requests
import json
from datetime import datetime, timedelta
from typing import List, Dict, Any


WIKIDATA_API = "https://www.wikidata.org/w/api.php"

PATHWAY_PROPERTIES = ['P3937', 'P2410', 'P686', 'P682', 'P680', 'P681']


def get_recent_changes(
    hours: int = 24,
    properties: List[str] = None
) -> List[Dict[str, Any]]:
    """Get recent changes to items with specific properties."""

    if properties is None:
        properties = PATHWAY_PROPERTIES

    # Calculate timestamp
    start_time = datetime.utcnow() - timedelta(hours=hours)
    timestamp = start_time.strftime("%Y%m%d%H%M%S")

    changes = []

    for prop in properties:
        params = {
            'action': 'query',
            'list': 'recentchanges',
            'rcprop': 'title|timestamp|user|comment',
            'rctype': 'edit',
            'rcnamespace': 0,  # Main namespace (items)
            'rclimit': 500,
            'rcstart': timestamp,
            'format': 'json',
        }

        response = requests.get(WIKIDATA_API, params=params)
        data = response.json()

        for change in data.get('query', {}).get('recentchanges', []):
            # Check if change involves our property
            if prop.lower() in change.get('comment', '').lower():
                changes.append({
                    'qid': change['title'],
                    'timestamp': change['timestamp'],
                    'property': prop,
                    'user': change.get('user', ''),
                    'comment': change.get('comment', '')
                })

    return changes


def get_item_details(qid: str) -> Dict[str, Any]:
    """Get current state of a Wikidata item."""

    params = {
        'action': 'wbgetentities',
        'ids': qid,
        'format': 'json',
        'props': 'labels|descriptions|claims',
        'languages': 'en'
    }

    response = requests.get(WIKIDATA_API, params=params)
    data = response.json()

    return data.get('entities', {}).get(qid, {})


if __name__ == '__main__':
    print("Checking for recent pathway-related changes...")

    changes = get_recent_changes(hours=24)

    print(f"\nFound {len(changes)} changes in the last 24 hours:")

    for change in changes[:20]:
        print(f"  {change['qid']}: {change['property']} - {change['comment'][:60]}...")
```

---

## 8. Coverage Assessment

### 8.1 Compare Wikidata Pathways to Reactome

```sparql
# Count Wikidata pathways with Reactome IDs vs total human Reactome pathways
SELECT
  (COUNT(DISTINCT ?pathway) AS ?wikidata_with_reactome)
WHERE {
  ?pathway wdt:P3937 ?reactomeId .
  FILTER(STRSTARTS(?reactomeId, "R-HSA-"))
}
```

### 8.2 Python Script for Coverage Analysis

```python
#!/usr/bin/env python3
"""
Assess coverage of Wikidata pathway data compared to external databases.
"""

import requests
import json
from collections import defaultdict
from typing import Dict, Set, Tuple


def get_reactome_human_pathways() -> Set[str]:
    """Get all human pathway IDs from Reactome API."""

    url = "https://reactome.org/ContentService/data/pathways/top/9606"  # 9606 = human
    response = requests.get(url)

    pathways = set()
    for pathway in response.json():
        pathways.add(pathway['stId'])

    # Also get child pathways
    url = "https://reactome.org/ContentService/data/pathways/low/diagram/entity/9606"
    response = requests.get(url)

    for pathway in response.json():
        pathways.add(pathway['stId'])

    return pathways


def get_wikidata_reactome_ids() -> Set[str]:
    """Get Reactome IDs present in Wikidata."""

    query = """
    SELECT DISTINCT ?reactomeId
    WHERE {
      ?pathway wdt:P3937 ?reactomeId .
      FILTER(STRSTARTS(?reactomeId, "R-HSA-"))
    }
    """

    url = "https://query.wikidata.org/sparql"
    response = requests.get(url, params={
        'query': query,
        'format': 'json'
    })

    ids = set()
    for binding in response.json()['results']['bindings']:
        ids.add(binding['reactomeId']['value'])

    return ids


def get_wikipathways_human_pathways() -> Set[str]:
    """Get all human pathway IDs from WikiPathways SPARQL."""

    query = """
    PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
    PREFIX dcterms: <http://purl.org/dc/terms/>

    SELECT DISTINCT ?wpId
    WHERE {
      ?pathway a wp:Pathway ;
               wp:organism <http://purl.obolibrary.org/obo/NCBITaxon_9606> ;
               dcterms:identifier ?wpId .
    }
    """

    url = "https://sparql.wikipathways.org/sparql"
    response = requests.get(url, params={
        'query': query,
        'format': 'json'
    })

    ids = set()
    for binding in response.json()['results']['bindings']:
        ids.add(binding['wpId']['value'])

    return ids


def get_wikidata_wikipathways_ids() -> Set[str]:
    """Get WikiPathways IDs present in Wikidata."""

    query = """
    SELECT DISTINCT ?wpId
    WHERE {
      ?pathway wdt:P2410 ?wpId .
    }
    """

    url = "https://query.wikidata.org/sparql"
    response = requests.get(url, params={
        'query': query,
        'format': 'json'
    })

    ids = set()
    for binding in response.json()['results']['bindings']:
        ids.add(binding['wpId']['value'])

    return ids


def assess_coverage() -> Dict[str, any]:
    """Assess coverage of pathway databases in Wikidata."""

    print("Fetching Reactome pathways...")
    reactome_all = get_reactome_human_pathways()

    print("Fetching Wikidata Reactome IDs...")
    wikidata_reactome = get_wikidata_reactome_ids()

    print("Fetching WikiPathways pathways...")
    wikipathways_all = get_wikipathways_human_pathways()

    print("Fetching Wikidata WikiPathways IDs...")
    wikidata_wikipathways = get_wikidata_wikipathways_ids()

    # Calculate coverage
    reactome_coverage = len(wikidata_reactome & reactome_all) / len(reactome_all) * 100
    wikipathways_coverage = len(wikidata_wikipathways & wikipathways_all) / len(wikipathways_all) * 100

    # Find gaps
    missing_reactome = reactome_all - wikidata_reactome
    missing_wikipathways = wikipathways_all - wikidata_wikipathways

    results = {
        'reactome': {
            'total_in_reactome': len(reactome_all),
            'in_wikidata': len(wikidata_reactome & reactome_all),
            'coverage_percent': round(reactome_coverage, 2),
            'missing_count': len(missing_reactome),
            'missing_sample': list(missing_reactome)[:20]
        },
        'wikipathways': {
            'total_in_wikipathways': len(wikipathways_all),
            'in_wikidata': len(wikidata_wikipathways & wikipathways_all),
            'coverage_percent': round(wikipathways_coverage, 2),
            'missing_count': len(missing_wikipathways),
            'missing_sample': list(missing_wikipathways)[:20]
        }
    }

    return results


def assess_cross_reference_completeness() -> Dict[str, any]:
    """Check completeness of cross-references within Wikidata."""

    query = """
    SELECT
      (COUNT(DISTINCT ?pathway) AS ?total)
      (COUNT(DISTINCT ?withReactome) AS ?hasReactome)
      (COUNT(DISTINCT ?withWikiPathways) AS ?hasWikiPathways)
      (COUNT(DISTINCT ?withGO) AS ?hasGO)
      (COUNT(DISTINCT ?withMultiple) AS ?hasMultiple)
    WHERE {
      ?pathway wdt:P31/wdt:P279* wd:Q4915012 .

      OPTIONAL {
        ?pathway wdt:P3937 ?reactome .
        BIND(?pathway AS ?withReactome)
      }
      OPTIONAL {
        ?pathway wdt:P2410 ?wp .
        BIND(?pathway AS ?withWikiPathways)
      }
      OPTIONAL {
        ?pathway wdt:P686 ?go .
        BIND(?pathway AS ?withGO)
      }
      OPTIONAL {
        ?pathway wdt:P3937 ?r ; wdt:P2410 ?w .
        BIND(?pathway AS ?withMultiple)
      }
    }
    """

    url = "https://query.wikidata.org/sparql"
    response = requests.get(url, params={
        'query': query,
        'format': 'json'
    })

    binding = response.json()['results']['bindings'][0]

    total = int(binding['total']['value'])

    return {
        'total_pathways': total,
        'with_reactome': int(binding['hasReactome']['value']),
        'with_wikipathways': int(binding['hasWikiPathways']['value']),
        'with_go': int(binding['hasGO']['value']),
        'with_multiple_refs': int(binding['hasMultiple']['value']),
        'reactome_percent': round(int(binding['hasReactome']['value']) / total * 100, 2),
        'wikipathways_percent': round(int(binding['hasWikiPathways']['value']) / total * 100, 2),
        'go_percent': round(int(binding['hasGO']['value']) / total * 100, 2),
    }


if __name__ == '__main__':
    print("Assessing pathway coverage in Wikidata...\n")

    coverage = assess_coverage()

    print("=== Reactome Coverage ===")
    print(f"  Total Reactome pathways: {coverage['reactome']['total_in_reactome']}")
    print(f"  Present in Wikidata: {coverage['reactome']['in_wikidata']}")
    print(f"  Coverage: {coverage['reactome']['coverage_percent']}%")
    print(f"  Missing: {coverage['reactome']['missing_count']}")

    print("\n=== WikiPathways Coverage ===")
    print(f"  Total WikiPathways: {coverage['wikipathways']['total_in_wikipathways']}")
    print(f"  Present in Wikidata: {coverage['wikipathways']['in_wikidata']}")
    print(f"  Coverage: {coverage['wikipathways']['coverage_percent']}%")
    print(f"  Missing: {coverage['wikipathways']['missing_count']}")

    print("\n=== Cross-Reference Completeness ===")
    xref = assess_cross_reference_completeness()
    print(f"  Total pathways: {xref['total_pathways']}")
    print(f"  With Reactome ID: {xref['with_reactome']} ({xref['reactome_percent']}%)")
    print(f"  With WikiPathways ID: {xref['with_wikipathways']} ({xref['wikipathways_percent']}%)")
    print(f"  With GO ID: {xref['with_go']} ({xref['go_percent']}%)")
    print(f"  With multiple references: {xref['with_multiple_refs']}")

    # Save results
    with open('pathway_coverage_assessment.json', 'w') as f:
        json.dump({
            'coverage': coverage,
            'cross_references': xref
        }, f, indent=2)

    print("\nResults saved to pathway_coverage_assessment.json")
```

### 8.3 SPARQL Query for Gap Identification

```sparql
# Find pathways missing external identifiers
SELECT ?pathway ?pathwayLabel
       (BOUND(?reactome) AS ?hasReactome)
       (BOUND(?wp) AS ?hasWikiPathways)
       (BOUND(?go) AS ?hasGO)
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .

  # Filter for human-related pathways
  {
    ?pathway wdt:P703 wd:Q15978631 .
  }
  UNION
  {
    ?pathway wdt:P527 ?component .
    ?component wdt:P703 wd:Q15978631 .
  }

  OPTIONAL { ?pathway wdt:P3937 ?reactome . }
  OPTIONAL { ?pathway wdt:P2410 ?wp . }
  OPTIONAL { ?pathway wdt:P686 ?go . }

  # Only show items missing at least one reference
  FILTER(!BOUND(?reactome) || !BOUND(?wp) || !BOUND(?go))

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel
LIMIT 500
```

---

## 9. Integration Strategy

### 9.1 Using Wikidata as Pathway Index

Wikidata serves as an excellent central index for pathways due to its:
- Comprehensive cross-references to specialized databases
- Standardized identifiers and structure
- Regular community updates
- Open licensing (CC0)

#### Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                     Unified Pathway Database                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐      │
│  │   Wikidata   │    │   Reactome   │    │ WikiPathways │      │
│  │    Index     │───▶│   Details    │    │   Details    │      │
│  │              │    │              │    │              │      │
│  │  - Q-IDs     │    │  - Reactions │    │  - GPML      │      │
│  │  - Labels    │    │  - Proteins  │    │  - Diagrams  │      │
│  │  - XRefs     │    │  - Evidence  │    │  - Curation  │      │
│  └──────┬───────┘    └──────────────┘    └──────────────┘      │
│         │                                                        │
│         │            ┌──────────────┐    ┌──────────────┐      │
│         └───────────▶│    KEGG      │    │     GO       │      │
│                      │   Details    │    │   Details    │      │
│                      │              │    │              │      │
│                      │  - Maps      │    │  - Ontology  │      │
│                      │  - KO        │    │  - Evidence  │      │
│                      └──────────────┘    └──────────────┘      │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

### 9.2 Fetching Details from Specialized Databases

```python
#!/usr/bin/env python3
"""
Unified pathway data fetcher using Wikidata as index.
"""

import requests
import json
from typing import Dict, Any, Optional
from dataclasses import dataclass
from functools import lru_cache


@dataclass
class PathwayInfo:
    """Unified pathway information."""
    wikidata_id: str
    name: str
    description: str
    reactome_id: Optional[str] = None
    wikipathways_id: Optional[str] = None
    kegg_id: Optional[str] = None
    go_id: Optional[str] = None
    genes: list = None
    compounds: list = None
    reactions: list = None


class WikidataPathwayIndex:
    """Use Wikidata as central pathway index."""

    WIKIDATA_SPARQL = "https://query.wikidata.org/sparql"

    def __init__(self):
        self.cache = {}

    def search_pathways(self, query: str, limit: int = 20) -> list:
        """Search for pathways by name."""

        sparql = f"""
        SELECT ?pathway ?pathwayLabel ?reactome ?wikipathways ?go
        WHERE {{
          ?pathway wdt:P31/wdt:P279* wd:Q4915012 .
          ?pathway rdfs:label ?label .
          FILTER(CONTAINS(LCASE(?label), LCASE("{query}")))
          FILTER(LANG(?label) = "en")

          OPTIONAL {{ ?pathway wdt:P3937 ?reactome . }}
          OPTIONAL {{ ?pathway wdt:P2410 ?wikipathways . }}
          OPTIONAL {{ ?pathway wdt:P686 ?go . }}

          SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
        }}
        LIMIT {limit}
        """

        response = requests.get(
            self.WIKIDATA_SPARQL,
            params={'query': sparql, 'format': 'json'}
        )

        results = []
        for binding in response.json()['results']['bindings']:
            results.append({
                'wikidata_id': binding['pathway']['value'].split('/')[-1],
                'name': binding['pathwayLabel']['value'],
                'reactome_id': binding.get('reactome', {}).get('value'),
                'wikipathways_id': binding.get('wikipathways', {}).get('value'),
                'go_id': binding.get('go', {}).get('value'),
            })

        return results

    def get_pathway_xrefs(self, wikidata_id: str) -> Dict[str, str]:
        """Get all cross-references for a pathway."""

        sparql = f"""
        SELECT ?reactome ?wikipathways ?go ?kegg
        WHERE {{
          wd:{wikidata_id} wdt:P3937 ?reactome .
          OPTIONAL {{ wd:{wikidata_id} wdt:P2410 ?wikipathways . }}
          OPTIONAL {{ wd:{wikidata_id} wdt:P686 ?go . }}
          OPTIONAL {{
            wd:{wikidata_id} wdt:P2888 ?keggUri .
            FILTER(CONTAINS(STR(?keggUri), "kegg"))
            BIND(STR(?keggUri) AS ?kegg)
          }}
        }}
        """

        response = requests.get(
            self.WIKIDATA_SPARQL,
            params={'query': sparql, 'format': 'json'}
        )

        bindings = response.json()['results']['bindings']
        if not bindings:
            return {}

        binding = bindings[0]
        return {
            'reactome': binding.get('reactome', {}).get('value'),
            'wikipathways': binding.get('wikipathways', {}).get('value'),
            'go': binding.get('go', {}).get('value'),
            'kegg': binding.get('kegg', {}).get('value'),
        }


class ReactomeClient:
    """Fetch detailed pathway data from Reactome."""

    BASE_URL = "https://reactome.org/ContentService"

    @lru_cache(maxsize=1000)
    def get_pathway(self, reactome_id: str) -> Dict[str, Any]:
        """Get pathway details from Reactome."""

        url = f"{self.BASE_URL}/data/query/{reactome_id}"
        response = requests.get(url)

        if response.status_code != 200:
            return {}

        return response.json()

    def get_pathway_participants(self, reactome_id: str) -> Dict[str, list]:
        """Get genes and compounds participating in pathway."""

        url = f"{self.BASE_URL}/data/participants/{reactome_id}"
        response = requests.get(url)

        if response.status_code != 200:
            return {'genes': [], 'compounds': []}

        data = response.json()

        genes = []
        compounds = []

        for participant in data:
            if participant.get('schemaClass') == 'EntityWithAccessionedSequence':
                # This is a protein/gene
                ref_entity = participant.get('referenceEntity', {})
                if ref_entity.get('geneName'):
                    genes.append({
                        'symbol': ref_entity['geneName'][0] if ref_entity['geneName'] else '',
                        'uniprot': ref_entity.get('identifier', ''),
                        'name': participant.get('displayName', '')
                    })
            elif participant.get('schemaClass') == 'SimpleEntity':
                # This is a compound
                ref_entity = participant.get('referenceEntity', {})
                compounds.append({
                    'name': participant.get('displayName', ''),
                    'chebi': ref_entity.get('identifier', ''),
                })

        return {'genes': genes, 'compounds': compounds}


class WikiPathwaysClient:
    """Fetch detailed pathway data from WikiPathways."""

    SPARQL_ENDPOINT = "https://sparql.wikipathways.org/sparql"

    def get_pathway_genes(self, wp_id: str) -> list:
        """Get genes in a WikiPathways pathway."""

        query = f"""
        PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
        PREFIX dcterms: <http://purl.org/dc/terms/>

        SELECT DISTINCT ?gene ?label ?entrez ?ensembl
        WHERE {{
          ?pathway dcterms:identifier "{wp_id}" .
          ?gene a wp:GeneProduct ;
                dcterms:isPartOf ?pathway ;
                rdfs:label ?label .

          OPTIONAL {{ ?gene wp:bdbEntrezGene ?entrez . }}
          OPTIONAL {{ ?gene wp:bdbEnsembl ?ensembl . }}
        }}
        """

        response = requests.get(
            self.SPARQL_ENDPOINT,
            params={'query': query, 'format': 'json'}
        )

        genes = []
        for binding in response.json()['results']['bindings']:
            genes.append({
                'label': binding['label']['value'],
                'entrez': binding.get('entrez', {}).get('value', ''),
                'ensembl': binding.get('ensembl', {}).get('value', ''),
            })

        return genes

    def get_pathway_metabolites(self, wp_id: str) -> list:
        """Get metabolites in a WikiPathways pathway."""

        query = f"""
        PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
        PREFIX dcterms: <http://purl.org/dc/terms/>

        SELECT DISTINCT ?metabolite ?label ?chebi ?hmdb
        WHERE {{
          ?pathway dcterms:identifier "{wp_id}" .
          ?metabolite a wp:Metabolite ;
                      dcterms:isPartOf ?pathway ;
                      rdfs:label ?label .

          OPTIONAL {{ ?metabolite wp:bdbChEBI ?chebi . }}
          OPTIONAL {{ ?metabolite wp:bdbHmdb ?hmdb . }}
        }}
        """

        response = requests.get(
            self.SPARQL_ENDPOINT,
            params={'query': query, 'format': 'json'}
        )

        metabolites = []
        for binding in response.json()['results']['bindings']:
            metabolites.append({
                'label': binding['label']['value'],
                'chebi': binding.get('chebi', {}).get('value', ''),
                'hmdb': binding.get('hmdb', {}).get('value', ''),
            })

        return metabolites


class UnifiedPathwayDatabase:
    """Unified pathway database using Wikidata as index."""

    def __init__(self):
        self.wikidata = WikidataPathwayIndex()
        self.reactome = ReactomeClient()
        self.wikipathways = WikiPathwaysClient()

    def get_pathway(self, wikidata_id: str) -> PathwayInfo:
        """Get comprehensive pathway information from all sources."""

        # Get cross-references from Wikidata
        xrefs = self.wikidata.get_pathway_xrefs(wikidata_id)

        # Initialize pathway info
        pathway = PathwayInfo(
            wikidata_id=wikidata_id,
            name='',
            description='',
            reactome_id=xrefs.get('reactome'),
            wikipathways_id=xrefs.get('wikipathways'),
            kegg_id=xrefs.get('kegg'),
            go_id=xrefs.get('go'),
            genes=[],
            compounds=[],
            reactions=[]
        )

        # Fetch from Reactome if available
        if pathway.reactome_id:
            reactome_data = self.reactome.get_pathway(pathway.reactome_id)
            if reactome_data:
                pathway.name = reactome_data.get('displayName', '')
                pathway.description = reactome_data.get('summation', [{}])[0].get('text', '')

                participants = self.reactome.get_pathway_participants(pathway.reactome_id)
                pathway.genes.extend(participants['genes'])
                pathway.compounds.extend(participants['compounds'])

        # Fetch from WikiPathways if available
        if pathway.wikipathways_id:
            wp_genes = self.wikipathways.get_pathway_genes(pathway.wikipathways_id)
            wp_metabolites = self.wikipathways.get_pathway_metabolites(pathway.wikipathways_id)

            # Merge with existing (deduplicate)
            existing_genes = {g.get('symbol', g.get('label', '')) for g in pathway.genes}
            for gene in wp_genes:
                if gene['label'] not in existing_genes:
                    pathway.genes.append(gene)

            pathway.compounds.extend(wp_metabolites)

        return pathway

    def search(self, query: str) -> list:
        """Search for pathways across all sources."""
        return self.wikidata.search_pathways(query)

    def export_pathway(self, wikidata_id: str, format: str = 'json') -> str:
        """Export pathway data in specified format."""

        pathway = self.get_pathway(wikidata_id)

        if format == 'json':
            return json.dumps({
                'wikidata_id': pathway.wikidata_id,
                'name': pathway.name,
                'description': pathway.description,
                'identifiers': {
                    'reactome': pathway.reactome_id,
                    'wikipathways': pathway.wikipathways_id,
                    'kegg': pathway.kegg_id,
                    'go': pathway.go_id,
                },
                'genes': pathway.genes,
                'compounds': pathway.compounds,
            }, indent=2)

        elif format == 'tsv':
            lines = [
                f"# Pathway: {pathway.name}",
                f"# Wikidata: {pathway.wikidata_id}",
                f"# Reactome: {pathway.reactome_id}",
                f"# WikiPathways: {pathway.wikipathways_id}",
                "",
                "## Genes",
                "symbol\tuniprot\tname"
            ]
            for gene in pathway.genes:
                lines.append(f"{gene.get('symbol', '')}\t{gene.get('uniprot', '')}\t{gene.get('name', '')}")

            lines.extend(["", "## Compounds", "name\tchebi"])
            for compound in pathway.compounds:
                lines.append(f"{compound.get('name', '')}\t{compound.get('chebi', '')}")

            return '\n'.join(lines)

        return ''


# Example usage
if __name__ == '__main__':
    db = UnifiedPathwayDatabase()

    # Search for pathways
    print("Searching for 'glycolysis'...")
    results = db.search('glycolysis')

    for r in results[:5]:
        print(f"  {r['wikidata_id']}: {r['name']}")
        print(f"    Reactome: {r['reactome_id']}")
        print(f"    WikiPathways: {r['wikipathways_id']}")

    # Get detailed pathway info
    if results:
        print(f"\nGetting details for {results[0]['wikidata_id']}...")
        pathway = db.get_pathway(results[0]['wikidata_id'])

        print(f"Name: {pathway.name}")
        print(f"Genes: {len(pathway.genes)}")
        print(f"Compounds: {len(pathway.compounds)}")
```

### 9.3 Building Unified Pathway Database

```python
#!/usr/bin/env python3
"""
Build a unified pathway database from multiple sources.
"""

import sqlite3
import json
from pathlib import Path
from typing import Dict, List, Any


def create_database(db_path: str):
    """Create unified pathway database schema."""

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Pathways table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pathways (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            wikidata_id TEXT UNIQUE NOT NULL,
            name TEXT NOT NULL,
            description TEXT,
            pathway_type TEXT,
            organism TEXT DEFAULT 'Homo sapiens',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)

    # External identifiers
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pathway_xrefs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pathway_id INTEGER NOT NULL,
            database TEXT NOT NULL,
            identifier TEXT NOT NULL,
            url TEXT,
            FOREIGN KEY (pathway_id) REFERENCES pathways(id),
            UNIQUE(pathway_id, database, identifier)
        )
    """)

    # Genes table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS genes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            wikidata_id TEXT UNIQUE,
            symbol TEXT NOT NULL,
            name TEXT,
            entrez_id TEXT,
            ensembl_id TEXT,
            uniprot_id TEXT
        )
    """)

    # Pathway-gene relationships
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pathway_genes (
            pathway_id INTEGER NOT NULL,
            gene_id INTEGER NOT NULL,
            evidence_source TEXT,
            evidence_code TEXT,
            FOREIGN KEY (pathway_id) REFERENCES pathways(id),
            FOREIGN KEY (gene_id) REFERENCES genes(id),
            PRIMARY KEY (pathway_id, gene_id)
        )
    """)

    # Compounds table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS compounds (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            wikidata_id TEXT UNIQUE,
            name TEXT NOT NULL,
            chebi_id TEXT,
            pubchem_id TEXT,
            inchi TEXT,
            smiles TEXT
        )
    """)

    # Pathway-compound relationships
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pathway_compounds (
            pathway_id INTEGER NOT NULL,
            compound_id INTEGER NOT NULL,
            role TEXT,
            FOREIGN KEY (pathway_id) REFERENCES pathways(id),
            FOREIGN KEY (compound_id) REFERENCES compounds(id),
            PRIMARY KEY (pathway_id, compound_id)
        )
    """)

    # Pathway hierarchy
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pathway_hierarchy (
            parent_id INTEGER NOT NULL,
            child_id INTEGER NOT NULL,
            relationship_type TEXT DEFAULT 'part_of',
            FOREIGN KEY (parent_id) REFERENCES pathways(id),
            FOREIGN KEY (child_id) REFERENCES pathways(id),
            PRIMARY KEY (parent_id, child_id)
        )
    """)

    # Indexes
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_pathways_wikidata ON pathways(wikidata_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_xrefs_database ON pathway_xrefs(database, identifier)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(symbol)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_genes_entrez ON genes(entrez_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_compounds_chebi ON compounds(chebi_id)")

    conn.commit()
    return conn


def import_wikidata_pathways(conn: sqlite3.Connection, pathways_file: str):
    """Import pathways from Wikidata extraction."""

    cursor = conn.cursor()

    with open(pathways_file) as f:
        for line in f:
            item = json.loads(line)

            wikidata_id = item['qid']
            name = item['labels'].get('en', wikidata_id)
            description = item['descriptions'].get('en', '')

            # Insert pathway
            cursor.execute("""
                INSERT OR IGNORE INTO pathways (wikidata_id, name, description)
                VALUES (?, ?, ?)
            """, (wikidata_id, name, description))

            pathway_id = cursor.lastrowid or cursor.execute(
                "SELECT id FROM pathways WHERE wikidata_id = ?", (wikidata_id,)
            ).fetchone()[0]

            # Insert cross-references
            claims = item.get('claims', {})

            if 'P3937' in claims:
                for reactome_id in claims['P3937']:
                    cursor.execute("""
                        INSERT OR IGNORE INTO pathway_xrefs (pathway_id, database, identifier, url)
                        VALUES (?, 'reactome', ?, ?)
                    """, (pathway_id, reactome_id, f"https://reactome.org/content/detail/{reactome_id}"))

            if 'P2410' in claims:
                for wp_id in claims['P2410']:
                    cursor.execute("""
                        INSERT OR IGNORE INTO pathway_xrefs (pathway_id, database, identifier, url)
                        VALUES (?, 'wikipathways', ?, ?)
                    """, (pathway_id, wp_id, f"https://www.wikipathways.org/pathways/{wp_id}"))

            if 'P686' in claims:
                for go_id in claims['P686']:
                    cursor.execute("""
                        INSERT OR IGNORE INTO pathway_xrefs (pathway_id, database, identifier, url)
                        VALUES (?, 'go', ?, ?)
                    """, (pathway_id, go_id, f"http://amigo.geneontology.org/amigo/term/{go_id}"))

    conn.commit()


def import_pathway_genes(conn: sqlite3.Connection, genes_file: str, relationships_file: str):
    """Import genes and their pathway relationships."""

    cursor = conn.cursor()

    # Import genes
    with open(genes_file) as f:
        for line in f:
            item = json.loads(line)

            wikidata_id = item['qid']
            claims = item.get('claims', {})

            symbol = claims.get('P353', [''])[0]
            if not symbol:
                continue

            name = item['labels'].get('en', symbol)
            entrez = claims.get('P351', [''])[0]
            ensembl = claims.get('P594', [''])[0]
            uniprot = claims.get('P352', [''])[0]

            cursor.execute("""
                INSERT OR IGNORE INTO genes (wikidata_id, symbol, name, entrez_id, ensembl_id, uniprot_id)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (wikidata_id, symbol, name, entrez, ensembl, uniprot))

    # Import pathway-gene relationships
    with open(relationships_file) as f:
        for line in f:
            rel = json.loads(line)

            if rel['property'] == 'P527':  # has part
                pathway_qid = rel['source']
                gene_qid = rel['target']

                cursor.execute("""
                    INSERT OR IGNORE INTO pathway_genes (pathway_id, gene_id, evidence_source)
                    SELECT p.id, g.id, 'wikidata_haspart'
                    FROM pathways p, genes g
                    WHERE p.wikidata_id = ? AND g.wikidata_id = ?
                """, (pathway_qid, gene_qid))

    conn.commit()


def query_pathway_by_gene(conn: sqlite3.Connection, gene_symbol: str) -> List[Dict]:
    """Find all pathways containing a gene."""

    cursor = conn.cursor()
    cursor.execute("""
        SELECT DISTINCT
            p.wikidata_id,
            p.name,
            p.description,
            GROUP_CONCAT(DISTINCT x.database || ':' || x.identifier) as xrefs
        FROM pathways p
        JOIN pathway_genes pg ON p.id = pg.pathway_id
        JOIN genes g ON pg.gene_id = g.id
        LEFT JOIN pathway_xrefs x ON p.id = x.pathway_id
        WHERE UPPER(g.symbol) = UPPER(?)
        GROUP BY p.id
    """, (gene_symbol,))

    results = []
    for row in cursor.fetchall():
        results.append({
            'wikidata_id': row[0],
            'name': row[1],
            'description': row[2],
            'cross_references': row[3].split(',') if row[3] else []
        })

    return results


def query_genes_in_pathway(conn: sqlite3.Connection, pathway_id: str) -> List[Dict]:
    """Get all genes in a pathway."""

    cursor = conn.cursor()
    cursor.execute("""
        SELECT
            g.wikidata_id,
            g.symbol,
            g.name,
            g.entrez_id,
            g.ensembl_id,
            g.uniprot_id
        FROM genes g
        JOIN pathway_genes pg ON g.id = pg.gene_id
        JOIN pathways p ON pg.pathway_id = p.id
        WHERE p.wikidata_id = ?
        ORDER BY g.symbol
    """, (pathway_id,))

    return [
        {
            'wikidata_id': row[0],
            'symbol': row[1],
            'name': row[2],
            'entrez_id': row[3],
            'ensembl_id': row[4],
            'uniprot_id': row[5]
        }
        for row in cursor.fetchall()
    ]


if __name__ == '__main__':
    import sys

    db_path = 'unified_pathways.db'

    print(f"Creating database: {db_path}")
    conn = create_database(db_path)

    if len(sys.argv) > 1:
        data_dir = sys.argv[1]

        print("Importing Wikidata pathways...")
        import_wikidata_pathways(conn, f"{data_dir}/pathways.jsonl")

        print("Importing genes and relationships...")
        import_pathway_genes(
            conn,
            f"{data_dir}/genes_with_pathways.jsonl",
            f"{data_dir}/pathway_relationships.jsonl"
        )

    # Example queries
    print("\nExample: Pathways containing TP53")
    for pathway in query_pathway_by_gene(conn, 'TP53'):
        print(f"  {pathway['name']} ({pathway['wikidata_id']})")

    conn.close()
```

---

## 10. Best Practices and Recommendations

### 10.1 Query Optimization

1. **Use LIMIT**: Always include LIMIT clauses to prevent timeouts
2. **Avoid SELECT ***: Specify only needed variables
3. **Use property paths sparingly**: `wdt:P31/wdt:P279*` can be expensive
4. **Cache results**: Store frequently-used query results locally

### 10.2 Data Quality Considerations

1. **Verify cross-references**: Not all Wikidata cross-references are current
2. **Check for duplicates**: Same pathway may exist under different Q-IDs
3. **Validate labels**: Some items lack English labels
4. **Monitor updates**: Use Recent Changes API for incremental updates

### 10.3 Integration Recommendations

1. **Use Wikidata as index**: Don't rely on it for complete pathway details
2. **Fetch details from sources**: Use Reactome/WikiPathways APIs for current data
3. **Build local cache**: Store fetched data to reduce API calls
4. **Implement fallbacks**: Handle missing cross-references gracefully

### 10.4 Licensing and Attribution

- **Wikidata**: CC0 (public domain dedication)
- **Reactome**: CC BY 4.0
- **WikiPathways**: CC BY 3.0
- **KEGG**: Requires license for commercial use
- **Gene Ontology**: CC BY 4.0

Always include appropriate attribution when redistributing data.

---

## Appendix A: Property Reference

| Property | Label | Description |
|----------|-------|-------------|
| P31 | instance of | Class membership |
| P279 | subclass of | Class hierarchy |
| P361 | part of | Parent pathway |
| P527 | has part | Components |
| P682 | biological process | GO BP annotation |
| P680 | molecular function | GO MF annotation |
| P681 | cellular component | GO CC annotation |
| P686 | Gene Ontology ID | GO term ID |
| P702 | encoded by | Gene encoding protein |
| P703 | found in taxon | Species |
| P2410 | WikiPathways ID | WikiPathways cross-ref |
| P2868 | subject has role | Functional role |
| P2888 | exact match | KEGG and other URIs |
| P3937 | Reactome pathway ID | Reactome cross-ref |

## Appendix B: SPARQL Endpoint URLs

| Endpoint | URL | Rate Limits |
|----------|-----|-------------|
| Wikidata | https://query.wikidata.org/sparql | 1 req/sec recommended |
| WikiPathways | https://sparql.wikipathways.org/sparql | No official limit |
| UniProt | https://sparql.uniprot.org/sparql | 1 req/sec recommended |

## Appendix C: Related Resources

- [Wikidata Query Service](https://query.wikidata.org/)
- [WikiPathways SPARQL Tutorial](https://www.wikipathways.org/index.php/Help:WikiPathways_Sparql_queries)
- [Reactome Content Service API](https://reactome.org/ContentService/)
- [Gene Ontology SPARQL](http://sparql.geneontology.org/)
