---
id: download-wikidata
title: "Wikidata Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# Wikidata Download Instructions

## Quick Start

```bash
# Download latest JSON dump (large!)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2
```

## Prerequisites

- **wget** or **curl** for downloads
- **bzip2** for decompression
- **jq** for JSON processing
- 100GB-2TB disk space (full dump is ~100GB compressed, ~1TB uncompressed)

## No Registration Required

Wikidata is freely available under CC0 (public domain).

## Download Methods

### Method 1: Full Dump (JSON)

```bash
# Complete Wikidata dump (JSON)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2

# Truthy statements only (smaller)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-truthy.nt.bz2
```

### Method 2: RDF/TTL Dumps

```bash
# N-Triples format (all statements)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.nt.bz2

# Truthy triples only
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-truthy.nt.bz2

# Lexemes (words/language)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-lexemes.json.bz2
```

### Method 3: SPARQL Queries (Recommended for Specific Data)

```bash
# Query biological databases cross-references
curl -G https://query.wikidata.org/sparql \
  --data-urlencode 'query=
    SELECT ?item ?itemLabel ?uniprot ?ensembl WHERE {
      ?item wdt:P31 wd:Q7187.  # instance of gene
      ?item wdt:P703 wd:Q15978631.  # found in Homo sapiens
      OPTIONAL { ?item wdt:P352 ?uniprot. }
      OPTIONAL { ?item wdt:P594 ?ensembl. }
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    LIMIT 1000
  ' \
  -H "Accept: text/csv" \
  -o human_genes.csv

# Query diseases with OMIM IDs
curl -G https://query.wikidata.org/sparql \
  --data-urlencode 'query=
    SELECT ?disease ?diseaseLabel ?omim WHERE {
      ?disease wdt:P31 wd:Q12136.  # instance of disease
      ?disease wdt:P492 ?omim.  # OMIM ID
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    LIMIT 5000
  ' \
  -H "Accept: text/csv" \
  -o diseases_omim.csv

# Query drugs with ChEMBL IDs
curl -G https://query.wikidata.org/sparql \
  --data-urlencode 'query=
    SELECT ?drug ?drugLabel ?chembl ?drugbank WHERE {
      ?drug wdt:P31 wd:Q12140.  # instance of medication
      OPTIONAL { ?drug wdt:P592 ?chembl. }
      OPTIONAL { ?drug wdt:P715 ?drugbank. }
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    LIMIT 5000
  ' \
  -H "Accept: text/csv" \
  -o drugs.csv
```

### Method 4: Wikidata API

```bash
# Get specific entity
curl "https://www.wikidata.org/wiki/Special:EntityData/Q7187.json" \
  -o gene_entity.json

# Search entities
curl "https://www.wikidata.org/w/api.php?action=wbsearchentities&search=BRCA1&language=en&format=json" \
  -o brca1_search.json

# Get entity with claims
curl "https://www.wikidata.org/w/api.php?action=wbgetentities&ids=Q17853272&format=json" \
  -o brca1_entity.json
```

### Method 5: WikidataIntegrator (Python)

```python
# pip install wikidataintegrator

from wikidataintegrator import wdi_core

# Get entity by QID
entity = wdi_core.WDItemEngine(wd_item_id='Q17853272')  # BRCA1

# Get all claims
for claim in entity.get_json_representation()['claims']:
    print(claim)

# Query with SPARQL
results = wdi_core.WDItemEngine.execute_sparql_query('''
    SELECT ?gene ?geneLabel ?uniprot WHERE {
      ?gene wdt:P31 wd:Q7187.
      ?gene wdt:P352 ?uniprot.
      FILTER(CONTAINS(?uniprot, "P38398"))
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
''')
```

### Method 6: Subset Dumps

```bash
# Download specific dated dump
DUMP_DATE="20240101"
wget "https://dumps.wikimedia.org/wikidatawiki/entities/${DUMP_DATE}/wikidata-${DUMP_DATE}-all.json.bz2"

# List available dumps
curl -s https://dumps.wikimedia.org/wikidatawiki/entities/ | grep -o '20[0-9]*'
```

## File Inventory

### Full Dumps

| File | Size (compressed) | Description |
|------|-------------------|-------------|
| latest-all.json.bz2 | ~90 GB | Complete JSON dump |
| latest-all.nt.bz2 | ~70 GB | All N-Triples |
| latest-truthy.nt.bz2 | ~30 GB | Truthy triples only |
| latest-lexemes.json.bz2 | ~1 GB | Lexeme data |

### Derived Subsets

| Data Type | Estimated Records |
|-----------|-------------------|
| Genes (Q7187) | ~2 million |
| Diseases (Q12136) | ~50,000 |
| Drugs (Q12140) | ~30,000 |
| Proteins (Q8054) | ~5 million |
| Chemical compounds | ~100 million |

## Post-Download Processing

```bash
# Stream process JSON dump (don't decompress fully!)
bzcat latest-all.json.bz2 | jq -c 'select(.type=="item") |
  select(.claims.P352?) |
  {qid: .id, label: .labels.en.value, uniprot: .claims.P352[0].mainsnak.datavalue.value}' \
  > proteins_uniprot.jsonl

# Extract specific entity type
bzcat latest-all.json.bz2 | python3 << 'EOF' > genes.jsonl
import sys
import json

for line in sys.stdin:
    line = line.strip()
    if line.startswith('[') or line.startswith(']'):
        continue
    if line.endswith(','):
        line = line[:-1]
    try:
        entity = json.loads(line)
        # Check if it's a gene (P31 = Q7187)
        claims = entity.get('claims', {})
        p31 = claims.get('P31', [])
        for claim in p31:
            if claim.get('mainsnak', {}).get('datavalue', {}).get('value', {}).get('id') == 'Q7187':
                print(json.dumps({
                    'qid': entity['id'],
                    'label': entity.get('labels', {}).get('en', {}).get('value'),
                    'claims': claims
                }))
                break
    except:
        pass
EOF

# Load into DuckDB
duckdb << 'EOF'
CREATE TABLE wikidata_genes AS
SELECT * FROM read_json_auto('genes.jsonl');
SELECT COUNT(*) FROM wikidata_genes;
EOF
```

## Verification

```bash
# Check dump structure
bzcat latest-all.json.bz2 | head -c 10000

# Count entities (streaming)
bzcat latest-all.json.bz2 | grep -c '"type":"item"'

# Verify SPARQL endpoint
curl -s "https://query.wikidata.org/sparql?query=SELECT%20(COUNT(*)%20as%20%3Fcount)%20WHERE%20%7B%20%3Fs%20%3Fp%20%3Fo%20%7D" \
  -H "Accept: application/json"
```

## Update Schedule

| Dump Type | Frequency |
|-----------|-----------|
| Full dumps | Weekly (Monday) |
| Incremental | Every few minutes |
| SPARQL | Real-time |

## Common Issues

- **Dump size**: Full dump is enormous; use SPARQL for specific data
- **JSON streaming**: Don't load entire dump in memory; stream process
- **Rate limits**: SPARQL endpoint has query limits; batch requests
- **Property IDs**: Learn common properties (P31=instance of, P352=UniProt, etc.)
- **Label language**: Always specify language for labels

## Key Properties for Bioinformatics

| Property | ID | Description |
|----------|-----|-------------|
| instance of | P31 | Entity type |
| UniProt ID | P352 | Protein ID |
| Ensembl gene ID | P594 | Gene ID |
| OMIM ID | P492 | Disease ID |
| ChEMBL ID | P592 | Compound ID |
| DrugBank ID | P715 | Drug ID |
| PubChem CID | P662 | Chemical ID |
| NCBI Gene ID | P351 | Gene ID |
| found in taxon | P703 | Organism |

## SPARQL Query Examples

```sparql
# Human proteins with structures
SELECT ?protein ?proteinLabel ?uniprot ?pdb WHERE {
  ?protein wdt:P31 wd:Q8054.  # instance of protein
  ?protein wdt:P703 wd:Q15978631.  # found in human
  ?protein wdt:P352 ?uniprot.
  OPTIONAL { ?protein wdt:P638 ?pdb. }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000

# Disease-gene associations
SELECT ?disease ?diseaseLabel ?gene ?geneLabel WHERE {
  ?disease wdt:P31 wd:Q12136.  # disease
  ?disease wdt:P2293 ?gene.  # genetic association
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

## Related Resources

- [UniProt](../../07.proteins.molecular.biology/7.1.protein.sequences.annotations/uniprot/) - Protein database
- [PubMed](../8.1.scientific.literature/pubmed/) - Literature
- [ChEMBL](../../02.compounds.molecules/2.2.pharmaceuticals/chembl/) - Bioactivity
