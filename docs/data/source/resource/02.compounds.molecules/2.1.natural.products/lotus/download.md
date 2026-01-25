---
id: download-lotus
title: "LOTUS Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# LOTUS Download Instructions

## Quick Start

```bash
# Download from Zenodo archive
wget https://zenodo.org/record/5794106/files/lotus.json.gz
gunzip lotus.json.gz

# Or query via Wikidata SPARQL
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=SELECT * WHERE { ?compound wdt:P703 ?taxon } LIMIT 10" \
  -H "Accept: application/json"
```

## Prerequisites

- **wget** or **curl** for downloads
- ~2 GB disk space for full dataset
- SPARQL knowledge for Wikidata queries
- JSON parser (jq recommended)

## No Registration Required

Data is CC0 (Public Domain) - completely unrestricted use.

## Download Methods

### Method 1: Zenodo Archive (Recommended for Bulk)

```bash
# Download complete LOTUS dataset
wget https://zenodo.org/record/5794106/files/lotus.json.gz
gunzip lotus.json.gz

# Also available in TSV format
wget https://zenodo.org/record/5794106/files/lotus.tsv.gz
gunzip lotus.tsv.gz

# Check latest version at:
# https://doi.org/10.5281/zenodo.5794106
```

### Method 2: Wikidata SPARQL (Recommended for Queries)

```bash
# Basic query: compounds found in taxon
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=
    SELECT ?compound ?compoundLabel ?taxon ?taxonLabel WHERE {
      ?compound wdt:P703 ?taxon .
      ?compound wdt:P235 ?inchikey .
      SERVICE wikibase:label { bd:serviceParam wikibase:language 'en' }
    } LIMIT 100" \
  -H "Accept: application/json" | jq

# Query compounds from specific organism
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=
    SELECT ?compound ?compoundLabel ?inchikey WHERE {
      ?compound wdt:P703 wd:Q158695 .  # Artemisia annua
      ?compound wdt:P235 ?inchikey .
      SERVICE wikibase:label { bd:serviceParam wikibase:language 'en' }
    }" \
  -H "Accept: application/json" | jq

# Get structure-organism-reference triplets
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=
    SELECT ?compound ?inchikey ?taxon ?taxonLabel ?reference WHERE {
      ?compound wdt:P235 ?inchikey .
      ?compound p:P703 ?stmt .
      ?stmt ps:P703 ?taxon .
      OPTIONAL { ?stmt prov:wasDerivedFrom/pr:P248 ?reference }
      SERVICE wikibase:label { bd:serviceParam wikibase:language 'en' }
    } LIMIT 1000" \
  -H "Accept: application/json" | jq
```

### Method 3: LNPN Web Interface

```bash
# Interactive search at LOTUS Natural Products Navigator
# https://lotus.naturalproducts.net

# Features:
# - Structure search
# - Organism browser
# - Export search results
```

### Method 4: Substructure Search (IDSM)

```bash
# Substructure search via IDSM/Sachem endpoint
# https://idsm.elixir-czech.cz/sparql/endpoint/wikidata

# Find compounds containing a specific scaffold
curl -G "https://idsm.elixir-czech.cz/sparql/endpoint/wikidata" \
  --data-urlencode "query=
    PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
    SELECT ?compound WHERE {
      ?compound sachem:substructureSearch [
        sachem:query 'c1ccc2c(c1)cccc2'
      ]
    } LIMIT 100" \
  -H "Accept: application/json"
```

## File Inventory

| File | Size | Description |
|------|------|-------------|
| lotus.json.gz | ~400 MB | Complete JSON export |
| lotus.tsv.gz | ~300 MB | Tab-separated format |
| lotus.sdf.gz | ~500 MB | Structure-Data File |

## Data Fields in Zenodo Export

| Field | Description |
|-------|-------------|
| structure_wikidata | Compound QID |
| structure_inchikey | InChI Key |
| structure_smiles | SMILES |
| organism_wikidata | Taxon QID |
| organism_name | Scientific name |
| organism_taxonomy_ncbi | NCBI Taxon ID |
| reference_wikidata | Reference QID |
| reference_doi | DOI if available |

## Post-Download Processing

```bash
# Preview JSON data
zcat lotus.json.gz | head -100 | jq '.[0]'

# Count total records
zcat lotus.json.gz | jq 'length'

# Extract unique organisms
zcat lotus.json.gz | jq -r '.[].organism_name' | sort | uniq | wc -l

# Extract unique compounds
zcat lotus.json.gz | jq -r '.[].structure_inchikey' | sort | uniq | wc -l

# Filter for specific organism
zcat lotus.json.gz | jq '[.[] | select(.organism_name | contains("Artemisia"))]' > artemisia_compounds.json

# Convert to CSV
zcat lotus.json.gz | jq -r '.[] | [.structure_inchikey, .organism_name, .reference_doi] | @csv' > lotus.csv

# Load into SQLite
sqlite3 lotus.db << 'EOF'
CREATE TABLE lotus (
  structure_qid TEXT,
  inchikey TEXT,
  smiles TEXT,
  organism_qid TEXT,
  organism_name TEXT,
  ncbi_taxon_id INTEGER,
  reference_qid TEXT,
  reference_doi TEXT
);
CREATE INDEX idx_inchikey ON lotus(inchikey);
CREATE INDEX idx_organism ON lotus(organism_name);
CREATE INDEX idx_taxon ON lotus(ncbi_taxon_id);
EOF

# Import from TSV
zcat lotus.tsv.gz | sqlite3 lotus.db '.mode tabs' '.import /dev/stdin lotus'
```

## SPARQL Query Examples

```sparql
# Count compounds per kingdom
SELECT ?kingdom (COUNT(DISTINCT ?compound) as ?count) WHERE {
  ?compound wdt:P703 ?taxon .
  ?taxon wdt:P105 wd:Q36732 ;  # taxonomic rank = species
         wdt:P171* ?kingdom_taxon .
  ?kingdom_taxon wdt:P105 wd:Q36732 ;  # kingdom rank
                  rdfs:label ?kingdom .
  FILTER(LANG(?kingdom) = "en")
}
GROUP BY ?kingdom
ORDER BY DESC(?count)

# Find all quinoline-containing natural products
PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
SELECT ?compound ?taxonLabel WHERE {
  SERVICE <https://idsm.elixir-czech.cz/sparql/endpoint/wikidata> {
    ?compound sachem:substructureSearch [
      sachem:query "c1ccc2ncccc2c1"  # quinoline core
    ]
  }
  ?compound wdt:P703 ?taxon .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
```

## Verification

```bash
# Check archive integrity
gzip -t lotus.json.gz && echo "OK"

# Verify JSON validity
zcat lotus.json.gz | jq 'length'

# Check expected fields
zcat lotus.json.gz | jq '.[0] | keys'

# Sample random entries
zcat lotus.json.gz | jq '.[1000:1010]'
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| LOTUS 2025 | 2025-11-20 | ~1.2 GB | Current |
| LOTUS 1.0 | 2021-02-28 | ~800 MB | Archived |

### Version Notes

LOTUS contains 750,000+ referenced structure-organism pairs:
- Data hosted on Wikidata with community curation
- Mirrored at lotus.naturalproducts.net (being phased out)
- Available under CC0 license for unrestricted use

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://query.wikidata.org/sparql` |
| Rate Limit | Subject to Wikidata limits |
| Auth Required | No |
| Documentation | https://www.wikidata.org/wiki/Wikidata:SPARQL_query_service |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Zenodo archive | Periodic snapshots |
| Wikidata | Continuous updates |
| LNPN interface | Daily sync |

## Common Issues

- **SPARQL timeout**: Large queries may timeout; use LIMIT and pagination
- **Wikidata updates**: Data is continuously updated; archive is snapshot
- **InChI Key matching**: Use for structure deduplication
- **Missing references**: Not all entries have DOI citations

## Integration with Wikidata

```bash
# Get additional compound properties from Wikidata
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=
    SELECT ?compound ?mass ?formula WHERE {
      VALUES ?compound { wd:Q312879 }  # Artemisinin
      OPTIONAL { ?compound wdt:P2067 ?mass }
      OPTIONAL { ?compound wdt:P274 ?formula }
    }" \
  -H "Accept: application/json" | jq
```

## License

CC0 (Public Domain) - No restrictions, no attribution required (but appreciated)
