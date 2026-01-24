---
id: guides-wikidata-wikipedia-integration
title: "Wikipedia & Wikidata Integration Guide"
type: guide
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [wikipedia, wikidata, api, integration, cross-references, guide]
---

**Parent:** [Wikidata Guides](./_index.md)

# Wikipedia & Wikidata Integration Guide

**Created:** January 2026
**Purpose:** Guide for integrating Wikipedia/Wikidata with medicinal plant, compound, and genetic data
**Status:** Active reference document

---

## Executive Summary

Wikipedia and Wikidata provide valuable structured data that can:

1. **Fill identifier gaps** - Map common names to scientific names to database IDs
2. **Cross-validate data** - Verify compound-plant relationships
3. **Provide metadata** - Get images, descriptions, common names in multiple languages
4. **Link to external databases** - Wikidata has extensive external ID mappings

**Key Advantage:** Wikidata is CC0 (public domain) and has a robust SPARQL endpoint for bulk queries.

---

## Part 1: Wikipedia API

### 1.1 Wikipedia REST API

**Base URL:** `https://en.wikipedia.org/api/rest_v1/`

#### Get Page Summary

```bash
# Get summary of a medicinal plant
curl "https://en.wikipedia.org/api/rest_v1/page/summary/Curcuma_longa"

# Response includes:
# - Title, description, extract (first paragraph)
# - Thumbnail image URL
# - Wikidata ID (Q-number)
# - Coordinates (if applicable)
```

**Example Response:**
```json
{
  "title": "Curcuma longa",
  "extract": "Curcuma longa, commonly known as turmeric, is a flowering plant...",
  "thumbnail": {
    "source": "https://upload.wikimedia.org/...",
    "width": 320,
    "height": 213
  },
  "wikibase_item": "Q42562",  // ← Wikidata ID
  "coordinates": {
    "lat": 10.0,
    "lon": 76.0
  }
}
```

#### Get Full Page Content

```bash
# Get full HTML content
curl "https://en.wikipedia.org/api/rest_v1/page/html/Curcuma_longa"

# Get mobile-optimized HTML
curl "https://en.wikipedia.org/api/rest_v1/page/mobile-html/Curcuma_longa"
```

---

### 1.2 Wikipedia MediaWiki API

**Base URL:** `https://en.wikipedia.org/w/api.php`

#### Search for Pages

```bash
# Search for a plant by common name
curl "https://en.wikipedia.org/w/api.php?action=query&list=search&srsearch=turmeric&format=json"

# Get suggestions (autocomplete)
curl "https://en.wikipedia.org/w/api.php?action=opensearch&search=turm&format=json"
```

#### Get Page Info and Wikidata ID

```bash
# Get Wikidata ID for a page
curl "https://en.wikipedia.org/w/api.php?action=query&titles=Curcuma_longa&prop=pageprops&format=json"

# Response:
{
  "query": {
    "pages": {
      "123456": {
        "pageid": 123456,
        "title": "Curcuma longa",
        "pageprops": {
          "wikibase_item": "Q42562"  // ← Wikidata ID
        }
      }
    }
  }
}
```

#### Get Infobox Data

```bash
# Get page properties (including infobox data)
curl "https://en.wikipedia.org/w/api.php?action=query&titles=Curcuma_longa&prop=revisions&rvprop=content&format=json"

# Parse the wikitext to extract infobox fields
# (This requires parsing MediaWiki markup)
```

---

## Part 2: Wikidata Query Service (SPARQL)

### 2.1 SPARQL Endpoint

**URL:** `https://query.wikidata.org/sparql`

**Query Interface:** https://query.wikidata.org/

#### Basic Query Structure

```sparql
SELECT ?item ?itemLabel WHERE {
  # Your query here

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

---

### 2.2 Example Queries for Traditional Medicine

#### Query 1: Get All Medicinal Plants with PubChem IDs

```sparql
SELECT ?plant ?plantLabel ?pubchemCID WHERE {
  ?plant wdt:P31 wd:Q756.           # instance of "plant"
  ?plant wdt:P2067 ?mass.           # has mass (indicates chemical data)
  ?plant wdt:P662 ?pubchemCID.      # has PubChem Compound ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

**Execute via API:**
```bash
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=SELECT ?plant ?plantLabel ?pubchemCID WHERE { ?plant wdt:P31 wd:Q756. ?plant wdt:P662 ?pubchemCID. SERVICE wikibase:label { bd:serviceParam wikibase:language \"en\". } } LIMIT 1000" \
  -H "Accept: application/sparql-results+json"
```

---

#### Query 2: Get Compounds from a Specific Plant

```sparql
# Get all compounds found in Turmeric (Curcuma longa)
SELECT ?compound ?compoundLabel ?pubchemCID ?chemicalFormula WHERE {
  wd:Q42562 wdt:P527 ?compound.       # Curcuma longa "has part" compound
  OPTIONAL { ?compound wdt:P662 ?pubchemCID. }
  OPTIONAL { ?compound wdt:P274 ?chemicalFormula. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

**Result includes:** Curcumin (Q421789), demethoxycurcumin, bisdemethoxycurcumin, etc.

---

#### Query 3: Map Common Names to Scientific Names

```sparql
# Find plants by common name "turmeric"
SELECT ?plant ?plantLabel ?scientificName ?commonName WHERE {
  ?plant wdt:P31 wd:Q756.             # instance of plant
  ?plant wdt:P225 ?scientificName.    # taxon name
  ?plant rdfs:label ?commonName.

  FILTER(CONTAINS(LCASE(?commonName), "turmeric"))
  FILTER(LANG(?commonName) = "en")

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

#### Query 4: Get All External Database IDs for a Compound

```sparql
# Get all external IDs for Curcumin (Q421789)
SELECT ?property ?propertyLabel ?value WHERE {
  wd:Q421789 ?prop ?value.
  ?property wikibase:directClaim ?prop.
  ?property wikibase:propertyType wikibase:ExternalId.

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

**Returns:**
- PubChem CID (P662)
- ChEMBL ID (P592)
- ChEBI ID (P683)
- KEGG ID (P665)
- CAS Registry Number (P231)
- DrugBank ID (P715)
- etc.

---

#### Query 5: Get Targets/Proteins for Compounds

```sparql
# Get protein targets for compounds from medicinal plants
SELECT ?compound ?compoundLabel ?protein ?proteinLabel ?uniprotID WHERE {
  ?compound wdt:P31 wd:Q11173.        # instance of "chemical compound"
  ?compound wdt:P703 wd:Q756.         # found in taxon "plant"
  ?compound wdt:P129 ?protein.        # physically interacts with protein

  OPTIONAL { ?protein wdt:P352 ?uniprotID. }  # UniProt ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

---

#### Query 6: Get Pathways for Compounds

```sparql
# Get biochemical pathways involving compounds
SELECT ?compound ?compoundLabel ?pathway ?pathwayLabel WHERE {
  ?compound wdt:P31 wd:Q11173.        # instance of chemical compound
  ?compound wdt:P2868 ?pathway.       # subject has role (pathway)

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

### 2.3 Bulk Data Download (Wikidata Dumps)

For large-scale integration, download Wikidata dumps:

**URL:** https://dumps.wikimedia.org/wikidatawiki/entities/

```bash
# Download latest JSON dump (warning: very large ~100GB compressed)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.gz

# Or download specific entity dumps
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2
```

**Processing Tools:**
- `wikidata-taxonomy` - Python tool for parsing dumps
- `wikidata-cli` - Command-line tool for querying
- `wikidata-sdk` - JavaScript library

---

## Part 3: Key Wikidata Properties for Our Use Case

### 3.1 Plant/Organism Properties

| Property | Description | Example |
|----------|-------------|---------|
| **P31** | instance of | Q756 (plant) |
| **P225** | taxon name | "Curcuma longa" |
| **P105** | taxon rank | Q7432 (species) |
| **P171** | parent taxon | Q34687 (Curcuma genus) |
| **P18** | image | Commons filename |
| **P527** | has part (compound) | Q421789 (curcumin) |
| **P2275** | World Flora Online ID | wfo-0000123456 |

### 3.2 Chemical Compound Properties

| Property | Description | Example |
|----------|-------------|---------|
| **P31** | instance of | Q11173 (chemical compound) |
| **P662** | PubChem CID | 969516 |
| **P683** | ChEBI ID | CHEBI:3962 |
| **P592** | ChEMBL ID | CHEMBL308 |
| **P231** | CAS Registry Number | 458-37-7 |
| **P274** | chemical formula | C21H20O6 |
| **P2017** | isomeric SMILES | ... |
| **P703** | found in taxon | Q42562 (Curcuma longa) |
| **P129** | physically interacts with | Q protein |

### 3.3 Protein/Gene Properties

| Property | Description | Example |
|----------|-------------|---------|
| **P31** | instance of | Q8054 (protein) |
| **P352** | UniProt ID | P04637 |
| **P353** | HGNC gene symbol | TP53 |
| **P594** | Ensembl gene ID | ENSG00000141510 |
| **P354** | HGNC ID | 11998 |
| **P638** | PDB structure ID | 1TUP |

### 3.4 Pathway Properties

| Property | Description | Example |
|----------|-------------|---------|
| **P31** | instance of | Q4915012 (biological pathway) |
| **P2888** | exact match (external URI) | Reactome URL |
| **P665** | KEGG pathway ID | hsa00010 |
| **P2868** | subject has role | Q... |

---

## Part 4: Integration Pipeline Examples

### 4.1 Example 1: Plant Name → Compounds → Targets

```python
import requests

# Step 1: Get Wikidata ID from common name via Wikipedia
def get_wikidata_id(plant_name):
    url = f"https://en.wikipedia.org/api/rest_v1/page/summary/{plant_name}"
    response = requests.get(url)
    data = response.json()
    return data.get('wikibase_item')  # Returns Q-number

# Step 2: Query Wikidata for compounds in the plant
def get_compounds(wikidata_id):
    query = f"""
    SELECT ?compound ?compoundLabel ?pubchemCID WHERE {{
      wd:{wikidata_id} wdt:P527 ?compound.
      OPTIONAL {{ ?compound wdt:P662 ?pubchemCID. }}
      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
    }}
    """

    url = "https://query.wikidata.org/sparql"
    response = requests.get(url, params={'query': query, 'format': 'json'})
    return response.json()['results']['bindings']

# Step 3: For each compound, get protein targets
def get_targets(compound_wikidata_id):
    query = f"""
    SELECT ?protein ?proteinLabel ?uniprotID WHERE {{
      wd:{compound_wikidata_id} wdt:P129 ?protein.
      OPTIONAL {{ ?protein wdt:P352 ?uniprotID. }}
      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
    }}
    """

    url = "https://query.wikidata.org/sparql"
    response = requests.get(url, params={'query': query, 'format': 'json'})
    return response.json()['results']['bindings']

# Usage
plant_name = "Curcuma_longa"
wikidata_id = get_wikidata_id(plant_name)
print(f"Wikidata ID: {wikidata_id}")

compounds = get_compounds(wikidata_id)
for compound in compounds:
    print(f"Compound: {compound['compoundLabel']['value']}")
    print(f"  PubChem CID: {compound.get('pubchemCID', {}).get('value', 'N/A')}")

    # Get targets for each compound
    compound_id = compound['compound']['value'].split('/')[-1]
    targets = get_targets(compound_id)
    for target in targets:
        print(f"    Target: {target['proteinLabel']['value']}")
        print(f"      UniProt: {target.get('uniprotID', {}).get('value', 'N/A')}")
```

---

### 4.2 Example 2: PubChem CID → Wikidata → External IDs

```python
def pubchem_to_wikidata(pubchem_cid):
    """Map PubChem CID to Wikidata ID"""
    query = f"""
    SELECT ?compound WHERE {{
      ?compound wdt:P662 "{pubchem_cid}".
    }}
    """

    url = "https://query.wikidata.org/sparql"
    response = requests.get(url, params={'query': query, 'format': 'json'})
    results = response.json()['results']['bindings']

    if results:
        return results[0]['compound']['value'].split('/')[-1]
    return None

def get_all_external_ids(wikidata_id):
    """Get all external database IDs for a compound"""
    query = f"""
    SELECT ?property ?propertyLabel ?value WHERE {{
      wd:{wikidata_id} ?prop ?value.
      ?property wikibase:directClaim ?prop.
      ?property wikibase:propertyType wikibase:ExternalId.
      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
    }}
    """

    url = "https://query.wikidata.org/sparql"
    response = requests.get(url, params={'query': query, 'format': 'json'})
    return response.json()['results']['bindings']

# Usage
pubchem_cid = "969516"  # Curcumin
wikidata_id = pubchem_to_wikidata(pubchem_cid)
print(f"Wikidata ID: {wikidata_id}")

external_ids = get_all_external_ids(wikidata_id)
for ext_id in external_ids:
    print(f"{ext_id['propertyLabel']['value']}: {ext_id['value']['value']}")
```

---

## Part 5: Data Quality Considerations

### 5.1 Completeness

**Wikidata Coverage:**
- ✅ Good: Major medicinal plants, common compounds
- ⚠️ Partial: Protein-compound interactions (not comprehensive)
- ❌ Limited: Obscure traditional medicine compounds, detailed pathway data

**Recommendation:** Use Wikidata for identifier mapping and metadata, NOT as primary source for compound-target interactions.

### 5.2 Validation

Always cross-validate Wikidata claims:
- Check references (`P248` property - "stated in")
- Verify against primary databases (PubChem, UniProt, KEGG)
- Prioritize claims with external database IDs

### 5.3 Update Frequency

- Wikipedia: Updated continuously (real-time)
- Wikidata: Updated continuously (real-time)
- Dumps: Weekly (for bulk processing)

---

## Part 6: Integration Priority

### High Priority Use Cases

1. **Common Name → Scientific Name Mapping**
   - Use Wikipedia API to resolve common names
   - Get Wikidata ID for further queries

2. **External ID Mapping**
   - Use Wikidata to map between PubChem, ChEBI, KEGG, CAS
   - Fill identifier gaps in primary databases

3. **Metadata Enrichment**
   - Get images, descriptions, taxonomy for user-facing features
   - Multi-language support (labels in 300+ languages)

### Medium Priority Use Cases

4. **Compound-Plant Associations**
   - Cross-validate against TCM/Ayurveda databases
   - Discover additional compound-plant links

5. **Protein/Gene Lookup**
   - Map gene symbols to UniProt IDs
   - Get basic protein metadata

### Low Priority Use Cases

6. **Compound-Target Interactions**
   - Coverage is limited compared to specialized databases
   - Better to use ChEMBL, BindingDB, BATMAN-TCM

7. **Pathway Data**
   - Wikidata pathway coverage is minimal
   - Use Reactome, KEGG, WikiPathways instead

---

## Part 7: License & Attribution

- **Wikipedia Content:** CC BY-SA 4.0 (share-alike)
- **Wikidata Data:** CC0 (public domain - no attribution required)
- **Wikipedia API:** Free, no API key required
- **Wikidata Query Service:** Free, no API key required

**Rate Limits:**
- Wikipedia API: ~200 requests/second (respect User-Agent)
- Wikidata Query Service: 60 seconds timeout per query

---

## Part 8: Code Examples

### Python Library: `wikidataintegrator`

```bash
pip install wikidataintegrator
```

```python
from wikidataintegrator import wdi_core

# Get all claims for an entity
item = wdi_core.WDItemEngine(wd_item_id='Q421789')  # Curcumin
print(item.get_wd_json_representation())

# Get specific property
pubchem_id = item.get_property('P662')  # PubChem CID
print(f"PubChem CID: {pubchem_id}")
```

### JavaScript Library: `wikidata-sdk`

```bash
npm install wikidata-sdk
```

```javascript
const wdk = require('wikidata-sdk');

// Search for entities
const url = wdk.searchEntities('curcumin');
console.log(url);

// Get entity data
const entityUrl = wdk.getEntities({
  ids: ['Q421789'],
  languages: ['en'],
  props: ['claims', 'labels']
});
```

---

## Part 9: Recommended Integration Workflow

```
User Input (common name) → Wikipedia API (get Wikidata ID)
                                ↓
                         Wikidata SPARQL
                                ↓
                    ┌───────────┴───────────┐
                    ↓                       ↓
            Get Compounds              Get External IDs
             (P527, P662)              (P352, P665, etc.)
                    ↓                       ↓
            PubChem API              UniProt/KEGG APIs
                    ↓                       ↓
              Full Compound Data      Full Protein/Pathway Data
```

**Key Insight:** Use Wikidata as a "router" to connect to authoritative databases, not as the final data source.

---

## Download

### Wikipedia & Wikidata Integration Sources

| Source | URL | Format | Access | Size |
|--------|-----|--------|--------|------|
| Wikidata dump | https://dumps.wikimedia.org/wikidatawiki/entities/ | JSON Lines (.bz2) | Bulk | ~90 GB |
| Wikipedia API | https://en.wikipedia.org/api/rest_v1/ | REST JSON | Query | Real-time |
| Wikidata SPARQL | https://query.wikidata.org/sparql | SPARQL/JSON | Query | Real-time |
| Wikipedia dump | https://dumps.wikimedia.org/enwiki/ | XML (.bz2) | Bulk | ~20 GB |

### SPARQL Query for Wikipedia-Wikidata Cross-Links

```bash
# Query for plants with Wikipedia articles and compounds
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=
    SELECT ?item ?wikipedia ?label ?compound WHERE {
      ?item wdt:P31 wd:Q756 .
      ?item wdt:P366 wd:Q9690 .
      ?wikipedia schema:about ?item .
      ?wikipedia schema:inLanguage \"en\" .
      OPTIONAL { ?item wdt:P527 ?compound }
      SERVICE wikibase:label { bd:serviceParam wikibase:language \"en\" }
    }
    LIMIT 100" \
  -H "Accept: application/sparql-results+json"
```

---

## Data Format

### Wikipedia REST API Response

**Get page summary:**

```json
{
  "type": "standard",
  "title": "Turmeric",
  "displaytitle": "Turmeric",
  "namespace": { "id": 0, "case": "first-letter", "content": "Main" },
  "wikibase_item": "Q170593",
  "titles": { "canonical": "Turmeric", "normalized": "Turmeric", "display": "Turmeric" },
  "extract": "Turmeric is a flowering plant...",
  "thumbnail": {
    "source": "https://upload.wikimedia.org/wikipedia/commons/...",
    "width": 200,
    "height": 250
  },
  "originalimage": {
    "source": "https://upload.wikimedia.org/wikipedia/commons/...",
    "width": 2448,
    "height": 3264
  },
  "lang": "en",
  "dir": "ltr",
  "revision": "1234567890",
  "lastmodified": "2024-01-15T10:30:45Z"
}
```

### Wikidata & Wikipedia Linked Format

```json
{
  "type": "item",
  "id": "Q170593",
  "labels": { "en": { "value": "turmeric" } },
  "sitelinks": {
    "enwiki": { "site": "enwiki", "title": "Turmeric", "url": "https://en.wikipedia.org/wiki/Turmeric" },
    "dewiki": { "site": "dewiki", "title": "Kurkuma", "url": "https://de.wikipedia.org/wiki/Kurkuma" }
  },
  "claims": {
    "P31": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P31", "datavalue": { "value": { "entity-type": "item", "numeric-id": 756 }, "type": "wikibase-entityid" } } }]
  }
}
```

---

## Schema

### Wikipedia-Wikidata Integration Schema

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| Wikidata Q-ID | String | Item identifier | Q170593 |
| Wikipedia title | String | Article title | Turmeric |
| Language code | String | Wikipedia language | en, de, fr |
| Wikipedia URL | URL | Full article link | https://en.wikipedia.org/wiki/Turmeric |
| Extract | Text | First paragraph summary | "Turmeric is a flowering plant..." |
| Thumbnail | URL | Article image | Wikipedia Commons image URL |
| Last modified | ISO 8601 | Article update time | 2024-01-15T10:30:45Z |

### Cross-Reference Properties

| P-code | Property | Type | Example | Purpose |
|--------|----------|------|---------|---------|
| P31 | Instance of | Item | Q756 (plant) | Type classification |
| P225 | Taxon name | String | Curcuma longa | Scientific name |
| P527 | Has part | Item | Q421789 (curcumin) | Composition |
| P18 | Image | File | Commons file | Visual reference |
| P373 | Commons category | String | Category name | File organization |

---

## Sample Data

### Sample Query: Turmeric with Wikipedia Link

**SPARQL Query:**
```sparql
SELECT ?item ?wikipedia ?label ?extract ?compound WHERE {
  ?item wdt:P31 wd:Q756 .
  ?item rdfs:label "turmeric"@en .
  ?wikipedia schema:about ?item .
  ?wikipedia schema:inLanguage "en" .
  OPTIONAL { ?item wdt:P527 ?compound }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
```

**Sample Results:**
```json
{
  "item": { "value": "http://www.wikidata.org/entity/Q170593" },
  "wikipedia": { "value": "https://en.wikipedia.org/wiki/Turmeric" },
  "label": { "value": "turmeric", "xml:lang": "en" },
  "compound": { "value": "http://www.wikidata.org/entity/Q421789" }
}
```

### Wikipedia REST API Result

**GET /page/summary/Turmeric:**

```json
{
  "title": "Turmeric",
  "wikibase_item": "Q170593",
  "extract": "Turmeric is a flowering plant in the ginger family, Zingiberaceae, the roots of which contain a bright yellow volatile essential oil and a polyphenolic compound called curcumin...",
  "thumbnail": {
    "source": "https://upload.wikimedia.org/wikipedia/commons/b/b6/Curcuma_longa_-_Köhler%E2%80%93s_Medizinal-Pflanzen-150.jpg",
    "width": 200,
    "height": 250
  },
  "lastmodified": "2024-01-15T10:30:45Z"
}
```

### Sample Medicinal Plants with Wikipedia Links

| Q-ID | Plant Name | Wikipedia Title | Language Articles | Main Compounds |
|------|-----------|-----------------|------------------|-----------------|
| Q170593 | Turmeric | Turmeric | 50+ languages | Curcumin |
| Q18208 | Ginseng | Ginseng | 40+ languages | Ginsenosides |
| Q163236 | Echinacea | Echinacea | 30+ languages | Alkylamides |
| Q1455 | Ginkgo biloba | Ginkgo biloba | 45+ languages | Flavonoids |
| Q39647 | Garlic | Garlic | 60+ languages | Allicin |

---

## License

### Wikipedia & Wikidata Licensing

- **Wikidata License:** CC0 1.0 Universal (Public Domain)
- **Wikipedia License:** CC-BY-SA 3.0 / CC-BY-SA 4.0 (depends on language)
- **Requirement for Wikipedia:** Attribution and ShareAlike when redistributing
- **Requirement for Wikidata:** None (but appreciated)
- **Images:** Vary - check each file on Wikimedia Commons

### Citation Format

```
Wikipedia contributors. "[Article Title]." Wikipedia, The Free Encyclopedia.
Wikimedia Foundation, [DATE].
https://en.wikipedia.org/wiki/[Title]

Wikidata contributors. "Wikidata: The Free Knowledge Base That Anyone Can Edit."
Wikimedia Foundation, [DATE].
https://www.wikidata.org/

For images: See Wikimedia Commons file attribution
```

---

## Data Set Size

### Wikipedia-Wikidata Integration Statistics

| Metric | Value | Notes |
|--------|-------|-------|
| Wikidata items with Wikipedia links | ~20 million | All languages combined |
| English Wikipedia articles | ~6.8 million | Main English language edition |
| Items with sitelinks to enwiki | ~6 million | Documented Wikidata items |
| Medicinal plants (with enwiki links) | ~8,000 | Plants with medical use + article |
| Languages with Wikipedia | ~300+ | Multilingual coverage |
| Average article length | ~5-10 KB | Compressed text |

### Storage Requirements

| Component | Size | Format |
|-----------|------|--------|
| Wikidata complete dump | ~90 GB | Compressed JSON Lines |
| Wikipedia dump (enwiki) | ~20 GB | Compressed XML |
| Filtered medicinal plants | ~500 MB | JSON Lines |
| Extracted summaries | ~100 MB | JSON |
| Cross-reference index | ~1 GB | SQLite |

### Update Frequency

- **Wikidata:** Weekly dumps (Monday UTC)
- **Wikipedia:** Monthly dumps (begin of month)
- **Wikipedia REST API:** Real-time (lag <1 second)
- **SPARQL endpoint:** Real-time (lag <1 minute)

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Q-ID | Unique identifier for an item in Wikidata | Q421789 (curcumin) |
| P-ID | Unique identifier for a property in Wikidata | P352 (UniProt ID) |
| SPARQL | Query language for RDF databases | SELECT ?item WHERE {...} |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Wikidata | Free structured knowledge base linked to Wikipedia | Wikimedia Foundation |
| Wikipedia API | RESTful interface to query Wikipedia content | MediaWiki |
| Entity | An item in Wikidata representing a concept, object, or topic | Q-IDs |
| Property | A relationship type connecting entities in Wikidata | P-IDs |
| Claim | A statement asserting a property value for an entity | P661: ChEMBL ID |
| wikidata-sdk | JavaScript library for working with Wikidata | npm package |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST, SPARQL |
| CC0 | Creative Commons Zero (Public Domain) | Wikidata license |
| JSON | JavaScript Object Notation | Data format |
| RDF | Resource Description Framework | Semantic web standard |
| SPARQL | SPARQL Protocol and RDF Query Language | Query language |
| URI | Uniform Resource Identifier | Web address format |
| WD | Wikidata (prefix) | Entity namespace |
| WDT | Wikidata Truthy (prefix) | Property namespace |

---

*Guide compiled January 2026. For related integration guides, see [integration-guide.md](./integration-guide.md)*
