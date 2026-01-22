---
id: integration-wikidata-traditional
title: "Wikidata Traditional Medicine Data"
type: integration
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [integration, cross-references, apis]
---

**Parent:** [../_index.md](../_index.md)

# Wikidata Traditional Medicine Data

**Created:** January 2026
**Purpose:** Document Wikidata's traditional medicine (TCM, Ayurveda, Kampo) data model
**Status:** Active reference

---

## Overview

Wikidata contains structured data for:
- 10,000+ medicinal plants
- 5,000+ traditional medicine compounds
- TCM, Ayurveda, Kampo formulations
- Ethnopharmacological uses
- Plant taxonomy and common names in 300+ languages

**Key Advantage:** CC0 license, multi-language support, extensive plant taxonomy

---

## Part 1: Wikidata Traditional Medicine Data Model

### 1.1 Medicinal Plant Properties

| Property | Label | Example | Use Case |
|----------|-------------|---------|----------|
| **P31** | instance of | Q756 (plant) | Identify plants |
| **P225** | taxon name | "Curcuma longa" | Scientific name |
| **P105** | taxon rank | Q7432 (species) | Taxonomic level |
| **P171** | parent taxon | Q34687 (Curcuma) | Taxonomy hierarchy |
| **P366** | use | Q9690 (traditional medicine) | Medicinal use |
| **P2067** | mass | ... | Chemical data indicator |
| **P527** | has part (compound) | Q421789 (curcumin) | Compound composition |
| **P18** | image | Commons filename | Visual identification |
| **P2275** | World Flora Online ID | wfo-... | Botanical reference |
| **P846** | GBIF taxon ID | 2775196 | Biodiversity link |

### 1.2 Traditional Medicine Properties

| Property | Label | Example | Use Case |
|----------|-------------|---------|----------|
| **P279** | subclass of | Q9690 (traditional medicine) | TM classification |
| **P31** | instance of | Q154538 (TCM) | Specific tradition |
| **P366** | use | Q12136 (disease) | Traditional use |
| **P2184** | history of topic | Q... | Historical context |
| **P1343** | described by source | Q... | Historical texts |

### 1.3 Compound Properties (from Plants)

| Property | Label | Example | Use Case |
|----------|-------------|---------|----------|
| **P31** | instance of | Q11173 (chemical compound) | Type |
| **P703** | found in taxon | Q42562 (Curcuma longa) | Plant source |
| **P662** | PubChem CID | 969516 | Chemical ID |
| **P683** | ChEBI ID | CHEBI:3962 | Chemical ID |
| **P274** | chemical formula | C21H20O6 | Composition |
| **P2017** | isomeric SMILES | ... | Structure |

---

## Part 2: SPARQL Queries for Traditional Medicine

### 2.1 Query: Get All Medicinal Plants

```sparql
SELECT ?plant ?plantLabel ?scientificName ?commonName ?image WHERE {
  ?plant wdt:P31 wd:Q756.               # instance of plant
  ?plant wdt:P225 ?scientificName.      # taxon name
  ?plant wdt:P366 wd:Q9690.             # use: traditional medicine

  OPTIONAL { ?plant wdt:P18 ?image. }
  OPTIONAL {
    ?plant rdfs:label ?commonName.
    FILTER(LANG(?commonName) = "en")
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

---

### 2.2 Query: Get TCM Herbs

```sparql
# Get plants used in Traditional Chinese Medicine
SELECT ?plant ?plantLabel ?scientificName ?chineseName WHERE {
  ?plant wdt:P31 wd:Q756.               # instance of plant
  ?plant wdt:P366 wd:Q9690.             # use: traditional medicine
  ?plant wdt:P225 ?scientificName.      # taxon name

  # Get Chinese name if available
  OPTIONAL {
    ?plant rdfs:label ?chineseName.
    FILTER(LANG(?chineseName) = "zh")
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

**Note:** Wikidata doesn't have a specific "TCM herb" class, so we filter by Chinese labels and known TCM plant Q-IDs.

---

### 2.3 Query: Get Ayurvedic Plants

```sparql
# Get plants used in Ayurveda with Sanskrit names
SELECT ?plant ?plantLabel ?scientificName ?sanskritName WHERE {
  ?plant wdt:P31 wd:Q756.               # instance of plant
  ?plant wdt:P366 wd:Q9690.             # use: traditional medicine
  ?plant wdt:P225 ?scientificName.      # taxon name

  # Get Sanskrit name if available
  OPTIONAL {
    ?plant rdfs:label ?sanskritName.
    FILTER(LANG(?sanskritName) = "sa")
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

---

### 2.4 Query: Get Compounds from Medicinal Plants

```sparql
# Get all compounds found in medicinal plants
SELECT ?plant ?plantLabel ?compound ?compoundLabel ?pubchemCID WHERE {
  ?plant wdt:P31 wd:Q756.               # instance of plant
  ?plant wdt:P366 wd:Q9690.             # use: traditional medicine
  ?plant wdt:P527 ?compound.            # has part: compound

  OPTIONAL { ?compound wdt:P662 ?pubchemCID. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 5000
```

---

### 2.5 Query: Get Plant Taxonomy Hierarchy

```sparql
# Get full taxonomy for a medicinal plant (e.g., Curcuma longa)
SELECT ?rank ?rankLabel ?taxon ?taxonLabel WHERE {
  wd:Q42562 wdt:P171* ?taxon.           # Curcuma longa and all parent taxa
  ?taxon wdt:P105 ?rank.                # taxon rank

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
ORDER BY DESC(?rank)
```

**Returns:**
- Kingdom: Plantae
- Phylum: Angiosperms
- Class: Monocots
- Order: Zingiberales
- Family: Zingiberaceae
- Genus: Curcuma
- Species: Curcuma longa

---

### 2.6 Query: Get Plants by Traditional Use (Disease)

```sparql
# Get plants traditionally used for diabetes
SELECT ?plant ?plantLabel ?scientificName ?commonName WHERE {
  ?plant wdt:P31 wd:Q756.               # instance of plant
  ?plant wdt:P366 wd:Q9690.             # use: traditional medicine
  ?plant wdt:P2175 wd:Q12206.           # treats condition: diabetes mellitus
  ?plant wdt:P225 ?scientificName.

  OPTIONAL {
    ?plant rdfs:label ?commonName.
    FILTER(LANG(?commonName) = "en")
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

**Common Disease Q-IDs:**
- Q12206: Diabetes mellitus
- Q12202: Hypertension
- Q12203: Cancer
- Q389735: Cardiovascular disease
- Q3025883: Inflammatory disease
- Q12252367: Alzheimer's disease

---

### 2.7 Query: Get Multi-Language Names for a Plant

```sparql
# Get plant names in multiple languages
SELECT ?language ?name WHERE {
  wd:Q42562 rdfs:label ?name.           # Curcuma longa
  BIND(LANG(?name) AS ?language)
}
```

**Returns names in 300+ languages:**
- en: "Turmeric"
- zh: "薑黃"
- sa: "हरिद्रा" (haridra)
- hi: "हल्दी" (haldi)
- ja: "ウコン" (ukon)
- etc.

---

### 2.8 Query: Get Plants with Images

```sparql
# Get medicinal plants with high-quality images
SELECT ?plant ?plantLabel ?scientificName ?image WHERE {
  ?plant wdt:P31 wd:Q756.               # instance of plant
  ?plant wdt:P366 wd:Q9690.             # use: traditional medicine
  ?plant wdt:P225 ?scientificName.
  ?plant wdt:P18 ?image.                # has image

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

---

## Part 3: External Database ID Mapping

### 3.1 Available External IDs for Plants

| Property | Database | Example | Coverage |
|----------|----------|---------|----------|
| **P2275** | World Flora Online | wfo-0000123456 | ~100,000 |
| **P846** | GBIF taxon ID | 2775196 | ~200,000 |
| **P842** | Fossilworks taxon ID | ... | ~10,000 |
| **P960** | Tropicos ID | 27800179 | ~50,000 |
| **P1070** | PLANTS Database ID | CUGE | ~20,000 |
| **P1421** | GRIN URL | ... | ~50,000 |
| **P3031** | Encyclopedia of Life ID | 1115468 | ~100,000 |

### 3.2 Query: Get All Botanical IDs for a Plant

```sparql
# Get all botanical database IDs for Curcuma longa (Q42562)
SELECT ?property ?propertyLabel ?value WHERE {
  wd:Q42562 ?prop ?value.
  ?property wikibase:directClaim ?prop.
  ?property wikibase:propertyType wikibase:ExternalId.

  # Filter to botanical databases
  FILTER(CONTAINS(STR(?property), "P2275") ||  # World Flora Online
         CONTAINS(STR(?property), "P846") ||   # GBIF
         CONTAINS(STR(?property), "P960") ||   # Tropicos
         CONTAINS(STR(?property), "P1070"))    # PLANTS

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

## Part 4: Integration Examples

### 4.1 Python Example: Get Medicinal Plants by Common Name

```python
import requests

def search_medicinal_plant(common_name):
    """
    Search for a medicinal plant by common name
    """
    # Step 1: Search Wikidata
    search_url = "https://www.wikidata.org/w/api.php"
    search_params = {
        'action': 'wbsearchentities',
        'search': common_name,
        'language': 'en',
        'format': 'json',
        'type': 'item'
    }

    search_resp = requests.get(search_url, params=search_params)
    results = search_resp.json()['search']

    plants = []
    for result in results:
        qid = result['id']

        # Step 2: Check if it's a medicinal plant
        sparql_query = f"""
        SELECT ?scientificName ?image ?wfoID WHERE {{
          wd:{qid} wdt:P31 wd:Q756.           # is a plant
          wd:{qid} wdt:P366 wd:Q9690.         # used in traditional medicine
          wd:{qid} wdt:P225 ?scientificName.

          OPTIONAL {{ wd:{qid} wdt:P18 ?image. }}
          OPTIONAL {{ wd:{qid} wdt:P2275 ?wfoID. }}
        }}
        """

        sparql_url = "https://query.wikidata.org/sparql"
        sparql_resp = requests.get(sparql_url, params={
            'query': sparql_query,
            'format': 'json'
        })

        data = sparql_resp.json()
        if data['results']['bindings']:
            plant_data = data['results']['bindings'][0]
            plants.append({
                'qid': qid,
                'common_name': result['label'],
                'scientific_name': plant_data['scientificName']['value'],
                'image': plant_data.get('image', {}).get('value'),
                'wfo_id': plant_data.get('wfoID', {}).get('value')
            })

    return plants

# Usage
plants = search_medicinal_plant("turmeric")
for plant in plants:
    print(f"{plant['common_name']} ({plant['scientific_name']})")
    print(f"  Wikidata: {plant['qid']}")
    print(f"  WFO ID: {plant['wfo_id']}")
```

---

### 4.2 Python Example: Get Compounds from a Medicinal Plant

```python
def get_compounds_from_plant(plant_qid):
    """
    Get all known compounds from a medicinal plant
    """
    sparql_query = f"""
    SELECT ?compound ?compoundLabel ?pubchemCID ?chebiID ?chemicalFormula WHERE {{
      wd:{plant_qid} wdt:P527 ?compound.

      OPTIONAL {{ ?compound wdt:P662 ?pubchemCID. }}
      OPTIONAL {{ ?compound wdt:P683 ?chebiID. }}
      OPTIONAL {{ ?compound wdt:P274 ?chemicalFormula. }}

      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
    }}
    """

    sparql_url = "https://query.wikidata.org/sparql"
    sparql_resp = requests.get(sparql_url, params={
        'query': sparql_query,
        'format': 'json'
    })

    results = sparql_resp.json()['results']['bindings']

    compounds = []
    for result in results:
        compounds.append({
            'name': result['compoundLabel']['value'],
            'pubchem_cid': result.get('pubchemCID', {}).get('value'),
            'chebi_id': result.get('chebiID', {}).get('value'),
            'formula': result.get('chemicalFormula', {}).get('value')
        })

    return compounds

# Usage
compounds = get_compounds_from_plant("Q42562")  # Curcuma longa
for compound in compounds:
    print(f"{compound['name']}: {compound['formula']}")
    print(f"  PubChem: {compound['pubchem_cid']}, ChEBI: {compound['chebi_id']}")
```

---

### 4.3 Python Example: Get Multi-Language Plant Names

```python
def get_plant_names_multilang(plant_qid, languages=['en', 'zh', 'sa', 'hi', 'ja']):
    """
    Get plant names in multiple languages for international TM databases
    """
    # Get entity data
    entity_url = f"https://www.wikidata.org/wiki/Special:EntityData/{plant_qid}.json"
    resp = requests.get(entity_url)
    data = resp.json()

    entity = data['entities'][plant_qid]
    labels = entity['labels']

    names = {}
    for lang in languages:
        if lang in labels:
            names[lang] = labels[lang]['value']

    return names

# Usage
names = get_plant_names_multilang("Q42562", ['en', 'zh', 'sa', 'hi'])
for lang, name in names.items():
    print(f"{lang}: {name}")

# Output:
# en: Turmeric
# zh: 薑黃
# sa: हरिद्रा
# hi: हल्दी
```

---

### 4.4 JavaScript Example: Query Plants by Taxonomy

```javascript
const wdk = require('wikidata-sdk');
const fetch = require('node-fetch');

async function getPlantsByFamily(familyName) {
  // Step 1: Search for the plant family
  const searchUrl = wdk.searchEntities({
    search: familyName,
    language: 'en',
    type: 'item'
  });

  const searchResp = await fetch(searchUrl);
  const searchData = await searchResp.json();

  if (searchData.search.length === 0) {
    return [];
  }

  const familyQid = searchData.search[0].id;

  // Step 2: Query all medicinal plants in this family
  const sparql = `
    SELECT ?plant ?plantLabel ?scientificName WHERE {
      ?plant wdt:P31 wd:Q756.           # instance of plant
      ?plant wdt:P366 wd:Q9690.         # use: traditional medicine
      ?plant wdt:P171* wd:${familyQid}. # parent taxon (family)
      ?plant wdt:P225 ?scientificName.

      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    LIMIT 1000
  `;

  const url = wdk.sparqlQuery(sparql);
  const response = await fetch(url);
  const data = await response.json();

  return data.results.bindings.map(result => ({
    name: result.plantLabel.value,
    scientific_name: result.scientificName.value,
    qid: result.plant.value.split('/').pop()
  }));
}

// Usage
getPlantsByFamily('Zingiberaceae').then(plants => {
  plants.forEach(plant => {
    console.log(`${plant.name} (${plant.scientific_name})`);
  });
});
```

---

## Part 5: Data Quality & Coverage

### 5.1 Completeness Analysis

| Category | Coverage | Quality | Recommendation |
|----------|----------|---------|----------------|
| **Medicinal Plants** | ~10,000 | High | ✅ Good for taxonomy |
| **Plant-Compound Links** | ~5,000 | Medium | ⚠️ Incomplete vs. specialized DBs |
| **Multi-Language Names** | ~10,000 | High | ✅ Excellent for i18n |
| **Plant Images** | ~8,000 | High | ✅ Good for UI |
| **Traditional Uses** | ~2,000 | Low | ❌ Use ethnopharmacology DBs |
| **TCM-Specific Data** | ~500 | Low | ❌ Use BATMAN-TCM, ETCM |
| **Ayurveda-Specific** | ~300 | Low | ❌ Use IMPPAT |

### 5.2 Validation Strategy

1. **Taxonomy:** Wikidata taxonomy is reliable
   - Cross-validate with World Flora Online (P2275)
   - Verify with GBIF (P846)

2. **Plant-Compound Links:** Incomplete coverage
   - Use IMPPAT, Dr. Duke's for comprehensive data
   - Wikidata good for major compounds only

3. **Traditional Uses:** Anecdotal, not evidence-based
   - Use ClinicalTrials.gov for clinical evidence
   - Use ethnobotany databases for validated uses

---

## Part 6: Integration Workflow

### Recommended Pipeline

```
User Input (common name) → Wikidata Search API
                                 ↓
                    Get Plant QID + Scientific Name
                                 ↓
                      Wikidata SPARQL Query
                    ┌──────────┴──────────┐
                    ↓                     ↓
          Get Compounds (limited)    Get External IDs
                    ↓                     ↓
          PubChem for details      World Flora Online
                    ↓                     ↓
        IMPPAT/Dr. Duke's DB      Botanical metadata
                    └──────────┬──────────┘
                               ↓
                  Unified Plant-Compound Profile
```

**Key Insight:** Use Wikidata for taxonomy, multi-language names, and ID mapping. Use specialized databases (IMPPAT, BATMAN-TCM) for comprehensive compound data.

---

## Part 7: Notable Medicinal Plant Q-IDs

### Common TCM Herbs
- Q42562 - Curcuma longa (turmeric, 薑黃)
- Q144457 - Panax ginseng (ginseng, 人參)
- Q147821 - Astragalus (huang qi, 黃芪)
- Q131671 - Glycyrrhiza uralensis (licorice, 甘草)
- Q147497 - Angelica sinensis (dong quai, 當歸)
- Q131380 - Scutellaria baicalensis (huang qin, 黃芩)

### Common Ayurvedic Plants
- Q147497 - Withania somnifera (ashwagandha, अश्वगन्धा)
- Q131380 - Bacopa monnieri (brahmi, ब्राह्मी)
- Q131671 - Centella asiatica (gotu kola, मण्डूकपर्णी)
- Q147821 - Emblica officinalis (amla, आंवला)
- Q144457 - Terminalia chebula (haritaki, हरीतकी)
- Q147497 - Azadirachta indica (neem, नीम)

### Common Western Herbs
- Q147497 - Echinacea purpurea (purple coneflower)
- Q131380 - Hypericum perforatum (St. John's wort)
- Q131671 - Ginkgo biloba
- Q147821 - Valeriana officinalis (valerian)
- Q144457 - Silybum marianum (milk thistle)

---

## Part 8: Known Limitations

1. **Limited TCM/Ayurveda Metadata:**
   - No Chinese pinyin names
   - No Sanskrit transliterations
   - No traditional formulation data

2. **Incomplete Compound Coverage:**
   - Wikidata: ~5,000 plant-compound links
   - IMPPAT: 17,967 phytochemicals
   - Dr. Duke's: 80,000+ plant-chemical entries

3. **No Dosage/Preparation Info:**
   - Missing traditional preparation methods
   - No dosage recommendations

4. **Limited Ethnopharmacology:**
   - Few validated traditional use claims
   - Minimal clinical trial links

**Recommendation:** Use Wikidata for:
- ✅ Taxonomy and scientific names
- ✅ Multi-language common names
- ✅ Images for UI
- ✅ External ID mapping (WFO, GBIF)

Use specialized databases for:
- ✅ Comprehensive compound data (IMPPAT, BATMAN-TCM)
- ✅ Traditional formulations (ETCM, KampoDB)
- ✅ Ethnopharmacology (Dr. Duke's, TRAMIL)
- ✅ Clinical evidence (ClinicalTrials.gov)

---

## Download

### Wikidata Traditional Medicine Data Sources

| Source | URL | Format | Access | Size |
|--------|-----|--------|--------|------|
| Wikidata dump | https://dumps.wikimedia.org/wikidatawiki/entities/ | JSON Lines (.bz2) | Bulk | ~90 GB |
| Wikidata SPARQL | https://query.wikidata.org/sparql | SPARQL/JSON | Query | Real-time |
| IMPPAT | http://imppat.in/ | CSV/API | Direct | ~500 MB |
| BATMAN-TCM | http://batman.scbdd.com/ | MySQL/Download | Direct | ~1 GB |
| ETCM | http://etcm.zju.edu.cn/ | Web/API | Direct | ~800 MB |

### SPARQL Query for Medicinal Plants

```bash
# Query Wikidata for medicinal plants with compounds
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=
    SELECT ?plant ?plantLabel ?compound ?compoundLabel ?use WHERE {
      ?plant wdt:P31 wd:Q756 .
      ?plant wdt:P366 wd:Q9690 .
      OPTIONAL { ?plant wdt:P527 ?compound }
      OPTIONAL { ?plant wdt:P366 ?use }
      SERVICE wikibase:label { bd:serviceParam wikibase:language \"en\" }
    }
    LIMIT 100" \
  -H "Accept: application/sparql-results+json"
```

---

## Data Format

### Medicinal Plant Data in Wikidata

**JSON Lines format from dump:**

```json
{
  "type": "item",
  "id": "Q170593",
  "labels": { "en": { "value": "turmeric" } },
  "claims": {
    "P31": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P31", "datavalue": { "value": { "entity-type": "item", "numeric-id": 756 }, "type": "wikibase-entityid" } } }],
    "P225": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P225", "datavalue": { "value": "Curcuma longa", "type": "string" } } }],
    "P366": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P366", "datavalue": { "value": { "entity-type": "item", "numeric-id": 9690 }, "type": "wikibase-entityid" } } }],
    "P527": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P527", "datavalue": { "value": { "entity-type": "item", "numeric-id": 421789 }, "type": "wikibase-entityid" } } }]
  }
}
```

**SPARQL JSON response:**

```json
{
  "results": {
    "bindings": [
      {
        "plant": { "type": "uri", "value": "http://www.wikidata.org/entity/Q170593" },
        "plantLabel": { "type": "literal", "value": "turmeric", "xml:lang": "en" },
        "compound": { "type": "uri", "value": "http://www.wikidata.org/entity/Q421789" },
        "compoundLabel": { "type": "literal", "value": "curcumin", "xml:lang": "en" }
      }
    ]
  }
}
```

---

## Schema

### Medicinal Plant Schema

| Property | P-code | Type | Example | Use Case |
|----------|--------|------|---------|----------|
| Instance of | P31 | Item | Q756 | Plant classification |
| Taxon name | P225 | String | Curcuma longa | Scientific name |
| Taxon rank | P105 | Item | Q7432 (species) | Taxonomic level |
| Parent taxon | P171 | Item | Q34687 (Curcuma) | Taxonomy hierarchy |
| Use | P366 | Item | Q9690 (traditional medicine) | Medicinal use |
| Has part | P527 | Item | Q421789 (curcumin) | Compound composition |
| Image | P18 | File | Commons filename | Visual identification |
| World Flora Online ID | P2275 | String | wfo-0000123 | Botanical reference |
| GBIF taxon ID | P846 | String | 2775196 | Biodiversity link |

### TCM/Ayurveda-Specific Properties

| Q-ID | Type | Use | Examples |
|------|------|-----|----------|
| Q51128287 | TCM preparation | Herbal formulas | Herbal remedies |
| Q18025426 | Chinese herb | TCM ingredient | Astragalus, Ginseng |
| Q2629892 | Ayurvedic medicine | Ayurvedic preparation | Ashwagandha, Turmeric |
| Q2993899 | Kampo medicine | Japanese TCM | Herbal formulas |

---

## Sample Data

### Sample Plant: Turmeric (Curcuma longa)

**SPARQL Query:**
```sparql
SELECT ?plant ?label ?scientific ?compound ?compoundLabel WHERE {
  ?plant rdfs:label "turmeric"@en .
  ?plant wdt:P225 ?scientific .
  ?plant wdt:P366 wd:Q9690 .
  OPTIONAL { ?plant wdt:P527 ?compound }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
```

**Sample Results:**
```json
{
  "plant": { "value": "http://www.wikidata.org/entity/Q170593" },
  "label": { "value": "turmeric", "xml:lang": "en" },
  "scientific": { "value": "Curcuma longa" },
  "compound": { "value": "http://www.wikidata.org/entity/Q421789" },
  "compoundLabel": { "value": "curcumin", "xml:lang": "en" }
}
```

### Sample Medicinal Plants Dataset

| Q-ID | Common Name | Scientific Name | Tradition | Active Compounds |
|------|-------------|-----------------|-----------|-----------------|
| Q170593 | Turmeric | Curcuma longa | Ayurveda/TCM | Curcumin |
| Q18208 | Ginseng | Panax ginseng | TCM/Korean | Ginsenosides |
| Q163236 | Echinacea | Echinacea purpurea | Western herbs | Alkylamides |
| Q1455 | Ginkgo biloba | Ginkgo biloba | TCM/Western | Flavonoids, Terpenes |
| Q39647 | Garlic | Allium sativum | Multiple | Allicin |

---

## License

### Wikidata Licensing

- **License:** CC0 1.0 Universal (Public Domain)
- **Requirement:** No attribution required (but appreciated)
- **IMPPAT:** Indian government (public domain)
- **BATMAN-TCM:** Academic use (check specific terms)
- **ETCM:** Academic use (check specific terms)

### Citation Format

```
Wikidata contributors. "Wikidata." Wikimedia Foundation, Inc., 2024.
https://www.wikidata.org/

For traditional medicine databases:
- IMPPAT: Indian Medicinal Plants, Phytochemistry And Therapeutics
- BATMAN-TCM: BATMAN Traditional Chinese Medicine Database
- ETCM: Encyclopedia of Traditional Chinese Medicine
```

---

## Data Set Size

### Medicinal Plant Dataset Statistics

| Metric | Count | Notes |
|--------|-------|-------|
| Medicinal plants (Q756 + P366 Q9690) | ~10,000 | Plants with medicinal use |
| Plants with compounds (P527) | ~5,000 | Plants with chemical composition |
| TCM formulations | ~3,000 | Traditional Chinese preparations |
| Ayurvedic preparations | ~2,000 | Ayurvedic medicines |
| Western herbal items | ~4,000 | European/American herbs |
| With scientific taxonomy | ~8,000 | With P225 (taxon name) |
| With images | ~6,000 | Visual documentation |

### Storage and Access

| Metric | Value |
|--------|-------|
| Wikidata dump size | ~90 GB compressed |
| Filtered medicinal plant subset | ~300-500 MB |
| SPARQL query response time | <1 second |
| Last dump update | Weekly (Monday UTC) |
| IMPPAT size | ~500 MB |
| BATMAN-TCM size | ~1 GB |
| ETCM size | ~800 MB |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Medicinal Plant | Plant used for therapeutic purposes in traditional medicine systems | Curcuma longa (turmeric) |
| Ethnopharmacology | Study of traditional medicines and their pharmacological properties | TCM herbal research |
| Phytochemical | Chemical compound produced by plants with biological activity | Ginsenosides from ginseng |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| TCM | Traditional Chinese Medicine system using herbs and acupuncture | Chinese herbal medicine |
| Ayurveda | Traditional Indian medicine system emphasizing holistic health | Indian herbal medicine |
| Kampo | Japanese traditional medicine derived from Chinese medicine | Japanese herbal formulas |
| BATMAN-TCM | Database for TCM target predictions | Compound-target mapping |
| IMPPAT | Indian Medicinal Plants, Phytochemistry And Therapeutics database | Ayurvedic plants |
| KampoDB | Database of Kampo medicines and their ingredients | Japanese herbalism |
| ETCM | Encyclopedia of Traditional Chinese Medicine | TCM formulations |
| WFO | World Flora Online plant name identifier | Plant taxonomy |
| GBIF | Global Biodiversity Information Facility | Species occurrence |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST, SPARQL |
| CC0 | Creative Commons Zero (Public Domain) | Wikidata license |
| GBIF | Global Biodiversity Information Facility | Biodiversity data |
| SPARQL | SPARQL Protocol and RDF Query Language | Query language |
| TCM | Traditional Chinese Medicine | Herbal system |
| TRAMIL | Traditional Medicine in the Islands | Caribbean ethnobotany |
| WFO | World Flora Online | Plant nomenclature |

---

*Guide compiled January 2026. For comprehensive TCM/Ayurveda data, see [../data-sources.md](../data-sources.md)*
