---
id: integration-wikidata-supplements
title: "Wikidata Dietary Supplements & Nutraceuticals"
type: integration
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [integration, cross-references, apis]
---

**Parent:** [../_index.md](../_index.md)

# Wikidata Dietary Supplements & Nutraceuticals

**Created:** January 2026
**Purpose:** Document Wikidata's dietary supplement and nutraceutical data model
**Status:** Active reference

---

## Overview

Wikidata contains structured data for:
- 5,000+ dietary supplements
- Vitamins, minerals, amino acids
- Herbal supplements
- Probiotic strains
- Omega-3 fatty acids
- Protein powders
- Sports nutrition

**Key Advantage:** CC0 license (public domain), multi-language support, linked to regulatory databases

---

## Part 1: Wikidata Supplement Data Model

### 1.1 Core Supplement Properties

| Property | Label | Example | Use Case |
|----------|-------|---------|----------|
| **P31** | instance of | Q28885102 (dietary supplement) | Identify supplements |
| **P3781** | has active ingredient | Q vitamin | Ingredient mapping |
| **P2175** | medical condition treated | Q disease | Health claims |
| **P1051** | NutritionFacts.org ID | ... | Nutrition data |
| **P715** | DrugBank ID | DB00116 (vitamin) | Cross-reference |
| **P3345** | Health Canada LNHPD | 12345 | Canadian approval |
| **P486** | MeSH ID | D014810 (vitamins) | Medical terminology |
| **P2892** | UMLS CUI | C0042890 (vitamin) | Medical concept |

### 1.2 Ingredient Properties

| Property | Label | Example | Use Case |
|----------|-------|---------|----------|
| **P31** | instance of | Q34956 (vitamin) | Type classification |
| **P527** | has part | Q compound | Multi-ingredient products |
| **P662** | PubChem CID | 5280795 (vitamin D3) | Chemical ID |
| **P231** | CAS Registry Number | 67-97-0 | Chemical ID |
| **P274** | chemical formula | C27H44O | Chemical composition |
| **P2114** | recommended dietary allowance | Q... | Dosage info |

### 1.3 Taxonomic Properties (Herbal Supplements)

| Property | Label | Example | Use Case |
|----------|-------|---------|----------|
| **P703** | found in taxon | Q42562 (Curcuma longa) | Plant source |
| **P225** | taxon name | "Curcuma longa" | Scientific name |
| **P18** | image | Commons filename | Visual ID |
| **P2275** | World Flora Online ID | wfo-... | Botanical reference |

---

## Part 2: SPARQL Queries for Supplements

### 2.1 Query: Get All Dietary Supplements

```sparql
SELECT ?supplement ?supplementLabel ?activeIngredient ?activeIngredientLabel WHERE {
  ?supplement wdt:P31 wd:Q28885102.     # instance of dietary supplement

  OPTIONAL { ?supplement wdt:P3781 ?activeIngredient. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

---

### 2.2 Query: Get Vitamin Supplements

```sparql
# Get all vitamin supplements with chemical IDs
SELECT ?vitamin ?vitaminLabel ?pubchemCID ?casNumber ?chemicalFormula WHERE {
  ?vitamin wdt:P31 wd:Q34956.           # instance of vitamin

  OPTIONAL { ?vitamin wdt:P662 ?pubchemCID. }
  OPTIONAL { ?vitamin wdt:P231 ?casNumber. }
  OPTIONAL { ?vitamin wdt:P274 ?chemicalFormula. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

**Common Vitamin Q-IDs:**
- Q18225 - Vitamin A
- Q192111 - Thiamine (B1)
- Q192423 - Riboflavin (B2)
- Q192423 - Niacin (B3)
- Q164787 - Pantothenic acid (B5)
- Q192760 - Pyridoxine (B6)
- Q18228 - Biotin (B7)
- Q165606 - Folate (B9)
- Q18225 - Cobalamin (B12)
- Q18225 - Ascorbic acid (C)
- Q18225 - Cholecalciferol (D3)
- Q18225 - Tocopherol (E)
- Q18225 - Phylloquinone (K1)

---

### 2.3 Query: Get Herbal Supplements by Plant Source

```sparql
# Get supplements derived from Turmeric (Curcuma longa)
SELECT ?supplement ?supplementLabel ?activeIngredient ?activeIngredientLabel ?pubchemCID WHERE {
  ?supplement wdt:P31 wd:Q28885102.     # instance of dietary supplement
  ?supplement wdt:P3781 ?activeIngredient.
  ?activeIngredient wdt:P703 wd:Q42562. # found in taxon: Curcuma longa

  OPTIONAL { ?activeIngredient wdt:P662 ?pubchemCID. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

### 2.4 Query: Get Supplements with Health Claims

```sparql
# Get supplements claimed to treat cardiovascular disease
SELECT ?supplement ?supplementLabel ?condition ?conditionLabel ?activeIngredient WHERE {
  ?supplement wdt:P31 wd:Q28885102.     # instance of dietary supplement
  ?supplement wdt:P2175 ?condition.     # medical condition treated
  ?condition wdt:P279* wd:Q389735.      # subclass of cardiovascular disease

  OPTIONAL { ?supplement wdt:P3781 ?activeIngredient. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

**Common Health Condition Q-IDs:**
- Q389735 - Cardiovascular disease
- Q12206 - Diabetes mellitus
- Q3025883 - Inflammatory disease
- Q2648537 - Osteoporosis
- Q12252367 - Alzheimer's disease
- Q18123741 - Metabolic syndrome

---

### 2.5 Query: Get Probiotic Strains

```sparql
# Get probiotic bacterial strains used in supplements
SELECT ?strain ?strainLabel ?species ?speciesLabel WHERE {
  ?strain wdt:P31 wd:Q855769.           # instance of bacterial strain
  ?strain wdt:P703 ?species.            # found in taxon (species)
  ?strain wdt:P366 wd:Q28885102.        # use: dietary supplement

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

**Common Probiotic Species Q-IDs:**
- Q131140 - Lactobacillus
- Q131140 - Bifidobacterium
- Q131140 - Streptococcus thermophilus

---

### 2.6 Query: Get Omega-3 Supplements

```sparql
# Get omega-3 fatty acid supplements
SELECT ?supplement ?supplementLabel ?omega3 ?omega3Label ?pubchemCID WHERE {
  ?supplement wdt:P31 wd:Q28885102.     # instance of dietary supplement
  ?supplement wdt:P3781 ?omega3.        # has active ingredient
  ?omega3 wdt:P31 wd:Q425213.           # instance of omega-3 fatty acid

  OPTIONAL { ?omega3 wdt:P662 ?pubchemCID. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

**Common Omega-3 Q-IDs:**
- Q425213 - Omega-3 fatty acid (class)
- Q425213 - EPA (eicosapentaenoic acid)
- Q425213 - DHA (docosahexaenoic acid)
- Q425213 - ALA (alpha-linolenic acid)

---

### 2.7 Query: Get Amino Acid Supplements

```sparql
# Get amino acid supplements
SELECT ?aminoAcid ?aminoAcidLabel ?pubchemCID ?chemicalFormula WHERE {
  ?aminoAcid wdt:P31 wd:Q8054.          # instance of amino acid
  ?aminoAcid wdt:P366 wd:Q28885102.     # use: dietary supplement

  OPTIONAL { ?aminoAcid wdt:P662 ?pubchemCID. }
  OPTIONAL { ?aminoAcid wdt:P274 ?chemicalFormula. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

### 2.8 Query: Get Supplements by Regulatory Status

```sparql
# Get supplements approved by Health Canada (LNHPD)
SELECT ?supplement ?supplementLabel ?lnhpdID ?ingredient WHERE {
  ?supplement wdt:P31 wd:Q28885102.     # instance of dietary supplement
  ?supplement wdt:P3345 ?lnhpdID.       # Health Canada LNHPD ID

  OPTIONAL { ?supplement wdt:P3781 ?ingredient. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

## Part 3: External Database ID Mapping

### 3.1 Available External IDs for Supplements

| Property | Database | Example | Coverage |
|----------|----------|---------|----------|
| **P3345** | Health Canada LNHPD | 12345 | ~2,000 |
| **P1051** | NutritionFacts.org | ... | ~500 |
| **P715** | DrugBank | DB00116 (vitamins) | ~200 |
| **P662** | PubChem CID | 5280795 | ~5,000 |
| **P231** | CAS Registry Number | 67-97-0 | ~3,000 |
| **P486** | MeSH ID | D014810 | ~1,000 |
| **P2892** | UMLS CUI | C0042890 | ~2,000 |

### 3.2 Query: Get All External IDs for a Supplement

```sparql
# Get all external IDs for Vitamin D (Q18225)
SELECT ?property ?propertyLabel ?value WHERE {
  wd:Q18225 ?prop ?value.
  ?property wikibase:directClaim ?prop.
  ?property wikibase:propertyType wikibase:ExternalId.

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

## Part 4: Integration Examples

### 4.1 Python Example: Get Herbal Supplements for a Plant

```python
import requests

def get_supplements_from_plant(plant_name):
    """
    Get dietary supplements derived from a specific plant
    """
    # Step 1: Search Wikidata for the plant
    search_url = "https://www.wikidata.org/w/api.php"
    search_params = {
        'action': 'wbsearchentities',
        'search': plant_name,
        'language': 'en',
        'format': 'json',
        'type': 'item'
    }

    search_resp = requests.get(search_url, params=search_params)
    plant_results = search_resp.json()['search']

    if not plant_results:
        return None

    plant_qid = plant_results[0]['id']
    print(f"Plant QID: {plant_qid}")

    # Step 2: Query supplements from this plant
    sparql_query = f"""
    SELECT ?supplement ?supplementLabel ?activeIngredient ?activeIngredientLabel ?pubchemCID WHERE {{
      ?supplement wdt:P31 wd:Q28885102.
      ?supplement wdt:P3781 ?activeIngredient.
      ?activeIngredient wdt:P703 wd:{plant_qid}.
      OPTIONAL {{ ?activeIngredient wdt:P662 ?pubchemCID. }}
      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
    }}
    """

    sparql_url = "https://query.wikidata.org/sparql"
    sparql_resp = requests.get(sparql_url, params={
        'query': sparql_query,
        'format': 'json'
    })

    results = sparql_resp.json()['results']['bindings']

    supplements = []
    for result in results:
        supplements.append({
            'supplement_name': result['supplementLabel']['value'],
            'active_ingredient': result['activeIngredientLabel']['value'],
            'pubchem_cid': result.get('pubchemCID', {}).get('value')
        })

    return supplements

# Usage
supplements = get_supplements_from_plant("Curcuma longa")
for supp in supplements:
    print(f"{supp['supplement_name']}: {supp['active_ingredient']} (PubChem: {supp['pubchem_cid']})")
```

---

### 4.2 Python Example: Get Vitamins with Dosage Info

```python
def get_vitamin_info(vitamin_name):
    """
    Get vitamin information including RDA (recommended dietary allowance)
    """
    # Step 1: Search for the vitamin
    search_url = "https://www.wikidata.org/w/api.php"
    search_params = {
        'action': 'wbsearchentities',
        'search': vitamin_name,
        'language': 'en',
        'format': 'json',
        'type': 'item'
    }

    search_resp = requests.get(search_url, params=search_params)
    vitamin_results = search_resp.json()['search']

    if not vitamin_results:
        return None

    vitamin_qid = vitamin_results[0]['id']

    # Step 2: Query vitamin data
    sparql_query = f"""
    SELECT ?vitaminLabel ?pubchemCID ?casNumber ?chemicalFormula WHERE {{
      wd:{vitamin_qid} rdfs:label ?vitaminLabel.
      OPTIONAL {{ wd:{vitamin_qid} wdt:P662 ?pubchemCID. }}
      OPTIONAL {{ wd:{vitamin_qid} wdt:P231 ?casNumber. }}
      OPTIONAL {{ wd:{vitamin_qid} wdt:P274 ?chemicalFormula. }}
      FILTER(LANG(?vitaminLabel) = "en")
    }}
    """

    sparql_url = "https://query.wikidata.org/sparql"
    sparql_resp = requests.get(sparql_url, params={
        'query': sparql_query,
        'format': 'json'
    })

    results = sparql_resp.json()['results']['bindings']

    if results:
        result = results[0]
        return {
            'name': result['vitaminLabel']['value'],
            'pubchem_cid': result.get('pubchemCID', {}).get('value'),
            'cas_number': result.get('casNumber', {}).get('value'),
            'formula': result.get('chemicalFormula', {}).get('value')
        }

    return None

# Usage
vitamin_d = get_vitamin_info("vitamin D3")
print(vitamin_d)
```

---

### 4.3 JavaScript Example: Query Supplements by Health Claim

```javascript
const wdk = require('wikidata-sdk');
const fetch = require('node-fetch');

async function getSupplementsByHealthClaim(conditionName) {
  // Step 1: Search for the health condition
  const searchUrl = wdk.searchEntities({
    search: conditionName,
    language: 'en',
    type: 'item'
  });

  const searchResp = await fetch(searchUrl);
  const searchData = await searchResp.json();

  if (searchData.search.length === 0) {
    return [];
  }

  const conditionQid = searchData.search[0].id;

  // Step 2: Query supplements for this condition
  const sparql = `
    SELECT ?supplement ?supplementLabel ?activeIngredient ?activeIngredientLabel WHERE {
      ?supplement wdt:P31 wd:Q28885102.
      ?supplement wdt:P2175 wd:${conditionQid}.
      OPTIONAL { ?supplement wdt:P3781 ?activeIngredient. }
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    LIMIT 100
  `;

  const url = wdk.sparqlQuery(sparql);
  const response = await fetch(url);
  const data = await response.json();

  return data.results.bindings.map(result => ({
    supplement: result.supplementLabel.value,
    ingredient: result.activeIngredientLabel?.value
  }));
}

// Usage
getSupplementsByHealthClaim('cardiovascular disease').then(supplements => {
  supplements.forEach(supp => {
    console.log(`${supp.supplement}: ${supp.ingredient}`);
  });
});
```

---

## Part 5: Data Quality & Coverage

### 5.1 Completeness Analysis

| Category | Coverage | Quality | Recommendation |
|----------|----------|---------|----------------|
| **Vitamins** | ~50 entries | High | ✅ Good for vitamin data |
| **Minerals** | ~30 entries | High | ✅ Good for mineral data |
| **Herbal Supplements** | ~2,000 entries | Medium | ⚠️ Cross-validate plant links |
| **Amino Acids** | ~20 entries | High | ✅ Good for amino acid data |
| **Probiotics** | ~100 strains | Low | ❌ Use specialized databases |
| **Omega-3** | ~10 entries | Medium | ⚠️ Limited coverage |
| **Health Claims** | ~500 claims | Low | ❌ Not evidence-based |

### 5.2 Validation Strategy

1. **Vitamins & Minerals:** Wikidata has good coverage, but:
   - Cross-validate dosages with NIH Office of Dietary Supplements
   - Verify chemical IDs with PubChem

2. **Herbal Supplements:** Wikidata links are incomplete:
   - Use IMPPAT, Dr. Duke's for plant-compound mapping
   - Verify taxonomic names with World Flora Online

3. **Health Claims:** Wikidata health claims are NOT evidence-based:
   - Use ClinicalTrials.gov for clinical evidence
   - Use Cochrane Library for systematic reviews
   - Use Health Canada LNHPD for approved claims

---

## Part 6: Integration Workflow

### Recommended Pipeline

```
User Query (supplement name) → Wikidata Search API
                                      ↓
                         Wikidata SPARQL (get supplement data)
                                      ↓
                            Extract External IDs
                    ┌─────────┴─────────┬──────────┐
                    ↓                   ↓          ↓
          Health Canada LNHPD     PubChem API   Dr. Duke's DB
                    ↓                   ↓          ↓
          Approved claims       Chemical data   Plant source
                    └─────────┬─────────┘
                              ↓
                    Unified Supplement Profile
```

---

## Part 7: Known Limitations

1. **Limited Dosage Data:** Wikidata rarely includes:
   - Recommended Dietary Allowance (RDA)
   - Upper Intake Levels (UL)
   - Dosage forms

2. **Incomplete Probiotic Data:** Only ~100 strains vs. thousands in specialized databases

3. **No Adverse Event Data:** Wikidata lacks:
   - Side effects
   - Drug-supplement interactions
   - Safety warnings

4. **Unreliable Health Claims:** Many claims lack scientific evidence

**Recommendation:** Use Wikidata for:
- ✅ External ID mapping (especially Health Canada LNHPD)
- ✅ Chemical structure data
- ✅ Plant-compound linking (herbal supplements)

Use specialized databases for:
- ✅ Dosage recommendations (NIH ODS)
- ✅ Safety data (Natural Medicines Database)
- ✅ Clinical evidence (ClinicalTrials.gov)
- ✅ Drug-supplement interactions (DrugBank)

---

## Part 8: Notable Supplement Q-IDs

### Vitamins
- Q18225 - Vitamin A
- Q18229 - Vitamin C
- Q18225 - Vitamin D
- Q18225 - Vitamin E
- Q18225 - Vitamin K
- Q192111 - Vitamin B1 (Thiamine)
- Q192423 - Vitamin B2 (Riboflavin)
- Q192760 - Vitamin B3 (Niacin)
- Q164787 - Vitamin B5 (Pantothenic acid)
- Q18228 - Vitamin B7 (Biotin)
- Q165606 - Vitamin B9 (Folate)
- Q18225 - Vitamin B12 (Cobalamin)

### Minerals
- Q706 - Iron
- Q706 - Calcium
- Q706 - Magnesium
- Q706 - Zinc
- Q706 - Selenium
- Q706 - Copper
- Q706 - Manganese

### Popular Herbal Supplements
- Q421789 - Curcumin (from turmeric)
- Q18228 - Ginseng
- Q18228 - Echinacea
- Q18228 - St. John's Wort
- Q18228 - Ginkgo biloba
- Q18228 - Garlic extract

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Dietary Supplement | Product intended to supplement the diet containing vitamins, minerals, herbs, or other substances | Vitamin D capsules |
| Q-ID | Unique identifier for an item in Wikidata | Q18225 (Vitamin A) |
| Bioavailability | Proportion of a substance that enters circulation and has an active effect | 10-15% for curcumin |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| DSLD | Dietary Supplement Label Database (NIH) | Supplement labels |
| LNHPD | Licensed Natural Health Products Database (Health Canada) | Canadian supplements |
| Phytochemical | Bioactive compound produced by plants | Curcumin, quercetin |
| RDA | Recommended Dietary Allowance | Daily intake |
| UL | Tolerable Upper Intake Level | Safety threshold |
| Bioactive | Having a biological effect on living tissue | Active compounds |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST, SPARQL |
| CC0 | Creative Commons Zero (Public Domain) | Wikidata license |
| DSLD | Dietary Supplement Label Database | NIH resource |
| LNHPD | Licensed Natural Health Products Database | Health Canada |
| NIH | National Institutes of Health | US agency |
| ODS | Office of Dietary Supplements | Part of NIH |
| RDA | Recommended Dietary Allowance | Nutrition standard |
| SPARQL | SPARQL Protocol and RDF Query Language | Query language |
| UL | Tolerable Upper Intake Level | Safety limit |

---

*Guide compiled January 2026. For related supplement databases, see [alt-sources.md](./alt-sources.md)*
