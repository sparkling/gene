---
id: integration-wikidata-pharma
title: "Wikidata Pharmaceutical & Drug Data"
type: integration
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [integration, cross-references, apis]
---

**Parent:** [../_index.md](../_index.md)

# Wikidata Pharmaceutical & Drug Data

**Created:** January 2026
**Purpose:** Document Wikidata's pharmaceutical and drug-related data model and queries
**Status:** Active reference

---

## Overview

Wikidata contains structured data for:
- 15,000+ pharmaceutical drugs
- Drug-target relationships
- Drug-disease relationships
- Clinical trial links
- Regulatory approval data
- Drug-drug interactions
- Pharmacokinetics data

**Key Advantage:** CC0 license (public domain), SPARQL endpoint for complex queries

---

## Part 1: Wikidata Drug Data Model

### 1.1 Core Drug Properties

| Property | Label | Example | Use Case |
|----------|-------|---------|----------|
| **P31** | instance of | Q12140 (medication) | Identify drugs |
| **P661** | ChEMBL ID | CHEMBL308 | Link to ChEMBL |
| **P715** | DrugBank ID | DB00945 | Link to DrugBank |
| **P592** | ChEMBL target ID | CHEMBL1827 | Target mapping |
| **P2175** | medical condition treated | Q12203 (cancer) | Indications |
| **P3780** | active ingredient | Q compound | Formulation data |
| **P769** | significant drug interaction | Q other_drug | Safety |
| **P2789** | FDA approval date | 2020-01-15 | Regulatory |
| **P3345** | EMA number | EMEA/H/C/004463 | EU approval |

### 1.2 Chemical & Pharmacological Properties

| Property | Label | Example | Use Case |
|----------|-------|---------|----------|
| **P662** | PubChem CID | 969516 | Chemical ID |
| **P231** | CAS Registry Number | 458-37-7 | Chemical ID |
| **P234** | InChI | InChI=1S/... | Structure |
| **P235** | InChIKey | VFLDPWHFBUODDF-... | Structure hash |
| **P274** | chemical formula | C21H20O6 | Composition |
| **P2114** | pharmacokinetics | Q... | ADME data |
| **P3780** | active ingredient in | Q drug | Formulations |
| **P2868** | subject has role | Q12379712 (ATC code) | Classification |

### 1.3 Target & Mechanism Properties

| Property | Label | Example | Use Case |
|----------|-------|---------|----------|
| **P129** | physically interacts with | Q protein | Direct targets |
| **P3489** | protein complex | Q complex | Multi-target |
| **P2888** | exact match (URI) | UniProt URL | External link |
| **P680** | molecular function | GO term | Mechanism |

---

## Part 2: SPARQL Queries for Drug Data

### 2.1 Query: Get All FDA-Approved Drugs

```sparql
SELECT ?drug ?drugLabel ?drugbankID ?approvalDate WHERE {
  ?drug wdt:P31 wd:Q12140.              # instance of medication
  ?drug wdt:P715 ?drugbankID.           # has DrugBank ID

  OPTIONAL { ?drug wdt:P2789 ?approvalDate. }  # FDA approval date

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 10000
```

**Execute via API:**
```bash
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=SELECT ?drug ?drugLabel ?drugbankID WHERE { ?drug wdt:P31 wd:Q12140. ?drug wdt:P715 ?drugbankID. SERVICE wikibase:label { bd:serviceParam wikibase:language \"en\". } } LIMIT 10000" \
  -H "Accept: application/sparql-results+json"
```

---

### 2.2 Query: Get Drugs for a Specific Disease

```sparql
# Get all drugs approved for treating diabetes (Q12206)
SELECT ?drug ?drugLabel ?drugbankID ?chemblID WHERE {
  ?drug wdt:P31 wd:Q12140.              # instance of medication
  ?drug wdt:P2175 wd:Q12206.            # medical condition treated: diabetes

  OPTIONAL { ?drug wdt:P715 ?drugbankID. }
  OPTIONAL { ?drug wdt:P661 ?chemblID. }

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

### 2.3 Query: Get Drug Targets with UniProt IDs

```sparql
# Get protein targets for FDA-approved drugs
SELECT ?drug ?drugLabel ?target ?targetLabel ?uniprotID WHERE {
  ?drug wdt:P31 wd:Q12140.              # instance of medication
  ?drug wdt:P715 ?drugbankID.           # has DrugBank ID (ensures FDA approved)
  ?drug wdt:P129 ?target.               # physically interacts with target

  OPTIONAL { ?target wdt:P352 ?uniprotID. }  # UniProt ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

---

### 2.4 Query: Get Drug-Drug Interactions

```sparql
# Get known drug-drug interactions for a specific drug
SELECT ?drug1 ?drug1Label ?drug2 ?drug2Label ?interactionType WHERE {
  ?drug1 wdt:P31 wd:Q12140.             # instance of medication
  ?drug1 wdt:P769 ?drug2.               # significant drug interaction with drug2

  OPTIONAL { ?drug1 wdt:P31 ?interactionType. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

---

### 2.5 Query: Get Drugs with Specific Mechanism of Action

```sparql
# Get all kinase inhibitors
SELECT ?drug ?drugLabel ?drugbankID ?target ?targetLabel WHERE {
  ?drug wdt:P31 wd:Q12140.              # instance of medication
  ?drug wdt:P129 ?target.               # physically interacts with target
  ?target wdt:P31 wd:Q8047.             # target is an enzyme
  ?target rdfs:label ?targetLabel.

  FILTER(CONTAINS(LCASE(?targetLabel), "kinase"))
  FILTER(LANG(?targetLabel) = "en")

  OPTIONAL { ?drug wdt:P715 ?drugbankID. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

### 2.6 Query: Get Natural Product-Derived Drugs

```sparql
# Get drugs derived from natural products (plants, fungi, bacteria)
SELECT ?drug ?drugLabel ?drugbankID ?derivedFrom ?derivedFromLabel WHERE {
  ?drug wdt:P31 wd:Q12140.              # instance of medication
  ?drug wdt:P1343 ?source.              # described by source (often indicates natural origin)
  ?drug wdt:P703 ?derivedFrom.          # found in taxon (organism source)

  OPTIONAL { ?drug wdt:P715 ?drugbankID. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

---

### 2.7 Query: Get Drugs with Clinical Trial Links

```sparql
# Get drugs with ClinicalTrials.gov identifiers
SELECT ?drug ?drugLabel ?drugbankID ?trialID WHERE {
  ?drug wdt:P31 wd:Q12140.              # instance of medication
  ?drug wdt:P3098 ?trialID.             # ClinicalTrials.gov ID

  OPTIONAL { ?drug wdt:P715 ?drugbankID. }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

### 2.8 Query: Get ATC Classification

```sparql
# Get Anatomical Therapeutic Chemical (ATC) classification
SELECT ?drug ?drugLabel ?atcCode WHERE {
  ?drug wdt:P31 wd:Q12140.              # instance of medication
  ?drug wdt:P267 ?atcCode.              # ATC code

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 1000
```

**ATC Code Structure:**
- Level 1: Anatomical main group (A-V)
- Level 2: Therapeutic subgroup
- Level 3: Pharmacological subgroup
- Level 4: Chemical subgroup
- Level 5: Chemical substance

---

## Part 3: External Database ID Mapping

### 3.1 Available External IDs for Drugs

| Property | Database | Example | Coverage |
|----------|----------|---------|----------|
| **P715** | DrugBank | DB00945 | ~10,000 |
| **P661** | ChEMBL | CHEMBL308 | ~15,000 |
| **P662** | PubChem CID | 969516 | ~50,000 |
| **P683** | ChEBI | CHEBI:3962 | ~30,000 |
| **P592** | ChEMBL target ID | CHEMBL1827 | ~5,000 |
| **P3345** | EMA number | EMEA/H/C/004463 | ~1,000 |
| **P3098** | ClinicalTrials.gov ID | NCT01234567 | ~5,000 |
| **P2115** | RxNorm ID | 197361 | ~8,000 |
| **P3345** | ATC code | A10BA02 | ~5,000 |

### 3.2 Query: Get All External IDs for a Drug

```sparql
# Get all external database IDs for Metformin (Q19484)
SELECT ?property ?propertyLabel ?value WHERE {
  wd:Q19484 ?prop ?value.
  ?property wikibase:directClaim ?prop.
  ?property wikibase:propertyType wikibase:ExternalId.

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
```

---

## Part 4: Integration Examples

### 4.1 Python Example: Query Drugs for a Disease

```python
import requests

def get_drugs_for_disease(disease_name):
    """
    Get FDA-approved drugs for a specific disease via Wikidata
    """
    # Step 1: Search Wikidata for the disease
    search_url = "https://www.wikidata.org/w/api.php"
    search_params = {
        'action': 'wbsearchentities',
        'search': disease_name,
        'language': 'en',
        'format': 'json',
        'type': 'item'
    }

    search_resp = requests.get(search_url, params=search_params)
    disease_results = search_resp.json()['search']

    if not disease_results:
        return None

    disease_qid = disease_results[0]['id']
    print(f"Disease QID: {disease_qid}")

    # Step 2: Query drugs for this disease
    sparql_query = f"""
    SELECT ?drug ?drugLabel ?drugbankID ?chemblID WHERE {{
      ?drug wdt:P31 wd:Q12140.
      ?drug wdt:P2175 wd:{disease_qid}.
      OPTIONAL {{ ?drug wdt:P715 ?drugbankID. }}
      OPTIONAL {{ ?drug wdt:P661 ?chemblID. }}
      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
    }}
    """

    sparql_url = "https://query.wikidata.org/sparql"
    sparql_resp = requests.get(sparql_url, params={
        'query': sparql_query,
        'format': 'json'
    })

    results = sparql_resp.json()['results']['bindings']

    drugs = []
    for result in results:
        drugs.append({
            'name': result['drugLabel']['value'],
            'wikidata_id': result['drug']['value'].split('/')[-1],
            'drugbank_id': result.get('drugbankID', {}).get('value'),
            'chembl_id': result.get('chemblID', {}).get('value')
        })

    return drugs

# Usage
drugs = get_drugs_for_disease("diabetes mellitus")
for drug in drugs:
    print(f"{drug['name']}: DrugBank {drug['drugbank_id']}, ChEMBL {drug['chembl_id']}")
```

---

### 4.2 Python Example: Get Drug Targets

```python
def get_drug_targets(drug_name):
    """
    Get protein targets for a drug via Wikidata
    """
    # Step 1: Search for the drug
    search_url = "https://www.wikidata.org/w/api.php"
    search_params = {
        'action': 'wbsearchentities',
        'search': drug_name,
        'language': 'en',
        'format': 'json',
        'type': 'item'
    }

    search_resp = requests.get(search_url, params=search_params)
    drug_results = search_resp.json()['search']

    if not drug_results:
        return None

    drug_qid = drug_results[0]['id']

    # Step 2: Query targets
    sparql_query = f"""
    SELECT ?target ?targetLabel ?uniprotID ?geneSymbol WHERE {{
      wd:{drug_qid} wdt:P129 ?target.
      OPTIONAL {{ ?target wdt:P352 ?uniprotID. }}
      OPTIONAL {{ ?target wdt:P353 ?geneSymbol. }}
      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
    }}
    """

    sparql_url = "https://query.wikidata.org/sparql"
    sparql_resp = requests.get(sparql_url, params={
        'query': sparql_query,
        'format': 'json'
    })

    results = sparql_resp.json()['results']['bindings']

    targets = []
    for result in results:
        targets.append({
            'name': result['targetLabel']['value'],
            'uniprot_id': result.get('uniprotID', {}).get('value'),
            'gene_symbol': result.get('geneSymbol', {}).get('value')
        })

    return targets

# Usage
targets = get_drug_targets("metformin")
for target in targets:
    print(f"{target['name']} ({target['gene_symbol']}): {target['uniprot_id']}")
```

---

### 4.3 JavaScript Example: Query via wikidata-sdk

```javascript
const wdk = require('wikidata-sdk');
const fetch = require('node-fetch');

async function getDrugsByATC(atcCode) {
  // SPARQL query
  const sparql = `
    SELECT ?drug ?drugLabel ?drugbankID WHERE {
      ?drug wdt:P31 wd:Q12140.
      ?drug wdt:P267 "${atcCode}".
      OPTIONAL { ?drug wdt:P715 ?drugbankID. }
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
  `;

  const url = wdk.sparqlQuery(sparql);
  const response = await fetch(url);
  const data = await response.json();

  return data.results.bindings.map(result => ({
    name: result.drugLabel.value,
    wikidata_id: result.drug.value.split('/').pop(),
    drugbank_id: result.drugbankID?.value
  }));
}

// Usage
getDrugsByATC('A10BA02').then(drugs => {
  drugs.forEach(drug => {
    console.log(`${drug.name}: ${drug.drugbank_id}`);
  });
});
```

---

## Part 5: Data Quality & Coverage

### 5.1 Completeness Analysis

| Data Type | Coverage | Quality | Recommendation |
|-----------|----------|---------|----------------|
| **DrugBank IDs** | ~10,000 drugs | High | ✅ Use for mapping |
| **ChEMBL IDs** | ~15,000 drugs | High | ✅ Use for mapping |
| **Drug-Target Links** | ~5,000 drugs | Medium | ⚠️ Cross-validate |
| **Drug-Disease Links** | ~8,000 drugs | Medium | ⚠️ Cross-validate |
| **Drug-Drug Interactions** | ~2,000 pairs | Low | ❌ Use DrugBank instead |
| **Clinical Trial IDs** | ~5,000 drugs | Low | ⚠️ Partial coverage |
| **ATC Codes** | ~5,000 drugs | High | ✅ Use for classification |

### 5.2 Validation Strategy

1. **Always cross-reference** Wikidata drug-target claims with:
   - ChEMBL (definitive source)
   - DrugBank (curated data)
   - BindingDB (measured affinities)

2. **Check references** in Wikidata claims:
   - Look for `P248` (stated in) property
   - Prioritize claims with PubMed references

3. **Use external IDs** to fetch detailed data:
   - DrugBank ID → Full drug profile
   - ChEMBL ID → Target interactions + bioassays
   - PubChem CID → Chemical structure + bioactivity

---

## Part 6: Integration Workflow

### Recommended Pipeline

```
User Query (disease name) → Wikidata Search API (get disease QID)
                                      ↓
                         Wikidata SPARQL (get drugs for disease)
                                      ↓
                            Extract External IDs
                    ┌─────────┴─────────┬─────────┐
                    ↓                   ↓         ↓
              DrugBank API         ChEMBL API   PubChem API
                    ↓                   ↓         ↓
             Full drug data      Target data   Chemical data
                    └─────────┬─────────┘
                              ↓
                    Unified Drug-Target-Disease Model
```

**Key Insight:** Use Wikidata as initial discovery layer, then fetch detailed data from authoritative sources.

---

## Part 7: License & Usage

- **Wikidata Data:** CC0 (public domain)
- **SPARQL Endpoint:** Free, no API key required
- **Rate Limits:** 60-second timeout per query
- **Attribution:** Not required (CC0), but recommended

---

## Part 8: Known Limitations

1. **Incomplete Target Data:** Wikidata has ~5,000 drug-target relationships, compared to:
   - ChEMBL: ~2.3 million
   - DrugBank: ~30,000

2. **Limited Pharmacokinetics:** Wikidata has minimal ADME data

3. **Few Drug-Drug Interactions:** ~2,000 pairs vs. DrugBank's ~1 million

4. **Outdated Clinical Trials:** ClinicalTrials.gov IDs may not be current

5. **No Quantitative Binding Data:** Use BindingDB or ChEMBL for Ki/Kd values

**Recommendation:** Use Wikidata for:
- ✅ External ID mapping
- ✅ Disease-drug discovery
- ✅ ATC classification
- ✅ Regulatory status

Use ChEMBL/DrugBank for:
- ✅ Comprehensive target data
- ✅ Binding affinities
- ✅ Pharmacokinetics
- ✅ Drug-drug interactions

---

## Download

### Wikidata Pharmaceutical Data Sources

| Source | URL | Format | Access | Size |
|--------|-----|--------|--------|------|
| Complete Wikidata dump | https://dumps.wikimedia.org/wikidatawiki/entities/ | JSON Lines (.bz2) | Bulk | ~90 GB |
| Wikidata SPARQL endpoint | https://query.wikidata.org/sparql | SPARQL/JSON | Query | Real-time |
| ChEMBL data | https://www.ebi.ac.uk/chembl/ | PostgreSQL/RDF | API/Bulk | ~2 GB |
| DrugBank | https://www.drugbank.com/ | XML/RDF | Registration | ~500 MB |

### SPARQL Query for Pharmaceutical Data

```bash
# Query Wikidata for drugs with ChEMBL IDs
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=
    SELECT ?drug ?drugLabel ?chembl ?drugbank WHERE {
      ?drug wdt:P31 wd:Q12140 .
      OPTIONAL { ?drug wdt:P661 ?chembl }
      OPTIONAL { ?drug wdt:P715 ?drugbank }
      SERVICE wikibase:label { bd:serviceParam wikibase:language \"en\" }
    }
    LIMIT 100" \
  -H "Accept: application/sparql-results+json"
```

---

## Data Format

### Pharmaceutical Data in Wikidata

**JSON Lines format from dump:**

```json
{
  "type": "item",
  "id": "Q12261",
  "labels": { "en": { "value": "ibuprofen" } },
  "claims": {
    "P31": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P31", "datavalue": { "value": { "entity-type": "item", "numeric-id": 12140 }, "type": "wikibase-entityid" } } }],
    "P715": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P715", "datavalue": { "value": "DB01050", "type": "string" } } }],
    "P661": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P661", "datavalue": { "value": "CHEMBL364", "type": "string" } } }],
    "P2175": [{ "rank": "normal", "mainsnak": { "snaktype": "value", "property": "P2175", "datavalue": { "value": { "entity-type": "item", "numeric-id": 4084 }, "type": "wikibase-entityid" } } }]
  }
}
```

**SPARQL JSON response:**

```json
{
  "results": {
    "bindings": [
      {
        "drug": { "type": "uri", "value": "http://www.wikidata.org/entity/Q12261" },
        "drugLabel": { "type": "literal", "value": "ibuprofen", "xml:lang": "en" },
        "chembl": { "type": "literal", "value": "CHEMBL364" },
        "drugbank": { "type": "literal", "value": "DB01050" }
      }
    ]
  }
}
```

---

## Schema

### Pharmaceutical Drug Schema

| Property | P-code | Type | Example | Use Case |
|----------|--------|------|---------|----------|
| Instance of | P31 | Item | Q12140 | Drug classification |
| ChEMBL ID | P661 | String | CHEMBL364 | Chemical structure mapping |
| DrugBank ID | P715 | String | DB01050 | Pharmacology data |
| Medical condition treated | P2175 | Item | Q4084 (pain) | Therapeutic indication |
| Active ingredient | P3780 | Item | Q compound | Formulation |
| Drug interaction | P769 | Item | Q other_drug | Safety interaction |
| FDA approval date | P2789 | Date | 1974-04-12 | Regulatory timeline |
| Route of administration | P636 | Item | Q12322 (oral) | Dosage form |
| ATC code | P3781 | String | M01AE01 | Anatomical classification |

### Drug Properties Available in Wikidata

| Property | Description | Typical Values |
|----------|-------------|-----------------|
| **P661** | ChEMBL ID | CHEMBL prefix + numeric ID |
| **P715** | DrugBank ID | DB followed by 5 digits |
| **P592** | ChEMBL target | Target molecule identifier |
| **P2175** | Treats condition | Q-ID of disease/symptom |
| **P3345** | EMA number | EMEA/H/C/XXXXX format |
| **P3781** | ATC code | 5-character alphanumeric |

---

## Sample Data

### Sample Drug Record: Ibuprofen

**SPARQL Query:**
```sparql
SELECT ?drug ?label ?chembl ?drugbank ?indication ?atp WHERE {
  ?drug wdt:P31 wd:Q12140 .
  ?drug rdfs:label "ibuprofen"@en .
  OPTIONAL { ?drug wdt:P661 ?chembl }
  OPTIONAL { ?drug wdt:P715 ?drugbank }
  OPTIONAL { ?drug wdt:P2175 ?indication }
  OPTIONAL { ?drug wdt:P3781 ?atp }
}
```

**Sample Results:**
```json
{
  "drug": { "value": "http://www.wikidata.org/entity/Q12261" },
  "label": { "value": "ibuprofen", "xml:lang": "en" },
  "chembl": { "value": "CHEMBL364" },
  "drugbank": { "value": "DB01050" },
  "indication": { "value": "http://www.wikidata.org/entity/Q4084" },
  "atp": { "value": "M01AE01" }
}
```

### Sample Drug List (Top 5 by Claims)

| Q-ID | Drug Name | ChEMBL | DrugBank | Indications |
|------|-----------|--------|----------|-------------|
| Q12261 | Ibuprofen | CHEMBL364 | DB01050 | Pain, fever, inflammation |
| Q18216 | Aspirin | CHEMBL25 | DB00945 | Pain, anticoagulation |
| Q27074 | Metformin | CHEMBL1642 | DB00331 | Type 2 diabetes |
| Q12312 | Paracetamol | CHEMBL112 | DB00316 | Pain, fever |
| Q12374 | Naproxen | CHEMBL1232 | DB00788 | Pain, inflammation |

---

## License

### Wikidata Licensing

- **License:** CC0 1.0 Universal (Public Domain)
- **Requirement:** No attribution required; optional but encouraged
- **DrugBank Link:** DrugBank data is CC-BY-NC 4.0 (not public domain)
- **ChEMBL Link:** ChEMBL data is CC-BY 4.0 (requires attribution)
- **Caveat:** When using linked external databases, follow their respective licenses

### Citation Format

```
Wikidata contributors. "Wikidata." Wikimedia Foundation, Inc., 2024.
https://www.wikidata.org/

For linked data, also cite:
- ChEMBL: Gaulton et al. (2017) Nucleic Acids Res 45(D1):D945-D954
- DrugBank: Wishart et al. (2018) Nucleic Acids Res 46(D1):D1074-D1082
```

---

## Data Set Size

### Pharmaceutical Dataset Statistics

| Metric | Count | Notes |
|--------|-------|-------|
| Medications (Q12140) | ~45,000 | Primary drug class |
| Drugs with ChEMBL ID | ~15,000 | Linked to ChEMBL |
| Drugs with DrugBank ID | ~10,000 | Linked to DrugBank |
| Drug-disease relationships | ~50,000 | Via P2175 property |
| Drug interactions documented | ~5,000 | Via P769 property |
| ATC classifications | ~30,000 | Via P3781 property |

### Storage and Access

| Metric | Value |
|--------|-------|
| Wikidata dump size | ~90 GB compressed |
| Filtered pharma subset | ~500 MB - 1 GB |
| SPARQL query response time | <1 second |
| Last dump update | Weekly (Monday UTC) |
| SPARQL endpoint lag | <1 minute |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Drug Target | Molecular entity (usually protein) that a drug binds to exert its effect | P04637 (p53) |
| ATC Code | Anatomical Therapeutic Chemical classification for drugs | A10BA02 (metformin) |
| Pharmacokinetics | Study of drug absorption, distribution, metabolism, excretion (ADME) | Half-life, bioavailability |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| ChEMBL ID | Identifier for compounds in the ChEMBL database | CHEMBL941 |
| DrugBank ID | Identifier for drugs in the DrugBank database | DB00331 |
| PubChem CID | Compound identifier in PubChem | CID 2244 |
| InChIKey | Hashed chemical structure identifier | 27-character string |
| ATC Classification | WHO drug classification system | 5-level hierarchy |
| DDI | Drug-Drug Interaction | Safety concern |
| ADME | Absorption, Distribution, Metabolism, Excretion | Pharmacokinetics |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST, SPARQL |
| ATC | Anatomical Therapeutic Chemical | WHO drug classification |
| CC0 | Creative Commons Zero (Public Domain) | Wikidata license |
| CID | Compound Identifier | PubChem ID type |
| DDI | Drug-Drug Interaction | Safety data |
| SPARQL | SPARQL Protocol and RDF Query Language | Query language |

---

*Guide compiled January 2026. For related pharmaceutical data sources, see [../data-sources.md](../data-sources.md)*
