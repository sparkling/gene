---
id: schemas-lotus-schema
title: "LOTUS Natural Products Database Schema"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, natural-products, wikidata, sparql, structure-organism, cc0]
---

**Parent:** [Schema Documentation](./_index.md)

# LOTUS Natural Products Database Schema

**Document ID:** LOTUS-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** LOTUS via Wikidata (Continuous Updates)

---

## TL;DR

LOTUS (Linked Open Total Unified Structure-organism) contains 750K+ referenced structure-organism pairs linking natural product chemical structures to their biological source organisms with literature citations. Data is fully integrated into Wikidata using SPARQL queries. Core properties: InChIKey (P235), found in taxon (P703), stated in (P248). The dataset represents 290K+ unique structures from 40K+ organisms with 75K+ references. License: CC0 (Creative Commons Zero).

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **Referenced Structure-Organism Pairs** | 750,000+ |
| **Unique Chemical Structures** | 290,000+ |
| **Distinct Organisms (Taxa)** | 40,000+ |
| **Literature References** | 75,000+ |
| **Initial Entries (Pre-validation)** | 2,500,000+ |
| **Validation Retention Rate** | ~30% |
| **Manual Validation Accuracy** | 97% true positives |

---

## Available Downloads

### Wikidata Export Tables

| File | Description |
|------|-------------|
| compounds.tsv | Chemical structures metadata (wikidataId, SMILES, InChI, InChIKey) |
| references.tsv | Bibliographical references (wikidataId, DOI, title, PMID) |
| taxa.tsv | Biological organisms (wikidataId, name, rank, taxonomy) |
| compound_reference_taxon.tsv | Documented structure-organism pairs |

**Download Sources:**
- Zenodo Archive: https://doi.org/10.5281/zenodo.5794106
- Wikidata Exporter: https://github.com/lotusnprod/lotus-wikidata-exporter
- Query Results: https://zenodo.org/communities/the-lotus-initiative/

**License:** CC0 (Creative Commons Zero) - completely unrestricted

---

## Data Architecture

### Dual Hosting Model

```
Wikidata (Primary)                    LNPN Interface
    |                                      |
    +-- SPARQL Query Service               +-- Simple Search
    |                                      |
    +-- RDF Exports                        +-- Structure Search
    |                                      |
    +-- Federated Queries                  +-- Advanced Filters
                                           |
                                           +-- API Access
```

### Single Source of Truth (SSOT)

```
PostgreSQL (Canonical)
       |
       +---> Wikidata (Knowledge Graph)
       |
       +---> LNPN (lotus.naturalproducts.net)
```

---

## Core Entity Relationships

### Structure-Organism-Reference Triplet

```
CHEMICAL_STRUCTURE
       |
       | (P703: found in taxon)
       v
    ORGANISM <-------- REFERENCE
                       (P248: stated in)
```

### Three Minimal Sufficient Objects

1. **Chemical Structure Object:**
   - SMILES (canonical and isomeric)
   - InChI string
   - InChIKey (27-character hash)

2. **Biological Organism Object:**
   - Taxon name (scientific)
   - Taxonomic database reference
   - Taxon ID (NCBI, OTT, GBIF)

3. **Reference Object:**
   - Article title
   - DOI / PMID / PMCID
   - Publication metadata

---

## Wikidata Property Mappings

### Core Properties

| Property ID | Property Name | Description |
|-------------|---------------|-------------|
| **P235** | InChIKey | 27-character hashed chemical identifier |
| **P234** | InChI | Full InChI string |
| **P233** | Canonical SMILES | Canonical SMILES notation |
| **P2017** | Isomeric SMILES | SMILES with stereochemistry |
| **P703** | found in taxon | Links compound to source organism |
| **P248** | stated in | Bibliographic reference for claim |
| **P31** | instance of | Type classification |

### Organism Properties

| Property ID | Property Name | Description |
|-------------|---------------|-------------|
| **P171** | parent taxon | Taxonomic parent |
| **P105** | taxon rank | Species, genus, family, etc. |
| **P225** | taxon name | Scientific name |
| **P685** | NCBI taxon ID | NCBI Taxonomy identifier |
| **P9157** | Open Tree of Life ID | OTT identifier (created for LOTUS) |

### Reference Properties

| Property ID | Property Name | Description |
|-------------|---------------|-------------|
| **P356** | DOI | Digital Object Identifier |
| **P698** | PubMed ID | PubMed article identifier |
| **P1476** | title | Publication title |
| **P577** | publication date | Date published |
| **P2093** | author name string | Author names |

---

## Data Schema Tables

### compounds.tsv

| Column | Type | Description |
|--------|------|-------------|
| wikidataId | VARCHAR | Wikidata Q-identifier (e.g., Q312879) |
| canonicalSmiles | TEXT | Canonical SMILES notation |
| isomericSmiles | TEXT | Isomeric SMILES with stereochemistry |
| inchi | TEXT | Standard InChI string |
| inchiKey | VARCHAR(27) | InChI key |

---

### taxa.tsv

| Column | Type | Description |
|--------|------|-------------|
| wikidataId | VARCHAR | Wikidata Q-identifier for organism |
| name | VARCHAR | Scientific name |
| rank | VARCHAR | Taxonomic rank (species, genus, etc.) |
| ncbiId | INTEGER | NCBI Taxonomy ID |
| ottId | INTEGER | Open Tree of Life ID |

---

### references.tsv

| Column | Type | Description |
|--------|------|-------------|
| wikidataId | VARCHAR | Wikidata Q-identifier for reference |
| doi | VARCHAR | Digital Object Identifier |
| pmid | INTEGER | PubMed ID |
| pmcid | VARCHAR | PubMed Central ID |
| title | TEXT | Publication title |

---

### compound_reference_taxon.tsv

| Column | Type | Description |
|--------|------|-------------|
| compoundWikidataId | VARCHAR | FK to compounds |
| taxonWikidataId | VARCHAR | FK to taxa |
| referenceWikidataId | VARCHAR | FK to references |

---

## SPARQL Query Examples

### Basic: All Referenced Structure-Organism Pairs

```sparql
SELECT DISTINCT ?compound ?compoundLabel ?taxon ?taxonLabel ?reference ?referenceLabel
WHERE {
  ?compound p:P235 [];           # Has InChIKey
           p:P703 [              # Found in taxon statement
             ps:P703 ?taxon;
             prov:wasDerivedFrom/pr:P248 ?reference;
             wikibase:rank wikibase:NormalRank;
           ] .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 1000
```

### Compounds in a Specific Organism

```sparql
# Find all compounds in Arabidopsis thaliana
SELECT DISTINCT ?compound ?compoundLabel ?inchikey
WHERE {
  ?compound wdt:P235 ?inchikey;
            wdt:P703 wd:Q158695.  # Arabidopsis thaliana
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
```

**Live Query:** https://w.wiki/4Vcv

### Organisms Containing a Specific Compound

```sparql
# Find all organisms containing beta-sitosterol
SELECT DISTINCT ?taxon ?taxonLabel
WHERE {
  wd:Q312879 wdt:P703 ?taxon.    # beta-sitosterol (Q312879)
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
```

**Live Query:** https://w.wiki/4VFn

### Compounds by Scaffold/Substructure

```sparql
# Organisms producing compounds with indolic scaffold
SELECT ?parentTaxon ?parentTaxonLabel (COUNT(DISTINCT ?compound) AS ?count)
WHERE {
  ?compound wdt:P235 ?ik;
            wdt:P703 ?taxon;
            wdt:P2017 ?smiles.
  ?taxon wdt:P171 ?parentTaxon.
  FILTER(CONTAINS(?smiles, "c1ccc2[nH]ccc2c1"))  # Indole pattern
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
GROUP BY ?parentTaxon ?parentTaxonLabel
ORDER BY DESC(?count)
```

**Live Query:** https://w.wiki/4VG9

### Bioactive Compounds from Specific Taxa

```sparql
# Bioactive compounds from Actinobacteria (2014-2019)
SELECT DISTINCT ?compound ?compoundLabel ?taxon ?taxonLabel ?reference ?date
WHERE {
  ?compound wdt:P235 ?ik;
            wdt:P703 ?taxon;
            wdt:P2868 ?bioactivity.  # subject has role
  ?taxon wdt:P171* wd:Q26262.        # Actinobacteria
  ?compound p:P703 [
    prov:wasDerivedFrom [
      pr:P248 ?reference;
      pr:P577 ?date.
    ]
  ].
  FILTER(YEAR(?date) >= 2014 && YEAR(?date) <= 2019)
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 500
```

**Live Query:** https://w.wiki/4VGC

### Terpenoids from a Genus

```sparql
# Terpenoids from Aspergillus (2010-2020)
SELECT DISTINCT ?compound ?compoundLabel ?taxon ?taxonLabel
WHERE {
  ?compound wdt:P235 ?ik;
            wdt:P703 ?taxon;
            wdt:P279* wd:Q422292.    # subclass of terpenoid
  ?taxon wdt:P171* wd:Q335130.       # Aspergillus
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
```

**Live Query:** https://w.wiki/4VGD

### Cross-Reference with ChEMBL

```sparql
# LOTUS compounds with ChEMBL IDs
SELECT ?compound ?compoundLabel ?inchikey ?chemblId
WHERE {
  ?compound wdt:P235 ?inchikey;
            wdt:P703 [];            # Has organism link (LOTUS)
            wdt:P592 ?chemblId.     # Has ChEMBL ID
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 500
```

---

## IDSM/Sachem Substructure Search

For substructure and similarity searches, use the IDSM SPARQL endpoint:

**Endpoint:** https://idsm.elixir-czech.cz/sparql/endpoint/wikidata

### Substructure Search Example

```sparql
PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
PREFIX idsm: <https://idsm.elixir-czech.cz/sparql/endpoint/wikidata>

SELECT ?compound ?inchikey
WHERE {
  SERVICE idsm: {
    ?compound sachem:substructureSearch [
      sachem:query "c1ccc2[nH]ccc2c1"  # Indole SMILES
    ].
  }
  ?compound wdt:P235 ?inchikey;
            wdt:P703 [].  # LOTUS compounds only
}
LIMIT 100
```

---

## MongoDB Schema (LNPN Instance)

### Collection: lotusUniqueNaturalProduct

| Field | Type | Description | Indexed |
|-------|------|-------------|---------|
| lotus_id | String | Internal identifier | Yes |
| inchi | String | Standard InChI | Yes |
| inchikey | String | InChI key | Yes |
| smiles | String | SMILES notation | Yes |
| inchi2D | String | 2D InChI (no stereochemistry) | Yes |
| inchikey2D | String | 2D InChI key | Yes |
| smiles2D | String | 2D SMILES | Yes |
| molecular_formula | String | Molecular formula | Yes |
| molecular_weight | Float | Molecular weight | Yes |
| heavy_atom_number | Integer | Heavy atom count | Yes |
| npl_score | Float | Natural product likeness | Yes |
| fsp3 | Float | Fraction sp3 carbons | Yes |
| lipinskiRuleOf5Failures | Integer | Lipinski violations | Yes |
| deep_smiles | String | DeepSMILES representation | Yes |
| iupac_name | String | IUPAC systematic name | Yes (text) |
| traditional_name | String | Common/traditional name | Yes (text) |
| allTaxa | Array[String] | All associated organisms | Yes (text) |
| allChemClassifications | Array[String] | Chemical classifications | Yes (text) |
| allWikidataIds | Array[String] | All Wikidata Q-IDs | Yes (text) |
| fragments | Array | Fragment analysis | Yes |
| fragmentsWithSugar | Array | Fragments with sugars | Yes |
| pfCounts | Object | Physicochemical counts | Yes |
| pubchemBits | Array | PubChem fingerprint bits | Yes |

---

## Sample Data

### Structure-Organism Pair (TSV Format)

```
compoundWikidataId	taxonWikidataId	referenceWikidataId
Q312879	Q158695	Q56636231
Q312879	Q147098	Q28189571
Q27104933	Q147098	Q28189571
```

### Compound Record (JSON)

```json
{
  "wikidataId": "Q312879",
  "name": "beta-Sitosterol",
  "canonicalSmiles": "CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C",
  "isomericSmiles": "CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C",
  "inchi": "InChI=1S/C29H50O/c1-7-21(19(2)3)8-9-22(20(4)5)25-12-13-26-24-11-10-23-18-27(30)14-16-28(23,6)25-15-17-29(25,26)6/h10,19-22,24-27,30H,7-9,11-18H2,1-6H3/t21-,22-,24+,25-,26+,27+,28+,29-/m1/s1",
  "inchiKey": "KZJWDPNRJALLNS-VJSFXXLFSA-N",
  "molecularFormula": "C29H50O",
  "molecularWeight": 414.71,
  "foundInTaxa": [
    {"wikidataId": "Q158695", "name": "Arabidopsis thaliana"},
    {"wikidataId": "Q147098", "name": "Allium cepa"}
  ],
  "references": [
    {"wikidataId": "Q56636231", "doi": "10.1016/j.phytochem.2018.01.001"}
  ]
}
```

### Organism Record (JSON)

```json
{
  "wikidataId": "Q158695",
  "name": "Arabidopsis thaliana",
  "rank": "species",
  "parentTaxon": "Arabidopsis",
  "ncbiId": 3702,
  "ottId": 1054861,
  "compoundCount": 1247
}
```

---

## Cross-References to Other Databases

### Chemical Identifier Mapping

| LOTUS Field | Cross-Reference | Database |
|-------------|-----------------|----------|
| inchiKey | standard_inchi_key | COCONUT |
| inchiKey | standard_inchi_key | ChEMBL |
| inchiKey | InChIKey | PubChem |
| wikidataId | wikidata_iri | COCONUT Organisms |
| ncbiId | tax_id | ChEMBL Targets |
| ncbiId | ncbi_id | COCONUT Organisms |

### COCONUT Cross-Reference

COCONUT and LOTUS share many compounds. To find overlapping records:

```sql
-- COCONUT PostgreSQL
SELECT c.identifier, c.standard_inchi_key, o.name AS organism
FROM molecules c
JOIN organism_molecule om ON c.id = om.molecule_id
JOIN organisms o ON om.organism_id = o.id
WHERE c.standard_inchi_key IN (
  -- InChI keys from LOTUS query
);
```

### ChEMBL Cross-Reference

```sql
-- ChEMBL: Find bioactivity for LOTUS compounds
SELECT md.chembl_id, cs.standard_inchi_key,
       a.standard_type, a.standard_value, a.standard_units
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities a ON md.molregno = a.molregno
WHERE md.natural_product = 1
AND cs.standard_inchi_key IN (
  -- InChI keys from LOTUS
);
```

### Federated Wikidata Query (ChEMBL + LOTUS)

```sparql
SELECT ?compound ?inchikey ?chemblId ?taxon ?taxonLabel
WHERE {
  ?compound wdt:P235 ?inchikey;
            wdt:P703 ?taxon;        # LOTUS: found in taxon
            wdt:P592 ?chemblId.     # ChEMBL ID
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 500
```

---

## Data Quality Indicators

### Validation Pipeline

| Stage | Description | Outcome |
|-------|-------------|---------|
| **Harmonization** | Standardize formats | Unified identifiers |
| **Processing** | Compute properties | InChI, SMILES normalization |
| **Validation** | Check against sources | 97% true positive rate |
| **Dissemination** | Publish to Wikidata | CC0 licensed |

### Quality Metrics

- **Structure Validation:** Must resolve to valid InChIKey
- **Organism Validation:** Must map to canonical taxonomic name
- **Reference Validation:** Verified against original publications
- **Retention Rate:** ~30% after validation (2.5M -> 750K)

---

## API Access

### LNPN REST API

**Base URL:** https://lotus.naturalproducts.net/api

| Endpoint | Method | Description |
|----------|--------|-------------|
| /search | GET | Simple text search |
| /structure | POST | Structure search (exact/similarity/substructure) |
| /compound/{id} | GET | Get compound details |
| /organism/{id} | GET | Get organism details |

### Wikidata Query Service

**SPARQL Endpoint:** https://query.wikidata.org/sparql

**Query Interface:** https://query.wikidata.org/

**Rate Limits:** 60 requests/minute for automated queries

---

## Integration Notes for Knowledge Base

### Key Identifiers for Cross-Referencing

| Identifier | Format | Example | Cross-ref to |
|------------|--------|---------|--------------|
| Wikidata QID | Q + digits | Q312879 | All Wikidata-linked DBs |
| InChI Key | 27 characters | KZJWDPNRJALLNS-VJSFXXLFSA-N | PubChem, ChEMBL, COCONUT |
| NCBI Taxon ID | Integer | 3702 | UniProt, ChEMBL, COCONUT |
| DOI | doi:... | 10.1016/j.phytochem.2018.01.001 | CrossRef, PubMed |
| PubMed ID | Integer | 29371045 | PubMed, Europe PMC |
| Open Tree of Life ID | Integer | 1054861 | OTL |

### Licensing

**License:** CC0 (Creative Commons Zero)
**Attribution Required:** No (but appreciated)
**Commercial Use:** Unrestricted
**Derivative Works:** Unrestricted

---

## Useful Links

| Resource | URL |
|----------|-----|
| LOTUS Web Interface | https://lotus.naturalproducts.net/ |
| Wikidata Project Page | https://www.wikidata.org/wiki/Wikidata:WikiProject_Chemistry/Natural_products |
| SPARQL Query Interface | https://query.wikidata.org/ |
| GitHub: lotus-web | https://github.com/lotusnprod/lotus-web |
| GitHub: lotus-wikidata-interact | https://github.com/lotusnprod/lotus-wikidata-interact |
| Zenodo Archive | https://doi.org/10.5281/zenodo.5794106 |
| IDSM SPARQL Endpoint | https://idsm.elixir-czech.cz/sparql/endpoint/wikidata |
| Publication (eLife) | https://elifesciences.org/articles/70780 |

---

## Download

### Data Access Methods

| Method | URL | Format |
|--------|-----|--------|
| **Web Interface** | https://lotus.naturalproducts.net/ | Interactive search |
| **SPARQL Query** | https://query.wikidata.org/ | SPARQL queries |
| **Zenodo Archive** | https://doi.org/10.5281/zenodo.5794106 | TSV, JSON export |
| **GitHub Downloads** | https://github.com/lotusnprod/lotus-web | Source code/exports |
| **Wikidata Dump** | https://www.wikidata.org/wiki/Wikidata:Data_access | RDF, JSON |

### SPARQL Example

```bash
# Query natural products from LOTUS
curl "https://query.wikidata.org/sparql?query=SELECT%20...&format=json"
```

---

## Data Format

| Format | Description | Use Case |
|--------|-------------|----------|
| Primary | RDF/Linked Data (Wikidata) | Semantic web queries |
| Alternative | TSV (Bulk export) | Data analysis |
| Alternative | JSON (API responses) | Programmatic access |
| Query | SPARQL | Complex queries |
| Encoding | UTF-8 | All formats |

---

## Data Set Size

| Component | Count | Size (Est.) |
|-----------|-------|------------|
| **Natural Products** | 750,000+ | ~100-150 MB |
| **Organisms** | 100,000+ | ~50 MB |
| **References (DOIs)** | 50,000+ | ~10 MB |
| **Structure Data** | All linked | Included |
| **Total Zenodo Archive** | Full export | ~500-800 MB |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "LOTUS0000001" |
| `name` | string | Entity name | "Artemisinin" |
| `type` | string | Record type | "natural_product" / "compound" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `found_in` | Organism | N:M |
| `has_structure` | Chemical Structure | 1:1 |
| `associated_with` | Activity | N:M |

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| LOTUS | CC0 (Public Domain) | Yes |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `wikidataId` | Wikidata Q-identifier for compounds, organisms, or references | Q312879 (beta-sitosterol) |
| `canonicalSmiles` | Simplified molecular-input line-entry system notation without stereochemistry | CC(C)CCCC(C)C |
| `isomericSmiles` | SMILES notation with stereochemistry information preserved | CC[C@H](CC[C@@H](C)... |
| `inchiKey` | 27-character hashed InChI for fast structure lookup | KZJWDPNRJALLNS-VJSFXXLFSA-N |
| `P703` | Wikidata property "found in taxon" linking compound to organism | Compound P703 -> Organism |
| `P248` | Wikidata property "stated in" providing bibliographic reference | Claim P248 -> Reference |
| `P235` | Wikidata property for InChIKey identifier | Compound P235 -> InChIKey |
| `ncbiId` | NCBI Taxonomy identifier for organisms | 3702 (Arabidopsis thaliana) |
| `ottId` | Open Tree of Life identifier for organisms | 1054861 |
| `npl_score` | Natural product likeness score indicating compound similarity to natural products | 0.85 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Natural Product | Secondary metabolite produced by living organisms | LOTUS scope |
| Structure-Organism Pair | Documented association between chemical structure and source organism | Core data model |
| Taxon | Taxonomic unit (species, genus, family) from which compound originates | P703 property |
| InChI | International Chemical Identifier, IUPAC-standard for chemical structures | inchiKey field |
| SPARQL | Query language for RDF databases like Wikidata | API access method |
| Federated Query | SPARQL query spanning multiple endpoints | Wikidata + IDSM |
| Substructure Search | Finding compounds containing a specific molecular fragment | IDSM endpoint |
| Murcko Framework | Core molecular scaffold with side chains removed | LNPN schema |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| LOTUS | Linked Open Total Unified Structure-organism | Database name |
| LNPN | LOTUS Natural Products Network | Web interface |
| SSOT | Single Source of Truth | PostgreSQL canonical |
| CC0 | Creative Commons Zero | License (public domain) |
| SMILES | Simplified Molecular Input Line Entry System | Structure notation |
| InChI | International Chemical Identifier | IUPAC standard |
| RDF | Resource Description Framework | Wikidata format |
| SPARQL | SPARQL Protocol and RDF Query Language | Query language |
| IDSM | Integrated Database for Structure and Metabolism | ELIXIR substructure search |
| OTT | Open Tree of Life Taxonomy | Taxonomy system |
| NCBI | National Center for Biotechnology Information | Taxonomy source |
| DOI | Digital Object Identifier | Publication identifier |
| PMID | PubMed Identifier | Literature reference |

---

## References

1. Rutz A, et al. (2022) "The LOTUS initiative for open knowledge management in natural products research." eLife 11:e70780.

2. Wikidata Project Page: https://www.wikidata.org/wiki/Wikidata:WikiProject_Chemistry/Natural_products

3. LOTUS Manuscript: https://lotus.nprod.net/lotus-manuscript/

4. Zenodo Data Repository: https://zenodo.org/communities/the-lotus-initiative/

5. IDSM/Sachem Documentation: https://idsm.elixir-czech.cz/
