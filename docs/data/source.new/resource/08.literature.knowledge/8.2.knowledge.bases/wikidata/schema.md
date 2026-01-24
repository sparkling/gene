---
id: schemas-wikidata-schema
title: "Wikidata Biomedical Schema"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, wikidata, sparql, knowledge-graph, genes, diseases, cc0]
---

**Parent:** [Schema Documentation](./_index.md)

# Wikidata Biomedical Schema

**Document ID:** WIKIDATA-SCHEMA
**Status:** Final
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Wikidata is a free and open knowledge base that serves as the central structured data repository for Wikimedia projects. For biomedical applications, it contains comprehensive coverage of human genes (~59,721), proteins (~27,306), diseases (~200,000), and medications (~45,000) under CC0 (public domain) license.

### SPARQL Endpoint

| Attribute | Value |
|-----------|-------|
| **Query Interface** | https://query.wikidata.org/ |
| **SPARQL Endpoint** | https://query.wikidata.org/sparql |
| **Max Results** | 1,000,000 (with streaming) |
| **Timeout** | 60 seconds (default) |
| **Rate Limit** | Not strictly enforced, but be considerate |

---

## Gene Properties

### Primary Gene Identifiers

| Property | Name | Format | Example | Usage Count (Human) |
|----------|------|--------|---------|---------------------|
| P351 | Entrez Gene ID | Numeric string | "7157" | ~59,721 |
| P352 | UniProt Protein ID | [A-Z][0-9][A-Z0-9]{3}[0-9] | "P04637" | ~27,306 |
| P353 | HGNC Gene Symbol | Alphanumeric | "TP53" | ~20,000+ |
| P354 | HGNC ID | Numeric with prefix | "HGNC:11998" | ~20,000+ |
| P594 | Ensembl Gene ID | ENSG[0-9]{11} | "ENSG00000141510" | ~60,000+ |
| P639 | RefSeq Protein ID | [NX]P_[0-9]+\.[0-9]+ | "NP_000537.3" | Variable |
| P637 | RefSeq RNA ID | [NX][MR]_[0-9]+\.[0-9]+ | "NM_000546.6" | Variable |

---

## Disease Properties

### Primary Disease Identifiers

| Property | Name | Format | Example |
|----------|------|--------|---------|
| P699 | Disease Ontology ID | DOID:[0-9]+ | "DOID:162" |
| P492 | OMIM ID | 6 digits | "191170" |
| P1395 | Orphanet ID | [0-9]+ | "319182" |
| P493 | ICD-9-CM | [0-9]{3}(.[0-9]+)? | "250.00" |
| P494 | ICD-10 | [A-Z][0-9]{2}(.[0-9]+)? | "E11.9" |
| P5270 | MONDO ID | MONDO:[0-9]+ | "MONDO:0005148" |
| P486 | MeSH ID | [A-Z][0-9]{6,9} | "D003920" |
| P2892 | UMLS CUI | C[0-9]+ | "C0011849" |

---

## Compound Properties

### Chemical Identifiers

| Property | Name | Format | Example |
|----------|------|--------|---------|
| P231 | CAS Registry Number | [0-9]+-[0-9]+-[0-9] | "50-00-0" |
| P662 | PubChem CID | [0-9]+ | "2244" |
| P592 | ChEMBL ID | CHEMBL[0-9]+ | "CHEMBL25" |
| P652 | ChEBI ID | [0-9]+ | "15377" |
| P715 | DrugBank ID | DB[0-9]+ | "DB00945" |
| P234 | InChI | InChI=... | Full InChI string |
| P235 | InChIKey | [A-Z]{14}-[A-Z]{10}-[A-Z] | "BSYNRYMUTXBXSQ-UHFFFAOYSA-N" |
| P233 | SMILES | Chemical notation | "CC(=O)Oc1ccccc1C(=O)O" |

---

## Download

### Data Access Methods

| Format | Size (approx) | URL | Use Case |
|--------|---------------|-----|----------|
| JSON Dumps | ~60 GB (compressed) | https://dumps.wikimedia.org/wikidatawiki/ | Full knowledge base export |
| RDF Dumps | ~80 GB (compressed) | https://dumps.wikimedia.org/wikidatawiki/entities/ | Semantic web/RDF applications |
| SPARQL Endpoint | - | https://query.wikidata.org/sparql | Real-time queries |
| REST API | - | https://www.wikidata.org/w/api.php | Programmatic access |
| Quickstatements | - | https://www.wikidata.org/wiki/Wikidata:QuickStatements | Batch data retrieval |

### Update Schedule

- JSON/RDF dumps: Daily (around 01:00 UTC)
- Real-time data: Available instantly via SPARQL
- API updates: Continuous (changes reflected immediately)

---

## Data Format

| Aspect | Details |
|--------|---------|
| **Primary Formats** | JSON-LD (API), Turtle/N-Triples (RDF), SPARQL responses |
| **Alternative Formats** | RDF/XML, CSV (custom exports) |
| **Encoding** | UTF-8 |
| **API Response** | JSON with RDF structure |
| **Query Language** | SPARQL 1.1 |
| **Compression** | bzip2 for dump files |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| **Total Items** | 100,000,000+ (all entities) |
| **Human Genes** | ~59,721 entries |
| **Human Proteins** | ~27,306 entries |
| **Human Diseases** | ~200,000+ entries |
| **Medications/Drugs** | ~45,000+ entries |
| **Gene-Disease Links** | 1,000,000+ associations |
| **Protein Properties** | 286+ cross-database mappings via properties |
| **Total Properties** | 10,000+ predicates |
| **Labels in Languages** | 300+ languages |
| **Total Statements** | 1,000,000,000+ RDF triples |
| **JSON Dump Size** | ~60 GB (compressed bzip2) |
| **RDF Dump Size** | ~80 GB (compressed bzip2) |
| **Growth Rate** | ~1 million new statements daily |
| **Update Frequency** | Continuous (real-time); dumps daily |

---

## Sample SPARQL Queries

### Get Human Genes with Multiple Identifiers

```sparql
SELECT ?gene ?geneLabel ?entrezId ?ensembl ?hgncSymbol ?uniprot WHERE {
  ?gene wdt:P31 wd:Q7187 .                    # instance of: gene
  ?gene wdt:P703 wd:Q15978631 .               # found in taxon: Homo sapiens
  ?gene wdt:P351 ?entrezId .                  # Entrez Gene ID
  OPTIONAL { ?gene wdt:P594 ?ensembl . }      # Ensembl gene ID
  OPTIONAL { ?gene wdt:P353 ?hgncSymbol . }   # HGNC symbol
  OPTIONAL { ?gene wdt:P352 ?uniprot . }      # UniProt ID
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

### Get Diseases with Cross-References

```sparql
SELECT ?disease ?diseaseLabel ?mondo ?omim ?orphanet ?icd10 WHERE {
  ?disease wdt:P31 wd:Q12136 .               # instance of: disease
  OPTIONAL { ?disease wdt:P5270 ?mondo . }    # MONDO ID
  OPTIONAL { ?disease wdt:P492 ?omim . }      # OMIM ID
  OPTIONAL { ?disease wdt:P1395 ?orphanet . } # Orphanet ID
  OPTIONAL { ?disease wdt:P494 ?icd10 . }     # ICD-10
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

### Get Gene-Disease Associations

```sparql
SELECT ?gene ?geneLabel ?disease ?diseaseLabel WHERE {
  ?gene wdt:P31 wd:Q7187 .                   # instance of: gene
  ?gene wdt:P703 wd:Q15978631 .              # found in taxon: human
  ?gene wdt:P2293 ?disease .                 # genetic association
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

---

## Coverage Statistics

| Entity Type | Approximate Count | Key Properties |
|-------------|-------------------|----------------|
| Human Genes | ~59,721 | P351 (Entrez), P353 (HGNC) |
| Proteins | ~27,306 | P352 (UniProt) |
| Diseases | ~200,000 | P699 (DOID), P5270 (MONDO) |
| Medications | ~45,000 | P715 (DrugBank), P592 (ChEMBL) |
| Chemical Compounds | ~1,000,000+ | P662 (PubChem), P235 (InChIKey) |

---

## License

**License:** CC0 1.0 Universal (Public Domain)
**Attribution:** Not required
**Commercial Use:** Unrestricted

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "Q123456" |
| `name` | string | Entity name | "Item Label" |
| `type` | string | Record type | "entity" / "property" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `is_instance_of` | Class | N:M |
| `connected_to` | Entity | N:M |
| `described_by` | Language | N:M |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `P351` | Property ID for Entrez Gene ID in Wikidata | `7157` (TP53) |
| `P352` | Property ID for UniProt Protein ID in Wikidata | `P04637` |
| `P699` | Property ID for Disease Ontology ID in Wikidata | `DOID:162` |
| `P662` | Property ID for PubChem CID in Wikidata | `2244` |
| `wdt:` | Wikidata truthy statement prefix in SPARQL | `wdt:P351` |
| `wd:` | Wikidata entity prefix in SPARQL | `wd:Q7187` (gene) |
| `Q7187` | Wikidata item ID for "gene" class | Instance type |
| `Q12136` | Wikidata item ID for "disease" class | Instance type |
| `Q15978631` | Wikidata item ID for "Homo sapiens" | Species filter |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Property | Wikidata predicate describing relationships | P351, P699 |
| Item | Wikidata entity representing a concept | Q7187, Q12136 |
| Statement | Fact about an item using a property and value | Gene has Entrez ID |
| Qualifier | Additional context for a statement | References, dates |
| SPARQL | Query language for RDF knowledge graphs | Query endpoint |
| Truthy Statements | Best-ranked statements for a property | wdt: prefix |
| Entity | Any Wikidata object (item, property, or lexeme) | Q-IDs, P-IDs |
| Knowledge Graph | Structured database of interlinked entities | Wikidata architecture |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| SPARQL | SPARQL Protocol and RDF Query Language | Knowledge graph query language |
| RDF | Resource Description Framework | Linked data standard |
| HGNC | HUGO Gene Nomenclature Committee | Human gene naming authority |
| OMIM | Online Mendelian Inheritance in Man | Genetic disease database |
| DOID | Disease Ontology ID | Disease classification system |
| MONDO | Mondo Disease Ontology | Unified disease ontology |
| MeSH | Medical Subject Headings | NLM controlled vocabulary |
| UMLS | Unified Medical Language System | NLM metathesaurus |
| CUI | Concept Unique Identifier | UMLS concept ID format |
| CAS | Chemical Abstracts Service | Chemical registry system |
| InChI | International Chemical Identifier | IUPAC chemical notation |
| InChIKey | InChI Hash Key | 27-character structure hash |
| CC0 | Creative Commons Zero | Public domain dedication |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
