---
id: schemas-mondo
title: MONDO Disease Ontology Schema
category: schemas
subcategory: diseases
tier: 1
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, diseases, ontology, rare-diseases, cross-references]
---

**Parent:** [Schema Documentation](./_index.md)

# MONDO Disease Ontology Schema

**Document ID:** MONDO-SCHEMA
**Status:** Final
**Last Updated:** January 2026
**Version:** 1.0

---

## Overview

MONDO (Monarch Disease Ontology) is a unified disease ontology that integrates multiple disease terminologies into a coherent, logic-based structure. Unlike loose cross-references, MONDO provides **precise 1:1 equivalence axioms** validated by OWL reasoning, enabling safe data propagation across OMIM, Orphanet, EFO, DOID, and NCIt.

### Key Resources

| Resource | URL |
|----------|-----|
| **Homepage** | https://mondo.monarchinitiative.org/ |
| **OBO Foundry** | https://obofoundry.org/ontology/mondo.html |
| **GitHub** | https://github.com/monarch-initiative/mondo |
| **BioPortal** | https://bioportal.bioontology.org/ontologies/MONDO |

---

## Statistics (v2026-01-06)

| Metric | Count |
|--------|-------|
| **Total diseases** | 25,938 |
| **Human diseases** | 22,977 |
| **Cancer classes** | 4,728 |
| **Mendelian conditions** | 11,639 |
| **Rare diseases** | 15,901 |
| **Database cross-references** | 129,914 |

---

## Identifier Format

### CURIE Format (Compact)
```
MONDO:0005015
```

### Pattern
```
MONDO:[0-9]{7}
```

**Example:**
- MONDO:0005015 = diabetes mellitus
- MONDO:0007947 = Marfan syndrome

---

## Cross-Reference Types

### Equivalence Axioms

| Source | Format | Example |
|--------|--------|---------|
| OMIM | OMIM:[0-9]{6} | OMIM:154700 |
| Orphanet | Orphanet:[0-9]+ | Orphanet:558 |
| DOID | DOID:[0-9]+ | DOID:14323 |
| NCIt | NCIT:C[0-9]+ | NCIT:C34807 |
| MeSH | MESH:D[0-9]{6} | MESH:D008382 |

### Mapping Predicates

| Predicate | Meaning |
|-----------|---------|
| `exactMatch` | 1:1 equivalence |
| `closeMatch` | Very similar |
| `narrowMatch` | More specific |
| `broadMatch` | More general |

---

## Download

### Access Methods

**Source:** https://github.com/monarch-initiative/mondo/releases

| Method | URL | Format |
|--------|-----|--------|
| **GitHub Releases** | https://github.com/monarch-initiative/mondo/releases | All formats |
| **Stable PURL** | http://purl.obolibrary.org/obo/mondo.owl | OWL, OBO, JSON |
| **OBO Foundry** | https://obofoundry.org/ontology/mondo.html | Multiple formats |
| **Zenodo** | https://zenodo.org/ | Archived versions |

### Download URLs

```bash
# Stable PURL (redirects to latest)
http://purl.obolibrary.org/obo/mondo.owl
http://purl.obolibrary.org/obo/mondo.obo
http://purl.obolibrary.org/obo/mondo.json
```

---

## Data Format

| Format | Description | Use Case |
|--------|-------------|----------|
| Primary | OWL (Web Ontology Language) | Full semantics, reasoning |
| Alternative | OBO (Open Biological Ontology) | Lightweight, readable |
| Alternative | JSON-LD | Linked data, web integration |
| Equivalence Axioms | MONDO:X exactMatch OMIM:Y | Cross-database mapping |
| Encoding | UTF-8 | All formats |

---

## Integration Notes

### Recommended Usage

1. **Use MONDO as pivot** - Map all disease IDs through MONDO
2. **Trust equivalence axioms** - Safe for data propagation
3. **Check subset membership** - Use `subset: mondo_rare` for filtering

### License

- **CC BY 4.0** - Free for commercial use with attribution
- Attribution: "Mondo Disease Ontology, Monarch Initiative"

---

## Data Set Size

| Component | Count | Size |
|-----------|-------|------|
| **Disease Terms** | 27,000+ | ~240 MB (OWL) |
| **Equivalence Axioms** | 20,000+ | Included in OWL |
| **Relationships** | 50,000+ | Included in OWL |
| **Subsets** | 10+ (mondo_rare, mondo_core, etc.) | Included |
| **mondo.owl** | All | 239.3 MB |
| **mondo.obo** | All | 50.7 MB |
| **mondo.json** | All | 102.7 MB |

---

## Sample Data

### Example Record
```json
{
  "mondo_id": "MONDO:0001816",
  "label": "Sickle cell anemia",
  "definition": "Autosomal recessive hemoglobin disorder",
  "parent_id": "MONDO:0001816",
  "xref": ["OMIM:603903"]
}
```

### Sample Query Result
| mondo_id | label | definition | parent_id | xref |
|---------|--------|-----------|-----------|------|
| MONDO:0001816 | Sickle cell anemia | Autosomal recessive hemoglobin disorder | MONDO:0007522 | OMIM:603903 |
| MONDO:0001816 | Cystic fibrosis | Autosomal recessive lung disease | MONDO:0008168 | OMIM:219700 |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "MONDO:0000001" |
| `name` | string | Entity name | "disease" |
| `type` | string | Record type | "disease" / "disorder" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `is_a` | MONDO Term | N:1 |
| `similar_to` | MONDO Term | N:M |
| `cross_referenced` | Database ID | N:M |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `MONDO:XXXXXXX` | MONDO disease identifier in 7-digit format | MONDO:0005015 (diabetes mellitus) |
| `exactMatch` | Equivalence axiom indicating 1:1 correspondence between disease IDs | MONDO:0007947 exactMatch OMIM:154700 |
| `closeMatch` | Mapping indicating high similarity but not exact equivalence | Approximate match |
| `narrowMatch` | Mapping to a more specific disease concept | Subtype mapping |
| `broadMatch` | Mapping to a more general disease concept | Parent mapping |
| `xref` | Cross-reference to external disease identifier | OMIM:154700, Orphanet:558 |
| `subset` | Filtering category for disease classification | mondo_rare |
| `equivalence axiom` | OWL-reasoned 1:1 mapping enabling safe data propagation | Validated equivalence |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Disease Ontology | Hierarchical classification of human diseases with logical definitions | MONDO structure |
| Mendelian Condition | Single-gene disorder following Mendelian inheritance patterns | OMIM linkage |
| Rare Disease | Disease affecting fewer than 1 in 2000 individuals | mondo_rare subset |
| Cancer Class | Neoplastic disease classification | NCIt alignment |
| Equivalence Propagation | Using 1:1 mappings to transfer annotations safely between databases | Integration strategy |
| Pivot Ontology | Central resource for mapping between multiple disease databases | MONDO use case |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| MONDO | Monarch Disease Ontology | Unified disease ontology |
| OMIM | Online Mendelian Inheritance in Man | Mendelian disease database |
| DOID | Disease Ontology Identifier | Alternative disease ontology |
| NCIt | NCI Thesaurus | Cancer terminology |
| EFO | Experimental Factor Ontology | Trait/disease ontology |
| MeSH | Medical Subject Headings | NLM disease vocabulary |
| OWL | Web Ontology Language | Ontology format |
| OBO | Open Biological Ontologies | Ontology format |
| CURIE | Compact URI Expression | Identifier format |
| CC BY | Creative Commons Attribution | License type |

---

## References

- Shefchek et al. (2020) "The Monarch Initiative in 2019" Nucleic Acids Research
- MONDO Documentation: https://mondo.readthedocs.io/
