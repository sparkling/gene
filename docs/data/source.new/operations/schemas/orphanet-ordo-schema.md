---
id: schemas-orphanet-ordo
title: Orphanet/ORDO Rare Disease Schema
category: schemas
subcategory: diseases
tier: 1
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, rare-diseases, orphanet, ontology, gene-disease]
---

**Parent:** [Schema Documentation](./_index.md)

# Orphanet/ORDO Rare Disease Schema

**Document ID:** ORPHANET-ORDO-SCHEMA
**Source:** Orphanet (https://www.orpha.net/) & ORDO Ontology (EBI OLS)
**Last Updated:** January 2026
**ORDO Version:** 4.7

---

## Database Statistics (January 2026)

### Orphanet Database

| Metric | Count |
|--------|-------|
| Rare Diseases | 6,528 |
| Linked Genes | 4,512 |
| Diagnostic Tests | 36,595 |
| Expert Centres | 8,722 |
| Languages | 8 (EN, CZ, NL, FR, DE, IT, ES, PL) |

### ORDO Ontology (via EBI OLS)

| Metric | Count |
|--------|-------|
| Total Terms | 15,788 |
| Properties | 35 |
| Version | 4.7 |

---

## Data Sources

| Resource | URL | Format | License |
|----------|-----|--------|---------|
| Orphanet Portal | https://www.orpha.net/ | Web | CC BY 4.0 |
| Orphadata Science | https://sciences.orphadata.com/ | XML, JSON | CC BY 4.0 |
| ORDO Ontology | https://www.ebi.ac.uk/ols4/ontologies/ordo | OWL | CC BY 4.0 |

---

## Available Datasets

| Product | File | Entries | Description |
|---------|------|---------|-------------|
| Genes | en_product6.xml | 4,128 | Gene-disease associations |
| Phenotypes | en_product4.xml | 4,337 | HPO-linked disease phenotypes |
| Disorders | en_product1.xml | 6,528 | Disease classifications |
| ORDO | ordo.owl | 15,788 | Full ontology |

---

## Gene-Disease Association Schema

```json
{
  "orpha_code": "90033",
  "disease_name": "Warm autoimmune hemolytic anemia",
  "gene_associations": [
    {
      "gene_symbol": "TNFRSF13B",
      "hgnc_id": "HGNC:11924",
      "entrez_id": "23495",
      "association_type": "Disease-causing germline mutation(s) in",
      "association_status": "Assessed"
    }
  ]
}
```

### Association Types

| Type | Description |
|------|-------------|
| Disease-causing germline mutation(s) in | Causal germline variants |
| Disease-causing somatic mutation(s) in | Causal somatic variants |
| Major susceptibility factor in | Increased risk factor |
| Modifying germline mutation in | Disease modifier |

---

## Cross-Reference Systems

| System | Prefix | Example |
|--------|--------|---------|
| OMIM | OMIM: | OMIM:614470 |
| ICD-10 | ICD-10: | ICD-10:D59.1 |
| MeSH | MSH: | MSH:D000744 |
| UMLS | UMLS: | UMLS:C0272118 |
| MONDO | MONDO: | MONDO:0015893 |

---

## License

**CC BY 4.0 (Creative Commons Attribution 4.0 International)**

- Attribution to Orphanet required
- Commercial use permitted

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `ORPHA:XXXXX` | Orphanet rare disease identifier | ORPHA:558 (Marfan syndrome) |
| `orpha_code` | Numeric Orphanet disease code | 90033 |
| `association_type` | Type of gene-disease relationship | Disease-causing germline mutation(s) in |
| `association_status` | Curation status of gene-disease link | Assessed |
| `hgnc_id` | HUGO Gene Nomenclature Committee identifier | HGNC:11924 |
| `entrez_id` | NCBI Gene identifier | 23495 |
| `gene_symbol` | Official gene symbol | TNFRSF13B |
| `diagnostic_test` | Laboratory test for disease diagnosis | Molecular genetic testing |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Rare Disease | Condition affecting fewer than 1 in 2000 individuals in Europe | Orphanet scope |
| ORDO | Orphanet Rare Disease Ontology, OWL representation of Orphanet | Ontology format |
| Disease-causing mutation | Genetic variant that directly causes the disease phenotype | Association type |
| Susceptibility factor | Genetic variant increasing disease risk | Association type |
| Modifying mutation | Variant affecting disease course or severity | Association type |
| Expert Centre | Specialized healthcare center for rare disease care | Healthcare resource |
| Diagnostic Prevalence | How commonly a disease is diagnosed in a population | Epidemiology |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ORDO | Orphanet Rare Disease Ontology | Ontology version |
| ORPHA | Orphanet | Disease identifier prefix |
| OLS | Ontology Lookup Service | EBI ontology browser |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming authority |
| OMIM | Online Mendelian Inheritance in Man | Cross-reference source |
| ICD-10 | International Classification of Diseases, 10th Edition | Cross-reference |
| MeSH | Medical Subject Headings | Cross-reference |
| UMLS | Unified Medical Language System | Cross-reference |
| MONDO | Monarch Disease Ontology | Cross-reference |
| HPO | Human Phenotype Ontology | Phenotype linkage |
| CC BY | Creative Commons Attribution | License type |
| XML | Extensible Markup Language | Data format |
| OWL | Web Ontology Language | Ontology format |

---

## References

1. Orphanet: https://www.orpha.net/
2. Orphadata Science: https://sciences.orphadata.com/
3. ORDO in OLS: https://www.ebi.ac.uk/ols4/ontologies/ordo
