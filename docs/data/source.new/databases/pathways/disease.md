---
id: pathways-disease
title: Disease and Phenotype Databases
category: shared
subcategory: disease
tier: 1
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
databases:
  - disgenet
  - omim
  - hpo
  - mondo
  - orphanet
  - ordo
  - clinvar
  - clingen
  - disease-ontology
  - malacards
  - monarch
  - decipher
  - gwas-catalog
tags: [disease, phenotypes, rare-diseases, gene-disease, variants]
---

# Disease and Phenotype Databases

**Document ID:** 43-43-PATHWAYS-DISEASE
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

Disease and phenotype databases link genes, variants, and biological pathways to human diseases. Primary sources include DisGeNET (1.1M+ gene-disease associations), OMIM (25K+ Mendelian disorders), and Human Phenotype Ontology (HPO, 13K+ phenotype terms). These databases enable genotype-phenotype correlation, disease mechanism understanding, and clinical variant interpretation.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary gene-disease source | DisGeNET | Largest curated collection, integrates multiple sources | Jan 2026 |
| Mendelian disease reference | OMIM | Gold standard for genetic disorders, comprehensive | Jan 2026 |
| Phenotype ontology | HPO | Standard vocabulary, cross-species integration via Monarch | Jan 2026 |
| Disease ontology | MONDO | Unified ontology with precise equivalence axioms | Jan 2026 |
| Rare disease integration | Orphanet/ORDO | ELIXIR Core Resource, comprehensive rare disease coverage | Jan 2026 |
| Clinical variant interpretation | ClinGen + ClinVar | NIH-funded authoritative sources | Jan 2026 |

---

[Full content from disease.md original file]

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **DisGeNET** | Download | `https://www.disgenet.org/downloads` |
| **OMIM** | API | `https://api.omim.org/` (registration required) |
| **HPO** | Download | `https://hpo.jax.org/app/data/ontology` |
| **MONDO** | Download | `https://mondo.monarchinitiative.org/` |
| **ClinVar** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/` |
| **Orphanet** | Download | `https://www.orphadata.com/` |

**Access Requirements:** Most are freely accessible; OMIM requires registration; DisGeNET has academic and commercial licenses.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | TSV, OBO, OWL |
| Alternative | JSON, XML, RDF |
| Ontology | OBO, OWL |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `disease_id` | string | Disease identifier | "MONDO:0005148" |
| `gene_symbol` | string | Gene symbol | "MTHFR" |
| `association_type` | string | Evidence type | "causal" |
| `score` | float | Association confidence | 0.85 |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Disease | N:M |
| `has_phenotype` | HPO Term | N:M |
| `variant_in` | ClinVar Entry | N:M |

## Sample Data

### Example Gene-Disease Association
```json
{
  "gene_symbol": "MTHFR",
  "disease": "Neural tube defects",
  "mondo_id": "MONDO:0005148",
  "score": 0.85,
  "evidence": ["GWAS", "literature"],
  "variants": ["rs1801133", "rs1801131"]
}
```

### Sample Query Result
| gene | disease | score | evidence |
|------|---------|-------|----------|
| MTHFR | Neural tube defects | 0.85 | GWAS |
| BRCA1 | Breast cancer | 0.99 | Clinical |

## License

| Source | License | Commercial Use |
|--------|---------|----------------|
| DisGeNET | CC BY-NC-SA (free) | Commercial license available |
| OMIM | Terms of use | Requires license |
| HPO | Custom (free) | Yes with attribution |
| MONDO | CC BY 4.0 | Yes |
| ClinVar | Public domain | Yes |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| DisGeNET | 1.1M+ gene-disease associations |
| OMIM | 25K+ Mendelian disorders |
| HPO phenotype terms | 13K+ terms |
| MONDO disease concepts | Cross-referenced ontology |
| ClinVar variants | 3M+ variant submissions |
| Total storage estimate | ~10-15 GB (combined sources) |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `gene-disease association` | A relationship between a gene and a disease phenotype, with supporting evidence | MTHFR associated with neural tube defects |
| `phenotype` | An observable characteristic or trait of an organism | Elevated homocysteine levels |
| `ontology` | A structured vocabulary defining concepts and relationships in a domain | Human Phenotype Ontology |
| `Mendelian disorder` | A genetic disorder caused by mutation in a single gene | Cystic fibrosis (CFTR gene) |
| `complex disease` | A condition influenced by multiple genes and environmental factors | Type 2 diabetes |
| `penetrance` | The proportion of individuals with a genotype who express the phenotype | 80% penetrance for BRCA1 mutations |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `DisGeNET` | Database of gene-disease associations with 1.1M+ curated relationships | Gene-disease links |
| `OMIM` | Online Mendelian Inheritance in Man - catalog of human genes and genetic disorders | Mendelian diseases |
| `HPO` | Human Phenotype Ontology - standardized vocabulary of 13K+ phenotypic abnormalities | Phenotype annotation |
| `MONDO` | Mondo Disease Ontology - unified disease classification with cross-references | Disease ontology |
| `Orphanet` | Portal for rare diseases with ORDO ontology | Rare disease resource |
| `ClinVar` | NCBI database of clinically relevant genetic variants and phenotypes | Variant interpretation |
| `ClinGen` | Clinical Genome Resource providing gene-disease validity curation | Clinical validity |
| `GWAS Catalog` | Curated collection of genome-wide association study results | Common disease genetics |
| `Monarch Initiative` | Platform integrating genotype-phenotype data across species | Cross-species integration |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| OMIM | Online Mendelian Inheritance in Man | Gold standard for genetic disorders |
| HPO | Human Phenotype Ontology | Standard phenotype vocabulary |
| MONDO | Mondo Disease Ontology | Unified disease classification |
| ORDO | Orphanet Rare Disease Ontology | Rare disease classification |
| GDA | Gene-Disease Association | DisGeNET relationship type |
| VUS | Variant of Uncertain Significance | ClinVar classification |
| ACMG | American College of Medical Genetics | Variant classification guidelines |
| P/LP | Pathogenic/Likely Pathogenic | ClinVar clinical significance |
| GWAS | Genome-Wide Association Study | Common variant discovery method |
| EFO | Experimental Factor Ontology | GWAS Catalog trait terms |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [primary.md](./primary.md) | Pathway-disease connections |
| [../../genetics/primary.md](../../genetics/primary.md) | Gene-disease variant data |
| [../../diseases/rare.md](../../diseases/rare.md) | Detailed rare disease sources |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial disease database catalog |
