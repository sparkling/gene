---
id: xrefs-disgenet
title: "DisGeNET Cross-References"
type: xrefs
parent: README.md
last_updated: 2026-01-23
---

# DisGeNET Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `03.diseases.phenotypes/3.3.disease.gene.associations/disgenet` |
| Secondary | Symlink | `01.genetics.genomics/1.5.expression.regulation/disgenet` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| NCBI Gene | geneid | High |
| UniProt | uniprot | High |
| UMLS | umls_cui | High |
| DO | doid | High |
| HPO | hpo | Medium |
| EFO | efo | Medium |
| MONDO | mondo | Medium |
| MeSH | mesh | High |
| OMIM | omim | High |
| ICD | icd9/icd10 | Medium |

## Integration Notes

DisGeNET integrates gene-disease associations from multiple sources, bridging disease genetics with expression and regulatory data.

**Primary use (Disease-Gene Associations):** Curated and text-mined gene-disease relationships, association scores, evidence types, variant-disease associations.

**Secondary use (Expression Regulation):** Disease-associated genes for expression analysis, regulatory context of disease genes, pathway enrichment for disease gene sets.

DisGeNET aggregates data from GWAS Catalog, ClinVar, CTD, UniProt, and text mining, making it valuable for both disease genetics and functional genomics applications.
