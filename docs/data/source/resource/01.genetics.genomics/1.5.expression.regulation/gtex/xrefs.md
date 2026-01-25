---
id: xrefs-gtex
title: "GTEx Cross-References"
type: xrefs
parent: README.md
last_updated: 2026-01-23
---

# GTEx Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `01.genetics.genomics/1.5.expression.regulation/gtex` |
| Secondary | Symlink | `03.diseases.phenotypes/3.7.mental.health.neurological/gtex` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| Ensembl Gene | ensembl_gene_id | High |
| Ensembl Transcript | ensembl_transcript_id | High |
| HGNC | hgnc_symbol | High |
| dbSNP | variant_id | High |
| Entrez Gene | entrez_id | High |
| UBERON | tissue_ontology | High |
| dbGaP | dbgap_study | High |

## Integration Notes

GTEx (Genotype-Tissue Expression) provides expression data across human tissues with genetic variant associations.

**Primary use (Expression Regulation):** Tissue-specific gene expression, eQTLs (expression quantitative trait loci), sQTLs (splicing QTLs), allele-specific expression.

**Secondary use (Mental Health/Neurological):** Brain region-specific expression, neurological tissue eQTLs, brain transcriptomics for psychiatric genetics.

GTEx's brain tissue data (13 regions) makes it invaluable for interpreting GWAS hits for neurological and psychiatric conditions in tissue-specific context.
