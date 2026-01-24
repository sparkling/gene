---
id: xrefs-clinvar
title: "ClinVar Cross-References"
type: xrefs
parent: _index.md
last_updated: 2026-01-23
---

# ClinVar Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `01.genetics.genomics/1.1.variant.repositories/clinvar` |
| Secondary | Symlink | `03.diseases.phenotypes/3.3.disease.gene.associations/clinvar` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| dbSNP | rsid | High |
| NCBI Gene | gene_id | High |
| OMIM | omim | High |
| MedGen | medgen_cui | High |
| Orphanet | orpha | Medium |
| UniProt | uniprot_variant | Medium |
| RefSeq | refseq | High |
| Ensembl | ensembl | High |
| PubMed | pubmed | High |
| HPO | hpo | Medium |

## Integration Notes

ClinVar aggregates variant-disease interpretations from clinical laboratories, bridging genetic variants with clinical significance.

**Primary use (Variant Repositories):** Clinical variant interpretations, pathogenicity classifications, review status, submitter information.

**Secondary use (Disease-Gene Associations):** Variant-disease relationships, clinical conditions, molecular consequences.

ClinVar is essential for clinical variant interpretation in genetic testing, connecting genomic findings to clinical phenotypes with standardized terminology.
