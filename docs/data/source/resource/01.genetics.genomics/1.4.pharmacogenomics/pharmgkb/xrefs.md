---
id: xrefs-pharmgkb
title: "PharmGKB Cross-References"
type: xrefs
parent: README.md
last_updated: 2026-01-23
---

# PharmGKB Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `01.genetics.genomics/1.4.pharmacogenomics/pharmgkb` |
| Secondary | Symlink | `02.compounds.molecules/2.5.drug.metabolism.pharmacokinetics/pharmgkb` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| NCBI Gene | gene_id | High |
| UniProt | uniprot_id | High |
| Ensembl | ensembl_id | High |
| dbSNP | rsid | High |
| DrugBank | drugbank_id | High |
| PubChem | pubchem_cid | High |
| ChEMBL | chembl_id | High |
| RxNorm | rxnorm_id | High |
| ATC | atc_code | High |
| HGNC | hgnc_id | High |

## Integration Notes

PharmGKB curates pharmacogenomic relationships, bridging genetic variation with drug response and metabolism.

**Primary use (Pharmacogenomics):** Gene-drug relationships, variant-drug associations, clinical annotations, dosing guidelines, drug labels with PGx information.

**Secondary use (Drug Metabolism/PK):** Drug metabolizing enzyme variants, transporter polymorphisms, drug-drug interactions via PGx mechanisms.

PharmGKB integrates with CPIC for clinical guidelines and provides the knowledge base connecting genetic variants to drug efficacy, toxicity, and dosing recommendations.
