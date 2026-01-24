---
id: xrefs-dgidb
title: "DGIdb Cross-References"
type: xrefs
parent: _index.md
last_updated: 2026-01-23
---

# DGIdb Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `02.compounds.molecules/2.7.compound.target.interactions/dgidb` |
| Secondary | Symlink | `04.pathways.networks/4.4.drug.target.interactions/dgidb` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| Entrez Gene | entrez_id | High |
| UniProt | uniprot_id | High |
| Ensembl | ensembl_id | High |
| ChEMBL | chembl_id | High |
| DrugBank | drugbank_id | High |
| PubChem | pubchem_cid | High |
| PharmGKB | pharmgkb_id | Medium |
| TTD | ttd_id | Medium |
| PubMed | pmid | High |

## Integration Notes

DGIdb (Drug Gene Interaction Database) aggregates drug-gene interaction data from multiple sources for druggable genome analysis.

**Primary use (Compound-Target Interactions):** Drug-gene interactions, interaction types, drug categories, source databases.

**Secondary use (Drug-Target Interactions):** Druggable genome annotations, gene categories (kinases, ion channels, etc.), clinical actionability.

DGIdb integrates interaction data from 30+ sources, enabling identification of druggable targets and known drugs for genes of interest.
