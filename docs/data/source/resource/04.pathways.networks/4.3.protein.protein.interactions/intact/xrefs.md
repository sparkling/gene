---
id: xrefs-intact
title: "IntAct Cross-References"
type: xrefs
parent: README.md
last_updated: 2026-01-23
---

# IntAct Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `04.pathways.networks/4.3.protein.protein.interactions/intact` |
| Secondary | Symlink | `07.proteins.molecular.biology/7.3.molecular.interactions/intact` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| UniProt | uniprotkb | High |
| ChEBI | chebi | High (small molecules) |
| Ensembl | ensembl | Medium |
| PubMed | pubmed | High |
| DOI | doi | High |
| PSI-MI | psi-mi | High (ontology) |
| IMEx | imex | High |

## Integration Notes

IntAct is the primary EMBL-EBI database for molecular interaction data, providing curated binary interactions with full experimental details.

**Primary use (Protein-Protein Interactions):** Experimentally validated interactions, interaction detection methods, stoichiometry, binding domains.

**Secondary use (Molecular Interactions):** Protein-small molecule interactions, protein-nucleic acid interactions, interaction network construction.

IntAct is part of the IMEx consortium, ensuring data exchange with other major interaction databases (MINT, DIP, BioGRID) through standardized formats (PSI-MI XML, MITAB).
