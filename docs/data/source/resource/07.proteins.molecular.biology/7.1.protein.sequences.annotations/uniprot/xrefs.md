---
id: xrefs-uniprot
title: "UniProt Cross-References"
type: xrefs
parent: README.md
last_updated: 2026-01-23
---

# UniProt Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `07.proteins.molecular.biology/7.1.protein.sequences.annotations/uniprot` |
| Secondary | Symlink | `08.literature.knowledge/8.3.identifier.mapping/uniprot` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| NCBI Gene | GeneID | High |
| Ensembl | Ensembl | High |
| RefSeq | RefSeq | High |
| PDB | PDB | High |
| GO | GO | High |
| InterPro | InterPro | High |
| Pfam | Pfam | High |
| KEGG | KEGG | High |
| Reactome | Reactome | High |
| ChEMBL | ChEMBL | High |
| DrugBank | DrugBank | Medium |
| OMIM | MIM | High |
| STRING | STRING | High |

## Integration Notes

UniProt is the universal protein resource, serving as both the authoritative protein sequence/annotation database and the primary hub for cross-database identifier mapping.

**Primary use (Protein Sequences/Annotations):** Protein sequences, functional annotations, post-translational modifications, subcellular localization, tissue expression, disease associations.

**Secondary use (Identifier Mapping):** ID mapping service connecting 200+ databases, accession history, isoform mapping, proteome reference sets.

UniProt's comprehensive cross-references make it essential for any integration pipeline requiring protein-centric identifier resolution across molecular biology databases.
