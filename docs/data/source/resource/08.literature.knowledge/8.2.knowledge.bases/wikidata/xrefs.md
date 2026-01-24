---
id: xrefs-wikidata
title: "Wikidata Cross-References"
type: xrefs
parent: _index.md
last_updated: 2026-01-23
---

# Wikidata Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary | Canonical | `08.literature.knowledge/8.2.knowledge.bases/wikidata` |
| Secondary | Symlink | `05.traditional.medicine/5.4.multi.system.integration/wikidata` |

## External ID Mappings

| External DB | ID Field | Coverage |
|-------------|----------|----------|
| UniProt | P352 | High (proteins) |
| PubChem | P662 | High (compounds) |
| ChEMBL | P592 | High (drugs) |
| OMIM | P492 | High (diseases) |
| Ensembl | P594 | High (genes) |
| PubMed | P698 | High (publications) |
| MeSH | P486 | Medium |
| HGNC | P354 | High (gene symbols) |
| DrugBank | P715 | Medium |
| Gene Ontology | P686 | High |

## Integration Notes

Wikidata appears in Traditional Medicine integration because it serves as a universal identifier bridge across biomedical domains. The knowledge graph connects traditional medicine concepts (herbs, formulas, symptoms) to modern biomedical entities through structured statements.

**Primary use (Knowledge Bases):** Entity resolution, cross-database linking, knowledge graph construction.

**Secondary use (Traditional Medicine Integration):** Linking traditional medicine terms to standardized identifiers, multilingual terminology mapping, connecting ethnobotanical knowledge to molecular data.

Use SPARQL queries to traverse relationships between traditional medicine and biomedical entities.
