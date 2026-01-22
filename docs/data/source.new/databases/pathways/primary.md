---
id: pathways-primary
title: Primary Pathway Databases
world: null
category: pathways
subcategory: null
tier: 1
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
databases:
  - reactome
  - wikipathways
  - kegg-pathway
  - metacyc
  - biocyc
  - smpdb
  - pharmgkb-pathways
  - panther
  - pathway-commons
  - ndex
  - pathbank
tags: [pathways, biological-processes, mechanisms, networks]
---

# Pathways Primary Databases

**Document ID:** 43-41-PATHWAYS-PRIMARY
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

Primary pathway databases provide curated biological pathway data essential for understanding gene-compound-disease relationships. Reactome (CC BY 4.0, 2,712 human pathways) and WikiPathways (CC0, 3,100+ pathways) are recommended as open-access priorities, with KEGG, MetaCyc, SMPDB, PharmGKB Pathways, PANTHER, Pathway Commons, and NDEx providing specialized coverage for metabolism, pharmacogenomics, and network analysis.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary pathway source | Reactome | CC BY 4.0 license allows commercial use; highest curation quality; Neo4j + REST API | Jan 2026 |
| Secondary pathway source | WikiPathways | CC0 public domain; community contributions capture emerging biology | Jan 2026 |
| Drug metabolism pathways | SMPDB | Best HMDB/DrugBank integration; comprehensive PK pathways | Jan 2026 |
| Pharmacogenomics pathways | PharmGKB Pathways | Clinical variant-drug relationships; CPIC guideline integration | Jan 2026 |
| KEGG approach | Reference only | Bulk FTP requires subscription; use Reactome/WikiPathways for open alternatives | Jan 2026 |
| Integration layer | Pathway Commons | Pre-integrated data from multiple sources; BioPAX standard format | Jan 2026 |
| Network repository | NDEx | Access specialized networks; NCI-PID archive preserved | Jan 2026 |
| Primary compound identifier | ChEBI | Used by Reactome, WikiPathways, Pathway Commons; open ontology | Jan 2026 |

---

[Full content from primary.md original file for pathways]

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `pathway` | A series of molecular interactions and reactions leading to a biological outcome | Glycolysis pathway converting glucose to pyruvate |
| `metabolic pathway` | A linked series of enzyme-catalyzed biochemical reactions | Folate metabolism pathway |
| `signaling pathway` | A cascade of protein interactions transmitting cellular signals | MAPK signaling pathway |
| `BioPAX` | Biological Pathway Exchange - standard format for pathway data exchange | BioPAX Level 3 OWL files |
| `SBML` | Systems Biology Markup Language - XML format for computational models | SBML Level 3 Version 2 |
| `gene set` | A defined group of genes sharing a biological function or pathway membership | KEGG pathway gene sets |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `Reactome` | Curated pathway database with 2,712 human pathways and CC BY 4.0 license | Primary pathway source |
| `WikiPathways` | Community-curated pathway database with CC0 license | Secondary source |
| `KEGG` | Kyoto Encyclopedia of Genes and Genomes - pathway and genome database | Metabolism pathways |
| `MetaCyc` | Curated database of metabolic pathways from all domains of life | Metabolic coverage |
| `SMPDB` | Small Molecule Pathway Database integrated with HMDB and DrugBank | Drug metabolism |
| `ChEBI` | Chemical Entities of Biological Interest - chemical ontology used across pathways | Compound identifiers |
| `NDEx` | Network Data Exchange - repository for biological networks | Network sharing |
| `Pathway Commons` | Integrated pathway data from multiple sources in BioPAX format | Data aggregation |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway + genome database |
| BioPAX | Biological Pathway Exchange | Standard pathway format |
| SBML | Systems Biology Markup Language | Computational model format |
| SMPDB | Small Molecule Pathway Database | Drug/metabolite pathways |
| GO | Gene Ontology | Functional annotation standard |
| PPI | Protein-Protein Interaction | Network interaction type |
| TF | Transcription Factor | Gene regulation proteins |
| GSEA | Gene Set Enrichment Analysis | Pathway analysis method |
| OWL | Web Ontology Language | BioPAX file format |
| RDF | Resource Description Framework | Linked data format |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [../_index.md](../_index.md) | Parent index |
| [disease.md](./disease.md) | Disease pathway databases |
| [../compounds/pharmaceuticals.md](../compounds/pharmaceuticals.md) | Drug database integration |

---

## References

1. Jassal B, et al. (2020). The reactome pathway knowledgebase. Nucleic Acids Res. 48(D1):D498-D503. https://doi.org/10.1093/nar/gkz1031

2. Martens M, et al. (2021). WikiPathways: connecting communities. Nucleic Acids Res. 49(D1):D613-D621. https://doi.org/10.1093/nar/gkaa1024

3. Kanehisa M, et al. (2023). KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Res. 51(D1):D587-D592. https://doi.org/10.1093/nar/gkac963

4. Caspi R, et al. (2020). The MetaCyc database of metabolic pathways and enzymes. Nucleic Acids Res. 48(D1):D445-D453. https://doi.org/10.1093/nar/gkz862

5. Jewison T, et al. (2014). SMPDB 2.0: Big Improvements to the Small Molecule Pathway Database. Nucleic Acids Res. 42:D478-D484. https://doi.org/10.1093/nar/gkt1067

6. Whirl-Carrillo M, et al. (2021). An Evidence-Based Framework for Evaluating Pharmacogenomics Knowledge. Clin Pharmacol Ther. 110(3):563-572. https://doi.org/10.1002/cpt.2350

7. Mi H, et al. (2021). PANTHER version 16. Nucleic Acids Res. 49(D1):D419-D426. https://doi.org/10.1093/nar/gkaa1106

8. Rodchenkov I, et al. (2020). Pathway Commons 2019 Update. Nucleic Acids Res. 48(D1):D489-D497. https://doi.org/10.1093/nar/gkz946

9. Pratt D, et al. (2015). NDEx, the Network Data Exchange. Cell Syst. 1(4):302-305. https://doi.org/10.1016/j.cels.2015.10.001

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial catalog: 10 pathway databases migrated from research.old/pathways-databases.md |
