# MVP Data Source Selection Plan

Comprehensive analysis of data sources for the biomedical knowledge graph MVP, scored against 3-tier criteria.

## Overview

- **Total Sources Selected**: 55
- **Categories Covered**: 9 (all first-level subcategories)
- **Analysis Date**: 2026-01-25
- **Methodology**: Parallel agent swarm analysis (9 agents)

## Selection Criteria

### Tier 1: Must-Have (9 points max)
| Criterion | Points | Description |
|-----------|--------|-------------|
| Identifier Overlap | 0-3 | Has 2+ shared ID types with other domains |
| Open License | 0-3 | CC0/CC-BY/Public Domain (commercial OK) |
| Bulk Download | 0-3 | Available without API-only restrictions |

### Tier 2: Strong Preference (9 points max)
| Criterion | Points | Description |
|-----------|--------|-------------|
| Canonical Status | 0-3 | THE go-to source in its domain |
| Active Maintenance | 0-3 | Updated within 12 months |
| Rich Cross-References | 0-3 | Links to many other databases |

### Tier 3: Nice-to-Have (9 points max)
| Criterion | Points | Description |
|-----------|--------|-------------|
| RDF/Semantic Format | 0-3 | OWL, RDF, SPARQL available |
| Versioned Releases | 0-3 | Clear version tracking |
| Reasonable Size | 0-3 | Manageable for MVP scope |

**Maximum Score**: 27 points

## Summary by Category

| Category | Sources | Top Score | Key Hub IDs |
|----------|---------|-----------|-------------|
| [01.genetics](./01-genetics.md) | 9 | 26/27 | rsID, Gene ID, OMIM |
| [02.compounds](./02-compounds.md) | 10 | 26/27 | PubChem CID, InChIKey, ChEBI |
| [03.diseases](./03-diseases.md) | 6 | 26/27 | MONDO, HPO, EFO |
| [04.pathways](./04-pathways.md) | 8 | 27/27 | GO, UniProt, Reactome |
| [05.traditional-medicine](./05-traditional-medicine.md) | 4 | 21/27 | PubChem CID, UniProt |
| [06.nutrition](./06-nutrition.md) | 4 | 22/27 | HMDB, PubChem CID |
| [07.proteins](./07-proteins.md) | 3 | 26/27 | UniProt (central hub) |
| [08.literature](./08-literature.md) | 8 | 26/27 | PMID, DOI, Wikidata QID |
| [09.microbiome](./09-microbiome.md) | 3 | 24/27 | NCBI Taxonomy, UniProt |

## Top Scoring Sources (25+ points)

| Score | Source | Category | License |
|-------|--------|----------|---------|
| **27/27** | Gene Ontology | Pathways | CC BY 4.0 |
| 26/27 | ClinVar | Genetics | CC0 |
| 26/27 | CIViC | Genetics | CC0 |
| 26/27 | LOTUS | Compounds | CC0 |
| 26/27 | Open Targets | Diseases | CC BY 4.0 |
| 26/27 | UniProt | Proteins | CC BY 4.0 |
| 26/27 | PDB | Proteins | CC0 |
| 26/27 | Wikidata | Literature | CC0 |
| 26/27 | OpenAlex | Literature | CC0 |
| 25/27 | dbSNP | Genetics | CC0 |
| 25/27 | GWAS Catalog | Genetics | CC BY 4.0 |
| 25/27 | ChEBI | Compounds | CC BY 4.0 |
| 25/27 | MONDO | Diseases | CC BY 4.0 |
| 25/27 | Orphanet | Diseases | CC BY 4.0 |
| 25/27 | Reactome | Pathways | CC BY 4.0 |
| 25/27 | STRING | Pathways | CC BY 4.0 |
| 25/27 | IntAct | Pathways | CC BY 4.0 |
| 25/27 | JASPAR | Pathways | CC BY 4.0 |
| 25/27 | PubMed | Literature | CC0 |

## License Distribution

| License | Count | Commercial Use |
|---------|-------|----------------|
| CC0/Public Domain | 15 | Yes |
| CC BY 4.0 | 25 | Yes |
| CC BY-SA | 5 | Yes (ShareAlike) |
| Custom/Academic | 10 | Check terms |
| **CC BY-NC** | **0** | **Excluded** |

All selected sources permit commercial use.

## Quick Reference

```
01.GENETICS (9)
├── ClinVar, dbSNP, gnomAD, GWAS Catalog
├── dbNSFP, AlphaMissense, PharmGKB
└── GTEx, CIViC

02.COMPOUNDS (10)
├── ChEBI, PubChem, ChEMBL, LOTUS, COCONUT
├── DGIdb, BindingDB, GtoPdb
└── RxNorm, USDA FoodData

03.DISEASES (6)
├── MONDO, HPO, Open Targets
└── Orphanet, MeSH, GDC/TCGA

04.PATHWAYS (8)
├── Gene Ontology, Reactome, WikiPathways
├── STRING, IntAct, STITCH
└── JASPAR, MSigDB

05.TRAD.MEDICINE (4)
├── BATMAN-TCM 2.0, IMPPAT 2.0
└── KampoDB, HIT 2.0

06.NUTRITION (4)
├── FooDB, HMDB
└── DSLD, Exposome-Explorer

07.PROTEINS (3)
├── UniProt
└── PDB, AlphaFold DB

08.LITERATURE (8)
├── PubMed, Wikidata, OpenAlex
├── UniProt ID Mapping, ClinicalTrials.gov
└── FDA OpenFDA, PMC ID Converter, NCBI E-Link

09.MICROBIOME (3)
├── HMP, VMH
└── HOMD
```

## Related Documents

- [Integration Map](./integration-map.md) - How sources link via shared identifiers
- [Exclusions](./exclusions.md) - Sources excluded and rationale
- [Hub Sources](./hub-sources.md) - Critical identity bridges
