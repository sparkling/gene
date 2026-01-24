# Category 04: Pathways & Networks - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| Category ID | 04 |
| Category Name | Pathways & Networks |
| Total Subcategories | 6 |
| Total Data Sources | 12 |
| Extraction Date | 2026-01-24 |

## Subcategories

| ID | Name | Data Sources |
|----|------|--------------|
| 4.1 | Metabolic Pathways | KEGG, Reactome, WikiPathways |
| 4.2 | Signaling Pathways | Pathway Commons |
| 4.3 | Protein-Protein Interactions | BioGRID, IntAct, STRING |
| 4.4 | Drug-Target Interactions | STITCH |
| 4.5 | Gene Function & Ontology | Gene Ontology, MSigDB |
| 4.6 | Regulatory Networks | JASPAR, Roadmap Epigenomics |

## Cross-Category Shared Identifiers

### Gene Identifiers

| Type | Format | Sources |
|------|--------|---------|
| UniProt | `[A-Z][0-9]{5}` or `[A-Z][0-9][A-Z0-9]{3}[0-9]` | Reactome, IntAct, Gene Ontology, JASPAR |
| Ensembl Gene | `ENSG[0-9]{11}` | STRING, WikiPathways, Gene Ontology |
| Ensembl Protein | `ENSP[0-9]{11}` | STRING, STITCH, Reactome |
| NCBI Gene ID | Integer | BioGRID, WikiPathways, Gene Ontology, MSigDB |
| HGNC Symbol | Text | BioGRID, IntAct, STRING, MSigDB |

### Pathway Identifiers

| Type | Format | Sources |
|------|--------|---------|
| KEGG Pathway | `hsa[0-9]{5}` or `map[0-9]{5}` | KEGG, Reactome, WikiPathways, STRING |
| Reactome | `R-[A-Z]{3}-[0-9]+` | Reactome, WikiPathways, Pathway Commons, STRING |
| WikiPathways | `WP[0-9]+` | WikiPathways, Pathway Commons, STRING |

### Ontology Terms

| Type | Format | Sources |
|------|--------|---------|
| GO Term | `GO:[0-9]{7}` | Gene Ontology, Reactome, BioGRID, IntAct, STRING |
| ChEBI | `CHEBI:[0-9]+` | Reactome, WikiPathways, STITCH |
| PSI-MI | `MI:[0-9]{4}` | IntAct, BioGRID, STRING |

### Taxonomy Identifiers

| Tax ID | Species |
|--------|---------|
| 9606 | Homo sapiens |
| 10090 | Mus musculus |
| 10116 | Rattus norvegicus |
| 7227 | Drosophila melanogaster |
| 6239 | Caenorhabditis elegans |
| 4932 | Saccharomyces cerevisiae |

## Data Format Summary

| Subcategory | Formats |
|-------------|---------|
| 4.1 Metabolic Pathways | XML (KGML, GPML), JSON, BioPAX (OWL), SBML |
| 4.2 Signaling Pathways | BioPAX Level 3 (OWL), SIF (TSV), JSON-LD, SBGN-ML |
| 4.3 Protein-Protein Interactions | PSI-MITAB (TSV), PSI-MI XML, JSON, TSV |
| 4.4 Drug-Target Interactions | TSV, JSON |
| 4.5 Gene Function & Ontology | OBO, OWL, GAF 2.2 (TSV), GMT, JSON |
| 4.6 Regulatory Networks | JASPAR/TRANSFAC/MEME (matrix), BED, bigWig, narrowPeak |

## Related Data Dictionaries

- [4.1 Metabolic Pathways](.//4.1.metabolic.pathways/_data-dictionary.md)
- [4.2 Signaling Pathways](./4.2.signaling.pathways/_data-dictionary.md)
- [4.3 Protein-Protein Interactions](./4.3.protein.protein.interactions/_data-dictionary.md)
- [4.4 Drug-Target Interactions](./4.4.drug.target.interactions/_data-dictionary.md)
- [4.5 Gene Function & Ontology](./4.5.gene.function.ontology/_data-dictionary.md)
- [4.6 Regulatory Networks](./4.6.regulatory.networks/_data-dictionary.md)
