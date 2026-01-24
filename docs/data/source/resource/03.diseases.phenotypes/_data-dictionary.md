# Category 03: Diseases & Phenotypes - Data Dictionary

## Overview

This data dictionary documents all fields used across the Diseases & Phenotypes category, covering 22 data sources organized into 7 subcategories.

**Extraction Date:** 2026-01-24
**Total Sources:** 22

## Subcategories

| ID | Name | Data Sources |
|----|------|--------------|
| 3.1 | Disease Ontologies | EFO, HPO, ICD, MeSH, MONDO, Orphanet/ORDO |
| 3.2 | Phenotype Databases | HPO Annotations, OMIM, Orphanet Phenotypes |
| 3.3 | Disease-Gene Associations & Pharmacogenomics | DisGeNET, Open Targets, PharmGKB |
| 3.4 | Cancer/Oncology & Drug-Disease Targets | GDC/TCGA, DGIdb |
| 3.5 | Rare/Orphan Diseases | DECIPHER, Orphanet, PanelApp |
| 3.6 | Autoimmune/Inflammatory Diseases | ImmunoBase, IPD-IMGT/HLA |
| 3.7 | Mental Health/Neurological | Allen Brain Atlas, PGC, SynGO |

---

## Cross-Category Identifiers

### Disease Identifiers

| Identifier | Pattern | Description | Used By |
|------------|---------|-------------|---------|
| MONDO | `MONDO:[0-9]{7}` | Unified disease ontology (pivot identifier) | MONDO, HPO, DisGeNET, Orphanet |
| EFO | `EFO_[0-9]{7}` | Experimental Factor Ontology (GWAS/Open Targets) | EFO, Open Targets, GWAS Catalog |
| OMIM | `[*#%+]?[0-9]{6}` | Mendelian disease identifier | OMIM, HPO, MONDO, Orphanet, DisGeNET |
| Orphanet | `Orphanet:[0-9]+` | Rare disease identifier | Orphanet, MONDO, HPO, DECIPHER |
| ICD-10 | `[A-Z][0-9]{2}(\.[0-9]{1,4})?` | WHO disease classification | ICD, Orphanet, MONDO |
| UMLS CUI | `C[0-9]{7}` | UMLS Concept Unique Identifier | DisGeNET, MeSH, HPO |

### Phenotype Identifiers

| Identifier | Pattern | Description | Used By |
|------------|---------|-------------|---------|
| HPO | `HP:[0-9]{7}` | Human Phenotype Ontology term | HPO, Orphanet, DECIPHER, DisGeNET |

### Gene Identifiers

| Identifier | Pattern | Description | Used By |
|------------|---------|-------------|---------|
| NCBI Gene ID | `[0-9]+` | NCBI Entrez Gene identifier | DisGeNET, HPO, PharmGKB |
| Ensembl Gene ID | `ENSG[0-9]{11}` | Ensembl gene identifier | Open Targets, GDC/TCGA, DECIPHER |
| HGNC Symbol | `[A-Z0-9]+` | HUGO gene symbol | All sources |
| HGNC ID | `HGNC:[0-9]+` | HUGO Gene Nomenclature Committee ID | Orphanet, PanelApp, PharmGKB |

### Variant Identifiers

| Identifier | Pattern | Description | Used By |
|------------|---------|-------------|---------|
| dbSNP rsID | `rs[0-9]+` | dbSNP reference SNP identifier | DisGeNET, ImmunoBase, PGC, PharmGKB |

---

## Common Field Patterns

### Evidence Scores

Numeric confidence/evidence scores typically ranging 0-1, used by:
- DisGeNET (GDA score)
- Open Targets (overall score)
- PharmGKB (evidence level)

### Cross-References

Links to external databases and ontologies, common targets include:
- OMIM
- MONDO
- HPO
- UMLS
- MeSH
- SNOMED CT
- ICD-10

### Publication References

Literature citations supporting annotations, formats:
- PMID integers
- Full citations

### Frequency Annotations

How often a phenotype/variant occurs, formats:
- HPO frequency terms
- Percentages
- Fractions (n/N)

---

## Subcategory Data Dictionaries

- [3.1 Disease Ontologies](./_data-dictionary.md)
- [3.2 Phenotype Databases](./3.2.phenotype.databases/_data-dictionary.md)
- [3.3 Disease-Gene Associations](./3.3.disease.gene.associations/_data-dictionary.md)
- [3.4 Cancer/Oncology](./3.4.cancer.oncology/_data-dictionary.md)
- [3.5 Rare/Orphan Diseases](./3.5.rare.orphan.diseases/_data-dictionary.md)
- [3.6 Autoimmune/Inflammatory](./3.6.autoimmune.inflammatory/_data-dictionary.md)
- [3.7 Mental Health/Neurological](./3.7.mental.health.neurological/_data-dictionary.md)
