# Category 08: Literature & Knowledge - Data Dictionary

**Category ID:** 08
**Category Name:** Literature & Knowledge
**Schema Version:** 1.0
**Last Updated:** 2026-01-24

## Overview

This category encompasses data sources related to scientific publications, knowledge bases, identifier mapping services, and regulatory/legal information. It provides unified access to literature databases, structured knowledge graphs, cross-reference mapping tools, and clinical/regulatory data.

## Subcategories

| ID | Name | Data Sources |
|----|------|--------------|
| 8.1 | Scientific Literature | PubMed, PubMed Central (PMC), Europe PMC, OpenAlex, Semantic Scholar |
| 8.2 | Knowledge Bases | Wikidata, Wikipedia |
| 8.3 | Identifier Mapping | NCBI E-Link, PMC ID Converter, UniProt ID Mapping |
| 8.4 | Regulatory & Legal | ClinicalTrials.gov, FDA OpenFDA |

---

## Cross-Reference Identifier Types

### Publication Identifiers

| Identifier | Format | Description | Sources |
|------------|--------|-------------|---------|
| PMID | `[0-9]+` | PubMed identifier | PubMed, PMC, Europe PMC, Semantic Scholar |
| PMCID | `PMC[0-9]+` | PubMed Central identifier | PMC, PubMed, Europe PMC, Semantic Scholar |
| DOI | `10\.[0-9]+/.+` | Digital Object Identifier | All literature sources, PMC ID Converter |
| ArXiv | `[0-9]+\.[0-9]+` | ArXiv preprint identifier | Semantic Scholar, OpenAlex |
| MAG ID | `[0-9]+` | Microsoft Academic Graph ID | OpenAlex, Semantic Scholar |

### Gene Identifiers

| Identifier | Format | Wikidata Property | Description | Sources |
|------------|--------|-------------------|-------------|---------|
| Entrez Gene ID | `[0-9]+` | P351 | NCBI Gene identifier | Wikidata, NCBI E-Link, UniProt ID Mapping |
| HGNC Symbol | `[A-Z0-9]+` | P353 | HUGO gene symbol | Wikidata, UniProt ID Mapping |
| HGNC ID | `HGNC:[0-9]+` | P354 | HGNC identifier | Wikidata, UniProt ID Mapping |
| Ensembl Gene ID | `ENSG[0-9]{11}` | P594 | Ensembl gene identifier | Wikidata, UniProt ID Mapping, OpenAlex |

### Protein Identifiers

| Identifier | Format | Wikidata Property | Description | Sources |
|------------|--------|-------------------|-------------|---------|
| UniProt Accession | `[A-Z][0-9][A-Z0-9]{3}[0-9]` | P352 | UniProtKB accession | Wikidata, UniProt ID Mapping |
| RefSeq Protein | `[NX]P_[0-9]+\.[0-9]+` | P639 | NCBI RefSeq protein | Wikidata, UniProt ID Mapping |
| PDB ID | `[0-9][A-Z0-9]{3}` | - | Protein Data Bank structure | UniProt ID Mapping, Wikidata |

### Disease Identifiers

| Identifier | Format | Wikidata Property | Description | Sources |
|------------|--------|-------------------|-------------|---------|
| DOID | `DOID:[0-9]+` | P699 | Disease Ontology ID | Wikidata |
| OMIM | `[0-9]{6}` | P492 | OMIM ID | Wikidata, UniProt ID Mapping |
| Orphanet | `[0-9]+` | P1395 | Orphanet ID | Wikidata |
| MONDO | `MONDO:[0-9]+` | P5270 | MONDO disease ontology | Wikidata |
| MeSH | `[A-Z][0-9]{6,9}` | P486 | MeSH descriptor | Wikidata, PubMed |

### Drug Identifiers

| Identifier | Format | Wikidata Property | Description | Sources |
|------------|--------|-------------------|-------------|---------|
| RxCUI | `[0-9]+` | - | RxNorm Concept Unique Identifier | FDA OpenFDA, UniProt ID Mapping |
| NDC | `[0-9]{4,5}-[0-9]{3,4}-[0-9]{1,2}` | - | National Drug Code | FDA OpenFDA |
| DrugBank ID | `DB[0-9]+` | P715 | DrugBank identifier | Wikidata, UniProt ID Mapping |
| ChEMBL ID | `CHEMBL[0-9]+` | P592 | ChEMBL compound ID | Wikidata, UniProt ID Mapping |
| PubChem CID | `[0-9]+` | P662 | PubChem compound ID | Wikidata |
| UNII | `[A-Z0-9]{10}` | - | Unique Ingredient Identifier | FDA OpenFDA |

### Clinical Trial Identifiers

| Identifier | Format | Description | Sources |
|------------|--------|-------------|---------|
| NCT ID | `NCT[0-9]{8}` | ClinicalTrials.gov identifier | ClinicalTrials.gov |
| EudraCT | `[0-9]{4}-[0-9]{6}-[0-9]{2}` | EU Clinical Trials Register ID | ClinicalTrials.gov |
| NDA/ANDA/BLA | `(NDA\|ANDA\|BLA)[0-9]+` | FDA application number | FDA OpenFDA |

### Knowledge Base Identifiers

| Identifier | Format | Description | Sources |
|------------|--------|-------------|---------|
| Wikidata QID | `Q[0-9]+` | Wikidata item ID | Wikidata, Wikipedia, OpenAlex |
| Wikipedia Page ID | `[0-9]+` | Wikipedia page identifier | Wikipedia |
| DBpedia URI | `http://dbpedia.org/resource/.+` | DBpedia resource URI | Wikipedia |

---

## Relationship Types

### Literature Relationships

| Relationship | From | To | Cardinality | Sources |
|--------------|------|----|-----------:|---------|
| cites | Article | Article | N:M | PubMed, Semantic Scholar, OpenAlex |
| related_to | Article | Article | N:M | Semantic Scholar, OpenAlex, NCBI E-Link |
| authored_by | Article | Author | N:M | All literature sources |
| published_in | Article | Journal | N:1 | All literature sources |
| indexed_with | Article | MeSH Term | N:M | PubMed, Europe PMC |
| annotates | Annotation | Article | N:1 | Europe PMC |

### ID Mapping Relationships

| Relationship | From | To | Cardinality | Sources |
|--------------|------|----|-----------:|---------|
| maps_to | Identifier | Identifier | N:M | All mapping services |
| links_to | Database Record | Database Record | N:M | NCBI E-Link |

### Knowledge Graph Relationships

| Relationship | From | To | Cardinality | Sources |
|--------------|------|----|-----------:|---------|
| instance_of | Entity | Class | N:M | Wikidata |
| subclass_of | Class | Class | N:M | Wikidata |
| found_in_taxon | Gene/Protein | Taxon | N:M | Wikidata |
| genetic_association | Gene | Disease | N:M | Wikidata |

### Regulatory Relationships

| Relationship | From | To | Cardinality | Sources |
|--------------|------|----|-----------:|---------|
| tests_drug | Clinical Trial | Drug | N:M | ClinicalTrials.gov |
| studies_condition | Clinical Trial | Condition | N:M | ClinicalTrials.gov |
| caused_by | Adverse Event | Drug | N:M | FDA OpenFDA |
| recalled | Recall | Product | N:M | FDA OpenFDA |

---

## Subcategory Data Dictionaries

- [8.1 Scientific Literature Data Dictionary](./8.1.scientific.literature/_data-dictionary.md)
- [8.2 Knowledge Bases Data Dictionary](./8.2.knowledge.bases/_data-dictionary.md)
- [8.3 Identifier Mapping Data Dictionary](./8.3.identifier.mapping/_data-dictionary.md)
- [8.4 Regulatory & Legal Data Dictionary](./8.4.regulatory.legal/_data-dictionary.md)
