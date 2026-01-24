---
id: schema-metahit
title: "MetaHIT Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: final
tags: [schema, database, metagenomics, gene-catalog, reference]
---

# MetaHIT Schema Documentation

**Document ID:** SCHEMA-METAHIT
**Version:** Legacy (Project Completed)
**Source Version:** IGC Gene Catalog v1

---

## TL;DR

MetaHIT (Metagenomics of the Human Intestinal Tract) established foundational gut microbiome resources including the Integrated Gene Catalog (IGC) with 10M+ microbial genes. Data includes reference gene catalogs, sample abundance profiles, enterotype classifications, and phenotype associations. Project completed but data remains a key reference.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| IGC genes | 9,879,896 | MetaHIT Publications |
| Original samples | 1,267 | MetaHIT Consortium |
| Countries | 8 European | MetaHIT Consortium |
| Enterotypes | 3 (Bacteroides, Prevotella, Ruminococcus) | Nature 2011 |
| Publications | 60+ | PubMed |

---

## Entity Relationship Overview

```
MetaHIT Data
  ├── Gene Catalog (IGC)
  │     ├── Gene ID
  │     ├── Sequence
  │     ├── Taxonomy assignment
  │     └── Functional annotation
  ├── Sample Profiles
  │     ├── Sample ID
  │     ├── Gene abundances
  │     └── Metadata
  └── Phenotype Associations
        ├── Enterotype
        ├── BMI/Obesity
        └── IBD status
```

---

## Core Tables/Entities

### IGC Gene Entry

**Description:** Integrated Gene Catalog entry.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| gene_id | string | Yes | IGC gene identifier |
| sequence | string | Yes | Nucleotide sequence |
| length | integer | Yes | Gene length (bp) |
| taxonomy | string | No | Taxonomic assignment |
| kegg_ko | string | No | KEGG Ortholog ID |
| eggnog | string | No | eggNOG annotation |
| cog | string | No | COG category |
| sample_origin | string | No | Source sample ID |

### Sample Profile

**Description:** Per-sample gene abundance data.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| sample_id | string | Yes | MetaHIT sample ID |
| gene_id | string | Yes | IGC gene identifier |
| abundance | float | Yes | Relative abundance |
| coverage | float | No | Read coverage |
| read_count | integer | No | Mapped reads |

### Sample Metadata

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| sample_id | string | Yes | MetaHIT sample ID |
| country | string | Yes | Country of origin |
| age | integer | No | Subject age |
| sex | string | No | male/female |
| bmi | float | No | Body mass index |
| ibd_status | string | No | CD/UC/healthy |
| enterotype | string | No | Enterotype classification |
| sequencing_depth | integer | No | Total reads |

---

## Gene Catalog Structure

### IGC Gene ID Format

```
Format: MH{sample_number}_GL{gene_number}
Example: MH0001_GL0000001

Components:
- MH: MetaHIT prefix
- {sample_number}: 4-digit sample origin
- GL: Gene Locus
- {gene_number}: 7-digit gene index
```

### Taxonomic Assignment

| Field | Description |
|-------|-------------|
| kingdom | Domain (Bacteria, Archaea) |
| phylum | Major phylum |
| class | Taxonomic class |
| order | Taxonomic order |
| family | Taxonomic family |
| genus | Genus assignment |
| species | Species (if confident) |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | FASTA (gene sequences) |
| Alternative | GFF, TSV (annotations) |
| Abundance | TSV matrix |
| Encoding | UTF-8 |
| Compression | gzip |

---

## Sample Record

### Gene Catalog Entry (FASTA)

```
>MH0001_GL0000001 length=843 taxonomy=Bacteroides kegg=K00001 eggnog=COG0001
ATGAAACGTGCAGCTGATCGTGATCGATCGATCGATCGATCGATCGATCG...
```

### Abundance Table (TSV)

```
gene_id	MH0001	MH0002	MH0003	MH0004
MH0001_GL0000001	0.00012	0.00008	0.00015	0.00003
MH0001_GL0000002	0.00005	0.00012	0.00002	0.00018
MH0001_GL0000003	0.00000	0.00001	0.00000	0.00005
```

### Metadata (TSV)

```
sample_id	country	age	sex	bmi	ibd_status	enterotype
MH0001	Denmark	45	male	24.5	healthy	Bacteroides
MH0002	Spain	38	female	22.1	healthy	Prevotella
MH0003	France	52	male	31.2	CD	Ruminococcus
```

---

## Enterotype Classification

| Enterotype | Dominant Genus | Characteristics |
|------------|---------------|-----------------|
| Type 1 | Bacteroides | High protein/fat diet |
| Type 2 | Prevotella | High fiber/carbohydrate |
| Type 3 | Ruminococcus | Variable, often Firmicutes |

---

## Functional Annotations

### KEGG Ortholog Groups

| Field | Description |
|-------|-------------|
| ko_id | KEGG Ortholog ID (K00001) |
| ko_name | Ortholog name |
| pathway | Associated pathways |
| module | KEGG module |

### eggNOG Categories

| Category | Description |
|----------|-------------|
| J | Translation |
| K | Transcription |
| L | Replication |
| C | Energy metabolism |
| G | Carbohydrate metabolism |
| E | Amino acid metabolism |
| S | Function unknown |

---

## Data Access Points

| Resource | URL | Content |
|----------|-----|---------|
| IGC Catalog | http://meta.genomics.cn | Gene sequences |
| ENA | https://www.ebi.ac.uk/ena | Raw reads |
| Publications | PubMed | Processed data |

---

## Cross-References

| Database | ID Type | Usage |
|----------|---------|-------|
| NCBI SRA | ERR/SRR accessions | Raw sequences |
| KEGG | KO IDs | Functional annotation |
| eggNOG | COG/NOG | Ortholog groups |
| NCBI Taxonomy | TaxID | Species assignment |

---

## Glossary

| Term | Definition |
|------|------------|
| IGC | Integrated Gene Catalog |
| Enterotype | Microbiome community type |
| Gene richness | Number of genes detected |
| Metagenome | Collective genome of community |
| Core microbiome | Genes shared across samples |

---

## Key Publications

| Year | Title | DOI |
|------|-------|-----|
| 2010 | A human gut microbial gene catalogue | 10.1038/nature08821 |
| 2011 | Enterotypes of the human gut microbiome | 10.1038/nature09944 |
| 2013 | Richness of human gut microbiome correlates with metabolic markers | 10.1038/nature12506 |
| 2014 | An integrated catalog of reference genes in the human gut microbiome | 10.1038/nbt.2942 |

---

## References

1. Qin J, et al. (2010). A human gut microbial gene catalogue established by metagenomic sequencing. Nature. https://doi.org/10.1038/nature08821
2. Li J, et al. (2014). An integrated catalog of reference genes in the human gut microbiome. Nature Biotechnology. https://doi.org/10.1038/nbt.2942
3. MetaHIT Gene Catalog: http://meta.genomics.cn
