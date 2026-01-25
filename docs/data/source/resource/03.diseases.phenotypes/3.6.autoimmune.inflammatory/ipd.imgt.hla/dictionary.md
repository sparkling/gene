# IPD-IMGT/HLA - Data Dictionary

## Overview

This data dictionary documents the schema for IPD-IMGT/HLA (Immuno Polymorphism Database - HLA).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | ipd.imgt.hla |
| **Name** | IPD-IMGT/HLA |
| **Parent** | 3.6.autoimmune.inflammatory |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### HLA Allele

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| allele_name | string | 1:1 | Yes | Full allele name | HLA-A*02:01:01:01 |
| gene | string | 1:1 | Yes | HLA gene | HLA-A |
| allele_group | string | 1:1 | Yes | Allele group | *02 |
| protein | string | 1:1 | No | Specific protein | :01 |
| synonymous | string | 1:1 | No | Coding synonymous | :01 |
| noncoding | string | 1:1 | No | Non-coding diff | :01 |
| ipd_accession | string | 1:1 | Yes | IPD accession | HLA00001 |
| date_assigned | date | 1:1 | Yes | Assignment date | 2005-01-15 |
| status | enum | 1:1 | Yes | Allele status | current |

### Sequence Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ipd_accession | string | 1:1 | Yes | IPD accession | HLA00001 |
| nucleotide_seq | string | 1:1 | Yes | DNA sequence | ATGCGTCA... |
| protein_seq | string | 1:1 | No | Amino acid sequence | MAVMAPRTL... |
| exon_boundaries | array | 1:N | No | Exon positions | [{start: 1, end: 73}] |
| cds_length | integer | 1:1 | No | Coding sequence length | 1098 |
| genomic_length | integer | 1:1 | No | Full genomic length | 3503 |

### Gene Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_symbol | string | 1:1 | Yes | HLA gene name | HLA-A |
| locus | string | 1:1 | Yes | MHC locus | Class I |
| chromosome | string | 1:1 | Yes | Chromosome | 6 |
| region | string | 1:1 | Yes | MHC region | 6p21.3 |
| total_alleles | integer | 1:1 | No | Number of alleles | 7000 |
| function | string | 1:1 | No | Gene function | Antigen presentation |

### Serological Equivalents

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| allele_name | string | 1:1 | Yes | HLA allele | HLA-A*02:01 |
| serological_type | string | 1:1 | No | Serological antigen | A2 |
| split | string | 1:1 | No | Split antigen | A2 |
| broad | string | 1:1 | No | Broad antigen | A2 |
| public | array | 1:N | No | Public epitopes | [Bw4, Bw6] |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| HLA Allele | HLA-[A-Z]+[0-9]*\*[0-9:]+[NLSQCA]? | HLA-A*02:01:01:01 | Full nomenclature |
| Allele Group | HLA-[A-Z]+\*[0-9]+ | HLA-A*02 | Serological equiv |
| IPD Accession | HLA[0-9]+ | HLA00001 | Database accession |
| Gene Symbol | HLA-[A-Z]+[0-9]* | HLA-DRB1 | Gene name |
| IMGT Accession | [A-Z]{1,2}[0-9]{5,6} | AJ567890 | EMBL accession |

---

## Enumerations

### HLA Nomenclature Fields

| Field | Example | Meaning |
|-------|---------|---------|
| Gene | HLA-A | Gene locus |
| Allele Group | *02 | Serological equivalent |
| Specific Protein | :01 | Unique protein sequence |
| Synonymous | :01 | Synonymous coding change |
| Non-coding | :01 | Non-coding difference |
| Suffix | N, L, S, Q, C, A | Expression status |

### Expression Suffixes

| Suffix | Meaning | Description |
|--------|---------|-------------|
| N | Null | No cell surface expression |
| L | Low | Low expression |
| S | Secreted | Soluble, secreted |
| Q | Questionable | Expression uncertain |
| C | Cytoplasm | Cytoplasmic only |
| A | Aberrant | Aberrant expression |

### HLA Gene Classes

| Class | Genes | Function |
|-------|-------|----------|
| Class I | HLA-A, HLA-B, HLA-C | Cytotoxic T cell response |
| Class II | HLA-DR, HLA-DQ, HLA-DP | Helper T cell response |
| Class III | C4, BF | Complement components |
| Non-classical | HLA-E, HLA-F, HLA-G | Specialized functions |

### MHC Loci

| Locus | Genes | Description |
|-------|-------|-------------|
| HLA-A | HLA-A | Class I heavy chain |
| HLA-B | HLA-B | Class I heavy chain |
| HLA-C | HLA-C | Class I heavy chain |
| HLA-DRA | HLA-DRA | DR alpha chain |
| HLA-DRB1 | HLA-DRB1 | DR beta 1 chain |
| HLA-DRB3/4/5 | HLA-DRB3, DRB4, DRB5 | DR beta secondary |
| HLA-DQA1 | HLA-DQA1 | DQ alpha chain |
| HLA-DQB1 | HLA-DQB1 | DQ beta chain |
| HLA-DPA1 | HLA-DPA1 | DP alpha chain |
| HLA-DPB1 | HLA-DPB1 | DP beta chain |

### Allele Status

| Status | Description |
|--------|-------------|
| current | Active allele |
| deleted | Removed allele |
| renamed | Renamed to new name |
| merged | Merged with another |

### Sequence Types

| Type | Description |
|------|-------------|
| genomic | Full genomic sequence |
| cds | Coding sequence only |
| protein | Translated protein |
| exon | Individual exon |

---

## Entity Relationships

### Gene to Alleles
- **Cardinality:** 1:N
- **Description:** Each gene has many alleles
- **Key Fields:** gene_symbol, allele_name

### Allele to Sequences
- **Cardinality:** 1:1
- **Description:** Each allele has one sequence
- **Key Fields:** ipd_accession, nucleotide_seq

### Allele to Serological
- **Cardinality:** N:1
- **Description:** Many alleles map to one serotype
- **Key Fields:** allele_name, serological_type

### Gene to Locus
- **Cardinality:** N:1
- **Description:** Genes belong to MHC class
- **Key Fields:** gene_symbol, locus

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HLA | Human Leukocyte Antigen | Gene system |
| MHC | Major Histocompatibility Complex | Genomic region |
| IPD | Immuno Polymorphism Database | Database name |
| IMGT | ImMunoGeneTics | Database system |
| WHO | World Health Organization | Nomenclature authority |
| EBI | European Bioinformatics Institute | Host institution |
| EMBL | European Molecular Biology Laboratory | Sequence database |
| CTL | Cytotoxic T Lymphocyte | Immune cell |
| TCR | T Cell Receptor | Immune receptor |
| APC | Antigen Presenting Cell | Immune cell |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| EMBL/ENA | Accession | Sequence source |
| UniProt | Accession | Protein data |
| RefSeq | Accession | Reference sequences |
| Ensembl | Gene ID | Gene annotation |
| HGNC | Gene symbol | Gene nomenclature |
| PDB | Structure ID | 3D structures |
| dbSNP | rsID | Variants |

---

## Data Quality Notes

1. **Total Alleles:** 30,000+ HLA alleles
2. **HLA-A Alleles:** 7,000+ alleles
3. **HLA-B Alleles:** 8,500+ alleles
4. **HLA-C Alleles:** 7,000+ alleles
5. **Class II Alleles:** Extensive DRB1, DQB1 coverage
6. **WHO Nomenclature:** Official naming authority
7. **Quarterly Updates:** Regular releases
8. **REST API:** Programmatic access available

