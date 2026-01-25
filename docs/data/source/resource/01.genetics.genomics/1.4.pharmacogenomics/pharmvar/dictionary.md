# PharmVar - Data Dictionary

## Overview

This data dictionary documents the schema for PharmVar (Pharmacogene Variation Consortium) allele database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | pharmvar |
| **Name** | PharmVar |
| **Parent** | 1.4.pharmacogenomics |
| **Total Fields** | 14 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Gene

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_symbol | string | 1:1 | Yes | HGNC symbol | CYP2D6 |
| gene_name | string | 1:1 | Yes | Full gene name | cytochrome P450 2D6 |
| chromosome | string | 1:1 | Yes | Chromosome location | 22q13.2 |
| refseq_gene | string | 1:1 | Yes | RefSeq gene ID | NG_008376.4 |
| refseq_mrna | string | 1:1 | Yes | RefSeq mRNA ID | NM_000106.6 |
| refseq_protein | string | 1:1 | Yes | RefSeq protein ID | NP_000097.3 |

### Allele (Haplotype)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| allele_name | string | 1:1 | Yes | Star allele name | *4 |
| full_name | string | 1:1 | Yes | Full designation | CYP2D6*4 |
| function | string | 1:1 | Yes | Functional status | No function |
| activity_value | float | 1:1 | No | Activity score | 0 |
| core_variants | array | 1:N | Yes | Defining variants | ["c.506-1G>A", "c.1846G>A"] |
| suballeles | array | 1:N | No | Suballele variants | ["*4.001", "*4.002"] |

### Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| hgvs_c | string | 1:1 | Yes | Coding HGVS | c.506-1G>A |
| hgvs_p | string | 1:1 | No | Protein HGVS | p.Arg296Cys |
| rsid | string | 1:1 | No | dbSNP identifier | rs3892097 |
| position_grch38 | string | 1:1 | Yes | GRCh38 coordinates | chr22:42128174 |
| reference_allele | string | 1:1 | Yes | Reference nucleotide | G |
| variant_allele | string | 1:1 | Yes | Variant nucleotide | A |
| variant_type | string | 1:1 | Yes | Variant classification | SNV, deletion |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Star Allele | Gene*number | CYP2D6*4 | Haplotype designation |
| Suballele | Gene*number.number | CYP2D6*4.001 | Suballele variant |
| rsID | rs + digits | rs3892097 | dbSNP identifier |
| RefSeq Gene | NG_ + digits | NG_008376.4 | Gene reference |
| RefSeq mRNA | NM_ + digits | NM_000106.6 | Transcript reference |
| RefSeq Protein | NP_ + digits | NP_000097.3 | Protein reference |
| HGVS | transcript:change | NM_000106.6:c.506-1G>A | Nomenclature |

---

## Enumerations

### Functional Status

| Status | Description | Activity Value |
|--------|-------------|----------------|
| No function | Complete loss of activity | 0 |
| Decreased function | Reduced activity | 0.5 |
| Normal function | Wild-type activity | 1 |
| Increased function | Enhanced activity | >1 |
| Uncertain function | Unknown/variable | null |

### Variant Type

| Type | Description |
|------|-------------|
| SNV | Single nucleotide variant |
| insertion | Nucleotide insertion |
| deletion | Nucleotide deletion |
| indel | Insertion-deletion |
| duplication | Gene/exon duplication |
| deletion/conversion | Hybrid/deletion allele |

### Allele Evidence

| Level | Meaning |
|-------|---------|
| Definitive | Well-established function |
| Moderate | Good evidence for function |
| Limited | Limited functional data |
| In vitro only | Only cell-based evidence |

---

## Entity Relationships

### Gene to Allele
- **Cardinality:** 1:N
- **Description:** Multiple star alleles per gene
- **Key Fields:** gene_symbol, allele_name

### Allele to Suballele
- **Cardinality:** 1:N
- **Description:** Core allele has suballele variants
- **Key Fields:** allele_name, suballele_name

### Allele to Variant
- **Cardinality:** 1:N
- **Description:** Multiple variants define one allele
- **Key Fields:** allele_name, variant_id

### Variant to Position
- **Cardinality:** 1:N
- **Description:** Variant mapped to multiple assemblies
- **Key Fields:** variant_id, assembly

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PharmVar | Pharmacogene Variation Consortium | Database name |
| PGx | Pharmacogenomics | Field abbreviation |
| HGVS | Human Genome Variation Society | Nomenclature standard |
| RefSeq | Reference Sequence | NCBI sequence database |
| SNV | Single Nucleotide Variant | Variant type |
| CNV | Copy Number Variant | Gene duplications/deletions |
| CYP | Cytochrome P450 | Enzyme family |
| GRCh38 | Genome Reference Consortium human 38 | Current assembly |
| GRCh37 | Genome Reference Consortium human 37 | Previous assembly |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant identifier |
| RefSeq | NG/NM/NP | Sequence reference |
| PharmGKB | PA accession | Clinical annotations |
| CPIC | Allele | Guideline reference |
| ClinVar | VCV | Clinical significance |
| gnomAD | Variant ID | Population frequency |

---

## Covered Genes

| Gene | Allele Count | Primary Function |
|------|--------------|------------------|
| CYP2D6 | 150+ | Drug metabolism |
| CYP2C19 | 40+ | Drug metabolism |
| CYP2C9 | 70+ | Drug metabolism |
| CYP3A4 | 50+ | Drug metabolism |
| CYP3A5 | 15+ | Drug metabolism |
| DPYD | 30+ | Fluoropyrimidine metabolism |
| TPMT | 40+ | Thiopurine metabolism |
| UGT1A1 | 130+ | Glucuronidation |
| NAT2 | 100+ | Acetylation |
| SLCO1B1 | 45+ | Drug transport |

---

## Data Quality Notes

1. **Cardinality:** One allele definition per star allele designation
2. **Nomenclature:** Star allele system is gene-specific standard
3. **Suballeles:** Represent minor sequence variations within core allele
4. **Evidence Levels:** Functional assignments based on published evidence
5. **Updates:** Regular curation as new alleles discovered
6. **Predecessor:** Replaces individual gene nomenclature committees
