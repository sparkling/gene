# gnomAD - Data Dictionary

## Overview

This data dictionary documents the schema for gnomAD (Genome Aggregation Database).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | gnomad |
| **Name** | gnomAD |
| **Parent** | 1.3.population.genetics |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Variant

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| variantId | string | 1:1 | Yes | chr-pos-ref-alt | 1-55516888-G-A |
| chrom | string | 1:1 | Yes | Chromosome | 1, X |
| pos | integer | 1:1 | Yes | Genomic position | 55516888 |
| ref | string | 1:1 | Yes | Reference allele | G |
| alt | string | 1:1 | Yes | Alternate allele | A |
| rsid | string | 1:1 | No | dbSNP identifier | rs12345 |

### Allele Counts (Global)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ac | integer | 1:1 | Yes | Allele count | 1234 |
| an | integer | 1:1 | Yes | Allele number | 1461870 |
| af | float | 1:1 | Yes | Allele frequency | 0.000844 |
| ac_hom | integer | 1:1 | Yes | Homozygote count | 2 |
| ac_hemi | integer | 1:1 | No | Hemizygote count (X/Y) | 5 |

### Population Frequencies

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| populations | array | 1:N | Yes | Per-population data | [{id: "nfe", ac: 800}] |
| pop.id | string | 1:1 | Yes | Population code | nfe, afr, eas |
| pop.ac | integer | 1:1 | Yes | Population allele count | 800 |
| pop.an | integer | 1:1 | Yes | Population allele number | 900000 |
| pop.af | float | 1:1 | Yes | Population frequency | 0.00089 |

### Quality Flags

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| flags | array | 1:N | No | QC flags | ["lcr", "segdup"] |
| filters | array | 1:N | Yes | Filter status | ["PASS"], ["RF"] |

### Gene Constraint

| Field Name | Data Type | Cardinality | Required | Description | Range |
|------------|-----------|-------------|----------|-------------|-------|
| pLI | float | 1:1 | No | Loss-of-function intolerance | 0-1 |
| oe_lof | float | 1:1 | No | Observed/expected LoF ratio | 0-2+ |
| oe_lof_lower | float | 1:1 | No | LOEUF lower bound | 0-2 |
| oe_lof_upper | float | 1:1 | No | LOEUF upper bound | 0-2+ |
| oe_mis | float | 1:1 | No | O/E missense ratio | 0-2+ |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Variant ID | chr-pos-ref-alt | 1-55516888-G-A | Primary key |
| Gene ID | ENSG + 11 digits | ENSG00000141510 | Ensembl gene |
| rsID | rs + digits | rs28934576 | dbSNP cross-ref |
| Transcript | ENST + 11 digits | ENST00000269305 | Ensembl transcript |

---

## Enumerations

### Genetic Ancestry Groups

| Code | Name | v4 Samples |
|------|------|------------|
| afr | African/African American | 63,000 |
| amr | Admixed American | 48,000 |
| asj | Ashkenazi Jewish | 15,000 |
| eas | East Asian | 27,000 |
| fin | Finnish | 22,000 |
| mid | Middle Eastern | 1,600 |
| nfe | Non-Finnish European | 418,000 |
| sas | South Asian | 36,000 |

### Quality Flags

| Flag | Meaning |
|------|---------|
| lcr | Low complexity region |
| segdup | Segmental duplication |
| decoy | Decoy sequence |
| par | Pseudoautosomal region |

### Filter Values

| Value | Meaning |
|-------|---------|
| PASS | Passed all filters |
| RF | Failed random forest |
| InbreedingCoeff | Failed inbreeding coefficient |
| AC0 | No alternate alleles |

### Constraint Categories

| pLI Range | Interpretation |
|-----------|----------------|
| >= 0.9 | Highly intolerant to LoF |
| 0.5 - 0.9 | Moderately intolerant |
| < 0.5 | Tolerant to LoF |

---

## Entity Relationships

### Variant to Population
- **Cardinality:** 1:N
- **Description:** Each variant has frequencies per population
- **Key Fields:** variantId, population.id

### Variant to Transcript Consequence
- **Cardinality:** 1:N
- **Description:** One variant may affect multiple transcripts
- **Key Fields:** variantId, transcript_id

### Gene to Constraint
- **Cardinality:** 1:1
- **Description:** One constraint metric set per gene
- **Key Fields:** gene_id

### Transcript to Gene
- **Cardinality:** N:1
- **Description:** Multiple transcripts per gene
- **Key Fields:** transcript_id, gene_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| gnomAD | Genome Aggregation Database | Database name |
| AF | Allele Frequency | Population metric |
| AC | Allele Count | Observed count |
| AN | Allele Number | Total chromosomes |
| LoF | Loss of Function | Variant category |
| pLI | Probability of LoF Intolerance | Constraint metric |
| LOEUF | LoF Observed/Expected Upper Fraction | Constraint metric |
| O/E | Observed/Expected ratio | Constraint metric |
| VEP | Variant Effect Predictor | Annotation tool |
| QC | Quality Control | Data filtering |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant identifier |
| ClinVar | VCV | Clinical significance |
| Ensembl | ENSG/ENST | Gene/transcript |
| OMIM | MIM | Disease associations |
| dbNSFP | Variant | Prediction scores |

---

## Data Subsets

| Dataset | Description | Size |
|---------|-------------|------|
| Exomes | Exome sequencing | 730,000 |
| Genomes | Whole genomes | 76,000 |
| Joint | Combined call set | 807,000 |
| v2 | GRCh37 legacy | 141,000 |
| v3 | GRCh38 genomes | 76,000 |
| v4 | GRCh38 exomes+genomes | 807,000 |

---

## Data Quality Notes

1. **Cardinality:** One aggregate frequency per variant per population
2. **Exclusions:** Pediatric disease cohorts excluded
3. **Healthy Controls:** Generally healthy adults only
4. **Ancestry Bias:** v3 ~60% European; v4 more diverse
5. **Constraint Metrics:** Available for protein-coding genes only
6. **Coverage:** Exomes have higher depth than genomes in coding regions
