# ImmunoBase - Data Dictionary

## Overview

This data dictionary documents the schema for ImmunoBase immune disease genetics database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | immunobase |
| **Name** | ImmunoBase |
| **Parent** | 3.6.autoimmune.inflammatory |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Disease Record

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| disease_id | string | 1:1 | Yes | Disease identifier | T1D |
| disease_name | string | 1:1 | Yes | Full disease name | Type 1 Diabetes |
| associated_regions | integer | 1:1 | No | Number of associated regions | 60 |
| gwas_snps | integer | 1:1 | No | GWAS-significant SNPs | 150 |
| credible_sets | integer | 1:1 | No | Fine-mapped credible sets | 45 |

### Association Record

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| snp_id | string | 1:1 | Yes | dbSNP identifier | rs2476601 |
| chromosome | string | 1:1 | Yes | Chromosome | 1 |
| position | integer | 1:1 | Yes | Genomic position | 114377568 |
| disease | string | 1:1 | Yes | Disease code | T1D |
| p_value | float | 1:1 | Yes | Association p-value | 1.2e-45 |
| odds_ratio | float | 1:1 | No | Effect size | 1.89 |
| risk_allele | string | 1:1 | No | Risk allele | A |
| gene | string | 1:1 | No | Nearest gene | PTPN22 |
| source | string | 1:1 | Yes | Data source | GWAS, ImmunoChip |

### Region Record

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| region_id | string | 1:1 | Yes | Region identifier | 1p13.2 |
| chromosome | string | 1:1 | Yes | Chromosome | 1 |
| start | integer | 1:1 | Yes | Start position | 114300000 |
| end | integer | 1:1 | Yes | End position | 114500000 |
| genes | array | 1:N | No | Genes in region | [PTPN22, RSBN1] |
| diseases | array | 1:N | Yes | Associated diseases | [T1D, RA, CEL] |
| lead_snp | string | 1:1 | No | Lead SNP | rs2476601 |
| shared | boolean | 1:1 | No | Shared across diseases | true |

### Credible Set

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| set_id | string | 1:1 | Yes | Credible set ID | CS_1p13.2_T1D |
| region | string | 1:1 | Yes | Region ID | 1p13.2 |
| disease | string | 1:1 | Yes | Disease code | T1D |
| variants | array | 1:N | Yes | Variants in set | [rs2476601, rs6679677] |
| posterior_probs | array | 1:N | No | Posterior probabilities | [0.85, 0.10] |
| credible_interval | float | 1:1 | No | Coverage (e.g., 95%) | 0.95 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| dbSNP ID | rs[0-9]+ | rs2476601 | Variant identifier |
| Region ID | [0-9]+[pq][0-9.]+ | 1p13.2 | Cytogenetic band |
| Gene Symbol | Text | PTPN22 | HGNC symbol |
| Disease Code | Text | T1D | Disease abbreviation |
| Credible Set | CS_[region]_[disease] | CS_1p13.2_T1D | Set identifier |

---

## Enumerations

### Diseases Covered

| Code | Name | GWAS Loci |
|------|------|-----------|
| T1D | Type 1 Diabetes | 60+ |
| RA | Rheumatoid Arthritis | 100+ |
| CEL | Celiac Disease | 40+ |
| MS | Multiple Sclerosis | 200+ |
| CD | Crohn's Disease | 140+ |
| UC | Ulcerative Colitis | 100+ |
| JIA | Juvenile Idiopathic Arthritis | 20+ |
| AS | Ankylosing Spondylitis | 40+ |
| PSO | Psoriasis | 60+ |
| PSC | Primary Sclerosing Cholangitis | 15+ |
| PBC | Primary Biliary Cholangitis | 30+ |
| VIT | Vitiligo | 50+ |

### Data Sources

| Source | Description |
|--------|-------------|
| GWAS | Genome-wide association studies |
| ImmunoChip | Custom genotyping array |
| Fine-mapping | Statistical fine-mapping |
| Meta-analysis | Combined analysis |
| Functional | Functional annotation |

### Association Categories

| Category | Description |
|----------|-------------|
| Genome-wide significant | P < 5e-8 |
| Suggestive | 5e-8 < P < 1e-5 |
| Fine-mapped | Credible set analysis |
| Replicated | Independent replication |

### Annotation Types

| Type | Description |
|------|-------------|
| Coding | Coding variants |
| Regulatory | Regulatory regions |
| eQTL | Expression QTL |
| pQTL | Protein QTL |
| Chromatin | Chromatin accessibility |

### Effect Directions

| Direction | Description |
|-----------|-------------|
| Risk | Increases disease risk |
| Protective | Decreases disease risk |
| Unknown | Direction not determined |

---

## Entity Relationships

### Disease to Associations
- **Cardinality:** 1:N
- **Description:** Diseases have multiple SNP associations
- **Key Fields:** disease_id, snp_id

### Region to Diseases
- **Cardinality:** N:M
- **Description:** Regions shared across diseases
- **Key Fields:** region_id, disease_id

### Credible Set to Variants
- **Cardinality:** 1:N
- **Description:** Credible sets contain multiple variants
- **Key Fields:** set_id, variants

### Region to Genes
- **Cardinality:** 1:N
- **Description:** Regions contain multiple genes
- **Key Fields:** region_id, genes

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GWAS | Genome-Wide Association Study | Study type |
| SNP | Single Nucleotide Polymorphism | Variant type |
| LD | Linkage Disequilibrium | Genetic correlation |
| OR | Odds Ratio | Effect size |
| CI | Credible Interval | Statistical interval |
| eQTL | Expression Quantitative Trait Locus | Functional annotation |
| T1D | Type 1 Diabetes | Disease |
| RA | Rheumatoid Arthritis | Disease |
| MS | Multiple Sclerosis | Disease |
| CD | Crohn's Disease | Disease |
| UC | Ulcerative Colitis | Disease |
| IBD | Inflammatory Bowel Disease | Disease group |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| dbSNP | rsID | Variant data |
| GWAS Catalog | Study ID | GWAS studies |
| Ensembl | Gene ID | Gene annotation |
| GTEx | eQTL ID | Expression data |
| Open Targets | Disease ID | Target-disease |
| PubMed | PMID | Literature |

---

## Data Quality Notes

1. **Disease Coverage:** 12 major immune diseases
2. **Associated Regions:** 350+ genetic regions
3. **SNP Associations:** 1,000+ GWAS variants
4. **Credible Sets:** 200+ fine-mapped sets
5. **Shared Loci:** Cross-disease sharing analysis
6. **ImmunoChip:** Custom array with 200K variants
7. **Expert Curation:** Curated associations
8. **Regular Updates:** Incorporated new GWAS findings

