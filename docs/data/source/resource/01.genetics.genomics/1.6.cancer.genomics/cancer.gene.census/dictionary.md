# Cancer Gene Census - Data Dictionary

## Overview

This data dictionary documents the schema for COSMIC Cancer Gene Census.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | cancer.gene.census |
| **Name** | Cancer Gene Census |
| **Parent** | 1.6.cancer.genomics |
| **Total Fields** | 18 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Gene Entry

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Gene_Symbol | string | 1:1 | Yes | HGNC gene symbol | TP53 |
| Name | string | 1:1 | Yes | Full gene name | tumor protein p53 |
| Entrez_GeneId | integer | 1:1 | Yes | NCBI Gene ID | 7157 |
| Genome_Location | string | 1:1 | Yes | Cytogenetic band | 17p13.1 |
| Tier | integer | 1:1 | Yes | Evidence tier (1 or 2) | 1 |
| Hallmark | string | 1:N | No | Cancer hallmarks | genome instability |
| Chr_Band | string | 1:1 | No | Chromosome band | 17p13 |
| Somatic | string | 1:1 | No | Somatic mutation evidence | yes |
| Germline | string | 1:1 | No | Germline mutation evidence | yes |
| Tumour_Types_Somatic | string | 1:N | No | Cancer types (somatic) | breast, lung, colorectal |
| Tumour_Types_Germline | string | 1:N | No | Cancer types (germline) | Li-Fraumeni syndrome |
| Cancer_Syndrome | string | 1:N | No | Associated syndromes | Li-Fraumeni syndrome |

### Functional Annotation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| Role_in_Cancer | string | 1:1 | Yes | Oncogene/TSG/fusion | TSG |
| Mutation_Types | string | 1:N | No | Mutation categories | Mis, N, F |
| Translocation_Partner | string | 1:N | No | Fusion partners | BCR, ABL1 |
| Molecular_Genetics | string | 1:1 | No | Dominant/recessive | Rec |
| Other_Germline_Mut | string | 1:1 | No | Other germline mutations | yes |
| Other_Syndrome | string | 1:N | No | Other syndromes | HNPCC |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Gene Symbol | HGNC symbol | TP53 | Gene identifier |
| Entrez Gene ID | Integer | 7157 | NCBI Gene ID |
| Ensembl ID | ENSG + digits | ENSG00000141510 | Ensembl gene |
| COSMIC Gene ID | Integer | 376 | COSMIC internal ID |

---

## Enumerations

### Tier Classification

| Tier | Description | Criteria |
|------|-------------|----------|
| 1 | Strong evidence | Documented in cancer, extensive literature |
| 2 | Moderate evidence | Strong indication, less literature |

### Role in Cancer

| Role | Description |
|------|-------------|
| oncogene | Promotes cancer when activated |
| TSG | Tumor suppressor gene |
| fusion | Involved in oncogenic fusions |
| oncogene, fusion | Both roles |
| TSG, fusion | Both roles |

### Mutation Types

| Code | Description |
|------|-------------|
| Mis | Missense |
| N | Nonsense |
| F | Frameshift |
| S | Splice site |
| D | Deletion |
| A | Amplification |
| T | Translocation |
| O | Other |

### Molecular Genetics

| Code | Description |
|------|-------------|
| Dom | Dominant (heterozygous mutation sufficient) |
| Rec | Recessive (biallelic inactivation required) |

### Cancer Hallmarks

| Hallmark | Description |
|----------|-------------|
| sustaining proliferative signaling | Uncontrolled growth signals |
| evading growth suppressors | Bypass cell cycle control |
| resisting cell death | Avoid apoptosis |
| enabling replicative immortality | Unlimited division |
| inducing angiogenesis | Blood vessel formation |
| activating invasion and metastasis | Spread capability |
| genome instability and mutation | DNA repair defects |
| tumor promoting inflammation | Inflammatory environment |
| deregulating cellular energetics | Metabolic changes |
| avoiding immune destruction | Immune evasion |

---

## Entity Relationships

### Gene to Cancer Type
- **Cardinality:** 1:N
- **Description:** One gene implicated in multiple cancers
- **Key Fields:** Gene_Symbol, Tumour_Types

### Gene to Syndrome
- **Cardinality:** 1:N
- **Description:** One gene causes multiple syndromes
- **Key Fields:** Gene_Symbol, Cancer_Syndrome

### Gene to Fusion Partner
- **Cardinality:** N:M
- **Description:** Genes form fusions with multiple partners
- **Key Fields:** Gene_Symbol, Translocation_Partner

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CGC | Cancer Gene Census | Database name |
| COSMIC | Catalogue of Somatic Mutations in Cancer | Parent database |
| TSG | Tumor Suppressor Gene | Gene category |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| Dom | Dominant | Inheritance pattern |
| Rec | Recessive | Inheritance pattern |
| LOF | Loss of Function | Mutation effect |
| GOF | Gain of Function | Mutation effect |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| COSMIC | COSMIC ID | Mutation data |
| HGNC | Symbol | Gene naming |
| NCBI Gene | Entrez ID | Gene reference |
| Ensembl | ENSG | Gene annotation |
| OMIM | MIM | Syndrome reference |
| ClinVar | VCV | Germline variants |

---

## Census Statistics

| Category | Count |
|----------|-------|
| Total genes | 730+ |
| Tier 1 genes | 580+ |
| Tier 2 genes | 150+ |
| Oncogenes | 250+ |
| TSGs | 300+ |
| Fusion genes | 300+ |

---

## Data Quality Notes

1. **Cardinality:** One entry per gene
2. **Curation:** Expert-curated from literature
3. **Evidence Tiers:** Tier 1 has strongest evidence
4. **Updates:** Regularly updated with new cancer genes
5. **Access:** Requires COSMIC license for commercial use
6. **Scope:** Focus on causally implicated cancer genes
