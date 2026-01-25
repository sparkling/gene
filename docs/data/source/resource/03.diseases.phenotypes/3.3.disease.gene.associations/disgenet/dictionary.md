# DisGeNET - Data Dictionary

## Overview

This data dictionary documents the schema for DisGeNET gene-disease association database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | disgenet |
| **Name** | DisGeNET |
| **Parent** | 3.3.disease.gene.associations |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Gene-Disease Association (GDA)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| geneId | integer | 1:1 | Yes | NCBI Gene ID | 3157 |
| geneSymbol | string | 1:1 | Yes | HGNC gene symbol | HMGCS1 |
| diseaseId | string | 1:1 | Yes | UMLS CUI | C0020443 |
| diseaseName | string | 1:1 | Yes | Disease name | Hypercholesterolemia |
| score | float | 1:1 | Yes | GDA score (0-1) | 0.72 |
| EI | float | 1:1 | No | Evidence Index (0-1) | 1.0 |
| associationType | string | 1:1 | No | Association ontology | GeneticVariation |
| source | string | 1:1 | Yes | Data source | CTD_human |
| pmid | array | 1:N | No | PubMed IDs | [12345678] |
| year | integer | 1:1 | No | First publication year | 2020 |
| NofPmids | integer | 1:1 | No | Publication count | 15 |
| NofSnps | integer | 1:1 | No | Associated SNP count | 3 |

### Variant-Disease Association (VDA)

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| snpId | string | 1:1 | Yes | dbSNP ID | rs2066843 |
| chromosome | string | 1:1 | Yes | Chromosome | 1 |
| position | integer | 1:1 | Yes | GRCh38 position | 67355291 |
| diseaseId | string | 1:1 | Yes | UMLS CUI | C0027051 |
| diseaseName | string | 1:1 | Yes | Disease name | Myocardial Infarction |
| score | float | 1:1 | Yes | VDA score (0-1) | 0.58 |
| EI | float | 1:1 | No | Evidence Index | 0.85 |
| source | string | 1:1 | Yes | Data source | GWAS Catalog |
| pmid | array | 1:N | No | PubMed IDs | [23456789] |
| consequence | string | 1:1 | No | VEP consequence | missense_variant |

### Gene Attributes

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| geneId | integer | 1:1 | Yes | NCBI Gene ID | 3157 |
| geneSymbol | string | 1:1 | Yes | HGNC symbol | HMGCS1 |
| geneName | string | 1:1 | Yes | Full gene name | HMG-CoA synthase 1 |
| uniprotId | array | 1:N | No | UniProt accessions | [Q01581] |
| proteinClass | string | 1:1 | No | Protein classification | Enzyme |
| pLI | float | 1:1 | No | LoF intolerance (0-1) | 0.85 |
| DSI | float | 1:1 | No | Disease Specificity Index | 0.45 |
| DPI | float | 1:1 | No | Disease Pleiotropy Index | 0.62 |
| NofDiseases | integer | 1:1 | No | Associated diseases | 12 |

### Disease Attributes

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| diseaseId | string | 1:1 | Yes | UMLS CUI | C0020443 |
| diseaseName | string | 1:1 | Yes | Preferred name | Hypercholesterolemia |
| diseaseType | string | 1:1 | No | Disease type | disease |
| diseaseClass | array | 1:N | No | MeSH disease classes | [Metabolic] |
| diseaseSemanticType | string | 1:1 | No | UMLS semantic type | Disease or Syndrome |
| DOID | string | 1:1 | No | Disease Ontology ID | DOID:1287 |
| HPO | array | 1:N | No | HPO phenotype terms | [HP:0003124] |
| NofGenes | integer | 1:1 | No | Associated gene count | 45 |
| NofVariants | integer | 1:1 | No | Associated variant count | 128 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| NCBI Gene ID | Integer | 3157 | Gene identifier |
| UMLS CUI | C[0-9]{7} | C0020443 | Disease identifier |
| dbSNP ID | rs[0-9]+ | rs2066843 | Variant identifier |
| HGNC ID | HGNC:[0-9]+ | HGNC:5007 | Gene nomenclature |
| Ensembl ID | ENSG[0-9]{11} | ENSG00000112972 | Gene annotation |
| UniProt ID | [A-Z0-9]{6,10} | Q01581 | Protein identifier |
| DOID | DOID:[0-9]+ | DOID:1287 | Disease ontology |
| HPO | HP:[0-9]{7} | HP:0003124 | Phenotype ontology |

---

## Enumerations

### Score Interpretation

| Range | Confidence | Description |
|-------|------------|-------------|
| > 0.7 | High | Strong evidence |
| 0.4 - 0.7 | Medium | Moderate evidence |
| 0.1 - 0.4 | Low | Weak evidence |
| < 0.1 | Very low | Predicted/inferred |

### Evidence Index (EI)

| Value | Meaning |
|-------|---------|
| 1.0 | All publications support |
| 0.5 - 1.0 | Majority support |
| < 0.5 | Significant contradiction |

### Association Types

| Type | Description |
|------|-------------|
| Therapeutic | Drug target relationship |
| Biomarker | Diagnostic/prognostic marker |
| GeneticVariation | Genetic association |
| Causal | Disease-causing mutation |
| SomaticMutation | Cancer somatic variant |
| GermlineMutation | Inherited variant |
| Susceptibility | Risk variant |
| Modifying | Disease modifier |
| AlteredExpression | Expression changes |

### Source Types

| Source | Type | Description |
|--------|------|-------------|
| UniProt | Curated | Manual curation |
| CTD_human | Curated | Comparative toxicogenomics |
| CGI | Curated | Cancer Genome Interpreter |
| ClinGen | Curated | Clinical Genome Resource |
| PanelApp | Curated | Genomics England |
| Orphanet | Curated | Rare diseases |
| GWAS Catalog | Curated | GWAS associations |
| ClinVar | Curated | Clinical variants |
| RGD | Animal | Rat Genome Database |
| MGD | Animal | Mouse Genome Database |
| BeFree | Literature | Text mining |
| LHGDN | Literature | Literature mining |

### Variant Consequences (VEP)

| Consequence | Impact | Percentage |
|-------------|--------|------------|
| missense_variant | Moderate | 28% |
| intron_variant | Modifier | 26% |
| frameshift_variant | High | 11% |
| intergenic_variant | Modifier | 11% |
| synonymous_variant | Low | 5% |
| stop_gained | High | 3% |
| splice_region_variant | Low-Moderate | 4% |

### Disease Semantic Types

| Type | Description |
|------|-------------|
| Disease or Syndrome | Clinical diseases |
| Neoplastic Process | Cancers |
| Mental or Behavioral Dysfunction | Psychiatric |
| Pathologic Function | Pathology |
| Sign or Symptom | Clinical manifestations |
| Congenital Abnormality | Birth defects |

---

## Entity Relationships

### Gene to Diseases (GDA)
- **Cardinality:** N:M
- **Description:** Genes associated with multiple diseases
- **Key Fields:** geneId, diseaseId

### Variant to Diseases (VDA)
- **Cardinality:** N:M
- **Description:** Variants associated with diseases
- **Key Fields:** snpId, diseaseId

### Gene to Variants
- **Cardinality:** 1:N
- **Description:** Genes have multiple variants
- **Key Fields:** geneId, snpId

### Disease to References
- **Cardinality:** 1:N
- **Description:** Literature support
- **Key Fields:** diseaseId, pmid

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GDA | Gene-Disease Association | Core data type |
| VDA | Variant-Disease Association | Core data type |
| EI | Evidence Index | Consistency metric |
| DSI | Disease Specificity Index | Specificity metric |
| DPI | Disease Pleiotropy Index | Pleiotropy metric |
| UMLS | Unified Medical Language System | Disease ontology |
| CUI | Concept Unique Identifier | UMLS ID format |
| VEP | Variant Effect Predictor | Consequence tool |
| pLI | Probability of LoF Intolerance | Gene constraint |
| GWAS | Genome-Wide Association Study | Study type |
| CTD | Comparative Toxicogenomics Database | Data source |
| CGI | Cancer Genome Interpreter | Data source |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| NCBI Gene | Gene ID | Gene information |
| UniProt | Accession | Protein data |
| dbSNP | rsID | Variant data |
| ClinVar | Variation ID | Clinical variants |
| GWAS Catalog | Study ID | GWAS data |
| HPO | HP ID | Phenotype terms |
| DOID | DOID | Disease ontology |
| MeSH | MeSH ID | Disease classification |
| PubMed | PMID | Literature |

---

## Data Quality Notes

1. **GDA Coverage:** 628,685+ gene-disease associations
2. **VDA Coverage:** 210,498+ variant-disease associations
3. **Gene Coverage:** 17,549 genes
4. **Disease Coverage:** 24,166 diseases/phenotypes
5. **Score System:** 0-1 confidence scoring
6. **Multi-Source:** 40+ integrated sources
7. **API Access:** REST API available
8. **CC BY-NC-SA:** Academic use; commercial license available

