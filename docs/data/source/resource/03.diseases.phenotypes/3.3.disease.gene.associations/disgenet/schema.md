---
id: schema-disgenet
title: "DisGeNET Database Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./README.md)

# DisGeNET Database Schema Documentation

**Source:** https://www.disgenet.org / https://disgenet.com
**Version:** 7.0 (as of January 2026)
**License:** CC BY-NC-SA 4.0 (Academic) / Commercial license available
**Reference:** Nucleic Acids Research, 2020 Jan 8;48(D1):D845-D855

---

## TL;DR

DisGeNET is a comprehensive knowledge platform integrating information on human disease-associated genes and variants. The core data model centers on Gene-Disease Associations (GDAs) and Variant-Disease Associations (VDAs).

---

## Database Statistics (v7.0)

| Entity Type | Count |
|-------------|-------|
| Gene-Disease Associations (GDAs) | 628,685+ |
| Variant-Disease Associations (VDAs) | 210,498+ |
| Genes | 17,549 |
| Variants | 117,337 |
| Diseases/Traits/Phenotypes | 24,166 |

### GDA Source Distribution

| Source Type | Associations |
|-------------|--------------|
| Curated | 81,746 |
| Animal Models | 11,517 |
| Inferred | 163,626 |
| Literature | 415,583 |

### VDA Source Distribution

| Source Type | Associations |
|-------------|--------------|
| Curated | 165,354 |
| Literature | 48,998 |

---

## Core Tables

### 1. gene_disease_association (GDA)

| Column | Type | Description |
|--------|------|-------------|
| `geneId` | Integer | NCBI Gene ID |
| `geneSymbol` | String | HGNC gene symbol |
| `diseaseId` | String | UMLS CUI (C#######) |
| `diseaseName` | String | Disease name |
| `score` | Float | DisGeNET GDA score (0-1) |
| `EI` | Float | Evidence Index (0-1) |
| `associationType` | String | Association type ontology |
| `source` | String | Data source |
| `pmid` | Integer[] | PubMed IDs |
| `year` | Integer | First publication year |
| `NofPmids` | Integer | Number of publications |
| `NofSnps` | Integer | Associated SNPs count |

### 2. variant_disease_association (VDA)

| Column | Type | Description |
|--------|------|-------------|
| `snpId` | String | dbSNP ID (rs#######) |
| `chromosome` | String | Chromosome (1-22, X, Y, MT) |
| `position` | Integer | GRCh38 position |
| `diseaseId` | String | UMLS CUI |
| `diseaseName` | String | Disease name |
| `score` | Float | DisGeNET VDA score (0-1) |
| `EI` | Float | Evidence Index |
| `source` | String | Data source |
| `pmid` | Integer[] | PubMed IDs |
| `associationType` | String | Association type |

### 3. gene

| Column | Type | Description |
|--------|------|-------------|
| `geneId` | Integer | NCBI Gene ID |
| `geneSymbol` | String | HGNC symbol |
| `geneName` | String | Full gene name |
| `uniprotId` | String[] | UniProt accessions |
| `proteinClass` | String | Protein classification |
| `pLI` | Float | Loss-of-function intolerance (0-1) |
| `DSI` | Float | Disease Specificity Index |
| `DPI` | Float | Disease Pleiotropy Index |
| `NofDiseases` | Integer | Number of associated diseases |

### 4. variant

| Column | Type | Description |
|--------|------|-------------|
| `snpId` | String | dbSNP ID |
| `chromosome` | String | Chromosome |
| `position` | Integer | GRCh38 position |
| `reference` | String | Reference allele |
| `alternative` | String | Alternative allele |
| `consequence` | String | VEP consequence |
| `geneSymbol` | String | Associated gene |
| `gnomAD_AF` | Float | gnomAD allele frequency |
| `DSI` | Float | Disease Specificity Index |
| `DPI` | Float | Disease Pleiotropy Index |

### 5. disease

| Column | Type | Description |
|--------|------|-------------|
| `diseaseId` | String | UMLS CUI |
| `diseaseName` | String | Preferred name |
| `diseaseType` | String | Disease, phenotype, group |
| `diseaseClass` | String[] | MeSH disease classes |
| `diseaseSemanticType` | String | UMLS semantic type |
| `DOID` | String | Human Disease Ontology ID |
| `HPO` | String[] | HPO phenotype terms |
| `NofGenes` | Integer | Associated genes count |
| `NofVariants` | Integer | Associated variants count |

---

## Score Metrics

### DisGeNET Score (GDA/VDA Score)

**Range:** 0.0 - 1.0

**Calculation factors:**
1. Number of sources supporting the association
2. Reliability weight of each source type:
   - Expert curated sources: Highest weight
   - Animal model evidence: Medium-high weight
   - Literature mining: Medium weight
   - Inferred associations: Lower weight
3. Number of supporting publications
4. Consistency across sources

**Interpretation:**
| Score Range | Confidence |
|-------------|------------|
| > 0.7 | High confidence |
| 0.4 - 0.7 | Medium confidence |
| 0.1 - 0.4 | Low confidence |
| < 0.1 | Very low / predicted |

### Evidence Index (EI)

**Range:** 0.0 - 1.0

**Definition:** Reflects consistency of evidence across publications.

**Calculation:**
```
EI = (Supporting publications) / (Total publications)
```

**Interpretation:**
| EI Value | Meaning |
|----------|---------|
| 1.0 | All publications support the association |
| 0.5 - 1.0 | Majority support, some contradictory |
| < 0.5 | Significant contradictory evidence |

### Disease Specificity Index (DSI)

**Range:** 0.0 - 1.0

**Definition:** Inversely proportional to number of diseases associated with a gene/variant.

**Formula:**
```
DSI = log2(N_diseases) / log2(N_total_diseases)
```
Where:
- N_diseases = diseases associated with the gene
- N_total_diseases = total diseases in DisGeNET

**Interpretation:**
| DSI Value | Meaning |
|-----------|---------|
| ~1.0 | Gene associated with single disease (highly specific) |
| ~0.5 | Gene associated with moderate number of diseases |
| ~0.0 | Gene associated with many diseases (e.g., TNF > 1500 diseases) |

### Disease Pleiotropy Index (DPI)

**Range:** 0.0 - 1.0

**Definition:** Measures whether associated diseases span multiple disease classes.

**Formula:**
```
DPI = 1 - (1 / N_disease_classes)
```

**Interpretation:**
| DPI Value | Meaning |
|-----------|---------|
| ~0.0 | All associated diseases in same MeSH class |
| ~0.5 | Diseases span ~2 classes |
| ~1.0 | Diseases span many different classes (pleiotropic) |

---

## Association Types (Ontology)

### Therapeutic Association Types

| Code | Type | Description |
|------|------|-------------|
| `Therapeutic` | Therapeutic relationship | Drug target for disease |
| `DrugTarget` | Drug target | Gene encodes drug target |

### Biomarker Association Types

| Code | Type | Description |
|------|------|-------------|
| `Biomarker` | Disease biomarker | Diagnostic/prognostic marker |
| `AlteredExpression` | Expression altered | Gene expression changes |
| `Diagnostic` | Diagnostic marker | Used in diagnosis |
| `Prognostic` | Prognostic marker | Predicts outcome |

### Genetic Association Types

| Code | Type | Description |
|------|------|-------------|
| `GeneticVariation` | Genetic association | Variant associated with disease |
| `Causal` | Causal mutation | Disease-causing mutation |
| `SomaticMutation` | Somatic mutation | Cancer somatic variant |
| `GermlineMutation` | Germline mutation | Inherited variant |
| `Susceptibility` | Susceptibility | Risk variant |
| `Modifying` | Modifier | Disease course modifier |

### Functional Association Types

| Code | Type | Description |
|------|------|-------------|
| `FunctionalAnnotation` | Functional role | Pathway involvement |
| `Posttranslational` | PTM association | Post-translational modification |

---

## Data Sources

### Curated Sources

| Source | Description | Association Types |
|--------|-------------|-------------------|
| UniProt | Manual curation | GDA |
| CTD (human) | Comparative Toxicogenomics | GDA |
| CGI | Cancer Genome Interpreter | GDA, VDA |
| ClinGen | Clinical Genome Resource | GDA |
| Genomics England PanelApp | Gene panels | GDA |
| Orphanet | Rare diseases | GDA |
| GWAS Catalog | GWAS associations | VDA |
| ClinVar | Clinical variants | VDA |

### Animal Model Sources

| Source | Description |
|--------|-------------|
| RGD | Rat Genome Database |
| MGD | Mouse Genome Database |
| CTD (mouse/rat) | Animal model toxicogenomics |

### Literature Sources

| Source | Description |
|--------|-------------|
| LHGDN | Literature Human Gene Disease Network |
| BeFree | Text mining pipeline |
| GAD | Genetic Association Database |

### Inferred Sources

| Source | Description |
|--------|-------------|
| Inferred from HPO | Phenotype-based inference |

---

## Variant Consequence Types (VEP)

| Consequence | Percentage | Impact |
|-------------|------------|--------|
| missense_variant | 28% | Moderate |
| intron_variant | 26% | Modifier |
| frameshift_variant | 11% | High |
| intergenic_variant | 11% | Modifier |
| synonymous_variant | 5% | Low |
| splice_region_variant | 4% | Low-Moderate |
| stop_gained | 3% | High |
| 3_prime_UTR_variant | 3% | Modifier |
| 5_prime_UTR_variant | 2% | Modifier |
| Other | 7% | Various |

---

## API Endpoints

**Base URL:** `https://www.disgenet.org/api`

### GDA Endpoints

| Endpoint | Description |
|----------|-------------|
| `GET /gda/gene/{gene_id}` | GDAs for a gene |
| `GET /gda/disease/{disease_id}` | GDAs for a disease |
| `GET /gda/source/{source}` | GDAs by source |
| `GET /gda/evidences/{type}` | GDAs by evidence type |

### VDA Endpoints

| Endpoint | Description |
|----------|-------------|
| `GET /vda/variant/{variant_id}` | VDAs for a variant |
| `GET /vda/disease/{disease_id}` | VDAs for a disease |
| `GET /vda/gene/{gene_id}` | VDAs for a gene |

### Attribute Endpoints

| Endpoint | Description |
|----------|-------------|
| `GET /gene/{gene_id}` | Gene attributes |
| `GET /variant/{variant_id}` | Variant attributes |
| `GET /disease/{disease_id}` | Disease attributes |

### Query Parameters

| Parameter | Description |
|-----------|-------------|
| `source` | Filter by source |
| `min_score` | Minimum score threshold |
| `min_ei` | Minimum evidence index |
| `format` | Response format (json/tsv/xml) |
| `limit` | Results limit |
| `offset` | Pagination offset |

---

## Download Files

### Tab-Separated Files

| File | Description |
|------|-------------|
| `curated_gene_disease_associations.tsv` | Expert-curated GDAs |
| `befree_gene_disease_associations.tsv` | Literature-mined GDAs |
| `all_gene_disease_associations.tsv` | Complete GDA dataset |
| `curated_variant_disease_associations.tsv` | Curated VDAs |
| `all_variant_disease_associations.tsv` | Complete VDA dataset |
| `disease_mappings.tsv` | Disease ID mappings |

### SQLite Database

Available for local deployment with full schema.

---

## Disease Semantic Types (UMLS)

| Semantic Type | Description |
|---------------|-------------|
| Disease or Syndrome | Clinical diseases |
| Neoplastic Process | Cancers |
| Mental or Behavioral Dysfunction | Psychiatric conditions |
| Pathologic Function | Pathological processes |
| Sign or Symptom | Clinical manifestations |
| Finding | Clinical findings |
| Congenital Abnormality | Birth defects |

---

## Cross-References

### Gene Identifiers

| Database | Format |
|----------|--------|
| NCBI Gene | Integer |
| HGNC | HGNC:##### |
| Ensembl | ENSG########### |
| UniProt | [A-Z0-9]{6,10} |

### Disease Identifiers

| Database | Format |
|----------|--------|
| UMLS CUI | C####### |
| MeSH | D##### |
| DOID | DOID:#### |
| HPO | HP:####### |
| Orphanet | ORPHA:##### |
| ICD-10 | [A-Z]##.# |
| EFO | EFO_####### |

### Variant Identifiers

| Database | Format |
|----------|--------|
| dbSNP | rs####### |
| ClinVar | RCV########### |
| COSMIC | COSM####### |

---

## Sample Data

### GDA Record

```json
{
  "geneId": 3157,
  "geneSymbol": "HMGCS1",
  "geneName": "3-hydroxy-3-methylglutaryl-CoA synthase 1",
  "diseaseId": "C0020443",
  "diseaseName": "Hypercholesterolemia",
  "score": 0.72,
  "EI": 1.0,
  "DSI": 0.45,
  "DPI": 0.62,
  "associationType": "GeneticVariation",
  "source": "CTD_human",
  "NofPmids": 15,
  "pmids": [12345678, 23456789, 34567890]
}
```

### VDA Record

```json
{
  "snpId": "rs2066843",
  "chromosome": "1",
  "position": 67355291,
  "reference": "G",
  "alternative": "C",
  "diseaseId": "C0027051",
  "diseaseName": "Myocardial Infarction",
  "score": 0.58,
  "EI": 0.85,
  "consequence": "missense_variant",
  "geneSymbol": "IL23R",
  "source": "GWAS Catalog",
  "NofPmids": 8
}
```

---

## Integration Notes

### Recommended ID Mappings

1. **Genes:** Use NCBI Gene ID as primary, map to Ensembl/UniProt
2. **Diseases:** Use UMLS CUI as primary, map to DOID/HPO
3. **Variants:** Use dbSNP rsID as primary, include coordinates

### Score Filtering Guidelines

| Use Case | Recommended Filters |
|----------|---------------------|
| High-confidence associations | score > 0.7, EI > 0.9 |
| Curated only | source IN (UniProt, ClinVar, CGI) |
| Drug targets | associationType = 'Therapeutic' |
| Diagnostic markers | associationType = 'Biomarker' |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | TSV (tab-separated values) |
| Alternative | JSON, XML (API), SQLite |
| Compression | gzip (.gz) |
| Encoding | UTF-8 |
| API Response | JSON, TSV, XML |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "GDA1234" |
| `name` | string | Entity name | "Gene Disease Association" |
| `type` | string | Record type | "association" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 500,000+ |
| Storage | Unknown |
| Last updated | January 2026 |

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| DisGeNET | CC BY-NC-SA 4.0 | No (requires separate license) |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `GDA` | Gene-Disease Association linking a gene to a disease with evidence | HMGCS1 - Hypercholesterolemia |
| `VDA` | Variant-Disease Association linking a genetic variant to a disease | rs2066843 - Myocardial Infarction |
| `score` | DisGeNET confidence score (0-1) based on sources, publications, and consistency | 0.72 (high confidence) |
| `EI` | Evidence Index measuring consistency of evidence across publications (0-1) | 1.0 (all support) |
| `DSI` | Disease Specificity Index, inverse of number of associated diseases (0-1) | 0.9 (disease-specific gene) |
| `DPI` | Disease Pleiotropy Index measuring spread across disease classes (0-1) | 0.8 (pleiotropic gene) |
| `diseaseId` | UMLS Concept Unique Identifier for diseases | C0020443 |
| `geneId` | NCBI Gene identifier for genes | 3157 |
| `snpId` | dbSNP rsID for genetic variants | rs2066843 |
| `associationType` | Ontology term describing the nature of the association | GeneticVariation, Biomarker |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Gene-Disease Association | Relationship between a gene and disease phenotype with supporting evidence | GDA records |
| Variant-Disease Association | Relationship between genetic variant and disease with evidence | VDA records |
| Evidence Index | Measure of evidence consistency (supporting/total publications) | EI metric |
| Disease Specificity | How specific a gene/variant is to particular diseases | DSI metric |
| Disease Pleiotropy | Whether gene affects multiple unrelated disease classes | DPI metric |
| Curated Source | Manually reviewed data from expert databases | UniProt, ClinVar |
| Literature Mining | Automated text extraction from publications | BeFree source |
| Association Type | Semantic category of gene-disease relationship | Therapeutic, Biomarker |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GDA | Gene-Disease Association | Core data type |
| VDA | Variant-Disease Association | Core data type |
| EI | Evidence Index | Consistency metric |
| DSI | Disease Specificity Index | Specificity metric |
| DPI | Disease Pleiotropy Index | Pleiotropy metric |
| UMLS | Unified Medical Language System | Disease ontology |
| CUI | Concept Unique Identifier | UMLS disease ID format |
| VEP | Variant Effect Predictor | Consequence annotation |
| CTD | Comparative Toxicogenomics Database | Data source |
| CGI | Cancer Genome Interpreter | Curated source |
| pLI | Probability of LoF Intolerance | Gene constraint |
| GWAS | Genome-Wide Association Study | VDA source |
| MeSH | Medical Subject Headings | Disease classification |
| HPO | Human Phenotype Ontology | Phenotype terms |
| DOID | Disease Ontology Identifier | Disease cross-reference |
