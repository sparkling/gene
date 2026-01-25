# OMIM - Data Dictionary

## Overview

This data dictionary documents the schema for OMIM (Online Mendelian Inheritance in Man).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | omim |
| **Name** | OMIM |
| **Parent** | 3.2.phenotype.databases |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### MIM Entry

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| mimNumber | integer | 1:1 | Yes | MIM number | 154700 |
| prefix | enum | 1:1 | Yes | Entry type prefix | # |
| preferredTitle | string | 1:1 | Yes | Official title | MARFAN SYNDROME |
| alternativeTitles | array | 1:N | No | Alternative names | [MFS] |
| geneSymbols | array | 1:N | No | Associated gene symbols | [FBN1] |
| textSections | array | 1:N | No | Content sections | Text, Clinical Features |
| created | date | 1:1 | Yes | Entry creation date | 1986-06-02 |
| updated | date | 1:1 | Yes | Last update date | 2024-01-15 |

### Clinical Synopsis

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| mimNumber | integer | 1:1 | Yes | MIM number | 154700 |
| clinicalFeatures | object | 1:1 | No | Features by system | Growth, Cardiovascular... |
| inheritance | string | 1:1 | No | Inheritance pattern | Autosomal dominant |
| molecularBasis | string | 1:1 | No | Molecular basis | FBN1 mutations |

### Allelic Variants

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| mimNumber | integer | 1:1 | Yes | Gene MIM number | 134797 |
| allelicVariantNumber | string | 1:1 | Yes | Variant suffix | .0001 |
| name | string | 1:1 | Yes | Variant name | CYS1039TYR |
| status | enum | 1:1 | No | Variant status | Pathogenic |
| mutations | array | 1:N | No | Mutation details | c.3116G>A |
| phenotypes | array | 1:N | No | Associated phenotypes | Marfan syndrome |

### Gene-Phenotype Relationships

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| geneMim | integer | 1:1 | Yes | Gene MIM number | 134797 |
| phenotypeMim | integer | 1:1 | Yes | Phenotype MIM number | 154700 |
| phenotypeMap | integer | 1:1 | Yes | Mapping key | 3 |
| inheritance | string | 1:1 | No | Inheritance pattern | AD |
| phenotypeName | string | 1:1 | Yes | Phenotype name | Marfan syndrome |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| MIM Number | [0-9]{6} | 154700 | Base identifier |
| Gene (*) | \*[0-9]{6} | *134797 | Gene with known sequence |
| Phenotype (#) | #[0-9]{6} | #154700 | Phenotype with known basis |
| Combined (+) | \+[0-9]{6} | +134610 | Gene and phenotype |
| Unknown (%) | %[0-9]{6} | %176300 | Phenotype, basis unknown |
| Moved (^) | ^[0-9]{6} | ^176400 | Entry moved/removed |
| Variant | [0-9]{6}\.[0-9]{4} | 134797.0001 | Allelic variant |

---

## Enumerations

### Entry Prefixes

| Prefix | Name | Description |
|--------|------|-------------|
| * | Asterisk | Gene with known sequence |
| # | Number sign | Phenotype, molecular basis known |
| + | Plus sign | Gene and phenotype, combined |
| % | Percent | Phenotype, basis unknown |
| None | No prefix | Other, mainly phenotypes |
| ^ | Caret | Entry moved or removed |

### Phenotype Mapping Keys

| Key | Meaning | Description |
|-----|---------|-------------|
| 1 | Disorder placed on map | Based on association |
| 2 | Disorder mapped by linkage | No mutation in gene |
| 3 | Molecular basis known | Mutation in gene |
| 4 | Chromosome deletion/duplication syndrome | Contiguous gene |

### Inheritance Patterns

| Code | Name | Description |
|------|------|-------------|
| AD | Autosomal dominant | One copy sufficient |
| AR | Autosomal recessive | Two copies required |
| XLD | X-linked dominant | X-linked, one copy |
| XLR | X-linked recessive | X-linked, males affected |
| YL | Y-linked | Y chromosome only |
| MT | Mitochondrial | Maternal inheritance |
| IC | Isolated cases | Sporadic |
| MU | Multifactorial | Complex inheritance |
| SMu | Somatic mutation | Cancer genes |
| SMo | Somatic mosaicism | Mosaic patterns |

### Clinical Synopsis Categories

| Category | Description |
|----------|-------------|
| Growth | Height, weight abnormalities |
| Head | Head/face features |
| Eyes | Ocular features |
| Ears | Hearing, ear features |
| Cardiovascular | Heart, vessel abnormalities |
| Respiratory | Lung features |
| Skeletal | Bone abnormalities |
| Muscle | Muscular features |
| Neurologic | Brain, nerve features |
| Metabolic | Metabolic abnormalities |
| Laboratory | Lab abnormalities |
| Inheritance | Inheritance pattern |
| Molecular Basis | Gene/mutation info |

### Text Section Types

| Section | Description |
|---------|-------------|
| Description | Disease description |
| Clinical Features | Clinical manifestations |
| Inheritance | Inheritance pattern |
| Mapping | Chromosomal location |
| Molecular Genetics | Gene information |
| Population Genetics | Prevalence |
| Diagnosis | Diagnostic criteria |
| Clinical Management | Treatment |
| History | Historical information |
| Animal Model | Animal studies |

---

## Entity Relationships

### Gene to Phenotypes
- **Cardinality:** 1:N
- **Description:** Genes associated with multiple phenotypes
- **Key Fields:** geneMim, phenotypeMim

### Phenotype to Allelic Variants
- **Cardinality:** N:M
- **Description:** Phenotypes caused by multiple variants
- **Key Fields:** phenotypeMim, allelicVariantNumber

### Entry to Text Sections
- **Cardinality:** 1:N
- **Description:** Entries contain multiple text sections
- **Key Fields:** mimNumber, textSections

### Entry to References
- **Cardinality:** 1:N
- **Description:** Literature citations
- **Key Fields:** mimNumber, referenceNumber

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| OMIM | Online Mendelian Inheritance in Man | Database name |
| MIM | Mendelian Inheritance in Man | Original name |
| AD | Autosomal Dominant | Inheritance pattern |
| AR | Autosomal Recessive | Inheritance pattern |
| XLR | X-Linked Recessive | Inheritance pattern |
| XLD | X-Linked Dominant | Inheritance pattern |
| MT | Mitochondrial | Inheritance pattern |
| JHU | Johns Hopkins University | Host institution |
| NCBI | National Center for Biotechnology Information | Data partner |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| PMID | PubMed Identifier | Literature reference |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| HGNC | HGNC ID | Gene nomenclature |
| Ensembl | Gene ID | Gene annotation |
| UniProt | Accession | Protein data |
| Orphanet | ORPHA ID | Rare diseases |
| HPO | HP ID | Phenotype terms |
| ClinVar | Variation ID | Clinical variants |
| MONDO | MONDO ID | Disease ontology |
| PubMed | PMID | Literature |

---

## Data Quality Notes

1. **Gene Entries:** 16,000+ genes documented
2. **Phenotype Entries:** 9,000+ phenotypes
3. **Gene-Phenotype Links:** 7,200+ relationships
4. **Allelic Variants:** 29,000+ documented variants
5. **Expert Authorship:** Expert-written entries
6. **Daily Updates:** Continuous updates
7. **API Access:** Requires registration
8. **Authoritative Source:** Gold standard for Mendelian genetics

