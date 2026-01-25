# NPClassifier - Data Dictionary

## Overview

This data dictionary documents the schema for NPClassifier natural product classification service.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | npclassifier |
| **Name** | NPClassifier |
| **Parent** | 2.6.chemical.ontology.classification |
| **Total Fields** | 10+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Classification Result

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| smiles | string | 1:1 | Yes | Input SMILES string | CC1=CC2=C(C=C1)NC=C2CCN |
| pathway_results | array | 1:N | Yes | Predicted pathways with scores | [{pathway, pathway_score}] |
| superclass_results | array | 1:N | Yes | Predicted superclasses | [{superclass, superclass_score}] |
| class_results | array | 1:N | Yes | Predicted classes | [{class, class_score}] |
| isglycoside | boolean | 1:1 | Yes | Glycoside moiety detected | true, false |

### Pathway Result Object

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| pathway | string | 1:1 | Yes | Pathway name | Alkaloids |
| pathway_score | decimal | 1:1 | Yes | Confidence score (0-1) | 0.9823 |

### Superclass Result Object

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| superclass | string | 1:1 | Yes | Superclass name | Indole alkaloids |
| superclass_score | decimal | 1:1 | Yes | Confidence score (0-1) | 0.9456 |

### Class Result Object

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| class | string | 1:1 | Yes | Class name | Simple tryptamines |
| class_score | decimal | 1:1 | Yes | Confidence score (0-1) | 0.8912 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| SMILES | varies | CC1=CC2=C(C=C1)NC=C2CCN | Input structure |
| InChIKey | 27 characters | APJYDQYYACXCRM-UHFFFAOYSA-N | Structure identifier |

---

## Enumerations

### Biosynthetic Pathway Classes (7 Total)

| Pathway | Description | Example Compounds |
|---------|-------------|-------------------|
| Terpenoids | Isoprene-derived (MVA/MEP pathway) | Steroids, carotenoids, monoterpenes |
| Polyketides | Acetyl-CoA derived | Macrolides, tetracyclines, anthraquinones |
| Alkaloids | Nitrogen-containing natural products | Indole alkaloids, isoquinolines |
| Amino acids and Peptides | Amino acid derived | NRPS products, cyclopeptides |
| Carbohydrates | Sugar derived | Glycosides, aminoglycosides, nucleosides |
| Fatty acids | Fatty acid derived | Lipids, prostaglandins |
| Shikimates and Phenylpropanoids | Shikimate pathway | Flavonoids, coumarins, lignans |

### Terpenoid Superclasses

| Superclass | Example Classes |
|------------|-----------------|
| Monoterpenoids | Acyclic, menthane monoterpenoids |
| Sesquiterpenoids | Eudesmane, germacrane sesquiterpenoids |
| Diterpenoids | Abietane, labdane diterpenoids |
| Triterpenoids | Oleanane, ursane triterpenoids |
| Steroids | Cholestane, androstane steroids |
| Carotenoids | C40 carotenoids, apocarotenoids |

### Alkaloid Superclasses

| Superclass | Example Classes |
|------------|-----------------|
| Indole alkaloids | Simple tryptamines, carbazole alkaloids |
| Isoquinoline alkaloids | Benzylisoquinolines, aporphines |
| Pyridine alkaloids | Simple pyridines, nicotinoids |
| Purine alkaloids | Xanthines, purines |
| Tropane alkaloids | Tropanes, ecgonine derivatives |

### Shikimates/Phenylpropanoid Superclasses

| Superclass | Example Classes |
|------------|-----------------|
| Flavonoids | Flavones, flavonols, isoflavones, anthocyanins |
| Coumarins | Simple coumarins, furanocoumarins |
| Lignans | Dibenzylbutane, arylnaphthalene lignans |
| Phenolic acids | Hydroxybenzoic, hydroxycinnamic acids |
| Stilbenoids | Monomeric, oligomeric stilbenoids |

### Confidence Score Interpretation

| Score Range | Interpretation | Recommendation |
|-------------|----------------|----------------|
| 0.9 - 1.0 | High confidence | Accept classification |
| 0.7 - 0.9 | Moderate confidence | Review structure |
| 0.5 - 0.7 | Low confidence | Manual verification needed |
| < 0.5 | Very low confidence | Likely not natural product |

### HTTP Status Codes

| Code | Meaning |
|------|---------|
| 200 | Success |
| 400 | Invalid SMILES |
| 500 | Server error |
| 503 | Service unavailable |

---

## Entity Relationships

### Input to Pathway Classification
- **Cardinality:** 1:N
- **Description:** One compound may have multiple pathway predictions
- **Key Fields:** smiles, pathway_results

### Pathway to Superclass
- **Cardinality:** 1:N
- **Description:** Each pathway contains multiple superclasses
- **Key Fields:** pathway, superclass

### Superclass to Class
- **Cardinality:** 1:N
- **Description:** Each superclass contains multiple classes
- **Key Fields:** superclass, class

### Compound to Glycoside Status
- **Cardinality:** 1:1
- **Description:** Binary glycoside detection
- **Key Fields:** smiles, isglycoside

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NP | Natural Product | Compound of biological origin |
| MVA | Mevalonate Pathway | Terpenoid biosynthesis |
| MEP | Methylerythritol Phosphate Pathway | Alternative terpenoid pathway |
| NRPS | Non-Ribosomal Peptide Synthetase | Peptide biosynthesis |
| BGC | Biosynthetic Gene Cluster | Genes for NP production |
| GNPS | Global Natural Products Social Molecular Networking | Integration platform |
| SMILES | Simplified Molecular-Input Line-Entry System | Structure notation |
| API | Application Programming Interface | Access method |
| REST | Representational State Transfer | API architecture |
| JSON | JavaScript Object Notation | Response format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| COCONUT | NP ID | Natural products database |
| LOTUS | Lotus ID | NP occurrence database |
| ClassyFire | ChemOnt ID | General chemical taxonomy |
| ChEBI | ChEBI ID | Chemical ontology |
| GNPS | Spectrum ID | Mass spectrometry |

---

## Data Quality Notes

1. **Deep Learning Model:** Neural network trained on 80,000+ compounds
2. **Accuracy:** >90% cross-validation accuracy
3. **Inference Speed:** <100ms per compound average
4. **Biosynthetic Focus:** Classification based on predicted biosynthetic origin
5. **Multi-Pathway:** Returns multiple predictions for hybrid natural products
6. **Glycoside Detection:** Automatic sugar moiety identification
7. **Open Access:** Free web interface and API
8. **GNPS Integration:** Part of GNPS molecular networking ecosystem

