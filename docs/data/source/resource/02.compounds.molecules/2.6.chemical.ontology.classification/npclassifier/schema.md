---
id: schema-npclassifier
title: "NPClassifier Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, natural-products, classification, biosynthesis, machine-learning]
---

# NPClassifier Schema Documentation

**Document ID:** SCHEMA-NPCLASSIFIER
**Version:** 1.0
**Source Version:** NPClassifier v1.0 (2021+)

---

## TL;DR

NPClassifier provides deep learning-based natural product classification via REST API and web interface. It classifies compounds into biosynthetic pathway classes, superclasses, and classes based on predicted biosynthetic origin, returning JSON responses with confidence scores and glycoside detection.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Pathway Classes | 7 | Biosynthetic origins |
| Superclasses | 35+ | Structural superclasses |
| Classes | 200+ | Fine-grained classes |
| Training Compounds | 80,000+ | COCONUT/LOTUS |
| Model Accuracy | >90% | Cross-validation |
| Average Inference | <100ms | Per compound |

---

## Entity Relationship Overview

```
                    ┌─────────────────────┐
                    │   Input Molecule    │
                    │      (SMILES)       │
                    └──────────┬──────────┘
                               │
                    ┌──────────▼──────────┐
                    │   Deep Learning     │
                    │      Model          │
                    └──────────┬──────────┘
                               │
        ┌──────────────────────┼──────────────────────┐
        │                      │                      │
        ▼                      ▼                      ▼
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│   Pathway     │     │  Superclass   │     │   Glycoside   │
│   (Level 1)   │     │  (Level 2)    │     │   Detection   │
│   7 classes   │     │  35+ classes  │     │   (Yes/No)    │
└───────┬───────┘     └───────┬───────┘     └───────────────┘
        │                     │
        │                     ▼
        │             ┌───────────────┐
        │             │    Class      │
        │             │  (Level 3)    │
        │             │  200+ classes │
        │             └───────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│                  Confidence Scores                         │
│  (Pathway conf, Superclass conf, Class conf)              │
└───────────────────────────────────────────────────────────┘
```

---

## Biosynthetic Pathway Classes

### 7 Major Pathway Types

| Pathway | Description | Example Compounds |
|---------|-------------|-------------------|
| Terpenoids | Isoprene-derived (MVA/MEP) | Steroids, carotenoids, monoterpenes, sesquiterpenes |
| Polyketides | Acetyl-CoA derived | Macrolides, tetracyclines, anthraquinones |
| Alkaloids | Nitrogen-containing | Indole alkaloids, isoquinolines, pyridines |
| Amino acids and Peptides | Amino acid derived | NRPS products, cyclopeptides, beta-lactams |
| Carbohydrates | Sugar derived | Glycosides, aminoglycosides, nucleosides |
| Fatty acids | Fatty acid derived | Lipids, acyl derivatives, prostaglandins |
| Shikimates and Phenylpropanoids | Shikimate pathway | Flavonoids, coumarins, lignans, tannins |

---

## Core API Response Entity

### Classification Result Object

**Description:** Complete classification response for a compound

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| smiles | String | Yes | Input SMILES string |
| pathway_results | Array | Yes | Predicted pathways with scores |
| superclass_results | Array | Yes | Predicted superclasses with scores |
| class_results | Array | Yes | Predicted classes with scores |
| isglycoside | Boolean | Yes | Glycoside moiety detected |

### Pathway Result Object

**Description:** Single pathway prediction

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| pathway | String | Yes | Pathway name |
| pathway_score | Float | Yes | Confidence score (0-1) |

### Superclass Result Object

**Description:** Single superclass prediction

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| superclass | String | Yes | Superclass name |
| superclass_score | Float | Yes | Confidence score (0-1) |

### Class Result Object

**Description:** Single class prediction

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| class | String | Yes | Class name |
| class_score | Float | Yes | Confidence score (0-1) |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| https://npclassifier.ucsd.edu | GET | Web interface |
| https://npclassifier.ucsd.edu/classify | GET | Single compound classification |
| https://npclassifier.ucsd.edu/classify | POST | Batch classification |

### Query Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| smiles | String | Yes | URL-encoded SMILES string |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | JSON (API responses) |
| Input | SMILES strings |
| Encoding | UTF-8 |
| Output | Nested JSON with arrays |

---

## Sample API Request

### Single Compound (GET)

```bash
# Classify curcumin
curl "https://npclassifier.ucsd.edu/classify?smiles=COC1%3DCC(%3DCC(%3DOC)C%3D1O)%2FC%3DC%2FC(%3DO)CC(%3DO)%2FC%3DC%2FC1%3DCC(%3DOC)C(%3DO)%3DCC%3D1"

# URL-decoded SMILES: COC1=CC(=CC(=OC)C=1O)/C=C/C(=O)CC(=O)/C=C/C1=CC(=OC)C(=O)=CC=1
```

### Single Compound (POST)

```bash
curl -X POST "https://npclassifier.ucsd.edu/classify" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC1=CC2=C(C=C1)NC=C2CCN"}'
```

---

## Sample Output Record

### Simple Natural Product (Tryptamine)

```json
{
  "smiles": "CC1=CC2=C(C=C1)NC=C2CCN",
  "pathway_results": [
    {"pathway": "Alkaloids", "pathway_score": 0.9823}
  ],
  "superclass_results": [
    {"superclass": "Indole alkaloids", "superclass_score": 0.9456}
  ],
  "class_results": [
    {"class": "Simple tryptamines and derivatives", "class_score": 0.8912}
  ],
  "isglycoside": false
}
```

### Complex Natural Product (Multi-pathway)

```json
{
  "smiles": "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
  "pathway_results": [
    {"pathway": "Terpenoids", "pathway_score": 0.9921},
    {"pathway": "Fatty acids", "pathway_score": 0.0234}
  ],
  "superclass_results": [
    {"superclass": "Steroids and steroid derivatives", "superclass_score": 0.9834},
    {"superclass": "Triterpenes", "superclass_score": 0.0089}
  ],
  "class_results": [
    {"class": "Cholestane steroids", "class_score": 0.9567}
  ],
  "isglycoside": false
}
```

### Glycoside Detection

```json
{
  "smiles": "OCC1OC(OC2=CC=C(CCO)C=C2)C(O)C(O)C1O",
  "pathway_results": [
    {"pathway": "Shikimates and Phenylpropanoids", "pathway_score": 0.8234},
    {"pathway": "Carbohydrates", "pathway_score": 0.6123}
  ],
  "superclass_results": [
    {"superclass": "Phenylpropanoids", "superclass_score": 0.7845}
  ],
  "class_results": [
    {"class": "Tyrosols and derivatives", "class_score": 0.7234}
  ],
  "isglycoside": true
}
```

---

## Pathway-Superclass-Class Hierarchy

### Terpenoids

| Superclass | Example Classes |
|------------|-----------------|
| Monoterpenoids | Acyclic monoterpenoids, Menthane monoterpenoids |
| Sesquiterpenoids | Eudesmane sesquiterpenoids, Germacrane sesquiterpenoids |
| Diterpenoids | Abietane diterpenoids, Labdane diterpenoids |
| Triterpenoids | Oleanane triterpenoids, Ursane triterpenoids |
| Steroids | Cholestane steroids, Androstane steroids |
| Carotenoids | C40 carotenoids, Apocarotenoids |

### Polyketides

| Superclass | Example Classes |
|------------|-----------------|
| Macrolides | 12-membered macrolides, 14-membered macrolides |
| Polyenes | Tetraenes, Pentaenes |
| Anthraquinones | Hydroxyanthraquinones |
| Tetracyclines | Classic tetracyclines |
| Naphthalenoids | Simple naphthalenoids |

### Alkaloids

| Superclass | Example Classes |
|------------|-----------------|
| Indole alkaloids | Simple tryptamines, Carbazole alkaloids |
| Isoquinoline alkaloids | Benzylisoquinolines, Aporphines |
| Pyridine alkaloids | Simple pyridines, Nicotinoids |
| Purine alkaloids | Xanthines, Purines |
| Tropane alkaloids | Tropanes, Ecgonine derivatives |

### Shikimates and Phenylpropanoids

| Superclass | Example Classes |
|------------|-----------------|
| Flavonoids | Flavones, Flavonols, Isoflavones, Anthocyanins |
| Coumarins | Simple coumarins, Furanocoumarins |
| Lignans | Dibenzylbutane lignans, Arylnaphthalene lignans |
| Phenolic acids | Hydroxybenzoic acids, Hydroxycinnamic acids |
| Stilbenoids | Monomeric stilbenoids, Oligomeric stilbenoids |

---

## Confidence Score Interpretation

| Score Range | Interpretation | Recommendation |
|-------------|----------------|----------------|
| 0.9 - 1.0 | High confidence | Accept classification |
| 0.7 - 0.9 | Moderate confidence | Review structure |
| 0.5 - 0.7 | Low confidence | Manual verification needed |
| < 0.5 | Very low confidence | Likely not natural product |

---

## Multiple Predictions

NPClassifier can return multiple predictions when structure fits multiple categories:

```json
{
  "pathway_results": [
    {"pathway": "Terpenoids", "pathway_score": 0.65},
    {"pathway": "Polyketides", "pathway_score": 0.45}
  ]
}
```

**Interpretation:** Mixed biosynthetic origin or hybrid natural product (meroterpenoid)

---

## Status Codes

| HTTP Code | Meaning |
|-----------|---------|
| 200 | Success |
| 400 | Invalid SMILES |
| 500 | Server error |
| 503 | Service unavailable |

---

## Glossary

| Term | Definition |
|------|------------|
| Pathway | Biosynthetic pathway of origin (terpenoid, polyketide, etc.) |
| Superclass | Broad structural/biosynthetic grouping |
| Class | Fine-grained structural classification |
| Glycoside | Compound with attached sugar moiety |
| NP | Natural Product |
| MVA | Mevalonate pathway (terpenoid biosynthesis) |
| MEP | Methylerythritol phosphate pathway |
| NRPS | Non-ribosomal peptide synthetase |
| BGC | Biosynthetic Gene Cluster |
| Meroterpenoid | Hybrid terpenoid-polyketide compound |
| Shikimate | Metabolic pathway producing aromatic amino acids |

---

## References

1. Kim HW, et al. (2021) "NPClassifier: A Deep Neural Network-Based Structural Classification Tool for Natural Products." J Nat Prod. 84(11):2795-2807. DOI: 10.1021/acs.jnatprod.1c00399
2. NPClassifier GitHub: https://github.com/mwang87/NP-Classifier
3. GNPS Integration: https://gnps.ucsd.edu

---

## Related Documents

- [NPClassifier Download Instructions](./download.md)
- [ClassyFire](../classyfire/README.md) - General chemical taxonomy
- [COCONUT](../../2.1.natural.products/coconut/README.md) - Natural products database
- [ChEBI](../chebi/README.md) - Chemical ontology
