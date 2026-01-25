---
id: npclassifier
title: "NPClassifier - Natural Product Classification"
type: source
parent: ../README.md
tier: 2
status: active
category: compounds.molecules
subcategory: chemical.ontology.classification
tags:
  - natural-products
  - classification
  - biosynthesis
  - pathway
  - taxonomy
---

# NPClassifier - Natural Product Classification

## Overview

NPClassifier is a deep learning-based tool for automated structural classification of natural products based on their predicted biosynthetic pathways. Unlike general chemical classifiers, NPClassifier specifically addresses the unique structural diversity and biosynthetic origins of natural products.

The tool classifies compounds into biosynthetic pathway classes (e.g., terpenoids, polyketides, alkaloids), superclasses, and fine-grained classes. This classification reflects how organisms synthesize these compounds, providing insights into their biological origin and potential bioactivities.

NPClassifier is particularly valuable for natural product dereplication, chemotaxonomy research, and connecting metabolomics data to biosynthetic gene cluster predictions.

## Key Statistics

| Metric | Value |
|--------|-------|
| Pathway Classes | 7 |
| Superclasses | 35+ |
| Classes | 200+ |
| Training Compounds | 80,000+ |
| Model Accuracy | >90% |

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| SMILES | String | CC1=CC2=C(C=C1)NC=C2CCN |
| InChI | InChI=... | InChI=1S/C11H14N2/... |
| InChIKey | 27-char | QEVHRUUCFGRFIF-UHFFFAOYSA-N |
| Pathway Class | Text | Terpenoids, Alkaloids |

## Primary Use Cases

1. **Natural Product Annotation** - Classify NP structures by biosynthetic origin
2. **Chemotaxonomy** - Study chemical patterns across taxa
3. **BGC Correlation** - Link structures to biosynthetic gene clusters
4. **Dereplication** - Rapid NP class identification
5. **Metabolomics** - Annotate unknown natural products

## Pathway Classes

| Pathway | Examples |
|---------|----------|
| Terpenoids | Steroids, carotenoids, monoterpenes |
| Polyketides | Macrolides, tetracyclines |
| Alkaloids | Indole alkaloids, pyridines |
| Amino acids/Peptides | NRPS products, cyclopeptides |
| Carbohydrates | Glycosides, aminoglycosides |
| Fatty acids | Lipids, acyl derivatives |
| Shikimates | Phenylpropanoids, flavonoids |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://npclassifier.ucsd.edu | Interactive tool |
| API | https://npclassifier.ucsd.edu/classify | REST endpoint |
| GNPS Integration | Available | Mass spectrometry workflows |
| Source Code | https://github.com/mwang87/NP-Classifier | Open source |

## Output Format

| Field | Description |
|-------|-------------|
| pathway | Biosynthetic pathway class |
| superclass | Structural superclass |
| class | Fine-grained class |
| isglycoside | Glycoside detection |
| confidence | Classification confidence |

## Limitations

- Deep learning model may misclassify unusual structures
- Synthetic compounds not designed for NP classification
- Biosynthetic pathway predictions are computational estimates
- Glycoside detection may miss non-standard sugar moieties

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [ClassyFire](../classyfire/_index.md) - General chemical taxonomy
- [COCONUT](../../2.1.natural.products/coconut/_index.md) - Natural products
- [ChEBI](../chebi/_index.md) - Chemical ontology

## References

1. Kim HW, et al. (2021) "NPClassifier: A Deep Neural Network-Based Structural Classification Tool for Natural Products." J Nat Prod. 84(11):2795-2807.
