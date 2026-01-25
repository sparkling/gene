---
id: classyfire
title: "ClassyFire - Chemical Taxonomy"
type: source
parent: ../README.md
tier: 2
status: active
category: compounds.molecules
subcategory: chemical.ontology.classification
tags:
  - classification
  - taxonomy
  - chemical-structure
  - ontology
  - cheminformatics
---

# ClassyFire - Chemical Taxonomy

## Overview

ClassyFire is a web-based tool and hierarchical chemical taxonomy for the automatic structural classification of chemical compounds. It classifies compounds into chemical "superclasses", "classes", and "subclasses" based on their structural features, providing standardized vocabulary for chemical description.

The ClassyFire taxonomy consists of over 4,800 chemical categories organized in a hierarchical tree structure with 11 kingdoms at the top level. The classification is rule-based and purely structural, making it applicable to any chemical compound regardless of origin or use.

ClassyFire is widely used for metabolomics data annotation, chemical database organization, and standardizing chemical descriptions across publications and databases.

## Key Statistics

| Metric | Value |
|--------|-------|
| Chemical Categories | 4,800+ |
| Kingdoms | 11 |
| Superclasses | 37 |
| Classes | 150+ |
| Subclasses | 400+ |
| Classification Rules | 7,000+ |

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| SMILES | String | CC(=O)OC1=CC=CC=C1C(=O)O |
| InChI | InChI=... | InChI=1S/C9H8O4/... |
| InChIKey | 27-char | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| ClassyFire ID | CHEMONTID:nnnnnnnn | CHEMONTID:0000311 |

## Primary Use Cases

1. **Metabolomics Annotation** - Classify unknown metabolites by structure
2. **Database Organization** - Standardize chemical categorization
3. **Literature Mining** - Consistent terminology for compounds
4. **Chemical Space Analysis** - Characterize compound libraries
5. **Data Integration** - Common classification across resources

## Classification Levels

| Level | Example (Aspirin) |
|-------|-------------------|
| Kingdom | Organic compounds |
| Superclass | Benzenoids |
| Class | Benzene and substituted derivatives |
| Subclass | Phenol esters |
| Direct Parent | Phenyl acetates |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://classyfire.wishartlab.com | Interactive tool |
| Web API | http://classyfire.wishartlab.com/api | REST classification |
| Batch Processing | Available | Multiple compounds |
| InChIKey Lookup | Available | Pre-computed results |

## Chemical Kingdoms

| Kingdom | Description |
|---------|-------------|
| Organic compounds | Carbon-based |
| Inorganic compounds | Non-carbon based |
| Organometallic compounds | Metal-carbon bonds |
| Organic salts | Ionic organic |
| Homogeneous mixtures | Single phase |
| Heterogeneous mixtures | Multiple phases |

## Limitations

- Commercial use requires contacting database maintainers
- Classification is purely structural; biological context not considered
- Complex mixtures and salts may have inconsistent classification
- API rate limits may affect high-throughput classification

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [ChEBI](../chebi/README.md) - Chemical ontology
- [NPClassifier](../npclassifier/README.md) - Natural product classification
- [PubChem](../pubchem/README.md) - Chemical repository

## References

1. Djoumbou Feunang Y, et al. (2016) "ClassyFire: automated chemical classification with a comprehensive, computable taxonomy." J Cheminform. 8:61.
