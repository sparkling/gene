---
id: phytohub
title: "PhytoHub - Dietary Phytochemical Metabolome"
type: source
parent: ../README.md
tier: 2
status: active
category: compounds.molecules
subcategory: food.compounds.nutrients
tags:
  - phytochemicals
  - metabolites
  - biomarkers
  - dietary-exposure
  - mass-spectrometry
---

# PhytoHub - Dietary Phytochemical Metabolome

## Overview

PhytoHub is a comprehensive database of dietary phytochemicals and their human and microbial metabolites. The database catalogs the full spectrum of compounds derived from plant foods as they are absorbed and transformed in the human body, providing essential reference data for metabolomics studies of dietary exposure.

The database includes parent phytochemicals found in foods along with Phase I and Phase II metabolites (human liver metabolism) and microbial metabolites (gut microbiota transformation). Each compound entry includes chemical structures, mass spectrometry data, and links to food sources and biological activities.

PhytoHub supports the identification of unknown metabolites in metabolomics studies and enables the development of dietary biomarkers for nutritional epidemiology research.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Compounds | 1,800+ |
| Parent Phytochemicals | 600+ |
| Human Metabolites | 400+ |
| Microbial Metabolites | 600+ |
| Dietary Sources | 350+ |
| MS/MS Spectra | 1,200+ |

## Primary Use Cases

1. **Metabolomics Analysis** - Identify phytochemical metabolites in biofluids
2. **Biomarker Discovery** - Find biomarkers of dietary intake
3. **Metabolic Pathway Mapping** - Trace phytochemical biotransformation
4. **Mass Spectrometry Reference** - MS/MS spectral matching
5. **Microbiome Research** - Understand gut microbial metabolism

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| PhytoHub ID | PHUB + digits | PHUB000001 |
| InChI Key | 27 characters | Standard format |
| PubChem CID | Integer | Cross-referenced |
| HMDB ID | HMDB + digits | Cross-referenced |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://phytohub.eu | Search and browse |
| Compound Search | http://phytohub.eu/compounds | Structure search |
| MS/MS Library | http://phytohub.eu/spectra | Spectral data |
| Download | Available | Contact for access |

## Compound Categories

| Category | Description |
|----------|-------------|
| Parent Compounds | Original food phytochemicals |
| Phase I Metabolites | Oxidation, reduction products |
| Phase II Metabolites | Conjugates (glucuronides, sulfates) |
| Microbial Metabolites | Gut bacteria products |

## Limitations

- Commercial use requires contacting database maintainers
- MS/MS spectra coverage not complete for all compounds
- Metabolite pathways are partially characterized
- Individual variation in metabolism not captured

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions

## Related Resources

- [Phenol-Explorer](../phenol.explorer/_index.md) - Polyphenol content
- [HMDB](../../_index.md) - Human metabolome
- [USDA FoodData](../usda.fooddata/_index.md) - Food composition

## References

1. Neveu V, et al. (2010) "Phenol-Explorer: an online comprehensive database on polyphenol contents in foods." Database (Oxford). 2010:bap024.
