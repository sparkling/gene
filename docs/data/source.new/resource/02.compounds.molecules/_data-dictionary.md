# Category 02: Compounds and Molecules - Data Dictionary

## Overview

This data dictionary documents the unified data schema for Category 02: Compounds and Molecules. This category encompasses chemical compounds from natural products, pharmaceuticals, food compounds, drug metabolism data, chemical classifications, and compound-target interactions.

**Category ID:** 02
**Category Name:** Compounds and Molecules
**Extraction Date:** 2026-01-24

## Subcategories

| ID | Name | Data Sources |
|----|------|--------------|
| 2.1 | Natural Products | COCONUT, Dr. Duke's, LOTUS, NPASS, NPAtlas |
| 2.2 | Pharmaceuticals | ChEMBL, DailyMed, DrugBank, Orange Book, RxNorm |
| 2.4 | Food Compounds and Nutrients | Phenol-Explorer, PhytoHub, USDA FoodData |
| 2.5 | Drug Metabolism and Pharmacokinetics | SuperCYP, SwissADME |
| 2.6 | Chemical Ontology and Classification | ChEBI, ClassyFire, NPClassifier, PubChem |
| 2.7 | Compound-Target Interactions | BindingDB, DGIdb, GtoPdb, TTD |

---

## Common Chemical Identifiers

These identifiers appear across multiple subcategories and enable cross-database linking.

| Field Name | Data Type | Format | Description | Found In |
|------------|-----------|--------|-------------|----------|
| inchi_key | string | 27 characters (XXXXXXXXXXXXXX-YYYYYYYYYY-Z) | InChI Key - primary cross-reference identifier for exact structure matching | 2.1, 2.2, 2.4, 2.5, 2.6, 2.7 |
| canonical_smiles | string | Variable | Canonical SMILES - standardized structure notation without stereochemistry | 2.1, 2.2, 2.4, 2.5, 2.6, 2.7 |
| inchi | string | Variable | International Chemical Identifier - IUPAC standard structure representation | 2.1, 2.2, 2.6 |
| cas_number | string | #-##-# or ##-##-# format | CAS Registry Number - Chemical Abstracts Service identifier | 2.1, 2.2, 2.4 |
| pubchem_cid | integer | Positive integer | PubChem Compound ID - links to PubChem database | 2.1, 2.2, 2.4, 2.6 |
| chembl_id | string | CHEMBL + digits | ChEMBL identifier - links to bioactivity data | 2.2, 2.7 |
| drugbank_id | string | DB + 5 digits | DrugBank identifier - links to comprehensive drug data | 2.2, 2.5 |
| uniprot_id | string | 6-10 characters | UniProt accession - links to protein target data | 2.2, 2.5, 2.7 |
| chebi_id | string | CHEBI:digits | ChEBI identifier - links to chemical ontology | 2.4, 2.6 |

---

## Common Molecular Properties

These property fields are found across multiple subcategories.

| Field Name | Data Type | Unit | Description | Found In |
|------------|-----------|------|-------------|----------|
| molecular_weight | decimal | Daltons (Da) | Average molecular mass considering isotope distribution | 2.1, 2.2, 2.4, 2.5, 2.6, 2.7 |
| molecular_formula | string | e.g., C9H8O4 | Elemental composition of the compound | 2.1, 2.2, 2.4, 2.5, 2.6 |
| logp | decimal | dimensionless | Partition coefficient (octanol/water) | 2.1, 2.2, 2.5, 2.6 |
| tpsa | decimal | Angstrom squared | Topological polar surface area | 2.1, 2.2, 2.5, 2.6 |
| hbd | integer | count | Hydrogen bond donors | 2.1, 2.2, 2.5, 2.7 |
| hba | integer | count | Hydrogen bond acceptors | 2.1, 2.2, 2.5, 2.7 |
| rotatable_bonds | integer | count | Number of rotatable bonds | 2.1, 2.2, 2.5, 2.7 |

---

## Cross-Reference Guidelines

### Structure Matching
- **Primary Key:** Use `inchi_key` for exact structure matching across databases
- **Secondary Key:** Use `canonical_smiles` for structure-based queries
- **Fallback:** Use `molecular_formula` + `molecular_weight` for approximate matching

### Database Linking
- **Drug Data:** Link via `drugbank_id` or `chembl_id`
- **Protein Targets:** Link via `uniprot_id`
- **Chemical Ontology:** Link via `chebi_id`
- **Literature:** Link via `pubchem_cid` for patent and publication data

### Identifier Precedence
1. InChI Key (most reliable for structure matching)
2. Database-specific IDs (DrugBank, ChEMBL, etc.)
3. CAS Number (industry standard but not always available)
4. PubChem CID (comprehensive but may have duplicates)

---

## Data Quality Notes

### Cardinality
- **1:1** - Single value per compound
- **1:N** - Multiple values allowed (arrays/lists)
- **Required** - Field must have a value
- **Optional** - Field may be null

### Source Coverage
- Not all sources provide all fields
- Check individual subcategory dictionaries for source-specific mappings
- Some fields require transformation during ETL

---

## Related Documentation

- [2.1 Natural Products Data Dictionary](./2.1.natural.products/_data-dictionary.md)
- [2.2 Pharmaceuticals Data Dictionary](./2.2.pharmaceuticals/_data-dictionary.md)
- [2.4 Food Compounds Data Dictionary](./2.4.food.compounds.nutrients/_data-dictionary.md)
- [2.5 Drug Metabolism Data Dictionary](./2.5.drug.metabolism.pharmacokinetics/_data-dictionary.md)
- [2.6 Chemical Ontology Data Dictionary](./2.6.chemical.ontology.classification/_data-dictionary.md)
- [2.7 Compound-Target Interactions Data Dictionary](./2.7.compound.target.interactions/_data-dictionary.md)
