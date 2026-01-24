# 2.6 Chemical Ontology and Classification - Data Dictionary

## Overview

This data dictionary documents the unified schema for chemical ontology and classification data from four major databases: ChEBI, ClassyFire, NPClassifier, and PubChem.

**Subcategory ID:** 2.6
**Subcategory Name:** Chemical Ontology and Classification
**Data Sources:** ChEBI, ClassyFire, NPClassifier, PubChem

---

## Unified Fields

### Core Identification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| compound_id | string | 1:1, Required | Primary identifier for compound/concept | `CHEBI:17234`, `CHEMONTID:0000300`, `2244` | ChEBI: CHEBI_ID, ClassyFire: chemont_id, PubChem: CID |
| name | string | 1:1, Required | Compound or classification name | `Glucose`, `Aspirin`, `Organic acids` | ChEBI, ClassyFire, NPClassifier, PubChem |

### Structure Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| canonical_smiles | string | 1:1, Optional | Canonical SMILES structure | `CC(=O)Oc1ccccc1C(=O)O` | ChEBI, ClassyFire, NPClassifier, PubChem |
| inchi | string | 1:1, Optional | InChI identifier | `InChI=1S/C9H8O4/c1-6(10)...` | ChEBI, PubChem |
| inchi_key | string | 1:1, Optional | InChI Key | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` | ChEBI, ClassyFire, PubChem |
| molecular_formula | string | 1:1, Optional | Molecular formula | `C9H8O4` | ChEBI, PubChem |

### Classification Fields

| Field Name | Data Type | Cardinality | Description | Example Values | Source Mappings |
|------------|-----------|-------------|-------------|----------------|-----------------|
| classification_hierarchy | object[] | 1:N, Optional | Hierarchical classification from kingdom to class | `[{"level": "kingdom", "name": "Organic compounds"}, {"level": "superclass", "name": "Benzenoids"}]` | ChEBI, ClassyFire, NPClassifier |
| cross_references | object[] | 1:N, Optional | Links to external databases | `[{"database": "KEGG", "id": "C00243"}]` | ChEBI, PubChem |

---

## Source-Specific Fields

### ChEBI

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| chebi_id | string | 1:1, Required | ChEBI identifier (CHEBI:##### format) | `CHEBI:17234`, `CHEBI:15377` |
| definition | string | 1:1, Optional | Textual definition of the entity | `A monosaccharide containing six carbon atoms...` |
| synonyms | string[] | 1:N, Optional | Alternative names | `["D-glucose", "dextrose", "blood sugar"]` |
| mass | decimal | 1:1, Optional | Monoisotopic mass | `180.0634` |
| charge | integer | 1:1, Optional | Formal charge | `-1`, `0`, `+1` |
| star | integer | 1:1, Optional | Curation level (1=low, 2=medium, 3=high) | `1`, `2`, `3` |
| is_a | string[] | 1:N, Optional | Parent classes in ontology hierarchy | `["CHEBI:35381", "CHEBI:33285"]` |
| has_role | string[] | 1:N, Optional | Biological/chemical roles | `["CHEBI:35610 (drug)", "CHEBI:27780 (metabolite)"]` |
| is_conjugate_acid_of | string | 1:1, Optional | Conjugate base reference | `CHEBI:28866` |
| is_enantiomer_of | string | 1:1, Optional | Mirror image isomer reference | `CHEBI:15904` |

### ClassyFire

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| chemont_id | string | 1:1, Required | ChemOnt identifier | `CHEMONTID:0000300` |
| cf_kingdom | object | 1:1, Optional | Top-level classification | `{"name": "Organic compounds", "chemont_id": "CHEMONTID:0000000"}` |
| cf_superclass | object | 1:1, Optional | Second-level classification | `{"name": "Benzenoids", "chemont_id": "CHEMONTID:0002448"}` |
| cf_class | object | 1:1, Optional | Third-level classification | `{"name": "Benzene and substituted derivatives", "chemont_id": "..."}` |
| cf_subclass | object | 1:1, Optional | Fourth-level classification | `{"name": "Benzoic acids and derivatives", "chemont_id": "..."}` |
| cf_direct_parent | object | 1:1, Optional | Most specific classification | `{"name": "Salicylic acids", "chemont_id": "..."}` |
| cf_alternative_parents | object[] | 1:N, Optional | Additional valid classification paths | `[{"name": "Phenol ethers", "chemont_id": "..."}]` |
| cf_substituents | string[] | 1:N, Optional | Functional group substituents | `["Carboxylic acid", "Phenol ether", "Benzene"]` |
| cf_molecular_framework | string | 1:1, Optional | Structural framework description | `Aromatic homomonocyclic compounds` |

### NPClassifier

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| npc_pathway | string | 1:1, Optional | Biosynthetic pathway (7 types) | `Shikimates and Phenylpropanoids`, `Terpenoids`, `Alkaloids`, `Polyketides`, `Fatty acids`, `Amino acids and Peptides`, `Carbohydrates` |
| npc_pathway_score | decimal | 1:1, Optional | Pathway confidence (0-1) | `0.95` |
| npc_superclass | string | 1:1, Optional | NP superclass (35+ types) | `Flavonoids`, `Terpene lactones`, `Indole alkaloids` |
| npc_superclass_score | decimal | 1:1, Optional | Superclass confidence (0-1) | `0.92` |
| npc_class | string | 1:1, Optional | NP class (200+ types) | `Flavonols`, `Sesquiterpene lactones`, `Beta-carbolines` |
| npc_class_score | decimal | 1:1, Optional | Class confidence (0-1) | `0.88` |
| is_glycoside | boolean | 1:1, Optional | Glycoside moiety detected | `true`, `false` |

### PubChem

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| pubchem_cid | integer | 1:1, Required | PubChem Compound ID | `2244`, `5988` |
| iupac_name | string | 1:1, Optional | IUPAC systematic name | `2-acetyloxybenzoic acid` |
| isomeric_smiles | string | 1:1, Optional | SMILES with stereochemistry | `C[C@@H](O)c1ccccc1` |
| xlogp | decimal | 1:1, Optional | Computed LogP | `1.2` |
| exact_mass | decimal | 1:1, Optional | Monoisotopic mass | `180.0423` |
| tpsa | decimal | 1:1, Optional | Topological polar surface area | `63.6` |
| complexity | decimal | 1:1, Optional | Molecular complexity score | `212` |
| bioactivity | object[] | 1:N, Optional | Bioassay results | `[{"assay_id": "AID1234", "activity": "active", "target": "COX-2"}]` |
| patents | object[] | 1:N, Optional | Patent references | `[{"patent_id": "US1234567", "title": "..."}]` |
| literature | object[] | 1:N, Optional | PubMed references | `[{"pmid": 12345678, "title": "..."}]` |
| fingerprint2d | string | 1:1, Optional | PubChem fingerprint (binary string) | `0000001100001010...` |

---

## Metadata Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| _source.database | string | 1:1, Required | Source database name | `ChEBI`, `ClassyFire`, `NPClassifier`, `PubChem` |
| _source.version | string | 1:1, Optional | Database version | `231`, `2.1`, `2024-01` |
| _source.access_date | date | 1:1, Optional | Date data was accessed | `2026-01-24` |
| _source.original_id | string | 1:1, Optional | Original identifier from source | - |

---

## Field Mappings by Source

### ChEBI to Unified Schema

| ChEBI Field | Unified Field |
|-------------|---------------|
| CHEBI_ID | chebi_id |
| NAME | name |
| SMILES | canonical_smiles |
| InChI | inchi |
| InChIKey | inchi_key |
| FORMULA | molecular_formula |
| DEFINITION | definition |
| SYNONYM | synonyms |
| MASS | mass |
| CHARGE | charge |
| STAR | star |
| IS_A | is_a |
| HAS_ROLE | has_role |
| IS_CONJUGATE_ACID_OF | is_conjugate_acid_of |
| IS_ENANTIOMER_OF | is_enantiomer_of |

### ClassyFire to Unified Schema

| ClassyFire Field | Unified Field |
|------------------|---------------|
| chemont_id | chemont_id |
| smiles | canonical_smiles |
| inchikey | inchi_key |
| kingdom | cf_kingdom |
| superclass | cf_superclass |
| class | cf_class |
| subclass | cf_subclass |
| direct_parent | cf_direct_parent |
| alternative_parents | cf_alternative_parents |
| substituents | cf_substituents |
| molecular_framework | cf_molecular_framework |

### NPClassifier to Unified Schema

| NPClassifier Field | Unified Field |
|--------------------|---------------|
| smiles | canonical_smiles |
| pathway | npc_pathway |
| pathway_score | npc_pathway_score |
| superclass | npc_superclass |
| superclass_score | npc_superclass_score |
| class | npc_class |
| class_score | npc_class_score |
| isglycoside | is_glycoside |

### PubChem to Unified Schema

| PubChem Field | Unified Field |
|---------------|---------------|
| CID | pubchem_cid |
| Title | name |
| CanonicalSMILES | canonical_smiles |
| IsomericSMILES | isomeric_smiles |
| InChI | inchi |
| InChIKey | inchi_key |
| MolecularFormula | molecular_formula |
| IUPACName | iupac_name |
| XLogP | xlogp |
| ExactMass | exact_mass |
| TPSA | tpsa |
| Complexity | complexity |
| bioactivity | bioactivity |
| patents | patents |
| literature | literature |
| Fingerprint2D | fingerprint2d |

---

## Data Quality Notes

- **Required Field:** Only `name` is strictly required across all sources
- **Ontology Data:** ChEBI provides formal ontology with parent-child relationships
- **Classification:** ClassyFire and NPClassifier provide hierarchical chemical classification
- **Cross-References:** ChEBI and PubChem provide extensive external database links
- **Bioactivity:** Only PubChem provides bioassay results

### ChEBI Star Levels

| Star | Curation Level | Description |
|------|----------------|-------------|
| 1 | Minimal | Automatically generated; basic data only |
| 2 | Intermediate | Partially curated; some manual review |
| 3 | Full | Fully manually curated; high quality |

### ChEBI Relationship Types

| Relationship | Description |
|--------------|-------------|
| is_a | Subsumption relationship (parent class) |
| has_part | Part-whole relationship |
| has_role | Biological or chemical role |
| is_conjugate_acid_of | Acid-base pair relationship |
| is_conjugate_base_of | Base-acid pair relationship |
| is_enantiomer_of | Stereoisomer relationship |
| is_tautomer_of | Tautomeric relationship |
| has_functional_parent | Functional group relationship |

### ClassyFire Classification Hierarchy

| Level | Description | Example |
|-------|-------------|---------|
| Kingdom | Top level | Organic compounds, Inorganic compounds |
| Superclass | Second level | Benzenoids, Lipids, Alkaloids |
| Class | Third level | Benzene derivatives, Fatty acids |
| Subclass | Fourth level | Benzoic acids, Saturated fatty acids |
| Direct Parent | Most specific | Salicylic acids, Palmitic acid |

### NPClassifier Biosynthetic Pathways

| Pathway | Description |
|---------|-------------|
| Shikimates and Phenylpropanoids | Aromatic compounds from shikimate pathway |
| Terpenoids | Isoprenoid-derived compounds |
| Alkaloids | Nitrogen-containing compounds |
| Polyketides | Compounds from polyketide synthases |
| Fatty acids | Lipid-derived compounds |
| Amino acids and Peptides | Protein building blocks and peptides |
| Carbohydrates | Sugar-based compounds |

### PubChem Data Sources

PubChem aggregates data from >800 sources including:
- Chemical vendors
- Government agencies (FDA, EPA, NIH)
- Academic databases
- Patent offices
- Published literature
