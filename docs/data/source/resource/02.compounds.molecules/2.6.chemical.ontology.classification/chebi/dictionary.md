# ChEBI - Data Dictionary

## Overview

This data dictionary documents the schema for ChEBI (Chemical Entities of Biological Interest) ontology database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | chebi |
| **Name** | ChEBI |
| **Parent** | 2.6.chemical.ontology.classification |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Entity Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| CHEBI_ID | string | 1:1 | Yes | Primary identifier | CHEBI:17234 |
| name | string | 1:1 | Yes | Primary name | D-glucose |
| definition | string | 1:1 | No | Textual definition | An aldohexose used as energy source |
| synonyms | array | 1:N | No | Alternative names | dextrose, grape sugar |
| status | string | 1:1 | Yes | Curation status | CHECKED |
| star | integer | 1:1 | Yes | Curation level (1-3) | 3 |

### Chemical Structure Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| formula | string | 1:1 | No | Molecular formula | C6H12O6 |
| mass | decimal | 1:1 | No | Average mass (Da) | 180.15588 |
| monoisotopic_mass | decimal | 1:1 | No | Monoisotopic mass | 180.06339 |
| charge | integer | 1:1 | No | Formal charge | 0 |
| InChI | string | 1:1 | No | IUPAC InChI | InChI=1S/C6H12O6/... |
| InChIKey | string | 1:1 | No | Hashed InChI (27 chars) | WQZGKKKJIJFFOK-GASJEMHNSA-N |
| SMILES | string | 1:1 | No | SMILES notation | OC[C@H]1OC(O)... |
| mol_file | string | 1:1 | No | MDL MOL block | V2000/V3000 format |

### Ontology Relationship Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| is_a | array | 1:N | No | Parent classes | CHEBI:4194 (D-aldohexose) |
| has_role | array | 1:N | No | Biological/chemical roles | CHEBI:77746 (human metabolite) |
| has_part | array | 1:N | No | Structural components | - |
| is_conjugate_acid_of | string | 1:1 | No | Conjugate base | CHEBI:18391 |
| is_conjugate_base_of | string | 1:1 | No | Conjugate acid | - |
| is_enantiomer_of | string | 1:1 | No | Mirror image isomer | CHEBI:37627 (L-glucose) |
| is_tautomer_of | string | 1:1 | No | Tautomeric form | - |

### Cross-Reference Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| KEGG_COMPOUND | string | 1:1 | No | KEGG Compound ID | C00031 |
| PubChem_CID | integer | 1:1 | No | PubChem Compound ID | 5793 |
| HMDB | string | 1:1 | No | Human Metabolome DB | HMDB0000122 |
| DrugBank | string | 1:1 | No | DrugBank ID | DB09341 |
| LIPID_MAPS | string | 1:1 | No | LIPID MAPS ID | - |
| CAS | string | 1:1 | No | CAS Registry Number | 50-99-7 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| ChEBI ID | CHEBI: + digits | CHEBI:17234 | Primary identifier |
| InChIKey | 27 characters | WQZGKKKJIJFFOK-GASJEMHNSA-N | Structure hash |
| CAS Number | varies | 50-99-7 | Chemical registry |
| KEGG ID | C + 5 digits | C00031 | Pathway database |
| HMDB ID | HMDB + 7 digits | HMDB0000122 | Metabolome database |
| PubChem CID | Integer | 5793 | Chemical database |

---

## Enumerations

### Curation Star Levels

| Star Rating | Description | Count |
|-------------|-------------|-------|
| 3-star | Fully manually curated with verified structure | 50,000+ |
| 2-star | Partially curated with some manual verification | 30,000+ |
| 1-star | Automatically generated from external sources | 95,000+ |

### Status Values

| Status | Description |
|--------|-------------|
| CHECKED | Manually verified entry |
| PRELIMINARY | Initial entry, not yet reviewed |
| DELETED | Removed entry |

### Chemical Ontology Hierarchy

| Level | Parent | Examples |
|-------|--------|----------|
| chemical entity | root | CHEBI:24431 |
| molecular entity | chemical entity | CHEBI:23367 |
| molecule | molecular entity | CHEBI:25367 |
| organic molecular entity | molecule | CHEBI:50860 |
| inorganic molecular entity | molecule | CHEBI:24835 |

### Major Chemical Classes

| Class | ChEBI ID | Description |
|-------|----------|-------------|
| carbohydrate | CHEBI:16646 | Sugars and derivatives |
| lipid | CHEBI:18059 | Fatty compounds |
| amino acid | CHEBI:33709 | Protein building blocks |
| nucleotide | CHEBI:36976 | DNA/RNA components |
| organic acid | CHEBI:64709 | Acidic organics |
| alkaloid | CHEBI:22315 | Nitrogen-containing naturals |

### Biological Role Categories

| Role Category | ChEBI ID | Examples |
|---------------|----------|----------|
| biochemical role | CHEBI:52206 | cofactor, substrate, metabolite |
| pharmacological role | CHEBI:52217 | drug, antibiotic, analgesic |
| biological role | CHEBI:24432 | hormone, neurotransmitter |

### Specific Roles

| Role | ChEBI ID | Description |
|------|----------|-------------|
| human metabolite | CHEBI:77746 | Found in human metabolism |
| plant metabolite | CHEBI:76924 | Found in plant metabolism |
| bacterial metabolite | CHEBI:76969 | Found in bacterial metabolism |
| cofactor | CHEBI:23357 | Enzyme helper molecule |
| enzyme inhibitor | CHEBI:23924 | Blocks enzyme activity |
| drug | CHEBI:23888 | Therapeutic compound |

### Chemical Role Categories

| Role | ChEBI ID | Description |
|------|----------|-------------|
| antioxidant | CHEBI:22586 | Prevents oxidation |
| buffer | CHEBI:22695 | pH stabilizer |
| solvent | CHEBI:46787 | Dissolving medium |
| reducing agent | CHEBI:63247 | Electron donor |

### Ontology Relationship Types

| Relationship | Description |
|--------------|-------------|
| is_a | Hierarchical subclass |
| has_role | Function assignment |
| has_part | Structural component |
| is_conjugate_acid_of | Proton donor pair |
| is_conjugate_base_of | Proton acceptor pair |
| is_enantiomer_of | Mirror image stereoisomer |
| is_tautomer_of | Structural isomer |

---

## Entity Relationships

### Entity to Parent Classes
- **Cardinality:** N:M
- **Description:** Entities may have multiple parent classes
- **Key Fields:** CHEBI_ID, is_a

### Entity to Roles
- **Cardinality:** 1:N
- **Description:** Each entity can have multiple biological/chemical roles
- **Key Fields:** CHEBI_ID, has_role

### Entity to Synonyms
- **Cardinality:** 1:N
- **Description:** Each entity has multiple synonyms
- **Key Fields:** CHEBI_ID, synonyms

### Entity to Cross-References
- **Cardinality:** 1:N
- **Description:** Links to external databases
- **Key Fields:** CHEBI_ID, external database IDs

### Conjugate Acid/Base Pairs
- **Cardinality:** 1:1
- **Description:** Proton donor/acceptor relationships
- **Key Fields:** is_conjugate_acid_of, is_conjugate_base_of

### Enantiomeric Relationships
- **Cardinality:** 1:1
- **Description:** Mirror image stereoisomers
- **Key Fields:** CHEBI_ID, is_enantiomer_of

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ChEBI | Chemical Entities of Biological Interest | EMBL-EBI chemical ontology |
| InChI | IUPAC International Chemical Identifier | Standard structure ID |
| SMILES | Simplified Molecular-Input Line-Entry System | Line notation |
| OBO | Open Biological Ontologies | Ontology format |
| OWL | Web Ontology Language | Semantic web format |
| SDF | Structure-Data File | MDL chemical format |
| OLS | Ontology Lookup Service | EBI ontology browser |
| SOAP | Simple Object Access Protocol | Web service protocol |
| CAS | Chemical Abstracts Service | Registry number source |
| HMDB | Human Metabolome Database | Cross-reference DB |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |
| CC BY | Creative Commons Attribution | License type |
| EBI | European Bioinformatics Institute | Database maintainer |
| EMBL | European Molecular Biology Laboratory | Parent organization |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PubChem | CID | Chemical data |
| KEGG | Compound ID | Pathway context |
| HMDB | HMDB ID | Metabolomics |
| DrugBank | DrugBank ID | Drug information |
| LIPID MAPS | LM ID | Lipid classification |
| Reactome | Entity ID | Pathway reactions |
| Rhea | Compound | Reaction database |
| UniProt | Annotation | Protein ligands |
| MetaCyc | Compound | Metabolic pathways |
| GO | Function terms | Molecular function |

---

## Data Quality Notes

1. **Comprehensive Coverage:** 175,000+ chemical entities
2. **Curation Quality:** 50,000+ fully curated (3-star) entries
3. **Role Ontology:** 2,500+ biological and chemical roles
4. **Structure Data:** Complete InChI, SMILES, formula, mass
5. **Multiple Formats:** OBO, OWL, SDF, TSV, JSON API
6. **Open Access:** CC BY 4.0 license
7. **Update Frequency:** Monthly updates
8. **API Access:** REST, SOAP, OLS interfaces

