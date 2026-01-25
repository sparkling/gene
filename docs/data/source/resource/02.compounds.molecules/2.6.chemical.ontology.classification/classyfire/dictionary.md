# ClassyFire - Data Dictionary

## Overview

This data dictionary documents the schema for ClassyFire chemical taxonomy classification service.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | classyfire |
| **Name** | ClassyFire |
| **Parent** | 2.6.chemical.ontology.classification |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Classification Result

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| id | integer | 1:1 | Yes | ClassyFire query ID | 123456 |
| identifier | string | 1:1 | Yes | Input identifier/label | aspirin |
| smiles | string | 1:1 | Yes | Canonical SMILES | CC(=O)Oc1ccccc1C(=O)O |
| inchikey | string | 1:1 | Yes | InChIKey | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |
| molecular_framework | string | 1:1 | No | Structural framework | Aromatic homomonocyclic compounds |
| description | string | 1:1 | No | Compound description | - |
| classification_version | string | 1:1 | Yes | ChemOnt version | 2.1 |

### Classification Level Object

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| name | string | 1:1 | Yes | Category name | Benzenoids |
| description | string | 1:1 | Yes | Category description | Aromatic compounds... |
| chemont_id | string | 1:1 | Yes | ChemOnt identifier | CHEMONTID:0000009 |
| url | string | 1:1 | Yes | Link to category page | http://classyfire.wishartlab.com/... |

### Hierarchy Fields

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| kingdom | object | 1:1 | Yes | Kingdom classification | Organic compounds |
| superclass | object | 1:1 | Yes | Superclass classification | Benzenoids |
| class | object | 1:1 | Yes | Class classification | Benzene derivatives |
| subclass | object | 1:1 | No | Subclass classification | Phenol esters |
| intermediate_nodes | array | 1:N | No | Intermediate categories | - |
| direct_parent | object | 1:1 | Yes | Most specific classification | Phenyl acetates |
| alternative_parents | array | 1:N | No | Alternative paths | Benzoic acids... |
| substituents | array | 1:N | No | Functional groups | Phenyl acetate, Carboxylic acid |

### Predicted Mappings

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| predicted_chebi_terms | array | 1:N | No | Predicted ChEBI mappings | CHEBI:15365 |
| predicted_lipidmaps_terms | array | 1:N | No | Predicted LIPID MAPS | - |
| external_descriptors | array | 1:N | No | External annotations | - |

### Query Status

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| status | string | 1:1 | Yes | Query status | in progress, Done, Error |
| classification_status | string | 1:1 | Yes | Classification state | classified |
| invalid_entities | array | 1:N | No | Failed compounds | - |
| entities | array | 1:N | No | Classified compounds | - |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| ChemOnt ID | CHEMONTID: + 7 digits | CHEMONTID:0000300 | Taxonomy node ID |
| Query ID | Integer | 123456 | Classification job ID |
| InChIKey | 27 characters | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | Structure identifier |
| URL ID | C + 7 digits | C0000300 | Web URL format |

---

## Enumerations

### Kingdoms (11 Total)

| Kingdom | ChemOnt ID | Description |
|---------|------------|-------------|
| Organic compounds | CHEMONTID:0000000 | Carbon-based compounds |
| Inorganic compounds | CHEMONTID:0000001 | Non-carbon based |
| Organometallic compounds | CHEMONTID:0000002 | Metal-carbon bonds |
| Organic salts | CHEMONTID:0000003 | Ionic organic compounds |
| Organic acids | CHEMONTID:0000004 | Acidic organics |
| Nucleosides, nucleotides, analogues | CHEMONTID:0000005 | Nucleic acid components |
| Lipids and lipid-like molecules | CHEMONTID:0000006 | Fatty compounds |
| Organic nitrogen compounds | CHEMONTID:0000007 | Nitrogen-containing |
| Organoheterocyclic compounds | CHEMONTID:0000008 | Heterocycles |
| Benzenoids | CHEMONTID:0000009 | Benzene derivatives |
| Phenylpropanoids and polyketides | CHEMONTID:0000010 | Natural products |

### Major Superclasses (37 Total)

| Superclass | Parent Kingdom | Examples |
|------------|----------------|----------|
| Benzenoids | Organic compounds | Benzene derivatives |
| Lipids and lipid-like | Organic compounds | Fatty acids, steroids |
| Organic acids | Organic compounds | Carboxylic acids |
| Organoheterocyclic | Organic compounds | Pyridines, furans |
| Organic nitrogen | Organic compounds | Amines, imines |
| Phenylpropanoids | Organic compounds | Flavonoids, stilbenes |
| Alkaloids | Organic compounds | Caffeine, morphine |

### Classification Hierarchy Levels

| Level | Count | Description |
|-------|-------|-------------|
| Kingdom | 11 | Broadest category |
| Superclass | 37 | Major structural class |
| Class | 150+ | Specific structure type |
| Subclass | 400+ | Detailed category |
| Intermediate | Variable | Additional levels |
| Direct Parent | 4,800+ | Most specific category |

### Query Status Values

| Status | Description |
|--------|-------------|
| in progress | Classification running |
| Done | Classification complete |
| Error | Classification failed |

### HTTP Status Codes

| Code | Meaning |
|------|---------|
| 200 | Success |
| 201 | Created (submission) |
| 404 | Not found |
| 429 | Rate limited |
| 500 | Server error |

---

## Entity Relationships

### Compound to Classification Hierarchy
- **Cardinality:** 1:1
- **Description:** Each compound has one primary classification path
- **Key Fields:** inchikey, kingdom through direct_parent

### Compound to Alternative Parents
- **Cardinality:** 1:N
- **Description:** Compounds may fit multiple classification paths
- **Key Fields:** inchikey, alternative_parents

### Classification Level to Parent
- **Cardinality:** N:1
- **Description:** Each level has one parent in primary hierarchy
- **Key Fields:** chemont_id

### Classification Level to Children
- **Cardinality:** 1:N
- **Description:** Each level may have multiple children
- **Key Fields:** chemont_id

### Compound to Substituents
- **Cardinality:** 1:N
- **Description:** Compounds have multiple functional groups
- **Key Fields:** inchikey, substituents

### Compound to Predicted Terms
- **Cardinality:** 1:N
- **Description:** Mapped to ChEBI/LIPID MAPS
- **Key Fields:** inchikey, predicted terms

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ChemOnt | Chemical Ontology | ClassyFire taxonomy system |
| CHEMONTID | ChemOnt Identifier | Taxonomy node ID |
| InChI | International Chemical Identifier | IUPAC structure standard |
| InChIKey | InChI Key | Hashed structure ID |
| SMILES | Simplified Molecular-Input Line-Entry System | Line notation |
| SDF | Structure-Data File | Chemical file format |
| REST | Representational State Transfer | API architecture |
| JSON | JavaScript Object Notation | Data format |
| XML | Extensible Markup Language | Alternative format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| ChEBI | CHEBI ID | Predicted mapping |
| LIPID MAPS | LM ID | Predicted mapping |
| PubChem | CID (via InChIKey) | Structure lookup |
| NPClassifier | Classification | Natural product taxonomy |

---

## Data Quality Notes

1. **Comprehensive Taxonomy:** 4,825+ chemical categories
2. **Rule-Based:** 7,000+ structural classification rules
3. **Pre-Classified:** 80M+ compounds with InChIKey lookup
4. **Multiple Paths:** Alternative parents capture structural diversity
5. **Hierarchical:** 6-level classification depth
6. **REST API:** Asynchronous classification with batch support
7. **Open Access:** Free for academic and research use
8. **Version Tracking:** Classification version recorded with results

