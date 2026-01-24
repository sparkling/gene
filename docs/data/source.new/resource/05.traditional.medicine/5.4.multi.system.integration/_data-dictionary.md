# 5.4 Multi-System Integration - Data Dictionary

## Overview

This data dictionary documents the unified schema for multi-system traditional medicine integration data from HIT 2.0, which provides experimentally validated herbal ingredient-target interactions across multiple traditional medicine systems.

**Subcategory ID:** 5.4
**Subcategory Name:** Multi-System Integration
**Data Sources:** HIT 2.0 (Herbal Ingredient-Target Database)

---

## Data Source Summary

| Database | Focus Area | Key Contributions |
|----------|------------|-------------------|
| HIT 2.0 | Herbal Ingredient-Target Interactions | Experimentally validated compound-target interactions with binding affinity data across TCM, Ayurveda, Kampo, and Western herbal systems |

**Key Characteristics:**
- Contains ONLY experimentally validated interactions (no computational predictions)
- Cross-references multiple traditional medicine systems
- Includes detailed binding affinity measurements
- Links to primary literature with full citation data

---

## Unified Fields

### Ingredient Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ingredient_id | string | 1:1 | Yes | Unique herbal ingredient identifier | HIT001234 |
| ingredient_name | string | 1:1 | No | Chemical compound name | - |

### Chemical Identifiers

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| pubchem_cid | integer | 1:1 | No | PubChem Compound Identifier | - |
| chembl_id | string | 1:1 | No | ChEMBL compound identifier | CHEMBL12345 |
| cas_number | string | 1:1 | No | CAS registry number | 123-45-6 |
| inchi_key | string | 1:1 | No | InChIKey hash (27 characters) | - |
| smiles | string | 1:1 | No | SMILES notation | - |
| molecular_formula | string | 1:1 | No | Chemical formula | - |
| molecular_weight | float | 1:1 | No | Molecular weight in Daltons | - |

### Target Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| target_id | string | 1:1 | No | Target protein identifier (UniProt or gene symbol) | - |
| gene_symbol | string | 1:1 | No | HGNC gene symbol | EGFR, MTOR, BCL2 |
| uniprot_id | string | 1:1 | No | UniProt protein accession | P00533 |

---

## Source-Specific Fields (HIT 2.0)

### Interaction Record

| Field Name | Data Type | Cardinality | Description | Examples/Values |
|------------|-----------|-------------|-------------|-----------------|
| hit_interaction_id | string | 1:1 | Unique interaction record identifier | HIT_INT_001234 |
| hit_interaction_type | enum | 1:1 | Type of compound-target interaction | inhibitor, agonist, antagonist, modulator, binder, substrate |
| hit_validation_status | enum | 1:1 | Experimental validation status | validated (HIT 2.0 only contains validated) |

### Experimental Method

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| hit_experimental_method | string | 1:1 | Experimental method used to validate interaction | binding assay, enzyme inhibition, cell-based assay |
| hit_assay_type | string | 1:1 | Specific assay type | radioligand binding, fluorescence polarization |

### Binding Affinity Data

| Field Name | Data Type | Cardinality | Description | Values |
|------------|-----------|-------------|-------------|--------|
| hit_binding_affinity_metric | enum | 1:1 | Type of binding affinity measurement | IC50, Ki, Kd, EC50, KD |
| hit_binding_affinity_value | float | 1:1 | Quantitative binding affinity value | Numeric value |
| hit_binding_affinity_unit | string | 1:1 | Unit for binding affinity | nM, uM, mM |
| hit_binding_affinity_relation | enum | 1:1 | Relation symbol for affinity value | =, <, >, ~ |

**Binding Affinity Interpretation:**
- **IC50:** Half maximal inhibitory concentration
- **Ki:** Inhibition constant
- **Kd:** Dissociation constant
- **EC50:** Half maximal effective concentration
- **KD:** Equilibrium dissociation constant

### Confidence Assessment

| Field Name | Data Type | Cardinality | Description | Values |
|------------|-----------|-------------|-------------|--------|
| hit_confidence_level | enum | 1:1 | Confidence level of experimental evidence | high, medium, low |

**Confidence Level Criteria:**
- **High:** Multiple independent validations, quantitative binding data
- **Medium:** Single validation study, semi-quantitative data
- **Low:** Limited evidence, qualitative data

### Literature Citation

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| hit_pubmed_id | integer | 1:1 | PubMed literature reference | 12345678 |
| hit_doi | string | 1:1 | Digital Object Identifier | 10.1000/xyz123 |
| hit_authors | string | 1:1 | Publication authors | Smith J, et al. |
| hit_journal | string | 1:1 | Publication journal | J Med Chem |
| hit_publication_year | integer | 1:1 | Publication year | 2023 |

### Target Classification

| Field Name | Data Type | Cardinality | Description | Values |
|------------|-----------|-------------|-------------|--------|
| hit_protein_class | string | 1:1 | Target protein functional class | Kinase, GPCR, Nuclear Receptor, Enzyme, Ion Channel, Transporter |
| hit_organism | enum | 1:1 | Target organism species | Human, Mouse, Rat |

### Herbal Source Information

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| hit_herbal_source | array[object] | 1:N | Herbal source information |

**Herbal Source Structure:**

| Sub-field | Data Type | Description | Values |
|-----------|-----------|-------------|--------|
| herb_name | string | Name of source herb | - |
| traditional_system | enum | Traditional medicine system | TCM, Ayurveda, Kampo, Western |
| plant_part | string | Part of plant containing compound | root, leaf, bark, etc. |

### Traditional System Classification

| Field Name | Data Type | Cardinality | Description | Values |
|------------|-----------|-------------|-------------|--------|
| hit_traditional_system | array[enum] | 1:N | Traditional medicine systems using this ingredient | TCM, Ayurveda, Kampo, Western, African Traditional |

### Biological Associations

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| hit_pathway_associations | array[string] | 1:N | Associated biological pathways | hsa04151, KEGG pathway IDs |
| hit_disease_associations | array[string] | 1:N | Associated diseases | MESH IDs, OMIM IDs |

---

## Field Mappings

### HIT 2.0 Mappings
| Original Field | Unified Field |
|----------------|---------------|
| ingredient_id | ingredient_id |
| ingredient_name | ingredient_name |
| pubchem_cid | pubchem_cid |
| chembl_id | chembl_id |
| cas_number | cas_number |
| inchi_key | inchi_key |
| smiles | smiles |
| molecular_formula | molecular_formula |
| molecular_weight | molecular_weight |
| target_id | target_id |
| gene_symbol | gene_symbol |
| uniprot_id | uniprot_id |
| interaction_id | hit_interaction_id |
| interaction_type | hit_interaction_type |
| experimental_method | hit_experimental_method |
| assay_type | hit_assay_type |
| binding_affinity_metric | hit_binding_affinity_metric |
| binding_affinity_value | hit_binding_affinity_value |
| binding_affinity_unit | hit_binding_affinity_unit |
| binding_affinity_relation | hit_binding_affinity_relation |
| confidence_level | hit_confidence_level |
| pubmed_id | hit_pubmed_id |
| doi | hit_doi |
| authors | hit_authors |
| journal | hit_journal |
| publication_year | hit_publication_year |
| protein_class | hit_protein_class |
| organism | hit_organism |
| herbal_source | hit_herbal_source |
| traditional_system | hit_traditional_system |
| pathway_associations | hit_pathway_associations |
| disease_associations | hit_disease_associations |
| validation_status | hit_validation_status |

---

## Metadata Field

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| _source.database | string | Source database name (HIT 2.0) |
| _source.version | string | Database version |
| _source.access_date | date | Date of data access |
| _source.original_id | string | Original identifier in source database |

---

## Data Quality and Validation

### Validation Criteria

All interactions in HIT 2.0 meet the following criteria:

1. **Experimental Validation:** Direct experimental evidence (no computational predictions)
2. **Literature Support:** Published in peer-reviewed literature with PubMed citation
3. **Quantitative Data:** Includes binding affinity measurements where available
4. **Species Information:** Target organism clearly specified

### Data Completeness

| Field Category | Completeness | Notes |
|----------------|--------------|-------|
| Ingredient Identifiers | High | Multiple cross-references (PubChem, ChEMBL, CAS) |
| Target Information | High | Gene symbols and UniProt IDs standardized |
| Binding Affinity | Medium | Available for most but not all interactions |
| Literature Citations | High | All entries have PubMed references |
| Traditional System | Medium | Source system identified where documented |

---

## Integration with Other Subcategories

### Cross-Reference Fields

| Field | Links To | Purpose |
|-------|----------|---------|
| pubchem_cid | 5.1, 5.2 | Compound standardization |
| gene_symbol | 5.1, 5.2 | Target identification |
| uniprot_id | 5.1, 5.2 | Protein sequence reference |
| traditional_system | 5.2 | System classification alignment |

### Integration Opportunities

1. **Compound Matching:** Use PubChem CID, InChIKey, or SMILES to link HIT 2.0 ingredients with compounds in TCM (5.1) and South/East Asian (5.2) databases
2. **Target Validation:** HIT 2.0 validated targets can confirm predicted targets from BATMAN-TCM (5.1) or computational studies
3. **Cross-System Analysis:** Compare traditional uses across systems for the same ingredient
4. **Binding Affinity Comparison:** Correlate HIT 2.0 affinity data with docking predictions from KampoDB (5.2)

---

## Example Records

### Complete Interaction Record

```json
{
  "ingredient_id": "HIT001234",
  "ingredient_name": "Curcumin",
  "pubchem_cid": 969516,
  "chembl_id": "CHEMBL116438",
  "cas_number": "458-37-7",
  "inchi_key": "VFLDPWHFBUODDF-FCXRPNKRSA-N",
  "smiles": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
  "molecular_formula": "C21H20O6",
  "molecular_weight": 368.38,
  "target_id": "P04637",
  "gene_symbol": "TP53",
  "uniprot_id": "P04637",
  "hit_interaction_id": "HIT_INT_001234",
  "hit_interaction_type": "modulator",
  "hit_experimental_method": "cell-based assay",
  "hit_assay_type": "reporter gene assay",
  "hit_binding_affinity_metric": "IC50",
  "hit_binding_affinity_value": 2.5,
  "hit_binding_affinity_unit": "uM",
  "hit_binding_affinity_relation": "=",
  "hit_confidence_level": "high",
  "hit_pubmed_id": 12345678,
  "hit_doi": "10.1000/example",
  "hit_authors": "Smith J, et al.",
  "hit_journal": "J Nat Prod",
  "hit_publication_year": 2023,
  "hit_protein_class": "Transcription Factor",
  "hit_organism": "Human",
  "hit_herbal_source": [
    {
      "herb_name": "Curcuma longa",
      "traditional_system": "Ayurveda",
      "plant_part": "rhizome"
    },
    {
      "herb_name": "Jiang Huang",
      "traditional_system": "TCM",
      "plant_part": "rhizome"
    }
  ],
  "hit_traditional_system": ["Ayurveda", "TCM"],
  "hit_pathway_associations": ["hsa04115", "hsa04210"],
  "hit_disease_associations": ["MESH:D009369"],
  "hit_validation_status": "validated",
  "_source": {
    "database": "HIT 2.0",
    "version": "2.0",
    "access_date": "2026-01-24",
    "original_id": "HIT001234"
  }
}
```

---

## Use Cases

### 1. Drug Repurposing
Identify herbal ingredients that target the same proteins as approved drugs for potential therapeutic applications.

### 2. Multi-Target Analysis
Find ingredients that modulate multiple targets within a disease pathway for polypharmacology approaches.

### 3. Traditional Knowledge Validation
Use validated interactions to support traditional therapeutic claims with modern experimental evidence.

### 4. Safety Assessment
Cross-reference ingredient-target interactions with known adverse effect pathways.

### 5. Cross-Cultural Comparison
Analyze how the same ingredient is used across different traditional medicine systems (TCM vs Ayurveda vs Kampo).
