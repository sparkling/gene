# Category 05: Traditional Medicine - Data Dictionary

## Overview

This data dictionary documents the unified schema for Category 05: Traditional Medicine, which integrates data from 13 distinct databases across four subcategories covering global traditional medicine systems.

**Category ID:** 05
**Category Name:** Traditional Medicine
**Total Data Sources:** 13
**Total Unified Fields:** 87
**Total Source-Specific Fields:** 156
**Extraction Date:** 2026-01-24

---

## Subcategories

| ID | Name | Data Sources | Description |
|----|------|--------------|-------------|
| 5.1 | Traditional Chinese Medicine | HERB, ETCM, TCMBank, TCMSID, SymMap, BATMAN-TCM | TCM herbs, formulas, compounds, and targets |
| 5.2 | South/East Asian Systems | NPACT, TM-MC, IMPPAT, KampoDB | Ayurveda, Kampo, and multi-system Asian medicine |
| 5.3 | Western/Global Herbal | NAPRALERT, EMA Herbal | Western herbal medicine and regulatory data |
| 5.4 | Multi-System Integration | HIT 2.0 | Cross-system validated herb-target interactions |

---

## Cross-Subcategory Common Fields

These fields appear across multiple subcategories and enable data integration.

| Field Name | Description | Present In | Semantic Importance |
|------------|-------------|------------|---------------------|
| pubchem_cid | PubChem Compound Identifier | 5.1, 5.2, 5.4 | Primary identifier for compound structure standardization |
| smiles | SMILES molecular notation | 5.1, 5.2, 5.3 | Enables structure-based searching and comparison |
| gene_symbol | HGNC gene symbol | 5.1, 5.2, 5.4 | Primary identifier for target protein gene |
| uniprot_id | UniProt protein accession | 5.1, 5.2, 5.4 | Cross-reference for protein sequence and annotation |
| botanical_name | Latin binomial scientific name | 5.1, 5.2, 5.3 | Taxonomic standardization for plant identification |
| molecular_weight | Molecular weight in Daltons | 5.1, 5.2, 5.3, 5.4 | Universal physical property for compound characterization |
| molecular_formula | Chemical elemental composition | 5.1, 5.2, 5.3, 5.4 | Basic chemical identity information |
| plant_part | Anatomical part of plant | 5.1, 5.2, 5.3, 5.4 | Tissue-level localization of compounds and therapeutic uses |
| traditional_system | Traditional medicine system | 5.2, 5.4 | Cultural/geographic origin of traditional knowledge |

---

## Entity Relationships

### Herb-Compound Relationship
- **Cardinality:** N:M (many-to-many)
- **Description:** Many herbs contain many compounds; many compounds found in multiple herbs
- **Key Fields:** herb_id, compound_id, source_herbs
- **Present In:** 5.1, 5.2

### Compound-Target Relationship
- **Cardinality:** N:M (many-to-many)
- **Description:** Many compounds interact with many targets; many targets have multiple ligands
- **Key Fields:** compound_id, target_id, predicted_targets, validated_targets
- **Present In:** 5.1, 5.2, 5.4

### Formula-Herb Relationship
- **Cardinality:** N:M (many-to-many)
- **Description:** Formulas contain multiple herbs; herbs appear in multiple formulas
- **Key Fields:** formula_id, herb_id, herb_role
- **Present In:** 5.1, 5.2

### Symptom-Herb Relationship
- **Cardinality:** N:M (many-to-many)
- **Description:** Symptoms treated by multiple herbs; herbs treat multiple symptoms
- **Key Fields:** symptom_id, herb_id, symptom_herb_associations
- **Present In:** 5.1

### Symptom-Gene Relationship
- **Cardinality:** N:M (many-to-many)
- **Description:** Symptoms associated with multiple genes via phenotype mapping
- **Key Fields:** symptom_id, gene_symbol, hpo_mappings
- **Present In:** 5.1

---

## Traditional Indication Fields

### TCM-Specific Fields

| Field Name | TCM Term | Description | Values |
|------------|----------|-------------|--------|
| tcm_nature | Xing (性) | Five-level thermal classification | cold, cool, neutral, warm, hot |
| tcm_flavor | Wei (味) | Taste-based therapeutic classification | sweet, bitter, sour, salty, pungent, astringent |
| meridian_tropism | Gui Jing (归经) | Target organ/meridian affinity | 12 primary meridians + pericardium |
| tcm_pattern | Zheng (证) | TCM syndrome/pattern diagnosis | Various pattern names |
| herb_role | Jun Chen Zuo Shi (君臣佐使) | Hierarchical role in formula | sovereign, minister, assistant, envoy |

### Ayurvedic Fields

| Field Name | Description | Sub-fields |
|------------|-------------|------------|
| ayurveda_properties | Ayurvedic pharmacological classification | rasa, guna, vipaka, virya, dosha_effect |

### Regulatory Fields

| Field Name | Description | Present In |
|------------|-------------|------------|
| use_category | EMA regulatory classification | 5.3 |
| pharmacopoeia_status | Chinese Pharmacopoeia listing | 5.1 |

---

## Compound-Target Interaction Fields

| Field Name | Description | Metrics/Values | Present In |
|------------|-------------|----------------|------------|
| binding_affinity | Quantitative interaction strength | IC50, Ki, Kd, EC50, KD (nM, uM, mM) | 5.1, 5.2, 5.4 |
| interaction_type | Classification of interaction | known/predicted (BATMAN-TCM); inhibitor/agonist/antagonist/modulator/binder/substrate (HIT 2.0) | 5.1, 5.4 |
| confidence_score | Prediction confidence or evidence quality | 0-1 scale | 5.1, 5.2, 5.4 |
| experimental_method | Assay method used for validation | Various assay types | 5.4 |

---

## Data Quality Notes

1. **Cardinality:** Fields marked 1:1 are single-valued; 1:N are multi-valued arrays
2. **N:M Relationships:** Indicate many-to-many associations between entities
3. **Semantic Definitions:** Capture the meaning and intended use of each field
4. **Cross-References:** PubChem, UniProt, KEGG enable multi-database integration
5. **TCM Terminology:** Includes Chinese terms (Pinyin) with definitions
6. **ADMET Fields:** Computational predictions, not experimental values

---

## See Also

- [5.1 Traditional Chinese Medicine Data Dictionary](./5.1.traditional.chinese.medicine/_data-dictionary.md)
- [5.2 South/East Asian Systems Data Dictionary](./5.2.south.east.asian.systems/_data-dictionary.md)
- [5.3 Western/Global Herbal Data Dictionary](./5.3.western.global.herbal/_data-dictionary.md)
- [5.4 Multi-System Integration Data Dictionary](./5.4.multi.system.integration/_data-dictionary.md)
