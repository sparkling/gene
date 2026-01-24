# 5.3 Western/Global Herbal - Data Dictionary

## Overview

This data dictionary documents the unified schema for Western and global herbal medicine data, integrating information from NAPRALERT (literature-based natural products) and EMA Herbal (European regulatory monographs).

**Subcategory ID:** 5.3
**Subcategory Name:** Western/Global Herbal
**Data Sources:** NAPRALERT, EMA Herbal

---

## Data Sources Summary

| Database | Focus Area | Key Contributions |
|----------|------------|-------------------|
| NAPRALERT | Natural Products Alert Database | Literature-curated ethnobotany, pharmacology, chemistry, toxicology |
| EMA Herbal | European Medicines Agency Herbal Monographs | Regulatory-approved indications, safety data, dosing |

---

## Unified Fields

### Organism Identification

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| organism_id | string | 1:1 | Yes | Unique identifier for biological organism/plant | NAPRALERT | - |
| scientific_name | string | 1:1 | No | Binomial scientific name | Both | - |
| family | string | 1:1 | No | Botanical family | Both | - |
| common_names | array[object] | 1:N | No | Vernacular names by language/region | Both | {"name": "...", "language": "..."} |
| plant_part | string | 1:1 | No | Plant part used (radix, folium, flos, etc.) | Both | - |

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| compound_id | string | 1:1 | No | Unique compound identifier | NAPRALERT | - |
| compound_name | string | 1:1 | No | Chemical compound name | NAPRALERT | - |
| cas_number | string | 1:1 | No | CAS registry number | NAPRALERT | - |
| molecular_formula | string | 1:1 | No | Chemical formula | NAPRALERT | - |
| molecular_weight | float | 1:1 | No | Molecular weight in Daltons | NAPRALERT | - |

---

## Source-Specific Fields

### NAPRALERT Database

#### Record Identification

| Field Name | Data Type | Cardinality | Description | Examples/Values |
|------------|-----------|-------------|-------------|-----------------|
| napralert_record_id | string | 1:1 | NAPRALERT citation record identifier | - |
| napralert_kingdom | enum | 1:1 | Biological kingdom of source organism | Plantae, Fungi, Bacteria, Animalia |
| napralert_geographic_distribution | array[string] | 1:N | Native geographic regions | - |
| napralert_compound_class | string | 1:1 | Chemical class classification | Diarylheptanoid, Alkaloid, Flavonoid |
| napralert_pubmed_id | integer | 1:1 | PubMed literature identifier | - |
| napralert_data_types | array[string] | 1:N | Types of data in record | Ethnobotany, Pharmacology, Chemistry, Toxicology |

#### Ethnobotanical Use Records

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| napralert_ethnobotanical_use | array[object] | 1:N | Traditional/indigenous use documentation |

**Ethnobotanical Use Structure:**

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| use_id | string | Unique use record identifier |
| traditional_use | string | Description of traditional use |
| preparation_method | string | How the plant is prepared |
| administration_route | string | Route of administration |
| dosage | string | Traditional dosage |
| geographic_origin | string | Region where use is documented |
| cultural_group | string | Cultural/ethnic group using the plant |
| plant_part | string | Plant part used |
| literature_ref | string | Literature citation |

#### Bioactivity Records

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| napralert_bioactivity | array[object] | 1:N | Pharmacological activity records |

**Bioactivity Structure:**

| Sub-field | Data Type | Description | Values |
|-----------|-----------|-------------|--------|
| activity_id | string | Unique activity record identifier | - |
| activity_type | string | Type of pharmacological activity | - |
| target | string | Biological target | - |
| test_system | string | Experimental test system | - |
| result | string | Activity result description | - |
| dose | string | Dose tested | - |
| unit | string | Dose unit | - |
| effect_type | enum | Classification of effect | Active, Inactive, Moderate |
| literature_ref | string | Literature citation | - |

#### Chemistry Records

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| napralert_chemistry_record | array[object] | 1:N | Chemical isolation/identification data |

**Chemistry Record Structure:**

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| chemistry_id | string | Unique chemistry record identifier |
| isolation_method | string | Method used to isolate compound |
| yield | string | Isolation yield |
| identification_method | string | Method used to identify compound |
| literature_ref | string | Literature citation |

#### Toxicology Records

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| napralert_toxicology_record | array[object] | 1:N | Safety and toxicity data |

**Toxicology Record Structure:**

| Sub-field | Data Type | Description | Values |
|-----------|-----------|-------------|--------|
| toxicology_id | string | Unique toxicology record identifier | - |
| toxicity_type | enum | Type of toxicity study | Acute, Chronic, Subacute, Reproductive, Developmental |
| test_system | string | Animal model or test system | - |
| endpoint | string | Toxicological endpoint | LD50, LC50, NOAEL, LOAEL |
| value | float | Numeric endpoint value | - |
| unit | string | Unit for endpoint value | - |
| route | string | Route of administration | - |
| adverse_effects | array[string] | Observed adverse effects | - |
| literature_ref | string | Literature citation | - |

---

### EMA Herbal Database

#### Monograph Identification

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| ema_substance_id | string | 1:1 | EMA herbal substance reference number | - |
| ema_latin_name | string | 1:1 | Botanical Latin name with plant part | - |
| ema_pharmacopoeia_reference | string | 1:1 | European Pharmacopoeia monograph reference | Ph. Eur. monograph number |
| ema_monograph_id | string | 1:1 | HMPC monograph identifier | EMA/HMPC/xxxxxx/xxxx |
| ema_adoption_date | date | 1:1 | HMPC monograph adoption date | - |
| ema_revision_date | date | 1:1 | Latest revision date | - |
| ema_status | enum | 1:1 | Monograph status | Final, Draft, Under revision |

#### Use Classification

| Field Name | Data Type | Cardinality | Description | Values |
|------------|-----------|-------------|-------------|--------|
| ema_use_category | enum | 1:1 | Regulatory use classification | Well-established, Traditional |
| ema_therapeutic_indication | array[string] | 1:N | Approved therapeutic indications | - |

#### Preparation Specifications

| Field Name | Data Type | Cardinality | Description | Values |
|------------|-----------|-------------|-------------|--------|
| ema_preparation_type | enum | 1:1 | Type of herbal preparation | Herbal substance, Powdered, Dry extract, Soft extract, Liquid extract, Tincture, Expressed juice, Essential oil |
| ema_drug_extract_ratio | string | 1:1 | DER specification for extracts | 4-7:1, 1:1, 1:5 |
| ema_standardization | string | 1:1 | Marker compound standardization | Valerenic acids 0.1-0.8%, Hypericin 0.1-0.3% |
| ema_preparation | array[object] | 1:N | Detailed preparation specifications | See structure below |

**Preparation Structure:**

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| preparation_type | string | Type of preparation |
| description | string | Preparation description |
| extraction_solvent | string | Solvent used for extraction |
| drug_extract_ratio | string | DER specification |
| standardization | string | Standardization requirements |

#### Posology (Dosing)

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| ema_posology | object | 1:1 | Dosage and administration details |

**Posology Structure:**

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| adults | string | Adult dosing |
| children | string | Pediatric dosing |
| elderly | string | Geriatric dosing |
| duration | string | Duration of use |

#### Safety Information

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| ema_contraindications | array[object] | 1:N | Conditions where use is contraindicated |
| ema_warnings | array[string] | 1:N | Special warnings and precautions |
| ema_interactions | array[object] | 1:N | Drug-herb interactions |
| ema_pregnancy_data | object | 1:1 | Pregnancy and lactation safety |
| ema_adverse_effects | array[object] | 1:N | Known adverse reactions |

**Contraindications Structure:**

| Sub-field | Data Type | Description | Values |
|-----------|-----------|-------------|--------|
| condition | string | Contraindicated condition | - |
| severity | enum | Severity of contraindication | Absolute, Relative |
| evidence_level | enum | Evidence supporting contraindication | Clinical, Preclinical, Theoretical |

**Interactions Structure:**

| Sub-field | Data Type | Description | Values |
|-----------|-----------|-------------|--------|
| interacting_substance | string | Drug or substance that interacts | - |
| effect | string | Effect of interaction | - |
| severity | enum | Severity of interaction | Major, Moderate, Minor |
| mechanism | string | Mechanism of interaction | - |
| evidence | string | Evidence supporting interaction | - |

**Pregnancy Data Structure:**

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| pregnancy_category | string | Pregnancy safety category |
| lactation_category | string | Lactation safety category |
| recommendation | string | Specific recommendation |

**Adverse Effects Structure:**

| Sub-field | Data Type | Description | Values |
|-----------|-----------|-------------|--------|
| effect | string | Adverse effect description | - |
| frequency | enum | Frequency of occurrence | Very common, Common, Uncommon, Rare, Very rare |
| meddra_term | string | MedDRA preferred term | - |
| system_organ_class | string | MedDRA system organ class | - |

#### Community List Entry

| Field Name | Data Type | Cardinality | Description |
|------------|-----------|-------------|-------------|
| ema_community_list_entry | object | 1:1 | EU Community list harmonized entry |

**Community List Entry Structure:**

| Sub-field | Data Type | Description |
|-----------|-----------|-------------|
| entry_id | string | Community list entry identifier |
| herbal_preparations | array[object] | Approved preparations |
| indication | string | Approved indication |
| route | string | Route of administration |
| posology | string | Dosing information |
| period_of_use | string | Maximum period of use |
| restrictions | array[string] | Use restrictions |

---

## Source Field Mappings

### NAPRALERT Mappings
| Original Field | Unified Field |
|----------------|---------------|
| organism_id | organism_id |
| scientific_name | scientific_name |
| family | family |
| common_names | common_names |
| compound_id | compound_id |
| compound_name | compound_name |
| cas_number | cas_number |
| molecular_formula | molecular_formula |
| molecular_weight | molecular_weight |
| plant_part | plant_part |
| record_id | napralert_record_id |
| kingdom | napralert_kingdom |
| geographic_distribution | napralert_geographic_distribution |
| ethnobotanical_use | napralert_ethnobotanical_use |
| bioactivity | napralert_bioactivity |
| chemistry_record | napralert_chemistry_record |
| toxicology_record | napralert_toxicology_record |
| compound_class | napralert_compound_class |
| pubmed_id | napralert_pubmed_id |
| data_types | napralert_data_types |

### EMA Herbal Mappings
| Original Field | Unified Field |
|----------------|---------------|
| scientific_name | scientific_name |
| family | family |
| common_names | common_names |
| plant_part | plant_part |
| substance_id | ema_substance_id |
| latin_name | ema_latin_name |
| pharmacopoeia_reference | ema_pharmacopoeia_reference |
| monograph_id | ema_monograph_id |
| adoption_date | ema_adoption_date |
| revision_date | ema_revision_date |
| status | ema_status |
| use_category | ema_use_category |
| therapeutic_indication | ema_therapeutic_indication |
| posology | ema_posology |
| preparation | ema_preparation |
| preparation_type | ema_preparation_type |
| drug_extract_ratio | ema_drug_extract_ratio |
| standardization | ema_standardization |
| contraindications | ema_contraindications |
| warnings | ema_warnings |
| interactions | ema_interactions |
| pregnancy_data | ema_pregnancy_data |
| adverse_effects | ema_adverse_effects |
| community_list_entry | ema_community_list_entry |

---

## Metadata Field

| Field Name | Data Type | Description |
|------------|-----------|-------------|
| _source.database | string | Source database name |
| _source.version | string | Database version |
| _source.access_date | date | Date of data access |
| _source.original_id | string | Original identifier in source database |

---

## Data Integration Notes

### NAPRALERT vs EMA Herbal Comparison

| Aspect | NAPRALERT | EMA Herbal |
|--------|-----------|------------|
| Focus | Literature-curated research data | Regulatory-approved clinical use |
| Coverage | Global natural products (plants, fungi, bacteria, animals) | European herbal medicines (plants only) |
| Evidence Type | Primary research literature | Systematic reviews, clinical evidence |
| Compound Data | Yes (chemistry records) | No (focuses on whole plant preparations) |
| Safety Data | Toxicology from literature | Regulatory safety assessments |
| Indications | Traditional and pharmacological | Approved therapeutic use |
| Geographic Scope | Global | European Union |

### Integration Opportunities

1. **Botanical Name Linking:** Both sources use scientific names for cross-referencing
2. **Safety Data Enrichment:** NAPRALERT toxicology can supplement EMA safety assessments
3. **Traditional Use Validation:** NAPRALERT ethnobotanical data supports EMA traditional use claims
4. **Compound Identification:** NAPRALERT provides chemical constituents for EMA-listed herbs
