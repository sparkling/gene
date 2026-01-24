# {{LEVEL_NAME}} - Data Dictionary

## Overview

This data dictionary documents the schema for {{LEVEL_NAME}}.

| Property | Value |
|----------|-------|
| **Level** | {{category/subcategory/source}} |
| **ID** | {{LEVEL_ID}} |
| **Name** | {{LEVEL_NAME}} |
| **Data Sources** | {{SOURCE_COUNT}} |
| **Total Fields** | {{FIELD_COUNT}} |
| **Last Updated** | {{DATE}} |

---

## Subcategories

| ID | Name | Data Sources | Description |
|----|------|--------------|-------------|
| {{SUBCAT_ID_1}} | {{SUBCAT_NAME_1}} | {{SOURCES_1}} | {{DESC_1}} |
| {{SUBCAT_ID_2}} | {{SUBCAT_NAME_2}} | {{SOURCES_2}} | {{DESC_2}} |

---

## Unified Fields

### {{FIELD_GROUP_1}}

| Field Name | Data Type | Cardinality | Required | Description | Sources | Examples |
|------------|-----------|-------------|----------|-------------|---------|----------|
| {{FIELD_1}} | {{TYPE}} | {{CARD}} | {{REQ}} | {{DESC}} | {{SOURCES}} | {{EXAMPLES}} |
| {{FIELD_2}} | {{TYPE}} | {{CARD}} | {{REQ}} | {{DESC}} | {{SOURCES}} | {{EXAMPLES}} |

### {{FIELD_GROUP_2}}

| Field Name | Data Type | Cardinality | Required | Description | Sources | Allowed Values |
|------------|-----------|-------------|----------|-------------|---------|----------------|
| {{FIELD_3}} | {{TYPE}} | {{CARD}} | {{REQ}} | {{DESC}} | {{SOURCES}} | {{VALUES}} |
| {{FIELD_4}} | {{TYPE}} | {{CARD}} | {{REQ}} | {{DESC}} | {{SOURCES}} | {{VALUES}} |

---

## Source-Specific Fields

### {{SOURCE_1}} Database

| Field Name | Data Type | Cardinality | Description | Examples |
|------------|-----------|-------------|-------------|----------|
| {{SRC_FIELD_1}} | {{TYPE}} | {{CARD}} | {{DESC}} | {{EXAMPLES}} |
| {{SRC_FIELD_2}} | {{TYPE}} | {{CARD}} | {{DESC}} | {{EXAMPLES}} |

---

## Entity Relationships

### {{RELATIONSHIP_1}}
- **Cardinality:** {{CARDINALITY}}
- **Description:** {{DESCRIPTION}}
- **Key Fields:** {{KEY_FIELDS}}
- **Present In:** {{SOURCES}}

### {{RELATIONSHIP_2}}
- **Cardinality:** {{CARDINALITY}}
- **Description:** {{DESCRIPTION}}
- **Key Fields:** {{KEY_FIELDS}}
- **Present In:** {{SOURCES}}

---

## Cross-Subcategory Common Fields

| Field Name | Description | Present In | Semantic Importance |
|------------|-------------|------------|---------------------|
| {{COMMON_FIELD_1}} | {{DESC}} | {{SUBCATS}} | {{IMPORTANCE}} |
| {{COMMON_FIELD_2}} | {{DESC}} | {{SUBCATS}} | {{IMPORTANCE}} |

---

## Data Quality Notes

1. **Cardinality:** Fields marked 1:1 are single-valued; 1:N are multi-valued arrays
2. **N:M Relationships:** Indicate many-to-many associations between entities
3. **Cross-References:** {{XREF_IDS}} enable multi-database integration
4. **Source-Specific:** Fields prefixed with source name are not unified

---

## See Also

{{SUBCATEGORY_DICTIONARY_LINKS}}
