# Data Source Documentation Template

This template provides a standardized format for documenting data sources in the Gene Platform catalog. It is designed for **human readers** (researchers, data engineers, integration specialists) rather than AI optimization.

---

## Template Usage

1. Copy this template to create a new data source document
2. Replace all `{{PLACEHOLDER}}` variables with actual values
3. Remove sections that don't apply (mark as "N/A" if partially applicable)
4. Add additional sections as needed for source-specific details
5. Save in the appropriate directory under `docs/data/source.new/`

---

## File Naming Convention

```
{{source-id}}-schema.md        # For schema documentation
{{source-id}}.md               # For general source documentation
{{source-id}}.yaml             # For machine-readable metadata
```

Examples: `clinvar-schema.md`, `foodb.md`, `chembl.yaml`

---

## Template Starts Below

Copy everything from the `---` line below to use as your template:

---

```markdown
---
id: {{SOURCE_ID}}
title: "{{SOURCE_NAME}} - Schema Documentation"
type: schema
parent: _index.md
last_updated: {{YYYY-MM-DD}}
status: {{draft|active|deprecated|migrated}}
tags: [{{tag1}}, {{tag2}}, {{tag3}}]
tier: {{1|2|3}}
category: {{genetics|traditional|nutrition|shared}}
subcategory: {{optional-subcategory}}
---

**Parent:** [Schema Documentation](./_index.md)

# {{SOURCE_NAME}}

**Document ID:** {{DOCUMENT_ID}}
**Status:** {{Draft|Final|Review}}
**Owner:** {{Team/Person}}
**Last Updated:** {{Month Year}}
**Version:** {{X.Y}}
**Source Version:** {{Database version if applicable}}

---

## 1. Overview

{{Brief description of the data source - 2-3 paragraphs covering:
- What the database/source is
- Who maintains it
- Primary use cases
- Why it's relevant to the Gene Platform}}

**Primary URL:** {{https://example.com}}
**Organization:** {{Maintaining organization}}
**Update Frequency:** {{Daily|Weekly|Monthly|Quarterly|Annually|Ad-hoc}}

---

## 2. Data Content

### Summary Statistics

| Metric | Value | Notes |
|--------|-------|-------|
| **Total Records** | {{X,XXX,XXX}} | {{As of date}} |
| **{{Entity Type 1}}** | {{count}} | {{description}} |
| **{{Entity Type 2}}** | {{count}} | {{description}} |
| **{{Entity Type 3}}** | {{count}} | {{description}} |
| **Date Range** | {{YYYY - YYYY}} | {{Coverage period}} |
| **Last Release** | {{YYYY-MM-DD}} | {{Version/release name}} |

### Data Categories

| Category | Description | Record Count |
|----------|-------------|--------------|
| {{Category 1}} | {{What it contains}} | {{count}} |
| {{Category 2}} | {{What it contains}} | {{count}} |
| {{Category 3}} | {{What it contains}} | {{count}} |

### Coverage

{{Describe what populations, organisms, conditions, or domains the data covers}}

- **Geographic Coverage:** {{Global|Regional|Country-specific}}
- **Temporal Coverage:** {{Date range of data}}
- **Taxonomic Coverage:** {{Species covered, if applicable}}

---

## 3. Access Methods

### 3.1 Web Interface

| Feature | URL | Description |
|---------|-----|-------------|
| **Browse** | {{https://example.com/browse}} | {{Browse records}} |
| **Search** | {{https://example.com/search}} | {{Search interface}} |
| **Download** | {{https://example.com/download}} | {{Download portal}} |

### 3.2 Programmatic Access

#### REST API

**Base URL:** `{{https://api.example.com/v1/}}`

| Endpoint | Method | Description | Auth Required |
|----------|--------|-------------|---------------|
| `/{{endpoint1}}` | GET | {{Description}} | {{Yes|No}} |
| `/{{endpoint2}}` | GET | {{Description}} | {{Yes|No}} |
| `/{{endpoint3}}` | POST | {{Description}} | {{Yes|No}} |

**Authentication:**
- **Type:** {{API Key|OAuth|Basic Auth|None}}
- **How to obtain:** {{Registration URL or instructions}}

**Rate Limits:**
- {{X requests per second/minute/hour}}
- {{Batch request limits}}

#### Example API Request

```bash
# {{Description of what this request does}}
curl -X GET "{{https://api.example.com/v1/endpoint}}" \
  -H "Authorization: Bearer {{YOUR_API_KEY}}" \
  -H "Accept: application/json"
```

**Example Response:**
```json
{
  "{{field1}}": "{{value1}}",
  "{{field2}}": {{value2}},
  "{{field3}}": ["{{value3a}}", "{{value3b}}"]
}
```

### 3.3 Bulk Download

**Download Portal:** {{https://example.com/downloads}}

| File | Format | Size | Description |
|------|--------|------|-------------|
| `{{filename1}}` | {{format}} | {{X.X GB}} | {{Description}} |
| `{{filename2}}` | {{format}} | {{X.X MB}} | {{Description}} |
| `{{filename3}}` | {{format}} | {{X.X MB}} | {{Description}} |

**Download Commands:**

```bash
# Download primary data file
wget {{https://example.com/downloads/filename.gz}}

# Download with authentication (if required)
curl -u "{{username}}:{{password}}" -O {{https://example.com/downloads/filename.gz}}

# Using rsync (if available)
rsync -avz {{rsync://example.com/data/}} ./local_directory/
```

### 3.4 Database Connection

{{If direct database access is available}}

| Parameter | Value |
|-----------|-------|
| **Protocol** | {{MySQL|PostgreSQL|MongoDB|etc.}} |
| **Host** | {{hostname}} |
| **Port** | {{port}} |
| **Database** | {{database_name}} |
| **Access** | {{Public|Registration required|Restricted}} |

---

## 4. Schema / Fields

### 4.1 Core Entities

```
{{ENTITY_1}}
    |
    | ({{relationship}})
    v
{{ENTITY_2}} ----> {{ENTITY_3}}
    |
    | ({{relationship}})
    v
{{ENTITY_4}}
```

### 4.2 Primary Table/Collection: `{{table_name}}`

**Description:** {{What this table contains}}

| Field | Type | Required | Description | Example |
|-------|------|----------|-------------|---------|
| `{{field1}}` | {{INTEGER|VARCHAR|TEXT|etc.}} | {{Yes|No}} | {{Description}} | `{{example}}` |
| `{{field2}}` | {{type}} | {{Yes|No}} | {{Description}} | `{{example}}` |
| `{{field3}}` | {{type}} | {{Yes|No}} | {{Description}} | `{{example}}` |
| `{{field4}}` | {{type}} | {{Yes|No}} | {{Description}} | `{{example}}` |
| `{{field5}}` | {{type}} | {{Yes|No}} | {{Description}} | `{{example}}` |

**Primary Key:** `{{field_name}}`
**Indexes:** `{{index1}}`, `{{index2}}`

### 4.3 Secondary Table: `{{table_name_2}}`

**Description:** {{What this table contains}}

| Field | Type | Required | Description | Example |
|-------|------|----------|-------------|---------|
| `{{field1}}` | {{type}} | {{Yes|No}} | {{Description}} | `{{example}}` |
| `{{field2}}` | {{type}} | {{Yes|No}} | {{Description}} | `{{example}}` |

### 4.4 Relationships

| From Table | From Column | To Table | To Column | Cardinality |
|------------|-------------|----------|-----------|-------------|
| `{{table1}}` | `{{column}}` | `{{table2}}` | `{{column}}` | {{1:1|1:N|N:M}} |
| `{{table1}}` | `{{column}}` | `{{table3}}` | `{{column}}` | {{1:1|1:N|N:M}} |

### 4.5 Enumerated Values

#### {{Field Name}} Values

| Value | Meaning | Notes |
|-------|---------|-------|
| `{{value1}}` | {{Description}} | {{Notes}} |
| `{{value2}}` | {{Description}} | {{Notes}} |
| `{{value3}}` | {{Description}} | {{Notes}} |

---

## 5. Usage Examples

### 5.1 Basic Query

**Goal:** {{What this query achieves}}

```sql
-- {{Description}}
SELECT {{fields}}
FROM {{table}}
WHERE {{conditions}}
ORDER BY {{field}}
LIMIT {{n}};
```

### 5.2 Join Query

**Goal:** {{What this query achieves}}

```sql
SELECT {{fields}}
FROM {{table1}} t1
JOIN {{table2}} t2 ON t1.{{field}} = t2.{{field}}
WHERE {{conditions}};
```

### 5.3 Python Example

```python
import {{library}}

def {{function_name}}({{parameters}}):
    """
    {{Docstring describing what this function does}}

    Args:
        {{param}}: {{Description}}

    Returns:
        {{Return type and description}}
    """
    # {{Implementation}}
    {{code}}
    return {{result}}

# Usage example
result = {{function_name}}({{arguments}})
print(result)
```

### 5.4 API Integration Example

```python
import requests

BASE_URL = "{{https://api.example.com/v1}}"
API_KEY = "{{your_api_key}}"

def fetch_{{entity}}({{parameters}}):
    """Fetch {{entity}} from {{source_name}}."""
    response = requests.get(
        f"{BASE_URL}/{{endpoint}}",
        headers={"Authorization": f"Bearer {API_KEY}"},
        params={{params_dict}}
    )
    response.raise_for_status()
    return response.json()

# Example usage
data = fetch_{{entity}}({{arguments}})
```

---

## 6. Integration

### 6.1 Cross-References

This data source links to other databases via the following identifiers:

| External Database | Identifier Type | Field/Column | Coverage |
|-------------------|-----------------|--------------|----------|
| {{Database 1}} | {{ID type}} | `{{field}}` | {{X%}} |
| {{Database 2}} | {{ID type}} | `{{field}}` | {{X%}} |
| {{Database 3}} | {{ID type}} | `{{field}}` | {{X%}} |

### 6.2 Hub Identifiers

| Entity Type | Hub Identifier | This Source's Field | Mapping Notes |
|-------------|----------------|---------------------|---------------|
| Proteins | UniProt ID | `{{field}}` | {{Direct/Requires mapping}} |
| Genes | NCBI Gene ID | `{{field}}` | {{Direct/Requires mapping}} |
| Compounds | InChIKey | `{{field}}` | {{Direct/Requires mapping}} |
| Diseases | MONDO ID | `{{field}}` | {{Direct/Requires mapping}} |

### 6.3 Related Gene Platform Sources

| Source | Relationship | Integration Point |
|--------|--------------|-------------------|
| [{{Related Source 1}}]({{link}}) | {{Description of relationship}} | {{Shared identifier or concept}} |
| [{{Related Source 2}}]({{link}}) | {{Description of relationship}} | {{Shared identifier or concept}} |

### 6.4 Integration Workflow

```
{{SOURCE_NAME}}
    |
    | Extract: {{what fields}}
    v
[Transformation Layer]
    |
    | Map: {{identifier mapping}}
    v
[Gene Platform Knowledge Graph]
    |
    | Link via: {{hub identifiers}}
    v
[Other Data Sources]
```

---

## 7. References

### Primary Documentation

1. **Official Website:** {{https://example.com}}
2. **API Documentation:** {{https://example.com/docs/api}}
3. **Schema Documentation:** {{https://example.com/docs/schema}}
4. **FAQ:** {{https://example.com/faq}}

### Publications

1. {{Author(s)}} ({{Year}}). "{{Title}}." *{{Journal}}*, {{Volume}}({{Issue}}):{{Pages}}. DOI: {{doi}}

2. {{Author(s)}} ({{Year}}). "{{Title}}." *{{Journal}}*, {{Volume}}({{Issue}}):{{Pages}}. DOI: {{doi}}

### Related Resources

- {{Resource 1}}: {{URL}}
- {{Resource 2}}: {{URL}}

---

## 8. License & Terms

### Data License

| Aspect | Details |
|--------|---------|
| **License Type** | {{CC BY|CC BY-SA|CC BY-NC|Proprietary|Custom}} |
| **License URL** | {{https://example.com/license}} |
| **Commercial Use** | {{Allowed|Restricted|Prohibited}} |
| **Attribution Required** | {{Yes|No}} |
| **Share-Alike** | {{Yes|No}} |
| **Modification Allowed** | {{Yes|No}} |

### Terms of Use

- {{Key term 1}}
- {{Key term 2}}
- {{Restrictions or requirements}}

### Citation Requirements

**Preferred Citation:**
```
{{Citation format required by the data source}}
```

**BibTeX:**
```bibtex
@article{{{citation_key}},
  author = {{{Authors}}},
  title = {{{Title}}},
  journal = {{{Journal}}},
  year = {{{Year}}},
  volume = {{{Volume}}},
  pages = {{{Pages}}},
  doi = {{{DOI}}}
}
```

---

## 9. Update History

### Version History

| Version | Date | Changes |
|---------|------|---------|
| {{X.Y}} | {{YYYY-MM-DD}} | {{Description of changes}} |
| {{X.Y-1}} | {{YYYY-MM-DD}} | {{Description of changes}} |
| {{X.Y-2}} | {{YYYY-MM-DD}} | {{Description of changes}} |

### Release Schedule

- **Update Frequency:** {{Daily|Weekly|Monthly|Quarterly|Annually}}
- **Next Expected Release:** {{Date or "TBD"}}
- **Release Notes URL:** {{https://example.com/releases}}

### Change Notifications

- **Mailing List:** {{email or subscription URL}}
- **RSS Feed:** {{feed URL}}
- **Twitter/Social:** {{@handle}}

---

## 10. Data Format

| Aspect | Value |
|--------|-------|
| **Primary Format** | {{VCF|JSON|XML|CSV|TSV|SQL|etc.}} |
| **Alternative Formats** | {{list of other formats}} |
| **Compression** | {{gzip|bzip2|zip|none}} |
| **Encoding** | {{UTF-8|ASCII|ISO-8859-1}} |
| **API Response Format** | {{JSON|XML|Both}} |

### Format-Specific Notes

{{Any important notes about data formats, such as:
- Field delimiters
- Quote handling
- Null value representation
- Multi-value field handling
- Special characters}}

---

## 11. Data Set Size

| Component | Records | Storage Size | Notes |
|-----------|---------|--------------|-------|
| **{{Component 1}}** | {{X,XXX,XXX}} | {{X.X GB}} | {{Notes}} |
| **{{Component 2}}** | {{X,XXX,XXX}} | {{X.X GB}} | {{Notes}} |
| **{{Component 3}}** | {{X,XXX,XXX}} | {{X.X MB}} | {{Notes}} |
| **Total (Compressed)** | - | {{X.X GB}} | {{Compression type}} |
| **Total (Uncompressed)** | - | {{X.X GB}} | - |

### Storage Requirements

| Deployment | Recommended Storage | Notes |
|------------|---------------------|-------|
| **Download** | {{X GB}} | Raw downloaded files |
| **Working** | {{X GB}} | Processing/staging |
| **Production** | {{X GB}} | Final indexed data |

---

## 12. Sample Data

### Example Record

```json
{
  "{{field1}}": "{{value1}}",
  "{{field2}}": {{value2}},
  "{{field3}}": "{{value3}}",
  "{{nested_object}}": {
    "{{subfield1}}": "{{subvalue1}}",
    "{{subfield2}}": {{subvalue2}}
  },
  "{{array_field}}": [
    "{{item1}}",
    "{{item2}}"
  ]
}
```

### Sample Query Results

| {{field1}} | {{field2}} | {{field3}} | {{field4}} |
|------------|------------|------------|------------|
| {{value}} | {{value}} | {{value}} | {{value}} |
| {{value}} | {{value}} | {{value}} | {{value}} |
| {{value}} | {{value}} | {{value}} | {{value}} |

---

## 13. Glossary

### Key Terms

| Term | Definition | Example |
|------|------------|---------|
| `{{term1}}` | {{Definition of this term in context of this data source}} | `{{example value}}` |
| `{{term2}}` | {{Definition}} | `{{example value}}` |
| `{{term3}}` | {{Definition}} | `{{example value}}` |
| `{{term4}}` | {{Definition}} | `{{example value}}` |
| `{{term5}}` | {{Definition}} | `{{example value}}` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| {{Domain term 1}} | {{Definition in context of this domain}} | {{Related field or concept}} |
| {{Domain term 2}} | {{Definition}} | {{Related field or concept}} |
| {{Domain term 3}} | {{Definition}} | {{Related field or concept}} |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| {{ACRONYM1}} | {{Full expansion}} | {{Additional context}} |
| {{ACRONYM2}} | {{Full expansion}} | {{Additional context}} |
| {{ACRONYM3}} | {{Full expansion}} | {{Additional context}} |

---

## 14. Troubleshooting & FAQ

### Common Issues

**Q: {{Common question 1}}**
A: {{Answer}}

**Q: {{Common question 2}}**
A: {{Answer}}

**Q: {{Common question 3}}**
A: {{Answer}}

### Known Limitations

- {{Limitation 1}}
- {{Limitation 2}}
- {{Limitation 3}}

### Support

- **Documentation:** {{URL}}
- **Help Desk:** {{email or URL}}
- **Community Forum:** {{URL}}
- **Issue Tracker:** {{URL}}

---

## Change Log (This Document)

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| {{1.0}} | {{YYYY-MM-DD}} | {{Author}} | {{Initial documentation}} |

---

*Document generated from template. Last template update: January 2026.*
```

---

## Section Requirements

### Required Sections (must include)

1. **Overview** - What the source is
2. **Data Content** - What data it contains, record counts
3. **Access Methods** - How to get the data
4. **Schema/Fields** - Data structure
5. **License & Terms** - Usage rights
6. **Data Set Size** - Storage estimates
7. **Data Format** - File formats
8. **Sample Data** - Example records
9. **Glossary** - Term definitions

### Recommended Sections (include when applicable)

10. **Usage Examples** - Query/code examples
11. **Integration** - Cross-references to other sources
12. **References** - Citations and documentation links
13. **Update History** - Version and release information
14. **Troubleshooting & FAQ** - Common questions

---

## Best Practices

### Writing Style

1. **Be specific** - Use exact values, not approximations ("24.3M records" not "millions of records")
2. **Be current** - Include dates for statistics ("As of January 2026")
3. **Be practical** - Include working code examples
4. **Be complete** - Document all fields and values
5. **Be consistent** - Follow the format for all sources

### Tables

- Use tables for structured information
- Include headers with clear labels
- Align columns for readability
- Keep rows concise

### Code Examples

- Test all code before including
- Include necessary imports
- Add comments explaining key steps
- Show expected output

### Cross-References

- Link to related documentation within the catalog
- Use relative paths for internal links
- Verify all links work

---

## Frontmatter Reference

```yaml
---
id: unique-identifier           # Required: lowercase-with-dashes
title: "Human Readable Title"   # Required: quoted string
type: schema                    # Required: schema|overview|guide
parent: _index.md               # Required: parent document
last_updated: YYYY-MM-DD        # Required: ISO date
status: draft|active|deprecated # Required: document status
tags: [tag1, tag2]              # Required: array of tags
tier: 1|2|3                     # Optional: implementation priority
category: genetics|traditional|nutrition|shared  # Required: data category
subcategory: tcm|ayurveda       # Optional: for traditional medicine only
---
```

---

## Directory Placement

| Category | Directory | Example |
|----------|-----------|---------|
| Genetics | `databases/genetics/` | `dbsnp-schema.md` |
| Traditional Medicine | `databases/traditional/` | `batman-tcm-schema.md` |
| Nutrition | `databases/nutrition/` | `foodb-schema.md` |
| Pathways | `databases/pathways/` | `reactome-schema.md` |
| Literature | `databases/literature/` | `pubmed-schema.md` |
| Compounds | `databases/compounds/` | `chembl-schema.md` |
| All Schemas | `operations/schemas/` | `*-schema.md` |

---

*Template maintained by Gene Platform Data Team. For questions, contact: data-team@geneplatform.org*
