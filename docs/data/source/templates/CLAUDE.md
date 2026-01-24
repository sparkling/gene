# Data Source Templates - AI Guide

Instructions for AI creating new data source documentation.

## Creating a New Data Source

### Step 1: Create Folder

```
docs/data/source/resource/{category}/{subcategory}/{source}/
```

### Step 2: Generate Required Files (7)

| File | Template | Purpose |
|------|----------|---------|
| `README.md` | [readme.template.md](./readme.template.md) | Overview + frontmatter |
| `schema.json` | [schema.template.json](./schema.template.json) | Data structure |
| `dictionary.md` | [dictionary.template.md](./dictionary.template.md) | Field semantics |
| `sample.json` | [sample.template.json](./sample.template.json) | Example records |
| `download.md` | [download.template.md](./download.template.md) | Access methods |
| `license.md` | [license.template.md](./license.template.md) | License terms |
| `mapping.xslt` | [mapping.template.xslt](./mapping.template.xslt) | Transform to parent |

### Step 3: Optional Files

| File | Template | When Needed |
|------|----------|-------------|
| `xrefs.md` | [xrefs.template.md](./xrefs.template.md) | Cross-references exist |

## Download Sub-Templates

Compose `download.md` by including applicable sub-templates:

| Access Method | Include Template |
|---------------|------------------|
| REST API | [download/rest-api.template.md](./download/rest-api.template.md) |
| SPARQL | [download/sparql.template.md](./download/sparql.template.md) |
| GraphQL | [download/graphql.template.md](./download/graphql.template.md) |
| FTP | [download/ftp.template.md](./download/ftp.template.md) |
| AWS S3 | [download/s3.template.md](./download/s3.template.md) |
| Rsync | [download/rsync.template.md](./download/rsync.template.md) |
| Python | [download/python.template.md](./download/python.template.md) |
| R | [download/r.template.md](./download/r.template.md) |
| Auth Required | [download/auth.template.md](./download/auth.template.md) |

## Template Workflow

1. **Read the base template** from this directory
2. **Replace all `{{PLACEHOLDER}}` values** with actual content
3. **Delete inapplicable sections** (marked with comments)
4. **Include relevant sub-templates** for download.md
5. **Validate against file size targets** (see README.md)

## File Size Targets

| File | Target | Max |
|------|--------|-----|
| README.md | 50-80 | 100 |
| schema.json | 100-150 | 300 |
| dictionary.md | 80-150 | 200 |
| sample.json | 30-50 | 100 |
| download.md | 100-200 | 300 |
| license.md | 30-50 | 80 |
| mapping.xslt | 100-200 | 300 |

## Key Principles

1. **7 required files** - README, schema, dictionary, sample, download, license, mapping
2. **xrefs.md is optional** - Only if cross-references exist
3. **license.md is ALWAYS required** - Not optional
4. **mapping.xslt uses XSLT 3.0** - Native JSON support, complex transforms
5. **Schemas are standalone** - No parent awareness in schema.json
6. **Mappings define relationships** - XSLT transforms to `../schema.json`

## Mapping Chain

```
category/schema.json              # Root unified schema
    ↑
subcategory/mapping.xslt          # Transforms subcategory → category
subcategory/schema.json           # Standalone subcategory schema
    ↑
source/mapping.xslt               # Transforms source → subcategory
source/schema.json                # Standalone source schema
```

## Quick Checklist

- [ ] Folder at correct hierarchy level
- [ ] README.md with YAML frontmatter
- [ ] schema.json with all fields
- [ ] dictionary.md with field semantics
- [ ] sample.json with 2-5 real records
- [ ] download.md with versions and access methods
- [ ] license.md with license terms
- [ ] mapping.xslt with XSLT 3.0 transformation
- [ ] xrefs.md if cross-references exist
