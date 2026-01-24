# Data Source Documentation Templates

Instructions for AI to create consistent data source documentation.

## Adding a New Data Source

### Step 1: Create Folder

Create the source folder at the appropriate hierarchy level:

```
docs/data/source/resource/{category}/{subcategory}/[{sub-subcategory}/]{source-name}/
```

**Hierarchy depth:** 3+ levels supported
- Category -> Subcategory -> Source (3 levels)
- Category -> Subcategory -> Sub-subcategory -> Source (4+ levels)

**Naming convention:** lowercase, dots for spaces (e.g., `gwas.catalog`, `batman.tcm`)

---

### Step 2: Create Required Files (7 files)

| Order | File | Template | Purpose |
|-------|------|----------|---------|
| 1 | `README.md` | [readme.template.md](./readme.template.md) | Overview + frontmatter |
| 2 | `schema.json` | [schema.template.json](./schema.template.json) | JSON Schema definition |
| 3 | `dictionary.md` | [dictionary.template.md](./dictionary.template.md) | Field semantics for this source |
| 4 | `sample.json` | [sample.template.json](./sample.template.json) | 2-5 real example records |
| 5 | `download.md` | [download.template.md](./download.template.md) | Access methods, API, versions |
| 6 | `license.md` | [license.template.md](./license.template.md) | License terms (always required) |
| 7 | `mapping.xslt` | [mapping.template.xslt](./mapping.template.xslt) | XSLT 3.0 transformation to parent |

---

### Step 3: Create Optional Files

| File | Template | When Required |
|------|----------|---------------|
| `xrefs.md` | [xrefs.template.md](./xrefs.template.md) | Cross-references to other databases exist |

---

### Step 4: Update Parent Dictionary

Add source-specific terms to `{parent}/dictionary.md` at the nearest category or subcategory level.

---

## AI File Selection Guide

| Need | Read This File |
|------|----------------|
| Overview, context, statistics | `README.md` |
| Data structure, field types | `schema.json` |
| Field definitions, semantics | `dictionary.md` |
| Example records for code generation | `sample.json` |
| Download/API instructions, versions | `download.md` |
| License, commercial use | `license.md` |
| Transform to unified schema | `mapping.xslt` |
| Cross-references to other DBs | `xrefs.md` |

---

## File Size Targets

| File | Target Lines | Max Lines |
|------|--------------|-----------|
| `README.md` | 50-80 | 100 |
| `schema.json` | 100-150 | 300 |
| `dictionary.md` | 80-150 | 200 |
| `sample.json` | 30-50 | 100 |
| `download.md` | 100-200 | 300 |
| `license.md` | 30-50 | 80 |
| `mapping.xslt` | 100-200 | 300 |
| `xrefs.md` | 50-100 | 150 |

---

## Schema Mapping Chain

**Key principle:** Schemas are standalone; mappings define relationships.

```
{category}/
├── schema.json          # Unified schema for domain (standalone)
├── dictionary.md        # Field definitions
│
└── {subcategory}/
    ├── schema.json      # Unified schema for subcategory (standalone)
    ├── dictionary.md    # Field definitions
    ├── mapping.xslt     # XSLT 3.0: subcategory → ../schema.json
    │
    └── {source}/
        ├── schema.json  # Source schema (standalone)
        ├── dictionary.md
        ├── mapping.xslt # XSLT 3.0: source → ../schema.json
        └── ...
```

**Mappings always transform UP to `../schema.json` (immediate parent level).**

---

## Download Sub-Templates

For sources with specific access methods, compose `download.md` from:

| Sub-Template | When to Include |
|--------------|-----------------|
| [rest-api.template.md](./download/rest-api.template.md) | Source has REST API |
| [sparql.template.md](./download/sparql.template.md) | Source has SPARQL endpoint |
| [graphql.template.md](./download/graphql.template.md) | Source has GraphQL API |
| [ftp.template.md](./download/ftp.template.md) | Source offers FTP bulk download |
| [s3.template.md](./download/s3.template.md) | Source offers AWS S3 access |
| [rsync.template.md](./download/rsync.template.md) | Source supports rsync mirroring |
| [python.template.md](./download/python.template.md) | Python package available |
| [r.template.md](./download/r.template.md) | R package available |
| [auth.template.md](./download/auth.template.md) | Authentication required |

---

## Quick Checklist

- [ ] Folder created at correct hierarchy level
- [ ] `README.md` - overview with frontmatter
- [ ] `schema.json` - all fields defined
- [ ] `dictionary.md` - field semantics documented
- [ ] `sample.json` - 2-5 real records
- [ ] `download.md` - acquisition instructions with versions
- [ ] `license.md` - license terms documented
- [ ] `mapping.xslt` - XSLT 3.0 transformation to parent schema
- [ ] `xrefs.md` - cross-references (if applicable)
- [ ] Parent `dictionary.md` updated
