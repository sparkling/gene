# Data Source Documentation - AI Guide

Instructions for AI querying existing data sources.

## How to Find Data Sources

```
docs/data/source/resource/{category}/{subcategory}/{source}/
```

## File Routing

| Query Type | Read This File |
|------------|----------------|
| "What is [source]?" | `README.md` |
| "What data is available?" | `README.md` |
| "Show me the schema" | `schema.json` |
| "What does field X mean?" | `dictionary.md` |
| "What are the possible values?" | `dictionary.md` |
| "Show example data" | `sample.json` |
| "How do I download?" | `download.md` |
| "What's the API?" | `download.md` |
| "What versions are available?" | `download.md` |
| "Can I use commercially?" | `license.md` |
| "What are the restrictions?" | `license.md` |
| "How to transform to unified?" | `mapping.xslt` |
| "What other DBs link to this?" | `xrefs.md` |

## Navigation Strategy

1. **Start at category level** - `README.md` for domain overview
2. **Drill to subcategory** - For unified schema across sources
3. **Drill to source** - For specific data details

## File Selection by Task

### For Code Generation
1. Read `sample.json` - Get real example records
2. Read `schema.json` - Get field types and constraints
3. Read `dictionary.md` - Understand field semantics

### For Data Integration
1. Read `mapping.xslt` - Understand transformation to unified schema
2. Read `xrefs.md` - Find cross-references to other databases
3. Read parent's `schema.json` - Get unified schema target

### For Data Access
1. Read `download.md` - Get API endpoints, rate limits, versions
2. Read `license.md` - Verify usage rights

### For Understanding Data
1. Read `README.md` - Get overview, limitations, use cases
2. Read `dictionary.md` - Get field definitions and semantics

## Hierarchy Levels

```
{category}/                    # Domain (e.g., genetics.genomics)
├── README.md                  # Category overview
├── schema.json                # Root unified schema
├── dictionary.md              # Domain field definitions
│
└── {subcategory}/             # Sub-domain (e.g., variant.repositories)
    ├── README.md              # Subcategory overview + source listing
    ├── schema.json            # Unified schema (standalone)
    ├── dictionary.md          # Field definitions
    ├── mapping.xslt           # Transform → ../schema.json
    │
    └── {source}/              # Individual data source
        ├── README.md          # Source overview + limitations
        ├── schema.json        # Source schema
        ├── dictionary.md      # Source field semantics
        ├── sample.json        # 2-5 real example records
        ├── download.md        # Access methods + versions
        ├── license.md         # License terms
        ├── mapping.xslt       # Transform → ../schema.json
        └── xrefs.md           # Cross-references (optional)
```

## Key Principles

1. **Schemas are standalone** - They document structure, not relationships
2. **Mappings define relationships** - XSLT 3.0 transforms to parent schema
3. **Each file has ONE purpose** - Load only what you need
4. **dictionary.md at every level** - Source-specific semantics in source folder
