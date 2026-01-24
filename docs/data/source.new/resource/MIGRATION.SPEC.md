# Data Source Migration Specification

**Version:** 1.0.0
**Created:** 2026-01-23
**Status:** Draft
**Purpose:** Define folder and file organization templates for data source migration

---

## Overview

This specification defines the standard folder structure and file organization pattern for migrating data source documentation from `/docs/data/source.new/` into the hierarchical resource taxonomy at `/docs/data/source.new/resource/`.

The taxonomy consists of:
- **9 main categories** (e.g., `01.genetics.genomics/`)
- **44 subcategories** (e.g., `1.1.variant.repositories/`)
- **126+ data source folders** (e.g., `dbsnp/`)

---

## 1. Standard Folder Structure Per Data Source

Each data source folder follows this structure:

```
resource/{category}/{subcategory}/{source.name}/
    _index.md              # Primary documentation file (REQUIRED)
    schema.md              # Technical schema documentation (REQUIRED)
    download.md            # Data acquisition procedures (CONDITIONAL)
    sample.data.json       # Sample records (RECOMMENDED)
    xrefs.md               # Cross-reference mappings (CONDITIONAL)
    changelog.md           # Version history (OPTIONAL)
```

### File Descriptions

| File | Status | Purpose |
|------|--------|---------|
| `_index.md` | REQUIRED | Primary entry point with overview, key features, use cases |
| `schema.md` | REQUIRED | Technical schema: tables, fields, relationships, API endpoints |
| `download.md` | CONDITIONAL | Required if download is non-trivial (FTP, auth, processing) |
| `sample.data.json` | RECOMMENDED | Example records for validation and testing |
| `xrefs.md` | CONDITIONAL | Required for polyhierarchical sources with symlinks |
| `changelog.md` | OPTIONAL | Version history for frequently updated sources |

### Minimal Structure (Simple Sources)

For simple data sources, the minimal structure is:

```
resource/{category}/{subcategory}/{source.name}/
    _index.md              # Combined overview + schema + download
```

---

## 2. File Naming Conventions

### 2.1 Folder Names

**Pattern:** `{name}` or `{name.subname}` using lowercase with dots as separators

**Rules:**
1. Use dots (`.`) to separate words, never hyphens or underscores
2. All lowercase
3. Omit common suffixes (db, database)
4. Use standard abbreviations where universal

**Examples:**
```
GOOD:
- dbsnp/
- clinvar/
- gwas.catalog/
- open.food.facts/
- batman.tcm/
- uk.biobank/
- gene.ontology/

BAD:
- dbSNP/           (mixed case)
- clin-var/        (hyphen)
- gwas_catalog/    (underscore)
- gwascatalog/     (no separator)
- dbsnp-database/  (redundant suffix + hyphen)
```

### 2.2 File Names

**Pattern:** `{descriptor}.{extension}` using lowercase with dots

**Standard Files:**
| File | Pattern | Description |
|------|---------|-------------|
| Index | `_index.md` | Entry point (underscore prefix for sorting) |
| Schema | `schema.md` | Technical documentation |
| Download | `download.md` | Acquisition procedures |
| Sample Data | `sample.data.json` | Example records |
| Cross-refs | `xrefs.md` | External mappings |
| Changelog | `changelog.md` | Version history |

**Custom Files (if needed):**
```
api.endpoints.md       # API documentation
etl.pipeline.md        # Processing instructions
query.examples.md      # Query recipes
integration.notes.md   # Integration guidance
```

### 2.3 Placeholder Folders

Sources not yet documented use curly brace notation:

```
{hit.2.0}/            # Placeholder - awaiting documentation
{pathwaycommons}/     # Placeholder - planned addition
{ebasis}/             # Placeholder - future source
```

---

## 3. Required File Content Sections

### 3.1 _index.md Structure

```markdown
---
id: {source.id}
title: "{Source Full Name}"
type: data-source
category: {category}
subcategory: {subcategory}
parent: ../_index.md
tier: 1 | 2 | 3
last_updated: YYYY-MM-DD
status: active | planned | deprecated
tags: [category-tags]
---

# {Source Name}

**Category:** [{Category Name}](../_index.md) > [{Subcategory Name}](./_index.md)

## Overview

{2-3 paragraph description of the data source, its purpose, and significance}

## Key Statistics

| Metric | Value |
|--------|-------|
| Records | {count} |
| Last Update | {date} |
| Coverage | {description} |
| Storage | {size estimate} |

## Primary Use Cases

1. {Use case 1}
2. {Use case 2}
3. {Use case 3}

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| {ID Type} | {pattern} | {example} |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| {API/FTP/Web} | {URL} | {rate limits, auth} |

## License

| Aspect | Value |
|--------|-------|
| License | {license name} |
| Commercial Use | {Yes/No/Restricted} |
| Attribution | {Required/Optional} |

## See Also

- [Schema Documentation](./schema.md)
- [Download Instructions](./download.md)
- [{Related Source}](../../../{path}/_index.md)
```

### 3.2 schema.md Structure

```markdown
---
id: schema-{source.id}
title: "{Source Name} Schema Documentation"
type: schema
parent: _index.md
last_updated: YYYY-MM-DD
status: migrated | draft | final
tags: [schema, database, {domain-tags}]
---

# {Source Name} Schema Documentation

**Document ID:** SCHEMA-{SOURCE-ID}
**Version:** {version}
**Source Version:** {upstream version/date}

---

## TL;DR

{2-3 sentence summary of what this schema contains and key facts}

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| {Metric 1} | {Value} | {Release/Date} |

---

## Entity Relationship Overview

```
{ASCII diagram showing core entity relationships}
```

---

## Core Tables/Entities

### {Entity 1}

**Description:** {purpose}

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| {field} | {type} | {Yes/No} | {description} |

---

## API Endpoints (if applicable)

| Endpoint | Method | Description |
|----------|--------|-------------|
| {path} | {GET/POST} | {description} |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | {format, compression} |
| Alternative | {other formats} |
| Encoding | UTF-8 |

---

## Sample Query/Record

```json
{
  "example": "record"
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| {term} | {definition} |

---

## References

1. {Official documentation URL}
2. {Key publication}
```

### 3.3 download.md Structure

```markdown
---
id: download-{source.id}
title: "{Source Name} Download Instructions"
type: download
parent: _index.md
last_updated: YYYY-MM-DD
---

# {Source Name} Download Instructions

## Quick Start

```bash
# {Minimal download command}
```

## Prerequisites

- {Requirement 1}
- {Requirement 2}

## Download Methods

### Method 1: {Primary Method}

{Step-by-step instructions}

### Method 2: {Alternative Method}

{Alternative instructions}

## File Inventory

| File | Size | Description |
|------|------|-------------|
| {filename} | {size} | {description} |

## Post-Download Processing

```bash
# {Processing commands}
```

## Verification

```bash
# {Checksum or validation commands}
```

## Update Schedule

| Release | Frequency | Notes |
|---------|-----------|-------|
| {Release type} | {frequency} | {notes} |
```

---

## 4. Cross-Reference Handling

### 4.1 Symlink Pattern for Polyhierarchical Sources

Sources that belong to multiple categories use symlinks pointing to the canonical location.

**Canonical Location:** The primary category where the source naturally belongs
**Symlinks:** Additional categories where the source is relevant

**Example: IMPPAT (Ayurvedic Database)**

```
# Canonical location
05.traditional.medicine/5.2.south.east.asian.systems/imppat/
    _index.md
    schema.md

# Symlink in compounds category
02.compounds.molecules/2.1.natural.products/imppat -> ../../../05.traditional.medicine/5.2.south.east.asian.systems/imppat

# Symlink in traditional medicine compounds
02.compounds.molecules/2.3.traditional.medicine.compounds/imppat -> ../../../05.traditional.medicine/5.2.south.east.asian.systems/imppat
```

### 4.2 xrefs.md for Cross-References

Create `xrefs.md` in the canonical location documenting symlink relationships:

```markdown
---
id: xrefs-{source.id}
title: "{Source Name} Cross-References"
type: xrefs
parent: _index.md
---

# {Source Name} Cross-References

## Taxonomy Locations

| Location | Relationship | Path |
|----------|--------------|------|
| Primary (canonical) | Native | `{path}` |
| Secondary | Symlink | `{symlink-path}` |

## External Cross-References

| Database | ID Mapping | Coverage |
|----------|------------|----------|
| UniProt | {mapping} | {%} |
| PubChem | {mapping} | {%} |

## ID Mapping Table

| Source ID | External ID | Notes |
|-----------|-------------|-------|
| {pattern} | {pattern} | {notes} |
```

### 4.3 Current Symlink Inventory

Existing polyhierarchical sources with symlinks:

| Source | Canonical Location | Symlink Locations |
|--------|-------------------|-------------------|
| imppat | 05.traditional.medicine/5.2.south.east.asian.systems/ | 02.compounds.molecules/2.1.natural.products/ |
| batman.tcm | 05.traditional.medicine/5.1.traditional.chinese.medicine/ | 02.compounds.molecules/2.3.traditional.medicine.compounds/ |
| kampodb | 05.traditional.medicine/5.2.south.east.asian.systems/ | 02.compounds.molecules/2.3.traditional.medicine.compounds/ |
| herb | 05.traditional.medicine/5.1.traditional.chinese.medicine/ | 02.compounds.molecules/2.3.traditional.medicine.compounds/ |
| foodb | 06.nutrition.food/6.1.food.composition/ | 02.compounds.molecules/2.4.food.compounds.nutrients/ |
| usda.fooddata | (symlink target) | 06.nutrition.food/6.1.food.composition/, 02.compounds.molecules/2.4.food.compounds.nutrients/ |
| phenol.explorer | (symlink target) | 06.nutrition.food/6.3.bioactive.food.compounds/ |
| phytohub | (symlink target) | 06.nutrition.food/6.3.bioactive.food.compounds/ |
| string | 04.pathways.networks/4.3.protein.protein.interactions/ | 07.proteins.molecular.biology/7.3.molecular.interactions/ |
| intact | 04.pathways.networks/4.3.protein.protein.interactions/ | 07.proteins.molecular.biology/7.3.molecular.interactions/ |
| reactome | 04.pathways.networks/4.1.metabolic.pathways/ | 07.proteins.molecular.biology/7.3.molecular.interactions/ |
| uniprot | 07.proteins.molecular.biology/7.1.protein.sequences.annotations/ | 08.literature.knowledge/8.3.identifier.mapping/ |
| disgenet | 03.diseases.phenotypes/3.3.disease.gene.associations/ | 01.genetics.genomics/1.5.expression.regulation/ |
| pharmgkb | 01.genetics.genomics/1.4.pharmacogenomics/ | 02.compounds.molecules/2.5.drug.metabolism.pharmacokinetics/ |
| cpic | 01.genetics.genomics/1.4.pharmacogenomics/ | 02.compounds.molecules/2.5.drug.metabolism.pharmacokinetics/ |
| chembl | 02.compounds.molecules/2.2.pharmaceuticals/ | 02.compounds.molecules/2.7.compound.target.interactions/ |
| dailymed | 02.compounds.molecules/2.2.pharmaceuticals/ | 08.literature.knowledge/8.4.regulatory.legal/ |
| omim | 03.diseases.phenotypes/3.2.phenotype.databases/ | 03.diseases.phenotypes/3.5.rare.orphan.diseases/ |
| wikidata | 08.literature.knowledge/8.2.knowledge.bases/ | 05.traditional.medicine/5.4.multi.system.integration/ |
| dr.dukes | 02.compounds.molecules/2.1.natural.products/ | 05.traditional.medicine/5.3.western.global.herbal/ |

---

## 5. Content Migration Mapping

### 5.1 Source Content Locations

| Content Type | Current Location | Destination |
|--------------|------------------|-------------|
| Schema docs | `operations/schemas/{source}-schema.md` | `resource/{path}/schema.md` |
| Download guides | `operations/downloads/` | `resource/{path}/download.md` |
| Integration notes | `operations/integration/` | `resource/{path}/xrefs.md` |
| Database overviews | `databases/{category}/` | `resource/{path}/_index.md` |
| Domain views | `domains/` | Reference via cross-links (no move) |

### 5.2 Migration Workflow

1. **Identify source** in existing schema documentation
2. **Locate canonical folder** in resource hierarchy
3. **Create _index.md** with frontmatter from existing content
4. **Move/transform schema.md** from operations/schemas/
5. **Create download.md** if non-trivial acquisition
6. **Create xrefs.md** if symlinks exist
7. **Update cross-references** in related files
8. **Validate symlinks** point to canonical location

---

## 6. Frontmatter Standards

### 6.1 Required Frontmatter Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | Yes | Unique identifier (lowercase, dots) |
| `title` | string | Yes | Human-readable title |
| `type` | enum | Yes | `data-source`, `schema`, `download`, `xrefs` |
| `parent` | path | Yes | Relative path to parent _index.md |
| `last_updated` | date | Yes | YYYY-MM-DD format |
| `status` | enum | Yes | `active`, `planned`, `deprecated`, `migrated` |

### 6.2 Optional Frontmatter Fields

| Field | Type | Description |
|-------|------|-------------|
| `category` | string | Primary category (`genetics`, `traditional`, `nutrition`, `shared`) |
| `subcategory` | string | Subcategory for traditional medicine (`tcm`, `ayurveda`, etc.) |
| `tier` | integer | Implementation priority (1, 2, 3) |
| `tags` | array | Searchable keywords |
| `version` | string | Schema/API version documented |
| `license` | string | Data license type |

### 6.3 ID Naming Pattern

```
# Data sources
id: {source.name}              # e.g., "dbsnp", "clinvar", "batman.tcm"

# Schemas
id: schema-{source.name}       # e.g., "schema-dbsnp"

# Downloads
id: download-{source.name}     # e.g., "download-gnomad"

# Cross-references
id: xrefs-{source.name}        # e.g., "xrefs-imppat"
```

---

## 7. Validation Checklist

### Per Data Source Folder

- [ ] `_index.md` exists with required frontmatter
- [ ] `schema.md` exists with technical documentation
- [ ] Frontmatter `id` matches folder name pattern
- [ ] Frontmatter `parent` points to valid file
- [ ] All internal links resolve
- [ ] Symlinks (if any) point to valid canonical location
- [ ] `xrefs.md` exists if source has symlinks

### Per Category

- [ ] `_index.md` exists at subcategory level
- [ ] All sources listed in subcategory index
- [ ] No orphan folders (empty or placeholder only)

### Global

- [ ] No duplicate source IDs across taxonomy
- [ ] All symlinks have corresponding canonical folders
- [ ] Cross-category references use correct relative paths

---

## 8. Example: Complete dbSNP Migration

### Folder Structure

```
resource/01.genetics.genomics/1.1.variant.repositories/dbsnp/
    _index.md
    schema.md
    download.md
    sample.data.json
```

### _index.md

```markdown
---
id: dbsnp
title: "dbSNP - Database of Single Nucleotide Polymorphisms"
type: data-source
category: genetics
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [variants, snp, ncbi, population-genetics]
---

# dbSNP

**Category:** [Genetics & Genomics](../../_index.md) > [Variant Repositories](../_index.md)

## Overview

dbSNP (Database of Single Nucleotide Polymorphisms) is NCBI's public archive for genetic variation data. It contains over 1 billion variants with population frequencies from the ALFA project covering 12 populations and 400,000+ subjects.

## Key Statistics

| Metric | Value |
|--------|-------|
| Variants | 1B+ |
| ALFA Subjects | 400K+ |
| Populations | 12 |
| Reference Build | GRCh38.p14 |

## Primary Use Cases

1. Variant annotation and lookup by rsID
2. Population frequency queries via ALFA
3. Clinical variant contextualization
4. VCF annotation pipelines

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| RefSNP ID | rs[0-9]+ | rs334 |
| SPDI | NC_*:[pos]:[del]:[ins] | NC_000011.10:5227001:T:A |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| REST API | api.ncbi.nlm.nih.gov/variation/v0 | 1 req/sec recommended |
| FTP | ftp.ncbi.nlm.nih.gov/snp/ | Bulk downloads |
| VCF Files | latest_release/ | Pre-built annotations |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (NCBI) |
| Commercial Use | Yes |
| Attribution | Optional |

## See Also

- [Schema Documentation](./schema.md)
- [Download Instructions](./download.md)
- [ClinVar](../clinvar/_index.md) - Clinical variant interpretations
- [gnomAD](../../1.3.population.genetics/gnomad/_index.md) - Complementary population data
```

---

## 9. Category Index Templates

Each subcategory folder requires an `_index.md`:

```markdown
---
id: {subcategory.id}
title: "{Subcategory Name}"
type: subcategory
parent: ../_index.md
last_updated: YYYY-MM-DD
---

# {Subcategory Name}

**Parent:** [{Category Name}](../_index.md)

## Overview

{Description of this subcategory and what data sources it contains}

## Data Sources

| Source | Tier | Description |
|--------|------|-------------|
| [{Source 1}](./{source1}/_index.md) | 1 | {brief description} |
| [{Source 2}](./{source2}/_index.md) | 2 | {brief description} |

## Integration Notes

{How sources in this subcategory relate to each other and to other categories}
```

---

## Appendix A: Complete Category Hierarchy

```
resource/
├── 01.genetics.genomics/
│   ├── 1.1.variant.repositories/
│   ├── 1.2.functional.prediction/
│   ├── 1.3.population.genetics/
│   ├── 1.4.pharmacogenomics/
│   ├── 1.5.expression.regulation/
│   └── 1.6.cancer.genomics/
├── 02.compounds.molecules/
│   ├── 2.1.natural.products/
│   ├── 2.2.pharmaceuticals/
│   ├── 2.3.traditional.medicine.compounds/
│   ├── 2.4.food.compounds.nutrients/
│   ├── 2.5.drug.metabolism.pharmacokinetics/
│   ├── 2.6.chemical.ontology.classification/
│   └── 2.7.compound.target.interactions/
├── 03.diseases.phenotypes/
│   ├── 3.1.disease.ontologies/
│   ├── 3.2.phenotype.databases/
│   ├── 3.3.disease.gene.associations/
│   ├── 3.4.cancer.oncology/
│   ├── 3.5.rare.orphan.diseases/
│   ├── 3.6.autoimmune.inflammatory/
│   └── 3.7.mental.health.neurological/
├── 04.pathways.networks/
│   ├── 4.1.metabolic.pathways/
│   ├── 4.2.signaling.pathways/
│   ├── 4.3.protein.protein.interactions/
│   ├── 4.4.drug.target.interactions/
│   ├── 4.5.gene.function.ontology/
│   └── 4.6.regulatory.networks/
├── 05.traditional.medicine/
│   ├── 5.1.traditional.chinese.medicine/
│   ├── 5.2.south.east.asian.systems/
│   ├── 5.3.western.global.herbal/
│   └── 5.4.multi.system.integration/
├── 06.nutrition.food/
│   ├── 6.1.food.composition/
│   ├── 6.2.dietary.supplements/
│   ├── 6.3.bioactive.food.compounds/
│   └── 6.4.metabolomics/
├── 07.proteins.molecular.biology/
│   ├── 7.1.protein.sequences.annotations/
│   ├── 7.2.protein.structures/
│   └── 7.3.molecular.interactions/
├── 08.literature.knowledge/
│   ├── 8.1.scientific.literature/
│   ├── 8.2.knowledge.bases/
│   ├── 8.3.identifier.mapping/
│   └── 8.4.regulatory.legal/
└── 09.microbiome/
    ├── 9.1.gut.microbiome/
    ├── 9.2.body.site.microbiomes/
    └── 9.3.microbe.host.interactions/
```

---

## Appendix B: Migration Priority

| Priority | Criteria | Sources |
|----------|----------|---------|
| P1 | Tier 1 + existing schema | dbSNP, ClinVar, gnomAD, ChEMBL, Reactome |
| P2 | Tier 1 + no schema | COCONUT, LOTUS, PharmGKB |
| P3 | Tier 2 + polyhierarchical | IMPPAT, BATMAN-TCM, FooDB |
| P4 | Tier 2 + simple | All remaining Tier 2 |
| P5 | Tier 3 + placeholders | Future sources |

---

*Specification maintained by Gene Platform Data Team*
*Last updated: 2026-01-23*
