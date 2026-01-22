# Data Sources Structure Improvement Plan

**Date:** January 2026
**Status:** Ready for Implementation
**Swarm Analysis:** 4 parallel agents completed

---

## Executive Summary

Analysis of 102 markdown files in `/home/claude/src/gene/docs/data/source/` revealed strong content consistency but significant navigation gaps. The structure is 67% incomplete for navigation links.

---

## File & Folder Organization

### Directory Structure

```
docs/data/source/                    # Root (102 files total)
├── index.md                         # Master navigation index
├── schema-gaps.md                   # Gap analysis (orphaned)
├── world1-schema-research.md        # Research doc (orphaned)
│
├── genetics/          (2 files)     # World 1: Modern Genetics
│   ├── primary.md                   # dbSNP, ClinVar, PharmGKB, etc.
│   └── population.md                # Population genetics
│
├── traditional/       (5 files)     # World 2: Traditional Medicine
│   ├── tcm.md                       # Chinese medicine (20 DBs)
│   ├── ayurveda.md                  # Indian medicine (16 DBs)
│   ├── kampo.md                     # Japanese medicine (15 DBs)
│   ├── western-herbal.md            # Western herbs (27 DBs)
│   └── global.md                    # African, Latin American (15 DBs)
│
├── compounds/         (3 files)     # Interventions
│   ├── pharmaceuticals.md           # PharmGKB, CPIC, DrugBank
│   ├── natural-products.md          # COCONUT, LOTUS, NPASS
│   └── drug-metabolism.md           # CYP450, transporters
│
├── diseases/          (10 files)    # Health Domains
│   ├── mental-cognitive.md
│   ├── cardio-metabolic.md
│   ├── cancer-oncology.md
│   ├── autoimmune.md
│   ├── rare.md
│   ├── womens-pediatric.md
│   ├── microbiome.md
│   ├── allergy-pain.md
│   ├── sleep-longevity-nutri.md
│   └── oral-skin-sensory.md
│
├── pathways/          (3 files)     # Biological Mechanisms
│   ├── primary.md                   # Reactome, WikiPathways, KEGG
│   ├── disease.md                   # DisGeNET, OMIM, HPO
│   └── processes.md                 # Process catalogs
│
├── literature/        (6 files)     # Research Papers
│   ├── sources.md                   # PubMed, PMC, OpenAlex
│   ├── public-sources.md
│   ├── coverage-analysis.md
│   ├── abstracts-vs-fulltext.md
│   ├── pipeline-design.md
│   └── data-structures.md
│
├── clinical/          (2 files)     # Biomarkers & Labs
│   ├── biomarkers-labs.md
│   └── environmental-mito.md
│
├── integration/       (11 files)    # Cross-references & APIs
│   ├── xrefs.md                     # ID mapping
│   ├── alt-sources.md
│   ├── wikipedia-wikidata.md
│   ├── wikidata-*.md                # (4 wikidata files)
│   └── ...
│
├── schemas/           (47 files)    # Database Schemas
│   ├── index.md                     # Schema navigation
│   ├── schemas-index.md             # Detailed catalog
│   ├── *-schema.md                  # 41 database schemas
│   └── ...
│
├── research/          (2 files)     # Priority Recommendations
│   ├── interventions-priority.md    # 6-agent swarm synthesis
│   └── literature-priority.md       # 5-agent swarm synthesis
│
├── downloads/         (5 files)     # Bulk Download Guides
│   ├── processing-pipeline.md
│   ├── traditional-medicine.md
│   ├── pharmaceuticals.md
│   ├── pathways-targets.md
│   └── wikidata-bulk.md
│
├── governance/        (2 files)     # Legal & Curation
│   ├── curation-framework.md
│   └── data-access-legal.md
│
└── community/         (1 file)      # Patient Networks
    └── patient-networks.md
```

### Organization Principles

| Principle | Implementation |
|-----------|---------------|
| **Conceptual Framework** | THREE WORLDS (Genetics, Traditional Medicine, Nutrition) |
| **Naming Convention** | Lowercase with hyphens (`mental-cognitive.md`) |
| **File Format** | Markdown with standard headers |
| **Tier System** | Tier 1 (MVP), Tier 2 (Post-MVP), Tier 3 (Future) |
| **Schemas Separate** | All database schemas in `/schemas/` subdirectory |

### File Counts by Directory

| Directory | Files | Purpose |
|-----------|-------|---------|
| schemas | 47 | Database field definitions |
| integration | 11 | Cross-reference APIs |
| diseases | 10 | Health domain databases |
| literature | 6 | Paper sources |
| traditional | 5 | Traditional medicine systems |
| downloads | 5 | Bulk data guides |
| compounds | 3 | Pharmaceuticals & natural products |
| pathways | 3 | Biological pathways |
| clinical | 2 | Biomarkers & labs |
| genetics | 2 | SNP & variant databases |
| governance | 2 | Legal framework |
| research | 2 | Priority syntheses |
| community | 1 | Patient networks |
| **root** | 3 | Index + 2 orphaned files |

---

## Current State Assessment

### Scores

| Category | Score | Status |
|----------|-------|--------|
| Content Formatting | 9/10 | Excellent |
| Folder Structure | 7/10 | Good foundation, missing indexes |
| Navigation/Links | 4/10 | Significant gaps |
| Schema Completeness | 8/10 | 14 missing from priorities |

### What's Working Well

- 100% lowercase naming consistency
- 100% hyphen separator consistency
- Standardized document headers (TL;DR, Key Decisions, etc.)
- Uniform database entry table format
- 95% schema content quality (fields, types, samples, APIs)
- Consistent Tier 1/2/3 priority system
- Thorough license documentation
- 13 logical domain categories

---

## Issues Identified

### Critical (Navigation)

| Issue | Scope | Files |
|-------|-------|-------|
| Files missing parent links | 67% of docs | 69 of 102 |
| Orphaned files (no links in/out) | Undiscoverable | 19 files |
| Subdirectories missing index.md | 92% of folders | 12 of 13 |
| index.md missing World 3 section | Main navigation | 1 file |

### Orphaned Files List

```
Root:
- world1-schema-research.md

Directories without links from index.md:
- community/patient-networks.md
- governance/curation-framework.md
- governance/data-access-legal.md
- genetics/population.md
- diseases/oral-skin-sensory.md
- pathways/processes.md
- integration/wikidata-master-reference.md
- integration/wikidata-supplements.md
- integration/wikidata-traditional.md
- integration/pathway-target-mapping.md

Files in index but no parent link:
- literature/public-sources.md
- literature/coverage-analysis.md
- literature/abstracts-vs-fulltext.md
- literature/pipeline-design.md
- literature/data-structures.md
- research/interventions-priority.md
- research/literature-priority.md
- downloads/*.md (5 files)
```

### Medium Priority

| Issue | Files |
|-------|-------|
| Missing schemas from priority lists | 14 databases |
| Schema files without `-schema` suffix | 4 files |
| Duplicate schema index files | 2 files |
| Orphaned root files to relocate | 2 files |
| Owner field typo | 1 file |

### Missing Schemas (from priority files)

**High Priority:**
1. PubMed/PMC - Literature source (critical for RAG)
2. DrugBank - Pharmaceutical database
3. HERB 2.0 - TCM gene expression
4. NPAtlas - Microbial natural products
5. NPASS - Quantitative bioactivity

**Medium Priority:**
6. TCMBank, SymMap, TM-MC 2.0, OpenAlex, Europe PMC, ODS API, OSADHI

---

## Implementation Plan

### Phase 1: Navigation Fixes (High Impact)

**1.1 Create subdirectory index files (12 new files)**

Create `index.md` in each directory:
```
docs/data/source/clinical/index.md
docs/data/source/community/index.md
docs/data/source/compounds/index.md
docs/data/source/diseases/index.md
docs/data/source/downloads/index.md
docs/data/source/genetics/index.md
docs/data/source/governance/index.md
docs/data/source/integration/index.md
docs/data/source/literature/index.md
docs/data/source/pathways/index.md
docs/data/source/research/index.md
docs/data/source/traditional/index.md
```

Template for each:
```markdown
# [Domain] Data Sources

**Parent:** [../index.md](../index.md)

## Documents

| Document | Content | Status |
|----------|---------|--------|
| [file.md](./file.md) | Description | Final |
```

**1.2 Add parent links to 69 files**

Add to header of each file:
```markdown
**Parent:** [../index.md](../index.md)
```

**1.3 Update main index.md**

- Add "World 3: Nutritional Science" section
- Add links to 19 orphaned files
- Add links to new subdirectory indexes

**1.4 Add missing links to docs/data/source/index.md**

Add these sections/links:
```markdown
## Community
| [community/patient-networks.md](./community/patient-networks.md) | Patient networks | Final |

## Governance
| [governance/curation-framework.md](./governance/curation-framework.md) | Curation framework | Final |
| [governance/data-access-legal.md](./governance/data-access-legal.md) | Legal framework | Final |
```

### Phase 2: Structure Cleanup (Medium Impact)

**2.1 Relocate orphaned root files**
```bash
mv docs/data/source/schema-gaps.md docs/data/source/schemas/gaps.md
mv docs/data/source/world1-schema-research.md docs/data/source/research/world1-schema-research.md
```

**2.2 Consolidate schema indexes**
- Merge `schemas/schemas-index.md` content into `schemas/index.md`
- Delete `schemas/schemas-index.md`

**2.3 Rename schema files for consistency**
```bash
mv schemas/foodb.md schemas/foodb-schema.md
mv schemas/dsld-nih.md schemas/dsld-nih-schema.md
mv schemas/usda-fooddata-central.md schemas/usda-fooddata-central-schema.md
mv schemas/pathway-formats.md schemas/pathway-formats-schema.md
```

**2.4 Fix owner field typo**
- File: `clinical/biomarkers-labs.md`
- Change: `**Owner:** Engineering` → `**Owner:** Data Engineering`

### Phase 3: Content Gaps (Future)

**3.1 Create missing critical schemas**
- `schemas/pubmed-schema.md`
- `schemas/drugbank-schema.md` (when license obtained)
- `schemas/herb2-schema.md`

**3.2 Add World context to headers (optional)**
Add to each document header:
```markdown
**World:** 1 - Modern Genetics
```

---

## Files to Modify

### Phase 1 Files
- `docs/data/source/index.md` - Add World 3, orphan links
- 69 files - Add parent links
- 12 new files - Create subdirectory indexes

### Phase 2 Files
- `docs/data/source/schema-gaps.md` - Move
- `docs/data/source/world1-schema-research.md` - Move
- `docs/data/source/schemas/schemas-index.md` - Merge & delete
- `docs/data/source/schemas/foodb.md` - Rename
- `docs/data/source/schemas/dsld-nih.md` - Rename
- `docs/data/source/schemas/usda-fooddata-central.md` - Rename
- `docs/data/source/schemas/pathway-formats.md` - Rename
- `docs/data/source/clinical/biomarkers-labs.md` - Fix owner

---

## Verification

1. **Link validation**: Run markdown link checker on `docs/data/source/`
2. **Structure check**: Verify each subdirectory has `index.md`
3. **Parent link check**: Grep for `**Parent:**` - should find 100+ matches
4. **Orphan check**: Verify all files linked from main index
5. **Navigation test**: Manually traverse from index → subdirectory → document → back

---

## Effort Summary

| Phase | Scope | Effort |
|-------|-------|--------|
| Phase 1 | 80+ files | 3-4 hours |
| Phase 2 | 8 files | 1 hour |
| Phase 3 | 15+ files | 4-6 hours (future) |
| **Total** | **100+ files** | **8-11 hours** |

---

## Target State Architecture (AI-Optimized)

### Research Summary (3-Agent Swarm)

**Enterprise Best Practices (Google, Netflix, Airbnb, LinkedIn):**
- Data catalogs with machine-readable metadata (W3C DCAT v3, Schema.org)
- Entity-first organization over function-first
- Centralized ontology with distributed ownership (Data Mesh pattern)
- Max 2-level nesting for discoverability (Dataplex)
- Automated lineage tracking (Apache Atlas pattern)
- URN identifiers for documents: `urn:gene:datasource:[domain]:[name]`
- Tag taxonomy with hierarchical tags: `domain:genetics`, `priority:tier-1`

**AI/LLM Documentation Patterns:**
- YAML frontmatter for structured metadata
- Semantic chunking: 256-512 tokens optimal for RAG
- Hierarchical embedding indexing (category → document → section)
- Explicit cross-references for graph traversal
- Consistent entity naming (kebab-case)

**Key Insight:** Current structure already follows Data Mesh principles (domain-driven). Enhance with DCAT-compliant frontmatter, machine-readable catalog, and explicit lineage relationships.

---

### Target State Directory Structure

```
docs/data/
├── _index.md                        # Master navigation (human + AI)
├── _catalog.yaml                    # Machine-readable database catalog
├── _ontology.yaml                   # Entity relationships & hierarchy
│
├── databases/                       # PRIMARY: One doc per database
│   ├── _index.md                    # Database catalog with tier matrix
│   │
│   ├── genetics/                    # WORLD 1: Modern Genetics
│   │   ├── _index.md
│   │   ├── clinvar.md               # One file per database
│   │   ├── dbsnp.md
│   │   ├── pharmgkb.md
│   │   ├── cpic.md
│   │   └── ...
│   │
│   ├── traditional/                 # WORLD 2: Traditional Medicine
│   │   ├── _index.md
│   │   ├── tcm/                     # Sub-grouped by system
│   │   │   ├── batman-tcm.md
│   │   │   ├── herb.md
│   │   │   └── tcmbank.md
│   │   ├── kampo/
│   │   │   ├── kampodb.md
│   │   │   └── stork.md
│   │   ├── ayurveda/
│   │   │   ├── imppat.md
│   │   │   └── osadhi.md
│   │   └── western/
│   │       ├── dsld.md
│   │       └── dr-dukes.md
│   │
│   ├── nutrition/                   # WORLD 3: Nutritional Science
│   │   ├── _index.md
│   │   ├── foodb.md
│   │   ├── usda-fdc.md
│   │   └── ...
│   │
│   ├── pathways/                    # Biological Mechanisms
│   │   ├── _index.md
│   │   ├── reactome.md
│   │   ├── kegg.md
│   │   └── wikipathways.md
│   │
│   ├── literature/                  # Research Papers
│   │   ├── _index.md
│   │   ├── pubmed.md
│   │   ├── pmc.md
│   │   └── openalex.md
│   │
│   └── compounds/                   # Natural Products & Drugs
│       ├── _index.md
│       ├── coconut.md
│       ├── lotus.md
│       ├── drugbank.md
│       └── chembl.md
│
├── domains/                         # Health domain VIEWS (not primary)
│   ├── _index.md                    # Links to databases by domain
│   ├── mental-cognitive.md          # References databases/*/...
│   ├── cardio-metabolic.md
│   ├── cancer-oncology.md
│   └── ...                          # Other health domains
│
├── operations/                      # HOW we work with data
│   ├── _index.md
│   ├── downloads/                   # Bulk download guides
│   │   ├── _index.md
│   │   └── *.md
│   ├── integration/                 # Cross-references & APIs
│   │   ├── _index.md
│   │   ├── id-mapping.md
│   │   └── *.md
│   ├── governance/                  # Legal & curation
│   │   ├── _index.md
│   │   └── *.md
│   └── schemas/                     # Database field definitions
│       ├── _index.md
│       └── *.md
│
└── research/                        # Analysis & decisions
    ├── _index.md
    ├── priorities/                  # Swarm synthesis docs
    │   ├── interventions.md
    │   └── literature.md
    └── analysis/                    # Gap analyses
        └── *.md
```

---

### AI-Optimized Document Template

Each database document follows this structure:

```yaml
---
# YAML FRONTMATTER (Machine-readable)
id: batman-tcm-2
name: BATMAN-TCM 2.0
world: 2
category: traditional-medicine
subcategory: tcm
tier: 1
license: CC-BY-NC
access_methods: [rest-api, bulk-download]
record_counts:
  formulas: 54832
  herbs: 8404
  compounds: 39171
  targets: 2300000
identifiers:
  primary: compound_id
  cross_refs: [pubchem_cid, chembl_id, inchikey]
relationships:
  - type: contains
    target: databases/compounds/pubchem.md
  - type: overlaps
    target: databases/traditional/tcm/herb.md
tags: [tcm, target-prediction, compound-herb, formula]
last_verified: 2026-01-15
---

# BATMAN-TCM 2.0

**Parent:** [../_index.md](../_index.md) | **World:** 2 - Traditional Medicine | **Tier:** 1

## TL;DR
One paragraph summary optimized for RAG retrieval (~100 words).

## Key Facts
| Attribute | Value |
|-----------|-------|
| Records | 54,832 formulas, 39,171 compounds |
| Access | REST API, Bulk Download |
| License | CC-BY-NC |
| Last Update | 2024 |

## Schema
[Link to operations/schemas/batman-tcm-schema.md]

## Integration Notes
Cross-references, ID mapping strategies.

## Download Guide
[Link to operations/downloads/tcm.md#batman-tcm]
```

---

### Machine-Readable Catalog (_catalog.yaml)

```yaml
# Auto-generated from document frontmatter
version: "2.0"
generated: "2026-01-22"
databases:
  - id: batman-tcm-2
    path: databases/traditional/tcm/batman-tcm.md
    world: 2
    tier: 1
    categories: [traditional-medicine, tcm]

  - id: pharmgkb
    path: databases/genetics/pharmgkb.md
    world: 1
    tier: 1
    categories: [genetics, pharmacogenomics]

# ... all 256+ databases
```

---

### Benefits of Target State

| Benefit | Current | Target |
|---------|---------|--------|
| **Database Discovery** | Scattered across 13 folders | Single `databases/` tree |
| **AI Retrieval** | No structured metadata | YAML frontmatter + _catalog.yaml |
| **RAG Chunking** | Variable structure | Consistent sections = predictable chunks |
| **Cross-references** | Manual links | frontmatter.relationships |
| **Tier Filtering** | Requires reading content | `tier: 1` in frontmatter |
| **World Navigation** | Implicit in folder | Explicit `world: 2` field |
| **Schema Lookup** | Separate /schemas folder | Link in each database doc |
| **Stale Detection** | None | `last_verified` field |

---

### Migration Strategy

**Approach:** Incremental restructure (not big-bang)

**Phase A: Add Frontmatter (Non-Breaking)**
1. Add YAML frontmatter to all 47 schema files
2. Add YAML frontmatter to top 20 Tier 1 database docs
3. Generate `_catalog.yaml` from frontmatter
4. No folder changes yet

**Phase B: Create Target Structure (Parallel)**
1. Create `databases/` hierarchy (empty)
2. Move schema files to `operations/schemas/`
3. Migrate database docs one World at a time
4. Update links incrementally

**Phase C: Domain Views (Final)**
1. Convert `diseases/` files to reference-only views
2. Create cross-reference indexes
3. Archive old structure

---

### Recommended Immediate Actions (Phase 1 + A Combined)

Instead of just fixing navigation links, also add frontmatter:

1. **Create 12 subdirectory _index.md files** with YAML frontmatter
2. **Add frontmatter to 20 Tier 1 database docs** (PharmGKB, CPIC, ClinVar, BATMAN-TCM, etc.)
3. **Generate _catalog.yaml** script
4. **Fix 69 missing parent links**
5. **Create _ontology.yaml** (THREE WORLDS structure)

This positions the codebase for full migration while delivering immediate navigation improvements.
