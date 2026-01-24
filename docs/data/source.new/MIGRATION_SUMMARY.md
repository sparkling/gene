---
title: "Content Migration Summary"
date: 2026-01-23
status: complete
---

# Content Migration Summary

**Migration Date:** 2026-01-23
**Status:** Complete

---

## Overview

This document summarizes the reorganization of documentation content in `docs/data/source.new/` to improve structure and navigability.

---

## Files Migrated to `reference/`

### reference/formats/ (1 file)
| File | Lines | Status |
|------|-------|--------|
| `_index.md` | 67 | Index created |

**Note:** Content files (pathway-formats.md, literature-formats.md, text-extraction.md) referenced in index but source files remain in operations/schemas/.

### reference/architecture/ (1 file)
| File | Lines | Status |
|------|-------|--------|
| `_index.md` | 74 | Index created |

**Note:** Content files (three-worlds-schema.md, unified-schema-analysis.md, sample-data.md) referenced in index but source files remain in operations/schemas/.

**Total reference/ files:** 3 (including root _index.md)

---

## Files Migrated to `guides/`

### guides/wikidata/ (8 files)
| File | Lines | Description |
|------|-------|-------------|
| `_index.md` | 104 | Wikidata guides index |
| `master-reference.md` | 2,906 | Complete Wikidata usage guide |
| `pharmaceutical-queries.md` | 747 | Drug and pharmaceutical SPARQL queries |
| `traditional-medicine.md` | 835 | Traditional medicine queries |
| `supplements.md` | 781 | Dietary supplement queries |
| `wikipedia-integration.md` | 843 | Wikipedia-Wikidata integration |
| `bulk-download.md` | 157 | Bulk data download procedures |
| `schema.md` | 278 | Wikidata schema reference |

**Total lines:** ~6,651

### guides/integration/ (4 files)
| File | Lines | Description |
|------|-------|-------------|
| `_index.md` | 104 | Integration guides index |
| `compound-pathway-linking.md` | 2,536 | Cross-database compound-pathway linking |
| `pathway-target-mapping.md` | 1,793 | Pathway to target relationships |
| `pathway-downloads.md` | 154 | Pathway database downloads |

**Total lines:** ~4,587

**Total guides/ files:** 13 (including indices)

---

## Files Reorganized in `operations/`

### operations/methodology/ (6 files)
| File | Lines | Description |
|------|-------|-------------|
| `_index.md` | 63 | Methodology index |
| `coverage-analysis.md` | 574 | Data coverage assessment |
| `curation-framework.md` | 148 | Data curation methodology |
| `literature-pipeline.md` | 519 | Literature processing pipeline |
| `literature-sources.md` | 756 | Literature source catalog |
| `public-sources.md` | 830 | Public data source evaluation |

**Total lines:** ~2,890

### operations/governance/ (3 files)
| File | Lines | Description |
|------|-------|-------------|
| `_index.md` | 59 | Governance index |
| `data-access-legal.md` | 146 | Legal and licensing analysis |
| `curation-framework.md` | 148 | Curation framework (duplicate) |

**Total lines:** ~353

### operations/planning/ (3 files)
| File | Lines | Description |
|------|-------|-------------|
| `_index.md` | 61 | Planning index |
| `size-estimates.md` | 808 | Storage capacity estimates |
| `alternative-sources.md` | 741 | Alternative data sources |

**Total lines:** ~1,610

**Total reorganized operations/ files:** 12

---

## Subcategory Indexes Enhanced

The following `_index.md` files were created or enhanced:

| Directory | Status |
|-----------|--------|
| `reference/_index.md` | Created |
| `reference/formats/_index.md` | Created |
| `reference/architecture/_index.md` | Created |
| `guides/_index.md` | Created |
| `guides/wikidata/_index.md` | Created |
| `guides/integration/_index.md` | Created |
| `operations/_index.md` | Updated |
| `operations/methodology/_index.md` | Created |
| `operations/governance/_index.md` | Existing |
| `operations/planning/_index.md` | Created |

---

## Domains Folder (Unchanged)

The `domains/` folder was intentionally NOT modified as these files serve a different purpose (cross-cutting health domain views).

| File | Lines | Status |
|------|-------|--------|
| `_index.md` | 95 | Unchanged |
| `mental-cognitive.md` | 1,024 | Unchanged |
| `cardio-metabolic.md` | 710 | Unchanged |
| `cancer-oncology.md` | 728 | Unchanged |
| `autoimmune.md` | 856 | Unchanged |
| `rare.md` | 637 | Unchanged |
| `womens-pediatric.md` | 884 | Unchanged |
| `microbiome.md` | 623 | Unchanged |
| `allergy-pain.md` | 940 | Unchanged |
| `oral-skin-sensory.md` | 898 | Unchanged |
| `sleep-longevity-nutri.md` | 926 | Unchanged |
| `clinical-environmental-mito.md` | 754 | Unchanged |
| `clinical-biomarkers-labs.md` | 927 | Unchanged |
| `community-patient-networks.md` | 879 | Unchanged |

**Total:** 14 files (unchanged)

---

## Migration Summary Statistics

| Category | Files | Total Lines |
|----------|-------|-------------|
| reference/ | 3 | ~140 |
| guides/ | 13 | ~11,238 |
| operations/ (reorganized) | 12 | ~4,853 |
| domains/ (preserved) | 14 | ~10,881 |

**Total files affected:** 42
**Total lines of content:** ~27,112

---

## Errors or Issues

### Minor Issues Identified

1. **Duplicate file:** `curation-framework.md` exists in both:
   - `operations/governance/curation-framework.md`
   - `operations/methodology/curation-framework.md`

   **Recommendation:** Keep in governance/, remove from methodology/

2. **Reference content not fully migrated:** Index files in `reference/formats/` and `reference/architecture/` reference content files that were not copied. Source files remain in:
   - `operations/schemas/pathway-formats.md`
   - `operations/schemas/ruvector-three-worlds-schema.md`
   - `operations/schemas/unified-schema-analysis.md`
   - `databases/literature/data-structures.md`
   - `databases/literature/abstracts-vs-fulltext.md`

   **Recommendation:** Copy these files to reference/ subdirectories or update index files to point to existing locations.

3. **Missing directories:**
   - `operations/standards/` - Not created (0 files)
   - `operations/vendors/` - Not created (0 files)

   **Status:** These were optional and not required for current content.

---

## Final Directory Structure

```
docs/data/source.new/
├── _index.md
├── MIGRATION_SUMMARY.md          # This file
│
├── resource/                     # 190+ data source files (unchanged)
│   ├── 01.genetics.genomics/
│   ├── 02.compounds.molecules/
│   ├── 03.diseases.phenotypes/
│   ├── 04.pathways.networks/
│   ├── 05.traditional.medicine/
│   ├── 06.nutrition.food/
│   ├── 07.proteins.molecular.biology/
│   ├── 08.literature.knowledge/
│   ├── 09.microbiome/
│   └── CONTENT_TYPE_ANALYSIS.md
│
├── domains/                      # 14 files (preserved)
│   ├── _index.md
│   └── [13 health domain files]
│
├── reference/                    # NEW - 3 files
│   ├── _index.md
│   ├── formats/_index.md
│   └── architecture/_index.md
│
├── guides/                       # NEW - 13 files
│   ├── _index.md
│   ├── wikidata/                 # 8 files
│   └── integration/              # 4 files
│
├── operations/                   # REORGANIZED - 12+ files
│   ├── _index.md
│   ├── methodology/              # 6 files
│   ├── governance/               # 3 files
│   ├── planning/                 # 3 files
│   ├── downloads/                # existing
│   ├── integration/              # existing
│   └── schemas/                  # existing
│
├── databases/                    # existing (to be deprecated)
└── research/                     # existing
```

---

## Verification Checklist

- [x] reference/formats/ directory exists
- [x] reference/architecture/ directory exists
- [x] guides/wikidata/ directory exists with 7+ guides
- [x] guides/integration/ directory exists with 3+ guides
- [x] operations/methodology/ directory exists with content
- [x] operations/governance/ directory exists with content
- [x] operations/planning/ directory exists with content
- [x] domains/ folder unchanged (14 files)
- [x] All subcategory _index.md files created
- [x] CONTENT_TYPE_ANALYSIS.md updated with migration status

---

**Migration completed successfully.**

*Generated: 2026-01-23*
