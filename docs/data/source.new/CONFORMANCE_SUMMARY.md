# Template Conformance Verification Report

**Generated:** 2026-01-23
**Report Type:** Re-verification of Template Conformance
**Base Directory:** `/home/claude/src/gene/docs/data/source.new`

## Executive Summary

A comprehensive audit of all 215 data source directories in the source.new dataset reveals **partial compliance** with the 3-file template standard:

- **Total Data Sources:** 215
- **Fully Compliant (all 3 files):** 101 sources (46.98%)
- **Partially Compliant (1-2 files):** 114 sources (53.02%)
- **Non-Compliant (0 files):** 0 sources (0%)

### Template Requirements

Each data source directory should contain:

| File | Purpose | Coverage |
|------|---------|----------|
| `_index.md` | Directory index and metadata | 215/215 (100%) ✓ |
| `schema.md` | Data schema and structure documentation | 102/215 (47.44%) ⚠️ |
| `download.md` | Download instructions and sources | 120/215 (55.81%) ⚠️ |

## Nutrition.food Resource Analysis

### Category Performance

The **06.nutrition.food** category shows **above-average compliance** with a rate of **61.54%**:

```
Category: 06.nutrition.food
├── Total sources: 13
├── Fully compliant: 8
├── Compliance rate: 61.54%
└── Status: GOOD (exceeds global average)
```

### Subdirectory Breakdown

#### 6.1 Food Composition (100% Data Source Compliance)
- **Status:** Category incomplete, but child sources compliant
- **Files Missing at Category Level:** schema.md, download.md
- **Compliant Data Sources:** 2/2
  - foodb ✓
  - open.food.facts ✓

#### 6.2 Dietary Supplements (100% Data Source Compliance)
- **Status:** Category incomplete, but child sources compliant
- **Files Missing at Category Level:** schema.md, download.md
- **Compliant Data Sources:** 3/3
  - consumerlab ✓
  - dsld ✓
  - natural.medicines ✓

#### 6.3 Bioactive Food Compounds (100% Data Source Compliance)
- **Status:** Category incomplete, but child sources compliant
- **Files Missing at Category Level:** schema.md, download.md
- **Compliant Data Sources:** 1/1
  - ebasis ✓

#### 6.4 Metabolomics (100% Data Source Compliance)
- **Status:** Category incomplete, but child sources compliant
- **Files Missing at Category Level:** schema.md, download.md
- **Compliant Data Sources:** 2/2
  - exposome.explorer ✓
  - hmdb ✓

### Key Findings

✓ **Strengths:**
- All 8 individual data sources in nutrition.food are fully compliant (100%)
- All 13 directories in the category have _index.md files
- Strong performance in schema.md coverage (8/13 = 61.54%)
- Strong performance in download.md coverage (8/13 = 61.54%)

⚠️ **Gaps:**
- Root nutrition.food directory: missing schema.md, download.md
- All 4 category directories (6.1-6.4): missing schema.md, download.md
- 5 missing files needed to achieve 100% nutrition.food compliance

## Global Compliance by Major Category

| Category | Total | Compliant | Rate | Status |
|----------|-------|-----------|------|--------|
| **resource** | 190 | 101 | 53.16% | GOOD |
| **reference** | 10 | 0 | 0% | NON-COMPLIANT |
| **ontology** | 4 | 0 | 0% | NON-COMPLIANT |
| **operations** | 4 | 0 | 0% | NON-COMPLIANT |
| **guides** | 3 | 0 | 0% | NON-COMPLIANT |
| **taxonomy** | 1 | 0 | 0% | NON-COMPLIANT |
| **clinical-views** | 1 | 0 | 0% | NON-COMPLIANT |

## Compliance Tier Analysis

```
Tier 0 (Non-indexed):           0 sources (0%)
├─ No _index.md file

Tier 1 (Indexed only):          114 sources (53.02%)
├─ Has _index.md
├─ Missing schema.md and/or download.md

Tier 2 (Partial):               0 sources (0%)
├─ Has _index.md
├─ Has one of: schema.md, download.md

Tier 3 (Fully compliant):       101 sources (46.98%)
├─ Has all three files
├─ _index.md ✓
├─ schema.md ✓
└─ download.md ✓
```

## Recommendations

### Immediate Priorities

#### Priority 1: Nutrition.food Category Directories
**Impact:** Quick path to 100% nutrition.food compliance

Required actions:
```
Create 4 schema.md files:
  /resource/06.nutrition.food/6.1.food.composition/schema.md
  /resource/06.nutrition.food/6.2.dietary.supplements/schema.md
  /resource/06.nutrition.food/6.3.bioactive.food.compounds/schema.md
  /resource/06.nutrition.food/6.4.metabolomics/schema.md

Create 4 download.md files:
  /resource/06.nutrition.food/6.1.food.composition/download.md
  /resource/06.nutrition.food/6.2.dietary.supplements/download.md
  /resource/06.nutrition.food/6.3.bioactive.food.compounds/download.md
  /resource/06.nutrition.food/6.4.metabolomics/download.md

Create 1 root schema.md file:
  /resource/06.nutrition.food/schema.md

Create 1 root download.md file:
  /resource/06.nutrition.food/download.md
```

**Expected Result:** Nutrition.food category compliance: 61.54% → 100%

#### Priority 2: Complete Global Data Source Coverage
**Impact:** Increase global compliance from 46.98% to 100%

Required actions:
- Create 57 additional schema.md files (102 exist, 215 needed)
- Create 95 additional download.md files (120 exist, 215 needed)

**Expected Result:** Global compliance: 46.98% → 100%

#### Priority 3: Clarify Directory Structure
**Impact:** Improved transparency about true compliance metrics

The non-resource categories (reference, ontology, operations, etc.) appear to be organizational rather than data sources. Consider:
- Marking these as non-data-source directories
- Or making them compliant with the template
- Updated metrics would show true data source compliance separately

## File Gap Summary

### Schema.md Gaps
- **Currently exist:** 102 files
- **Needed for full compliance:** 215 files
- **Gap:** 113 files (52.56% gap)
- **Nutrition.food gap:** 5 files (6.1, 6.2, 6.3, 6.4, root)

### Download.md Gaps
- **Currently exist:** 120 files
- **Needed for full compliance:** 215 files
- **Gap:** 95 files (44.19% gap)
- **Nutrition.food gap:** 5 files (6.1, 6.2, 6.3, 6.4, root)

## Detailed File Locations

### Nutrition.food Compliant Data Sources
All leaf node data sources are fully compliant:

```
✓ /resource/06.nutrition.food/6.1.food.composition/foodb
✓ /resource/06.nutrition.food/6.1.food.composition/open.food.facts
✓ /resource/06.nutrition.food/6.2.dietary.supplements/consumerlab
✓ /resource/06.nutrition.food/6.2.dietary.supplements/dsld
✓ /resource/06.nutrition.food/6.2.dietary.supplements/natural.medicines
✓ /resource/06.nutrition.food/6.3.bioactive.food.compounds/ebasis
✓ /resource/06.nutrition.food/6.4.metabolomics/exposome.explorer
✓ /resource/06.nutrition.food/6.4.metabolomics/hmdb
```

### Nutrition.food Non-Compliant Directories
Category and root directories missing template files:

```
✗ /resource/06.nutrition.food (missing: schema.md, download.md)
✗ /resource/06.nutrition.food/6.1.food.composition (missing: schema.md, download.md)
✗ /resource/06.nutrition.food/6.2.dietary.supplements (missing: schema.md, download.md)
✗ /resource/06.nutrition.food/6.3.bioactive.food.compounds (missing: schema.md, download.md)
✗ /resource/06.nutrition.food/6.4.metabolomics (missing: schema.md, download.md)
```

## Conclusion

The source.new dataset shows **good compliance at the data source level** (leaf nodes), with 100% of individual data sources containing _index.md files. However, **schema.md and download.md coverage is incomplete**, especially at category levels.

The nutrition.food resource demonstrates **above-average performance** with individual data sources at 100% compliance. Achieving complete nutrition.food compliance requires creating 5 additional files at category levels.

A coordinated effort to fill the identified gaps would bring the entire dataset to full template compliance, improving documentation completeness and data accessibility.

### Report Artifacts

- **JSON Report:** `/home/claude/src/gene/docs/data/source.new/CONFORMANCE_REPORT.json`
- **Summary:** `/home/claude/src/gene/docs/data/source.new/CONFORMANCE_SUMMARY.md`

---

*For detailed metrics and file lists, see the accompanying CONFORMANCE_REPORT.json file.*
