# Legacy Directory Cleanup Plan

**Document ID:** CLEANUP-PLAN-001
**Created:** 2026-01-24
**Updated:** 2026-01-24
**Status:** COMPLETED

---

## Executive Summary

**Root Cause:** Taxonomy reorganization left behind old directories with duplicate subcategory numbers.

| Subcategory # | CURRENT Directory | LEGACY Directory | Action |
|---------------|-------------------|------------------|--------|
| 2.3 | *(none - intentional gap)* | `2.3.traditional.medicine.compounds` | Keep (cross-ref) |
| 3.3 | `3.3.disease.gene.associations` | `3.3.pharmacogenomics` | **DELETE legacy** |
| 3.4 | `3.4.cancer.oncology` | `3.4.drug.disease.targets` | **DELETE legacy** |

**The CURRENT directories are complete** with unified schemas, data dictionaries, and proper symlinks.
**The LEGACY directories are duplicates** that should be removed after backup.

---

## Detailed Analysis

### 1. Directory: `2.3.traditional.medicine.compounds`

**Location:** `02.compounds.molecules/2.3.traditional.medicine.compounds/`

**Current State:**
```
2.3.traditional.medicine.compounds/
├── _index.md              (1.5KB - explains cross-reference)
├── batman.tcm -> ../../05.traditional.medicine/5.1.traditional.chinese.medicine/batman.tcm
├── herb -> ../../05.traditional.medicine/5.1.traditional.chinese.medicine/herb
└── kampodb -> ../../05.traditional.medicine/5.2.south.east.asian.systems/kampodb
```

**Analysis:** This is an INTENTIONAL cross-reference directory. It allows users browsing Category 02 (Compounds) to discover traditional medicine compounds that are primarily housed in Category 05 (Traditional Medicine).

**Action:** ✅ NO CHANGE REQUIRED - This is proper design, not legacy.

---

### 2. Directory: `3.3.pharmacogenomics`

**Location:** `03.diseases.phenotypes/3.3.pharmacogenomics/`

**Current State:**
```
3.3.pharmacogenomics/
├── _index.md              (890 bytes)
└── pharmgkb/
    └── schema.md          (11,162 bytes - MORE DETAILED)
```

**Primary Location:** `01.genetics.genomics/1.4.pharmacogenomics/pharmgkb/`
```
pharmgkb/
├── _index.md              (2,306 bytes)
├── download.md            (6,853 bytes)
├── schema.md              (7,549 bytes)
└── xrefs.md               (1,430 bytes)
```

**Analysis:**
- PharmGKB exists in TWO locations
- The LEGACY location (3.3) has a MORE COMPREHENSIVE schema.md (11KB vs 7.5KB)
- The legacy schema includes additional ClinPGx API details, more field descriptions
- Primary location has complete documentation set (download, xrefs)

**Data to Preserve:**
- Additional content from `3.3.pharmacogenomics/pharmgkb/schema.md`:
  - ClinPGx API rebranding information
  - S3 download URLs
  - Rate limit information
  - Extended field descriptions

**Action Required:**
1. MERGE: Integrate unique content from legacy schema.md into primary location
2. BACKUP: Create backup of legacy files before deletion
3. CONVERT: Replace directory with symlink to primary location
4. UPDATE: Update `_index.md` to explain cross-reference

---

### 3. Directory: `3.4.drug.disease.targets`

**Location:** `03.diseases.phenotypes/3.4.drug.disease.targets/`

**Current State:**
```
3.4.drug.disease.targets/
├── _index.md              (1,004 bytes)
└── dgidb/
    └── schema.md          (15,794 bytes - MUCH MORE DETAILED)
```

**Primary Location:** `02.compounds.molecules/2.7.compound.target.interactions/dgidb/`
```
dgidb/
├── _index.md              (3,488 bytes)
├── download.md            (6,329 bytes)
├── schema.md              (6,623 bytes)
└── xrefs.md               (1,355 bytes)
```

**Analysis:**
- DGIdb exists in TWO locations
- The LEGACY location (3.4) has a MUCH MORE COMPREHENSIVE schema.md (15.8KB vs 6.6KB)
- Legacy schema includes BOTH DGIdb AND Open Targets Platform documentation
- Open Targets Platform content should be extracted to its own data source

**Data to Preserve:**
- Additional DGIdb content from legacy schema.md:
  - Extended GraphQL schema documentation
  - More detailed field descriptions
  - Additional API examples
- Open Targets Platform content:
  - Should be extracted to separate data source in appropriate category

**Action Required:**
1. MERGE: Integrate unique DGIdb content from legacy schema.md into primary location
2. EXTRACT: Create new Open Targets Platform data source (likely in 3.3 or 4.4)
3. BACKUP: Create backup of legacy files before deletion
4. CONVERT: Replace directory with symlink to primary location
5. UPDATE: Update `_index.md` to explain cross-reference

---

## Execution Plan

### Phase 1: Backup (CRITICAL - Do First)

```bash
# Create backup directory
mkdir -p /home/claude/src/gene/docs/data/source.new/resource/_legacy_backup/

# Backup all legacy content with full structure
cp -r 03.diseases.phenotypes/3.3.pharmacogenomics/ _legacy_backup/
cp -r 03.diseases.phenotypes/3.4.drug.disease.targets/ _legacy_backup/

# Verify backups
ls -laR _legacy_backup/
```

### Phase 2: Content Merge (if needed)

Check if legacy schema.md files have unique content not in primary locations:

#### 2a. PharmGKB Schema Comparison
- Legacy: `3.3.pharmacogenomics/pharmgkb/schema.md` (11KB)
- Primary: `01.genetics.genomics/1.4.pharmacogenomics/pharmgkb/schema.md` (7.5KB)
- **Action:** If legacy has unique content (ClinPGx API, S3 URLs), merge into primary

#### 2b. DGIdb Schema Comparison
- Legacy: `3.4.drug.disease.targets/dgidb/schema.md` (16KB - includes Open Targets!)
- Primary: `02.compounds.molecules/2.7.compound.target.interactions/dgidb/schema.md` (6.6KB)
- **Action:**
  - Merge unique DGIdb content into primary
  - Note: Open Targets already exists at `3.3.disease.gene.associations/open.targets/`

### Phase 3: Delete Legacy Directories

```bash
# Delete legacy directories (after backup confirmed)
rm -rf 03.diseases.phenotypes/3.3.pharmacogenomics/
rm -rf 03.diseases.phenotypes/3.4.drug.disease.targets/
```

**Note:** We do NOT create replacement directories because:
- `3.3.disease.gene.associations` already exists with proper structure
- `3.4.cancer.oncology` already exists with proper structure
- These have unified schemas and data dictionaries
- They have symlinks to data in other categories

### Phase 4: Verification

1. Confirm only 7 subcategories remain in Category 03 (3.1-3.7)
2. Verify no broken references in documentation
3. Verify unified schemas still valid
4. Update verification reports

---

## Risk Assessment

| Risk | Mitigation |
|------|------------|
| Data loss | Backup created before any deletions |
| Broken links | Symlinks tested before removing originals |
| Schema validation failure | Re-run verification after cleanup |
| Missing Open Targets data | Extract to new location before cleanup |

---

## Rollback Plan

If issues are discovered after cleanup:

```bash
# Restore from backup
cp -r _legacy_backup/3.3.pharmacogenomics/ 03.diseases.phenotypes/
cp -r _legacy_backup/3.4.drug.disease.targets/ 03.diseases.phenotypes/
```

---

## Success Criteria

- [ ] All unique content from legacy directories preserved
- [ ] Cross-reference directories created with proper symlinks
- [ ] Primary data sources have complete, merged documentation
- [ ] Open Targets Platform has its own data source location
- [ ] Schema validation passes
- [ ] No broken links in documentation

---

## Files Affected

### Files to Merge
| Source (Legacy) | Target (Primary) |
|-----------------|------------------|
| `3.3.pharmacogenomics/pharmgkb/schema.md` | `1.4.pharmacogenomics/pharmgkb/schema.md` |
| `3.4.drug.disease.targets/dgidb/schema.md` | `2.7.compound.target.interactions/dgidb/schema.md` |

### Files to Create
| File | Purpose |
|------|---------|
| `3.3.pharmacogenomics/_index.md` | Cross-reference index |
| `3.4.drug.disease.targets/_index.md` | Cross-reference index |
| `open.targets/_index.md` | New data source (extracted) |
| `open.targets/schema.md` | Open Targets schema |
| `open.targets/download.md` | Open Targets download info |

### Directories to Convert
| Directory | Conversion |
|-----------|------------|
| `3.3.pharmacogenomics/` | Data dir → Cross-reference with symlinks |
| `3.4.drug.disease.targets/` | Data dir → Cross-reference with symlinks |

---

*Plan created: 2026-01-24*
*Awaiting approval before execution*
