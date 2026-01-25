# Compounds and Molecules - Template Conformance Verification

**Date:** 2026-01-23
**Category:** 02.compounds.molecules
**Total Data Sources:** 23
**Conformance Rate:** 13.0% (3 fully compliant)

## Executive Summary

The compounds.molecules resource category shows **significant conformance gaps** across template requirements. Only 3 of 23 data sources are fully compliant with the required template structure.

### Key Metrics

| Metric | Value | Status |
|--------|-------|--------|
| **Fully Compliant Sources** | 3/23 | Critical Gap |
| **schema.md Coverage** | 5/23 (21.7%) | Major Gap |
| **download.md Coverage** | 4/23 (17.4%) | Major Gap |
| **xrefs.md Coverage** | 7/23 (30.4%) | Moderate Gap |
| **README.md Section Gaps** | 3 sources | Minor Gap |

---

## Template Requirements

### 1. README.md (Data Source Files)

**Required Frontmatter:**
- `id` - Unique identifier
- `title` - Full name
- `type` - Must be "data-source"
- `category` - Primary category
- `subcategory` - Secondary category
- `parent` - Parent category link
- `tier` - Data tier (1-3)
- `last_updated` - Update date
- `status` - Current status
- `tags` - Searchable tags

**Required Sections (## Headers):**
1. Overview (2-3 paragraphs)
2. Primary Use Cases (numbered list)
3. Key Identifiers (table format)
4. Access Methods (table with URL/notes)
5. License (table with usage rights)

### 2. schema.md (Database Schema)

**Required Sections:**
1. Database Statistics (metrics table)
2. Entity Relationship Overview (ER diagram or description)
3. Core Tables/Entities (entity descriptions)
4. Data Formats (format specifications)
5. Sample Record (example with all fields)

### 3. download.md (Download Instructions)

**Required Sections:**
1. Quick Start (2-3 step guide)
2. Prerequisites (system requirements)
3. Download Methods (multiple options)
4. File Inventory (list of available files)
5. Verification (checksums, validation)

### 4. xrefs.md (Optional but Recommended)

Cross-reference mapping to related data sources and integration points.

---

## Conformance by Category

### Natural Products (2.1)

| Source | Files | Compliance | Issues |
|--------|-------|-----------|--------|
| LOTUS | README.md, schema.md | 83% | Missing download.md |
| COCONUT | README.md, schema.md, xrefs.md | 83% | Missing download.md |
| Dr. Duke's | README.md, xrefs.md | 40% | Missing schema.md, download.md |
| NPASS | README.md | 20% | Missing schema.md, download.md, xrefs.md |
| NPAtlas | README.md | 20% | Missing schema.md, download.md, xrefs.md |

### Pharmaceuticals (2.2)

| Source | Files | Compliance | Issues |
|--------|-------|-----------|--------|
| ChEMBL | README.md, schema.md, download.md, xrefs.md | 100% | None |
| DrugBank | README.md, download.md | 67% | Missing schema.md |
| RxNorm | README.md | 20% | Missing schema.md, download.md, xrefs.md |
| DailyMed | README.md, xrefs.md | 40% | Missing schema.md, download.md |
| Orange Book | README.md | 20% | Missing schema.md, download.md, xrefs.md |

### Food Compounds & Nutrients (2.4)

| Source | Files | Compliance | Issues |
|--------|-------|-----------|--------|
| Phenol-Explorer | README.md, xrefs.md | 40% | Missing schema.md, download.md |
| PhytoHub | README.md, xrefs.md | 40% | Missing schema.md, download.md |
| USDA FoodData | README.md, xrefs.md | 40% | Missing schema.md, download.md |

### Drug Metabolism & Pharmacokinetics (2.5)

| Source | Files | Compliance | Issues |
|--------|-------|-----------|--------|
| SuperCYP | README.md | 20% | Missing schema.md, download.md, xrefs.md |
| SwissADME | README.md | 0% | Missing Key Identifiers section; Missing schema.md, download.md, xrefs.md |

### Chemical Ontology & Classification (2.6)

| Source | Files | Compliance | Issues |
|--------|-------|-----------|--------|
| ChEBI | README.md, schema.md, download.md | 100% | xrefs.md optional |
| PubChem | README.md, schema.md, download.md | 100% | xrefs.md optional |
| ClassyFire | README.md | 0% | Missing Key Identifiers section; Missing schema.md, download.md, xrefs.md |
| NPClassifier | README.md | 0% | Missing Key Identifiers section; Missing schema.md, download.md, xrefs.md |

### Compound-Target Interactions (2.7)

| Source | Files | Compliance | Issues |
|--------|-------|-----------|--------|
| BindingDB | README.md | 20% | Missing schema.md, download.md, xrefs.md |
| DGIdb | README.md, xrefs.md | 40% | Missing schema.md, download.md |
| GtoPdb | README.md | 20% | Missing schema.md, download.md, xrefs.md |
| TTD | README.md | 20% | Missing schema.md, download.md, xrefs.md |

---

## Fully Compliant Sources (3)

### 1. ChEMBL (2.2.pharmaceuticals/chembl)
- **Files:** README.md, schema.md, download.md, xrefs.md
- **Status:** 100% Compliant
- **All Sections:** Complete
- **Notes:** Exemplary documentation template

### 2. ChEBI (2.6.chemical.ontology.classification/chebi)
- **Files:** README.md, schema.md, download.md
- **Status:** 100% Compliant (xrefs.md optional)
- **All Sections:** Complete
- **Notes:** Schema documentation is thorough

### 3. PubChem (2.6.chemical.ontology.classification/pubchem)
- **Files:** README.md, schema.md, download.md
- **Status:** 100% Compliant (xrefs.md optional)
- **All Sections:** Complete
- **Notes:** Well-structured download instructions

---

## Critical Issues (0% Compliance)

### 3 Sources with Zero Compliance

1. **SwissADME** (2.5.drug.metabolism.pharmacokinetics/swissadme)
   - Missing: Key Identifiers section in README.md
   - Missing: schema.md, download.md, xrefs.md
   - Action: Complete README.md + create 3 supporting files

2. **ClassyFire** (2.6.chemical.ontology.classification/classyfire)
   - Missing: Key Identifiers section in README.md
   - Missing: schema.md, download.md, xrefs.md
   - Action: Complete README.md + create 3 supporting files

3. **NPClassifier** (2.6.chemical.ontology.classification/npclassifier)
   - Missing: Key Identifiers section in README.md
   - Missing: schema.md, download.md, xrefs.md
   - Action: Complete README.md + create 3 supporting files

---

## High Priority Issues (20 Sources)

### Missing Both schema.md AND download.md (15 sources)

Sources missing critical supplementary documentation:

- 2.1.natural.products/dr.dukes
- 2.1.natural.products/npass
- 2.1.natural.products/npatlas
- 2.2.pharmaceuticals/dailymed
- 2.2.pharmaceuticals/orange.book
- 2.2.pharmaceuticals/rxnorm
- 2.4.food.compounds.nutrients/phenol.explorer
- 2.4.food.compounds.nutrients/phytohub
- 2.4.food.compounds.nutrients/usda.fooddata
- 2.5.drug.metabolism.pharmacokinetics/supercyp
- 2.7.compound.target.interactions/bindingdb
- 2.7.compound.target.interactions/dgidb
- 2.7.compound.target.interactions/gtopdb
- 2.7.compound.target.interactions/ttd

**Total Files Needed:** 30 (15 schema.md + 15 download.md)

### Missing Only schema.md (3 sources)

- 2.2.pharmaceuticals/drugbank
- 2.1.natural.products/coconut (has schema, missing download)
- 2.1.natural.products/lotus (has schema, missing download)

---

## Medium Priority Issues (2 Sources)

### Missing Only xrefs.md

These sources are otherwise compliant but lack cross-reference mapping:

- 2.2.pharmaceuticals/chembl
- 2.6.chemical.ontology.classification/pubchem

---

## Recommendations

### Immediate Actions (Priority 1)

1. **Create missing README.md sections (3 files)**
   - Add "Key Identifiers" section to:
     - SwissADME
     - ClassyFire
     - NPClassifier
   - Include table with identifier patterns and examples

2. **Create 18 missing schema.md files**
   - Templates available: See ChEMBL, ChEBI, PubChem examples
   - Minimum requirements:
     - Database Statistics table
     - Core table/entity descriptions
     - Sample record with field documentation

3. **Create 19 missing download.md files**
   - Templates available: See ChEMBL and PubChem examples
   - Minimum requirements:
     - Download URLs/methods
     - File size and format information
     - Verification checksums

### Medium-term Actions (Priority 2)

4. **Add 16 missing xrefs.md files**
   - Cross-reference mappings to related sources
   - Integration point documentation
   - Use existing xrefs.md files as templates

5. **Establish CI/CD validation**
   - Add pre-commit hooks to validate template conformance
   - Require all data sources to pass conformance check
   - Generate automated reports on each commit

### Long-term Actions (Priority 3)

6. **Documentation standardization**
   - Create detailed template guidelines
   - Establish naming conventions
   - Document schema design patterns

7. **Automated schema generation**
   - Build tools to extract schema from source documentation
   - Generate schema.md from database introspection
   - Semi-automated download.md generation from API docs

---

## File Locations

- **Report Location:** `/home/claude/src/gene/docs/data/source.new/resource/02.compounds.molecules/`
- **JSON Report:** `CONFORMANCE_REPORT.json`
- **This Report:** `CONFORMANCE_SUMMARY.md`

---

## Template Examples

Reference these fully compliant sources for template structure:

1. **ChEMBL** (complete example)
   - Path: `2.2.pharmaceuticals/chembl/`
   - All files present and properly formatted

2. **ChEBI** (schema + download example)
   - Path: `2.6.chemical.ontology.classification/chebi/`
   - Well-documented database structure

3. **PubChem** (large dataset example)
   - Path: `2.6.chemical.ontology.classification/pubchem/`
   - Comprehensive download documentation

---

## Contact & Issues

For template conformance questions or to report issues, refer to the project documentation guidelines.
