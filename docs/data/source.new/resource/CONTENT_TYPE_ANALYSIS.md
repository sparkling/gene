---
title: "Content Type Analysis - Template Fit Assessment"
date: 2026-01-23
status: "Migration Complete"
version: "1.1"
migration_date: 2026-01-23
---

# Content Type Analysis

## Migration Status: COMPLETE

**Migration completed:** 2026-01-23

### Summary of Changes

| Action | Status | Details |
|--------|--------|---------|
| Create `reference/` folder | Complete | Formats and architecture subdirectories created |
| Create `guides/` folder | Complete | Wikidata (7 guides) and Integration (3 guides) |
| Reorganize `operations/` | Complete | methodology/, governance/, planning/ populated |
| Preserve `domains/` | Complete | 13 health domain files unchanged |
| Update index files | Complete | All subcategory _index.md files created |

---

## Executive Summary

Analysis of remaining content in `databases/`, `domains/`, and `operations/` reveals **6 distinct document types**. Only **1 type** fully fits the current template framework. The others represent valuable documentation that needs either new templates or a separate documentation structure.

| Content Type | Fits Templates? | Files | Recommendation |
|--------------|-----------------|-------|----------------|
| Data Source Documentation | ✅ YES | ~40 | Already migrated |
| Multi-Database Overviews | ⚠️ PARTIAL | 12 | Adapt to subcategory _index.md |
| Health Domain Views | ❌ NO | 13 | New template or separate folder |
| Technical Format References | ❌ NO | 8 | New template or operations/ |
| Integration Guides | ❌ NO | 10 | New template or guides/ folder |
| Methodology Guides | ❌ NO | 7 | Separate operations/ folder |

---

## Current Template Framework

The `resource/` templates support **individual data source documentation**:

| Template | Purpose | Structure |
|----------|---------|-----------|
| `_index.md` | Primary data source overview | One source, statistics, access methods |
| `schema.md` | Technical schema details | Tables, fields, API endpoints |
| `download.md` | Acquisition guide | Download commands, file inventory |
| `xrefs.md` | Cross-references | Taxonomy locations, ID mappings |
| Subcategory `_index.md` | Category index | List of sources in subcategory |

**Key limitation:** Templates assume **1 document = 1 data source**

---

## Content Type 1: Multi-Database Overviews

### Description
Documents that aggregate multiple databases by functional area, providing comparative analysis and selection guidance.

### Examples (12 files)

| File | Databases Covered |
|------|-------------------|
| `databases/genetics/primary.md` | dbSNP, ClinVar, gnomAD, dbNSFP, AlphaMissense, CADD, SpliceAI (30+) |
| `databases/genetics/population.md` | TOPMed, All of Us, ALFA, UK Biobank, GenomeAsia (6) |
| `databases/compounds/pharmaceuticals.md` | PharmGKB, CPIC, DrugBank, ChEMBL, PubChem (13) |
| `databases/compounds/natural-products.md` | COCONUT, LOTUS, NPASS, NPAtlas, Dr. Duke's (5) |
| `databases/compounds/drug-metabolism.md` | PharmVar, SuperCYP, KEGG Drug, HMDB (18) |
| `databases/pathways/primary.md` | Reactome, WikiPathways, KEGG, MetaCyc (10) |
| `databases/pathways/disease.md` | DisGeNET, OMIM, HPO, MONDO, Orphanet (12) |
| `databases/pathways/processes.md` | Gene Ontology, Rhea, BRENDA, BioModels (8) |
| `databases/traditional/tcm.md` | BATMAN-TCM, HERB, SymMap, ETCM (23) |
| `databases/traditional/ayurveda.md` | IMPPAT, Ayurvedic databases |
| `databases/traditional/kampo.md` | KampoDB, Japanese herbal |
| `databases/traditional/western-herbal.md` | DSLD, Dr. Duke's, EMA Herbal (25) |

### Template Fit: ⚠️ PARTIAL

**Problem:** These are **1-to-many** documents (one doc → many databases), but templates are **1-to-1**.

**Solution Options:**
1. **Adapt as enhanced subcategory _index.md** - Add comparative tables, selection guidance
2. **Create new "overview.md" template** - Multi-database comparison document
3. **Keep separate** - Retain in databases/ as reference views

### Recommendation: **Option 1 - Enhance subcategory _index.md**

Extend the subcategory _index.md template to include:
- Comparative decision tables
- Database selection criteria
- Combined statistics
- Integration recommendations

---

## Content Type 2: Health Domain Views

### Description
Cross-cutting documents that aggregate databases by **clinical application area** rather than data type. These don't map to the taxonomy structure.

### Examples (13 files)

| File | Clinical Focus | Databases Referenced |
|------|----------------|---------------------|
| `domains/mental-cognitive.md` | Psychiatric, cognitive, nootropics | PGC, Allen Brain, GTEx, ChEMBL (35+) |
| `domains/cardio-metabolic.md` | Cardiovascular, metabolic | CARDIoGRAMplusC4D, GLGC, DIAGRAM (10) |
| `domains/cancer-oncology.md` | Cancer genetics, oncology | COSMIC, TCGA, OncoKB, CIViC (8) |
| `domains/autoimmune.md` | Autoimmune, hormonal | IPD-IMGT/HLA, ImmunoBase (23) |
| `domains/rare.md` | Rare diseases, biobanks | Orphanet, DECIPHER, UK Biobank (9) |
| `domains/womens-pediatric.md` | Women's health, pediatric | Reproductive, pregnancy, pediatric PGx (30+) |
| `domains/microbiome.md` | Gut, oral, skin microbiome | HMP, GMrepo, HOMD, MASI (12) |
| `domains/allergy-pain.md` | Allergy, histamine, pain | GWAS Catalog allergy, pain genetics (26) |
| `domains/oral-skin-sensory.md` | Oral, skin, vision, hearing | HOMD, RetNet, BitterDB (20+) |
| `domains/sleep-longevity-nutri.md` | Sleep, aging, nutrigenomics | CircaDB, HAGR, NutriGenomeDB (30+) |
| `domains/clinical-environmental-mito.md` | Toxicology, mitochondrial | CTD, T3DB, MITOMAP (10+) |
| `domains/clinical-biomarkers-labs.md` | Biomarkers, lab references | MarkerDB, LOINC, CALIPER (21) |
| `domains/community-patient-networks.md` | Patient communities, biohacking | PatientsLikeMe, Open Humans (30+) |

### Template Fit: ❌ NO

**Problem:** These cross-cut the taxonomy. A single domain file references databases from 5+ different resource/ categories.

**Why they're valuable:**
- Provide clinical use-case perspective
- Guide database selection for specific health applications
- Document databases NOT YET in resource/ (longevity, toxicology, biomarkers)

### Recommendation: **Keep separate as `domains/` folder**

These serve a different purpose than data source documentation:
- **resource/** = "What is database X?"
- **domains/** = "What databases do I need for clinical application Y?"

Create symlinks or cross-references between them.

---

## Content Type 3: Technical Format References

### Description
Multi-format technical specifications that document data interchange standards, not individual databases.

### Examples (8 files)

| File | Content |
|------|---------|
| `operations/schemas/pathway-formats.md` | BioPAX, SBML, PSI-MI, CX format specs |
| `databases/literature/data-structures.md` | PubMed XML, Europe PMC JSON, OpenAlex |
| `databases/literature/abstracts-vs-fulltext.md` | Text extraction methodologies |
| `operations/schemas/ruvector-three-worlds-schema.md` | Architecture design document |
| `operations/schemas/unified-schema-analysis.md` | Cross-database schema comparison |
| `operations/schemas/sample-data.md` | Example data records |
| `operations/schemas/schemas-navigation.md` | Schema index/navigation |
| `operations/schemas/schemas-index.md` | Schema catalog |

### Template Fit: ❌ NO

**Problem:** These document **formats and standards**, not data sources.

### Recommendation: **Create `reference/` folder**

New structure:
```
docs/data/source.new/reference/
├── formats/
│   ├── pathway-formats.md      (BioPAX, SBML, GPML)
│   ├── variant-formats.md      (VCF, GFF, BED)
│   └── literature-formats.md   (PubMed XML, OpenAlex JSON)
├── architecture/
│   ├── unified-schema.md
│   └── three-worlds.md
└── examples/
    └── sample-data.md
```

---

## Content Type 4: Integration Guides

### Description
Step-by-step guides for using a data source for specific purposes (e.g., "How to use Wikidata for pharmaceutical data").

### Examples (10 files)

| File | Purpose |
|------|---------|
| `operations/wikidata/wikidata-master-reference.md` | Complete Wikidata usage guide |
| `operations/wikidata/wikidata-pharma.md` | Pharmaceutical queries |
| `operations/wikidata/wikidata-traditional.md` | Traditional medicine queries |
| `operations/wikidata/wikidata-supplements.md` | Supplement queries |
| `operations/wikidata/wikidata-bulk.md` | Bulk data download |
| `operations/wikidata/wikipedia-wikidata.md` | Wikipedia integration |
| `operations/integration-guide.md` | General integration patterns |
| `operations/compound-pathway-linking.md` | Cross-database linking |
| `operations/pathway-target-mapping.md` | Pathway-target relationships |
| `operations/pathways-targets.md` | Target identification |

### Template Fit: ❌ NO

**Problem:** These are **how-to guides**, not data source documentation.

### Recommendation: **Create `guides/` folder or integrate into resource/**

Option A: Separate guides folder
```
docs/data/source.new/guides/
├── wikidata/
│   ├── overview.md
│   ├── pharma-queries.md
│   └── bulk-download.md
└── integration/
    ├── compound-pathway-linking.md
    └── cross-database-mapping.md
```

Option B: Add `guide.md` template to resource/ data sources
```
resource/08.literature.knowledge/8.2.knowledge.bases/wikidata/
├── _index.md
├── schema.md
├── download.md
├── xrefs.md
└── guide.md  ← NEW: Usage guides
```

---

## Content Type 5: Methodology Guides

### Description
Operational methodology documents covering data curation, pipeline design, and governance.

### Examples (7 files)

| File | Content |
|------|---------|
| `databases/literature/pipeline-design.md` | Literature pipeline architecture |
| `databases/literature/coverage-analysis.md` | Coverage assessment methodology |
| `databases/literature/public-sources.md` | Public source evaluation |
| `databases/literature/sources.md` | Source catalog |
| `operations/curation-framework.md` | Data curation methodology |
| `operations/data-access-legal.md` | Legal/licensing analysis |
| `operations/size-estimates.md` | Storage estimates |

### Template Fit: ❌ NO

**Problem:** These are **operational documents**, not data source documentation.

### Recommendation: **Keep in `operations/` folder**

These belong in a separate operational documentation structure:
```
docs/data/source.new/operations/
├── methodology/
│   ├── curation-framework.md
│   ├── pipeline-design.md
│   └── coverage-analysis.md
├── governance/
│   ├── data-access-legal.md
│   └── licensing.md
└── planning/
    ├── size-estimates.md
    └── roadmap.md
```

---

## Content Type 6: Navigation/Index Files

### Description
Internal navigation files that provide directory indexes.

### Examples
- `databases/_index.md`
- `domains/_index.md`
- `operations/_index.md`
- `operations/schemas/_index.md`
- etc.

### Template Fit: N/A

**Status:** These are replaced by the resource/ folder structure and don't need migration.

---

## Summary: What Fits vs What Doesn't

### ✅ FITS Current Templates (Already Migrated)

| Content | Template | Status |
|---------|----------|--------|
| Individual database overviews | `_index.md` | ✅ 190 files in resource/ |
| Technical schemas | `schema.md` | ✅ 40 files migrated |
| Download instructions | `download.md` | ✅ 30 files created |
| Cross-references | `xrefs.md` | ✅ 30 files created |

### ⚠️ PARTIALLY FITS (Needs Template Extension)

| Content | Current Template | Enhancement Needed |
|---------|------------------|-------------------|
| Multi-database overviews | Subcategory `_index.md` | Add comparison tables, decision guidance |
| Integration guides | None | Add `guide.md` template |

### ❌ DOES NOT FIT (Needs Separate Structure)

| Content | Recommended Location | Reason |
|---------|---------------------|--------|
| Health domain views | `domains/` (keep) | Cross-cuts taxonomy |
| Technical format references | `reference/formats/` | Documents standards, not sources |
| Architecture documents | `reference/architecture/` | System design docs |
| Methodology guides | `operations/methodology/` | Operational docs |
| Governance/legal | `operations/governance/` | Policy docs |

---

## Recommended Final Structure

```
docs/data/source.new/
├── resource/                    # Data source documentation (current)
│   ├── 01.genetics.genomics/
│   ├── 02.compounds.molecules/
│   ├── ...
│   └── 09.microbiome/
│
├── domains/                     # Health domain views (keep as-is)
│   ├── mental-cognitive.md
│   ├── cardio-metabolic.md
│   └── ...
│
├── reference/                   # Technical references (new)
│   ├── formats/
│   │   ├── pathway-formats.md
│   │   └── variant-formats.md
│   └── architecture/
│       └── unified-schema.md
│
├── guides/                      # Integration guides (new)
│   ├── wikidata/
│   └── cross-database/
│
└── operations/                  # Operational docs (reorganize)
    ├── methodology/
    ├── governance/
    └── planning/
```

---

## Action Items

### Immediate (Fits Templates)
1. Enhance subcategory `_index.md` files with content from multi-database overviews
2. Consider adding `guide.md` template for integration documentation

### Retain Separately
3. Keep `domains/` folder - provides valuable clinical perspective
4. Reorganize `operations/` into methodology/governance/planning

### Create New Structures
5. Create `reference/` folder for format specifications
6. Create `guides/` folder for integration guides

### Document Gaps
7. Many databases in domain files don't have resource/ entries yet (~100 databases)
8. Consider creating entries for high-priority missing databases

---

## Conclusion

The current template framework is well-suited for **individual data source documentation** but doesn't accommodate:
- Cross-cutting health domain views
- Multi-database comparison documents
- Technical format specifications
- Integration how-to guides
- Operational methodology

**Recommendation:** Keep the template framework focused on data sources, and maintain separate folder structures for other content types. This preserves the clarity of resource/ while retaining valuable documentation in appropriate locations.

---

## Generated: 2026-01-23
