# Content Type Taxonomy

**Document ID:** CONTENT-TYPE-TAXONOMY
**Status:** Final Analysis
**Date:** January 2026
**Purpose:** Classification of document types in source.new/ to guide migration decisions

---

## Executive Summary

Analysis of 423 markdown files across `/home/claude/src/gene/docs/data/source.new/` reveals **15 distinct document types**. The current template system in `resource/` supports 4 primary types well. This taxonomy identifies which documents fit templates, which need adaptation, and which should remain in their current locations.

### Document Type Distribution

| Type | Count | Template Support | Migration Status |
|------|-------|------------------|------------------|
| `data-source` | 127 | Full | Migrate to resource/ |
| `subcategory` | 54 | Full | Already in resource/ |
| `schema` | 53 | Full | Keep in operations/schemas/ |
| `xrefs` | 30 | Full | Migrate to resource/ |
| `download` | 30 | Full | Migrate to resource/ |
| `health-domain` | 10 | None | Keep in domains/ |
| `integration` | 10 | None | Keep in operations/integration/ |
| `category` | 9 | Full | Already in resource/ |
| `download-guide` | 5 | Partial | Keep in operations/downloads/ |
| `research` | 4 | None | Keep in research/ |
| `research-priority` | 2 | None | Keep in research/ |
| `governance` | 2 | None | Keep in operations/governance/ |
| `clinical` | 2 | None | Evaluate for new template |
| `community` | 1 | None | Keep in domains/ |

---

## Category 1: Data Source Documentation (Fits Templates)

### Type 1.1: data-source (_index.md)

**Definition:** Primary documentation for an individual database or data source.

**Example Files:**
- `/resource/01.genetics.genomics/1.1.variant.repositories/clinvar/_index.md`
- `/resource/05.traditional.medicine/5.1.tcm.databases/batman-tcm/_index.md`
- `/resource/06.nutrition.food/6.1.food.composition/foodb/_index.md`

**Template:** `TEMPLATES.md` Template 1 (_index.md)

**Required Sections:**
- Frontmatter with id, title, type, category, subcategory, tier, status, tags
- Overview (2-3 paragraphs)
- Key Statistics table
- Primary Use Cases
- Key Identifiers table
- Access Methods table
- License table
- See Also links

**Migration Recommendation:** **MIGRATE TO resource/** - These are the core content type, fully supported by templates.

---

### Type 1.2: schema (Technical Schema Documentation)

**Definition:** Detailed technical documentation of database schemas, file formats, and data structures.

**Example Files:**
- `/operations/schemas/clinvar-schema.md` (878 lines, comprehensive)
- `/operations/schemas/chembl-schema.md`
- `/operations/schemas/pathway-formats.md`

**Template:** `TEMPLATES.md` Template 2 (schema.md)

**Required Sections:**
- Database Statistics
- Entity Relationship Overview (ASCII diagram)
- Core Tables/Entities with field definitions
- API Endpoints
- Data Formats
- Sample Record (JSON)
- Glossary
- References

**Current Location:** `/operations/schemas/`

**Migration Recommendation:** **KEEP IN operations/schemas/** - These are operational reference documents, not data source descriptions. The centralized schema library serves as a technical reference that spans multiple data sources.

**Note:** Some data sources in resource/ have inline schema documentation. Consider creating symlinks or cross-references rather than duplicating content.

---

### Type 1.3: download (Download Instructions)

**Definition:** Step-by-step instructions for acquiring data from a source.

**Example Files:**
- `/resource/01.genetics.genomics/1.3.population.genetics/gnomad/download.md`
- `/resource/06.nutrition.food/6.1.food.composition/foodb/download.md`
- `/resource/07.proteins.molecular.biology/7.2.protein.structures/pdb/download.md`

**Template:** `TEMPLATES.md` Template 3 (download.md)

**Required Sections:**
- Quick Start (minimal command)
- Prerequisites
- Download Methods (Primary/Alternative)
- File Inventory table
- Post-Download Processing
- Verification commands
- Update Schedule

**Migration Recommendation:** **MIGRATE TO resource/** - Place alongside corresponding _index.md files in the source directory.

---

### Type 1.4: xrefs (Cross-References)

**Definition:** Documentation of how a source integrates with other databases, ID mappings, and cross-reference capabilities.

**Example Files:**
- `/resource/01.genetics.genomics/1.1.variant.repositories/clinvar/xrefs.md`
- `/resource/05.traditional.medicine/5.1.traditional.chinese.medicine/batman.tcm/xrefs.md`
- `/resource/06.nutrition.food/6.1.food.composition/foodb/xrefs.md`

**Template:** `TEMPLATES.md` Template 4 (xrefs.md)

**Required Sections:**
- Taxonomy Locations (if polyhierarchical)
- External ID Mappings table
- Integration Notes

**Migration Recommendation:** **MIGRATE TO resource/** - Place alongside corresponding _index.md files.

---

### Type 1.5: subcategory (_index.md for folders)

**Definition:** Navigation index for a subcategory containing multiple data sources.

**Example Files:**
- `/resource/01.genetics.genomics/1.1.variant.repositories/_index.md`
- `/resource/05.traditional.medicine/5.1.traditional.chinese.medicine/_index.md`

**Template:** `TEMPLATES.md` Template 5 (Subcategory _index.md)

**Required Sections:**
- Parent link
- Overview
- Data Sources table (Source | Tier | Description)
- Integration Notes

**Migration Recommendation:** **ALREADY IN resource/** - These are correctly placed.

---

### Type 1.6: category (_index.md for top-level)

**Definition:** Navigation index for a top-level category.

**Example Files:**
- `/resource/01.genetics.genomics/_index.md`
- `/resource/05.traditional.medicine/_index.md`

**Migration Recommendation:** **ALREADY IN resource/** - These are correctly placed.

---

## Category 2: Health Domain Views (Does NOT Fit Templates)

### Type 2.1: health-domain

**Definition:** Cross-cutting views of data sources organized by health/disease domain rather than by database. These aggregate information from multiple sources around a clinical theme.

**Example Files:**
- `/domains/mental-cognitive.md` (985 lines)
- `/domains/cardio-metabolic.md`
- `/domains/autoimmune.md`
- `/domains/cancer-oncology.md`
- `/domains/microbiome.md`
- `/domains/womens-pediatric.md`
- `/domains/allergy-pain.md`
- `/domains/rare.md`

**Characteristics:**
- Multi-database aggregation (35+ databases in mental-cognitive.md)
- Domain-specific organization (not by data type)
- Storage estimates for domain-specific data
- Priority tiering for each domain
- Download instructions for multiple sources
- Integration recommendations specific to domain

**Why Templates Don't Fit:**
1. **Scope:** Templates are per-database; these span dozens of databases
2. **Organization:** Templates follow database taxonomy; these follow clinical taxonomy
3. **Purpose:** Templates describe "what is this database"; these describe "what data exists for this health domain"
4. **Audience:** Templates serve data engineers; these serve domain scientists and product managers

**Migration Recommendation:** **KEEP IN domains/** - Create new template type for health-domain documents.

**Proposed New Template:** Create `TEMPLATE-HEALTH-DOMAIN.md` with sections:
- Domain Overview
- Key Decisions table
- Database Catalog (organized by subdomain)
- Integration Recommendations (Priority 1/2/3)
- API Integration Summary
- Download section
- Schema (Core Entity Fields, Relationships)
- Data Set Size
- License Summary
- Glossary (Domain-Specific Terms, Acronyms)

---

### Type 2.2: community (Patient Networks)

**Definition:** Documentation of patient community and citizen science data sources.

**Example Files:**
- `/domains/community-patient-networks.md`

**Migration Recommendation:** **KEEP IN domains/** - Similar structure to health-domain but focused on community resources.

---

## Category 3: Integration Methodologies (Does NOT Fit Templates)

### Type 3.1: integration

**Definition:** Technical guides for connecting multiple data sources, including API chaining, identifier mapping, and ETL pipelines.

**Example Files:**
- `/operations/integration/integration-guide.md` (700 lines)
- `/operations/integration/compound-pathway-linking.md`
- `/operations/integration/pathway-target-mapping.md`
- `/operations/integration/wikidata-master-reference.md`
- `/operations/integration/xrefs.md`

**Characteristics:**
- Multi-database focus
- API endpoint documentation
- Sample code (Python, SQL, cURL)
- Data flow diagrams
- Schema recommendations
- Licensing summary for multiple sources

**Why Templates Don't Fit:**
1. **Scope:** Integration spans multiple databases
2. **Content:** Heavy on code examples and API patterns
3. **Purpose:** How to combine data, not how to access it

**Migration Recommendation:** **KEEP IN operations/integration/** - These are operational guides, not data source documentation.

---

### Type 3.2: download-guide

**Definition:** Processing pipeline documentation for bulk data operations.

**Example Files:**
- `/operations/downloads/processing-pipeline.md`
- `/operations/downloads/traditional-medicine.md`
- `/operations/downloads/pharmaceuticals.md`
- `/operations/downloads/pathways-targets.md`

**Characteristics:**
- Storage requirement calculations
- Streaming parser code
- Memory-efficient processing patterns
- Parallel processing strategies
- Database loading scripts

**Why Templates Don't Fit:**
1. **Scope:** Cross-database processing concerns
2. **Content:** Infrastructure and pipeline code
3. **Audience:** DevOps and data engineers

**Migration Recommendation:** **KEEP IN operations/downloads/** - These are operational infrastructure guides.

---

## Category 4: Format Specifications (Partial Template Fit)

### Type 4.1: Format Schema Documents

**Definition:** Technical specifications for data exchange formats used across multiple databases.

**Example Files:**
- `/operations/schemas/pathway-formats.md` (GPML, BioPAX, SBML, PSI-MI, CX)
- `/operations/schemas/kgml-schema.md`
- `/operations/schemas/wikipathways-gpml-schema.md`

**Characteristics:**
- Format-centric rather than database-centric
- XML/JSON schema definitions
- Parsing code examples
- Format comparison tables

**Why Partial Fit:**
- Schema template works for structure
- Missing: format comparison, interoperability notes

**Migration Recommendation:** **KEEP IN operations/schemas/** - Add format-specific template variant if needed.

---

## Category 5: Research and Analysis (Does NOT Fit Templates)

### Type 5.1: research

**Definition:** Research analysis documents, often synthesizing findings from swarm agents or deep analysis sessions.

**Example Files:**
- `/research/literature-priority.md`
- `/research/interventions-priority.md`
- `/research/schema-gaps.md`
- `/research/world1-schema-research.md`
- `/research/genetics-schema-research.md`

**Characteristics:**
- Decision recommendations
- Swarm synthesis outputs
- Cost-benefit analysis
- Timeline and phasing plans
- Risk assessments

**Why Templates Don't Fit:**
1. **Purpose:** Strategic decision documents, not reference documentation
2. **Lifecycle:** Point-in-time analysis, not maintained reference
3. **Structure:** Narrative and analysis, not structured reference

**Migration Recommendation:** **KEEP IN research/** - These are project artifacts, not operational documentation.

---

### Type 5.2: research-priority

**Definition:** Final recommendations synthesized from research.

**Example Files:**
- `/research/literature-priority.md`
- `/research/interventions-priority.md`

**Migration Recommendation:** **KEEP IN research/** - Decision documents for project planning.

---

## Category 6: Governance (Does NOT Fit Templates)

### Type 6.1: governance

**Definition:** Legal, licensing, and data access governance documentation.

**Example Files:**
- `/operations/governance/data-access-legal.md`
- `/operations/governance/curation-framework.md`

**Characteristics:**
- License compatibility analysis
- Data use agreement requirements
- Ethical considerations
- Curation workflows

**Migration Recommendation:** **KEEP IN operations/governance/** - Operational policy documents.

---

## Category 7: Literature Pipeline (Does NOT Fit Templates)

### Type 7.1: literature (Pipeline Documentation)

**Definition:** Technical documentation for the literature/research paper processing pipeline.

**Example Files:**
- `/databases/literature/pipeline-design.md`
- `/databases/literature/data-structures.md`
- `/databases/literature/public-sources.md`
- `/databases/literature/abstracts-vs-fulltext.md`
- `/databases/literature/coverage-analysis.md`
- `/databases/literature/sources.md`

**Characteristics:**
- ETL pipeline specifications
- Embedding strategy documentation
- Storage calculations
- MeSH filtering logic
- Processing cost analysis

**Why Templates Don't Fit:**
1. **Scope:** Pipeline architecture, not database description
2. **Content:** System design documents
3. **Lifecycle:** Evolving implementation specs

**Migration Recommendation:** **EVALUATE** - Consider whether to:
- Keep in `/databases/literature/` as implementation docs
- Create data source entries for PubMed, PMC, etc. in resource/
- Both (reference docs + operational docs)

---

## Migration Decision Matrix

| Document Type | Count | Current Location | Action | Destination |
|---------------|-------|------------------|--------|-------------|
| data-source | 127 | Mixed | **Migrate** | resource/{category}/{subcategory}/{source}/ |
| schema | 53 | operations/schemas/ | **Keep** | operations/schemas/ |
| subcategory | 54 | resource/ | **Keep** | resource/ |
| xrefs | 30 | resource/ | **Keep** | resource/ |
| download | 30 | resource/ | **Keep** | resource/ |
| health-domain | 10 | domains/ | **Keep** | domains/ (create template) |
| integration | 10 | operations/integration/ | **Keep** | operations/integration/ |
| category | 9 | resource/ | **Keep** | resource/ |
| download-guide | 5 | operations/downloads/ | **Keep** | operations/downloads/ |
| research | 4 | research/ | **Keep** | research/ |
| research-priority | 2 | research/ | **Keep** | research/ |
| governance | 2 | operations/governance/ | **Keep** | operations/governance/ |
| clinical | 2 | databases/ | **Evaluate** | May need new template |
| community | 1 | domains/ | **Keep** | domains/ |

---

## Template Gap Analysis

### Templates That Exist and Work Well

1. **_index.md (data-source)** - Comprehensive, well-structured
2. **schema.md** - Detailed technical reference
3. **download.md** - Clear acquisition guide
4. **xrefs.md** - Cross-reference documentation
5. **subcategory _index.md** - Navigation structure

### Templates Needed

1. **health-domain.md** - For domain-view documents
   - Sections: Domain Overview, Database Catalog, Integration Recommendations, Combined Schema, Glossary

2. **format-spec.md** - For data format specifications
   - Sections: Format Overview, Schema Definition, Parsing Examples, Interoperability Notes

3. **integration-guide.md** - For multi-database integration
   - Sections: Use Case, Data Flow, API Chaining, Code Examples, Schema Recommendations

4. **CLAUDE.md (LLM-optimized)** - Already defined in TEMPLATE-GUIDE.md
   - Compact 50-100 line summaries for AI context windows

---

## Final Recommendations

### Immediate Actions

1. **Keep Current Structure for:**
   - `operations/schemas/` - Technical reference library
   - `operations/integration/` - Integration methodology
   - `operations/downloads/` - Processing pipelines
   - `operations/governance/` - Policy documents
   - `research/` - Analysis and decisions
   - `domains/` - Health domain views

2. **Complete Migration to resource/ for:**
   - All remaining data-source _index.md files
   - All download.md files paired with data sources
   - All xrefs.md files paired with data sources

3. **Create New Templates for:**
   - Health domain documents
   - Integration guides
   - Format specifications

### Directory Structure (Final)

```
source.new/
├── resource/                    # CANONICAL DATA SOURCE TAXONOMY
│   ├── 01.genetics.genomics/
│   │   ├── _index.md            # Category
│   │   └── 1.1.variant.repositories/
│   │       ├── _index.md        # Subcategory
│   │       └── clinvar/
│   │           ├── _index.md    # Data source
│   │           ├── download.md  # Download instructions
│   │           └── xrefs.md     # Cross-references
│   └── ...
├── domains/                     # HEALTH DOMAIN VIEWS (cross-cutting)
│   ├── mental-cognitive.md
│   ├── cardio-metabolic.md
│   └── ...
├── operations/                  # OPERATIONAL DOCUMENTATION
│   ├── schemas/                 # Technical schema library
│   ├── integration/             # Integration guides
│   ├── downloads/               # Processing pipelines
│   └── governance/              # Policy documents
├── research/                    # RESEARCH & ANALYSIS
│   ├── literature-priority.md
│   └── ...
└── databases/                   # LEGACY (to be deprecated)
    └── literature/              # Move PubMed/PMC to resource/
```

---

## Appendix: Document Type Signatures

### How to Identify Document Types

| Type | Frontmatter `type:` | Key Sections | File Location Pattern |
|------|---------------------|--------------|----------------------|
| data-source | `data-source` | Overview, Key Statistics, Key Identifiers | `*/_index.md` |
| schema | `schema` | Entity Relationship, Core Tables, Sample Record | `operations/schemas/*.md` |
| download | `download` | Quick Start, Download Methods, File Inventory | `*/download.md` |
| xrefs | `xrefs` | External ID Mappings, Integration Notes | `*/xrefs.md` |
| health-domain | `health-domain` | Database Catalog, Integration Recommendations | `domains/*.md` |
| integration | `integration` | Data Flow, API Examples, Schema Recommendations | `operations/integration/*.md` |
| research | `research` | Executive Summary, Recommendations, Risk Assessment | `research/*.md` |

---

*Taxonomy generated January 2026 for Gene platform data documentation migration.*
