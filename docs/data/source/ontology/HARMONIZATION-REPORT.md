# Ontology Harmonization Report

**Author:** Kurt Cagle, Lead Architect
**Date:** 2026-01-25
**Status:** Phase 1 Complete

---

## Executive Summary

This report documents the analysis and remediation of discrepancies between the Gene Data Source ontology definitions (OWL, SHACL, SKOS) and the instance data. The harmonization effort identified **7 major categories of issues** and implemented fixes for the most critical ones.

---

## Issues Identified

### 1. Category/Subcategory Title Formatting (FIXED)

**Problem:** Category titles in `datasources.ttl` used CamelCase without spaces (e.g., "GeneticsGenomics") instead of matching SKOS prefLabels ("Genetics and Genomics").

**Solution:** Updated all 9 category and 38 subcategory titles to use proper spacing and match SKOS taxonomy labels.

**Files Modified:**
- `/home/claude/src/gene/docs/data/source/ontology/instances/datasources.ttl`

### 2. Missing Notation Properties (FIXED)

**Problem:** SHACL shapes require `ds:notation` for categories (pattern "^[0-9]{2}$") and subcategories (pattern "^[0-9]+\.[0-9]+$"), but these were missing from all instances.

**Solution:** Added `ds:notation` to all categories ("01" through "09") and subcategories ("1.1", "1.2", etc.).

**Files Modified:**
- `/home/claude/src/gene/docs/data/source/ontology/instances/datasources.ttl`

### 3. OWL Ontology Property Gaps (FIXED)

**Problem:** Category-specific instance files use properties not defined in OWL:
- `ds:approximateRecordCount` - for DataSource approximate record counts
- `ds:updateFrequencyLabel` - for human-readable update frequency

**Solution:** Added missing datatype properties to OWL ontology.

**Files Modified:**
- `/home/claude/src/gene/docs/data/source/ontology/owl/datasource-ontology.ttl`

### 4. Property Namespace Conflicts (NEEDS PHASE 2)

**Problem:** Category-specific files (`01-genetics-genomics.ttl`, etc.) use properties differently than the OWL ontology specifies:

| Instance Usage | OWL Definition | Issue |
|---------------|----------------|-------|
| `ds:maintainer` (string) | `ds:maintainedBy` (object property) | Should reference `ds:Maintainer` instance |
| `ds:category tax:xxx` | `ds:belongsToCategory cat:xxx` | Wrong property and range |
| `ds:tier "1"^^xsd:integer` | `ds:tier tax:Tier1` | Should be object reference |
| `ds:status "active"` | `ds:status tax:Active` | Should be object reference |
| `ds:forSource` | Not defined | Inverse pattern needed |

**Recommendation:** Phase 2 should update category-specific instance files to use correct OWL properties.

### 5. Missing Required Properties in Main File (NEEDS PHASE 2)

**Problem:** DataSource instances in `datasources.ttl` are missing SHACL-required properties:

| Property | Severity | Present |
|----------|----------|---------|
| `ds:description` | Warning | NO (only `rdfs:comment`) |
| `ds:website` | Warning | NO |
| `ds:hasLicense` | Violation | NO |
| `ds:hasVersion` | Warning | NO |
| `ds:hasAccessMethod` | Warning | NO |
| `ds:maintainedBy` | Warning | NO |

**Recommendation:** Phase 2 should enrich main file with required properties from category-specific files.

### 6. Dual Instance Representation (ARCHITECTURAL)

**Problem:** Two parallel representations exist:
1. `datasources.ttl` - Minimal, follows OWL structure
2. `01-genetics-genomics.ttl` etc. - Detailed, uses non-standard properties

**Recommendation:** Consolidate into single canonical representation or formally document relationship.

### 7. ID Pattern Inconsistencies (LOW PRIORITY)

**Problem:** Inconsistent use of hyphens vs dots in IDs:
- Hyphens: `1000-genomes`, `brca-exchange`
- Dots: `dr.dukes`, `hit.2.0`

**Recommendation:** Establish ID naming convention for new sources.

---

## Changes Made (Phase 1)

### 1. datasources.ttl Updates

**Categories (9 total):**
```turtle
# Before
cat:GeneticsGenomics a ds:Category ;
    ds:title "GeneticsGenomics" .

# After
cat:GeneticsGenomics a ds:Category ;
    ds:title "Genetics and Genomics" ;
    ds:notation "01" .
```

**Subcategories (38 total):**
```turtle
# Before
subcat:VariantRepositories a ds:Subcategory ;
    ds:title "VariantRepositories" .

# After
subcat:VariantRepositories a ds:Subcategory ;
    ds:title "Variant Repositories" ;
    ds:notation "1.1" .
```

### 2. OWL Ontology Updates

**New Properties Added:**
```turtle
ds:approximateRecordCount a owl:DatatypeProperty, owl:FunctionalProperty ;
    rdfs:label "approximate record count" ;
    rdfs:domain ds:DataSource ;
    rdfs:range xsd:string .

ds:updateFrequencyLabel a owl:DatatypeProperty, owl:FunctionalProperty ;
    rdfs:label "update frequency label" ;
    rdfs:domain ds:DataSource ;
    rdfs:range xsd:string .
```

---

## Validation Status

After Phase 1 fixes:

| Component | Status | Notes |
|-----------|--------|-------|
| Category titles | PASS | Match SKOS prefLabels |
| Category notation | PASS | Format "XX" |
| Subcategory titles | PASS | Proper spacing |
| Subcategory notation | PASS | Format "X.Y" |
| OWL property coverage | PARTIAL | Core properties defined |
| SHACL validation | PARTIAL | Some warnings remain |

---

## Phase 2 Recommendations

### Priority 1: Harmonize Category-Specific Files
1. Update `ds:maintainer` to `ds:maintainedBy` with object reference
2. Update `ds:category` to `ds:belongsToCategory`
3. Update `ds:tier` and `ds:status` to use SKOS concept references
4. Add inverse properties or document pattern

### Priority 2: Enrich Main Instance File
1. Add `ds:description` to DataSources (copy from category files)
2. Add `ds:website` links
3. Add `ds:hasLicense` references
4. Add `ds:hasVersion` references
5. Add `ds:hasAccessMethod` references

### Priority 3: Consolidation Strategy
1. Define canonical representation
2. Create validation pipeline
3. Add automated consistency checks

---

## Files Modified

| File | Changes |
|------|---------|
| `ontology/instances/datasources.ttl` | Fixed 9 categories, 38 subcategories |
| `ontology/owl/datasource-ontology.ttl` | Added 2 properties |

---

## Appendix: Category Notation Mapping

| Notation | Category ID | Title |
|----------|-------------|-------|
| 01 | genetics.genomics | Genetics and Genomics |
| 02 | compounds.molecules | Compounds and Molecules |
| 03 | diseases.phenotypes | Diseases and Phenotypes |
| 04 | pathways.networks | Pathways and Networks |
| 05 | traditional.medicine | Traditional Medicine |
| 06 | nutrition.food | Nutrition and Food |
| 07 | proteins.molecular.biology | Proteins and Molecular Biology |
| 08 | literature.knowledge | Literature and Knowledge |
| 09 | microbiome | Microbiome |

---

## Appendix: Subcategory Notation Mapping

| Category | Notation | Subcategory ID | Title |
|----------|----------|----------------|-------|
| 01 | 1.1 | variant.repositories | Variant Repositories |
| 01 | 1.2 | functional.prediction | Functional Prediction |
| 01 | 1.3 | population.genetics | Population Genetics |
| 01 | 1.4 | pharmacogenomics | Pharmacogenomics |
| 01 | 1.5 | expression.regulation | Expression and Regulation |
| 01 | 1.6 | cancer.genomics | Cancer Genomics |
| 02 | 2.1 | natural.products | Natural Products |
| 02 | 2.2 | pharmaceuticals | Pharmaceuticals |
| 02 | 2.3 | traditional.medicine.compounds | Traditional Medicine Compounds |
| 02 | 2.4 | food.compounds.nutrients | Food Compounds and Nutrients |
| 02 | 2.5 | drug.metabolism.pharmacokinetics | Drug Metabolism and Pharmacokinetics |
| 02 | 2.6 | chemical.ontology.classification | Chemical Ontology and Classification |
| 02 | 2.7 | compound.target.interactions | Compound-Target Interactions |
| 03 | 3.1 | disease.ontologies | Disease Ontologies |
| 03 | 3.2 | phenotype.databases | Phenotype Databases |
| 03 | 3.3 | disease.gene.associations | Disease-Gene Associations |
| 03 | 3.4 | cancer.oncology | Cancer and Oncology |
| 03 | 3.5 | rare.orphan.diseases | Rare and Orphan Diseases |
| 03 | 3.6 | autoimmune.inflammatory | Autoimmune and Inflammatory |
| 03 | 3.7 | mental.health.neurological | Mental Health and Neurological |
| 04 | 4.1 | metabolic.pathways | Metabolic Pathways |
| 04 | 4.2 | signaling.pathways | Signaling Pathways |
| 04 | 4.3 | protein.protein.interactions | Protein-Protein Interactions |
| 04 | 4.4 | drug.target.interactions | Drug-Target Interactions |
| 04 | 4.5 | gene.function.ontology | Gene Function Ontology |
| 04 | 4.6 | regulatory.networks | Regulatory Networks |
| 05 | 5.1 | traditional.chinese.medicine | Traditional Chinese Medicine |
| 05 | 5.2 | south.east.asian.systems | South and East Asian Systems |
| 05 | 5.3 | western.global.herbal | Western and Global Herbal |
| 05 | 5.4 | multi.system.integration | Multi-System Integration |
| 06 | 6.1 | food.composition | Food Composition |
| 06 | 6.2 | dietary.supplements | Dietary Supplements |
| 06 | 6.3 | bioactive.food.compounds | Bioactive Food Compounds |
| 06 | 6.4 | metabolomics | Metabolomics |
| 07 | 7.1 | protein.sequences.annotations | Protein Sequences and Annotations |
| 07 | 7.2 | protein.structures | Protein Structures |
| 07 | 7.3 | molecular.interactions | Molecular Interactions |
| 08 | 8.1 | scientific.literature | Scientific Literature |
| 08 | 8.2 | knowledge.bases | Knowledge Bases |
| 08 | 8.3 | identifier.mapping | Identifier Mapping |
| 08 | 8.4 | regulatory.legal | Regulatory and Legal |
| 09 | 9.1 | gut.microbiome | Gut Microbiome |
| 09 | 9.2 | body.site.microbiomes | Body Site Microbiomes |
| 09 | 9.3 | microbe.host.interactions | Microbe-Host Interactions |
