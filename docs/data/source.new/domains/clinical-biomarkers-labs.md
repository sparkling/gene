---
id: clinical-biomarkers-labs
title: Biomarker and Laboratory Reference Data Sources
type: clinical
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [clinical, biomarkers, labs, reference-ranges, loinc, markerdb, caliper]
---

# Biomarker and Laboratory Reference Data Sources

**Document ID:** 43-81-BIOMARKERS-LABS
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [Health Domains](./_index.md)

---

## TL;DR

Biomarker and laboratory reference databases provide reference ranges, standardized test codes, and clinical chemistry data for comprehensive lab result interpretation. MarkerDB 2.0 and HMDB 5.0 serve as primary sources for biomarker reference intervals, LOINC provides universal test identification, and CALIPER delivers pediatric-specific ranges. Integration of 21+ databases enables complete lab interpretation from test identification through age/sex-stratified reference ranges to functional/optimal interpretation.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary biomarker database | MarkerDB 2.0 | REST API, age/sex/BMI stratification, free | Jan 2026 |
| Metabolite concentrations | HMDB 5.0 | 220K metabolites, comprehensive coverage | Jan 2026 |
| Test code standardization | LOINC | Universal standard, free with attribution | Jan 2026 |
| Pediatric reference | CALIPER | Gold standard for birth-18 years | Jan 2026 |
| Cancer biomarkers | OncoKB | FDA-recognized, MSK maintained | Jan 2026 |
| Clinical terminology | SNOMED CT + UMLS | Cross-vocabulary integration | Jan 2026 |

---

## Database Catalog

### Biomarker Databases

#### 1. MarkerDB 2.0

| Field | Value |
|-------|-------|
| **URL** | https://markerdb.ca/ |
| **Maintainer** | University of Alberta (HMDB team) |
| **Content** | 218 protein biomarkers, 1,664 chemical biomarkers, 154 karyotype biomarkers, 32,447 genetic markers |
| **Records** | 34,000+ total biomarkers |
| **License** | Free for academic use |
| **API** | REST API for programmatic access |
| **Bulk Download** | Yes (downloadable tables) |
| **Data Formats** | JSON, downloadable tables |
| **Update Frequency** | Periodic (Version 2.0, 2025) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |

**Key Features:**
- Reference concentration intervals for plasma/serum/blood (2.5th-97.5th percentiles)
- CALIPER pediatric intervals integrated
- HMDB 5.0 adult reference values
- Age/sex/BMI/ethnicity-stratified ranges

**Biomarker Types:**
- Diagnostic biomarkers
- Predictive biomarkers
- Prognostic biomarkers
- Exposure biomarkers
- Condition-specific biomarkers

**Publication:** [NAR 2025](https://academic.oup.com/nar/article/53/D1/D1415/7899528)

---

#### 2. HMDB 5.0 (Human Metabolome Database)

| Field | Value |
|-------|-------|
| **URL** | https://hmdb.ca/ |
| **Maintainer** | University of Alberta |
| **Content** | 220,945 metabolite entries with concentration data. 19,715+ compound concentrations in serum/plasma/urine/CSF. |
| **Records** | 220,945 metabolites |
| **License** | Free for academic; commercial requires permission |
| **API** | REST API (contact required: eponine@ualberta.ca) |
| **Bulk Download** | https://hmdb.ca/downloads |
| **Data Formats** | XML, SDF, CSV |
| **Update Frequency** | Periodic (Version 5.0, 2022) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~10 GB (full), ~2 GB (human metabolites) |

**Key Features:**
- Normal/abnormal ranges
- Vitamin levels
- SRM-1950 reference values
- Spectral data (NMR, MS/MS, GC-MS)

**Publication:** [NAR 2022](https://academic.oup.com/nar/article/50/D1/D622/6431815)

---

#### 3. TheMarker (Therapeutic Biomarkers)

| Field | Value |
|-------|-------|
| **URL** | https://idrblab.org/themarker/ |
| **Maintainer** | Zhejiang University |
| **Content** | FDA-defined therapeutic biomarkers (5 types) |
| **Records** | Comprehensive therapeutic biomarker coverage |
| **License** | Free, no login required |
| **API** | Web interface only |
| **Bulk Download** | Limited |
| **Data Formats** | Web interface |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB |

**Biomarker Types:**
- Pharmacodynamic (PDY)
- Safety (SAF)
- Monitoring (MOI)
- Predictive (PRD)
- Surrogate Endpoint (SUR)

**Publication:** [NAR 2024](https://academic.oup.com/nar/article/52/D1/D1450/7321069)

---

#### 4. OncoKB (Cancer Biomarkers)

| Field | Value |
|-------|-------|
| **URL** | https://www.oncokb.org/ |
| **Maintainer** | Memorial Sloan Kettering Cancer Center |
| **Content** | Cancer genomic alterations, biomarker-drug associations |
| **Records** | Comprehensive cancer biomarker coverage |
| **License** | Free for academic; commercial license required |
| **API** | REST API (token required) - https://api.oncokb.org/oncokb-website/api |
| **Bulk Download** | Cancer Gene List, Biomarker-Drug Association List (free) |
| **Data Formats** | JSON |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 (MVP) for oncology |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- First FDA-recognized somatic tumor mutation database
- R Package: oncoKBData (Bioconductor)
- GitHub: https://github.com/oncokb/oncokb

---

#### 5. BIONDA (Text Mining-Based Biomarkers)

| Field | Value |
|-------|-------|
| **URL** | http://bionda.mpc.ruhr-uni-bochum.de/ |
| **Maintainer** | Ruhr University Bochum |
| **Content** | Molecular biomarkers (genes, proteins, miRNAs) for diseases |
| **Records** | Text-mined from Europe PMC |
| **License** | Free |
| **API** | Europe PMC API integration |
| **Bulk Download** | Via GitHub |
| **Data Formats** | Various |
| **Update Frequency** | Periodic |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~100 MB |

**GitHub:** https://github.com/mpc-bioinformatics/bionda

**Publication:** [PMC 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9710600/)

---

#### 6. GOBIOM Database

| Field | Value |
|-------|-------|
| **URL** | https://gobiomdbplus.com/ |
| **Maintainer** | Commercial provider |
| **Content** | Clinically evaluated and preclinical biomarkers from global clinical trials |
| **Records** | Comprehensive (340 data points per marker) |
| **License** | Commercial subscription required |
| **API** | Subscription-based |
| **Bulk Download** | Subscription-based |
| **Data Formats** | Multiple |
| **Update Frequency** | Regular |
| **Priority** | Tier 3 |
| **Storage Estimate** | N/A (commercial) |

**Biomarker Types:**
- Biochemical
- Genomic
- Proteomic
- Imaging
- Metabolite markers

---

### Laboratory Test Standards and Terminologies

#### 7. LOINC (Logical Observation Identifiers Names and Codes)

| Field | Value |
|-------|-------|
| **URL** | https://loinc.org/ |
| **Maintainer** | Regenstrief Institute |
| **Content** | Universal codes for lab tests, clinical observations |
| **Records** | 100,000+ codes |
| **License** | Free with attribution |
| **API** | FHIR-based access |
| **Bulk Download** | https://loinc.org/downloads/ (free registration) |
| **Data Formats** | CSV, FHIR |
| **Update Frequency** | Regular (Version 2.70, August 2025) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~78 MB |

**Key Features:**
- Standard identifiers for lab tests
- RELMA mapping assistant
- SearchLOINC web search

**Important Note:** LOINC provides universal CODES for lab tests but does NOT provide standardized reference ranges. Reference intervals must be determined locally by each laboratory.

---

#### 8. SNOMED CT (Laboratory Content)

| Field | Value |
|-------|-------|
| **URL** | https://www.nlm.nih.gov/healthit/snomedct/index.html |
| **Maintainer** | SNOMED International / NLM |
| **Content** | Clinical terminology including laboratory test concepts |
| **Records** | ~2 GB (full SNOMED CT) |
| **License** | Free UMLS license required |
| **API** | SNOMED CT Browser, terminology servers |
| **Bulk Download** | Via UMLS Terminology Services |
| **Data Formats** | RF2 format |
| **Update Frequency** | Biannual |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 GB |

**Key Features:**
- Pathology Bounded Code List (PBCL)
- PaLM observable entity reference sets
- Cross-maps to LOINC, ICD-10, ICD-O-3

---

#### 9. UMLS Metathesaurus

| Field | Value |
|-------|-------|
| **URL** | https://www.nlm.nih.gov/research/umls/index.html |
| **Maintainer** | NIH/NLM |
| **Content** | Integrates vocabularies including LOINC, SNOMED CT, CPT, ICD-10-CM |
| **Records** | ~30 GB (full Metathesaurus) |
| **License** | Free UMLS license (individual registration) |
| **API** | REST API (https://documentation.uts.nlm.nih.gov/) |
| **Bulk Download** | https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html |
| **Data Formats** | RRF format |
| **Update Frequency** | Regular |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~30 GB |

**Browser:** https://uts.nlm.nih.gov/

---

#### 10. MedlinePlus Lab Tests

| Field | Value |
|-------|-------|
| **URL** | https://medlineplus.gov/ |
| **Maintainer** | NIH/NLM |
| **Content** | 125+ lab test information articles (English/Spanish) |
| **Records** | 125+ test explanations |
| **License** | Public domain (NLM) |
| **API** | MedlinePlus Connect (https://connect.medlineplus.gov/service) |
| **Bulk Download** | No (API access) |
| **Data Formats** | XML |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- LOINC-mapped
- Consumer-friendly explanations
- Rate Limit: 85 requests/minute per IP

**Documentation:** https://medlineplus.gov/medlineplus-connect/technical-information/

---

### Clinical Chemistry Reference Databases

#### 11. CALIPER (Pediatric Reference Intervals)

| Field | Value |
|-------|-------|
| **URL** | https://caliperproject.ca/ |
| **Maintainer** | SickKids Hospital, Toronto |
| **Content** | 200+ age- and sex-specific reference intervals (birth to 18 years) |
| **Records** | 200+ test intervals |
| **License** | Free for healthcare professionals |
| **API** | None (web interface and mobile apps) |
| **Bulk Download** | No (web only) |
| **Data Formats** | Web interface |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 MB |

**Database:** https://caliperproject.ca/caliper/database/

**Test Types:**
- Biochemical
- Immunological
- Hematological
- Nutritional
- Endocrine
- Fertility

**Platforms Supported:**
- Abbott ARCHITECT
- Beckman Coulter AU and DxC
- Ortho VITROS
- Roche Cobas
- Siemens ADVIA

**Mobile Apps:** iOS and Android (CALIPER Reference App)

**Publication:** [Critical Reviews in Clinical Laboratory Sciences 2017](https://www.tandfonline.com/doi/full/10.1080/10408363.2017.1379945)

---

#### 12. ABIM Laboratory Reference Ranges

| Field | Value |
|-------|-------|
| **URL** | https://www.abim.org/Media/bfijryql/laboratory-reference-ranges.pdf |
| **Maintainer** | American Board of Internal Medicine |
| **Content** | Standard reference ranges for board examinations |
| **Records** | Comprehensive clinical chemistry reference |
| **License** | Public |
| **API** | None |
| **Bulk Download** | PDF |
| **Data Formats** | PDF |
| **Update Frequency** | Annual (January 2025) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~1 MB |

**Key Features:**
- Widely used as clinical baseline reference
- Board examination standard

---

#### 13. Metabolomics Workbench

| Field | Value |
|-------|-------|
| **URL** | https://www.metabolomicsworkbench.org/ |
| **Maintainer** | NIH Common Fund |
| **Content** | 4,215+ metabolomics studies, 175,466 metabolite structures |
| **Records** | 175,466 structures |
| **License** | Free (NIH Common Fund) |
| **API** | REST service (https://www.metabolomicsworkbench.org/tools/mw_rest.php) |
| **Bulk Download** | Yes |
| **Data Formats** | JSON, various |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 GB |

**Key Features:**
- RefMet standardized nomenclature
- Experimental data and protocols

**Database:** https://www.metabolomicsworkbench.org/databases/metabolitedatabase.php

---

#### 14. BioGPS (Gene Expression Biomarkers)

| Field | Value |
|-------|-------|
| **URL** | http://biogps.org/ |
| **Maintainer** | Scripps Research |
| **Content** | Gene expression profiles across tissues/cell types |
| **Records** | ~6,000 datasets from ArrayExpress + BioGPS-specific |
| **License** | Free |
| **API** | REST API (http://biogps.org/api/) |
| **Bulk Download** | Yes |
| **Data Formats** | JSON |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~1 GB |

**Related Service:** MyGene.info (http://mygene.info/)

**Publication:** [NAR 2016](https://academic.oup.com/nar/article/44/D1/D313/2502613)

---

### Functional Medicine Reference Ranges

#### 15. Optimal DX

| Field | Value |
|-------|-------|
| **URL** | https://www.optimaldx.com/ |
| **Maintainer** | Optimal DX |
| **Content** | Functional/optimal ranges for common biomarkers |
| **Records** | 500+ optimal ranges |
| **License** | Free educational content |
| **API** | None |
| **Bulk Download** | No |
| **Data Formats** | Web content |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~10 MB |

**Key Features:**
- Narrower "optimal" ranges vs. standard lab ranges
- Functional Blood Chemistry Analysis methodology

**Example Optimal vs Standard Ranges:**

| Biomarker | Standard Range | Functional/Optimal Range |
|-----------|---------------|-------------------------|
| TSH | 0.5-4.5 mIU/L | 1.0-2.5 mIU/L |
| Vitamin D | 30-100 ng/mL | 50-80 ng/mL |
| HbA1c | <5.7% | 5.0-5.4% |
| CRP | <3.0 mg/L | <0.3 mg/L |
| Fasting Glucose | 70-100 mg/dL | 82-88 mg/dL |
| Ferritin (men) | 12-300 ng/mL | 50-150 ng/mL |

---

#### 16. Functional Reference Limits (Research)

| Field | Value |
|-------|-------|
| **URL** | https://pmc.ncbi.nlm.nih.gov/articles/PMC10151278/ |
| **Maintainer** | Academic publication |
| **Content** | Academic framework for physiological relationships and enhanced interpretation |
| **Records** | Research framework |
| **License** | Open access publication |
| **API** | None |
| **Bulk Download** | No |
| **Data Formats** | PDF |
| **Update Frequency** | N/A |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~5 MB |

---

### Age/Sex-Specific Reference Databases

#### 17. Chinese Population Reference Intervals Study

| Field | Value |
|-------|-------|
| **URL** | https://pmc.ncbi.nlm.nih.gov/articles/PMC7956001/ |
| **Maintainer** | Academic publication |
| **Content** | Gender and age-specific reference intervals for common biochemical analytes |
| **Records** | Population-specific intervals |
| **License** | Open access |
| **API** | None |
| **Bulk Download** | No |
| **Data Formats** | PDF/tables |
| **Update Frequency** | N/A |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~5 MB |

**Key Findings:**
- Gender-specific RIs needed for: ALT, AST, GGT, DBil, TBil, UA, Cr
- Age/gender-specific RIs for urea and ALP
- No stratification needed for: K, Na, Cl, Ca, IP

---

#### 18. Transgender Reference Intervals

| Field | Value |
|-------|-------|
| **URL** | https://academic.oup.com/jalm/article/7/5/1131/6588198 |
| **Maintainer** | Journal of Applied Laboratory Medicine |
| **Content** | Reference intervals for transgender men and women on stable hormone therapy |
| **Records** | Clinical chemistry analytes |
| **License** | Academic publication |
| **API** | None |
| **Bulk Download** | No |
| **Data Formats** | PDF |
| **Update Frequency** | N/A |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~5 MB |

---

#### 19. CLSI EP28 Guidelines

| Field | Value |
|-------|-------|
| **URL** | https://clsi.org/shop/standards/ep28/ |
| **Maintainer** | Clinical and Laboratory Standards Institute |
| **Content** | Guidelines for defining, establishing, verifying reference intervals |
| **Records** | Framework document |
| **License** | Paid standard document |
| **API** | None |
| **Bulk Download** | Purchase required |
| **Data Formats** | PDF |
| **Update Frequency** | Periodic |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~10 MB |

**Key Features:**
- Gold standard methodology for RI establishment
- Minimum sample: 120 apparently healthy individuals
- Method: 2.5th and 97.5th percentiles (parametric or non-parametric)

---

### Commercial Laboratory APIs

#### 20. Quest Diagnostics

| Field | Value |
|-------|-------|
| **URL** | https://testdirectory.questdiagnostics.com/ |
| **Maintainer** | Quest Diagnostics |
| **Content** | Lab test directory with reference ranges |
| **Records** | Comprehensive test catalog |
| **License** | Business agreement required |
| **API** | Via Health Gorilla, HL7/FHIR integrations |
| **Bulk Download** | No |
| **Data Formats** | FHIR Observation resources with LOINC codes |
| **Update Frequency** | Continuous |
| **Priority** | Tier 3 |
| **Storage Estimate** | N/A (integration) |

---

#### 21. Health Gorilla (Aggregator)

| Field | Value |
|-------|-------|
| **URL** | https://developer.healthgorilla.com/docs/diagnostic-network |
| **Maintainer** | Health Gorilla |
| **Content** | Unified access to LabCorp, Quest, BioReference |
| **Records** | Aggregated lab network |
| **License** | Commercial |
| **API** | REST API |
| **Bulk Download** | No |
| **Data Formats** | FHIR, JSON |
| **Update Frequency** | Real-time |
| **Priority** | Tier 3 |
| **Storage Estimate** | N/A (integration) |

**Note:** Reference ranges from commercial labs are institution-specific and included with test results, not available as standalone databases.

---

## Comparison Matrix

### Data Coverage by Category

| Database | Biomarkers | Reference Ranges | Test Codes | Pediatric | Metabolites |
|----------|------------|------------------|------------|-----------|-------------|
| MarkerDB 2.0 | PRIMARY | YES | Cross-ref | CALIPER | YES |
| HMDB 5.0 | Metabolites | YES | No | Partial | PRIMARY |
| OncoKB | Cancer | No | No | No | No |
| LOINC | No | No | PRIMARY | No | No |
| CALIPER | No | PRIMARY (peds) | Cross-ref | PRIMARY | No |
| MedlinePlus | No | Partial | LOINC | No | No |
| BioGPS | Gene expr | No | No | No | No |
| Metabolomics WB | No | Research | No | No | YES |

### Access and Licensing

| Database | API | Bulk Download | Academic | Commercial | Priority |
|----------|-----|---------------|----------|------------|----------|
| MarkerDB 2.0 | REST | Yes | Free | Contact | Tier 1 |
| HMDB 5.0 | Contact | Yes | Free | Permission | Tier 1 |
| OncoKB | REST | Partial | Free | License | Tier 1 |
| LOINC | FHIR | Yes | Free | Free | Tier 1 |
| SNOMED CT | Various | Yes | Free (UMLS) | National | Tier 2 |
| MedlinePlus | REST | No | Free | Free | Tier 1 |
| CALIPER | No | No | Free | Unknown | Tier 1 |
| BioGPS | REST | Yes | Free | Free | Tier 2 |

---

## Integration Priority

### Tier 1: Core Implementation (Weeks 1-4)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| LOINC | Test Codes | Universal test code standardization | Free |
| MarkerDB 2.0 | Biomarkers | Reference intervals with API | Free |
| CALIPER | Pediatric | Birth-18 years reference intervals | Free |
| MedlinePlus Connect | Consumer | Consumer-friendly test explanations | Public Domain |
| OncoKB | Oncology | Cancer biomarkers, FDA-recognized | Free (academic) |
| HMDB 5.0 | Metabolites | Metabolite concentrations | Free |

### Tier 2: Enrichment (Weeks 5-8)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| SNOMED CT | Terminology | Laboratory terminology | UMLS Free |
| BioGPS | Gene Expression | Tissue expression biomarkers | Free |
| TheMarker | Therapeutic | FDA therapeutic biomarker types | Free |
| Metabolomics Workbench | Research | Metabolomics study data | Free |
| ABIM Reference | Clinical | Standard clinical reference | Public |

### Tier 3: Specialized (Weeks 9-12)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| Functional ranges | Optimal | User preference option | Various |
| Age/sex studies | Population | Population-specific ranges | Open |
| GOBIOM | Commercial | Comprehensive clinical trials | Subscription |
| Quest/LabCorp | Integration | Lab result integration | Commercial |
| CLSI Guidelines | Framework | RI establishment methodology | Paid |

---

## Data Pipeline Architecture

```
Lab Result Input
         |
         v
+------------------+
| Normalize Codes  | <-- LOINC mapping
+------------------+
         |
         v
+------------------+
| Reference Lookup | --> MarkerDB 2.0, CALIPER (by age/sex)
+------------------+
         |
         v
+------------------+
| Interpretation   | --> MedlinePlus (consumer explanation)
+------------------+
         |
         v
+------------------+
| Biomarker Context| --> OncoKB (cancer), HMDB (metabolites)
+------------------+
         |
         v
+------------------+
| Optional Optimal | --> Functional ranges (user preference)
+------------------+
         |
         v
    Interpreted Result
```

---

## Data Quality Considerations

### Reference Interval Variability Factors

1. **Analytical method** - Different assays produce different ranges
2. **Population** - Age, sex, race, ethnicity, pregnancy
3. **Pre-analytical** - Fasting state, time of day, posture
4. **Geographic** - Regional population differences
5. **Instrument platform** - Manufacturer-specific calibration

### Recommendations for Gene Platform

1. **Use LOINC** for standardized test identification
2. **Implement MarkerDB 2.0** for biomarker reference intervals
3. **Integrate CALIPER** for pediatric populations
4. **Consider HMDB 5.0** for metabolite concentrations
5. **Display ranges with appropriate disclaimers** about laboratory variation
6. **Support multiple reference range sources** for user customization

---

## Storage Estimates Summary

| Tier | Databases | Estimated Storage |
|------|-----------|-------------------|
| Tier 1 | LOINC, MarkerDB, CALIPER, MedlinePlus, OncoKB, HMDB | ~13 GB |
| Tier 2 | SNOMED CT, BioGPS, TheMarker, Metabolomics WB, ABIM | ~8 GB |
| Tier 3 | Functional, Population studies, GOBIOM, CLSI | ~1 GB |
| **Total** | All | **~22 GB** |

*Note: UMLS Metathesaurus adds ~30 GB if full vocabulary needed.*

---

## API Summary Table

| Database | API Available | Authentication | Rate Limits | Format |
|----------|---------------|----------------|-------------|--------|
| MarkerDB 2.0 | Yes | None | Unknown | JSON |
| HMDB 5.0 | Yes | Contact required | Unknown | XML, JSON |
| OncoKB | Yes | Token required | Unknown | JSON |
| LOINC | FHIR | Registration | Unknown | FHIR |
| MedlinePlus Connect | Yes | None | 85/min/IP | XML |
| UMLS | REST | License required | Unknown | JSON |
| BioGPS | Yes | None | Unknown | JSON |
| Metabolomics Workbench | REST | None | Unknown | JSON |
| CALIPER | No | N/A | N/A | Web only |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [index.md](./../index.md) | Parent index |
| [drug-metabolism.md](./../compounds/drug-metabolism.md) | Drug metabolism integration |
| [pharmaceuticals.md](./../compounds/pharmaceuticals.md) | Pharmacogenomics data |
| [primary.md](./../genetics/primary.md) | Genetic biomarkers |

---

## Open Questions

- [ ] CALIPER data extraction - manual curation or partnership?
- [ ] Functional ranges - how to present vs. standard ranges (user setting)?
- [ ] Commercial lab integration - Health Gorilla partnership?
- [ ] HMDB API access - formal agreement required?
- [ ] Pediatric vs. adult cutoff - 18 or 21 years?

---

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `biomarker` | Measurable indicator of a biological state or condition | HbA1c for diabetes |
| `reference_interval` | Range of values that includes 95% of a healthy population (2.5th-97.5th percentile) | Glucose: 70-100 mg/dL |
| `reference_range` | Synonym for reference interval; lab-specific normal values | TSH: 0.5-4.5 mIU/L |
| `LOINC_code` | Universal identifier for laboratory tests and clinical observations | 2345-7 (Glucose) |
| `analyte` | Substance being measured in a laboratory test | Creatinine, ALT |
| `percentile` | Statistical measure indicating percentage of population below a value | 97.5th percentile |
| `sensitivity` | Ability of a test to correctly identify positive cases | 95% sensitivity |
| `specificity` | Ability of a test to correctly identify negative cases | 90% specificity |
| `metabolite` | Small molecule intermediate or product of metabolism | Lactate, pyruvate |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| CALIPER | Canadian Laboratory Initiative on Pediatric Reference Intervals | Pediatric ranges |
| HMDB | Human Metabolome Database containing metabolite concentrations | Metabolomics |
| MarkerDB | Database of biomarkers with reference intervals | Biomarker lookup |
| OncoKB | Memorial Sloan Kettering cancer biomarker database | Oncology |
| functional_range | Narrower "optimal" range beyond standard lab reference | Optimal health |
| SRM-1950 | Standard Reference Material for human serum metabolite concentrations | NIST standard |
| RefMet | Standardized nomenclature for metabolites | Naming convention |
| FHIR | Fast Healthcare Interoperability Resources standard | Data exchange |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| LOINC | Logical Observation Identifiers Names and Codes | Universal test codes |
| SNOMED CT | Systematized Nomenclature of Medicine Clinical Terms | Clinical terminology |
| UMLS | Unified Medical Language System | Vocabulary integration |
| CALIPER | Canadian Laboratory Initiative on Pediatric Reference Intervals | Pediatric reference |
| HMDB | Human Metabolome Database | Metabolite data |
| CRP | C-Reactive Protein | Inflammation marker |
| HbA1c | Glycated Hemoglobin A1c | Diabetes marker |
| TSH | Thyroid Stimulating Hormone | Thyroid function |
| ALT | Alanine Aminotransferase | Liver enzyme |
| AST | Aspartate Aminotransferase | Liver enzyme |
| GGT | Gamma-Glutamyl Transferase | Liver enzyme |
| ALP | Alkaline Phosphatase | Bone/liver enzyme |
| PDY | Pharmacodynamic biomarker | Drug response |
| SAF | Safety biomarker | Adverse effects |
| MOI | Monitoring biomarker | Disease tracking |
| PRD | Predictive biomarker | Treatment response |
| SUR | Surrogate endpoint biomarker | Clinical trials |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Initial catalog from data-sources-biomarkers-labs.md |
