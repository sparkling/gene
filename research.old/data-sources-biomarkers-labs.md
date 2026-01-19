# Biomarker and Lab Reference Databases Inventory

**Last Updated:** January 2026
**Purpose:** Reference ranges, biomarker databases, clinical chemistry references for Gene Platform
**Constraint:** Focus on publicly available or academically accessible sources

---

## Summary

| Category | Sources | Total Records | Est. Size | Primary Use |
|----------|---------|---------------|-----------|-------------|
| Biomarker Databases | 6 | 250K+ markers | ~5 GB | Biomarker lookup & reference |
| Lab Reference Standards | 4 | 100K+ codes | ~500 MB | Lab test standardization |
| Clinical Chemistry Refs | 5 | 50K+ analytes | ~2 GB | Reference intervals |
| Functional Medicine Refs | 3 | 500+ optimal ranges | ~10 MB | Optimal range interpretation |
| Age/Sex-Specific Refs | 3 | 10K+ intervals | ~200 MB | Population-specific ranges |

---

## 1. Biomarker Databases

### 1.1 MarkerDB 2.0 - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://markerdb.ca/ |
| **Content** | 218 protein biomarkers, 1,664 chemical biomarkers, 154 karyotype biomarkers, 32,447 genetic markers |
| **Version** | 2.0 (2025) |
| **API** | REST API for programmatic access |
| **Formats** | JSON, downloadable tables |
| **License** | Free for academic use |
| **Size** | ~500 MB |
| **Key Features** | Reference concentration intervals for plasma/serum/blood (2.5th-97.5th percentiles), CALIPER pediatric intervals integrated, HMDB 5.0 adult reference values, age/sex/BMI/ethnicity-stratified ranges |
| **Publication** | [NAR 2025](https://academic.oup.com/nar/article/53/D1/D1415/7899528) |

**Biomarker Types:**
- Diagnostic biomarkers
- Predictive biomarkers
- Prognostic biomarkers
- Exposure biomarkers
- Condition-specific biomarkers

### 1.2 Human Metabolome Database (HMDB) 5.0 - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://hmdb.ca/ |
| **Content** | 220,945 metabolite entries with concentration data |
| **Version** | 5.0 (2022) |
| **API** | REST API (contact required: eponine@ualberta.ca) |
| **Formats** | XML, SDF, CSV |
| **Download** | https://hmdb.ca/downloads |
| **License** | Free for academic, commercial requires permission |
| **Size** | ~10 GB (full), ~2 GB (human metabolites) |
| **Key Features** | 19,715+ compound concentrations in serum/plasma/urine/CSF, normal/abnormal ranges, vitamin levels, SRM-1950 reference values |
| **Publication** | [NAR 2022](https://academic.oup.com/nar/article/50/D1/D622/6431815) |

### 1.3 TheMarker (Therapeutic Biomarkers)
| Field | Value |
|-------|-------|
| **URL** | https://idrblab.org/themarker/ |
| **Content** | FDA-defined therapeutic biomarkers (5 types) |
| **Types** | Pharmacodynamic (PDY), Safety (SAF), Monitoring (MOI), Predictive (PRD), Surrogate Endpoint (SUR) |
| **API** | Web interface, no formal API |
| **License** | Free, no login required |
| **Size** | ~100 MB |
| **Publication** | [NAR 2024](https://academic.oup.com/nar/article/52/D1/D1450/7321069) |

### 1.4 OncoKB (Cancer Biomarkers) - PRIMARY for Oncology
| Field | Value |
|-------|-------|
| **URL** | https://www.oncokb.org/ |
| **Content** | Cancer genomic alterations, biomarker-drug associations |
| **Organization** | Memorial Sloan Kettering Cancer Center |
| **API** | REST API (token required) |
| **API Docs** | https://api.oncokb.org/oncokb-website/api |
| **Download** | Cancer Gene List, Biomarker-Drug Association List (free) |
| **License** | Free for academic research, commercial license required |
| **Size** | ~50 MB |
| **GitHub** | https://github.com/oncokb/oncokb |
| **R Package** | oncoKBData (Bioconductor) |
| **Note** | First FDA-recognized somatic tumor mutation database |

### 1.5 BIONDA (Text Mining-Based Biomarkers)
| Field | Value |
|-------|-------|
| **URL** | http://bionda.mpc.ruhr-uni-bochum.de/ |
| **Content** | Molecular biomarkers (genes, proteins, miRNAs) for diseases |
| **Method** | Text mining from Europe PMC |
| **API** | Europe PMC API integration |
| **License** | Free |
| **GitHub** | https://github.com/mpc-bioinformatics/bionda |
| **Size** | ~100 MB |
| **Publication** | [PMC 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9710600/) |

### 1.6 GOBIOM Database (Commercial)
| Field | Value |
|-------|-------|
| **URL** | https://gobiomdbplus.com/ |
| **Content** | Clinically evaluated and preclinical biomarkers from global clinical trials |
| **Types** | Biochemical, Genomic, Proteomic, Imaging, Metabolite markers |
| **Data Points** | 340 data points per marker |
| **License** | Commercial subscription required |
| **Note** | Comprehensive but not freely accessible |

---

## 2. Laboratory Test Standards & Terminologies

### 2.1 LOINC (Logical Observation Identifiers Names and Codes) - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://loinc.org/ |
| **Content** | Universal codes for lab tests, clinical observations |
| **Version** | 2.70 (August 2025) |
| **Download** | https://loinc.org/downloads/ (free registration) |
| **API** | FHIR-based access |
| **Tools** | RELMA (mapping assistant), SearchLOINC (web search) |
| **License** | Free with attribution |
| **Size** | 78 MB |
| **Key Features** | Standard identifiers for lab tests, not reference ranges (ranges are institution-specific) |

**Important Note:** LOINC provides universal **codes** for lab tests but does NOT provide standardized **reference ranges**. Reference intervals must be determined locally by each laboratory.

### 2.2 SNOMED CT (Laboratory Content)
| Field | Value |
|-------|-------|
| **URL** | https://www.nlm.nih.gov/healthit/snomedct/index.html |
| **Content** | Clinical terminology including laboratory test concepts |
| **Distribution** | Via UMLS Terminology Services (NLM) |
| **API** | SNOMED CT Browser, various terminology servers |
| **License** | Free UMLS license required |
| **Size** | ~2 GB (full SNOMED CT) |
| **Lab Content** | Pathology Bounded Code List (PBCL), PaLM observable entity reference sets |
| **Cross-maps** | LOINC, ICD-10, ICD-O-3 |
| **GitHub** | SNOMED International open source tools |

### 2.3 UMLS Metathesaurus
| Field | Value |
|-------|-------|
| **URL** | https://www.nlm.nih.gov/research/umls/index.html |
| **Content** | Integrates vocabularies including LOINC, SNOMED CT, CPT, ICD-10-CM |
| **API** | REST API (https://documentation.uts.nlm.nih.gov/) |
| **Browser** | https://uts.nlm.nih.gov/ |
| **Download** | https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html |
| **License** | Free UMLS license (individual registration) |
| **Size** | ~30 GB (full Metathesaurus) |

### 2.4 MedlinePlus Lab Tests (Consumer Information)
| Field | Value |
|-------|-------|
| **URL** | https://medlineplus.gov/ |
| **Content** | 125+ lab test information articles (English/Spanish) |
| **API** | MedlinePlus Connect (https://connect.medlineplus.gov/service) |
| **Formats** | XML |
| **License** | Public domain (NLM) |
| **Key Features** | LOINC-mapped, consumer-friendly explanations, what tests mean |
| **Rate Limit** | 85 requests/minute per IP |
| **Documentation** | https://medlineplus.gov/medlineplus-connect/technical-information/ |

---

## 3. Clinical Chemistry Reference Databases

### 3.1 CALIPER (Pediatric Reference Intervals) - PRIMARY for Pediatrics
| Field | Value |
|-------|-------|
| **URL** | https://caliperproject.ca/ |
| **Database** | https://caliperproject.ca/caliper/database/ |
| **Content** | 200+ age- and sex-specific reference intervals (birth to 18 years) |
| **Test Types** | Biochemical, immunological, hematological, nutritional, endocrine, fertility |
| **API** | No formal API; web interface and mobile apps |
| **Mobile Apps** | iOS and Android (CALIPER Reference App) |
| **License** | Free for healthcare professionals |
| **Size** | ~50 MB |
| **Organization** | SickKids Hospital, Toronto, Canada |
| **Key Features** | Multi-platform assay support, continuous updates |
| **Publication** | [Critical Reviews in Clinical Laboratory Sciences 2017](https://www.tandfonline.com/doi/full/10.1080/10408363.2017.1379945) |

**Platforms Supported:**
- Abbott ARCHITECT
- Beckman Coulter AU and DxC
- Ortho VITROS
- Roche Cobas
- Siemens ADVIA

### 3.2 ABIM Laboratory Reference Ranges
| Field | Value |
|-------|-------|
| **URL** | https://www.abim.org/Media/bfijryql/laboratory-reference-ranges.pdf |
| **Content** | Standard reference ranges for board examinations |
| **Version** | January 2025 |
| **Format** | PDF |
| **License** | Public |
| **Size** | ~1 MB |
| **Note** | Widely used as clinical baseline reference |

### 3.3 Metabolomics Workbench
| Field | Value |
|-------|-------|
| **URL** | https://www.metabolomicsworkbench.org/ |
| **Content** | 4,215+ metabolomics studies, 175,466 metabolite structures |
| **API** | REST service (https://www.metabolomicsworkbench.org/tools/mw_rest.php) |
| **Database** | https://www.metabolomicsworkbench.org/databases/metabolitedatabase.php |
| **License** | NIH Common Fund - Free |
| **Size** | ~5 GB |
| **Key Features** | RefMet standardized nomenclature, experimental data, protocols |
| **Funding** | NIH Common Fund Metabolomics Program |

### 3.4 Clinical Laboratory Reference (CLR Online)
| Field | Value |
|-------|-------|
| **URL** | Medical Laboratory Observer catalog |
| **Content** | Critical values, reference interval tables |
| **Format** | PDFs |
| **License** | Subscription-based |
| **Note** | Industry reference for laboratory professionals |

### 3.5 BioGPS (Gene Expression Biomarkers)
| Field | Value |
|-------|-------|
| **URL** | http://biogps.org/ |
| **Content** | Gene expression profiles across tissues/cell types |
| **Datasets** | ~6,000 from ArrayExpress + BioGPS-specific |
| **API** | REST API (http://biogps.org/api/) |
| **Related Service** | MyGene.info (http://mygene.info/) |
| **License** | Free |
| **Size** | ~1 GB (metadata) |
| **Publication** | [NAR 2016](https://academic.oup.com/nar/article/44/D1/D313/2502613) |

---

## 4. Functional Medicine Reference Ranges

### 4.1 Optimal DX
| Field | Value |
|-------|-------|
| **URL** | https://www.optimaldx.com/ |
| **Content** | Functional/optimal ranges for common biomarkers |
| **Approach** | Narrower "optimal" ranges vs. standard lab ranges |
| **Format** | Web content, educational resources |
| **License** | Free educational content |
| **Key Features** | Functional Blood Chemistry Analysis methodology |

**Example Optimal vs Standard Ranges:**

| Biomarker | Standard Range | Functional/Optimal Range |
|-----------|---------------|-------------------------|
| TSH | 0.5-4.5 mIU/L | 1.0-2.5 mIU/L |
| Vitamin D | 30-100 ng/mL | 50-80 ng/mL |
| HbA1c | <5.7% | 5.0-5.4% |
| CRP | <3.0 mg/L | <0.3 mg/L |
| Fasting Glucose | 70-100 mg/dL | 82-88 mg/dL |
| Ferritin (men) | 12-300 ng/mL | 50-150 ng/mL |

### 4.2 Functional Reference Limits (Research)
| Field | Value |
|-------|-------|
| **Publication** | [PMC 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10151278/) |
| **Content** | Academic framework for physiological relationships and enhanced interpretation |
| **Focus** | Describing physiological relationships beyond disease detection |
| **License** | Open access publication |

### 4.3 FullScript Lab Interpretation Guide
| Field | Value |
|-------|-------|
| **URL** | https://fullscript.com/blog/lab-interpretation-in-functional-medicine |
| **Content** | Optimal ranges developed by Medical Advisory Team |
| **Sources** | ScienceDirect, PubMed |
| **Format** | Educational content |
| **License** | Free content |

---

## 5. Age/Sex-Specific Reference Databases

### 5.1 CALIPER (Birth to 18 Years) - See Section 3.1
Primary resource for pediatric age-stratified reference intervals.

### 5.2 Chinese Population Reference Intervals Study
| Field | Value |
|-------|-------|
| **Publication** | [PMC 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7956001/) |
| **Content** | Gender and age-specific reference intervals for common biochemical analytes |
| **Method** | Indirect sampling from real laboratory data |
| **Key Findings** | Gender-specific RIs needed for: ALT, AST, GGT, DBil, TBil, UA, Cr |
| **Note** | Age/gender-specific RIs for urea and ALP |
| **No stratification needed for** | K, Na, Cl, Ca, IP |

### 5.3 Transgender Reference Intervals
| Field | Value |
|-------|-------|
| **Publication** | [JALM 2022](https://academic.oup.com/jalm/article/7/5/1131/6588198) |
| **Content** | Reference intervals for transgender men and women on stable hormone therapy |
| **Tests Covered** | Clinical chemistry analytes |
| **Note** | Important for inclusive reference ranges |

### 5.4 CLSI EP28 Guidelines (Framework)
| Field | Value |
|-------|-------|
| **URL** | https://clsi.org/shop/standards/ep28/ |
| **Content** | Guidelines for defining, establishing, verifying reference intervals |
| **Organization** | Clinical and Laboratory Standards Institute + IFCC |
| **Minimum Sample** | 120 apparently healthy individuals |
| **Method** | 2.5th and 97.5th percentiles (parametric or non-parametric) |
| **License** | Paid standard document |
| **Note** | Gold standard methodology for RI establishment |

---

## 6. Commercial Laboratory APIs (Integration Only)

### 6.1 Quest Diagnostics
| Field | Value |
|-------|-------|
| **Test Directory** | https://testdirectory.questdiagnostics.com/ |
| **API Access** | Via Health Gorilla, HL7/FHIR integrations |
| **Data** | FHIR Observation resources with LOINC codes, reference ranges |
| **License** | Business agreement required |

### 6.2 LabCorp
| Field | Value |
|-------|-------|
| **URL** | https://www.labcorp.com/organizations/data |
| **API Access** | Via Human API, Health Gorilla, direct HL7 |
| **License** | Business agreement required |

### 6.3 Health Gorilla (Aggregator)
| Field | Value |
|-------|-------|
| **URL** | https://developer.healthgorilla.com/docs/diagnostic-network |
| **Content** | Unified access to LabCorp, Quest, BioReference |
| **API** | REST API |
| **License** | Commercial |

**Note:** Reference ranges from commercial labs are institution-specific and included with test results, not available as standalone databases.

---

## 7. Data Quality Considerations

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

## 8. Integration Priority

### Tier 1 - Core Implementation (Weeks 1-4)
1. **LOINC** - Test code standardization
2. **MarkerDB 2.0** - Biomarker reference intervals with API
3. **CALIPER** - Pediatric reference intervals
4. **MedlinePlus Connect** - Consumer-friendly test explanations

### Tier 2 - Enrichment (Weeks 5-8)
1. **HMDB 5.0** - Metabolite concentrations
2. **OncoKB** - Cancer biomarkers
3. **SNOMED CT** - Laboratory terminology
4. **BioGPS** - Gene expression biomarkers

### Tier 3 - Specialized (Weeks 9-12)
1. **Functional medicine optimal ranges** - User preference option
2. **TheMarker** - Therapeutic biomarkers
3. **Metabolomics Workbench** - Research metabolomics data
4. **Age/sex stratification** from population studies

---

## 9. API Summary Table

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

## 10. License Summary

| Database | Academic Use | Commercial Use | Download |
|----------|--------------|----------------|----------|
| MarkerDB 2.0 | Free | Contact required | Yes |
| HMDB 5.0 | Free | Permission required | Yes |
| OncoKB | Free | License required | Partial |
| LOINC | Free | Free with attribution | Yes |
| SNOMED CT | Free (UMLS license) | National license | Yes |
| CALIPER | Free | Unknown | No (web only) |
| MedlinePlus | Public domain | Public domain | N/A |
| BioGPS | Free | Free | Yes |

---

## Sources

### Primary Publications
- [MarkerDB 2.0 - NAR 2025](https://academic.oup.com/nar/article/53/D1/D1415/7899528)
- [HMDB 5.0 - NAR 2022](https://academic.oup.com/nar/article/50/D1/D622/6431815)
- [CALIPER White Paper - CRCLS 2017](https://www.tandfonline.com/doi/full/10.1080/10408363.2017.1379945)
- [OncoKB - JCO Precision Oncology](https://pmc.ncbi.nlm.nih.gov/articles/mid/NIHMS897314/)
- [TheMarker - NAR 2024](https://academic.oup.com/nar/article/52/D1/D1450/7321069)
- [BioGPS - NAR 2016](https://academic.oup.com/nar/article/44/D1/D313/2502613)
- [Metabolomics Workbench - NAR 2016](https://academic.oup.com/nar/article/44/D1/D463/2502588)
- [CLSI EP28 Guidelines](https://clsi.org/shop/standards/ep28/)

### Web Resources
- [LOINC Downloads](https://loinc.org/downloads/)
- [CALIPER Database](https://caliperproject.ca/caliper/database/)
- [MedlinePlus Connect API](https://medlineplus.gov/medlineplus-connect/web-service/)
- [UMLS API Documentation](https://documentation.uts.nlm.nih.gov/)
- [Optimal DX](https://www.optimaldx.com/optimal-range)
- [ABIM Reference Ranges](https://www.abim.org/Media/bfijryql/laboratory-reference-ranges.pdf)
