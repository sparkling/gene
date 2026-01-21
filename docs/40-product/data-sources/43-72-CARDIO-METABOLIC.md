# Cardiovascular and Metabolic Data Sources

**Document ID:** 43-72-CARDIO-METABOLIC
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-00-INDEX.md](./43-00-INDEX.md)

---

## TL;DR

Cardiovascular and metabolic genetics databases provide GWAS summary statistics for heart disease, diabetes, obesity, lipids, and blood pressure phenotypes. Total storage ~88GB across 9 primary consortia with ~500+ summary statistic files covering 1.65M+ individuals. Most require bulk download; only GWAS Catalog, Metabolomics Workbench, and OpenGWAS provide REST APIs.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary CAD source | CARDIoGRAMplusC4D | Largest CAD GWAS meta-analysis, research-standard | Jan 2026 |
| Primary lipid source | GLGC | 923 loci across 1.65M individuals, multi-ancestry | Jan 2026 |
| Primary diabetes source | DIAGRAM/T2DGGI | Latest multi-ancestry T2D GWAS (2.5M participants) | Jan 2026 |
| Primary obesity source | GIANT | Largest BMI/height GWAS consortium | Jan 2026 |
| Metabolite reference | HMDB | 220K metabolites, comprehensive annotations | Jan 2026 |
| API-first integration | GWAS Catalog + OpenGWAS | REST APIs available, comprehensive coverage | Jan 2026 |

---

## Database Catalog

### 1. CARDIoGRAMplusC4D (Cardiovascular Genetics)

| Field | Value |
|-------|-------|
| **URL** | https://cardiogramplusc4d.org/ |
| **Content** | Coronary artery disease and myocardial infarction GWAS |
| **Records** | 9 major datasets; 60K+ cases, 123K+ controls |
| **License** | Research-only; citation required; no re-identification |
| **API** | None (bulk download only) |
| **Update Frequency** | Periodic (study-based releases) |
| **Priority** | Tier 2 (Medium) |
| **Storage Estimate** | ~1.8 GB |

**Available Datasets:**

| Dataset | Cases | Controls | Size |
|---------|-------|----------|------|
| CARDIoGRAM GWAS | 22,233 | 64,762 | 77 MB |
| C4D GWAS | 15,420 | 15,062 | 14 MB |
| CARDIoGRAMplusC4D Metabochip | 63,746 | 130,681 | 4.3 MB |
| 1000 Genomes (Additive) | 60,801 | 123,504 | 260 MB |
| 1000 Genomes (Recessive) | - | - | 190 MB |
| MI Additive | - | - | 389 MB |
| Exome Chip | 42,335 | 78,240 | 2.3 MB |
| Chromosome X-CAD | 43,000+ | 58,000+ | 584 MB |
| UK Biobank Meta-analysis | - | - | 255 MB |

**Data Fields:** SNP ID (rsID), chromosomal position, effect/reference alleles, P-values, odds ratios, sample sizes, heterogeneity (I-squared)

**Access Methods:**
1. Direct file downloads from consortium website
2. HeartBioPortal2.0 web interface (https://heartbioportal.org/)
3. Knowledge Portal Network (https://www.kp4cd.org/node/620)

---

### 2. Global Lipids Genetics Consortium (GLGC)

| Field | Value |
|-------|-------|
| **URL** | https://www.lipidgenetics.org/ |
| **Content** | GWAS for lipid traits (HDL-C, LDL-C, TC, TG, nonHDL-C) |
| **Records** | 923 loci; 1.65M individuals; 500+ summary statistic files |
| **License** | Open access; citation required (Graham et al. 2021 Nature) |
| **API** | None (HTTP downloads from UMich server) |
| **Update Frequency** | Study-based releases |
| **Priority** | Tier 2 (Medium) |
| **Storage Estimate** | ~2 GB |

**Available Dataset Categories:**

| Category | Traits | Files | URL Path |
|----------|--------|-------|----------|
| Ancestry-Specific GWAS | HDL-C, LDL-C, nonHDL-C, TC, TG | 90+ | ancestry_specific |
| Trans-Ancestry GWAS | HDL-C, LDL-C, nonHDL-C, TC | 20+ | trans_ancestry |
| LDL-C Polygenic Scores | All ancestries (6 populations) | 12+ | prs_weights |
| Trans-Ancestry Credible Sets | All 5 traits | 5+ | credible_sets |
| ChrX Summary Statistics | All lipid traits | 180 | chrx_summary_stats |
| Sex-Specific GWAS | All lipid traits | 20 | sex_specific |
| Sex & Ancestry-Specific | All lipid traits | 200 | sex_ancestry |

**Lipid Traits:**
- HDL-C: High-density lipoprotein cholesterol
- LDL-C: Low-density lipoprotein cholesterol
- nonHDL-C: Non-HDL cholesterol
- TC: Total cholesterol
- TG: Triglycerides

**Data Downloads:** https://csg.sph.umich.edu/willer/public/glgc-lipids2021/

---

### 3. DIAGRAM/DIAMANTE/T2DGGI (Diabetes Genetics)

| Field | Value |
|-------|-------|
| **URL** | https://diagram-consortium.org/ |
| **Content** | Type 2 diabetes GWAS summary statistics |
| **Records** | 2.5M+ participants across multi-ancestry studies |
| **License** | Research-only; no re-identification; citation required |
| **API** | None (bulk download); T2D Knowledge Portal for browsing |
| **Update Frequency** | Study-based releases (latest: 2024) |
| **Priority** | Tier 2 (Medium) |
| **Storage Estimate** | ~4 GB |

**Early DIAGRAM Analyses:**

| Dataset | Size |
|---------|------|
| Stage 1 GWAS Summary Statistics | ~32 MB |
| Stage 1 + Stage 2 Metabochip | ~3.3 MB |
| Sex-Specific Summary Statistics | ~7 MB |
| Stage 2 Metabochip Only | ~3.5 MB |

**Expanded Analyses (2015-2018):**

| Dataset | Size |
|---------|------|
| Trans-ethnic T2D GWAS meta-analysis | ~46 MB |
| Gaulton et al. (2015) Metabochip | ~1 MB |
| Gaulton et al. (2015) Credible Sets | ~400 KB |
| Fuchsberger et al. (2016) GoT2D/T2D-GENES | ~145 MB |
| Scott et al. (2017) 1000G GWAS | ~155 MB |
| Scott et al. (2017) BMI-adjusted | ~138 MB |
| Mahajan et al. (2018a) ExomeChip European | ~6 MB |
| Mahajan et al. (2018a) Trans-Ethnic | ~9 MB |
| Mahajan et al. (2018b) Various analyses | 175-347 MB each |

**Multi-Ancestry Data (2022-2024):**

| Dataset | Size |
|---------|------|
| Mahajan et al. (2022) Multi-ancestry | ~465 MB |
| Mahajan et al. (2022) South Asian | ~196 MB |
| Mahajan et al. (2022) East Asian | ~169 MB |
| Mahajan et al. (2022) European | ~191 MB |
| **Suzuki et al. (2024) T2DGGI All-ancestry** | ~680 MB |
| Suzuki et al. (2024) African | ~360 MB |
| Suzuki et al. (2024) East Asian | ~320 MB |
| Suzuki et al. (2024) European | ~365 MB |
| Suzuki et al. (2024) Hispanic | ~350 MB |
| Suzuki et al. (2024) South Asian | ~270 MB |

**Access Methods:**
1. Direct downloads: https://diagram-consortium.org/downloads.html
2. T2D Knowledge Portal: https://t2d.hugeamp.org/

---

### 4. GIANT Consortium (Obesity Genetics)

| Field | Value |
|-------|-------|
| **URL** | https://giant-consortium.web.broadinstitute.org/ |
| **Content** | Anthropometric traits: BMI, height, WHR, waist/hip circumference |
| **Records** | ~5M individuals (2022 height GWAS); 27M SNP sites |
| **License** | Open access; citation required; allele frequencies excluded for privacy |
| **API** | None (bulk download only) |
| **Update Frequency** | Study-based releases |
| **Priority** | Tier 2 (Medium) |
| **Storage Estimate** | ~10 GB (compressed) |

**Traits Covered:**
- BMI: Body Mass Index
- Height
- WHR: Waist-Hip Ratio
- WC: Waist Circumference
- HIP: Hip Circumference
- Body shape PCs: Principal components

**Key Datasets:**

| Dataset | Sample Size | Description |
|---------|-------------|-------------|
| 2022 Height GWAS | ~5M | Multi-ancestry polygenic scores |
| 2018 Meta-analysis (Yengo) | ~700,000 | BMI and Height combined GWAS |
| 2018 Exome Array | Various | Protein-altering variants |
| 2017 Gene-Environment | Various | Smoking/activity interactions |
| 2015-2016 GIANT+UKB | Various | Sex/age-stratified analyses |
| 2010-2014 Historical | Various | Earlier GWAS releases |

**Data Downloads:**
- Broad Institute: https://giant-consortium.web.broadinstitute.org/index.php/GIANT_consortium_data_files
- Zenodo Archive: https://zenodo.org/record/1251813 (DOI: 10.5281/zenodo.1251813)

---

### 5. ICBP (Blood Pressure Genetics)

| Field | Value |
|-------|-------|
| **URL** | http://www.bloodpressuregenetics.org/ |
| **Content** | Systolic and diastolic blood pressure GWAS |
| **Records** | 535+ loci associated with blood pressure traits |
| **License** | Public data limited to p-values; full statistics require dbGaP application |
| **API** | None; LDHub platform for cross-study queries |
| **Update Frequency** | Study-based releases |
| **Priority** | Tier 3 (Lower - controlled access) |
| **Storage Estimate** | ~500 MB |

**Available Data:**

| Dataset | Description | Access |
|---------|-------------|--------|
| Cardio-MetaboChip (2016) | SBP and DBP summary statistics | Controlled (dbGaP) |
| UKB + ICBP Discovery | Combined UK Biobank + ICBP GWAS | LDHub/Request |
| Summary P-values | Genome-wide meta-analysis | Public download |

**Public Summary File (ICBP-summary-Nature.csv):**
- `rsid`: SNP ID (rs number)
- `chr.hg18`: Chromosome
- `pos.hg18`: Physical position (hg18)
- `pval.GC.SBP`: SBP p-values (genomic control corrected)
- `pval.GC.DBP`: DBP p-values (genomic control corrected)

**Access Methods:**
1. Controlled access via dbGaP (phs000585.v2.p1)
2. ICBP Steering Committee request
3. LDHub platform: http://ldsc.broadinstitute.org/ldhub/

---

### 6. HMDB (Human Metabolome Database)

| Field | Value |
|-------|-------|
| **URL** | https://hmdb.ca/ |
| **Content** | Small molecule metabolites in human body |
| **Records** | 220,945 metabolite entries |
| **License** | Free academic use; commercial requires agreement |
| **API** | Contact-based (email for credentials) |
| **Update Frequency** | Periodic (current: v5.0) |
| **Priority** | Tier 3 (Lower - contact required for API) |
| **Storage Estimate** | ~70 GB (full database with spectra) |

**Database Content:**
- Chemical structures, names, identifiers
- Biological roles and pathways
- Physiological concentrations
- Tissue/biofluid locations
- Disease associations
- Genetic associations
- MS/MS, GC-MS, and NMR spectra

**Available Downloads:**

| File Type | Size | Description |
|-----------|------|-------------|
| Metabolite Metabolizing Enzymes (Protein FASTA) | 1.87 MB | Protein sequences |
| Metabolite Metabolizing Enzymes (Gene FASTA) | 2.88 MB | Gene sequences |
| Metabolite Structures (SDF) | 92 MB | Chemical structures |
| All Metabolites (XML) | 910 MB | Complete metabolite data |
| All Proteins (XML) | 34.7 MB | Protein annotations |
| Urine Metabolites (XML) | 28.1 MB | Biofluid-specific |
| Serum Metabolites (XML) | 202 MB | Biofluid-specific |
| CSF Metabolites (XML) | 8.49 MB | Biofluid-specific |
| Saliva Metabolites (XML) | 16.5 MB | Biofluid-specific |
| Feces Metabolites (XML) | 61.2 MB | Biofluid-specific |
| MS/MS Spectra Images | 163 MB | Spectral data |
| NMR FID Files | 1.91 GB | Spectral data |
| All Peaklists Combined | 19.8 GB | Spectral data |
| All XML Spectra Combined | 42.8 GB | Spectral data |

**API Access Contacts:**
- Academic/Research: eponine@ualberta.ca or samackay@ualberta.ca
- Commercial: samackay@ualberta.ca

---

### 7. Metabolomics Workbench

| Field | Value |
|-------|-------|
| **URL** | https://www.metabolomicsworkbench.org/ |
| **Content** | Metabolomics study repository with compound data |
| **Records** | Thousands of studies with metabolite profiles |
| **License** | Open access |
| **API** | **Full REST API available** |
| **Update Frequency** | Continuous (study submissions) |
| **Priority** | Tier 1 (High - API available) |
| **Storage Estimate** | Variable (study-dependent) |

**REST API Base URL:** `https://www.metabolomicsworkbench.org/rest/`

**API Pattern:** `/[context]/[input_item]/[input_value]/[output_item]/[output_format]`

**API Contexts:**
- `compound`: Metabolite structures/identifiers
- `study`: Study metadata and results
- `refmet`: Standardized nomenclature
- `gene`: Gene information
- `protein`: Protein data
- `moverz`: Mass spectrometry searches

**Example Queries:**
```
/compound/pubchem_cid/5281365/smiles
/study/study_id/ST000001/summary
/refmet/name/Cholesterol/all
/moverz/REFMET/255.2/M+H/0.2/txt
```

**Output Formats:** JSON (default), TXT

**Client Libraries:**
- Python: `mwtab` package (PyPI)
- R: `metabolomicsWorkbenchR` (Bioconductor)

**Documentation:** https://www.metabolomicsworkbench.org/tools/MWRestAPIv1.2.pdf

---

### 8. GWAS Catalog (EBI)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Content** | Curated GWAS publications and associations |
| **Records** | 7,100+ studies; all cardiovascular/metabolic traits |
| **License** | Open (EMBL-EBI terms) |
| **API** | **Full REST API available** |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (High - comprehensive, API) |
| **Storage Estimate** | ~5 GB (full catalog) |

**REST API Base URL:** `https://www.ebi.ac.uk/gwas/rest/api`

**API Endpoints:**
- `/studies` - GWAS studies
- `/associations` - SNP-trait associations
- `/variants` - Variant information
- `/traits` - Trait ontology (EFO)

**Client Libraries:**
- R: `gwasrapidd` (CRAN)
- Python: `pandasGWAS` (PyPI)

**Documentation:** https://www.ebi.ac.uk/gwas/docs/programmatic-access

---

### 9. OpenGWAS (IEU)

| Field | Value |
|-------|-------|
| **URL** | https://opengwas.io/ |
| **Content** | Aggregated GWAS summary statistics |
| **Records** | 50,000+ GWAS datasets |
| **License** | Open access |
| **API** | **Full REST API available** |
| **Update Frequency** | Continuous |
| **Priority** | Tier 1 (High - fast batch queries) |
| **Storage Estimate** | Hosted (API access) |

**Features:**
- Fast batch API queries
- MR (Mendelian Randomization) analysis integration
- Cross-study comparisons

**Client Libraries:**
- R: `TwoSampleMR`, `gwasglue2`

---

### 10. MetaboLights (EBI)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/metabolights/ |
| **Content** | Cross-species metabolomics experiments |
| **Records** | Thousands of metabolomics studies |
| **License** | Open access |
| **API** | FTP access + Python package |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 (Medium) |
| **Storage Estimate** | Variable (study-dependent) |

**FTP Access:** `ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/`

**Python Package:** `metabolights-utils` (PyPI)

**Features:**
- Metabolite structures and spectra
- Study metadata as XML
- Cross-species coverage

---

## Integration Summary

### Data Volume Summary

| Database | Estimated Size | Primary Format |
|----------|---------------|----------------|
| CARDIoGRAMplusC4D | ~1.8 GB | Text/CSV |
| GLGC | ~2 GB | Text/CSV |
| DIAGRAM/T2DGGI | ~4 GB | Text/ZIP |
| GIANT | ~10 GB | Text/GZ |
| ICBP | ~500 MB | CSV |
| HMDB | ~70 GB (full) | XML/SDF/FASTA |
| GWAS Catalog | ~5 GB | JSON/TSV |
| **Total** | **~93 GB** | Mixed |

### API Availability Matrix

| Database | REST API | Bulk Download | Client Libraries |
|----------|----------|---------------|------------------|
| CARDIoGRAMplusC4D | No | Yes | No |
| GLGC | No | Yes | No |
| DIAGRAM | No | Yes | No |
| GIANT | No | Yes | No |
| ICBP | No | Partial | No |
| HMDB | Contact-based | Yes | No |
| Metabolomics Workbench | **Yes** | Yes | Python, R |
| GWAS Catalog | **Yes** | Yes | Python, R |
| OpenGWAS | **Yes** | Yes | R |
| MetaboLights | FTP | Yes | Python |

### License Summary

| Database | License Type | Commercial Use |
|----------|--------------|----------------|
| CARDIoGRAMplusC4D | Research-only | Restricted |
| GLGC | Open (citation required) | Likely allowed |
| DIAGRAM | Research-only | Restricted |
| GIANT | Open (citation required) | Likely allowed |
| ICBP | Controlled access | Application required |
| HMDB | Free academic | Separate agreement |
| Metabolomics Workbench | Open | Yes |
| GWAS Catalog | Open (EMBL-EBI terms) | Yes |
| OpenGWAS | Open | Yes |

---

## Integration Recommendations

### Priority 1: High (API available, open access)

| Source | Rationale |
|--------|-----------|
| GWAS Catalog | Comprehensive trait associations, REST API, curated |
| Metabolomics Workbench | Metabolite data with full REST API |
| OpenGWAS | 50K+ GWAS datasets, fast batch queries |

### Priority 2: Medium (bulk download, open)

| Source | Rationale |
|--------|-----------|
| GLGC | Essential lipid genetics, 923 loci |
| GIANT | Largest obesity genetics consortium |
| DIAGRAM/T2DGGI | Latest multi-ancestry T2D GWAS (2024) |
| CARDIoGRAMplusC4D | Primary CAD genetics resource |

### Priority 3: Lower (controlled access or contact-required)

| Source | Rationale |
|--------|-----------|
| ICBP | Blood pressure requires dbGaP application |
| HMDB | Full API requires contact |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-00-INDEX.md](./43-00-INDEX.md) | Parent index document |
| [43-DATA-SOURCES.md](../43-DATA-SOURCES.md) | Master data sources overview |
| [43-11-GENETICS-PRIMARY.md](./43-11-GENETICS-PRIMARY.md) | Related genetics sources |
| [43-42-PATHWAYS-METABOLISM.md](./43-42-PATHWAYS-METABOLISM.md) | Metabolic pathway databases |
| [44-ARCHITECTURE.md](../44-ARCHITECTURE.md) | Informs database design |
| [45-DATA-MODEL.md](../45-DATA-MODEL.md) | Informs entity structure |

---

## Key Citations

1. **CARDIoGRAMplusC4D**: Nikpay et al. (2015) Nat Genet 47:1121-1130
2. **GLGC**: Graham et al. (2021) Nature 600:675-679
3. **DIAGRAM/T2DGGI**: Suzuki et al. (2024) Nature
4. **GIANT**: Yengo et al. (2018) Hum Mol Genet 27:3641-3649
5. **ICBP**: Evangelou et al. (2018) Nat Genet 50:1412-1425
6. **HMDB**: Wishart et al. (2022) NAR 50:D622-D631

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial migration from research.old/data-sources-cardio-metabolic.md |
