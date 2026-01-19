# Cardiovascular and Metabolic Genetics Databases

Research document for Gene Platform integration of cardiovascular and metabolic genetics data sources.

---

## Table of Contents

1. [CARDIoGRAMplusC4D (Cardiovascular Genetics)](#1-cardiogramplusc4d-cardiovascular-genetics)
2. [Global Lipids Genetics Consortium (GLGC)](#2-global-lipids-genetics-consortium-glgc)
3. [DIAGRAM/DIAMANTE/T2DGGI (Diabetes Genetics)](#3-diagramdiamantet2dggi-diabetes-genetics)
4. [GIANT Consortium (Obesity Genetics)](#4-giant-consortium-obesity-genetics)
5. [ICBP (Blood Pressure Genetics)](#5-icbp-blood-pressure-genetics)
6. [HMDB (Human Metabolome Database)](#6-hmdb-human-metabolome-database)
7. [Additional Resources](#7-additional-resources)
8. [Integration Summary](#8-integration-summary)

---

## 1. CARDIoGRAMplusC4D (Cardiovascular Genetics)

### Overview

The CARDIoGRAMplusC4D (Coronary ARtery DIsease Genome wide Replication and Meta-analysis plus Coronary Artery Disease Genetics) consortium is a global collaborative effort combining data from large-scale genetic studies to identify risk loci for coronary artery disease (CAD) and myocardial infarction (MI).

### URLs

| Resource | URL |
|----------|-----|
| Main Website | https://cardiogramplusc4d.org/ |
| Data Downloads | https://cardiogramplusc4d.org/data-downloads/ |
| Publications | https://cardiogramplusc4d.org/publications/ |
| Knowledge Portal | https://www.kp4cd.org/node/620 |
| HeartBioPortal2.0 | https://heartbioportal.org/ |

### Available Datasets

| Dataset | Cases | Controls | Format | Size |
|---------|-------|----------|--------|------|
| CARDIoGRAM GWAS | 22,233 | 64,762 | Text delimited | 77 MB |
| C4D GWAS | 15,420 | 15,062 | Text delimited | 14 MB |
| CARDIoGRAMplusC4D Metabochip | 63,746 | 130,681 | Text delimited | 4.3 MB |
| 1000 Genomes (Additive) | 60,801 | 123,504 | Text delimited | 260 MB |
| 1000 Genomes (Recessive) | - | - | Text delimited | 190 MB |
| MI Additive | - | - | Text delimited | 389 MB |
| Exome Chip | 42,335 | 78,240 | Text delimited | 2.3 MB |
| Chromosome X-CAD | 43,000+ | 58,000+ | ZIP | 584 MB |
| UK Biobank Meta-analysis | - | - | Text delimited | 255 MB |

**Total Estimated Size**: ~1.8 GB

### Data Content

Each dataset includes:
- SNP ID (rsID)
- Chromosomal position
- Effect/reference alleles
- P-values
- Odds ratios
- Sample sizes
- Heterogeneity measures (I-squared)

### API Access

**No dedicated REST API available.**

Data access methods:
1. Direct file downloads from consortium website
2. Integration via HeartBioPortal2.0 (web-based query interface)
3. Knowledge Portal Network browsing

### License/Terms

- **Usage**: Research purposes only
- **Restriction**: Users must not attempt to de-identify individual subjects
- **Attribution**: Required citation of relevant publication and CARDIoGRAMplusC4D acknowledgment
- **License Type**: Custom research-use agreement (not standard open license)

### Contacts

- Jemma Hopewell (University of Oxford)
- Tim Assimes (Stanford University)

---

## 2. Global Lipids Genetics Consortium (GLGC)

### Overview

The GLGC is a worldwide collaboration investigating the genetic etiology of quantitative lipid traits. The consortium has identified over 923 genomic loci associated with lipid traits through GWAS involving more than 1.65 million individuals from globally diverse populations.

### URLs

| Resource | URL |
|----------|-----|
| Main Website | https://www.lipidgenetics.org/ |
| Data Downloads (2021) | https://csg.sph.umich.edu/willer/public/glgc-lipids2021/ |
| GitHub | https://github.com/Global-Lipids-Genetics |
| Knowledge Portal | https://www.kp4cd.org/node/784 |

### Available Datasets

| Dataset Category | Traits | Files | Link |
|------------------|--------|-------|------|
| Ancestry-Specific GWAS | HDL-C, LDL-C, nonHDL-C, TC, TG | 90 + README | [ancestry_specific](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific) |
| Trans-Ancestry GWAS | HDL-C, LDL-C, nonHDL-C, TC | 20 + README | [trans_ancestry](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry) |
| LDL-C Polygenic Scores | ALL, African, East Asian, European, Hispanic, South Asian | 12 + README | [prs_weights](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/prs_weights) |
| Trans-Ancestry Credible Sets | HDL-C, LDL-C, nonHDL-C, TC, TG | 5 + README | [credible_sets](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/credible_sets) |
| ChrX Summary Statistics | All lipid traits | 180 | [chrx_summary_stats](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/chrx_summary_stats) |
| Sex-Specific GWAS | All lipid traits | 20 | [sex_specific](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/sex_specific_summary_stats) |
| Sex & Ancestry-Specific | All lipid traits | 200 | [sex_ancestry](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/sex_and_ancestry_specific_summary_stats) |

**Total Files**: 500+ summary statistic files

### Lipid Traits Covered

- **HDL-C**: High-density lipoprotein cholesterol
- **LDL-C**: Low-density lipoprotein cholesterol
- **nonHDL-C**: Non-HDL cholesterol
- **TC**: Total cholesterol
- **TG**: Triglycerides

### API Access

**No dedicated REST API available.**

Data access:
1. Direct HTTP downloads from University of Michigan server
2. GitHub repository for analysis scripts

### License/Terms

- **Access**: Publicly available without explicit restrictions
- **Attribution**: Citation of Graham et al. (2021) Nature required
- **Usage**: Research purposes

### Contact

- Gina Peloso (gpeloso@bu.edu) - for joining consortium

---

## 3. DIAGRAM/DIAMANTE/T2DGGI (Diabetes Genetics)

### Overview

The DIAGRAM (DIAbetes Genetics Replication And Meta-analysis) consortium performs large-scale studies to characterize the genetic basis of type 2 diabetes. The consortium has evolved into DIAMANTE and T2DGGI (Type 2 Diabetes Global Genomics Initiative), assembling over 2.5 million participants for the largest multi-ancestry T2D GWAS.

### URLs

| Resource | URL |
|----------|-----|
| Main Website | https://diagram-consortium.org/ |
| Data Downloads | https://diagram-consortium.org/downloads.html |
| T2D Knowledge Portal | https://t2d.hugeamp.org/ |
| Knowledge Portal Network | https://kp4cd.org/node/167 |

### Available Datasets

#### Early DIAGRAM Analyses

| Dataset | Size |
|---------|------|
| Stage 1 GWAS Summary Statistics | ~32 MB |
| Stage 1 + Stage 2 Metabochip | ~3.3 MB |
| Sex-Specific Summary Statistics | ~7 MB |
| Stage 2 Metabochip Only | ~3.5 MB |

#### Expanded Analyses (2015-2018)

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

#### Multi-Ancestry Data (2022-2024)

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

**Total Estimated Size**: ~4 GB

### Data Format

- Tab-separated text files in ZIP archives
- PDF documentation guides for each dataset
- SHA1 and MD5 checksums provided

### API Access

**No dedicated REST API available.**

Access methods:
1. Direct file downloads from consortium website
2. Type 2 Diabetes Knowledge Portal (web interface)

### License/Terms

- **Restriction**: Users must not attempt to de-identify individual subjects
- **Rationale**: Sample size and data precision could enable re-identification
- **Attribution**: Required citation of relevant publications

---

## 4. GIANT Consortium (Obesity Genetics)

### Overview

The Genetic Investigation of ANthropometric Traits (GIANT) consortium is an international collaboration identifying genetic loci that modulate human body size and shape, including height, BMI, and waist-related measures.

### URLs

| Resource | URL |
|----------|-----|
| Main Website | https://giant-consortium.web.broadinstitute.org/ |
| Data Files | https://giant-consortium.web.broadinstitute.org/index.php/GIANT_consortium_data_files |
| Zenodo Archive | https://zenodo.org/record/1251813 |
| Knowledge Portal | https://kp4cd.org/node/180 |

### Available Traits

- **BMI**: Body Mass Index
- **Height**
- **WHR**: Waist-Hip Ratio
- **WC**: Waist Circumference
- **HIP**: Hip Circumference
- **Body shape PCs**: Principal components

### Key Datasets

| Dataset | Sample Size | Description |
|---------|-------------|-------------|
| 2022 Height GWAS | ~5M | Multi-ancestry polygenic scores |
| 2018 Meta-analysis (Yengo) | ~700,000 | BMI and Height combined GWAS |
| 2018 Exome Array | Various | Protein-altering variants |
| 2017 Gene-Environment | Various | Smoking/activity interactions |
| 2015-2016 GIANT+UKB | Various | Sex/age-stratified analyses |
| 2010-2014 Historical | Various | Earlier GWAS releases |

**Total Size**: All 27M sites available via Zenodo (~10+ GB compressed)

### Data Format

- Gzipped text files (.txt.gz)
- P-values and effect directions for ~2 million SNPs
- **Note**: Allele frequency data excluded to protect privacy

### API Access

**No dedicated REST API available.**

Access methods:
1. Direct downloads from Broad Institute portal
2. Zenodo repository (DOI: 10.5281/zenodo.1251813)
3. Joel Hirschhorn Lab archives

### License/Terms

- **Access**: Publicly available
- **Citation**: Required for each dataset (Nature, Nature Genetics publications)
- **Privacy**: Allele frequencies excluded

---

## 5. ICBP (Blood Pressure Genetics)

### Overview

The International Consortium for Blood Pressure (ICBP) investigates blood pressure genetics, formed by CHARGE-BP and GBPGEN parent consortia. Major findings include 535+ loci associated with blood pressure traits.

### URLs

| Resource | URL |
|----------|-----|
| Blood Pressure Genetics | http://www.bloodpressuregenetics.org/ |
| ICBP Study (dbGaP) | https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000585.v2.p1 |
| EGA Archive | https://ega-archive.org/studies/phs000585 |
| Knowledge Portal | https://www.kp4cd.org/node/679 |
| OmicsDI | https://www.omicsdi.org/dataset/dbgap/phs000585 |

### Available Data

| Dataset | Description | Access |
|---------|-------------|--------|
| Cardio-MetaboChip (2016) | SBP and DBP summary statistics | Controlled (dbGaP) |
| UKB + ICBP Discovery | Combined UK Biobank + ICBP GWAS | LDHub/Request |
| Summary P-values | Genome-wide meta-analysis | Public download |

### Data Content

The summary file (ICBP-summary-Nature.csv) contains:
- `rsid`: SNP ID (rs number)
- `chr.hg18`: Chromosome
- `pos.hg18`: Physical position (hg18)
- `pval.GC.SBP`: SBP p-values (genomic control corrected)
- `pval.GC.DBP`: DBP p-values (genomic control corrected)

### API Access

**No dedicated REST API available.**

Access methods:
1. **Controlled access via dbGaP** (phs000585.v2.p1) - requires application
2. **ICBP Steering Committee request** - contact corresponding authors
3. **LDHub platform** (http://ldsc.broadinstitute.org/ldhub/)

### License/Terms

- **Public data**: Limited to p-values only
- **Full statistics**: Controlled access via dbGaP
- **Effect sizes/directions**: Requires ICBP approval

---

## 6. HMDB (Human Metabolome Database)

### Overview

The Human Metabolome Database (HMDB) is a freely available database containing detailed information about small molecule metabolites found in the human body. Current version: HMDB 5.0 with 220,945 metabolite entries.

### URLs

| Resource | URL |
|----------|-----|
| Main Website | https://hmdb.ca/ |
| Downloads | https://hmdb.ca/downloads |
| API Information | https://hmdb.ca/simple/api |
| About | https://hmdb.ca/about |

### Database Content

- **220,945** metabolite entries
- Chemical structures, names, identifiers
- Biological roles and pathways
- Physiological concentrations
- Tissue/biofluid locations
- Disease associations
- Genetic associations
- MS/MS, GC-MS, and NMR spectra

### Available Downloads

#### Protein/Gene Sequences (FASTA)

| File | Size | Released |
|------|------|----------|
| Metabolite Metabolizing Enzymes (Protein) | 1.87 MB | 2021-11-02 |
| Metabolite Metabolizing Enzymes (Gene) | 2.88 MB | 2021-11-02 |

#### Metabolite Structures (SDF)

| File | Size | Released |
|------|------|----------|
| Metabolite Structures | 92 MB | 2021-11-02 |

#### Metabolite/Protein Data (XML)

| File | Size | Released |
|------|------|----------|
| All Metabolites | 910 MB | 2021-11-17 |
| All Proteins | 34.7 MB | 2021-11-09 |
| Urine Metabolites | 28.1 MB | 2021-10-24 |
| Serum Metabolites | 202 MB | 2021-10-24 |
| CSF Metabolites | 8.49 MB | 2021-10-24 |
| Saliva Metabolites | 16.5 MB | 2021-10-24 |
| Feces Metabolites | 61.2 MB | 2021-10-24 |
| Sweat Metabolites | 3.27 MB | 2021-10-24 |

#### Spectra Files

| Type | Size | Released |
|------|------|----------|
| MS/MS Spectra Images | 163 MB | 2023-07-01 |
| NMR FID Files | 1.91 GB | 2023-07-01 |
| GC-MS Peaklists (Predicted) | 323 MB | 2023-07-01 |
| MS-MS Peaklists (Predicted) | 186 MB | 2023-07-01 |
| All Peaklists Combined | 19.8 GB | 2023-07-01 |
| All XML Spectra Combined | 42.8 GB | 2023-07-01 |

**Total Estimated Size**: ~70 GB (full database with spectra)

### API Access

**Contact-based API access:**

| User Type | Contact |
|-----------|---------|
| Academic/Research | eponine@ualberta.ca or samackay@ualberta.ca |
| Commercial | samackay@ualberta.ca |

**Note**: No public REST API documentation. API access requires contacting the HMDB team for credentials and specifications.

### License/Terms

- **Access**: Freely available for academic use
- **Citation**: HMDB 5.0 publication (NAR 2022)
- **Commercial**: Requires separate agreement

---

## 7. Additional Resources

### 7.1 Metabolomics Workbench

**URL**: https://www.metabolomicsworkbench.org/

**REST API**: Full programmatic access available

**Base URL**: `https://www.metabolomicsworkbench.org/rest/`

**API Pattern**: `/[context]/[input_item]/[input_value]/[output_item]/[output_format]`

**Contexts**:
- `compound`: Metabolite structures/identifiers
- `study`: Study metadata and results
- `refmet`: Standardized nomenclature
- `gene`: Gene information
- `protein`: Protein data
- `moverz`: Mass spectrometry searches

**Example Queries**:
```
/compound/pubchem_cid/5281365/smiles
/study/study_id/ST000001/summary
/refmet/name/Cholesterol/all
/moverz/REFMET/255.2/M+H/0.2/txt
```

**Output Formats**: JSON (default), TXT

**Client Libraries**:
- Python: `mwtab` package (PyPI)
- R: `metabolomicsWorkbenchR` (Bioconductor)

**Documentation**: https://www.metabolomicsworkbench.org/tools/MWRestAPIv1.2.pdf

### 7.2 MetaboLights (EBI)

**URL**: https://www.ebi.ac.uk/metabolights/

**FTP Access**: `ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/`

**Python Package**: `metabolights-utils` (PyPI)

**Features**:
- Cross-species metabolomics experiments
- Metabolite structures and spectra
- Study metadata as XML

### 7.3 GWAS Catalog

**URL**: https://www.ebi.ac.uk/gwas/

**REST API**: https://www.ebi.ac.uk/gwas/rest/api

**Programmatic Access**: https://www.ebi.ac.uk/gwas/docs/programmatic-access

**Client Libraries**:
- R: `gwasrapidd` (CRAN)
- Python: `pandasGWAS` (PyPI)

**API Endpoints**:
- `/studies` - GWAS studies
- `/associations` - SNP-trait associations
- `/variants` - Variant information
- `/traits` - Trait ontology (EFO)

**Coverage**: 7,100+ studies, all cardiovascular/metabolic traits included

### 7.4 OpenGWAS

**URL**: https://opengwas.io/

**Features**:
- 50,000+ GWAS summary datasets
- Fast batch API queries
- MR analysis integration

**Client Libraries**:
- R: `TwoSampleMR`, `gwasglue2`

---

## 8. Integration Summary

### Data Volume Summary

| Database | Estimated Size | Primary Format |
|----------|---------------|----------------|
| CARDIoGRAMplusC4D | ~1.8 GB | Text/CSV |
| GLGC | ~2 GB | Text/CSV |
| DIAGRAM/T2DGGI | ~4 GB | Text/ZIP |
| GIANT | ~10 GB | Text/GZ |
| ICBP | ~500 MB | CSV |
| HMDB | ~70 GB (full) | XML/SDF/FASTA |
| **Total** | **~88 GB** | Mixed |

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

### License Summary

| Database | License Type | Commercial Use |
|----------|--------------|----------------|
| CARDIoGRAMplusC4D | Research-only | Restricted |
| GLGC | Open (citation required) | Likely allowed |
| DIAGRAM | Research-only | Restricted |
| GIANT | Open (citation required) | Likely allowed |
| ICBP | Controlled access | Application required |
| HMDB | Free academic, commercial contact | Separate agreement |
| Metabolomics Workbench | Open | Yes |
| GWAS Catalog | Open (EMBL-EBI terms) | Yes |

### Recommended Integration Priority

1. **High Priority** (API available, open access):
   - GWAS Catalog - comprehensive trait associations
   - Metabolomics Workbench - metabolite data with REST API
   - OpenGWAS - summary statistics access

2. **Medium Priority** (bulk download, open):
   - GLGC - lipid genetics
   - GIANT - obesity genetics
   - DIAGRAM - diabetes genetics (latest T2DGGI)

3. **Lower Priority** (controlled access or contact-required):
   - ICBP - blood pressure (requires dbGaP application)
   - HMDB - full API (requires contact)
   - CARDIoGRAMplusC4D - cardiovascular (download only)

### Key Citations

1. **CARDIoGRAMplusC4D**: Nikpay et al. (2015) Nat Genet 47:1121-1130
2. **GLGC**: Graham et al. (2021) Nature 600:675-679
3. **DIAGRAM**: Suzuki et al. (2024) Nature (T2DGGI)
4. **GIANT**: Yengo et al. (2018) Hum Mol Genet 27:3641-3649
5. **ICBP**: Evangelou et al. (2018) Nat Genet 50:1412-1425
6. **HMDB**: Wishart et al. (2022) NAR 50:D622-D631

---

*Document generated: 2026-01-19*
*For Gene Platform integration planning*
