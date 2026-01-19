# Cancer and Oncology Genetics Databases

Research compilation for Gene Platform integration covering somatic mutation databases, clinical interpretation resources, and hereditary cancer variant databases.

---

## 1. COSMIC (Catalogue Of Somatic Mutations In Cancer)

### Overview
COSMIC is the world's largest and most comprehensive resource for exploring the impact of somatic mutations in human cancer. Maintained by the Wellcome Sanger Institute.

### URL
- **Main Portal**: https://cancer.sanger.ac.uk/cosmic
- **Downloads**: https://cancer.sanger.ac.uk/cosmic/download
- **Knowledge Base**: https://www.cosmickb.org/

### Content
- **Somatic coding mutations**: >38 million mutations
- **Cancer samples**: >1.4 million samples
- **Genomic variants**: >25 million across >6,800 cancer types
- **Gene fusions**: Comprehensive fusion database
- **Mutational signatures**: Catalog of mutational processes
- **Drug resistance data**: Clinically relevant resistance mutations
- **Copy number variants**: Amplifications and deletions
- **Non-coding mutations**: Regulatory region variants

### Data Products
| Product | Description |
|---------|-------------|
| COSMIC Core | Complete somatic mutation dataset |
| Cell Lines Project | Cancer cell line mutations |
| Cancer Mutation Census (CMC) | Clinically actionable mutations |
| Actionability | Drug-variant associations |
| Cancer Gene Census (CGC) | Curated cancer driver genes |

### API Access
- **GA4GH Beacon**: Standard variant query interface
- **NIH Clinical Tables API**: Limited access (v89, GRCh37 only)
  - URL: https://clinicaltables.nlm.nih.gov/apidoc/cosmic/v3/doc.html
  - Note: No pagination support due to licensing
- **Download Files**: Structured flat files (TSV/CSV)

### Download Formats
- Tab-separated files (TSV)
- VCF files (GRCh37 and GRCh38)
- Complete datasets by genome version
- 3 previous release versions available

### License
| User Type | Access |
|-----------|--------|
| Academic (research) | Free registration required |
| Academic (patient testing) | License fee required |
| Commercial | License via QIAGEN required |

- **Commercial Access**: https://digitalinsights.qiagen.com/products-overview/cosmic/
- Registration required for all downloads
- Contact: cosmic@sanger.ac.uk

### Data Size
- ~1 GB RAM required to load database
- Quarterly updates (current: v103, November 2025)
- Supports GRCh37 and GRCh38

### Integration Notes
- Python access via `gget cosmic` package
- Google Cloud Life Sciences: Public dataset available
- Commercial licensing may limit redistribution

---

## 2. TCGA / GDC (Genomic Data Commons)

### Overview
The NCI Genomic Data Commons (GDC) is the primary repository for TCGA (The Cancer Genome Atlas) data, providing harmonized genomic and clinical data from major cancer research programs.

### URL
- **GDC Portal**: https://portal.gdc.cancer.gov/
- **GDC API**: https://gdc.cancer.gov/access-data
- **Documentation**: https://gdc.cancer.gov/developers/gdc-application-programming-interface-api

### Content
- **TCGA**: 33 cancer types, >11,000 patients
- **TARGET**: Pediatric cancers
- **CPTAC**: Clinical proteomic tumor analysis
- **CGCI**: Cancer Genome Characterization Initiative
- **HCMI**: Human Cancer Model Initiative

### Data Types
| Data Type | Description |
|-----------|-------------|
| WGS/WXS | Whole genome/exome sequencing |
| RNA-Seq | Gene expression data |
| miRNA-Seq | MicroRNA expression |
| DNA Methylation | Epigenetic profiles (SeSAMe pipeline) |
| SNP Arrays | Copy number data |
| RPPA | Protein expression |
| Clinical | Patient demographics, outcomes |

### API Access
- **GDC REST API**: Full programmatic access
  - Search, filter, download, upload
  - JSON responses
  - Authentication for controlled data
- **GraphQL**: Query interface
- **R Package**: `GenomicDataCommons` (Bioconductor)
- **Python**: `gdclient` and various packages

### Download Methods
1. **GDC Data Portal**: Web-based browsing and cart system
2. **GDC Data Transfer Tool (DTT)**: Command-line for large files
   - Version: 2.3.0
   - Supports: Linux, macOS (Intel/Silicon), Windows
   - Features: Resume capability, high-throughput
3. **API Downloads**: Programmatic access
4. **TCGADownloadHelper**: Simplified extraction tool (2025)
   - GitHub: https://github.com/alex-baumann-ur/TCGADownloadHelper

### License
| Data Type | Access |
|-----------|--------|
| Open Access | No authentication required |
| Controlled Access | dbGaP authorization + eRA Commons |

- Data downloads are free
- CC0 or custom NCI terms depending on dataset
- May require data use agreements for controlled data

### Data Size
- **Total**: ~2.5 petabytes
- Current release: Data Release 43.0 (May 2025)
- Individual files: 50 MB to 200-300 GB

### Integration Notes
- Standard file formats (BAM, VCF, MAF, TSV)
- GDC harmonization pipelines ensure consistency
- Metadata-rich with standardized ontologies

---

## 3. Cancer Gene Census (CGC)

### Overview
A high-confidence, expertly curated catalog of genes causally implicated in cancer, maintained as part of COSMIC by the Wellcome Sanger Institute.

### URL
- **Main Page**: https://cancer.sanger.ac.uk/census
- **Sanger Data Page**: https://www.sanger.ac.uk/data/cancer-gene-census/

### Content
- **Tier 1 Genes**: Strong functional + mutational evidence
- **Tier 2 Genes**: Mutational patterns without full functional characterization
- Documented mutation mechanisms
- Associated cancer types
- Literature references

### Data Fields
| Field | Description |
|-------|-------------|
| Gene Symbol | HGNC symbol |
| Tier | 1 (high confidence) or 2 |
| Hallmark | Cancer hallmarks affected |
| Mutation Types | Point, amplification, translocation, etc. |
| Role in Cancer | Oncogene, TSG, or both |
| Cancer Types | Associated malignancies |

### Acceptance Criteria
- Literature-based evidence required
- Independent confirmation from 2+ groups
- Clear mutation patterns in specified diseases
- Statistical interpretations excluded (conservative approach)

### API Access
- Part of COSMIC data downloads
- Same access restrictions as COSMIC
- Can export as CSV/TSV from web interface

### Download Formats
- CSV export from web interface
- TSV in COSMIC downloads
- Searchable/sortable web tables

### License
- Same as COSMIC (see above)
- Academic free with registration
- Commercial via QIAGEN

### Data Size
- Hundreds of genes (curated subset of genome)
- Regular updates by COSMIC curation team

### Integration Notes
- Often used as filter/annotation for variant analysis
- Compare with OncoKB gene list (1,216 genes as of Dec 2025)
- Foundation for many cancer gene panels

---

## 4. OncoKB

### Overview
Memorial Sloan Kettering Cancer Center's (MSK) Precision Oncology Knowledge Base. FDA-recognized human genetic variant database for clinical variant interpretation.

### URL
- **Main Portal**: https://www.oncokb.org/
- **API Documentation**: https://api.oncokb.org/oncokb-website/api
- **GitHub**: https://github.com/oncokb/oncokb

### Content (as of Dec 2025)
| Metric | Count |
|--------|-------|
| Genes | 978 |
| Alterations | 10,623 |
| Cancer Types | 144 |
| Drugs | 159 |
| Cancer Gene List | 1,216 genes |

### Evidence Levels
| Level | Description |
|-------|-------------|
| Level 1 | FDA-recognized biomarker |
| Level 2 | Standard care biomarker |
| Level 3A | Compelling clinical evidence |
| Level 3B | Standard care in another cancer type |
| Level 4 | Compelling biological evidence |
| Level R1 | Standard care resistance |
| Level R2 | Compelling resistance evidence |

### API Access
- **Authentication**: Token required (via license)
- **Demo Instance**: https://demo.oncokb.org (BRAF, TP53, ROS1 only)
- **Production**: https://www.oncokb.org (licensed users)

### API Endpoints
| Endpoint Type | Description |
|---------------|-------------|
| Protein Change | Gene + alteration + tumor type |
| Copy Number | Amplification, deletion, gain, loss |
| Structural Variants | Fusions, rearrangements |
| Genomic Coordinates | TCGA MAF or HGVS format |

### Download Options
- **Public Downloads** (no registration):
  - Cancer Gene List
  - All Curated Genes List
  - Biomarker-Drug Association List
- **API**: Recommended for variant annotation (always current)
- **Full Database**: Last resort for licensed users only

### License
| Use Case | Cost |
|----------|------|
| Academic research | Free (with registration) |
| Commercial research | Fee based on company size |
| Clinical/hospital use | Fee based on report volume |
| AI/ML training | Prohibited |
| AI/ML benchmarking | Requires explicit permission |

- Contact: contact@oncokb.org

### Data Size
- Focused curation (quality over quantity)
- Updated frequently (17 novel biomarkers added in 2024)
- Germline variants expected in 2025

### Integration Notes
- Strongly recommends API over downloads
- Downloads can become outdated quickly
- Umbrella terms use logic that requires API

---

## 5. CIViC (Clinical Interpretation of Variants in Cancer)

### Overview
An open-access, community-driven knowledgebase for expert crowdsourcing of clinical variant interpretations. Maintained by The McDonnell Genome Institute at Washington University School of Medicine.

### URL
- **Main Portal**: https://civicdb.org/
- **Documentation**: https://civic.readthedocs.io/
- **GitHub**: https://github.com/griffithlab/civic-v2
- **AWS Open Data**: https://registry.opendata.aws/civic/

### Content
| Metric | Count |
|--------|-------|
| Variants | >3,200 |
| Genes | >470 |
| Publications | >3,100 |
| Contributors | >300 |

### Evidence Types
| Type | Description |
|------|-------------|
| Predictive | Drug response predictions |
| Diagnostic | Supports diagnosis |
| Prognostic | Outcome predictions |
| Predisposing | Germline risk factors |
| Oncogenic | Functional oncogenicity |

### API Access
- **GraphQL API**: Primary interface
- **GraphiQL**: Interactive query interface
- **REST API v2**: Current version
- **REST API v1**: Deprecated
- **CIViCpy**: Python SDK and analysis toolkit

### Download Formats
- **Nightly TSV/CSV dumps**: Available on AWS
- **Monthly releases**: Full data exports
- **API queries**: Real-time data access
- Formats: JSON, TSV, VCF

### License
- **Data**: CC0 1.0 Universal (Public Domain)
- **Source Code**: MIT License
- Freely queryable, downloadable, reusable, redistributable

### Data Size
- Focused expert-curated resource
- Nightly updates available
- Growing through community contributions

### Integration Notes
- Completely open with no restrictions
- MCP integration available (CIViC MCP, 2025)
- Excellent for research and clinical decision support
- Community contribution model

---

## 6. Hereditary Cancer Databases

### 6.1 ClinVar (Hereditary Cancer Variants)

#### Overview
NCBI's free public database of human genetic variants and disease relationships. Primary resource for hereditary cancer variant classification.

#### URL
- **Main Portal**: https://www.ncbi.nlm.nih.gov/clinvar/
- **FTP**: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/
- **Simple ClinVar**: https://simple-clinvar.broadinstitute.org/

#### Content
- >3 million variants
- >2,800 submitting organizations
- Germline and somatic classifications
- Oncogenicity classifications (new 2025)

#### Hereditary Cancer Genes Covered
- BRCA1, BRCA2 (breast/ovarian)
- MLH1, MSH2, MSH6, PMS2, EPCAM (Lynch syndrome)
- APC, MUTYH (polyposis syndromes)
- TP53 (Li-Fraumeni)
- PTEN (Cowden syndrome)
- CDH1 (hereditary diffuse gastric cancer)
- PALB2, CHEK2, ATM (moderate penetrance)

#### API Access
- **E-utilities**: NCBI programmatic access
- **Clinical Tables API**: https://clinicaltables.nlm.nih.gov/api/variants/v4/search
- FTP downloads in multiple formats

#### Download Formats
- VCF (GRCh37, GRCh38)
- XML (VCV, RCV formats)
- Tab-delimited summaries
- Weekly and monthly releases

#### License
- **Public domain** (US Government work)
- Free to use, download, redistribute

#### Data Size
- Multi-gigabyte download files
- Monthly release cycle

---

### 6.2 InSiGHT Database

#### Overview
International Society for Gastrointestinal Hereditary Tumours database focusing on Lynch syndrome and related hereditary GI cancer syndromes.

#### URL
- **Main Portal**: https://www.insight-group.org/mutations
- **Organization**: https://www.insight-group.org/

#### Content
- **Primary Focus**: MMR genes (MLH1, MSH2, MSH6, PMS2)
- ~3,000 unique germline variants
- Gene distribution: MLH1 (40%), MSH2 (34%), MSH6 (18%), PMS2 (8%)
- Also includes: APC, MUTYH variants

#### Classification System
- 5-tiered IARC classification
- Linked to clinical recommendations
- Based on:
  - Variant characteristics
  - Family data
  - Functional assay results

#### API Access
- Web-based query interface
- Classifications mirrored to ClinVar
- Export capabilities

#### Related Resources
- **PLSD (Prospective Lynch Syndrome Database)**
  - Cancer risk estimates by age, gender, gene
  - Published: https://hccpjournal.biomedcentral.com/

#### License
- Academic access (registration may be required)
- Data sharing with ClinVar

---

### 6.3 BRCA Exchange

#### Overview
Global resource for BRCA1 and BRCA2 variant data aggregation and expert classification, part of the GA4GH BRCA Challenge.

#### URL
- **Main Portal**: https://brcaexchange.org/
- **Variants**: https://brcaexchange.org/variants
- **GitHub**: https://github.com/BRCAChallenge/brca-exchange

#### Content
- >20,000 unique BRCA1/BRCA2 variants
- >6,100 expert-classified variants
- ~3,700 pathogenic variants

#### Data Sources Aggregated
| Source | Type |
|--------|------|
| ClinVar | Clinical submissions |
| LOVD | Locus-specific database |
| BIC | Breast Cancer Information Core |
| ExAC/gnomAD | Population frequencies |
| 1000 Genomes | Population data |
| ENIGMA | Expert classifications |

#### API Access
- REST API available
- PostgreSQL database backend
- Version history tracking
- VR (Variation Representation) specification support

#### Download Options
- Full variant data via API or web portal
- Historical data releases (tarballs)
- Pipeline intermediate outputs included
- Monthly updates with change tracking

#### License
- **Open access**
- Freely queryable, downloadable, redistributable
- Data use policy permits broad reuse

---

### 6.4 Additional Hereditary Cancer Resources

#### NCCN Guidelines
- **URL**: https://www.nccn.org/guidelines/
- Current: Version 2.2026
- Covers: BRCA-Related, Li-Fraumeni, Cowden/PHTS, Lynch syndrome
- Testing criteria and management recommendations

#### GeneReviews
- **URL**: https://www.ncbi.nlm.nih.gov/books/NBK1247/
- BRCA1/BRCA2-Associated HBOC (updated March 2025)
- Comprehensive clinical descriptions
- Free access

---

## Summary Comparison Table

| Database | Focus | Variants | License | API | Best For |
|----------|-------|----------|---------|-----|----------|
| COSMIC | Somatic mutations | 38M+ | Academic free, Commercial paid | Limited | Somatic variant annotation |
| GDC/TCGA | Multi-omics cancer | 2.5 PB | Open/Controlled | Full REST | Research data access |
| Cancer Gene Census | Cancer driver genes | ~700 genes | Same as COSMIC | Via COSMIC | Gene filtering |
| OncoKB | Clinical actionability | 10K+ alterations | Academic free, Clinical paid | Full REST | Treatment decisions |
| CIViC | Clinical interpretation | 3K+ variants | CC0 (open) | GraphQL | Research, clinical support |
| ClinVar | Germline classification | 3M+ | Public domain | E-utilities | Variant classification |
| InSiGHT | Lynch syndrome | 3K+ variants | Academic | Web | Lynch syndrome |
| BRCA Exchange | BRCA1/2 | 20K+ variants | Open | REST | BRCA interpretation |

---

## Recommended Integration Priority

### Tier 1 (Essential)
1. **ClinVar** - Foundation for all variant classification
2. **CIViC** - Open license, clinical interpretations
3. **GDC/TCGA** - Research data backbone

### Tier 2 (High Value)
4. **OncoKB** - Clinical actionability (licensing required for clinical use)
5. **COSMIC** - Somatic mutation reference (licensing considerations)
6. **BRCA Exchange** - BRCA-specific expertise

### Tier 3 (Specialized)
7. **Cancer Gene Census** - Gene-level filtering
8. **InSiGHT** - Lynch syndrome focus

---

## Technical Integration Considerations

### Data Harmonization Challenges
- Variant nomenclature differences (HGVS, protein, genomic)
- Reference genome versions (GRCh37 vs GRCh38)
- Classification system differences (ACMG, AMP/ASCO/CAP)
- Update frequency variations

### Recommended Tools
- **VEP** (Variant Effect Predictor) - Ensembl annotation
- **ANNOVAR** - Multi-database annotation
- **OpenCRAVAT** - Modular annotation system
- **ClinGen Allele Registry** - Variant normalization

### API Rate Limits
| Database | Limit Notes |
|----------|-------------|
| COSMIC | Registration required |
| GDC | Generally unrestricted |
| OncoKB | Token-based, contact for limits |
| CIViC | Open, reasonable use |
| ClinVar | E-utilities guidelines apply |

---

*Research compiled: January 2026*
*For Gene Platform data integration planning*
