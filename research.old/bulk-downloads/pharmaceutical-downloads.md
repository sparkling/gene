# Pharmaceutical and Drug Database Bulk Downloads

This document provides comprehensive information on bulk download methods for major pharmaceutical and drug databases, including direct URLs, file sizes, licensing requirements, and recommended processing approaches.

---

## Table of Contents

1. [DrugBank](#1-drugbank)
2. [ChEMBL](#2-chembl)
3. [PubChem](#3-pubchem)
4. [PharmGKB](#4-pharmgkb)
5. [DGIdb](#5-dgidb)
6. [Open Targets](#6-open-targets)
7. [BindingDB](#7-bindingdb)
8. [OpenFDA](#8-openfda)
9. [DailyMed](#9-dailymed)
10. [RxNorm](#10-rxnorm)

---

## 1. DrugBank

### Overview
DrugBank is the gold standard knowledge resource for drug, drug-target, and pharmaceutical information. DrugBank 6.0 (2024) contains 4,563 FDA-approved drugs, 6,231 investigational drugs, 1,413,413 drug-drug interactions, and 2,475 drug-food interactions.

### Download URLs
- **Releases Page**: https://go.drugbank.com/releases
- **Latest Release**: https://go.drugbank.com/releases/latest
- **Data Packages**: https://go.drugbank.com/data_packages
- **Download Help**: https://go.drugbank.com/releases/help

### Available Files
| File | Description | Estimated Size |
|------|-------------|----------------|
| `drugbank_all_full_database.xml.zip` | Complete database in XML format | ~500MB compressed |
| `drugbank_all_structures.sdf.zip` | Chemical structures in SDF format | ~100MB |
| `drugbank_all_target_sequences.fasta.zip` | Target protein sequences | ~10MB |
| `drugbank_all_drug_sequences.fasta.zip` | Drug sequences (biologics) | ~5MB |
| SQL dumps (MySQL/PostgreSQL/MS SQL) | Relational database format | ~200MB |

### Licensing
- **Academic License**: Free for researchers at academic/non-profit institutions
- **Commercial License**: Required for commercial use - contact DrugBank
- **Registration**: Required for all downloads

### Command Line Download
```bash
# Requires account credentials
curl -Lfv -o drugbank_full.zip \
  -u YOUR_EMAIL:YOUR_PASSWORD \
  https://go.drugbank.com/releases/latest/downloads/all-full-database
```

### Recommended Processing
- **Python**: Use `dbparser` R package or custom XML parsing with `lxml`
- **Database**: Load CSV exports directly into SQL databases
- **Memory**: XML file is large - use streaming parsers (SAX) for low-memory processing
- **Format**: Structure data available in SMILES, SDF, MOL, PDB, InChI, and InChIKey formats

---

## 2. ChEMBL

### Overview
ChEMBL is a manually curated database of bioactive molecules with drug-like properties. ChEMBL 35 (December 2024) contains 2,496,335 compounds, 1,740,546 assays, 21,123,501 activities, and 16,003 targets.

### Download URLs
- **FTP (Latest)**: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
- **FTP (Releases)**: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/
- **ChEMBL 35**: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/
- **Documentation**: https://chembl.gitbook.io/chembl-interface-documentation/downloads

### Available Files
| File | Description | Size |
|------|-------------|------|
| `chembl_35_sqlite.tar.gz` | SQLite database | ~5GB compressed, ~35GB extracted |
| `chembl_35_mysql.tar.gz` | MySQL dump | ~4GB compressed |
| `chembl_35_postgresql.tar.gz` | PostgreSQL dump | ~4GB compressed |
| `chembl_35.sdf.gz` | All structures in SDF | ~1.5GB |
| `chembl_35.fa.gz` | Target sequences FASTA | ~10MB |
| `chembl_35_release_notes.txt` | Release documentation | ~100KB |

### Licensing
- **License**: Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0)
- **Registration**: Not required for download

### Command Line Download
```bash
# Download SQLite database
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/chembl_35_sqlite.tar.gz

# Extract
tar -xzf chembl_35_sqlite.tar.gz
```

### Python Package (Recommended)
```bash
pip install chembl-downloader
```

```python
from chembl_downloader import download_extract_sqlite

# Automatically downloads and extracts the latest SQLite database
path = download_extract_sqlite()
```

### Recommended Processing
- **Disk Space**: Ensure ~35GB free space for SQLite import
- **Database**: SQLite version is easiest for local analysis
- **Python**: Use `chembl_downloader` for reproducible workflows
- **RDKit**: Use with RDKit for chemical structure analysis

---

## 3. PubChem

### Overview
PubChem is the world's largest collection of freely accessible chemical information. As of 2025, it contains 119 million compounds, 322 million substances, and 295 million bioactivities from over 1,000 data sources.

### Download URLs
- **FTP Root**: https://ftp.ncbi.nlm.nih.gov/pubchem/
- **Compound**: https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/
- **Substance**: https://ftp.ncbi.nlm.nih.gov/pubchem/Substance/
- **BioAssay**: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/
- **RDF**: https://ftp.ncbi.nlm.nih.gov/pubchem/RDF/
- **Documentation**: https://pubchem.ncbi.nlm.nih.gov/docs/downloads

### Directory Structure
```
/pubchem/
  /Compound/           # Full compound data dump
    /CURRENT-Full/     # Complete current release
    /Daily/            # Daily incremental updates
    /Extras/           # Simplified data files
  /Substance/          # Substance records from submitters
  /Compound3D/         # 3D conformer structures
  /Bioassay/           # Bioactivity data
  /RDF/                # RDF format data
  /specifications/     # File format specifications
```

### Key Data Files
| File | Location | Description |
|------|----------|-------------|
| `CID-SMILES.gz` | `/Compound/Extras/` | All compound SMILES (~5GB) |
| `CID-InChI-Key.gz` | `/Compound/Extras/` | InChI keys for all compounds |
| `CID-Synonym-filtered.gz` | `/Compound/Extras/` | Compound synonyms/names |
| `CID-Title.gz` | `/Compound/Extras/` | Preferred compound names |
| SDF files | `/Compound/CURRENT-Full/SDF/` | Full structure data (many files) |

### Total Database Size
- **Full Database**: Multiple terabytes (TB)
- **Compound SDF (full)**: ~500GB+
- **Extras files**: Varies, most are 1-10GB compressed

### Licensing
- **License**: Public domain (no restrictions)
- **Registration**: Not required

### Filtering Strategies
```bash
# Download just SMILES for all compounds
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz

# Download specific compound range (SDF)
# Files are split by CID ranges (e.g., Compound_000000001_000500000.sdf.gz)
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_000000001_000500000.sdf.gz
```

### Recommended Processing
- **Start Small**: Begin with Extras files (SMILES, InChI) before full SDF
- **Incremental Updates**: Use Daily directory for updates after initial download
- **PUG-REST**: For targeted queries, use the PUG-REST API instead of bulk download
- **Storage**: Plan for TB-scale storage if downloading full database
- **Streaming**: Process SDF files in streaming mode due to size

---

## 4. PharmGKB

### Overview
PharmGKB (Pharmacogenomics Knowledge Base) is a comprehensive resource for pharmacogenomic data, including gene-drug associations, clinical annotations, and drug dosing guidelines.

### Download URLs
- **Main Site**: https://www.pharmgkb.org/
- **Downloads**: https://www.pharmgkb.org/downloads
- **Clinical Annotations Help**: https://www.pharmgkb.org/page/downloadClinicalAnnotationsHelp

### Available Files
| File | Description | Format |
|------|-------------|--------|
| `relationships.zip` | Gene-drug-disease relationships | TSV |
| `clinical_annotations.zip` | Variant-drug clinical annotations | TSV |
| `var_drug_ann.zip` | Variant annotations for drugs | TSV |
| `var_fa_ann.zip` | Variant functional annotations | TSV |
| `var_pheno_ann.zip` | Variant phenotype annotations | TSV |
| `pathways-tsv.zip` | Pathway data | TSV |
| `drugs.zip` | Drug information | TSV |
| `genes.zip` | Gene information | TSV |
| `automated_annotations.zip` | Automated literature annotations | TSV |
| `guideline_annotations.zip` | CPIC guideline annotations | TSV |
| `pathway diagrams` | Visual pathway data | BioPAX, GPML |

### File Sizes
- Individual files typically range from 1-50MB
- Complete download package: ~200-500MB

### Licensing
- **License**: Creative Commons Attribution-ShareAlike 4.0 (CC BY-SA 4.0)
- **Registration**: Required (free)
- **Terms**: Research use only, no redistribution without permission

### Account Registration
1. Visit https://www.pharmgkb.org/
2. Create free account
3. Agree to license terms
4. Access Downloads section

### Recommended Processing
- **Format**: TSV files load easily into pandas, R, or spreadsheets
- **Integration**: Cross-reference with DrugBank using drug identifiers
- **Updates**: Check monthly for new annotations
- **CPIC Guidelines**: Focus on guideline_annotations for clinical implementation

---

## 5. DGIdb

### Overview
DGIdb (Drug Gene Interaction Database) aggregates drug-gene interaction data from multiple sources. DGIdb 5.0 is the current version, containing comprehensive drug-gene interactions for precision medicine applications.

### Download URLs
- **Downloads Page**: https://dgidb.org/downloads
- **Legacy Downloads**: https://dgidb.genome.wustl.edu/downloads
- **API Documentation**: https://dgidb.org/api

### Available Files
| File | Description | Format |
|------|-------------|--------|
| `interactions.tsv` | All drug-gene interaction claims | TSV |
| `genes.tsv` | All gene claims | TSV |
| `drugs.tsv` | All drug claims | TSV |
| `categories.tsv` | Genes in druggable categories | TSV |
| SQL database dump | Complete database | PostgreSQL |

### File Sizes
- TSV files: 10-100MB each
- SQL dump: ~500MB

### Licensing
- **License**: Open access
- **Registration**: Not required for downloads

### Command Line Download
```bash
# Download interactions
wget https://dgidb.org/data/monthly_tsvs/interactions.tsv

# Download genes
wget https://dgidb.org/data/monthly_tsvs/genes.tsv

# Download drugs
wget https://dgidb.org/data/monthly_tsvs/drugs.tsv
```

### Alternative Access
```python
# GraphQL API for custom queries
import requests

query = """
{
  genes(names: ["BRAF"]) {
    nodes {
      name
      interactions {
        nodes {
          drug {
            name
          }
          interactionScore
        }
      }
    }
  }
}
"""

response = requests.post('https://dgidb.org/api/graphql', json={'query': query})
```

### Recommended Processing
- **TSV Files**: Easiest for quick analysis
- **GraphQL API**: Best for targeted queries
- **SQL Dump**: Best for comprehensive local analysis
- **Integration**: Map to Entrez gene IDs for cross-database linking

---

## 6. Open Targets

### Overview
Open Targets Platform integrates public data to enable systematic drug target identification and prioritization. Data is released quarterly in Parquet format.

### Download URLs
- **Downloads Page**: https://platform.opentargets.org/downloads
- **Documentation**: https://platform-docs.opentargets.org/data-access/datasets
- **FTP**: https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/
- **Latest Release**: https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/

### Available Datasets (Release 25.03+)
| Dataset | Description | Approximate Size |
|---------|-------------|------------------|
| `targets/` | Target (gene/protein) information | ~500MB |
| `diseases/` | Disease/phenotype data | ~100MB |
| `drugs/` | Drug/molecule data | ~200MB |
| `evidence/` | Evidence linking targets to diseases | ~5GB |
| `association_by_overall_direct/` | Direct associations | ~1GB |
| `association_by_overall_indirect/` | Indirect associations | ~2GB |
| `association_by_datasource_direct/` | Per-source associations | ~3GB |
| `interaction/` | Molecular interactions | ~500MB |
| `mouse_phenotypes/` | Mouse phenotype data | ~200MB |

### Total Size
- **Core Platform Data**: ~10-15GB
- **With Genetics Data**: ~350GB (Azure mirror)

### Licensing
- **License**: Open access (mostly CC0 or CC BY)
- **Registration**: Not required

### Download Methods
```bash
# Using wget (recommended)
wget -r -np -nH --cut-dirs=4 \
  https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.03/output/targets/

# Using lftp (for large downloads)
lftp -c "mirror https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.03/output/"

# Using rsync
rsync -avz rsync://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.03/ ./opentargets/
```

### Google Cloud Access
```bash
# Using gsutil
gsutil -m cp -r gs://open-targets-data-releases/25.03/output/ ./
```

### Reading Parquet Files
```python
import pandas as pd

# Read partitioned Parquet dataset
targets = pd.read_parquet('targets/')

# Or use PyArrow for large datasets
import pyarrow.parquet as pq
dataset = pq.ParquetDataset('evidence/')
table = dataset.read()
```

### Recommended Processing
- **Format**: Parquet files - use pandas, PyArrow, or Spark
- **Tool**: Use `p2j` for Parquet to JSON conversion if needed
- **Storage**: Plan for 15-50GB depending on datasets needed
- **Updates**: Quarterly releases (March, June, September, December)

---

## 7. BindingDB

### Overview
BindingDB is a public database of measured binding affinities between proteins and drug-like molecules. As of 2024, it contains 3,173,561 binding measurements for 1,390,080 compounds and 11,382 targets.

### Download URLs
- **Main Downloads**: https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp
- **SDF Downloads**: http://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFdownload.jsp

### Available Files - Full Database
| File | Description | Size |
|------|-------------|------|
| `BindingDB_All_202601_tsv.zip` | All data in TSV format | 545.23 MB |
| `BindingDB_All_2D_202601_sdf.zip` | All 2D structures (SDF) | 1.45 GB |
| `BindingDB_All_3D_202601_sdf.zip` | All 3D structures (SDF) | 831.86 MB |
| `BDB-mySQL_All_202508_dmp.zip` | MySQL database dump | 267.80 MB |

### Subset Downloads
| Source/Type | TSV | 2D SDF | 3D SDF |
|-------------|-----|--------|--------|
| BindingDB Curated | 18.03 MB | 35.46 MB | 86.03 MB |
| ChEMBL Source | 310.22 MB | 885.20 MB | 1.60 GB |
| Patent Data | 149.65 MB | 437.24 MB | 0.98 GB |
| PubChem Source | 25.40 MB | 62.53 MB | 107.41 MB |
| COVID-19 | 9.56 MB | 19.09 MB | 28.90 MB |
| PDSP Ki Data | 5.32 MB | 15.40 MB | 31.70 MB |

### Reference Files
| File | Description | Size |
|------|-------------|------|
| `BindingDBTargetSequences.fasta` | Target protein sequences | 5.96 MB |
| `BindingDB_CID.txt` | PubChem CID mappings | 21.33 MB |
| `BindingDB_SID.txt` | PubChem SID mappings | 22.11 MB |
| `BindingDB_UniProt_mapping.txt` | UniProt ID mappings | 473.80 KB |

### Licensing
- **License**: Creative Commons Attribution-ShareAlike 3.0 (CC BY-SA 3.0)
- **Registration**: Not required
- **Citation**: Required for publications

### Command Line Download
```bash
# Download full TSV
wget https://www.bindingdb.org/bind/downloads/BindingDB_All_202601_tsv.zip

# Download full 2D SDF
wget https://www.bindingdb.org/bind/downloads/BindingDB_All_2D_202601_sdf.zip
```

### Recommended Processing
- **TSV Format**: Easiest to work with - one row per binding measurement
- **SMILES**: Included in TSV for structure analysis
- **Spreadsheets**: TSV files work with Excel/LibreOffice (watch for size limits)
- **Updates**: Monthly releases - check download page for latest version

---

## 8. OpenFDA

### Overview
OpenFDA provides open access to FDA data including drug adverse events (FAERS), drug labeling, recalls, and more. Data is available via API and bulk downloads.

### Download URLs
- **Downloads Page**: https://open.fda.gov/apis/downloads/
- **Download Index (JSON)**: https://api.fda.gov/download.json
- **File Hosting**: https://download.open.fda.gov/

### Available Endpoints
| Endpoint | Description | Data Volume |
|----------|-------------|-------------|
| `/drug/event/` | FAERS adverse events | ~130GB total (~1000 files) |
| `/drug/label/` | Drug labeling/SPL | ~10GB |
| `/drug/ndc/` | NDC directory | ~100MB |
| `/drug/enforcement/` | Drug recalls | ~500MB |
| `/drug/drugsfda/` | Drugs@FDA data | ~100MB |
| `/device/event/` | Device adverse events | ~20GB |
| `/food/event/` | Food adverse events | ~5GB |

### FAERS Data Details
- **Format**: Zipped JSON files
- **Structure**: Split into chunks up to 150MB each
- **Per Quarter**: Up to 2GB per quarter
- **Total Files**: ~1000 files (as of 2024)
- **Total Size**: ~130GB

### Licensing
- **License**: Public domain
- **Registration**: Optional (API key increases rate limits)
- **API Limits**: 40 requests/minute without key, 240/minute with key

### Download Methods
```bash
# Get download manifest
curl -s https://api.fda.gov/download.json | jq '.results.drug.event.partitions[].file'

# Download all FAERS files
curl -s https://api.fda.gov/download.json | \
  jq -r '.results.drug.event.partitions[].file' | \
  xargs -I {} wget {}

# Download specific endpoint
wget https://download.open.fda.gov/drug/label/drug-label-0001-of-0001.json.zip
```

### Important Notes
- **Full Re-download Required**: Updates may affect all records - cannot download only "new" files
- **Quarterly Updates**: FAERS data updated quarterly
- **JSON Format**: Same structure as API responses

### Recommended Processing
- **Storage**: Plan for 150GB+ for complete FAERS data
- **Incremental**: Not supported - must re-download everything for updates
- **Parsing**: Use streaming JSON parsers for large files
- **Integration**: Cross-reference with RxNorm for drug normalization

---

## 9. DailyMed

### Overview
DailyMed provides FDA-approved drug labeling (Structured Product Labeling - SPL) for prescription and OTC drugs, as well as other product types.

### Download URLs
- **All Drug Labels**: https://dailymed.nlm.nih.gov/dailymed/spl-resources-all-drug-labels.cfm
- **SPL Resources**: https://dailymed.nlm.nih.gov/dailymed/spl-resources.cfm
- **Mapping Files**: https://dailymed.nlm.nih.gov/dailymed/spl-resources-all-mapping-files.cfm
- **Indexing Files**: https://dailymed.nlm.nih.gov/dailymed/spl-resources-all-indexing-files.cfm

### Full Release Files
| Category | Parts | Total Files | Size Per Part |
|----------|-------|-------------|---------------|
| Human Prescription | 5 parts | ~60,000 | 2.93-3.00 GB each |
| Human OTC | 11 parts | ~73,655 | 0.37-3.00 GB each |
| Homeopathic | 1 file | ~15,529 | 5.29 GB |
| Animal | 1 file | ~3,241 | 1.02 GB |
| Remainder | 1 file | ~2,367 | 571 MB |

### Periodic Updates
| Type | Frequency | Typical Size |
|------|-----------|--------------|
| Daily | Every day | 100-500 MB |
| Weekly | Every week | 500 MB - 1 GB |
| Monthly | Every month | 2-4 GB |

### Recent Monthly Updates (Example)
- December 2025: 8,200 labels, 2.96 GB
- November 2025: 9,479 labels, 3.53 GB

### Licensing
- **License**: Public domain (US Government work)
- **Registration**: Not required

### Download Methods
```bash
# Download full Human Rx (all parts)
wget https://dailymed.nlm.nih.gov/dailymed/downloads/dm_spl_release_human_rx_part1.zip
wget https://dailymed.nlm.nih.gov/dailymed/downloads/dm_spl_release_human_rx_part2.zip
# ... continue for all parts

# Download by SET ID (API)
wget "https://dailymed.nlm.nih.gov/dailymed/downloadzipfile.cfm?setId={SET_ID}"

# Download specific version
wget "https://dailymed.nlm.nih.gov/dailymed/getFile.cfm?type=zip&setid={SET_ID}&version={VERSION}"
```

### Web Services
```bash
# Get SPL history
curl "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{SET_ID}/history.json"

# Search SPLs
curl "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name=aspirin"
```

### Recommended Processing
- **Initial Download**: Get full release files (split into parts due to size)
- **Updates**: Use daily/weekly incremental downloads after initial load
- **Format**: SPL files are XML - use XML parsers
- **Mapping**: Download RxNorm mapping files for cross-referencing

---

## 10. RxNorm

### Overview
RxNorm provides normalized names and identifiers for clinical drugs. It is part of the UMLS (Unified Medical Language System) and serves as a drug vocabulary standard.

### Download URLs
- **Main Page**: https://www.nlm.nih.gov/research/umls/rxnorm/index.html
- **Files Documentation**: https://www.nlm.nih.gov/research/umls/rxnorm/docs/rxnormfiles.html
- **Automation Guide**: https://documentation.uts.nlm.nih.gov/automating-downloads.html

### Generic URLs (Always Current)
```
https://download.nlm.nih.gov/umls/kss/rxnorm/RxNorm_full_current.zip
https://download.nlm.nih.gov/umls/kss/rxnorm/RxNorm_weekly_current.zip
```

### Available Files
| File | Description | Release Frequency |
|------|-------------|-------------------|
| `RxNorm_full_MMDDYYYY.zip` | Complete RxNorm release | Monthly (1st Monday) |
| `RxNorm_weekly_MMDDYYYY.zip` | Weekly incremental update | Weekly (Wed/Thu) |

### File Contents
RxNorm files are pipe-delimited Rich Release Format (RRF) files:
- `RXNCONSO.RRF` - Concept names and sources
- `RXNREL.RRF` - Relationships between concepts
- `RXNSAT.RRF` - Simple concept attributes
- `RXNSTY.RRF` - Semantic types
- `RXNDOC.RRF` - Documentation
- `RXNATOMARCHIVE.RRF` - Archive of retired atoms

### File Sizes
- **Full Release**: ~300-400 MB compressed
- **Weekly Update**: ~50-100 MB compressed

### Licensing
- **License**: UMLS Metathesaurus License Agreement (free)
- **Registration**: Required - UMLS Terminology Services (UTS) account
- **Proprietary Sources**: Some included data requires additional licensing

### Account Registration
1. Visit https://uts.nlm.nih.gov/
2. Create UTS account (free)
3. Accept UMLS License Agreement
4. Generate API key for automated downloads

### Automated Download
```bash
# Using curl with API key
curl -L -o RxNorm_full_current.zip \
  "https://download.nlm.nih.gov/umls/kss/rxnorm/RxNorm_full_current.zip?apiKey=YOUR_API_KEY"

# Using the UTS API
curl "https://uts-ws.nlm.nih.gov/releases" \
  -H "Authorization: Bearer YOUR_API_KEY"
```

### Python Package
```bash
pip install umls_downloader
```

```python
from umls_downloader import download_rxnorm

# Automatically handles authentication and download
path = download_rxnorm()
```

### Recommended Processing
- **Format**: RRF files load into relational databases
- **MySQL/PostgreSQL**: Load scripts available from NLM
- **Updates**: Use weekly releases for incremental updates
- **Integration**: Primary vocabulary for drug name normalization

---

## Summary Comparison

| Database | Size | License | Account Required | Update Frequency |
|----------|------|---------|------------------|------------------|
| DrugBank | ~500MB | Academic free / Commercial paid | Yes | Periodic |
| ChEMBL | ~5GB (SQLite) | CC BY-SA 3.0 | No | Quarterly |
| PubChem | TBs | Public Domain | No | Daily |
| PharmGKB | ~500MB | CC BY-SA 4.0 | Yes | Monthly |
| DGIdb | ~100MB | Open Access | No | Monthly |
| Open Targets | ~10-15GB | CC0/CC BY | No | Quarterly |
| BindingDB | ~1.5GB (TSV) | CC BY-SA 3.0 | No | Monthly |
| OpenFDA | ~130GB (FAERS) | Public Domain | Optional | Quarterly |
| DailyMed | ~30GB | Public Domain | No | Daily |
| RxNorm | ~400MB | UMLS License | Yes | Monthly/Weekly |

---

## Best Practices

### Storage Planning
```
Minimal Setup (core data):        ~20GB
Standard Setup (main databases):  ~50GB
Comprehensive Setup:              ~200GB+
Full PubChem + OpenFDA:           ~500GB+
```

### Download Priority
1. **Start with**: DrugBank, ChEMBL, DGIdb (smaller, high value)
2. **Add next**: PharmGKB, BindingDB, RxNorm (specialized data)
3. **Large datasets**: Open Targets, DailyMed (comprehensive)
4. **Massive**: PubChem, OpenFDA (plan storage carefully)

### Update Strategy
- **Daily**: PubChem (if needed), DailyMed
- **Weekly**: RxNorm
- **Monthly**: DrugBank, PharmGKB, DGIdb, BindingDB
- **Quarterly**: ChEMBL, Open Targets, OpenFDA FAERS

### Data Integration Tips
1. Use RxNorm as the primary drug vocabulary linker
2. Map DrugBank IDs to ChEMBL IDs via InChI keys
3. Use PubChem CIDs as universal structure identifiers
4. Cross-reference targets using UniProt IDs

---

## References

- DrugBank: https://go.drugbank.com/
- ChEMBL: https://www.ebi.ac.uk/chembl/
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
- PharmGKB: https://www.pharmgkb.org/
- DGIdb: https://dgidb.org/
- Open Targets: https://platform.opentargets.org/
- BindingDB: https://www.bindingdb.org/
- OpenFDA: https://open.fda.gov/
- DailyMed: https://dailymed.nlm.nih.gov/
- RxNorm: https://www.nlm.nih.gov/research/umls/rxnorm/

---

*Last Updated: January 2026*
