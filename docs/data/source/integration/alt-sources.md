# Alternative Data Sources

**Document ID:** 43-85-ALT-DATA-SOURCES
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [index.md](../index.md)

---

## TL;DR

Alternative data sources provide backup acquisition methods when primary APIs/FTP are limited. **Hugging Face datasets** (clearly licensed, easy download) are the safest option with 5 genomics-relevant datasets including PubMed (36M+ citations, public domain). **AWS/GCS public data** provides large genomics datasets (gnomAD 30TB, 1000 Genomes 260TB) with free egress from same region. **Zenodo** offers research supplements and archived databases like OpenSNP. **Shadow libraries (Anna's Archive, Sci-Hub) should NOT be used** for commercial applications due to extreme legal risk ($1B+ lawsuits against AI companies). Total alternative source storage: 70GB minimum to 1.75TB full.

---

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **Primary Alternative** | Hugging Face | Clear licenses, preprocessed, easy download |
| **Large Genomics Data** | AWS/GCS Public | Free egress from same region, summary files |
| **Research Supplements** | Zenodo | DOI-based, archived databases |
| **Shadow Libraries** | DO NOT USE | Extreme legal risk for commercial use |
| **Storage Strategy** | Summary files first | Download subsets, process in cloud |
| **Torrent Sources** | Academic Torrents only | BioTorrents defunct |

---

## Database Catalog

### 1. Hugging Face Datasets

#### 1.1 Overview

| Attribute | Value |
|-----------|-------|
| **Provider** | Hugging Face Inc. |
| **URL** | https://huggingface.co/datasets |
| **License** | Per-dataset (mostly open) |
| **Commercial Use** | Per-dataset license |
| **Access Method** | `datasets` Python library |

**Installation:**
```bash
pip install datasets
```

#### 1.2 ncbi/pubmed

| Attribute | Value |
|-----------|-------|
| **URL** | https://huggingface.co/datasets/ncbi/pubmed |
| **Content** | Full PubMed baseline - 36M+ citations |
| **Size** | ~100 GB |
| **License** | Public domain |
| **Format** | Parquet |
| **Commercial Use** | YES |

**Use Case:** RAG system for research retrieval, literature corpus

**Download:**
```python
from datasets import load_dataset
ds = load_dataset('ncbi/pubmed')
```

---

#### 1.3 InstaDeepAI/genomics-long-range-benchmark

| Attribute | Value |
|-----------|-------|
| **URL** | https://huggingface.co/datasets/InstaDeepAI/genomics-long-range-benchmark |
| **Content** | SNP pathogenicity data from OMIM/gnomAD |
| **Size** | ~1 GB |
| **License** | Apache 2.0 |
| **Format** | Parquet |
| **Commercial Use** | YES |

**Use Case:** Variant pathogenicity scoring, ML model training

---

#### 1.4 microsoft/msr_genomics_kbcomp

| Attribute | Value |
|-----------|-------|
| **URL** | https://huggingface.co/datasets/microsoft/msr_genomics_kbcomp |
| **Content** | Gene regulation + PubMed mentions |
| **Size** | ~500 MB |
| **License** | MIT |
| **Format** | Parquet |
| **Commercial Use** | YES |

**Use Case:** Gene-research linking, knowledge base completion

---

#### 1.5 katielink/genomic-benchmarks

| Attribute | Value |
|-----------|-------|
| **URL** | https://huggingface.co/datasets/katielink/genomic-benchmarks |
| **Content** | Genomic sequence classification |
| **Size** | ~2 GB |
| **License** | MIT |
| **Format** | Parquet |
| **Commercial Use** | YES |

**Use Case:** ML model training, sequence analysis

---

#### 1.6 InstaDeepAI/multi_species_genomes

| Attribute | Value |
|-----------|-------|
| **URL** | https://huggingface.co/datasets/InstaDeepAI/multi_species_genomes |
| **Content** | 850 species, 174B nucleotides |
| **Size** | ~100 GB |
| **License** | Apache 2.0 |
| **Format** | Parquet |
| **Commercial Use** | YES |

**Use Case:** Comparative genomics, cross-species analysis

---

### 2. Cloud Public Datasets

#### 2.1 AWS Open Data Registry

| Attribute | Value |
|-----------|-------|
| **Provider** | Amazon Web Services |
| **URL** | https://registry.opendata.aws |
| **Access** | S3 buckets (free from same region) |
| **Authentication** | None required for public data |
| **Cost Optimization** | Run in us-east-1 for zero egress |

**Key Genomics Datasets:**

| Dataset | S3 Path | Size | Description |
|---------|---------|------|-------------|
| gnomAD | s3://gnomad-public-us-east-1/ | ~30 TB | Population variant frequencies |
| 1000 Genomes | s3://1000genomes/ | ~260 TB | Population sequencing data |
| TCGA | Via GDC portal | ~2.5 PB | Cancer genomics |
| ClinVar | s3://aws-roda-hcls-datalake/clinvar/ | ~500 MB | Clinical variants |
| dbSNP | s3://aws-roda-hcls-datalake/dbsnp/ | ~15 GB | SNP database |

**Download Example:**
```bash
# No credentials needed for public data
aws s3 cp s3://gnomad-public-us-east-1/release/4.1/vcf/joint/ ./gnomad/ --recursive --no-sign-request
```

---

#### 2.2 Google Cloud Public Datasets

| Attribute | Value |
|-----------|-------|
| **Provider** | Google Cloud Platform |
| **URL** | https://cloud.google.com/public-datasets |
| **Access** | GCS buckets, BigQuery |
| **Cost Optimization** | Query via BigQuery, download summaries |

**Relevant Datasets:**

| Dataset | Access | Description |
|---------|--------|-------------|
| gnomAD | BigQuery | Query-based access |
| 1000 Genomes | GCS | Population data |
| ClinVar | GCS | Clinical variants |
| NCBI datasets | GCS | Various NCBI data |

---

### 3. Research Archives

#### 3.1 Zenodo

| Attribute | Value |
|-----------|-------|
| **Provider** | CERN / European Commission |
| **URL** | https://zenodo.org |
| **Content** | Research data, supplementary files, datasets |
| **Access** | DOI-based download, REST API |
| **License** | Per-dataset |
| **Commercial Use** | Per-dataset license |

**CLI Tool:**
```bash
pip install zenodo_get
zenodo_get 10.5281/zenodo.XXXXXXX
```

**Relevant Searches:**
- "genomics datasets"
- "SNP data"
- "genetic variation"
- "pharmacogenomics"

**Key Archived Data:**
- OpenSNP archived data
- Dosage sensitivity maps
- 3D genomics datasets
- Paper supplementary data

---

#### 3.2 Academic Torrents

| Attribute | Value |
|-----------|-------|
| **Provider** | Academic Torrents |
| **URL** | https://academictorrents.com |
| **Content** | Research datasets, course materials |
| **Access** | BitTorrent protocol |
| **License** | Varies by dataset |

**Relevant Collections:**
- Search "genomics" for available datasets
- Educational genetics courses
- Research dataset supplements

**Note:** BioTorrents (dedicated biological data) is now defunct. Academic Torrents is the primary remaining academic torrent tracker.

---

### 4. Archived/Defunct Database Dumps

#### 4.1 OpenSNP Archive

| Attribute | Value |
|-----------|-------|
| **Status** | Shut down April 2025 |
| **Reason** | Privacy concerns (authoritarian governments, 23andMe bankruptcy) |
| **Archive Location** | Zenodo, GenomePrep |
| **Content** | User-contributed genotype data |
| **Note** | Historical data only, no new contributions |

**Background:** Founder Bastian Greshake Tzovaras is now Director of Research at Open Humans Foundation. OpenSNP was replaced by privacy-preserving alternatives.

---

#### 4.2 GenomePrep

| Attribute | Value |
|-----------|-------|
| **Content** | Archived OpenSNP data |
| **Access** | Check Zenodo for availability |
| **Status** | Archive only |

---

### 5. Shadow Libraries (DO NOT USE)

#### 5.1 Legal Risk Assessment

| Source | Status | Legal Risk | Recommendation |
|--------|--------|------------|----------------|
| Anna's Archive | Active | EXTREME | AVOID |
| Sci-Hub | Active | EXTREME | AVOID |
| Library Genesis | Active | EXTREME | AVOID |
| Z-Library | Active | EXTREME | AVOID |

#### 5.2 Legal Precedents (2024-2026)

| Case | Year | Outcome |
|------|------|---------|
| Bartz v. Anthropic | 2025 | Judge ruled LibGen downloads = willful infringement. Settlement: $1.5B |
| Meta Litigation | 2025 | Unsealed emails revealed 81+ TB from Anna's Archive torrents |
| OpenAI | 2024-2025 | Court allowed plaintiffs to pursue shadow library claim |

**Bottom Line:** For commercial applications, using shadow library content carries **existential legal risk**. All recommended sources above provide legal alternatives with sufficient coverage.

---

## Integration Priority

| Source | Priority | License | Size | Key Use Case |
|--------|----------|---------|------|--------------|
| HF: ncbi/pubmed | HIGH | Public domain | 100 GB | Literature RAG |
| HF: genomics-benchmark | HIGH | Apache 2.0 | 1 GB | Pathogenicity |
| AWS: ClinVar | HIGH | Open | 500 MB | Clinical variants |
| AWS: dbSNP | HIGH | Open | 15 GB | SNP reference |
| Zenodo | MEDIUM | Varies | Varies | Supplements |
| HF: msr_genomics | MEDIUM | MIT | 500 MB | Gene linking |
| AWS: gnomAD (summary) | MEDIUM | Open | ~5 GB | Frequencies |
| Academic Torrents | LOW | Varies | Varies | Educational |
| Shadow Libraries | NEVER | Illegal | N/A | DO NOT USE |

---

## Data Acquisition Strategy

### Priority Order

1. **Hugging Face** (First choice)
   - Clean, preprocessed datasets
   - Clear licensing
   - Easy download via `datasets` library

2. **AWS/GCS Public Data** (For large datasets)
   - gnomAD, 1000 Genomes
   - Process in cloud to avoid download
   - Use summary files when possible

3. **Zenodo** (For supplementary data)
   - Research dataset supplements
   - Archived databases

4. **Academic Torrents** (For educational content)
   - Course materials
   - Historical datasets

---

## Download Tools

### Hugging Face

```bash
pip install datasets

# Python usage
from datasets import load_dataset
ds = load_dataset('ncbi/pubmed')
```

### Zenodo

```bash
pip install zenodo_get
zenodo_get 10.5281/zenodo.XXXXXXX
```

### AWS S3

```bash
# No credentials needed for public data
aws s3 cp s3://aws-roda-hcls-datalake/clinvar/ ./clinvar/ --recursive --no-sign-request

# List bucket contents
aws s3 ls s3://gnomad-public-us-east-1/ --no-sign-request
```

### Academic Torrents

```bash
# Use any BitTorrent client with magnet links from academictorrents.com
# Recommended: qBittorrent, Transmission
```

---

## Storage Requirements

### Minimum Configuration (70 GB)

| Source Category | Size | Description |
|-----------------|------|-------------|
| Hugging Face datasets | 5 GB | Core genomics datasets |
| Research papers (subset) | 10 GB | Selected literature |
| Cloud genomics (summaries) | 50 GB | ClinVar, dbSNP, gnomAD summary |
| Archived databases | 5 GB | OpenSNP, supplements |

### Recommended Configuration (370 GB)

| Source Category | Size | Description |
|-----------------|------|-------------|
| Hugging Face datasets | 50 GB | Full genomics suite |
| Research papers | 100 GB | Comprehensive literature |
| Cloud genomics | 200 GB | Extended variant data |
| Archived databases | 20 GB | Full archives |

### Full Configuration (1.75 TB)

| Source Category | Size | Description |
|-----------------|------|-------------|
| Hugging Face datasets | 200 GB | All relevant datasets |
| Research papers | 1 TB | Complete literature corpus |
| Cloud genomics | 500 GB | Full frequency data |
| Archived databases | 50 GB | Complete archives |

---

## Legal Considerations

| Source | Legal Status | Commercial Use | Recommendation |
|--------|--------------|----------------|----------------|
| Hugging Face | Clear licenses | Per-dataset | Safe to use |
| AWS/GCS Public | Open access | YES | Safe to use |
| Zenodo | Per-dataset license | Check each | Verify license |
| Academic Torrents | Varies | Verify each | Check license |
| Anna's Archive | Gray area | NO | Personal research only |
| Sci-Hub | Copyright concerns | NO | Personal research only |

**Recommendation:** Prioritize clearly licensed sources (Hugging Face, AWS, official FTP) for any commercial or public-facing use.

---

## Dependencies

### Upstream Dependencies

| Dependency | Purpose | Risk if Unavailable |
|------------|---------|---------------------|
| Hugging Face Hub | Dataset downloads | LOW - direct URLs available |
| AWS S3 | Large genomics data | LOW - GCS alternative |
| Zenodo | Research archives | LOW - not critical |
| BitTorrent | Academic torrents | LOW - direct URLs |

### Downstream Dependents

| Dependent | Usage |
|-----------|-------|
| Literature RAG | PubMed corpus for retrieval |
| Variant Analysis | gnomAD, ClinVar, dbSNP |
| ML Training | Genomics benchmark datasets |
| Research Search | Paper corpus |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Initial document from research.old/data-sources-alt.md |
