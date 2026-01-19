# Alternative Data Sources

**Last Updated:** January 2026
**Purpose:** Alternative acquisition methods for research data when primary APIs/FTP are limited.

---

## Summary

| Category | Sources | Content Type | Access Method |
|----------|---------|--------------|---------------|
| Dataset Repositories | 3 | ML/genomics datasets | Direct download |
| Research Archives | 2 | Scientific papers | Torrents, DOI lookup |
| Cloud Public Data | 2 | Large genomics datasets | S3/GCS buckets |
| Archived Projects | 2 | Defunct database dumps | Zenodo, archives |

---

## 1. Hugging Face Datasets (Free, Direct Download)

Hugging Face hosts numerous genomics and biomedical datasets under open licenses.

### 1.1 ncbi/pubmed
| Field | Value |
|-------|-------|
| **URL** | https://huggingface.co/datasets/ncbi/pubmed |
| **Content** | Full PubMed baseline - 36M+ citations |
| **Size** | ~100 GB |
| **License** | Public domain |
| **Format** | Parquet |
| **Use Case** | RAG system for research retrieval |

### 1.2 InstaDeepAI/genomics-long-range-benchmark
| Field | Value |
|-------|-------|
| **URL** | https://huggingface.co/datasets/InstaDeepAI/genomics-long-range-benchmark |
| **Content** | SNP pathogenicity data from OMIM/gnomAD |
| **Size** | ~1 GB |
| **License** | Apache 2.0 |
| **Use Case** | Variant pathogenicity scoring |

### 1.3 microsoft/msr_genomics_kbcomp
| Field | Value |
|-------|-------|
| **URL** | https://huggingface.co/datasets/microsoft/msr_genomics_kbcomp |
| **Content** | Gene regulation + PubMed mentions |
| **Size** | ~500 MB |
| **License** | MIT |
| **Use Case** | Gene-research linking |

### 1.4 katielink/genomic-benchmarks
| Field | Value |
|-------|-------|
| **URL** | https://huggingface.co/datasets/katielink/genomic-benchmarks |
| **Content** | Genomic sequence classification |
| **Size** | ~2 GB |
| **License** | MIT |
| **Use Case** | ML model training |

### 1.5 InstaDeepAI/multi_species_genomes
| Field | Value |
|-------|-------|
| **URL** | https://huggingface.co/datasets/InstaDeepAI/multi_species_genomes |
| **Content** | 850 species, 174B nucleotides |
| **Size** | ~100 GB |
| **License** | Apache 2.0 |
| **Use Case** | Comparative genomics |

---

## 2. Research Paper Archives

### 2.1 Anna's Archive
| Field | Value |
|-------|-------|
| **URL** | https://annas-archive.org |
| **Content** | 61M books, 95M+ scientific papers |
| **Total Size** | ~1.1 PB (full archive) |
| **Access** | DOI lookup, torrent bulk downloads |
| **Metadata** | Available as MariaDB/ElasticSearch dumps |

**Access Methods:**
- **Single Paper:** DOI-based lookup via web interface
- **Bulk Download:** Torrent files for metadata and content
- **SciDB Subset:** Scientific papers at annas-archive.li/scidb

**Use Case:** Building local research corpus for RAG, accessing papers behind paywalls.

### 2.2 Sci-Hub
| Field | Value |
|-------|-------|
| **URL** | https://sci-hub.se (mirrors vary) |
| **Content** | 85M+ scientific papers |
| **Total Size** | ~77 TB |
| **Access** | DOI-based lookup |

**Access Methods:**
- **Single Paper:** Enter DOI to retrieve PDF
- **Bulk Download:** Torrent collections available (large)

**Legal Note:** Copyright status varies by jurisdiction. Use for personal research.

---

## 3. Academic Torrents

| Field | Value |
|-------|-------|
| **URL** | https://academictorrents.com |
| **Content** | Research datasets, course materials |
| **Categories** | Datasets, papers, courses |

**Relevant Collections:**
- Search "genomics" for available datasets
- Educational course materials on genetics
- Research dataset supplements

**Note:** BioTorrents (dedicated biological data torrents) is now defunct. Academic Torrents is the primary remaining academic torrent tracker.

---

## 4. Zenodo (Open Research Data)

| Field | Value |
|-------|-------|
| **URL** | https://zenodo.org |
| **Content** | Research data, supplementary files, datasets |
| **Access** | DOI-based download, REST API |
| **Tool** | `zenodo_get` CLI for bulk download |

**Relevant Searches:**
- "genomics datasets"
- "SNP data"
- "genetic variation"
- "pharmacogenomics"

**Key Datasets:**
- OpenSNP archived data
- Dosage sensitivity maps
- 3D genomics datasets
- Supplementary data from published papers

---

## 5. Cloud Public Datasets (Free Egress from Same Region)

### 5.1 AWS Open Data Registry
| Field | Value |
|-------|-------|
| **URL** | https://registry.opendata.aws |
| **Access** | S3 buckets (free from same region) |

**Key Genomics Datasets:**

| Dataset | S3 Path | Size |
|---------|---------|------|
| **gnomAD** | s3://gnomad-public-us-east-1/ | ~30 TB |
| **1000 Genomes** | s3://1000genomes/ | ~260 TB |
| **TCGA** | Via GDC portal | ~2.5 PB |
| **ClinVar** | s3://aws-roda-hcls-datalake/clinvar/ | ~500 MB |
| **dbSNP** | s3://aws-roda-hcls-datalake/dbsnp/ | ~15 GB |

**Cost Optimization:**
- Run processing in us-east-1 to avoid egress fees
- Use AWS Free Tier EC2 for initial exploration
- Download only summary/subset files for MVP

### 5.2 Google Cloud Public Datasets
| Field | Value |
|-------|-------|
| **URL** | https://cloud.google.com/public-datasets |
| **Access** | GCS buckets, BigQuery |

**Relevant Datasets:**
- gnomAD (BigQuery)
- 1000 Genomes
- ClinVar
- Various NCBI datasets

---

## 6. Archived/Defunct Database Dumps

### 6.1 OpenSNP Archive
| Field | Value |
|-------|-------|
| **Status** | Shut down April 2025 |
| **Reason** | Privacy concerns (authoritarian governments, 23andMe bankruptcy) |
| **Archive Location** | Zenodo, GenomePrep |
| **Content** | User-contributed genotype data |
| **Note** | Historical data only, no new contributions |

**Background:** Founder Bastian Greshake Tzovaras is now Director of Research at Open Humans Foundation. OpenSNP was replaced by privacy-preserving alternatives rather than direct successors.

### 6.2 GenomePrep
| Field | Value |
|-------|-------|
| **Content** | Archived OpenSNP data |
| **Access** | Check Zenodo for availability |

---

## 7. Data Acquisition Strategy

### Priority Order for Alternative Sources

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

4. **Anna's Archive / Sci-Hub** (For research papers)
   - Build local corpus for RAG
   - Use for papers not in PubMed Open Access

### Download Tools

```bash
# Hugging Face datasets
pip install datasets
python -c "from datasets import load_dataset; ds = load_dataset('ncbi/pubmed')"

# Zenodo
pip install zenodo_get
zenodo_get 10.5281/zenodo.XXXXXXX

# AWS S3 (no credentials needed for public data)
aws s3 cp s3://gnomad-public-us-east-1/release/4.1/vcf/joint/ ./gnomad/ --recursive --no-sign-request

# Academic Torrents
# Use any BitTorrent client with magnet links from academictorrents.com
```

---

## 8. Legal Considerations

| Source | Legal Status | Recommendation |
|--------|--------------|----------------|
| Hugging Face | Clear licenses | Safe to use per license terms |
| AWS/GCS Public | Open access | Safe to use |
| Zenodo | Per-dataset license | Check each dataset |
| Academic Torrents | Varies | Verify dataset license |
| Anna's Archive | Gray area | Personal research use |
| Sci-Hub | Copyright concerns | Personal research use |

**Recommendation:** Prioritize clearly licensed sources (Hugging Face, AWS, official FTP) for any commercial or public-facing use. Use Anna's Archive/Sci-Hub only for personal research and RAG training where papers are otherwise inaccessible.

---

## 9. Storage Requirements

| Source Category | Minimum | Recommended | Full |
|-----------------|---------|-------------|------|
| Hugging Face datasets | 5 GB | 50 GB | 200 GB |
| Research papers (subset) | 10 GB | 100 GB | 1 TB |
| Cloud genomics (summaries) | 50 GB | 200 GB | 500 GB |
| Archived databases | 5 GB | 20 GB | 50 GB |
| **Total Alternative Sources** | **70 GB** | **370 GB** | **1.75 TB** |

---

## References

- [Hugging Face Datasets](https://huggingface.co/datasets)
- [AWS Open Data Registry](https://registry.opendata.aws)
- [Academic Torrents](https://academictorrents.com)
- [Zenodo](https://zenodo.org)
- [OpenSNP Shutdown Announcement](https://opensnp.org/shutdown)
