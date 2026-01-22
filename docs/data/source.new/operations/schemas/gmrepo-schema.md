---
id: schema-gmrepo
title: "GMrepo REST API Schema"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./_index.md)

# GMrepo REST API Schema

**Document ID:** GMREPO-SCHEMA
**Source:** GMrepo - Human Gut Microbiome Repository (https://gmrepo.humangut.info)
**Last Updated:** January 2026
**Version:** v3 (current)

---

## Database Statistics (January 2026)

### GMrepo v3 (Current)

| Metric | Count |
|--------|-------|
| Projects | 890 |
| Runs/Samples | 118,965 |
| 16S rRNA Samples | 87,048 |
| Metagenomic (WGS) Samples | 31,917 |
| Annotated Diseases | 302 |
| Marker Taxa | 1,299 |

### Historical Versions

| Version | Projects | Samples | Diseases | Year |
|---------|----------|---------|----------|------|
| v1 | ~300 | ~71,000 | ~100 | 2020 |
| v2 | 353 | 71,642 | 133 | 2022 |
| v3 | 890 | 118,965 | 302 | 2025 |

---

## Database Information

| Attribute | Value |
|-----------|-------|
| URL | https://gmrepo.humangut.info |
| v2 Archive | https://gmrepo2022.humangut.info |
| GitHub | https://github.com/evolgeniusteam/GMrepoProgrammableAccess |
| Documentation | https://evolgeniusteam.github.io/gmrepodocumentation/ |
| Institution | Huazhong University of Science and Technology, Wuhan, China |
| Contact | Wei-Hua Chen (weihuachen@hust.edu.cn) |
| License | CC BY-NC 3.0 |
| Founded | 2020 |

---

## REST API Overview

### Base URL
```
https://gmrepo.humangut.info/api/
```

### Key API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/getRunsByProjectID` | GET | Fetch runs/samples by project ID |
| `/getRunsByPhenotype` | GET | Fetch samples by disease phenotype |
| `/getSpeciesAbundance` | GET | Get species relative abundance data |
| `/getMarkerTaxa` | GET | Retrieve disease marker taxa |
| `/getProjectList` | GET | List all available projects |
| `/getPhenotypeList` | GET | List all annotated phenotypes |
| `/getTaxonInfo` | GET | Get taxonomic information |
| `/getStatistics` | GET | Database statistics |

---

## Data Schema

### Project Schema

```json
{
  "project_id": "PRJEB6070",
  "project_name": "Metagenome of healthy individuals",
  "description": "Human gut microbiome study...",
  "sequencing_type": "16S",
  "sample_count": 156,
  "phenotypes": ["Health"],
  "publication_pmid": "25567908",
  "submission_date": "2014-06-15"
}
```

### Field Specifications

| Field | Type | Description |
|-------|------|-------------|
| `project_id` | string | Project accession (PRJNA/PRJEB format) |
| `project_name` | string | Project title |
| `description` | string | Project description |
| `sequencing_type` | string | "16S" or "WGS" |
| `sample_count` | integer | Number of samples |
| `phenotypes` | array | Associated phenotype list |
| `publication_pmid` | string | PubMed ID if published |
| `submission_date` | string | Submission date |

---

### Run/Sample Schema

```json
{
  "run_id": "SRR1234567",
  "project_id": "PRJNA12345",
  "sample_id": "SAMN1234567",
  "phenotype": "Type 2 diabetes",
  "phenotype_mesh_id": "D003924",
  "host_age": 45,
  "host_sex": "male",
  "host_bmi": 28.5,
  "country": "China",
  "sequencing_type": "16S",
  "antibiotics_usage": "no",
  "collection_date": "2018-06-15"
}
```

### Field Specifications

| Field | Type | Description |
|-------|------|-------------|
| `run_id` | string | Sequencing run accession |
| `project_id` | string | Parent project ID |
| `sample_id` | string | Sample accession |
| `phenotype` | string | Disease/health phenotype |
| `phenotype_mesh_id` | string | MeSH ID for phenotype |
| `host_age` | integer/null | Subject age |
| `host_sex` | string | "male", "female", or null |
| `host_bmi` | float/null | Body mass index |
| `country` | string | Country of sample origin |
| `sequencing_type` | string | Sequencing method |
| `antibiotics_usage` | string | Recent antibiotic use |
| `collection_date` | string | Sample collection date |

---

### Species Abundance Schema

```json
{
  "run_id": "SRR1234567",
  "species": [
    {
      "taxon_id": 816,
      "taxon_name": "Bacteroides",
      "rank": "genus",
      "relative_abundance": 0.2534,
      "read_count": 15234
    },
    {
      "taxon_id": 1263,
      "taxon_name": "Ruminococcus",
      "rank": "genus",
      "relative_abundance": 0.1234,
      "read_count": 7421
    }
  ]
}
```

### Abundance Field Specifications

| Field | Type | Description |
|-------|------|-------------|
| `taxon_id` | integer | NCBI Taxonomy ID |
| `taxon_name` | string | Taxonomic name |
| `rank` | string | Taxonomic rank (genus, species, etc.) |
| `relative_abundance` | float | Relative abundance (0-1) |
| `read_count` | integer | Number of sequencing reads |

---

### Marker Taxa Schema

```json
{
  "phenotype": "Type 2 diabetes",
  "phenotype_mesh_id": "D003924",
  "markers": [
    {
      "taxon_id": 816,
      "taxon_name": "Bacteroides",
      "rank": "genus",
      "direction": "increased",
      "effect_size": 0.45,
      "p_value": 0.001,
      "study_count": 12,
      "sample_count": 1500
    }
  ]
}
```

### Marker Field Specifications

| Field | Type | Description |
|-------|------|-------------|
| `phenotype` | string | Disease phenotype |
| `phenotype_mesh_id` | string | MeSH identifier |
| `taxon_id` | integer | NCBI Taxonomy ID |
| `taxon_name` | string | Taxonomic name |
| `direction` | string | "increased" or "decreased" |
| `effect_size` | float | Effect magnitude |
| `p_value` | float | Statistical significance |
| `study_count` | integer | Number of supporting studies |
| `sample_count` | integer | Total samples analyzed |

---

### Phenotype/Disease Schema

```json
{
  "phenotype": "Type 2 diabetes",
  "mesh_id": "D003924",
  "mesh_term": "Diabetes Mellitus, Type 2",
  "sample_count": 5234,
  "project_count": 45,
  "marker_count": 89,
  "category": "Metabolic Diseases"
}
```

---

## Curated Metadata Categories

GMrepo provides manually curated metadata for each sample:

| Category | Description | Data Type |
|----------|-------------|-----------|
| Phenotype | Disease/health status | MeSH standardized |
| Age | Subject age | Integer |
| Sex | Subject sex | Male/Female |
| BMI | Body mass index | Float |
| Country | Geographic origin | ISO country |
| Antibiotics | Recent antibiotic use | Yes/No/Unknown |
| Sequencing | 16S or WGS | Enum |

---

## Sample Phenotype Categories (v3)

### Top Disease Categories

| Category | Sample Count (approx) |
|----------|----------------------|
| Metabolic Disorders | 25,000+ |
| Gastrointestinal Diseases | 20,000+ |
| Autoimmune Conditions | 8,000+ |
| Neurological Disorders | 5,000+ |
| Cardiovascular Diseases | 4,000+ |
| Cancer | 3,000+ |
| Healthy Controls | 40,000+ |

---

## API Usage Examples

### Fetch Runs by Project
```bash
curl "https://gmrepo.humangut.info/api/getRunsByProjectID?projectID=PRJEB6070"
```

### Fetch Samples by Phenotype
```bash
curl "https://gmrepo.humangut.info/api/getRunsByPhenotype?phenotype=Type%202%20diabetes"
```

### Get Species Abundance
```bash
curl "https://gmrepo.humangut.info/api/getSpeciesAbundance?runID=SRR1234567"
```

### Get Marker Taxa
```bash
curl "https://gmrepo.humangut.info/api/getMarkerTaxa?phenotype=Type%202%20diabetes"
```

### Get All Phenotypes
```bash
curl "https://gmrepo.humangut.info/api/getPhenotypeList"
```

---

## Programmatic Access

### Python Example

```python
import requests

base_url = "https://gmrepo.humangut.info/api"

# Fetch samples by phenotype
response = requests.get(f"{base_url}/getRunsByPhenotype",
                       params={"phenotype": "Type 2 diabetes"})
samples = response.json()

# Get species abundance for a sample
response = requests.get(f"{base_url}/getSpeciesAbundance",
                       params={"runID": samples[0]["run_id"]})
abundance = response.json()
```

### R Example

```r
library(httr)
library(jsonlite)

base_url <- "https://gmrepo.humangut.info/api"

# Fetch marker taxa
response <- GET(paste0(base_url, "/getMarkerTaxa"),
                query = list(phenotype = "Type 2 diabetes"))
markers <- fromJSON(content(response, "text"))
```

### Perl Example

```perl
use LWP::UserAgent;
use JSON;

my $ua = LWP::UserAgent->new;
my $base_url = "https://gmrepo.humangut.info/api";

my $response = $ua->get("$base_url/getRunsByPhenotype?phenotype=Type%202%20diabetes");
my $data = decode_json($response->content);
```

---

## Key Features

### Disease Marker Identification
- Pre-computed disease-associated taxa
- Cross-dataset statistical validation
- Effect size and significance metrics

### Cross-Dataset Comparison
- Standardized phenotype annotation (MeSH)
- Consistent taxonomic assignment
- Meta-analysis ready

### Graphical Query Builder
- Complex sample filtering
- Multi-criteria searches
- Export capabilities

---

## Data Quality

### Curation Process
1. Manual phenotype annotation from original publications
2. MeSH standardization for cross-study comparison
3. Quality control on sequencing metadata
4. Consistent taxonomic classification using standard pipelines

### Taxonomic Assignment
- 16S: SILVA/Greengenes databases
- WGS: MetaPhlAn/Kraken2 pipelines
- Relative abundance normalization

---

## Cross-References

| Database | Identifier Type | Usage |
|----------|-----------------|-------|
| NCBI Taxonomy | taxon_id | Microbial species identification |
| MeSH | mesh_id | Disease/phenotype standardization |
| SRA/ENA | run_id, project_id | Raw sequence data access |
| BioSample | sample_id | Sample metadata |

---

## Integration Notes

### Commercial Use
- License: CC BY-NC 3.0
- Commercial use requires separate permission
- Contact: weihuachen@hust.edu.cn

### Data Download
- Processed abundance data available
- Raw sequences via SRA/ENA links
- Bulk download options available

### Update Frequency
- Active curation and updates
- New projects added regularly
- Systematic literature review for data collection

---

## Recognition

| Designation | Year |
|-------------|------|
| ELIXIR Recommended Resource | 2022 |
| Database Commons Listed | 2020 |
| NAR Database Issue | 2020, 2022, 2025 |

### Database Ranking
- 405th of 6,933 databases overall (94.17th percentile)
- 41st of 723 metadata databases

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `run_id` | Sequencing run accession identifier from SRA/ENA archives | `SRR1234567` |
| `project_id` | BioProject accession identifying a collection of related samples | `PRJEB6070` |
| `sample_id` | BioSample accession for individual biological samples | `SAMN1234567` |
| `phenotype` | Disease or health status annotation for a sample | `Type 2 diabetes` |
| `phenotype_mesh_id` | Medical Subject Headings identifier for standardized phenotype | `D003924` |
| `relative_abundance` | Proportion of a taxon in a sample (0-1 scale) | `0.2534` |
| `taxon_id` | NCBI Taxonomy identifier for microbial species or genus | `816` |
| `direction` | Whether a marker taxon is increased or decreased in disease | `increased` |
| `effect_size` | Magnitude of association between taxon and phenotype | `0.45` |
| `sequencing_type` | Method used for microbiome profiling | `16S` or `WGS` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Marker Taxa | Microbial taxa statistically associated with specific diseases or phenotypes | phenotype, effect_size |
| 16S rRNA | Ribosomal RNA gene used for bacterial identification and taxonomy | sequencing_type |
| Metagenomic (WGS) | Whole-genome shotgun sequencing of all DNA in a sample | sequencing_type |
| Gut Microbiome | Community of microorganisms inhabiting the human gastrointestinal tract | phenotype |
| Cross-Dataset Comparison | Analysis methods for comparing results across different studies | phenotype_mesh_id |
| Taxonomic Assignment | Process of identifying microbial species from sequence data | taxon_id |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GMrepo | Gut Microbiome Repository | Human gut microbiome database |
| MeSH | Medical Subject Headings | NIH controlled vocabulary for indexing |
| SRA | Sequence Read Archive | NCBI raw sequence data repository |
| ENA | European Nucleotide Archive | EBI sequence data repository |
| WGS | Whole-Genome Shotgun | Metagenomic sequencing method |
| 16S | 16S Ribosomal RNA | Bacterial identification gene |
| NCBI | National Center for Biotechnology Information | US bioinformatics resource |
| SILVA | - | 16S/18S rRNA sequence database |
| MetaPhlAn | Metagenomic Phylogenetic Analysis | Taxonomic profiling tool |
| ELIXIR | European Life Sciences Infrastructure | European bioinformatics network |

---

## References

1. Wu S, et al. (2025). GMrepo v3: a curated human gut microbiome database with expanded disease coverage. Nucleic Acids Research. https://doi.org/10.1093/nar/gkaf1190

2. Dai D, et al. (2022). GMrepo v2: a curated human gut microbiome database with special focus on disease markers. Nucleic Acids Research. 50(D1):D777-D784. https://doi.org/10.1093/nar/gkab1019

3. Wu S, et al. (2020). GMrepo: a database of curated and consistently annotated human gut metagenomes. Nucleic Acids Research. 48(D1):D545-D553. https://doi.org/10.1093/nar/gkz764
