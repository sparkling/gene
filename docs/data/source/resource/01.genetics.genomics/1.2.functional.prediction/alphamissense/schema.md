---
id: schema-alphamissense
title: AlphaMissense Schema Documentation
type: schema
parent: README.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, variant, pathogenicity]
---

**Parent:** [Schema Documentation](./README.md)

# AlphaMissense Schema Documentation

**Document ID:** SCHEMA-ALPHAMISSENSE
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** Google DeepMind, Zenodo Repository, VEP Plugin Documentation

---

## TL;DR

AlphaMissense is Google DeepMind's AI-powered pathogenicity predictor for all possible human missense variants. Pre-computed predictions cover 71 million single nucleotide missense variants across 19,000 canonical transcripts, with calibrated pathogenicity scores (0-1) and three-class classifications. The model achieved state-of-the-art performance, correctly classifying 89% of ClinVar pathogenic and 90% of benign variants.

**Key Facts:**
- **Coverage:** 71M missense variants across 19K canonical transcripts
- **Score Range:** 0.0 (benign) to 1.0 (pathogenic)
- **Thresholds:** <0.34 = likely_benign, >0.564 = likely_pathogenic
- **License:** CC BY 4.0 (commercial-friendly for predictions)
- **File Size:** ~643 MB (hg38), ~622 MB (hg19) compressed

---

## Overview

AlphaMissense predicts the pathogenicity of all possible single amino acid substitutions in the human proteome. Built on AlphaFold2's protein structure prediction architecture, it incorporates evolutionary conservation, structural context, and population frequency data to generate calibrated pathogenicity scores.

**Publication:** Cheng J, et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. DOI: 10.1126/science.adg7492

**Repository:** https://github.com/google-deepmind/alphamissense (archived May 2025)

**Data Access:** Google Cloud Storage: `gs://dm_alphamissense`

---

## License Information

| Component | License | Commercial Use |
|-----------|---------|----------------|
| **Pre-computed Predictions** | CC BY 4.0 | Yes - with attribution |
| **Source Code** | Apache 2.0 | Yes |
| **Model Weights** | Not released | N/A |

**Attribution Requirement:** Cite the Science publication when using predictions.

**Note:** The Zenodo-hosted version uses CC BY-NC-SA 4.0 (non-commercial), but the official Google Cloud release uses CC BY 4.0.

---

## Database Statistics

| Metric | Value |
|--------|-------|
| **Total Variants (hg38)** | 71,000,000+ |
| **Canonical Transcripts** | 19,233 |
| **Isoform Transcripts** | ~60,000 |
| **Amino Acid Substitutions** | 216,000,000+ |
| **GENCODE Version (hg38)** | V32 |
| **GENCODE Version (hg19)** | V27 |
| **UniProt Release** | 2021_02 |
| **Reference Genome (hg38)** | GRCh38.p13 |
| **Reference Genome (hg19)** | GRCh37.p13 |

---

## File Inventory

| File | Size | Description |
|------|------|-------------|
| `AlphaMissense_hg38.tsv.gz` | 643 MB | All SNM predictions, hg38 coordinates |
| `AlphaMissense_hg19.tsv.gz` | 622 MB | All SNM predictions, hg19 coordinates |
| `AlphaMissense_gene_hg38.tsv.gz` | 254 KB | Gene-level averages, hg38 |
| `AlphaMissense_gene_hg19.tsv.gz` | 244 KB | Gene-level averages, hg19 |
| `AlphaMissense_aa_substitutions.tsv.gz` | 1.2 GB | All AA substitutions (216M) |
| `AlphaMissense_isoforms_hg38.tsv.gz` | 1.2 GB | Non-canonical isoforms, hg38 |
| `AlphaMissense_isoforms_aa_substitutions.tsv.gz` | 2.5 GB | Isoform AA substitutions |

**Total Dataset Size:** ~6.1 GB compressed

---

## Genomic Variant File Schema

### File: AlphaMissense_hg38.tsv.gz / AlphaMissense_hg19.tsv.gz

Tab-separated values, bgzip compressed, sorted by genomic coordinates.

**Header Line:**
```
#CHROM	POS	REF	ALT	genome	uniprot_id	transcript_id	protein_variant	am_pathogenicity	am_class
```

### Column Definitions

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | string | Chromosome identifier: `chr<N>` where N in {1-22, X, Y, M} |
| `POS` | integer | Genomic position (1-based coordinate) |
| `REF` | string | Reference nucleotide allele (+ strand) |
| `ALT` | string | Alternative nucleotide allele (+ strand) |
| `genome` | string | Genome build: `hg38` or `hg19` |
| `uniprot_id` | string | UniProtKB accession (UniProt release 2021_02) |
| `transcript_id` | string | Ensembl transcript ID (GENCODE annotation) |
| `protein_variant` | string | Amino acid change: `<Ref_AA><Position><Alt_AA>` (e.g., `V2L`) |
| `am_pathogenicity` | float | Pathogenicity score [0.0-1.0] |
| `am_class` | string | Classification: `likely_benign`, `ambiguous`, `likely_pathogenic` |

### Field Details

#### CHROM
- **Format:** `chr1`, `chr2`, ..., `chr22`, `chrX`, `chrY`, `chrM`
- **Note:** Includes mitochondrial variants on `chrM`

#### POS
- **Coordinate System:** 1-based (VCF-style)
- **Reference:** GRCh38.p13 (hg38) or GRCh37.p13 (hg19)

#### REF / ALT
- **Values:** Single nucleotide (A, C, G, T)
- **Strand:** Positive strand orientation
- **Scope:** Only SNVs causing missense changes

#### uniprot_id
- **Format:** UniProt accession (e.g., `P00533`, `Q9Y6K9`)
- **Source:** UniProt release 2021_02
- **Note:** Canonical isoform only in main files

#### transcript_id
- **Format:** Ensembl transcript ID (e.g., `ENST00000275493.7`)
- **Source:** GENCODE V32 (hg38) or V27 (hg19)
- **Note:** Includes version number

#### protein_variant
- **Format:** `<Reference_AA><Position><Alternate_AA>`
- **Example:** `V600E` (Valine at position 600 to Glutamic acid)
- **Position:** 1-based within protein sequence
- **Amino Acids:** Single-letter codes (A-Y, excluding J, O, U)

#### am_pathogenicity
- **Range:** 0.0 to 1.0 (continuous)
- **Interpretation:** Predicted probability of clinical pathogenicity
- **Calibration:** Trained on ClinVar pathogenic/benign variants
- **Missing:** No missing values in canonical files

#### am_class
Classification derived from pathogenicity score thresholds:

| Class | Score Range | Interpretation |
|-------|-------------|----------------|
| `likely_benign` | < 0.34 | Low probability of pathogenicity |
| `ambiguous` | 0.34 - 0.564 | Uncertain classification |
| `likely_pathogenic` | > 0.564 | High probability of pathogenicity |

**Threshold Selection:** Optimized using ClinVar validation set to minimize false positives/negatives.

---

## Gene-Level File Schema

### File: AlphaMissense_gene_hg38.tsv.gz / AlphaMissense_gene_hg19.tsv.gz

| Column | Type | Description |
|--------|------|-------------|
| `transcript_id` | string | Ensembl transcript ID |
| `mean_am_pathogenicity` | float | Average pathogenicity across all missense variants |

**Example:**
```
transcript_id	mean_am_pathogenicity
ENST00000000233.5	0.74227
ENST00000000412.3	0.378343
ENST00000001008.4	0.42229
```

**Use Case:** Identify genes with high overall pathogenicity burden.

---

## Amino Acid Substitution File Schema

### File: AlphaMissense_aa_substitutions.tsv.gz

Protein-centric view of all possible amino acid substitutions.

| Column | Type | Description |
|--------|------|-------------|
| `uniprot_id` | string | UniProtKB accession |
| `protein_variant` | string | Amino acid change |
| `am_pathogenicity` | float | Pathogenicity score |
| `am_class` | string | Classification |

**Coverage:** 216M substitutions across 20K UniProt canonical isoforms

---

## Sample Data

### Genomic Variants (hg38)
```tsv
#CHROM	POS	REF	ALT	genome	uniprot_id	transcript_id	protein_variant	am_pathogenicity	am_class
chr1	69094	G	A	hg38	Q8NH21	ENST00000641515.2	G2R	0.0897	likely_benign
chr1	69094	G	C	hg38	Q8NH21	ENST00000641515.2	G2A	0.1034	likely_benign
chr1	69094	G	T	hg38	Q8NH21	ENST00000641515.2	G2V	0.1156	likely_benign
chr7	140753336	A	T	hg38	P15056	ENST00000288602.11	V600E	0.9972	likely_pathogenic
chr17	7673803	G	A	hg38	P04637	ENST00000269305.9	R175H	0.9856	likely_pathogenic
```

### Gene-Level Averages
```tsv
transcript_id	mean_am_pathogenicity
ENST00000366794.4	0.8234
ENST00000380152.8	0.4521
ENST00000288602.11	0.6892
```

---

## Classification Performance

### ClinVar Validation

| Metric | Value |
|--------|-------|
| **Pathogenic Recall** | 89% classified as likely_pathogenic |
| **Benign Recall** | 90% classified as likely_benign |
| **Ambiguous Rate** | ~11-10% |

### Comparison with Other Predictors

| Predictor | AUROC (ClinVar) |
|-----------|-----------------|
| **AlphaMissense** | 0.940 |
| REVEL | 0.908 |
| CADD | 0.866 |
| PolyPhen-2 | 0.823 |
| SIFT | 0.789 |

---

## Download Instructions

### Google Cloud Storage (Official)

```bash
# List available files
gsutil ls gs://dm_alphamissense/

# Download hg38 predictions
gsutil cp gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz .

# Download all files
gsutil -m cp -r gs://dm_alphamissense/* ./alphamissense/
```

### Zenodo (Alternative)

**Record:** https://zenodo.org/records/8360242

```bash
# Download via wget
wget https://zenodo.org/records/8360242/files/AlphaMissense_hg38.tsv.gz
```

### Tabix Indexing

For efficient querying, create tabix index:

```bash
# Index the file
tabix -s 1 -b 2 -e 2 -S 1 AlphaMissense_hg38.tsv.gz

# Query specific region
tabix AlphaMissense_hg38.tsv.gz chr7:140753300-140753400
```

---

## Integration Examples

### Python: Parse TSV File

```python
import gzip
import csv

def parse_alphamissense(filepath: str):
    """Parse AlphaMissense TSV file."""
    with gzip.open(filepath, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            yield {
                'chrom': row['#CHROM'],
                'pos': int(row['POS']),
                'ref': row['REF'],
                'alt': row['ALT'],
                'uniprot_id': row['uniprot_id'],
                'transcript_id': row['transcript_id'],
                'protein_variant': row['protein_variant'],
                'pathogenicity': float(row['am_pathogenicity']),
                'classification': row['am_class']
            }

# Filter pathogenic variants
pathogenic = [v for v in parse_alphamissense('AlphaMissense_hg38.tsv.gz')
              if v['classification'] == 'likely_pathogenic']
```

### Python: Query with Pysam

```python
import pysam

def query_alphamissense(filepath: str, chrom: str, start: int, end: int):
    """Query AlphaMissense with tabix index."""
    tbx = pysam.TabixFile(filepath)

    for row in tbx.fetch(chrom, start, end):
        fields = row.split('\t')
        yield {
            'chrom': fields[0],
            'pos': int(fields[1]),
            'ref': fields[2],
            'alt': fields[3],
            'pathogenicity': float(fields[8]),
            'classification': fields[9]
        }

# Query BRAF V600 region
variants = list(query_alphamissense(
    'AlphaMissense_hg38.tsv.gz',
    'chr7', 140753300, 140753400
))
```

### VEP Plugin Integration

```bash
# Run VEP with AlphaMissense plugin
vep -i input.vcf \
    --plugin AlphaMissense,file=AlphaMissense_hg38.tsv.gz \
    -o output.vcf
```

---

## Known Limitations

1. **Canonical Transcripts Only:** Main files cover only canonical isoforms; use isoform files for complete coverage

2. **Missense Only:** Does not predict frameshift, nonsense, or splice variants

3. **Model Weights Unavailable:** Cannot generate new predictions; limited to pre-computed data

4. **UniProt Version:** Based on 2021_02 release; newer proteins not covered

5. **No Structural Variants:** Only single nucleotide substitutions

6. **Archived Repository:** Code repository archived May 2025; no future updates planned

---

## Download

### Google Cloud Storage (Official)

The official and recommended source for AlphaMissense data:

```bash
# List available files
gsutil ls gs://dm_alphamissense/

# Download hg38 predictions (main file)
gsutil cp gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz .

# Download all files
gsutil -m cp -r gs://dm_alphamissense/* ./alphamissense/
```

**URL:** `gs://dm_alphamissense/`

### Zenodo (Alternative)

Alternative access for users without Google Cloud access:

**Repository:** https://zenodo.org/records/8360242

Files available:
- AlphaMissense_hg38.tsv.gz
- AlphaMissense_hg19.tsv.gz
- AlphaMissense_gene_hg38.tsv.gz
- AlphaMissense_gene_hg19.tsv.gz
- AlphaMissense_aa_substitutions.tsv.gz
- AlphaMissense_isoforms_hg38.tsv.gz
- AlphaMissense_isoforms_aa_substitutions.tsv.gz

```bash
# Download via wget
wget https://zenodo.org/records/8360242/files/AlphaMissense_hg38.tsv.gz
```

### File Sizes

| File | Compressed Size | Uncompressed |
|------|-----------------|--------------|
| AlphaMissense_hg38.tsv.gz | 643 MB | ~2.1 GB |
| AlphaMissense_hg19.tsv.gz | 622 MB | ~2.0 GB |
| AlphaMissense_gene_hg38.tsv.gz | 254 KB | ~1.2 MB |
| AlphaMissense_aa_substitutions.tsv.gz | 1.2 GB | ~4.0 GB |
| AlphaMissense_isoforms_hg38.tsv.gz | 1.2 GB | ~3.8 GB |

**Total Dataset:** ~6.1 GB compressed (~13 GB uncompressed)

---

## Data Format

| Format | Description | Availability |
|--------|-------------|--------------|
| **Primary** | TSV (Tab-Separated Values), gzip compressed (.tsv.gz) | All files |
| **Compression** | gzip (.gz) | All files |
| **Encoding** | UTF-8 | All files |
| **Variant Format** | Genomic coordinates (1-based) with VCF-style chromosome names | Main files |
| **Protein Format** | UniProt accession + protein coordinates | AA substitution files |

### TSV Structure

Files are tab-separated with header line (prefixed with #) and sorted by genomic position:
- **Tab separator** for all fields
- **Header row** with column names
- **No quoted fields** - raw tab-delimited format
- **One variant per row** for genomic files
- **One substitution per row** for amino acid files

### Coordinate System

- **Genomic positions:** 1-based (VCF-style)
- **Protein positions:** 1-based within protein sequence
- **Chromosome format:** chr1, chr2, ..., chr22, chrX, chrY, chrM

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "CHROM:POS:REF:ALT" |
| `name` | string | Entity name | "chr7:140753336:A:T" |
| `type` | string | Record type | "variant" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 71,000,000+ |
| Storage | 643 MB (compressed, hg38) |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `CHROM` | Chromosome identifier in VCF-style format | `chr7` |
| `POS` | Genomic position using 1-based coordinate system | `140753336` |
| `REF` | Reference nucleotide allele on the positive strand | `A` |
| `ALT` | Alternative nucleotide allele on the positive strand | `T` |
| `am_pathogenicity` | Pathogenicity score from 0.0 (benign) to 1.0 (pathogenic) | `0.9972` |
| `am_class` | Classification based on pathogenicity thresholds | `likely_pathogenic` |
| `protein_variant` | Amino acid change in format RefAA-Position-AltAA | `V600E` |
| `uniprot_id` | UniProtKB accession for the canonical protein isoform | `P15056` |
| `transcript_id` | Ensembl transcript identifier with version number | `ENST00000288602.11` |
| `mean_am_pathogenicity` | Average pathogenicity score across all missense variants for a gene | `0.74227` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Missense Variant | Single nucleotide variant that changes one amino acid to another | protein_variant |
| Pathogenicity Score | Predicted probability that a variant causes clinical disease | am_pathogenicity |
| Canonical Transcript | Primary transcript isoform selected for a gene | transcript_id |
| Calibration | Training process using known pathogenic/benign variants from ClinVar | am_pathogenicity |
| AUROC | Area Under Receiver Operating Characteristic curve, measures prediction accuracy | Classification Performance |
| Tabix Index | Index format enabling rapid genomic coordinate queries | POS, CHROM |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| SNM | Single Nucleotide Missense | Variant type covered by AlphaMissense |
| SNV | Single Nucleotide Variant | General term for single-base changes |
| GENCODE | GENCODE gene annotation project | Source of transcript annotations |
| GRCh38 | Genome Reference Consortium Human Build 38 | Primary reference genome (hg38) |
| GRCh37 | Genome Reference Consortium Human Build 37 | Legacy reference genome (hg19) |
| VEP | Variant Effect Predictor | Ensembl tool with AlphaMissense plugin |
| ClinVar | Clinical Variant database | Validation dataset for pathogenicity |
| REVEL | Rare Exome Variant Ensemble Learner | Competing pathogenicity predictor |
| CADD | Combined Annotation Dependent Depletion | Competing variant scoring method |
| dbNSFP | Database for Nonsynonymous SNPs Functional Predictions | Integrates AlphaMissense scores |

---

## Related Resources

- **Ensembl VEP Plugin:** https://github.com/Ensembl/VEP_plugins/blob/release/115/AlphaMissense.pm
- **AlphaMissense Browser:** https://alphamissense.hegelab.org/
- **dbNSFP Integration:** AlphaMissense scores included in dbNSFP v4.4+
- **ClinVar:** Validation dataset source

---

## References

1. Cheng J, Novati G, Pan J, et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 381(6664):eadg7492. DOI: 10.1126/science.adg7492

2. Jumper J, Evans R, Pritzel A, et al. (2021). Highly accurate protein structure prediction with AlphaFold. Nature. 596:583-589.

3. GENCODE Consortium. Human Release 32 (GRCh38.p13). https://www.gencodegenes.org/human/release_32.html

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
