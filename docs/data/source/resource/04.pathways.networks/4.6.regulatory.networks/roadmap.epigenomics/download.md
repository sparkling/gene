---
id: download-roadmap-epigenomics
title: "Roadmap Epigenomics Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# Roadmap Epigenomics - Download Documentation

## Overview

Roadmap Epigenomics data is available through multiple portals including the WashU Epigenome Browser, AWS, and NCBI GEO. Data volumes exceed 150 TB total.

## Primary Data Portal

### Base URL

```
https://egg2.wustl.edu/roadmap/data/
```

### Directory Structure

```
roadmap/data/
├── byFileType/
│   ├── chromhmmSegmentations/
│   ├── signal/
│   ├── peaks/
│   ├── alignments/
│   └── dnamethylation/
├── byDataType/
│   ├── chipseq/
│   ├── dnase/
│   ├── rnaseq/
│   └── dnamethylation/
└── byReference/
    ├── hg19/
    └── GRCh38/
```

## Chromatin State Files

### 15-State Model

```bash
# Single epigenome
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E003_15_coreMarks_mnemonics.bed.gz

# All epigenomes
for eid in E{001..127}; do
  wget "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/${eid}_15_coreMarks_mnemonics.bed.gz"
done
```

### 18-State Model (Extended)

```bash
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E003_18_core_K27ac_mnemonics.bed.gz
```

## Signal Tracks (bigWig)

### Histone Modifications

```bash
# H3K4me3 signal (-log10 p-value)
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K4me3.pval.signal.bigwig

# H3K27ac signal
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K27ac.pval.signal.bigwig

# H3K4me1 signal
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K4me1.pval.signal.bigwig

# H3K27me3 signal
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K27me3.pval.signal.bigwig

# H3K36me3 signal
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K36me3.pval.signal.bigwig

# Fold change over control
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-H3K4me3.fc.signal.bigwig
```

### DNase-seq / ATAC-seq

```bash
# DNase signal
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-DNase.pval.signal.bigwig
```

## Peak Files

```bash
# Narrow peaks
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E003-H3K4me3.narrowPeak.gz

# Broad peaks (for broad marks like H3K27me3)
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E003-H3K27me3.broadPeak.gz

# Gapped peaks
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/gappedPeak/E003-H3K27me3.gappedPeak.gz
```

## DNA Methylation

### Whole Genome Bisulfite Sequencing

```bash
# Fractional methylation bigWig
wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/E003_WGBS_FractionalMethylation.bigwig

# Read coverage bigWig
wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/ReadCoverage_bigwig/E003_WGBS_ReadCoverage.bigwig

# BED format
wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/bed/E003.WGBS.bed.gz
```

### RRBS

```bash
wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/RRBS/E003_RRBS_FractionalMethylation.bigwig
```

## RNA-seq Expression

```bash
# RPKM values (protein coding)
wget https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.RPKM.pc.gz

# All genes RPKM
wget https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.RPKM.all.gz

# N (read counts)
wget https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.N.pc.gz
```

## Metadata Files

```bash
# Epigenome metadata
wget https://egg2.wustl.edu/roadmap/data/byFileType/metadata/EID_metadata.tab

# Sample-experiment mapping
wget https://egg2.wustl.edu/roadmap/data/byFileType/metadata/SampleID_to_EID.tab

# Mark availability matrix
wget https://egg2.wustl.edu/roadmap/data/byFileType/metadata/EID_markAvailability.tab
```

## Python Examples

### Download Multiple Files

```python
import requests
import os
from concurrent.futures import ThreadPoolExecutor

def download_file(url, output_dir):
    """Download a file."""
    filename = url.split('/')[-1]
    filepath = os.path.join(output_dir, filename)

    response = requests.get(url, stream=True)
    with open(filepath, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    return filepath

def download_chromatin_states(epigenomes, output_dir):
    """Download chromatin state files for multiple epigenomes."""
    base_url = "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final"

    urls = [
        f"{base_url}/{eid}_15_coreMarks_mnemonics.bed.gz"
        for eid in epigenomes
    ]

    os.makedirs(output_dir, exist_ok=True)

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(
            lambda url: download_file(url, output_dir),
            urls
        ))

    return results

# Example
epigenomes = ['E003', 'E004', 'E005', 'E006']
download_chromatin_states(epigenomes, './chromatin_states/')
```

### Process Chromatin States

```python
import pandas as pd
import pyBigWig

def load_chromatin_states(filepath):
    """Load chromatin state BED file."""
    return pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'state'],
        compression='gzip'
    )

def get_state_distribution(states_df):
    """Calculate genome coverage by state."""
    states_df['length'] = states_df['end'] - states_df['start']
    distribution = states_df.groupby('state')['length'].sum()
    distribution = distribution / distribution.sum() * 100
    return distribution.sort_values(ascending=False)

def get_signal_at_regions(bigwig_file, bed_df):
    """Get mean signal values at BED regions."""
    bw = pyBigWig.open(bigwig_file)

    signals = []
    for _, row in bed_df.iterrows():
        try:
            val = bw.stats(row['chrom'], row['start'], row['end'], type='mean')[0]
            signals.append(val if val is not None else 0)
        except:
            signals.append(0)

    bw.close()
    return signals
```

### Query Specific Regions

```python
def get_enhancers_near_gene(states_df, gene_tss, window=100000):
    """Find enhancer states near a gene TSS."""
    chrom, pos = gene_tss

    # Filter to window
    nearby = states_df[
        (states_df['chrom'] == chrom) &
        (states_df['start'] >= pos - window) &
        (states_df['end'] <= pos + window)
    ]

    # Filter to enhancer states
    enhancer_states = ['6_EnhG', '7_Enh', '12_EnhBiv']
    enhancers = nearby[nearby['state'].isin(enhancer_states)]

    return enhancers
```

## UCSC Genome Browser

### Track Hub Access

```
hub https://egg2.wustl.edu/roadmap/data/byFileType/trackhub/hub.txt
```

### Direct Track URLs

```bash
# bigBed chromatin states
https://egg2.wustl.edu/roadmap/data/byFileType/trackhub/hg19/E003.15_coreMarks.bb

# bigWig signal
https://egg2.wustl.edu/roadmap/data/byFileType/trackhub/hg19/E003-H3K4me3.pval.signal.bigwig
```

## AWS Access

```bash
# S3 bucket (may require AWS credentials)
aws s3 ls s3://encode-public/2012/07/01/
```

## Data Size Estimates

| Data Type | Size per Epigenome | Total |
|-----------|-------------------|-------|
| Chromatin states | ~5 MB | ~600 MB |
| Signal tracks (5 marks) | ~2 GB | ~250 GB |
| DNA methylation | ~500 MB | ~60 GB |
| Alignments | ~10 GB | >1 TB |
| All processed data | - | ~150 TB |

---

## Dataset Versions

### Current Release: Roadmap Epigenomics (Final)

| Property | Value |
|----------|-------|
| Version | Final Release |
| Release Date | 2015-02-18 |
| Total Size | ~150 TB |
| Epigenomes | 127 |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| chromhmmSegmentations (15-state) | ~600 MB | 127 | Chromatin states |
| signal (bigWig) | ~250 GB | 127 x 5 marks | Histone signals |
| dnamethylation (WGBS) | ~60 GB | ~40 | DNA methylation |
| rna/expression | ~1 GB | 57 | RNA-seq RPKM |

### Data Freeze Information

| Aspect | Value |
|--------|-------|
| Freeze Date | 2015-02-18 |
| Reference Genome | hg19 (primary), GRCh38 (converted) |
| Publication | Nature 2015 |
| Status | Final (no updates planned) |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| Base URL | `https://egg2.wustl.edu/roadmap/data/` |
| Authentication | None required |
| Rate Limit | No limit |
| Response Format | bigWig, BED, TSV |

### UCSC Track Hub

| Property | Value |
|----------|-------|
| Hub URL | `https://egg2.wustl.edu/roadmap/data/byFileType/trackhub/hub.txt` |
| Genome | hg19 |

### Directory Structure

| Path | Content |
|------|---------|
| `/byFileType/chromhmmSegmentations/` | Chromatin state BED files |
| `/byFileType/signal/` | Signal bigWig files |
| `/byDataType/dnamethylation/` | DNA methylation data |
| `/byDataType/rna/expression/` | RNA-seq expression |

---

## See Also

- [Schema Documentation](./schema.md)
- [Roadmap Analysis Pipeline](http://egg2.wustl.edu/roadmap/web_portal/processed_data.html)
