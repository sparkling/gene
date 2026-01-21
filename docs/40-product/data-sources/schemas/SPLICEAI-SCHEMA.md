# SpliceAI Schema Documentation

**Document ID:** SCHEMA-SPLICEAI
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** SpliceAI GitHub, Broad Institute Lookup API

---

## TL;DR

SpliceAI is a deep learning tool predicting splice-altering effects of genetic variants. It outputs delta scores (0-1) for acceptor/donor gain/loss events. The Broad Institute provides a free lookup API for academic use; commercial use requires an Illumina license. Key integration points: VCF INFO field annotation, HGVS query support, and multi-transcript analysis.

---

## Overview

SpliceAI uses a 32-layer deep neural network trained on RNA-seq data to predict cryptic splicing effects from primary sequence. It annotates SNVs and simple indels with probability scores for four splice-altering mechanisms: acceptor gain, acceptor loss, donor gain, and donor loss.

**Primary Use Cases:**
- Variant interpretation in clinical genetics
- Filtering VUS (variants of uncertain significance)
- Prioritizing candidates in rare disease diagnosis
- Splicing research and annotation pipelines

---

## License & Access

| Aspect | Details |
|--------|---------|
| **Source Code** | PolyForm Strict License 1.0.0 |
| **Trained Models** | CC BY-NC 4.0 (non-commercial) |
| **Academic Use** | Free via Broad Institute Lookup |
| **Commercial Use** | Requires Illumina license |
| **Lookup API** | Rate-limited (handful of queries per-user per-minute) |

**Recommendation:** For batch processing (>100 variants), use local installation rather than the web API.

---

## API Endpoints

### Broad Institute SpliceAI Lookup

**Base URLs by Genome Version:**

| Service | Endpoint |
|---------|----------|
| SpliceAI hg37 | `https://spliceai-37-xwkwwwxdwq-uc.a.run.app` |
| SpliceAI hg38 | `https://spliceai-38-xwkwwwxdwq-uc.a.run.app` |
| Pangolin hg37 | `https://pangolin-37-xwkwwwxdwq-uc.a.run.app` |
| Pangolin hg38 | `https://pangolin-38-xwkwwwxdwq-uc.a.run.app` |
| Normalization | `https://liftover-xwkwwwxdwq-uc.a.run.app` |

**Web Interface:** https://spliceailookup.broadinstitute.org/

### Query Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `hg` | string | Yes | Genome version: `37` or `38` |
| `variant` | string | Yes | Normalized variant (e.g., `chr8-140300616-T-G`) |
| `distance` | integer | No | Maximum distance to splice sites in bp (default: 50, max: 10000) |
| `mask` | integer | No | Score masking: `0` = unmasked, `1` = masked |
| `bc` | string | No | Gencode set: `basic` or `comprehensive` |
| `raw` | string | No | Original input before normalization |
| `variant_consequence` | string | No | VEP consequence type |

### Example API Requests

```bash
# Simple variant query
curl "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/?hg=38&variant=chr17-43106533-C-T&distance=500"

# With comprehensive transcripts
curl "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/?hg=38&variant=chr17-43106533-C-T&distance=500&bc=comprehensive"
```

### Supported Input Formats

| Format | Example |
|--------|---------|
| Standard VCF | `chr8-140300616-T-G` |
| No chr prefix | `8-140300616-T-G` |
| Colon-separated | `1:930130:C:G` |
| HGVS coding | `NM_001089.3(ABCA3):c.875A>T` |
| HGVS genomic | `NC_000017.11:g.43106533C>T` |
| Insertion | `1-1042601-A-AGAGAG` |
| Deletion | `1-1042601-AGAGAG-A` |

**Supported Chromosomes:** 1-22, X, Y, M (mitochondrial)

---

## VCF INFO Field Format

### Field Structure

```
SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
```

### Field Dictionary

| Field | Type | Range | Description |
|-------|------|-------|-------------|
| `ALLELE` | string | - | Alternate allele sequence |
| `SYMBOL` | string | - | Gene symbol (HGNC) |
| `DS_AG` | float | 0.00-1.00 | Delta score: acceptor gain |
| `DS_AL` | float | 0.00-1.00 | Delta score: acceptor loss |
| `DS_DG` | float | 0.00-1.00 | Delta score: donor gain |
| `DS_DL` | float | 0.00-1.00 | Delta score: donor loss |
| `DP_AG` | integer | -10000 to +10000 | Delta position: acceptor gain (bp from variant) |
| `DP_AL` | integer | -10000 to +10000 | Delta position: acceptor loss (bp from variant) |
| `DP_DG` | integer | -10000 to +10000 | Delta position: donor gain (bp from variant) |
| `DP_DL` | integer | -10000 to +10000 | Delta position: donor loss (bp from variant) |

### Delta Score Interpretation

| Score Range | Confidence | Interpretation |
|-------------|------------|----------------|
| >= 0.80 | High precision | Strong evidence of splice-altering effect |
| >= 0.50 | Recommended | Likely splice-altering |
| >= 0.20 | High recall | Possible splice-altering effect |
| < 0.20 | Low | Unlikely to alter splicing |

### Delta Position Interpretation

- **Negative values:** Splice site change occurs upstream (5') of the variant
- **Positive values:** Splice site change occurs downstream (3') of the variant
- **Value magnitude:** Distance in base pairs from variant position

---

## Extended Output Fields

### Transcript Annotations

| Field | Description |
|-------|-------------|
| `transcript_id` | Ensembl transcript accession |
| `transcript_priority` | MANE Select, MANE Plus Clinical, Canonical, or empty |
| `biotype` | protein_coding, lncRNA, etc. |
| `strand` | + or - |

### Score Types

| Type | Description |
|------|-------------|
| `ref_scores` | Splice site strengths for reference haplotype |
| `alt_scores` | Splice site strengths for alternate haplotype |
| `delta_scores` | Difference (alt - ref) for each position |

### Gene Links

| Resource | Format |
|----------|--------|
| Ensembl | `ENSG00000141510` |
| OMIM | `191170` |
| GTEx | Link to expression data |

---

## Sample Data

### VCF Output Example

```vcf
##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI variant annotation. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
17	43106533	.	C	T	.	.	SpliceAI=T|BRCA1|0.00|0.00|0.91|0.08|-28|-31|2|-31
2	179446218	.	G	A	.	.	SpliceAI=A|TTN|0.07|0.00|0.00|0.00|35|1|-5|1
11	5227002	.	T	A	.	.	SpliceAI=A|HBB|0.00|0.02|0.00|0.00|-26|-26|-15|-15
```

### JSON API Response Example

```json
{
  "variant": "chr17-43106533-C-T",
  "scores": [
    {
      "gene": "BRCA1",
      "transcript": "ENST00000357654.9",
      "transcript_priority": "MANE Select",
      "biotype": "protein_coding",
      "strand": "-",
      "DS_AG": 0.00,
      "DS_AL": 0.00,
      "DS_DG": 0.91,
      "DS_DL": 0.08,
      "DP_AG": -28,
      "DP_AL": -31,
      "DP_DG": 2,
      "DP_DL": -31,
      "max_delta": 0.91,
      "max_delta_type": "donor_gain"
    }
  ],
  "assembly": "GRCh38",
  "distance": 500
}
```

### Multi-Transcript Example

```json
{
  "variant": "chr2-179446218-G-A",
  "scores": [
    {
      "gene": "TTN",
      "transcript": "ENST00000589042.6",
      "transcript_priority": "MANE Select",
      "DS_AG": 0.07,
      "DS_AL": 0.00,
      "DS_DG": 0.00,
      "DS_DL": 0.00
    },
    {
      "gene": "TTN",
      "transcript": "ENST00000460472.6",
      "transcript_priority": "",
      "DS_AG": 0.05,
      "DS_AL": 0.00,
      "DS_DG": 0.00,
      "DS_DL": 0.00
    }
  ]
}
```

---

## Command-Line Usage

### Installation

```bash
# Via pip
pip install spliceai

# Via conda
conda install -c bioconda spliceai
```

**Dependency:** TensorFlow >= 1.2.0

### Basic Usage

```bash
spliceai -I input.vcf -O output.vcf -R genome.fa -A grch38
```

### Parameters

| Flag | Required | Description |
|------|----------|-------------|
| `-I` | Yes | Input VCF file path |
| `-O` | Yes | Output VCF file path |
| `-R` | Yes | Reference genome FASTA |
| `-A` | Yes | Annotation: `grch37` or `grch38` (uses GENCODE V24) |
| `-D` | No | Maximum distance to splice sites (default: 50, max: 10000) |
| `-M` | No | Mask scores: `0` = unmasked, `1` = masked (default: 0) |

### Custom Annotation Files

For non-human species or custom transcripts:

```bash
spliceai -I input.vcf -O output.vcf -R genome.fa -A custom_annotations.txt
```

**Annotation file format (tab-delimited):**
```
#NAME	CHROM	STRAND	TX_START	TX_END	EXON_START	EXON_END
BRCA1	chr17	-	43044295	43170245	43044295,43047643,...	43045802,43049194,...
```

---

## Variant Scope & Limitations

### Supported Variants

| Type | Supported | Notes |
|------|-----------|-------|
| SNVs | Yes | Full support |
| Simple insertions | Yes | Most cases |
| Simple deletions | Yes | Length <= 2 x distance parameter |
| Complex indels | Partial | May have reduced accuracy |
| Structural variants | No | Not annotated |

### Excluded Regions

- Variants within 5kb of chromosome ends
- Intergenic variants (outside gene boundaries)
- Deletions longer than 2x the `-D` parameter

### Masking Option

When `-M 1` is set:
- Scores are masked for positions where the variant overlaps a reference splice site
- Useful for filtering out loss-of-function variants at canonical splice sites

---

## Integration Patterns

### Python: Parse VCF Annotations

```python
import re
from typing import Dict, List, Optional

def parse_spliceai_info(info_field: str) -> List[Dict]:
    """Parse SpliceAI INFO field from VCF."""
    results = []

    match = re.search(r'SpliceAI=([^;]+)', info_field)
    if not match:
        return results

    for annotation in match.group(1).split(','):
        fields = annotation.split('|')
        if len(fields) >= 10:
            results.append({
                'allele': fields[0],
                'gene': fields[1],
                'ds_ag': float(fields[2]),
                'ds_al': float(fields[3]),
                'ds_dg': float(fields[4]),
                'ds_dl': float(fields[5]),
                'dp_ag': int(fields[6]),
                'dp_al': int(fields[7]),
                'dp_dg': int(fields[8]),
                'dp_dl': int(fields[9]),
                'max_delta': max(float(fields[2]), float(fields[3]),
                                float(fields[4]), float(fields[5]))
            })

    return results

def get_max_score(annotations: List[Dict]) -> Optional[float]:
    """Get maximum delta score across all transcripts."""
    if not annotations:
        return None
    return max(a['max_delta'] for a in annotations)
```

### Python: Query Lookup API

```python
import requests
from typing import Dict, Optional

def query_spliceai(
    variant: str,
    genome: str = "38",
    distance: int = 500,
    gencode: str = "basic"
) -> Optional[Dict]:
    """Query SpliceAI Lookup API."""
    base_url = f"https://spliceai-{genome}-xwkwwwxdwq-uc.a.run.app"

    params = {
        "hg": genome,
        "variant": variant,
        "distance": distance,
        "bc": gencode,
        "mask": 0
    }

    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        print(f"API error: {e}")
        return None

# Example usage
result = query_spliceai("chr17-43106533-C-T", genome="38", distance=500)
```

### Batch Processing Script

```python
import subprocess
from pathlib import Path

def run_spliceai_batch(
    input_vcf: Path,
    output_vcf: Path,
    reference: Path,
    assembly: str = "grch38",
    distance: int = 500
) -> bool:
    """Run SpliceAI on a VCF file."""
    cmd = [
        "spliceai",
        "-I", str(input_vcf),
        "-O", str(output_vcf),
        "-R", str(reference),
        "-A", assembly,
        "-D", str(distance)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"SpliceAI error: {e.stderr}")
        return False
```

---

## Performance Considerations

| Scenario | Recommendation |
|----------|----------------|
| < 100 variants | Use Broad Institute Lookup API |
| 100-10,000 variants | Local installation recommended |
| > 10,000 variants | Local installation required, consider GPU |
| Production pipeline | Local installation with GPU support |

### GPU Acceleration

SpliceAI supports GPU via TensorFlow:

```bash
# Install TensorFlow with GPU support
pip install tensorflow-gpu

# Run with GPU
CUDA_VISIBLE_DEVICES=0 spliceai -I input.vcf -O output.vcf -R genome.fa -A grch38
```

---

## Related Resources

| Resource | URL |
|----------|-----|
| GitHub Repository | https://github.com/Illumina/SpliceAI |
| Broad Lookup | https://spliceailookup.broadinstitute.org/ |
| Publication | Jaganathan et al., Cell 2019 |
| GENCODE | https://www.gencodegenes.org/ |
| Pangolin (related) | https://github.com/tkzeng/Pangolin |

---

## References

1. Jaganathan K, et al. (2019). Predicting Splicing from Primary Sequence with Deep Learning. Cell. 176(3):535-548.e24. doi:10.1016/j.cell.2018.12.015
2. SpliceAI GitHub: https://github.com/Illumina/SpliceAI
3. Broad Institute Lookup: https://spliceailookup.broadinstitute.org/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
