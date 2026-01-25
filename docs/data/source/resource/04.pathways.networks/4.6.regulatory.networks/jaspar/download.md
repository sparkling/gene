---
id: download-jaspar
title: "JASPAR Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
---

# JASPAR - Download Documentation

## Overview

JASPAR provides transcription factor binding profiles via REST API and bulk downloads. All data is available under CC BY 4.0 license.

## REST API

### Base URL

```
https://jaspar.genereg.net/api/v1
```

### Get Matrix by ID

```bash
# Get matrix details
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106.3/"

# Get specific version
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106.1/"

# Get latest version
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106/"
```

### Get Matrix in Different Formats

```bash
# PFM format
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106.3/?format=pfm"

# MEME format
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106.3/?format=meme"

# TRANSFAC format
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106.3/?format=transfac"

# JASPAR format
curl "https://jaspar.genereg.net/api/v1/matrix/MA0106.3/?format=jaspar"
```

### Search Matrices

```bash
# Search by TF name
curl "https://jaspar.genereg.net/api/v1/matrix/?name=TP53"

# Search by collection
curl "https://jaspar.genereg.net/api/v1/matrix/?collection=CORE"

# Search by taxonomy
curl "https://jaspar.genereg.net/api/v1/matrix/?tax_group=vertebrates"

# Search by species
curl "https://jaspar.genereg.net/api/v1/matrix/?species=9606"

# Search by TF class
curl "https://jaspar.genereg.net/api/v1/matrix/?class=Zinc-coordinating"

# Combined search
curl "https://jaspar.genereg.net/api/v1/matrix/?collection=CORE&tax_group=vertebrates&class=Basic%20domains"
```

### Pagination

```bash
# With pagination
curl "https://jaspar.genereg.net/api/v1/matrix/?collection=CORE&page=1&page_size=100"
```

### Get Binding Sites

```bash
# Get sites for a matrix
curl "https://jaspar.genereg.net/api/v1/sites/MA0106.3/"
```

### List Available Options

```bash
# List collections
curl "https://jaspar.genereg.net/api/v1/collections/"

# List tax groups
curl "https://jaspar.genereg.net/api/v1/taxon/"

# List TF classes
curl "https://jaspar.genereg.net/api/v1/tffm/"
```

## Bulk Downloads

### Download Portal

```
https://jaspar.genereg.net/download/
```

### CORE Collection

```bash
# All vertebrate PFMs
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt

# All vertebrate PFMs (MEME format)
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

# All vertebrate PFMs (TRANSFAC format)
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_transfac.txt

# Plant PFMs
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt

# All taxa PFMs
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt
```

### Redundant Collections

```bash
# All redundant matrices (includes all versions)
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_redundant_pfms_jaspar.zip
```

## Output Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| jaspar | JASPAR text format | Native JASPAR |
| pfm | Simple PFM | Generic |
| meme | MEME format | MEME Suite tools |
| transfac | TRANSFAC format | Legacy tools |

## Python Examples

### API Access

```python
import requests

class JASPARClient:
    def __init__(self):
        self.base_url = "https://jaspar.genereg.net/api/v1"

    def get_matrix(self, matrix_id, format="json"):
        """Get matrix by ID."""
        url = f"{self.base_url}/matrix/{matrix_id}/"
        params = {"format": format} if format != "json" else {}

        response = requests.get(url, params=params)
        return response.json() if format == "json" else response.text

    def search_matrices(self, name=None, collection="CORE", tax_group=None):
        """Search for matrices."""
        params = {"collection": collection}
        if name:
            params["name"] = name
        if tax_group:
            params["tax_group"] = tax_group

        response = requests.get(f"{self.base_url}/matrix/", params=params)
        return response.json()

    def get_pfm(self, matrix_id):
        """Get PFM as dictionary."""
        data = self.get_matrix(matrix_id)
        return data.get("pfm", {})

# Example usage
client = JASPARClient()
tp53 = client.get_matrix("MA0106.3")
pfm = client.get_pfm("MA0106.3")
```

### Convert to PWM

```python
import numpy as np

def pfm_to_pwm(pfm, pseudocount=0.1, background=None):
    """Convert PFM to PWM (log-odds scores)."""
    if background is None:
        background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    # Convert to numpy array
    matrix = np.array([pfm["A"], pfm["C"], pfm["G"], pfm["T"]])

    # Add pseudocount
    matrix = matrix + pseudocount

    # Normalize to probabilities
    matrix = matrix / matrix.sum(axis=0)

    # Calculate log-odds
    bg = np.array([background["A"], background["C"], background["G"], background["T"]])
    pwm = np.log2(matrix / bg[:, np.newaxis])

    return {
        "A": pwm[0].tolist(),
        "C": pwm[1].tolist(),
        "G": pwm[2].tolist(),
        "T": pwm[3].tolist()
    }
```

### Score Sequence

```python
def score_sequence(sequence, pwm):
    """Score a DNA sequence using PWM."""
    nuc_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    matrix = np.array([pwm["A"], pwm["C"], pwm["G"], pwm["T"]])

    if len(sequence) != matrix.shape[1]:
        raise ValueError("Sequence length must match motif length")

    score = 0
    for i, nuc in enumerate(sequence.upper()):
        if nuc in nuc_to_idx:
            score += matrix[nuc_to_idx[nuc], i]

    return score

def scan_sequence(sequence, pwm, threshold=0.8):
    """Scan sequence for motif matches."""
    motif_length = len(pwm["A"])
    max_score = sum(max(pwm[n][i] for n in "ACGT") for i in range(motif_length))
    min_score = sum(min(pwm[n][i] for n in "ACGT") for i in range(motif_length))

    threshold_score = min_score + threshold * (max_score - min_score)

    hits = []
    for i in range(len(sequence) - motif_length + 1):
        subseq = sequence[i:i+motif_length]
        score = score_sequence(subseq, pwm)

        if score >= threshold_score:
            hits.append({
                "start": i,
                "end": i + motif_length,
                "sequence": subseq,
                "score": score,
                "relative_score": (score - min_score) / (max_score - min_score)
            })

    return hits
```

### Using Biopython

```python
from Bio import motifs
from Bio.motifs import jaspar

# Read JASPAR format file
with open("JASPAR2024_CORE_vertebrates.txt") as f:
    for m in jaspar.read(f, "jaspar"):
        print(f"{m.matrix_id}: {m.name}")
        print(m.counts)
        pwm = m.counts.normalize(pseudocounts=0.5)
        pssm = pwm.log_odds()
```

## MEME Suite Integration

```bash
# Scan sequences for motif
fimo --thresh 1e-4 JASPAR2024_CORE_vertebrates_pfms_meme.txt sequences.fasta

# Find motifs de novo and compare to JASPAR
meme sequences.fasta -dna -mod zoops -nmotifs 10
tomtom meme_out/meme.txt JASPAR2024_CORE_vertebrates_pfms_meme.txt
```

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| API | No strict limit |
| Recommended | 10 requests/second |
| Bulk download | Unlimited |

---

## Dataset Versions

### Current Release: JASPAR 2024

| Property | Value |
|----------|-------|
| Version | JASPAR 2024 |
| Release Date | 2024-01-01 |
| Total Profiles | ~3,700 |
| CORE Vertebrates | ~900 |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt | ~2 MB | ~900 | Vertebrate PFMs |
| JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt | ~500 KB | ~250 | Plant PFMs |
| JASPAR2024_CORE_non-redundant_pfms_jaspar.txt | ~3 MB | ~1,200 | All taxa PFMs |
| JASPAR2024_CORE_redundant_pfms_jaspar.zip | ~10 MB | ~3,700 | All versions |

### Previous Versions

| Version | Release | CORE Profiles | Status |
|---------|---------|---------------|--------|
| JASPAR 2022 | 2022-01-01 | ~1,000 | Archived |
| JASPAR 2020 | 2020-01-01 | ~750 | Archived |
| JASPAR 2018 | 2018-01-01 | ~590 | Archived |

---

## API Access

### Configuration

| Property | Value |
|----------|-------|
| Base URL | `https://jaspar.genereg.net/api/v1` |
| Authentication | None required |
| Rate Limit | 10 requests/second recommended |
| Response Format | JSON, PFM, MEME, TRANSFAC |

### API Endpoints

| Operation | Endpoint | Example |
|-----------|----------|---------|
| Get Matrix | `/matrix/{id}/` | `/matrix/MA0106.3/` |
| Get Format | `/matrix/{id}/` | `?format=meme` |
| Search | `/matrix/` | `?name=TP53&collection=CORE` |
| By Species | `/matrix/` | `?species=9606&tax_group=vertebrates` |
| Sites | `/sites/{id}/` | `/sites/MA0106.3/` |
| Collections | `/collections/` | List all collections |

---

## See Also

- [Schema Documentation](./schema.md)
- [JASPAR Database](https://jaspar.genereg.net/)
