---
id: download-uniprot-id-mapping
title: "UniProt ID Mapping Service Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
---

# UniProt ID Mapping Service - Download Documentation

## Overview

The UniProt ID Mapping service converts between UniProt accessions and 100+ external database identifiers. It supports both interactive web use and programmatic REST API access.

## REST API (New)

### Base URL

```
https://rest.uniprot.org/idmapping
```

### Step 1: Submit Job

```bash
# Submit mapping request
curl --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'ids=P04637,P53_HUMAN,Q9Y6K9' \
  --form 'from=UniProtKB_AC-ID' \
  --form 'to=Ensembl'
```

Response:
```json
{
  "jobId": "abc123def456"
}
```

### Step 2: Check Status

```bash
curl 'https://rest.uniprot.org/idmapping/status/abc123def456'
```

Response (while running):
```json
{
  "jobStatus": "RUNNING"
}
```

Response (complete):
```json
{
  "jobStatus": "FINISHED",
  "results": "https://rest.uniprot.org/idmapping/results/abc123def456"
}
```

### Step 3: Retrieve Results

```bash
# Get results (JSON)
curl 'https://rest.uniprot.org/idmapping/results/abc123def456'

# Get results (TSV)
curl 'https://rest.uniprot.org/idmapping/results/abc123def456' \
  -H 'Accept: text/plain'
```

## Database Codes (From/To)

### Gene Databases

| Code | Database |
|------|----------|
| Ensembl | Ensembl gene |
| Ensembl_Genomes | Ensembl Genomes |
| GeneID | NCBI Gene |
| KEGG | KEGG gene |
| HGNC | HGNC ID |
| MGI | Mouse Genome Informatics |
| RGD | Rat Genome Database |
| SGD | Saccharomyces Genome Database |
| WormBase | WormBase |
| FlyBase | FlyBase |
| ZFIN | Zebrafish Information Network |

### Protein Databases

| Code | Database |
|------|----------|
| UniProtKB_AC-ID | UniProt accession or ID |
| UniProtKB | UniProt with entry info |
| UniRef100 | UniRef100 cluster |
| UniRef90 | UniRef90 cluster |
| UniRef50 | UniRef50 cluster |
| UniParc | UniParc ID |
| RefSeq_Protein | RefSeq protein |
| PDB | Protein Data Bank |
| AlphaFoldDB | AlphaFold DB |

### Sequence Databases

| Code | Database |
|------|----------|
| EMBL | EMBL nucleotide |
| EMBL-GenBank-DDBJ | Combined nucleotide |
| RefSeq_Nucleotide | RefSeq nucleotide |
| PIR | Protein Information Resource |

### Pathway/Ontology Databases

| Code | Database |
|------|----------|
| Reactome | Reactome pathway |
| KEGG | KEGG pathway |
| GO | Gene Ontology |
| InterPro | InterPro domain |
| Pfam | Pfam family |

### Other Databases

| Code | Database |
|------|----------|
| ChEMBL | ChEMBL target |
| DrugBank | DrugBank target |
| GlyConnect | Glycoprotein |
| OpenTargets | Open Targets |
| Proteomes | Proteome ID |

## Mapping Examples

### UniProt to Gene IDs

```bash
# UniProt to Ensembl Gene
curl --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'ids=P04637,P38398' \
  --form 'from=UniProtKB_AC-ID' \
  --form 'to=Ensembl'

# UniProt to NCBI Gene
curl --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'ids=P04637,P38398' \
  --form 'from=UniProtKB_AC-ID' \
  --form 'to=GeneID'
```

### Gene IDs to UniProt

```bash
# Ensembl to UniProt
curl --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'ids=ENSG00000141510,ENSG00000012048' \
  --form 'from=Ensembl' \
  --form 'to=UniProtKB'

# NCBI Gene to UniProt
curl --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'ids=7157,672' \
  --form 'from=GeneID' \
  --form 'to=UniProtKB'
```

### UniProt to Structure/Pathway

```bash
# UniProt to PDB
curl --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'ids=P04637' \
  --form 'from=UniProtKB_AC-ID' \
  --form 'to=PDB'

# UniProt to Reactome
curl --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'ids=P04637' \
  --form 'from=UniProtKB_AC-ID' \
  --form 'to=Reactome'
```

## Output Formats

| Format | Accept Header | Description |
|--------|---------------|-------------|
| JSON | application/json | Default |
| TSV | text/plain | Tab-separated |
| FASTA | text/fasta | Sequences |
| XML | application/xml | Full records |

## Streaming Large Results

```bash
# Stream results with pagination
curl 'https://rest.uniprot.org/idmapping/results/stream/abc123def456?format=tsv'
```

## Python Examples

### Asynchronous Mapping

```python
import requests
import time

def map_ids(ids, from_db, to_db):
    """Map identifiers between databases."""
    api_url = "https://rest.uniprot.org/idmapping"

    # Submit job
    response = requests.post(
        f"{api_url}/run",
        data={
            "ids": ",".join(ids),
            "from": from_db,
            "to": to_db
        }
    )
    job_id = response.json()["jobId"]

    # Poll for completion
    while True:
        status = requests.get(f"{api_url}/status/{job_id}").json()

        if status["jobStatus"] == "FINISHED":
            break
        elif status["jobStatus"] == "ERROR":
            raise Exception("Mapping job failed")

        time.sleep(1)

    # Get results
    results = requests.get(f"{api_url}/results/{job_id}").json()

    return {r["from"]: r["to"] for r in results["results"]}


# Example usage
mapping = map_ids(
    ids=["P04637", "P38398", "P42336"],
    from_db="UniProtKB_AC-ID",
    to_db="Ensembl"
)
print(mapping)
```

### Batch Processing

```python
def batch_map_ids(ids, from_db, to_db, batch_size=10000):
    """Map large lists of IDs in batches."""
    all_results = {}

    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        results = map_ids(batch, from_db, to_db)
        all_results.update(results)
        print(f"Processed {min(i+batch_size, len(ids))}/{len(ids)}")

    return all_results
```

## Rate Limits

| Access Type | Limit |
|-------------|-------|
| Default | 100 requests/second |
| Large batches | Up to 100,000 IDs per job |

## Web Interface

For interactive use: https://www.uniprot.org/id-mapping

1. Paste or upload IDs
2. Select source database
3. Select target database
4. Submit and download results

## Pre-computed Mappings (Bulk)

Download pre-computed mapping files:

```bash
# ID mapping file
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

# Selected mappings
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
```

### idmapping_selected.tab columns

| Column | Content |
|--------|---------|
| 1 | UniProtKB-AC |
| 2 | UniProtKB-ID |
| 3 | GeneID (EntrezGene) |
| 4 | RefSeq |
| 5 | GI |
| 6 | PDB |
| 7 | GO |
| 8 | UniRef100 |
| 9 | UniRef90 |
| 10 | UniRef50 |
| 11 | UniParc |
| 12 | PIR |
| 13 | NCBI-taxon |
| 14 | MIM |
| 15 | UniGene |
| 16 | PubMed |
| 17 | EMBL |
| 18 | EMBL-CDS |
| 19 | Ensembl |
| 20 | Ensembl_TRS |
| 21 | Ensembl_PRO |
| 22 | Additional PubMed |

## See Also

- [Schema Documentation](./schema.md)
- [UniProt ID Mapping Help](https://www.uniprot.org/help/id_mapping)
