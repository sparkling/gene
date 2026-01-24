---
id: download-uniprot
title: "UniProt Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# UniProt Download Instructions

## Quick Start

```bash
# Download UniProt Swiss-Prot (curated, ~570K entries)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for decompression
- **rsync** for efficient mirroring
- 50-500GB disk space depending on data scope

## No Registration Required

UniProt data is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: FTP Complete Downloads

```bash
# Swiss-Prot (curated)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

# TrEMBL (automated annotation, ~250M entries)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

# Combined (Swiss-Prot + TrEMBL)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz
```

### Method 2: Taxonomic Subsets

```bash
# Human proteome
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz

# Reference proteomes (all species)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README

# Archaea, Bacteria, Eukaryota, Viruses
rsync -av rsync://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/ ./eukaryota/
```

### Method 3: ID Mapping Files

```bash
# UniProt ID mapping (all external IDs)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

# Selected ID types
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz

# Specific organism mapping (human)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
```

### Method 4: REST API Bulk Download

```bash
# Download human reviewed proteins (Swiss-Prot)
curl -o human_swissprot.fasta.gz \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(organism_id:9606)%20AND%20(reviewed:true)&compressed=true"

# Download with specific fields (TSV)
curl -o human_proteins.tsv.gz \
  "https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=(organism_id:9606)%20AND%20(reviewed:true)&fields=accession,id,gene_names,protein_name,length,go_p,go_c,go_f&compressed=true"

# Download specific proteins by accession
curl -o specific_proteins.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=accession:(P04637+P53_HUMAN)"
```

### Method 5: UniRef Clusters

```bash
# UniRef100 (100% identity clusters)
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz

# UniRef90 (90% identity)
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz

# UniRef50 (50% identity)
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
```

### Method 6: rsync Mirror

```bash
# Mirror Swiss-Prot
rsync -av rsync://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot* ./

# Mirror entire release
rsync -av rsync://ftp.uniprot.org/pub/databases/uniprot/current_release/ ./uniprot_current/
```

## File Inventory

### Knowledgebase Files

| File | Size | Description |
|------|------|-------------|
| uniprot_sprot.fasta.gz | ~90 MB | Swiss-Prot sequences |
| uniprot_sprot.xml.gz | ~1 GB | Swiss-Prot full XML |
| uniprot_sprot.dat.gz | ~600 MB | Swiss-Prot flat file |
| uniprot_trembl.fasta.gz | ~50 GB | TrEMBL sequences |

### ID Mapping Files

| File | Size | Description |
|------|------|-------------|
| idmapping.dat.gz | ~15 GB | All ID mappings |
| idmapping_selected.tab.gz | ~5 GB | Selected mappings |
| HUMAN_9606_idmapping.dat.gz | ~200 MB | Human mappings |

### Reference Proteomes

| File | Size | Description |
|------|------|-------------|
| UP000005640_9606.fasta.gz | ~25 MB | Human reference proteome |
| UP000000589_10090.fasta.gz | ~20 MB | Mouse reference proteome |

### UniRef Files

| File | Size | Description |
|------|------|-------------|
| uniref100.fasta.gz | ~70 GB | 100% identity clusters |
| uniref90.fasta.gz | ~30 GB | 90% identity clusters |
| uniref50.fasta.gz | ~10 GB | 50% identity clusters |

## Post-Download Processing

```bash
# Decompress
gunzip uniprot_sprot.fasta.gz

# Index FASTA for quick access
samtools faidx uniprot_sprot.fasta

# Create BLAST database
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot

# Parse ID mapping file
zcat idmapping_selected.tab.gz | \
  awk -F'\t' '$2=="GeneID" {print $1"\t"$3}' > uniprot_to_ncbi.tsv

# Extract human proteins
zcat uniprot_sprot.fasta.gz | \
  awk '/^>/ {p = /OS=Homo sapiens/} p' > human_proteins.fasta

# Convert to TSV
python3 << 'EOF'
from Bio import SeqIO
import gzip

with gzip.open('uniprot_sprot.fasta.gz', 'rt') as f:
    for record in SeqIO.parse(f, 'fasta'):
        acc = record.id.split('|')[1]
        print(f"{acc}\t{len(record.seq)}\t{record.description}")
EOF
```

## Verification

```bash
# Check FASTA integrity
zcat uniprot_sprot.fasta.gz | grep -c "^>"

# Verify release version
head uniprot_sprot.fasta

# Check ID mapping format
zcat idmapping_selected.tab.gz | head -5
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Swiss-Prot | Every 4 weeks |
| TrEMBL | Every 4 weeks |
| Reference proteomes | Quarterly |
| ID mapping | With each release |

## Common Issues

- **Large file downloads**: Use rsync for resume capability
- **Memory issues**: Use streaming parsers for TrEMBL
- **ID versioning**: UniProt accessions are stable; entry versions change
- **Isoforms**: Use varsplic files for isoform sequences
- **Deprecated accessions**: Check secondary accessions in XML/dat files

## File Formats Comparison

| Format | Best For |
|--------|----------|
| FASTA | Sequence similarity search |
| XML | Complete annotation parsing |
| TSV | Quick tabular analysis |
| DAT | Legacy pipelines |
| RDF | Semantic web applications |

## REST API Examples

```bash
# Search and download
curl "https://rest.uniprot.org/uniprotkb/search?query=kinase+AND+organism_id:9606&format=json&size=100"

# Get single entry
curl "https://rest.uniprot.org/uniprotkb/P04637.json"

# ID mapping service
curl -X POST "https://rest.uniprot.org/idmapping/run" \
  --data "from=UniProtKB_AC-ID&to=GeneID&ids=P04637,P01308"
```

## Related Resources

- [RefSeq](../refseq/_index.md) - NCBI protein sequences
- [PDB](../../7.2.protein.structures/pdb/) - Protein structures
- [InterPro](../../7.3.molecular.interactions/) - Protein families
