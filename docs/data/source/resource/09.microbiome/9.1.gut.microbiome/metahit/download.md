---
id: download-metahit
title: "MetaHIT Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# MetaHIT Download Instructions

## Quick Start

```bash
# Download IGC gene catalog from BGI
wget http://meta.genomics.cn/meta/dataTools/download/IGC.fa.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for decompression
- **DIAMOND** or **BLAST** for searching (optional)
- 50-200GB storage for gene catalog
- 5-50TB for raw sequences

## No Registration Required

MetaHIT data is publicly available. Note: Project completed; data is legacy but foundational.

## Download Methods

### Method 1: IGC Gene Catalog (Primary Resource)

```bash
# Download Integrated Gene Catalog
# From BGI Meta server
wget http://meta.genomics.cn/meta/dataTools/download/IGC.fa.gz

# Alternative mirror (if available)
# wget ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/IGC.fa.gz

# Decompress
gunzip IGC.fa.gz

# Note: File is ~4GB compressed, ~20GB uncompressed
```

### Method 2: Gene Annotations

```bash
# Download gene annotations (functional, taxonomic)
wget http://meta.genomics.cn/meta/dataTools/download/IGC.annotation.txt.gz

# KEGG annotations
wget http://meta.genomics.cn/meta/dataTools/download/IGC.kegg.txt.gz

# eggNOG annotations
wget http://meta.genomics.cn/meta/dataTools/download/IGC.eggnog.txt.gz

# Taxonomic assignments
wget http://meta.genomics.cn/meta/dataTools/download/IGC.taxonomy.txt.gz
```

### Method 3: Raw Sequences from ENA

```bash
# MetaHIT raw sequences are deposited in ENA
# Project accession: PRJEB1220

# List available files
curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB1220&result=read_run&fields=run_accession,fastq_ftp" \
  > metahit_files.tsv

# Download specific samples
while IFS=$'\t' read -r run_accession ftp_url; do
  wget "$ftp_url"
done < metahit_files.tsv

# Using ENA browser
# https://www.ebi.ac.uk/ena/browser/view/PRJEB1220
```

### Method 4: Processed Abundance Data

```bash
# Sample abundance profiles (if available)
# Check publication supplementary materials

# Example from Nature 2010 paper
# wget [supplementary_data_url]/abundance_table.tsv

# Note: Abundance matrices are often in publication supplements
# Check DOI: 10.1038/nature08821 for supplementary files
```

### Method 5: Build Local Search Database

```bash
# Create DIAMOND database for fast searching
diamond makedb --in IGC.fa -d IGC

# Or create BLAST database
makeblastdb -in IGC.fa -dbtype nucl -out IGC

# Search your sequences against IGC
diamond blastx -d IGC -q query.fasta -o results.tsv
```

## File Inventory

### Gene Catalog Files

| File | Size | Description |
|------|------|-------------|
| IGC.fa.gz | ~4 GB | Gene sequences (FASTA) |
| IGC.annotation.txt.gz | ~1 GB | Combined annotations |
| IGC.kegg.txt.gz | ~500 MB | KEGG assignments |
| IGC.taxonomy.txt.gz | ~500 MB | Taxonomic assignments |

### Raw Sequence Data (ENA)

| Study | Samples | Size (est.) |
|-------|---------|-------------|
| PRJEB1220 | 124 | ~1 TB |
| PRJEB2054 | 396 | ~3 TB |
| PRJEB5224 | 292 | ~2 TB |

## Post-Download Processing

```bash
# Index gene catalog for fast access
python3 << 'EOF'
from Bio import SeqIO
import gzip

# Create gene index
gene_index = {}
with gzip.open("IGC.fa.gz", "rt") as f:
    for record in SeqIO.parse(f, "fasta"):
        gene_index[record.id] = {
            'length': len(record.seq),
            'description': record.description
        }

print(f"Indexed {len(gene_index)} genes")
EOF

# Parse annotations
python3 << 'EOF'
import pandas as pd
import gzip

# Load annotations
annotations = pd.read_csv(
    "IGC.annotation.txt.gz",
    sep="\t",
    compression="gzip"
)

print(f"Annotations: {len(annotations)} genes")
print(f"Columns: {list(annotations.columns)}")

# Summary statistics
print("\nTaxonomic distribution:")
print(annotations['phylum'].value_counts().head(10))
EOF

# Map reads to IGC
bowtie2-build IGC.fa IGC_index
bowtie2 -x IGC_index -U sample.fastq.gz -S sample.sam

# Calculate gene abundances
# Using samtools and custom script
samtools view -bS sample.sam | samtools sort > sample.sorted.bam
samtools idxstats sample.sorted.bam > gene_counts.tsv
```

### Create Subset Catalogs

```python
import gzip
from Bio import SeqIO

# Create human-gut-specific subset
def filter_catalog(input_fa, output_fa, taxonomy_filter):
    """Filter gene catalog by taxonomy."""
    count = 0
    with gzip.open(input_fa, "rt") as fin, gzip.open(output_fa, "wt") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            if taxonomy_filter in record.description:
                SeqIO.write(record, fout, "fasta")
                count += 1
    return count

# Filter for Bacteroides genes
n = filter_catalog("IGC.fa.gz", "IGC.Bacteroides.fa.gz", "Bacteroides")
print(f"Extracted {n} Bacteroides genes")
```

## Verification

```bash
# Count genes in catalog
zcat IGC.fa.gz | grep -c "^>"

# Expected: ~9.9 million genes

# Check annotation completeness
zcat IGC.annotation.txt.gz | wc -l

# Verify KEGG assignments
zcat IGC.kegg.txt.gz | awk -F'\t' '{print $2}' | sort | uniq -c | head

# Check sequence integrity
python3 << 'EOF'
from Bio import SeqIO
import gzip

valid = 0
invalid = 0
with gzip.open("IGC.fa.gz", "rt") as f:
    for record in SeqIO.parse(f, "fasta"):
        if len(record.seq) > 0 and set(str(record.seq)).issubset("ATGCN"):
            valid += 1
        else:
            invalid += 1

print(f"Valid sequences: {valid}")
print(f"Invalid sequences: {invalid}")
EOF
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| IGC updates | Legacy (no new updates) |
| Raw data | Archived in ENA |

Note: MetaHIT project is completed. IGC remains a reference resource but is not actively updated.

## Common Issues

- **Large file size**: IGC is ~20GB uncompressed
- **Server availability**: BGI server may be slow; use mirrors
- **Legacy format**: Older annotation formats may need parsing
- **Newer alternatives**: Consider HMP, GMrepo for current data
- **Version tracking**: Multiple IGC versions exist

## Integration with Current Resources

```bash
# Map IGC genes to current databases
# Cross-reference with UniProt
diamond blastp -d uniprot -q IGC_proteins.fa -o igc_uniprot.tsv

# Link to KEGG current
# Use KEGG API for updated pathway information

# Compare with HMP gene catalog
# Newer catalogs may supersede IGC for some purposes
```

## Related Resources

- [HMP](../../9.2.body.site.microbiomes/hmp/README.md) - US microbiome project
- [GMrepo](../gmrepo/README.md) - Curated repository
- [KEGG](../../../../04.pathways.networks/4.1.metabolic.pathways/kegg/) - Pathway annotations
