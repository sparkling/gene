---
id: download-homd
title: "HOMD Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# HOMD (Human Oral Microbiome Database) Download Instructions

## Quick Start

```bash
# Download all 16S reference sequences
wget https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq_V15.23.fasta
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for decompression
- **BLAST** or **VSEARCH** for sequence searching
- 5-50GB storage depending on data scope

## No Registration Required

HOMD data is freely available for academic use.

## Download Methods

### Method 1: 16S rRNA Reference Sequences

```bash
# Download 16S reference database (main resource)
wget https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq_V15.23.fasta

# Download with taxonomy
wget https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq_V15.23.taxonomy.txt

# Download QIIME-formatted database
wget https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq_V15.23.qiime.fasta
wget https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq_V15.23.qiime.taxonomy
```

### Method 2: Genome Sequences

```bash
# Download all annotated genomes
wget -r -np -nH --cut-dirs=2 \
  https://www.homd.org/ftp/genomes/PROKKA/

# Download specific genome by HOMT ID
wget https://www.homd.org/ftp/genomes/PROKKA/HMT-096/

# Download genome sequence (FASTA)
wget https://www.homd.org/ftp/genomes/all/GCF_000007465.2_ASM746v2/GCF_000007465.2_ASM746v2_genomic.fna.gz

# Download annotations (GFF)
wget https://www.homd.org/ftp/genomes/all/GCF_000007465.2_ASM746v2/GCF_000007465.2_ASM746v2_genomic.gff.gz
```

### Method 3: Taxonomic Data

```bash
# Download complete taxonomy table
wget https://www.homd.org/ftp/taxonomy/HOMD_taxonomy.txt

# Download NCBI-linked taxonomy
wget https://www.homd.org/ftp/taxonomy/HOMD_NCBI_taxonomy.txt

# Download phylogenetic tree
wget https://www.homd.org/ftp/trees/HOMD_16S_tree.nwk
```

### Method 4: Web Interface Export

1. Visit https://www.homd.org
2. Use Taxon Table browser
3. Select taxa of interest
4. Export to CSV/Excel

### Method 5: BLAST Database

```bash
# Download pre-built BLAST database (if available)
wget https://www.homd.org/ftp/BLAST/HOMD_16S.nhr
wget https://www.homd.org/ftp/BLAST/HOMD_16S.nin
wget https://www.homd.org/ftp/BLAST/HOMD_16S.nsq

# Or build your own
makeblastdb -in HOMD_16S_rRNA_RefSeq_V15.23.fasta -dbtype nucl -out HOMD_16S
```

### Method 6: FTP Bulk Download

```bash
# List available files
curl -l ftp://ftp.homd.org/

# Download entire FTP site
wget -r -np ftp://ftp.homd.org/

# Selective download
wget -r -np -A "*.fasta,*.gff,*.txt" ftp://ftp.homd.org/
```

## File Inventory

### 16S Reference Database

| File | Size | Description |
|------|------|-------------|
| HOMD_16S_rRNA_RefSeq_V15.23.fasta | ~15 MB | All 16S sequences |
| HOMD_16S_rRNA_RefSeq_V15.23.taxonomy.txt | ~500 KB | Taxonomy mapping |
| HOMD_16S_rRNA_RefSeq_V15.23.qiime.fasta | ~15 MB | QIIME format |

### Genome Data

| Category | Size (est.) | Description |
|----------|-------------|-------------|
| All genomes | ~10 GB | Complete + draft |
| Complete genomes | ~2 GB | Finished assemblies |
| Annotations (GFF) | ~3 GB | Gene annotations |
| Protein sequences | ~4 GB | Predicted proteins |

### Taxonomy Files

| File | Size | Description |
|------|------|-------------|
| HOMD_taxonomy.txt | ~200 KB | Full taxonomy table |
| HOMD_phenotype.txt | ~100 KB | Phenotypic data |
| HOMD_tree.nwk | ~500 KB | Phylogenetic tree |

## Post-Download Processing

```bash
# Create BLAST database
makeblastdb -in HOMD_16S_rRNA_RefSeq_V15.23.fasta \
  -dbtype nucl -out HOMD_16S -title "HOMD 16S"

# Create VSEARCH database (for DADA2/QIIME)
vsearch --makeudb_usearch HOMD_16S_rRNA_RefSeq_V15.23.fasta \
  --output HOMD_16S.udb

# Parse taxonomy
python3 << 'EOF'
import pandas as pd

# Load taxonomy
tax = pd.read_csv("HOMD_16S_rRNA_RefSeq_V15.23.taxonomy.txt",
                  sep="\t", header=None,
                  names=["seq_id", "taxonomy"])

# Parse into levels
tax[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = \
    tax['taxonomy'].str.split(';', expand=True)

# Save parsed taxonomy
tax.to_csv("HOMD_taxonomy_parsed.tsv", sep="\t", index=False)

# Summary
print(f"Total sequences: {len(tax)}")
print(f"Unique genera: {tax['Genus'].nunique()}")
print(f"Unique species: {tax['Species'].nunique()}")
EOF

# Format for QIIME2
python3 << 'EOF'
from Bio import SeqIO

# Convert to QIIME2 format
with open("HOMD_16S_rRNA_RefSeq_V15.23.fasta") as fin:
    with open("homd_seqs.fasta", "w") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            # Extract HOMT ID from header
            homt_id = record.description.split()[0]
            record.id = homt_id
            record.description = ""
            SeqIO.write(record, fout, "fasta")
EOF
```

### Create Oral-Specific Classifier

```bash
# Train QIIME2 classifier on HOMD
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path HOMD_16S_rRNA_RefSeq_V15.23.qiime.fasta \
  --output-path homd-seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path HOMD_16S_rRNA_RefSeq_V15.23.qiime.taxonomy \
  --output-path homd-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads homd-seqs.qza \
  --i-reference-taxonomy homd-taxonomy.qza \
  --o-classifier homd-classifier.qza
```

## Verification

```bash
# Count sequences
grep -c "^>" HOMD_16S_rRNA_RefSeq_V15.23.fasta

# Check FASTA format
head -10 HOMD_16S_rRNA_RefSeq_V15.23.fasta

# Verify taxonomy
head HOMD_16S_rRNA_RefSeq_V15.23.taxonomy.txt

# Test BLAST database
blastn -db HOMD_16S -query test.fasta -outfmt 6 | head
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| 16S database | Quarterly |
| New genomes | Continuous |
| Taxonomy updates | Semi-annual |

## Common Issues

- **Naming conventions**: HOMT IDs may update between versions
- **NCBI sync**: Some NCBI taxonomy may differ
- **Provisional taxa**: May lack full phenotype data
- **Sequence variants**: Multiple 16S copies per species
- **URL changes**: FTP paths may change between versions

## Integration with Oral Studies

```bash
# Classify oral microbiome samples
vsearch --usearch_global sample_16S.fasta \
  --db HOMD_16S_rRNA_RefSeq_V15.23.fasta \
  --id 0.97 \
  --blast6out classification.tsv

# Compare to whole-mouth reference
# Link HMT IDs to specific oral sites
python3 << 'EOF'
import pandas as pd

# Map classifications to oral sites
oral_sites = pd.read_csv("HOMD_oral_sites.txt", sep="\t")
classifications = pd.read_csv("classification.tsv", sep="\t",
                              names=["query", "target", "identity", ...])

# Add site information
merged = classifications.merge(oral_sites, left_on="target", right_on="homt_id")
print(merged.groupby("oral_site").size())
EOF
```

## Related Resources

- [HMP](../hmp/README.md) - Multi-site microbiome data
- [mBodyMap](../mbodymap/README.md) - Body-wide atlas
- [GMrepo](../../9.1.gut.microbiome/gmrepo/README.md) - Gut microbiome data
