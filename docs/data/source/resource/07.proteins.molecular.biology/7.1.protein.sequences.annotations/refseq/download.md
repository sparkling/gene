---
id: download-refseq
title: "NCBI RefSeq Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# NCBI RefSeq Download Instructions

## Quick Start

```bash
# Download human protein sequences (RefSeq)
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gunzip** for decompression
- **Aspera** for high-speed transfers (optional)
- **E-utilities** account for API access (optional)
- 50GB-5TB storage depending on scope

## No Registration Required

RefSeq data is public domain (US Government work). API key recommended for heavy use.

## Download Methods

### Method 1: Species-Specific Downloads

```bash
# Human
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.gbff.gz

# Mouse
wget https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.faa.gz

# Common model organisms
for org in H_sapiens M_musculus D_melanogaster C_elegans S_cerevisiae; do
  wget "https://ftp.ncbi.nlm.nih.gov/refseq/${org}/mRNA_Prot/*.protein.faa.gz"
done
```

### Method 2: RefSeq Release (Complete)

```bash
# Download release catalog
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/

# Download specific division
# Options: archaea, bacteria, fungi, invertebrate, plant, protozoa, vertebrate_mammalian, vertebrate_other, viral
wget -r -np -nH --cut-dirs=2 \
  https://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/

# Download protein sequences only
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete.*.protein.faa.gz
```

### Method 3: RefSeq Select (Curated Subset)

```bash
# RefSeq Select for human (high-quality subset)
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/refseqgene_human.protein.faa.gz

# MANE transcripts (NCBI-Ensembl matched annotations)
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.3.summary.txt.gz
```

### Method 4: E-utilities API

```bash
# Search for specific proteins
esearch -db protein -query "TP53[gene] AND human[organism] AND refseq[filter]" | \
  efetch -format fasta > tp53_proteins.fasta

# Download by accession list
echo -e "NP_000537.3\nNP_001119584.1" > accessions.txt
cat accessions.txt | epost -db protein | efetch -format fasta > proteins.fasta

# Batch download with history
esearch -db protein -query "kinase[title] AND human[organism] AND refseq[filter]" | \
  efetch -format fasta > human_kinases.fasta

# Get GenPept format
efetch -db protein -id NP_000537.3 -format gp > tp53.gp
```

### Method 5: Aspera High-Speed Download

```bash
# Install Aspera Connect
# Download link: https://www.ibm.com/aspera/connect/

# High-speed transfer
ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  -k1 -T -l 300m \
  anonftp@ftp.ncbi.nlm.nih.gov:/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz ./
```

### Method 6: NCBI Datasets CLI

```bash
# Install datasets CLI
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
chmod +x datasets

# Download human RefSeq proteins
./datasets download gene taxon human --include protein --filename human_proteins.zip

# Download specific gene
./datasets download gene symbol TP53 --taxon human --include protein
```

## File Inventory

### Species-Specific Files

| File | Size | Description |
|------|------|-------------|
| human.protein.faa.gz | ~25 MB | Human proteins (FASTA) |
| human.rna.fna.gz | ~100 MB | Human transcripts (FASTA) |
| human.rna.gbff.gz | ~1 GB | Human transcripts (GenBank) |

### Release Files

| Directory | Size | Description |
|-----------|------|-------------|
| complete/ | ~500 GB | All RefSeq sequences |
| vertebrate_mammalian/ | ~100 GB | Mammalian sequences |
| bacteria/ | ~200 GB | Bacterial genomes |
| viral/ | ~10 GB | Viral sequences |

### Catalog Files

| File | Size | Description |
|------|------|-------------|
| release-catalog/ | ~500 MB | Release documentation |
| RefSeq-release*.txt | ~1 GB | Accession catalog |

## Post-Download Processing

```bash
# Decompress
gunzip human.protein.faa.gz

# Index FASTA
samtools faidx human.protein.faa

# Create BLAST database
makeblastdb -in human.protein.faa -dbtype prot -out human_refseq

# Extract specific proteins
grep -A1 "NP_000537" human.protein.faa > tp53.fasta

# Convert GenBank to FASTA
python3 << 'EOF'
from Bio import SeqIO
for record in SeqIO.parse("human.rna.gbff", "genbank"):
    SeqIO.write(record, "individual_genes/" + record.id + ".gb", "genbank")
EOF

# Parse feature annotations
python3 << 'EOF'
from Bio import SeqIO
for record in SeqIO.parse("human.protein.gpff", "genbank"):
    for feature in record.features:
        if feature.type == "Region":
            print(f"{record.id}\t{feature.qualifiers.get('region_name', [''])[0]}\t{feature.location}")
EOF
```

## Verification

```bash
# Check FASTA integrity
zcat human.protein.faa.gz | grep -c "^>"

# Verify MD5 checksums
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/md5checksums.txt
md5sum -c md5checksums.txt

# Check file format
head human.protein.faa
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| RefSeq Release | Bimonthly (6x/year) |
| Species updates | Continuous |
| MANE updates | Quarterly |

## Common Issues

- **Large files**: Use Aspera for faster downloads
- **Accession versions**: Always include version number (NP_000537.3 not NP_000537)
- **Deprecated GI**: Use accession.version instead
- **Redundancy**: WP_ proteins are non-redundant across prokaryotes
- **API limits**: Use API key for >3 requests/second

## RefSeq vs GenBank

| Aspect | RefSeq | GenBank |
|--------|--------|---------|
| Redundancy | Non-redundant | Contains duplicates |
| Curation | NCBI-curated | Submitter-provided |
| Prefix | NP_, NM_, NC_ | Various |
| Updates | Centrally managed | Submitter updates |

## Filtering Examples

```bash
# Get only curated proteins (NP_)
zcat complete.protein.faa.gz | awk '/^>NP_/{p=1} /^>/{if(!/^>NP_/)p=0} p' > curated_proteins.faa

# Filter by organism
zcat complete.protein.faa.gz | awk '/^>.*\[Homo sapiens\]/{p=1} /^>/{if(!/\[Homo sapiens\]/)p=0} p' > human.faa

# Extract proteins with specific domain
grep -B1 "kinase" human.protein.faa | grep "^>" > kinase_headers.txt
```

## Related Resources

- [UniProt](../uniprot/) - Cross-referenced protein database
- [NCBI Gene](../../../../01.genetics.genomics/1.2.gene.databases/ncbi.gene/) - Gene records
- [Ensembl](../../../../01.genetics.genomics/1.2.gene.databases/ensembl/) - Parallel annotation
