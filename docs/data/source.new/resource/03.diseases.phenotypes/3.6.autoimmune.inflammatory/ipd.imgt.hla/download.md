---
id: download-ipd-imgt-hla
title: "IPD-IMGT/HLA Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# IPD-IMGT/HLA Download Instructions

## Quick Start

```bash
# Download nucleotide sequences
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_nuc.fasta

# Download protein sequences
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_prot.fasta
```

## Prerequisites

- **wget** or **curl** for downloads
- **FASTA parser** (Biopython recommended)
- Approximately 500MB disk space

## No Registration Required

IPD-IMGT/HLA data is freely available with attribution.

## Download Methods

### Method 1: FTP Downloads - FASTA Sequences

```bash
# Nucleotide sequences
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_nuc.fasta

# Protein sequences
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_prot.fasta

# Gene-specific nucleotide sequences
for gene in A B C DRB1 DQB1 DPB1; do
    wget "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/${gene}_nuc.fasta"
done

# Gene-specific protein sequences
for gene in A B C DRB1 DQB1 DPB1; do
    wget "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/${gene}_prot.fasta"
done

# CDS (coding sequence) files
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_gen.fasta
```

### Method 2: XML/DAT Format

```bash
# Full database in XML
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.gz

# EMBL format (DAT)
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat.gz

# Alignments (MSF format)
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/A_nuc.txt
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/A_prot.txt
```

### Method 3: REST API

```bash
# Get allele information
curl "https://www.ebi.ac.uk/ipd/imgt/hla/api/allele/HLA-A*02:01:01:01" \
  -H "Accept: application/json" \
  -o hla_a0201.json

# Search alleles
curl "https://www.ebi.ac.uk/ipd/imgt/hla/api/search?query=HLA-A*02" \
  -H "Accept: application/json" \
  -o hla_a02_search.json

# Get allele sequence
curl "https://www.ebi.ac.uk/ipd/imgt/hla/api/allele/HLA-A*02:01:01:01/sequence" \
  -H "Accept: application/json" \
  -o hla_a0201_sequence.json

# Get allele groups
curl "https://www.ebi.ac.uk/ipd/imgt/hla/api/groups" \
  -H "Accept: application/json" \
  -o hla_groups.json
```

### Method 4: Alignment Files

```bash
# Pre-computed alignments for each gene
for gene in A B C DRB1 DQB1 DPB1 DQA1 DPA1; do
    # Nucleotide alignments
    wget "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/${gene}_nuc.txt"
    # Protein alignments
    wget "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/${gene}_prot.txt"
done

# Full genomic alignments
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/Alignments_Full.zip
```

### Method 5: Reference Data

```bash
# Allele list (current)
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allelelist.txt

# Allele list with version history
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allelelist_history.txt

# Deleted alleles
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Deleted_alleles.txt

# Nomenclature equivalents (two-field to full)
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/wmda/hla_nom.txt
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/wmda/hla_nom_g.txt
wget https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/wmda/hla_nom_p.txt
```

## File Inventory

### FASTA Files

| File | Size | Description |
|------|------|-------------|
| hla_nuc.fasta | ~100 MB | All nucleotide sequences |
| hla_prot.fasta | ~30 MB | All protein sequences |
| hla_gen.fasta | ~150 MB | Genomic sequences |
| {gene}_nuc.fasta | Variable | Gene-specific nucleotide |
| {gene}_prot.fasta | Variable | Gene-specific protein |

### Reference Files

| File | Size | Description |
|------|------|-------------|
| Allelelist.txt | ~5 MB | Current allele list |
| hla.xml.gz | ~200 MB | Full database XML |
| hla.dat.gz | ~100 MB | EMBL format |

### Alignment Files

| File | Description |
|------|-------------|
| {gene}_nuc.txt | Nucleotide alignment |
| {gene}_prot.txt | Protein alignment |
| Alignments_Full.zip | Complete set |

## Post-Download Processing

```bash
# Parse FASTA sequences
python3 << 'EOF'
from Bio import SeqIO
import pandas as pd

alleles = []
for record in SeqIO.parse('hla_prot.fasta', 'fasta'):
    # Parse header: >HLA:HLA00001 A*01:01:01:01 366 bp
    parts = record.description.split()
    allele_name = parts[1] if len(parts) > 1 else ''

    alleles.append({
        'id': record.id,
        'allele': allele_name,
        'length': len(record.seq),
        'gene': allele_name.split('*')[0] if '*' in allele_name else ''
    })

df = pd.DataFrame(alleles)
print(f"Total alleles: {len(df)}")
print(f"\nAlleles per gene:")
print(df['gene'].value_counts())

df.to_csv('hla_alleles.tsv', sep='\t', index=False)
EOF

# Extract sequences for specific gene
python3 << 'EOF'
from Bio import SeqIO

gene = 'A'
output_records = []

for record in SeqIO.parse('hla_prot.fasta', 'fasta'):
    if f'{gene}*' in record.description:
        output_records.append(record)

SeqIO.write(output_records, f'HLA-{gene}_prot.fasta', 'fasta')
print(f"HLA-{gene} alleles: {len(output_records)}")
EOF

# Parse allele nomenclature
python3 << 'EOF'
import pandas as pd

# Parse allele list
alleles = pd.read_csv('Allelelist.txt', sep=',', comment='#',
                      names=['allele_id', 'allele_name'])

print(f"Total alleles: {len(alleles)}")

# Parse into components
def parse_allele(name):
    if '*' not in name:
        return {}
    gene, rest = name.split('*')
    fields = rest.split(':')
    return {
        'gene': gene,
        'allele_group': fields[0] if len(fields) > 0 else '',
        'protein': fields[1] if len(fields) > 1 else '',
        'synonymous': fields[2] if len(fields) > 2 else '',
        'noncoding': fields[3] if len(fields) > 3 else ''
    }

parsed = pd.DataFrame([parse_allele(a) for a in alleles['allele_name']])
alleles = pd.concat([alleles, parsed], axis=1)
alleles.to_csv('hla_alleles_parsed.tsv', sep='\t', index=False)

# Count unique at each resolution
print(f"\nUnique genes: {alleles['gene'].nunique()}")
print(f"Unique allele groups: {len(alleles.groupby(['gene', 'allele_group']))}")
EOF

# Parse alignment files
python3 << 'EOF'
def parse_alignment(filename):
    """Parse IMGT alignment format"""
    sequences = {}
    current_pos = 0

    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith(' '):
                continue
            if line.startswith('*'):
                continue

            parts = line.split()
            if len(parts) >= 2:
                allele = parts[0]
                seq = ''.join(parts[1:])
                if allele not in sequences:
                    sequences[allele] = ''
                sequences[allele] += seq

    return sequences

# Parse HLA-A alignment
seqs = parse_alignment('A_prot.txt')
print(f"Aligned alleles: {len(seqs)}")
print(f"Alignment length: {len(list(seqs.values())[0]) if seqs else 0}")
EOF
```

## Verification

```bash
# Check FASTA format
head -20 hla_prot.fasta

# Count sequences
grep "^>" hla_prot.fasta | wc -l

# Check allele list
head -20 Allelelist.txt

# Verify specific allele
grep "A\*02:01:01:01" hla_prot.fasta

# Check gene distribution
grep "^>" hla_prot.fasta | cut -d'*' -f1 | cut -d' ' -f2 | sort | uniq -c | sort -rn
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database releases | Quarterly |
| New alleles | Continuous submissions |
| Nomenclature updates | As needed |

## Common Issues

- **Allele naming**: Follow WHO nomenclature (HLA-A*02:01:01:01)
- **Resolution levels**: Two-field vs four-field naming
- **Null alleles**: Denoted by 'N' suffix (e.g., A*02:01:01:01N)
- **Ambiguity codes**: G and P groups for typing ambiguity
- **Deleted alleles**: Check Deleted_alleles.txt for removed entries

## HLA Nomenclature

| Field | Example | Meaning |
|-------|---------|---------|
| Gene | HLA-A | Locus |
| Allele group | *02 | Serological equivalent |
| Specific protein | :01 | Protein sequence |
| Synonymous | :01 | Coding synonymous |
| Non-coding | :01 | Non-coding differences |
| Suffix | N, L, S, Q | Expression level |

### Expression Suffixes

| Suffix | Meaning |
|--------|---------|
| N | Null (not expressed) |
| L | Low expression |
| S | Secreted only |
| Q | Questionable expression |

## Integration Examples

```bash
# Map HLA to drug hypersensitivity (PharmGKB)
python3 << 'EOF'
# Example: HLA-B*57:01 and Abacavir
# HLA-B*15:02 and Carbamazepine

important_alleles = {
    'HLA-B*57:01': 'Abacavir hypersensitivity',
    'HLA-B*15:02': 'Carbamazepine SJS/TEN',
    'HLA-B*58:01': 'Allopurinol SCAR',
    'HLA-A*31:01': 'Carbamazepine hypersensitivity',
    'HLA-B*35:01': 'Nevirapine hypersensitivity'
}

# Extract these from FASTA
from Bio import SeqIO

for record in SeqIO.parse('hla_prot.fasta', 'fasta'):
    for allele, drug in important_alleles.items():
        if allele.replace('HLA-', '') in record.description:
            print(f"{allele}: {len(record.seq)} aa - {drug}")
            break
EOF

# Create typing reference panel
python3 << 'EOF'
from Bio import SeqIO

# Select common alleles for reference panel
common_alleles = ['A*01:01', 'A*02:01', 'A*03:01', 'A*24:02',
                  'B*07:02', 'B*08:01', 'B*44:02', 'B*51:01']

reference = []
for record in SeqIO.parse('hla_prot.fasta', 'fasta'):
    for allele in common_alleles:
        if f'{allele}:' in record.description or record.description.endswith(allele):
            reference.append(record)
            break

SeqIO.write(reference, 'hla_reference_panel.fasta', 'fasta')
print(f"Reference panel: {len(reference)} alleles")
EOF
```

## Related Resources

- [ImmunoBase](../immunobase/) - Autoimmune genetics
- [PharmGKB](../../../../04.drugs.compounds/4.3.pharmacogenomics/pharmgkb/) - HLA pharmacogenomics
- [IEDB](../../../../02.proteins.structures/2.4.immune.proteomics/iedb/) - Immune epitope data
