---
id: download-pgc
title: "Psychiatric Genomics Consortium (PGC) Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# Psychiatric Genomics Consortium (PGC) Download Instructions

## Quick Start

```bash
# Download schizophrenia GWAS summary statistics
wget https://pgc.unc.edu/for-researchers/download-results/pgc3_scz_wave3.2_public.tsv.gz

# Download depression summary statistics
wget https://pgc.unc.edu/for-researchers/download-results/PGC_MDD_2023.tsv.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **gzip** for decompression
- Approximately 10GB for full summary statistics
- **dbGaP application** for individual-level data

## Registration Requirements

| Data Type | Access Level | Requirements |
|-----------|--------------|--------------|
| Summary statistics | Open | None |
| LD scores | Open | None |
| PRS weights | Open | Citation agreement |
| Individual genotypes | Controlled | dbGaP approval |

## Download Methods

### Method 1: PGC Website Downloads

```bash
# Navigate to: https://pgc.unc.edu/for-researchers/download-results/

# Schizophrenia (SCZ3)
wget https://pgc.unc.edu/for-researchers/download-results/pgc3_scz_wave3.2_public.tsv.gz

# Major Depression (MDD3)
wget https://pgc.unc.edu/for-researchers/download-results/PGC_MDD_2023.tsv.gz

# Bipolar Disorder (BIP3)
wget https://pgc.unc.edu/for-researchers/download-results/pgc-bip2021-all.tsv.gz

# ADHD
wget https://pgc.unc.edu/for-researchers/download-results/adhd_eur_jun2017.gz

# Autism Spectrum Disorder
wget https://pgc.unc.edu/for-researchers/download-results/iPSYCH-PGC_ASD_Nov2017.gz

# Anorexia Nervosa
wget https://pgc.unc.edu/for-researchers/download-results/pgc.an2019.tsv.gz

# PTSD
wget https://pgc.unc.edu/for-researchers/download-results/pts_eur_freeze2_overall.results.gz
```

### Method 2: GWAS Catalog Links

```bash
# PGC studies are also in GWAS Catalog
# Search by study accession

# Schizophrenia (GCST90308669)
curl "https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/GCST90308669/associations" \
  -o scz_gwas_catalog.json

# Download via GWAS Catalog FTP
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90308669/
```

### Method 3: Cross-Disorder Analysis

```bash
# Cross-disorder summary statistics
wget https://pgc.unc.edu/for-researchers/download-results/pgc_cdg2019.tsv.gz

# Genetic correlation matrices
wget https://pgc.unc.edu/for-researchers/download-results/rg_matrix_pgc.tsv
```

### Method 4: LD Hub for Analysis

```bash
# LD Hub provides pre-computed LD scores
# https://ldsc.broadinstitute.org/

# Download LDSC-formatted files
# PGC provides files pre-formatted for LDSC analysis
```

### Method 5: PRS Weights

```bash
# Polygenic Risk Score weights
# Often provided as supplementary materials

# Example PRS weights format
# Download from specific study supplementary data
```

## File Inventory

### Major Study Files

| Study | File | Size | Samples |
|-------|------|------|---------|
| SCZ3 | pgc3_scz_wave3.2_public.tsv.gz | ~500 MB | 320,000 |
| MDD3 | PGC_MDD_2023.tsv.gz | ~800 MB | 1,200,000 |
| BIP3 | pgc-bip2021-all.tsv.gz | ~400 MB | 413,000 |
| ADHD | adhd_eur_jun2017.gz | ~200 MB | 225,000 |
| ASD | iPSYCH-PGC_ASD_Nov2017.gz | ~150 MB | 54,000 |
| AN | pgc.an2019.tsv.gz | ~300 MB | 72,500 |

### Summary Statistics Columns

| Column | Description |
|--------|-------------|
| SNP | rs ID |
| CHR | Chromosome |
| BP | Base pair position |
| A1 | Effect allele |
| A2 | Other allele |
| OR/BETA | Effect size |
| SE | Standard error |
| P | P-value |
| FRQ | Allele frequency |
| INFO | Imputation quality |

## Post-Download Processing

```bash
# Decompress and inspect
gunzip -k pgc3_scz_wave3.2_public.tsv.gz
head pgc3_scz_wave3.2_public.tsv

# Parse summary statistics
python3 << 'EOF'
import pandas as pd

# Read summary statistics
scz = pd.read_csv('pgc3_scz_wave3.2_public.tsv.gz', sep='\t', compression='gzip')

print(f"Total variants: {len(scz)}")
print(f"\nColumns: {list(scz.columns)}")

# Find genome-wide significant hits
sig = scz[scz['P'] < 5e-8]
print(f"\nGenome-wide significant SNPs: {len(sig)}")

# Save significant hits
sig.to_csv('scz_significant_snps.tsv', sep='\t', index=False)
EOF

# Manhattan plot preparation
python3 << 'EOF'
import pandas as pd
import numpy as np

scz = pd.read_csv('pgc3_scz_wave3.2_public.tsv.gz', sep='\t', compression='gzip')

# Filter for plotting (reduce points)
scz['log10p'] = -np.log10(scz['P'])

# Sample for visualization
scz_plot = scz.sample(n=min(500000, len(scz)))
scz_plot.to_csv('scz_manhattan_data.tsv', sep='\t', index=False)
print(f"Points for plotting: {len(scz_plot)}")
EOF

# LDSC format conversion
python3 << 'EOF'
import pandas as pd

scz = pd.read_csv('pgc3_scz_wave3.2_public.tsv.gz', sep='\t', compression='gzip')

# Convert to LDSC format
# Required columns: SNP, A1, A2, N, Z or BETA/SE
ldsc = scz[['SNP', 'A1', 'A2', 'P']].copy()
# Calculate Z from OR and SE if available

ldsc.to_csv('scz_ldsc.sumstats.gz', sep='\t', index=False, compression='gzip')
EOF

# Gene-based analysis preparation
python3 << 'EOF'
import pandas as pd

scz = pd.read_csv('pgc3_scz_wave3.2_public.tsv.gz', sep='\t', compression='gzip')

# Format for MAGMA
magma_input = scz[['SNP', 'CHR', 'BP', 'P']].copy()
magma_input.columns = ['SNP', 'CHR', 'POS', 'P']
magma_input.to_csv('scz_magma.tsv', sep='\t', index=False)
print("Prepared for MAGMA gene-based analysis")
EOF
```

## Verification

```bash
# Check file integrity
zcat pgc3_scz_wave3.2_public.tsv.gz | head -5

# Count variants
zcat pgc3_scz_wave3.2_public.tsv.gz | wc -l

# Check P-value range
zcat pgc3_scz_wave3.2_public.tsv.gz | awk -F'\t' 'NR>1 {print $NF}' | sort -g | head

# Verify chromosome coverage
zcat pgc3_scz_wave3.2_public.tsv.gz | cut -f2 | sort | uniq -c
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| SCZ3 | 2022-04 | ~500 MB | Current |
| MDD3 | 2023-05 | ~800 MB | Current |
| BIP3 | 2021-05 | ~400 MB | Current |
| ADHD3 | 2023-06 | ~300 MB | Current |

### Version Notes

PGC current releases statistics:
- Schizophrenia (SCZ3): 320,404 samples, 287 loci
- Depression (MDD3): 1.2M+ samples, 243 loci
- Bipolar (BIP3): 413,466 samples, 64 loci
- ADHD: 225,534 samples, 27 loci
- Cross-disorder: 8 psychiatric disorders analyzed jointly

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://pgc.unc.edu/for-researchers/download-results/` |
| Rate Limit | N/A (file downloads) |
| Auth Required | No (summary stats), Yes (individual data via dbGaP) |
| Documentation | https://pgc.unc.edu |

## Update Schedule

| Study Type | Frequency |
|------------|-----------|
| Major GWAS | As published (years) |
| Meta-analyses | Periodic |
| Cross-disorder | Periodic |

## Common Issues

- **Genome build**: Check if GRCh37 or GRCh38
- **Effect allele**: Confirm which allele A1 represents
- **Sample overlap**: Different studies may share samples
- **Population**: Most PGC studies are European ancestry
- **File formats**: Various formats (tsv, txt, gz)

## Major PGC Studies

| Disorder | Study | Year | Key Findings |
|----------|-------|------|--------------|
| Schizophrenia | SCZ3 | 2022 | 287 loci |
| Depression | MDD3 | 2019 | 102 loci |
| Bipolar | BIP3 | 2021 | 64 loci |
| ADHD | ADHD2 | 2023 | 27 loci |
| Autism | ASD | 2019 | 5 loci |
| Anorexia | AN | 2019 | 8 loci |
| OCD | OCD | 2018 | 1 locus |
| Tourette | TS | 2019 | 1 locus |

## Analysis Tools

| Tool | Purpose | Link |
|------|---------|------|
| LDSC | Genetic correlation | https://github.com/bulik/ldsc |
| MAGMA | Gene-based tests | https://ctg.cncr.nl/software/magma |
| PRSice | Polygenic scores | https://prsice.info/ |
| FUMA | Functional mapping | https://fuma.ctglab.nl/ |
| LocusZoom | Regional plots | http://locuszoom.org/ |

## Integration Examples

```bash
# Calculate genetic correlation with LDSC
# Requires LDSC installation

# Prepare sumstats
python3 munge_sumstats.py \
  --sumstats pgc3_scz_wave3.2_public.tsv.gz \
  --out scz_ldsc \
  --merge-alleles w_hm3.snplist

# Calculate heritability
python3 ldsc.py \
  --h2 scz_ldsc.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out scz_h2

# Genetic correlation
python3 ldsc.py \
  --rg scz_ldsc.sumstats.gz,mdd_ldsc.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out scz_mdd_rg
EOF

# Gene-based analysis with MAGMA
# magma --bfile g1000_eur \
#       --gene-annot genes.annot \
#       --pval scz_magma.tsv ncol=N \
#       --out scz_genes

# Cross-reference with drug targets
python3 << 'EOF'
import pandas as pd

# Load significant genes (from MAGMA output)
# genes = pd.read_csv('scz_genes.genes.out', sep='\t')

# Load drug target database
# druggable = pd.read_csv('druggable_genome.tsv', sep='\t')

# Find druggable SCZ genes
# overlap = genes.merge(druggable, on='GENE')
# print(f"Druggable SCZ genes: {len(overlap)}")
EOF
```

## Related Resources

- [GWAS Catalog](../../../../01.genetics.genomics/1.5.expression.regulation/gwas.catalog/) - Full GWAS database
- [Allen Brain Atlas](../allen.brain.atlas/) - Brain expression
- [SynGO](../syngo/) - Synaptic genes
- [Open Targets](../../3.3.disease.gene.associations/open.targets/) - Drug targets
