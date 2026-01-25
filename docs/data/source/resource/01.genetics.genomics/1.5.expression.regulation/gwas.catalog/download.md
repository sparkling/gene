---
id: download-gwas-catalog
title: "GWAS Catalog Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# GWAS Catalog Download Instructions

## Quick Start

```bash
# Download all associations
wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -O gwas_catalog_associations.tsv
```

## Prerequisites

- **wget** or **curl** for downloads
- Approximately 500MB-2GB disk space

## No Registration Required

GWAS Catalog data is freely available under EMBL-EBI terms of use.

## Download Methods

### Method 1: Web Downloads (Recommended)

```bash
# All associations (main file)
wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative \
  -O gwas_catalog_associations.tsv

# All studies
wget https://www.ebi.ac.uk/gwas/api/search/downloads/studies \
  -O gwas_catalog_studies.tsv

# All ancestry data
wget https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry \
  -O gwas_catalog_ancestry.tsv
```

### Method 2: FTP Downloads

```bash
# Full release (TSV)
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv

# With mapped traits
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv

# Studies with publications
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-studies_ontology-annotated.tsv

# Ancestry information
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv
```

### Method 3: Harmonized Summary Statistics

```bash
# List available studies with summary stats
curl "https://www.ebi.ac.uk/gwas/summary-statistics/api/studies" | \
  jq '.["_embedded"]["studies"][] | {accession, trait}'

# Download specific study summary stats
STUDY_ID="GCST000001"
wget "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/${STUDY_ID}/" -r -np
```

### Method 4: REST API

```bash
# Get associations for a trait
curl "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/EFO_0000384/associations?size=100" \
  -o breast_cancer_associations.json

# Get SNP information
curl "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs7329174" \
  -o rs7329174.json

# Get study details
curl "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000001" \
  -o study.json

# Search by gene
curl "https://www.ebi.ac.uk/gwas/rest/api/associations/search/findByGene?geneName=BRCA1" \
  -o brca1_associations.json
```

### Method 5: Specific Data Subsets

```bash
# Download by chromosome
for chr in {1..22} X Y; do
  wget "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative?chromosome=${chr}" \
    -O "gwas_chr${chr}.tsv"
done

# Download by trait
curl "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative?efoTrait=EFO_0000384" \
  -o breast_cancer_gwas.tsv
```

### Method 6: Knowledge Graph Format

```bash
# Download knowledge graph triples
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-kg/gwas-catalog-kg.ttl.gz

# Download SPARQL endpoint data
# Access via: https://www.ebi.ac.uk/gwas/sparql
```

## File Inventory

### Main Catalog Files

| File | Size | Description |
|------|------|-------------|
| gwas-catalog-associations.tsv | ~200 MB | All associations |
| gwas-catalog-associations_ontology-annotated.tsv | ~250 MB | With EFO mapping |
| gwas-catalog-studies.tsv | ~50 MB | All studies |
| gwas-catalog-ancestry.tsv | ~20 MB | Ancestry data |

### Summary Statistics

| Data Type | Size | Description |
|-----------|------|-------------|
| Per-study summary stats | 10MB-10GB each | Full GWAS results |
| Harmonized format | Variable | Standardized coordinates |

### Supplementary Files

| File | Size | Description |
|------|------|-------------|
| gwas-diagram | ~10 MB | SNP-trait diagrams |
| gwas-kg.ttl.gz | ~100 MB | Knowledge graph |

## Post-Download Processing

```bash
# Parse associations file
python3 << 'EOF'
import pandas as pd

# Load associations
assoc = pd.read_csv('gwas_catalog_associations.tsv', sep='\t', low_memory=False)
print(f"Total associations: {len(assoc)}")
print(assoc.columns.tolist())

# Filter genome-wide significant
gw_sig = assoc[assoc['P-VALUE'] < 5e-8]
print(f"Genome-wide significant: {len(gw_sig)}")

# Group by trait
trait_counts = assoc['DISEASE/TRAIT'].value_counts().head(20)
print(trait_counts)
EOF

# Extract specific SNPs
grep -E "(rs7329174|rs1234567)" gwas_catalog_associations.tsv > selected_snps.tsv

# Create BED file
python3 << 'EOF'
import pandas as pd

assoc = pd.read_csv('gwas_catalog_associations.tsv', sep='\t', low_memory=False)

# Filter valid coordinates
bed_data = assoc[['CHR_ID', 'CHR_POS', 'SNPS', 'P-VALUE']].dropna()
bed_data['CHR_ID'] = 'chr' + bed_data['CHR_ID'].astype(str)
bed_data['START'] = bed_data['CHR_POS'].astype(int) - 1
bed_data['END'] = bed_data['CHR_POS'].astype(int)

bed_data[['CHR_ID', 'START', 'END', 'SNPS', 'P-VALUE']].to_csv(
    'gwas_catalog.bed', sep='\t', index=False, header=False
)
EOF

# Link to rsIDs
awk -F'\t' 'NR==1 || $22 ~ /rs[0-9]+/' gwas_catalog_associations.tsv > gwas_with_rsids.tsv
```

## Verification

```bash
# Check file structure
head -5 gwas_catalog_associations.tsv

# Count associations
wc -l gwas_catalog_associations.tsv

# Check for specific trait
grep -i "diabetes" gwas_catalog_associations.tsv | wc -l

# Verify column names
head -1 gwas_catalog_associations.tsv | tr '\t' '\n' | nl
```

## Dataset Versions

### Current Release

| Property | Value |
|----------|-------|
| Version | 2026-01 (Weekly) |
| Release Date | 2026-01-20 |
| Total Size | ~500 MB |
| Studies | 7,000+ |
| Associations | 500,000+ |
| Publications | 6,000+ |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| gwas-catalog-associations.tsv | ~200 MB | 500K+ | All associations |
| gwas-catalog-associations_ontology.tsv | ~250 MB | 500K+ | With EFO mapping |
| gwas-catalog-studies.tsv | ~50 MB | 7K+ | All studies |
| gwas-catalog-ancestry.tsv | ~20 MB | 15K+ | Ancestry data |

### Summary Statistics

| Category | Studies | Size Range |
|----------|---------|------------|
| Full summary stats | 5,000+ | 10 MB - 10 GB |
| Harmonized | 4,000+ | Standardized |
| Available traits | 3,000+ | Various diseases |

### Trait Categories

| Category | Associations | Description |
|----------|--------------|-------------|
| Cancer | 50,000+ | Various cancer types |
| Cardiovascular | 40,000+ | Heart, blood pressure |
| Metabolic | 35,000+ | Diabetes, obesity |
| Neurological | 30,000+ | Brain disorders |
| Autoimmune | 25,000+ | Immune conditions |

---

## API Access

| Property | Value |
|----------|-------|
| Base URL | https://www.ebi.ac.uk/gwas/rest/api/ |
| Rate Limit | 100 req/min |
| Auth Required | No |
| Response Format | JSON |

---

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Catalog updates | Weekly |
| Summary statistics | As submitted |
| Major releases | Quarterly |

## Common Issues

- **Multiple SNPs per row**: Some associations report multiple linked SNPs
- **P-value format**: May be in scientific notation; parse carefully
- **Coordinate versions**: Check genome build (GRCh38 standard)
- **Trait mapping**: Use ontology-annotated file for standardized traits
- **Missing data**: Not all fields populated; handle NaN appropriately

## Key Column Descriptions

| Column | Description |
|--------|-------------|
| SNPS | Reported SNP ID(s) |
| CHR_ID | Chromosome |
| CHR_POS | Position (GRCh38) |
| P-VALUE | Association p-value |
| PVALUE_MLOG | -log10(p-value) |
| OR or BETA | Effect size |
| 95% CI | Confidence interval |
| RISK ALLELE FREQUENCY | Risk allele frequency |
| DISEASE/TRAIT | Reported trait |
| MAPPED_TRAIT | EFO-mapped trait |
| MAPPED_TRAIT_URI | EFO URI |
| STUDY ACCESSION | GCST identifier |
| PUBMEDID | Publication reference |

## EFO Trait Mapping

```bash
# Extract EFO mappings
awk -F'\t' '{print $36"\t"$35}' gwas_catalog_associations.tsv | \
  sort -u > trait_efo_mapping.tsv

# Download EFO ontology for full hierarchy
wget https://github.com/EBISPOT/efo/releases/download/current/efo.obo
```

## API Endpoints Reference

| Endpoint | Description |
|----------|-------------|
| /associations | All associations |
| /studies | All studies |
| /singleNucleotidePolymorphisms | SNP data |
| /efoTraits | Traits |
| /genes | Gene-based queries |

## Related Resources

- [dbSNP](../../1.1.variant.repositories/dbsnp/) - SNP identifiers
- [Open Targets](../../03.diseases.phenotypes/3.3.disease.gene.associations/open.targets/) - Target validation
- [ClinVar](../../1.1.variant.repositories/clinvar/) - Clinical variants
