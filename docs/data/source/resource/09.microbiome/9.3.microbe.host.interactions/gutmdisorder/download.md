---
id: download-gutmdisorder
title: "gutMDisorder Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# gutMDisorder Download Instructions

## Quick Start

```bash
# Download via web interface
# Visit http://bio-annotation.cn/gutMDisorder and use Download section
```

## Prerequisites

- Web browser for download interface
- **Python** or **R** for data processing
- 100MB-500MB storage for full dataset

## No Registration Required

gutMDisorder data is freely available for academic use.

## Download Methods

### Method 1: Web Interface Download

1. Visit http://bio-annotation.cn/gutMDisorder
2. Navigate to "Download" section
3. Select data tables:
   - Microbe-disease associations
   - Microbe information
   - Disease information
4. Download TSV/Excel files

### Method 2: Query-Based Export

1. Visit http://bio-annotation.cn/gutMDisorder/search
2. Search by:
   - Disease name
   - Microbe name/taxon
   - Direction (increased/decreased)
3. Export search results to TSV

### Method 3: Programmatic Access

```python
import pandas as pd

# Note: Check website for current download URLs
# These are example patterns - verify actual URLs

# Download main association table (if direct URL available)
# url = "http://bio-annotation.cn/gutMDisorder/data/associations.tsv"
# associations = pd.read_csv(url, sep="\t")

# Alternative: Parse from web interface
# May require browser automation or manual download
```

### Method 4: By Disease Category

```bash
# If categorized downloads are available:
# wget http://bio-annotation.cn/gutMDisorder/data/metabolic_disorders.tsv
# wget http://bio-annotation.cn/gutMDisorder/data/gastrointestinal.tsv
# wget http://bio-annotation.cn/gutMDisorder/data/neurological.tsv
```

## File Inventory

### Expected Download Files

| File | Size (est.) | Description |
|------|-------------|-------------|
| associations.tsv | ~5 MB | All microbe-disease links |
| microbes.tsv | ~500 KB | Microbe taxonomy info |
| diseases.tsv | ~100 KB | Disease annotations |
| evidence.tsv | ~2 MB | Literature evidence |

### Data Fields (Associations Table)

| Column | Description |
|--------|-------------|
| microbe_taxid | NCBI Taxonomy ID |
| microbe_name | Species/genus name |
| disease_name | Disease name |
| disease_icd | ICD-10 code |
| direction | increased/decreased/altered |
| pmid | PubMed reference |
| sample_size | Study sample size |

## Post-Download Processing

```python
import pandas as pd

# Load associations
df = pd.read_csv("associations.tsv", sep="\t")

# Summary statistics
print(f"Total associations: {len(df)}")
print(f"Unique microbes: {df['microbe_taxid'].nunique()}")
print(f"Unique diseases: {df['disease_name'].nunique()}")

# Direction distribution
print("\nDirection distribution:")
print(df['direction'].value_counts())

# Top diseases by association count
print("\nTop 10 diseases:")
print(df['disease_name'].value_counts().head(10))

# Top microbes by association count
print("\nTop 10 microbes:")
print(df['microbe_name'].value_counts().head(10))

# Filter by disease
diabetes = df[df['disease_name'].str.contains('Diabetes', case=False)]
print(f"\nDiabetes associations: {len(diabetes)}")

# Filter by direction
decreased = df[df['direction'] == 'decreased']
increased = df[df['direction'] == 'increased']
print(f"Decreased in disease: {len(decreased)}")
print(f"Increased in disease: {len(increased)}")
```

### Create Disease-Centric Summary

```python
import pandas as pd

# Load data
df = pd.read_csv("associations.tsv", sep="\t")

# Create disease-centric summary
disease_summary = df.groupby('disease_name').agg({
    'microbe_name': lambda x: len(set(x)),
    'direction': lambda x: {
        'increased': sum(x == 'increased'),
        'decreased': sum(x == 'decreased'),
        'altered': sum(x == 'altered')
    },
    'pmid': 'nunique'
}).reset_index()

disease_summary.columns = ['disease', 'microbe_count', 'directions', 'publication_count']

# Expand directions
disease_summary['increased_count'] = disease_summary['directions'].apply(lambda x: x['increased'])
disease_summary['decreased_count'] = disease_summary['directions'].apply(lambda x: x['decreased'])

disease_summary.to_csv("disease_summary.tsv", sep="\t", index=False)
print("Created disease_summary.tsv")
```

### Create Microbe-Centric Summary

```python
import pandas as pd

# Create microbe-centric summary
microbe_summary = df.groupby(['microbe_taxid', 'microbe_name']).agg({
    'disease_name': lambda x: list(set(x)),
    'direction': lambda x: {d: sum(x == d) for d in ['increased', 'decreased', 'altered']},
    'pmid': 'nunique'
}).reset_index()

# Save
microbe_summary.to_csv("microbe_summary.tsv", sep="\t", index=False)
```

### Generate Dysbiosis Patterns

```python
import pandas as pd

# Find consistent patterns
def get_consistent_microbes(df, min_studies=3):
    """Find microbes with consistent direction across studies."""
    results = []

    for taxid, group in df.groupby('microbe_taxid'):
        name = group['microbe_name'].iloc[0]

        for disease in group['disease_name'].unique():
            disease_group = group[group['disease_name'] == disease]

            if len(disease_group) >= min_studies:
                directions = disease_group['direction'].value_counts()
                total = len(disease_group)
                dominant = directions.index[0]
                consistency = directions.iloc[0] / total

                if consistency >= 0.7:  # 70% consistent
                    results.append({
                        'microbe_taxid': taxid,
                        'microbe_name': name,
                        'disease': disease,
                        'direction': dominant,
                        'studies': total,
                        'consistency': consistency
                    })

    return pd.DataFrame(results)

consistent = get_consistent_microbes(df)
consistent.to_csv("consistent_patterns.tsv", sep="\t", index=False)
print(f"Found {len(consistent)} consistent microbe-disease patterns")
```

## Verification

```bash
# Check file integrity
wc -l associations.tsv

# Check column headers
head -1 associations.tsv

# Count unique values
python3 << 'EOF'
import pandas as pd
df = pd.read_csv("associations.tsv", sep="\t")
print(f"Total rows: {len(df)}")
print(f"Unique microbe-disease pairs: {len(df.groupby(['microbe_taxid', 'disease_name']))}")
print(f"Unique PMIDs: {df['pmid'].nunique()}")
EOF
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Major updates | Annually |
| Literature curation | Continuous |

## Common Issues

- **No direct download link**: Use web interface
- **Duplicate entries**: Same pair from different studies
- **Direction conflicts**: Different studies may show opposite results
- **Missing metadata**: Some fields may be empty
- **Commercial use**: Contact maintainers

## Integration with Other Databases

```python
# Link to GMrepo for sample data
import pandas as pd

gutmd = pd.read_csv("associations.tsv", sep="\t")

# Get unique taxids for GMrepo query
taxids = gutmd['microbe_taxid'].unique()
print(f"Query GMrepo for {len(taxids)} taxa")

# Link to gutMGene for mechanism
# Use microbe_taxid to find gene expression effects

# Link to VMH for metabolic modeling
# Use microbe names to find AGORA models
```

## Related Resources

- [GMrepo](../../9.1.gut.microbiome/gmrepo/README.md) - Sample-level data
- [gutMGene](../../9.1.gut.microbiome/gutmgene/README.md) - Gene expression effects
- [VMH](../vmh/README.md) - Metabolic models
