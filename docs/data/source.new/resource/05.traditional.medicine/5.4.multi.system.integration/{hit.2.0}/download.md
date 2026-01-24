---
id: download-hit2
title: "HIT 2.0 Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# HIT 2.0 Download Instructions

## Quick Start

```bash
# Visit the web interface for data access
# URL: http://hit2.badd-cao.net/
```

## Prerequisites

- Web browser for data access
- Python 3.8+ for data processing (optional)
- Spreadsheet software for TSV/Excel analysis

## Download Methods

### Primary: Web Interface

Access HIT 2.0 through the main portal:
- **URL**: http://hit2.badd-cao.net/
- **Data Export**: Bulk TSV/Excel files available
- **Search Options**: By ingredient, target, herb

## File Inventory

| Category | Records | Description |
|----------|---------|-------------|
| Herbal Ingredients | 10,000+ | Validated compounds |
| Validated Targets | 5,000+ | Experimentally confirmed |
| Herb-Target Interactions | 100,000+ | Literature-validated |
| Literature References | 50,000+ | PubMed citations |
| Binding Affinity Data | 30,000+ | Quantitative measurements |
| Traditional Systems | Multiple | TCM, Ayurveda, Western |

## Evidence Types

| Evidence Type | Description | Typical Count |
|---------------|-------------|---------------|
| Binding Assay | Direct binding measurement | High |
| Enzyme Inhibition | IC50, Ki values | High |
| Cell-Based Assay | Functional validation | Medium |
| In Vivo | Animal model confirmation | Medium |
| Clinical | Human studies | Limited |

## Data Structure

### Interaction Record Fields
- Ingredient ID (HIT internal)
- Ingredient name
- PubChem CID
- Target ID (UniProt / Gene Symbol)
- Target name
- Interaction type (inhibitor, agonist, etc.)
- Binding affinity (IC50, Ki, Kd, EC50)
- Affinity value and unit
- Experimental method
- PubMed reference

### Multi-System Coverage

| System | Coverage Level |
|--------|----------------|
| Traditional Chinese Medicine | Comprehensive |
| Ayurveda | Major herbs |
| Japanese Kampo | Selected |
| Western Herbal | Common herbs |
| African Traditional | Limited |

## Download Steps

1. Navigate to http://hit2.badd-cao.net/
2. Access Download section
3. Select data tables:
   - Ingredients
   - Targets
   - Interactions
   - References
4. Download TSV/Excel files
5. Import for analysis

## Data Processing Example

```python
import pandas as pd

# Load HIT 2.0 interactions
interactions = pd.read_csv('hit2_interactions.tsv', sep='\t')

# Filter by evidence quality
binding_data = interactions[
    interactions['affinity_type'].isin(['IC50', 'Ki', 'Kd'])
]
print(f"Interactions with binding data: {len(binding_data)}")

# Filter potent interactions
potent = binding_data[
    (binding_data['affinity_value'] <= 1000) &
    (binding_data['affinity_unit'] == 'nM')
]
print(f"Potent interactions (< 1uM): {len(potent)}")

# Cross-reference with predicted databases
# HIT provides validation for BATMAN-TCM, KampoDB predictions
```

## Validation Workflow

Use HIT 2.0 to validate predictions from other databases:

```python
# Load predicted interactions (e.g., from BATMAN-TCM)
predicted = pd.read_csv('batman_predicted.tsv', sep='\t')

# Load HIT validated interactions
validated = pd.read_csv('hit2_interactions.tsv', sep='\t')

# Find overlap (validated predictions)
validated_ids = set(zip(validated['compound_id'], validated['target_id']))
predicted['validated'] = predicted.apply(
    lambda x: (x['compound_id'], x['target_id']) in validated_ids,
    axis=1
)
print(f"Validated predictions: {predicted['validated'].sum()}")
```

## Integration Value

HIT 2.0 serves as validation for other databases:

| Database | Integration Purpose |
|----------|---------------------|
| BATMAN-TCM | Confirm predicted TTIs |
| KampoDB | Validate docking predictions |
| IMPPAT | Cross-reference STITCH predictions |
| SymMap | Ground symptom-target links |

## Verification

Check against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Ingredients | 10,000+ |
| Targets | 5,000+ |
| Interactions | 100,000+ |
| References | 50,000+ |
| Binding Affinity | 30,000+ |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Update Frequency | Periodic literature curation |
| License | Academic use free |
| Commercial Use | Contact maintainers |

## Notes

- Experimental validation only (no predictions)
- Multi-system traditional medicine coverage
- Quantitative binding affinity data included
- Full literature citations
- Gold standard for herbal compound-target validation
- Literature curation may lag publications
