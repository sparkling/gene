---
id: download-supercyp
title: "SuperCYP Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# SuperCYP Download Instructions

## Quick Start

```bash
# Access via web interface
# http://bioinformatics.charite.de/supercyp

# Browse and search CYP interactions
# Request bulk data from maintainers
```

## Prerequisites

- Web browser for interface access
- Contact maintainers for bulk data access

## Registration

Free for academic use. Contact for commercial access.

## Download Methods

### Method 1: Web Interface (Primary)

```bash
# 1. Navigate to SuperCYP
#    http://bioinformatics.charite.de/supercyp

# 2. Search by drug name
#    - Enter drug name in search box
#    - View CYP interaction profile
#    - Export individual results

# 3. Browse by CYP enzyme
#    - Select CYP isoform (CYP3A4, CYP2D6, etc.)
#    - View all substrates/inhibitors/inducers
#    - Export lists
```

### Method 2: Drug-Specific Query

```bash
# 1. Go to drug search page
# 2. Enter drug name (e.g., "ketoconazole")
# 3. View results showing:
#    - CYP interactions
#    - Interaction types
#    - Ki/IC50 values
#    - Literature references
```

### Method 3: CYP-Specific Query

```bash
# 1. Select CYP enzyme from list
#    - CYP3A4, CYP3A5
#    - CYP2D6
#    - CYP2C9, CYP2C19
#    - CYP1A2
#    - CYP2B6, CYP2E1

# 2. Filter by interaction type:
#    - Substrates
#    - Inhibitors
#    - Inducers

# 3. Export filtered results
```

### Method 4: Bulk Data Request

```bash
# Contact SuperCYP team for bulk data:
# Email: supercyp@charite.de

# Request may include:
# - Complete drug-CYP interaction table
# - CYP enzyme reference data
# - Literature citations
```

## Data Categories

| Category | Content | Access |
|----------|---------|--------|
| Substrates | Drugs metabolized by CYPs | Web browse/export |
| Inhibitors | CYP inhibiting compounds | Web browse/export |
| Inducers | CYP inducing compounds | Web browse/export |
| Ki/IC50 | Quantitative inhibition | Web display |
| References | PubMed citations | Web links |

## Expected Data Fields

### Drug-CYP Interaction Export

| Field | Description |
|-------|-------------|
| Drug Name | Compound name |
| CAS Number | CAS registry number |
| CYP Enzyme | CYP isoform |
| Interaction Type | Substrate/Inhibitor/Inducer |
| Strength | Strong/Moderate/Weak |
| Ki (uM) | Inhibition constant |
| IC50 (uM) | Half-maximal inhibition |
| Evidence | In vitro/Clinical |
| Reference | PubMed ID |

## Post-Download Processing

```bash
# If you obtain bulk data export (CSV/TSV format):

# Preview data
head -20 supercyp_data.csv

# Count interactions by CYP
cut -d',' -f3 supercyp_data.csv | sort | uniq -c | sort -rn

# Filter strong inhibitors
awk -F',' '$5=="strong" && $4=="inhibitor"' supercyp_data.csv > strong_inhibitors.csv

# Load into SQLite
sqlite3 supercyp.db << 'EOF'
CREATE TABLE interactions (
  drug_name TEXT,
  cas_number TEXT,
  cyp_enzyme TEXT,
  interaction_type TEXT,
  strength TEXT,
  ki_value REAL,
  ic50_value REAL,
  evidence TEXT,
  reference TEXT
);
.mode csv
.import supercyp_data.csv interactions
CREATE INDEX idx_cyp ON interactions(cyp_enzyme);
CREATE INDEX idx_drug ON interactions(drug_name);
EOF

# Query CYP3A4 inhibitors
sqlite3 supercyp.db "SELECT drug_name, strength, ki_value FROM interactions WHERE cyp_enzyme='CYP3A4' AND interaction_type='inhibitor' ORDER BY ki_value;"
```

## Verification

```bash
# Check data completeness
wc -l supercyp_data.csv

# Verify CYP coverage
cut -d',' -f3 supercyp_data.csv | sort | uniq

# Check for required fields
head -1 supercyp_data.csv
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| SuperCYP 2.0 | 2013-06 | ~20 MB | Current |
| SuperCYP 1.0 | 2010 | ~10 MB | Archived |

### Version Notes

SuperCYP 2.0 contains comprehensive CYP450 interaction data:
- 3,300 drugs with CYP interaction profiles
- 2,000+ CYP substrates, inhibitors, inducers
- Ki/IC50 values where available
- Links to primary literature
- Major CYP isoforms: 3A4, 2D6, 2C9, 2C19, 1A2, 2B6, 2E1

## API Access

| Property | Value |
|----------|-------|
| Base URL | `http://bioinformatics.charite.de/supercyp` |
| Rate Limit | N/A (web interface) |
| Auth Required | No |
| Documentation | Contact supercyp@charite.de |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database updates | Periodic |
| Literature additions | Ongoing |

## Common Issues

- **Access restrictions**: Academic use free; contact for commercial
- **Export limitations**: Web interface may limit export size
- **Data currency**: Check publication dates for recent drug additions
- **Multiple entries**: Same drug may have multiple CYP interactions

## Integration Notes

Cross-reference with:
- **DrugBank**: Via drug name or CAS number
- **ChEMBL**: Via compound structure
- **PharmGKB**: For clinical annotations
- **SwissADME**: For ADMET predictions

## Clinical DDI Assessment

Use SuperCYP data for:

```
1. Identify victim drugs (substrates of major CYPs)
2. Identify perpetrator drugs (inhibitors/inducers)
3. Assess interaction potential based on:
   - Shared CYP pathways
   - Inhibitor strength (Ki values)
   - Clinical evidence level
```

## CYP Contribution Reference

| Drug | Primary CYP | Fm (fraction metabolized) |
|------|-------------|---------------------------|
| Use SuperCYP to identify primary metabolizing enzymes |
| Helps predict DDI magnitude based on Fm |

## License

Free for academic use. Contact for commercial licensing.
