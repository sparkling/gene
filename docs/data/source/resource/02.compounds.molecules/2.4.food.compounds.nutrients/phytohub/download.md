---
id: download-phytohub
title: "PhytoHub Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# PhytoHub Download Instructions

## Quick Start

```bash
# Access via web interface
# http://phytohub.eu

# Browse and search compounds
# Export search results or request bulk data
```

## Prerequisites

- Web browser for interface access
- MS data analysis software for spectra
- Contact maintainers for programmatic access

## Registration

Free for academic use. Contact for commercial access.

## Download Methods

### Method 1: Web Interface Browse

```bash
# 1. Navigate to PhytoHub
#    http://phytohub.eu

# 2. Compound Browser
#    http://phytohub.eu/compounds
#    - Browse by compound class
#    - View parent-metabolite relationships
#    - Export individual compound data

# 3. Spectral Library
#    http://phytohub.eu/spectra
#    - Browse MS/MS spectra
#    - Download in MSP format
```

### Method 2: Search-Based Export

```bash
# 1. Use compound search
#    http://phytohub.eu/compounds/search

# 2. Search by:
#    - Compound name
#    - InChI Key
#    - Molecular formula
#    - Parent compound (for metabolites)

# 3. Export search results
```

### Method 3: Spectral Library Download

```bash
# 1. Navigate to spectral library
#    http://phytohub.eu/spectra

# 2. Filter by ionization mode (positive/negative)

# 3. Download spectra in:
#    - MSP format (NIST compatible)
#    - MGF format (Mascot Generic)
#    - JSON format

# Use in metabolomics software:
# - MS-DIAL
# - MZmine
# - GNPS
```

### Method 4: Bulk Data Request

```bash
# For complete database:
# 1. Contact PhytoHub team via website
# 2. Describe research purpose
# 3. Request appropriate data files

# Available bulk exports:
# - Complete compound list
# - Parent-metabolite mappings
# - Food source associations
# - Full spectral library
```

## Data Categories

| Category | Content | Format |
|----------|---------|--------|
| Compounds | Structures, properties | CSV, SDF |
| Metabolites | Transformation data | CSV |
| Spectra | MS/MS reference | MSP, MGF |
| Food sources | Dietary associations | CSV |

## Spectral Library Format

### MSP Format

```
Name: Quercetin
PRECURSORMZ: 303.049
PRECURSORTYPE: [M+H]+
IONMODE: Positive
COLLISIONENERGY: 30
Num Peaks: 3
153.018 100
137.023 45
121.028 30
```

### MGF Format

```
BEGIN IONS
TITLE=Quercetin
PEPMASS=303.049
CHARGE=1+
153.018 100
137.023 45
121.028 30
END IONS
```

## Post-Download Processing

```bash
# Preview MSP file
head -20 phytohub_spectra.msp

# Count spectra in library
grep -c "^Name:" phytohub_spectra.msp

# Extract compound names
grep "^Name:" phytohub_spectra.msp | cut -d':' -f2 | sort | uniq

# Convert MSP to MGF (using ms-convert or custom script)
# Many metabolomics tools accept both formats

# Load into MS-DIAL
# 1. Open MS-DIAL
# 2. Go to Settings > Library
# 3. Import MSP file as reference library
```

## Data Fields

### Compound Export

| Field | Description |
|-------|-------------|
| phytohub_id | Database identifier |
| name | Compound name |
| type | parent, phase1, phase2, microbial |
| parent_name | Parent compound (if metabolite) |
| smiles | Structure |
| inchi_key | Structure hash |
| mol_weight | Molecular weight |
| formula | Molecular formula |

### Spectrum Export

| Field | Description |
|-------|-------------|
| compound_id | PhytoHub ID |
| name | Compound name |
| precursor_mz | Precursor ion |
| adduct | Ion adduct type |
| collision_energy | Fragmentation energy |
| ion_mode | Positive/Negative |
| peaks | Fragment list |

## Verification

```bash
# Check MSP file validity
grep -c "^Name:" phytohub_spectra.msp
grep -c "^Num Peaks:" phytohub_spectra.msp
# Counts should match

# Verify spectral format
head -50 phytohub_spectra.msp

# Check for required fields
grep -c "PRECURSORMZ:" phytohub_spectra.msp
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| PhytoHub 1.4 | 2015 | ~100 MB | Current |
| Earlier versions | 2013-2014 | Varies | Archived |

### Version Notes

PhytoHub contains dietary phytochemical metabolite data:
- 1,150+ compounds and 67,000+ spectral annotations
- MS/MS spectral library for metabolomics
- Parent-metabolite relationships
- Phase I, Phase II, and microbial metabolites

## API Access

| Property | Value |
|----------|-------|
| Base URL | `http://phytohub.eu` |
| Rate Limit | N/A (web interface) |
| Auth Required | No |
| Documentation | Contact maintainers |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Database updates | Periodic |
| Spectral additions | Ongoing |
| Last major update | 2015 |

## Common Issues

- **Spectral matching**: Use appropriate mass tolerance (5-10 ppm)
- **Ionization mode**: Match library mode to experimental data
- **Collision energy**: Library spectra at specific CE values
- **Missing spectra**: Not all compounds have reference spectra

## Integration with Metabolomics Workflows

```bash
# Use with MS-DIAL
# - Import as custom library
# - Use for annotation of dietary metabolites

# Use with GNPS
# - Upload as reference library
# - Match against experimental spectra

# Use with MZmine
# - Import as spectral library
# - Compound annotation
```

## Citation Requirements

When using data, cite:
```
PhytoHub: A database of dietary phytochemicals and their metabolites
for metabolomics research. http://phytohub.eu
```

## License

Free for academic use. Contact for commercial licensing.
