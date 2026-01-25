---
id: download-tmmc
title: "TM-MC Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# TM-MC Download Instructions

## Quick Start

```bash
# Visit the web interface for data access
# URL: http://tm-mc.org/
```

## Prerequisites

- Web browser for data access
- Python 3.8+ for data processing (optional)
- Spreadsheet software for CSV/Excel analysis

## Download Methods

### Primary: Web Interface

Access TM-MC through the main portal:
- **URL**: http://tm-mc.org/
- **Data Export**: CSV/Excel via web interface
- **Browse Options**: By traditional system, plant, compound

## File Inventory

| Category | Records | Description |
|----------|---------|-------------|
| Medicinal Plants | 2,500+ | Multi-system coverage |
| Traditional Systems | 5+ | Ayurveda, TCM, Unani, etc. |
| Compounds | 8,000+ | Chemical constituents |
| Molecular Targets | 1,500+ | Protein targets |
| Therapeutic Indications | 500+ | Traditional uses |
| Plant-Compound Links | 25,000+ | Associations |

## Traditional Systems Covered

| System | Region | Coverage Level |
|--------|--------|----------------|
| Ayurveda | India | Comprehensive |
| TCM | China | Major herbs |
| Unani | Middle East/India | Selected |
| Siddha | South India | Selected |
| Kampo | Japan | Selected |
| Jamu | Indonesia | Limited |

## Data Structure

### Plant Entity Fields
- Plant ID (TM-MC internal)
- Botanical name
- Common names
- Traditional systems
- Associated compounds
- Therapeutic uses

### Cross-System Features
- Shared plants across systems
- Common compounds linking traditions
- Convergent therapeutic uses
- Overlapping molecular targets

## Download Steps

1. Navigate to http://tm-mc.org/
2. Browse by traditional system, plant, or compound
3. Search for specific entities
4. Export search results as CSV/Excel
5. Combine exports for cross-system analysis

## Cross-System Analysis Example

```python
import pandas as pd

# Load plant data
plants = pd.read_csv('tmmc_plants.csv')

# Find plants used in multiple systems
multi_system = plants[plants['systems'].str.count(',') >= 1]
print(f"Plants in multiple systems: {len(multi_system)}")

# Find Ayurveda-TCM overlap
ayurveda_plants = set(plants[plants['systems'].str.contains('Ayurveda')]['plant_id'])
tcm_plants = set(plants[plants['systems'].str.contains('TCM')]['plant_id'])
overlap = ayurveda_plants & tcm_plants
print(f"Ayurveda-TCM shared plants: {len(overlap)}")

# Load compounds
compounds = pd.read_csv('tmmc_compounds.csv')

# Find compounds in plants from multiple systems
# This enables convergent validation of traditional uses
```

## Unique Value Proposition

TM-MC enables cross-cultural validation:

| Comparison | Insight |
|------------|---------|
| Ayurveda vs TCM | Identify overlapping herbs and uses |
| Multi-system convergence | Validate traditional claims |
| Molecular mechanism | Link tradition to modern biology |

## Verification

Check against published statistics:

| Entity | Expected Count |
|--------|----------------|
| Medicinal Plants | 2,500+ |
| Traditional Systems | 5+ |
| Compounds | 8,000+ |
| Targets | 1,500+ |
| Indications | 500+ |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Update Frequency | Periodic |
| Notification | Check website |
| License | Academic use free |

## Notes

- Coverage varies by traditional system
- Not all plants have molecular data
- Therapeutic term standardization challenges
- Updates depend on literature curation
- Unique cross-system perspective
- Contact maintainers for commercial use
