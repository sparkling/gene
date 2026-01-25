---
id: download-ema-herbal
title: "EMA Herbal Medicines Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# EMA Herbal Medicines Download Instructions

## Quick Start

```bash
# Public access - regulatory documents freely available
# URL: https://www.ema.europa.eu/en/human-regulatory-overview/herbal-medicinal-products
```

## Prerequisites

- Web browser for document access
- PDF reader for monographs
- No subscription required

## Access Model

EMA herbal monographs are public regulatory documents:

| Access Type | Cost | Description |
|-------------|------|-------------|
| Web Interface | Free | Browse and search |
| PDF Download | Free | Full monograph documents |
| Bulk Access | Free | All published monographs |

## Document Statistics

| Category | Count | Description |
|----------|-------|-------------|
| Herbal Monographs | 200+ | Comprehensive assessments |
| Community Lists | 80+ | Harmonized substances |
| Assessment Reports | 200+ | Scientific background |
| Herbal Substances | 250+ | Evaluated materials |
| Traditional Use | 150+ | Registrations supported |
| Well-Established Use | 50+ | Clinical evidence |

## Document Types

### Herbal Monograph
Complete regulatory assessment including:
- Herbal substance identification
- Pharmaceutical preparations
- Therapeutic indications
- Posology and administration
- Contraindications
- Special warnings
- Interactions
- Pregnancy/lactation
- Effects on driving
- Undesirable effects
- Overdose information

### Assessment Report
Scientific background for monograph decisions:
- Historical use data
- Clinical trial summaries
- Pharmacological studies
- Toxicological data
- Risk assessment

### Community List Entry
Simplified entry for substances with:
- Long tradition of safe use
- Harmonized conditions across EU

## Download Steps

1. Navigate to EMA herbal medicines portal
2. Search by herbal substance name (Latin binomial)
3. Select the relevant monograph
4. Download PDF documents:
   - Final Monograph
   - Assessment Report
   - List of References
5. Review for regulatory requirements

## Bulk Download

EMA provides access to all published documents:

```
Navigation Path:
1. EMA Website > Human Regulatory > Herbal Medicinal Products
2. Herbal Medicines > Community Monographs
3. Browse alphabetically by substance
4. Download individual PDFs or batch download
```

## Monograph Structure

| Section | Content |
|---------|---------|
| Introduction | Botanical identity, quality standards |
| Pharmaceutical | Preparations, extracts, standardization |
| Clinical Particulars | Indications, dosage, duration |
| Clinical Safety | Contraindications, warnings, interactions |
| Pharmacological | Mechanism (when known) |
| Non-clinical Safety | Toxicology data |

## Regulatory Categories

| Category | Evidence Level | Example |
|----------|----------------|---------|
| Well-Established Use | Clinical trials, systematic reviews | Ginkgo for cognitive function |
| Traditional Use | >=30 years use (>=15 in EU) | Valerian for sleep |
| Community List | Harmonized across EU | Specific substances |

## Data Extraction Example

For integration purposes:

```python
# Extract key data from EMA monographs
# (Manual or PDF parsing required)

monograph_data = {
    'substance': 'Valeriana officinalis L., radix',
    'category': 'Traditional Use',
    'indication': 'Relief of mild nervous tension and sleep disorders',
    'posology': '0.3-3g of dried root, 1-3 times daily',
    'duration': '2-4 weeks',
    'contraindications': ['Children under 12'],
    'ema_reference': 'EMA/HMPC/150848/2017'
}
```

## Verification

Check published monograph count:

| Document Type | Expected Count |
|---------------|----------------|
| Final Monographs | 200+ |
| Assessment Reports | 200+ |
| Community List Entries | 80+ |

## Update Schedule

| Aspect | Value |
|--------|-------|
| Update Frequency | Regulatory review cycles |
| New Monographs | Ongoing HMPC work program |
| Revisions | Based on new safety data |

## Notes

- Official EU regulatory documents
- Public domain (regulatory documents)
- Evidence-based safety assessments
- Quality standards linked to European Pharmacopoeia
- Updates follow regulatory timelines
- Citation recommended but not required
- Focused on European market requirements
