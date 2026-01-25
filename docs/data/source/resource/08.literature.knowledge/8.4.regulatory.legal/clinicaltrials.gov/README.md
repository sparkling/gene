---
id: clinicaltrials.gov
title: "ClinicalTrials.gov"
type: data-source
category: literature
subcategory: regulatory.legal
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [clinical-trials, regulatory, research, interventions, fda]
---

# ClinicalTrials.gov

**Category:** [Literature](../../../README.md) > [Regulatory & Legal](../README.md)

## Overview

ClinicalTrials.gov is the world's largest clinical trials registry, maintained by the National Library of Medicine (NLM). It provides information about publicly and privately supported clinical studies conducted around the world.

The database contains detailed information about trial design, eligibility criteria, interventions, outcomes, and results. Registration is required by law for many trials involving FDA-regulated products.

ClinicalTrials.gov is essential for drug development research, competitive intelligence, patient recruitment, and systematic reviews of clinical evidence.

## Key Statistics

| Metric | Value |
|--------|-------|
| Registered Studies | 500,000+ |
| Studies with Results | 70,000+ |
| Countries | 220+ |
| Sponsors/Collaborators | 400,000+ |
| Last Update | Daily |

## Primary Use Cases

1. Clinical trial discovery and monitoring
2. Drug development competitive intelligence
3. Systematic review of interventions
4. Patient recruitment support
5. Regulatory compliance tracking

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| NCT Number | `NCT[0-9]{8}` | NCT00000001 |
| Secondary IDs | Sponsor-specific | 2020-001234-12 |
| Drug Names | Text | Pembrolizumab |
| Condition | MeSH-based | Breast Cancer |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://clinicaltrials.gov | Search and browse |
| REST API | https://clinicaltrials.gov/api/v2 | New version (2024) |
| Downloads | https://clinicaltrials.gov/ct2/resources/download | Full database |
| RSS Feeds | Available | Study updates |

## License

| Aspect | Value |
|--------|-------|
| License | Public Domain (US Government) |
| Commercial Use | Yes |
| Attribution | Recommended |

## See Also

- [Schema Documentation](./schema.md)
- [FDA openFDA](../fda.openfda/README.md) - Drug approvals
- [PubMed](../../8.1.scientific.literature/pubmed/README.md) - Trial publications
