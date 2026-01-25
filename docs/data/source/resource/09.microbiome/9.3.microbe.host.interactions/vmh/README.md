---
id: vmh
title: "VMH - Virtual Metabolic Human"
type: data-source
category: microbiome
subcategory: microbe.host.interactions
parent: ../README.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [metabolism, microbiome, models, reconstruction, flux-balance]
---

# VMH - Virtual Metabolic Human

**Category:** [Microbiome](../../../README.md) > [Microbe-Host Interactions](../README.md)

## Overview

The Virtual Metabolic Human (VMH) is a comprehensive resource connecting human and gut microbial metabolism. It provides genome-scale metabolic reconstructions (GEMs) for human cells and thousands of gut microbes, enabling computational modeling of microbiome-host metabolic interactions.

VMH integrates data on reactions, metabolites, genes, microbes, and diseases to support constraint-based modeling approaches like flux balance analysis. It serves as a platform for personalized computational modeling of human-microbiome metabolism.

VMH is essential for systems biology approaches to microbiome research, metabolic modeling, and understanding diet-microbiome-host interactions.

## Key Statistics

| Metric | Value |
|--------|-------|
| Microbial GEMs | 7,200+ |
| Reactions | 6,000+ (human) |
| Metabolites | 4,000+ |
| Microbe-Disease Links | 15,000+ |
| Last Update | 2024 |

## Primary Use Cases

1. Genome-scale metabolic modeling
2. Microbiome-host metabolic analysis
3. Personalized nutrition modeling
4. Drug-microbiome interactions
5. Metabolic engineering design

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| VMH Reaction ID | Text | EX_glc(e) |
| VMH Metabolite ID | Text | glc_D[e] |
| Microbe ID | AGORA naming | Bacteroides_thetaiotaomicron |
| NCBI Taxon ID | Numeric | 818 |
| KEGG Reaction | `R[0-9]+` | R00200 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.vmh.life | Browse and search |
| Downloads | https://www.vmh.life/#download | SBML models |
| API | https://www.vmh.life/api | REST API |
| AGORA | Integrated | Microbial models |

## License

| Aspect | Value |
|--------|-------|
| License | Creative Commons Attribution 4.0 |
| Commercial Use | Yes |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [gutMDisorder](../gutmdisorder/README.md) - Disease associations
- [HMP](../../9.2.body.site.microbiomes/hmp/README.md) - Reference microbiomes
- [KEGG](../../../../04.pathways.networks/4.1.metabolic.pathways/kegg/) - Pathway data
