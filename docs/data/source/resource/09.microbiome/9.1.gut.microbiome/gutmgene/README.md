---
id: gutmgene
title: "gutMGene"
type: data-source
category: microbiome
subcategory: gut.microbiome
parent: ../README.md
tier: 2
last_updated: 2026-01-23
status: active
tags: [gut, microbiome, gene-expression, host-microbe, regulation]
---

# gutMGene

**Category:** [Microbiome](../../../README.md) > [Gut Microbiome](../README.md)

## Overview

gutMGene is a comprehensive database linking gut microbiota to host gene expression changes. It catalogs experimentally validated associations between specific gut bacteria and changes in host gene expression across various tissues and conditions.

The database integrates data from literature curation of microbiota manipulation experiments (germ-free colonization, antibiotic treatment, probiotic supplementation) that measured host transcriptional responses.

gutMGene is valuable for understanding the molecular mechanisms of microbiota-host crosstalk and identifying microbial targets for therapeutic intervention.

## Key Statistics

| Metric | Value |
|--------|-------|
| Microbe-Gene Links | 50,000+ |
| Microbe Species | 500+ |
| Host Genes | 10,000+ |
| Tissues | 30+ |
| Last Update | 2023 |

## Primary Use Cases

1. Microbe-host gene associations
2. Mechanism of microbiome effects
3. Therapeutic target identification
4. Transcriptional response prediction
5. Probiotic mechanism research

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| NCBI Taxon ID | Numeric | 853 (Faecalibacterium) |
| Gene Symbol | Text | TNF, IL6 |
| Entrez Gene ID | Numeric | 7124 |
| PMID | Numeric | 12345678 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://bio-annotation.cn/gutmgene | Search and browse |
| Downloads | Available | TSV format |
| API | N/A | No public API |

## License

| Aspect | Value |
|--------|-------|
| License | Free for academic use |
| Commercial Use | Contact maintainers |
| Attribution | Required |

## See Also

- [Schema Documentation](./schema.md)
- [GMrepo](../gmrepo/README.md) - Microbiome repository
- [gutMDisorder](../../9.3.microbe.host.interactions/gutmdisorder/README.md) - Disease links
