---
id: batman.tcm
title: "BATMAN-TCM 2.0"
type: data-source
category: traditional-medicine
subcategory: traditional-chinese-medicine
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [tcm, traditional-chinese-medicine, herb-target, network-pharmacology, api]
---

# BATMAN-TCM 2.0

**Category:** [Traditional Medicine](../../_index.md) > [Traditional Chinese Medicine](../_index.md)

## Overview

BATMAN-TCM 2.0 (Bioinformatics Analysis Tool for Molecular mechANism of Traditional Chinese Medicine) is the most comprehensive Traditional Chinese Medicine database with programmatic API access. It provides bidirectional query support linking TCM ingredients, herbs, and formulas to molecular targets, enabling network pharmacology analysis of traditional medicines.

The database integrates curated known interactions from ChEMBL, CTD, STITCH, and DGIdb with machine learning-predicted interactions achieving ROC AUC of 0.9663. Version 2.0 represents a significant expansion with 3.16x more ingredients and 62x more known target-TCM interactions compared to version 1.0.

BATMAN-TCM supports both traditional queries (TCM to targets) and reverse queries (disease genes to candidate herbs), making it valuable for drug discovery and mechanism elucidation of traditional formulas.

## Key Statistics

| Metric | Value |
|--------|-------|
| TCM Formulas | 54,832 |
| Herbs | 8,404 |
| Ingredients/Compounds | 39,171 |
| Known TTIs | 17,068 |
| Predicted TTIs | 2,319,272 |
| Target Proteins | ~15,000+ |
| Prediction ROC AUC | 0.9663 |

## Primary Use Cases

1. Network pharmacology analysis of TCM formulas and herbs
2. Target prediction for natural compounds from Chinese medicine
3. Reverse drug discovery: finding herbs for disease gene signatures
4. Pathway and disease enrichment analysis of TCM mechanisms
5. Cross-validation of compound-target interactions

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Formula ID | Internal string | Formula-specific |
| Herb ID | Internal string | Herb-specific |
| Ingredient ID | Internal string | BATMAN-TCM-ING-5679 |
| Target ID | UniProt/Gene Symbol | P12345, TP53 |
| PubChem CID | Integer | 5280343 |
| ChEMBL ID | CHEMBL prefix | CHEMBL25 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://bionet.ncpsb.org.cn/batman-tcm/ | Interactive browsing |
| REST API | http://bionet.ncpsb.org.cn/batman-tcm/api | JSON output |
| Bulk Download | Via web interface | TSV files |

## Data Model

```
Formula (54,832) -> Herbs (8,404) -> Ingredients (39,171) -> Targets (~15,000)
                                            |
                    Known TTI: 17,068       |       Predicted TTI: 2,319,272
```

## API Usage

```bash
# Query herb targets
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=herb&query=Ginseng&output=json"

# Query compound targets
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=compound&query=ginsenoside_Rg1&output=json"

# Reverse query: find herbs for gene set
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=genes&genes=TP53,EGFR,MTOR&output=json"
```

## License

| Aspect | Value |
|--------|-------|
| License | CC BY-NC 4.0 |
| Commercial Use | Requires separate agreement |
| Attribution | Required |

## Limitations

- Server may have connectivity issues; implement retry logic
- Rate limits not documented; use 3+ second delays
- Commercial restrictions under CC BY-NC 4.0

## See Also

- [Schema Documentation](./schema.md)
- [HERB Database](../herb/_index.md)
- [SymMap](../symmap/_index.md)
- [HIT 2.0](../../5.4.multi.system.integration/{hit.2.0}/_index.md)
