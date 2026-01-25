# MASI - Data Dictionary

## Overview

This data dictionary documents MASI microbiome-associated signaling interaction entries.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | masi |
| **Name** | MASI |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| interaction_id | integer | Yes | Unique interaction identifier |
| metabolite.name | string | Yes | Metabolite name |
| target.symbol | string | Yes | Gene/protein symbol |
| target.type | string | Yes | Target type (GPCR, enzyme, etc.) |
| interaction.effect | string | Yes | activate, inhibit, modulate, etc. |
| evidence.pmid | integer | Yes | PubMed reference |

---

## Effect Types

| Effect | Description | Example |
|--------|-------------|---------|
| activate | Increases activity | Butyrate activates GPR43 |
| inhibit | Decreases activity | Butyrate inhibits HDACs |
| modulate | Context-dependent | Bile acids modulate FXR |
| substrate | Acts as substrate | TMA processed by FMO3 |
| induce | Increases expression | LPS induces cytokines |
| suppress | Decreases expression | SCFAs suppress inflammation |

---

## Metabolite Classes

| Class | Examples | Host Effects |
|-------|----------|--------------|
| Short-chain fatty acids | Butyrate, Propionate | Anti-inflammatory |
| Indoles | Indole, IAA | AHR activation |
| Bile acid derivatives | Secondary bile acids | FXR, TGR5 signaling |
| Polyamines | Putrescine, Spermidine | Cell proliferation |
| Vitamins | B12, K2, Folate | Cofactors |
| Amino acid derivatives | GABA, Serotonin precursors | Neuroactive |

---

## Host Target Types

| Type | Examples | Function |
|------|----------|----------|
| GPCR | GPR41, GPR43 | SCFA sensing |
| Nuclear receptor | AHR, FXR, PXR | Transcription |
| Enzyme | FMO3, CYP | Metabolism |
| Pattern recognition | TLR4, NOD2 | Immune activation |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| MASI | Microbiome-Associated Signaling Interactions |
| SCFA | Short-Chain Fatty Acid |
| GPCR | G Protein-Coupled Receptor |
| AHR | Aryl Hydrocarbon Receptor |
| FXR | Farnesoid X Receptor |
| HDAC | Histone Deacetylase |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
