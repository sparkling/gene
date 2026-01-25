# VMH - Data Dictionary

## Overview

This data dictionary documents Virtual Metabolic Human database entries.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | vmh |
| **Name** | Virtual Metabolic Human |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

| Field Name | Data Type | Required | Description |
|------------|-----------|----------|-------------|
| model_id | string | Yes | Model identifier |
| model_type | string | No | human or microbe |
| organism_name | string | No | Species name |
| reactions_count | integer | No | Number of reactions |
| metabolites_count | integer | No | Number of metabolites |
| genes_count | integer | No | Number of genes |

---

## Compartment Codes

| Code | Compartment |
|------|-------------|
| [c] | Cytosol |
| [e] | Extracellular |
| [m] | Mitochondria |
| [n] | Nucleus |
| [r] | Endoplasmic reticulum |
| [l] | Lysosome |
| [x] | Peroxisome |
| [g] | Golgi |

---

## Reaction ID Prefixes

| Prefix | Type | Description |
|--------|------|-------------|
| EX_ | Exchange | Boundary transport reaction |
| DM_ | Demand | Demand reaction |
| sink_ | Sink | Sink reaction |
| (none) | Internal | Metabolic reaction |

---

## Model Statistics

| Model | Reactions | Metabolites | Use Case |
|-------|-----------|-------------|----------|
| Recon3D | 13,543 | 4,140 | Human metabolism |
| AGORA2 | Variable | Variable | Gut microbes (7,206 models) |

---

## Subsystem Categories

| Category | Subsystems |
|----------|------------|
| Carbohydrate | Glycolysis, TCA cycle, Pentose phosphate |
| Amino acid | Biosynthesis, Degradation |
| Lipid | Fatty acid metabolism, Steroid |
| Nucleotide | Purine, Pyrimidine |
| Energy | Oxidative phosphorylation |
| Transport | Exchange, ABC transporters |

---

## Acronyms

| Acronym | Expansion |
|---------|-----------|
| VMH | Virtual Metabolic Human |
| GEM | Genome-scale Metabolic Reconstruction |
| FBA | Flux Balance Analysis |
| GPR | Gene-Protein-Reaction |
| AGORA | Assembly of Gut Organisms through Reconstruction and Analysis |
| SBML | Systems Biology Markup Language |

---

## See Also

- [schema.json](./schema.json) - JSON Schema definition
- [sample.json](./sample.json) - Example records
