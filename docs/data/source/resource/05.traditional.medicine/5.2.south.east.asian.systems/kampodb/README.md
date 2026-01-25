---
id: kampodb
title: "KampoDB - Kampo Medicine Database"
type: source
parent: ../README.md
tier: 1
status: active
category: traditional.medicine
subcategory: south.east.asian.systems
tags:
  - kampo
  - japanese-medicine
  - rest-api
  - docking
  - network-pharmacology
---

# KampoDB - Kampo Medicine Database

## Overview

KampoDB provides a comprehensive REST API for exploring Japanese Kampo medicine. The database implements a 4-layer hierarchical structure (Formula -> Crude Drug -> Compound -> Protein) that is fully traversable via clean JSON endpoints, making it one of the most programmatically accessible traditional medicine resources.

A key differentiator is the extensive molecular docking data: over 3 million docking simulations predict binding affinities between natural compounds and human proteins. This enables quantitative assessment of potential drug-target interactions beyond simple association predictions.

KampoDB includes both ligand-based and structure-based target prediction methods, GO/KEGG pathway annotations, and disease associations for systematic analysis of Kampo formula mechanisms.

## Key Statistics

| Metric | Value |
|--------|-------|
| Kampo Formulas | 298 |
| Crude Drugs | 180 |
| Natural Compounds | 3,002 |
| Proteins/Genes | 62,906 |
| Docking Simulations | 3,063,505 |
| Known Target Proteins | 460 |
| Predicted Target Proteins | 1,369 |

## Primary Use Cases

1. Exploring composition of Kampo formulas
2. Target prediction with docking scores
3. Pathway enrichment analysis for formulas
4. Compound-target network construction
5. Comparative analysis with TCM databases

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Formula ID | String code | KT (Kakkonto) |
| Crude Drug ID | Integer 0-179 | 40 (Cinnamon Bark) |
| Compound ID | PubChem CID | 969516 (Curcumin) |
| Protein ID | NCBI Gene ID | 2475 (MTOR) |
| Pathway ID | KEGG hsa prefix | hsa04150 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| REST API | https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/ | Full JSON access |
| Web Interface | https://wakanmoview.inm.u-toyama.ac.jp/kampo/ | Interactive browsing |
| Documentation | On-site | API examples |

## Data Model (4-Layer Hierarchy)

```
Layer 1: Formulas (298) - Traditional Kampo prescriptions
         |
         v (contains)
Layer 2: Crude Drugs (180) - Raw medicinal materials
         |
         v (contains)
Layer 3: Compounds (3,002) - Natural chemicals (PubChem CIDs)
         |
         v (targets via docking)
Layer 4: Proteins (62,906) - Human targets (NCBI Gene IDs)
```

## API Endpoints

| Category | Endpoint | Example |
|----------|----------|---------|
| List all | `/api/{entity}/` | `/api/formula/` |
| Entity info | `/api/{entity}/{id}/info` | `/api/formula/KT/info` |
| Relations | `/api/{entity}/{id}/{target}` | `/api/formula/KT/compound` |
| Docking | `/api/docking/compound/{id}/protein/{id}` | Binding affinity |
| Enrichment | `/api/{entity}/{id}/pathway` | KEGG pathways |

## API Examples

```bash
# List all formulas
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/"

# Get compounds in Kakkonto formula
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/formula/KT/compound"

# Get docking score
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/docking/compound/969516/protein/2475"

# Get protein pathways
curl "https://wakanmoview.inm.u-toyama.ac.jp/kampo/api/protein/2475/pathway"
```

## Docking Data

| Field | Description |
|-------|-------------|
| affinity_kcal_mol | Binding energy (more negative = stronger) |
| domain_id | Protein domain identifier |
| target_name | Domain with residue range |
| ligand_atoms | Compound atom count |
| protein_atoms | Target domain atom count |

## Limitations

- No bulk download endpoint (query individual entities)
- Rate limiting not documented
- No confidence scores on compound-target relationships
- Gene IDs require manual UniProt mapping

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
