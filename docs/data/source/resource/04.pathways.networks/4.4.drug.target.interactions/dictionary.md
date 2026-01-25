# 4.4 Drug-Target Interactions - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| Subcategory ID | 4.4 |
| Subcategory Name | Drug-Target Interactions |
| Data Sources | STITCH |
| Schema ID | `https://gene.ai/schemas/4.4-drug-target-interactions.json` |

## Unified Fields

These fields are harmonized across all data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `chemical_id` | string | Required (1:1) | STITCH chemical identifier (PubChem CID-based with prefix) | STITCH | `CIDm00002244` |
| `protein_id` | string | Required (1:1) | STRING protein identifier (taxid.ENSP format) | STITCH | `9606.ENSP00000269305` |
| `chemical_name` | string | Optional (1:1) | Preferred chemical name | STITCH | `aspirin` |
| `combined_score` | integer | Optional (1:1) | Probabilistic confidence score (0-1000) | STITCH | `976` |

## STITCH-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `experimental` | integer | Optional (1:1) | Experimental evidence score (binding assays, structures) | `900` |
| `prediction` | integer | Optional (1:1) | Computational prediction score | `0` |
| `database` | integer | Optional (1:1) | Curated database annotation score | `700` |
| `textmining` | integer | Optional (1:1) | Literature co-occurrence score | `650` |
| `experimental_direct` | integer | Optional (1:1) | Direct experimental evidence | `900` |
| `experimental_transferred` | integer | Optional (1:1) | Evidence transferred from homologs | `0` |
| `molecular_weight` | float | Optional (1:1) | Molecular weight in Daltons | `180.157` |
| `smiles` | string | Optional (1:1) | Canonical SMILES representation | `CC(=O)OC1=CC=CC=C1C(=O)O` |
| `inchi_key` | string | Optional (1:1) | International Chemical Identifier key | `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` |
| `mode` | enum | Optional (1:1) | Action type/mechanism | `activation`, `inhibition`, `binding`, `catalysis`, `reaction`, `expression`, `ptmod` |
| `action` | string | Optional (1:1) | Specific action annotation | `inhibitor` |
| `is_directional` | boolean | Optional (1:1) | Whether interaction has directionality | `true` |
| `a_is_acting` | boolean | Optional (1:1) | Whether A (chemical) acts on B (protein) | `true` |

## Chemical ID Formats

STITCH uses PubChem CID-based identifiers with prefixes to indicate stereochemistry handling.

| Prefix | Meaning | Example |
|--------|---------|---------|
| `CIDm` | Merged stereo (stereo-agnostic) | `CIDm00002244` |
| `CIDs` | Stereospecific | `CIDs00002244` |
| `CID0` | PubChem exact | `CID000002244` |

## Action Types

| Action | Description |
|--------|-------------|
| `activation` | Chemical activates protein (agonist) |
| `inhibition` | Chemical inhibits protein (antagonist, inhibitor) |
| `binding` | Physical binding without clear functional effect |
| `catalysis` | Enzyme-substrate relationship |
| `reaction` | Metabolic reaction (product/substrate) |
| `expression` | Affects protein expression (inducer/repressor) |
| `ptmod` | Post-translational modification |

## Evidence Score Channels

STITCH integrates evidence from multiple channels, each scored 0-1000.

| Channel | Description | Data Sources |
|---------|-------------|--------------|
| `experimental` | Laboratory evidence | Binding assays, structural data, kinetic measurements |
| `prediction` | Computational predictions | Structure-based predictions, chemogenomics |
| `database` | Curated annotations | DrugBank, ChEMBL, KEGG Drug, TTD |
| `textmining` | Literature mining | PubMed co-occurrence analysis |

The `combined_score` integrates all channels using a Bayesian approach.

## Source Metadata

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `_source.database` | string | Required | Source database name | `STITCH` |
| `_source.version` | string | Optional | Database version | `5.0` |
| `_source.access_date` | date | Optional | Date data was retrieved | `2026-01-24` |
| `_source.original_id` | string | Optional | Original identifier in source | `CIDm00002244` |

## Field Mappings by Source

### STITCH Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `chemical` | `chemical_id` |
| `protein` | `protein_id` |
| `chemical_name` | `chemical_name` |
| `combined_score` | `combined_score` |
| `experimental` | `experimental` |
| `prediction` | `prediction` |
| `database` | `database` |
| `textmining` | `textmining` |
| `experimental_direct` | `experimental_direct` |
| `experimental_transferred` | `experimental_transferred` |
| `molecular_weight` | `molecular_weight` |
| `SMILES_string` | `smiles` |
| `InChIKey` | `inchi_key` |
| `mode` | `mode` |
| `action` | `action` |
| `is_directional` | `is_directional` |
| `a_is_acting` | `a_is_acting` |

## Score Interpretation

| Score Range | Interpretation |
|-------------|----------------|
| 900-1000 | Highest confidence |
| 700-899 | High confidence |
| 400-699 | Medium confidence |
| 150-399 | Low confidence |
| 0-149 | Very low confidence |
