# 4.2 Signaling Pathways - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| Subcategory ID | 4.2 |
| Subcategory Name | Signaling Pathways |
| Data Sources | Pathway Commons |
| Schema ID | `https://gene.ai/schemas/4.2-signaling-pathways.json` |

## Unified Fields

These fields are harmonized across all data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `uri` | URI | Required (1:1) | Unique resource identifier for the entity (BioPAX URI) | Pathway Commons | `http://identifiers.org/reactome/R-HSA-109581` |
| `display_name` | string | Required (1:1) | Human-readable entity name | Pathway Commons | `Apoptosis` |
| `organism` | BioSource | Optional (1:1) | Species for the pathway | Pathway Commons | `{"display_name": "Homo sapiens", "tax_id": 9606}` |
| `data_source` | Provenance | Optional (1:1) | Source database of the record | Pathway Commons | `Reactome` |

## Pathway Commons-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `pathway_components` | array&lt;Process&gt; | Optional (1:N) | Contained reactions and controls | `["BiochemicalReaction69541"]` |
| `pathway_order` | array&lt;PathwayStep&gt; | Optional (1:N) | Ordered steps in the pathway | `[{"step": 1}, {"step": 2}]` |
| `left` | array&lt;PhysicalEntity&gt; | Optional (1:N) | Reaction substrates (BiochemicalReaction) | `["Protein69523"]` |
| `right` | array&lt;PhysicalEntity&gt; | Optional (1:N) | Reaction products (BiochemicalReaction) | `["Complex69543"]` |
| `conversion_direction` | enum | Optional (1:1) | Direction of biochemical conversion | `LEFT-TO-RIGHT`, `RIGHT-TO-LEFT`, `REVERSIBLE` |
| `ec_numbers` | array&lt;string&gt; | Optional (1:N) | EC enzyme classification numbers | `["1.1.1.1"]` |
| `control_type` | enum | Optional (1:1) | Type of regulatory control | `ACTIVATION`, `INHIBITION` |
| `controller` | array&lt;Entity&gt; | Optional (1:N) | Entity performing the control | `["TP53 protein"]` |
| `controlled` | Process | Optional (1:1) | Process being controlled | `BiochemicalReaction69541` |
| `entity_reference` | Reference | Optional (1:1) | Canonical sequence or structure reference | `{"db": "UniProt", "id": "P04637"}` |
| `features` | array&lt;EntityFeature&gt; | Optional (1:N) | PTMs, domains, binding sites | `[{"type": "phosphorylation", "position": "S15"}]` |
| `cellular_location` | CellularLocation | Optional (1:1) | Subcellular location from GO | `{"term": "cytosol", "go": "GO:0005829"}` |
| `components` | array&lt;PhysicalEntity&gt; | Optional (1:N) | Members of a complex | `["APAF1", "CYCS"]` |
| `component_stoichiometry` | array&lt;Stoichiometry&gt; | Optional (1:N) | Stoichiometric ratios of complex components | `[{"entity": "APAF1", "coefficient": 7}]` |

## SIF Interaction Types

The Simple Interaction Format (SIF) provides binary interaction representations.

| Type | Description | BioPAX Source |
|------|-------------|---------------|
| `INTERACTS_WITH` | Physical binding | MolecularInteraction |
| `IN_COMPLEX_WITH` | Shared complex membership | Complex |
| `CONTROLS-STATE-CHANGE-OF` | State modification | Control |
| `CONTROLS-TRANSPORT-OF` | Localization change | Control (Transport) |
| `CONTROLS-PHOSPHORYLATION-OF` | Phosphorylation | Control (Modification) |
| `CONTROLS-EXPRESSION-OF` | Transcriptional regulation | TemplateReactionRegulation |
| `CATALYSIS-PRECEDES` | Enzyme cascade ordering | Catalysis ordering |
| `NEIGHBOR_OF` | Same reaction | Neighborhood |
| `REACTS-WITH` | Co-substrates | Shared reaction |
| `USED-TO-PRODUCE` | Substrate to product | Reaction flow |
| `CHEMICAL-AFFECTS` | Chemical regulation | Control |

## Source Metadata

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `_source.database` | string | Required | Source database name | `Pathway Commons` |
| `_source.version` | string | Optional | Database version | `12.0` |
| `_source.access_date` | date | Optional | Date data was retrieved | `2026-01-24` |
| `_source.original_id` | string | Optional | Original identifier in source | `http://identifiers.org/reactome/R-HSA-109581` |

## Field Mappings by Source

### Pathway Commons Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `uri` | `uri` |
| `displayName` | `display_name` |
| `organism` | `organism` |
| `dataSource` | `data_source` |
| `pathwayComponent` | `pathway_components` |
| `pathwayOrder` | `pathway_order` |
| `left` | `left` |
| `right` | `right` |
| `conversionDirection` | `conversion_direction` |
| `eCNumber` | `ec_numbers` |
| `controlType` | `control_type` |
| `controller` | `controller` |
| `controlled` | `controlled` |
| `entityReference` | `entity_reference` |
| `feature` | `features` |
| `cellularLocation` | `cellular_location` |
| `component` | `components` |
| `componentStoichiometry` | `component_stoichiometry` |

## BioPAX Entity Types

Pathway Commons uses BioPAX Level 3 ontology. Key entity types include:

| Entity Type | Description |
|-------------|-------------|
| `Pathway` | Collection of related biological processes |
| `BiochemicalReaction` | Conversion of substrates to products |
| `Control` | Regulatory relationship |
| `Catalysis` | Enzymatic activity |
| `Modulation` | Non-catalytic regulation |
| `TemplateReaction` | Transcription/translation |
| `Complex` | Multi-protein assembly |
| `Protein` | Protein entity |
| `SmallMolecule` | Chemical compound |
| `Gene` | Genetic element |
| `Rna` | RNA molecule |
