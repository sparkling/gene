# 4.1 Metabolic Pathways - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| Subcategory ID | 4.1 |
| Subcategory Name | Metabolic Pathways |
| Data Sources | KEGG, Reactome, WikiPathways |
| Schema ID | `https://gene.ai/schemas/4.1-metabolic-pathways.json` |

## Unified Fields

These fields are harmonized across all data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `pathway_id` | string | Required (1:1) | Primary identifier for the pathway | KEGG, Reactome, WikiPathways | KEGG: `hsa00010`, Reactome: `R-HSA-1430728`, WikiPathways: `WP554` |
| `pathway_name` | string | Required (1:1) | Human-readable name of the pathway | KEGG, Reactome, WikiPathways | `Glycolysis / Gluconeogenesis`, `Metabolism`, `Cholesterol biosynthesis pathway` |
| `organism` | string | Optional (1:1) | Species/organism the pathway applies to | KEGG, Reactome, WikiPathways | KEGG: `hsa`, Reactome: `Homo sapiens` |
| `gene_entries` | array&lt;entity&gt; | Optional (1:N) | Genes/gene products participating in the pathway | KEGG, Reactome, WikiPathways | Array of objects with `id`, `symbol`, `name` |
| `compound_entries` | array&lt;entity&gt; | Optional (1:N) | Metabolites/compounds in the pathway | KEGG, Reactome, WikiPathways | Array of objects with `id`, `name` |
| `reactions` | array&lt;reaction&gt; | Optional (1:N) | Biochemical reactions comprising the pathway | KEGG, Reactome, WikiPathways | Array of objects with `id`, `substrates`, `products`, `enzyme` |
| `cross_references` | array&lt;xref&gt; | Optional (1:N) | Links to external databases | KEGG, Reactome, WikiPathways | Array of objects with `database`, `id` |

## KEGG-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `kegg_entry_type` | enum | Optional (1:1) | Type of KGML entry | `ortholog`, `enzyme`, `gene`, `group`, `compound`, `map`, `brite`, `other` |
| `kegg_relation_type` | enum | Optional (1:1) | Type of molecular relation | `ECrel`, `PPrel`, `GErel`, `PCrel`, `maplink` |
| `kegg_relation_subtype` | array&lt;string&gt; | Optional (1:N) | Specific interaction type | `phosphorylation`, `activation`, `inhibition` |
| `kegg_graphics` | object | Optional (1:1) | Visual rendering properties (x, y coordinates, colors) | `{"x": 547, "y": 308, "type": "rectangle"}` |
| `kegg_reaction_type` | enum | Optional (1:1) | Reversibility of reaction | `reversible`, `irreversible` |

## Reactome-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `reactome_db_id` | integer | Optional (1:1) | Internal Reactome database identifier | `1430728` |
| `reactome_st_id` | string | Optional (1:1) | Stable identifier persistent across releases | `R-HSA-1430728` |
| `reactome_st_id_version` | string | Optional (1:1) | Versioned stable identifier | `R-HSA-1430728.16` |
| `reactome_schema_class` | string | Optional (1:1) | Neo4j node type classification | `TopLevelPathway`, `Pathway`, `Reaction`, `Complex`, `EntityWithAccessionedSequence`, `SimpleEntity` |
| `compartment` | array&lt;string&gt; | Optional (1:N) | Cellular locations where reaction occurs | `["cytosol", "nucleoplasm"]` |
| `go_biological_process` | object | Optional (1:1) | Associated GO Biological Process annotation | `{"accession": "GO:0008152", "displayName": "metabolic process"}` |
| `literature_references` | array&lt;reference&gt; | Optional (1:N) | PubMed citations supporting the annotation | `[{"pmid": 10583946, "title": "..."}]` |
| `input` | array&lt;string&gt; | Optional (1:N) | Reaction input entities (substrates) | `["ATP", "DHCR7"]` |
| `output` | array&lt;string&gt; | Optional (1:N) | Reaction output entities (products) | `["ADP", "cholesterol"]` |
| `catalyst_activity` | array&lt;object&gt; | Optional (1:N) | Enzymatic activities catalyzing the reaction | `[{"activity": "GO:0016491", "physicalEntity": "DHCR7"}]` |
| `has_event` | array&lt;string&gt; | Optional (1:N) | Child events within a pathway | `["R-HSA-191273"]` |
| `reference_entity` | object | Optional (1:1) | External database reference (UniProt, ChEBI) | `{"databaseName": "UniProt", "identifier": "Q9UBM7"}` |
| `is_in_disease` | boolean | Optional (1:1) | Flag indicating disease-related pathway | `false` |
| `is_inferred` | boolean | Optional (1:1) | Flag indicating orthology inference | `false` |

## WikiPathways-Specific Fields

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `wikipathways_graph_id` | string | Optional (1:1) | Unique identifier within GPML document | `abc123` |
| `wikipathways_datanode_type` | enum | Optional (1:1) | Classification of pathway node | `GeneProduct`, `Metabolite`, `Protein`, `Rna`, `Complex`, `Pathway`, `Unknown` |
| `wikipathways_text_label` | string | Optional (1:1) | Display text on node | `HMGCR` |
| `wikipathways_arrow_head` | enum | Optional (1:1) | Interaction type using MIM notation | `Arrow`, `mim-catalysis`, `mim-stimulation`, `mim-inhibition`, `mim-conversion`, `mim-modification`, `mim-binding`, `mim-cleavage`, `TBar` |
| `wikipathways_group_style` | enum | Optional (1:1) | Visual grouping style for complexes | `None`, `Group`, `Complex`, `Pathway` |
| `wikipathways_state` | array&lt;object&gt; | Optional (1:N) | Post-translational modification markers | `[{"TextLabel": "P", "GraphRef": "tp53node"}]` |
| `wikipathways_shape` | object | Optional (1:N) | Cellular compartment shapes | `Rectangle`, `RoundedRectangle`, `Oval`, `Mitochondria`, `Endoplasmic Reticulum`, `Golgi Apparatus`, `Nucleus`, `Cell` |
| `wikipathways_graphics` | object | Optional (1:1) | Visual rendering properties | `{"CenterX": 218.39, "CenterY": 190.95, "Width": 66.667, "Height": 20.0, "FillColor": "ffcccc"}` |
| `wikipathways_xref_database` | string | Optional (1:1) | External database for cross-reference | `Entrez Gene` |

## Source Metadata

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `_source.database` | string | Required | Source database name | `KEGG`, `Reactome`, `WikiPathways` |
| `_source.version` | string | Optional | Database version | `2024.1` |
| `_source.access_date` | date | Optional | Date data was retrieved | `2026-01-24` |
| `_source.original_id` | string | Optional | Original identifier in source | `hsa00010` |

## Field Mappings by Source

### KEGG Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `pathway_id` | `pathway_id` |
| `name` | `pathway_name` |
| `org` | `organism` |
| `entry` | `gene_entries` |
| `compound` | `compound_entries` |
| `reaction` | `reactions` |
| `link` | `cross_references` |
| `entry_type` | `kegg_entry_type` |
| `relation_type` | `kegg_relation_type` |
| `relation_subtype` | `kegg_relation_subtype` |
| `graphics` | `kegg_graphics` |
| `reaction_type` | `kegg_reaction_type` |

### Reactome Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `stId` | `reactome_st_id` |
| `dbId` | `reactome_db_id` |
| `stIdVersion` | `reactome_st_id_version` |
| `displayName` | `pathway_name` |
| `speciesName` | `organism` |
| `schemaClass` | `reactome_schema_class` |
| `compartment` | `compartment` |
| `goBiologicalProcess` | `go_biological_process` |
| `literatureReference` | `literature_references` |
| `input` | `input` |
| `output` | `output` |
| `catalystActivity` | `catalyst_activity` |
| `hasEvent` | `has_event` |
| `referenceEntity` | `reference_entity` |
| `isInDisease` | `is_in_disease` |
| `isInferred` | `is_inferred` |

### WikiPathways Field Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `wpid` | `pathway_id` |
| `Name` | `pathway_name` |
| `Organism` | `organism` |
| `DataNode` | `gene_entries` |
| `Metabolite` | `compound_entries` |
| `Interaction` | `reactions` |
| `Xref` | `cross_references` |
| `GraphId` | `wikipathways_graph_id` |
| `DataNode_Type` | `wikipathways_datanode_type` |
| `TextLabel` | `wikipathways_text_label` |
| `ArrowHead` | `wikipathways_arrow_head` |
| `Group_Style` | `wikipathways_group_style` |
| `State` | `wikipathways_state` |
| `Shape` | `wikipathways_shape` |
| `Graphics` | `wikipathways_graphics` |
| `Xref_Database` | `wikipathways_xref_database` |
