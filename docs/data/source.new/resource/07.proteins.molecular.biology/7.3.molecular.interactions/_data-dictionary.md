# Data Dictionary: 7.3 Molecular Interactions

## Overview

This data dictionary documents the unified schema for molecular interaction data. Note that the primary data sources (IntAct, Reactome, STRING) are maintained in Category 04: Pathways & Networks. This subcategory serves as a reference point for protein-centric interaction queries.

**Subcategory ID:** 7.3
**Data Sources:** IntAct, Reactome, STRING (via Category 04)
**Schema ID:** https://gene.taxonomy/schemas/7.3-molecular-interactions

---

## Cross-Reference to Category 04

The molecular interaction data is maintained in Category 04 to avoid duplication. See the following locations for complete data:

| Source | Location | Description |
|--------|----------|-------------|
| IntAct | `../../04.pathways.networks/intact/` | Molecular interaction database from EMBL-EBI |
| Reactome | `../../04.pathways.networks/reactome/` | Pathway database with interaction data |
| STRING | `../../04.pathways.networks/string/` | Protein-protein interaction networks |

---

## Unified Fields

These fields provide a standardized representation of molecular interactions for cross-category integration.

| Field | Data Type | Cardinality | Description | Sources | Example |
|-------|-----------|-------------|-------------|---------|---------|
| `interaction_id` | string | Required (1:1) | Unique identifier for the molecular interaction | IntAct, STRING, Reactome | `EBI-1234567` |
| `interaction_type` | string (enum) | Optional (1:1) | Type of molecular interaction | All | `physical`, `genetic`, `functional`, `regulatory` |
| `participants` | array[object] | Required (1:N) | Molecules participating in the interaction | All | See Participants Object Structure |
| `detection_method` | string | Optional (1:1) | Experimental method used to detect the interaction | IntAct, Reactome | `two hybrid`, `co-immunoprecipitation` |
| `confidence_score` | number | Optional (1:1) | Confidence score for the interaction (0-1) | STRING, IntAct | `0.9`, `0.7` |
| `pubmed_ids` | array[integer] | Optional (1:N) | PubMed IDs supporting the interaction | All | `[12345678, 23456789]` |

### Interaction Type Values

| Value | Description | Common Sources |
|-------|-------------|----------------|
| `physical` | Direct physical contact between molecules | IntAct, STRING |
| `genetic` | Genetic interaction (epistasis, synthetic lethality) | STRING |
| `functional` | Functional association without direct contact | STRING |
| `regulatory` | Regulatory relationship (activation, inhibition) | Reactome |

---

## Participants Object Structure

The `participants` array contains objects describing each molecule in the interaction.

| Field | Data Type | Cardinality | Description | Example |
|-------|-----------|-------------|-------------|---------|
| `participant_id` | string | Required (1:1) | Unique identifier for this participant | `EBI-P04637` |
| `participant_type` | string | Optional (1:1) | Type of molecule (protein, small molecule, etc.) | `protein`, `small_molecule`, `nucleic_acid` |
| `uniprot_accession` | string | Optional (1:1) | UniProt accession if participant is a protein | `P04637` |
| `gene_name` | string | Optional (1:1) | Gene symbol for protein-coding participants | `TP53` |
| `organism` | string | Optional (1:1) | Source organism scientific name | `Homo sapiens` |
| `tax_id` | integer | Optional (1:1) | NCBI taxonomy ID | `9606` |

### Participant Type Values

| Value | Description |
|-------|-------------|
| `protein` | Protein molecule |
| `small_molecule` | Small molecule or chemical compound |
| `nucleic_acid` | DNA or RNA molecule |
| `complex` | Protein complex |

---

## Source-Specific Field Availability

### IntAct

IntAct provides experimentally validated molecular interactions with detailed experimental evidence.

| Field | Availability | Notes |
|-------|--------------|-------|
| `interaction_id` | Always | IntAct interaction accession (EBI-XXXXXXX) |
| `interaction_type` | Always | Primarily physical interactions |
| `participants` | Always | With UniProt accessions |
| `detection_method` | Always | MI ontology term for detection method |
| `confidence_score` | Sometimes | MI score when available |
| `pubmed_ids` | Always | Primary literature citations |

### Reactome

Reactome provides pathway-centric interactions with biological context.

| Field | Availability | Notes |
|-------|--------------|-------|
| `interaction_id` | Always | Reactome stable identifier |
| `interaction_type` | Always | Includes regulatory interactions |
| `participants` | Always | With role annotations (enzyme, substrate) |
| `detection_method` | Sometimes | When experimentally validated |
| `confidence_score` | Rarely | Not typically provided |
| `pubmed_ids` | Always | Literature supporting pathway annotations |

### STRING

STRING provides functional protein associations with confidence scores.

| Field | Availability | Notes |
|-------|--------------|-------|
| `interaction_id` | Always | STRING interaction identifier |
| `interaction_type` | Always | Physical, genetic, or functional |
| `participants` | Always | STRING protein identifiers with UniProt mapping |
| `detection_method` | Sometimes | For experimentally validated interactions |
| `confidence_score` | Always | Combined score (0-1) from multiple evidence types |
| `pubmed_ids` | Sometimes | When available from text mining |

---

## Source Metadata

The `_source` object provides metadata about data provenance.

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `primary_source` | string | Name of the primary data source | `IntAct`, `Reactome`, `STRING` |
| `source_id` | string | Original identifier in the source | `EBI-1234567` |
| `extraction_date` | string (date) | Date data was extracted | `2026-01-24` |
| `source_version` | string | Version of the source database | `2026_01` |

---

## Required Fields

The following fields are required for a valid record:

- `interaction_id`
- `participants` (at least 2 participants)

---

## Integration with Protein Data (7.1, 7.2)

Molecular interactions can be linked to protein sequence (7.1) and structure (7.2) data through the following fields:

| This Schema | Links To | Via Field | Notes |
|-------------|----------|-----------|-------|
| `participants[].uniprot_accession` | 7.1 Protein Sequences | `accession` | Links interaction participants to sequence data |
| `participants[].uniprot_accession` | 7.2 Protein Structures | `uniprot_accession` | Links interaction participants to structural data |
| `participants[].gene_name` | 7.1 Protein Sequences | `gene` | Alternative gene-based linking |
| `participants[].tax_id` | 7.1 Protein Sequences | (via organism lookup) | Organism-based filtering |

---

## Usage Notes

1. **For complete interaction data**, refer to Category 04: Pathways & Networks data dictionaries.

2. **Protein-centric queries**: Use the `participants[].uniprot_accession` field to retrieve all interactions involving a specific protein.

3. **Cross-species comparisons**: Use `participants[].tax_id` to filter interactions by organism.

4. **Evidence quality**: Use `confidence_score` and `detection_method` to filter for high-confidence experimental interactions.

5. **Structure-interaction integration**: Link to 7.2 Protein Structures to visualize interaction interfaces on 3D structures.

---

## Related Documentation

- [Category 04: Pathways & Networks Data Dictionary](../../04.pathways.networks/_data-dictionary.md)
- [7.1 Protein Sequences & Annotations Data Dictionary](./../7.1.protein.sequences.annotations/_data-dictionary.md)
- [7.2 Protein Structures Data Dictionary](./../7.2.protein.structures/_data-dictionary.md)
