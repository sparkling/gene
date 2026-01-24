# Data Dictionary: 7.2 Protein Structures

## Overview

This data dictionary documents all fields for protein structure data, integrating PDB (experimental structures), AlphaFold DB (AI predictions), and SWISS-MODEL (homology models) into a unified schema.

**Subcategory ID:** 7.2
**Data Sources:** PDB, AlphaFold DB, SWISS-MODEL
**Schema ID:** https://gene.taxonomy/schemas/7.2-protein-structures

---

## Unified Fields

These fields are harmonized across all data sources with consistent semantics.

| Field | Data Type | Cardinality | Description | Sources | Example |
|-------|-----------|-------------|-------------|---------|---------|
| `entry_id` | string | Required (1:1) | Unique identifier for the structure entry. PDB uses 4-character codes, AlphaFold uses AF-{UniProt}-F{n} format, SWISS-MODEL uses internal IDs | PDB, AlphaFold DB, SWISS-MODEL | `1TUP`, `AF-P04637-F1` |
| `uniprot_accession` | string | Optional (1:1) | Link to UniProt entry for the modeled/determined protein | PDB, AlphaFold DB, SWISS-MODEL | `P04637`, `Q9Y6K9` |
| `sequence` | string | Required (1:1) | Amino acid sequence of the protein chain(s) in the structure | PDB, AlphaFold DB, SWISS-MODEL | `MEEPQSDPSVEPPLSQETFSDLWKLL...` |
| `organism` | string | Optional (1:1) | Scientific name of the organism from which the protein originates | PDB, AlphaFold DB, SWISS-MODEL | `Homo sapiens`, `Escherichia coli` |
| `tax_id` | integer | Optional (1:1) | Numeric identifier from NCBI taxonomy database | PDB, AlphaFold DB, SWISS-MODEL | `9606`, `10090`, `562` |
| `structure_type` | string (enum) | Required (1:1) | Whether structure is experimental or predicted | All | `experimental`, `predicted_ai`, `predicted_homology` |
| `chains` | array[object] | Optional (1:N) | Polymer chains in the structure | PDB, AlphaFold DB, SWISS-MODEL | See Chain Object Structure |
| `residues` | array[object] | Optional (1:N) | Per-residue information including positions, names, B-factors | PDB, AlphaFold DB, SWISS-MODEL | See Residue Object Structure |

### Structure Type Values

| Value | Description | Sources |
|-------|-------------|---------|
| `experimental` | Experimentally determined structure (X-ray, EM, NMR) | PDB |
| `predicted_ai` | AI/deep learning predicted structure | AlphaFold DB |
| `predicted_homology` | Homology/comparative modeling predicted structure | SWISS-MODEL |

---

## Chain Object Structure

The `chains` array contains objects with the following structure:

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `chain_id` | string | Single character or short string identifying the polymer chain | `A`, `B`, `L` |
| `sequence` | string | Amino acid sequence of this chain | `MEEPQ...` |
| `entity_type` | string | Type of polymer entity | `polypeptide(L)`, `polyribonucleotide` |

---

## Residue Object Structure

The `residues` array contains objects with the following structure:

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `residue_number` | integer | Position of the residue in the sequence (1-based) | `1`, `102`, `393` |
| `residue_name` | string | Standard three-letter amino acid code | `MET`, `CYS`, `ALA` |
| `chain_id` | string | Chain this residue belongs to | `A` |
| `b_factor` | number | B-factor (experimental) or confidence score (predicted) | `15.23`, `92.5` |
| `occupancy` | number | Fraction of molecules with atom at this position (0-1) | `1.0`, `0.5` |

---

## PDB-Specific Fields (Experimental Structures)

These fields are only available when the source is PDB.

| Field | Data Type | Cardinality | Description | Example |
|-------|-----------|-------------|-------------|---------|
| `title` | string | Optional (1:1) | Descriptive title of the deposited structure | `CRYSTAL STRUCTURE OF THE TETRAMERIZATION DOMAIN OF P53` |
| `experimental_method` | string (enum) | Optional (1:1) | Primary experimental technique used to determine the structure | `X-RAY DIFFRACTION`, `ELECTRON MICROSCOPY`, `SOLUTION NMR` |
| `resolution` | number | Optional (1:1) | Highest resolution of the diffraction/reconstruction data in Angstroms | `1.7`, `2.5`, `3.2` |
| `r_factor_work` | number | Optional (1:1) | Crystallographic R-factor measuring model-data agreement | `0.194`, `0.21` |
| `r_factor_free` | number | Optional (1:1) | Cross-validation R-factor using test set reflections | `0.256`, `0.27` |
| `deposition_date` | string (date) | Optional (1:1) | Date structure was deposited to the PDB | `1994-05-15`, `2020-03-01` |
| `release_status` | string (enum) | Optional (1:1) | Release status of the entry | `REL`, `HOLD`, `OBS` |
| `em_resolution` | number | Optional (1:1) | Resolution of cryo-EM reconstruction in Angstroms | `2.8`, `3.5`, `4.0` |
| `nmr_models` | integer | Optional (1:1) | Number of conformational models in NMR ensemble | `20`, `10` |

### Experimental Method Values

| Value | Description |
|-------|-------------|
| `X-RAY DIFFRACTION` | X-ray crystallography |
| `ELECTRON MICROSCOPY` | Cryo-electron microscopy |
| `SOLUTION NMR` | Solution nuclear magnetic resonance |
| `SOLID-STATE NMR` | Solid-state nuclear magnetic resonance |

### Release Status Values

| Value | Description |
|-------|-------------|
| `REL` | Released and publicly available |
| `HOLD` | On hold, not yet released |
| `OBS` | Obsolete, superseded by another entry |

---

## AlphaFold DB-Specific Fields (AI Predictions)

These fields are only available when the source is AlphaFold DB.

| Field | Data Type | Cardinality | Description | Example |
|-------|-----------|-------------|-------------|---------|
| `alphafold_entry_id` | string | Optional (1:1) | Unique AlphaFold identifier in format AF-{UniProt}-F{fragment} | `AF-P04637-F1`, `AF-Q9Y6K9-F1` |
| `uniprot_id` | string | Optional (1:1) | UniProt mnemonic identifier for the protein | `P53_HUMAN`, `BRCA1_HUMAN` |
| `gene` | string | Optional (1:1) | Gene symbol encoding the protein | `TP53`, `BRCA1` |
| `model_created_date` | string (date) | Optional (1:1) | Date the AlphaFold model was generated | `2022-07-01`, `2023-01-15` |
| `model_version` | integer | Optional (1:1) | AlphaFold model version (v4 is most recent) | `4`, `3` |
| `sequence_length` | integer | Optional (1:1) | Number of residues in the full UniProt sequence | `393`, `1863` |
| `sequence_checksum` | string | Optional (1:1) | CRC64 checksum of the source sequence | `AD5C149FD8106131` |
| `plddt` | array[number] | Optional (1:N) | Per-residue pLDDT confidence scores (0-100) | `[92.5, 45.3, 78.9]` |
| `plddt_mean` | number | Optional (1:1) | Mean pLDDT score across all residues | `85.2` |
| `predicted_aligned_error` | array[array[number]] | Optional (1:N) | Matrix of predicted aligned errors between all residue pairs in Angstroms | `[[0.5, 2.3], [2.3, 0.4]]` |
| `max_predicted_aligned_error` | number | Optional (1:1) | Maximum value of the PAE scale (typically 31.75 Angstroms) | `31.75` |
| `ptm` | number | Optional (1:1) | Global fold quality metric predicting TM-score to true structure | `0.85`, `0.92` |
| `fragment_number` | integer | Optional (1:1) | Fragment number for proteins longer than 1400 residues | `1`, `2` |

### pLDDT Confidence Score Interpretation

| Score Range | Confidence Level | Interpretation |
|-------------|------------------|----------------|
| 90-100 | Very high | High accuracy, suitable for analysis |
| 70-90 | Confident | Good backbone prediction |
| 50-70 | Low | Treat with caution |
| 0-50 | Very low | May be unstructured or incorrectly predicted |

---

## SWISS-MODEL-Specific Fields (Homology Models)

These fields are only available when the source is SWISS-MODEL.

| Field | Data Type | Cardinality | Description | Example |
|-------|-----------|-------------|-------------|---------|
| `template` | string | Optional (1:1) | PDB template used for homology modeling | `1tup.1.A`, `3sak.1.A` |
| `template_identity` | number | Optional (1:1) | Fraction of identical residues in target-template alignment (0-1) | `0.95`, `0.35` |
| `template_similarity` | number | Optional (1:1) | Fraction of similar residues including conservative substitutions (0-1) | `0.98`, `0.55` |
| `coverage` | number | Optional (1:1) | Fraction of target sequence covered by the model (0-1) | `0.50`, `0.85` |
| `model_range` | object | Optional (1:1) | Residue range covered by model | `{"from": 94, "to": 289}` |
| `qmean` | number | Optional (1:1) | Global model quality Z-score (-4 to 0, higher is better, > -1 is excellent) | `-0.85`, `-1.5`, `-2.3` |
| `qmean_disco` | number | Optional (1:1) | Distance constraint-based quality score (0-1, higher is better) | `0.78`, `0.65` |
| `qmean_local` | array[number] | Optional (1:N) | Local quality estimate per residue stored in B-factor column (0-1) | `[0.72, 0.85]` |
| `oligo_state` | string | Optional (1:1) | Biological assembly state of the model | `monomer`, `tetramer`, `dimer` |
| `model_created_date` | string (date) | Optional (1:1) | Date the homology model was computed | `2024-06-15` |
| `resolution` | number | Optional (1:1) | Resolution of the experimental template structure in Angstroms | `1.7`, `2.5` |
| `template_method` | string | Optional (1:1) | Experimental method used to determine the template structure | `X-ray`, `Cryo-EM` |

### Model Range Object Structure

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `from` | integer | First residue position of the modeled region | `94` |
| `to` | integer | Last residue position of the modeled region | `289` |

### QMEAN Score Interpretation

| Score Range | Quality Level | Interpretation |
|-------------|---------------|----------------|
| > -1.0 | Excellent | High-quality model |
| -1.0 to -2.0 | Good | Reliable for most purposes |
| -2.0 to -3.0 | Fair | Use with caution |
| < -3.0 | Poor | Significant quality issues |

---

## Source Metadata

The `_source` object provides metadata about data provenance.

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `primary_source` | string | Name of the primary data source | `PDB`, `AlphaFold DB`, `SWISS-MODEL` |
| `source_id` | string | Original identifier in the source | `1TUP`, `AF-P04637-F1` |
| `extraction_date` | string (date) | Date data was extracted | `2026-01-24` |
| `source_version` | string | Version of the source database | `2026_01` |

---

## Field Mappings

### PDB to Unified Schema

| PDB Field | Unified Field |
|-----------|---------------|
| `entry.id` | `entry_id` |
| `entity_poly.pdbx_seq_one_letter_code` | `sequence` |
| `entity_src_gen.pdbx_gene_src_scientific_name` | `organism` |
| `entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id` | `tax_id` |
| `entity_poly.entity_id` | `chains[].chain_id` |
| `entity.type` | `chains[].entity_type` |
| `struct.title` | `title` |
| `exptl.method` | `experimental_method` |
| `refine.ls_d_res_high` | `resolution` |
| `refine.ls_R_factor_R_work` | `r_factor_work` |
| `refine.ls_R_factor_R_free` | `r_factor_free` |
| `pdbx_database_status.recvd_initial_deposition_date` | `deposition_date` |
| `pdbx_database_status.status_code` | `release_status` |
| `em_3d_reconstruction.resolution` | `em_resolution` |
| `pdbx_nmr_ensemble.conformers_submitted_total_number` | `nmr_models` |

### AlphaFold DB to Unified Schema

| AlphaFold DB Field | Unified Field |
|--------------------|---------------|
| `entryId` | `entry_id` |
| `uniprotAccession` | `uniprot_accession` |
| `uniprotSequence` | `sequence` |
| `organismScientificName` | `organism` |
| `taxId` | `tax_id` |
| `uniprotId` | `uniprot_id` |
| `gene` | `gene` |
| `modelCreatedDate` | `model_created_date` |
| `latestVersion` | `model_version` |
| `sequenceLength` | `sequence_length` |
| `sequenceChecksum` | `sequence_checksum` |
| `confidenceScore` / `pLDDT` | `plddt` |
| `predicted_aligned_error` | `predicted_aligned_error` |
| `max_predicted_aligned_error` | `max_predicted_aligned_error` |
| `pTM` | `ptm` |
| `fragmentNumber` | `fragment_number` |

### SWISS-MODEL to Unified Schema

| SWISS-MODEL Field | Unified Field |
|-------------------|---------------|
| `model_id` | `entry_id` |
| `uniprot_ac` | `uniprot_accession` |
| `sequence` | `sequence` |
| `organism` | `organism` |
| `tax_id` | `tax_id` |
| `template` | `template` |
| `identity` | `template_identity` |
| `similarity` | `template_similarity` |
| `coverage` | `coverage` |
| `from` | `model_range.from` |
| `to` | `model_range.to` |
| `qmean` | `qmean` |
| `qmean_disco` | `qmean_disco` |
| `qmean_local` | `qmean_local` |
| `oligo_state` | `oligo_state` |
| `created_date` | `model_created_date` |
| `template_resolution` | `resolution` |
| `template_method` | `template_method` |

---

## Required Fields

The following fields are required for a valid record:

- `entry_id`
- `sequence`
- `structure_type`

---

## Cross-Reference Mappings

| From Source | To Source | From Field | To Field | Mapping Type | Notes |
|-------------|-----------|------------|----------|--------------|-------|
| PDB | UniProt | entity_poly.pdbx_seq_one_letter_code | sequence.value | N:1 | Multiple PDB entries may represent same UniProt protein |
| AlphaFold DB | UniProt | uniprotAccession | accession | 1:1 | Direct mapping via UniProt accession |
| SWISS-MODEL | UniProt | uniprot_ac | accession | N:1 | Multiple models may exist per UniProt entry |
| SWISS-MODEL | PDB | template | entry.id | N:1 | Models reference PDB templates |
