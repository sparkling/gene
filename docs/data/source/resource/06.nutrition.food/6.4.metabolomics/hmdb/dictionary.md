# HMDB - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | hmdb |
| **Name** | Human Metabolome Database |
| **Total Fields** | 100+ |
| **Last Updated** | 2026-01 |

---

## Core Fields

### Identification

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| accession | string | 1:1 | Yes | HMDB ID | `HMDB0000001` |
| name | string | 1:1 | Yes | Primary name | `1-Methylhistidine` |
| version | integer | 1:1 | No | Record version | `5` |

### Chemical Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| chemical_formula | string | 1:1 | No | Molecular formula | `C7H11N3O2` |
| average_molecular_weight | number | 1:1 | No | Molecular weight (Da) | `169.181` |
| smiles | string | 1:1 | No | SMILES notation | `CN1C=NC...` |
| inchikey | string | 1:1 | No | InChI key | `BRMWTNUJHUMWMS...` |
| cas_registry_number | string | 1:1 | No | CAS number | `332-80-9` |

### Biological Properties

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| biofluid_locations | array | 1:N | No | Biofluids present | `Blood, Urine` |
| tissue_locations | array | 1:N | No | Tissues present | `Muscle, Liver` |
| cellular_locations | array | 1:N | No | Cell compartments | `Cytoplasm` |

### Concentration Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| concentration_value | number | 1:1 | No | Measured concentration | `5.2` |
| concentration_units | string | 1:1 | No | Units | `uM` |
| condition | string | 1:1 | No | Health status | `Normal` |
| biofluid | string | 1:1 | No | Sample type | `Blood` |

---

## Enumerations

### Biofluids

| Value | Description |
|-------|-------------|
| Blood | Whole blood |
| Serum | Blood serum |
| Plasma | Blood plasma |
| Urine | Urinary metabolites |
| CSF | Cerebrospinal fluid |
| Saliva | Salivary metabolites |
| Feces | Fecal metabolites |
| Breast_milk | Maternal milk |

### Physical State

| Value | Description |
|-------|-------------|
| Solid | Solid at room temperature |
| Liquid | Liquid at room temperature |
| Gas | Gaseous at room temperature |

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HMDB | Human Metabolome Database | Database name |
| MS | Mass Spectrometry | Analytical method |
| MS/MS | Tandem Mass Spectrometry | Fragmentation |
| NMR | Nuclear Magnetic Resonance | Spectroscopy |
| CSF | Cerebrospinal Fluid | Biofluid |
| OMIM | Online Mendelian Inheritance in Man | Disease database |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |
| SMPDB | Small Molecule Pathway Database | Pathway database |
| InChI | International Chemical Identifier | Structure identifier |

---

## Data Quality Notes

1. **Concentration references**: All values from peer-reviewed literature
2. **Spectra included**: MS/MS and NMR reference spectra available
3. **Version control**: Records versioned for tracking updates
4. **CC BY-NC 4.0**: Non-commercial use freely available

---

## See Also

- [Schema Definition](./schema.json)
- [Sample Data](./sample.json)
- [License Terms](./license.md)
