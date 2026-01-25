# PhytoHub - Data Dictionary

## Overview

This data dictionary documents the schema for PhytoHub dietary phytochemical metabolome database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | phytohub |
| **Name** | PhytoHub |
| **Parent** | 2.4.food.compounds.nutrients |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| phytohub_id | string | 1:1 | Yes | PhytoHub identifier | PHUB000001 |
| name | string | 1:1 | Yes | Compound name | Quercetin |
| compound_type | string | 1:1 | Yes | Compound category | parent, phase1, phase2, microbial |
| parent_id | string | 1:1 | No | Parent compound ID | PHUB000001 |
| smiles | string | 1:1 | No | SMILES structure | O=c1c(O)c(-c2ccc(O)c(O)c2)... |
| inchi | string | 1:1 | No | InChI identifier | InChI=1S/... |
| inchi_key | string | 1:1 | Yes | InChI Key | REFJWTPEDVJJIY-UHFFFAOYSA-N |
| molecular_formula | string | 1:1 | No | Chemical formula | C15H10O7 |
| molecular_weight | decimal | 1:1 | No | Exact mass | 302.236 |
| chemical_class | string | 1:1 | No | Chemical classification | Flavonoid, Phenolic acid |

### Metabolite Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| metabolite_id | string | 1:1 | Yes | Primary identifier | PHUB000002 |
| metabolite_name | string | 1:1 | Yes | Metabolite name | Quercetin-3-O-glucuronide |
| transformation | string | 1:1 | Yes | Chemical transformation | Glucuronidation |
| enzyme | string | 1:1 | No | Metabolizing enzyme | UGT1A9 |
| location | string | 1:1 | No | Metabolic location | liver, gut |
| biofluid | string | 1:1 | No | Detection matrix | plasma, urine |

### MS Spectrum Data

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| spectrum_id | integer | 1:1 | Yes | Primary identifier | 1001 |
| precursor_mz | decimal | 1:1 | Yes | Precursor ion m/z | 303.049 |
| adduct | string | 1:1 | Yes | Ion adduct | [M+H]+, [M-H]- |
| collision_energy | decimal | 1:1 | No | CE in eV | 30 |
| ionization_mode | string | 1:1 | Yes | Ionization mode | positive, negative |
| instrument_type | string | 1:1 | No | MS instrument | Q-TOF, Orbitrap |
| peaks | array | 1:N | Yes | m/z, intensity pairs | [[153.018, 100], [137.023, 45]] |

### Food Source Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| source_id | integer | 1:1 | Yes | Primary identifier | 201 |
| food_name | string | 1:1 | Yes | Food name | Onion |
| food_group | string | 1:1 | No | Food category | Vegetables |
| typical_content | decimal | 1:1 | No | Amount in food | 39.21 |
| content_unit | string | 1:1 | No | Unit | mg/100g |

### External IDs

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| database | string | 1:1 | Yes | External database | PubChem, HMDB, ChEBI |
| external_id | string | 1:1 | Yes | External identifier | 5280343 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| PhytoHub ID | PHUB + 6 digits | PHUB000001 | Primary identifier |
| InChI Key | 27 characters | REFJWTPEDVJJIY-UHFFFAOYSA-N | Structure hash |
| PubChem CID | Integer | 5280343 | Chemical database |
| HMDB ID | HMDB + digits | HMDB0005794 | Metabolite database |
| ChEBI ID | CHEBI: + digits | CHEBI:16243 | Chemical ontology |

---

## Enumerations

### Compound Types

| Type | Description |
|------|-------------|
| parent | Dietary phytochemical (food-derived) |
| phase1 | Phase I metabolite (oxidation/reduction) |
| phase2 | Phase II metabolite (conjugation) |
| microbial | Gut microbial metabolite |

### Transformation Types

| Transformation | Description | Location |
|----------------|-------------|----------|
| Glucuronidation | UGT conjugation | Liver |
| Sulfation | SULT conjugation | Liver |
| Methylation | COMT conjugation | Liver |
| Oxidation | CYP450 oxidation | Liver |
| Reduction | Reductase reaction | Liver/gut |
| Ring fission | Aromatic ring cleavage | Gut |
| Dehydroxylation | Hydroxyl removal | Gut |
| Demethylation | Methyl group removal | Gut |

### Ion Adducts

| Adduct | Mode | Mass Shift |
|--------|------|------------|
| [M+H]+ | Positive | +1.008 |
| [M+Na]+ | Positive | +22.990 |
| [M+NH4]+ | Positive | +18.034 |
| [M-H]- | Negative | -1.008 |
| [M+FA-H]- | Negative | +44.998 |
| [M-H2O-H]- | Negative | -19.018 |

### Chemical Classes

| Class | Examples |
|-------|----------|
| Flavonoid | Quercetin, Catechin |
| Phenolic acid | Caffeic acid, Ferulic acid |
| Stilbene | Resveratrol |
| Lignan | Enterolactone |
| Terpene | Limonene |
| Carotenoid | Beta-carotene |
| Alkaloid | Caffeine |

### Instrument Types

| Type | Description |
|------|-------------|
| Q-TOF | Quadrupole time-of-flight |
| Orbitrap | Orbital ion trap |
| QQQ | Triple quadrupole |
| QTRAP | Quadrupole ion trap |
| FT-ICR | Fourier transform ion cyclotron |

---

## Entity Relationships

### Parent to Metabolites
- **Cardinality:** 1:N
- **Description:** One parent produces multiple metabolites
- **Key Fields:** parent_id, metabolite_id

### Compound to MS Spectra
- **Cardinality:** 1:N
- **Description:** Multiple spectra per compound (different conditions)
- **Key Fields:** phytohub_id, spectrum_id

### Compound to Food Sources
- **Cardinality:** N:M
- **Description:** Compounds in multiple foods; foods contain multiple compounds
- **Key Fields:** phytohub_id, source_id

### Compound to External IDs
- **Cardinality:** 1:N
- **Description:** One compound mapped to multiple databases
- **Key Fields:** phytohub_id, database, external_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| MS | Mass Spectrometry | Analysis method |
| MS/MS | Tandem Mass Spectrometry | Fragmentation |
| m/z | Mass-to-charge ratio | MS measurement |
| CE | Collision Energy | Fragmentation parameter |
| Q-TOF | Quadrupole Time-Of-Flight | MS instrument |
| UGT | UDP-Glucuronosyltransferase | Phase II enzyme |
| SULT | Sulfotransferase | Phase II enzyme |
| COMT | Catechol-O-Methyltransferase | Phase II enzyme |
| CYP450 | Cytochrome P450 | Phase I enzyme |
| HMDB | Human Metabolome Database | Reference database |
| MSP | Mass Spectra in NIST format | Spectrum format |
| MGF | Mascot Generic Format | Spectrum format |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PubChem | CID | Chemical data |
| HMDB | HMDB ID | Metabolite reference |
| ChEBI | ChEBI ID | Chemical ontology |
| Phenol-Explorer | Internal | Polyphenol composition |
| MassBank | Accession | Reference spectra |
| GNPS | Spectrum ID | Spectral library |
| MetaCyc | Reaction ID | Metabolic reactions |

---

## Data Quality Notes

1. **Metabolomics Focus:** Designed for dietary biomarker identification
2. **MS Reference Library:** 1,200+ MS/MS spectra for compound identification
3. **Metabolite Coverage:** 600+ parent compounds, 1,000+ metabolites
4. **Gut Microbiome:** Includes microbial transformation products
5. **Multi-Database Links:** Cross-referenced to PubChem, HMDB, ChEBI
6. **Open Access:** Free for research use
7. **Spectrum Quality:** Curated reference spectra with known conditions
