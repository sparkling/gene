# Dr. Duke's Phytochemical Database - Data Dictionary

## Overview

This data dictionary documents the schema for Dr. Duke's Phytochemical and Ethnobotanical Databases.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | dr.dukes |
| **Name** | Dr. Duke's Phytochemical Database |
| **Parent** | 2.1.natural.products |
| **Total Fields** | 30+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Plant Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| plant_id | integer | 1:1 | Yes | Primary key identifier | 42 |
| scientific_name | string | 1:1 | Yes | Binomial nomenclature | Curcuma longa |
| common_name | string | 1:1 | No | Vernacular name | Turmeric |
| family | string | 1:1 | No | Botanical family | Zingiberaceae |
| genus | string | 1:1 | No | Genus name | Curcuma |
| species | string | 1:1 | No | Species epithet | longa |
| part | string | 1:1 | No | Plant part (leaf, root, etc.) | Rhizome |
| synonyms | string | 1:N | No | Alternative scientific names | Listed synonyms |

### Chemical Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| chemical_id | integer | 1:1 | Yes | Primary key identifier | 458 |
| chemical_name | string | 1:1 | Yes | Compound name | Curcumin |
| cas_number | string | 1:1 | No | CAS Registry Number | 458-37-7 |
| formula | string | 1:1 | No | Molecular formula | C21H20O6 |
| molecular_weight | decimal | 1:1 | No | Molecular weight | 368.38 |
| smiles | string | 1:1 | No | SMILES notation | (if available) |
| inchikey | string | 1:1 | No | InChI Key | (if available) |
| chemical_class | string | 1:1 | No | Chemical classification | Polyphenol |

### Plant-Chemical Relationship

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| plant_chemical_id | integer | 1:1 | Yes | Primary key | 12345 |
| plant_part | string | 1:1 | No | Part containing compound | Root, Leaf |
| concentration_low | decimal | 1:1 | No | Minimum concentration | 1000 |
| concentration_high | decimal | 1:1 | No | Maximum concentration | 6000 |
| concentration_units | string | 1:1 | No | Unit of measurement | ppm |
| reference | string | 1:1 | No | Literature reference | Duke, 1992 |

### Activity Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| activity_id | integer | 1:1 | Yes | Primary key | 101 |
| activity_name | string | 1:1 | Yes | Biological/pharmacological effect | Anti-inflammatory |
| description | string | 1:1 | No | Detailed description | Reduces inflammation |
| category | string | 1:1 | No | Activity category | Pharmacological |
| dosage | string | 1:1 | No | Effective dosage | 50 mg/kg |

### Ethnobotanical Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| ethno_id | integer | 1:1 | Yes | Primary key | 201 |
| use | string | 1:1 | Yes | Traditional use | Digestive aid |
| culture | string | 1:1 | No | Cultural origin | Ayurveda |
| region | string | 1:1 | No | Geographic region | India |
| preparation | string | 1:1 | No | Preparation method | Decoction |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| Plant ID | Integer | 42 | Internal plant identifier |
| Chemical ID | Integer | 458 | Internal chemical identifier |
| CAS Number | XX-XX-X pattern | 458-37-7 | CAS Registry Number |
| BINOMIAL | Genus species | Curcuma longa | Scientific name |

---

## Enumerations

### Activity Types (1,900+ documented)

| Activity | Description |
|----------|-------------|
| ACE-Inhibitor | Angiotensin-converting enzyme inhibition |
| AChE-Inhibitor | Acetylcholinesterase inhibition |
| Analgesic | Pain relief |
| Antibacterial | Against bacteria |
| Anticancer | Antitumor activity |
| Antidiabetic | Blood sugar regulation |
| Antifungal | Against fungi |
| Anti-inflammatory | Reduces inflammation |
| Antioxidant | Reduces oxidative stress |
| Antiviral | Against viruses |
| Cardioprotective | Heart protection |
| Hepatoprotective | Liver protection |
| Hypotensive | Lowers blood pressure |
| Immunostimulant | Boosts immune system |
| Neuroprotective | Brain protection |

### Traditional Uses

| Use | Description |
|-----|-------------|
| Analgesic | Pain relief |
| Antirheumatic | Arthritis/rheumatism treatment |
| Carminative | Gas/bloating relief |
| Digestive | Digestive aid |
| Febrifuge | Fever reduction |
| Tonic | General health tonic |
| Vulnerary | Wound healing |

### Concentration Units

| Unit | Description |
|------|-------------|
| ppm | Parts per million |
| % | Percentage by weight |
| mg/kg | Milligrams per kilogram |
| ug/g | Micrograms per gram |

---

## Entity Relationships

### Plant to Chemicals
- **Cardinality:** N:M
- **Description:** Plants contain multiple chemicals; chemicals found in multiple plants
- **Key Fields:** plant_id, chemical_id

### Chemical to Activities
- **Cardinality:** N:M
- **Description:** Chemicals have multiple activities; activities shared by chemicals
- **Key Fields:** chemical_id, activity_id

### Plant to Ethnobotany
- **Cardinality:** 1:N
- **Description:** Plants have multiple documented traditional uses
- **Key Fields:** plant_id, ethno_id

### Chemical to Toxicity
- **Cardinality:** 1:N
- **Description:** Chemicals may have multiple toxicity measurements
- **Key Fields:** chemical_id, tox_id

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| USDA | United States Department of Agriculture | Database maintainer |
| ARS | Agricultural Research Service | USDA research division |
| CAS | Chemical Abstracts Service | Chemical registry system |
| NAL | National Agricultural Library | Data host |
| CC0 | Creative Commons Zero | Public domain license |
| ppm | Parts Per Million | Concentration unit |
| SMILES | Simplified Molecular-Input Line-Entry System | Structure notation |
| InChI | International Chemical Identifier | IUPAC identifier |
| TCM | Traditional Chinese Medicine | Related practice |
| ACE | Angiotensin-Converting Enzyme | Common activity target |
| AChE | Acetylcholinesterase | Common activity target |
| LD50 | Lethal Dose 50% | Toxicity measure |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| NCBI Taxonomy | Taxon ID | Plant species lookup |
| PubChem | CAS Number | Chemical cross-reference |
| ChEMBL | InChI Key | Bioactivity data |
| COCONUT | InChI Key | Natural products |
| ChEBI | ChEBI ID | Chemical ontology |
| MeSH | MeSH ID | Activity terminology |

---

## Data Quality Notes

1. **Coverage Period:** 1992-2016 (may miss recent discoveries)
2. **Focus:** Traditional medicinal plants (Western, Asian, Global)
3. **Strength:** Ethnobotanical context often missing from other databases
4. **Limitation:** Quantitative data may be sparse for some entries
5. **License:** CC0 (Public Domain) - unrestricted use
6. **Data Dictionary Status:** Preliminary - some fields undocumented
