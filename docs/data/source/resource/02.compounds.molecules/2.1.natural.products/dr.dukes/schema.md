---
id: schema-dr-dukes
title: "Dr. Duke's Phytochemical and Ethnobotanical Databases Schema"
type: schema
parent: README.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./README.md)

# Dr. Duke's Phytochemical and Ethnobotanical Databases Schema

**Document ID:** DR-DUKES-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** USDA Agricultural Research Service (1992-2016)

---

## TL;DR

Dr. Duke's Phytochemical and Ethnobotanical Databases contain approximately 49,788 indexed entries covering plants, chemicals, biological activities, and ethnobotanical uses. Core tables include PLANTS, CHEMICALS, ACTIVITIES, and ETHNOBOTANY with relationships linking phytochemicals to their plant sources and documented biological/traditional uses. Data is available as CSV downloads under CC0 license.

---

## Database Overview

| Attribute | Value |
|-----------|-------|
| **Full Name** | Dr. Duke's Phytochemical and Ethnobotanical Databases |
| **Maintainer** | USDA Agricultural Research Service / National Agricultural Library |
| **Creator** | James A. Duke, Ph.D. |
| **Coverage Period** | 1992-2016 |
| **Total Indexed Entries** | ~49,788 |
| **License** | CC0 (Creative Commons Zero - Public Domain) |
| **Interactive Search** | https://phytochem.nal.usda.gov/ |
| **Data Downloads** | https://agdatacommons.nal.usda.gov/articles/dataset/Dr_Duke_s_Phytochemical_and_Ethnobotanical_Databases/24660351 |

---

## Available Downloads

| File | Description |
|------|-------------|
| Duke-Source-CSV.zip | Raw database tables in CSV format |
| DrDukesDatabaseDataDictionary-prelim.csv | Data dictionary (preliminary) |

**Download URLs:**
- CSV Data: https://ndownloader.figshare.com/files/43363335
- Data Dictionary: https://ndownloader.figshare.com/files/43363338

**Note:** The data dictionary is preliminary and some variables are not yet defined. Contact nal-adc-curator@ars.usda.gov for updates.

---

## Entity Categories

The database organizes information into six main searchable entity types:

| Entity Type | Description |
|-------------|-------------|
| **Plant** | Scientific and common plant names |
| **Chemical** | Phytochemical compounds |
| **Biological Activity** | Documented biological/pharmacological activities |
| **Ethnobotanical Activity** | Traditional uses of plants |
| **Ethnobotanical Plant** | Plants with documented traditional uses |
| **Syndrome** | Health conditions and syndromes |

---

## Core Tables Structure

### PLANTS Table

**Description:** Plant species with taxonomic information

| Column | Type | Description |
|--------|------|-------------|
| plant_id | INTEGER | Primary key |
| scientific_name | VARCHAR | Binomial nomenclature |
| common_name | VARCHAR | Common English name |
| family | VARCHAR | Botanical family |
| genus | VARCHAR | Genus name |
| species | VARCHAR | Species epithet |
| part | VARCHAR | Plant part (leaf, root, etc.) |
| synonyms | TEXT | Alternative scientific names |

**Example Records:**
- Curcuma longa (Turmeric)
- Zingiber officinale (Ginger)
- Allium sativum (Garlic)
- Camellia sinensis (Tea)

---

### CHEMICALS Table

**Description:** Phytochemical compounds found in plants

| Column | Type | Description |
|--------|------|-------------|
| chemical_id | INTEGER | Primary key |
| chemical_name | VARCHAR | Chemical/compound name |
| cas_number | VARCHAR | CAS Registry Number |
| formula | VARCHAR | Molecular formula |
| molecular_weight | DECIMAL | Molecular weight |
| smiles | VARCHAR | SMILES notation (if available) |
| inchikey | VARCHAR | InChI Key (if available) |
| chemical_class | VARCHAR | Chemical classification |

**Example Records:**
- Curcumin (CAS: 458-37-7)
- Quercetin (CAS: 117-39-5)
- Allicin (CAS: 539-86-6)
- Gingerol (CAS: 23513-14-6)

---

### PLANT_CHEMICAL Table (Relationship)

**Description:** Links plants to their chemical constituents

| Column | Type | Description |
|--------|------|-------------|
| plant_chemical_id | INTEGER | Primary key |
| plant_id | INTEGER | Foreign key to PLANTS |
| chemical_id | INTEGER | Foreign key to CHEMICALS |
| plant_part | VARCHAR | Part of plant (root, leaf, etc.) |
| concentration_low | DECIMAL | Minimum concentration |
| concentration_high | DECIMAL | Maximum concentration |
| concentration_units | VARCHAR | Units (ppm, %, mg/kg, etc.) |
| reference | VARCHAR | Literature reference |

---

### ACTIVITIES Table

**Description:** Biological activities of chemicals

| Column | Type | Description |
|--------|------|-------------|
| activity_id | INTEGER | Primary key |
| activity_name | VARCHAR | Activity type |
| description | TEXT | Detailed description |
| category | VARCHAR | Activity category |

**Sample Activity Types (1,900+ documented):**

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

---

### CHEMICAL_ACTIVITY Table (Relationship)

**Description:** Links chemicals to their biological activities

| Column | Type | Description |
|--------|------|-------------|
| chem_activity_id | INTEGER | Primary key |
| chemical_id | INTEGER | Foreign key to CHEMICALS |
| activity_id | INTEGER | Foreign key to ACTIVITIES |
| dosage | VARCHAR | Effective dosage (if known) |
| reference | VARCHAR | Literature citation |
| confidence | VARCHAR | Evidence level |

---

### ETHNOBOTANY Table

**Description:** Traditional/ethnobotanical uses of plants

| Column | Type | Description |
|--------|------|-------------|
| ethno_id | INTEGER | Primary key |
| plant_id | INTEGER | Foreign key to PLANTS |
| use | VARCHAR | Traditional use |
| culture | VARCHAR | Cultural origin |
| region | VARCHAR | Geographic region |
| preparation | TEXT | How prepared/administered |
| reference | VARCHAR | Literature/ethnographic source |

**Sample Uses:**

| Use | Description |
|-----|-------------|
| Analgesic | Pain relief |
| Antirheumatic | Arthritis/rheumatism treatment |
| Carminative | Gas/bloating relief |
| Digestive | Digestive aid |
| Febrifuge | Fever reduction |
| Tonic | General health tonic |
| Vulnerary | Wound healing |

---

### TOXICITY Table

**Description:** Toxicity data for chemicals

| Column | Type | Description |
|--------|------|-------------|
| tox_id | INTEGER | Primary key |
| chemical_id | INTEGER | Foreign key to CHEMICALS |
| organism | VARCHAR | Test organism |
| route | VARCHAR | Administration route |
| ld_type | VARCHAR | LD50, LD100, etc. |
| ld_value | DECIMAL | Lethal dose value |
| ld_units | VARCHAR | Units (mg/kg, etc.) |
| reference | VARCHAR | Source citation |

---

### SYNDROMES Table

**Description:** Health conditions and syndromes

| Column | Type | Description |
|--------|------|-------------|
| syndrome_id | INTEGER | Primary key |
| syndrome_name | VARCHAR | Condition name |
| description | TEXT | Description |
| icd_code | VARCHAR | ICD code (if mapped) |

---

## Entity Relationship Diagram

```
PLANTS
   |
   +---> PLANT_CHEMICAL <---+---> CHEMICALS
   |                        |        |
   |                        |        +---> CHEMICAL_ACTIVITY <---> ACTIVITIES
   |                        |        |
   |                        |        +---> TOXICITY
   |
   +---> ETHNOBOTANY
   |
   +---> (Syndromes via activities)
```

---

## Interactive Search Interface

**URL:** https://phytochem.nal.usda.gov/

### Search Capabilities

1. **Chemical Search**
   - Find chemicals by name
   - Get plants containing the chemical
   - Get biological activities

2. **Plant Search**
   - Find plants by scientific or common name
   - Get all chemicals in the plant
   - Get ethnobotanical uses

3. **Activity Search**
   - Find chemicals with specific activity
   - Get plants containing active chemicals

4. **Ethnobotanical Search**
   - Search traditional uses
   - Find plants used for specific conditions

### Export Options

- Search results downloadable as spreadsheets
- PDF export available

---

## Data Dictionary Fields (Preliminary)

Based on the available data dictionary, here are documented fields:

### Plant Fields

| Field Name | Description | Data Type |
|------------|-------------|-----------|
| BINOMIAL | Scientific binomial name | VARCHAR |
| COMMON | Common name | VARCHAR |
| FAMILY | Plant family | VARCHAR |
| PART | Plant part | VARCHAR |
| SYNONYM | Alternative names | TEXT |

### Chemical Fields

| Field Name | Description | Data Type |
|------------|-------------|-----------|
| CHEMICAL | Chemical name | VARCHAR |
| CAS | CAS Registry Number | VARCHAR |
| MW | Molecular weight | DECIMAL |
| FORMULA | Molecular formula | VARCHAR |

### Activity Fields

| Field Name | Description | Data Type |
|------------|-------------|-----------|
| ACTIVITY | Activity name | VARCHAR |
| CHEMICAL | Associated chemical | VARCHAR |
| DOSAGE | Effective dose | VARCHAR |
| REFERENCE | Citation | VARCHAR |

### Ethnobotany Fields

| Field Name | Description | Data Type |
|------------|-------------|-----------|
| USE | Traditional use | VARCHAR |
| PLANT | Plant name | VARCHAR |
| CULTURE | Cultural source | VARCHAR |
| REGION | Geographic region | VARCHAR |

**Note:** Data dictionary is preliminary. Some fields may be undocumented or have modified meanings in current implementation.

---

## Sample Data Records

### Sample Plant-Chemical Relationship

| Plant | Part | Chemical | Concentration | Units |
|-------|------|----------|---------------|-------|
| Curcuma longa | Rhizome | Curcumin | 1,000-6,000 | ppm |
| Curcuma longa | Rhizome | Turmerone | 500-2,000 | ppm |
| Zingiber officinale | Rhizome | Gingerol | 3,000-5,000 | ppm |
| Allium sativum | Bulb | Allicin | 2,000-6,000 | ppm |

### Sample Chemical-Activity Relationship

| Chemical | Activity | Reference |
|----------|----------|-----------|
| Curcumin | Anti-inflammatory | Duke, 1992 |
| Curcumin | Antioxidant | Duke, 1992 |
| Curcumin | Anticancer | Multiple |
| Quercetin | ACE-Inhibitor | Duke, 1992 |
| Quercetin | Antiviral | Duke, 1992 |
| Allicin | Antibacterial | Duke, 1992 |

### Sample Ethnobotanical Record

| Plant | Use | Culture | Region |
|-------|-----|---------|--------|
| Curcuma longa | Anti-inflammatory | Ayurveda | India |
| Curcuma longa | Digestive | Traditional | Asia |
| Zingiber officinale | Carminative | TCM | China |
| Zingiber officinale | Antiemetic | Folk | Global |

---

## Integration Notes for Knowledge Base

### Key Identifiers for Cross-Referencing

| Identifier | Format | Cross-ref to |
|------------|--------|--------------|
| Scientific Name | Binomial | NCBI Taxonomy, Wikidata, COCONUT |
| CAS Number | Variable | PubChem, ChEMBL, SciFinder |
| Chemical Name | Text | PubChem (name search), ChEBI |

### Mapping to Other Databases

1. **Plants to NCBI Taxonomy**
   - Use scientific name for lookup
   - May require fuzzy matching for synonyms

2. **Chemicals to PubChem**
   - CAS number is most reliable
   - Chemical name as fallback (may have multiple matches)

3. **Chemicals to COCONUT/ChEMBL**
   - Use CAS number or InChI Key
   - Verify with molecular formula

4. **Activities to Ontologies**
   - Map to MeSH terms
   - Map to ChEBI roles
   - Map to Gene Ontology (biological process)

### Data Quality Notes

- **Coverage Period:** 1992-2016 (may miss recent discoveries)
- **Focus:** Traditional medicinal plants (Western, Asian, Global)
- **Strength:** Ethnobotanical context often missing from other databases
- **Limitation:** Quantitative data may be sparse for some entries

---

## Query Patterns

### Find All Chemicals in a Plant

```sql
SELECT p.scientific_name, p.common_name, p.part,
       c.chemical_name, c.cas_number,
       pc.concentration_low, pc.concentration_high, pc.concentration_units
FROM plants p
JOIN plant_chemical pc ON p.plant_id = pc.plant_id
JOIN chemicals c ON pc.chemical_id = c.chemical_id
WHERE p.scientific_name LIKE 'Curcuma%'
ORDER BY c.chemical_name;
```

### Find Plants Containing a Specific Chemical

```sql
SELECT p.scientific_name, p.common_name, p.family, pc.plant_part
FROM plants p
JOIN plant_chemical pc ON p.plant_id = pc.plant_id
JOIN chemicals c ON pc.chemical_id = c.chemical_id
WHERE c.chemical_name = 'Quercetin'
ORDER BY p.scientific_name;
```

### Find Chemicals with Specific Activity

```sql
SELECT c.chemical_name, c.cas_number, a.activity_name, ca.dosage
FROM chemicals c
JOIN chemical_activity ca ON c.chemical_id = ca.chemical_id
JOIN activities a ON ca.activity_id = a.activity_id
WHERE a.activity_name = 'Anti-inflammatory'
ORDER BY c.chemical_name;
```

### Find Ethnobotanical Uses of a Plant

```sql
SELECT p.scientific_name, e.use, e.culture, e.region, e.preparation
FROM plants p
JOIN ethnobotany e ON p.plant_id = e.plant_id
WHERE p.scientific_name LIKE 'Zingiber officinale%'
ORDER BY e.use;
```

### Find Plants Used for Specific Condition

```sql
SELECT p.scientific_name, p.common_name, e.use, e.culture
FROM plants p
JOIN ethnobotany e ON p.plant_id = e.plant_id
WHERE e.use LIKE '%digestive%'
   OR e.use LIKE '%Digestive%'
ORDER BY p.scientific_name;
```

---

## Licensing

**License:** CC0 (Creative Commons Zero)
- Public domain dedication
- No attribution required
- Unrestricted commercial use
- Can be modified and redistributed

---

## Contact and Support

- **Data Questions:** nal-adc-curator@ars.usda.gov
- **Technical Support:** National Agricultural Library
- **Reference:** agref@usda.gov

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Records | 20,000+ |
| Storage | Unknown |
| Last updated | January 2026 |

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | JSON (API) |
| Alternative | CSV |
| Encoding | UTF-8 |

---

## Download

| Source | Method | URL |
|--------|--------|-----|
| DrugBank DUKEs | HTTP | https://www.drugbank.ca/download |
| See main database | API | See main database |

**Access Requirements:** Registration required (free account)

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| DrugBank DUKEs | CC BY-NC 4.0 | No (non-commercial only) |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "42" |
| `name` | string | Entity name | "Curcuma longa" |
| `type` | string | Record type | "plant" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `plant_id` | Primary key identifier for plant species | `42` |
| `scientific_name` | Binomial nomenclature for plant species | `Curcuma longa` |
| `common_name` | Vernacular name for plant species | `Turmeric` |
| `chemical_id` | Identifier for phytochemical compound | `458` |
| `cas_number` | Chemical Abstracts Service registry number | `458-37-7` |
| `activity_name` | Documented biological or pharmacological effect | `Anti-inflammatory` |
| `concentration_low` | Minimum measured concentration of compound in plant | `1000` |
| `concentration_high` | Maximum measured concentration of compound in plant | `6000` |
| `concentration_units` | Unit of measurement for compound concentration | `ppm` |
| `ld_value` | Lethal dose value from toxicity testing | `3.5` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Phytochemical | Chemical compound produced naturally by plants | chemical_id |
| Ethnobotanical Use | Traditional use of plant in cultural medicine practice | ETHNOBOTANY table |
| Biological Activity | Documented pharmacological or physiological effect | activity_name |
| Plant Part | Specific portion of plant containing compound (root, leaf, etc.) | plant_part |
| Binomial Nomenclature | Two-part scientific naming system for species | scientific_name |
| Lethal Dose (LD50) | Dose causing death in 50% of test subjects | ld_value |

### Acronyms

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

---

## References

1. Duke, J.A. (1992) Handbook of Phytochemical Constituents of GRAS Herbs and Other Economic Plants. CRC Press.

2. USDA National Agricultural Library: https://www.nal.usda.gov/

3. Ag Data Commons: https://agdatacommons.nal.usda.gov/

4. Interactive Database: https://phytochem.nal.usda.gov/
