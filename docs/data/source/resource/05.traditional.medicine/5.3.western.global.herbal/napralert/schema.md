---
id: schema-napralert
title: "NAPRALERT Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, natural-products, ethnobotany, pharmacognosy, literature]
---

# NAPRALERT Schema Documentation

**Document ID:** SCHEMA-NAPRALERT
**Version:** 1.0
**Source Version:** NAPRALERT (Current - 1975-present)

---

## TL;DR

NAPRALERT is the world's largest relational database of natural products with 50+ years of curated literature. It contains 200,000+ literature records covering 60,000+ plant species and 180,000+ natural products. Data spans ethnobotany, pharmacology, chemistry, and toxicology from global sources. Access is subscription-based with web interface for searching.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Literature Records | 200,000+ | Curated citations |
| Plant Species | 60,000+ | Global coverage |
| Natural Products | 180,000+ | Isolated compounds |
| Ethnobotanical Records | 100,000+ | Traditional uses |
| Pharmacological Records | 150,000+ | Bioactivity data |
| Years of Coverage | 1975-present | 50+ years |
| Geographic Regions | Global | All continents |

---

## Entity Relationship Overview

```
                    +---------------------+
                    |  Literature Source  |
                    |    (200,000+)       |
                    +----------+----------+
                               |
              +----------------+----------------+
              |                                 |
              v                                 v
    +---------+---------+           +---------+---------+
    |     Organisms     |           | Natural Products  |
    |     (60,000+)     |           |    (180,000+)     |
    +---------+---------+           +---------+---------+
              |                               |
    +---------+---------+           +---------+---------+
    |         |         |           |         |         |
    v         v         v           v         v         v
+---+---+ +---+---+ +---+---+ +-----+---+ +---+---+ +---+---+
|Ethno- | |Pharma-| |Toxi-  | |Bioactiv-| |Chem-  | |Struct-|
|botany | |cology | |cology | |ity      | |istry  | |ure    |
+-------+ +-------+ +-------+ +---------+ +-------+ +-------+
```

---

## Core Entities

### Organism

**Description:** Biological source of natural products

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| organism_id | String | Yes | NAPRALERT internal ID |
| scientific_name | String | Yes | Binomial nomenclature |
| family | String | Yes | Taxonomic family |
| kingdom | Enum | Yes | Plantae, Fungi, Bacteria, Animalia |
| common_names | Array[CommonName] | No | Names by language/region |
| synonyms | Array[String] | No | Taxonomic synonyms |
| geographic_distribution | Array[String] | No | Native regions |
| parts_used | Array[String] | No | Organs studied |
| ethnobotanical_uses | Array[EthnobotanicalUse] | No | Traditional uses |
| compounds | Array[CompoundRef] | No | Isolated chemicals |

### NaturalProduct (Compound)

**Description:** Chemical compound from natural source

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| compound_id | String | Yes | NAPRALERT internal ID |
| compound_name | String | Yes | Primary name |
| cas_number | String | No | CAS Registry Number |
| molecular_formula | String | No | Chemical formula |
| molecular_weight | Float | No | MW in Daltons |
| compound_class | String | No | Chemical class |
| sources | Array[OrganismRef] | Yes | Biological sources |
| bioactivities | Array[Bioactivity] | No | Pharmacological data |
| structural_data | StructuralData | No | Structure information |

### LiteratureRecord

**Description:** Published literature reference

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| record_id | String | Yes | NAPRALERT citation ID |
| pubmed_id | Integer | No | PubMed ID if available |
| authors | Array[String] | Yes | Author list |
| title | String | Yes | Article title |
| journal | String | Yes | Journal name |
| year | Integer | Yes | Publication year |
| volume | String | No | Volume number |
| pages | String | No | Page range |
| language | String | No | Original language |
| data_types | Array[Enum] | Yes | Ethno, Pharma, Chem, Tox |

---

## Data Category Entities

### EthnobotanicalUse

**Description:** Traditional/indigenous use documentation

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| use_id | String | Yes | Internal ID |
| organism_id | String | Yes | Source organism |
| traditional_use | String | Yes | Therapeutic application |
| preparation_method | String | No | How prepared |
| administration_route | String | No | How administered |
| dosage | String | No | Traditional dosage |
| geographic_origin | String | Yes | Region/country |
| cultural_group | String | No | Indigenous group |
| plant_part | String | No | Part used |
| literature_ref | String | Yes | Source citation |

### Bioactivity

**Description:** Pharmacological activity record

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| activity_id | String | Yes | Internal ID |
| compound_id | String | No | If compound-specific |
| organism_id | String | No | If whole organism |
| activity_type | String | Yes | Type of bioactivity |
| target | String | No | Biological target |
| test_system | String | Yes | In vitro/vivo model |
| result | String | Yes | Activity result |
| dose | String | No | Test concentration |
| unit | String | No | Concentration unit |
| effect_type | Enum | Yes | Active/Inactive/Moderate |
| literature_ref | String | Yes | Source citation |

### ChemistryRecord

**Description:** Chemical isolation/identification data

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| chemistry_id | String | Yes | Internal ID |
| compound_id | String | Yes | Compound reference |
| organism_id | String | Yes | Source organism |
| plant_part | String | No | Tissue source |
| isolation_method | String | No | Extraction method |
| yield | String | No | Isolation yield |
| identification_method | String | No | Structure elucidation |
| literature_ref | String | Yes | Source citation |

### ToxicologyRecord

**Description:** Safety and toxicity data

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| toxicology_id | String | Yes | Internal ID |
| compound_id | String | No | If compound-specific |
| organism_id | String | No | If whole organism |
| toxicity_type | Enum | Yes | Acute, chronic, etc. |
| test_system | String | Yes | Model organism |
| endpoint | String | Yes | LD50, LC50, etc. |
| value | Float | No | Numeric result |
| unit | String | No | Units |
| route | String | No | Exposure route |
| adverse_effects | Array[String] | No | Observed effects |
| literature_ref | String | Yes | Source citation |

---

## Activity Classifications

### Bioactivity Categories

| Category | Examples |
|----------|----------|
| Antimicrobial | Antibacterial, antifungal, antiviral |
| Anti-inflammatory | COX inhibition, cytokine modulation |
| Anticancer | Cytotoxicity, apoptosis induction |
| Cardiovascular | Vasodilation, cardiotonic |
| CNS | Sedative, analgesic, anxiolytic |
| Metabolic | Antidiabetic, lipid-lowering |
| Immunological | Immunostimulant, immunosuppressant |
| Antioxidant | Free radical scavenging, metal chelation |
| Gastrointestinal | Antiulcer, antidiarrheal |
| Endocrine | Hormone modulation |

### Test System Types

| System | Description |
|--------|-------------|
| In vitro enzyme | Purified enzyme assays |
| In vitro cell | Cell culture models |
| In vitro receptor | Receptor binding |
| In vivo animal | Whole animal models |
| Ex vivo tissue | Isolated tissue studies |
| Clinical | Human studies |
| Computational | In silico predictions |

---

## Geographic Coverage

### Regional Distribution

| Region | Coverage Level |
|--------|----------------|
| Asia | Comprehensive (TCM, Ayurveda, Kampo) |
| Africa | Extensive (traditional healers) |
| Americas | Strong (ethnobotany focus) |
| Europe | Good (herbal medicine) |
| Oceania | Moderate (Pacific islands) |
| Global Marine | Growing (marine natural products) |

### Ethnobotanical Systems

| System | Region | Records |
|--------|--------|---------|
| Traditional Chinese Medicine | China | 30,000+ |
| Ayurveda | India | 25,000+ |
| African Traditional Medicine | Africa | 20,000+ |
| Native American Medicine | Americas | 15,000+ |
| European Herbal | Europe | 10,000+ |
| Other Systems | Global | 10,000+ |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | Relational database |
| Access | Web interface queries |
| Export | Custom reports (subscription) |
| Encoding | UTF-8 |

---

## Sample Records

### Organism Record

```json
{
  "organism_id": "NAP_ORG_045678",
  "scientific_name": "Curcuma longa L.",
  "family": "Zingiberaceae",
  "kingdom": "Plantae",
  "common_names": [
    {"name": "Turmeric", "language": "English"},
    {"name": "Haldi", "language": "Hindi"},
    {"name": "Jiang Huang", "language": "Chinese"}
  ],
  "geographic_distribution": ["India", "Southeast Asia"],
  "parts_used": ["rhizome"],
  "compound_count": 235,
  "ethnobotanical_record_count": 450,
  "pharmacological_record_count": 2100
}
```

### Natural Product Record

```json
{
  "compound_id": "NAP_CPD_012345",
  "compound_name": "Curcumin",
  "cas_number": "458-37-7",
  "molecular_formula": "C21H20O6",
  "molecular_weight": 368.38,
  "compound_class": "Diarylheptanoid",
  "sources": [
    {
      "organism_id": "NAP_ORG_045678",
      "organism_name": "Curcuma longa",
      "plant_part": "rhizome",
      "yield": "2-5%"
    }
  ],
  "bioactivity_count": 3500,
  "literature_count": 1200
}
```

### Ethnobotanical Use Record

```json
{
  "use_id": "NAP_ETH_098765",
  "organism": {
    "id": "NAP_ORG_045678",
    "name": "Curcuma longa"
  },
  "traditional_use": "Anti-inflammatory, wound healing",
  "preparation_method": "Fresh rhizome paste",
  "administration_route": "Topical",
  "geographic_origin": "India",
  "cultural_group": "Ayurvedic tradition",
  "plant_part": "rhizome",
  "literature": {
    "id": "NAP_LIT_543210",
    "authors": ["Sharma RA", "Gescher AJ"],
    "journal": "Eur J Cancer",
    "year": 2005
  }
}
```

### Bioactivity Record

```json
{
  "activity_id": "NAP_ACT_234567",
  "compound": {
    "id": "NAP_CPD_012345",
    "name": "Curcumin"
  },
  "activity_type": "Anti-inflammatory",
  "target": "COX-2",
  "test_system": "In vitro enzyme assay",
  "result": "IC50 = 5.1 uM",
  "dose": 5.1,
  "unit": "uM",
  "effect_type": "Active",
  "literature": {
    "id": "NAP_LIT_654321",
    "pubmed_id": 15489888,
    "year": 2004
  }
}
```

---

## Search Capabilities

### Query Types

| Query Type | Description | Fields |
|------------|-------------|--------|
| Organism | By scientific/common name | Name, family, region |
| Compound | By name or CAS | Name, CAS, class |
| Bioactivity | By activity type | Activity, target, result |
| Ethnobotany | By use or region | Use, region, culture |
| Literature | By author or year | Author, journal, year |
| Combined | Multi-field searches | Any combination |

### Search Operators

| Operator | Function | Example |
|----------|----------|---------|
| AND | Both terms | turmeric AND anti-inflammatory |
| OR | Either term | curcumin OR curcuminoid |
| NOT | Exclude | medicinal NOT food |
| NEAR | Proximity | cancer NEAR treatment |
| * | Wildcard | curcu* |

---

## Cross-References

| External DB | ID Type | Mapping |
|-------------|---------|---------|
| PubMed | PMID | Literature citations |
| CAS | Registry Number | Compound identifiers |
| PubChem | CID | Via CAS mapping |
| NCBI Taxonomy | TaxID | Organism identifiers |
| UniProt | Accession | Target proteins |
| ChEMBL | CHEMBL ID | Via structure |

---

## Glossary

| Term | Definition |
|------|------------|
| Ethnobotany | Study of plant-human cultural relationships |
| Pharmacognosy | Study of medicinal substances from natural sources |
| Natural Product | Chemical compound from biological source |
| Bioactivity | Biological effect of substance |
| CAS Number | Chemical Abstracts Service registry number |
| DER | Drug-Extract Ratio |
| Phytochemistry | Chemistry of plant compounds |
| Secondary Metabolite | Non-essential plant compound |

---

## Subscription Tiers

| Tier | Access Level | Use Case |
|------|--------------|----------|
| Free | Limited queries (5-10/day) | Preliminary searches |
| Academic | Unlimited institutional | University research |
| Commercial | Full access + reports | Industry R&D |
| Custom | Tailored data exports | Specialized needs |

---

## References

1. Loub WD, et al. "NAPRALERT: computer handling of natural product research data." J Chem Inf Comput Sci. 1985.
2. University of Illinois at Chicago: https://napralert.pharmacy.uic.edu/
3. Farnsworth NR. "Ethnopharmacology and drug development." Ciba Found Symp. 1994.
