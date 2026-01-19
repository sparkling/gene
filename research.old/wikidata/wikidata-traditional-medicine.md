# Wikidata Coverage of Traditional Medicine Systems

## Overview

This document provides comprehensive research on Wikidata's coverage of three major traditional medicine systems:
- **Traditional Chinese Medicine (TCM)**
- **Kampo (Japanese Traditional Medicine)**
- **Ayurveda (Indian Traditional Medicine)**

Wikidata serves as a structured, machine-readable knowledge base that can complement specialized databases for traditional medicine research.

---

## 1. Traditional Chinese Medicine (TCM) in Wikidata

### 1.1 Core Q-IDs

| Concept | Q-ID | Description |
|---------|------|-------------|
| Traditional Chinese Medicine | Q200253 | Alternative medical practice drawn from traditional medicine in China |
| Chinese herbology | Q3271191 | Study of Chinese herbal medicine |
| Traditional medicine | Q771035 | Parent category for all traditional medicine systems |
| Medicinal plant | Q188840 | Plants or derivatives used to treat medical conditions |
| Crude drug | Q735160 | Unrefined plant/animal/mineral used in traditional medicine |
| Qi (vital energy) | Q844401 | Physical life-force concept in Chinese philosophy |
| Meridian (Chinese medicine) | Q178999 | Pathways through which qi flows |
| Acupuncture | Q131149 | TCM technique using needles at specific points |
| Five elements (Wu Xing) | Q188593 | Five phases theory (Wood, Fire, Earth, Metal, Water) |
| Yin and yang | Q83064 | Dualistic concept of complementary forces |

### 1.2 TCM Formulas and Prescriptions

TCM formulas (prescriptions/decoctions) have limited but growing coverage in Wikidata. Key patterns:

**Formula Classification Q-IDs:**
| Type | Q-ID | Notes |
|------|------|-------|
| Chinese herbal formula | Q5102267 | General class for TCM formulas |
| Decoction | Q1180573 | Liquid preparation method |
| Herbal tea | Q1194458 | Includes medicinal tisanes |

**Example TCM Formula Items:**
- Six Gentlemen decoction (Liu Jun Zi Tang)
- Four Substances Decoction (Si Wu Tang)
- Rehmannia Six Formula (Liu Wei Di Huang Wan)

### 1.3 Chinese Medicinal Plants Coverage

Plants used in TCM can be identified through:
1. Items with `P2175` (medical condition treated) statements
2. Items categorized under Chinese herbology (`P361` part of Q3271191)
3. Items with `P1535` (used by) traditional Chinese medicine

**Notable Plant Q-IDs:**
| Plant | Q-ID | Chinese Name |
|-------|------|--------------|
| Ginseng | Q83371 | Ren Shen (人参) |
| Astragalus | Q311174 | Huang Qi (黄芪) |
| Licorice root | Q218848 | Gan Cao (甘草) |
| Ginger | Q35625 | Sheng Jiang (生姜) |
| Cinnamon | Q28165 | Rou Gui (肉桂) |
| Rehmannia | Q2451052 | Di Huang (地黄) |
| Angelica sinensis | Q1131442 | Dang Gui (当归) |
| Atractylodes | Q2869855 | Bai Zhu (白术) |
| Poria cocos | Q1130961 | Fu Ling (茯苓) |
| Ephedra | Q156903 | Ma Huang (麻黄) |

### 1.4 TCM Concepts

**Core Theoretical Concepts:**
| Concept | Q-ID | Description |
|---------|------|-------------|
| Qi | Q844401 | Vital energy/life force |
| Meridian system | Q178999 | Energy pathways |
| Acupuncture point | Q199831 | Specific points on meridians |
| Zang-fu | Q2461549 | Organ theory system |
| Blood (Xue) | Q7873 | Xue concept (broader) |
| Jing (essence) | Q6204633 | Vital essence concept |
| Shen (spirit) | Q49657 | Spirit/mind concept |

### 1.5 Key Properties for TCM

| Property | P-code | Description | Use Case |
|----------|--------|-------------|----------|
| medical condition treated | P2175 | Disease treated by drug/therapy | Link herbs to conditions |
| has active ingredient | P3781 | Active component of substance | Phytochemical content |
| part of | P361 | Item is part of larger system | Categorize TCM items |
| used by | P1535 | System that uses this item | Tag as TCM |
| instance of | P31 | Item is instance of class | Classify items |
| subclass of | P279 | Item is subclass of | Taxonomic hierarchy |
| medical treatment | P2176 | Treatment for condition | Link conditions to therapies |
| interaction | P769 | Drug interactions | Safety information |
| found in taxon | P703 | Organism containing compound | Source species |

### 1.6 SPARQL Queries for TCM

#### Query 1: All Items Related to Traditional Chinese Medicine
```sparql
# Find all items that are part of or used in Traditional Chinese Medicine
SELECT DISTINCT ?item ?itemLabel ?itemDescription WHERE {
  {
    ?item wdt:P361 wd:Q200253 .  # part of TCM
  } UNION {
    ?item wdt:P1535 wd:Q200253 .  # used by TCM
  } UNION {
    ?item wdt:P361 wd:Q3271191 .  # part of Chinese herbology
  } UNION {
    ?item wdt:P1535 wd:Q3271191 .  # used by Chinese herbology
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en,zh". }
}
ORDER BY ?itemLabel
```

#### Query 2: TCM Medicinal Plants with Treated Conditions
```sparql
# Find medicinal plants used in TCM with conditions they treat
SELECT DISTINCT ?plant ?plantLabel ?condition ?conditionLabel WHERE {
  ?plant wdt:P31/wdt:P279* wd:Q188840 .  # instance of medicinal plant
  ?plant wdt:P2175 ?condition .           # medical condition treated
  {
    ?plant wdt:P1535 wd:Q200253 .         # used by TCM
  } UNION {
    ?plant wdt:P361 wd:Q3271191 .         # part of Chinese herbology
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en,zh". }
}
ORDER BY ?plantLabel
```

#### Query 3: Chinese Herbal Formulas
```sparql
# Find Chinese herbal formulas/prescriptions
SELECT DISTINCT ?formula ?formulaLabel ?formulaDescription WHERE {
  {
    ?formula wdt:P31/wdt:P279* wd:Q5102267 .  # instance of herbal formula
  } UNION {
    ?formula wdt:P31 wd:Q12140 .              # pharmaceutical drug
    ?formula wdt:P1535 wd:Q200253 .           # used by TCM
  } UNION {
    ?formula wdt:P361 wd:Q200253 .            # part of TCM
    ?formula wdt:P31 wd:Q12140 .              # is a drug
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en,zh". }
}
ORDER BY ?formulaLabel
```

#### Query 4: TCM Herbs with Active Ingredients
```sparql
# Find TCM herbs and their active ingredients
SELECT DISTINCT ?herb ?herbLabel ?ingredient ?ingredientLabel WHERE {
  ?herb wdt:P31/wdt:P279* wd:Q188840 .   # medicinal plant
  ?herb wdt:P3781 ?ingredient .           # has active ingredient
  {
    ?herb wdt:P1535 wd:Q200253 .
  } UNION {
    ?herb wdt:P1535 wd:Q3271191 .
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?herbLabel
```

#### Query 5: Count of TCM-Related Items
```sparql
# Count items related to TCM
SELECT ?type (COUNT(DISTINCT ?item) AS ?count) WHERE {
  {
    SELECT ?item ("TCM items" AS ?type) WHERE {
      ?item wdt:P361|wdt:P1535 wd:Q200253 .
    }
  } UNION {
    SELECT ?item ("Chinese herbology items" AS ?type) WHERE {
      ?item wdt:P361|wdt:P1535 wd:Q3271191 .
    }
  } UNION {
    SELECT ?item ("Medicinal plants" AS ?type) WHERE {
      ?item wdt:P31/wdt:P279* wd:Q188840 .
    }
  }
}
GROUP BY ?type
ORDER BY DESC(?count)
```

---

## 2. Kampo (Japanese Traditional Medicine) in Wikidata

### 2.1 Core Q-IDs

| Concept | Q-ID | Description |
|---------|------|-------------|
| Kampo | Q1723498 | Traditional Japanese herbal medicine |
| Japanese Pharmacopoeia | Q6150186 | Official pharmacopoeia of Japan |
| Crude drug | Q735160 | Raw medicinal materials |

### 2.2 Kampo Formulas

Kampo formulas (prescriptions) that exist in Wikidata:

| Formula | Q-ID | Japanese Name | Romaji |
|---------|------|---------------|--------|
| Kakkonto | Q11253809 | 葛根湯 | Kakkonto |
| Hochuekkito | (varies) | 補中益気湯 | Hochuekkito |
| Yokukansan | (varies) | 抑肝散 | Yokukansan |
| Daikenchuto | (varies) | 大建中湯 | Daikenchuto |
| Sho-saiko-to | (varies) | 小柴胡湯 | Shosaikoto |

**Example: Kakkonto (Q11253809)**
- Used for: Initial stages of common cold
- Components: Pueraria root, Ephedra, Jujube, Cinnamon, Peony, Licorice, Ginger
- Listed in: Japanese Pharmacopoeia

### 2.3 Japanese Crude Drugs

The Japanese Pharmacopoeia (17th edition) lists:
- 324 herbal medicines total
- 176 crude drugs
- 35 Kampo extract preparations
- 148 formulations covered by national health insurance

**Key Crude Drug Q-IDs:**
| Drug | Q-ID | Japanese Name |
|------|------|---------------|
| Ephedrae Herba | (search needed) | 麻黄 (Mao) |
| Puerariae Radix | (search needed) | 葛根 (Kakkon) |
| Glycyrrhizae Radix | (search needed) | 甘草 (Kanzo) |
| Zingiberis Rhizoma | (search needed) | 生姜 (Shokyo) |
| Cinnamomi Cortex | (search needed) | 桂皮 (Keihi) |

### 2.4 Links to Specialized Databases

**KampoDB:**
- URL: http://wakanmoview.inm.u-toyama.ac.jp/kampo/
- Contains: 42 Kampo medicines, 54 crude drugs, 1,230 compounds
- Maintained by: University of Toyama, Institute of Natural Medicine
- Wikidata property: Not yet established (as of research date)

**STORK (Standards of Reporting Kampo Products):**
- URL: http://mpdb.nibiohn.go.jp/stork/
- Purpose: Standardized reporting for Kampo in clinical research
- Maintained by: National Institutes of Biomedical Innovation, Health and Nutrition (NIBN)
- Wikidata property: Not yet established

### 2.5 SPARQL Queries for Kampo

#### Query 1: All Kampo-Related Items
```sparql
# Find all items related to Kampo medicine
SELECT DISTINCT ?item ?itemLabel ?itemDescription WHERE {
  {
    ?item wdt:P361 wd:Q1723498 .  # part of Kampo
  } UNION {
    ?item wdt:P1535 wd:Q1723498 .  # used by Kampo
  } UNION {
    ?item wdt:P31 wd:Q1723498 .    # instance of Kampo
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],ja,en". }
}
ORDER BY ?itemLabel
```

#### Query 2: Kampo Formulas
```sparql
# Find Kampo formulas/prescriptions
SELECT DISTINCT ?formula ?formulaLabel ?formulaDescription ?japaneseLabel WHERE {
  {
    ?formula wdt:P31/wdt:P279* wd:Q12140 .  # pharmaceutical drug
    ?formula wdt:P1535 wd:Q1723498 .         # used by Kampo
  } UNION {
    ?formula wdt:P361 wd:Q1723498 .          # part of Kampo
  }
  OPTIONAL { ?formula rdfs:label ?japaneseLabel . FILTER(LANG(?japaneseLabel) = "ja") }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],ja,en". }
}
ORDER BY ?formulaLabel
```

#### Query 3: Japanese Crude Drugs
```sparql
# Find crude drugs in Japanese Pharmacopoeia context
SELECT DISTINCT ?drug ?drugLabel ?drugDescription WHERE {
  ?drug wdt:P31/wdt:P279* wd:Q735160 .  # instance of crude drug
  {
    ?drug wdt:P1535 wd:Q1723498 .       # used by Kampo
  } UNION {
    ?drug wdt:P17 wd:Q17 .              # country: Japan
    ?drug wdt:P31/wdt:P279* wd:Q188840 . # medicinal plant
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],ja,en". }
}
ORDER BY ?drugLabel
```

#### Query 4: Kampo Items with Treated Conditions
```sparql
# Find Kampo medicines and conditions they treat
SELECT DISTINCT ?medicine ?medicineLabel ?condition ?conditionLabel WHERE {
  ?medicine wdt:P1535 wd:Q1723498 .   # used by Kampo
  ?medicine wdt:P2175 ?condition .     # medical condition treated
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],ja,en". }
}
ORDER BY ?medicineLabel
```

---

## 3. Ayurveda in Wikidata

### 3.1 Core Q-IDs

| Concept | Q-ID | Description |
|---------|------|-------------|
| Ayurveda | Q132325 | Traditional medicine system from Indian subcontinent |
| Vata (dosha) | Q1124583 | One of three doshas - air/space element |
| Pitta (dosha) | Q1124576 | One of three doshas - fire/water element |
| Kapha (dosha) | Q1124572 | One of three doshas - earth/water element |
| Dosha | Q766998 | Ayurvedic concept of bodily humor |
| Siddha medicine | Q740784 | South Indian traditional medicine system |
| Unani medicine | Q252255 | Greco-Arabic traditional medicine |

### 3.2 Dosha Concepts

The three doshas form the foundation of Ayurvedic diagnosis and treatment:

| Dosha | Q-ID | Elements | Characteristics |
|-------|------|----------|-----------------|
| Vata | Q1124583 | Air + Space | Dry, cold, light, mobile |
| Pitta | Q1124576 | Fire + Water | Hot, moist, liquid, sharp |
| Kapha | Q1124572 | Earth + Water | Heavy, cold, tender, soft |

### 3.3 Ayurvedic Preparations and Formulations

**Types of Ayurvedic Preparations:**
| Type | Description |
|------|-------------|
| Churna | Powder preparations |
| Kashayam | Decoctions |
| Arishtam/Asavam | Fermented preparations |
| Ghritam | Ghee-based preparations |
| Tailam | Oil-based preparations |
| Rasayana | Rejuvenation formulas |
| Bhasma | Calcined metal/mineral preparations |
| Guggulu | Resin-based formulations |

### 3.4 Indian Medicinal Plants

India has documented extensive traditional knowledge of medicinal plants:

**Notable Ayurvedic Plants:**
| Plant | Q-ID | Sanskrit/Hindi Name |
|-------|------|---------------------|
| Ashwagandha | Q190224 | Withania somnifera |
| Tulsi (Holy Basil) | Q183268 | Ocimum tenuiflorum |
| Turmeric | Q42562 | Curcuma longa |
| Neem | Q213558 | Azadirachta indica |
| Amla | Q3358008 | Phyllanthus emblica |
| Brahmi | Q132490 | Bacopa monnieri |
| Guggul | Q13588216 | Commiphora wightii |
| Shatavari | Q310780 | Asparagus racemosus |
| Triphala | (compound) | Three-fruit formula |

### 3.5 Links to Specialized Databases

**IMPPAT (Indian Medicinal Plants, Phytochemistry And Therapeutics):**
- Version 2.0 contains: 4,010 Indian medicinal plants
- Phytochemicals: 17,967 compounds
- Therapeutic uses: 1,095 documented
- Plants in Ayurveda: 1,328
- Plants in Siddha: 1,151
- Wikidata property: Not yet established

**TKDL (Traditional Knowledge Digital Library):**
- Ayurveda formulations: 80,000+
- Unani formulations: 1,000,000+
- Siddha formulations: 12,000+
- Purpose: Protection against biopiracy/improper patents
- Access: Restricted to patent offices
- Wikidata property: Not yet established

### 3.6 SPARQL Queries for Ayurveda

#### Query 1: All Ayurveda-Related Items
```sparql
# Find all items related to Ayurveda
SELECT DISTINCT ?item ?itemLabel ?itemDescription WHERE {
  {
    ?item wdt:P361 wd:Q132325 .   # part of Ayurveda
  } UNION {
    ?item wdt:P1535 wd:Q132325 .  # used by Ayurveda
  } UNION {
    ?item wdt:P31 wd:Q132325 .    # instance of Ayurveda
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en,hi,sa". }
}
ORDER BY ?itemLabel
```

#### Query 2: Dosha-Related Items
```sparql
# Find items related to the three doshas
SELECT DISTINCT ?item ?itemLabel ?dosha ?doshaLabel WHERE {
  VALUES ?dosha { wd:Q1124583 wd:Q1124576 wd:Q1124572 }  # Vata, Pitta, Kapha
  {
    ?item wdt:P361 ?dosha .
  } UNION {
    ?item wdt:P31 ?dosha .
  } UNION {
    ?item wdt:P1535 ?dosha .
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en,hi,sa". }
}
ORDER BY ?doshaLabel ?itemLabel
```

#### Query 3: Ayurvedic Medicinal Plants
```sparql
# Find medicinal plants used in Ayurveda
SELECT DISTINCT ?plant ?plantLabel ?condition ?conditionLabel WHERE {
  ?plant wdt:P31/wdt:P279* wd:Q188840 .  # medicinal plant
  {
    ?plant wdt:P1535 wd:Q132325 .        # used by Ayurveda
  } UNION {
    ?plant wdt:P361 wd:Q132325 .         # part of Ayurveda
  }
  OPTIONAL { ?plant wdt:P2175 ?condition . }  # medical condition treated
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en,hi,sa". }
}
ORDER BY ?plantLabel
```

#### Query 4: Indian Medicinal Plants with Therapeutic Uses
```sparql
# Find Indian medicinal plants with their uses
SELECT DISTINCT ?plant ?plantLabel ?use ?useLabel WHERE {
  ?plant wdt:P31/wdt:P279* wd:Q188840 .  # medicinal plant
  ?plant wdt:P17 wd:Q668 .                # country: India
  ?plant wdt:P2175 ?use .                 # medical condition treated
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en,hi". }
}
ORDER BY ?plantLabel
LIMIT 500
```

#### Query 5: Ayurvedic Preparations
```sparql
# Find Ayurvedic medicinal preparations
SELECT DISTINCT ?preparation ?preparationLabel ?preparationDescription WHERE {
  {
    ?preparation wdt:P31/wdt:P279* wd:Q12140 .  # pharmaceutical drug
    ?preparation wdt:P1535 wd:Q132325 .         # used by Ayurveda
  } UNION {
    ?preparation wdt:P361 wd:Q132325 .          # part of Ayurveda
    ?preparation wdt:P31/wdt:P279* wd:Q12140 .  # is a drug
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en,hi,sa". }
}
ORDER BY ?preparationLabel
```

---

## 4. Coverage Assessment

### 4.1 Estimated Item Counts

| Category | Estimated Wikidata Items | Specialized Database Coverage |
|----------|--------------------------|------------------------------|
| **TCM** | | |
| TCM herbs | ~200-500 | TCMSP: 499 herbs |
| TCM formulas | ~50-100 | TCMID: 49,000+ prescriptions |
| TCM compounds | ~500-1000 | TCMSP: 29,384 ingredients |
| **Kampo** | | |
| Kampo formulas | ~50-100 | KampoDB: 42 formulas |
| Japanese crude drugs | ~100-200 | JP 17th ed: 176 crude drugs |
| **Ayurveda** | | |
| Ayurvedic plants | ~200-500 | IMPPAT 2.0: 4,010 plants |
| Ayurvedic preparations | ~50-100 | TKDL: 80,000+ formulations |

### 4.2 Completeness vs Specialized Databases

**TCM Coverage Gaps:**
| Specialized Database | Items | Wikidata Coverage |
|---------------------|-------|-------------------|
| TCMSP | 499 herbs, 29,384 ingredients | ~10-20% |
| TCMID | 49,000+ prescriptions | <1% |
| ETCM | 403 herbs, 3,962 prescriptions | ~10-15% |
| SymMap | TCM symptoms mapped | Limited |

**Kampo Coverage Gaps:**
| Specialized Database | Items | Wikidata Coverage |
|---------------------|-------|-------------------|
| KampoDB | 42 Kampo medicines, 54 crude drugs | ~20-30% |
| Japanese Pharmacopoeia | 148 Kampo formulations | ~10-20% |
| STORK | Standardized product info | Not linked |

**Ayurveda Coverage Gaps:**
| Specialized Database | Items | Wikidata Coverage |
|---------------------|-------|-------------------|
| IMPPAT 2.0 | 4,010 plants, 17,967 phytochemicals | ~5-10% |
| TKDL | 80,000+ Ayurveda formulations | <0.1% |

### 4.3 What's Missing in Wikidata

**General Gaps:**
1. **Standardized Properties**: No dedicated properties for traditional medicine systems
2. **Formula Composition**: Lacking structured ingredient lists for most formulas
3. **Dosage Information**: Rarely included
4. **Preparation Methods**: Not systematically captured
5. **Traditional Indications**: Often missing or incomplete
6. **Cross-system Equivalents**: Limited linking between TCM/Kampo/Ayurveda equivalents

**TCM-Specific Gaps:**
- Systematic coverage of classic formulas (e.g., Shang Han Lun formulas)
- Channel tropism (meridian affinity) data
- Temperature/nature properties (hot/cold/warm/cool)
- Taste properties (five flavors)

**Kampo-Specific Gaps:**
- Sho (pattern) classifications
- Links to modern clinical research
- Manufacturer-specific product variants

**Ayurveda-Specific Gaps:**
- Rasa (taste), Guna (quality), Virya (potency), Vipaka (post-digestive effect)
- Prabhava (special properties)
- Yoga (formulation) details

### 4.4 Cross-References to Specialized Databases

**Existing Wikidata Properties (General Medicine):**
| Property | P-code | Database |
|----------|--------|----------|
| DrugBank ID | P715 | DrugBank |
| ChEMBL ID | P592 | ChEMBL |
| PubChem CID | P662 | PubChem |
| ChEBI ID | P683 | ChEBI |
| KEGG ID | P665 | KEGG |
| IPNI plant ID | P961 | Int'l Plant Names Index |
| Plant Finder ID | P6034 | Missouri Botanical Garden |

**Missing External ID Properties (Traditional Medicine):**
| Database | Status | Proposed Use |
|----------|--------|--------------|
| TCMSP | No property | TCM herbs, compounds |
| TCMID | No property | TCM prescriptions |
| KampoDB | No property | Kampo formulas |
| STORK | No property | Kampo standardization |
| IMPPAT | No property | Indian medicinal plants |
| TKDL | No property (restricted) | Traditional formulations |
| ETCM | No property | TCM encyclopedia |

---

## 5. Bulk Extraction Strategy

### 5.1 Filtering Wikidata Dump

The Wikidata JSON dump (compressed ~100GB) can be filtered for traditional medicine items.

**Filter Criteria:**
```python
# Q-ID patterns to match for traditional medicine
TCM_QIDS = {
    'Q200253',   # Traditional Chinese Medicine
    'Q3271191',  # Chinese herbology
    'Q844401',   # Qi
    'Q178999',   # Meridian
    'Q5102267',  # Chinese herbal formula
}

KAMPO_QIDS = {
    'Q1723498',  # Kampo
    'Q6150186',  # Japanese Pharmacopoeia
}

AYURVEDA_QIDS = {
    'Q132325',   # Ayurveda
    'Q1124583',  # Vata
    'Q1124576',  # Pitta
    'Q1124572',  # Kapha
    'Q766998',   # Dosha
    'Q740784',   # Siddha medicine
}

GENERAL_QIDS = {
    'Q188840',   # Medicinal plant
    'Q735160',   # Crude drug
    'Q771035',   # Traditional medicine
    'Q207123',   # Herbal medicine
    'Q12140',    # Pharmaceutical drug
}
```

### 5.2 Properties to Extract

**Core Properties:**
```python
EXTRACT_PROPERTIES = {
    # Classification
    'P31',    # instance of
    'P279',   # subclass of
    'P361',   # part of
    'P1535',  # used by

    # Medical
    'P2175',  # medical condition treated
    'P2176',  # drug or therapy used for treatment
    'P3781',  # has active ingredient
    'P769',   # significant drug interaction

    # Chemical
    'P274',   # chemical formula
    'P231',   # CAS Registry Number
    'P662',   # PubChem CID
    'P683',   # ChEBI ID

    # Biological
    'P703',   # found in taxon
    'P171',   # parent taxon
    'P225',   # taxon name

    # External IDs
    'P715',   # DrugBank ID
    'P592',   # ChEMBL ID
    'P665',   # KEGG ID
    'P961',   # IPNI plant ID
}
```

### 5.3 Extraction Script Approach

```python
#!/usr/bin/env python3
"""
Filter Wikidata JSON dump for traditional medicine items.
"""

import json
import gzip
from pathlib import Path

# Q-IDs to track
TRADITIONAL_MEDICINE_QIDS = {
    # TCM
    'Q200253', 'Q3271191', 'Q844401', 'Q178999', 'Q5102267',
    # Kampo
    'Q1723498', 'Q6150186', 'Q11253809',
    # Ayurveda
    'Q132325', 'Q1124583', 'Q1124576', 'Q1124572', 'Q766998',
    # General
    'Q188840', 'Q735160', 'Q771035', 'Q207123', 'Q12140',
}

# Properties indicating traditional medicine relevance
RELEVANT_PROPERTIES = {'P31', 'P279', 'P361', 'P1535'}

def is_traditional_medicine_item(item):
    """Check if item is related to traditional medicine."""
    claims = item.get('claims', {})

    for prop in RELEVANT_PROPERTIES:
        if prop in claims:
            for claim in claims[prop]:
                mainsnak = claim.get('mainsnak', {})
                datavalue = mainsnak.get('datavalue', {})
                if datavalue.get('type') == 'wikibase-entityid':
                    value = datavalue.get('value', {})
                    qid = value.get('id', '')
                    if qid in TRADITIONAL_MEDICINE_QIDS:
                        return True
    return False

def extract_traditional_medicine(dump_path, output_path):
    """Extract traditional medicine items from Wikidata dump."""
    with gzip.open(dump_path, 'rt', encoding='utf-8') as f_in:
        with open(output_path, 'w', encoding='utf-8') as f_out:
            for line in f_in:
                line = line.strip().rstrip(',')
                if line in ['[', ']']:
                    continue
                try:
                    item = json.loads(line)
                    if is_traditional_medicine_item(item):
                        f_out.write(json.dumps(item) + '\n')
                except json.JSONDecodeError:
                    continue

if __name__ == '__main__':
    extract_traditional_medicine(
        'latest-all.json.gz',
        'traditional_medicine_items.jsonl'
    )
```

### 5.4 SPARQL-Based Bulk Extraction

For smaller extractions, SPARQL queries with pagination:

```sparql
# Paginated extraction of all traditional medicine items
SELECT ?item ?itemLabel ?itemDescription ?type ?typeLabel WHERE {
  {
    ?item wdt:P31|wdt:P279|wdt:P361|wdt:P1535 wd:Q200253 .
    BIND("TCM" AS ?type)
  } UNION {
    ?item wdt:P31|wdt:P279|wdt:P361|wdt:P1535 wd:Q1723498 .
    BIND("Kampo" AS ?type)
  } UNION {
    ?item wdt:P31|wdt:P279|wdt:P361|wdt:P1535 wd:Q132325 .
    BIND("Ayurveda" AS ?type)
  } UNION {
    ?item wdt:P31/wdt:P279* wd:Q188840 .
    BIND("Medicinal Plant" AS ?type)
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?item
LIMIT 10000
OFFSET 0
```

### 5.5 Recommended Extraction Workflow

1. **Initial Survey**: Run count queries to estimate volume
2. **Identify Q-IDs**: Collect all relevant Q-IDs through exploratory queries
3. **Extract Core Items**: Use SPARQL for items <100k, dump filtering for larger
4. **Enrich with Relationships**: Second pass to capture relationships
5. **Validate Coverage**: Compare against specialized databases
6. **Regular Updates**: Set up incremental extraction from recent changes

---

## 6. Recommendations

### 6.1 For Data Integration Projects

1. **Use Wikidata as a Hub**: Cross-reference with specialized databases
2. **Supplement with Specialized Sources**:
   - TCM: TCMSP, TCMID, ETCM
   - Kampo: KampoDB, STORK
   - Ayurveda: IMPPAT, TKDL (if accessible)
3. **Contribute Back**: Add missing items and properties to Wikidata

### 6.2 Proposed New Wikidata Properties

| Proposed Property | Description | Domain |
|-------------------|-------------|--------|
| TCM nature | Hot/cold/warm/cool property | TCM herbs |
| TCM flavor | Five flavors classification | TCM herbs |
| Channel tropism | Meridian affinity | TCM herbs |
| Dosha effect | Increase/decrease effect on doshas | Ayurvedic herbs |
| Kampo pattern | Sho classification | Kampo formulas |
| TCMSP ID | Identifier in TCMSP database | Chemical/herb |
| KampoDB ID | Identifier in KampoDB | Kampo items |
| IMPPAT ID | Identifier in IMPPAT database | Indian plants |

### 6.3 Query Best Practices

1. **Use Language Preferences**: Include native scripts (zh, ja, hi, sa)
2. **Apply Property Paths**: Use `/wdt:P279*` for subclass traversal
3. **Limit Results**: Always use LIMIT for exploratory queries
4. **Cache Results**: Wikidata content changes; cache with timestamps
5. **Respect Rate Limits**: Use delays between SPARQL queries

---

## References

### External Databases
- [TCMSP](https://old.tcmsp-e.com/tcmsp.php) - Traditional Chinese Medicine Systems Pharmacology
- [TCMID](http://www.megabionet.org/tcmid/) - Traditional Chinese Medicine Integrated Database
- [ETCM](http://www.tcmip.cn/ETCM/) - Encyclopedia of Traditional Chinese Medicine
- [KampoDB](http://wakanmoview.inm.u-toyama.ac.jp/kampo/) - Kampo Medicine Database
- [STORK](http://mpdb.nibiohn.go.jp/stork/) - Standards of Reporting Kampo Products
- [IMPPAT](https://cb.imsc.res.in/imppat/) - Indian Medicinal Plants Database

### Wikidata Resources
- [Wikidata SPARQL Query Service](https://query.wikidata.org/)
- [Wikidata SPARQL Tutorial](https://www.wikidata.org/wiki/Wikidata:SPARQL_tutorial)
- [Wikidata SPARQL Examples](https://www.wikidata.org/wiki/Wikidata:SPARQL_query_service/queries/examples)
- [WikiProject Medicine](https://www.wikidata.org/wiki/Wikidata:WikiProject_Medicine)

### Research Literature
- KampoDB paper: Sawada et al., Scientific Reports 8:11216 (2018)
- IMPPAT paper: Mohanraj et al., Scientific Reports 8:4329 (2018)
- Wikidata for biomedical data: Turki et al., J Biomed Inform (2019)
