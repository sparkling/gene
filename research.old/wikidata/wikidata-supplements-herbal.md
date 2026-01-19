# Wikidata Coverage: Western Herbal Medicine, Supplements, Vitamins, and Minerals

**Research Date:** 2026-01-18
**Data Source:** Wikidata SPARQL Endpoint (query.wikidata.org)

## Executive Summary

This document analyzes Wikidata's coverage of Western herbal medicine, dietary supplements, vitamins, minerals, and natural products. While Wikidata has extensive chemical compound data (1.3M+ compounds with structural identifiers), specialized classification for supplements and herbal medicine categories is relatively sparse.

---

## 1. Western Herbal Medicine

### Key Q-IDs
| Entity | Q-ID | Description |
|--------|------|-------------|
| Herbal medicine | Q188748 | The use of plants for medicinal purposes |
| Medicinal plant | Q2140699 | Plant used for medicinal purposes |

### Coverage Statistics

| Category | Count | Notes |
|----------|-------|-------|
| Items classified as herbal medicine (P31:Q188748) | 18 | Very limited direct instances |
| Medicinal plants (P31:Q2140699) | 296 | Direct instances only |
| Medicinal plants (with subclasses) | ~546 | Including inherited classifications |
| Medicinal plants with P2175 (treats) property | 0 | Gap - no treatment links |

### Relevant Properties

| Property | ID | Usage Count | Description |
|----------|-----|-------------|-------------|
| medical condition treated | P2175 | 7,038 | Indicates therapeutic use |
| has active ingredient | P3781 | 2,905 | Links products to compounds |
| found in taxon | P703 | 2,572,326 | Links compounds to source organisms |

### SPARQL Query: All Medicinal Plants

```sparql
# All medicinal plants with their scientific names and therapeutic uses
SELECT DISTINCT ?plant ?plantLabel ?scientificName ?treats ?treatsLabel
WHERE {
  ?plant wdt:P31/wdt:P279* wd:Q2140699 .  # Instance of medicinal plant (or subclass)

  OPTIONAL { ?plant wdt:P225 ?scientificName . }  # Scientific name
  OPTIONAL { ?plant wdt:P2175 ?treats . }         # Medical condition treated

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
ORDER BY ?plantLabel
LIMIT 1000
```

### SPARQL Query: Herbal Preparations

```sparql
# Herbal medicine preparations and products
SELECT ?item ?itemLabel ?ingredientLabel ?treatsLabel
WHERE {
  ?item wdt:P31/wdt:P279* wd:Q188748 .  # Instance of herbal medicine

  OPTIONAL { ?item wdt:P3781 ?ingredient . }  # Active ingredient
  OPTIONAL { ?item wdt:P2175 ?treats . }      # Treats condition

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

### Gap Analysis - Herbal Medicine
- **Major Gap:** Only 18 items directly classified as herbal medicine
- **Missing:** Links between medicinal plants and conditions treated (P2175)
- **Missing:** Herbal preparation formulations
- **Opportunity:** Many compounds have P703 (found in taxon) which could link to medicinal plants

---

## 2. Dietary Supplements

### Key Q-IDs
| Entity | Q-ID | Description |
|--------|------|-------------|
| Dietary supplement | Q179478 | Product intended to supplement the diet |

### Coverage Statistics

| Category | Count | Notes |
|----------|-------|-------|
| Direct instances of dietary supplement | 0 | No items directly classified |
| Items with DSLD ID (P8743) | 18 | Links to NIH DSLD database |
| Items with LNHPD ID (P6955) | 88 | Links to Health Canada database |

### External Database Links

| Database | Property | Count | Database Size |
|----------|----------|-------|---------------|
| DSLD (Dietary Supplement Label Database) | P8743 | 18 | ~100,000+ products |
| Health Canada LNHPD | P6955 | 88 | ~90,000+ products |

### SPARQL Query: Supplements with External IDs

```sparql
# Items with dietary supplement database identifiers
SELECT ?item ?itemLabel ?dsld ?lnhpd
WHERE {
  { ?item wdt:P8743 ?dsld . }  # DSLD ID
  UNION
  { ?item wdt:P6955 ?lnhpd . } # LNHPD ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

### Gap Analysis - Dietary Supplements
- **Critical Gap:** No items directly classified as dietary supplement
- **Minimal Coverage:** Only 18 items linked to DSLD, 88 to LNHPD
- **Comparison:** DSLD contains 100,000+ products; Wikidata links <0.02%
- **Opportunity:** Import structured data from DSLD/LNHPD

---

## 3. Vitamins

### Key Q-IDs
| Entity | Q-ID | Description |
|--------|------|-------------|
| Vitamin | Q34956 | Organic compound essential in small quantities |
| Vitamin A | Q18225 | Retinoids and carotenoids |
| Vitamin B | Q183206 | B-complex vitamins |
| Vitamin C | Q140829 | Ascorbic acid |
| Vitamin D | Q175621 | Calciferols |
| Vitamin E | Q141180 | Tocopherols and tocotrienols |
| Vitamin K | Q182338 | Phylloquinone and menaquinones |
| Provitamin | Q423960 | Vitamin precursor |

### Vitamin Classification Hierarchy

```
Vitamin (Q34956)
├── Fat-soluble vitamin (Q4389811)
│   ├── Vitamin A (Q18225)
│   ├── Vitamin D (Q175621)
│   ├── Vitamin E (Q141180)
│   └── Vitamin K (Q182338)
├── Water-soluble vitamin (Q11548920)
│   └── Vitamin B (Q183206)
└── Provitamin (Q423960)
```

### Vitamin Forms Coverage

| Vitamin | Forms in Wikidata | Examples |
|---------|-------------------|----------|
| Vitamin A | 15+ | Retinol, retinal, retinoic acid, beta-carotene |
| Vitamin B | 7+ | Biotin, B12, thiamine, pantothenate, folic acid |
| Vitamin C | ~3 | Ascorbic acid, vitamin C derivatives |
| Vitamin D | 10+ | D1-D7, ergocalciferol, cholecalciferol, calcitriol |
| Vitamin E | 8+ | Alpha/beta/gamma/delta-tocopherol, tocotrienols |
| Vitamin K | 10+ | K1 (phytonadione), K2 (menaquinones), K3-K7 |

### Vitamins with Treatment Properties (P2175)

| Vitamin | Conditions Treated (per Wikidata) |
|---------|-----------------------------------|
| Cholecalciferol (D3) | Cancer, osteoporosis, cardiovascular disease, chronic renal insufficiency, skin diseases |
| Vitamin K | Vitamin K deficiency |
| Vitamin A | (none linked) |
| Vitamin E | (none linked) |
| Vitamin C | (none linked) |

### SPARQL Query: All Vitamins with Biological Data

```sparql
# All vitamins with targets, processes, and therapeutic uses
SELECT ?vitamin ?vitaminLabel ?treatsLabel ?processLabel ?targetLabel
WHERE {
  ?vitamin wdt:P279*/wdt:P31* wd:Q34956 .  # Vitamin or subclass

  OPTIONAL { ?vitamin wdt:P2175 ?treats . }   # Medical condition treated
  OPTIONAL { ?vitamin wdt:P682 ?process . }   # Biological process
  OPTIONAL { ?vitamin wdt:P129 ?target . }    # Physically interacts with

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
ORDER BY ?vitaminLabel
```

### SPARQL Query: B Vitamins Complete

```sparql
# All B vitamins with their forms
SELECT ?vitamin ?vitaminLabel ?form ?formLabel
WHERE {
  ?vitamin wdt:P279* wd:Q183206 .  # Subclass of vitamin B
  OPTIONAL { ?form wdt:P279 ?vitamin . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

### Gap Analysis - Vitamins
- **Good:** Basic vitamin hierarchy is well-structured
- **Gap:** Only cholecalciferol and vitamin K have P2175 treatment links
- **Gap:** Most vitamins lack P682 (biological process) annotations
- **Gap:** Vitamin-drug interactions not systematically captured

---

## 4. Minerals (Nutritional)

### Key Q-IDs
| Entity | Q-ID | Description |
|--------|------|-------------|
| Dietary mineral | Q15123592 | Mineral required for nutrition |

### Essential Minerals (Expected Coverage)

| Mineral | Element Q-ID | Dietary Role |
|---------|--------------|--------------|
| Calcium | Q706 | Bone health, nerve function |
| Magnesium | Q660 | Enzyme function, muscle |
| Potassium | Q703 | Electrolyte balance |
| Sodium | Q658 | Electrolyte balance |
| Iron | Q677 | Oxygen transport |
| Zinc | Q758 | Immune function, enzymes |
| Copper | Q753 | Enzyme cofactor |
| Manganese | Q731 | Enzyme function |
| Selenium | Q876 | Antioxidant |
| Iodine | Q1103 | Thyroid function |
| Chromium | Q725 | Glucose metabolism |
| Molybdenum | Q1098 | Enzyme cofactor |

### Coverage Statistics

| Category | Count | Notes |
|----------|-------|-------|
| Direct instances of dietary mineral (Q15123592) | 0 | No items classified |
| Chemical elements in Wikidata | 118 | All elements present |
| Elements with P2175 (treats) | ~0 | Minerals lack treatment links |

### SPARQL Query: Nutritional Minerals

```sparql
# Essential minerals with biological roles
SELECT ?element ?elementLabel ?deficiency ?deficiencyLabel ?process ?processLabel
WHERE {
  VALUES ?element {
    wd:Q706   # Calcium
    wd:Q660   # Magnesium
    wd:Q703   # Potassium
    wd:Q658   # Sodium
    wd:Q677   # Iron
    wd:Q758   # Zinc
    wd:Q753   # Copper
    wd:Q731   # Manganese
    wd:Q876   # Selenium
    wd:Q1103  # Iodine
    wd:Q725   # Chromium
    wd:Q1098  # Molybdenum
  }

  OPTIONAL { ?element wdt:P2175 ?deficiency . }  # Treats (deficiency diseases)
  OPTIONAL { ?element wdt:P682 ?process . }      # Biological process

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

### Gap Analysis - Minerals
- **Critical Gap:** No items classified as dietary mineral (Q15123592)
- **Gap:** Elements lack nutritional role annotations
- **Gap:** No links between minerals and deficiency diseases
- **Opportunity:** Add P2868 (subject has role) for nutritional roles

---

## 5. Natural Products

### Key Q-IDs
| Entity | Q-ID | Description |
|--------|------|-------------|
| Natural product | Q193166 | Chemical compound produced by living organism |
| Secondary metabolite | Q81915 | Organic compound not directly involved in growth |
| Alkaloid | Q422204 | Nitrogen-containing organic compound |
| Flavonoid | Q182990 | Polyphenolic secondary metabolite |
| Terpenoid | Q407550 | Terpene-derived compound |

### Coverage Statistics

| Category | Count | Notes |
|----------|-------|-------|
| Direct instances of natural product | 0 | Class not used for instances |
| Compounds with P703 (found in taxon) | 2,572,326 | Rich source organism data |
| Compounds in plants with SMILES | 227,315 | Plant-derived with structures |
| Items with Natural Products Atlas ID | 932 | Links to NPA database |
| Items with LOTUS ID (P11802) | 51 | Links to LOTUS database |
| Items with ChEBI ID | 157,189 | Links to ChEBI ontology |
| Secondary metabolites | 1 | Extremely sparse |

### Phytochemical Classes

| Class | Q-ID | Count in Hierarchy | Notes |
|-------|------|-------------------|-------|
| Alkaloid | Q422204 | 1 | Sparse hierarchy |
| Flavonoid | Q182990 | 6 | Limited coverage |
| Terpenoid | Q407550 | 1 | Sparse hierarchy |

### Chemical Identifier Coverage

| Property | ID | Count | Description |
|----------|-----|-------|-------------|
| Canonical SMILES | P233 | 1,364,830 | Structural representation |
| InChI | P234 | 1,364,612 | Standard identifier |
| InChIKey | P235 | 1,365,873 | Hashed identifier |
| CAS Registry Number | P231 | 945,081 | Chemical registry |
| PubChem CID | P662 | 1,329,508 | PubChem compound ID |

### SPARQL Query: Plant Compounds with Structures

```sparql
# Plant-derived compounds with chemical structures
SELECT ?compound ?compoundLabel ?taxon ?taxonLabel ?smiles ?inchikey
WHERE {
  ?compound wdt:P703 ?taxon .           # Found in taxon
  ?taxon wdt:P31/wdt:P279* wd:Q756 .   # Taxon is a plant

  OPTIONAL { ?compound wdt:P233 ?smiles . }
  OPTIONAL { ?compound wdt:P235 ?inchikey . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
LIMIT 1000
```

### SPARQL Query: Natural Products with Database Links

```sparql
# Natural products with external database identifiers
SELECT ?compound ?compoundLabel ?npa ?lotus ?chebi
WHERE {
  { ?compound wdt:P8121 ?npa . }   # Natural Products Atlas
  UNION
  { ?compound wdt:P11802 ?lotus . } # LOTUS
  UNION
  { ?compound wdt:P683 ?chebi . }   # ChEBI

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
LIMIT 500
```

### Gap Analysis - Natural Products
- **Good:** 2.5M+ compounds with source organism (P703)
- **Good:** 227K plant compounds with structural data
- **Gap:** Q193166 (natural product) not used as instance class
- **Gap:** Secondary metabolite classification nearly empty
- **Gap:** Only 932 links to Natural Products Atlas (~33K compounds)
- **Gap:** Only 51 links to LOTUS (~800K compound-taxon pairs)

---

## 6. Coverage Statistics Summary

### Overall Counts

| Category | Wikidata Count | Reference Database | Coverage % |
|----------|----------------|-------------------|------------|
| Medicinal plants | 296-546 | ~28,000 (MPNS) | 1-2% |
| Herbal medicines | 18 | ~10,000+ (traditional) | <1% |
| Dietary supplements | 0 | ~100,000 (DSLD) | 0% |
| Vitamins (distinct) | ~50 | ~13 main vitamins | Good |
| Vitamin forms | ~60 | ~100+ forms | 60% |
| Dietary minerals | 0 (as class) | 12 essential | Class unused |
| Natural products (NPA) | 932 | ~33,000 | 3% |
| LOTUS compounds | 51 | ~800,000 | <0.01% |
| Plant compounds (P703) | 227,315 | - | Rich data |

### Property Usage Statistics

| Property | Count | Description |
|----------|-------|-------------|
| P2175 (treats) | 7,038 | Therapeutic relationships |
| P3781 (active ingredient) | 2,905 | Product-compound links |
| P703 (found in taxon) | 2,572,326 | Source organism |
| P682 (biological process) | 365,262 | Process involvement |
| P129 (interacts with) | 4,264 | Molecular interactions |
| P233 (SMILES) | 1,364,830 | Chemical structures |
| P234 (InChI) | 1,364,612 | Standard identifiers |
| P235 (InChIKey) | 1,365,873 | Hashed identifiers |
| P231 (CAS) | 945,081 | Registry numbers |
| P662 (PubChem CID) | 1,329,508 | PubChem links |
| P683 (ChEBI) | 157,189 | ChEBI ontology |
| P2868 (has role) | 221,377 | Chemical roles |

---

## 7. Key Findings and Recommendations

### Strengths
1. **Chemical Structure Data:** 1.3M+ compounds with SMILES/InChI/InChIKey
2. **Source Organism Links:** 2.5M+ P703 (found in taxon) statements
3. **Vitamin Hierarchy:** Well-structured classification of major vitamins
4. **External Identifiers:** Strong links to PubChem, ChEBI, CAS

### Critical Gaps
1. **Dietary Supplement Class (Q179478):** Not used - 0 instances
2. **Dietary Mineral Class (Q15123592):** Not used - 0 instances
3. **Natural Product Class (Q193166):** Not used for instances
4. **Medicinal Plant Therapeutic Links:** 0 plants have P2175 (treats)
5. **LOTUS Integration:** Only 51 of 800K+ compound-taxon pairs

### Recommendations for Data Improvement

1. **Import LOTUS Data**
   - LOTUS has 800K+ compound-taxon pairs
   - Add P11802 (LOTUS ID) systematically
   - Enriches P703 (found in taxon) significantly

2. **Add Dietary Supplement Products**
   - Use Q179478 for supplement products
   - Link to DSLD (P8743) and LNHPD (P6955)
   - Add P3781 (active ingredient) relationships

3. **Classify Dietary Minerals**
   - Add P31:Q15123592 to essential mineral elements
   - Add P2868 (has role) for nutritional roles
   - Link to deficiency diseases via P2175

4. **Enrich Vitamin Data**
   - Add P2175 (treats) for vitamin deficiency diseases
   - Add P682 (biological process) for metabolic roles
   - Add P129 (interacts with) for molecular targets

5. **Link Medicinal Plants to Conditions**
   - Many plants have traditional uses documented
   - Add P2175 (treats) based on pharmacological evidence
   - Link plant compounds to therapeutic effects

---

## 8. Useful SPARQL Queries

### Count All Items with Therapeutic Uses

```sparql
SELECT (COUNT(DISTINCT ?item) AS ?count)
WHERE {
  ?item wdt:P2175 ?condition .
}
```

### Find Compounds in Specific Plant

```sparql
# Example: Compounds in Ginkgo biloba
SELECT ?compound ?compoundLabel ?smiles
WHERE {
  ?compound wdt:P703 wd:Q43284 .  # Found in Ginkgo biloba
  OPTIONAL { ?compound wdt:P233 ?smiles . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

### Vitamins Missing Treatment Data

```sparql
# Vitamins without P2175 (treats) property
SELECT ?vitamin ?vitaminLabel
WHERE {
  ?vitamin wdt:P279* wd:Q34956 .  # Subclass of vitamin
  FILTER NOT EXISTS { ?vitamin wdt:P2175 ?treats . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

### Plants with Bioactive Compounds

```sparql
# Plants containing compounds with known therapeutic uses
SELECT DISTINCT ?plant ?plantLabel ?compound ?compoundLabel ?conditionLabel
WHERE {
  ?compound wdt:P703 ?plant .
  ?plant wdt:P31/wdt:P279* wd:Q756 .  # Plant
  ?compound wdt:P2175 ?condition .

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
LIMIT 500
```

---

## 9. External Database Comparison

| Database | Focus | Size | Wikidata Links |
|----------|-------|------|----------------|
| DSLD | US Supplement Labels | ~100K products | P8743 (18 items) |
| LNHPD | Canadian Natural Health | ~90K products | P6955 (88 items) |
| LOTUS | Natural Product-Taxon | 800K+ pairs | P11802 (51 items) |
| Natural Products Atlas | Microbial Natural Products | ~33K compounds | P8121 (932 items) |
| ChEBI | Chemical Ontology | ~60K entities | P683 (157K items) |
| PubChem | Chemical Database | ~117M compounds | P662 (1.3M items) |
| MPNS | Medicinal Plant Names | ~28K species | No direct property |

---

## 10. Conclusion

Wikidata has **excellent chemical compound coverage** with 1.3M+ structures and strong links to PubChem and ChEBI. However, **specialized classification for supplements, herbal medicine, and dietary minerals is severely lacking**. The key opportunity is leveraging the existing compound-taxon relationships (2.5M+ P703 statements) to build out therapeutic and nutritional annotations.

Priority improvements:
1. Integrate LOTUS for comprehensive plant compound coverage
2. Populate Q179478 (dietary supplement) with products
3. Add nutritional role annotations to mineral elements
4. Link medicinal plants to conditions they treat

The foundation exists for a comprehensive natural products knowledge base; systematic curation efforts could transform Wikidata into a valuable resource for herbal medicine and supplement research.
