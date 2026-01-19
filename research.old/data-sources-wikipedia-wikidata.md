# Wikipedia and Wikidata Data Sources for Gene Platform

**Research Date:** 2026-01-19
**Purpose:** Comprehensive reference for Wikipedia/Wikidata/DBpedia data sources for genetics, health conditions, and supplements

---

## Table of Contents

1. [Wikipedia Categories for Genetics, Health, and Supplements](#1-wikipedia-categories)
2. [Wikidata SPARQL Queries for Biomedical Entities](#2-wikidata-sparql-queries)
3. [DBpedia Structured Extracts](#3-dbpedia-structured-extracts)
4. [Missing Category Coverage Gaps](#4-coverage-gap-analysis)
5. [Infobox Data Extraction Opportunities](#5-infobox-data-extraction)
6. [Cross-Linking with Scientific Databases](#6-cross-linking-scientific-databases)
7. [Recommended Extraction Strategy](#7-recommended-extraction-strategy)

---

## 1. Wikipedia Categories

### 1.1 Genetics Categories

**Primary Category:** [Category:Genetics](https://en.wikipedia.org/wiki/Category:Genetics)

| Subcategory | Est. Pages | Relevance | Description |
|-------------|------------|-----------|-------------|
| Genetic disorders | 1,000+ | Critical | All genetic conditions |
| Genes | 500+ | Critical | Gene articles by organism |
| Genomics | 300+ | High | Genome-level analysis |
| DNA | 200+ | High | DNA structure, replication |
| Chromosomes | 150+ | High | Chromosome biology |
| Human genetics | 400+ | Critical | Human-specific genetics |
| Medical genetics | 600+ | Critical | Clinical genetics |
| Molecular genetics | 250+ | High | Molecular mechanisms |
| Population genetics | 100+ | Medium | Population-level variation |
| Pharmacogenomics | 50+ | Critical | Drug-gene interactions |
| Epigenetics | 100+ | High | Epigenetic modifications |
| Gene expression | 150+ | High | Transcription, regulation |
| Genetic mapping | 75+ | Medium | Linkage, association |
| Genetic engineering | 200+ | Medium | Gene editing, GMOs |

**Category Hierarchy:**
```
Category:Genetics
├── Category:Applied genetics
│   ├── Category:Genetic engineering
│   └── Category:Pharmacogenomics
├── Category:Classical genetics
│   └── Category:Mendelian genetics
├── Category:Molecular genetics
│   ├── Category:DNA
│   ├── Category:Gene expression
│   └── Category:RNA
├── Category:Human genetics
│   ├── Category:Medical genetics
│   └── Category:Genetic disorders
├── Category:Genomics
│   ├── Category:Comparative genomics
│   └── Category:Functional genomics
└── Category:Genes
    ├── Category:Human genes by chromosome
    └── Category:Genes by protein encoded
```

### 1.2 Health Conditions Categories

**Primary Category:** [Category:Diseases and disorders](https://en.wikipedia.org/wiki/Category:Diseases_and_disorders)

| Category | Subcategories | Pages | Description |
|----------|---------------|-------|-------------|
| [Rare diseases](https://en.wikipedia.org/wiki/Category:Rare_diseases) | 19 | ~823 | Orphan diseases |
| [Genetic disorders](https://en.wikipedia.org/wiki/Category:Genetic_disorders) | 25+ | ~1,200 | Inherited conditions |
| [Rare syndromes](https://en.wikipedia.org/wiki/Category:Rare_syndromes) | 3 | ~380 | Syndromic conditions |
| [Rare genetic syndromes](https://en.wikipedia.org/wiki/Category:Rare_genetic_syndromes) | 1 | ~185 | Genetic syndromes |
| Autoimmune diseases | 10+ | ~300 | Immune disorders |
| Metabolic disorders | 15+ | ~500 | Metabolic conditions |
| Neurological disorders | 20+ | ~800 | Nervous system |
| Cardiovascular diseases | 15+ | ~600 | Heart/vascular |
| Cancer | 30+ | ~1,500 | Neoplasms |
| Infectious diseases | 25+ | ~1,000 | Pathogen-caused |

**Medical Condition Category Structure:**
```
Category:Diseases and disorders
├── Category:Diseases and disorders by body system
│   ├── Category:Neurological disorders
│   ├── Category:Cardiovascular diseases
│   ├── Category:Respiratory diseases
│   └── Category:Digestive diseases
├── Category:Genetic disorders
│   ├── Category:Autosomal dominant disorders
│   ├── Category:Autosomal recessive disorders
│   ├── Category:X-linked recessive disorders
│   └── Category:Mitochondrial diseases
├── Category:Rare diseases
│   ├── Category:Rare syndromes
│   └── Category:Orphan diseases
└── Category:Cancer
    ├── Category:Carcinomas
    └── Category:Sarcomas
```

### 1.3 Supplements Categories

**Primary Category:** [Category:Dietary supplements](https://en.wikipedia.org/wiki/Category:Dietary_supplements)

| Subcategory | Pages | Description |
|-------------|-------|-------------|
| [Vitamins](https://en.wikipedia.org/wiki/Category:Vitamins) | ~26 | Vitamin articles |
| [B vitamins](https://en.wikipedia.org/wiki/Category:B_vitamins) | ~21 | B-complex vitamins |
| Minerals (nutrition) | ~30 | Dietary minerals |
| Amino acid supplements | ~20 | Amino acids |
| Bodybuilding supplements | ~40 | Sports nutrition |
| Effervescent supplements | ~10 | Effervescent forms |
| Nutritional supplement companies | ~50 | Industry companies |

**Vitamin Category Hierarchy:**
```
Category:Vitamins (10 subcategories, 26 pages)
├── Category:B vitamins (2 subcategories, 21 pages)
│   ├── Thiamine (B1)
│   ├── Riboflavin (B2)
│   ├── Niacin (B3)
│   ├── Pantothenic acid (B5)
│   ├── Pyridoxine (B6)
│   ├── Biotin (B7)
│   ├── Folate (B9)
│   └── Cobalamin (B12)
├── Category:Vitamin A
├── Category:Vitamin C
├── Category:Vitamin D
├── Category:Vitamin E
└── Category:Vitamin K
```

### 1.4 Alternative Medicine Categories

**Primary Category:** [Category:Alternative medicine](https://en.wikipedia.org/wiki/Category:Alternative_medicine)

| Subcategory | Pages | Description |
|-------------|-------|-------------|
| Herbal medicine | ~150 | Plant-based medicine |
| Traditional Chinese medicine | ~200 | TCM practices |
| Ayurveda | ~100 | Indian medicine |
| Homeopathy | ~75 | Homeopathic remedies |
| Naturopathy | ~50 | Natural therapies |
| Acupuncture | ~40 | Acupuncture points |

**Source:** [Wikipedia:WikiProject Alternative medicine](https://en.wikipedia.org/wiki/Wikipedia:WikiProject_Alternative_medicine)

---

## 2. Wikidata SPARQL Queries

### 2.1 SPARQL Endpoint

**Primary Endpoint:** https://query.wikidata.org/sparql
**GUI Interface:** https://query.wikidata.org/

### 2.2 Biomedical Entity Coverage Statistics

| Entity Type | Wikidata Count | Source Database | Coverage |
|-------------|----------------|-----------------|----------|
| Human genes | ~59,721 | NCBI Gene | ~100% |
| Mouse genes | ~73,355 | NCBI Gene | ~100% |
| Human proteins | ~27,306 | UniProt SwissProt | ~100% |
| Mouse proteins | ~16,728 | UniProt SwissProt | ~100% |
| Chemical compounds | ~1,300,000+ | PubChem | ~1% |
| Medications | ~45,000 | DrugBank/ChEMBL | ~50% |
| Diseases | ~200,000 | Disease Ontology | ~90% |
| Genetic associations | ~12,593 | Various | Partial |
| Biological pathways | ~3,000 | Reactome/WikiPathways | ~80% |

**Sources:**
- [Wikidata as a knowledge graph for the life sciences](https://pmc.ncbi.nlm.nih.gov/articles/PMC7077981/)
- [Wikidata as a semantic framework for the Gene Wiki initiative](https://pmc.ncbi.nlm.nih.gov/articles/PMC4795929/)

### 2.3 Essential SPARQL Queries for Gene Platform

#### Query 1: All Human Genes with Identifiers

```sparql
# All human protein-coding genes with major identifiers
SELECT DISTINCT ?gene ?geneLabel ?entrezId ?hgncSymbol ?hgncId ?ensemblId
WHERE {
  ?gene wdt:P31 wd:Q20747295 .              # instance of protein-coding gene
  ?gene wdt:P703 wd:Q15978631 .             # found in taxon: Homo sapiens

  OPTIONAL { ?gene wdt:P351 ?entrezId . }    # Entrez Gene ID
  OPTIONAL { ?gene wdt:P353 ?hgncSymbol . }  # HGNC gene symbol
  OPTIONAL { ?gene wdt:P354 ?hgncId . }      # HGNC ID
  OPTIONAL { ?gene wdt:P594 ?ensemblId . }   # Ensembl gene ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?hgncSymbol
```

#### Query 2: Gene-Disease Associations

```sparql
# Genes associated with diseases (using P2293 genetic association)
SELECT DISTINCT ?gene ?geneLabel ?hgncSymbol ?disease ?diseaseLabel ?omimId
WHERE {
  ?gene wdt:P31/wdt:P279* wd:Q7187 .        # instance of gene
  ?gene wdt:P703 wd:Q15978631 .              # Homo sapiens
  ?gene wdt:P353 ?hgncSymbol .               # HGNC symbol

  ?disease wdt:P2293 ?gene .                 # genetic association to gene
  OPTIONAL { ?disease wdt:P492 ?omimId . }   # OMIM ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?diseaseLabel ?geneLabel
LIMIT 1000
```

#### Query 3: Drugs and Their Protein Targets

```sparql
# Drugs with their protein targets and pathways
SELECT DISTINCT ?drug ?drugLabel ?target ?targetLabel ?pathway ?pathwayLabel
WHERE {
  ?drug wdt:P31 wd:Q12140 .                  # instance of medication
  ?drug wdt:P129 ?target .                   # drug interacts with target
  ?target wdt:P31 wd:Q8054 .                 # target is protein
  ?target wdt:P703 wd:Q15978631 .            # human protein

  OPTIONAL {
    ?pathway wdt:P31/wdt:P279* wd:Q4915012 . # biological pathway
    ?pathway wdt:P527 ?target .               # pathway has part target
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?drugLabel
LIMIT 500
```

#### Query 4: Vitamins with Biological Data

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

#### Query 5: Plant Compounds with Structures

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

#### Query 6: CYP450 Enzymes (Pharmacogenomics)

```sparql
# CYP450 enzymes in humans (pharmacogenomics relevant)
SELECT DISTINCT ?enzyme ?enzymeLabel ?ecNumber ?uniprotId ?geneLabel ?hgncSymbol
WHERE {
  ?enzyme wdt:P31/wdt:P279* wd:Q423111 .    # subclass of cytochrome P450
  ?enzyme wdt:P703 wd:Q15978631 .            # Homo sapiens

  OPTIONAL { ?enzyme wdt:P591 ?ecNumber . }  # EC number
  OPTIONAL { ?enzyme wdt:P352 ?uniprotId . } # UniProt ID
  OPTIONAL {
    ?enzyme wdt:P702 ?gene .                 # encoded by
    ?gene wdt:P353 ?hgncSymbol .             # HGNC symbol
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
ORDER BY ?enzymeLabel
```

#### Query 7: Dietary Supplements with Database Links

```sparql
# Items with dietary supplement database identifiers
SELECT ?item ?itemLabel ?dsld ?lnhpd ?chebi
WHERE {
  { ?item wdt:P8743 ?dsld . }   # DSLD ID
  UNION
  { ?item wdt:P6955 ?lnhpd . }  # LNHPD ID
  UNION
  {
    ?item wdt:P31/wdt:P279* wd:Q324546 .  # dietary supplement
    ?item wdt:P683 ?chebi .               # ChEBI ID
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

#### Query 8: Rare Genetic Diseases

```sparql
# Rare genetic diseases with OMIM and Orphanet identifiers
SELECT DISTINCT ?disease ?diseaseLabel ?omimId ?orphanetId ?geneLabel
WHERE {
  ?disease wdt:P31/wdt:P279* wd:Q929833 .   # rare disease
  ?disease wdt:P31/wdt:P279* wd:Q18553442 . # genetic disorder

  OPTIONAL { ?disease wdt:P492 ?omimId . }     # OMIM ID
  OPTIONAL { ?disease wdt:P1550 ?orphanetId . } # Orphanet ID
  OPTIONAL { ?disease wdt:P2293 ?gene . ?gene rdfs:label ?geneLabel . FILTER(LANG(?geneLabel) = "en") }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
ORDER BY ?diseaseLabel
```

### 2.4 Federated Queries

#### Federated Query: Wikidata + WikiPathways

```sparql
# Federated query combining Wikidata metabolites with WikiPathways
# Run on WikiPathways endpoint: https://sparql.wikipathways.org/sparql
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX wd: <http://www.wikidata.org/entity/>

SELECT DISTINCT ?metabolite ?metaboliteLabel ?pathway ?pathwayTitle
WHERE {
  # WikiPathways part
  ?pathway a wp:Pathway .
  ?pathway wp:organismName "Homo sapiens" .
  ?pathway dc:title ?pathwayTitle .
  ?pathway wp:isAbout ?pwDataNode .
  ?pwDataNode wp:bdbWikidata ?metabolite .

  # Wikidata part for labels
  SERVICE <https://query.wikidata.org/sparql> {
    ?metabolite rdfs:label ?metaboliteLabel .
    FILTER(LANG(?metaboliteLabel) = "en")
  }
}
LIMIT 200
```

**Available Federated Endpoints:**

| Endpoint | URL | Description |
|----------|-----|-------------|
| Wikidata | https://query.wikidata.org/sparql | Main Wikidata |
| WikiPathways | https://sparql.wikipathways.org/sparql | Pathway data |
| UniProt | https://sparql.uniprot.org/sparql | Protein data |
| ChEMBL Mirror | https://chemblmirror.rdf.bigcat-bioinformatics.org/sparql | Drug bioactivity |
| Rhea | https://sparql.rhea-db.org/sparql | Biochemical reactions |
| neXtProt | https://sparql.nextprot.org/sparql | Human protein data |

**Source:** [WikiPathways Federated SPARQL queries](https://www.wikipathways.org/federated.html)

---

## 3. DBpedia Structured Extracts

### 3.1 DBpedia SPARQL Endpoint

**Primary Endpoint:** http://dbpedia.org/sparql
**Interface:** [DBpedia SPARQL Query Editor](https://dbpedia.org/sparql)

### 3.2 Ontology Classes for Biomedical Data

| Class | URI | Est. Instances | Description |
|-------|-----|----------------|-------------|
| Drug | dbo:Drug | ~10,000 | Pharmaceutical compounds |
| Disease | dbo:Disease | ~5,000 | Medical conditions |
| Gene | dbo:Gene | ~500 | Gene entities |
| HumanGene | dbo:HumanGene | ~200 | Human-specific genes |
| Protein | dbo:Protein | ~1,000 | Protein entities |
| ChemicalCompound | dbo:ChemicalCompound | ~50,000 | Chemical substances |
| AnatomicalStructure | dbo:AnatomicalStructure | ~3,000 | Body structures |
| Enzyme | dbo:Enzyme | ~500 | Enzyme proteins |

**Source:** [DBpedia Ontology Classes](https://dief.tools.dbpedia.org/server/ontology/classes/)

### 3.3 DBpedia SPARQL Queries

#### Query 1: All Drugs with Wikipedia Links

```sparql
PREFIX dbo: <http://dbpedia.org/ontology/>
PREFIX dbr: <http://dbpedia.org/resource/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT DISTINCT ?drug ?drugLabel ?drugClass ?wikipedia
WHERE {
  ?drug a dbo:Drug .
  ?drug rdfs:label ?drugLabel .
  FILTER(LANG(?drugLabel) = "en")

  OPTIONAL { ?drug dbo:drugClass ?drugClass . }
  OPTIONAL { ?drug foaf:isPrimaryTopicOf ?wikipedia . }
}
ORDER BY ?drugLabel
LIMIT 500
```

#### Query 2: Diseases with Genes and Symptoms

```sparql
PREFIX dbo: <http://dbpedia.org/ontology/>
PREFIX dbp: <http://dbpedia.org/property/>

SELECT DISTINCT ?disease ?diseaseLabel ?gene ?symptom
WHERE {
  ?disease a dbo:Disease .
  ?disease rdfs:label ?diseaseLabel .
  FILTER(LANG(?diseaseLabel) = "en")

  OPTIONAL { ?disease dbo:gene ?gene . }
  OPTIONAL { ?disease dbo:symptom ?symptom . }
}
ORDER BY ?diseaseLabel
LIMIT 500
```

#### Query 3: Chemical Compounds with Identifiers

```sparql
PREFIX dbo: <http://dbpedia.org/ontology/>
PREFIX dbp: <http://dbpedia.org/property/>

SELECT DISTINCT ?compound ?compoundLabel ?casNumber ?pubchem ?drugbank
WHERE {
  ?compound a dbo:ChemicalCompound .
  ?compound rdfs:label ?compoundLabel .
  FILTER(LANG(?compoundLabel) = "en")

  OPTIONAL { ?compound dbp:casNumber ?casNumber . }
  OPTIONAL { ?compound dbo:pubchem ?pubchem . }
  OPTIONAL { ?compound dbo:drugbank ?drugbank . }
}
LIMIT 500
```

### 3.4 DBpedia Coverage Statistics

| Domain | Entities | Notes |
|--------|----------|-------|
| Total entities | 6.0 million | 2016-04 release |
| Classified entities | 5.2 million | With ontology types |
| Species | 301,000 | Biological species |
| Diseases | ~5,000 | Medical conditions |
| Drugs | ~10,000 | Medications |
| Total triples | 850+ million | As of June 2021 |

**Source:** [DBpedia Wikipedia](https://en.wikipedia.org/wiki/DBpedia)

---

## 4. Coverage Gap Analysis

### 4.1 Critical Gaps in Wikidata

| Gap Area | Current State | Reference Size | Gap % | Priority |
|----------|---------------|----------------|-------|----------|
| Dietary supplements (Q179478) | 0 instances | ~100,000 DSLD | 100% | Critical |
| Dietary minerals (Q15123592) | 0 instances | 12 essential | 100% | Critical |
| Medicinal plants with P2175 | 0 plants | ~28,000 MPNS | 100% | High |
| LOTUS compound links | 51 items | ~800,000 pairs | 99.99% | High |
| Natural Products Atlas | 932 items | ~33,000 | 97% | High |
| DSLD supplement links | 18 items | ~100,000 | 99.98% | High |
| LNHPD links | 88 items | ~90,000 | 99.9% | Medium |
| Gene-GO alignment | 0 links | Variable | 100% | Medium |
| EntitySchema coverage | Limited | All biomedical | 80%+ | Medium |

**Source:** Existing research in `/home/claude/src/gene/research.old/wikidata/wikidata-supplements-herbal.md`

### 4.2 Ontology Integration Gaps

| Ontology | Wikidata Items | Valid Links | Gap |
|----------|----------------|-------------|-----|
| OBO Ontologies represented | 182 | 15 valid | 92% lack URI |
| Disease Ontology | 10,550 | Linked | Good |
| Gene Ontology | Limited | Not aligned | Major gap |
| MeSH integration | Ongoing | Partial | Active work |

**Key Finding:** "Wikidata covers 12,593 genetic association relations. 10,550 of these associations involve items from Disease Ontology. However, none of them are linked to items from Gene Ontology."

**Source:** [Framework for integrating biomedical knowledge in Wikidata](https://pmc.ncbi.nlm.nih.gov/articles/PMC11471508/)

### 4.3 Wikipedia Category Gaps

| Category Area | Current Coverage | Gap Description |
|---------------|------------------|-----------------|
| Rare diseases | ~823 pages | ~10,000+ exist globally |
| Genetic disorders | ~1,200 pages | ~6,000+ known |
| Supplements | ~226 pages | ~100,000+ products |
| Medicinal plants | ~500 pages | ~28,000 species |
| Traditional medicine | ~500 pages | Thousands of formulas |

### 4.4 Missing Biomedical Relation Types

Research identified that 69.3% (579,412) of semantic relations from PMI-generated new associations are missing relation types in Wikidata.

**Source:** [A framework for integrating biomedical knowledge](https://pmc.ncbi.nlm.nih.gov/articles/PMC11471508/)

---

## 5. Infobox Data Extraction

### 5.1 Wikipedia Infobox Templates

| Template | Fields | Wikidata Sync | Extraction Value |
|----------|--------|---------------|------------------|
| [Infobox drug](https://en.wikipedia.org/wiki/Template:Infobox_drug) | 50+ | Yes | High |
| [Infobox gene](https://en.wikipedia.org/wiki/Template:Infobox_gene) | 30+ | Yes (auto) | Critical |
| [Infobox protein](https://en.wikipedia.org/wiki/Template:Infobox_protein) | 25+ | Yes | High |
| [Infobox medical condition](https://en.wikipedia.org/wiki/Template:Infobox_medical_condition) | 25+ | Partial | High |

### 5.2 Infobox Drug Fields

| Field | Type | Description | Wikidata Property |
|-------|------|-------------|-------------------|
| tradename | Text | Brand names | - |
| Drugs.com | Link | Drug reference | P8769 |
| MedlinePlus | ID | NIH reference | P604 |
| routes | Text | Administration routes | P636 |
| CAS_number | ID | CAS Registry | P231 |
| PubChem | ID | PubChem CID | P662 |
| DrugBank | ID | DrugBank ID | P715 |
| ChemSpider | ID | ChemSpider ID | P661 |
| UNII | ID | FDA UNII | P652 |
| KEGG | ID | KEGG Drug | P665 |
| ChEBI | ID | ChEBI ID | P683 |
| ChEMBL | ID | ChEMBL ID | P592 |
| ATC_prefix | Code | ATC class | P267 |
| pregnancy_category | Code | Pregnancy risk | P3489 |
| legal_status | Text | Regulation | P6680 |

**Source:** [Template:Infobox drug](https://en.wikipedia.org/wiki/Template:Infobox_drug)

### 5.3 Infobox Gene Fields (Auto-populated from Wikidata)

| Field | Source | Wikidata Property |
|-------|--------|-------------------|
| Symbol | Wikidata | P353 (HGNC symbol) |
| HGNC | Wikidata | P354 (HGNC ID) |
| Entrez | Wikidata | P351 (Entrez Gene ID) |
| Ensembl | Wikidata | P594 (Ensembl gene ID) |
| OMIM | Wikidata | P492 (OMIM ID) |
| RefSeq | Wikidata | P639/P637 (RefSeq IDs) |
| UniProt | Wikidata | P352 (UniProt ID) |
| Chr | Wikidata | P1057 (Chromosome) |
| Start | Wikidata | P644 (Genomic start) |
| End | Wikidata | P645 (Genomic end) |

**Key Note:** "The data in the Infobox gene template is sourced from Wikidata. The Lua implementation is located at Module:Infobox gene."

**Source:** [Template:Infobox gene](https://en.wikipedia.org/wiki/Template:Infobox_gene)

### 5.4 Infobox Medical Condition Fields

| Field | Description | Potential Wikidata Mapping |
|-------|-------------|---------------------------|
| name | Condition name | rdfs:label |
| synonyms | Alternative names | P2561 (alias) |
| specialty | Medical specialty | P1995 |
| symptoms | Signs and symptoms | P780 |
| complications | Possible complications | - |
| onset | Usual onset age | - |
| duration | Duration | - |
| types | Subtypes | P279 |
| causes | Etiology | P828 |
| risks | Risk factors | P1554 |
| diagnosis | Diagnostic methods | - |
| differential | Differential diagnosis | - |
| prevention | Prevention | - |
| treatment | Treatment options | P2176 |
| medication | Medications used | P2176 |
| prognosis | Expected outcome | - |
| frequency | Prevalence | - |

**Source:** [Template:Infobox medical condition](https://en.wikipedia.org/wiki/Template:Infobox_medical_condition)

### 5.5 Extraction Precision

Research has achieved an average extraction precision of 91% for 1,727 distinct infobox templates using the DBpedia extraction framework.

**Source:** [Semantic Data Extraction from Infobox Wikipedia Template](https://www.researchgate.net/publication/255716376_Semantic_Data_Extraction_from_Infobox_Wikipedia_Template)

---

## 6. Cross-Linking Scientific Databases

### 6.1 Wikidata External Identifier Properties

| Property | Name | Database | Usage Count (Est.) |
|----------|------|----------|-------------------|
| P662 | PubChem CID | PubChem | ~1,329,508 |
| P233 | SMILES | Chemical | ~1,364,830 |
| P234 | InChI | Chemical | ~1,364,612 |
| P235 | InChIKey | Chemical | ~1,365,873 |
| P231 | CAS Number | CAS | ~945,081 |
| P352 | UniProt ID | UniProt | ~27,306 (human) |
| P351 | Entrez Gene ID | NCBI Gene | ~59,721 (human) |
| P353 | HGNC Symbol | HGNC | ~20,000+ |
| P594 | Ensembl Gene ID | Ensembl | ~60,000+ |
| P715 | DrugBank ID | DrugBank | ~10,000+ |
| P592 | ChEMBL ID | ChEMBL | ~50,000+ |
| P683 | ChEBI ID | ChEBI | ~157,189 |
| P486 | MeSH ID | NLM | ~50,000+ |
| P492 | OMIM ID | OMIM | ~25,000+ |
| P3937 | Reactome ID | Reactome | ~3,000+ |
| P2410 | WikiPathways ID | WikiPathways | ~1,500+ |
| P665 | KEGG ID | KEGG | ~30,000+ |
| P1550 | Orphanet ID | Orphanet | ~6,000+ |
| P8121 | NPA ID | Natural Products Atlas | ~932 |
| P11802 | LOTUS ID | LOTUS | ~51 |
| P8743 | DSLD ID | NIH DSLD | ~18 |
| P6955 | LNHPD ID | Health Canada | ~88 |

**Sources:**
- [Wikidata Property:P715](https://www.wikidata.org/wiki/Property:P715)
- [Wikidata Property:P662](https://www.wikidata.org/wiki/Property:P662)
- [Wikidata Property:P352](https://www.wikidata.org/wiki/Property:P352)

### 6.2 Database Interconnection Map

```
                    ┌─────────────────┐
                    │    Wikidata     │
                    │  (Central Hub)  │
                    └────────┬────────┘
                             │
       ┌─────────────────────┼─────────────────────┐
       │                     │                     │
       ▼                     ▼                     ▼
┌─────────────┐       ┌─────────────┐       ┌─────────────┐
│   PubChem   │◄─────►│   ChEMBL    │◄─────►│  DrugBank   │
│  (P662)     │       │  (P592)     │       │  (P715)     │
└──────┬──────┘       └──────┬──────┘       └──────┬──────┘
       │                     │                     │
       └─────────────┬───────┴───────────────┬─────┘
                     │                       │
              ┌──────▼──────┐         ┌──────▼──────┐
              │   UniProt   │◄───────►│    KEGG     │
              │   (P352)    │         │   (P665)    │
              └──────┬──────┘         └──────┬──────┘
                     │                       │
       ┌─────────────┼─────────────┐         │
       │             │             │         │
       ▼             ▼             ▼         ▼
┌─────────────┐ ┌─────────────┐ ┌─────────────┐
│  Reactome   │ │WikiPathways │ │    OMIM     │
│  (P3937)    │ │  (P2410)    │ │   (P492)    │
└─────────────┘ └─────────────┘ └─────────────┘
```

### 6.3 Cross-Reference Coverage Analysis

| Source DB | Target DB | Wikidata Mapping | Direct Links |
|-----------|-----------|------------------|--------------|
| DrugBank | PubChem | Via P662, P715 | Yes |
| DrugBank | UniProt | Via P352, P715 | Yes |
| DrugBank | ChEMBL | Via P592, P715 | Yes |
| ChEMBL | UniProt | Via P352, P8189 | Yes |
| ChEMBL | PubChem | Via P662, P592 | Yes |
| NCBI Gene | Ensembl | Via P351, P594 | Yes |
| NCBI Gene | HGNC | Via P351, P354 | Yes |
| UniProt | PDB | Via P352, P638 | Yes |
| KEGG | Reactome | Via P665, P3937 | Partial |
| OMIM | Orphanet | Via P492, P1550 | Yes |

**Source:** [Mapping Between Databases of Compounds and Protein Targets](https://pmc.ncbi.nlm.nih.gov/articles/PMC7449375/)

### 6.4 Integration Opportunities

| Database | Current State | Opportunity | Estimated Items |
|----------|---------------|-------------|-----------------|
| LOTUS | 51 items | Full integration | ~800,000 pairs |
| Natural Products Atlas | 932 items | Expand coverage | ~33,000 compounds |
| DSLD | 18 items | Product-level data | ~100,000 products |
| LNHPD | 88 items | Canadian supplements | ~90,000 products |
| MPNS | No property | Medicinal plant names | ~28,000 species |
| HMDB | Limited | Metabolome data | ~114,000 metabolites |
| PharmGKB | Limited | PGx relationships | ~40,000 annotations |
| ClinVar | Limited | Clinical variants | ~1,000,000 variants |
| DisGeNET | Limited | Gene-disease | ~600,000 associations |

---

## 7. Recommended Extraction Strategy

### 7.1 Priority Data Sources

| Priority | Source | Data Type | Method |
|----------|--------|-----------|--------|
| 1 | Wikidata SPARQL | Genes, Proteins, Pathways | Direct query |
| 2 | Wikipedia Infoboxes | Drugs, Conditions | DBpedia/custom |
| 3 | Wikidata Dumps | Full entities | Filtered extraction |
| 4 | DBpedia | Drugs, Diseases | SPARQL queries |
| 5 | Federated endpoints | Cross-DB data | SPARQL federation |

### 7.2 Extraction Pipeline

```
┌────────────────────────────────────────────────────────────┐
│                    DATA EXTRACTION PIPELINE                 │
└────────────────────────────────────────────────────────────┘

Step 1: Wikidata Core Entities
├── Query all human genes (Q20747295, ~60K)
├── Query all human proteins (Q8054, ~27K)
├── Query all diseases (Q12136, ~200K)
├── Query all medications (Q12140, ~45K)
└── Query all pathways (Q4915012, ~3K)

Step 2: Wikipedia Infobox Extraction
├── Extract Infobox:drug data via DBpedia
├── Extract Infobox:medical_condition via DBpedia
├── Extract Infobox:gene (supplementary)
└── Extract Infobox:protein (supplementary)

Step 3: Cross-Reference Enrichment
├── Add PubChem links (P662)
├── Add ChEMBL links (P592)
├── Add DrugBank links (P715)
├── Add UniProt links (P352)
└── Add OMIM links (P492)

Step 4: Relationship Extraction
├── Gene-Disease associations (P2293)
├── Drug-Target interactions (P129)
├── Pathway components (P527)
└── Condition treatments (P2175)

Step 5: Gap Filling
├── Supplement data from DSLD/LNHPD
├── Plant compounds from LOTUS
├── Natural products from NPA
└── Traditional medicine formulas
```

### 7.3 SPARQL Query Templates for Bulk Extraction

#### Template A: Complete Gene Export

```sparql
SELECT DISTINCT ?gene ?geneLabel ?entrez ?hgnc ?ensembl ?uniprot
       ?chromosome ?start ?end ?strand
WHERE {
  ?gene wdt:P31/wdt:P279* wd:Q7187 .
  ?gene wdt:P703 wd:Q15978631 .

  OPTIONAL { ?gene wdt:P351 ?entrez . }
  OPTIONAL { ?gene wdt:P353 ?hgnc . }
  OPTIONAL { ?gene wdt:P594 ?ensembl . }
  OPTIONAL { ?gene wdt:P688 ?protein . ?protein wdt:P352 ?uniprot . }
  OPTIONAL { ?gene wdt:P1057 ?chromItem . ?chromItem rdfs:label ?chromosome . FILTER(LANG(?chromosome)="en") }
  OPTIONAL { ?gene wdt:P644 ?start . }
  OPTIONAL { ?gene wdt:P645 ?end . }
  OPTIONAL { ?gene wdt:P2548 ?strandItem . ?strandItem rdfs:label ?strand . FILTER(LANG(?strand)="en") }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

#### Template B: Complete Drug Export

```sparql
SELECT DISTINCT ?drug ?drugLabel ?drugbank ?pubchem ?chembl ?chebi
       ?atc ?smiles ?inchikey ?targetLabel
WHERE {
  ?drug wdt:P31/wdt:P279* wd:Q12140 .

  OPTIONAL { ?drug wdt:P715 ?drugbank . }
  OPTIONAL { ?drug wdt:P662 ?pubchem . }
  OPTIONAL { ?drug wdt:P592 ?chembl . }
  OPTIONAL { ?drug wdt:P683 ?chebi . }
  OPTIONAL { ?drug wdt:P267 ?atc . }
  OPTIONAL { ?drug wdt:P233 ?smiles . }
  OPTIONAL { ?drug wdt:P235 ?inchikey . }
  OPTIONAL { ?drug wdt:P129 ?target . ?target rdfs:label ?targetLabel . FILTER(LANG(?targetLabel)="en") }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

#### Template C: Complete Disease Export

```sparql
SELECT DISTINCT ?disease ?diseaseLabel ?omim ?orphanet ?mesh ?icd10
       ?geneLabel ?treatmentLabel
WHERE {
  ?disease wdt:P31/wdt:P279* wd:Q12136 .

  OPTIONAL { ?disease wdt:P492 ?omim . }
  OPTIONAL { ?disease wdt:P1550 ?orphanet . }
  OPTIONAL { ?disease wdt:P486 ?mesh . }
  OPTIONAL { ?disease wdt:P494 ?icd10 . }
  OPTIONAL { ?disease wdt:P2293 ?gene . ?gene rdfs:label ?geneLabel . FILTER(LANG(?geneLabel)="en") }
  OPTIONAL { ?disease wdt:P2176 ?treatment . ?treatment rdfs:label ?treatmentLabel . FILTER(LANG(?treatmentLabel)="en") }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
```

### 7.4 Recommended Tools

| Tool | Purpose | URL |
|------|---------|-----|
| Wikidata Query Service | Interactive SPARQL | https://query.wikidata.org/ |
| DBpedia SPARQL Editor | DBpedia queries | https://dbpedia.org/sparql |
| wdtk (Wikidata Toolkit) | Dump processing | https://github.com/Wikidata/Wikidata-Toolkit |
| pywikibot | Bot operations | https://www.mediawiki.org/wiki/Manual:Pywikibot |
| SPARQLWrapper | Python SPARQL | https://github.com/RDFLib/sparqlwrapper |

---

## References

### Academic Sources
1. [Wikidata as a knowledge graph for the life sciences](https://pmc.ncbi.nlm.nih.gov/articles/PMC7077981/) - PMC7077981
2. [Wikidata as a semantic framework for the Gene Wiki initiative](https://pmc.ncbi.nlm.nih.gov/articles/PMC4795929/) - PMC4795929
3. [A framework for integrating biomedical knowledge in Wikidata](https://pmc.ncbi.nlm.nih.gov/articles/PMC11471508/) - PMC11471508
4. [Comparing the Chemical Structure and Protein Content of ChEMBL, DrugBank, HMDB and TTD](https://pmc.ncbi.nlm.nih.gov/articles/PMC3916886/) - PMC3916886
5. [Mapping Between Databases of Compounds and Protein Targets](https://pmc.ncbi.nlm.nih.gov/articles/PMC7449375/) - PMC7449375

### Wikidata Resources
- [Wikidata SPARQL Examples](https://www.wikidata.org/wiki/Wikidata:SPARQL_query_service/queries/examples)
- [WikiProject Molecular Biology Properties](https://www.wikidata.org/wiki/Wikidata:WikiProject_Molecular_biology/Properties)
- [Wikidata Federated Queries](https://www.wikidata.org/wiki/Wikidata:SPARQL_query_service/Federated_queries)
- [Sebotic SPARQL Collection](https://github.com/sebotic/SPARQL)

### Wikipedia Resources
- [Portal:Medicine/Categories](https://en.wikipedia.org/wiki/Portal:Medicine/Categories)
- [Template:Infobox drug](https://en.wikipedia.org/wiki/Template:Infobox_drug)
- [Template:Infobox gene](https://en.wikipedia.org/wiki/Template:Infobox_gene)
- [Template:Infobox medical condition](https://en.wikipedia.org/wiki/Template:Infobox_medical_condition)
- [WikiProject Medicine](https://en.wikipedia.org/wiki/Wikipedia:WikiProject_Medicine)
- [WikiProject Molecular Biology](https://en.wikipedia.org/wiki/Wikipedia:WikiProject_Molecular_Biology)

### DBpedia Resources
- [DBpedia SPARQL Endpoint](https://dbpedia.org/sparql)
- [DBpedia Ontology Classes](https://dief.tools.dbpedia.org/server/ontology/classes/)
- [DBpedia Extraction Framework](https://github.com/dbpedia/extraction-framework)
- [DBpedia Mappings Wiki](http://mappings.dbpedia.org/)

### External Databases
- [GARD - Genetic and Rare Diseases](https://rarediseases.info.nih.gov/diseases)
- [NORD Rare Diseases](https://rarediseases.org/rare-diseases/)
- [WikiPathways SPARQL](https://www.wikipathways.org/sparql.html)
- [WikiPathways Federated Queries](https://www.wikipathways.org/federated.html)

---

*Document generated: 2026-01-19*
*For Gene Platform data integration*
