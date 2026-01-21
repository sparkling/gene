# Wikipedia & Wikidata Data Sources

**Document ID:** 43-84-WIKIPEDIA-WIKIDATA
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-DATA-SOURCES.md](../43-DATA-SOURCES.md)

---

## TL;DR

Wikidata provides **59,721 human genes with 100% NCBI Gene coverage** and **200K+ diseases** under CC0 license, making it an ideal free knowledge graph backbone. Use **Wikidata SPARQL** as primary extraction method for genes, proteins, diseases, and cross-database identifiers. **DBpedia** supplements with ~10K drugs and ~5K diseases extracted from Wikipedia infoboxes. Critical coverage gaps exist for dietary supplements (100% missing), medicinal plants (100% missing), and Gene Ontology alignment (0 links). Wikipedia categories provide curated genetics content (~3K pages) but require DBpedia/custom extraction for structured data.

---

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **Primary Knowledge Graph** | Wikidata SPARQL | CC0 license, comprehensive gene/protein coverage |
| **Supplementary Extraction** | DBpedia | Wikipedia infobox structured data |
| **Gene Data Source** | Wikidata + NCBI | ~60K human genes, 100% coverage |
| **Disease Taxonomy** | Wikidata Disease Ontology | 10,550 DO-linked diseases |
| **Drug Metadata** | DBpedia + Wikidata | ~10K drugs with identifiers |
| **Federated Queries** | WikiPathways, UniProt | Cross-database enrichment |
| **Extraction Method** | SPARQL queries | Direct, no web scraping needed |

---

## Database Catalog

### 1. Wikidata Knowledge Graph

#### 1.1 Overview

| Attribute | Value |
|-----------|-------|
| **Provider** | Wikimedia Foundation |
| **URL** | https://www.wikidata.org |
| **SPARQL Endpoint** | https://query.wikidata.org/sparql |
| **GUI Interface** | https://query.wikidata.org/ |
| **Total Items** | 100+ million |
| **License** | CC0 (Public Domain) |
| **Commercial Use** | YES - fully open |

#### 1.2 Biomedical Entity Coverage

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

#### 1.3 External Identifier Properties

| Property | Name | Database | Usage Count |
|----------|------|----------|-------------|
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

#### 1.4 Supplement/Natural Product Identifiers (Sparse)

| Property | Name | Database | Usage Count | Gap |
|----------|------|----------|-------------|-----|
| P8121 | NPA ID | Natural Products Atlas | ~932 | 97% |
| P11802 | LOTUS ID | LOTUS | ~51 | 99.99% |
| P8743 | DSLD ID | NIH DSLD | ~18 | 99.98% |
| P6955 | LNHPD ID | Health Canada | ~88 | 99.9% |

---

### 2. Wikipedia Categories

#### 2.1 Genetics Categories

**Primary Category:** [Category:Genetics](https://en.wikipedia.org/wiki/Category:Genetics)

| Subcategory | Est. Pages | Relevance | Description |
|-------------|------------|-----------|-------------|
| Genetic disorders | 1,000+ | Critical | All genetic conditions |
| Genes | 500+ | Critical | Gene articles by organism |
| Human genetics | 400+ | Critical | Human-specific genetics |
| Medical genetics | 600+ | Critical | Clinical genetics |
| Pharmacogenomics | 50+ | Critical | Drug-gene interactions |
| Genomics | 300+ | High | Genome-level analysis |
| DNA | 200+ | High | DNA structure, replication |
| Chromosomes | 150+ | High | Chromosome biology |
| Molecular genetics | 250+ | High | Molecular mechanisms |
| Epigenetics | 100+ | High | Epigenetic modifications |
| Gene expression | 150+ | High | Transcription, regulation |
| Population genetics | 100+ | Medium | Population-level variation |
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

#### 2.2 Health Conditions Categories

**Primary Category:** [Category:Diseases and disorders](https://en.wikipedia.org/wiki/Category:Diseases_and_disorders)

| Category | Subcategories | Pages | Description |
|----------|---------------|-------|-------------|
| Rare diseases | 19 | ~823 | Orphan diseases |
| Genetic disorders | 25+ | ~1,200 | Inherited conditions |
| Rare syndromes | 3 | ~380 | Syndromic conditions |
| Rare genetic syndromes | 1 | ~185 | Genetic syndromes |
| Autoimmune diseases | 10+ | ~300 | Immune disorders |
| Metabolic disorders | 15+ | ~500 | Metabolic conditions |
| Neurological disorders | 20+ | ~800 | Nervous system |
| Cardiovascular diseases | 15+ | ~600 | Heart/vascular |
| Cancer | 30+ | ~1,500 | Neoplasms |
| Infectious diseases | 25+ | ~1,000 | Pathogen-caused |

#### 2.3 Supplements Categories

**Primary Category:** [Category:Dietary supplements](https://en.wikipedia.org/wiki/Category:Dietary_supplements)

| Subcategory | Pages | Description |
|-------------|-------|-------------|
| Vitamins | ~26 | Vitamin articles |
| B vitamins | ~21 | B-complex vitamins |
| Minerals (nutrition) | ~30 | Dietary minerals |
| Amino acid supplements | ~20 | Amino acids |
| Bodybuilding supplements | ~40 | Sports nutrition |
| Effervescent supplements | ~10 | Effervescent forms |

**Vitamin Hierarchy:**
```
Category:Vitamins (10 subcategories, 26 pages)
├── Category:B vitamins (21 pages)
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

#### 2.4 Alternative Medicine Categories

**Primary Category:** [Category:Alternative medicine](https://en.wikipedia.org/wiki/Category:Alternative_medicine)

| Subcategory | Pages | Description |
|-------------|-------|-------------|
| Herbal medicine | ~150 | Plant-based medicine |
| Traditional Chinese medicine | ~200 | TCM practices |
| Ayurveda | ~100 | Indian medicine |
| Homeopathy | ~75 | Homeopathic remedies |
| Naturopathy | ~50 | Natural therapies |
| Acupuncture | ~40 | Acupuncture points |

---

### 3. DBpedia

#### 3.1 Overview

| Attribute | Value |
|-----------|-------|
| **Provider** | DBpedia Association |
| **URL** | https://dbpedia.org |
| **SPARQL Endpoint** | http://dbpedia.org/sparql |
| **Total Entities** | 6.0 million (classified: 5.2M) |
| **Total Triples** | 850+ million |
| **License** | CC BY-SA / GNU FDL |
| **Commercial Use** | YES (with attribution) |

#### 3.2 Biomedical Ontology Classes

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
| Species | N/A | 301,000 | Biological species |

---

### 4. Federated SPARQL Endpoints

| Endpoint | URL | Description |
|----------|-----|-------------|
| Wikidata | https://query.wikidata.org/sparql | Main knowledge graph |
| WikiPathways | https://sparql.wikipathways.org/sparql | Pathway data |
| UniProt | https://sparql.uniprot.org/sparql | Protein data |
| ChEMBL Mirror | https://chemblmirror.rdf.bigcat-bioinformatics.org/sparql | Drug bioactivity |
| Rhea | https://sparql.rhea-db.org/sparql | Biochemical reactions |
| neXtProt | https://sparql.nextprot.org/sparql | Human protein data |

---

## Coverage Gap Analysis

### Critical Gaps in Wikidata

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

### Ontology Integration Gaps

| Ontology | Wikidata Items | Valid Links | Gap |
|----------|----------------|-------------|-----|
| OBO Ontologies represented | 182 | 15 valid | 92% lack URI |
| Disease Ontology | 10,550 | Linked | Good |
| Gene Ontology | Limited | Not aligned | Major gap |
| MeSH integration | Ongoing | Partial | Active work |

### Wikipedia Category Gaps

| Category Area | Current Coverage | Global Estimate |
|---------------|------------------|-----------------|
| Rare diseases | ~823 pages | ~10,000+ exist |
| Genetic disorders | ~1,200 pages | ~6,000+ known |
| Supplements | ~226 pages | ~100,000+ products |
| Medicinal plants | ~500 pages | ~28,000 species |
| Traditional medicine | ~500 pages | Thousands of formulas |

---

## Essential SPARQL Queries

### Query 1: All Human Genes with Identifiers

```sparql
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

### Query 2: Gene-Disease Associations

```sparql
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

### Query 3: Drugs and Their Protein Targets

```sparql
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

### Query 4: Vitamins with Biological Data

```sparql
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

### Query 5: CYP450 Enzymes (Pharmacogenomics)

```sparql
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

### Query 6: Rare Genetic Diseases

```sparql
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

---

## Wikipedia Infobox Extraction

### Key Infobox Templates

| Template | Fields | Wikidata Sync | Extraction Value |
|----------|--------|---------------|------------------|
| Infobox drug | 50+ | Yes | High |
| Infobox gene | 30+ | Yes (auto) | Critical |
| Infobox protein | 25+ | Yes | High |
| Infobox medical condition | 25+ | Partial | High |

### Infobox Gene Fields (Auto from Wikidata)

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

### Infobox Drug Fields

| Field | Type | Wikidata Property |
|-------|------|-------------------|
| tradename | Text | - |
| Drugs.com | Link | P8769 |
| MedlinePlus | ID | P604 |
| routes | Text | P636 |
| CAS_number | ID | P231 |
| PubChem | ID | P662 |
| DrugBank | ID | P715 |
| ChemSpider | ID | P661 |
| UNII | ID | P652 |
| KEGG | ID | P665 |
| ChEBI | ID | P683 |
| ChEMBL | ID | P592 |
| ATC_prefix | Code | P267 |

### Extraction Precision

Research achieves **91% average extraction precision** for 1,727 distinct infobox templates using DBpedia extraction framework.

---

## DBpedia SPARQL Queries

### Query 1: All Drugs with Wikipedia Links

```sparql
PREFIX dbo: <http://dbpedia.org/ontology/>
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

### Query 2: Diseases with Genes and Symptoms

```sparql
PREFIX dbo: <http://dbpedia.org/ontology/>

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

### Query 3: Chemical Compounds with Identifiers

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

---

## Integration Priority

| Source | Priority | License | Method | Key Use Case |
|--------|----------|---------|--------|--------------|
| Wikidata Genes | HIGH | CC0 | SPARQL | Gene identifiers, cross-refs |
| Wikidata Diseases | HIGH | CC0 | SPARQL | Disease taxonomy |
| Wikidata Proteins | HIGH | CC0 | SPARQL | Protein data |
| Wikidata Drugs | MEDIUM | CC0 | SPARQL | Drug metadata |
| DBpedia Drugs | MEDIUM | CC BY-SA | SPARQL | Wikipedia drug info |
| DBpedia Diseases | MEDIUM | CC BY-SA | SPARQL | Disease descriptions |
| WikiPathways (fed) | MEDIUM | CC0 | Federated SPARQL | Pathway data |
| Wikipedia Categories | LOW | CC BY-SA | API/Scraping | Curated lists |

---

## Recommended Extraction Strategy

### Phase 1: Wikidata Core Entities

| Entity | Query | Est. Records | Priority |
|--------|-------|--------------|----------|
| Human genes (Q20747295) | Gene query | ~60K | Week 1 |
| Human proteins (Q8054) | Protein query | ~27K | Week 1 |
| Diseases (Q12136) | Disease query | ~200K | Week 2 |
| Medications (Q12140) | Drug query | ~45K | Week 2 |
| Pathways (Q4915012) | Pathway query | ~3K | Week 3 |

### Phase 2: Cross-Reference Enrichment

| Identifier | Property | Records |
|------------|----------|---------|
| PubChem | P662 | ~1.3M |
| ChEMBL | P592 | ~50K |
| DrugBank | P715 | ~10K |
| UniProt | P352 | ~27K |
| OMIM | P492 | ~25K |

### Phase 3: DBpedia Supplementation

| Entity | Class | Est. Records |
|--------|-------|--------------|
| Drugs | dbo:Drug | ~10K |
| Diseases | dbo:Disease | ~5K |
| Compounds | dbo:ChemicalCompound | ~50K |

### Phase 4: Gap Filling (Future)

| Gap | Source | Est. Records |
|-----|--------|--------------|
| Supplements | DSLD/LNHPD import | ~100K |
| Plant compounds | LOTUS import | ~800K |
| Natural products | NPA import | ~33K |
| Traditional medicine | Manual curation | ~5K |

---

## Recommended Tools

| Tool | Purpose | URL |
|------|---------|-----|
| Wikidata Query Service | Interactive SPARQL | https://query.wikidata.org/ |
| DBpedia SPARQL Editor | DBpedia queries | https://dbpedia.org/sparql |
| wdtk (Wikidata Toolkit) | Dump processing | https://github.com/Wikidata/Wikidata-Toolkit |
| pywikibot | Bot operations | https://www.mediawiki.org/wiki/Manual:Pywikibot |
| SPARQLWrapper | Python SPARQL | https://github.com/RDFLib/sparqlwrapper |

---

## Dependencies

### Upstream Dependencies

| Dependency | Purpose | Risk if Unavailable |
|------------|---------|---------------------|
| Wikidata SPARQL | Primary data source | HIGH - no equivalent |
| DBpedia SPARQL | Supplementary data | LOW - Wikidata backup |
| WikiPathways SPARQL | Pathway federation | LOW - Reactome alternative |

### Downstream Dependents

| Dependent | Usage |
|-----------|-------|
| Gene Identifier Mapping | Cross-database ID resolution |
| Disease Taxonomy | Disease classification and hierarchy |
| Drug Metadata | Drug information display |
| Knowledge Graph | Entity relationships |
| Search Autocomplete | Entity labels and aliases |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Initial document from research.old/data-sources-wikipedia-wikidata.md |
