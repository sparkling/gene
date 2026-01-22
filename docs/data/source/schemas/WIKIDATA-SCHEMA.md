# Wikidata Biomedical Schema

**Document ID:** WIKIDATA-SCHEMA
**Status:** Final
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-84-WIKIPEDIA-WIKIDATA](../../integration/WIKIPEDIA-WIKIDATA.md)

---

## Overview

Wikidata is a free and open knowledge base that serves as the central structured data repository for Wikimedia projects. For biomedical applications, it contains comprehensive coverage of human genes (~59,721), proteins (~27,306), diseases (~200,000), and medications (~45,000) under CC0 (public domain) license.

### SPARQL Endpoint

| Attribute | Value |
|-----------|-------|
| **Query Interface** | https://query.wikidata.org/ |
| **SPARQL Endpoint** | https://query.wikidata.org/sparql |
| **Max Results** | 1,000,000 (with streaming) |
| **Timeout** | 60 seconds (default) |
| **Rate Limit** | Not strictly enforced, but be considerate |

---

## Gene Properties

### Primary Gene Identifiers

| Property | Name | Format | Example | Usage Count (Human) |
|----------|------|--------|---------|---------------------|
| P351 | Entrez Gene ID | Numeric string | "7157" | ~59,721 |
| P352 | UniProt Protein ID | [A-Z][0-9][A-Z0-9]{3}[0-9] | "P04637" | ~27,306 |
| P353 | HGNC Gene Symbol | Alphanumeric | "TP53" | ~20,000+ |
| P354 | HGNC ID | Numeric with prefix | "HGNC:11998" | ~20,000+ |
| P594 | Ensembl Gene ID | ENSG[0-9]{11} | "ENSG00000141510" | ~60,000+ |
| P639 | RefSeq Protein ID | [NX]P_[0-9]+\.[0-9]+ | "NP_000537.3" | Variable |
| P637 | RefSeq RNA ID | [NX][MR]_[0-9]+\.[0-9]+ | "NM_000546.6" | Variable |

### Gene Property Details

#### P351 - Entrez Gene ID (NCBI Gene)

```yaml
property: P351
label: Entrez Gene ID
description: Identifier for a gene per the NCBI Gene database
datatype: external-id
formatter_url: https://www.ncbi.nlm.nih.gov/gene/$1
constraints:
  format: "^[1-9][0-9]*$"
  single_value: preferred
  item_requires_statement:
    P31: Q7187 (gene) or subclass
    OR P703: any (found in taxon)
example_query: |
  SELECT ?gene ?geneLabel ?entrezId WHERE {
    ?gene wdt:P351 ?entrezId .
    ?gene wdt:P703 wd:Q15978631 .  # Homo sapiens
    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
  } LIMIT 10
```

#### P352 - UniProt Protein ID

```yaml
property: P352
label: UniProt protein ID
description: Identifier in UniProtKB
datatype: external-id
formatter_url: https://www.uniprot.org/uniprot/$1
constraints:
  format: "^([A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})|([OPQ][0-9][A-Z0-9]{3}[0-9])$"
  item_requires_statement:
    P31: Q8054 (protein) or Q7187 (gene)
inverse_property: P638 (encoded by)
```

#### P353 - HGNC Gene Symbol

```yaml
property: P353
label: HGNC gene symbol
description: Approved symbol for a human gene per HGNC
datatype: external-id
formatter_url: https://www.genenames.org/data/gene-symbol-report/#!/symbol/$1
constraints:
  format: "^[A-Z][A-Z0-9-]*$"
  single_value: recommended
  qualifier_required: P703 wd:Q15978631 (human)
example: "TP53", "BRCA1", "APOE"
```

#### P594 - Ensembl Gene ID

```yaml
property: P594
label: Ensembl gene ID
description: Identifier from Ensembl database
datatype: external-id
formatter_url: https://www.ensembl.org/id/$1
constraints:
  format: "^ENS[A-Z]*G[0-9]{11}$"
  item_requires_statement:
    P31: Q7187 (gene) or subclass
note: |
  Mapping complexity exists: some Ensembl IDs map to multiple NCBI Gene IDs.
  Distribution: 1:2 (486 cases), 1:3 (48), 1:4 (12), 1:5+ (rare)
```

### Gene Entity Classes

| Wikidata Class | QID | Description | Human Count |
|----------------|-----|-------------|-------------|
| gene | Q7187 | Generic gene class | N/A |
| protein-coding gene | Q20747295 | Gene encoding protein | ~19,000 |
| non-coding RNA gene | Q427087 | Gene encoding ncRNA | ~22,000 |
| pseudogene | Q277338 | Non-functional gene copy | ~15,000 |
| human gene | Q54800823 | Human-specific gene | ~59,000 |

### Sample Gene Entity (TP53)

```json
{
  "id": "Q14860489",
  "type": "item",
  "labels": {
    "en": { "value": "TP53", "language": "en" }
  },
  "descriptions": {
    "en": { "value": "human gene", "language": "en" }
  },
  "claims": {
    "P31": [
      { "value": { "id": "Q20747295" } }  // instance of: protein-coding gene
    ],
    "P703": [
      { "value": { "id": "Q15978631" } }  // found in taxon: Homo sapiens
    ],
    "P351": [
      { "value": "7157" }  // Entrez Gene ID
    ],
    "P352": [
      { "value": "P04637" }  // UniProt ID
    ],
    "P353": [
      { "value": "TP53" }  // HGNC symbol
    ],
    "P354": [
      { "value": "11998" }  // HGNC ID
    ],
    "P594": [
      { "value": "ENSG00000141510" }  // Ensembl gene ID
    ],
    "P1057": [
      { "value": { "id": "Q840741" } }  // chromosome: 17
    ],
    "P644": [
      { "value": "7661779" }  // genomic start
    ],
    "P645": [
      { "value": "7687538" }  // genomic end
    ]
  }
}
```

---

## Disease Properties

### Primary Disease Identifiers

| Property | Name | Format | Example | Usage Count |
|----------|------|--------|---------|-------------|
| P492 | OMIM ID | 6-digit number | "191170" | ~25,000+ |
| P1550 | Orphanet ID | Numeric | "558" | ~6,000+ |
| P486 | MeSH ID | D[0-9]{6} | "D009369" | ~50,000+ |
| P699 | Disease Ontology ID | DOID:[0-9]+ | "DOID:162" | ~10,000+ |
| P5270 | MONDO ID | MONDO:[0-9]{7} | "MONDO:0005015" | ~26,000 |

### Disease Property Details

#### P492 - OMIM ID

```yaml
property: P492
label: OMIM ID
description: Online Mendelian Inheritance in Man identifier
datatype: external-id
formatter_url: https://www.omim.org/entry/$1
constraints:
  format: "^[1-9][0-9]{5}$"
  single_value: preferred
qualifier_properties:
  - P1480 (sourcing circumstances)
  - P518 (applies to part)
```

#### P1550 - Orphanet ID

```yaml
property: P1550
label: Orphanet ID
description: Identifier in Orphanet rare disease database
datatype: external-id
formatter_url: https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=EN&Expert=$1
constraints:
  format: "^[1-9][0-9]*$"
  item_requires_statement:
    P31: Q112193867 (rare disease) or subclass
```

### Disease Entity Classes

| Wikidata Class | QID | Description | Count |
|----------------|-----|-------------|-------|
| disease | Q12136 | Human disease | ~200,000 |
| genetic disorder | Q18553442 | Inherited condition | ~6,000 |
| rare disease | Q929833 | Orphan disease | ~7,000 |
| syndrome | Q179630 | Syndrome | ~5,000 |
| infectious disease | Q18123741 | Pathogen-caused | ~2,000 |

### Gene-Disease Association Properties

| Property | Name | Description | Example |
|----------|------|-------------|---------|
| P2293 | genetic association | Gene linked to disease | disease → gene |
| P780 | symptoms and signs | Clinical manifestations | disease → symptom |
| P2176 | drug or therapy used for treatment | Treatment | disease → drug |

---

## Chemical/Drug Properties

### Primary Chemical Identifiers

| Property | Name | Format | Example | Usage Count |
|----------|------|--------|---------|-------------|
| P662 | PubChem CID | Numeric | "5284373" | ~1,329,508 |
| P683 | ChEBI ID | Numeric | "15377" | ~157,189 |
| P231 | CAS Registry Number | [0-9]+-[0-9]{2}-[0-9] | "50-00-0" | ~945,081 |
| P715 | DrugBank ID | DB[0-9]{5} | "DB00945" | ~10,000+ |
| P592 | ChEMBL ID | CHEMBL[0-9]+ | "CHEMBL25" | ~50,000+ |
| P235 | InChIKey | [A-Z]{14}-[A-Z]{10}-[A-Z] | "BSYNRYMUTXBXSQ-..." | ~1,365,873 |

### Chemical Structure Properties

| Property | Name | Description |
|----------|------|-------------|
| P233 | SMILES | Simplified molecular notation |
| P234 | InChI | International Chemical Identifier |
| P235 | InChIKey | Hashed InChI |
| P274 | chemical formula | Molecular formula |
| P2067 | mass | Molecular weight |

### Drug-Target Properties

| Property | Name | Description |
|----------|------|-------------|
| P129 | physically interacts with | Drug-protein interaction |
| P2175 | medical condition treated | Therapeutic indication |
| P3780 | active ingredient in | Formulation component |
| P636 | route of administration | Delivery method |

---

## SPARQL Query Patterns

### Count Human Genes

```sparql
SELECT (COUNT(DISTINCT ?gene) AS ?count) WHERE {
  ?gene wdt:P351 ?entrezId .
  ?gene wdt:P703 wd:Q15978631 .  # Homo sapiens
}
# Result: ~59,721 genes
```

### Count Diseases with OMIM

```sparql
SELECT (COUNT(DISTINCT ?disease) AS ?count) WHERE {
  ?disease wdt:P492 ?omimId .
}
# Result: ~25,000+ diseases
```

### Gene with All Identifiers

```sparql
SELECT ?gene ?geneLabel ?entrez ?uniprot ?hgnc ?ensembl ?refseq WHERE {
  ?gene wdt:P31/wdt:P279* wd:Q7187 .  # gene or subclass
  ?gene wdt:P703 wd:Q15978631 .        # human

  OPTIONAL { ?gene wdt:P351 ?entrez . }
  OPTIONAL { ?gene wdt:P352 ?uniprot . }
  OPTIONAL { ?gene wdt:P353 ?hgnc . }
  OPTIONAL { ?gene wdt:P594 ?ensembl . }
  OPTIONAL { ?gene wdt:P639 ?refseq . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
ORDER BY ?hgnc
LIMIT 100
```

### Disease-Gene Associations

```sparql
SELECT ?disease ?diseaseLabel ?gene ?geneLabel ?omim ?orphanet WHERE {
  ?disease wdt:P2293 ?gene .           # genetic association
  ?gene wdt:P703 wd:Q15978631 .        # human gene

  OPTIONAL { ?disease wdt:P492 ?omim . }
  OPTIONAL { ?disease wdt:P1550 ?orphanet . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
ORDER BY ?diseaseLabel
LIMIT 500
```

### Drugs and Targets

```sparql
SELECT ?drug ?drugLabel ?target ?targetLabel ?pubchem ?drugbank WHERE {
  ?drug wdt:P31 wd:Q12140 .            # medication
  ?drug wdt:P129 ?target .              # interacts with
  ?target wdt:P31/wdt:P279* wd:Q8054 . # protein
  ?target wdt:P703 wd:Q15978631 .       # human

  OPTIONAL { ?drug wdt:P662 ?pubchem . }
  OPTIONAL { ?drug wdt:P715 ?drugbank . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 500
```

---

## Entity Statistics (Estimated)

| Entity Type | Wikidata QID | Count | Coverage |
|-------------|--------------|-------|----------|
| Human genes | wdt:P703 Q15978631 | ~59,721 | ~100% of NCBI |
| Mouse genes | wdt:P703 Q83310 | ~73,355 | ~100% of NCBI |
| Human proteins | Q8054 + Q15978631 | ~27,306 | ~100% SwissProt |
| Diseases | Q12136 | ~200,000 | ~90% DO |
| Medications | Q12140 | ~45,000 | ~50% DrugBank |
| Rare diseases | Q929833 | ~7,000 | Variable |
| Genetic disorders | Q18553442 | ~6,000 | Variable |
| Chemical compounds | Q11173 | ~1,300,000+ | ~1% PubChem |

---

## Cross-Reference Property Summary

### Genomic Databases
- P351: NCBI Gene (Entrez)
- P594: Ensembl
- P354: HGNC
- P639/P637: RefSeq

### Protein Databases
- P352: UniProt
- P638: Encodes (inverse)
- P705: UniRef50

### Disease/Phenotype
- P492: OMIM
- P1550: Orphanet
- P5270: MONDO
- P699: Disease Ontology
- P3841: HPO

### Chemical/Drug
- P662: PubChem
- P683: ChEBI
- P231: CAS
- P715: DrugBank
- P592: ChEMBL

### Pathway/Function
- P3937: Reactome
- P2410: WikiPathways
- P686: Gene Ontology

---

## Data Quality Notes

1. **Gene Coverage**: Human gene coverage is nearly complete thanks to Gene Wiki initiative importing from NCBI Gene

2. **Mapping Complexity**: Some 1:N mappings exist between Ensembl and NCBI Gene IDs (486 cases of 1:2, 48 of 1:3)

3. **Missing Links**: Gene Ontology alignment is minimal (no systematic GO links)

4. **Supplement Gap**: Dietary supplements have nearly 0% coverage in Wikidata

5. **Update Frequency**: Biomedical data is maintained by bots and community contributors; updates are continuous

---

## References

- Gene Wiki: https://en.wikipedia.org/wiki/Gene_Wiki
- Wikidata:WikiProject Molecular Biology
- Malone et al. (2016) "Wikidata as a semantic framework for the Gene Wiki initiative" Database
