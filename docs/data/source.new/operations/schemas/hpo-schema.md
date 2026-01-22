---
id: schema-hpo
title: Human Phenotype Ontology (HPO) Schema
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, ontology, phenotype, rare-disease]
---

# Human Phenotype Ontology (HPO) Schema

**Document ID:** HPO-SCHEMA
**Status:** Final
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [Schema Documentation](./_index.md)

---

## Overview

The Human Phenotype Ontology (HPO) provides a standardized vocabulary of phenotypic abnormalities and clinical features encountered in human disease. It enables precise clinical phenotyping for rare disease diagnosis, genomic interpretation, and translational research.

### Key Resources

| Resource | URL |
|----------|-----|
| **Homepage** | https://hpo.jax.org/ |
| **OBO Foundry** | https://obofoundry.org/ontology/hp.html |
| **GitHub** | https://github.com/obophenotype/human-phenotype-ontology |
| **Latest Release** | https://github.com/obophenotype/human-phenotype-ontology/releases/latest |
| **Documentation** | https://hpo.jax.org/app/help/introduction |

---

## Statistics (v2026-01-08)

| Metric | Count |
|--------|-------|
| **Total terms** | 13,000+ |
| **Disease annotations** | 156,000+ |
| **New terms (recent)** | 2,239 |
| **New annotations (recent)** | 49,235 |
| **PubMed citations** | 9,573 PMIDs |
| **GitHub commits** | 8,308 |
| **Releases** | 61 |
| **Contributors** | 24+ |
| **Languages** | 10+ |

### Term Distribution by Category

| Branch | Description | Est. Terms |
|--------|-------------|------------|
| Phenotypic abnormality | Clinical features | ~11,000 |
| Mode of inheritance | Inheritance patterns | ~20 |
| Clinical modifier | Severity/modifiers | ~200 |
| Clinical course | Temporal patterns | ~100 |
| Past medical history | History terms | ~50 |

---

## Identifier Format

### Standard Format
```
HP:0000001
```

### Pattern
```
HP:[0-9]{7}
```

### URI Format
```
http://purl.obolibrary.org/obo/HP_0000001
```

**Examples:**
- HP:0000118 = Phenotypic abnormality (root)
- HP:0001250 = Seizure
- HP:0001627 = Abnormal heart morphology
- HP:0002719 = Recurrent infections

---

## File Formats

### Available Downloads (v2026-01-08)

| File | Size | Description |
|------|------|-------------|
| hp.obo | 10.4 MB | Standard OBO format |
| hp.owl | 75.6 MB | OWL format |
| hp.json | 21.9 MB | JSON format |
| hp-full.obo | 18.8 MB | With imports |
| hp-full.owl | 84.7 MB | Full OWL with imports |
| hp-base.obo | 10.9 MB | Base ontology |
| hp-international.obo | 21.8 MB | Multi-language |
| hp-international.owl | 207.1 MB | Multi-language OWL |
| phenotype.hpoa | 35.0 MB | Disease annotations |
| genes_to_phenotype.txt | 20.1 MB | Gene-phenotype links |
| phenotype_to_genes.txt | 64.6 MB | Phenotype-gene links |

### Download URLs

```bash
# Stable PURL (redirects to latest)
http://purl.obolibrary.org/obo/hp.obo
http://purl.obolibrary.org/obo/hp.owl
http://purl.obolibrary.org/obo/hp.json

# Direct GitHub release
https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-01-08/hp.obo
```

---

## OBO Format Structure

### Header Section

```obo
format-version: 1.2
data-version: hp/releases/2026-01-08
date: 08:01:2026 12:00
saved-by: Peter Robinson
auto-generated-by: OWL API
default-namespace: human_phenotype
ontology: hp
property_value: http://purl.org/dc/elements/1.1/creator "Human Phenotype Ontology Consortium"
property_value: http://purl.org/dc/elements/1.1/description "A standardized vocabulary of phenotypic abnormalities and clinical features encountered in human disease."
property_value: http://purl.org/dc/elements/1.1/rights "https://hpo.jax.org/app/license"
property_value: http://purl.org/dc/elements/1.1/title "Human Phenotype Ontology"
```

### Subset Definitions

```obo
subsetdef: hpo_slim "Core clinical terminology"
subsetdef: hposlim_core "HPO slim core"
subsetdef: secondary_consequence "Secondary consequences"
```

### Synonym Types

```obo
synonymtypedef: abbreviation "abbreviation"
synonymtypedef: layperson "layperson term"
synonymtypedef: obsolete_synonym "obsolete_synonym"
synonymtypedef: plural_form "plural form"
synonymtypedef: uk_spelling "UK spelling"
```

### Term Structure

```obo
[Term]
id: HP:0001250
name: Seizure
def: "A seizure is an event due to abnormal synchronous neuronal activity in the brain." [HPO:probinson, PMID:22983531]
comment: Seizures are sudden, transient symptoms resulting from excessive or synchronous neuronal activity in the brain.
synonym: "Epileptic seizure" EXACT []
synonym: "Seizures" EXACT plural_form []
synonym: "Convulsion" RELATED layperson []
synonym: "Fit" RELATED layperson []
xref: MSH:D012640
xref: SNOMEDCT_US:91175000
xref: UMLS:C0036572
is_a: HP:0012638 ! Abnormal nervous system physiology
created_by: peter
creation_date: 2008-02-27T02:20:00Z
```

### Sample Phenotype Entry (Cardiomegaly)

```obo
[Term]
id: HP:0001640
name: Cardiomegaly
alt_id: HP:0003911
def: "Increased size of the heart, clinically defined as >2 standard deviations above the mean for age, sex, and body surface area." [PMID:19793835]
subset: hposlim_core
synonym: "Cardiac hypertrophy" NARROW []
synonym: "Enlarged heart" EXACT layperson []
synonym: "Increased heart size" EXACT []
xref: MSH:D006332
xref: SNOMEDCT_US:8186001
xref: UMLS:C0018800
is_a: HP:0001627 ! Abnormal heart morphology
property_value: http://purl.org/dc/terms/creator "HPO:probinson" xsd:string
created_by: peter
creation_date: 2008-02-27T02:20:00Z
```

---

## Ontology Hierarchy

### Top-Level Structure

```
HP:0000001 All
├── HP:0000005 Mode of inheritance
│   ├── HP:0000006 Autosomal dominant inheritance
│   ├── HP:0000007 Autosomal recessive inheritance
│   ├── HP:0001417 X-linked inheritance
│   ├── HP:0001419 X-linked recessive inheritance
│   ├── HP:0001423 X-linked dominant inheritance
│   ├── HP:0001427 Mitochondrial inheritance
│   └── HP:0010982 Polygenic inheritance
├── HP:0000118 Phenotypic abnormality
│   ├── HP:0000152 Abnormality of head or neck
│   ├── HP:0000478 Abnormality of the eye
│   ├── HP:0000598 Abnormality of the ear
│   ├── HP:0000707 Abnormality of the nervous system
│   ├── HP:0000769 Abnormality of the breast
│   ├── HP:0000818 Abnormality of the endocrine system
│   ├── HP:0000924 Abnormality of the skeletal system
│   ├── HP:0001197 Abnormality of prenatal development
│   ├── HP:0001507 Growth abnormality
│   ├── HP:0001574 Abnormality of the integument
│   ├── HP:0001608 Abnormality of the voice
│   ├── HP:0001626 Abnormality of the cardiovascular system
│   ├── HP:0001871 Abnormality of blood and blood-forming tissues
│   ├── HP:0001939 Abnormality of metabolism/homeostasis
│   ├── HP:0002011 Abnormality of the musculature
│   ├── HP:0002012 Abnormality of the abdominal organs
│   ├── HP:0002086 Abnormality of the respiratory system
│   ├── HP:0002715 Abnormality of the immune system
│   ├── HP:0003011 Abnormality of the musculature of the limbs
│   ├── HP:0025031 Abnormality of the digestive system
│   ├── HP:0025142 Constitutional symptom
│   └── HP:0040064 Abnormality of limbs
├── HP:0012823 Clinical modifier
│   ├── HP:0012824 Severity
│   ├── HP:0012830 Position
│   ├── HP:0025280 Pain characteristic
│   └── HP:0031797 Clinical course
├── HP:0032223 Blood group
├── HP:0032443 Past medical history
└── HP:0040279 Frequency
```

### Phenotypic Abnormality Subcategories

| HP ID | Name | Description |
|-------|------|-------------|
| HP:0000707 | Abnormality of the nervous system | Neurological |
| HP:0001626 | Abnormality of the cardiovascular system | Cardiac |
| HP:0001939 | Abnormality of metabolism/homeostasis | Metabolic |
| HP:0002086 | Abnormality of the respiratory system | Pulmonary |
| HP:0000924 | Abnormality of the skeletal system | Skeletal |
| HP:0001574 | Abnormality of the integument | Skin/hair |
| HP:0000478 | Abnormality of the eye | Ocular |
| HP:0000598 | Abnormality of the ear | Hearing |
| HP:0002715 | Abnormality of the immune system | Immunological |
| HP:0001871 | Abnormality of blood | Hematological |

---

## Phenotype Annotation Format (phenotype.hpoa)

### File Structure

The phenotype.hpoa file is a tab-separated file containing disease-phenotype associations.

### Column Specification (12 columns)

| # | Column | Required | Description | Example |
|---|--------|----------|-------------|---------|
| 1 | database_id | Yes | Disease database ID | OMIM:154700 |
| 2 | disease_name | Yes | Disease name | Marfan syndrome |
| 3 | qualifier | No | NOT if negated | NOT (or empty) |
| 4 | hpo_id | Yes | HPO term ID | HP:0001166 |
| 5 | reference | Yes | Evidence source | OMIM:154700 |
| 6 | evidence | Yes | Evidence code | IEA, TAS, or PCS |
| 7 | onset | No | Age of onset | HP:0003577 |
| 8 | frequency | No | Frequency info | HP:0040280 or 5/13 |
| 9 | sex | No | Sex specificity | MALE or FEMALE |
| 10 | modifier | No | Clinical modifiers | HP:0012828 (semicolon-sep) |
| 11 | aspect | Yes | HPO branch | P, I, C, M, or H |
| 12 | biocuration | Yes | Curator/date | HPO:skoehler[2024-01-15] |

### Evidence Codes

| Code | Name | Description |
|------|------|-------------|
| IEA | Inferred from Electronic Annotation | Automated extraction (OMIM text) |
| TAS | Traceable Author Statement | Expert assertion with citation |
| PCS | Published Clinical Study | From clinical research |

### Aspect Codes

| Code | Branch | Description |
|------|--------|-------------|
| P | Phenotypic abnormality | Clinical features |
| I | Inheritance | Mode of inheritance |
| C | Clinical course | Onset, progression |
| M | Clinical modifier | Severity, laterality |
| H | Past medical history | Historical features |

### Sample phenotype.hpoa Entries

```tsv
#description: HPO annotations for rare diseases. See https://hpo.jax.org/
#date: 2026-01-08
#tracker: https://github.com/obophenotype/human-phenotype-ontology/issues
#HPO-version: hp/releases/2026-01-08/hp.obo
database_id	disease_name	qualifier	hpo_id	reference	evidence	onset	frequency	sex	modifier	aspect	biocuration
OMIM:154700	Marfan syndrome		HP:0001166	OMIM:154700	TAS	HP:0003577	HP:0040283		HP:0012828	P	HPO:skoehler[2024-01-15]
OMIM:154700	Marfan syndrome		HP:0001519	OMIM:154700	TAS		HP:0040281			P	HPO:skoehler[2024-01-15]
OMIM:154700	Marfan syndrome		HP:0000483	PMID:20301565	PCS	HP:0011463	12/45			P	HPO:probinson[2024-03-20]
OMIM:154700	Marfan syndrome		HP:0000006	OMIM:154700	TAS					I	HPO:skoehler[2024-01-15]
OMIM:154700	Marfan syndrome	NOT	HP:0001249	PMID:18451852	PCS					P	HPO:ivasilevsky[2024-02-10]
```

### Database ID Prefixes

| Prefix | Source | Example |
|--------|--------|---------|
| OMIM: | OMIM | OMIM:154700 |
| ORPHA: | Orphanet | ORPHA:558 |
| DECIPHER: | DECIPHER | DECIPHER:16 |
| MONDO: | MONDO | MONDO:0007947 |

---

## Gene Annotation Files

### genes_to_phenotype.txt

Maps genes to their associated phenotypes.

```tsv
#Format: gene_id<tab>gene_symbol<tab>hpo_id<tab>hpo_name<tab>frequency<tab>disease_id
7157	TP53	HP:0002664	Neoplasm	-	OMIM:191170
7157	TP53	HP:0001909	Leukemia	HP:0040283	OMIM:191170
7157	TP53	HP:0100526	Neoplasm of the lung	-	OMIM:614327
```

### phenotype_to_genes.txt

Maps phenotypes to associated genes.

```tsv
#Format: hpo_id<tab>hpo_name<tab>ncbi_gene_id<tab>gene_symbol<tab>disease_id
HP:0001250	Seizure	2890	GRIA1	OMIM:620003
HP:0001250	Seizure	6323	SCN1A	OMIM:607208
HP:0001250	Seizure	6326	SCN2A	OMIM:613721
HP:0001250	Seizure	6329	SCN8A	OMIM:614558
```

---

## Frequency Annotations

### HPO Frequency Terms

| HP ID | Name | Meaning |
|-------|------|---------|
| HP:0040280 | Obligate | 100% |
| HP:0040281 | Very frequent | 80-99% |
| HP:0040282 | Frequent | 30-79% |
| HP:0040283 | Occasional | 5-29% |
| HP:0040284 | Very rare | 1-4% |
| HP:0040285 | Excluded | 0% |

### Numeric Frequency

Can also be expressed as:
- Fraction: `12/45` (12 of 45 patients)
- Percentage: `27%`

---

## Cross-References

### External Mappings

| Source | Format | Example |
|--------|--------|---------|
| UMLS | CUI | UMLS:C0036572 |
| SNOMED CT | SCTID | SNOMEDCT_US:91175000 |
| MeSH | D-number | MSH:D012640 |
| MedDRA | PT code | MedDRA:10039906 |
| ICD-10 | Code | ICD10CM:G40 |

### Ontology Alignments

| Ontology | Relationship |
|----------|--------------|
| MONDO | Disease associations |
| Gene Ontology | Some mappings |
| Uberon | Anatomy alignment |
| ChEBI | Chemical/metabolite |
| PATO | Quality ontology |

---

## API and Programmatic Access

### OBO Foundry PURL

```bash
curl -L http://purl.obolibrary.org/obo/hp.obo -o hp.obo
```

### Python - Parse HPO

```python
import obonet

# Load HPO
graph = obonet.read_obo('http://purl.obolibrary.org/obo/hp.obo')

# Get term info
term = graph.nodes['HP:0001250']
print(f"Name: {term['name']}")
print(f"Definition: {term.get('def', 'N/A')}")
print(f"Synonyms: {term.get('synonym', [])}")

# Get ancestors
import networkx as nx
ancestors = nx.ancestors(graph, 'HP:0001250')
print(f"Ancestors: {len(ancestors)}")

# Get children
children = list(graph.predecessors('HP:0001250'))
print(f"Children: {children}")
```

### Parse phenotype.hpoa

```python
import pandas as pd

columns = [
    'database_id', 'disease_name', 'qualifier', 'hpo_id',
    'reference', 'evidence', 'onset', 'frequency', 'sex',
    'modifier', 'aspect', 'biocuration'
]

# Skip comment lines
df = pd.read_csv(
    'phenotype.hpoa',
    sep='\t',
    comment='#',
    names=columns,
    header=0
)

# Filter for Marfan syndrome
marfan = df[df['database_id'] == 'OMIM:154700']
print(f"Marfan annotations: {len(marfan)}")

# Get phenotypes by evidence
experimental = df[df['evidence'].isin(['TAS', 'PCS'])]
print(f"Experimental annotations: {len(experimental)}")
```

### Build Disease-Phenotype Matrix

```python
import numpy as np
from collections import defaultdict

# Load annotations
disease_phenotypes = defaultdict(set)
phenotype_diseases = defaultdict(set)

with open('phenotype.hpoa') as f:
    for line in f:
        if line.startswith('#') or line.startswith('database_id'):
            continue
        parts = line.strip().split('\t')
        disease_id = parts[0]
        hpo_id = parts[3]
        qualifier = parts[2]

        if qualifier != 'NOT':
            disease_phenotypes[disease_id].add(hpo_id)
            phenotype_diseases[hpo_id].add(disease_id)

print(f"Diseases: {len(disease_phenotypes)}")
print(f"Phenotypes: {len(phenotype_diseases)}")
```

---

## Phenopackets Integration

HPO is the foundation for the **GA4GH Phenopacket Schema**, a standard for sharing phenotypic data.

### Phenopacket Structure

```json
{
  "id": "patient-001",
  "subject": {
    "id": "patient-001",
    "timeAtLastEncounter": {
      "age": { "iso8601duration": "P25Y" }
    },
    "sex": "MALE"
  },
  "phenotypicFeatures": [
    {
      "type": {
        "id": "HP:0001250",
        "label": "Seizure"
      },
      "onset": {
        "age": { "iso8601duration": "P6M" }
      }
    },
    {
      "type": {
        "id": "HP:0001249",
        "label": "Intellectual disability"
      },
      "excluded": false
    }
  ],
  "diseases": [
    {
      "term": {
        "id": "OMIM:607208",
        "label": "Dravet syndrome"
      }
    }
  ]
}
```

---

## Multi-Language Support

HPO is available in 10+ languages:

| Language | File | Coverage |
|----------|------|----------|
| English | hp.obo | 100% |
| German | hp-international.obo | High |
| Spanish | hp-international.obo | High |
| French | hp-international.obo | High |
| Italian | hp-international.obo | Medium |
| Dutch | hp-international.obo | Medium |
| Portuguese | hp-international.obo | Medium |
| Turkish | hp-international.obo | Medium |
| Japanese | hp-international.obo | Medium |
| Chinese | hp-international.obo | Growing |

---

## Integration with Disease Databases

### OMIM Integration

- ~7,000 OMIM diseases annotated
- IEA evidence from text mining
- TAS evidence from expert curation

### Orphanet Integration

- ~4,000 Orphanet diseases
- Annotations imported from HOOM module
- Regular synchronization

### MONDO Alignment

- HPO terms link to MONDO diseases
- Shared gene associations
- Consistent identifier mapping

---

## License

- **Open Access** - https://hpo.jax.org/app/license
- Free for academic and commercial use
- Attribution recommended

---

## Download

### Ontology Files

**Source:** https://github.com/obophenotype/human-phenotype-ontology/releases/latest

| File | Format | Size | Use Case |
|------|--------|------|----------|
| hp.obo | OBO (text) | 10.4 MB | Standard use (recommended) |
| hp.owl | OWL (XML) | 75.6 MB | Semantic web, reasoners |
| hp.json | JSON | 21.9 MB | Programmatic parsing |
| hp-full.obo | OBO + imports | 18.8 MB | With merged imports |
| phenotype.hpoa | Annotation TSV | 35.0 MB | Disease annotations |
| genes_to_phenotype.txt | TSV | 20.1 MB | Gene-phenotype links |

### PURL Access (Always Current)

```bash
# Stable URLs that redirect to latest
curl -L http://purl.obolibrary.org/obo/hp.obo
curl -L http://purl.obolibrary.org/obo/hp.owl
curl -L http://purl.obolibrary.org/obo/hp.json
```

---

## Data Format

| Format | Description | Encoding |
|--------|-------------|----------|
| Primary | OBO (Open Biological Ontology) | UTF-8 |
| Alternative | OWL (Web Ontology Language) | XML |
| Alternative | JSON | UTF-8 |
| Annotation | GAF 2.2 / HPOA (Tab-separated) | UTF-8 |
| Gene Maps | TSV (Tab-separated values) | UTF-8 |

---

## Data Set Size

| Component | Records | Size |
|-----------|---------|------|
| **HPO Terms** | 13,000+ | 10-150 MB (by format) |
| **Phenotype Annotations** | 156,000+ | 35 MB |
| **Diseases Annotated** | ~9,500 | Included in annotations |
| **Gene Associations** | ~300,000+ | 20-65 MB |
| **Synonyms & Cross-refs** | 50,000+ | Included in OBO/OWL |
| **Full OBO Release** | All content | 45 MB (compressed: 10 MB) |
| **Full OWL Release** | With axioms | 150 MB (compressed: 75 MB) |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "HP:0001250" |
| `name` | string | Entity name | "Seizure" |
| `type` | string | Record type | "phenotype" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `is_a` | HPO Term | N:M |
| `part_of` | HPO Term | N:M |
| `related_to` | Disease / Gene | N:M |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `HP:XXXXXXX` | Human Phenotype Ontology term identifier in 7-digit format | HP:0001250 (Seizure) |
| `phenotype.hpoa` | Annotation file containing disease-phenotype associations | OMIM:154700 - HP:0001166 |
| `evidence` | Code indicating type of evidence supporting annotation (IEA, TAS, PCS) | TAS (expert assertion) |
| `aspect` | HPO branch classification (P=Phenotype, I=Inheritance, C=Course, M=Modifier) | P |
| `qualifier` | Indicates if phenotype is negated (NOT) or observed | NOT (excluded phenotype) |
| `frequency` | How often phenotype occurs in affected individuals | HP:0040281 (Very frequent) |
| `onset` | Age at which phenotype first appears | HP:0003577 (Congenital onset) |
| `is_a` | Ontology relationship indicating subclass/parent relationship | HP:0001250 is_a HP:0012638 |
| `xref` | Cross-reference to external database identifier | UMLS:C0036572 |
| `synonym` | Alternative name for a phenotype term | "Convulsion" for Seizure |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Phenotypic Abnormality | Any observable clinical feature or trait deviating from normal | HP:0000118 branch |
| Mode of Inheritance | Pattern by which a trait is passed through generations | HP:0000005 branch |
| Clinical Modifier | Terms describing severity, laterality, or other qualifications | HP:0012823 branch |
| Clinical Course | Temporal aspects of disease including onset and progression | HP:0031797 branch |
| IEA (Inferred from Electronic Annotation) | Automated annotation from text mining OMIM entries | Evidence code |
| TAS (Traceable Author Statement) | Expert curation with literature citation | Evidence code |
| PCS (Published Clinical Study) | Annotation from published research study | Evidence code |
| Phenopacket | GA4GH standard format for sharing phenotypic data using HPO | Data exchange |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HPO | Human Phenotype Ontology | Phenotype vocabulary |
| HP | Human Phenotype | Term prefix |
| OBO | Open Biological Ontologies | Ontology format |
| OWL | Web Ontology Language | Ontology format |
| IEA | Inferred from Electronic Annotation | Text mining evidence |
| TAS | Traceable Author Statement | Expert curation evidence |
| PCS | Published Clinical Study | Research evidence |
| HPOA | Human Phenotype Ontology Annotation | Annotation file format |
| OMIM | Online Mendelian Inheritance in Man | Disease database |
| ORPHA | Orphanet | Rare disease database |
| CUI | Concept Unique Identifier | UMLS term |
| GA4GH | Global Alliance for Genomics and Health | Standards organization |
| SNOMED CT | Systematized Nomenclature of Medicine Clinical Terms | Clinical terminology |
| MeSH | Medical Subject Headings | NLM vocabulary |
| PURL | Persistent Uniform Resource Locator | Stable URL |

---

## References

- Kohler et al. (2024) "The Human Phenotype Ontology in 2024: phenotypes around the world" Nucleic Acids Research
- HPO GitHub: https://github.com/obophenotype/human-phenotype-ontology
- Phenopackets: https://phenopackets-schema.readthedocs.io/
