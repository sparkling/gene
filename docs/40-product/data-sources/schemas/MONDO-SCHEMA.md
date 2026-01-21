# MONDO Disease Ontology Schema

**Document ID:** MONDO-SCHEMA
**Status:** Final
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-75-RARE-DISEASES](../43-75-RARE-DISEASES.md)

---

## Overview

MONDO (Monarch Disease Ontology) is a unified disease ontology that integrates multiple disease terminologies into a coherent, logic-based structure. Unlike loose cross-references, MONDO provides **precise 1:1 equivalence axioms** validated by OWL reasoning, enabling safe data propagation across OMIM, Orphanet, EFO, DOID, and NCIt.

### Key Resources

| Resource | URL |
|----------|-----|
| **Homepage** | https://mondo.monarchinitiative.org/ |
| **OBO Foundry** | https://obofoundry.org/ontology/mondo.html |
| **GitHub** | https://github.com/monarch-initiative/mondo |
| **Latest Release** | https://github.com/monarch-initiative/mondo/releases/latest |
| **Documentation** | https://mondo.readthedocs.io/ |
| **BioPortal** | https://bioportal.bioontology.org/ontologies/MONDO |

---

## Statistics (v2026-01-06)

| Metric | Count |
|--------|-------|
| **Total diseases** | 25,938 |
| **Human diseases** | 22,977 |
| **Cancer classes** | 4,728 |
| **Infectious diseases** | 1,075 |
| **Mendelian conditions** | 11,639 |
| **Rare diseases** | 15,901 |
| **Database cross-references** | 129,914 |
| **Term definitions** | 17,952 |
| **Exact synonyms** | 73,886 |
| **Narrow synonyms** | 2,554 |
| **Broad synonyms** | 1,420 |
| **Related synonyms** | 30,220 |

---

## Identifier Format

### URI Format (Full)
```
http://purl.obolibrary.org/obo/MONDO_0005015
```

### CURIE Format (Compact)
```
MONDO:0005015
```

### Pattern
```
MONDO:[0-9]{7}
```

**Example:**
- MONDO:0005015 = diabetes mellitus
- MONDO:0007947 = Marfan syndrome
- MONDO:0005148 = type 2 diabetes mellitus

---

## File Formats

### Available Downloads (v2026-01-06)

| File | Size | Description |
|------|------|-------------|
| mondo.owl | 239.3 MB | Full OWL with equivalence axioms |
| mondo.obo | 50.7 MB | OBO format using xrefs |
| mondo.json | 102.7 MB | JSON-LD equivalent to OWL |
| mondo-base.owl | 224.6 MB | Base ontology |
| mondo-base.obo | 47.1 MB | Base in OBO format |
| mondo-simple.owl | 213.0 MB | Simplified hierarchy |
| mondo-rare.owl | 157.5 MB | Rare diseases subset |
| mondo-clingen.owl | 217.4 MB | ClinGen subset |
| mondo-international.owl | 264.3 MB | With translations |
| equivalencies.json | 19.7 MB | Equivalence mappings only |
| equivalencies.obo | 3.8 MB | OBO equivalences |

### Download URLs

```bash
# Stable PURL (redirects to latest)
http://purl.obolibrary.org/obo/mondo.owl
http://purl.obolibrary.org/obo/mondo.obo
http://purl.obolibrary.org/obo/mondo.json

# Direct GitHub release
https://github.com/monarch-initiative/mondo/releases/download/v2026-01-06/mondo.obo
```

---

## OBO Format Structure

### Header Section

```obo
format-version: 1.2
data-version: mondo/releases/2026-01-06/mondo.obo
date: 06:01:2026 14:30
saved-by: nicole
auto-generated-by: OWL API
ontology: mondo
property_value: http://purl.org/dc/elements/1.1/description "Mondo Disease Ontology"
property_value: http://purl.org/dc/elements/1.1/title "Mondo Disease Ontology"
property_value: http://purl.org/dc/terms/license https://creativecommons.org/licenses/by/4.0/
```

### Subset Definitions

```obo
subsetdef: clingen "ClinGen disease subset"
subsetdef: do_inheritance_inconsistent "DO inheritance inconsistent"
subsetdef: gard_rare "GARD rare disease"
subsetdef: mondo_rare "Mondo rare disease"
subsetdef: ordo_clinical_syndrome "Orphanet clinical syndrome"
subsetdef: ordo_disorder "Orphanet disorder"
subsetdef: ordo_etiological_subtype "Orphanet etiological subtype"
subsetdef: ordo_group_of_disorders "Orphanet group of disorders"
subsetdef: ordo_morphological_anomaly "Orphanet morphological anomaly"
```

### Term Structure

```obo
[Term]
id: MONDO:0005015
name: diabetes mellitus
def: "A metabolic disorder characterized by abnormally high blood sugar levels due to diminished production of insulin or resistance to insulin's effects." [NCIT:C2985]
subset: clingen
subset: ncit
synonym: "diabetes" EXACT []
synonym: "DM" EXACT ABBREVIATION []
synonym: "diabetes mellitus syndrome" EXACT []
xref: DOID:9351
xref: EFO:0000400
xref: ICD10CM:E08-E13
xref: ICD9:250
xref: MedDRA:10012601
xref: MESH:D003920
xref: NCIT:C2985
xref: OMIM:125853
xref: Orphanet:101941
xref: SCTID:73211009
xref: UMLS:C0011849
is_a: MONDO:0005066 ! metabolic disease
property_value: exactMatch DOID:9351
property_value: exactMatch Orphanet:101941
property_value: closeMatch ICD10CM:E08-E13
```

### Sample Disease Entry (Marfan Syndrome)

```obo
[Term]
id: MONDO:0007947
name: Marfan syndrome
def: "Marfan syndrome is a systemic disorder of connective tissue with a high degree of clinical variability." [Orphanet:558]
subset: clingen
subset: mondo_rare
subset: ordo_disorder
synonym: "MFS" EXACT ABBREVIATION []
synonym: "Marfan syndrome type 1" EXACT []
synonym: "Marfan's syndrome" EXACT []
xref: DOID:14323
xref: GARD:6975
xref: ICD10CM:Q87.4
xref: MESH:D008382
xref: MedDRA:10026859
xref: NCIT:C34807
xref: OMIM:154700
xref: Orphanet:558
xref: SCTID:19346006
xref: UMLS:C0024796
is_a: MONDO:0019065 ! heritable thoracic aortic aneurysm
is_a: MONDO:0019295 ! rare disease with thoracic involvement
relationship: RO:0004020 HGNC:3603 ! has_material_basis_in_mutation_in FBN1
property_value: exactMatch DOID:14323
property_value: exactMatch OMIM:154700
property_value: exactMatch Orphanet:558
property_value: closeMatch GARD:6975
```

---

## Cross-Reference Types

### Equivalence Axioms

MONDO provides **validated 1:1 equivalences** using OWL reasoning:

| Source | Format | Example |
|--------|--------|---------|
| OMIM | OMIM:[0-9]{6} | OMIM:154700 |
| Orphanet | Orphanet:[0-9]+ | Orphanet:558 |
| DOID | DOID:[0-9]+ | DOID:14323 |
| EFO | EFO:[0-9]{7} | EFO:0000400 |
| NCIt | NCIT:C[0-9]+ | NCIT:C34807 |
| MeSH | MESH:D[0-9]{6} | MESH:D008382 |

### Mapping Predicates

| Predicate | Meaning | Example |
|-----------|---------|---------|
| `exactMatch` | 1:1 equivalence | Orphanet:558 |
| `closeMatch` | Very similar | GARD:6975 |
| `narrowMatch` | More specific | - |
| `broadMatch` | More general | - |
| `relatedMatch` | Related concept | - |

### Cross-Reference Coverage

| Source | Cross-refs | Status |
|--------|------------|--------|
| OMIM | ~16,000 | Validated |
| Orphanet | ~8,000 | Validated |
| DOID | ~10,000 | Validated |
| NCIt | ~8,000 | Validated |
| MeSH | ~12,000 | Curated |
| ICD-10 | ~5,000 | Mapped |
| UMLS | ~15,000 | Mapped |
| GARD | ~6,000 | Curated |
| SCTID (SNOMED) | ~10,000 | Mapped |

---

## Integrated Data Sources

### Primary Sources (with Equivalences)

| Source | Type | Coverage |
|--------|------|----------|
| **OMIM** | Phenotype database | Mendelian diseases |
| **Orphanet (ORDO)** | Rare disease info | Rare diseases |
| **Disease Ontology** | Disease classifications | General diseases |
| **NCIt** | Cancer terminology | Neoplasms |
| **EFO** | Experimental factors | Research diseases |
| **GARD** | Rare disease info | NIH rare diseases |

### Secondary Sources (Cross-References)

| Source | Description |
|--------|-------------|
| MeSH | Medical Subject Headings |
| UMLS | Unified Medical Language System |
| MedDRA | Medical Dictionary for Regulatory Activities |
| ICD-10-CM | International Classification of Diseases |
| ICD-11 Foundation | Latest ICD version |
| SNOMED CT | Clinical terminology |
| MedGen | NCBI medical genetics |
| OncoTree | Cancer classification |
| NANDO | Japanese rare diseases |

---

## OWL Structure

### Core Classes

```
MONDO:0000001 disease
├── MONDO:0021178 disease or disorder
│   ├── MONDO:0005066 metabolic disease
│   │   └── MONDO:0005015 diabetes mellitus
│   ├── MONDO:0005071 nervous system disorder
│   ├── MONDO:0005070 neoplasm
│   │   └── MONDO:0004992 cancer
│   ├── MONDO:0002025 hereditary disease
│   │   └── MONDO:0003847 Mendelian disease
│   └── MONDO:0021145 rare disease
│       └── MONDO:0015286 rare genetic disease
```

### Key Relationship Properties

| Property | URI | Description |
|----------|-----|-------------|
| `is_a` | rdfs:subClassOf | Hierarchical parent |
| `has_material_basis_in_mutation_in` | RO:0004020 | Causal gene |
| `disease_has_feature` | RO:0000053 | Phenotype relationship |
| `characterized_by` | RO:0000057 | Clinical feature |
| `has_onset` | MONDO:0000001 | Age of onset |

### Gene-Disease Relationships

```owl
Class: MONDO:0007947  # Marfan syndrome
  SubClassOf:
    RO:0004020 some HGNC:3603  # has_material_basis_in_mutation_in FBN1
```

---

## SSSOM Mappings

MONDO provides mappings in **SSSOM (Simple Standard for Sharing Ontology Mappings)** format:

### SSSOM File Structure

```tsv
subject_id	predicate_id	object_id	mapping_justification	subject_label	object_label
MONDO:0007947	skos:exactMatch	OMIM:154700	semapv:ManualMappingCuration	Marfan syndrome	MARFAN SYNDROME
MONDO:0007947	skos:exactMatch	Orphanet:558	semapv:ManualMappingCuration	Marfan syndrome	Marfan syndrome
MONDO:0007947	skos:exactMatch	DOID:14323	semapv:ManualMappingCuration	Marfan syndrome	Marfan syndrome
```

### Download SSSOM Mappings

```bash
https://github.com/monarch-initiative/mondo/releases/download/v2026-01-06/mondo_mappings.sssom.tsv
```

---

## Query Examples

### Python - Parse OBO File

```python
import obonet

# Load MONDO
graph = obonet.read_obo('http://purl.obolibrary.org/obo/mondo.obo')

# Get disease term
disease = graph.nodes['MONDO:0007947']
print(f"Name: {disease['name']}")
print(f"Definition: {disease.get('def', 'N/A')}")
print(f"Xrefs: {disease.get('xref', [])}")

# Find OMIM equivalent
for xref in disease.get('xref', []):
    if xref.startswith('OMIM:'):
        print(f"OMIM: {xref}")
```

### SPARQL - Find Rare Genetic Diseases

```sparql
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

SELECT ?disease ?label ?omim ?orphanet WHERE {
  ?disease rdfs:subClassOf* obo:MONDO_0015286 .  # rare genetic disease
  ?disease rdfs:label ?label .

  OPTIONAL {
    ?disease oboInOwl:hasDbXref ?omim .
    FILTER(STRSTARTS(?omim, "OMIM:"))
  }
  OPTIONAL {
    ?disease oboInOwl:hasDbXref ?orphanet .
    FILTER(STRSTARTS(?orphanet, "Orphanet:"))
  }
}
ORDER BY ?label
LIMIT 100
```

### Build OMIM-Orphanet Mapping

```python
import obonet

graph = obonet.read_obo('mondo.obo')

omim_to_orphanet = {}

for node_id, data in graph.nodes(data=True):
    if not node_id.startswith('MONDO:'):
        continue

    xrefs = data.get('xref', [])
    omim_ids = [x for x in xrefs if x.startswith('OMIM:')]
    orphanet_ids = [x for x in xrefs if x.startswith('Orphanet:')]

    if omim_ids and orphanet_ids:
        for omim in omim_ids:
            omim_to_orphanet[omim] = {
                'mondo': node_id,
                'orphanet': orphanet_ids,
                'name': data.get('name', '')
            }

# Result: OMIM → Orphanet mapping via MONDO
print(f"Total mappings: {len(omim_to_orphanet)}")
```

---

## Disease Categories

### By Type

| Category | MONDO Class | Est. Count |
|----------|-------------|------------|
| Neoplasm/Cancer | MONDO:0005070 | 4,728 |
| Infectious disease | MONDO:0005550 | 1,075 |
| Hereditary disease | MONDO:0003847 | 11,639 |
| Rare disease | MONDO:0021145 | 15,901 |
| Congenital disease | MONDO:0000839 | ~2,000 |
| Metabolic disease | MONDO:0005066 | ~2,500 |
| Autoimmune disease | MONDO:0007179 | ~800 |

### By Organ System

| System | MONDO Class | Description |
|--------|-------------|-------------|
| Nervous system | MONDO:0005071 | Neurological |
| Cardiovascular | MONDO:0004995 | Heart/vascular |
| Respiratory | MONDO:0005087 | Lung diseases |
| Musculoskeletal | MONDO:0100070 | Bone/muscle |
| Skin | MONDO:0024458 | Dermatological |
| Immune | MONDO:0005046 | Immunological |
| Endocrine | MONDO:0005151 | Hormonal |
| Hematological | MONDO:0005570 | Blood disorders |

---

## Integration Notes

### Recommended Usage

1. **Use MONDO as pivot** - Map all disease IDs through MONDO
2. **Trust equivalence axioms** - Safe for data propagation
3. **Check subset membership** - Use `subset: mondo_rare` for filtering
4. **Validate xrefs** - Some are equivalences, some are looser matches

### Version Compatibility

| MONDO Version | Breaking Changes |
|---------------|------------------|
| v2024+ | Uses MONDO native IDs exclusively |
| v2022-v2023 | Transition from source IDs |
| Pre-2022 | May use source database IDs |

### License

- **CC BY 4.0** - Free for commercial use with attribution
- Attribution: "Mondo Disease Ontology, Monarch Initiative"

---

## References

- Shefchek et al. (2020) "The Monarch Initiative in 2019" Nucleic Acids Research
- MONDO Documentation: https://mondo.readthedocs.io/
- OBO Foundry: https://obofoundry.org/ontology/mondo.html
