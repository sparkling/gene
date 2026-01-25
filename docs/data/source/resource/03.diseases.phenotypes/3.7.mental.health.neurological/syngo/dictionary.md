# SynGO - Data Dictionary

## Overview

This data dictionary documents the schema for SynGO (Synaptic Gene Ontologies).

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | syngo |
| **Name** | SynGO |
| **Parent** | 3.7.mental.health.neurological |
| **Total Fields** | 20+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Gene Annotation

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| gene_symbol | string | 1:1 | Yes | HGNC gene symbol | SYN1 |
| hgnc_id | string | 1:1 | Yes | HGNC identifier | HGNC:11494 |
| uniprot_id | string | 1:1 | No | UniProt accession | P17600 |
| syngo_terms | array | 1:N | Yes | SynGO term annotations | [SYNGO:synapse] |
| evidence_codes | array | 1:N | Yes | Evidence for annotation | [ECO:0005593] |
| pubmed_ids | array | 1:N | No | Supporting literature | [12345678] |
| annotated_by | string | 1:1 | Yes | Curator | Expert curator |
| annotation_date | date | 1:1 | Yes | Annotation date | 2024-01-15 |

### SynGO Term

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| term_id | string | 1:1 | Yes | SynGO term ID | SYNGO:synapse |
| name | string | 1:1 | Yes | Term name | synapse |
| namespace | string | 1:1 | Yes | Ontology branch | cellular_component |
| definition | string | 1:1 | No | Term definition | The junction between... |
| parents | array | 1:N | No | Parent terms | [SYNGO:presynapse] |
| children | array | 1:N | No | Child terms | [SYNGO:postsynapse] |
| go_term | string | 1:1 | No | GO equivalent | GO:0045202 |

### Annotation Evidence

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| annotation_id | string | 1:1 | Yes | Annotation identifier | ANT_0001234 |
| gene_symbol | string | 1:1 | Yes | Gene annotated | SYN1 |
| syngo_term | string | 1:1 | Yes | SynGO term | SYNGO:presynapse |
| evidence_code | string | 1:1 | Yes | Evidence type | ECO:0005593 |
| species | string | 1:1 | Yes | Species | Homo sapiens |
| brain_region | string | 1:1 | No | Brain region | Hippocampus |
| synapse_type | string | 1:1 | No | Synapse type | glutamatergic |
| reference | string | 1:1 | Yes | PMID reference | PMID:12345678 |
| model_system | string | 1:1 | No | Experimental system | Primary neurons |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| SynGO Term | SYNGO:[a-z_]+ | SYNGO:synapse | Ontology term |
| Gene Symbol | Text | SYN1 | HGNC symbol |
| HGNC ID | HGNC:[0-9]+ | HGNC:11494 | Gene identifier |
| UniProt ID | [A-Z0-9]{6,10} | P17600 | Protein identifier |
| GO Term | GO:[0-9]{7} | GO:0045202 | Gene Ontology mapping |
| ECO | ECO:[0-9]{7} | ECO:0005593 | Evidence code |
| PMID | PMID:[0-9]+ | PMID:12345678 | Literature reference |

---

## Enumerations

### Ontology Branches

| Branch | Description | Term Count |
|--------|-------------|------------|
| Cellular Component | Synaptic localization | 90+ |
| Biological Process | Synaptic function | 90+ |
| Presynaptic | Axon terminal | 50+ |
| Postsynaptic | Dendritic spine | 60+ |

### Cellular Component Categories

| Category | Description |
|----------|-------------|
| synapse | Generic synapse |
| presynapse | Presynaptic terminal |
| postsynapse | Postsynaptic compartment |
| synaptic vesicle | Vesicle structure |
| active zone | Release site |
| postsynaptic density | PSD structure |
| synaptic cleft | Extracellular space |
| presynaptic cytosol | Presynaptic cytoplasm |
| postsynaptic cytosol | Postsynaptic cytoplasm |

### Biological Process Categories

| Category | Description |
|----------|-------------|
| synapse organization | Synapse formation/maintenance |
| synaptic signaling | Signal transmission |
| synaptic vesicle cycle | Vesicle trafficking |
| neurotransmitter release | NT exocytosis |
| postsynaptic modulation | Postsynaptic changes |
| synaptic plasticity | LTP/LTD |
| trans-synaptic signaling | Retrograde signaling |

### Evidence Codes

| Code | Name | Description |
|------|------|-------------|
| ECO:0005593 | Biological assay | Functional assay |
| ECO:0005589 | Immunolocalization | Antibody staining |
| ECO:0005644 | Electron microscopy | EM localization |
| ECO:0006063 | Biochemical assay | In vitro biochemistry |
| ECO:0007322 | Imaging | Live imaging |
| ECO:0001255 | Phenotype | Knockout phenotype |

### Synapse Types

| Type | Description |
|------|-------------|
| glutamatergic | Excitatory (glutamate) |
| GABAergic | Inhibitory (GABA) |
| cholinergic | Acetylcholine |
| dopaminergic | Dopamine |
| serotonergic | Serotonin |
| glycinergic | Glycine |
| peptidergic | Neuropeptides |

### Model Systems

| System | Description |
|--------|-------------|
| Primary neurons | Primary cultured neurons |
| Brain slices | Acute brain slices |
| In vivo | Intact animal |
| Cell line | Immortalized cells |
| Organoids | Brain organoids |
| Synaptosomes | Isolated synaptosomes |

### Species

| Species | Description |
|---------|-------------|
| Homo sapiens | Human |
| Mus musculus | Mouse |
| Rattus norvegicus | Rat |
| Drosophila melanogaster | Fruit fly |
| Caenorhabditis elegans | Nematode |

---

## Entity Relationships

### Gene to Annotations
- **Cardinality:** 1:N
- **Description:** Genes have multiple SynGO annotations
- **Key Fields:** gene_symbol, annotation_id

### Term to Hierarchy
- **Cardinality:** N:M
- **Description:** Ontology parent-child relationships
- **Key Fields:** term_id, parents

### Annotation to Evidence
- **Cardinality:** 1:N
- **Description:** Annotations have multiple evidence types
- **Key Fields:** annotation_id, evidence_code

### Gene to GO Terms
- **Cardinality:** N:M
- **Description:** SynGO maps to GO terms
- **Key Fields:** gene_symbol, go_term

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| SynGO | Synaptic Gene Ontologies | Database name |
| GO | Gene Ontology | Parent ontology |
| ECO | Evidence and Conclusion Ontology | Evidence codes |
| PSD | Postsynaptic Density | Structure |
| LTP | Long-Term Potentiation | Plasticity |
| LTD | Long-Term Depression | Plasticity |
| GABA | Gamma-Aminobutyric Acid | Neurotransmitter |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| EM | Electron Microscopy | Technique |
| ICC | Immunocytochemistry | Technique |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| Gene Ontology | GO ID | Parent ontology |
| UniProt | Accession | Protein data |
| HGNC | HGNC ID | Gene nomenclature |
| PubMed | PMID | Literature |
| ECO | ECO ID | Evidence ontology |
| Ensembl | Gene ID | Gene annotation |
| Allen Brain Atlas | Structure ID | Brain regions |

---

## Data Quality Notes

1. **Annotated Genes:** 1,200+ synaptic genes
2. **Annotations:** 5,800+ gene annotations
3. **Supporting Papers:** 2,800+ publications
4. **Ontology Terms:** 180+ synaptic terms
5. **Expert Curation:** Neuroscience expert review
6. **Evidence-Based:** Experimental evidence required
7. **Public API:** Web service access
8. **CC BY 4.0:** Free for research and commercial use

