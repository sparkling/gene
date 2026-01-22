---
id: schema-gene-ontology
title: "Gene Ontology (GO) Schema"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./_index.md)

# Gene Ontology (GO) Schema

**Source:** https://geneontology.org
**Maintained by:** Gene Ontology Consortium
**License:** CC BY 4.0
**Formats:** OBO, OWL, GAF 2.2, GO-CAM JSON

---

## TL;DR

Gene Ontology provides the most widely used structured vocabulary for describing gene product functions across all species. It comprises three ontologies (Molecular Function, Biological Process, Cellular Component), standardized annotation files (GAF 2.2 with 15 columns), evidence codes for annotation provenance, and GO-CAM for causal activity models. Access via REST API, SPARQL endpoint, or bulk downloads from current.geneontology.org.

---

## Overview

The Gene Ontology (GO) is a computational representation of current scientific knowledge about the functions of genes across all species. It provides:

1. **Ontology:** Structured vocabulary of ~47,000 terms organized in three aspects
2. **Annotations:** Gene-to-term associations with evidence codes
3. **GO-CAM:** Causal Activity Models describing how gene products work together

---

## Three Ontologies

| Ontology | Abbreviation | Aspect Code | Root Term ID | Description |
|----------|--------------|-------------|--------------|-------------|
| **Molecular Function** | MF | F | GO:0003674 | Activities performed by gene products at molecular level |
| **Biological Process** | BP | P | GO:0008150 | Larger processes accomplished by multiple molecular activities |
| **Cellular Component** | CC | C | GO:0005575 | Locations where gene products perform functions |

### Molecular Function Examples

| GO ID | Term Name | Definition |
|-------|-----------|------------|
| GO:0003824 | catalytic activity | Catalysis of a biochemical reaction |
| GO:0005488 | binding | Selective, non-covalent interaction with another entity |
| GO:0005215 | transporter activity | Movement of substances across a membrane |
| GO:0140110 | transcription regulator activity | Interacting with DNA or transcription factors |

### Biological Process Examples

| GO ID | Term Name | Definition |
|-------|-----------|------------|
| GO:0008152 | metabolic process | Chemical reactions and pathways |
| GO:0009987 | cellular process | Process carried out at cellular level |
| GO:0023052 | signaling | Cell communication via signals |
| GO:0050896 | response to stimulus | Change in activity triggered by stimulus |

### Cellular Component Examples

| GO ID | Term Name | Definition |
|-------|-----------|------------|
| GO:0016020 | membrane | Lipid bilayer boundary |
| GO:0005634 | nucleus | Membrane-bounded organelle containing chromosomes |
| GO:0005737 | cytoplasm | Contents of cell excluding nucleus |
| GO:0005576 | extracellular region | Space external to plasma membrane |

---

## GAF 2.2 Format (Gene Association File)

### Field Dictionary

| Column | Field Name | Required | Cardinality | Description |
|--------|------------|----------|-------------|-------------|
| 1 | DB | Yes | 1 | Database contributing the annotation (e.g., UniProtKB, MGI) |
| 2 | DB_Object_ID | Yes | 1 | Unique identifier in source database |
| 3 | DB_Object_Symbol | Yes | 1 | Gene/protein symbol |
| 4 | Qualifier | No | 0+ | Annotation qualifier: NOT, contributes_to, colocalizes_with |
| 5 | GO_ID | Yes | 1 | GO term identifier (GO:0000000 format) |
| 6 | DB_Reference | Yes | 1+ | Reference for annotation (PMID, GO_REF, etc.) |
| 7 | Evidence_Code | Yes | 1 | Evidence code (see Evidence Codes table) |
| 8 | With_From | No | 0+ | Additional identifiers supporting evidence |
| 9 | Aspect | Yes | 1 | F (function), P (process), or C (component) |
| 10 | DB_Object_Name | No | 0-1 | Full gene/protein name |
| 11 | DB_Object_Synonym | No | 0+ | Pipe-separated synonyms |
| 12 | DB_Object_Type | Yes | 1 | gene, protein, transcript, complex, etc. |
| 13 | Taxon | Yes | 1-2 | NCBI Taxonomy ID(s) in taxon:##### format |
| 14 | Date | Yes | 1 | Annotation date in YYYYMMDD format |
| 15 | Assigned_By | Yes | 1 | Group that made annotation |

### Sample GAF Record

```tsv
UniProtKB	P04637	TP53	enables	GO:0003700	PMID:8657117	IDA		F	Cellular tumor antigen p53	p53|LFS1	protein	taxon:9606	20230415	UniProt
UniProtKB	P04637	TP53	involved_in	GO:0006915	PMID:9054499	IMP		P	Cellular tumor antigen p53	p53|LFS1	protein	taxon:9606	20230415	UniProt
UniProtKB	P04637	TP53	located_in	GO:0005634	PMID:2649976	IDA		C	Cellular tumor antigen p53	p53|LFS1	protein	taxon:9606	20230415	UniProt
```

---

## Evidence Codes

### Experimental Evidence

| Code | Name | Description | Weight |
|------|------|-------------|--------|
| EXP | Inferred from Experiment | General experimental evidence | High |
| IDA | Inferred from Direct Assay | Direct biochemical assay | High |
| IPI | Inferred from Physical Interaction | Physical interaction evidence | High |
| IMP | Inferred from Mutant Phenotype | Mutant analysis | High |
| IGI | Inferred from Genetic Interaction | Genetic interaction data | High |
| IEP | Inferred from Expression Pattern | Expression evidence | Medium |
| HTP | High Throughput Experiment | Large-scale experimental | Medium |
| HDA | High Throughput Direct Assay | Large-scale direct assay | Medium |
| HMP | High Throughput Mutant Phenotype | Large-scale mutant analysis | Medium |
| HGI | High Throughput Genetic Interaction | Large-scale genetic | Medium |
| HEP | High Throughput Expression Pattern | Large-scale expression | Medium |

### Computational Evidence

| Code | Name | Description | Weight |
|------|------|-------------|--------|
| ISS | Inferred from Sequence Similarity | Sequence-based prediction | Medium |
| ISO | Inferred from Sequence Orthology | Ortholog evidence | Medium |
| ISA | Inferred from Sequence Alignment | Alignment-based | Medium |
| ISM | Inferred from Sequence Model | HMM or profile match | Medium |
| IGC | Inferred from Genomic Context | Genomic context | Low |
| IBA | Inferred from Biological Ancestor | Phylogenetic inference | Medium |
| IBD | Inferred from Biological Descendant | Phylogenetic inference | Medium |
| IKR | Inferred from Key Residues | Key residue analysis | Medium |
| IRD | Inferred from Rapid Divergence | Rapid divergence | Low |
| RCA | Reviewed Computational Analysis | Reviewed computational | Medium |

### Author/Curator Evidence

| Code | Name | Description | Weight |
|------|------|-------------|--------|
| TAS | Traceable Author Statement | Author statement with reference | Medium |
| NAS | Non-traceable Author Statement | Author statement without reference | Low |
| IC | Inferred by Curator | Curator inference | Medium |
| ND | No biological Data available | No data exists | N/A |

### Electronic Evidence

| Code | Name | Description | Weight |
|------|------|-------------|--------|
| IEA | Inferred from Electronic Annotation | Automated annotation | Low |

---

## GO-CAM (Causal Activity Models)

GO-CAM extends traditional GO annotations by providing causal models describing how gene products work together.

### Model Structure

```json
{
  "id": "gomodel:5fce9b7300002411",
  "title": "Insulin signaling pathway",
  "state": "production",
  "date": "2023-04-15",
  "activities": [
    {
      "id": "activity_1",
      "gene_product": "UniProtKB:P06213",
      "molecular_function": "GO:0004716",
      "biological_process": "GO:0008286",
      "cellular_component": "GO:0005886",
      "enabled_by": "UniProtKB:P06213"
    },
    {
      "id": "activity_2",
      "gene_product": "UniProtKB:P42336",
      "molecular_function": "GO:0016303",
      "biological_process": "GO:0048015",
      "cellular_component": "GO:0005829"
    }
  ],
  "causal_relations": [
    {
      "subject": "activity_1",
      "predicate": "RO:0002413",
      "predicate_label": "directly_provides_input_for",
      "object": "activity_2"
    }
  ]
}
```

### GO-CAM Causal Relations (RO Predicates)

| RO ID | Relation Label | Description |
|-------|---------------|-------------|
| RO:0002413 | directly_provides_input_for | Output of one activity is input to another |
| RO:0002629 | directly_positively_regulates | Direct positive regulation |
| RO:0002630 | directly_negatively_regulates | Direct negative regulation |
| RO:0012001 | is_small_molecule_activator_of | Small molecule activates activity |
| RO:0012002 | is_small_molecule_inhibitor_of | Small molecule inhibits activity |
| RO:0002418 | causally_upstream_of_positive_effect | Upstream with positive effect |
| RO:0002419 | causally_upstream_of_negative_effect | Upstream with negative effect |

---

## OBO Format

### Term Structure

```obo
[Term]
id: GO:0008152
name: metabolic process
namespace: biological_process
def: "The chemical reactions and pathways, including anabolism and catabolism, by which living organisms transform chemical substances. Metabolic processes typically transform small molecules, but also include macromolecular processes such as DNA repair and replication, and protein synthesis and degradation." [GOC:go_curators, ISBN:0198547684]
comment: Note that metabolic processes do not include single functions or processes such as protein-protein interactions, protein-nucleic acid interactions etc. If you want to capture these basic activities, please apply the relevant terms from the GO molecular function branch.
subset: goslim_chembl
subset: goslim_generic
subset: goslim_metagenomics
subset: goslim_pir
subset: goslim_plant
synonym: "metabolic process resulting in cell growth" NARROW []
synonym: "metabolism" EXACT []
synonym: "metabolism resulting in cell growth" NARROW []
xref: Wikipedia:Metabolism
is_a: GO:0008150 ! biological_process
```

### Relationship Types in OBO

| Relationship | Description |
|--------------|-------------|
| is_a | Subclass relationship |
| part_of | Component relationship |
| has_part | Contains as component |
| regulates | Regulates another process |
| positively_regulates | Increases rate/extent |
| negatively_regulates | Decreases rate/extent |

---

## OWL Format

```xml
<owl:Class rdf:about="http://purl.obolibrary.org/obo/GO_0008152">
    <rdfs:label rdf:datatype="http://www.w3.org/2001/XMLSchema#string">metabolic process</rdfs:label>
    <rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/GO_0008150"/>
    <obo:IAO_0000115 rdf:datatype="http://www.w3.org/2001/XMLSchema#string">The chemical reactions and pathways, including anabolism and catabolism, by which living organisms transform chemical substances.</obo:IAO_0000115>
    <oboInOwl:hasOBONamespace rdf:datatype="http://www.w3.org/2001/XMLSchema#string">biological_process</oboInOwl:hasOBONamespace>
    <oboInOwl:hasExactSynonym rdf:datatype="http://www.w3.org/2001/XMLSchema#string">metabolism</oboInOwl:hasExactSynonym>
    <oboInOwl:hasDbXref rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Wikipedia:Metabolism</oboInOwl:hasDbXref>
</owl:Class>
```

---

## API Endpoints

### Base URLs

| Service | URL |
|---------|-----|
| GO API | https://api.geneontology.org/api |
| GOlr (Solr) | https://golr-aux.geneontology.org/solr |
| AmiGO | https://amigo.geneontology.org |
| QuickGO | https://www.ebi.ac.uk/QuickGO |

### GO REST API

#### Get GO Term

```bash
curl "https://api.geneontology.org/api/ontology/term/GO:0008150"
```

**Response:**
```json
{
  "goid": "GO:0008150",
  "label": "biological_process",
  "definition": "A biological process represents a specific objective that the organism is genetically programmed to achieve.",
  "comment": "Note that, in addition to having a single inherent discrete beginning and end, a biological process must have more than one distinct step.",
  "synonyms": [],
  "namespace": "biological_process",
  "alt_ids": [],
  "xrefs": [],
  "subsets": ["goslim_generic", "goslim_pir"],
  "replaced_by": [],
  "consider": []
}
```

#### Search Terms

```bash
curl "https://api.geneontology.org/api/search/entity?q=apoptosis&category=ontology_class"
```

**Response:**
```json
{
  "numFound": 156,
  "docs": [
    {
      "id": "GO:0006915",
      "label": "apoptotic process",
      "category": "ontology_class",
      "description": "A programmed cell death process which begins when a cell receives an internal (e.g. DNA damage) or external signal (e.g. an extracellular death ligand)..."
    }
  ]
}
```

#### Get Annotations for Gene

```bash
curl "https://api.geneontology.org/api/bioentity/gene/HGNC:11998/function"
```

**Response:**
```json
{
  "associations": [
    {
      "subject": {
        "id": "HGNC:11998",
        "label": "TP53",
        "taxon": {
          "id": "NCBITaxon:9606",
          "label": "Homo sapiens"
        }
      },
      "object": {
        "id": "GO:0003700",
        "label": "DNA-binding transcription factor activity",
        "taxon": null
      },
      "evidence": {
        "type": "IDA",
        "has_supporting_reference": ["PMID:8657117"]
      },
      "provided_by": "UniProt"
    }
  ]
}
```

### GOlr (Solr) Queries

```bash
# Search annotations
curl "https://golr-aux.geneontology.org/solr/select?q=TP53&fq=document_category:annotation&wt=json"

# Get term ancestors
curl "https://golr-aux.geneontology.org/solr/select?q=GO:0006915&fq=document_category:ontology_class&fl=isa_closure,isa_closure_label&wt=json"

# Filter by evidence code
curl "https://golr-aux.geneontology.org/solr/select?q=*:*&fq=document_category:annotation&fq=evidence_type:IDA&fq=taxon:NCBITaxon:9606&wt=json"
```

---

## SPARQL Endpoint

**Endpoint:** https://sparql.geneontology.org/sparql

### Example Queries

#### Get All GO Terms in Biological Process

```sparql
PREFIX GO: <http://purl.obolibrary.org/obo/GO_>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX obo: <http://purl.obolibrary.org/obo/>

SELECT ?term ?label ?definition
WHERE {
  ?term rdfs:subClassOf* GO:0008150 .
  ?term rdfs:label ?label .
  OPTIONAL { ?term obo:IAO_0000115 ?definition }
}
LIMIT 100
```

#### Get Annotations for Human Gene

```sparql
PREFIX GO: <http://purl.obolibrary.org/obo/GO_>
PREFIX enabled_by: <http://purl.obolibrary.org/obo/RO_0002333>
PREFIX UniProtKB: <http://identifiers.org/uniprot/>

SELECT ?model ?activity ?goterm ?gotermLabel
WHERE {
  ?model a <http://purl.obolibrary.org/obo/GO_0003674> .
  ?activity enabled_by: UniProtKB:P04637 .
  ?activity a ?goterm .
  ?goterm rdfs:label ?gotermLabel .
}
```

#### Count Annotations by Evidence Code

```sparql
PREFIX eco: <http://purl.obolibrary.org/obo/ECO_>
PREFIX goa: <http://geneontology.org/lego/>

SELECT ?evidenceCode (COUNT(?annotation) AS ?count)
WHERE {
  ?annotation goa:evidence ?evidence .
  ?evidence a ?evidenceCode .
}
GROUP BY ?evidenceCode
ORDER BY DESC(?count)
```

---

## Bulk Downloads

**FTP:** https://current.geneontology.org/

### Ontology Files

| File | Description | Size (approx) |
|------|-------------|---------------|
| go-basic.obo | Filtered GO (no obsolete, relationships) | 35 MB |
| go.obo | Full GO ontology | 45 MB |
| go.owl | OWL format | 150 MB |
| go-plus.owl | GO with external axioms | 500 MB |
| goslim_generic.obo | Generic GO slim | 50 KB |
| goslim_chembl.obo | ChEMBL GO slim | 25 KB |

### Annotation Files (GAF)

| File | Species | Annotations (approx) |
|------|---------|---------------------|
| goa_human.gaf.gz | Homo sapiens | 800,000 |
| goa_mouse.gaf.gz | Mus musculus | 600,000 |
| goa_rat.gaf.gz | Rattus norvegicus | 400,000 |
| goa_fly.gaf.gz | Drosophila melanogaster | 300,000 |
| goa_worm.gaf.gz | Caenorhabditis elegans | 200,000 |
| goa_yeast.gaf.gz | Saccharomyces cerevisiae | 100,000 |
| goa_zebrafish.gaf.gz | Danio rerio | 200,000 |
| goa_uniprot_all.gaf.gz | All UniProt species | 500,000,000+ |

### GO-CAM Models

| File | Description |
|------|-------------|
| go-cams.json | All GO-CAM models in JSON |
| go-cams.owl | All GO-CAM models in OWL |
| noctua-models/ | Individual model files |

---

## Cross-References

| Database | ID Format | Relationship |
|----------|-----------|--------------|
| UniProt | P#####, Q##### | Primary protein identifiers in annotations |
| NCBI Gene | Integer | Gene identifiers |
| Ensembl | ENSG########### | Transcript/protein IDs |
| ChEBI | CHEBI:##### | Chemical entity cross-references |
| PubMed | Integer | Literature references |
| Reactome | R-HSA-###### | Pathway cross-references |
| EC | #.#.#.# | Enzyme classification |
| InterPro | IPR###### | Protein domain annotations |
| KEGG | K##### | Pathway/orthology mappings |
| PDB | #### | Structure-based annotations |

---

## AmiGO Browser

**URL:** https://amigo.geneontology.org

### Search Types

| Search Type | Description | Example Query |
|-------------|-------------|---------------|
| Ontology | Browse/search GO terms | "apoptosis" |
| Annotations | Search gene-to-GO associations | "TP53 human" |
| Gene Products | Search specific genes/proteins | "UniProtKB:P04637" |
| Bioentities | Combined gene/protein search | "BRCA1" |

### URL Patterns

```
# View GO term
https://amigo.geneontology.org/amigo/term/GO:0006915

# View gene product
https://amigo.geneontology.org/amigo/gene_product/UniProtKB:P04637

# Search annotations
https://amigo.geneontology.org/amigo/search/annotation?q=TP53

# Term enrichment (via PANTHER)
https://amigo.geneontology.org/amigo/tools/enrichment
```

---

## Statistics (as of 2025)

| Metric | Count |
|--------|-------|
| GO Terms (total) | ~47,000 |
| Molecular Function terms | ~12,000 |
| Biological Process terms | ~30,000 |
| Cellular Component terms | ~4,500 |
| Annotations (all species) | ~500,000,000 |
| Human annotations | ~800,000 |
| Contributing groups | 40+ |
| Species with annotations | 4,000+ |
| GO-CAM models | 1,500+ |

---

## Data Quality Tiers

| Tier | Evidence Codes | Description |
|------|---------------|-------------|
| Gold | IDA, IPI, IMP, IGI, IEP, EXP | Direct experimental evidence |
| Silver | ISS, ISO, ISA, IBA, IBD, TAS | Inferred/author evidence |
| Bronze | IEA | Electronic annotation |

---

## GO Slims

Pre-defined subsets for high-level analysis:

| Slim | Terms | Use Case |
|------|-------|----------|
| goslim_generic | ~150 | General-purpose overview |
| goslim_chembl | ~100 | Drug discovery |
| goslim_metagenomics | ~150 | Metagenomics studies |
| goslim_plant | ~200 | Plant biology |
| goslim_yeast | ~180 | Yeast research |
| goslim_pir | ~300 | Protein Information Resource |

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 2024-01 | Jan 2024 | Monthly ontology release |
| GAF 2.2 | 2019 | Current annotation format |
| GAF 2.1 | 2014 | Previous annotation format |
| OBO 1.4 | 2017 | Current OBO format |
| GO-CAM | 2019 | Causal activity models introduced |

---

## Download

### Ontology Files

**Source:** https://current.geneontology.org/

| File | Format | Use Case |
|------|--------|----------|
| go-basic.obo | OBO (35 MB) | Standard use (recommended) |
| go.obo | OBO (45 MB) | Full ontology with relationships |
| go.owl | OWL (150 MB) | Semantic web, reasoner-compatible |
| go-plus.owl | OWL (500 MB) | Full with external axioms |
| go.json | JSON (21.9 MB) | Programmatic parsing |

### Annotation Files

| Species | File | Records |
|---------|------|---------|
| Human | goa_human.gaf.gz | ~800,000 |
| Mouse | goa_mouse.gaf.gz | ~600,000 |
| All UniProt | goa_uniprot_all.gaf.gz | 500,000,000+ |
| GO-CAM Models | go-cams.json, go-cams.owl | 1,500+ |

---

## Data Format

| Format | Description | Encoding |
|--------|-------------|----------|
| Primary | OBO (text-based, human-readable) | UTF-8 |
| Alternative | OWL (XML, semantic web) | XML |
| Alternative | JSON | UTF-8 |
| Annotation | GAF 2.2 (tab-separated) | UTF-8 |
| Query | SPARQL | Standard RDF |
| API Response | JSON | application/json |

---

## Data Set Size

| Component | Size | Records |
|-----------|------|---------|
| **go-basic.obo** | 35 MB | ~47,000 terms |
| **go.obo (full)** | 45 MB | ~47,000 + relationships |
| **go.owl** | 150 MB | ~47,000 terms + axioms |
| **go-plus.owl** | 500 MB | With external imports |
| **GAF (human)** | ~50 MB (gzipped) | ~800,000 annotations |
| **GAF (all species)** | ~5 GB (gzipped) | 500,000,000+ annotations |
| **GO-CAM models** | ~100 MB | 1,500+ causal models |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "GO:0008150" |
| `name` | string | Entity name | "Biological Process" |
| `type` | string | Record type | "ontology_term" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| Gene Ontology | CC BY 4.0 | Yes |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `GO:XXXXXXX` | Gene Ontology term identifier in 7-digit format | GO:0008152 |
| `DB_Object_ID` | Unique identifier in source database | P04637 |
| `Evidence_Code` | Code indicating type of evidence for annotation | IDA, IMP, IEA |
| `Aspect` | Ontology branch indicator (F, P, or C) | F (Molecular Function) |
| `Qualifier` | Annotation modifier (NOT, enables, part_of) | enables |
| `With_From` | Supporting identifiers for evidence | PMID:12345 |
| `Taxon` | NCBI Taxonomy identifier | taxon:9606 |
| `Assigned_By` | Annotation source organization | UniProt, MGI |
| `gomodel` | GO-CAM model identifier | gomodel:5fce9b7300002411 |
| `RO:XXXXXXX` | Relations Ontology identifier for causal relations | RO:0002413 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Molecular Function | Activity performed by gene product at molecular level | Aspect F |
| Biological Process | Larger processes accomplished by multiple activities | Aspect P |
| Cellular Component | Locations where gene products perform functions | Aspect C |
| GO-CAM | Gene Ontology Causal Activity Model | Extended annotation format |
| Evidence Code | Standardized code indicating annotation provenance | IDA, IMP, IEA |
| GO Slim | Reduced subset of GO terms for high-level analysis | goslim_generic |
| Annotation | Association between gene product and GO term | GAF record |
| Direct Assay | Biochemical experiment directly measuring function | IDA evidence |
| Electronic Annotation | Computationally derived annotation | IEA evidence |
| Causal Relation | Directed relationship between activities in GO-CAM | directly_provides_input_for |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GO | Gene Ontology | Functional annotation standard |
| GAF | Gene Association File | Annotation file format |
| MF | Molecular Function | GO ontology aspect |
| BP | Biological Process | GO ontology aspect |
| CC | Cellular Component | GO ontology aspect |
| GO-CAM | GO Causal Activity Model | Pathway-level annotation |
| IDA | Inferred from Direct Assay | Experimental evidence |
| IMP | Inferred from Mutant Phenotype | Experimental evidence |
| IEA | Inferred from Electronic Annotation | Computational evidence |
| TAS | Traceable Author Statement | Author-based evidence |
| RO | Relations Ontology | Predicate vocabulary |
| OBO | Open Biological Ontologies | Ontology format |
| OWL | Web Ontology Language | Semantic web format |
| GOlr | GO Lucene/Solr | Search index |
| AmiGO | GO Browser | Web interface |
| PANTHER | Protein Analysis Through Evolutionary Relationships | Classification system |

---

## References

1. Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. Genetics, 224(1), iyad031. PMID: 36866529

2. Thomas, P.D., et al. (2022). PANTHER: Making genome-scale phylogenetics accessible to all. Protein Science, 31(1), 8-22. PMID: 34717010

3. Carbon, S., et al. (2021). The Gene Ontology resource: enriching a GOld mine. Nucleic Acids Research, 49(D1), D325-D334. PMID: 33290552

4. Mungall, C.J., et al. (2019). GOCAM: A Gene Ontology Causal Activity Model. Journal of Biomedical Semantics, 10(1), 1-13.

---

*Schema version: 1.0*
*Last updated: 2025*
*For current information, visit https://geneontology.org*
