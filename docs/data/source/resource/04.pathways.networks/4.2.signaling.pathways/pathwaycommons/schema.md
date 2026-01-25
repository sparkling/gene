---
id: schema-pathwaycommons
title: "Pathway Commons Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, biopax, pathways, network, sif, owl]
---

# Pathway Commons Schema Documentation

**Document ID:** SCHEMA-PATHWAYCOMMONS
**Version:** 1.0
**Source Version:** PC2 v12 (2023)

---

## TL;DR

Pathway Commons aggregates pathway data from 22+ sources into BioPAX Level 3 OWL format. Provides both rich BioPAX for full pathway detail and simplified SIF (Simple Interaction Format) for network analysis. The PC2 API enables queries by gene, pathway, and graph traversal operations.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Integrated Databases | 22+ | PC website |
| Pathways | 5,772 | PC12 release |
| Interactions | 2,394,814 | Combined |
| Physical Entities | 1,475,993 | Proteins, chemicals, complexes |
| Species | 20+ | Multi-organism |
| Publications | 41,182 | PubMed references |
| UniProt Proteins | 93,710 | Human focus |

---

## Entity Relationship Overview

```
                      ┌─────────────────┐
                      │     Pathway     │
                      │   (BioPAX OWL)  │
                      └────────┬────────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                ▼
    ┌─────────────────┐ ┌──────────┐ ┌─────────────────┐
    │  Biochemical    │ │ Control  │ │   Catalysis     │
    │    Reaction     │ │          │ │                 │
    └────────┬────────┘ └────┬─────┘ └────────┬────────┘
             │               │                │
             ▼               ▼                ▼
    ┌─────────────────┐ ┌──────────────────────────────┐
    │  Physical       │ │        Simple Interaction    │
    │  Entity         │ │        (SIF Export)          │
    └────────┬────────┘ └──────────────────────────────┘
             │
    ┌────────┴────────┬──────────────┐
    │                 │              │
    ▼                 ▼              ▼
┌─────────┐   ┌───────────┐   ┌──────────┐
│ Protein │   │  Small    │   │ Complex  │
│         │   │ Molecule  │   │          │
└─────────┘   └───────────┘   └──────────┘
```

---

## Core BioPAX Entities

### Pathway

**Description:** Container for biological pathway information including reactions and controls.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| uri | URI | Yes | Unique identifier |
| displayName | String | Yes | Human-readable name |
| name | String[] | No | Alternative names |
| organism | BioSource | No | Species |
| pathwayComponent | Process[] | No | Contained reactions/controls |
| pathwayOrder | PathwayStep[] | No | Ordered steps |
| dataSource | Provenance | Yes | Source database |
| xref | Xref[] | No | External references |

### BiochemicalReaction

**Description:** Chemical transformation where one set of molecules converts to another.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| uri | URI | Yes | Unique identifier |
| displayName | String | Yes | Reaction name |
| left | PhysicalEntity[] | Yes | Substrates |
| right | PhysicalEntity[] | Yes | Products |
| conversionDirection | String | No | LEFT-TO-RIGHT, RIGHT-TO-LEFT, REVERSIBLE |
| deltaG | DeltaG[] | No | Gibbs free energy change |
| eCNumber | String[] | No | EC enzyme classification |
| spontaneous | Boolean | No | Is spontaneous |

### Protein

**Description:** Protein entity with sequence and modification information.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| uri | URI | Yes | Unique identifier |
| displayName | String | Yes | Gene symbol or name |
| name | String[] | No | Alternative names |
| entityReference | ProteinReference | No | Canonical sequence |
| feature | EntityFeature[] | No | PTMs, domains |
| cellularLocation | CellularLocation | No | Subcellular location |
| memberPhysicalEntity | Protein[] | No | Generic members |

### SmallMolecule

**Description:** Small molecule or metabolite entity.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| uri | URI | Yes | Unique identifier |
| displayName | String | Yes | Chemical name |
| entityReference | SmallMoleculeReference | No | Chemical structure |
| cellularLocation | CellularLocation | No | Subcellular location |
| chemicalFormula | String | No | Molecular formula |
| molecularWeight | Float | No | MW in Daltons |

### Complex

**Description:** Multi-component protein/molecule assembly.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| uri | URI | Yes | Unique identifier |
| displayName | String | Yes | Complex name |
| component | PhysicalEntity[] | Yes | Member entities |
| componentStoichiometry | Stoichiometry[] | No | Component ratios |

### Control

**Description:** Regulatory relationship between controller and controlled process.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| uri | URI | Yes | Unique identifier |
| controlType | String | Yes | ACTIVATION, INHIBITION |
| controller | Entity[] | Yes | Regulating entity |
| controlled | Process | Yes | Regulated process |
| catalysis | Boolean | No | Is enzymatic |

---

## SIF Format Schema

### Simple Interaction Format

**Description:** Tab-delimited network format for tools like Cytoscape.

| Column | Type | Description |
|--------|------|-------------|
| Source | String | Source entity (gene symbol) |
| Interaction Type | String | Relationship type |
| Target | String | Target entity (gene symbol) |

### SIF Interaction Types

| Type | BioPAX Source | Description |
|------|---------------|-------------|
| INTERACTS_WITH | MolecularInteraction | Physical binding |
| IN_COMPLEX_WITH | Complex | Shared complex membership |
| CONTROLS-STATE-CHANGE-OF | Control | State modification |
| CONTROLS-TRANSPORT-OF | Control (Transport) | Localization change |
| CONTROLS-PHOSPHORYLATION-OF | Control (Modification) | Phosphorylation |
| CONTROLS-EXPRESSION-OF | TemplateReactionRegulation | Transcriptional |
| CATALYSIS-PRECEDES | Catalysis ordering | Enzyme cascade |
| NEIGHBOR_OF | Neighborhood | Same reaction |
| CONSUMPTION-CONTROLLED-BY | Control (Consumption) | Substrate depletion |
| CONTROLS-PRODUCTION-OF | Control (Production) | Product formation |
| REACTS-WITH | Shared reaction | Co-substrates |
| USED-TO-PRODUCE | Reaction flow | Substrate to product |
| CHEMICAL-AFFECTS | Control | Chemical regulation |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /pc2/search | GET | Full-text search |
| /pc2/get | GET | Retrieve by URI |
| /pc2/graph | GET | Graph queries |
| /pc2/traverse | GET | Property traversal |
| /pc2/top_pathways | GET | Top-ranked pathways |

### Search Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| q | String | Search query |
| type | String | BioPAX class filter |
| datasource | String | Source database filter |
| organism | String | Species filter |
| page | Integer | Results page (0-based) |

### Graph Query Types

| Kind | Description | Parameters |
|------|-------------|------------|
| neighborhood | Connected subnetwork | source |
| pathsbetween | Shortest paths between sources | source (multiple) |
| pathsfromto | Directed paths | source, target |
| commonstream | Common regulators/targets | source (multiple), direction |

---

## Data Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| BioPAX Level 3 | .owl | Full RDF/OWL pathway data |
| SIF | .sif/.txt | Simple interaction network |
| Extended SIF | .txt | SIF + pubmed, mediators |
| GSEA GMT | .gmt | Gene set format |
| SBGN-ML | .sbgn | Graphical notation |
| JSON-LD | .json | Linked data |

---

## Sample Record

### BioPAX Pathway (OWL/XML)

```xml
<bp:Pathway rdf:about="http://identifiers.org/reactome/R-HSA-109581">
  <bp:displayName rdf:datatype="xsd:string">Apoptosis</bp:displayName>
  <bp:name rdf:datatype="xsd:string">Programmed Cell Death</bp:name>
  <bp:organism>
    <bp:BioSource rdf:about="http://identifiers.org/taxonomy/9606">
      <bp:displayName>Homo sapiens</bp:displayName>
    </bp:BioSource>
  </bp:organism>
  <bp:dataSource>
    <bp:Provenance rdf:about="http://pathwaycommons.org/pc12/reactome">
      <bp:displayName>Reactome</bp:displayName>
    </bp:Provenance>
  </bp:dataSource>
  <bp:xref>
    <bp:PublicationXref rdf:about="pubmed:12345678">
      <bp:db>PubMed</bp:db>
      <bp:id>12345678</bp:id>
    </bp:PublicationXref>
  </bp:xref>
</bp:Pathway>
```

### SIF Output

```
TP53	INTERACTS_WITH	MDM2
TP53	CONTROLS-PHOSPHORYLATION-OF	BCL2
BRCA1	IN_COMPLEX_WITH	BARD1
BRCA1	IN_COMPLEX_WITH	RAD51
ATM	CONTROLS-STATE-CHANGE-OF	TP53
ATM	CONTROLS-PHOSPHORYLATION-OF	CHEK2
EGFR	CONTROLS-STATE-CHANGE-OF	AKT1
PI3K	CATALYSIS-PRECEDES	AKT1
```

### Extended SIF

```
PARTICIPANT_A	INTERACTION_TYPE	PARTICIPANT_B	INTERACTION_DATA_SOURCE	INTERACTION_PUBMED_ID	PATHWAY_NAMES	MEDIATOR_IDS
TP53	controls-phosphorylation-of	MDM2	reactome;kegg	11900253;16596261	Regulation of TP53 Activity	CHEK2;ATM
EGFR	controls-state-change-of	AKT1	reactome	10508618	PI3K/AKT Signaling	PIK3CA;PIK3R1
```

### Search Response (JSON)

```json
{
  "numHits": 157,
  "maxHitsPerPage": 100,
  "searchHit": [
    {
      "uri": "http://identifiers.org/reactome/R-HSA-109581",
      "biopaxClass": "Pathway",
      "name": ["Apoptosis", "Programmed Cell Death"],
      "dataSource": ["Reactome"],
      "organism": ["Homo sapiens"],
      "pathway": ["Signal Transduction"],
      "numProcesses": 456,
      "numParticipants": 789,
      "excerpt": "Apoptosis is a form of <b>programmed cell death</b>..."
    }
  ]
}
```

### Graph Query Response (JSON)

```json
{
  "source": ["TP53"],
  "target": ["MDM2"],
  "kind": "pathsbetween",
  "format": "EXTENDED_BINARY_SIF",
  "providers": ["reactome", "kegg"],
  "graph": {
    "nodes": [
      {"id": "TP53", "type": "Protein"},
      {"id": "MDM2", "type": "Protein"},
      {"id": "ATM", "type": "Protein"}
    ],
    "edges": [
      {"source": "ATM", "target": "TP53", "type": "controls-phosphorylation-of"},
      {"source": "TP53", "target": "MDM2", "type": "controls-expression-of"}
    ]
  }
}
```

---

## Glossary

| Term | Definition |
|------|------------|
| BioPAX | Biological Pathway Exchange - OWL ontology for pathway data |
| BioPAX Level 3 | Current BioPAX specification with full pathway semantics |
| SIF | Simple Interaction Format - tab-delimited network format |
| OWL | Web Ontology Language - RDF-based semantic format |
| URI | Uniform Resource Identifier - unique entity reference |
| PhysicalEntity | Abstract BioPAX class for molecules/complexes |
| EntityReference | Canonical (generic) form of an entity |
| Control | BioPAX class for regulatory relationships |
| Xref | External database cross-reference |
| Provenance | Data source attribution record |
| PC2 | Pathway Commons version 2 API |
| SBGN | Systems Biology Graphical Notation |

---

## Integrated Data Sources

| Source | ID | Type | Records |
|--------|-------|------|---------|
| Reactome | reactome | Pathway | 2,547 |
| KEGG | kegg | Pathway | 538 |
| WikiPathways | wp | Pathway | 896 |
| PANTHER | panther | Pathway | 177 |
| HumanCyc | humancyc | Pathway | 332 |
| PhosphoSitePlus | psp | PTM | 15,000+ |
| BioGRID | biogrid | PPI | 1.9M |
| IntAct | intact | PPI | 1.1M |
| HPRD | hprd | PPI | 40,000 |
| DIP | dip | PPI | 80,000 |
| BIND | bind | PPI | 170,000 |
| CTD | ctd | Chem-Gene | 500,000+ |

---

## References

1. https://www.pathwaycommons.org/
2. https://www.pathwaycommons.org/pc2/
3. https://www.biopax.org/release/biopax-level3.owl
4. Cerami EG et al. (2011) Pathway Commons, a web resource for biological pathway data. Nucleic Acids Res. 39:D845-850
5. Rodchenkov I et al. (2020) Pathway Commons 2019 Update. Nucleic Acids Res. 48:D489-D497
