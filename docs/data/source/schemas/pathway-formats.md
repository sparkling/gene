# Pathway Database Formats - Detailed Technical Reference

This document provides comprehensive technical specifications for pathway database formats, extending the existing pathway-target-data-models.md with detailed schema definitions, example data structures, and parsing code snippets.

---

## Table of Contents

1. [Reactome Data Model](#1-reactome-data-model)
2. [KEGG Data Model](#2-kegg-data-model)
3. [WikiPathways GPML Format](#3-wikipathways-gpml-format)
4. [Gene Ontology Structure](#4-gene-ontology-structure)
5. [BioPAX (Biological Pathway Exchange)](#5-biopax-biological-pathway-exchange)
6. [SBML (Systems Biology Markup Language)](#6-sbml-systems-biology-markup-language)
7. [PSI-MI (Molecular Interactions)](#7-psi-mi-molecular-interactions)
8. [CX Format (NDEx)](#8-cx-format-ndex)

---

## 1. Reactome Data Model

### 1.1 Neo4j Graph Schema Overview

Reactome uses Neo4j as its graph database backend. The schema follows an object-oriented design with inheritance hierarchies.

#### Node Type Hierarchy

```
DatabaseObject (base class)
├── Event
│   ├── Pathway
│   ├── ReactionLikeEvent
│   │   ├── Reaction
│   │   ├── BlackBoxEvent
│   │   ├── Polymerisation
│   │   ├── Depolymerisation
│   │   └── FailedReaction
│   └── TopLevelPathway
├── PhysicalEntity
│   ├── Complex
│   ├── EntitySet
│   │   ├── CandidateSet
│   │   ├── DefinedSet
│   │   └── OpenSet
│   ├── GenomeEncodedEntity
│   │   └── EntityWithAccessionedSequence
│   ├── Drug
│   ├── SimpleEntity
│   ├── Polymer
│   └── OtherEntity
├── ReferenceEntity
│   ├── ReferenceSequence
│   │   ├── ReferenceGeneProduct
│   │   ├── ReferenceIsoform
│   │   └── ReferenceDNASequence
│   ├── ReferenceMolecule
│   ├── ReferenceGroup
│   └── ReferenceTherapeutic
├── Regulation
│   ├── PositiveRegulation
│   │   ├── PositiveGeneExpressionRegulation
│   │   └── Requirement
│   └── NegativeRegulation
│       └── NegativeGeneExpressionRegulation
├── CatalystActivity
├── GO_Term
│   ├── GO_BiologicalProcess
│   ├── GO_CellularComponent
│   └── GO_MolecularFunction
└── Species
```

### 1.2 Core Node Types

#### Pathway Node

```cypher
// Cypher query to retrieve Pathway node structure
MATCH (p:Pathway {stId: 'R-HSA-109581'})
RETURN p
```

**Properties:**

| Property | Type | Description | Example |
|----------|------|-------------|---------|
| `dbId` | Long | Internal database ID | `109581` |
| `stId` | String | Stable identifier | `R-HSA-109581` |
| `displayName` | String | Human-readable name | `Apoptosis` |
| `name` | List[String] | All names including synonyms | `["Apoptosis", "Programmed cell death"]` |
| `speciesName` | String | Species | `Homo sapiens` |
| `isInDisease` | Boolean | Disease-related pathway | `false` |
| `isInferred` | Boolean | Computationally inferred | `false` |
| `releaseDate` | String | First release date | `2004-09-20` |
| `schemaClass` | String | Node type | `Pathway` |
| `hasDiagram` | Boolean | Has pathway diagram | `true` |
| `doi` | String | DOI identifier | `10.3180/R-HSA-109581` |

#### Reaction Node

```cypher
MATCH (r:Reaction {stId: 'R-HSA-109606'})
RETURN r
```

**Properties:**

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Long | Internal database ID |
| `stId` | String | Stable identifier |
| `displayName` | String | Reaction name |
| `category` | String | Reaction category (transition, binding, dissociation, etc.) |
| `isChimeric` | Boolean | Contains entities from multiple species |
| `systematicName` | String | Systematic biochemical name |
| `releaseStatus` | String | RELEASED, NEW, UPDATED |

#### PhysicalEntity Node

```cypher
MATCH (pe:PhysicalEntity {stId: 'R-HSA-109606'})
RETURN pe
```

**Properties:**

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Long | Internal database ID |
| `stId` | String | Stable identifier |
| `displayName` | String | Entity name |
| `schemaClass` | String | Subtype (Complex, EntityWithAccessionedSequence, etc.) |
| `speciesName` | String | Species |
| `startCoordinate` | Integer | Protein start position (for EWAS) |
| `endCoordinate` | Integer | Protein end position (for EWAS) |

#### ReferenceEntity Node

```cypher
MATCH (re:ReferenceEntity {identifier: 'P04637'})
RETURN re
```

**Properties:**

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Long | Internal database ID |
| `identifier` | String | External database ID (e.g., UniProt) |
| `databaseName` | String | Source database name |
| `url` | String | External URL |
| `name` | List[String] | Names from reference database |
| `geneName` | List[String] | Gene names |
| `secondaryIdentifier` | List[String] | Secondary IDs |
| `description` | String | Entity description |

### 1.3 Relationship Types

#### Event Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `hasEvent` | Pathway | Event | Sub-events within pathway |
| `precedingEvent` | Event | Event | Temporal ordering |
| `followingEvent` | Event | Event | Reverse temporal ordering |
| `inferredTo` | Event | Event | Orthology inference |

#### Reaction Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `input` | ReactionLikeEvent | PhysicalEntity | Reaction inputs |
| `output` | ReactionLikeEvent | PhysicalEntity | Reaction outputs |
| `catalystActivity` | ReactionLikeEvent | CatalystActivity | Enzymatic catalysis |
| `regulatedBy` | ReactionLikeEvent | Regulation | Regulatory relationships |
| `normalReaction` | FailedReaction | Reaction | Normal version of failed reaction |

#### PhysicalEntity Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `hasComponent` | Complex | PhysicalEntity | Complex components |
| `hasMember` | EntitySet | PhysicalEntity | Set members |
| `hasCandidate` | CandidateSet | PhysicalEntity | Candidate members |
| `referenceEntity` | PhysicalEntity | ReferenceEntity | Link to reference DB |
| `hasModifiedResidue` | EWAS | AbstractModifiedResidue | PTM annotations |
| `compartment` | PhysicalEntity | Compartment | Cellular location |

#### CatalystActivity Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `physicalEntity` | CatalystActivity | PhysicalEntity | Catalyzing entity |
| `activity` | CatalystActivity | GO_MolecularFunction | GO activity term |
| `activeUnit` | CatalystActivity | PhysicalEntity | Active subunit |

### 1.4 Complete Cypher Query Examples

```cypher
// Get pathway with all sub-events
MATCH (p:Pathway {stId: 'R-HSA-109581'})-[:hasEvent*]->(e:Event)
RETURN p.displayName AS Pathway, collect(DISTINCT e.displayName) AS Events

// Get reaction with inputs, outputs, and catalysts
MATCH (r:Reaction {stId: 'R-HSA-69541'})
OPTIONAL MATCH (r)-[:input]->(i:PhysicalEntity)
OPTIONAL MATCH (r)-[:output]->(o:PhysicalEntity)
OPTIONAL MATCH (r)-[:catalystActivity]->(ca:CatalystActivity)-[:physicalEntity]->(c:PhysicalEntity)
RETURN r.displayName AS Reaction,
       collect(DISTINCT i.displayName) AS Inputs,
       collect(DISTINCT o.displayName) AS Outputs,
       collect(DISTINCT c.displayName) AS Catalysts

// Get protein with all pathways
MATCH (re:ReferenceEntity {identifier: 'P04637'})<-[:referenceEntity]-(pe:PhysicalEntity)
MATCH (pe)<-[:hasComponent|hasMember|hasCandidate|input|output*]-(e:Event)
MATCH (e)<-[:hasEvent*]-(p:TopLevelPathway)
RETURN DISTINCT p.displayName AS Pathway, p.stId AS PathwayId

// Get all proteins in a pathway
MATCH (p:Pathway {stId: 'R-HSA-109581'})-[:hasEvent*]->(r:ReactionLikeEvent)
MATCH (r)-[:input|output]->(pe:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
WHERE re.databaseName = 'UniProt'
RETURN DISTINCT re.identifier AS UniProtID, re.geneName AS GeneName
```

### 1.5 JSON API Response Structures

#### Pathway Detail Response

```json
{
  "dbId": 109581,
  "stId": "R-HSA-109581",
  "displayName": "Apoptosis",
  "speciesName": "Homo sapiens",
  "schemaClass": "Pathway",
  "releaseDate": "2004-09-20",
  "isInDisease": false,
  "isInferred": false,
  "hasDiagram": true,
  "hasEHLD": true,
  "diagramWidth": 1500,
  "diagramHeight": 1200,
  "stIdVersion": "R-HSA-109581.4",
  "doi": "10.3180/R-HSA-109581",
  "summation": [
    {
      "dbId": 76018,
      "displayName": "Summation",
      "text": "Apoptosis is a distinct form of cell death...",
      "schemaClass": "Summation"
    }
  ],
  "hasEvent": [
    {
      "dbId": 5357769,
      "stId": "R-HSA-5357769",
      "displayName": "Caspase activation via extrinsic apoptotic signalling pathway",
      "schemaClass": "Pathway"
    },
    {
      "dbId": 109606,
      "stId": "R-HSA-109606",
      "displayName": "Intrinsic Pathway for Apoptosis",
      "schemaClass": "Pathway"
    }
  ],
  "literatureReference": [
    {
      "dbId": 76019,
      "displayName": "Apoptosis--the p53 network.",
      "pubMedIdentifier": 12505355,
      "title": "Apoptosis--the p53 network.",
      "journal": "J. Cell. Sci.",
      "year": 2003,
      "volume": "116",
      "pages": "4077-85",
      "authors": ["Vousden KH", "Lu X"]
    }
  ],
  "compartment": [
    {
      "dbId": 984,
      "stId": "R-HSA-70101",
      "displayName": "cytosol",
      "accession": "GO:0005829"
    }
  ],
  "goBiologicalProcess": {
    "dbId": 8868,
    "displayName": "apoptotic process",
    "accession": "GO:0006915",
    "databaseName": "GO"
  },
  "created": {
    "dateTime": "2004-09-20 00:00:00",
    "author": [
      {
        "dbId": 111,
        "displayName": "Bhalla, US"
      }
    ]
  },
  "modified": {
    "dateTime": "2024-01-15 12:30:00",
    "author": [
      {
        "dbId": 222,
        "displayName": "Jassal, B"
      }
    ]
  }
}
```

#### Reaction Detail Response

```json
{
  "dbId": 69541,
  "stId": "R-HSA-69541",
  "displayName": "Cytochrome c binds APAF1",
  "speciesName": "Homo sapiens",
  "schemaClass": "Reaction",
  "category": "binding",
  "isChimeric": false,
  "input": [
    {
      "dbId": 69523,
      "stId": "R-HSA-69523",
      "displayName": "CYCS [cytosol]",
      "schemaClass": "EntityWithAccessionedSequence",
      "stoichiometry": 1
    },
    {
      "dbId": 69542,
      "stId": "R-HSA-69542",
      "displayName": "APAF1 [cytosol]",
      "schemaClass": "EntityWithAccessionedSequence",
      "stoichiometry": 1
    }
  ],
  "output": [
    {
      "dbId": 69543,
      "stId": "R-HSA-69543",
      "displayName": "CYCS:APAF1 [cytosol]",
      "schemaClass": "Complex",
      "stoichiometry": 1
    }
  ],
  "catalystActivity": [],
  "positivelyRegulatedBy": [],
  "negativelyRegulatedBy": [],
  "literatureReference": [
    {
      "pubMedIdentifier": 9358771
    }
  ]
}
```

#### PhysicalEntity Detail Response

```json
{
  "dbId": 69523,
  "stId": "R-HSA-69523",
  "displayName": "CYCS [cytosol]",
  "speciesName": "Homo sapiens",
  "schemaClass": "EntityWithAccessionedSequence",
  "startCoordinate": 1,
  "endCoordinate": 105,
  "referenceEntity": {
    "dbId": 69520,
    "identifier": "P99999",
    "databaseName": "UniProt",
    "name": ["Cytochrome c"],
    "geneName": ["CYCS"],
    "url": "https://www.uniprot.org/uniprot/P99999",
    "secondaryIdentifier": ["CYCS_HUMAN"],
    "description": "Cytochrome c, somatic"
  },
  "compartment": [
    {
      "displayName": "cytosol",
      "accession": "GO:0005829"
    }
  ],
  "hasModifiedResidue": [
    {
      "dbId": 69521,
      "displayName": "AcLys-Lys1",
      "coordinate": 1,
      "psiMod": {
        "identifier": "00723",
        "databaseName": "PSI-MOD"
      }
    }
  ],
  "crossReference": [
    {
      "databaseName": "NCBI Gene",
      "identifier": "54205"
    },
    {
      "databaseName": "OMIM",
      "identifier": "123970"
    }
  ]
}
```

### 1.6 BioPAX Export Format

Reactome exports data in BioPAX Level 3 OWL format. Example structure:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
         xmlns:bp="http://www.biopax.org/release/biopax-level3.owl#"
         xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
         xml:base="http://www.reactome.org/biopax/">

  <bp:Pathway rdf:ID="Pathway109581">
    <bp:displayName rdf:datatype="xsd:string">Apoptosis</bp:displayName>
    <bp:name rdf:datatype="xsd:string">Apoptosis</bp:name>
    <bp:organism rdf:resource="#BioSource1"/>
    <bp:pathwayComponent rdf:resource="#BiochemicalReaction69541"/>
    <bp:pathwayComponent rdf:resource="#Pathway109606"/>
    <bp:xref rdf:resource="#UnificationXref_R-HSA-109581"/>
    <bp:comment rdf:datatype="xsd:string">
      Apoptosis is a distinct form of cell death...
    </bp:comment>
  </bp:Pathway>

  <bp:BiochemicalReaction rdf:ID="BiochemicalReaction69541">
    <bp:displayName rdf:datatype="xsd:string">Cytochrome c binds APAF1</bp:displayName>
    <bp:left rdf:resource="#Protein69523"/>
    <bp:left rdf:resource="#Protein69542"/>
    <bp:right rdf:resource="#Complex69543"/>
    <bp:conversionDirection rdf:datatype="xsd:string">LEFT-TO-RIGHT</bp:conversionDirection>
  </bp:BiochemicalReaction>

  <bp:Protein rdf:ID="Protein69523">
    <bp:displayName rdf:datatype="xsd:string">CYCS</bp:displayName>
    <bp:entityReference rdf:resource="#ProteinReference_P99999"/>
    <bp:cellularLocation rdf:resource="#CellularLocationVocabulary_cytosol"/>
  </bp:Protein>

  <bp:ProteinReference rdf:ID="ProteinReference_P99999">
    <bp:displayName rdf:datatype="xsd:string">Cytochrome c</bp:displayName>
    <bp:organism rdf:resource="#BioSource1"/>
    <bp:xref rdf:resource="#UnificationXref_P99999"/>
    <bp:name rdf:datatype="xsd:string">CYCS</bp:name>
  </bp:ProteinReference>

  <bp:UnificationXref rdf:ID="UnificationXref_P99999">
    <bp:db rdf:datatype="xsd:string">UniProt</bp:db>
    <bp:id rdf:datatype="xsd:string">P99999</bp:id>
  </bp:UnificationXref>

  <bp:BioSource rdf:ID="BioSource1">
    <bp:displayName rdf:datatype="xsd:string">Homo sapiens</bp:displayName>
    <bp:xref rdf:resource="#UnificationXref_9606"/>
  </bp:BioSource>

</rdf:RDF>
```

### 1.7 Python Parsing Code

```python
import requests
from typing import Dict, List, Optional

class ReactomeParser:
    """Parser for Reactome API responses."""

    BASE_URL = "https://reactome.org/ContentService"

    def get_pathway(self, pathway_id: str) -> Dict:
        """Fetch pathway details."""
        url = f"{self.BASE_URL}/data/query/{pathway_id}"
        response = requests.get(url, headers={"Accept": "application/json"})
        return response.json()

    def get_pathway_events(self, pathway_id: str) -> List[Dict]:
        """Get all events in a pathway."""
        url = f"{self.BASE_URL}/data/pathway/{pathway_id}/containedEvents"
        response = requests.get(url, headers={"Accept": "application/json"})
        return response.json()

    def get_reaction_participants(self, reaction_id: str) -> Dict:
        """Get inputs, outputs, and catalysts for a reaction."""
        url = f"{self.BASE_URL}/data/query/{reaction_id}"
        response = requests.get(url, headers={"Accept": "application/json"})
        data = response.json()
        return {
            "inputs": data.get("input", []),
            "outputs": data.get("output", []),
            "catalysts": [
                ca.get("physicalEntity", {})
                for ca in data.get("catalystActivity", [])
            ]
        }

    def map_uniprot_to_pathways(self, uniprot_id: str) -> List[Dict]:
        """Map UniProt ID to pathways."""
        url = f"{self.BASE_URL}/data/mapping/UniProt/{uniprot_id}/pathways"
        response = requests.get(url, headers={"Accept": "application/json"})
        return response.json()

    def parse_biopax(self, biopax_content: str) -> Dict:
        """Parse BioPAX XML content."""
        from lxml import etree

        namespaces = {
            'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
            'bp': 'http://www.biopax.org/release/biopax-level3.owl#'
        }

        root = etree.fromstring(biopax_content.encode())

        pathways = []
        for pathway in root.findall('.//bp:Pathway', namespaces):
            pw_data = {
                'id': pathway.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}ID'),
                'name': pathway.findtext('bp:displayName', namespaces=namespaces),
                'components': [
                    comp.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource')
                    for comp in pathway.findall('bp:pathwayComponent', namespaces)
                ]
            }
            pathways.append(pw_data)

        return {'pathways': pathways}

# Usage example
parser = ReactomeParser()
pathway = parser.get_pathway("R-HSA-109581")
print(f"Pathway: {pathway['displayName']}")
print(f"Events: {len(pathway.get('hasEvent', []))}")
```

---

## 2. KEGG Data Model

### 2.1 KGML (KEGG Markup Language) XML Schema

KGML is the XML format used to represent KEGG pathway maps with graphical and relational information.

#### Complete KGML Schema Definition

```xml
<!-- KGML DTD (Document Type Definition) -->
<!ELEMENT pathway (entry*, relation*, reaction*)>
<!ATTLIST pathway
    name        CDATA   #REQUIRED
    org         CDATA   #REQUIRED
    number      CDATA   #REQUIRED
    title       CDATA   #IMPLIED
    image       CDATA   #IMPLIED
    link        CDATA   #IMPLIED
>

<!ELEMENT entry (graphics*, component*)>
<!ATTLIST entry
    id          ID      #REQUIRED
    name        CDATA   #REQUIRED
    type        (ortholog|enzyme|gene|group|compound|map|brite|other) #REQUIRED
    reaction    CDATA   #IMPLIED
    link        CDATA   #IMPLIED
>

<!ELEMENT graphics EMPTY>
<!ATTLIST graphics
    name        CDATA   #IMPLIED
    x           CDATA   #REQUIRED
    y           CDATA   #REQUIRED
    coords      CDATA   #IMPLIED
    type        (rectangle|circle|roundrectangle|line) #REQUIRED
    width       CDATA   #IMPLIED
    height      CDATA   #IMPLIED
    fgcolor     CDATA   #IMPLIED
    bgcolor     CDATA   #IMPLIED
>

<!ELEMENT component EMPTY>
<!ATTLIST component
    id          IDREF   #REQUIRED
>

<!ELEMENT relation (subtype*)>
<!ATTLIST relation
    entry1      IDREF   #REQUIRED
    entry2      IDREF   #REQUIRED
    type        (ECrel|PPrel|GErel|PCrel|maplink) #REQUIRED
>

<!ELEMENT subtype EMPTY>
<!ATTLIST subtype
    name        CDATA   #REQUIRED
    value       CDATA   #REQUIRED
>

<!ELEMENT reaction (substrate*, product*)>
<!ATTLIST reaction
    id          IDREF   #REQUIRED
    name        CDATA   #REQUIRED
    type        (reversible|irreversible) #REQUIRED
>

<!ELEMENT substrate EMPTY>
<!ATTLIST substrate
    id          IDREF   #REQUIRED
    name        CDATA   #REQUIRED
>

<!ELEMENT product EMPTY>
<!ATTLIST product
    id          IDREF   #REQUIRED
    name        CDATA   #REQUIRED
>
```

### 2.2 Entry Types

| Type | Description | Name Format | Example |
|------|-------------|-------------|---------|
| `ortholog` | KO (KEGG Orthology) | `ko:K00001` | K number |
| `enzyme` | Enzyme | `ec:1.1.1.1` | EC number |
| `gene` | Gene product | `hsa:7157` | organism:gene_id |
| `group` | Complex/group | Entry IDs | Multiple components |
| `compound` | Chemical compound | `cpd:C00001` | C number |
| `map` | Linked pathway | `path:hsa00010` | Pathway reference |
| `brite` | BRITE hierarchy | `br:hsa00001` | BRITE ID |
| `other` | Other types | Various | Various |

### 2.3 Relation Types with Subtypes

#### ECrel (Enzyme-Enzyme Relation)

Represents metabolic relationships through shared compounds.

| Subtype | Value | Description |
|---------|-------|-------------|
| `compound` | Entry ID | Shared compound between enzymes |

#### PPrel (Protein-Protein Relation)

| Subtype | Value | Description |
|---------|-------|-------------|
| `activation` | `-->` | Activation |
| `inhibition` | `--\|` | Inhibition |
| `binding/association` | `---` | Binding/association |
| `dissociation` | `-+-` | Dissociation |
| `indirect effect` | `..>` | Indirect effect |
| `state change` | `...` | State change |
| `missing interaction` | `-/-` | Missing interaction |
| `phosphorylation` | `+p` | Phosphorylation |
| `dephosphorylation` | `-p` | Dephosphorylation |
| `glycosylation` | `+g` | Glycosylation |
| `ubiquitination` | `+u` | Ubiquitination |
| `methylation` | `+m` | Methylation |

#### GErel (Gene Expression Relation)

| Subtype | Value | Description |
|---------|-------|-------------|
| `expression` | `-->` | Gene expression activation |
| `repression` | `--\|` | Gene expression repression |

#### PCrel (Protein-Compound Relation)

| Subtype | Value | Description |
|---------|-------|-------------|
| `activation` | `-->` | Compound activates protein |
| `inhibition` | `--\|` | Compound inhibits protein |

### 2.4 Complete KGML Example

```xml
<?xml version="1.0"?>
<!DOCTYPE pathway SYSTEM "https://www.kegg.jp/kegg/xml/KGML_v0.7.2_.dtd">
<pathway name="path:hsa04115" org="hsa" number="04115"
         title="p53 signaling pathway"
         image="https://www.kegg.jp/kegg/pathway/hsa/hsa04115.png"
         link="https://www.kegg.jp/kegg-bin/show_pathway?hsa04115">

  <!-- Gene entries -->
  <entry id="1" name="hsa:7157" type="gene"
         link="https://www.kegg.jp/dbget-bin/www_bget?hsa:7157"
         reaction="">
    <graphics name="TP53, BCC7, LFS1, P53, TRP53" type="rectangle"
              x="547" y="308" width="46" height="17"
              fgcolor="#000000" bgcolor="#BFFFBF"/>
  </entry>

  <entry id="2" name="hsa:1026" type="gene"
         link="https://www.kegg.jp/dbget-bin/www_bget?hsa:1026">
    <graphics name="CDKN1A, CAP20, CDKN1, CIP1, MDA-6, P21, SDI1, WAF1, p21CIP1"
              type="rectangle" x="716" y="239" width="46" height="17"
              fgcolor="#000000" bgcolor="#BFFFBF"/>
  </entry>

  <entry id="3" name="hsa:472" type="gene"
         link="https://www.kegg.jp/dbget-bin/www_bget?hsa:472">
    <graphics name="ATM, AT1, ATA, ATC, ATD, ATDC, ATE, TEL1, TELO1"
              type="rectangle" x="321" y="308" width="46" height="17"
              fgcolor="#000000" bgcolor="#BFFFBF"/>
  </entry>

  <!-- Compound entry -->
  <entry id="50" name="cpd:C00027" type="compound"
         link="https://www.kegg.jp/dbget-bin/www_bget?C00027">
    <graphics name="C00027" type="circle"
              x="200" y="400" width="8" height="8"
              fgcolor="#000000" bgcolor="#FFFFFF"/>
  </entry>

  <!-- Pathway map reference -->
  <entry id="100" name="path:hsa04110" type="map"
         link="https://www.kegg.jp/dbget-bin/www_bget?hsa04110">
    <graphics name="Cell cycle" type="roundrectangle"
              x="716" y="171" width="100" height="34"
              fgcolor="#000000" bgcolor="#FFFFFF"/>
  </entry>

  <!-- Group entry (complex) -->
  <entry id="200" name="undefined" type="group">
    <graphics fgcolor="#000000" bgcolor="#FFFFFF"
              type="rectangle" x="400" y="350" width="60" height="40"/>
    <component id="1"/>
    <component id="3"/>
  </entry>

  <!-- Relations -->
  <relation entry1="3" entry2="1" type="PPrel">
    <subtype name="activation" value="-->"/>
    <subtype name="phosphorylation" value="+p"/>
  </relation>

  <relation entry1="1" entry2="2" type="GErel">
    <subtype name="expression" value="-->"/>
  </relation>

  <relation entry1="2" entry2="100" type="maplink">
    <subtype name="compound" value="4"/>
  </relation>

  <!-- Reactions (for metabolic pathways) -->
  <reaction id="300" name="rn:R00001" type="irreversible">
    <substrate id="50" name="cpd:C00027"/>
    <product id="51" name="cpd:C00001"/>
  </reaction>

</pathway>
```

### 2.5 KEGG Flat File Format

#### Gene Entry Format

```
ENTRY       7157              CDS       T01001
NAME        TP53, P53, BCC7, LFS1, TRP53
DEFINITION  tumor protein p53 (EC:2.7.11.1)
ORTHOLOGY   K04451  tumor protein p53 [EC:2.7.11.1]
ORGANISM    hsa  Homo sapiens (human)
PATHWAY     hsa04010  MAPK signaling pathway
            hsa04110  Cell cycle
            hsa04115  p53 signaling pathway
            hsa04151  PI3K-Akt signaling pathway
            hsa04210  Apoptosis
            hsa04218  Cellular senescence
            hsa04310  Wnt signaling pathway
            hsa05160  Hepatitis C
            hsa05200  Pathways in cancer
            hsa05202  Transcriptional misregulation in cancer
            hsa05206  MicroRNAs in cancer
            hsa05210  Colorectal cancer
            hsa05212  Pancreatic cancer
            hsa05213  Endometrial cancer
            hsa05215  Prostate cancer
            hsa05218  Melanoma
            hsa05219  Bladder cancer
            hsa05220  Chronic myeloid leukemia
            hsa05222  Small cell lung cancer
            hsa05223  Non-small cell lung cancer
            hsa05224  Breast cancer
            hsa05225  Hepatocellular carcinoma
            hsa05226  Gastric cancer
NETWORK     nt06110  MAPK signaling (virus)
            nt06160  Human T-cell leukemia virus 1 (HTLV-1)
            nt06162  Hepatitis B virus (HBV)
            nt06163  Hepatitis C virus (HCV)
            nt06164  Kaposi sarcoma-associated herpesvirus (KSHV)
            nt06165  Epstein-Barr virus (EBV)
            nt06166  Human papillomavirus (HPV)
            nt06167  Human cytomegalovirus (HCMV)
            nt06170  Influenza A virus (IAV)
            nt06210  ERK signaling
            nt06214  PI3K signaling
            nt06260  Colorectal cancer
            nt06261  Gastric cancer
            nt06262  Pancreatic cancer
            nt06263  Hepatocellular carcinoma
            nt06264  Renal cell carcinoma
            nt06265  Bladder cancer
            nt06266  Non-small cell lung cancer
            nt06267  Small cell lung cancer
            nt06268  Melanoma
            nt06269  Basal cell carcinoma
            nt06270  Breast cancer
            nt06271  Endometrial cancer
            nt06272  Prostate cancer
            nt06273  Glioma
            nt06274  Thyroid cancer
            nt06275  Acute myeloid leukemia
            nt06276  Chronic myeloid leukemia
MODULE      M00691  DNA damage-induced cell cycle checkpoints
DISEASE     H00004  Chronic myeloid leukemia
            H00005  Chronic lymphocytic leukemia
            H00006  Hairy-cell leukemia
            H00007  Hodgkin lymphoma
            H00008  Burkitt lymphoma
            H00010  Multiple myeloma
            H00013  Small cell lung cancer
            H00014  Non-small cell lung cancer
            H00015  Malignant pleural mesothelioma
            H00016  Oral cancer
            H00017  Esophageal cancer
            H00018  Gastric cancer
            H00019  Pancreatic cancer
            H00020  Colorectal cancer
            H00022  Bladder cancer
            H00025  Penile cancer
            H00026  Endometrial cancer
            H00027  Ovarian cancer
            H00028  Choriocarcinoma
            H00029  Vulvar cancer
            H00031  Breast cancer
            H00036  Osteosarcoma
            H00038  Malignant melanoma
            H00039  Basal cell carcinoma
            H00040  Squamous cell carcinoma
            H00041  Kaposi sarcoma
            H00042  Glioma
            H00044  Cancer of the adrenal cortex
            H00046  Cholangiocarcinoma
            H00047  Gallbladder cancer
            H00048  Hepatocellular carcinoma
            H01667  Medulloblastoma
            H01738  Noonan syndrome
            H02434  Li-Fraumeni syndrome
DRUG_TARGET D06272  Cenersen sodium (USAN)
            D12742  Rezatapopt (USAN/INN)
BRITE       KEGG Orthology (KO) [BR:hsa00001]
             09180 Brite Hierarchies
              09182 Protein families: genetic information processing
               03000 Transcription factors [BR:hsa03000]
                7157 (TP53)
            Protein kinases [BR:hsa01001]
             Serine/threonine kinases: Atypical group
              RIO family
               7157 (TP53)
            Transcription factors [BR:hsa03000]
             Eukaryotic type
              beta-Scaffold factors with minor groove contacts
               p53
                7157 (TP53)
POSITION    17:7661779..7687538
MOTIF       Pfam: P53 P53_tetramer P53_TAD
            PROSITE: P53_TAD
DBLINKS     NCBI-GeneID: 7157
            NCBI-ProteinID: NP_000537 NP_001119584 NP_001119585 NP_001119586
            OMIM: 191170
            HGNC: 11998
            Ensembl: ENSG00000141510
            Vega: OTTHUMG00000163020
            Pharos: P04637(Tchem)
            UniProt: P04637 A4GQV2 E9PR17 K7PPA8 Q53GA5
STRUCTURE   PDB: 1A1U 1AIE 1C26 1DT7 1GZH 1HS5 1JSP 1KZY 1MA3 1OLG 1OLH
                 1PES 1PET 1SAE 1SAF 1SAK 1SAL 1TSR 1TUP 1UOL 1XQH 1YC5
                 1YCQ 1YCR 1YCS 2AC0 2ADY 2AHI 2ATA 2B3G 2BIM 2BIN 2BIO
                 2BIP 2BIQ 2FOJ 2FOO 2FOS 2GE0 2GS0 2H1L 2H2D 2H2F 2H4F
                 2H4H 2H4J 2H59 2J0Z 2J10 2J11 2J1W 2J1X 2J1Y 2J1Z 2J20
                 2J21 2K8F 2L14 2LY4 2MEJ 2MWO 2MWY 2MZD 2N4T 2OCJ 2PCX
                 2Q1I 2VUK 2WGX 2X0U 2X0V 2X0W 2XWR 2YDR 2YGJ 3DAB 3DAC
                 3D06 3D07 3D08 3D09 3D0A 3IGL 3KMD 3KZ8 3LW1 3OFV 3PDH
                 3Q01 3Q05 3Q06 3TS8 3TU1 3VD1 3ZME 4AGL 4AGM 4AGN 4AGO
                 4AGP 4AGQ 4AIP 4BUZ 4BV2 4CZJ 4CZK 4CZL 4HJE 4IBB 4IBC
                 4IBD 4IBE 4IBF 4IBG 4IBQ 4IBR 4IBS 4IJT 4KVP 4LO9 4LOA
                 4LOE 4LOF 4MZI 4MZR 4QO1 4RP6 4RP7 4RP8 4WCE 4X34 4XR8
                 4ZYF 5A7B 5AB9 5ABA 5BN1 5BNA 5BUE 5G4N 5G4O 5HP0 5HPD
                 5HOU 5HOV 5LGY 5MCT 5MF7 5OLN 5OLO 5OLP 5OLQ 5OLD 5OLE
AASEQ       393
            MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
            DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
            SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
            RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
            SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
            PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
            GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
NTSEQ       1182
            atggaggagccgcagtcagatcctagcgtcgagccccctctgagtcaggaaacattttca
            gacctatggaaactacttcctgaaaacaacgttctgtcccccttgccgtcccaagcaatg
            gatgatttgatgctgtccccggacgatattgaacaatggttcactgaagacccaggtcca
            gatgaagctcccagaatgccagaggctgctccccccgtggcccctgcaccagcagctcct
            acaccggcggcccctgcaccagccccctcctggcccctgtcatcttctgtcccttcccag
            aaaacctaccagggcagctacggtttccgtctgggcttcttgcattctgggacagccaag
            tctgtgacttgcacgtactcccctgccctcaacaagatgttttgccaactggccaagacc
            tgccctgtgcagctgtgggttgattccacacccccgcccggcacccgcgtccgcgccatg
            gccatctacaagcagtcacagcacatgacggaggttgtgaggcgctgcccccaccatgag
            cgctgctcagatagcgatggtctggcccctcctcagcatcttatccgagtggaaggaaat
            ttgcgtgtggagtatttggatgacagaaacacttttcgacatagtgtggtggtgccctat
            gagccgcctgaggttggctctgactgtaccaccatccactacaactacatgtgtaacagt
            tcctgcatgggcggcatgaaccggaggcccatcctcaccatcatcacactggaagactcc
            agtggtaatctactgggacggaacagctttgaggtgcgtgtttgtgcctgtcctgggaga
            gaccggcgcacagaggaagagaatctccgcaagaaaggggagcctcaccacgagctgccc
            ccagggagcactaagcgagcactgcccaacaacaccagctcctctccccagccaaagaag
            aaaccactggatggagaatatttcacccttcagatccgtgggcgtgagcgcttcgagatg
            ttccgagagctgaatgaggccttggaactcaaggatgcccaggctgggaaggagccaggg
            gggagcagggctcactccagccacctgaagtccaaaaagggtcagtctacctcccgccat
            aaaaaactcatgttcaagacagaagggcctgactcagactga
///
```

#### Field Reference for Gene Entry

| Field | Description | Cardinality |
|-------|-------------|-------------|
| `ENTRY` | Entry ID, type, and organism code | Required |
| `NAME` | Gene names and symbols | Required |
| `DEFINITION` | Gene function description | Required |
| `ORTHOLOGY` | KO assignment | Optional |
| `ORGANISM` | Species code and name | Required |
| `PATHWAY` | Pathway memberships | Optional |
| `NETWORK` | Disease network memberships | Optional |
| `MODULE` | KEGG module memberships | Optional |
| `DISEASE` | Disease associations | Optional |
| `DRUG_TARGET` | Drugs targeting this gene | Optional |
| `BRITE` | BRITE hierarchy classifications | Optional |
| `POSITION` | Chromosomal position | Optional |
| `MOTIF` | Protein domain/motif annotations | Optional |
| `DBLINKS` | External database cross-references | Optional |
| `STRUCTURE` | PDB structure IDs | Optional |
| `AASEQ` | Amino acid sequence with length | Optional |
| `NTSEQ` | Nucleotide sequence with length | Optional |

### 2.6 Python Parsing Code

```python
import requests
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

@dataclass
class KGMLEntry:
    id: str
    name: str
    entry_type: str
    link: Optional[str] = None
    graphics_name: Optional[str] = None
    x: Optional[float] = None
    y: Optional[float] = None
    components: List[str] = None

@dataclass
class KGMLRelation:
    entry1: str
    entry2: str
    relation_type: str
    subtypes: List[Tuple[str, str]] = None

@dataclass
class KGMLReaction:
    id: str
    name: str
    reaction_type: str
    substrates: List[Tuple[str, str]] = None
    products: List[Tuple[str, str]] = None

class KEGGParser:
    """Parser for KEGG KGML and flat file formats."""

    BASE_URL = "https://rest.kegg.jp"

    def parse_kgml(self, kgml_content: str) -> Dict:
        """Parse KGML XML content."""
        root = ET.fromstring(kgml_content)

        pathway = {
            "name": root.get("name"),
            "org": root.get("org"),
            "number": root.get("number"),
            "title": root.get("title"),
            "entries": [],
            "relations": [],
            "reactions": []
        }

        # Parse entries
        for entry in root.findall("entry"):
            graphics = entry.find("graphics")
            components = [c.get("id") for c in entry.findall("component")]

            entry_data = KGMLEntry(
                id=entry.get("id"),
                name=entry.get("name"),
                entry_type=entry.get("type"),
                link=entry.get("link"),
                graphics_name=graphics.get("name") if graphics is not None else None,
                x=float(graphics.get("x")) if graphics is not None else None,
                y=float(graphics.get("y")) if graphics is not None else None,
                components=components if components else None
            )
            pathway["entries"].append(entry_data)

        # Parse relations
        for relation in root.findall("relation"):
            subtypes = [
                (st.get("name"), st.get("value"))
                for st in relation.findall("subtype")
            ]

            rel_data = KGMLRelation(
                entry1=relation.get("entry1"),
                entry2=relation.get("entry2"),
                relation_type=relation.get("type"),
                subtypes=subtypes
            )
            pathway["relations"].append(rel_data)

        # Parse reactions
        for reaction in root.findall("reaction"):
            substrates = [
                (s.get("id"), s.get("name"))
                for s in reaction.findall("substrate")
            ]
            products = [
                (p.get("id"), p.get("name"))
                for p in reaction.findall("product")
            ]

            rxn_data = KGMLReaction(
                id=reaction.get("id"),
                name=reaction.get("name"),
                reaction_type=reaction.get("type"),
                substrates=substrates,
                products=products
            )
            pathway["reactions"].append(rxn_data)

        return pathway

    def parse_flat_file(self, content: str) -> Dict:
        """Parse KEGG flat file format."""
        entry = {}
        current_field = None
        current_value = []

        for line in content.split('\n'):
            if line.startswith('///'):
                break

            # Check if new field
            if line and line[0] != ' ':
                # Save previous field
                if current_field:
                    if current_field in ['PATHWAY', 'DISEASE', 'DBLINKS',
                                         'DRUG_TARGET', 'STRUCTURE', 'BRITE']:
                        entry[current_field] = current_value
                    else:
                        entry[current_field] = '\n'.join(current_value)

                # Parse new field
                parts = line.split(None, 1)
                current_field = parts[0]
                current_value = [parts[1]] if len(parts) > 1 else []
            else:
                # Continuation of current field
                current_value.append(line.strip())

        # Handle final field
        if current_field:
            if current_field in ['PATHWAY', 'DISEASE', 'DBLINKS',
                                 'DRUG_TARGET', 'STRUCTURE', 'BRITE']:
                entry[current_field] = current_value
            else:
                entry[current_field] = '\n'.join(current_value)

        return entry

    def get_pathway_kgml(self, pathway_id: str) -> str:
        """Fetch KGML for a pathway."""
        url = f"{self.BASE_URL}/get/{pathway_id}/kgml"
        response = requests.get(url)
        return response.text

    def get_gene_entry(self, gene_id: str) -> Dict:
        """Fetch and parse gene entry."""
        url = f"{self.BASE_URL}/get/{gene_id}"
        response = requests.get(url)
        return self.parse_flat_file(response.text)

    def extract_gene_pathway_relations(self, kgml_data: Dict) -> List[Dict]:
        """Extract gene-to-gene relations from parsed KGML."""
        relations = []

        # Build entry lookup
        entry_map = {e.id: e for e in kgml_data["entries"]}

        for relation in kgml_data["relations"]:
            entry1 = entry_map.get(relation.entry1)
            entry2 = entry_map.get(relation.entry2)

            if entry1 and entry2 and entry1.entry_type == "gene" and entry2.entry_type == "gene":
                for subtype_name, subtype_value in (relation.subtypes or []):
                    relations.append({
                        "source": entry1.name,
                        "target": entry2.name,
                        "relation_type": relation.relation_type,
                        "interaction": subtype_name,
                        "symbol": subtype_value
                    })

        return relations

# Usage example
parser = KEGGParser()

# Parse KGML
kgml_content = parser.get_pathway_kgml("hsa04115")
pathway_data = parser.parse_kgml(kgml_content)
print(f"Pathway: {pathway_data['title']}")
print(f"Entries: {len(pathway_data['entries'])}")
print(f"Relations: {len(pathway_data['relations'])}")

# Parse gene entry
gene_data = parser.get_gene_entry("hsa:7157")
print(f"Gene: {gene_data.get('NAME')}")
print(f"Pathways: {len(gene_data.get('PATHWAY', []))}")
```

---

## 3. WikiPathways GPML Format

### 3.1 Complete XML Schema

GPML (GenMAPP Pathway Markup Language) is the XML format used by WikiPathways and PathVisio.

#### Root Pathway Element

```xml
<?xml version="1.0" encoding="UTF-8"?>
<Pathway xmlns="http://pathvisio.org/GPML/2013a"
         Name="Apoptosis"
         Version="20231115"
         Organism="Homo sapiens"
         Data-Source="WikiPathways"
         Last-Modified="20231115120000"
         License="CC BY 4.0">

  <!-- Optional header elements -->
  <Comment Source="WikiPathways-description">
    Apoptosis is a form of programmed cell death...
  </Comment>
  <Comment Source="WikiPathways-category">Cellular Process</Comment>

  <!-- BiopaxRef for literature citations -->
  <BiopaxRef>abc123</BiopaxRef>

  <!-- Graphics settings -->
  <Graphics BoardWidth="1200.0" BoardHeight="800.0"/>

  <!-- DataNodes (genes, proteins, metabolites) -->
  <!-- Interactions (edges) -->
  <!-- Labels (text annotations) -->
  <!-- Shapes (compartments) -->
  <!-- Groups (complexes) -->
  <!-- InfoBox -->
  <!-- Legend -->
  <!-- Biopax section -->

</Pathway>
```

### 3.2 DataNode Types and Xref Format

#### DataNode Element Structure

```xml
<DataNode TextLabel="CASP3" GraphId="a1b2c" Type="GeneProduct" GroupRef="group1">
  <Comment>Caspase-3, executioner caspase involved in apoptosis</Comment>
  <BiopaxRef>ref1</BiopaxRef>
  <Graphics CenterX="450.0" CenterY="300.0"
            Width="80.0" Height="20.0"
            ZOrder="32768"
            FontSize="10"
            Valign="Middle"
            ShapeType="Rectangle"
            Color="000000"
            FillColor="ffffff"/>
  <Xref Database="Ensembl" ID="ENSG00000164305"/>
</DataNode>
```

#### DataNode Types

| Type | Description | Typical Databases |
|------|-------------|-------------------|
| `GeneProduct` | Gene or its protein product | Ensembl, Entrez Gene, HGNC |
| `Protein` | Protein entity | UniProt, Ensembl |
| `Rna` | RNA molecule | Ensembl, miRBase |
| `Metabolite` | Small molecule | ChEBI, KEGG Compound, HMDB, PubChem |
| `Pathway` | Reference to another pathway | WikiPathways, Reactome |
| `Complex` | Protein complex | Complex Portal |
| `Unknown` | Unknown/unspecified type | - |

#### Xref Database Codes

| Database | Code | ID Format | Example |
|----------|------|-----------|---------|
| Ensembl | `Ensembl` or `En` | `ENSG[0-9]{11}` | `ENSG00000164305` |
| Entrez Gene | `Entrez Gene` or `L` | Numeric | `836` |
| HGNC | `HGNC` or `H` | `HGNC:[0-9]+` or numeric | `1504` |
| UniProt | `Uniprot-TrEMBL` or `S` | `[A-Z][0-9][A-Z0-9]{3}[0-9]` | `P42574` |
| ChEBI | `ChEBI` or `Ce` | `CHEBI:[0-9]+` | `CHEBI:15377` |
| KEGG Compound | `KEGG Compound` or `Ck` | `C[0-9]{5}` | `C00001` |
| HMDB | `HMDB` | `HMDB[0-9]+` | `HMDB0000001` |
| PubChem Compound | `PubChem-compound` | Numeric | `2244` |
| WikiPathways | `WikiPathways` or `Wp` | `WP[0-9]+` | `WP254` |
| Reactome | `Reactome` or `Re` | `R-[A-Z]{3}-[0-9]+` | `R-HSA-109581` |
| NCBI Protein | `RefSeq` | `[NX]P_[0-9]+` | `NP_004337` |
| miRBase | `miRBase mature sequence` | `MIMAT[0-9]+` | `MIMAT0000062` |
| Complex Portal | `Complex Portal` | `CPX-[0-9]+` | `CPX-2158` |

### 3.3 Interaction and GraphicalLine Elements

#### Interaction Element (Biological Relationship)

```xml
<Interaction GraphId="id12345">
  <Comment>Caspase-9 activates Caspase-3</Comment>
  <BiopaxRef>ref2</BiopaxRef>
  <Graphics ZOrder="12288" LineThickness="1.0" Color="000000">
    <Point X="350.0" Y="300.0" GraphRef="nodeA" RelX="1.0" RelY="0.0"/>
    <Point X="410.0" Y="300.0" GraphRef="nodeB" RelX="-1.0" RelY="0.0" ArrowHead="mim-stimulation"/>
  </Graphics>
  <Xref Database="" ID=""/>
</Interaction>
```

#### GraphicalLine Element (Visual Only)

```xml
<GraphicalLine GraphId="gline1">
  <Graphics ZOrder="12290" LineThickness="2.0" Color="0000ff" LineStyle="Dashed">
    <Point X="100.0" Y="100.0"/>
    <Point X="200.0" Y="100.0"/>
    <Point X="200.0" Y="200.0"/>
  </Graphics>
</GraphicalLine>
```

#### ArrowHead Types (MIM Notation)

| ArrowHead | Symbol | Meaning |
|-----------|--------|---------|
| `Arrow` | Simple arrow | General direction/flow |
| `mim-catalysis` | Circle on line | Catalysis |
| `mim-stimulation` | Open arrowhead | Stimulation/Activation |
| `mim-inhibition` | T-bar | Inhibition |
| `mim-conversion` | Filled arrowhead | Conversion |
| `mim-modification` | Filled arrowhead | Covalent modification |
| `mim-binding` | Two lines | Non-covalent binding |
| `mim-cleavage` | Diagonal line | Covalent cleavage |
| `mim-branching-left` | Branch left | Branching point |
| `mim-branching-right` | Branch right | Branching point |
| `mim-transcription-translation` | Angled arrow | Transcription/Translation |
| `mim-gap` | Gap in line | Continuation |
| `TBar` | T-bar | Inhibition (alias) |

#### Line Styles

| Style | Description |
|-------|-------------|
| `Solid` | Continuous line (default) |
| `Dashed` | Dashed line |
| `Double` | Double line |

#### Anchor Points

```xml
<Interaction GraphId="main_interaction">
  <Graphics ZOrder="12288" LineThickness="1.0">
    <Point X="100.0" Y="200.0" GraphRef="nodeA" RelX="1.0" RelY="0.0"/>
    <Point X="300.0" Y="200.0" GraphRef="nodeB" RelX="-1.0" RelY="0.0" ArrowHead="Arrow"/>
    <Anchor Position="0.5" Shape="None" GraphId="anchor1"/>
  </Graphics>
</Interaction>

<!-- Another interaction connecting to the anchor -->
<Interaction GraphId="branch_interaction">
  <Graphics ZOrder="12289" LineThickness="1.0">
    <Point X="200.0" Y="100.0" GraphRef="nodeC" RelX="0.0" RelY="1.0"/>
    <Point X="200.0" Y="200.0" GraphRef="anchor1" RelX="0.0" RelY="0.0" ArrowHead="mim-catalysis"/>
  </Graphics>
</Interaction>
```

### 3.4 State Elements for Modifications

State elements represent post-translational modifications or other states on DataNodes.

```xml
<DataNode TextLabel="TP53" GraphId="tp53node" Type="GeneProduct">
  <Graphics CenterX="300.0" CenterY="200.0" Width="80.0" Height="20.0" ZOrder="32768"/>
  <Xref Database="Ensembl" ID="ENSG00000141510"/>
</DataNode>

<!-- Phosphorylation state -->
<State TextLabel="P" GraphId="state1" GraphRef="tp53node">
  <Comment>Phosphorylation at Ser15</Comment>
  <Graphics RelX="1.0" RelY="-1.0" Width="15.0" Height="15.0"
            ShapeType="Oval" FillColor="ff0000"/>
  <Xref Database="PSI-MOD" ID="MOD:00046"/>
</State>

<!-- Ubiquitination state -->
<State TextLabel="Ub" GraphId="state2" GraphRef="tp53node">
  <Comment>Ubiquitination</Comment>
  <Graphics RelX="1.0" RelY="1.0" Width="15.0" Height="15.0"
            ShapeType="Oval" FillColor="0000ff"/>
  <Xref Database="PSI-MOD" ID="MOD:01148"/>
</State>

<!-- Acetylation state -->
<State TextLabel="Ac" GraphId="state3" GraphRef="tp53node">
  <Comment>Acetylation at Lys382</Comment>
  <Graphics RelX="-1.0" RelY="-1.0" Width="15.0" Height="15.0"
            ShapeType="Oval" FillColor="00ff00"/>
  <Xref Database="PSI-MOD" ID="MOD:00394"/>
</State>
```

#### State ShapeTypes

| ShapeType | Description |
|-----------|-------------|
| `Oval` | Circle/oval (most common for PTMs) |
| `Rectangle` | Rectangular state |
| `RoundedRectangle` | Rounded rectangle |

### 3.5 Group Element (Complexes)

```xml
<!-- Group definition -->
<Group GroupId="apoptosome" Style="Complex" TextLabel="Apoptosome">
  <Comment>Apoptosome complex formed during apoptosis</Comment>
  <BiopaxRef>ref3</BiopaxRef>
</Group>

<!-- Group members -->
<DataNode TextLabel="APAF1" GraphId="apaf1" Type="GeneProduct" GroupRef="apoptosome">
  <Graphics CenterX="400.0" CenterY="400.0" Width="60.0" Height="20.0" ZOrder="32769"/>
  <Xref Database="Ensembl" ID="ENSG00000120868"/>
</DataNode>

<DataNode TextLabel="CYCS" GraphId="cycs" Type="GeneProduct" GroupRef="apoptosome">
  <Graphics CenterX="400.0" CenterY="420.0" Width="60.0" Height="20.0" ZOrder="32770"/>
  <Xref Database="Ensembl" ID="ENSG00000172115"/>
</DataNode>

<DataNode TextLabel="CASP9" GraphId="casp9" Type="GeneProduct" GroupRef="apoptosome">
  <Graphics CenterX="400.0" CenterY="440.0" Width="60.0" Height="20.0" ZOrder="32771"/>
  <Xref Database="Ensembl" ID="ENSG00000132906"/>
</DataNode>
```

#### Group Styles

| Style | Description |
|-------|-------------|
| `None` | No visual grouping |
| `Group` | Simple grouping |
| `Complex` | Protein complex |
| `Pathway` | Nested pathway |

### 3.6 Label and Shape Elements

#### Label Element

```xml
<Label TextLabel="Mitochondria" GraphId="label1">
  <Graphics CenterX="200.0" CenterY="100.0"
            Width="100.0" Height="25.0"
            ZOrder="28672"
            FontSize="12"
            FontWeight="Bold"
            FontStyle="Italic"
            Color="666666"/>
</Label>
```

#### Shape Element (Compartments)

```xml
<Shape GraphId="mitochondria_shape" TextLabel="Mitochondria">
  <Comment>Mitochondrial compartment</Comment>
  <Graphics CenterX="300.0" CenterY="400.0"
            Width="200.0" Height="300.0"
            ZOrder="16384"
            FontSize="12"
            Color="cccccc"
            FillColor="ffffcc"
            ShapeType="RoundedRectangle"
            Rotation="0.0"/>
  <Xref Database="Gene Ontology" ID="GO:0005739"/>
</Shape>
```

#### ShapeTypes

| ShapeType | Description |
|-----------|-------------|
| `Rectangle` | Rectangle |
| `RoundedRectangle` | Rectangle with rounded corners |
| `Oval` | Ellipse/Circle |
| `Triangle` | Triangle |
| `Pentagon` | Pentagon |
| `Hexagon` | Hexagon |
| `Octagon` | Octagon |
| `Arc` | Arc shape |
| `Brace` | Brace shape |
| `Mitochondria` | Mitochondria shape |
| `Sarcoplasmic Reticulum` | SR shape |
| `Endoplasmic Reticulum` | ER shape |
| `Golgi Apparatus` | Golgi shape |
| `Cell` | Cell outline |
| `Nucleus` | Nucleus shape |
| `Organelle` | Generic organelle |
| `Vesicle` | Vesicle shape |

### 3.7 Biopax Section for Citations

```xml
<Biopax>
  <bp:PublicationXref xmlns:bp="http://www.biopax.org/release/biopax-level3.owl#"
                      rdf:id="ref1" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
    <bp:ID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">12345678</bp:ID>
    <bp:DB rdf:datatype="http://www.w3.org/2001/XMLSchema#string">PubMed</bp:DB>
    <bp:TITLE rdf:datatype="http://www.w3.org/2001/XMLSchema#string">
      Caspase activation in apoptosis
    </bp:TITLE>
    <bp:SOURCE rdf:datatype="http://www.w3.org/2001/XMLSchema#string">
      Cell Death Differ.
    </bp:SOURCE>
    <bp:YEAR rdf:datatype="http://www.w3.org/2001/XMLSchema#string">2023</bp:YEAR>
    <bp:AUTHORS rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Smith J</bp:AUTHORS>
    <bp:AUTHORS rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Jones A</bp:AUTHORS>
  </bp:PublicationXref>

  <bp:openControlledVocabulary xmlns:bp="http://www.biopax.org/release/biopax-level3.owl#">
    <bp:TERM xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
             rdf:datatype="http://www.w3.org/2001/XMLSchema#string">apoptotic process</bp:TERM>
    <bp:ID xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
           rdf:datatype="http://www.w3.org/2001/XMLSchema#string">GO:0006915</bp:ID>
    <bp:Ontology xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
                 rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Gene Ontology</bp:Ontology>
  </bp:openControlledVocabulary>
</Biopax>
```

### 3.8 Python Parsing Code

```python
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional
from dataclasses import dataclass, field

@dataclass
class GPMLDataNode:
    graph_id: str
    text_label: str
    node_type: str
    database: Optional[str] = None
    identifier: Optional[str] = None
    center_x: float = 0.0
    center_y: float = 0.0
    width: float = 80.0
    height: float = 20.0
    group_ref: Optional[str] = None
    states: List['GPMLState'] = field(default_factory=list)

@dataclass
class GPMLState:
    graph_id: str
    text_label: str
    parent_ref: str
    database: Optional[str] = None
    identifier: Optional[str] = None

@dataclass
class GPMLInteraction:
    graph_id: str
    source_ref: str
    target_ref: str
    arrow_head: Optional[str] = None
    line_style: str = "Solid"

@dataclass
class GPMLGroup:
    group_id: str
    style: str
    text_label: Optional[str] = None
    members: List[str] = field(default_factory=list)

class GPMLParser:
    """Parser for WikiPathways GPML format."""

    NAMESPACE = {"gpml": "http://pathvisio.org/GPML/2013a"}

    def parse(self, gpml_content: str) -> Dict:
        """Parse GPML XML content."""
        root = ET.fromstring(gpml_content)

        # Handle namespace
        ns = self.NAMESPACE
        if not root.tag.startswith('{'):
            ns = {}

        pathway = {
            "name": root.get("Name"),
            "organism": root.get("Organism"),
            "version": root.get("Version"),
            "data_nodes": [],
            "interactions": [],
            "groups": [],
            "labels": [],
            "shapes": []
        }

        # Parse DataNodes
        for node in root.findall(".//gpml:DataNode", ns) or root.findall(".//DataNode"):
            xref = node.find("gpml:Xref", ns) or node.find("Xref")
            graphics = node.find("gpml:Graphics", ns) or node.find("Graphics")

            data_node = GPMLDataNode(
                graph_id=node.get("GraphId"),
                text_label=node.get("TextLabel"),
                node_type=node.get("Type", "Unknown"),
                database=xref.get("Database") if xref is not None else None,
                identifier=xref.get("ID") if xref is not None else None,
                center_x=float(graphics.get("CenterX", 0)) if graphics is not None else 0,
                center_y=float(graphics.get("CenterY", 0)) if graphics is not None else 0,
                width=float(graphics.get("Width", 80)) if graphics is not None else 80,
                height=float(graphics.get("Height", 20)) if graphics is not None else 20,
                group_ref=node.get("GroupRef")
            )
            pathway["data_nodes"].append(data_node)

        # Parse States
        for state in root.findall(".//gpml:State", ns) or root.findall(".//State"):
            xref = state.find("gpml:Xref", ns) or state.find("Xref")

            state_data = GPMLState(
                graph_id=state.get("GraphId"),
                text_label=state.get("TextLabel"),
                parent_ref=state.get("GraphRef"),
                database=xref.get("Database") if xref is not None else None,
                identifier=xref.get("ID") if xref is not None else None
            )

            # Find parent and add state
            for node in pathway["data_nodes"]:
                if node.graph_id == state_data.parent_ref:
                    node.states.append(state_data)
                    break

        # Parse Interactions
        for interaction in root.findall(".//gpml:Interaction", ns) or root.findall(".//Interaction"):
            graphics = interaction.find("gpml:Graphics", ns) or interaction.find("Graphics")

            if graphics is not None:
                points = graphics.findall("gpml:Point", ns) or graphics.findall("Point")

                if len(points) >= 2:
                    source_point = points[0]
                    target_point = points[-1]

                    interaction_data = GPMLInteraction(
                        graph_id=interaction.get("GraphId"),
                        source_ref=source_point.get("GraphRef"),
                        target_ref=target_point.get("GraphRef"),
                        arrow_head=target_point.get("ArrowHead"),
                        line_style=graphics.get("LineStyle", "Solid")
                    )
                    pathway["interactions"].append(interaction_data)

        # Parse Groups
        for group in root.findall(".//gpml:Group", ns) or root.findall(".//Group"):
            group_data = GPMLGroup(
                group_id=group.get("GroupId"),
                style=group.get("Style", "None"),
                text_label=group.get("TextLabel")
            )

            # Find members
            for node in pathway["data_nodes"]:
                if node.group_ref == group_data.group_id:
                    group_data.members.append(node.graph_id)

            pathway["groups"].append(group_data)

        return pathway

    def extract_gene_list(self, pathway_data: Dict) -> List[Dict]:
        """Extract gene identifiers from parsed GPML."""
        genes = []

        for node in pathway_data["data_nodes"]:
            if node.node_type in ["GeneProduct", "Protein"]:
                genes.append({
                    "label": node.text_label,
                    "database": node.database,
                    "identifier": node.identifier,
                    "modifications": [
                        {"label": s.text_label, "db": s.database, "id": s.identifier}
                        for s in node.states
                    ]
                })

        return genes

    def extract_interactions(self, pathway_data: Dict) -> List[Dict]:
        """Extract interactions as source-target pairs."""
        # Build node lookup
        node_map = {n.graph_id: n for n in pathway_data["data_nodes"]}

        interactions = []
        for inter in pathway_data["interactions"]:
            source = node_map.get(inter.source_ref)
            target = node_map.get(inter.target_ref)

            if source and target:
                interactions.append({
                    "source_label": source.text_label,
                    "source_id": source.identifier,
                    "target_label": target.text_label,
                    "target_id": target.identifier,
                    "interaction_type": inter.arrow_head or "unknown"
                })

        return interactions

# Usage example
import requests

def fetch_wikipathways_gpml(pathway_id: str) -> str:
    """Fetch GPML from WikiPathways."""
    url = f"https://webservice.wikipathways.org/getPathway?pwId={pathway_id}&format=json"
    response = requests.get(url)
    data = response.json()
    return data.get("pathway", {}).get("gpml", "")

# Parse a pathway
gpml_content = fetch_wikipathways_gpml("WP254")
parser = GPMLParser()
pathway = parser.parse(gpml_content)

print(f"Pathway: {pathway['name']}")
print(f"Organism: {pathway['organism']}")
print(f"DataNodes: {len(pathway['data_nodes'])}")
print(f"Interactions: {len(pathway['interactions'])}")
print(f"Groups: {len(pathway['groups'])}")

# Extract genes
genes = parser.extract_gene_list(pathway)
for gene in genes[:5]:
    print(f"  {gene['label']}: {gene['database']}:{gene['identifier']}")
```

---

## 4. Gene Ontology Structure

### 4.1 OBO Format Specification

The OBO (Open Biomedical Ontologies) format is a human-readable ontology format used by the Gene Ontology.

#### OBO File Structure

```obo
format-version: 1.4
data-version: releases/2024-01-17
subsetdef: goslim_generic "Generic GO slim"
default-namespace: gene_ontology
ontology: go

[Term]
id: GO:0006915
name: apoptotic process
namespace: biological_process
alt_id: GO:0006917
def: "A programmed cell death process which begins when a cell receives an internal or external signal..." [GOC:curators, PMID:18846107]
comment: This term should be used to annotate gene products involved in apoptosis.
subset: goslim_generic
synonym: "apoptosis" EXACT []
synonym: "programmed cell death by apoptosis" EXACT []
synonym: "type I programmed cell death" NARROW []
xref: Reactome:R-HSA-109581 "Apoptosis"
xref: Wikipedia:Apoptosis
is_a: GO:0012501 ! programmed cell death
relationship: part_of GO:0043067 ! regulation of programmed cell death

[Typedef]
id: part_of
name: part of
is_transitive: true
```

### 4.2 GO Term Structure Fields

| Field | Description | Cardinality | Example |
|-------|-------------|-------------|---------|
| `id` | GO identifier | Required, 1 | `GO:0006915` |
| `name` | Term name | Required, 1 | `apoptotic process` |
| `namespace` | Ontology branch (P/F/C) | Required, 1 | `biological_process` |
| `def` | Definition with dbxrefs | Optional, 1 | `"A programmed..." [GOC:curators]` |
| `alt_id` | Alternative/merged IDs | Optional, * | `GO:0006917` |
| `synonym` | Term synonyms | Optional, * | `"apoptosis" EXACT []` |
| `xref` | Cross-references | Optional, * | `Reactome:R-HSA-109581` |
| `is_a` | Parent term | Optional, * | `GO:0012501` |
| `relationship` | Other relationships | Optional, * | `part_of GO:0043067` |
| `is_obsolete` | Obsolete flag | Optional, 1 | `true` |

### 4.3 GO Namespaces

| Namespace | Code | Root Term |
|-----------|------|-----------|
| `biological_process` | P | GO:0008150 |
| `molecular_function` | F | GO:0003674 |
| `cellular_component` | C | GO:0005575 |

### 4.4 GAF (Gene Association File) 2.2 Format

| Col | Name | Required | Example |
|-----|------|----------|---------|
| 1 | DB | Yes | `UniProtKB` |
| 2 | DB Object ID | Yes | `P04637` |
| 3 | DB Object Symbol | Yes | `TP53` |
| 4 | Qualifier | Yes | `enables` |
| 5 | GO ID | Yes | `GO:0006915` |
| 6 | DB:Reference | Yes | `PMID:12505355` |
| 7 | Evidence Code | Yes | `IDA` |
| 8 | With/From | Optional | `UniProtKB:P04637` |
| 9 | Aspect | Yes | `P` |
| 10 | DB Object Name | Optional | `Cellular tumor antigen p53` |
| 11 | DB Object Synonym | Optional | `P53\|TRP53` |
| 12 | DB Object Type | Yes | `protein` |
| 13 | Taxon | Yes | `taxon:9606` |
| 14 | Date | Yes | `20240115` |
| 15 | Assigned By | Yes | `UniProt` |
| 16 | Annotation Extension | Optional | `occurs_in(CL:0000066)` |
| 17 | Gene Product Form ID | Optional | `UniProtKB:P04637-1` |

#### GAF Example

```tsv
!gaf-version: 2.2
UniProtKB	P04637	TP53	involved_in	GO:0006915	PMID:12505355	IDA		P	Cellular tumor antigen p53	P53|TRP53	protein	taxon:9606	20231115	UniProt
```

### 4.5 Evidence Codes

| Code | Category | Description |
|------|----------|-------------|
| `IDA` | Experimental | Inferred from Direct Assay |
| `IPI` | Experimental | Inferred from Physical Interaction |
| `IMP` | Experimental | Inferred from Mutant Phenotype |
| `IGI` | Experimental | Inferred from Genetic Interaction |
| `IEP` | Experimental | Inferred from Expression Pattern |
| `ISS` | Computational | Inferred from Sequence Similarity |
| `ISO` | Computational | Inferred from Sequence Orthology |
| `TAS` | Author | Traceable Author Statement |
| `IEA` | Electronic | Inferred from Electronic Annotation |

### 4.6 GPAD/GPI Formats

GPAD (Gene Product Association Data) and GPI (Gene Product Information) separate annotation from gene product metadata.

#### GPAD 2.0 Example

```tsv
!gpad-version: 2.0
UniProtKB	P04637	RO:0002331	GO:0006915	PMID:12505355	ECO:0000314		20240115	UniProt
```

#### GPI 2.0 Example

```tsv
!gpi-version: 2.0
UniProtKB	P04637	TP53	Cellular tumor antigen p53	P53|TRP53	PR:000000001	NCBITaxon:9606		Ensembl:ENSG00000141510
```

### 4.7 Python Parsing Code

```python
from typing import Dict, List, Set
from dataclasses import dataclass, field

@dataclass
class GOTerm:
    id: str
    name: str
    namespace: str
    definition: str = ""
    is_a: List[str] = field(default_factory=list)
    relationships: List[Dict] = field(default_factory=list)
    is_obsolete: bool = False

class OBOParser:
    """Parser for OBO format ontology files."""

    def parse(self, obo_content: str) -> Dict[str, GOTerm]:
        terms = {}
        current_term = None

        for line in obo_content.split('\n'):
            line = line.strip()

            if line == '[Term]':
                if current_term and current_term.id:
                    terms[current_term.id] = current_term
                current_term = GOTerm(id="", name="", namespace="")
                continue

            if not current_term:
                continue

            if ': ' in line:
                key, value = line.split(': ', 1)
                if key == 'id':
                    current_term.id = value
                elif key == 'name':
                    current_term.name = value
                elif key == 'namespace':
                    current_term.namespace = value
                elif key == 'is_a':
                    parent_id = value.split(' ! ')[0]
                    current_term.is_a.append(parent_id)

        if current_term and current_term.id:
            terms[current_term.id] = current_term

        return terms

    def get_ancestors(self, terms: Dict[str, GOTerm], term_id: str) -> Set[str]:
        ancestors = set()
        to_visit = [term_id]

        while to_visit:
            current_id = to_visit.pop()
            term = terms.get(current_id)
            if term:
                for parent_id in term.is_a:
                    if parent_id not in ancestors:
                        ancestors.add(parent_id)
                        to_visit.append(parent_id)

        return ancestors

# Usage
parser = OBOParser()
with open('go.obo', 'r') as f:
    terms = parser.parse(f.read())

term = terms.get('GO:0006915')
print(f"Term: {term.name}, Parents: {term.is_a}")
```

---

## 5. BioPAX (Biological Pathway Exchange)

### 5.1 Level 3 Specification Overview

BioPAX is an OWL-based ontology for biological pathway data exchange.

### 5.2 Core Class Hierarchy

```
Entity
├── PhysicalEntity
│   ├── Protein
│   ├── SmallMolecule
│   ├── Complex
│   ├── Dna / DnaRegion
│   └── Rna / RnaRegion
├── Gene
├── Pathway
└── Interaction
    ├── Conversion
    │   ├── BiochemicalReaction
    │   ├── Transport
    │   ├── ComplexAssembly
    │   └── Degradation
    ├── Control
    │   ├── Catalysis
    │   └── Modulation
    ├── MolecularInteraction
    └── TemplateReaction
```

### 5.3 Core Classes

#### Pathway

```xml
<bp:Pathway rdf:ID="Pathway_Apoptosis">
  <bp:displayName>Apoptosis</bp:displayName>
  <bp:organism rdf:resource="#BioSource_Homo_sapiens"/>
  <bp:pathwayComponent rdf:resource="#BiochemicalReaction_1"/>
  <bp:pathwayComponent rdf:resource="#Catalysis_1"/>
  <bp:xref rdf:resource="#UnificationXref_GO_0006915"/>
</bp:Pathway>
```

#### BiochemicalReaction

```xml
<bp:BiochemicalReaction rdf:ID="BiochemicalReaction_1">
  <bp:displayName>Cytochrome c release</bp:displayName>
  <bp:left rdf:resource="#Protein_CytC_mito"/>
  <bp:right rdf:resource="#Protein_CytC_cyto"/>
  <bp:conversionDirection>LEFT-TO-RIGHT</bp:conversionDirection>
</bp:BiochemicalReaction>
```

#### Protein

```xml
<bp:Protein rdf:ID="Protein_TP53">
  <bp:displayName>TP53</bp:displayName>
  <bp:entityReference rdf:resource="#ProteinReference_P04637"/>
  <bp:cellularLocation rdf:resource="#CellularLocationVocabulary_nucleus"/>
  <bp:feature rdf:resource="#ModificationFeature_phospho_Ser15"/>
</bp:Protein>
```

#### Catalysis

```xml
<bp:Catalysis rdf:ID="Catalysis_1">
  <bp:controller rdf:resource="#Protein_CASP9"/>
  <bp:controlled rdf:resource="#BiochemicalReaction_1"/>
  <bp:controlType>ACTIVATION</bp:controlType>
</bp:Catalysis>
```

### 5.4 Xref Types

| Type | Description |
|------|-------------|
| `UnificationXref` | Same entity in another database |
| `RelationshipXref` | Related entity |
| `PublicationXref` | Literature reference |

```xml
<bp:UnificationXref rdf:ID="UnificationXref_P04637">
  <bp:db>UniProt</bp:db>
  <bp:id>P04637</bp:id>
</bp:UnificationXref>
```

### 5.5 Python Parsing Code

```python
from rdflib import Graph, Namespace
from rdflib.namespace import RDF

BP = Namespace("http://www.biopax.org/release/biopax-level3.owl#")

class BioPAXParser:
    def __init__(self):
        self.graph = Graph()

    def parse_file(self, filepath: str):
        self.graph.parse(filepath, format="xml")

    def get_pathways(self):
        pathways = []
        for pw_uri in self.graph.subjects(RDF.type, BP.Pathway):
            name = str(self.graph.value(pw_uri, BP.displayName))
            pathways.append({'uri': str(pw_uri), 'name': name})
        return pathways

    def get_proteins(self):
        proteins = []
        for p_uri in self.graph.subjects(RDF.type, BP.Protein):
            name = str(self.graph.value(p_uri, BP.displayName))
            ref = self.graph.value(p_uri, BP.entityReference)
            proteins.append({'uri': str(p_uri), 'name': name, 'ref': str(ref) if ref else None})
        return proteins

# Usage
parser = BioPAXParser()
parser.parse_file("pathway.owl")
print(parser.get_pathways())
```

---

## 6. SBML (Systems Biology Markup Language)

### 6.1 Model Structure Overview

SBML is an XML format for computational models of biological systems, used by BioModels and other repositories.

### 6.2 Core Elements

```xml
<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="apoptosis_model" name="Apoptosis Model" timeUnits="second">

    <!-- Unit definitions -->
    <listOfUnitDefinitions>
      <unitDefinition id="per_second">
        <listOfUnits>
          <unit kind="second" exponent="-1"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>

    <!-- Compartments -->
    <listOfCompartments>
      <compartment id="cytosol" name="Cytosol" spatialDimensions="3"
                   size="1" units="litre" constant="true"/>
      <compartment id="mitochondria" name="Mitochondria" spatialDimensions="3"
                   size="0.1" units="litre" constant="true"/>
    </listOfCompartments>

    <!-- Species -->
    <listOfSpecies>
      <species id="CASP3" name="Caspase-3" compartment="cytosol"
               initialConcentration="100" substanceUnits="mole"
               hasOnlySubstanceUnits="false" boundaryCondition="false"
               constant="false">
        <annotation>
          <rdf:RDF>
            <rdf:Description rdf:about="#CASP3">
              <bqbiol:is>
                <rdf:Bag>
                  <rdf:li rdf:resource="http://identifiers.org/uniprot/P42574"/>
                </rdf:Bag>
              </bqbiol:is>
            </rdf:Description>
          </rdf:RDF>
        </annotation>
      </species>
      <species id="CASP9" name="Caspase-9" compartment="cytosol"
               initialConcentration="50" substanceUnits="mole"
               hasOnlySubstanceUnits="false" boundaryCondition="false"
               constant="false"/>
      <species id="CYCS" name="Cytochrome c" compartment="mitochondria"
               initialConcentration="200" substanceUnits="mole"
               hasOnlySubstanceUnits="false" boundaryCondition="false"
               constant="false"/>
    </listOfSpecies>

    <!-- Parameters -->
    <listOfParameters>
      <parameter id="k1" name="Activation rate" value="0.1" units="per_second" constant="true"/>
      <parameter id="k2" name="Degradation rate" value="0.01" units="per_second" constant="true"/>
      <parameter id="Km" name="Michaelis constant" value="10" units="mole" constant="true"/>
    </listOfParameters>

    <!-- Reactions -->
    <listOfReactions>
      <reaction id="CASP3_activation" name="Caspase-3 activation" reversible="false">
        <listOfReactants>
          <speciesReference species="CASP3" stoichiometry="1" constant="true"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="CASP3_active" stoichiometry="1" constant="true"/>
        </listOfProducts>
        <listOfModifiers>
          <modifierSpeciesReference species="CASP9"/>
        </listOfModifiers>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci>k1</ci>
              <ci>CASP9</ci>
              <ci>CASP3</ci>
            </apply>
          </math>
        </kineticLaw>
        <annotation>
          <rdf:RDF>
            <rdf:Description rdf:about="#CASP3_activation">
              <bqbiol:isVersionOf>
                <rdf:Bag>
                  <rdf:li rdf:resource="http://identifiers.org/go/GO:0006919"/>
                </rdf:Bag>
              </bqbiol:isVersionOf>
            </rdf:Description>
          </rdf:RDF>
        </annotation>
      </reaction>
    </listOfReactions>

  </model>
</sbml>
```

### 6.3 Element Reference

#### Compartment

| Attribute | Type | Description |
|-----------|------|-------------|
| `id` | SId | Unique identifier |
| `name` | String | Human-readable name |
| `spatialDimensions` | Double | 0, 1, 2, or 3 |
| `size` | Double | Volume/area/length |
| `units` | UnitSIdRef | Unit of size |
| `constant` | Boolean | Size is constant |

#### Species

| Attribute | Type | Description |
|-----------|------|-------------|
| `id` | SId | Unique identifier |
| `name` | String | Human-readable name |
| `compartment` | SIdRef | Containing compartment |
| `initialConcentration` | Double | Initial concentration |
| `initialAmount` | Double | Initial amount |
| `substanceUnits` | UnitSIdRef | Substance units |
| `hasOnlySubstanceUnits` | Boolean | Interpret as amount only |
| `boundaryCondition` | Boolean | Fixed boundary |
| `constant` | Boolean | Value is constant |

#### Reaction

| Attribute | Type | Description |
|-----------|------|-------------|
| `id` | SId | Unique identifier |
| `name` | String | Human-readable name |
| `reversible` | Boolean | Reversible reaction |
| `fast` | Boolean | Fast reaction (deprecated) |

#### SpeciesReference

| Attribute | Type | Description |
|-----------|------|-------------|
| `species` | SIdRef | Reference to species |
| `stoichiometry` | Double | Stoichiometric coefficient |
| `constant` | Boolean | Stoichiometry is constant |

### 6.4 MIRIAM Annotations

SBML uses MIRIAM (Minimum Information Requested in Annotation of Models) qualifiers:

| Qualifier | Description |
|-----------|-------------|
| `bqbiol:is` | Identity relationship |
| `bqbiol:isVersionOf` | More general term |
| `bqbiol:hasVersion` | More specific term |
| `bqbiol:isHomologTo` | Homolog |
| `bqbiol:isEncodedBy` | Encoded by gene |
| `bqbiol:encodes` | Encodes protein |
| `bqbiol:occursIn` | Cellular location |
| `bqmodel:is` | Model identity |
| `bqmodel:isDerivedFrom` | Derived from model |

### 6.5 Python Parsing Code

```python
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import List, Dict, Optional

SBML_NS = "{http://www.sbml.org/sbml/level3/version2/core}"

@dataclass
class SBMLSpecies:
    id: str
    name: str
    compartment: str
    initial_concentration: float = 0.0
    annotations: List[str] = field(default_factory=list)

@dataclass
class SBMLReaction:
    id: str
    name: str
    reversible: bool
    reactants: List[Dict] = field(default_factory=list)
    products: List[Dict] = field(default_factory=list)
    modifiers: List[str] = field(default_factory=list)
    kinetic_law: Optional[str] = None

class SBMLParser:
    """Parser for SBML Level 3 models."""

    def parse(self, sbml_content: str) -> Dict:
        root = ET.fromstring(sbml_content)
        model = root.find(f"{SBML_NS}model")

        result = {
            "id": model.get("id"),
            "name": model.get("name"),
            "compartments": [],
            "species": [],
            "parameters": [],
            "reactions": []
        }

        # Parse compartments
        compartments = model.find(f"{SBML_NS}listOfCompartments")
        if compartments:
            for comp in compartments.findall(f"{SBML_NS}compartment"):
                result["compartments"].append({
                    "id": comp.get("id"),
                    "name": comp.get("name"),
                    "size": float(comp.get("size", 1))
                })

        # Parse species
        species_list = model.find(f"{SBML_NS}listOfSpecies")
        if species_list:
            for species in species_list.findall(f"{SBML_NS}species"):
                sp = SBMLSpecies(
                    id=species.get("id"),
                    name=species.get("name", ""),
                    compartment=species.get("compartment"),
                    initial_concentration=float(species.get("initialConcentration", 0))
                )
                result["species"].append(sp)

        # Parse reactions
        reactions = model.find(f"{SBML_NS}listOfReactions")
        if reactions:
            for rxn in reactions.findall(f"{SBML_NS}reaction"):
                reaction = SBMLReaction(
                    id=rxn.get("id"),
                    name=rxn.get("name", ""),
                    reversible=rxn.get("reversible", "true") == "true"
                )

                # Reactants
                reactants = rxn.find(f"{SBML_NS}listOfReactants")
                if reactants:
                    for ref in reactants.findall(f"{SBML_NS}speciesReference"):
                        reaction.reactants.append({
                            "species": ref.get("species"),
                            "stoichiometry": float(ref.get("stoichiometry", 1))
                        })

                # Products
                products = rxn.find(f"{SBML_NS}listOfProducts")
                if products:
                    for ref in products.findall(f"{SBML_NS}speciesReference"):
                        reaction.products.append({
                            "species": ref.get("species"),
                            "stoichiometry": float(ref.get("stoichiometry", 1))
                        })

                result["reactions"].append(reaction)

        return result

# Usage
parser = SBMLParser()
with open("model.xml", "r") as f:
    model = parser.parse(f.read())

print(f"Model: {model['name']}")
print(f"Species: {len(model['species'])}")
print(f"Reactions: {len(model['reactions'])}")
```

---

## 7. PSI-MI (Molecular Interactions)

### 7.1 Overview

PSI-MI (Proteomics Standards Initiative - Molecular Interactions) defines XML and TAB formats for molecular interaction data. Used by IntAct, BioGRID, and other databases.

### 7.2 PSI-MI XML 3.0 Structure

```xml
<?xml version="1.0" encoding="UTF-8"?>
<entrySet xmlns="http://psi.hupo.org/mi/mif300"
          xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
          level="3" version="0" minorVersion="0">
  <entry>
    <source>
      <names>
        <shortLabel>IntAct</shortLabel>
        <fullName>IntAct molecular interaction database</fullName>
      </names>
      <bibref>
        <xref>
          <primaryRef db="pubmed" dbAc="MI:0446" id="24234451"/>
        </xref>
      </bibref>
    </source>

    <interactorList>
      <interactor id="1">
        <names>
          <shortLabel>TP53_HUMAN</shortLabel>
          <fullName>Cellular tumor antigen p53</fullName>
          <alias type="gene name" typeAc="MI:0301">TP53</alias>
        </names>
        <xref>
          <primaryRef db="uniprotkb" dbAc="MI:0486" id="P04637" version="312"/>
          <secondaryRef db="intact" dbAc="MI:0469" id="EBI-366083"/>
          <secondaryRef db="ensembl" dbAc="MI:0476" id="ENSG00000141510"/>
        </xref>
        <interactorType>
          <names>
            <shortLabel>protein</shortLabel>
          </names>
          <xref>
            <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0326"/>
          </xref>
        </interactorType>
        <organism ncbiTaxId="9606">
          <names>
            <shortLabel>human</shortLabel>
            <fullName>Homo sapiens</fullName>
          </names>
        </organism>
        <sequence>MEEPQSDPSVEPPLSQETFSDLWKLL...</sequence>
      </interactor>

      <interactor id="2">
        <names>
          <shortLabel>MDM2_HUMAN</shortLabel>
          <fullName>E3 ubiquitin-protein ligase Mdm2</fullName>
        </names>
        <xref>
          <primaryRef db="uniprotkb" dbAc="MI:0486" id="Q00987"/>
        </xref>
        <interactorType>
          <names><shortLabel>protein</shortLabel></names>
          <xref>
            <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0326"/>
          </xref>
        </interactorType>
        <organism ncbiTaxId="9606">
          <names><shortLabel>human</shortLabel></names>
        </organism>
      </interactor>
    </interactorList>

    <interactionList>
      <interaction id="100">
        <names>
          <shortLabel>tp53-mdm2</shortLabel>
          <fullName>TP53 interacts with MDM2</fullName>
        </names>
        <xref>
          <primaryRef db="intact" dbAc="MI:0469" id="EBI-12345"/>
          <secondaryRef db="imex" dbAc="MI:0670" id="IM-12345-1"/>
        </xref>
        <experimentList>
          <experimentRef>1</experimentRef>
        </experimentList>
        <participantList>
          <participant id="10">
            <interactorRef>1</interactorRef>
            <biologicalRole>
              <names><shortLabel>unspecified role</shortLabel></names>
              <xref>
                <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0499"/>
              </xref>
            </biologicalRole>
            <experimentalRoleList>
              <experimentalRole>
                <names><shortLabel>prey</shortLabel></names>
                <xref>
                  <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0498"/>
                </xref>
              </experimentalRole>
            </experimentalRoleList>
          </participant>
          <participant id="11">
            <interactorRef>2</interactorRef>
            <biologicalRole>
              <names><shortLabel>unspecified role</shortLabel></names>
              <xref>
                <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0499"/>
              </xref>
            </biologicalRole>
            <experimentalRoleList>
              <experimentalRole>
                <names><shortLabel>bait</shortLabel></names>
                <xref>
                  <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0496"/>
                </xref>
              </experimentalRole>
            </experimentalRoleList>
          </participant>
        </participantList>
        <interactionType>
          <names>
            <shortLabel>physical association</shortLabel>
          </names>
          <xref>
            <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0915"/>
          </xref>
        </interactionType>
        <confidenceList>
          <confidence>
            <unit>
              <names><shortLabel>intact-miscore</shortLabel></names>
            </unit>
            <value>0.85</value>
          </confidence>
        </confidenceList>
      </interaction>
    </interactionList>

    <experimentList>
      <experimentDescription id="1">
        <names>
          <shortLabel>smith-2020-1</shortLabel>
        </names>
        <bibref>
          <xref>
            <primaryRef db="pubmed" dbAc="MI:0446" id="32000000"/>
          </xref>
        </bibref>
        <interactionDetectionMethod>
          <names>
            <shortLabel>two hybrid</shortLabel>
          </names>
          <xref>
            <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0018"/>
          </xref>
        </interactionDetectionMethod>
        <hostOrganismList>
          <hostOrganism ncbiTaxId="4932">
            <names><shortLabel>yeast</shortLabel></names>
          </hostOrganism>
        </hostOrganismList>
      </experimentDescription>
    </experimentList>

  </entry>
</entrySet>
```

### 7.3 PSI-MITAB 2.8 Format

Tab-delimited format with 42 columns (extensible).

| Col | Name | Description | Example |
|-----|------|-------------|---------|
| 1 | ID(s) interactor A | Unique identifier | `uniprotkb:P04637` |
| 2 | ID(s) interactor B | Unique identifier | `uniprotkb:Q00987` |
| 3 | Alt. ID(s) A | Alternative identifiers | `intact:EBI-366083` |
| 4 | Alt. ID(s) B | Alternative identifiers | `intact:EBI-389668` |
| 5 | Alias(es) A | Aliases | `uniprotkb:TP53(gene name)` |
| 6 | Alias(es) B | Aliases | `uniprotkb:MDM2(gene name)` |
| 7 | Detection method | MI term | `psi-mi:"MI:0018"(two hybrid)` |
| 8 | Publication | Publication IDs | `pubmed:32000000` |
| 9 | Taxid A | NCBI taxonomy | `taxid:9606(human)` |
| 10 | Taxid B | NCBI taxonomy | `taxid:9606(human)` |
| 11 | Interaction type | MI term | `psi-mi:"MI:0915"(physical association)` |
| 12 | Source database | MI term | `psi-mi:"MI:0469"(IntAct)` |
| 13 | Interaction ID | Unique ID | `intact:EBI-12345` |
| 14 | Confidence score | Score value | `intact-miscore:0.85` |
| 15-42 | Extended fields | Various | See specification |

#### MITAB Example

```tsv
uniprotkb:P04637	uniprotkb:Q00987	intact:EBI-366083	intact:EBI-389668	uniprotkb:TP53(gene name)	uniprotkb:MDM2(gene name)	psi-mi:"MI:0018"(two hybrid)	pubmed:32000000	taxid:9606(human)	taxid:9606(human)	psi-mi:"MI:0915"(physical association)	psi-mi:"MI:0469"(IntAct)	intact:EBI-12345	intact-miscore:0.85	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
```

### 7.4 Python Parsing Code

```python
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import List, Dict, Optional

PSI_NS = "{http://psi.hupo.org/mi/mif300}"

@dataclass
class Interactor:
    id: str
    short_label: str
    full_name: str = ""
    uniprot_id: Optional[str] = None
    organism: Optional[int] = None
    interactor_type: str = "protein"

@dataclass
class Interaction:
    id: str
    participants: List[str] = field(default_factory=list)
    interaction_type: str = ""
    detection_method: str = ""
    confidence: Optional[float] = None
    pubmed_ids: List[str] = field(default_factory=list)

class PSIMIParser:
    """Parser for PSI-MI XML format."""

    def parse_xml(self, content: str) -> Dict:
        root = ET.fromstring(content)

        result = {
            "interactors": {},
            "interactions": []
        }

        for entry in root.findall(f".//{PSI_NS}entry"):
            # Parse interactors
            for interactor in entry.findall(f".//{PSI_NS}interactor"):
                iid = interactor.get("id")
                names = interactor.find(f"{PSI_NS}names")
                short_label = names.findtext(f"{PSI_NS}shortLabel", "")
                full_name = names.findtext(f"{PSI_NS}fullName", "")

                # Get UniProt ID
                uniprot_id = None
                xref = interactor.find(f"{PSI_NS}xref")
                if xref is not None:
                    primary = xref.find(f"{PSI_NS}primaryRef")
                    if primary is not None and primary.get("db") == "uniprotkb":
                        uniprot_id = primary.get("id")

                # Get organism
                organism = interactor.find(f"{PSI_NS}organism")
                taxid = int(organism.get("ncbiTaxId")) if organism is not None else None

                result["interactors"][iid] = Interactor(
                    id=iid,
                    short_label=short_label,
                    full_name=full_name,
                    uniprot_id=uniprot_id,
                    organism=taxid
                )

            # Parse interactions
            for interaction in entry.findall(f".//{PSI_NS}interaction"):
                iid = interaction.get("id")

                participants = []
                for participant in interaction.findall(f".//{PSI_NS}participant"):
                    ref = participant.find(f"{PSI_NS}interactorRef")
                    if ref is not None:
                        participants.append(ref.text)

                int_type = ""
                type_elem = interaction.find(f"{PSI_NS}interactionType")
                if type_elem is not None:
                    names = type_elem.find(f"{PSI_NS}names")
                    if names is not None:
                        int_type = names.findtext(f"{PSI_NS}shortLabel", "")

                result["interactions"].append(Interaction(
                    id=iid,
                    participants=participants,
                    interaction_type=int_type
                ))

        return result

    def parse_mitab(self, content: str) -> List[Dict]:
        interactions = []

        for line in content.strip().split('\n'):
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 14:
                continue

            interactions.append({
                "interactor_a": self._parse_id(fields[0]),
                "interactor_b": self._parse_id(fields[1]),
                "alias_a": fields[4] if len(fields) > 4 else "",
                "alias_b": fields[5] if len(fields) > 5 else "",
                "detection_method": self._parse_mi(fields[6]) if len(fields) > 6 else "",
                "publication": fields[7] if len(fields) > 7 else "",
                "taxid_a": self._parse_taxid(fields[8]) if len(fields) > 8 else None,
                "taxid_b": self._parse_taxid(fields[9]) if len(fields) > 9 else None,
                "interaction_type": self._parse_mi(fields[10]) if len(fields) > 10 else "",
                "source": self._parse_mi(fields[11]) if len(fields) > 11 else "",
                "interaction_id": fields[12] if len(fields) > 12 else "",
                "confidence": self._parse_confidence(fields[13]) if len(fields) > 13 else None
            })

        return interactions

    def _parse_id(self, field: str) -> Dict:
        if ':' in field:
            db, acc = field.split(':', 1)
            return {"db": db, "id": acc}
        return {"db": "", "id": field}

    def _parse_mi(self, field: str) -> str:
        # psi-mi:"MI:0018"(two hybrid) -> two hybrid
        if '(' in field and ')' in field:
            return field.split('(')[1].rstrip(')')
        return field

    def _parse_taxid(self, field: str) -> Optional[int]:
        # taxid:9606(human) -> 9606
        if field.startswith('taxid:'):
            taxid_str = field.replace('taxid:', '').split('(')[0]
            return int(taxid_str) if taxid_str.isdigit() else None
        return None

    def _parse_confidence(self, field: str) -> Optional[float]:
        # intact-miscore:0.85 -> 0.85
        if ':' in field:
            try:
                return float(field.split(':')[1])
            except ValueError:
                return None
        return None

# Usage
parser = PSIMIParser()

# Parse XML
with open("interactions.xml", "r") as f:
    data = parser.parse_xml(f.read())
print(f"Interactors: {len(data['interactors'])}")
print(f"Interactions: {len(data['interactions'])}")

# Parse MITAB
with open("interactions.txt", "r") as f:
    interactions = parser.parse_mitab(f.read())
print(f"MITAB interactions: {len(interactions)}")
```

---

## 8. CX Format (NDEx)

### 8.1 Overview

CX (Cytoscape Exchange) is a JSON-based format for biological networks used by NDEx (Network Data Exchange).

### 8.2 CX Structure

CX uses an "aspects" model where different types of data are stored in separate arrays.

```json
[
  {"numberVerification": [{"longNumber": 281474976710655}]},

  {"metaData": [
    {"name": "nodes", "elementCount": 3, "version": "1.0"},
    {"name": "edges", "elementCount": 2, "version": "1.0"},
    {"name": "nodeAttributes", "elementCount": 6, "version": "1.0"},
    {"name": "edgeAttributes", "elementCount": 4, "version": "1.0"},
    {"name": "networkAttributes", "elementCount": 3, "version": "1.0"}
  ]},

  {"networkAttributes": [
    {"n": "name", "v": "Apoptosis Pathway"},
    {"n": "description", "v": "Core apoptosis signaling network"},
    {"n": "organism", "v": "Homo sapiens"}
  ]},

  {"nodes": [
    {"@id": 0, "n": "TP53"},
    {"@id": 1, "n": "MDM2"},
    {"@id": 2, "n": "BAX"}
  ]},

  {"edges": [
    {"@id": 100, "s": 0, "t": 1, "i": "inhibits"},
    {"@id": 101, "s": 0, "t": 2, "i": "activates"}
  ]},

  {"nodeAttributes": [
    {"po": 0, "n": "uniprot", "v": "P04637"},
    {"po": 0, "n": "gene_name", "v": "TP53"},
    {"po": 1, "n": "uniprot", "v": "Q00987"},
    {"po": 1, "n": "gene_name", "v": "MDM2"},
    {"po": 2, "n": "uniprot", "v": "Q07812"},
    {"po": 2, "n": "gene_name", "v": "BAX"}
  ]},

  {"edgeAttributes": [
    {"po": 100, "n": "interaction_type", "v": "inhibition"},
    {"po": 100, "n": "pubmed", "v": "8875929"},
    {"po": 101, "n": "interaction_type", "v": "activation"},
    {"po": 101, "n": "pubmed", "v": "9926942"}
  ]},

  {"cyVisualProperties": [
    {
      "properties_of": "network",
      "properties": {
        "NETWORK_BACKGROUND_PAINT": "#FFFFFF"
      }
    },
    {
      "properties_of": "nodes:default",
      "properties": {
        "NODE_SHAPE": "ellipse",
        "NODE_SIZE": 50,
        "NODE_FILL_COLOR": "#CCCCCC"
      }
    }
  ]},

  {"status": [{"error": "", "success": true}]}
]
```

### 8.3 Core Aspects

#### nodes

```json
{"nodes": [
  {"@id": 0, "n": "TP53", "r": "uniprot:P04637"}
]}
```

| Field | Type | Description |
|-------|------|-------------|
| `@id` | Integer | Unique node ID |
| `n` | String | Node name |
| `r` | String | Node represents (external ID) |

#### edges

```json
{"edges": [
  {"@id": 100, "s": 0, "t": 1, "i": "inhibits"}
]}
```

| Field | Type | Description |
|-------|------|-------------|
| `@id` | Integer | Unique edge ID |
| `s` | Integer | Source node ID |
| `t` | Integer | Target node ID |
| `i` | String | Interaction type |

#### nodeAttributes

```json
{"nodeAttributes": [
  {"po": 0, "n": "uniprot", "v": "P04637", "d": "string"}
]}
```

| Field | Type | Description |
|-------|------|-------------|
| `po` | Integer | Node ID (property of) |
| `n` | String | Attribute name |
| `v` | Any | Attribute value |
| `d` | String | Data type (optional) |

#### edgeAttributes

```json
{"edgeAttributes": [
  {"po": 100, "n": "score", "v": 0.95, "d": "double"}
]}
```

| Field | Type | Description |
|-------|------|-------------|
| `po` | Integer | Edge ID (property of) |
| `n` | String | Attribute name |
| `v` | Any | Attribute value |
| `d` | String | Data type |

#### networkAttributes

```json
{"networkAttributes": [
  {"n": "name", "v": "My Network"},
  {"n": "description", "v": "Network description"},
  {"n": "version", "v": "1.0"}
]}
```

### 8.4 Data Types

| Type | Description | Example |
|------|-------------|---------|
| `string` | String value | `"TP53"` |
| `boolean` | Boolean value | `true` |
| `double` | Floating point | `0.95` |
| `integer` | Integer value | `42` |
| `long` | Long integer | `281474976710655` |
| `list_of_string` | String list | `["A", "B"]` |
| `list_of_double` | Double list | `[1.0, 2.0]` |

### 8.5 Python Parsing Code

```python
import json
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional

@dataclass
class CXNode:
    id: int
    name: str
    represents: Optional[str] = None
    attributes: Dict[str, Any] = field(default_factory=dict)

@dataclass
class CXEdge:
    id: int
    source: int
    target: int
    interaction: str = ""
    attributes: Dict[str, Any] = field(default_factory=dict)

@dataclass
class CXNetwork:
    name: str = ""
    description: str = ""
    nodes: Dict[int, CXNode] = field(default_factory=dict)
    edges: Dict[int, CXEdge] = field(default_factory=dict)
    network_attributes: Dict[str, Any] = field(default_factory=dict)

class CXParser:
    """Parser for NDEx CX format."""

    def parse(self, cx_content: str) -> CXNetwork:
        """Parse CX JSON content."""
        data = json.loads(cx_content)
        return self.parse_aspects(data)

    def parse_aspects(self, aspects: List[Dict]) -> CXNetwork:
        """Parse list of CX aspects."""
        network = CXNetwork()

        for aspect_wrapper in aspects:
            for aspect_name, aspect_data in aspect_wrapper.items():
                if aspect_name == "networkAttributes":
                    self._parse_network_attributes(network, aspect_data)
                elif aspect_name == "nodes":
                    self._parse_nodes(network, aspect_data)
                elif aspect_name == "edges":
                    self._parse_edges(network, aspect_data)
                elif aspect_name == "nodeAttributes":
                    self._parse_node_attributes(network, aspect_data)
                elif aspect_name == "edgeAttributes":
                    self._parse_edge_attributes(network, aspect_data)

        return network

    def _parse_network_attributes(self, network: CXNetwork, data: List[Dict]):
        for attr in data:
            name = attr.get("n")
            value = attr.get("v")
            if name == "name":
                network.name = value
            elif name == "description":
                network.description = value
            else:
                network.network_attributes[name] = value

    def _parse_nodes(self, network: CXNetwork, data: List[Dict]):
        for node_data in data:
            node_id = node_data.get("@id")
            node = CXNode(
                id=node_id,
                name=node_data.get("n", ""),
                represents=node_data.get("r")
            )
            network.nodes[node_id] = node

    def _parse_edges(self, network: CXNetwork, data: List[Dict]):
        for edge_data in data:
            edge_id = edge_data.get("@id")
            edge = CXEdge(
                id=edge_id,
                source=edge_data.get("s"),
                target=edge_data.get("t"),
                interaction=edge_data.get("i", "")
            )
            network.edges[edge_id] = edge

    def _parse_node_attributes(self, network: CXNetwork, data: List[Dict]):
        for attr in data:
            node_id = attr.get("po")
            name = attr.get("n")
            value = attr.get("v")
            if node_id in network.nodes:
                network.nodes[node_id].attributes[name] = value

    def _parse_edge_attributes(self, network: CXNetwork, data: List[Dict]):
        for attr in data:
            edge_id = attr.get("po")
            name = attr.get("n")
            value = attr.get("v")
            if edge_id in network.edges:
                network.edges[edge_id].attributes[name] = value

    def to_networkx(self, network: CXNetwork):
        """Convert CX network to NetworkX graph."""
        import networkx as nx

        G = nx.DiGraph()
        G.graph["name"] = network.name
        G.graph["description"] = network.description
        G.graph.update(network.network_attributes)

        for node_id, node in network.nodes.items():
            G.add_node(node_id, name=node.name, represents=node.represents, **node.attributes)

        for edge_id, edge in network.edges.items():
            G.add_edge(edge.source, edge.target, id=edge_id, interaction=edge.interaction, **edge.attributes)

        return G

    def get_node_by_attribute(self, network: CXNetwork, attr_name: str, attr_value: Any) -> Optional[CXNode]:
        """Find node by attribute value."""
        for node in network.nodes.values():
            if node.attributes.get(attr_name) == attr_value:
                return node
        return None

# Usage
parser = CXParser()

# Parse from file
with open("network.cx", "r") as f:
    network = parser.parse(f.read())

print(f"Network: {network.name}")
print(f"Nodes: {len(network.nodes)}")
print(f"Edges: {len(network.edges)}")

# Find node by UniProt ID
tp53 = parser.get_node_by_attribute(network, "uniprot", "P04637")
if tp53:
    print(f"Found TP53: {tp53.name}")

# Convert to NetworkX
import networkx as nx
G = parser.to_networkx(network)
print(f"NetworkX nodes: {G.number_of_nodes()}")
print(f"NetworkX edges: {G.number_of_edges()}")
```

### 8.6 NDEx API Integration

```python
import requests
from typing import Optional

class NDExClient:
    """Simple client for NDEx REST API."""

    BASE_URL = "https://www.ndexbio.org/v2"

    def __init__(self, username: Optional[str] = None, password: Optional[str] = None):
        self.auth = (username, password) if username and password else None

    def get_network(self, network_uuid: str) -> CXNetwork:
        """Fetch network by UUID."""
        url = f"{self.BASE_URL}/network/{network_uuid}"
        response = requests.get(url, auth=self.auth)
        response.raise_for_status()

        parser = CXParser()
        return parser.parse_aspects(response.json())

    def search_networks(self, search_string: str, start: int = 0, size: int = 25) -> List[Dict]:
        """Search for networks."""
        url = f"{self.BASE_URL}/search/network"
        params = {"searchString": search_string, "start": start, "size": size}
        response = requests.get(url, params=params, auth=self.auth)
        response.raise_for_status()
        return response.json().get("networks", [])

# Usage
client = NDExClient()

# Search for apoptosis networks
results = client.search_networks("apoptosis")
print(f"Found {len(results)} networks")

for result in results[:5]:
    print(f"  {result['name']}: {result['nodeCount']} nodes, {result['edgeCount']} edges")

# Get specific network
if results:
    uuid = results[0]['externalId']
    network = client.get_network(uuid)
    print(f"Loaded: {network.name}")
```

---

## Summary: Format Comparison

| Format | Type | Primary Use | Databases |
|--------|------|-------------|-----------|
| Reactome Neo4j | Graph DB | Pathway curation | Reactome |
| KGML | XML | Pathway visualization | KEGG |
| GPML | XML | Pathway editing | WikiPathways |
| OBO/GAF | Text | Ontology/Annotations | GO, other ontologies |
| BioPAX | OWL/RDF | Pathway exchange | Reactome, PathwayCommons |
| SBML | XML | Kinetic modeling | BioModels |
| PSI-MI | XML/TAB | Interactions | IntAct, BioGRID |
| CX | JSON | Network exchange | NDEx |

---

## Cross-References

- See also: [pathway-target-data-models.md](pathway-target-data-models.md) for API documentation
- See also: [drug-target-interaction-models.md](drug-target-interaction-models.md) for drug-target data

---

*Document generated: 2025*
*Last updated: January 2025*
