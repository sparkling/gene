# Reactome Neo4j Graph Database Schema

**Source:** https://reactome.org
**Version:** 95 (as of January 2026)
**License:** CC BY 4.0
**Format:** Neo4j Graph Database

---

## Overview

Reactome represents biological pathways as an interconnected graph where molecular reactions convert input physical entities to output entities. The database uses a hierarchical event structure with explicit representation of modifications, complexes, and cellular compartments.

---

## Core Entity Types (Node Classes)

### 1. Event Classes

Events represent biological processes that convert inputs to outputs.

#### TopLevelPathway
Top-level entry points in the pathway hierarchy.

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `stId` | String | Stable identifier (e.g., "R-HSA-1430728") |
| `stIdVersion` | String | Versioned stable ID |
| `displayName` | String | Human-readable name |
| `name` | String[] | Alternative names |
| `speciesName` | String | Species (e.g., "Homo sapiens") |
| `releaseDate` | Date | First release date |
| `lastUpdatedDate` | Date | Last modification date |
| `isInDisease` | Boolean | Disease pathway flag |
| `isInferred` | Boolean | Inferred from orthology |
| `doi` | String | Digital Object Identifier |
| `hasDiagram` | Boolean | Has pathway diagram |
| `hasEHLD` | Boolean | Has Enhanced High-Level Diagram |
| `maxDepth` | Integer | Maximum hierarchy depth |

**Sample JSON:**
```json
{
  "dbId": 1430728,
  "displayName": "Metabolism",
  "stId": "R-HSA-1430728",
  "stIdVersion": "R-HSA-1430728.16",
  "isInDisease": false,
  "isInferred": false,
  "maxDepth": 7,
  "speciesName": "Homo sapiens",
  "doi": "10.3180/R-HSA-1430728.15",
  "hasDiagram": true,
  "hasEHLD": true,
  "schemaClass": "TopLevelPathway"
}
```

#### Pathway
Groupings of related events.

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `stId` | String | Stable identifier |
| `displayName` | String | Pathway name |
| `compartment` | Compartment[] | Cellular locations |
| `goBiologicalProcess` | GO_BiologicalProcess | GO annotation |
| `literatureReference` | LiteratureReference[] | PubMed citations |
| `summation` | Summation[] | Pathway description |
| `reviewStatus` | ReviewStatus | Curation status |
| `hasEvent` | Event[] | Child events |

#### Reaction / BlackBoxEvent / Polymerisation / Depolymerisation
Direct input-to-output conversions.

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `stId` | String | Stable identifier |
| `displayName` | String | Reaction name |
| `isChimeric` | Boolean | Mixed evidence sources |
| `category` | String | "transition", "omitted", etc. |
| `input` | PhysicalEntity[] | Reaction inputs |
| `output` | PhysicalEntity[] | Reaction outputs |
| `catalystActivity` | CatalystActivity[] | Enzymes |

---

### 2. PhysicalEntity Classes

Physical entities are molecules, complexes, and molecular sets.

#### EntityWithAccessionedSequence
Proteins and nucleic acids with known sequences.

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `stId` | String | Stable identifier |
| `displayName` | String | Entity name |
| `compartment` | Compartment | Cellular location |
| `referenceEntity` | ReferenceGeneProduct | UniProt reference |

**Sample JSON:**
```json
{
  "peDbId": 196423,
  "displayName": "DHCR7 [endoplasmic reticulum membrane]",
  "schemaClass": "EntityWithAccessionedSequence",
  "refEntities": [{
    "dbId": 53616,
    "stId": "uniprot:Q9UBM7",
    "identifier": "Q9UBM7",
    "schemaClass": "ReferenceGeneProduct",
    "displayName": "UniProt:Q9UBM7 DHCR7",
    "url": "http://purl.uniprot.org/uniprot/Q9UBM7"
  }]
}
```

#### SimpleEntity (Small Molecules)
Fully characterized non-protein molecules.

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `stId` | String | Stable identifier |
| `displayName` | String | Molecule name with compartment |
| `name` | String[] | Alternative names |
| `referenceEntity` | ReferenceMolecule | ChEBI reference |
| `compartment` | Compartment | Cellular location |

**Sample JSON (ATP):**
```json
{
  "dbId": 29358,
  "displayName": "ATP [nucleoplasm]",
  "stId": "R-ALL-29358",
  "name": ["ATP", "Adenosine 5'-triphosphate", "ATP(4-)"],
  "referenceEntity": {
    "dbId": 8869364,
    "stId": "chebi:30616",
    "databaseName": "ChEBI",
    "identifier": "30616",
    "name": ["ATP(4-)", "ATP", "atp"],
    "formula": "C10H12N5O13P3",
    "schemaClass": "ReferenceMolecule"
  },
  "schemaClass": "SimpleEntity"
}
```

#### Complex
Multi-component assemblies.

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `displayName` | String | Complex name |
| `hasComponent` | PhysicalEntity[] | Component entities |
| `compartment` | Compartment | Cellular location |

#### EntitySet / DefinedSet / CandidateSet
Functionally interchangeable entities.

| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `displayName` | String | Set name |
| `hasMember` | PhysicalEntity[] | Member entities |

---

### 3. Reference Entity Classes

Invariant features shared across variant forms.

#### ReferenceGeneProduct (Proteins)
| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `stId` | String | UniProt accession |
| `identifier` | String | UniProt ID |
| `databaseName` | String | "UniProt" |
| `url` | String | UniProt URL |

#### ReferenceMolecule (Small Molecules)
| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `stId` | String | ChEBI ID |
| `identifier` | String | ChEBI number |
| `databaseName` | String | "ChEBI" |
| `formula` | String | Molecular formula |
| `moleculeType` | String | "Chemical", "Drug", etc. |
| `url` | String | ChEBI URL |

---

### 4. Annotation Classes

#### Compartment (GO Cellular Component)
| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `displayName` | String | Compartment name |
| `accession` | String | GO ID |
| `databaseName` | String | "GO" |
| `definition` | String | GO definition |

**Sample Compartments:**
- nucleoplasm (GO:0005654)
- cytosol (GO:0005829)
- endoplasmic reticulum membrane (GO:0005789)
- plasma membrane (GO:0005886)
- mitochondrial matrix (GO:0005759)

#### GO_BiologicalProcess
| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `displayName` | String | Process name |
| `accession` | String | GO ID |
| `definition` | String | GO definition |
| `url` | String | QuickGO URL |

#### LiteratureReference
| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `title` | String | Article title |
| `journal` | String | Journal name |
| `pubMedIdentifier` | Integer | PubMed ID |
| `year` | Integer | Publication year |
| `volume` | Integer | Journal volume |
| `pages` | String | Page numbers |

#### Species
| Property | Type | Description |
|----------|------|-------------|
| `dbId` | Integer | Internal database ID |
| `displayName` | String | Species name |
| `name` | String[] | All names |
| `taxId` | String | NCBI Taxonomy ID |
| `abbreviation` | String | 3-letter code |

**Supported Species (16 main):**
| Species | Tax ID | Abbreviation |
|---------|--------|--------------|
| Homo sapiens | 9606 | HSA |
| Mus musculus | 10090 | MMU |
| Rattus norvegicus | 10116 | RNO |
| Danio rerio | 7955 | DRE |
| Drosophila melanogaster | 7227 | DME |
| Caenorhabditis elegans | 6239 | CEL |
| Saccharomyces cerevisiae | 4932 | SCE |
| Bos taurus | 9913 | BTA |
| Canis familiaris | 9615 | CFA |
| Gallus gallus | 9031 | GGA |
| Sus scrofa | 9823 | SSC |
| Xenopus tropicalis | 8364 | XTR |
| Schizosaccharomyces pombe | 4896 | SPO |
| Dictyostelium discoideum | 44689 | DDI |
| Plasmodium falciparum | 5833 | PFA |
| Mycobacterium tuberculosis | 1773 | MTU |

---

## Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `hasEvent` | Pathway | Event | Pathway contains events |
| `eventOf` | Event | Pathway | Event belongs to pathway |
| `input` | Reaction | PhysicalEntity | Reaction input |
| `output` | Reaction | PhysicalEntity | Reaction output |
| `hasComponent` | Complex | PhysicalEntity | Complex component |
| `hasMember` | EntitySet | PhysicalEntity | Set member |
| `referenceEntity` | PhysicalEntity | ReferenceEntity | External reference |
| `compartment` | PhysicalEntity | Compartment | Localization |
| `orthologousEvent` | Event | Event | Ortholog mapping |
| `inferredFrom` | Event | Event | Inference source |

---

## Neo4j Node Type Hierarchy

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

---

## Cypher Query Examples

### Get Pathway with Sub-events
```cypher
MATCH (p:Pathway {stId: 'R-HSA-109581'})-[:hasEvent*]->(e:Event)
RETURN p.displayName AS Pathway, collect(DISTINCT e.displayName) AS Events
```

### Get Reaction with Inputs, Outputs, Catalysts
```cypher
MATCH (r:Reaction {stId: 'R-HSA-69541'})
OPTIONAL MATCH (r)-[:input]->(i:PhysicalEntity)
OPTIONAL MATCH (r)-[:output]->(o:PhysicalEntity)
OPTIONAL MATCH (r)-[:catalystActivity]->(ca:CatalystActivity)-[:physicalEntity]->(c:PhysicalEntity)
RETURN r.displayName AS Reaction,
       collect(DISTINCT i.displayName) AS Inputs,
       collect(DISTINCT o.displayName) AS Outputs,
       collect(DISTINCT c.displayName) AS Catalysts
```

### Get Protein with All Pathways
```cypher
MATCH (re:ReferenceEntity {identifier: 'P04637'})<-[:referenceEntity]-(pe:PhysicalEntity)
MATCH (pe)<-[:hasComponent|hasMember|hasCandidate|input|output*]-(e:Event)
MATCH (e)<-[:hasEvent*]-(p:TopLevelPathway)
RETURN DISTINCT p.displayName AS Pathway, p.stId AS PathwayId
```

### Get All Proteins in a Pathway
```cypher
MATCH (p:Pathway {stId: 'R-HSA-109581'})-[:hasEvent*]->(r:ReactionLikeEvent)
MATCH (r)-[:input|output]->(pe:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
WHERE re.databaseName = 'UniProt'
RETURN DISTINCT re.identifier AS UniProtID, re.geneName AS GeneName
```

---

## Extended Relationship Types

### Event Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `hasEvent` | Pathway | Event | Sub-events within pathway |
| `precedingEvent` | Event | Event | Temporal ordering |
| `followingEvent` | Event | Event | Reverse temporal ordering |
| `inferredTo` | Event | Event | Orthology inference |

### Reaction Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `input` | ReactionLikeEvent | PhysicalEntity | Reaction inputs |
| `output` | ReactionLikeEvent | PhysicalEntity | Reaction outputs |
| `catalystActivity` | ReactionLikeEvent | CatalystActivity | Enzymatic catalysis |
| `regulatedBy` | ReactionLikeEvent | Regulation | Regulatory relationships |
| `normalReaction` | FailedReaction | Reaction | Normal version of failed reaction |

### CatalystActivity Relationships

| Relationship | From | To | Description |
|--------------|------|-----|-------------|
| `physicalEntity` | CatalystActivity | PhysicalEntity | Catalyzing entity |
| `activity` | CatalystActivity | GO_MolecularFunction | GO activity term |
| `activeUnit` | CatalystActivity | PhysicalEntity | Active subunit |

---

## Statistics (Release 95, December 2024)

| Metric | Count |
|--------|-------|
| Human Pathways | 2,712 |
| Human Reactions | 13,872 |
| Human Proteins | 11,196 |
| Small Molecules | 1,925 |
| Species | 24 |

---

## API Endpoints

| Endpoint | Description |
|----------|-------------|
| `GET /data/query/{id}` | Fetch entity by ID |
| `GET /data/pathway/{id}/containedEvents` | Get pathway events |
| `GET /data/participants/{id}` | Get reaction participants |
| `GET /data/species/main` | List supported species |
| `GET /data/database/version` | Get database version |
| `GET /data/discover/{id}` | Schema.org metadata |

**Base URL:** `https://reactome.org/ContentService`

---

## Export Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| BioPAX Level 3 | OWL-based pathway exchange | Integration |
| SBML Level 3 | Systems biology models | Simulation |
| SBGN-ML | Graphical notation | Visualization |
| PSI-MITAB | Protein interactions | Network analysis |
| JSON | API responses | Programmatic access |

---

## Cross-References

| Database | ID Format | Example |
|----------|-----------|---------|
| UniProt | Q##### | Q9UBM7 |
| ChEBI | CHEBI:##### | CHEBI:30616 |
| GO | GO:####### | GO:0005654 |
| Ensembl | ENSG########### | ENSG00000113161 |
| KEGG | C#####, map##### | C00002 |
| PubMed | ######## | 10583946 |

---

## Python Integration

```python
import requests
from typing import Dict, List

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

# Usage
parser = ReactomeParser()
pathway = parser.get_pathway("R-HSA-109581")
print(f"Pathway: {pathway['displayName']}")
```

---

## BioPAX Level 3 Export

Reactome exports data in BioPAX Level 3 OWL format:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
         xmlns:bp="http://www.biopax.org/release/biopax-level3.owl#"
         xml:base="http://www.reactome.org/biopax/">

  <bp:Pathway rdf:ID="Pathway109581">
    <bp:displayName>Apoptosis</bp:displayName>
    <bp:organism rdf:resource="#BioSource1"/>
    <bp:pathwayComponent rdf:resource="#BiochemicalReaction69541"/>
  </bp:Pathway>

  <bp:BiochemicalReaction rdf:ID="BiochemicalReaction69541">
    <bp:displayName>Cytochrome c binds APAF1</bp:displayName>
    <bp:left rdf:resource="#Protein69523"/>
    <bp:right rdf:resource="#Complex69543"/>
    <bp:conversionDirection>LEFT-TO-RIGHT</bp:conversionDirection>
  </bp:BiochemicalReaction>

  <bp:Protein rdf:ID="Protein69523">
    <bp:displayName>CYCS</bp:displayName>
    <bp:entityReference rdf:resource="#ProteinReference_P99999"/>
  </bp:Protein>

</rdf:RDF>
```
