# WikiPathways GPML Schema Documentation

**Source:** https://www.wikipathways.org
**Format:** GPML (Graphical Pathway Markup Language)
**Version:** GPML2013a (current) / GPML2021 (latest)
**License:** CC0 1.0 Universal (Public Domain)

---

## Overview

GPML (Graphical Pathway Markup Language) is the native XML format for WikiPathways, designed to store both biological semantics and graphical layout information for pathway diagrams.

---

## Root Element: Pathway

```xml
<Pathway xmlns="http://pathvisio.org/GPML/2013a"
         Name="Cholesterol biosynthesis pathway"
         Version="20200316"
         Organism="Homo sapiens">
  <!-- Child elements -->
</Pathway>
```

| Attribute | Type | Required | Description |
|-----------|------|----------|-------------|
| `Name` | String | Yes | Pathway title |
| `Version` | String | No | Version date (YYYYMMDD) |
| `Organism` | String | Yes | Species name |
| `xmlns` | URI | Yes | GPML namespace |

---

## DataNode Elements

DataNodes represent biological entities (genes, proteins, metabolites, pathways).

### DataNode Structure

```xml
<DataNode TextLabel="HMGCR" GraphId="abc123" Type="GeneProduct">
  <Comment Source="GenMAPP remarks">Rate-limiting enzyme</Comment>
  <Attribute Key="org.pathvisio.model.BackpageHead" Value="HMGCR | ..." />
  <Graphics CenterX="218.39" CenterY="190.95" Width="66.667" Height="20.0"
            ZOrder="32768" FontSize="12" Valign="Middle" />
  <Xref Database="Entrez Gene" ID="3156" />
</DataNode>
```

| Attribute | Type | Required | Description |
|-----------|------|----------|-------------|
| `TextLabel` | String | Yes | Display text |
| `GraphId` | String | Yes | Unique node identifier |
| `Type` | String | Yes | Node type (see below) |
| `GroupRef` | String | No | Group membership |

### DataNode Types

| Type | Description | Typical Databases |
|------|-------------|-------------------|
| `GeneProduct` | Genes and proteins | Ensembl, Entrez Gene, UniProt |
| `Metabolite` | Small molecules | ChEBI, HMDB, KEGG, PubChem |
| `Protein` | Specific proteins | UniProt |
| `Rna` | RNA molecules | Ensembl, miRBase |
| `Complex` | Protein complexes | Complex Portal |
| `Pathway` | Linked pathways | WikiPathways, Reactome |
| `Unknown` | Uncharacterized | - |

---

## Graphics Element

Defines visual properties of nodes and edges.

### Node Graphics

```xml
<Graphics CenterX="497.109" CenterY="390.192"
          Width="130.0" Height="20.0"
          ZOrder="32768"
          FillColor="ffcccc"
          FontSize="12"
          FontWeight="Bold"
          Valign="Middle"
          ShapeType="RoundedRectangle" />
```

| Attribute | Type | Description |
|-----------|------|-------------|
| `CenterX` | Float | X coordinate |
| `CenterY` | Float | Y coordinate |
| `Width` | Float | Node width |
| `Height` | Float | Node height |
| `ZOrder` | Integer | Z-index (stacking) |
| `FillColor` | Hex | Background color |
| `Color` | Hex | Text/border color |
| `FontSize` | Integer | Font size |
| `FontWeight` | String | "Normal", "Bold" |
| `Valign` | String | Vertical alignment |
| `ShapeType` | String | Shape type |

### Shape Types

- `Rectangle`
- `RoundedRectangle`
- `Oval`
- `Hexagon`
- `Pentagon`
- `Triangle`
- `Organelle` (mitochondria, ER, etc.)
- `Brace`

---

## Xref Element (Cross-Reference)

Links nodes to external databases.

```xml
<Xref Database="Ensembl" ID="ENSG00000113161" />
<Xref Database="ChEBI" ID="CHEBI:16113" />
<Xref Database="Entrez Gene" ID="3156" />
```

| Attribute | Type | Required | Description |
|-----------|------|----------|-------------|
| `Database` | String | Yes | Database name |
| `ID` | String | Yes | Identifier |

### Supported Databases

| Database | ID Format | Entity Type |
|----------|-----------|-------------|
| `Ensembl` | ENSG########### | Genes |
| `Entrez Gene` | ####### | Genes |
| `HGNC` | HGNC:#### | Genes |
| `UniProt` | [A-Z0-9]{6,10} | Proteins |
| `Uniprot-TrEMBL` | [A-Z0-9]{6,10} | Proteins |
| `ChEBI` | CHEBI:##### | Metabolites |
| `HMDB` | HMDB####### | Metabolites |
| `KEGG Compound` | C##### | Metabolites |
| `PubChem-compound` | ####### | Metabolites |
| `CAS` | ###-##-# | Metabolites |
| `DrugBank` | DB##### | Drugs |
| `WikiPathways` | WP#### | Pathways |
| `Reactome` | R-HSA-####### | Pathways |
| `Complex Portal` | CPX-#### | Complexes |
| `Wikidata` | Q####### | Any |

---

## Interaction Elements

Define relationships between nodes.

### Interaction Structure

```xml
<Interaction GraphId="a76ac">
  <Graphics ConnectorType="Elbow" ZOrder="12288" LineThickness="1.0">
    <Point X="125.275" Y="488.127" GraphRef="fce72" RelX="0.5" RelY="1.0" />
    <Point X="260.0" Y="568.811" GraphRef="e993b" RelX="-1.0" RelY="0.0"
           ArrowHead="mim-conversion" />
    <Anchor Position="0.8765" Shape="None" GraphId="cd6c7" />
  </Graphics>
  <Xref Database="" ID="" />
</Interaction>
```

| Attribute | Type | Description |
|-----------|------|-------------|
| `GraphId` | String | Unique identifier |
| `ConnectorType` | String | Line routing |
| `LineThickness` | Float | Line width |
| `LineStyle` | String | "Solid", "Dashed", "Double" |

### Connector Types

- `Straight` - Direct line
- `Elbow` - Right-angle turns
- `Curved` - Bezier curves
- `Segmented` - Multiple waypoints

### ArrowHead Types (MIM Notation)

| ArrowHead | Meaning |
|-----------|---------|
| `mim-conversion` | Metabolic conversion |
| `mim-stimulation` | Activation/stimulation |
| `mim-inhibition` | Inhibition |
| `mim-catalysis` | Enzymatic catalysis |
| `mim-binding` | Binding interaction |
| `mim-modification` | Post-translational modification |
| `mim-transcription-translation` | Gene expression |
| `Arrow` | Generic arrow |
| `TBar` | Inhibition bar |

---

## Point Elements

Define line endpoints and waypoints.

```xml
<Point X="300.916" Y="90.142" GraphRef="fb3" RelX="0.0" RelY="1.0" />
```

| Attribute | Type | Description |
|-----------|------|-------------|
| `X` | Float | X coordinate |
| `Y` | Float | Y coordinate |
| `GraphRef` | String | Connected node ID |
| `RelX` | Float | Relative X (-1 to 1) |
| `RelY` | Float | Relative Y (-1 to 1) |
| `ArrowHead` | String | Arrow style |

---

## Anchor Elements

Define connection points along edges.

```xml
<Anchor Position="0.495" Shape="None" GraphId="cb45c" />
```

| Attribute | Type | Description |
|-----------|------|-------------|
| `Position` | Float | 0.0-1.0 along line |
| `Shape` | String | "None", "Circle" |
| `GraphId` | String | Unique identifier |

---

## Label Elements

Free-text labels on the diagram.

```xml
<Label TextLabel="Endoplasmic Reticulum" GraphId="xyz789">
  <Graphics CenterX="400.0" CenterY="300.0" Width="150.0" Height="25.0"
            ZOrder="28672" FillColor="ffffff" FontSize="14" FontWeight="Bold" />
</Label>
```

---

## Group Elements

Visual grouping of related nodes.

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
```

| Attribute | Type | Description |
|-----------|------|-------------|
| `GroupId` | String | Group identifier |
| `Style` | String | Visual style |
| `TextLabel` | String | Group label |

### Group Styles

| Style | Description |
|-------|-------------|
| `None` | No visual grouping |
| `Group` | Simple grouping |
| `Complex` | Protein complex |
| `Pathway` | Nested pathway |

---

## State Elements (Post-Translational Modifications)

State elements represent PTMs or other states on DataNodes.

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

### Common PTM State Labels

| Label | PSI-MOD ID | Description |
|-------|------------|-------------|
| `P` | MOD:00046 | Phosphorylation |
| `Ub` | MOD:01148 | Ubiquitination |
| `Ac` | MOD:00394 | Acetylation |
| `Me` | MOD:00599 | Methylation |
| `SUMO` | MOD:01149 | SUMOylation |
| `Glyc` | MOD:00693 | Glycosylation |

### State ShapeTypes

| ShapeType | Description |
|-----------|-------------|
| `Oval` | Circle/oval (most common for PTMs) |
| `Rectangle` | Rectangular state |
| `RoundedRectangle` | Rounded rectangle |

---

## Shape Elements

Background shapes and annotations.

```xml
<Shape GraphId="membrane" TextLabel="Plasma membrane">
  <Graphics CenterX="400.0" CenterY="500.0" Width="800.0" Height="50.0"
            ZOrder="16384" FillColor="cccccc" ShapeType="Rectangle" />
</Shape>
```

---

## Comment Elements

Metadata and descriptions.

```xml
<Comment Source="WikiPathways-description">
  Cholesterol is a waxy steroid metabolite found in cell membranes...
</Comment>
<Comment Source="WikiPathways-category">Metabolic Process</Comment>
```

| Source | Description |
|--------|-------------|
| `WikiPathways-description` | Pathway description |
| `WikiPathways-category` | Pathway category |
| `GenMAPP remarks` | Legacy annotations |
| `HomologyMapper` | Ortholog mapping |

---

## BiopaxRef Elements

Literature citations (BioPAX format).

```xml
<BiopaxRef>ec9</BiopaxRef>
```

Links to Biopax elements at the end of the file containing publication details.

---

## Sample Complete Pathway

```xml
<?xml version="1.0" encoding="UTF-8"?>
<Pathway xmlns="http://pathvisio.org/GPML/2013a"
         Name="Cholesterol biosynthesis pathway"
         Organism="Homo sapiens">

  <Comment Source="WikiPathways-description">
    Cholesterol biosynthesis from Acetyl-CoA...
  </Comment>

  <Graphics BoardWidth="642.628" BoardHeight="648.152" />

  <DataNode TextLabel="HMGCR" GraphId="d8b" Type="GeneProduct">
    <Graphics CenterX="218.39" CenterY="190.95" Width="66.667" Height="20.0"
              ZOrder="32768" FontSize="12" Valign="Middle" />
    <Xref Database="Entrez Gene" ID="3156" />
  </DataNode>

  <DataNode TextLabel="Cholesterol" GraphId="e9f" Type="Metabolite">
    <Graphics CenterX="497.109" CenterY="80.142" Width="66.667" Height="20.0"
              ZOrder="32768" FontSize="12" Valign="Middle" Color="0000ff" />
    <Xref Database="CAS" ID="57-88-5" />
  </DataNode>

  <DataNode TextLabel="HMG-CoA" GraphId="e89f4" Type="Metabolite">
    <Graphics CenterX="300.916" CenterY="155.899" Width="66.667" Height="20.0"
              ZOrder="32768" FontSize="12" Valign="Middle" Color="0000ff" />
    <Xref Database="HMDB" ID="HMDB0001375" />
  </DataNode>

  <Interaction GraphId="a94b2">
    <Graphics ZOrder="12288" LineThickness="1.0">
      <Point X="300.916" Y="90.142" GraphRef="fb3" RelX="0.0" RelY="1.0" />
      <Point X="300.916" Y="145.899" GraphRef="e89f4" RelX="0.0" RelY="-1.0"
             ArrowHead="mim-conversion" />
      <Anchor Position="0.495" Shape="None" GraphId="cb45c" />
    </Graphics>
    <Xref Database="" ID="" />
  </Interaction>

</Pathway>
```

---

## Statistics (January 2026)

| Metric | Count |
|--------|-------|
| Total Pathways | 3,100+ |
| Human Pathways | 955+ |
| Mouse Pathways | 259 |
| Rat Pathways | 155 |
| Organisms | 48 |

### Supported Organisms

```
Homo sapiens, Mus musculus, Rattus norvegicus, Danio rerio,
Drosophila melanogaster, Caenorhabditis elegans, Saccharomyces cerevisiae,
Bos taurus, Gallus gallus, Arabidopsis thaliana, Zea mays,
Plasmodium falciparum, Mycobacterium tuberculosis, Escherichia coli,
and 34+ others
```

---

## API Endpoints

| Endpoint | Description |
|----------|-------------|
| `getPathway` | Get full GPML content |
| `getPathwayInfo` | Get pathway metadata |
| `findPathwaysByText` | Search pathways |
| `listOrganisms` | List all organisms |
| `listPathways` | List pathways by organism |
| `getCurationTags` | Get curation status |

**Base URL:** `https://webservice.wikipathways.org`

**Response Formats:** JSON, XML

---

## Curation Tags

| Tag | Description |
|-----|-------------|
| `Curation:AnalysisCollection` | Approved for analysis |
| `Curation:FeaturedPathway` | Highlighted quality |
| `Curation:NeedsWork` | Requires improvements |
| `Curation:Stub` | Incomplete pathway |
| `Curation:Reactive` | Derived from Reactome |

---

## BridgeDb Integration

WikiPathways uses BridgeDb for identifier mapping:

```
https://bridgedb.github.io/
```

Enables conversion between:
- Ensembl <-> Entrez Gene <-> HGNC
- ChEBI <-> HMDB <-> KEGG <-> PubChem
- UniProt <-> RefSeq <-> Gene Symbol

---

## BioPAX Section for Citations

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
  </bp:PublicationXref>

  <bp:openControlledVocabulary xmlns:bp="http://www.biopax.org/release/biopax-level3.owl#">
    <bp:TERM rdf:datatype="http://www.w3.org/2001/XMLSchema#string">apoptotic process</bp:TERM>
    <bp:ID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">GO:0006915</bp:ID>
    <bp:Ontology rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Gene Ontology</bp:Ontology>
  </bp:openControlledVocabulary>
</Biopax>
```

---

## Python Integration

```python
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional
from dataclasses import dataclass, field
import requests

@dataclass
class GPMLDataNode:
    graph_id: str
    text_label: str
    node_type: str
    database: Optional[str] = None
    identifier: Optional[str] = None
    center_x: float = 0.0
    center_y: float = 0.0
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

class GPMLParser:
    """Parser for WikiPathways GPML format."""

    NAMESPACE = {"gpml": "http://pathvisio.org/GPML/2013a"}

    def parse(self, gpml_content: str) -> Dict:
        """Parse GPML XML content."""
        root = ET.fromstring(gpml_content)
        ns = self.NAMESPACE

        pathway = {
            "name": root.get("Name"),
            "organism": root.get("Organism"),
            "data_nodes": [],
            "interactions": [],
            "groups": []
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
                center_x=float(graphics.get("CenterX", 0)) if graphics else 0,
                center_y=float(graphics.get("CenterY", 0)) if graphics else 0,
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
                database=xref.get("Database") if xref else None,
                identifier=xref.get("ID") if xref else None
            )
            # Attach to parent node
            for node in pathway["data_nodes"]:
                if node.graph_id == state_data.parent_ref:
                    node.states.append(state_data)
                    break

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

def fetch_wikipathways_gpml(pathway_id: str) -> str:
    """Fetch GPML from WikiPathways."""
    url = f"https://webservice.wikipathways.org/getPathway?pwId={pathway_id}&format=json"
    response = requests.get(url)
    data = response.json()
    return data.get("pathway", {}).get("gpml", "")

# Usage
gpml_content = fetch_wikipathways_gpml("WP254")
parser = GPMLParser()
pathway = parser.parse(gpml_content)
print(f"Pathway: {pathway['name']}")
print(f"DataNodes: {len(pathway['data_nodes'])}")
```
