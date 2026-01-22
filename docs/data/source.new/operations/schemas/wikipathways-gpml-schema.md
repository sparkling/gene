---
id: schemas-wikipathways-gpml-schema
title: "WikiPathways GPML Schema Documentation"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, wikipathways, gpml, pathway, xml, visualization, cc0]
---

**Parent:** [Schema Documentation](./_index.md)

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
| `ShapeType` | String | Shape (Rectangle, Oval, etc.) |

---

## Xref Element

Cross-references to external databases.

```xml
<Xref Database="Entrez Gene" ID="3156" />
```

### Xref Database Codes

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

---

## Interaction Elements

### Interaction Element (Biological Relationship)

```xml
<Interaction GraphId="id12345">
  <Comment>Caspase-9 activates Caspase-3</Comment>
  <Graphics ZOrder="12288" LineThickness="1.0" Color="000000">
    <Point X="350.0" Y="300.0" GraphRef="nodeA" RelX="1.0" RelY="0.0"/>
    <Point X="410.0" Y="300.0" GraphRef="nodeB" RelX="-1.0" RelY="0.0" ArrowHead="mim-stimulation"/>
  </Graphics>
  <Xref Database="" ID=""/>
</Interaction>
```

### ArrowHead Types (MIM Notation)

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
| `TBar` | T-bar | Inhibition (alias) |

---

## Group Element (Complexes)

```xml
<!-- Group definition -->
<Group GroupId="apoptosome" Style="Complex" TextLabel="Apoptosome">
  <Comment>Apoptosome complex formed during apoptosis</Comment>
</Group>

<!-- Group members -->
<DataNode TextLabel="APAF1" GraphId="apaf1" Type="GeneProduct" GroupRef="apoptosome">
  <Graphics CenterX="400.0" CenterY="400.0" Width="60.0" Height="20.0" ZOrder="32769"/>
  <Xref Database="Ensembl" ID="ENSG00000120868"/>
</DataNode>
```

### Group Styles

| Style | Description |
|-------|-------------|
| `None` | No visual grouping |
| `Group` | Simple grouping |
| `Complex` | Protein complex |
| `Pathway` | Nested pathway |

---

## State Elements for Modifications

State elements represent post-translational modifications or other states on DataNodes.

```xml
<State TextLabel="P" GraphId="state1" GraphRef="tp53node">
  <Comment>Phosphorylation at Ser15</Comment>
  <Graphics RelX="1.0" RelY="-1.0" Width="15.0" Height="15.0"
            ShapeType="Oval" FillColor="ff0000"/>
  <Xref Database="PSI-MOD" ID="MOD:00046"/>
</State>
```

---

## Shape Element (Compartments)

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

### ShapeTypes

| ShapeType | Description |
|-----------|-------------|
| `Rectangle` | Rectangle |
| `RoundedRectangle` | Rectangle with rounded corners |
| `Oval` | Ellipse/Circle |
| `Mitochondria` | Mitochondria shape |
| `Endoplasmic Reticulum` | ER shape |
| `Golgi Apparatus` | Golgi shape |
| `Nucleus` | Nucleus shape |
| `Cell` | Cell outline |

---

## API Access

### Get Pathway GPML

```bash
curl "https://webservice.wikipathways.org/getPathway?pwId=WP254&format=json"
```

### Search Pathways

```bash
curl "https://webservice.wikipathways.org/findPathwaysByText?query=apoptosis&species=Homo%20sapiens&format=json"
```

---

## License

**License:** CC0 1.0 Universal (Public Domain for pathway data)
**GPML Spec:** CC BY 4.0

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `GraphId` | Unique identifier for any element within a GPML document | abc123 |
| `DataNode` | Visual element representing a biological entity (gene, protein, metabolite) | HMGCR gene node |
| `Interaction` | Line connecting two DataNodes representing a biological relationship | Activation arrow |
| `Xref` | Cross-reference linking a node to an external database identifier | Entrez Gene:3156 |
| `TextLabel` | Display text shown on a DataNode or other element | "HMGCR" |
| `Type` | Classification of DataNode (GeneProduct, Metabolite, Protein) | GeneProduct |
| `ArrowHead` | Symbol at the end of an interaction line indicating relationship type | mim-stimulation |
| `Group` | Container element for grouping related DataNodes (e.g., complexes) | Apoptosome complex |
| `State` | Modification marker attached to a DataNode (e.g., phosphorylation) | P (phosphate) |
| `Shape` | Visual element representing cellular compartments or regions | Mitochondria shape |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| MIM Notation | Molecular Interaction Map symbols for standardized pathway diagrams | ArrowHead types |
| GeneProduct | DataNode type for genes and their protein products | DataNode Type |
| Metabolite | DataNode type for small molecules and metabolites | DataNode Type |
| Complex | Group style for multi-protein assemblies | Group Style |
| Post-translational Modification | Covalent changes to proteins (phosphorylation, etc.) | State elements |
| Compartment | Cellular region where reactions occur | Shape elements |
| Pathway Reference | Link to another pathway within WikiPathways | DataNode Type=Pathway |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GPML | Graphical Pathway Markup Language | WikiPathways XML format |
| MIM | Molecular Interaction Map | Diagram notation standard |
| WP | WikiPathways | Pathway identifier prefix |
| CC0 | Creative Commons Zero | Public domain license |
| HMDB | Human Metabolome Database | Metabolite database |
| ChEBI | Chemical Entities of Biological Interest | Metabolite database |
| GO | Gene Ontology | Compartment annotations |
| PTM | Post-Translational Modification | Protein modifications |
| PSI-MOD | Proteomics Standards Initiative Modification | PTM ontology |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming authority |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
