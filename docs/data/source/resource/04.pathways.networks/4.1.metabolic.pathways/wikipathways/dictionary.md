# WikiPathways GPML - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | wikipathways |
| **Name** | WikiPathways GPML |
| **Total Fields** | 40 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| Pathway.Name | String | Yes | Pathway title | Cholesterol biosynthesis pathway |
| Pathway.Organism | String | Yes | Species name | Homo sapiens |
| Pathway.Version | String | No | Version date (YYYYMMDD) | 20200316 |
| DataNode.TextLabel | String | Yes | Display text | HMGCR |
| DataNode.GraphId | String | Yes | Unique node identifier | abc123 |
| DataNode.Type | String | Yes | Node classification | GeneProduct, Metabolite |
| Interaction.GraphId | String | Yes | Unique interaction identifier | id12345 |
| Xref.Database | String | Yes | External database name | Entrez Gene |
| Xref.ID | String | Yes | External identifier | 3156 |
| Graphics.CenterX | Float | Yes | X coordinate | 218.39 |
| Graphics.CenterY | Float | Yes | Y coordinate | 190.95 |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| WikiPathways | WP[0-9]+ | Pathway identifier | WP254 |
| Ensembl | ENSG[0-9]{11} | Gene identifier | ENSG00000164305 |
| Entrez Gene | Numeric | NCBI Gene ID | 836 |
| HGNC | HGNC:[0-9]+ | Gene naming | HGNC:1504 |
| UniProt | [A-Z][0-9][A-Z0-9]{3}[0-9] | Protein accession | P42574 |
| ChEBI | CHEBI:[0-9]+ | Metabolite identifier | CHEBI:15377 |
| HMDB | HMDB[0-9]+ | Metabolome database | HMDB0000001 |
| PubChem | Numeric | Compound identifier | 2244 |
| Reactome | R-[A-Z]{3}-[0-9]+ | Pathway reference | R-HSA-109581 |

---

## Enumerations

### DataNode Types

| Value | Description | Typical Databases |
|-------|-------------|-------------------|
| GeneProduct | Genes and proteins | Ensembl, Entrez Gene, UniProt |
| Metabolite | Small molecules | ChEBI, HMDB, KEGG, PubChem |
| Protein | Specific proteins | UniProt |
| Rna | RNA molecules | Ensembl, miRBase |
| Complex | Protein complexes | Complex Portal |
| Pathway | Linked pathways | WikiPathways, Reactome |
| Unknown | Uncharacterized | - |

### ArrowHead Types (MIM Notation)

| Value | Symbol | Meaning |
|-------|--------|---------|
| Arrow | Simple arrow | General direction/flow |
| mim-catalysis | Circle on line | Catalysis |
| mim-stimulation | Open arrowhead | Stimulation/Activation |
| mim-inhibition | T-bar | Inhibition |
| mim-conversion | Filled arrowhead | Conversion |
| mim-modification | Filled arrowhead | Covalent modification |
| mim-binding | Two lines | Non-covalent binding |
| mim-cleavage | Diagonal line | Covalent cleavage |
| TBar | T-bar | Inhibition (alias) |

### Group Styles

| Value | Description |
|-------|-------------|
| None | No visual grouping |
| Group | Simple grouping |
| Complex | Protein complex |
| Pathway | Nested pathway |

### ShapeTypes

| Value | Description |
|-------|-------------|
| Rectangle | Rectangle |
| RoundedRectangle | Rounded corners |
| Oval | Ellipse/Circle |
| Mitochondria | Mitochondria shape |
| Endoplasmic Reticulum | ER shape |
| Golgi Apparatus | Golgi shape |
| Nucleus | Nucleus shape |
| Cell | Cell outline |

---

## Entity Relationships

### Pathway to DataNode
- **Cardinality:** 1:N
- **Description:** Pathways contain data nodes
- **Key Fields:** DataNode elements within Pathway

### DataNode to Xref
- **Cardinality:** 1:N
- **Description:** Nodes have external references
- **Key Fields:** DataNode.GraphId, Xref.Database, Xref.ID

### Interaction to DataNode
- **Cardinality:** N:M
- **Description:** Interactions connect nodes
- **Key Fields:** Point.GraphRef

### DataNode to Group
- **Cardinality:** N:1
- **Description:** Nodes can belong to groups
- **Key Fields:** DataNode.GroupRef, Group.GroupId

---

## Acronyms

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

## See Also

- [schema.md](./schema.md) - Full schema documentation
- [sample.json](./sample.json) - Example records
