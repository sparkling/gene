---
id: schemas-kgml-schema
title: "KEGG KGML Schema Documentation"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, kegg, kgml, pathway, xml, licensed]
---

**Parent:** [Schema Documentation](./_index.md)

# KEGG KGML Schema Documentation

**Source:** https://www.kegg.jp
**Format:** KGML (KEGG Markup Language)
**License:** Academic/Commercial (Requires License)
**Status:** Reference Only

---

## License Notice

KEGG requires licensing for commercial use. This schema documentation is provided for reference purposes only. For data access:
- **Academic users:** Limited free access via REST API
- **Commercial users:** Requires paid license from Kanehisa Laboratories
- **Alternative:** Use Reactome (CC BY 4.0) for similar pathway data

---

## Overview

KGML (KEGG Markup Language) is the XML format used to represent KEGG pathway maps with graphical and relational information.

---

## DTD Schema Definition

```xml
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

<!ELEMENT relation (subtype*)>
<!ATTLIST relation
    entry1      IDREF   #REQUIRED
    entry2      IDREF   #REQUIRED
    type        (ECrel|PPrel|GErel|PCrel|maplink) #REQUIRED
>

<!ELEMENT reaction (substrate*, product*)>
<!ATTLIST reaction
    id          IDREF   #REQUIRED
    name        CDATA   #REQUIRED
    type        (reversible|irreversible) #REQUIRED
>
```

---

## Entry Types

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

---

## Relation Types

### ECrel (Enzyme-Enzyme Relation)

Represents metabolic relationships through shared compounds.

| Subtype | Value | Description |
|---------|-------|-------------|
| `compound` | Entry ID | Shared compound between enzymes |

### PPrel (Protein-Protein Relation)

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

### GErel (Gene Expression Relation)

| Subtype | Value | Description |
|---------|-------|-------------|
| `expression` | `-->` | Gene expression activation |
| `repression` | `--\|` | Gene expression repression |

### PCrel (Protein-Compound Relation)

| Subtype | Value | Description |
|---------|-------|-------------|
| `activation` | `-->` | Compound activates protein |
| `inhibition` | `--\|` | Compound inhibits protein |

---

## Sample KGML

```xml
<?xml version="1.0"?>
<!DOCTYPE pathway SYSTEM "https://www.kegg.jp/kegg/xml/KGML_v0.7.2_.dtd">
<pathway name="path:hsa04115" org="hsa" number="04115"
         title="p53 signaling pathway"
         image="https://www.kegg.jp/kegg/pathway/hsa/hsa04115.png"
         link="https://www.kegg.jp/kegg-bin/show_pathway?hsa04115">

  <!-- Gene entry -->
  <entry id="1" name="hsa:7157" type="gene"
         link="https://www.kegg.jp/dbget-bin/www_bget?hsa:7157">
    <graphics name="TP53" type="rectangle"
              x="547" y="308" width="46" height="17"
              fgcolor="#000000" bgcolor="#BFFFBF"/>
  </entry>

  <!-- Compound entry -->
  <entry id="50" name="cpd:C00027" type="compound">
    <graphics name="C00027" type="circle"
              x="200" y="400" width="8" height="8"/>
  </entry>

  <!-- Pathway map reference -->
  <entry id="100" name="path:hsa04110" type="map">
    <graphics name="Cell cycle" type="roundrectangle"
              x="716" y="171" width="100" height="34"/>
  </entry>

  <!-- Protein-protein relation with phosphorylation -->
  <relation entry1="3" entry2="1" type="PPrel">
    <subtype name="activation" value="-->"/>
    <subtype name="phosphorylation" value="+p"/>
  </relation>

  <!-- Gene expression relation -->
  <relation entry1="1" entry2="2" type="GErel">
    <subtype name="expression" value="-->"/>
  </relation>

  <!-- Metabolic reaction -->
  <reaction id="300" name="rn:R00001" type="irreversible">
    <substrate id="50" name="cpd:C00027"/>
    <product id="51" name="cpd:C00001"/>
  </reaction>

</pathway>
```

---

## Graphics Element

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | String | Display label |
| `x` | Integer | X coordinate |
| `y` | Integer | Y coordinate |
| `type` | Enum | rectangle, circle, roundrectangle, line |
| `width` | Integer | Node width |
| `height` | Integer | Node height |
| `fgcolor` | Hex | Foreground color |
| `bgcolor` | Hex | Background color |

---

## Cross-References

| Database | ID Format | Example |
|----------|-----------|---------|
| KEGG Gene | org:id | hsa:7157 |
| KEGG Compound | C##### | C00001 |
| KEGG Orthology | K##### | K00001 |
| KEGG Reaction | R##### | R00001 |
| EC Number | #.#.#.# | 1.1.1.1 |

---

## Open Alternatives

Due to KEGG licensing restrictions, consider these open alternatives:

| Database | License | Coverage |
|----------|---------|----------|
| **Reactome** | CC BY 4.0 | 2,712 human pathways |
| **WikiPathways** | CC0 | 3,100+ pathways |
| **Pathway Commons** | Open | Aggregated pathways |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `pathway` | Root XML element representing a KEGG pathway map | path:hsa04115 |
| `entry` | Node element representing a biological entity in the pathway | Gene, compound, map |
| `relation` | Edge element representing interaction between two entries | PPrel, GErel |
| `reaction` | Metabolic reaction with substrate and product entries | rn:R00001 |
| `graphics` | Visual rendering properties for pathway elements | x, y coordinates |
| `subtype` | Specific interaction type within a relation | activation, phosphorylation |
| `name` | KEGG identifier for an entry (organism:id or compound format) | hsa:7157, cpd:C00001 |
| `type` | Classification of entry (gene, compound, map) or relation (PPrel, GErel) | gene, PPrel |
| `id` | Unique identifier within a KGML document | 1, 50, 100 |
| `org` | Three-letter organism code for the pathway | hsa (human) |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| KEGG Orthology (KO) | Functional orthologs grouped across species | K##### entries |
| EC Number | Enzyme Commission classification for enzyme function | ec:1.1.1.1 |
| PPrel | Protein-Protein relation for signaling interactions | relation type |
| GErel | Gene Expression relation for transcriptional regulation | relation type |
| ECrel | Enzyme-Enzyme relation through shared metabolites | relation type |
| PCrel | Protein-Compound relation for compound effects | relation type |
| maplink | Reference to another pathway within current pathway | relation type |
| BRITE | KEGG hierarchical classification system | br:hsa00001 |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| KGML | KEGG Markup Language | XML pathway format |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |
| KO | KEGG Orthology | Functional classification |
| EC | Enzyme Commission | Enzyme classification |
| DTD | Document Type Definition | XML schema format |
| PPrel | Protein-Protein Relation | Signaling interaction |
| GErel | Gene Expression Relation | Transcriptional regulation |
| ECrel | Enzyme-Enzyme Relation | Metabolic connection |
| PCrel | Protein-Compound Relation | Small molecule effect |
| BRITE | Biomolecular Relations in Information Transmission and Expression | KEGG hierarchy |

---

## References

- KEGG KGML Documentation: https://www.kegg.jp/kegg/xml/docs/
- KEGG API: https://www.kegg.jp/kegg/rest/keggapi.html
- Kanehisa M, Goto S. (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28(1):27-30.
