# KEGG KGML - Data Dictionary

## Overview

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | kegg |
| **Name** | KEGG KGML |
| **Total Fields** | 35 |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

| Field Name | Data Type | Required | Description | Examples |
|------------|-----------|----------|-------------|----------|
| pathway.name | String | Yes | KEGG pathway identifier | path:hsa04115 |
| pathway.org | String | Yes | Three-letter organism code | hsa |
| pathway.number | String | Yes | Pathway number | 04115 |
| pathway.title | String | No | Human-readable pathway name | p53 signaling pathway |
| entry.id | ID | Yes | Unique entry identifier within KGML | 1, 50, 100 |
| entry.name | String | Yes | KEGG identifier (organism:id or compound) | hsa:7157, cpd:C00001 |
| entry.type | Enum | Yes | Entry classification | ortholog, enzyme, gene, compound, map |
| relation.entry1 | IDREF | Yes | Source entry ID | 3 |
| relation.entry2 | IDREF | Yes | Target entry ID | 1 |
| relation.type | Enum | Yes | Relationship type | ECrel, PPrel, GErel, PCrel, maplink |
| reaction.id | IDREF | Yes | Reaction entry reference | 300 |
| reaction.name | String | Yes | KEGG reaction identifier | rn:R00001 |
| reaction.type | Enum | Yes | Reversibility | reversible, irreversible |

---

## Identifiers

| ID Type | Format | Description | Examples |
|---------|--------|-------------|----------|
| KEGG Pathway | path:{org}{number} | Pathway identifier | path:hsa00010 |
| KEGG Gene | {org}:{id} | Organism-specific gene | hsa:7157 |
| KEGG Compound | C##### | Compound identifier | C00001 |
| KEGG Orthology | K##### | KO group identifier | K00001 |
| KEGG Reaction | R##### | Reaction identifier | R00001 |
| EC Number | #.#.#.# | Enzyme classification | 1.1.1.1 |
| BRITE | br:{org}##### | Hierarchy identifier | br:hsa00001 |

---

## Enumerations

### Entry Types

| Value | Description |
|-------|-------------|
| ortholog | KEGG Orthology (KO) group |
| enzyme | Enzyme by EC number |
| gene | Gene product |
| group | Complex or group |
| compound | Chemical compound |
| map | Linked pathway |
| brite | BRITE hierarchy |
| other | Other types |

### Relation Types

| Value | Description |
|-------|-------------|
| ECrel | Enzyme-Enzyme relation (metabolic) |
| PPrel | Protein-Protein relation (signaling) |
| GErel | Gene Expression relation (transcriptional) |
| PCrel | Protein-Compound relation |
| maplink | Reference to another pathway |

### PPrel Subtypes

| Value | Symbol | Description |
|-------|--------|-------------|
| activation | --> | Activation |
| inhibition | --\| | Inhibition |
| phosphorylation | +p | Phosphorylation |
| dephosphorylation | -p | Dephosphorylation |
| ubiquitination | +u | Ubiquitination |
| methylation | +m | Methylation |

---

## Entity Relationships

### Pathway to Entry
- **Cardinality:** 1:N
- **Description:** A pathway contains multiple entries
- **Key Fields:** pathway.name, entry.id

### Entry to Relation
- **Cardinality:** N:M
- **Description:** Entries connected via relations
- **Key Fields:** entry.id, relation.entry1, relation.entry2

### Entry to Reaction
- **Cardinality:** N:M
- **Description:** Entries participate in reactions
- **Key Fields:** entry.id, reaction.substrate, reaction.product

---

## Acronyms

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

## See Also

- [schema.md](./schema.md) - Full schema documentation
- [sample.json](./sample.json) - Example records
