---
id: guides-integration-pathway-target-mapping
title: "Pathway and Target Database Data Models"
type: guide
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [integration, pathways, targets, reactome, kegg, wikipathways, uniprot, string, guide]
---

**Parent:** [Integration Guides](./_index.md)


This document provides detailed data models, schemas, and field descriptions for pathway and target databases used to link compounds to genes.

## Table of Contents

1. [Reactome](#1-reactome)
2. [KEGG](#2-kegg)
3. [WikiPathways](#3-wikipathways)
4. [UniProt](#4-uniprot)
5. [PubChem](#5-pubchem)
6. [STRING/STITCH](#6-stringstitch)
7. [Cross-Database ID Mappings](#7-cross-database-id-mappings)

---

## 1. Reactome

**Base URL:** `https://reactome.org/ContentService/`

**Documentation:** [Reactome Content Service](https://reactome.org/dev/content-service)

### 1.1 API Overview

Reactome provides a REST API (version 1.2) following OpenAPI 3.0.1 standards. The API is organized into 16 functional endpoint categories.

### 1.2 Endpoint Categories

| Category | Description |
|----------|-------------|
| Events | Pathway and reaction queries |
| Entities | Physical entity operations (complexes, proteins, chemicals) |
| Pathways | Hierarchical pathway traversal |
| Participants | Event participant retrieval |
| Mapping | Cross-resource identifier mapping |
| Species | Organism-specific data filtering |
| Orthology | Computationally inferred events across species |
| Schema | Data model and class queries |
| Query | General object retrieval and search |
| Database | System metadata (version, name) |
| Person | Author and curator information |
| References | External database cross-references |
| Interactors | Protein-protein interaction data |
| Exporter | Format conversion (PDF, SBML, SBGN, images) |
| Discover | Schema.org dataset markup |
| Diseases | Disease annotation queries |

### 1.3 Key Endpoints

#### Pathway Hierarchy
```
GET /data/eventsHierarchy/{species}     # Full event trees
GET /data/pathways/top/{species}        # Top-level pathways
GET /data/pathway/{id}/containedEvents  # Nested event retrieval
```

#### Entity Operations
```
GET /data/complex/{id}/subunits         # Complex decomposition
GET /data/entity/{id}/otherForms        # Entity variant forms
GET /data/entity/{id}/componentOf       # Structural containment
GET /data/query/{id}                    # Complete entity information
GET /data/query/{id}/{attribute}        # Specific entity attributes
```

#### Cross-Mapping
```
GET /data/mapping/{resource}/{identifier}/pathways  # ID to pathway mapping
GET /data/orthology/{id}/species/{speciesId}        # Orthologous events
```

#### Export
```
GET /exporter/diagram/{identifier}.{ext}        # PNG, SVG, GIF, JPG
GET /exporter/document/event/{identifier}.pdf   # Full pathway documents
```

### 1.4 Data Schema

#### Core DatabaseObject (Base Class)

All Reactome entities inherit from `DatabaseObject`:

```json
{
  "dbId": 109581,
  "stId": "R-HSA-109581",
  "displayName": "Apoptosis",
  "schemaClass": "Pathway",
  "speciesName": "Homo sapiens",
  "created": { "instanceEdit": {} },
  "modified": { "instanceEdit": {} }
}
```

| Field | Type | Description |
|-------|------|-------------|
| `dbId` | Long | Unique database identifier |
| `stId` | String | Stable identifier (format: R-{species}-{number}) |
| `displayName` | String | Human-readable label |
| `schemaClass` | String | Entity type classification |
| `speciesName` | String | Organism designation |
| `created` | InstanceEdit | Creation timestamp |
| `modified` | InstanceEdit | Last modification timestamp |

#### Primary Entity Classes

| Class | Count | Description |
|-------|-------|-------------|
| PhysicalEntity | 406,613 | Biological molecules (Complexes, Drugs, Polymers) |
| Event | 117,945 | Pathways (23,290) and ReactionLikeEvents (94,655) |
| ReferenceEntity | 930,514 | Reference sequences and molecular definitions |
| Person | 171,942 | Curator and contributor information |
| Publication | 42,286 | Literature references |

#### PhysicalEntity Subtypes

- **Complex** - Multi-component molecular assemblies
- **EntityWithAccessionedSequence (EWAS)** - Proteins with UniProt IDs
- **SimpleEntity** - Small molecules (ChEBI references)
- **Drug** - Therapeutic compounds
- **Polymer** - Polymeric molecules
- **EntitySet** - Collections (DefinedSet, OpenSet, CandidateSet)

#### Event Subtypes

- **Pathway** - Biological processes with hierarchical structure
- **Reaction** - Biochemical transformations
- **BlackBoxEvent** - Events with unknown mechanism
- **FailedReaction** - Documented unsuccessful reactions

### 1.5 Pathway Response Structure

```json
{
  "dbId": 109581,
  "stId": "R-HSA-109581",
  "displayName": "Apoptosis",
  "speciesName": "Homo sapiens",
  "releaseDate": "2004-09-20",
  "hasEvent": [
    {
      "stId": "R-HSA-5357769",
      "displayName": "Caspase activation via extrinsic apoptotic signalling pathway"
    },
    {
      "stId": "R-HSA-109606",
      "displayName": "Intrinsic Pathway for Apoptosis"
    },
    {
      "stId": "R-HSA-75153",
      "displayName": "Apoptotic execution phase"
    },
    {
      "stId": "R-HSA-169911",
      "displayName": "Regulation of Apoptosis"
    }
  ],
  "summation": [
    {
      "text": "Programmed cell death which begins when a cell receives an internal or external signal..."
    }
  ],
  "literatureReference": [
    {
      "pubMedIdentifier": 12505355
    }
  ]
}
```

### 1.6 Graph Database Relationships

Reactome uses Neo4j with these key relationship types:

| Relationship | Description |
|--------------|-------------|
| `hasComponent` | Complex members |
| `hasMember` | Set participants |
| `hasCandidate` | Alternative set members |
| `referenceEntity` | Links to external IDs |
| `input` | Reaction inputs |
| `output` | Reaction outputs |
| `catalystActivity` | Enzymatic roles |
| `hasEvent` | Pathway sub-events |
| `regulatedBy` | Regulatory relationships |

### 1.7 Cross-References

Reactome integrates with:
- UniProt (protein sequences)
- ChEBI (small molecules)
- Ensembl (genes)
- NCBI Gene
- KEGG
- PubMed
- GO (Gene Ontology)

### 1.8 Response Formats

- **JSON** (Accept: application/json) - Structured data
- **TSV** (Accept: text/plain) - ID | displayName | schemaClass
- **PDF, PNG, SVG** - Diagram exports

---

## 2. KEGG

**Base URL:** `https://rest.kegg.jp/`

**Documentation:** [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html)

### 2.1 API Operations

| Operation | Description | Example |
|-----------|-------------|---------|
| `info` | Database release statistics | `/info/kegg` |
| `list` | Entry identifiers and names | `/list/pathway/hsa` |
| `find` | Keyword/property search | `/find/compound/C7H10O5/formula` |
| `get` | Retrieve database entries | `/get/hsa:10458/aaseq` |
| `conv` | Convert between KEGG and external IDs | `/conv/eco/ncbi-geneid` |
| `link` | Cross-database references | `/link/pathway/hsa` |
| `ddi` | Drug-drug interactions | `/ddi/D00564` |

### 2.2 Identifier Systems

#### KEGG Database Codes

| Code | Full Name | Description |
|------|-----------|-------------|
| `pathway` / `path` | Pathway | Biological pathways |
| `brite` / `br` | BRITE | Functional hierarchies |
| `module` / `md` | Module | Functional units |
| `orthology` / `ko` | Orthology | Cross-species gene groups |
| `genes` | Genes | Organism-specific genes |
| `genome` / `gn` | Genome | Complete genomes |
| `compound` / `cpd` | Compound | Small molecules |
| `glycan` / `gl` | Glycan | Carbohydrates |
| `reaction` / `rn` | Reaction | Biochemical reactions |
| `enzyme` / `ec` | Enzyme | EC classifications |
| `disease` / `ds` | Disease | Human diseases |
| `drug` / `dr` | Drug | Therapeutic drugs |

#### External Database Codes

| Code | Database |
|------|----------|
| `pubmed` | PubMed |
| `ncbi-geneid` | NCBI Gene |
| `ncbi-proteinid` | NCBI Protein |
| `uniprot` | UniProt |
| `pubchem` | PubChem |
| `chebi` | ChEBI |

#### ID Format Patterns

```
# KEGG pathway: <org prefix><5 digits>
hsa00010          # Human glycolysis pathway

# KEGG gene: <org>:<gene_id>
hsa:7157          # Human TP53 gene

# KEGG compound: C<5 digits>
C00027            # Hydrogen peroxide

# KEGG drug: D<5 digits>
D00564            # Aspirin

# K number (orthology): K<5 digits>
K04451            # p53 orthology group
```

### 2.3 Pathway Entry Structure

#### Pathway Prefix Classification

| Prefix | Description |
|--------|-------------|
| `map` | Manually drawn reference pathway |
| `ko` | Reference pathway highlighting KOs |
| `ec` | Metabolic pathway highlighting EC numbers |
| `rn` | Metabolic pathway highlighting reactions |
| `<org>` | Organism-specific pathway (e.g., hsa) |
| `vg/vx` | Virus-related pathways |

#### Pathway Number Classification

| Range | Type |
|-------|------|
| 011xx | Global maps linked to KOs |
| 012xx | Overview maps linked to KOs |
| 013xx | Multi-organism overview maps |
| 010xx | Chemical structure maps |
| 07xxx | Drug structure maps |
| Other | Regular maps linked to KOs |

### 2.4 Gene Entry Format

Example: `/get/hsa:7157` (TP53)

```
ENTRY       7157              CDS       T01001
SYMBOL      TP53, P53, TRP53
NAME        tumor protein p53
ORTHOLOGY   K04451  tumor protein p53
PATHWAY     hsa04010  MAPK signaling pathway
            hsa04110  Cell cycle
            hsa04115  p53 signaling pathway
            hsa05200  Pathways in cancer
            ...
NETWORK     nt06210  ERK signaling
            nt06214  PI3K signaling
DISEASE     H00004  Chronic myeloid leukemia
            H00005  Chronic lymphocytic leukemia
            H00006  Hairy-cell leukemia
            ...
DRUG_TARGET D10761  Cenersen sodium
            D12742  Rezatapopt
BRITE       Transcription factors [BR:hsa03000]
             Class: p53
POSITION    17:complement(7661779..7687538)
MOTIF       Pfam: P53 P53_tetramer P53_TAD
DBLINKS     NCBI-GeneID: 7157
            NCBI-ProteinID: NP_000537
            UniProt: P04637
            Ensembl: ENSG00000141510
AASEQ       393
            MEEPQSDPSV...
NTSEQ       1182
            atggaggagc...
```

| Field | Description |
|-------|-------------|
| ENTRY | Identifier, type (CDS), organism |
| SYMBOL | Gene symbols |
| NAME | Official gene name |
| ORTHOLOGY | K number and description |
| PATHWAY | Associated pathways |
| NETWORK | Disease-specific networks |
| DISEASE | Associated diseases |
| DRUG_TARGET | Therapeutic compounds |
| BRITE | Hierarchical classifications |
| POSITION | Chromosomal location |
| MOTIF | Protein domains |
| DBLINKS | External database cross-references |
| AASEQ | Amino acid sequence |
| NTSEQ | Nucleotide sequence |

### 2.5 Compound Entry Format

Example: `/get/cpd:C00027` (Hydrogen peroxide)

```
ENTRY       C00027                      Compound
NAME        Hydrogen peroxide;
            Oxydol;
            H2O2
FORMULA     H2O2
EXACT_MASS  34.0055
MOL_WEIGHT  34.01
REACTION    R00009 R00011 R00017 ...
PATHWAY     map00630  Glyoxylate and dicarboxylate metabolism
            map04210  Apoptosis
            map04216  Ferroptosis
            ...
ENZYME      1.1.3.4      1.1.3.9      1.1.3.13     ...
BRITE       Compounds with biological roles [BR:br08001]
DBLINKS     CAS: 7722-84-1
            PubChem: 3329
            ChEBI: 16240
ATOM        2
            1   O1a O    22.5385  -15.9592
            2   O1a O    23.7509  -15.2592
BOND        1
            1     1   2 1
```

### 2.6 KGML (KEGG Markup Language)

KGML is the XML format for pathway diagrams.

#### Structure Elements

```xml
<?xml version="1.0"?>
<pathway name="path:hsa05210" org="hsa" number="05210"
         title="Colorectal cancer">

  <!-- Gene Entry -->
  <entry id="39" name="hsa:7157" type="gene"
         link="https://www.kegg.jp/dbget-bin/www_bget?hsa:7157">
    <graphics name="TP53..." type="rectangle"
              x="786" y="516" width="46" height="17"/>
  </entry>

  <!-- Pathway Map Reference -->
  <entry id="76" name="path:hsa04310" type="map"
         link="https://www.kegg.jp/dbget-bin/www_bget?hsa04310">
    <graphics name="Wnt signaling pathway" type="roundrectangle"
              x="230" y="240" width="92" height="34"/>
  </entry>

  <!-- Relationships -->
  <relation entry1="33" entry2="34" type="PPrel">
    <subtype name="activation" value="-->"/>
    <subtype name="phosphorylation" value="+p"/>
  </relation>

  <relation entry1="39" entry2="40" type="GErel">
    <subtype name="expression" value="-->"/>
  </relation>

</pathway>
```

#### Relation Types

| Type | Description |
|------|-------------|
| `PPrel` | Protein-protein relation |
| `GErel` | Gene expression relation |
| `PCrel` | Protein-compound relation |
| `ECrel` | Enzyme-compound relation |
| `maplink` | Link to another map |

#### Subtype Values

| Subtype | Symbol | Description |
|---------|--------|-------------|
| `activation` | `-->` | Activation |
| `inhibition` | `--|` | Inhibition |
| `expression` | `-->` | Gene expression |
| `repression` | `--|` | Gene repression |
| `phosphorylation` | `+p` | Phosphorylation |
| `dephosphorylation` | `-p` | Dephosphorylation |
| `binding/association` | `---` | Binding |
| `dissociation` | `-+-` | Dissociation |

### 2.7 Link Response Format

Example: `/link/pathway/cpd:C00027`

```
cpd:C00027	path:map00630
cpd:C00027	path:map01110
cpd:C00027	path:map04022
cpd:C00027	path:map04210
...
```

Tab-separated format: `source_id<TAB>target_id`

### 2.8 Conversion Response Format

Example: `/conv/uniprot/hsa:7157`

```
hsa:7157	up:P04637
hsa:7157	up:K7PPA8
hsa:7157	up:Q53GA5
```

One KEGG gene can map to multiple UniProt IDs (isoforms).

---

## 3. WikiPathways

**Base URL:** `https://webservice.wikipathways.org/`

**API Specification:** `https://www.wikipathways.org/openapi/openapi.json`

### 3.1 Response Formats

All endpoints support multiple output formats via `format` parameter:
- JSON
- XML (default)
- HTML
- JPG
- PDF

### 3.2 Core Endpoints

#### Discovery

```
GET /listOrganisms              # Available organisms
GET /listPathways               # All pathways (optional organism filter)
GET /findPathwaysByText         # Keyword search
GET /findPathwaysByXref         # External reference search
GET /findPathwaysByLiterature   # PubMed ID, author, or title search
GET /findInteractions           # Entity interaction discovery
```

#### Pathway Access

```
GET /getPathway?pwId={id}       # Full pathway with GPML
GET /getPathwayInfo?pwId={id}   # Metadata only
GET /getPathwayAs?pwId={id}&fileType={type}  # Format conversion
GET /getColoredPathway          # Visualization (SVG, PDF, PNG)
GET /getXrefList?pwId={id}&code={db}  # Cross-references
```

#### History

```
GET /getPathwayHistory?pwId={id}&timestamp={ts}
GET /getRecentChanges?timestamp={ts}  # 2-month retention
```

#### Curation

```
GET /getCurationTags?pwId={id}
GET /getCurationTagsByName?tagName={name}
GET /getCurationTagHistory?pwId={id}
```

#### Ontology

```
GET /getOntologyTermsByPathway?pwId={id}
GET /getPathwaysByOntologyTerm?term={term}
GET /getPathwaysByParentOntologyTerm?term={term}
```

### 3.3 Pathway Response Structure

Example: `/getPathway?pwId=WP254&format=json`

```json
{
  "id": "WP254",
  "name": "Apoptosis",
  "organism": "Homo sapiens",
  "revision": "140926",
  "url": "https://www.wikipathways.org/pathways/WP254.html",
  "gpml": "<?xml version=\"1.0\" encoding=\"UTF-8\"?>..."
}
```

### 3.4 GPML (GenMAPP Pathway Markup Language)

GPML is the XML format for WikiPathways pathway models.

#### Root Element

```xml
<?xml version="1.0" encoding="UTF-8"?>
<Pathway xmlns="http://pathvisio.org/GPML/2021"
         Name="Apoptosis"
         Organism="Homo sapiens"
         Version="20231115">
  ...
</Pathway>
```

#### DataNode Element (Genes/Proteins)

```xml
<DataNode TextLabel="CASP3" GraphId="a1b2c" Type="GeneProduct">
  <Comment>Caspase-3, executioner caspase</Comment>
  <Graphics CenterX="450.0" CenterY="300.0"
            Width="80.0" Height="20.0" ZOrder="32768"/>
  <Xref Database="Ensembl" ID="ENSG00000164305"/>
  <Xref Database="Entrez Gene" ID="836"/>
  <Xref Database="HGNC" ID="1504"/>
</DataNode>
```

#### DataNode Types

| Type | Description |
|------|-------------|
| `GeneProduct` | Gene or protein |
| `Metabolite` | Small molecule |
| `Pathway` | Reference to another pathway |
| `Protein` | Protein entity |
| `Rna` | RNA molecule |
| `Complex` | Protein complex |

#### Interaction Element

```xml
<Interaction GraphId="id123">
  <Graphics ZOrder="12288" LineThickness="1.0">
    <Point X="450.0" Y="300.0" GraphRef="a1b2c" RelX="0.0" RelY="1.0"/>
    <Point X="450.0" Y="400.0" GraphRef="d4e5f" RelX="0.0" RelY="-1.0" ArrowHead="Arrow"/>
  </Graphics>
  <Xref Database="" ID=""/>
</Interaction>
```

#### Interaction ArrowHead Types

| ArrowHead | Meaning |
|-----------|---------|
| `Arrow` | Activation/stimulation |
| `TBar` | Inhibition |
| `mim-catalysis` | Catalysis |
| `mim-inhibition` | Inhibition |
| `mim-stimulation` | Stimulation |
| `mim-conversion` | Conversion |
| `mim-binding` | Binding |

#### Group Element (Complexes)

```xml
<Group GroupId="group1" Style="Complex">
  <Comment>Apoptosome complex</Comment>
</Group>

<DataNode TextLabel="APAF1" GroupRef="group1" .../>
<DataNode TextLabel="CYCS" GroupRef="group1" .../>
<DataNode TextLabel="CASP9" GroupRef="group1" .../>
```

#### Label Element

```xml
<Label TextLabel="Mitochondria" GraphId="label1">
  <Graphics CenterX="200.0" CenterY="100.0"
            Width="100.0" Height="20.0" ZOrder="28672"/>
</Label>
```

#### Shape Element (Compartments)

```xml
<Shape GraphId="shape1" Type="RoundedRectangle">
  <Graphics CenterX="300.0" CenterY="300.0"
            Width="200.0" Height="300.0" ZOrder="16384"
            Color="cccccc" FillColor="ffffff"/>
</Shape>
```

### 3.5 Cross-Reference Databases

WikiPathways uses these database codes:

| Code | Database |
|------|----------|
| `En` | Ensembl |
| `L` | Entrez Gene |
| `H` | HGNC |
| `S` | UniProt |
| `Ce` | ChEBI |
| `Ck` | KEGG Compound |
| `Wp` | WikiPathways |
| `Re` | Reactome |

---

## 4. UniProt

**Base URL:** `https://rest.uniprot.org/`

### 4.1 Endpoints

```
GET /uniprotkb/{accession}                    # Single entry
GET /uniprotkb/{accession}.json               # JSON format
GET /uniprotkb/search?query={query}           # Search
GET /uniprotkb/stream?query={query}           # Large result streaming
GET /uniref/{id}                              # UniRef clusters
GET /uniparc/{id}                             # UniParc archives
GET /proteomes/{id}                           # Proteomes
```

### 4.2 Query Parameters

| Parameter | Description |
|-----------|-------------|
| `query` | Search query (field:value syntax) |
| `fields` | Comma-separated return fields |
| `format` | Response format (json, tsv, fasta, xml) |
| `size` | Results per page |
| `cursor` | Pagination cursor |

### 4.3 Query Field Examples

```
accession:P04637                      # By accession
gene:TP53                             # By gene name
organism_id:9606                      # By organism (human)
keyword:Tumor suppressor              # By keyword
go:0006915                            # By GO term (apoptosis)
pathway:reactome                      # Has Reactome annotation
xref:kegg                             # Has KEGG cross-reference
```

### 4.4 Protein Entry JSON Schema

Example: `/uniprotkb/P04637.json`

```json
{
  "entryType": "UniProtKB reviewed (Swiss-Prot)",
  "primaryAccession": "P04637",
  "secondaryAccessions": ["Q15086", "Q16535", ...],
  "uniProtkbId": "P53_HUMAN",
  "entryAudit": {
    "firstPublicDate": "1987-08-13",
    "lastAnnotationUpdateDate": "2024-01-24",
    "lastSequenceUpdateDate": "1990-01-01",
    "entryVersion": 312,
    "sequenceVersion": 4
  },
  "annotationScore": 5.0,
  "organism": {
    "scientificName": "Homo sapiens",
    "commonName": "Human",
    "taxonId": 9606,
    "lineage": ["Eukaryota", "Metazoa", "Chordata", ...]
  },
  "proteinExistence": "1: Evidence at protein level",
  "proteinDescription": {
    "recommendedName": {
      "fullName": {
        "value": "Cellular tumor antigen p53"
      }
    },
    "alternativeNames": [...]
  },
  "genes": [
    {
      "geneName": {
        "value": "TP53"
      },
      "synonyms": [
        {"value": "P53"}
      ]
    }
  ],
  "comments": [...],
  "features": [...],
  "keywords": [...],
  "references": [...],
  "uniProtKBCrossReferences": [...],
  "sequence": {
    "value": "MEEPQSDPSV...",
    "length": 393,
    "molWeight": 43653,
    "crc64": "AD5C149FD8106131"
  }
}
```

### 4.5 Key JSON Fields

#### Comments Array

```json
{
  "comments": [
    {
      "commentType": "FUNCTION",
      "texts": [
        {
          "value": "Acts as a tumor suppressor in many tumor types..."
        }
      ]
    },
    {
      "commentType": "SUBUNIT",
      "texts": [
        {
          "value": "Binds DNA as a homotetramer..."
        }
      ]
    },
    {
      "commentType": "SUBCELLULAR LOCATION",
      "subcellularLocations": [
        {
          "location": {
            "value": "Nucleus"
          }
        },
        {
          "location": {
            "value": "Cytoplasm"
          }
        }
      ]
    },
    {
      "commentType": "DISEASE",
      "disease": {
        "diseaseId": "Li-Fraumeni syndrome",
        "diseaseAccession": "DI-00652",
        "acronym": "LFS",
        "description": "An autosomal dominant..."
      },
      "evidences": [
        {
          "evidenceCode": "ECO:0000269",
          "source": "PubMed",
          "id": "1565966"
        }
      ]
    }
  ]
}
```

#### Comment Types

| Type | Description |
|------|-------------|
| FUNCTION | Protein function |
| SUBUNIT | Quaternary structure |
| SUBCELLULAR LOCATION | Cellular localization |
| TISSUE SPECIFICITY | Expression patterns |
| DISEASE | Disease associations |
| PTM | Post-translational modifications |
| INTERACTION | Protein interactions |
| COFACTOR | Required cofactors |
| CATALYTIC ACTIVITY | Enzymatic activity |
| PATHWAY | Pathway involvement |

#### Features Array

```json
{
  "features": [
    {
      "type": "Chain",
      "location": {
        "start": {"value": 1},
        "end": {"value": 393}
      },
      "description": "Cellular tumor antigen p53",
      "featureId": "PRO_0000185703"
    },
    {
      "type": "DNA binding",
      "location": {
        "start": {"value": 102},
        "end": {"value": 292}
      },
      "description": "p53-type"
    },
    {
      "type": "Modified residue",
      "location": {
        "start": {"value": 15},
        "end": {"value": 15}
      },
      "description": "Phosphoserine"
    }
  ]
}
```

#### Feature Types

| Type | Description |
|------|-------------|
| Chain | Mature protein |
| Signal | Signal peptide |
| Propeptide | Propeptide |
| Domain | Protein domain |
| Region | Region of interest |
| Motif | Short sequence motif |
| Binding site | Binding site |
| Active site | Active site |
| Modified residue | PTM site |
| Disulfide bond | Disulfide bond |
| Natural variant | Sequence variant |
| Mutagenesis | Experimental mutation |

#### Cross-References

```json
{
  "uniProtKBCrossReferences": [
    {
      "database": "Reactome",
      "id": "R-HSA-109581",
      "properties": [
        {
          "key": "PathwayName",
          "value": "Apoptosis"
        }
      ]
    },
    {
      "database": "KEGG",
      "id": "hsa:7157"
    },
    {
      "database": "GO",
      "id": "GO:0006915",
      "properties": [
        {
          "key": "GoTerm",
          "value": "P:apoptotic process"
        },
        {
          "key": "GoEvidenceType",
          "value": "IDA:UniProtKB"
        }
      ]
    },
    {
      "database": "STRING",
      "id": "9606.ENSP00000269305"
    },
    {
      "database": "PDB",
      "id": "1A1U",
      "properties": [
        {
          "key": "Method",
          "value": "X-ray"
        },
        {
          "key": "Resolution",
          "value": "2.20 A"
        }
      ]
    }
  ]
}
```

### 4.6 TSV Response Format

```
Entry	Protein names	Gene Names	Organism	Reactome	KEGG
P04637	Cellular tumor antigen p53	TP53 P53	Homo sapiens (Human)	R-HSA-109581;R-HSA-111447;...	hsa:7157
```

### 4.7 Evidence Codes

| Code | Description |
|------|-------------|
| ECO:0000269 | Experimental evidence from publication |
| ECO:0000250 | Sequence similarity |
| ECO:0000255 | Sequence model |
| ECO:0000256 | Sequence analysis |
| ECO:0000305 | Curator inference |
| ECO:0000303 | Non-traceable author statement |

---

## 5. PubChem

**Base URL:** `https://pubchem.ncbi.nlm.nih.gov/rest/pug/`

### 5.1 URL Pattern

```
/rest/pug/<domain>/<namespace>/<identifiers>/<operation>/<output>
```

#### Domains

| Domain | Description |
|--------|-------------|
| `compound` | Chemical compounds |
| `substance` | Substance records |
| `assay` | Bioassay data |
| `gene` | Gene records |
| `protein` | Protein records |
| `pathway` | Pathway data |
| `cell` | Cell line data |
| `sources` | Data sources |

#### Namespaces (Compound)

| Namespace | Description |
|-----------|-------------|
| `cid` | Compound ID |
| `name` | Common name |
| `smiles` | SMILES string |
| `inchi` | InChI string |
| `inchikey` | InChI key |
| `formula` | Molecular formula |

#### Operations

| Operation | Description |
|-----------|-------------|
| `record` | Full record |
| `property` | Specific properties |
| `synonyms` | Name synonyms |
| `sids` | Substance IDs |
| `aids` | Assay IDs |
| `cids` | Compound IDs |
| `description` | Text descriptions |

### 5.2 CID Lookup Response

```
GET /compound/name/aspirin/cids/JSON
```

```json
{
  "IdentifierList": {
    "CID": [2244]
  }
}
```

### 5.3 Compound Full Record

```
GET /compound/cid/2244/JSON
```

```json
{
  "PC_Compounds": [
    {
      "id": {
        "id": {
          "cid": 2244
        }
      },
      "atoms": {
        "aid": [1, 2, 3, ...],
        "element": [8, 6, 6, ...]
      },
      "bonds": {
        "aid1": [1, 2, 3, ...],
        "aid2": [2, 3, 4, ...],
        "order": [2, 1, 2, ...]
      },
      "coords": [
        {
          "type": [1, 5, 255],
          "conformers": [
            {
              "x": [3.7320, 2.8660, ...],
              "y": [0.0600, 0.5600, ...],
              "style": {
                "annotation": [8, 8, ...],
                "aid1": [5, 5, ...],
                "aid2": [6, 7, ...]
              }
            }
          ]
        }
      ],
      "charge": 0,
      "props": [
        {
          "urn": {
            "label": "Compound",
            "name": "Canonicalized",
            "datatype": 5,
            "release": "2024.01.08"
          },
          "value": {
            "ival": 1
          }
        },
        {
          "urn": {
            "label": "IUPAC Name",
            "name": "Preferred"
          },
          "value": {
            "sval": "2-acetyloxybenzoic acid"
          }
        },
        {
          "urn": {
            "label": "Log P"
          },
          "value": {
            "fval": 1.2
          }
        },
        {
          "urn": {
            "label": "Molecular Formula"
          },
          "value": {
            "sval": "C9H8O4"
          }
        },
        {
          "urn": {
            "label": "SMILES",
            "name": "Canonical"
          },
          "value": {
            "sval": "CC(=O)OC1=CC=CC=C1C(=O)O"
          }
        },
        {
          "urn": {
            "label": "InChI"
          },
          "value": {
            "sval": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
          }
        },
        {
          "urn": {
            "label": "InChIKey"
          },
          "value": {
            "sval": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
          }
        }
      ],
      "count": {
        "heavy_atom": 13,
        "atom_chiral": 0,
        "atom_chiral_def": 0,
        "atom_chiral_undef": 0,
        "bond_chiral": 0,
        "bond_chiral_def": 0,
        "bond_chiral_undef": 0,
        "isotope_atom": 0,
        "covalent_unit": 1,
        "tautomers": 1
      }
    }
  ]
}
```

### 5.4 Property Request

```
GET /compound/cid/2244/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey/JSON
```

```json
{
  "PropertyTable": {
    "Properties": [
      {
        "CID": 2244,
        "MolecularFormula": "C9H8O4",
        "MolecularWeight": "180.16",
        "CanonicalSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
        "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
      }
    ]
  }
}
```

### 5.5 Available Properties

| Property | Description |
|----------|-------------|
| `MolecularFormula` | Chemical formula |
| `MolecularWeight` | Molecular weight |
| `CanonicalSMILES` | Canonical SMILES |
| `IsomericSMILES` | Isomeric SMILES |
| `InChI` | InChI identifier |
| `InChIKey` | Hashed InChI |
| `IUPACName` | IUPAC name |
| `XLogP` | Computed logP |
| `ExactMass` | Exact mass |
| `MonoisotopicMass` | Monoisotopic mass |
| `TPSA` | Topological polar surface area |
| `Complexity` | Molecular complexity |
| `Charge` | Formal charge |
| `HBondDonorCount` | H-bond donors |
| `HBondAcceptorCount` | H-bond acceptors |
| `RotatableBondCount` | Rotatable bonds |
| `HeavyAtomCount` | Heavy atom count |
| `AtomStereoCount` | Stereo atoms |
| `DefinedAtomStereoCount` | Defined stereo atoms |
| `UndefinedAtomStereoCount` | Undefined stereo atoms |
| `BondStereoCount` | Stereo bonds |
| `CovalentUnitCount` | Covalent units |

### 5.6 Bioassay Data Model

```
GET /assay/aid/1259411/JSON
```

```json
{
  "PC_AssaySubmit": {
    "assay": {
      "descr": {
        "aid": {
          "id": 1259411,
          "version": 1
        },
        "aid_source": {
          "db": {
            "name": "CCRIS",
            "source_id": {
              "str": "22070"
            }
          }
        },
        "name": "CCRIS carcinogenicity studies",
        "description": [...],
        "protocol": [...],
        "xref": [...],
        "results": [
          {
            "tid": 1,
            "name": "Species",
            "type": "string"
          },
          {
            "tid": 2,
            "name": "Strain/Sex",
            "type": "string"
          },
          {
            "tid": 3,
            "name": "Route",
            "type": "string"
          },
          {
            "tid": 4,
            "name": "Dose",
            "type": "string"
          },
          {
            "tid": 5,
            "name": "Tumor Site: Type of Lesion",
            "type": "string"
          },
          {
            "tid": 6,
            "name": "Results",
            "type": "string"
          },
          {
            "tid": 7,
            "name": "Reference",
            "type": "string"
          }
        ],
        "activity_outcome_method": "confirmatory"
      }
    },
    "data": [
      {
        "sid": 363897951,
        "outcome": 1,
        "data": [
          {"tid": 1, "value": {"sval": "rat"}},
          {"tid": 2, "value": {"sval": "Sprague-Dawley/male"}},
          {"tid": 3, "value": {"sval": "oral"}},
          {"tid": 6, "value": {"sval": "NEGATIVE"}}
        ]
      }
    ]
  }
}
```

#### Outcome Codes

| Code | Meaning |
|------|---------|
| 1 | Inactive/Negative |
| 2 | Active/Positive |
| 3 | Inconclusive |
| 4 | Unspecified |
| 5 | Probe |

### 5.7 ID Systems

| ID Type | Format | Example |
|---------|--------|---------|
| CID | Integer | 2244 |
| SID | Integer | 363897951 |
| AID | Integer | 1259411 |
| GI (Protein) | Integer | 4507879 |
| Gene ID | Integer | 7157 |

---

## 6. STRING/STITCH

**STRING Base URL:** `https://string-db.org/api/`

**STITCH Base URL:** `https://stitch.embl.de/api/`

### 6.1 API Endpoints

```
GET /api/{format}/{method}?{parameters}
```

#### Formats

| Format | Description |
|--------|-------------|
| `json` | JSON format |
| `tsv` | Tab-separated values |
| `tsv-no-header` | TSV without header |
| `psi-mi` | PSI-MI XML |
| `psi-mi-tab` | PSI-MI TAB |

#### Methods

| Method | Description |
|--------|-------------|
| `network` | Get interaction network |
| `interaction_partners` | Get interaction partners |
| `resolve` | Resolve protein identifiers |
| `enrichment` | Functional enrichment |
| `functional_annotation` | GO terms and pathways |
| `get_string_ids` | Map to STRING IDs |

### 6.2 Network Request

```
GET /api/json/network?identifiers=TP53&species=9606
```

### 6.3 Interaction Response Structure

```json
[
  {
    "stringId_A": "9606.ENSP00000269305",
    "stringId_B": "9606.ENSP00000391127",
    "preferredName_A": "TP53",
    "preferredName_B": "SFN",
    "ncbiTaxonId": 9606,
    "score": 0.999,
    "nscore": 0,
    "fscore": 0,
    "pscore": 0,
    "ascore": 0.088,
    "escore": 0.926,
    "dscore": 0.9,
    "tscore": 0.986
  }
]
```

### 6.4 Score Components

| Score | Full Name | Description |
|-------|-----------|-------------|
| `score` | Combined score | Overall confidence (0-1) |
| `nscore` | Neighborhood | Gene neighborhood evidence |
| `fscore` | Fusion | Gene fusion evidence |
| `pscore` | Phylogenetic | Co-occurrence across genomes |
| `ascore` | Coexpression | Gene co-expression evidence |
| `escore` | Experimental | Experimental evidence |
| `dscore` | Database | Curated database evidence |
| `tscore` | Textmining | Text mining evidence |

**Score Calculation:**
```
combined_score = 1 - (1-nscore)*(1-fscore)*(1-pscore)*(1-ascore)*(1-escore)*(1-dscore)*(1-tscore)
```

Scores are corrected for randomly expected interactions.

### 6.5 Protein Resolution

```
GET /api/json/resolve?identifier=TP53&species=9606
```

```json
[
  {
    "queryIndex": 0,
    "queryItem": "TP53",
    "stringId": "9606.ENSP00000269305",
    "ncbiTaxonId": 9606,
    "taxonName": "Homo sapiens",
    "preferredName": "TP53",
    "annotation": "tumor suppressor in many tumor types; induces growth arrest or apoptosis..."
  }
]
```

### 6.6 Functional Annotation

```
GET /api/json/functional_annotation?identifiers=TP53&species=9606
```

```json
[
  {
    "category": "Process",
    "term": "GO:0006915",
    "description": "apoptotic process",
    "number_of_genes": 1,
    "number_of_genes_in_background": 1489,
    "inputGenes": "TP53",
    "preferredNames": "TP53"
  },
  {
    "category": "KEGG",
    "term": "hsa04115",
    "description": "p53 signaling pathway",
    "number_of_genes": 1,
    "number_of_genes_in_background": 72,
    "inputGenes": "TP53",
    "preferredNames": "TP53"
  }
]
```

### 6.7 Enrichment Analysis

```
GET /api/json/enrichment?identifiers=TP53%0dBRCA1%0dATM&species=9606
```

```json
[
  {
    "category": "KEGG",
    "term": "hsa01524",
    "description": "Platinum drug resistance",
    "number_of_genes": 3,
    "number_of_genes_in_background": 70,
    "ncbiTaxonId": 9606,
    "inputGenes": "TP53,BRCA1,ATM",
    "preferredNames": "TP53,BRCA1,ATM",
    "p_value": 4.88e-08,
    "fdr": 1.64e-05
  },
  {
    "category": "Process",
    "term": "GO:0072331",
    "description": "signal transduction by p53 class mediator",
    "number_of_genes": 3,
    "number_of_genes_in_background": 155,
    "p_value": 5.97e-08,
    "fdr": 1.91e-05
  }
]
```

### 6.8 STRING Identifier Format

```
{taxon_id}.{Ensembl_protein_id}

Example: 9606.ENSP00000269305
- 9606 = NCBI Taxonomy ID for Homo sapiens
- ENSP00000269305 = Ensembl protein ID
```

### 6.9 STITCH (Compound-Protein Interactions)

STITCH extends STRING for chemical-protein interactions.

#### STITCH Compound ID Format

```
CIDm{pubchem_cid}    # Merged stereoisomers
CIDs{pubchem_cid}    # Specific stereoisomer

Example: CIDm000002244 (Aspirin)
```

#### STITCH Interaction Types

| Type | Description |
|------|-------------|
| `activation` | Compound activates protein |
| `inhibition` | Compound inhibits protein |
| `binding` | Direct binding |
| `reaction` | Metabolic reaction |
| `catalysis` | Enzyme catalysis |
| `expression` | Affects expression |

### 6.10 Request Parameters

| Parameter | Description |
|-----------|-------------|
| `identifiers` | Newline-separated protein names |
| `species` | NCBI Taxonomy ID |
| `required_score` | Minimum score threshold (0-1000) |
| `network_type` | `functional` or `physical` |
| `add_nodes` | Number of additional interactors |
| `limit` | Maximum results |
| `caller_identity` | Application identifier |

---

## 7. Cross-Database ID Mappings

### 7.1 Common Identifier Systems

| Database | ID Type | Format | Example |
|----------|---------|--------|---------|
| UniProt | Accession | `[A-Z][0-9][A-Z0-9]{3}[0-9]` | P04637 |
| Ensembl Gene | ENSG | `ENSG[0-9]{11}` | ENSG00000141510 |
| Ensembl Protein | ENSP | `ENSP[0-9]{11}` | ENSP00000269305 |
| NCBI Gene | Gene ID | Integer | 7157 |
| KEGG Gene | org:id | `[a-z]{3}:[0-9]+` | hsa:7157 |
| PubChem | CID | Integer | 2244 |
| ChEBI | CHEBI ID | `CHEBI:[0-9]+` | CHEBI:15365 |
| KEGG Compound | C number | `C[0-9]{5}` | C00027 |
| Reactome | Stable ID | `R-{species}-{number}` | R-HSA-109581 |
| GO | GO ID | `GO:[0-9]{7}` | GO:0006915 |
| STRING | taxon.protein | `{taxon}.ENSP{number}` | 9606.ENSP00000269305 |

### 7.2 ID Conversion Services

#### UniProt ID Mapping
```
POST https://rest.uniprot.org/idmapping/run
{
  "from": "UniProtKB_AC-ID",
  "to": "KEGG",
  "ids": "P04637"
}
```

#### KEGG Conversion
```
GET https://rest.kegg.jp/conv/uniprot/hsa:7157
# Returns: hsa:7157    up:P04637
```

#### STRING Resolution
```
GET https://string-db.org/api/json/get_string_ids?identifiers=TP53&species=9606
```

### 7.3 Mapping Table

| From | To | Service |
|------|-----|---------|
| UniProt | KEGG | UniProt ID Mapping |
| UniProt | Ensembl | UniProt cross-refs |
| UniProt | STRING | UniProt cross-refs |
| UniProt | Reactome | UniProt cross-refs |
| KEGG | UniProt | KEGG conv |
| KEGG | NCBI Gene | KEGG conv |
| KEGG | PubChem | KEGG conv |
| Gene Symbol | STRING | STRING resolve |
| Gene Symbol | UniProt | UniProt search |
| PubChem CID | KEGG | KEGG conv |
| Compound Name | PubChem | PubChem search |

### 7.4 Example Cross-Reference Flow

To link a compound to genes via pathways:

```
1. Compound (PubChem CID) -> KEGG Compound
   GET https://rest.kegg.jp/conv/compound/pubchem:2244

2. KEGG Compound -> Pathways
   GET https://rest.kegg.jp/link/pathway/cpd:C00027

3. Pathway -> Genes
   GET https://rest.kegg.jp/link/genes/path:hsa04210

4. KEGG Gene -> UniProt
   GET https://rest.kegg.jp/conv/uniprot/hsa:7157

5. UniProt -> Detailed Protein Data
   GET https://rest.uniprot.org/uniprotkb/P04637.json
```

---

## Appendix A: Response Format Summary

| Database | Formats | Default |
|----------|---------|---------|
| Reactome | JSON, TSV, PDF, PNG, SVG, SBML, SBGN | JSON |
| KEGG | Text, AASEQ, NTSEQ, MOL, KCF, IMAGE, KGML, JSON | Text |
| WikiPathways | JSON, XML, GPML, JPG, PDF, PNG, SVG | XML |
| UniProt | JSON, XML, TSV, FASTA, RDF, GFF | JSON |
| PubChem | JSON, XML, CSV, SDF, PNG | JSON |
| STRING | JSON, TSV, PSI-MI, PSI-MI-TAB | JSON |

## Appendix B: Rate Limits

| Database | Rate Limit | Notes |
|----------|------------|-------|
| Reactome | None specified | Reasonable use expected |
| KEGG | None specified | Bulk downloads restricted |
| WikiPathways | None specified | Use `caller_identity` |
| UniProt | Streaming recommended | Large result pagination |
| PubChem | 5 requests/second | Automatic throttling |
| STRING | Use `caller_identity` | Register application |

## Appendix C: Authentication

| Database | Auth Required | Method |
|----------|---------------|--------|
| Reactome | No | - |
| KEGG | No | - |
| WikiPathways | Write ops only | Session token |
| UniProt | No | - |
| PubChem | No | - |
| STRING | No | API key optional |

---

## Download

| Database | Method | Endpoint/URL |
|----------|--------|-------------|
| Reactome | REST API | https://reactome.org/ContentService/ |
| Reactome | SQL Dump | ftp://ftp.reactome.org/Reactome/ |
| KEGG | REST API | https://www.kegg.jp/kegg/rest/keggapi.html |
| KEGG | FTP Download | ftp://ftp.kegg.jp/pub/kegg/genes/ |
| WikiPathways | REST API | https://webservice.wikipathways.org/ |
| UniProt | REST API | https://rest.uniprot.org/ |
| UniProt | FASTA Download | https://www.uniprot.org/uniprotkb/ |
| PubChem | REST API | https://pubchem.ncbi.nlm.nih.gov/rest/pug/ |
| STRING | REST API | https://string-db.org/api/ |

---

## Data Format

| Format | Database | Description |
|--------|----------|-------------|
| JSON | Reactome, UniProt, KEGG | Structured REST responses |
| SPARQL | All RDF-compliant | Query language for linked data |
| KGML | KEGG | Pathway markup format |
| SBML | Reactome | Systems Biology Markup Language |
| XML | UniProt, PubChem | Hierarchical structure |
| TSV | All | Tab-separated for bulk data |
| FASTA | UniProt | Sequence format |
| RDF/OWL | WikiPathways | Semantic web representation |

---

## Schema

### Core Fields

| Field | Type | Description | Source |
|-------|------|-------------|--------|
| pathway_id | String | Unique pathway identifier | KEGG: hsa04260, Reactome: R-HSA-1234 |
| pathway_name | String | Descriptive pathway name | Chemokine signaling pathway |
| pathway_category | String | Pathway classification | Signal transduction, Metabolism |
| protein_id | String | UniProt or Ensembl protein ID | P35354 (UniProt) or ENSP00000001 |
| protein_name | String | Gene/protein symbol | COX2, TNF |
| protein_description | String | Full protein name | Cyclooxygenase-2, Tumor necrosis factor |
| gene_id | String | Ensembl gene ID | ENSG00000073756 |
| interaction_type | String | Type of relationship | "enzymatic_reaction", "binding", "phosphorylation" |
| stoichiometry | String | Ratio in reaction | "1:1", "2:3" |
| species | String | Organism taxonomy | "Homo sapiens" |
| evidence_score | Float | Confidence of annotation (0-1) | 0.95 |
| data_source | String | Originating database | "Reactome", "KEGG", "STRING" |

### Relationships

| From | To | Type | Via |
|------|-----|------|-----|
| Pathway | Protein | participates_in | Reactome Events |
| Protein | Gene | encodes | UniProt |
| Protein | Protein | interacts_with | STRING |
| Pathway | Disease | implicated_in | Reactome Diseases |

---

## Sample Data

### JSON Format

```json
{
  "pathway": {
    "id": "hsa04260",
    "name": "Chemokine signaling pathway",
    "category": "Signal transduction",
    "source": "KEGG",
    "proteins": [
      {
        "uniprot_id": "P35354",
        "gene_symbol": "COX2",
        "ensembl_id": "ENSG00000073756",
        "protein_name": "Cyclooxygenase-2",
        "interactions": [
          {
            "with_protein": "Q9BXM7",
            "interaction_type": "downstream_target",
            "evidence": 0.92,
            "source": "STRING"
          }
        ]
      }
    ]
  }
}
```

### Query Result Table

| Pathway | Gene Symbol | UniProt ID | Protein Name | Interaction Type | Evidence | Data Source |
|---------|------------|-----------|--------------|-----------------|----------|-----------|
| hsa04260 | CCL2 | P13500 | C-C motif chemokine ligand 2 | ligand | 0.98 | KEGG |
| hsa04260 | CCR2 | P41597 | C-C chemokine receptor 2 | receptor | 0.95 | STRING |
| R-HSA-109704 | ADORA2A | P29274 | Adenosine receptor A2a | signaling | 0.92 | Reactome |

---

## License

| Database | License | Commercial Use | Citation |
|----------|---------|-----------------|----------|
| Reactome | CC BY 4.0 | Yes | Required |
| KEGG | Academic/Commercial | Yes (paid tier) | Required |
| WikiPathways | CC BY 4.0 | Yes | Required |
| UniProt | CC BY 4.0 | Yes | Required |
| PubChem | CC0 1.0 | Yes | Recommended |
| STRING | CC BY 4.0 | Yes | Required |

---

## Data Set Size

| Resource | Records | Compressed Size | Uncompressed Size | Format |
|----------|---------|-----------------|-------------------|--------|
| Reactome Pathways | 2,700+ | 500 MB | 3 GB | JSON/RDF |
| Reactome Reactions | 14K+ | 600 MB | 4 GB | JSON |
| KEGG Pathways | 600+ | 200 MB | 1.5 GB | KGML |
| UniProt Human | 21K+ proteins | 2 GB | 12 GB | XML/FASTA |
| WikiPathways | 2,800+ | 300 MB | 2 GB | GPML |
| STRING (Human) | 19K+ proteins | 5 GB | 15 GB | TSV |
| **Integrated pathway data** | **~50K records** | **~8 GB** | **~37 GB** | - |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Pathway | Series of molecular interactions leading to a biological outcome | Apoptosis pathway |
| Target | Molecular entity (protein/gene) that a compound acts upon | COX-2 enzyme |
| Data Model | Structured representation of database entities and relationships | Schema design |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Reactome | Curated database of biological pathways and reactions | Pathway analysis |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Metabolic pathways |
| WikiPathways | Community-curated open pathway database | Collaborative curation |
| STRING | Search Tool for Retrieval of Interacting Genes/Proteins | PPI network |
| STITCH | Search Tool for Interacting Chemicals | Chemical-protein interactions |
| SBML | Systems Biology Markup Language | Pathway format |
| BioPAX | Biological Pathway Exchange format | Pathway interchange |
| GPML | Graphical Pathway Markup Language | WikiPathways format |
| KGML | KEGG Markup Language | KEGG pathway format |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST, GraphQL |
| BioPAX | Biological Pathway Exchange | OWL/RDF format |
| GPML | Graphical Pathway Markup Language | XML-based |
| JSON | JavaScript Object Notation | Data format |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway DB |
| KGML | KEGG Markup Language | Pathway format |
| PPI | Protein-Protein Interaction | Network edges |
| REST | Representational State Transfer | API style |
| SBML | Systems Biology Markup Language | Model format |
| TSV | Tab-Separated Values | Data format |
| XML | Extensible Markup Language | Data format |

---

*Document generated: 2025*
*Last updated: January 2025*
