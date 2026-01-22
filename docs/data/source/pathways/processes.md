# Biological Process and Ontology Databases

This catalog documents major databases for biological processes, ontologies, enzyme data, and biochemical reactions. These resources are essential for understanding gene function, metabolic pathways, and cellular processes.

---

## Table of Contents

1. [Gene Ontology (GO)](#1-gene-ontology-go)
2. [Rhea](#2-rhea)
3. [BRENDA](#3-brenda)
4. [IntEnz](#4-intenz)
5. [Reactome (Biological Processes)](#5-reactome-biological-processes)
6. [UniPathway](#6-unipathway)
7. [SABIO-RK](#7-sabio-rk)
8. [BioModels](#8-biomodels)
9. [ChEBI](#9-chebi)
10. [Human Protein Atlas](#10-human-protein-atlas)
11. [Cross-Database Integration](#cross-database-integration)
12. [Summary Comparison](#summary-comparison)

---

## 1. Gene Ontology (GO)

**URL:** http://geneontology.org
**Maintained by:** Gene Ontology Consortium

### Overview

The Gene Ontology (GO) is the most widely used structured vocabulary for describing gene product functions. It provides a computational representation of current scientific knowledge about the functions of genes across all species.

### Three Ontologies

| Ontology | Abbreviation | Description | Root Term |
|----------|--------------|-------------|-----------|
| **Molecular Function** | MF | Activities performed by gene products at molecular level | GO:0003674 |
| **Biological Process** | BP | Larger processes accomplished by multiple molecular activities | GO:0008150 |
| **Cellular Component** | CC | Locations where gene products perform functions | GO:0005575 |

#### Molecular Function Examples
- Catalytic activity (GO:0003824)
- Binding (GO:0005488)
- Transporter activity (GO:0005215)
- Transcription regulator activity (GO:0140110)

#### Biological Process Examples
- Metabolic process (GO:0008152)
- Cellular process (GO:0009987)
- Signaling (GO:0023052)
- Response to stimulus (GO:0050896)

#### Cellular Component Examples
- Membrane (GO:0016020)
- Nucleus (GO:0005634)
- Cytoplasm (GO:0005737)
- Extracellular region (GO:0005576)

### GO Annotations (GAF Files)

Gene Ontology Annotation Files (GAF) link gene products to GO terms.

#### GAF 2.2 Format (15 columns)

| Column | Field | Description |
|--------|-------|-------------|
| 1 | DB | Database contributing the annotation |
| 2 | DB Object ID | Unique identifier in DB |
| 3 | DB Object Symbol | Gene symbol |
| 4 | Qualifier | Annotation qualifier (NOT, contributes_to, etc.) |
| 5 | GO ID | GO term identifier |
| 6 | DB:Reference | Reference for annotation |
| 7 | Evidence Code | Evidence supporting annotation |
| 8 | With/From | Additional identifiers |
| 9 | Aspect | F (function), P (process), C (component) |
| 10 | DB Object Name | Full gene name |
| 11 | DB Object Synonym | Synonyms |
| 12 | DB Object Type | gene, protein, transcript, etc. |
| 13 | Taxon | NCBI taxonomy ID |
| 14 | Date | Annotation date (YYYYMMDD) |
| 15 | Assigned By | Group that made annotation |

#### Evidence Codes

| Code | Type | Description |
|------|------|-------------|
| **Experimental** | | |
| EXP | Inferred from Experiment | General experimental evidence |
| IDA | Inferred from Direct Assay | Direct biochemical assay |
| IPI | Inferred from Physical Interaction | Physical interaction evidence |
| IMP | Inferred from Mutant Phenotype | Mutant analysis |
| IGI | Inferred from Genetic Interaction | Genetic interaction data |
| IEP | Inferred from Expression Pattern | Expression evidence |
| **Computational** | | |
| ISS | Inferred from Sequence Similarity | Sequence-based prediction |
| ISO | Inferred from Sequence Orthology | Ortholog evidence |
| ISA | Inferred from Sequence Alignment | Alignment-based |
| ISM | Inferred from Sequence Model | HMM or other model |
| IGC | Inferred from Genomic Context | Genomic context |
| IBA | Inferred from Biological Ancestor | Phylogenetic inference |
| IBD | Inferred from Biological Descendant | Phylogenetic inference |
| IKR | Inferred from Key Residues | Key residue analysis |
| IRD | Inferred from Rapid Divergence | Rapid divergence |
| RCA | Reviewed Computational Analysis | Reviewed computational |
| **Author/Curator** | | |
| TAS | Traceable Author Statement | Author statement with reference |
| NAS | Non-traceable Author Statement | Author statement without reference |
| IC | Inferred by Curator | Curator inference |
| **Electronic** | | |
| IEA | Inferred from Electronic Annotation | Automated annotation |

### GO-CAM (Causal Activity Models)

GO-CAM extends traditional GO annotations by providing causal models that describe how gene products work together.

#### Key Features
- **Activities:** Gene products performing molecular functions
- **Causal Relations:** How activities influence each other
- **Biological Context:** Where and when activities occur
- **Standard Relations:** Uses RO (Relations Ontology) predicates

#### GO-CAM Relations
```
directly_provides_input_for
directly_positively_regulates
directly_negatively_regulates
is_small_molecule_activator_of
is_small_molecule_inhibitor_of
```

#### Model Structure
```json
{
  "id": "gomodel:5fce9b7300002411",
  "title": "Insulin signaling pathway",
  "activities": [
    {
      "gene_product": "UniProtKB:P06213",
      "molecular_function": "GO:0004716",
      "biological_process": "GO:0008286",
      "cellular_component": "GO:0005886"
    }
  ],
  "causal_relations": [
    {
      "subject": "activity_1",
      "predicate": "directly_positively_regulates",
      "object": "activity_2"
    }
  ]
}
```

### Data Formats

#### OBO Format
```obo
[Term]
id: GO:0008150
name: biological_process
namespace: biological_process
def: "A biological process represents a specific objective that the organism is genetically programmed to achieve."
subset: goslim_generic

[Term]
id: GO:0008152
name: metabolic process
namespace: biological_process
def: "The chemical reactions and pathways, including anabolism and catabolism."
is_a: GO:0008150 ! biological_process
```

#### OWL Format
```xml
<owl:Class rdf:about="http://purl.obolibrary.org/obo/GO_0008152">
    <rdfs:label>metabolic process</rdfs:label>
    <rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/GO_0008150"/>
    <obo:IAO_0000115>The chemical reactions and pathways...</obo:IAO_0000115>
</owl:Class>
```

### AmiGO Browser

**URL:** http://amigo.geneontology.org

AmiGO is the official web-based tool for searching and browsing the Gene Ontology.

#### Features
- Term search and browsing
- Gene product search
- Annotation statistics
- GO term enrichment analysis (via PANTHER)
- Visualization tools

#### Search Types
| Search Type | Description |
|-------------|-------------|
| Ontology | Browse/search GO terms |
| Annotations | Search gene-to-GO associations |
| Gene Products | Search specific genes |
| Bioentities | Combined gene/protein search |

### API Access

#### GO REST API
```bash
# Get GO term
curl "http://api.geneontology.org/api/ontology/term/GO:0008150"

# Search terms
curl "http://api.geneontology.org/api/search/entity?q=apoptosis"

# Get annotations for a gene
curl "http://api.geneontology.org/api/bioentity/gene/HGNC:11998/function"
```

#### GOlr (Solr-based search)
```bash
# Search annotations
curl "http://golr-aux.geneontology.org/solr/select?q=TP53&fq=document_category:annotation"
```

### Bulk Downloads

**FTP:** http://current.geneontology.org/

| File | Description |
|------|-------------|
| `go-basic.obo` | Basic GO ontology (filtered) |
| `go.obo` | Full GO ontology |
| `go.owl` | OWL format |
| `goa_human.gaf.gz` | Human annotations |
| `goa_uniprot_all.gaf.gz` | All UniProt annotations |
| `go-plus.owl` | GO with external axioms |

#### Species-Specific Downloads
```
annotations/
├── goa_human.gaf.gz
├── goa_mouse.gaf.gz
├── goa_rat.gaf.gz
├── goa_fly.gaf.gz
├── goa_worm.gaf.gz
├── goa_yeast.gaf.gz
└── goa_zebrafish.gaf.gz
```

### Links to Genes/Compounds

| External Resource | Link Type |
|-------------------|-----------|
| UniProt | Primary protein identifiers |
| NCBI Gene | Gene identifiers |
| Ensembl | Transcript/protein IDs |
| ChEBI | Chemical entity annotations |
| PubMed | Literature references |
| Reactome | Pathway cross-references |

---

## 2. Rhea

**URL:** https://www.rhea-db.org
**Maintained by:** Swiss Institute of Bioinformatics (SIB)

### Overview

Rhea is an expert-curated knowledgebase of biochemical reactions. All reactions are balanced, transport-aware, and described using the ChEBI (Chemical Entities of Biological Interest) ontology.

### Coverage

| Statistic | Count (approx.) |
|-----------|-----------------|
| Total reactions | 15,000+ |
| ChEBI compounds | 17,000+ |
| Enzyme-linked reactions | 11,000+ |
| Transport reactions | 1,500+ |

### Data Content

#### Reaction Types
- **Chemical reactions:** Substrate-product transformations
- **Transport reactions:** Movement across membranes
- **Spontaneous reactions:** Non-enzyme catalyzed
- **Enzyme-catalyzed:** Linked to EC numbers

#### Reaction Structure
```
RHEA:10000
Name: ATP + H2O = ADP + phosphate
Status: approved
EC: 3.6.1.3
Direction: bidirectional

Participants:
  Left: CHEBI:30616 (ATP), CHEBI:15377 (H2O)
  Right: CHEBI:456216 (ADP), CHEBI:43474 (phosphate)
```

### ChEBI Integration

Every compound in Rhea is referenced using ChEBI identifiers:

| ChEBI Role | Description |
|------------|-------------|
| Substrate | Input compound |
| Product | Output compound |
| Cofactor | Required helper molecule |
| Catalyst | Enzyme (protein) |

### Data Formats

#### TSV Format
```tsv
RHEA_ID	EQUATION	EC	STATUS
10000	ATP + H2O = ADP + phosphate	3.6.1.3	approved
```

#### RDF/Turtle Format
```turtle
rhea:10000 a rh:Reaction ;
    rh:equation "ATP + H2O = ADP + phosphate" ;
    rh:status rh:Approved ;
    rh:ec "3.6.1.3" ;
    rh:side [ rh:contains rhea:compound_30616 ] ;
    rh:side [ rh:contains rhea:compound_456216 ] .
```

### SPARQL Access

**Endpoint:** https://sparql.rhea-db.org/sparql

#### Example Queries
```sparql
# Find all reactions catalyzed by a specific EC number
PREFIX rh: <http://rdf.rhea-db.org/>
SELECT ?reaction ?equation
WHERE {
  ?reaction rh:ec "2.7.1.1" ;
            rh:equation ?equation .
}

# Find reactions involving ATP
PREFIX rh: <http://rdf.rhea-db.org/>
PREFIX chebi: <http://purl.obolibrary.org/obo/CHEBI_>
SELECT ?reaction ?equation
WHERE {
  ?reaction rh:side ?side ;
            rh:equation ?equation .
  ?side rh:contains ?participant .
  ?participant rh:compound chebi:30616 .
}
```

### REST API

```bash
# Get reaction by Rhea ID
curl "https://www.rhea-db.org/rhea/10000.json"

# Search reactions
curl "https://www.rhea-db.org/rhea?query=ATP"

# Get all reactions for an EC number
curl "https://www.rhea-db.org/rhea?query=ec:2.7.1.1"
```

### Downloads

**FTP:** https://ftp.expasy.org/databases/rhea/

| File | Description |
|------|-------------|
| `rhea-reactions.txt.gz` | All reactions (TSV) |
| `rhea-directions.tsv` | Reaction directionality |
| `rhea2ec.tsv` | Rhea to EC mapping |
| `rhea2uniprot.tsv` | Rhea to UniProt mapping |
| `rhea2go.tsv` | Rhea to GO mapping |
| `rhea.rdf.gz` | Full RDF dump |

### Links to Genes/Compounds

| Resource | Link Type |
|----------|-----------|
| UniProt | Enzyme annotations |
| ChEBI | All compound references |
| EC numbers | Enzyme classification |
| GO | Molecular function terms |
| MetaCyc | Pathway integration |
| KEGG | Reaction mappings |

---

## 3. BRENDA

**URL:** https://www.brenda-enzymes.org
**Maintained by:** Technische Universität Braunschweig

### Overview

BRENDA (BRaunschweig ENzyme DAtabase) is the most comprehensive enzyme database, containing functional and molecular information on enzymes extracted from primary literature.

### Coverage

| Statistic | Count (approx.) |
|-----------|-----------------|
| EC numbers covered | 7,500+ |
| Organisms | 100,000+ |
| Literature references | 170,000+ |
| Protein sequences | 480,000+ |
| Kinetic parameters | 250,000+ |

### Data Content

#### Enzyme Information
- **Nomenclature:** Systematic/recommended names, synonyms
- **Reaction:** Catalyzed reactions
- **Substrate specificity:** Natural and artificial substrates
- **Inhibitors:** Compounds that inhibit activity
- **Activators:** Compounds that enhance activity

#### Kinetic Data

| Parameter | Description | Unit |
|-----------|-------------|------|
| Km | Michaelis constant | mM |
| Vmax | Maximum velocity | µmol/min/mg |
| kcat | Turnover number | s⁻¹ |
| kcat/Km | Catalytic efficiency | M⁻¹s⁻¹ |
| Ki | Inhibition constant | mM |
| IC50 | Half-maximal inhibitory concentration | µM |
| pH optimum | Optimal pH | - |
| Temperature optimum | Optimal temperature | °C |

#### Example Entry
```
EC 2.7.1.1 (Hexokinase)

SUBSTRATE:
  D-glucose (natural)
  D-fructose (natural)
  D-mannose

KINETIC PARAMETERS:
  Km (D-glucose): 0.12 mM [Homo sapiens, hexokinase I]
  Km (ATP): 0.45 mM [Homo sapiens, hexokinase I]
  kcat: 180 s⁻¹

INHIBITORS:
  glucose-6-phosphate (product inhibition, Ki = 0.02 mM)

ORGANISM:
  Homo sapiens
  Saccharomyces cerevisiae
  Rattus norvegicus
```

### Substrate/Product Relationships

BRENDA catalogs detailed substrate-product relationships:

```
REACTION: ATP + D-glucose = ADP + D-glucose-6-phosphate

Substrates:
  - ATP (CHEBI:30616)
  - D-glucose (CHEBI:17634)

Products:
  - ADP (CHEBI:456216)
  - D-glucose-6-phosphate (CHEBI:4170)

Cofactors:
  - Mg2+ (required)
```

### Data Formats

#### BRENDA Flat File Format
```
ID	EC_NUMBER
DE	Hexokinase
AN	ATP:D-hexose 6-phosphotransferase
SY	hexokinase type I
CA	ATP + D-hexose = ADP + D-hexose 6-phosphate
SU	D-glucose
SU	D-fructose
KM	0.12	mM	D-glucose	Homo sapiens	HK1
```

#### JSON Format (via API)
```json
{
  "ec": "2.7.1.1",
  "name": "hexokinase",
  "systematic_name": "ATP:D-hexose 6-phosphotransferase",
  "reaction": "ATP + D-hexose = ADP + D-hexose 6-phosphate",
  "substrates": ["D-glucose", "D-fructose", "D-mannose"],
  "km_values": [
    {
      "substrate": "D-glucose",
      "value": 0.12,
      "unit": "mM",
      "organism": "Homo sapiens"
    }
  ]
}
```

### API Access

BRENDA provides a SOAP API (requires registration):

```python
# Python example using zeep
from zeep import Client

client = Client('https://www.brenda-enzymes.org/soap/brenda_server.php?wsdl')

# Get substrates for an EC number
result = client.service.getSubstrate(
    email="user@example.com",
    password="password",
    ecNumber="2.7.1.1",
    organism="Homo sapiens"
)
```

### Download Options

| Access | Requirements |
|--------|--------------|
| Web interface | Free, registration recommended |
| SOAP API | Registration required |
| Full download | Academic license required |
| Commercial use | License agreement |

#### Available Downloads (with license)
- Complete BRENDA flat file
- Enzyme-specific data dumps
- Cross-reference files

### Links to Genes/Compounds

| Resource | Link Type |
|----------|-----------|
| UniProt | Protein sequences |
| PDB | 3D structures |
| ChEBI | Compound identifiers |
| KEGG | Pathway integration |
| EC numbers | Primary classification |
| PubMed | Literature references |
| Gene names | HGNC, MGI, etc. |

---

## 4. IntEnz

**URL:** https://www.ebi.ac.uk/intenz/
**Maintained by:** European Bioinformatics Institute (EBI)

### Overview

IntEnz is the official repository of enzyme nomenclature, providing authoritative EC (Enzyme Commission) number assignments and enzyme classification.

### Coverage

| Statistic | Count (approx.) |
|-----------|-----------------|
| Total EC numbers | 8,400+ |
| Active entries | 6,500+ |
| Transferred entries | 500+ |
| Deleted entries | 1,400+ |

### EC Classification Hierarchy

```
EC 1.-.-.- Oxidoreductases
EC 2.-.-.- Transferases
    EC 2.1.-.- Transferring one-carbon groups
    EC 2.2.-.- Transferring aldehyde/ketone groups
    EC 2.3.-.- Acyltransferases
    EC 2.4.-.- Glycosyltransferases
    EC 2.5.-.- Transferring alkyl/aryl groups
    EC 2.6.-.- Transferring nitrogenous groups
    EC 2.7.-.- Transferring phosphorus-containing groups
        EC 2.7.1.- Phosphotransferases (alcohol acceptor)
            EC 2.7.1.1 Hexokinase
            EC 2.7.1.2 Glucokinase
            EC 2.7.1.11 6-phosphofructokinase
EC 3.-.-.- Hydrolases
EC 4.-.-.- Lyases
EC 5.-.-.- Isomerases
EC 6.-.-.- Ligases
EC 7.-.-.- Translocases
```

### Data Content

#### Entry Structure
```
EC 2.7.1.1

Accepted name: hexokinase
Systematic name: ATP:D-hexose 6-phosphotransferase
Reaction: ATP + D-hexose = ADP + D-hexose 6-phosphate

Other name(s):
  - hexokinase type I
  - hexokinase type II
  - brain hexokinase

Comments: D-Glucose, D-mannose, D-fructose, sorbitol and
D-glucosamine can act as acceptors...

History: EC 2.7.1.1 created 1961
```

### Data Formats

#### XML Format
```xml
<enzyme>
  <ec>2.7.1.1</ec>
  <accepted_name>hexokinase</accepted_name>
  <systematic_name>ATP:D-hexose 6-phosphotransferase</systematic_name>
  <reaction>ATP + D-hexose = ADP + D-hexose 6-phosphate</reaction>
  <cofactors>
    <cofactor>Mg2+</cofactor>
  </cofactors>
  <xrefs>
    <xref db="UniProt" id="P19367"/>
    <xref db="Rhea" id="RHEA:22968"/>
  </xrefs>
</enzyme>
```

### Integration with UniProt

IntEnz is tightly integrated with UniProt:

```
UniProt Entry: P19367 (HXK1_HUMAN)
EC: 2.7.1.1
Function: Hexokinase that phosphorylates glucose to glucose-6-phosphate
IntEnz link: https://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec=2.7.1.1
```

### API Access

```bash
# Get enzyme by EC number
curl "https://www.ebi.ac.uk/intenz/ws/enzyme/EC/2.7.1.1"

# Search enzymes
curl "https://www.ebi.ac.uk/intenz/ws/search?query=hexokinase"
```

### Downloads

**FTP:** https://ftp.ebi.ac.uk/pub/databases/intenz/

| File | Description |
|------|-------------|
| `enzyme.dat` | Full IntEnz database |
| `enzclass.txt` | EC classification tree |
| `enzyme.xml` | XML format |
| `enzyme.rdf` | RDF format |

### Links to Genes/Compounds

| Resource | Link Type |
|----------|-----------|
| UniProt | Primary protein link |
| Rhea | Reaction descriptions |
| GO | Molecular function terms |
| KEGG | Pathway mappings |
| BRENDA | Kinetic data |
| PDB | Structural data |

---

## 5. Reactome (Biological Processes)

**URL:** https://reactome.org
**Maintained by:** Reactome Project (EMBL-EBI, OICR, NYU)

### Overview

Reactome is a free, open-source, curated and peer-reviewed pathway database. This section focuses on its biological process organization.

### Event Hierarchy

Reactome organizes biological knowledge as a hierarchy of "events":

```
Top-Level Pathways
├── Cell Cycle
│   ├── Mitotic G1 phase and G1/S transition
│   │   ├── Cyclin D associated events in G1
│   │   └── G1/S Transition
│   ├── S Phase
│   └── Mitotic G2-G2/M phases
├── Metabolism
│   ├── Carbohydrate metabolism
│   ├── Lipid metabolism
│   └── Amino acid metabolism
├── Signal Transduction
│   ├── Signaling by EGFR
│   ├── Signaling by Insulin receptor
│   └── Signaling by Wnt
└── Immune System
    ├── Innate Immune System
    └── Adaptive Immune System
```

### Event Types

| Type | Description |
|------|-------------|
| Pathway | Collection of related events |
| Reaction | Biochemical transformation |
| BlackBoxEvent | Incompletely characterized event |
| Polymerisation | Chain-building reactions |
| Depolymerisation | Chain-breaking reactions |

### Reactions and Participants

#### Reaction Structure
```json
{
  "stId": "R-HSA-70171",
  "name": "Glycolysis",
  "schemaClass": "Pathway",
  "species": "Homo sapiens",
  "hasEvent": [
    {
      "stId": "R-HSA-70262",
      "name": "Hexokinase phosphorylates glucose",
      "schemaClass": "Reaction",
      "input": ["glucose", "ATP"],
      "output": ["glucose-6-phosphate", "ADP"],
      "catalyst": ["Hexokinase 1"]
    }
  ]
}
```

#### Participant Types
| Type | Description |
|------|-------------|
| Input | Consumed by reaction |
| Output | Produced by reaction |
| Catalyst | Enzyme facilitating reaction |
| Regulator | Positive/negative regulator |

### Biological Process Coverage

| Category | Pathways (approx.) |
|----------|-------------------|
| Metabolism | 400+ |
| Signal transduction | 300+ |
| Immune system | 250+ |
| Gene expression | 200+ |
| Cell cycle | 150+ |
| Apoptosis | 50+ |
| Autophagy | 40+ |

### API Access

```bash
# Get pathway details
curl "https://reactome.org/ContentService/data/pathway/R-HSA-70171/containedEvents"

# Get participants in a reaction
curl "https://reactome.org/ContentService/data/participants/R-HSA-70262"

# Search biological processes
curl "https://reactome.org/ContentService/search/query?query=apoptosis&types=Pathway"
```

### Downloads

| File | Description |
|------|-------------|
| `ReactomePathways.gmt` | Gene sets for pathway analysis |
| `ReactomePathwaysRelation.txt` | Pathway hierarchy |
| `Ensembl2Reactome.txt` | Gene to pathway mapping |
| `ChEBI2Reactome.txt` | Compound to pathway mapping |
| `UniProt2Reactome.txt` | Protein to pathway mapping |

### Links to Genes/Compounds

| Resource | Link Type |
|----------|-----------|
| UniProt | Protein participants |
| ChEBI | Small molecules |
| Ensembl | Gene identifiers |
| GO | Cross-references |
| KEGG | Pathway mappings |

---

## 6. UniPathway

**URL:** Archived (data integrated into UniProt)
**Status:** Discontinued, content preserved in UniProt

### Overview

UniPathway was a metabolic pathway database that provided a hierarchical classification of metabolic pathways. While no longer actively maintained, its data lives on in UniProt annotations.

### Historical Classification

```
Biosynthesis
├── Amino-acid biosynthesis
│   ├── Alanine biosynthesis
│   ├── Arginine biosynthesis
│   └── ...
├── Cofactor biosynthesis
│   ├── Biotin biosynthesis
│   ├── Cobalamin biosynthesis
│   └── ...
└── Lipid biosynthesis
    ├── Fatty acid biosynthesis
    └── ...

Degradation
├── Amino-acid degradation
├── Carbohydrate degradation
└── ...
```

### UniProt Integration

UniPathway identifiers (UPAxxxxx) are preserved in UniProt entries:

```
UniProt Entry: P00558 (PGK1_HUMAN)

Pathway:
  - Glycolysis; D-glyceraldehyde 3-phosphate from D-glucose: step 1/4.
  - UPA00109; D-glyceraldehyde 3-phosphate from D-glucose
```

### Data Access

UniPathway data can still be accessed through:
- UniProt pathway annotations
- Archived data files (limited availability)
- Cross-references in MetaCyc

### Pathway Structure

```
UPA00109: D-glyceraldehyde 3-phosphate from D-glucose
├── Step 1: Glucose → Glucose-6-phosphate (EC 2.7.1.1)
├── Step 2: Glucose-6-phosphate → Fructose-6-phosphate (EC 5.3.1.9)
├── Step 3: Fructose-6-phosphate → Fructose-1,6-bisphosphate (EC 2.7.1.11)
└── Step 4: Fructose-1,6-bisphosphate → 2 × D-glyceraldehyde-3-phosphate (EC 4.1.2.13)
```

---

## 7. SABIO-RK

**URL:** http://sabiork.h-its.org
**Maintained by:** Heidelberg Institute for Theoretical Studies (HITS)

### Overview

SABIO-RK (System for the Analysis of Biochemical Pathways - Reaction Kinetics) is a database for biochemical reaction kinetics. It provides standardized kinetic data for modeling.

### Coverage

| Statistic | Count (approx.) |
|-----------|-----------------|
| Kinetic records | 55,000+ |
| Reactions | 12,000+ |
| Enzymes | 3,500+ |
| Organisms | 900+ |
| Publications | 6,000+ |

### Data Content

#### Kinetic Parameters
| Parameter | Symbol | Description |
|-----------|--------|-------------|
| Michaelis constant | Km | Substrate concentration at half-Vmax |
| Maximum velocity | Vmax | Maximum reaction rate |
| Turnover number | kcat | Reactions per enzyme per second |
| Catalytic efficiency | kcat/Km | Efficiency measure |
| Inhibition constant | Ki | Inhibitor binding strength |
| Hill coefficient | nH | Cooperativity measure |
| Equilibrium constant | Keq | Reaction equilibrium |

#### Experimental Conditions
```json
{
  "kinetic_law_id": 12345,
  "reaction": "ATP + D-glucose = ADP + D-glucose-6-phosphate",
  "enzyme": "Hexokinase",
  "organism": "Homo sapiens",
  "tissue": "brain",
  "pH": 7.4,
  "temperature": 37,
  "temperature_unit": "°C",
  "parameters": [
    {
      "type": "Km",
      "value": 0.12,
      "unit": "mM",
      "substrate": "D-glucose"
    },
    {
      "type": "kcat",
      "value": 180,
      "unit": "s^(-1)"
    }
  ]
}
```

### Rate Constants

SABIO-RK stores various rate constants for kinetic modeling:

| Constant Type | Description |
|---------------|-------------|
| Forward rate (k+) | Rate of forward reaction |
| Reverse rate (k-) | Rate of reverse reaction |
| Association (kon) | Binding rate |
| Dissociation (koff) | Unbinding rate |

### Data Formats

#### SBML Export
```xml
<listOfReactions>
  <reaction id="R1" name="Hexokinase reaction">
    <kineticLaw>
      <math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply>
          <divide/>
          <apply><times/><ci>Vmax</ci><ci>glucose</ci></apply>
          <apply><plus/><ci>Km</ci><ci>glucose</ci></apply>
        </apply>
      </math>
      <listOfLocalParameters>
        <localParameter id="Km" value="0.12" units="mM"/>
        <localParameter id="Vmax" value="100" units="umol_per_min_per_mg"/>
      </listOfLocalParameters>
    </kineticLaw>
  </reaction>
</listOfReactions>
```

### API Access

```bash
# Search kinetic data
curl "http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/sbml?q=glucose%20AND%20hexokinase"

# Get specific entry
curl "http://sabiork.h-its.org/sabioRestWebServices/kineticLaws/12345"

# Export as SBML
curl "http://sabiork.h-its.org/sabioRestWebServices/kineticLaws/12345/sbml"
```

#### Query Parameters
| Parameter | Description |
|-----------|-------------|
| `Organism` | Filter by organism |
| `Tissue` | Filter by tissue |
| `ECNumber` | Filter by EC number |
| `Substrate` | Filter by substrate |
| `Product` | Filter by product |

### Downloads

| Format | Description |
|--------|-------------|
| SBML | Systems Biology Markup Language |
| CSV | Tabular data |
| JSON | Programmatic access |
| BioPAX | Biological Pathway Exchange |

### Links to Genes/Compounds

| Resource | Link Type |
|----------|-----------|
| UniProt | Enzyme identifiers |
| ChEBI | Compound identifiers |
| EC numbers | Enzyme classification |
| KEGG | Compound/reaction mapping |
| PubMed | Literature references |
| BioModels | Model integration |

---

## 8. BioModels

**URL:** https://www.ebi.ac.uk/biomodels/
**Maintained by:** European Bioinformatics Institute (EBI)

### Overview

BioModels is a repository of mathematical models of biological systems. Models are encoded in standard formats and can be used for simulation and analysis.

### Coverage

| Category | Count (approx.) |
|----------|-----------------|
| Curated models | 1,000+ |
| Non-curated models | 1,000+ |
| Path2Models | 140,000+ |
| Total entries | 150,000+ |

### Model Types

| Type | Description |
|------|-------------|
| Kinetic | ODE-based dynamic models |
| Logical | Boolean network models |
| Constraint-based | Flux balance models |
| Spatial | PDE-based spatial models |
| Stochastic | Stochastic simulation models |

### SBML Format

SBML (Systems Biology Markup Language) is the primary format:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="BIOMD0000000001" name="Glycolysis Model">
    <listOfCompartments>
      <compartment id="cytoplasm" size="1"/>
    </listOfCompartments>

    <listOfSpecies>
      <species id="glucose" compartment="cytoplasm" initialConcentration="5"/>
      <species id="G6P" compartment="cytoplasm" initialConcentration="0"/>
      <species id="ATP" compartment="cytoplasm" initialConcentration="2"/>
      <species id="ADP" compartment="cytoplasm" initialConcentration="0"/>
    </listOfSpecies>

    <listOfReactions>
      <reaction id="hexokinase" reversible="false">
        <listOfReactants>
          <speciesReference species="glucose" stoichiometry="1"/>
          <speciesReference species="ATP" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="G6P" stoichiometry="1"/>
          <speciesReference species="ADP" stoichiometry="1"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
              <ci>Vmax</ci>
              <apply><divide/>
                <ci>glucose</ci>
                <apply><plus/><ci>Km</ci><ci>glucose</ci></apply>
              </apply>
            </apply>
          </math>
        </kineticLaw>
      </reaction>
    </listOfReactions>

    <listOfParameters>
      <parameter id="Vmax" value="100"/>
      <parameter id="Km" value="0.12"/>
    </listOfParameters>
  </model>
</sbml>
```

### Additional Formats

| Format | Description |
|--------|-------------|
| SBML | Primary exchange format |
| CellML | Alternative model format |
| MATLAB | Simulation scripts |
| VCML | Virtual Cell format |
| SED-ML | Simulation descriptions |
| COMBINE Archive | Bundled model + simulations |

### Model Curation

| Level | Description |
|-------|-------------|
| Curated | Manually verified, reproduces publication |
| Non-curated | Community submitted, not verified |
| Path2Models | Automatically generated from pathways |

### API Access

```bash
# Get model information
curl "https://www.ebi.ac.uk/biomodels/BIOMD0000000001?format=json"

# Download SBML
curl "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000001?filename=BIOMD0000000001_url.xml"

# Search models
curl "https://www.ebi.ac.uk/biomodels/search?query=glycolysis&format=json"
```

### Downloads

**FTP:** https://ftp.ebi.ac.uk/pub/databases/biomodels/

| Directory | Content |
|-----------|---------|
| `releases/` | Release archives |
| `curated/` | Curated models |
| `non-curated/` | Community models |

### Links to Genes/Compounds

| Resource | Link Type |
|----------|-----------|
| UniProt | Enzyme identifiers in models |
| ChEBI | Compound identifiers |
| GO | Biological process annotations |
| Reactome | Pathway-derived models |
| KEGG | Pathway-derived models |
| PubMed | Publication references |

---

## 9. ChEBI

**URL:** https://www.ebi.ac.uk/chebi/
**Maintained by:** European Bioinformatics Institute (EBI)

### Overview

ChEBI (Chemical Entities of Biological Interest) is a freely available dictionary of molecular entities focused on small chemical compounds. It provides an ontology of chemical roles and relationships.

### Coverage

| Statistic | Count (approx.) |
|-----------|-----------------|
| Total entities | 175,000+ |
| 3-star (fully curated) | 50,000+ |
| 2-star (partially curated) | 30,000+ |
| 1-star (automated) | 95,000+ |
| Chemical roles | 2,500+ |

### Chemical Ontology

ChEBI organizes chemicals into a hierarchical ontology:

```
chemical entity (CHEBI:24431)
├── molecular entity (CHEBI:23367)
│   ├── atom (CHEBI:33250)
│   ├── molecular ion (CHEBI:25213)
│   └── molecule (CHEBI:25367)
│       ├── organic molecular entity (CHEBI:50860)
│       │   ├── carbohydrate (CHEBI:16646)
│       │   ├── lipid (CHEBI:18059)
│       │   ├── amino acid (CHEBI:33709)
│       │   └── nucleotide (CHEBI:36976)
│       └── inorganic molecular entity (CHEBI:24835)
└── chemical substance (CHEBI:59999)
```

### Role Ontology

ChEBI assigns biological and chemical roles:

```
role (CHEBI:50906)
├── biological role (CHEBI:24432)
│   ├── biochemical role (CHEBI:52206)
│   │   ├── cofactor (CHEBI:23357)
│   │   ├── coenzyme (CHEBI:23354)
│   │   ├── enzyme inhibitor (CHEBI:23924)
│   │   └── metabolite (CHEBI:25212)
│   │       ├── human metabolite (CHEBI:77746)
│   │       └── plant metabolite (CHEBI:76924)
│   └── pharmacological role
│       ├── drug (CHEBI:23888)
│       ├── antibiotic (CHEBI:22582)
│       └── anti-inflammatory agent (CHEBI:67079)
├── chemical role (CHEBI:51086)
│   ├── antioxidant (CHEBI:22586)
│   ├── buffer (CHEBI:22695)
│   └── solvent (CHEBI:46787)
└── application (CHEBI:33232)
    ├── food additive (CHEBI:64047)
    └── pharmaceutical (CHEBI:52217)
```

### Entry Structure

```
CHEBI:17234 (D-glucose)

Name: D-glucose
Definition: An aldohexose used as a source of energy and metabolic intermediate.

InChI: InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1
InChIKey: WQZGKKKJIJFFOK-GASJEMHNSA-N
SMILES: OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O

Formula: C6H12O6
Mass: 180.15588
Charge: 0

Roles:
  - is_a: monosaccharide (CHEBI:35381)
  - is_a: D-aldohexose (CHEBI:4194)
  - has_role: human metabolite (CHEBI:77746)
  - has_role: Saccharomyces cerevisiae metabolite

Relationships:
  - is_conjugate_acid_of: D-gluconate (CHEBI:18391)
  - is_enantiomer_of: L-glucose (CHEBI:37627)
```

### Data Formats

#### OBO Format
```obo
[Term]
id: CHEBI:17234
name: D-glucose
def: "An aldohexose used as a source of energy..."
synonym: "Dextrose" RELATED
synonym: "Grape sugar" RELATED
xref: KEGG:C00031
xref: PubChem:5793
is_a: CHEBI:4194 ! D-aldohexose
relationship: has_role CHEBI:77746 ! human metabolite
```

#### SDF Format
```sdf
D-glucose
     RDKit          3D

 24 24  0  0  0  0  0  0  0  0999 V2000
    1.0469    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ...
M  END
> <CHEBI_ID>
CHEBI:17234

> <CHEBI_NAME>
D-glucose

$$$$
```

### Links to Pathways

ChEBI links to pathway databases:

| Pathway Database | Link Type |
|-----------------|-----------|
| Reactome | Participant in reactions |
| Rhea | Reaction substrate/product |
| KEGG | Compound cross-reference |
| MetaCyc | Metabolic pathways |
| UniProt | Cofactor/ligand annotations |

### API Access

```bash
# Get compound by ChEBI ID
curl "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17234"

# Search compounds
curl "https://www.ebi.ac.uk/chebi/advancedSearchFT.do?searchString=glucose"

# Get ontology children
curl "https://www.ebi.ac.uk/ols/api/ontologies/chebi/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCHEBI_17234/children"
```

### Web Services

```bash
# SOAP service for programmatic access
# WSDL: https://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl

# OLS (Ontology Lookup Service)
curl "https://www.ebi.ac.uk/ols/api/ontologies/chebi"
```

### Downloads

**FTP:** https://ftp.ebi.ac.uk/pub/databases/chebi/

| File | Description |
|------|-------------|
| `chebi.obo.gz` | OBO format ontology |
| `chebi_lite.obo.gz` | Lite version (no structures) |
| `ChEBI_complete.sdf.gz` | Structure-data file |
| `ChEBI_complete_3star.sdf.gz` | 3-star entries only |
| `database_accession.tsv` | Cross-references |
| `relation.tsv` | Ontology relations |

### Links to Genes/Compounds

| Resource | Link Type |
|----------|-----------|
| UniProt | Cofactor/ligand annotations |
| Rhea | Reaction participants |
| Reactome | Pathway molecules |
| KEGG | Compound mappings |
| PubChem | Structure cross-reference |
| DrugBank | Drug compounds |
| HMDB | Human metabolites |

---

## 10. Human Protein Atlas

**URL:** https://www.proteinatlas.org
**Maintained by:** KTH Royal Institute of Technology, Uppsala University

### Overview

The Human Protein Atlas (HPA) provides comprehensive information about human protein expression at tissue and cellular levels, combining transcriptomics, proteomics, and antibody-based imaging.

### Major Sections

| Section | Description |
|---------|-------------|
| Tissue Atlas | Protein expression in normal tissues |
| Cell Atlas | Subcellular localization |
| Pathology Atlas | Cancer expression patterns |
| Blood Atlas | Proteins in blood |
| Brain Atlas | Brain region expression |
| Single Cell Atlas | Single-cell transcriptomics |
| Metabolic Atlas | Metabolic pathway mapping |

### Tissue Expression Data

#### Expression Categories
| Category | Description |
|----------|-------------|
| Tissue enriched | >4-fold higher in one tissue |
| Group enriched | >4-fold in 2-7 tissues |
| Tissue enhanced | >4-fold in one or more |
| Low tissue specificity | Broadly expressed |
| Not detected | No expression |

#### Expression Levels
```
Protein expression (IHC):
  - Not detected
  - Low
  - Medium
  - High

RNA expression (TPM):
  - Not detected (< 1 TPM)
  - Low (1-10 TPM)
  - Medium (10-100 TPM)
  - High (> 100 TPM)
```

### Subcellular Localization

The Cell Atlas maps protein localization to cellular compartments:

```
Cellular Compartments:
├── Nucleus
│   ├── Nucleoplasm
│   ├── Nuclear membrane
│   ├── Nucleoli
│   └── Nuclear speckles
├── Cytoplasm
│   ├── Cytosol
│   ├── Mitochondria
│   ├── Golgi apparatus
│   └── Endoplasmic reticulum
├── Plasma membrane
├── Cell junctions
└── Extracellular
```

#### Localization Categories
| Category | Description |
|----------|-------------|
| Main location | Primary localization |
| Additional location | Secondary localization |
| Uncertain | Low confidence |
| Approved | High confidence |
| Supported | Medium confidence |

### Pathology Atlas

Cancer-related protein expression:

#### Prognostic Data
```json
{
  "gene": "TP53",
  "cancer": "breast cancer",
  "prognostic_type": "unfavorable",
  "p_value": 0.0001,
  "expression_correlation": "high expression = poor survival"
}
```

#### Cancer Types Covered
- Breast cancer
- Colorectal cancer
- Lung cancer
- Prostate cancer
- Liver cancer
- Pancreatic cancer
- And many more (17+ types)

### Data Formats

#### TSV Format
```tsv
Gene	Gene name	Tissue	Cell type	Level	Reliability
ENSG00000141510	TP53	liver	hepatocytes	High	Approved
ENSG00000141510	TP53	brain	neurons	Medium	Supported
```

#### XML Format
```xml
<entry>
  <identifier db="Ensembl">ENSG00000141510</identifier>
  <name>TP53</name>
  <tissueExpression tissue="liver">
    <cellType name="hepatocytes" level="High" reliability="Approved"/>
  </tissueExpression>
  <subcellularLocation>
    <location>Nucleus</location>
    <location>Cytoplasm</location>
  </subcellularLocation>
</entry>
```

### API Access

```bash
# Get gene information
curl "https://www.proteinatlas.org/api/search_download.php?search=TP53&format=json"

# Get tissue expression
curl "https://www.proteinatlas.org/ENSG00000141510-TP53/tissue"

# Download expression data
curl "https://www.proteinatlas.org/download/normal_tissue.tsv.zip"
```

### Downloads

**URL:** https://www.proteinatlas.org/about/download

| File | Description |
|------|-------------|
| `normal_tissue.tsv.zip` | Tissue expression data |
| `subcellular_location.tsv.zip` | Localization data |
| `pathology.tsv.zip` | Cancer expression |
| `rna_tissue_consensus.tsv.zip` | RNA expression (tissues) |
| `rna_single_cell_type.tsv.zip` | Single-cell RNA data |

### Links to Genes/Compounds

| Resource | Link Type |
|----------|-----------|
| Ensembl | Primary gene identifier |
| UniProt | Protein identifiers |
| NCBI Gene | Gene cross-reference |
| Reactome | Pathway annotations |
| GO | Functional annotations |
| KEGG | Pathway mapping |
| STRING | Protein interactions |

---

## Cross-Database Integration

### Identifier Mapping

| Database | Primary ID | Cross-references |
|----------|------------|------------------|
| GO | GO:xxxxxxx | UniProt, EC, Reactome |
| Rhea | RHEA:xxxxx | ChEBI, EC, UniProt |
| BRENDA | EC x.x.x.x | UniProt, PDB, ChEBI |
| IntEnz | EC x.x.x.x | UniProt, Rhea, GO |
| Reactome | R-HSA-xxxxx | UniProt, ChEBI, GO |
| SABIO-RK | Entry ID | EC, ChEBI, UniProt |
| BioModels | BIOMDxxxxxxxx | UniProt, ChEBI, GO |
| ChEBI | CHEBI:xxxxx | KEGG, PubChem, Rhea |
| HPA | ENSG ID | UniProt, Ensembl |

### Common Integration Patterns

```
Gene → UniProt → GO (function)
                 → Reactome (pathways)
                 → ChEBI (ligands)
                 → HPA (expression)

Compound → ChEBI → Rhea (reactions)
                 → Reactome (pathways)
                 → SABIO-RK (kinetics)

Enzyme → EC → BRENDA (kinetics)
            → IntEnz (nomenclature)
            → Rhea (reactions)
            → UniProt (sequences)
```

### Integrated Query Example

```python
# Pseudocode for integrated query
def get_enzyme_comprehensive(gene_symbol):
    """Get comprehensive enzyme information across databases."""

    # Start with UniProt
    uniprot_id = map_gene_to_uniprot(gene_symbol)
    protein_info = fetch_uniprot(uniprot_id)

    # Get GO annotations
    go_terms = fetch_go_annotations(uniprot_id)

    # Get EC number
    ec_number = protein_info.ec_number

    # Get reaction information
    rhea_reactions = fetch_rhea_by_ec(ec_number)

    # Get kinetic data
    brenda_kinetics = fetch_brenda(ec_number, organism="Homo sapiens")
    sabio_kinetics = fetch_sabio_rk(ec_number)

    # Get pathway context
    reactome_pathways = fetch_reactome(uniprot_id)

    # Get expression data
    hpa_expression = fetch_hpa(gene_symbol)

    return {
        "protein": protein_info,
        "function": go_terms,
        "reactions": rhea_reactions,
        "kinetics": {
            "brenda": brenda_kinetics,
            "sabio": sabio_kinetics
        },
        "pathways": reactome_pathways,
        "expression": hpa_expression
    }
```

---

## Summary Comparison

### Database Focus Areas

| Database | Primary Focus | Data Type |
|----------|--------------|-----------|
| GO | Function vocabulary | Ontology, annotations |
| Rhea | Biochemical reactions | Reactions, equations |
| BRENDA | Enzyme properties | Kinetics, substrates |
| IntEnz | Enzyme classification | Nomenclature |
| Reactome | Biological pathways | Events, reactions |
| SABIO-RK | Reaction kinetics | Rate constants |
| BioModels | Mathematical models | SBML models |
| ChEBI | Chemical entities | Ontology, structures |
| HPA | Protein expression | Tissue/cell data |

### Access Methods

| Database | Web | REST API | SPARQL | Bulk Download | License |
|----------|-----|----------|--------|---------------|---------|
| GO | Yes | Yes | Yes | Yes | CC BY 4.0 |
| Rhea | Yes | Yes | Yes | Yes | CC BY 4.0 |
| BRENDA | Yes | SOAP | No | License | Academic free |
| IntEnz | Yes | Yes | No | Yes | Public domain |
| Reactome | Yes | Yes | No | Yes | CC0 |
| SABIO-RK | Yes | Yes | No | Yes | Free |
| BioModels | Yes | Yes | No | Yes | CC0 |
| ChEBI | Yes | Yes | Yes | Yes | CC BY 4.0 |
| HPA | Yes | Limited | No | Yes | CC BY-SA 3.0 |

### Data Format Support

| Database | TSV | XML | JSON | RDF | OBO | SBML | SDF |
|----------|-----|-----|------|-----|-----|------|-----|
| GO | Yes | Yes | Yes | Yes | Yes | - | - |
| Rhea | Yes | Yes | Yes | Yes | - | - | - |
| BRENDA | Yes | - | Yes | - | - | - | - |
| IntEnz | Yes | Yes | - | Yes | - | - | - |
| Reactome | Yes | Yes | Yes | Yes | - | Yes | - |
| SABIO-RK | Yes | Yes | Yes | - | - | Yes | - |
| BioModels | - | Yes | Yes | - | - | Yes | - |
| ChEBI | Yes | - | - | Yes | Yes | - | Yes |
| HPA | Yes | Yes | Yes | - | - | - | - |

### Update Frequency

| Database | Update Cycle |
|----------|--------------|
| GO | Monthly |
| Rhea | Quarterly |
| BRENDA | Annual |
| IntEnz | Continuous |
| Reactome | Quarterly |
| SABIO-RK | Continuous |
| BioModels | Continuous |
| ChEBI | Monthly |
| HPA | Annual |

---

## References and Further Reading

1. **Gene Ontology Consortium** (2021). The Gene Ontology resource: enriching a GOld mine. Nucleic Acids Research, 49(D1), D325-D334.

2. **Bansal, P., et al.** (2022). Rhea, the reaction knowledgebase in 2022. Nucleic Acids Research, 50(D1), D665-D670.

3. **Chang, A., et al.** (2021). BRENDA, the ELIXIR core data resource in 2021. Nucleic Acids Research, 49(D1), D498-D508.

4. **Jassal, B., et al.** (2020). The reactome pathway knowledgebase. Nucleic Acids Research, 48(D1), D498-D503.

5. **Wittig, U., et al.** (2018). SABIO-RK: an updated resource for manually curated biochemical reaction kinetics. Nucleic Acids Research, 46(D1), D471-D478.

6. **Malik-Sheriff, R.S., et al.** (2020). BioModels - 15 years of sharing computational models in life science. Nucleic Acids Research, 48(D1), D1259-D1268.

7. **Hastings, J., et al.** (2016). ChEBI in 2016: Improved services and an expanding collection of metabolites. Nucleic Acids Research, 44(D1), D1214-D1219.

8. **Uhlen, M., et al.** (2015). Tissue-based map of the human proteome. Science, 347(6220), 1260419.

---

*Document generated: 2025*
*For the latest information, please consult the respective database websites.*
