# Research Prompt History

This document captures all prompts used during the Gene Platform research phase (January 17-18, 2026).

---

## Session Overview

- **Total Research Output**: 50 files, 51,009 lines
- **Location**: `/home/ubuntu/src/gene/docs/research/`
- **Method**: Claude Code Task tool with parallel subagents (pre-claude-flow init)

---

## User Prompts (Chronological)

### 1. Traditional Medicine Research

**User Request:**
> "now we need to research sources of data, research papers, and databases of ayurvedic, tcm, kampo, and western herbal medicine, unless we have already done so? We need to know how we can download the data, and how to use it. We should be able to link compounds and formulations to genes and biochemical pathways. swarm"

**Output Files Created:**
- `interventions-tcm.md` (714 lines)
- `interventions-ayurveda.md` (833 lines)
- `interventions-kampo.md` (1,058 lines)
- `interventions-western-herbal.md` (813 lines)
- `DATA-INTEGRATION-GUIDE.md`

---

### 2. Data Models Research

**User Request:**
> "swarm and get the models and the model descriptions for all the data sources"

**Agent Prompts Used:**

**TCM Agent:**
```
Document the complete data models and schemas for TCM databases:
- BATMAN-TCM 2.0: API structure, JSON response format, field definitions
- TCMBank: Data structure, downloadable formats
- HERB 2.0: Database schema, compound/target/herb relationships
- TCMID: Entity relationships, data fields
- SymMap: Integration schema, mapping structures

For each, provide:
1. Entity-relationship diagrams (text-based)
2. Field definitions with data types
3. Example records/responses
4. How entities link to genes/pathways (field mappings)
```

**Ayurveda Agent:**
```
Document the complete data models and schemas for Ayurveda databases:
- IMPPAT 2.0: Database schema, phytochemical/plant/target structure
- NPACT: Data model, compound-activity relationships
- InDiaMed: Entity structure, formulation data
- Ayurveda Knowledge Base sources

For each, provide:
1. Entity-relationship diagrams (text-based)
2. Field definitions with data types
3. Example records
4. How to link to genes/pathways
```

**Kampo Agent:**
```
Document the complete data models and schemas for Kampo databases:
- KampoDB: 4-layer data model (Formula→Crude Drug→Constituent→Target)
- STORK: Schema structure
- Kampo-related entries in KEGG

For each, provide:
1. Entity-relationship diagrams (text-based)
2. Field definitions with data types
3. Example records
4. Gene/pathway linking fields
```

**Western Herbal Agent:**
```
Document the complete data models and schemas for Western herbal databases:
- DSLD (Dietary Supplement Label Database): API schema, product/ingredient structure
- Dr. Duke's Phytochemical Database: CSV structure, field definitions
- NIH ODS: Fact sheet data model
- Natural Medicines Database structure (if publicly documented)

For each, provide:
1. Entity-relationship diagrams (text-based)
2. Field definitions with data types
3. Example records
4. How to link compounds to genes/pathways
```

**Pathway/Target Agent:**
```
Document the data models for pathway and target databases that link to traditional medicine:
- Reactome: Pathway schema, entity types
- UniProt: Protein/target data model
- PubChem: Compound schema, bioactivity structure
- KEGG PATHWAY: KGML format, entity relationships
- STRING: Protein interaction schema

For each, provide:
1. Entity-relationship diagrams (text-based)
2. Field definitions with data types
3. ID mapping fields (how to link compounds→targets→genes→pathways)
```

**Output Files Created:**
- `models/tcm-data-models.md` (830 lines)
- `models/ayurveda-data-models.md` (736 lines)
- `models/kampo-data-models.md` (890 lines)
- `models/western-herbal-data-models.md` (1,214 lines)
- `models/pathway-target-data-models.md` (1,613 lines)

---

### 3. Pharmaceutical Research

**User Request:**
> "now we need to research sources of data, research papers, and databases of pharmaceuticals, unless we have already done so? We need to know how we can download the data, and how to use it. Wikipedia should also be a great source of information. We should be able to link compounds and formulations to genes and biochemical pathways. swarm"

**Agent Prompts Used:**

**Wikipedia/Wikidata Agent:**
```
Research Wikipedia and Wikidata as sources for pharmaceutical data:
- Wikidata pharmaceutical entities (drugs, compounds, mechanisms)
- Wikipedia drug infoboxes and their structured data
- SPARQL queries for drug-gene-pathway relationships
- Bulk download options for Wikidata
- How to extract and use this data programmatically

Document access methods, data structures, and linking to genes/pathways.
```

**PharmGKB/CPIC Agent:**
```
Document complete data models for pharmacogenomics databases:
- PharmGKB: API schema, variant annotations, drug-gene relationships
- CPIC: Guideline data structure, dosing recommendations
- Clinical annotation formats

For each provide entity relationships, field definitions, example records, and download methods.
```

**DrugBank/ChEMBL Agent:**
```
Document complete data models for drug databases:
- DrugBank: XML schema, drug-target-enzyme relationships (note: requires license for bulk)
- ChEMBL: SQL schema, bioactivity data model, compound-target-assay structure
- PubChem BioAssay integration

For each provide entity relationships, field definitions, and how to link to genes/pathways.
```

**Drug-Target Interaction Agent:**
```
Document data models for drug-target interaction databases:
- DGIdb: Interaction schema, source integration
- Open Targets: Evidence data model, association scores
- BindingDB: Binding affinity schema
- STRING: Protein-protein interaction model

For each provide entity relationships and pathway linking methods.
```

**FDA/Regulatory Agent:**
```
Document data models for regulatory drug databases:
- OpenFDA: API schema, adverse events, drug labels
- DailyMed: SPL format, drug label structure
- RxNorm: Concept schema, drug naming relationships
- NDC Directory: Product data model

For each provide entity relationships, access methods, and field definitions.
```

**Output Files Created:**
- `interventions-wikipedia-wikidata.md` (844 lines)
- `models/pharmgkb-cpic-data-models.md` (1,051 lines)
- `models/drugbank-chembl-data-models.md` (1,834 lines)
- `models/drug-target-interaction-models.md` (1,450 lines)
- `models/fda-regulatory-data-models.md` (1,709 lines)

---

### 4. Bulk Download Documentation

**User Request:**
> "we need to download wikipedia first, consider that iterating over all entries using their SPARQL API is considered bad form and bad API usage. This goes for any other data source, as well. Consider that we need to process large amount of information from each source"

**Key Constraint Established:** Bulk downloads required, not API iteration.

**Agent Prompts Used:**

**Wikidata/Wikipedia Downloads Agent:**
```
Document how to bulk download and process Wikidata/Wikipedia data:
- Wikidata JSON dumps (location, size, format)
- Wikipedia dump formats
- Streaming parsers for large files (ijson for JSON, lxml.iterparse for XML)
- Filtering strategies during parsing
- Storage requirements and processing pipeline

CRITICAL: Do NOT document SPARQL iteration - only bulk download methods.
Provide code examples for streaming large files.
```

**Traditional Medicine Downloads Agent:**
```
Document bulk download methods for traditional medicine databases:
- TCMBank bulk exports
- BATMAN-TCM data downloads
- LOTUS Natural Products
- IMPPAT downloads
- KampoDB exports

For each: download URLs, file formats, sizes, and processing code.
CRITICAL: Bulk downloads only, not API iteration.
```

**Pharmaceutical Downloads Agent:**
```
Document bulk download methods for pharmaceutical databases:
- ChEMBL FTP downloads
- DrugBank (if available without license)
- OpenFDA bulk downloads
- PubChem FTP structure

For each: download URLs, file formats, sizes, and processing code.
CRITICAL: Bulk downloads only, not API iteration.
```

**Pathway/Target Downloads Agent:**
```
Document bulk download methods for pathway and target databases:
- Reactome bulk downloads
- UniProt FTP
- STRING downloads
- GO annotations bulk
- KEGG (if bulk available)

For each: download URLs, file formats, sizes, and processing code.
CRITICAL: Bulk downloads only, not API iteration.
```

**Data Processing Pipeline Agent:**
```
Design a complete data processing pipeline for Gene Platform:
- Storage architecture for bulk data
- Streaming parser implementations
- ETL pipeline stages
- Database loading strategies
- Storage size estimates
- Hardware requirements

Create a comprehensive pipeline document with code examples.
```

**Output Files Created:**
- `bulk-downloads/wikidata-wikipedia-downloads.md` (1,220 lines)
- `bulk-downloads/traditional-medicine-downloads.md` (688 lines)
- `bulk-downloads/pharmaceutical-downloads.md` (696 lines)
- `bulk-downloads/pathway-target-downloads.md` (794 lines)
- `bulk-downloads/data-processing-pipeline.md` (2,059 lines)

---

### 5. Comprehensive Wikidata Coverage

**User Request:**
> "we should consider wikidata for all data source categories: Kampo, TCM, ayurveda, western herbal, supplements, minerals, vitamins, and biochemical processes. In short, everything"

**Agent Prompts Used:**

**Traditional Medicine in Wikidata Agent:**
```
Research and document Wikidata coverage for traditional medicine:
- TCM herbs, formulas, and concepts (Q-IDs, properties)
- Kampo formulas and crude drugs in Wikidata
- Ayurveda concepts, plants, formulations
- Property codes (P-codes) for traditional medicine relationships

Document Q-IDs, SPARQL queries for bulk extraction, and linking to genes/pathways.
```

**Supplements/Herbal in Wikidata Agent:**
```
Research and document Wikidata coverage for supplements and herbal:
- Vitamins (all forms, Q-IDs)
- Minerals and trace elements
- Amino acids
- Western herbs and botanicals
- Dietary supplements

Document Q-IDs, property relationships, and extraction methods.
```

**Biochemistry in Wikidata Agent:**
```
Research and document Wikidata coverage for biochemistry:
- Genes (human and model organisms)
- Proteins and enzymes
- Metabolic pathways
- Biological processes
- Biochemical compounds

Document Q-IDs, property codes, and how to link compounds→genes→pathways.
```

**Complete Extraction Guide Agent:**
```
Create a master Wikidata extraction guide for Gene Platform:
- Complete list of relevant Q-IDs across all categories
- Property codes (P-codes) for relationships
- Python code for bulk extraction from Wikidata dumps
- Filtering and processing strategies
- Output schema for Gene Platform database

This should be the definitive guide for extracting ALL relevant data from Wikidata.
```

**Output Files Created:**
- `wikidata/wikidata-traditional-medicine.md` (789 lines)
- `wikidata/wikidata-supplements-herbal.md` (518 lines)
- `wikidata/wikidata-biochemistry.md` (906 lines)
- `wikidata/wikidata-complete-extraction.md` (2,557 lines)

---

### 6. Biochemical Pathways Research

**User Request:**
> "now we need to research sources of data, research papers, and databases of biochemical pathways and biological processes, unless we have already done so? We need to know how we can download the data, and how to use it. Wikipedia should also be a great source of information. We should be able to link compounds and formulations to genes and biochemical pathways. swarm"

**Agent Prompts Used:**

**Pathway Databases Agent:**
```
Research and document pathway databases for Gene Platform:
- Reactome: Coverage, data model, bulk downloads, API
- KEGG PATHWAY: Structure, licensing, alternatives
- WikiPathways: CC0 pathways, GPML format
- PathBank: Human pathways
- BioCyc/MetaCyc: Metabolic pathways

For each: coverage statistics, download methods, data formats, and linking strategies.
```

**Biological Processes Agent:**
```
Research and document biological process databases:
- Gene Ontology (GO): Biological process terms, annotations
- Rhea: Biochemical reactions
- BRENDA: Enzyme data
- ChEBI: Chemical entities of biological interest
- IntEnz: Enzyme nomenclature

For each: coverage, download methods, and linking to genes/compounds.
```

**Compound-Pathway Linking Agent:**
```
Create a comprehensive guide for linking compounds to pathways:
- ID mapping strategies (PubChem CID → UniProt → Reactome)
- Cross-reference databases
- Mapping hubs and bridge databases
- Code examples for complete linking pipelines
- Handling ambiguous mappings

This should be the definitive linking guide for Gene Platform.
```

**Pathway Formats Agent:**
```
Document pathway data formats in detail:
- BioPAX: OWL structure, levels, parsing
- GPML: WikiPathways format, XML schema
- KGML: KEGG format
- SBML: Systems biology format
- OBO/GAF: GO annotation format
- CX/CX2: Cytoscape exchange format

For each: structure, parsing code, and conversion between formats.
```

**Wikidata Pathways Agent:**
```
Document Wikidata coverage for pathways and biological processes:
- Pathway entities (Q-IDs)
- Biological process properties
- Gene-pathway relationships in Wikidata
- SPARQL queries for pathway data
- Federated queries linking Wikidata to Reactome/UniProt

Focus on bulk extraction, not API iteration.
```

**Output Files Created:**
- `pathways-databases.md` (797 lines)
- `biological-processes-databases.md` (1,650 lines)
- `compound-pathway-linking.md` (2,360 lines)
- `models/pathway-formats-detailed.md` (3,044 lines)
- `wikidata/wikidata-pathways-detailed.md` (2,293 lines)

---

## Stored Context (Memory)

### Research Papers Context
```json
{
  "purpose": "Research papers needed for: 1) RAG system for AI chat about user genetics, 2) Evidence linking for SNP interpretations, 3) Supporting recommendations with citations, 4) Building knowledge base with evidence levels. Key requirement: Must be public/accessible data only - no commercial licenses or partnerships. Platform focus: genetics, SNPs, biochemical pathways, supplements, autoimmune conditions."
}
```

### Interventions Context
```json
{
  "purpose": "Research all data sources for interventions (treatments, compounds, herbs, drugs) across: 1) Kampo (Japanese traditional medicine), 2) TCM (Traditional Chinese Medicine), 3) Ayurveda (Indian traditional medicine), 4) Western herbal medicine, 5) Pharmaceuticals. For each source identify: coverage, access methods (API/bulk/scrape), data structures/schemas, licensing, and whether full content or abstracts only are available."
}
```

---

## Key Constraints Established

1. **Public Data Only**: No commercial licenses, no partnerships required
2. **Bulk Downloads**: Never iterate via APIs - always bulk download and process locally
3. **Compound-Gene-Pathway Linking**: All data must be linkable through ID mapping
4. **Evidence-Based**: Research papers needed for citations and evidence levels
5. **Multi-Modality**: Cover pharma + supplements + TCM + Ayurveda + Kampo + Western herbal

---

## File Inventory

### Root Research Files
| File | Lines | Description |
|------|-------|-------------|
| `DATA-INTEGRATION-GUIDE.md` | ~500 | Consolidated integration guide |
| `pathways-databases.md` | 797 | 10 pathway database sources |
| `biological-processes-databases.md` | 1,650 | GO, Rhea, BRENDA, ChEBI |
| `compound-pathway-linking.md` | 2,360 | Complete linking guide |
| `interventions-wikipedia-wikidata.md` | 844 | Wikipedia/Wikidata pharma |
| `interventions-tcm.md` | 714 | TCM databases |
| `interventions-ayurveda.md` | 833 | Ayurveda databases |
| `interventions-kampo.md` | 1,058 | Kampo databases |
| `interventions-western-herbal.md` | 813 | Western herbal databases |
| `interventions-pharma.md` | 1,277 | Pharmaceutical databases |

### Models Directory
| File | Lines | Description |
|------|-------|-------------|
| `tcm-data-models.md` | 830 | BATMAN-TCM, TCMBank schemas |
| `ayurveda-data-models.md` | 736 | IMPPAT, NPACT schemas |
| `kampo-data-models.md` | 890 | KampoDB 4-layer model |
| `western-herbal-data-models.md` | 1,214 | DSLD, Dr. Duke's schemas |
| `pathway-target-data-models.md` | 1,613 | Reactome, UniProt schemas |
| `pharmgkb-cpic-data-models.md` | 1,051 | PharmGKB, CPIC schemas |
| `drugbank-chembl-data-models.md` | 1,834 | DrugBank, ChEMBL schemas |
| `drug-target-interaction-models.md` | 1,450 | DGIdb, Open Targets |
| `fda-regulatory-data-models.md` | 1,709 | OpenFDA, DailyMed, RxNorm |
| `pathway-formats-detailed.md` | 3,044 | BioPAX, GPML, KGML, SBML |

### Bulk Downloads Directory
| File | Lines | Description |
|------|-------|-------------|
| `wikidata-wikipedia-downloads.md` | 1,220 | JSON dumps, streaming parsers |
| `traditional-medicine-downloads.md` | 688 | TCMBank, LOTUS downloads |
| `pharmaceutical-downloads.md` | 696 | ChEMBL, OpenFDA bulk |
| `pathway-target-downloads.md` | 794 | Reactome, UniProt FTP |
| `data-processing-pipeline.md` | 2,059 | Complete ETL pipeline |

### Wikidata Directory
| File | Lines | Description |
|------|-------|-------------|
| `wikidata-traditional-medicine.md` | 789 | TCM, Kampo, Ayurveda Q-IDs |
| `wikidata-supplements-herbal.md` | 518 | Vitamins, minerals, herbs |
| `wikidata-biochemistry.md` | 906 | Genes, proteins, pathways |
| `wikidata-complete-extraction.md` | 2,557 | Master extraction guide |
| `wikidata-pathways-detailed.md` | 2,293 | Pathway SPARQL, federation |

---

## Total Statistics

- **Files**: 50
- **Lines**: 51,009
- **Categories Covered**: 6 (Traditional Medicine, Pharmaceuticals, Supplements, Pathways, Biological Processes, Wikidata)
- **Databases Documented**: 40+
- **Research Date**: January 17-18, 2026

---

*Document created: January 18, 2026*
*Purpose: Historical record of research prompts for Gene Platform*
