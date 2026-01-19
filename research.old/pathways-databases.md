# Biochemical Pathway Databases - Comprehensive Research

**Research Date:** January 2026
**Purpose:** Data sources for genetics/health knowledge base platform
**Status:** Comprehensive review of all major biochemical pathway databases

---

## Executive Summary

Biochemical pathway databases are essential for understanding how genetic variants affect biological processes and how interventions (drugs, supplements, natural products) modulate cellular functions. This document catalogs **10 major pathway databases** covering metabolic reactions, signaling cascades, disease pathways, and drug mechanisms.

### Key Findings

1. **Largest Databases by Pathway Count:**
   - KEGG PATHWAY: 559 reference pathways across all organisms
   - WikiPathways: 3,000+ community-curated pathways
   - Reactome: 2,700+ expert-curated human pathways
   - MetaCyc: 2,937 pathways from 3,295 organisms

2. **Best for Open Access:**
   - Reactome (CC BY 4.0) - fully open, commercial use allowed
   - WikiPathways (CC0) - public domain, no restrictions
   - Pathway Commons (varies by source) - aggregates open data

3. **Best for Drug/Compound Integration:**
   - SMPDB: Drug metabolism and disease pathways linked to DrugBank/HMDB
   - PharmGKB Pathways: Pharmacokinetic/pharmacodynamic focus
   - KEGG: COMPOUND and DRUG modules (API restrictions apply)

4. **Best API Access:**
   - Reactome: REST API, GraphQL, Neo4j graph database
   - Pathway Commons: REST API, BioPAX export
   - NDEx: CX format API for network exchange
   - WikiPathways: REST API, SPARQL endpoint

5. **Licensing Considerations:**
   - KEGG requires academic subscription for bulk FTP access
   - MetaCyc/BioCyc: Free academic, paid commercial
   - Most others: Open access with attribution

---

## Database Catalog

### 1. Reactome

| Attribute | Details |
|-----------|---------|
| **URL** | https://reactome.org |
| **Maintainer** | Ontario Institute for Cancer Research (OICR), EMBL-EBI, NYU Langone Health |
| **Content** | Curated biological pathways including metabolic, signaling, disease, and cellular processes |
| **Coverage** | 2,712 human pathways, 13,872 reactions, 11,196 proteins, 1,925 small molecules (Release 90, December 2024) |
| **Species** | 24 species with human as primary, others via orthology projection |
| **Access Method** | Web interface; REST API; GraphQL API; Neo4j graph database; Bulk download |
| **Data Format** | BioPAX Level 3, SBML, SBGN-ML, PSI-MITAB, JSON |
| **API Documentation** | https://reactome.org/dev |
| **Download** | https://reactome.org/download-data |
| **Schema** | Event hierarchy: Pathway -> Reaction -> PhysicalEntity (with compartment localization) |
| **Licensing** | CC BY 4.0 - allows commercial use with attribution |
| **Full Content** | Yes - complete pathway diagrams, reaction mechanisms, literature references |
| **Gene/Protein Links** | Yes - UniProt IDs, Ensembl gene IDs, NCBI Gene IDs |
| **Compound-Pathway Links** | ChEBI IDs for small molecules; reactions annotated with substrates/products/catalysts |
| **Updates** | Quarterly releases |

**Key Features:**
- **Pathway Browser**: Interactive visualization with zoom to molecular detail
- **Analysis Service**: Pathway enrichment, expression analysis, species comparison
- **Icon Library**: Standardized pathway notation (SBGN)
- **Interactors**: Protein-protein interaction overlays from IntAct
- **Literature Evidence**: Each reaction annotated with PubMed citations

**Compound Integration:**
- Small molecules identified by ChEBI IDs
- Drug entities linked to their targets and mechanisms
- Metabolites positioned within reaction networks
- Example: "Metabolism of amino acids and derivatives" includes 500+ compounds

**API Example:**
```
GET https://reactome.org/ContentService/data/pathway/R-HSA-1430728
```

**Notes:** Gold standard for curated human pathways. Preferred over KEGG for commercial applications due to open licensing. Neo4j option excellent for network analysis.

---

### 2. KEGG PATHWAY

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Maintainer** | Kyoto University Bioinformatics Center (Kanehisa Laboratories) |
| **Content** | Molecular interaction, reaction, and relation networks |
| **Coverage** | 559 reference pathways; 700,000+ organism-specific pathways |
| **Categories** | Metabolism (138), Genetic Information (25), Environmental Information (29), Cellular Processes (43), Organismal Systems (102), Human Diseases (99), Drug Development (23) |
| **Access Method** | Web interface; REST API (limited); FTP (subscription required) |
| **Data Format** | KGML (KEGG Markup Language), PNG/SVG images |
| **API Documentation** | https://www.kegg.jp/kegg/rest/keggapi.html |
| **Schema** | Pathways (map) -> Entries (gene, compound, enzyme) -> Relations/Reactions |
| **Licensing** | Web access free; FTP bulk download requires academic subscription |
| **Full Content** | Limited - API provides pathway lists and basic info; bulk data requires subscription |
| **Gene/Protein Links** | KEGG Gene IDs cross-referenced to NCBI, UniProt |
| **Compound-Pathway Links** | KEGG COMPOUND IDs (C numbers), integrated with KEGG DRUG |
| **Updates** | Continuous updates |

**Key Features:**
- **KEGG COMPOUND**: 19,430 metabolites and small molecules (C numbers)
- **KEGG DRUG**: 12,265 drugs with therapeutic categories (D numbers)
- **KEGG ENZYME**: EC number-based enzyme classification
- **KEGG REACTION**: 12,049 biochemical reactions (R numbers)
- **KEGG GLYCAN**: Carbohydrate structures
- **BRITE Hierarchies**: Functional classification trees

**KGML Format:**
```xml
<pathway name="path:hsa00010" org="hsa" number="00010" title="Glycolysis / Gluconeogenesis">
  <entry id="1" name="hsa:2821" type="gene">
    <graphics name="GPI" .../>
  </entry>
  <relation entry1="1" entry2="2" type="PPrel">
    <subtype name="activation" value="-->"/>
  </relation>
</pathway>
```

**Compound Integration:**
- C numbers link metabolites to pathways
- D numbers link drugs to targets and pathways
- Enzyme-compound relationships via EC numbers
- Example: C00031 (D-Glucose) appears in 50+ pathways

**API Limitations:**
- 10 requests per second rate limit
- No bulk download without subscription
- Commercial use requires license agreement

**Notes:** Most comprehensive pathway-compound integration but licensing restrictions limit bulk data use. For open-access alternatives, use Reactome + WikiPathways.

---

### 3. WikiPathways

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.wikipathways.org |
| **Maintainer** | Gladstone Institutes, Maastricht University (community-driven) |
| **Content** | Community-curated biological pathways |
| **Coverage** | 3,100+ pathways across 33 species |
| **Species Focus** | Human (955), Mouse (259), Rat (155), Zebrafish (76), and 29 others |
| **Access Method** | Web interface; REST API; SPARQL endpoint; GitHub releases |
| **Data Format** | GPML (Graphical Pathway Markup Language), BioPAX, JSON |
| **API Documentation** | https://webservice.wikipathways.org/index |
| **Download** | Monthly releases on GitHub and Zenodo |
| **Schema** | Pathway -> DataNodes (GeneProduct, Metabolite, Protein) -> Interactions |
| **Licensing** | CC0 1.0 Universal (Public Domain Dedication) |
| **Full Content** | Yes - complete pathway graphics, data nodes, interactions, literature |
| **Gene/Protein Links** | Ensembl, Entrez Gene, UniProt, HGNC |
| **Compound-Pathway Links** | ChEBI, HMDB, PubChem, KEGG Compound IDs |
| **Updates** | Monthly releases; continuous community edits |

**Key Features:**
- **Community Curation**: Wikipedia-style open editing model
- **Pathway Widget**: Embeddable pathway viewer
- **PathVisio Integration**: Desktop editing tool
- **Curation Tags**: Quality indicators (Featured, Stub, etc.)
- **Analysis Portal**: Enrichment analysis integration with g:Profiler, Enrichr
- **Collections**: Curated pathway sets (COVID-19, Rare Diseases, etc.)

**GPML Format:**
```xml
<Pathway Name="Cholesterol biosynthesis" Organism="Homo sapiens">
  <DataNode TextLabel="HMGCR" GraphId="abc123" Type="GeneProduct">
    <Xref Database="Ensembl" ID="ENSG00000113161"/>
  </DataNode>
  <DataNode TextLabel="Cholesterol" GraphId="def456" Type="Metabolite">
    <Xref Database="ChEBI" ID="CHEBI:16113"/>
  </DataNode>
  <Interaction>
    <Graphics>...</Graphics>
    <Xref Database="" ID=""/>
  </Interaction>
</Pathway>
```

**Compound Integration:**
- Metabolites identified by ChEBI, HMDB, KEGG, or PubChem IDs
- Drug nodes can reference DrugBank
- Xref database links enable cross-database queries
- BridgeDb identifier mapping service

**Zenodo Archive:**
- DOI: 10.5281/zenodo.10651908 (example, version-specific)
- Complete GPML files for all pathways
- RDF/OWL exports available

**Notes:** Best open-access pathway database. CC0 license allows unrestricted use including commercial. Excellent for integration with other resources.

---

### 4. MetaCyc / BioCyc

| Attribute | Details |
|-----------|---------|
| **URL** | https://metacyc.org (MetaCyc), https://biocyc.org (BioCyc) |
| **Maintainer** | SRI International |
| **Content** | Experimentally elucidated metabolic pathways from scientific literature |
| **Coverage** | MetaCyc: 2,937 pathways, 18,124 reactions, 18,400 compounds from 3,295 organisms |
| **BioCyc Collection** | 20,000+ organism-specific Pathway/Genome Databases (PGDBs) including HumanCyc, EcoCyc |
| **Access Method** | Web interface; Pathway Tools software; SmartTables API; FTP (subscription tiers) |
| **Data Format** | BioPAX, SBML, Pathway Tools flat files |
| **API Documentation** | https://biocyc.org/web-services.shtml |
| **Schema** | Pathways -> Reactions -> Compounds + Enzymes -> Genes |
| **Licensing** | Tier 1 (academic free for MetaCyc/EcoCyc); Tier 2/3 (subscription); Commercial (paid license) |
| **Full Content** | Tier-dependent - full access requires subscription |
| **Gene/Protein Links** | UniProt, NCBI Gene, species-specific identifiers |
| **Compound-Pathway Links** | MetaCyc compound IDs, cross-referenced to ChEBI, KEGG, PubChem |
| **Updates** | Bi-annual releases (Version 28.0, June 2024) |

**Key Features:**
- **Pathway Tools Software**: Desktop application for pathway analysis and visualization
- **SmartTables**: Programmatic data access and analysis
- **Comparative Analysis**: Cross-species pathway comparisons
- **Atom Mappings**: Track individual atoms through reactions
- **Gibbs Energy**: Thermodynamic predictions for reactions
- **Transport Reactions**: Membrane transport pathways

**Organism-Specific Databases:**
- **HumanCyc**: Human metabolic pathways (Tier 2)
- **EcoCyc**: E. coli reference database (Tier 1, free)
- **PlantCyc**: Plant metabolism (Tier 2)
- **YeastCyc**: S. cerevisiae (Tier 2)

**Compound Integration:**
- Comprehensive metabolite database with structures
- Reaction stoichiometry and atom balancing
- Enzyme-compound kinetic parameters where available
- Cofactor and coenzyme tracking

**Licensing Tiers:**
| Tier | Access | Cost |
|------|--------|------|
| Tier 1 | MetaCyc, EcoCyc, Web services | Free (academic) |
| Tier 2 | Individual organism PGDBs | Subscription |
| Tier 3 | Full BioCyc collection | Subscription |
| Commercial | All databases | Paid license |

**Notes:** Highest quality metabolic pathway curation. MetaCyc free for academic use. Pathway Tools software powerful but has learning curve.

---

### 5. Pathway Commons

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pathwaycommons.org |
| **Maintainer** | University of Toronto, Memorial Sloan Kettering Cancer Center |
| **Content** | Integrated biological pathway and interaction data from multiple sources |
| **Coverage** | 5,772 pathways, 2.4 million interactions, 22 data sources |
| **Data Sources** | Reactome, KEGG, WikiPathways, HumanCyc, NCI-PID, PANTHER, PhosphoSitePlus, HPRD, and others |
| **Access Method** | Web interface; REST API; Bulk download |
| **Data Format** | BioPAX Level 3, SIF (Simple Interaction Format), SBGN-ML |
| **API Documentation** | https://www.pathwaycommons.org/pc2/home |
| **Download** | https://www.pathwaycommons.org/archives/PC2/ |
| **Schema** | Unified BioPAX model: Pathway -> Interaction -> PhysicalEntity |
| **Licensing** | Varies by source - check individual data source licenses |
| **Full Content** | Yes - complete integrated network |
| **Gene/Protein Links** | UniProt, HGNC, NCBI Gene |
| **Compound-Pathway Links** | ChEBI IDs, KEGG Compound |
| **Updates** | Version 12 (March 2020); periodic updates |

**Key Features:**
- **Data Integration**: Unified access to multiple pathway databases
- **Network Queries**: Neighborhood, paths-between, common-stream queries
- **Paxtools Library**: Java library for BioPAX manipulation
- **CyPath2 App**: Cytoscape integration
- **SIF Network**: Simplified binary interaction format

**Integrated Data Sources:**
| Source | Type | Interactions |
|--------|------|--------------|
| Reactome | Pathways | 261,460 |
| KEGG | Pathways | 142,536 |
| PhosphoSitePlus | PTM | 219,000 |
| HumanCyc | Metabolic | 16,000 |
| HPRD | PPI | 39,240 |
| BioGRID | PPI | 645,000 |
| PANTHER | Pathways | 14,800 |
| NCI-PID | Cancer | 23,475 |

**SIF Format:**
```
TP53	controls-state-change-of	MDM2
BRCA1	interacts-with	BRCA2
EGFR	controls-phosphorylation-of	MAPK1
```

**Compound Integration:**
- Small molecules from ChEBI ontology
- Drug-target interactions from source databases
- Metabolic reactions with compound substrates/products
- Chemical modification annotations

**Notes:** Excellent for cross-database queries. License varies by source - commercial users should verify individual source licenses.

---

### 6. SMPDB (Small Molecule Pathway Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://smpdb.ca |
| **Maintainer** | University of Alberta (Wishart Research Group) |
| **Content** | Small molecule pathways including metabolic, drug action, drug metabolism, disease pathways |
| **Coverage** | 99,900+ pathways total: 1,945 normal metabolic, 1,556 drug metabolism, 89,596 metabolic disease |
| **Species** | Human-focused |
| **Access Method** | Web interface; Bulk download |
| **Data Format** | SMPDB XML, SVG/PNG pathway images, BioPAX |
| **Download** | https://smpdb.ca/downloads |
| **Schema** | Pathway -> Small Molecules + Proteins -> Reactions |
| **Licensing** | CC BY-NC 4.0 (non-commercial with attribution) |
| **Full Content** | Yes - pathway diagrams, compound lists, protein participants |
| **Gene/Protein Links** | UniProt, HMDB protein IDs |
| **Compound-Pathway Links** | HMDB IDs, DrugBank IDs, PubChem CIDs |
| **Updates** | Version 2.0 (2020); periodic updates |

**Key Features:**
- **HMDB Integration**: Direct links to Human Metabolome Database
- **DrugBank Integration**: Drug mechanism pathways linked to DrugBank
- **Disease Pathways**: 89,000+ inborn error of metabolism pathways
- **Visualization**: Detailed KEGG-style pathway diagrams
- **Subcellular Localization**: Compartment annotations for reactions
- **Tissue Expression**: Pathway tissue distributions

**Pathway Categories:**
| Category | Count | Description |
|----------|-------|-------------|
| Metabolic | 1,945 | Normal human metabolism |
| Drug Metabolism | 1,556 | Pharmacokinetics of 800+ drugs |
| Drug Action | 6,803 | Pharmacodynamics and mechanisms |
| Disease | 89,596 | Metabolic disease pathways |

**Compound Integration:**
- Primary metabolites linked via HMDB IDs
- Drugs linked via DrugBank IDs
- Chemical structures available
- Concentration data from HMDB
- Example: Acetaminophen pathway includes all metabolites and enzymes involved

**Download Files:**
- smpdb_pathways.csv - pathway metadata
- smpdb_metabolites.csv - metabolite-pathway associations
- smpdb_proteins.csv - protein-pathway associations
- Pathway images (SVG/PNG)
- BioPAX exports

**Notes:** Best database for drug metabolism pathways. Essential for pharmacokinetic modeling. Excellent HMDB/DrugBank integration.

---

### 7. PharmGKB Pathways

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pharmgkb.org/pathways |
| **Maintainer** | Stanford University |
| **Content** | Pharmacokinetic (PK) and pharmacodynamic (PD) pathways for drugs |
| **Coverage** | 200+ curated drug pathways |
| **Focus** | Pharmacogenomics - genetic variants affecting drug response |
| **Access Method** | Web interface; Bulk download (registration required) |
| **Data Format** | GPML, PNG images, TSV |
| **Download** | https://www.pharmgkb.org/downloads |
| **Schema** | Drug -> Pathways -> Genes/Variants -> Clinical Outcomes |
| **Licensing** | CC BY-SA 4.0 |
| **Full Content** | Yes - pathway diagrams, gene annotations, variant effects |
| **Gene/Protein Links** | HGNC, Ensembl, UniProt |
| **Compound-Pathway Links** | PharmGKB drug IDs, cross-referenced to DrugBank, PubChem, ChEMBL |
| **Updates** | Continuous curation |

**Key Features:**
- **PK Pathways**: Absorption, distribution, metabolism, excretion
- **PD Pathways**: Drug target interactions and downstream effects
- **Variant Annotations**: Clinical significance of genetic variants
- **Clinical Guidelines**: CPIC and DPWG dosing guidelines
- **Drug Labels**: FDA pharmacogenomic labeling information

**Pathway Types:**
| Type | Description | Example |
|------|-------------|---------|
| PK | Drug metabolism enzymes | Warfarin VKORC1/CYP2C9 pathway |
| PD | Drug mechanism of action | Statins pharmacodynamics |
| Disease | Drug-disease-gene relationships | Cardiovascular pathways |

**Compound Integration:**
- Drugs annotated with:
  - Metabolizing enzymes (CYP450, UGT, etc.)
  - Transporters (OATP, P-gp, etc.)
  - Drug targets
- Variant-drug-response associations
- Example: Clopidogrel pathway shows CYP2C19 variants affecting efficacy

**Key Genes by Drug Class:**
| Drug Class | Key Genes |
|------------|-----------|
| Warfarin | VKORC1, CYP2C9, CYP4F2 |
| Clopidogrel | CYP2C19, PON1 |
| Statins | SLCO1B1, ABCG2 |
| SSRIs | CYP2D6, CYP2C19 |
| Opioids | CYP2D6, OPRM1 |

**Notes:** Essential for pharmacogenomics. Best source for clinically actionable variant-drug interactions. Excellent integration with clinical guidelines.

---

### 8. PANTHER Pathways

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.pantherdb.org/pathway/ |
| **Maintainer** | University of Southern California |
| **Content** | Signaling and metabolic pathways |
| **Coverage** | 177 curated pathways, 7,180 pathway components |
| **Species** | Human-focused with ortholog mappings to 142 organisms |
| **Access Method** | Web interface; FTP download |
| **Data Format** | SBML, BioPAX, PSI-MI XML 2.5 |
| **Download** | ftp://ftp.pantherdb.org/pathway/ |
| **Schema** | Pathway -> Pathway Components -> Genes -> GO Annotations |
| **Licensing** | CC BY 4.0 |
| **Full Content** | Yes - pathway diagrams, gene lists, functional annotations |
| **Gene/Protein Links** | UniProt, Ensembl, NCBI Gene |
| **Compound-Pathway Links** | Limited - primarily gene/protein focused |
| **Updates** | Version 18.0 (February 2023); annual updates |

**Key Features:**
- **Classification System**: Hierarchical protein family classification
- **GO Integration**: Gene Ontology annotations for all genes
- **PANTHER Slim**: Simplified GO subsets for analysis
- **Overrepresentation Analysis**: Statistical enrichment testing
- **Gene List Analysis**: Functional classification of gene sets

**Pathway Categories:**
| Category | Count | Examples |
|----------|-------|----------|
| Signaling | 75 | Wnt, Notch, TGF-beta |
| Metabolic | 40 | Glycolysis, Krebs cycle |
| Disease | 30 | Alzheimer, Huntington |
| Other | 32 | Apoptosis, cell cycle |

**PANTHER Analysis Tools:**
- Gene List Analysis: Input genes, get pathway enrichment
- Pathway Diagrams: Interactive visualization
- Species Comparison: Orthology-based pathway comparison
- Protein Classification: Family and subfamily assignments

**Compound Integration:**
- Metabolite placeholders in pathway diagrams
- No direct compound database links
- Focus is on enzyme/gene annotation

**Notes:** Strong for signaling pathways and GO integration. Use alongside compound-focused databases like SMPDB for comprehensive coverage.

---

### 9. NCI Pathway Interaction Database (Archived) / NDEx

| Attribute | Details |
|-----------|---------|
| **Original URL** | https://pid.nci.nih.gov (archived 2016) |
| **Current Archive** | NDEx (Network Data Exchange): https://www.ndexbio.org |
| **Maintainer** | Originally NCI/Nature; now UCSD (NDEx) |
| **Content** | Cancer-related signaling pathways and biomolecular interactions |
| **Coverage** | NCI-PID: 2,697 interactions, 1,273 pathways; NDEx total: 15,000+ networks |
| **Access Method** | NDEx web interface; REST API; CyNDEx Cytoscape app |
| **Data Format** | CX (Cytoscape Exchange), BioPAX (legacy) |
| **API Documentation** | https://home.ndexbio.org/using-the-ndex-server-api/ |
| **Download** | https://www.ndexbio.org/index.html#/ |
| **Schema** | Networks -> Nodes (genes, proteins, compounds) -> Edges (interactions) |
| **Licensing** | Public domain (NCI-PID); varies by network on NDEx |
| **Full Content** | Yes - all NCI-PID data preserved in NDEx |
| **Gene/Protein Links** | UniProt, HGNC, Entrez Gene |
| **Compound-Pathway Links** | ChEBI, PubChem (where annotated) |
| **Updates** | NCI-PID frozen; NDEx continuously updated |

**NCI-PID Content (Archived):**
| Source | Pathways | Interactions |
|--------|----------|--------------|
| NCI-Nature Curated | 212 | 10,896 |
| BioCarta Import | 254 | 6,292 |
| Reactome Import | 807 | 35,551 |

**NDEx Features:**
- **Network Repository**: Share and discover biological networks
- **CX Format**: JSON-based network exchange format
- **DOI Assignment**: Citable network publications
- **Access Control**: Public, private, or group sharing
- **Integrated Search**: Query across all public networks
- **IQuery**: Pathway enrichment analysis service

**CX Format Example:**
```json
{
  "nodes": [
    {"@id": 1, "n": "TP53"},
    {"@id": 2, "n": "MDM2"}
  ],
  "edges": [
    {"@id": 3, "s": 1, "t": 2, "i": "inhibits"}
  ],
  "nodeAttributes": [
    {"po": 1, "n": "type", "v": "protein"}
  ]
}
```

**Key Networks on NDEx:**
- NCI Pathway Interaction Database (complete archive)
- SIGNOR: Signaling Network Open Resource
- STRING PPI networks
- NDEx cancer hallmarks networks
- User-contributed networks

**Notes:** NCI-PID was authoritative for cancer pathways. Content preserved in NDEx along with many other networks. NDEx is emerging as standard for network sharing.

---

### 10. NDEx (Network Data Exchange)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ndexbio.org |
| **Maintainer** | UC San Diego, Cytoscape Consortium |
| **Content** | Repository for biological network models of any type |
| **Coverage** | 15,000+ public networks; millions of nodes and edges |
| **Network Types** | Pathways, PPI, regulatory, signaling, metabolic, disease, drug-target |
| **Access Method** | Web interface; REST API; Python client (ndex2); Cytoscape via CyNDEx |
| **Data Format** | CX (Cytoscape Exchange JSON), CX2 (next generation) |
| **API Documentation** | https://home.ndexbio.org/using-the-ndex-server-api/ |
| **Python Client** | pip install ndex2 |
| **Schema** | Network -> Nodes -> Edges -> Attributes (flexible schema per network) |
| **Licensing** | Per-network licensing (many CC0 or CC BY) |
| **Full Content** | Yes - complete network structures and metadata |
| **Gene/Protein Links** | Flexible - depends on network |
| **Compound-Pathway Links** | Network-dependent; many include ChEBI, DrugBank IDs |
| **Updates** | Continuous - user uploads and curated content |

**Key Features:**
- **Network Publication**: DOI assignment for citation
- **Version Control**: Track network revisions
- **IQuery Service**: Automated pathway enrichment
- **HiView Integration**: Hierarchical network visualization
- **Network Sets**: Organized collections
- **Provenance Tracking**: Network history and lineage

**Notable Network Collections:**
| Collection | Description | Networks |
|------------|-------------|----------|
| NCI-PID | Cancer pathways (archive) | 1,273 |
| SIGNOR | Signaling interactions | 200+ |
| STRING | PPI networks | Multiple |
| Drug-target | Drug mechanism networks | 500+ |
| Disease | Disease-specific networks | 1,000+ |

**Python API Example:**
```python
import ndex2

client = ndex2.client.Ndex2()
networks = client.search_networks('cancer signaling')
network = client.get_network_as_cx_stream(uuid)
```

**CX2 Format:**
- Next-generation format
- Better performance for large networks
- Standardized visual properties
- Enhanced attribute system

**Compound Integration:**
- Network-dependent
- Best networks include:
  - ChEBI IDs for metabolites
  - DrugBank IDs for drugs
  - PubChem CIDs
- Query for drug/compound networks specifically

**Notes:** Platform for sharing any biological network. Excellent for finding specialized pathway networks not in major databases. Strong Cytoscape integration.

---

## Comparative Analysis

### Coverage Comparison

| Database | Pathways | Reactions | Species | Compounds |
|----------|----------|-----------|---------|-----------|
| Reactome | 2,712 | 13,872 | 24 | 1,925 |
| KEGG | 559 ref | 12,049 | 8,000+ | 19,430 |
| WikiPathways | 3,100+ | N/A | 33 | Varies |
| MetaCyc | 2,937 | 18,124 | 3,295 | 18,400 |
| Pathway Commons | 5,772 | 2.4M int. | Human | Varies |
| SMPDB | 99,900 | N/A | Human | Linked |
| PharmGKB | 200+ | N/A | Human | Drug focus |
| PANTHER | 177 | N/A | 142 | Limited |
| NDEx | 15,000+ | Varies | Varies | Varies |

### Access Methods Comparison

| Database | Web UI | REST API | SPARQL | Bulk DL | Cytoscape |
|----------|--------|----------|--------|---------|-----------|
| Reactome | Yes | Yes | Yes | Yes | App |
| KEGG | Yes | Limited | No | Paid | Via API |
| WikiPathways | Yes | Yes | Yes | Yes | App |
| MetaCyc | Yes | Yes | No | Paid | Export |
| Pathway Commons | Yes | Yes | No | Yes | CyPath2 |
| SMPDB | Yes | No | No | Yes | Export |
| PharmGKB | Yes | Limited | No | Yes | Export |
| PANTHER | Yes | No | No | Yes | Export |
| NDEx | Yes | Yes | No | Yes | CyNDEx |

### Licensing Summary

| License Type | Databases |
|--------------|-----------|
| **CC0 (Public Domain)** | WikiPathways |
| **CC BY 4.0** (Commercial OK) | Reactome, PANTHER |
| **CC BY-SA 4.0** (ShareAlike) | PharmGKB |
| **CC BY-NC 4.0** (Non-commercial) | SMPDB |
| **Tiered/Subscription** | KEGG (bulk), MetaCyc/BioCyc |
| **Varies by Source** | Pathway Commons, NDEx |

### Data Format Comparison

| Format | Databases | Description |
|--------|-----------|-------------|
| **BioPAX** | Reactome, WikiPathways, Pathway Commons, MetaCyc | Standard pathway exchange |
| **SBML** | Reactome, PANTHER, MetaCyc | Systems biology models |
| **GPML** | WikiPathways, PharmGKB | Graphical pathway markup |
| **KGML** | KEGG | KEGG-specific XML |
| **CX/CX2** | NDEx | Cytoscape exchange JSON |
| **SIF** | Pathway Commons | Simple interaction format |
| **SBGN-ML** | Reactome | Systems biology graphical notation |

### Compound-Pathway Integration Quality

| Database | Compound IDs | Drug Links | Metabolite Data | Quality |
|----------|--------------|------------|-----------------|---------|
| SMPDB | HMDB, DrugBank | Excellent | Comprehensive | High |
| KEGG | KEGG C/D | Excellent | Comprehensive | High |
| Reactome | ChEBI | Good | Good | High |
| PharmGKB | PharmGKB | Excellent | Limited | High |
| MetaCyc | MetaCyc, ChEBI | Limited | Excellent | High |
| WikiPathways | Mixed | Variable | Variable | Medium |
| Pathway Commons | ChEBI | Aggregated | Aggregated | Medium |
| PANTHER | Limited | Limited | Limited | Low |
| NDEx | Network-dependent | Variable | Variable | Variable |

---

## Recommendations for Integration

### Priority 1: Core Databases (Essential)

1. **Reactome** - Primary pathway database
   - Best overall curation quality
   - CC BY 4.0 allows commercial use
   - Excellent API and Neo4j access
   - Strong compound-reaction annotations

2. **WikiPathways** - Community pathways
   - CC0 license (no restrictions)
   - Community contributions capture emerging biology
   - GPML format well-documented
   - Monthly releases on Zenodo

3. **SMPDB** - Drug metabolism focus
   - Best HMDB/DrugBank integration
   - Drug metabolism pathways essential for pharmacology
   - Small molecule-centric design

### Priority 2: Specialized Databases

4. **PharmGKB Pathways** - Pharmacogenomics
   - Clinical variant-drug relationships
   - Essential for personalized medicine
   - CPIC guideline integration

5. **MetaCyc** (Tier 1 access) - Metabolic reference
   - Highest quality metabolic curation
   - Free for academic use
   - Comprehensive reaction database

6. **NDEx** - Network repository
   - Access to specialized networks
   - NCI-PID archive
   - Community contributions

### Priority 3: Aggregated Access

7. **Pathway Commons** - Integration layer
   - Unified access to multiple sources
   - BioPAX standard format
   - Pre-integrated network

### KEGG Consideration

- Use for reference but not bulk data (license restrictions)
- Alternative: Map KEGG pathway IDs to Reactome equivalents
- Use WikiPathways for similar coverage without restrictions

### Data Harmonization Strategy

1. **Primary Identifiers:**
   - Genes: Ensembl Gene ID, HGNC Symbol
   - Proteins: UniProt ID
   - Compounds: ChEBI ID (preferred), HMDB ID, PubChem CID

2. **Pathway Mapping:**
   - Create cross-reference table: Reactome <-> WikiPathways <-> KEGG
   - Use Pathway Commons for pre-built mappings

3. **Compound-Pathway Links:**
   - Primary: SMPDB (drug metabolism)
   - Secondary: Reactome (general metabolism)
   - Reference: KEGG (via API for specific queries)

4. **Integration Approach:**
   ```
   Gene/Variant -> Protein -> Pathway -> Compound
        |             |           |          |
        v             v           v          v
     Ensembl      UniProt    Reactome    ChEBI/HMDB
   ```

---

## Technical Notes

### API Rate Limits

| Database | Rate Limit | Notes |
|----------|------------|-------|
| Reactome | 100/min | Generous for most uses |
| KEGG | 10/sec | Strict enforcement |
| WikiPathways | Reasonable | No published limit |
| NDEx | Varies | Account-based quotas |
| Pathway Commons | Open | No published limit |

### Recommended Downloads

| Database | Files | Size | Format |
|----------|-------|------|--------|
| Reactome | All pathways | ~500 MB | BioPAX, SBML |
| WikiPathways | Monthly release | ~200 MB | GPML, BioPAX |
| SMPDB | All pathways | ~1 GB | CSV, XML |
| Pathway Commons | Complete | ~2 GB | BioPAX, SIF |
| PharmGKB | Pathways | ~100 MB | TSV, GPML |

### Schema Mapping Example

```sql
-- Unified pathway-compound table
CREATE TABLE pathway_compound (
    pathway_id VARCHAR(50),
    pathway_source ENUM('reactome','wikipathways','smpdb','kegg'),
    compound_id VARCHAR(50),
    compound_id_type ENUM('chebi','hmdb','pubchem','kegg','drugbank'),
    role ENUM('substrate','product','catalyst','inhibitor','modifier'),
    evidence_code VARCHAR(20),
    pubmed_id INTEGER
);
```

---

## References

1. Jassal B, et al. (2020). The reactome pathway knowledgebase. Nucleic Acids Res. 48(D1):D498-D503. https://doi.org/10.1093/nar/gkz1031

2. Kanehisa M, et al. (2023). KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Res. 51(D1):D587-D592. https://doi.org/10.1093/nar/gkac963

3. Martens M, et al. (2021). WikiPathways: connecting communities. Nucleic Acids Res. 49(D1):D613-D621. https://doi.org/10.1093/nar/gkaa1024

4. Caspi R, et al. (2020). The MetaCyc database of metabolic pathways and enzymes. Nucleic Acids Res. 48(D1):D445-D453. https://doi.org/10.1093/nar/gkz862

5. Rodchenkov I, et al. (2020). Pathway Commons 2019 Update. Nucleic Acids Res. 48(D1):D489-D497. https://doi.org/10.1093/nar/gkz946

6. Jewison T, et al. (2014). SMPDB 2.0: Big Improvements to the Small Molecule Pathway Database. Nucleic Acids Res. 42:D478-D484. https://doi.org/10.1093/nar/gkt1067

7. Whirl-Carrillo M, et al. (2021). An Evidence-Based Framework for Evaluating Pharmacogenomics Knowledge. Clin Pharmacol Ther. 110(3):563-572. https://doi.org/10.1002/cpt.2350

8. Mi H, et al. (2021). PANTHER version 16. Nucleic Acids Res. 49(D1):D419-D426. https://doi.org/10.1093/nar/gkaa1106

9. Pratt D, et al. (2015). NDEx, the Network Data Exchange. Cell Syst. 1(4):302-305. https://doi.org/10.1016/j.cels.2015.10.001

---

*Last updated: January 2026*
