# Pathways Primary Databases

**Document ID:** 43-41-PATHWAYS-PRIMARY
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../index.md](./../index.md)

---

## TL;DR

Primary pathway databases provide curated biological pathway data essential for understanding gene-compound-disease relationships. Reactome (CC BY 4.0, 2,712 human pathways) and WikiPathways (CC0, 3,100+ pathways) are recommended as open-access priorities, with KEGG, MetaCyc, SMPDB, PharmGKB Pathways, PANTHER, Pathway Commons, and NDEx providing specialized coverage for metabolism, pharmacogenomics, and network analysis.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary pathway source | Reactome | CC BY 4.0 license allows commercial use; highest curation quality; Neo4j + REST API | Jan 2026 |
| Secondary pathway source | WikiPathways | CC0 public domain; community contributions capture emerging biology | Jan 2026 |
| Drug metabolism pathways | SMPDB | Best HMDB/DrugBank integration; comprehensive PK pathways | Jan 2026 |
| Pharmacogenomics pathways | PharmGKB Pathways | Clinical variant-drug relationships; CPIC guideline integration | Jan 2026 |
| KEGG approach | Reference only | Bulk FTP requires subscription; use Reactome/WikiPathways for open alternatives | Jan 2026 |
| Integration layer | Pathway Commons | Pre-integrated data from multiple sources; BioPAX standard format | Jan 2026 |
| Network repository | NDEx | Access specialized networks; NCI-PID archive preserved | Jan 2026 |
| Primary compound identifier | ChEBI | Used by Reactome, WikiPathways, Pathway Commons; open ontology | Jan 2026 |

---

## Database Catalog

### 1. Reactome

| Field | Value |
|-------|-------|
| **URL** | https://reactome.org |
| **Content** | Curated biological pathways including metabolic, signaling, disease, and cellular processes |
| **Records** | 2,712 human pathways, 13,872 reactions, 11,196 proteins, 1,925 small molecules (Release 90, December 2024) |
| **Species** | 24 species with human as primary; others via orthology projection |
| **License** | CC BY 4.0 (allows commercial use with attribution) |
| **API** | REST API, GraphQL API, Neo4j graph database |
| **API Documentation** | https://reactome.org/dev |
| **Download** | https://reactome.org/download-data |
| **Data Formats** | BioPAX Level 3, SBML, SBGN-ML, PSI-MITAB, JSON |
| **Update Frequency** | Quarterly releases |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB (BioPAX + SBML exports) |

**Key Features:**
- Event hierarchy: Pathway -> Reaction -> PhysicalEntity (with compartment localization)
- Each reaction annotated with PubMed citations
- ChEBI IDs for small molecules
- Protein-protein interaction overlays from IntAct
- Standardized pathway notation (SBGN)

**Compound Integration:**
- Small molecules identified by ChEBI IDs
- Drug entities linked to targets and mechanisms
- Metabolites positioned within reaction networks
- Example: "Metabolism of amino acids and derivatives" includes 500+ compounds

**API Example:**
```
GET https://reactome.org/ContentService/data/pathway/R-HSA-1430728
```

---

### 2. WikiPathways

| Field | Value |
|-------|-------|
| **URL** | https://www.wikipathways.org |
| **Content** | Community-curated biological pathways |
| **Records** | 3,100+ pathways across 33 species |
| **Species Focus** | Human (955), Mouse (259), Rat (155), Zebrafish (76), and 29 others |
| **License** | CC0 1.0 Universal (Public Domain Dedication) |
| **API** | REST API, SPARQL endpoint |
| **API Documentation** | https://webservice.wikipathways.org/index |
| **Download** | Monthly releases on GitHub and Zenodo |
| **Data Formats** | GPML (Graphical Pathway Markup Language), BioPAX, JSON |
| **Update Frequency** | Monthly releases; continuous community edits |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~200 MB (GPML + BioPAX exports) |

**Key Features:**
- Wikipedia-style open editing model
- Curation tags for quality indicators (Featured, Stub, etc.)
- Curated collections (COVID-19, Rare Diseases, etc.)
- BridgeDb identifier mapping service
- PathVisio desktop editing tool

**Compound Integration:**
- Metabolites identified by ChEBI, HMDB, KEGG, or PubChem IDs
- Drug nodes can reference DrugBank
- Xref database links enable cross-database queries

**GPML Format Example:**
```xml
<Pathway Name="Cholesterol biosynthesis" Organism="Homo sapiens">
  <DataNode TextLabel="HMGCR" GraphId="abc123" Type="GeneProduct">
    <Xref Database="Ensembl" ID="ENSG00000113161"/>
  </DataNode>
  <DataNode TextLabel="Cholesterol" GraphId="def456" Type="Metabolite">
    <Xref Database="ChEBI" ID="CHEBI:16113"/>
  </DataNode>
</Pathway>
```

---

### 3. KEGG PATHWAY

| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Content** | Molecular interaction, reaction, and relation networks |
| **Records** | 559 reference pathways; 700,000+ organism-specific pathways |
| **Categories** | Metabolism (138), Genetic Information (25), Environmental Information (29), Cellular Processes (43), Organismal Systems (102), Human Diseases (99), Drug Development (23) |
| **License** | Web access free; FTP bulk download requires academic subscription |
| **API** | REST API (limited, 10 requests/sec rate limit) |
| **API Documentation** | https://www.kegg.jp/kegg/rest/keggapi.html |
| **Data Formats** | KGML (KEGG Markup Language), PNG/SVG images |
| **Update Frequency** | Continuous updates |
| **Priority** | Tier 2 (Reference only due to license) |
| **Storage Estimate** | N/A (subscription required for bulk) |

**Key Features:**
- KEGG COMPOUND: 19,430 metabolites (C numbers)
- KEGG DRUG: 12,265 drugs (D numbers)
- KEGG ENZYME: EC number-based classification
- KEGG REACTION: 12,049 biochemical reactions (R numbers)
- BRITE hierarchies for functional classification

**Compound Integration:**
- C numbers link metabolites to pathways
- D numbers link drugs to targets and pathways
- Example: C00031 (D-Glucose) appears in 50+ pathways

**API Limitations:**
- 10 requests per second rate limit
- No bulk download without subscription
- Commercial use requires license agreement

**Recommendation:** Use for reference queries only; prefer Reactome + WikiPathways for bulk data.

---

### 4. MetaCyc / BioCyc

| Field | Value |
|-------|-------|
| **URL** | https://metacyc.org (MetaCyc), https://biocyc.org (BioCyc) |
| **Content** | Experimentally elucidated metabolic pathways from scientific literature |
| **Records** | MetaCyc: 2,937 pathways, 18,124 reactions, 18,400 compounds from 3,295 organisms |
| **BioCyc Collection** | 20,000+ organism-specific Pathway/Genome Databases (PGDBs) including HumanCyc, EcoCyc |
| **License** | Tier 1 (academic free for MetaCyc/EcoCyc); Tier 2/3 (subscription); Commercial (paid) |
| **API** | SmartTables API, Pathway Tools software |
| **API Documentation** | https://biocyc.org/web-services.shtml |
| **Data Formats** | BioPAX, SBML, Pathway Tools flat files |
| **Update Frequency** | Bi-annual releases (Version 28.0, June 2024) |
| **Priority** | Tier 2 (MetaCyc free tier) |
| **Storage Estimate** | ~300 MB (MetaCyc exports) |

**Licensing Tiers:**
| Tier | Access | Cost |
|------|--------|------|
| Tier 1 | MetaCyc, EcoCyc, Web services | Free (academic) |
| Tier 2 | Individual organism PGDBs | Subscription |
| Tier 3 | Full BioCyc collection | Subscription |
| Commercial | All databases | Paid license |

**Key Features:**
- Atom mappings for tracking atoms through reactions
- Gibbs energy predictions for reactions
- Transport reactions and membrane pathways
- HumanCyc for human-specific metabolism

**Compound Integration:**
- MetaCyc compound IDs cross-referenced to ChEBI, KEGG, PubChem
- Reaction stoichiometry and atom balancing
- Enzyme-compound kinetic parameters where available

---

### 5. SMPDB (Small Molecule Pathway Database)

| Field | Value |
|-------|-------|
| **URL** | https://smpdb.ca |
| **Content** | Small molecule pathways: metabolic, drug action, drug metabolism, disease |
| **Records** | 99,900+ pathways: 1,945 normal metabolic, 1,556 drug metabolism, 6,803 drug action, 89,596 metabolic disease |
| **Species** | Human-focused |
| **License** | CC BY-NC 4.0 (non-commercial with attribution) |
| **API** | Bulk download only (no REST API) |
| **Download** | https://smpdb.ca/downloads |
| **Data Formats** | SMPDB XML, SVG/PNG pathway images, BioPAX |
| **Update Frequency** | Version 2.0 (2020); periodic updates |
| **Priority** | Tier 1 (MVP) for drug metabolism |
| **Storage Estimate** | ~1 GB (CSV + XML + images) |

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

**Download Files:**
- smpdb_pathways.csv - pathway metadata
- smpdb_metabolites.csv - metabolite-pathway associations
- smpdb_proteins.csv - protein-pathway associations
- Pathway images (SVG/PNG)

---

### 6. PharmGKB Pathways

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org/pathways |
| **Content** | Pharmacokinetic (PK) and pharmacodynamic (PD) pathways for drugs |
| **Records** | 200+ curated drug pathways |
| **Focus** | Pharmacogenomics - genetic variants affecting drug response |
| **License** | CC BY-SA 4.0 |
| **API** | Bulk download (registration required) |
| **Download** | https://www.pharmgkb.org/downloads |
| **Data Formats** | GPML, PNG images, TSV |
| **Update Frequency** | Continuous curation |
| **Priority** | Tier 1 (MVP) for pharmacogenomics |
| **Storage Estimate** | ~100 MB (TSV + GPML) |

**Pathway Types:**
| Type | Description | Example |
|------|-------------|---------|
| PK | Drug metabolism enzymes | Warfarin VKORC1/CYP2C9 pathway |
| PD | Drug mechanism of action | Statins pharmacodynamics |
| Disease | Drug-disease-gene relationships | Cardiovascular pathways |

**Key Genes by Drug Class:**
| Drug Class | Key Genes |
|------------|-----------|
| Warfarin | VKORC1, CYP2C9, CYP4F2 |
| Clopidogrel | CYP2C19, PON1 |
| Statins | SLCO1B1, ABCG2 |
| SSRIs | CYP2D6, CYP2C19 |
| Opioids | CYP2D6, OPRM1 |

**Compound Integration:**
- Drugs annotated with metabolizing enzymes (CYP450, UGT, etc.)
- Transporters (OATP, P-gp, etc.)
- Variant-drug-response associations
- Links to CPIC and DPWG clinical guidelines

---

### 7. PANTHER Pathways

| Field | Value |
|-------|-------|
| **URL** | http://www.pantherdb.org/pathway/ |
| **Content** | Signaling and metabolic pathways |
| **Records** | 177 curated pathways, 7,180 pathway components |
| **Species** | Human-focused with ortholog mappings to 142 organisms |
| **License** | CC BY 4.0 |
| **API** | FTP download |
| **Download** | ftp://ftp.pantherdb.org/pathway/ |
| **Data Formats** | SBML, BioPAX, PSI-MI XML 2.5 |
| **Update Frequency** | Version 18.0 (February 2023); annual updates |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

**Pathway Categories:**
| Category | Count | Examples |
|----------|-------|----------|
| Signaling | 75 | Wnt, Notch, TGF-beta |
| Metabolic | 40 | Glycolysis, Krebs cycle |
| Disease | 30 | Alzheimer, Huntington |
| Other | 32 | Apoptosis, cell cycle |

**Key Features:**
- Hierarchical protein family classification
- Gene Ontology annotations for all genes
- PANTHER Slim simplified GO subsets
- Overrepresentation analysis tools

**Compound Integration:**
- Limited - primarily gene/protein focused
- Metabolite placeholders in pathway diagrams
- Use alongside compound-focused databases like SMPDB

---

### 8. Pathway Commons

| Field | Value |
|-------|-------|
| **URL** | https://www.pathwaycommons.org |
| **Content** | Integrated biological pathway and interaction data from multiple sources |
| **Records** | 5,772 pathways, 2.4 million interactions from 22 data sources |
| **Data Sources** | Reactome, KEGG, WikiPathways, HumanCyc, NCI-PID, PANTHER, PhosphoSitePlus, HPRD, BioGRID |
| **License** | Varies by source - check individual data source licenses |
| **API** | REST API |
| **API Documentation** | https://www.pathwaycommons.org/pc2/home |
| **Download** | https://www.pathwaycommons.org/archives/PC2/ |
| **Data Formats** | BioPAX Level 3, SIF (Simple Interaction Format), SBGN-ML |
| **Update Frequency** | Version 12 (March 2020); periodic updates |
| **Priority** | Tier 2 (Integration layer) |
| **Storage Estimate** | ~2 GB (complete BioPAX + SIF) |

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

**Key Features:**
- Unified access to multiple pathway databases
- Network queries: neighborhood, paths-between, common-stream
- Paxtools Java library for BioPAX manipulation
- CyPath2 Cytoscape integration

**SIF Format Example:**
```
TP53    controls-state-change-of    MDM2
BRCA1   interacts-with              BRCA2
EGFR    controls-phosphorylation-of MAPK1
```

---

### 9. NDEx (Network Data Exchange)

| Field | Value |
|-------|-------|
| **URL** | https://www.ndexbio.org |
| **Content** | Repository for biological network models of any type |
| **Records** | 15,000+ public networks; millions of nodes and edges |
| **Network Types** | Pathways, PPI, regulatory, signaling, metabolic, disease, drug-target |
| **License** | Per-network (many CC0 or CC BY) |
| **API** | REST API, Python client (ndex2) |
| **API Documentation** | https://home.ndexbio.org/using-the-ndex-server-api/ |
| **Python Client** | `pip install ndex2` |
| **Data Formats** | CX (Cytoscape Exchange JSON), CX2 (next generation) |
| **Update Frequency** | Continuous - user uploads and curated content |
| **Priority** | Tier 2 |
| **Storage Estimate** | Variable (selective download recommended) |

**Notable Network Collections:**
| Collection | Description | Networks |
|------------|-------------|----------|
| NCI-PID | Cancer pathways (archive) | 1,273 |
| SIGNOR | Signaling interactions | 200+ |
| STRING | PPI networks | Multiple |
| Drug-target | Drug mechanism networks | 500+ |
| Disease | Disease-specific networks | 1,000+ |

**Key Features:**
- DOI assignment for citable networks
- Version control for network revisions
- IQuery pathway enrichment service
- HiView hierarchical visualization
- CyNDEx Cytoscape integration

**Python API Example:**
```python
import ndex2

client = ndex2.client.Ndex2()
networks = client.search_networks('cancer signaling')
network = client.get_network_as_cx_stream(uuid)
```

---

### 10. PathBank

| Field | Value |
|-------|-------|
| **URL** | https://pathbank.org |
| **Content** | Human and microbial pathway database |
| **Records** | 110,000+ pathway maps covering primary metabolism, secondary metabolism, drug metabolism, signaling |
| **Species** | Human + 380+ microbial species |
| **License** | CC BY-NC 4.0 |
| **API** | Bulk download |
| **Download** | https://pathbank.org/downloads |
| **Data Formats** | CSV, SVG/PNG images, SBML |
| **Update Frequency** | Periodic |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~500 MB |

**Key Features:**
- Integrates with HMDB, DrugBank, SMPDB
- Detailed pathway diagrams with subcellular localization
- Microbial pathways for microbiome research
- Metabolite concentration data

---

## Comparative Analysis

### Coverage Comparison

| Database | Pathways | Reactions | Species | Compounds |
|----------|----------|-----------|---------|-----------|
| Reactome | 2,712 | 13,872 | 24 | 1,925 |
| WikiPathways | 3,100+ | N/A | 33 | Varies |
| KEGG | 559 ref | 12,049 | 8,000+ | 19,430 |
| MetaCyc | 2,937 | 18,124 | 3,295 | 18,400 |
| SMPDB | 99,900 | N/A | Human | Linked |
| PharmGKB | 200+ | N/A | Human | Drug focus |
| PANTHER | 177 | N/A | 142 | Limited |
| Pathway Commons | 5,772 | 2.4M int. | Human | Varies |
| NDEx | 15,000+ | Varies | Varies | Varies |
| PathBank | 110,000+ | N/A | 381 | Linked |

### Licensing Summary

| License Type | Databases | Commercial Use |
|--------------|-----------|----------------|
| **CC0 (Public Domain)** | WikiPathways | Yes |
| **CC BY 4.0** | Reactome, PANTHER | Yes (with attribution) |
| **CC BY-SA 4.0** | PharmGKB | Yes (with ShareAlike) |
| **CC BY-NC 4.0** | SMPDB, PathBank | No (academic only) |
| **Tiered/Subscription** | KEGG (bulk), MetaCyc/BioCyc | Varies |
| **Varies by Source** | Pathway Commons, NDEx | Check individual |

### API Availability

| Database | REST API | SPARQL | Bulk DL | Rate Limit |
|----------|----------|--------|---------|------------|
| Reactome | Yes | Yes | Yes | 100/min |
| WikiPathways | Yes | Yes | Yes | Reasonable |
| KEGG | Limited | No | Paid | 10/sec |
| MetaCyc | Yes | No | Paid | N/A |
| SMPDB | No | No | Yes | N/A |
| PharmGKB | Limited | No | Yes | N/A |
| PANTHER | No | No | Yes | N/A |
| Pathway Commons | Yes | No | Yes | Open |
| NDEx | Yes | No | Yes | Account-based |

### Data Format Support

| Format | Description | Databases |
|--------|-------------|-----------|
| **BioPAX** | Standard pathway exchange | Reactome, WikiPathways, Pathway Commons, MetaCyc |
| **SBML** | Systems biology models | Reactome, PANTHER, MetaCyc |
| **GPML** | Graphical pathway markup | WikiPathways, PharmGKB |
| **KGML** | KEGG-specific XML | KEGG |
| **CX/CX2** | Cytoscape exchange JSON | NDEx |
| **SIF** | Simple interaction format | Pathway Commons |
| **SBGN-ML** | Graphical notation | Reactome |

---

## Integration Recommendations

### Priority 1: Core (MVP)

| Database | Use Case | Integration Notes |
|----------|----------|-------------------|
| **Reactome** | Primary pathway source | Neo4j for network queries; BioPAX for integration |
| **WikiPathways** | Community pathways | Monthly Zenodo releases; GPML parsing |
| **SMPDB** | Drug metabolism | Essential for PK modeling; DrugBank/HMDB links |
| **PharmGKB Pathways** | Pharmacogenomics | Variant-drug-response for personalized medicine |

### Priority 2: Specialized

| Database | Use Case | Integration Notes |
|----------|----------|-------------------|
| **MetaCyc** (Tier 1) | Metabolic reference | Highest quality metabolic curation |
| **Pathway Commons** | Integration layer | Pre-integrated; BioPAX standard |
| **NDEx** | Specialized networks | NCI-PID archive; community networks |
| **PANTHER** | Signaling pathways | GO integration for enrichment |

### Priority 3: Reference

| Database | Use Case | Integration Notes |
|----------|----------|-------------------|
| **KEGG** | Reference queries | API for specific lookups; no bulk |
| **PathBank** | Microbiome pathways | Future integration for gut-brain axis |

### Data Harmonization Strategy

**Primary Identifiers:**
- Genes: Ensembl Gene ID, HGNC Symbol
- Proteins: UniProt ID
- Compounds: ChEBI ID (preferred), HMDB ID, PubChem CID

**Pathway Mapping:**
- Create cross-reference table: Reactome <-> WikiPathways <-> KEGG
- Use Pathway Commons for pre-built mappings

**Integration Flow:**
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
| PharmGKB | Pathways | ~100 MB | TSV, GPML |
| Pathway Commons | Complete | ~2 GB | BioPAX, SIF |

### Schema Mapping

```sql
-- Unified pathway-compound table
CREATE TABLE pathway_compound (
    pathway_id VARCHAR(50),
    pathway_source ENUM('reactome','wikipathways','smpdb','kegg','metacyc'),
    compound_id VARCHAR(50),
    compound_id_type ENUM('chebi','hmdb','pubchem','kegg','drugbank'),
    role ENUM('substrate','product','catalyst','inhibitor','modifier'),
    evidence_code VARCHAR(20),
    pubmed_id INTEGER
);
```

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [index.md](./../index.md) | Parent index |
| [disease.md](./disease.md) | Disease pathway databases |
| [pharmaceuticals.md](./../compounds/pharmaceuticals.md) | Drug database integration |

---

## References

1. Jassal B, et al. (2020). The reactome pathway knowledgebase. Nucleic Acids Res. 48(D1):D498-D503. https://doi.org/10.1093/nar/gkz1031

2. Martens M, et al. (2021). WikiPathways: connecting communities. Nucleic Acids Res. 49(D1):D613-D621. https://doi.org/10.1093/nar/gkaa1024

3. Kanehisa M, et al. (2023). KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Res. 51(D1):D587-D592. https://doi.org/10.1093/nar/gkac963

4. Caspi R, et al. (2020). The MetaCyc database of metabolic pathways and enzymes. Nucleic Acids Res. 48(D1):D445-D453. https://doi.org/10.1093/nar/gkz862

5. Jewison T, et al. (2014). SMPDB 2.0: Big Improvements to the Small Molecule Pathway Database. Nucleic Acids Res. 42:D478-D484. https://doi.org/10.1093/nar/gkt1067

6. Whirl-Carrillo M, et al. (2021). An Evidence-Based Framework for Evaluating Pharmacogenomics Knowledge. Clin Pharmacol Ther. 110(3):563-572. https://doi.org/10.1002/cpt.2350

7. Mi H, et al. (2021). PANTHER version 16. Nucleic Acids Res. 49(D1):D419-D426. https://doi.org/10.1093/nar/gkaa1106

8. Rodchenkov I, et al. (2020). Pathway Commons 2019 Update. Nucleic Acids Res. 48(D1):D489-D497. https://doi.org/10.1093/nar/gkz946

9. Pratt D, et al. (2015). NDEx, the Network Data Exchange. Cell Syst. 1(4):302-305. https://doi.org/10.1016/j.cels.2015.10.001

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial catalog: 10 pathway databases migrated from research.old/pathways-databases.md |
