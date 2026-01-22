---
title: "Pathway Databases"
parent: ../_index.md
category: shared
last_updated: 2026-01-22
status: draft
---

# Pathway Databases

Databases for biological pathways, disease mechanisms, and systems biology networks.

## Database Catalog

| Database | Type | Tier | Coverage | Access Method | Size |
|----------|------|------|----------|---------------|------|
| **Reactome** | Pathway Knowledge | 1 | 2,600+ pathways, 11,000+ proteins | REST API, GraphQL | Human-focused |
| **DisGeNET** | Gene-Disease Network | 1 | 1.1M associations, 30,000+ diseases | REST API, Cytoscape | Comprehensive |
| **KEGG Pathway** | Reference Pathways | 2 | 550+ pathways (human) | REST API, FTP | Multi-organism |
| **WikiPathways** | Community Curation | 2 | 3,000+ pathways | REST API, SPARQL | Open access |
| **PathwayCommons** | Aggregator | 2 | 4,800+ pathways, 9 sources | REST API, Download | Meta-database |
| **BioCyc** | Metabolic Pathways | 3 | 18,000+ pathways, 17,000+ organisms | Web, Subscription | Comprehensive |
| **PANTHER** | Protein Classification | 2 | Pathway diagrams, GO terms | Web, API | Phylogenetic |
| **MSigDB** | Gene Sets | 2 | 32,000+ gene sets | Download, API | Computational |
| **STRING** | Protein Interactions | 2 | 67M interactions, 14,000+ organisms | REST API, Download | Interaction-focused |
| **SIGNOR** | Signaling Network | 2 | 26,000+ relations, 8,000+ proteins | Web, Download | Signaling-focused |

## Primary Use Cases

### Tier 1 (MVP) Focus

#### Reactome
- **Purpose**: Curated pathway knowledge base
- **Content**:
  - 2,600+ pathways
  - 11,000+ proteins
  - Reactions, interactions, regulatory events
  - Cross-species orthology
- **Use**:
  - Gene → Pathway mapping
  - Pathway enrichment analysis
  - Disease mechanism visualization
- **Access**: https://reactome.org/
  - REST API: `https://reactome.org/ContentService/`
  - GraphQL: `https://reactome.org/ContentService/graph-query`
  - Bulk download: BioPAX, SBML, PSI-MI
- **Update**: Quarterly releases

#### DisGeNET
- **Purpose**: Gene-disease association network
- **Content**:
  - 1.1 million gene-disease associations
  - 30,000+ diseases
  - 21,000+ genes
  - Curated and text-mined data
- **Use**:
  - Link genes to diseases
  - Disease similarity analysis
  - Comorbidity networks
- **Access**: https://www.disgenet.org/
  - REST API: `https://www.disgenet.org/api/`
  - Cytoscape plugin
  - RDF SPARQL endpoint
- **Update**: Annual major releases, quarterly updates

## Tier 2 Databases

### KEGG Pathway
- **Purpose**: Reference metabolic and regulatory pathways
- **Content**:
  - 550+ human pathways
  - Drug pathways
  - Disease pathways
  - Metabolic maps
- **Use**: Metabolic pathway analysis
- **Access**: https://www.kegg.jp/
  - REST API: `https://rest.kegg.jp/`
  - Subscription for bulk access
- **Note**: Free for academic use with rate limits

### WikiPathways
- **Purpose**: Community-curated pathway database
- **Content**:
  - 3,000+ pathways
  - 30+ species
  - Open access, editable
- **Use**: Pathway visualization, enrichment analysis
- **Access**: https://www.wikipathways.org/
  - REST API: `https://webservice.wikipathways.org/`
  - SPARQL: `https://sparql.wikipathways.org/`
  - GPML format (pathway markup)

### PathwayCommons
- **Purpose**: Aggregates multiple pathway databases
- **Content**:
  - 4,800+ pathways
  - Integrates: Reactome, KEGG, HumanCyc, PID, etc.
  - 9 source databases
- **Use**: Comprehensive pathway queries
- **Access**: https://www.pathwaycommons.org/
  - REST API: `https://www.pathwaycommons.org/pc2/`
  - BioPAX, SIF formats

### STRING
- **Purpose**: Protein-protein interaction networks
- **Content**:
  - 67 million interactions
  - 14,000+ organisms
  - Physical and functional associations
- **Use**: Network analysis, pathway discovery
- **Access**: https://string-db.org/
  - REST API: `https://string-db.org/api/`
  - Bulk download: TSV

## Data Integration Workflow

```
Gene List
    ↓
┌───────────┬───────────┬───────────┐
│           │           │           │
Reactome    KEGG        WikiPathways
│           │           │
Pathways    Metabolic   Community
│           Maps        Pathways
│           │           │
└───────────┴───────────┴───────────┘
              ↓
      Enriched Pathways
              ↓
         DisGeNET
              ↓
    Disease Associations
```

## Access Methods

### Reactome
```bash
# REST API - Get pathway by ID
curl "https://reactome.org/ContentService/data/query/R-HSA-73857"

# Get pathways for gene
curl "https://reactome.org/ContentService/data/pathways/low/entity/APOE"

# Pathway enrichment analysis
curl -X POST "https://reactome.org/AnalysisService/identifiers/?pageSize=20&page=1" \
     -H "Content-Type: text/plain" \
     -d "APOE BRCA1 TP53"

# Bulk download
wget https://reactome.org/download/current/ReactomePathways.txt
```

### DisGeNET
```bash
# REST API (requires API key for high volume)
curl "https://www.disgenet.org/api/gda/gene/APOE"

# Download full database
wget https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz

# Cytoscape integration
# Install DisGeNET Cytoscape app
```

### KEGG
```bash
# REST API - Get pathway
curl "https://rest.kegg.jp/get/hsa00010/kgml"

# List human pathways
curl "https://rest.kegg.jp/list/pathway/hsa"

# Get genes in pathway
curl "https://rest.kegg.jp/link/hsa/pathway:hsa00010"

# Note: Rate limiting applies, bulk access requires subscription
```

### WikiPathways
```bash
# REST API - List all pathways
curl "https://webservice.wikipathways.org/listPathways?organism=Homo+sapiens&format=json"

# Get pathway info
curl "https://webservice.wikipathways.org/getPathway?pwId=WP4844&format=json"

# SPARQL query
curl -X POST "https://sparql.wikipathways.org/sparql" \
     -H "Accept: application/json" \
     --data-urlencode "query=SELECT ?pathway WHERE { ?pathway a wp:Pathway }"
```

## Pathway Enrichment Analysis

### Standard Workflow
```
Input: Gene List (from genetics, traditional medicine, or nutrition analysis)
    ↓
Step 1: Map to Pathway Databases
    - Reactome: Biological processes
    - KEGG: Metabolic pathways
    - WikiPathways: Disease pathways
    ↓
Step 2: Statistical Enrichment
    - Hypergeometric test
    - Fisher's exact test
    - False discovery rate (FDR) correction
    ↓
Step 3: Disease Mapping
    - DisGeNET: Gene-disease associations
    - Disease similarity scoring
    ↓
Output: Enriched Pathways + Disease Associations
```

### Example: Gene Set to Diseases
```
Genes: [APOE, BRCA1, TP53, MTHFR, CYP2D6]
    ↓
Reactome Pathways:
    - DNA Repair (BRCA1, TP53)
    - Cholesterol metabolism (APOE)
    - One-carbon metabolism (MTHFR)
    - Drug metabolism (CYP2D6)
    ↓
DisGeNET Diseases:
    - Alzheimer's disease (APOE)
    - Breast cancer (BRCA1, TP53)
    - Cardiovascular disease (APOE, MTHFR)
    - Drug response variation (CYP2D6)
```

## Data Schema Examples

### Reactome Pathway
```json
{
  "stId": "R-HSA-73857",
  "displayName": "Lipoprotein metabolism",
  "species": "Homo sapiens",
  "genes": ["APOE", "LDLR", "APOB"],
  "reactions": 42,
  "hasDiagram": true
}
```

### DisGeNET Association
```json
{
  "geneId": "348",
  "geneSymbol": "APOE",
  "diseaseId": "C0002395",
  "diseaseName": "Alzheimer Disease",
  "score": 0.7,
  "source": ["CURATED", "LITERATURE"]
}
```

## Integration with Gene Platform

### Genetics → Pathways
```
User Variants (VCF)
    ↓
Affected Genes
    ↓
Reactome/KEGG Pathways
    ↓
DisGeNET Diseases
```

### Traditional Medicine → Pathways
```
TCM Herbs (HERB, BATMAN-TCM)
    ↓
Compound Targets
    ↓
Reactome/KEGG Pathways
    ↓
Disease Mechanisms
```

### Nutrition → Pathways
```
Foods (FooDB, USDA)
    ↓
Compounds (Polyphenols, etc.)
    ↓
Gene Targets
    ↓
Reactome/KEGG Pathways
    ↓
Health Benefits
```

## Update Strategy

### Automated Updates (Tier 1)
- **Reactome**: Quarterly releases (check in March, June, September, December)
- **DisGeNET**: Annual major release + quarterly updates

### Manual Checks (Tier 2)
- **KEGG**: Monthly updates (check for new pathways)
- **WikiPathways**: Continuous updates (check monthly)
- **PathwayCommons**: Biannual updates

## Storage Estimates

| Database | Storage Required | Format |
|----------|------------------|--------|
| Reactome | 500 MB - 1 GB | BioPAX, SBML, JSON |
| DisGeNET | 1-2 GB | TSV, RDF |
| KEGG | 100-500 MB | KGML, JSON |
| WikiPathways | 500 MB - 1 GB | GPML, JSON |
| PathwayCommons | 5-10 GB | BioPAX, SIF |

## Visualization Tools

| Tool | Purpose | Integration |
|------|---------|-------------|
| Cytoscape | Network visualization | DisGeNET, Reactome, STRING |
| Pathway Browser | Reactome viewer | Built-in to Reactome |
| KEGG Mapper | Pathway mapping | KEGG API |
| PathVisio | WikiPathways editor | GPML format |
| Graphviz | Graph rendering | DOT format export |

## Cross-Category Pathway Analysis

### Example: Alzheimer's Disease
```
Genetics:
- APOE ε4 variant → Cholesterol metabolism pathway → Increased risk

Traditional Medicine:
- Ginkgo biloba → Flavonoids → Neuroprotection pathway

Nutrition:
- Omega-3 fatty acids → Anti-inflammatory pathway → Brain health

→ Integrated pathway analysis:
  - Cholesterol metabolism
  - Inflammation
  - Oxidative stress
  - Neuroprotection
```

## Navigation

- **Parent**: [Database Sources](../_index.md)
- **Related**: [Genetics](../genetics/_index.md), [Traditional Medicine](../traditional/_index.md), [Nutrition](../nutrition/_index.md)
