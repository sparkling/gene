---
title: "Compound Databases"
parent: ../_index.md
category: shared
last_updated: 2026-01-22
status: draft
---

# Compound Databases

Chemical compound databases for natural products, drugs, and bioactive molecules.

## Database Catalog

| Database | Type | Tier | Coverage | Access Method | Size |
|----------|------|------|----------|---------------|------|
| **COCONUT** | Natural Products | 1 | 400,000+ natural products | REST API, Download | Open access |
| **LOTUS** | Natural Products | 1 | 750,000+ compound-organism pairs | SPARQL, Download | Open access |
| **DrugBank** | Drugs | 2 | 15,000+ drugs | REST API (subscription) | Comprehensive |
| **ChEMBL** | Bioactivity | 2 | 2.3M compounds, 20M activities | REST API, SQL | Open access |
| **PubChem** | Chemical Repository | 2 | 110M+ compounds | REST API, FTP | Massive |
| **ChEBI** | Chemical Ontology | 2 | 60,000+ entities | REST API, FTP | Ontology |
| **ZINC** | Purchasable Compounds | 3 | 1.4B compounds | Web, Download | Virtual screening |
| **NuBBE DB** | Brazilian Natural Products | 3 | 2,500+ compounds | Web | Regional |
| **StreptomeDB** | Streptomyces Products | 3 | 6,500+ compounds | Web | Specialized |
| **Super Natural II** | Natural Products | 3 | 325,000+ compounds | Web | 3D structures |

## Primary Use Cases

### Tier 1 (MVP) Focus

#### COCONUT (COlleCtion of Open Natural ProdUcTs)
- **Purpose**: Comprehensive natural product database
- **Content**:
  - 400,000+ unique natural products
  - Chemical structures (SMILES, InChI)
  - Source organisms
  - Literature references
- **Use**:
  - Link herbs/foods to chemical structures
  - Natural product similarity search
  - Traditional medicine compound identification
- **Access**: https://coconut.naturalproducts.net/
  - REST API: `https://coconut.naturalproducts.net/api/`
  - Bulk download: JSON, SDF
  - Free, no authentication
- **Update**: Periodic (check quarterly)

#### LOTUS (naturaL prOducTs occUrrences database)
- **Purpose**: Natural product-organism associations
- **Content**:
  - 750,000+ structure-organism pairs
  - 300,000+ unique structures
  - 30,000+ organisms
  - Literature-validated
- **Use**:
  - Map organisms to compounds
  - Source verification for traditional herbs
  - Biosynthetic pathway inference
- **Access**: https://lotus.naturalproducts.net/
  - SPARQL endpoint: `https://lotus.naturalproducts.net/sparql`
  - Bulk download: RDF, TSV
  - Free, open access
- **Update**: Continuous (community contributions)

## Tier 2 Databases

### DrugBank
- **Purpose**: Comprehensive drug database
- **Content**:
  - 15,000+ drugs (approved, investigational, experimental)
  - Drug-target interactions
  - Pharmacokinetics (ADME)
  - Drug-drug interactions
- **Use**: Drug-gene interactions, repurposing, safety
- **Access**: https://go.drugbank.com/
  - REST API (subscription required)
  - XML download (free for academic)
- **Note**: Free academic license available

### ChEMBL
- **Purpose**: Bioactivity database
- **Content**:
  - 2.3 million compounds
  - 20 million bioactivity measurements
  - 1.9 million assays
  - Target information
- **Use**: Compound-target predictions, bioactivity analysis
- **Access**: https://www.ebi.ac.uk/chembl/
  - REST API: `https://www.ebi.ac.uk/chembl/api/data/`
  - PostgreSQL dump
  - Free, open access

### PubChem
- **Purpose**: Public chemical repository
- **Content**:
  - 110+ million compounds
  - Chemical properties
  - Bioassay data
  - Literature links
- **Use**: Chemical structure lookup, property prediction
- **Access**: https://pubchem.ncbi.nlm.nih.gov/
  - REST API: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/`
  - FTP download
  - Free, no authentication

### ChEBI (Chemical Entities of Biological Interest)
- **Purpose**: Chemical ontology
- **Content**:
  - 60,000+ chemical entities
  - Hierarchical classification
  - Biological roles
  - Relationships
- **Use**: Chemical classification, ontology mapping
- **Access**: https://www.ebi.ac.uk/chebi/
  - REST API
  - OBO format (ontology)
  - Free, open access

## Data Integration Workflow

```
Traditional Herb / Food
         ↓
COCONUT / LOTUS (natural product structures)
         ↓
    Compounds
         ↓
┌────────┴────────┬────────┐
│                 │        │
ChEMBL          DrugBank   PubChem
│                 │        │
Bioactivity    Drug Info  Properties
│                 │        │
└────────┬────────┴────────┘
         ↓
   Target Genes
         ↓
    Pathways
```

## Access Methods

### COCONUT
```bash
# REST API - Search by name
curl "https://coconut.naturalproducts.net/api/compounds/search?text=curcumin"

# Get compound details
curl "https://coconut.naturalproducts.net/api/compounds/CNP0123456"

# Bulk download
wget https://coconut.naturalproducts.net/download/COCONUT_DB.sdf
wget https://coconut.naturalproducts.net/download/COCONUT_DB.json
```

### LOTUS
```bash
# SPARQL query - Find compounds from organism
curl -X POST "https://lotus.naturalproducts.net/sparql" \
  -H "Accept: application/json" \
  --data-urlencode "query=
  PREFIX wd: <http://www.wikidata.org/entity/>
  SELECT ?compound ?organism WHERE {
    ?compound lotus:fromOrganism ?organism .
    FILTER(CONTAINS(?organism, 'Ginkgo'))
  }"

# Bulk download
wget https://lotus.naturalproducts.net/download/lotus.tsv.gz
```

### ChEMBL
```bash
# REST API - Search by name
curl "https://www.ebi.ac.uk/chembl/api/data/molecule?molecule_synonyms__icontains=curcumin&format=json"

# Get compound bioactivity
curl "https://www.ebi.ac.uk/chembl/api/data/activity?molecule_chembl_id=CHEMBL365467&format=json"

# Get target information
curl "https://www.ebi.ac.uk/chembl/api/data/target/CHEMBL1824&format=json"
```

### DrugBank
```bash
# REST API (requires subscription and API key)
API_KEY="your_api_key"
curl "https://api.drugbank.com/v1/drugs/DB00945" \
  -H "Authorization: Bearer $API_KEY"

# Academic XML download (free)
# Register at https://go.drugbank.com/releases/latest
wget --user=username --password=password https://go.drugbank.com/releases/latest/downloads/all-full-database
```

### PubChem
```bash
# REST API - Search by name
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/curcumin/JSON"

# Get compound properties
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/969516/property/MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"

# Search by structure similarity
curl -X POST "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/JSON" \
  --data "smiles=CCO&Threshold=95"
```

## Compound Structure Formats

### Common Formats
- **SMILES**: Simplified Molecular Input Line Entry System
  - Example: `O=C(C=Cc1ccc(O)c(OC)c1)CC(=O)C=Cc1ccc(O)c(OC)c1` (Curcumin)
- **InChI**: International Chemical Identifier
  - Example: `InChI=1S/C21H20O6/c1-26-20-11-14(7-9-18(20)24)3-5-16(22)13-17(23)6-4-15-8-10-19(25)21(12-15)27-2/h3-12,24-25H,13H2,1-2H3`
- **InChIKey**: Hashed InChI for faster lookup
  - Example: `VFLDPWHFBUODDF-UHFFFAOYSA-N`
- **SDF**: Structure-Data File (3D structures)

## Example Use Cases

### Use Case 1: Map TCM Herb to Compounds
```bash
# Step 1: Get herb compounds from COCONUT
HERB="Ginkgo biloba"
curl "https://coconut.naturalproducts.net/api/compounds/search?text=${HERB}" > ginkgo_compounds.json

# Step 2: Extract compound IDs and get details
cat ginkgo_compounds.json | jq -r '.results[].coconut_id' | while read id; do
  curl "https://coconut.naturalproducts.net/api/compounds/${id}" >> ginkgo_details.json
done

# Step 3: Extract SMILES and search ChEMBL for bioactivity
cat ginkgo_details.json | jq -r '.smiles' | while read smiles; do
  curl "https://www.ebi.ac.uk/chembl/api/data/molecule?molecule_structures__canonical_smiles=${smiles}&format=json"
done
```

### Use Case 2: Find Drug-Gene Interactions
```bash
# Search DrugBank for APOE-affecting drugs
curl "https://api.drugbank.com/v1/targets?gene_name=APOE" \
  -H "Authorization: Bearer $API_KEY"
```

### Use Case 3: Compound Similarity Search
```bash
# Find compounds similar to curcumin in PubChem
CURCUMIN_CID="969516"
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/${CURCUMIN_CID}/cids/JSON?Threshold=90"
```

## Data Schema Examples

### COCONUT Compound
```json
{
  "coconut_id": "CNP0123456",
  "name": "Curcumin",
  "smiles": "O=C(C=Cc1ccc(O)c(OC)c1)CC(=O)C=Cc1ccc(O)c(OC)c1",
  "inchi": "InChI=1S/C21H20O6/...",
  "molecular_formula": "C21H20O6",
  "molecular_weight": 368.38,
  "source_organism": ["Curcuma longa"],
  "references": ["DOI:10.1021/..."]
}
```

### ChEMBL Activity
```json
{
  "activity_id": "ACT123456",
  "molecule_chembl_id": "CHEMBL365467",
  "target_chembl_id": "CHEMBL1824",
  "target_name": "Cyclooxygenase-2",
  "standard_type": "IC50",
  "standard_value": "0.45",
  "standard_units": "uM",
  "pchembl_value": "6.35"
}
```

## Cross-Database Mapping

### Identifier Cross-References
```
COCONUT ID ←→ PubChem CID ←→ ChEMBL ID ←→ DrugBank ID
    ↓              ↓              ↓              ↓
Natural       Chemical       Bioactivity    Drug
Product       Properties     Data           Information
```

### Example Mapping Workflow
```bash
# Start with natural product name
NAME="Curcumin"

# 1. COCONUT → SMILES
SMILES=$(curl "https://coconut.naturalproducts.net/api/compounds/search?text=${NAME}" | jq -r '.results[0].smiles')

# 2. SMILES → PubChem CID
CID=$(curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${SMILES}/cids/JSON" | jq -r '.IdentifierList.CID[0]')

# 3. PubChem CID → ChEMBL ID
# (Use PubChem's cross-reference links)

# 4. ChEMBL → Bioactivity
curl "https://www.ebi.ac.uk/chembl/api/data/activity?molecule_chembl_id=CHEMBL365467&format=json"
```

## Integration with Gene Platform

### Genetics ← Compounds
```
Genetic Variant → Gene → Pathway
                           ↑
                      ChEMBL/DrugBank
                           ↑
                    Compound Target
```

### Traditional Medicine → Compounds
```
TCM Herb (HERB, BATMAN-TCM)
         ↓
    COCONUT/LOTUS
         ↓
Natural Product Compounds
         ↓
    ChEMBL (Bioactivity)
         ↓
    Target Genes
```

### Nutrition → Compounds
```
Food (FooDB, USDA)
         ↓
Food Compound (Phenol-Explorer, PhytoHub)
         ↓
    PubChem/ChEBI
         ↓
Chemical Properties
         ↓
    ChEMBL (Bioactivity)
```

## Update Strategy

### Automated Updates (Tier 1)
- **COCONUT**: Check quarterly for new releases
- **LOTUS**: Continuous updates, sync monthly

### Manual Updates (Tier 2)
- **ChEMBL**: Quarterly releases (check in March, June, September, December)
- **DrugBank**: Annual major releases + quarterly updates
- **PubChem**: Daily updates (incremental sync)

## Storage Estimates

| Database | Storage Required | Format |
|----------|------------------|--------|
| COCONUT | 1-2 GB | JSON, SDF |
| LOTUS | 500 MB - 1 GB | TSV, RDF |
| ChEMBL | 50-100 GB (SQL) | PostgreSQL, SDF, JSON |
| DrugBank (academic) | 1-2 GB | XML, JSON |
| PubChem (full) | 100+ TB | SDF, XML (use API instead) |

## Compound Property Prediction

### Lipinski's Rule of Five (Drug-likeness)
- Molecular weight < 500 Da
- LogP < 5
- H-bond donors < 5
- H-bond acceptors < 10

### ADME Properties
- **A**bsorption: Caco-2 permeability
- **D**istribution: Plasma protein binding
- **M**etabolism: CYP450 interactions
- **E**xcretion: Renal clearance

### Toxicity Prediction
- hERG inhibition (cardiac toxicity)
- Hepatotoxicity
- Mutagenicity (Ames test)

## Navigation

- **Parent**: [Database Sources](../_index.md)
- **Related**: [Traditional Medicine](../traditional/_index.md), [Nutrition](../nutrition/_index.md), [Pathways](../pathways/_index.md)
