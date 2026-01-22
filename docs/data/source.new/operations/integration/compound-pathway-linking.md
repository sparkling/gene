---
id: integration-compound-pathway-linking
title: "Comprehensive Guide: Linking Compounds to Genes and Biochemical Pathways"
type: integration
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [integration, cross-references, apis]
---

**Parent:** [../_index.md](../_index.md)


This guide provides a complete reference for integrating compound data (drugs, natural products, metabolites) with gene and pathway information. It covers identifier mapping strategies, database connections, API usage, and practical code examples.

---

## Table of Contents

1. [Identifier Mapping Hub](#1-identifier-mapping-hub)
2. [Compound to Target Protein Links](#2-compound-to-target-protein-links)
3. [Target Protein to Gene Links](#3-target-protein-to-gene-links)
4. [Gene/Protein to Pathway Links](#4-geneprotein-to-pathway-links)
5. [Compound to Pathway Direct Links](#5-compound-to-pathway-direct-links)
6. [Integration APIs](#6-integration-apis)
7. [Example Pipelines](#7-example-pipelines)
8. [SPARQL Federated Queries](#8-sparql-federated-queries)
9. [Bulk Data Integration](#9-bulk-data-integration)
10. [Code Examples](#10-code-examples)

---

## 1. Identifier Mapping Hub

Successful data integration requires standardized identifiers as "hubs" that connect disparate databases.

### 1.1 PubChem CID as Compound Hub

**PubChem Compound ID (CID)** is the most comprehensive compound identifier, linking to nearly all other chemical databases.

| Database | PubChem Mapping | Coverage |
|----------|-----------------|----------|
| ChEMBL | Direct via PubChem Substance | ~2M compounds |
| DrugBank | Via DrugBank external links | ~15K drugs |
| KEGG COMPOUND | Via PubChem cross-refs | ~18K compounds |
| ChEBI | Via PubChem cross-refs | ~150K entities |
| HMDB | Via PubChem cross-refs | ~220K metabolites |

**Key resources:**
- PubChem PUG REST: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/`
- PubChem Identifier Exchange: `https://pubchem.ncbi.nlm.nih.gov/idexchange/`

**Alternative hubs:**
- **InChIKey**: Structure-based, database-agnostic (use for exact structure matching)
- **ChEBI ID**: Preferred for ontological relationships
- **SMILES**: For structure-based queries (canonical SMILES recommended)

### 1.2 UniProt as Protein Hub

**UniProt Accession** is the definitive protein identifier, with comprehensive cross-references.

| Database | UniProt Mapping | Notes |
|----------|-----------------|-------|
| Ensembl | Direct cross-reference | ENSP (protein), ENSG (gene) |
| Entrez Gene | Direct cross-reference | GeneID |
| RefSeq | Direct cross-reference | NP_, XP_ accessions |
| PDB | Direct cross-reference | Structure data |
| Reactome | Direct cross-reference | Pathway participation |

**Key resource:**
- UniProt ID Mapping: `https://www.uniprot.org/id-mapping/`

### 1.3 Ensembl/Entrez as Gene Hub

| Hub | Best For | Cross-references |
|-----|----------|------------------|
| **Ensembl Gene ID** (ENSG) | Genomic context, variants | UniProt, RefSeq, Entrez |
| **Entrez Gene ID** | NCBI ecosystem, PubMed links | UniProt, Ensembl, KEGG |
| **HGNC ID** | Official human gene nomenclature | All major databases |

**Recommendation**: Use **HGNC ID** as the authoritative human gene identifier, with Ensembl for genomic data and Entrez for NCBI resources.

### 1.4 Reactome/GO as Pathway Hub

| Resource | ID Format | Best For |
|----------|-----------|----------|
| **Reactome** | R-HSA-XXXXXX | Curated biochemical pathways |
| **GO (Gene Ontology)** | GO:XXXXXXX | Biological process terms |
| **KEGG Pathway** | hsa:XXXXX | Metabolic and signaling |
| **WikiPathways** | WP:XXXX | Community-curated |

**Mapping strategy:**
```
Reactome <-> GO (via pathway2GO mappings)
KEGG <-> Reactome (partial, via gene membership)
WikiPathways <-> Reactome (via GPML annotations)
```

---

## 2. Compound to Target Protein Links

### 2.1 DrugBank: Drug to Target

**Coverage**: ~15,000 drug entries with ~5,000 protein targets

**Data structure:**
```
Drug (DrugBank ID) → Target (UniProt AC) + Action (inhibitor/agonist/etc.)
```

**Access methods:**

1. **XML Download** (requires academic license):
   ```
   https://go.drugbank.com/releases/latest
   ```

2. **DrugBank API** (subscription required):
   ```
   GET https://api.drugbank.com/v1/drugs/{drugbank_id}/targets
   ```

3. **Free alternatives**:
   - UniProt drug annotations
   - Open Targets Platform
   - DGIdb (aggregated)

**Key fields:**
- `target.uniprot_id`: UniProt accession
- `target.gene_name`: HGNC symbol
- `actions`: pharmacological action type
- `known_action`: yes/no/unknown

### 2.2 ChEMBL: Compound to Assay to Target

**Coverage**: ~2.4M compounds, ~15M activities, ~15K targets

**Data hierarchy:**
```
Compound (CHEMBL ID) → Assay (CHEMBL ID) → Target (CHEMBL ID → UniProt)
```

**REST API examples:**

```bash
# Get compound by ChEMBL ID
curl "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL25.json"

# Get activities for a compound
curl "https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=CHEMBL25&limit=100"

# Get target info
curl "https://www.ebi.ac.uk/chembl/api/data/target/CHEMBL2093872.json"
```

**Activity filtering criteria:**
- `pchembl_value >= 6` (IC50/Ki/Kd ≤ 1 μM)
- `standard_type` in ['IC50', 'Ki', 'Kd', 'EC50']
- `assay_type = 'B'` (binding assays)

**SQL query (ChEMBL database):**
```sql
SELECT DISTINCT
    md.chembl_id AS compound_id,
    cs.canonical_smiles,
    td.chembl_id AS target_id,
    cp.accession AS uniprot_id,
    act.pchembl_value
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
JOIN activities act ON md.molregno = act.molregno
JOIN assays ass ON act.assay_id = ass.assay_id
JOIN target_dictionary td ON ass.tid = td.tid
JOIN target_components tc ON td.tid = tc.tid
JOIN component_sequences cp ON tc.component_id = cp.component_id
WHERE act.pchembl_value >= 6
  AND ass.assay_type = 'B';
```

### 2.3 BindingDB: Compound to Binding to Protein

**Coverage**: ~2.8M binding measurements, ~1.2M compounds, ~9K targets

**Data structure:**
```
Compound (BindingDB ID/SMILES) → Binding Affinity → Protein (UniProt)
```

**Download**: `https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp`

**Key columns in TSV:**
- `Ligand SMILES`: Canonical SMILES
- `Ligand InChI Key`: For cross-referencing
- `Target Name Assigned by Curator`: Protein name
- `UniProt (SwissProt) Primary ID of Target Chain`: UniProt AC
- `Ki (nM)`, `IC50 (nM)`, `Kd (nM)`: Binding affinities

**Filtering:**
```python
# High-confidence interactions
df_filtered = df[
    (df['Ki (nM)'].notna() & (df['Ki (nM)'] < 1000)) |
    (df['IC50 (nM)'].notna() & (df['IC50 (nM)'] < 1000)) |
    (df['Kd (nM)'].notna() & (df['Kd (nM)'] < 1000))
]
```

### 2.4 STITCH: Compound to Protein Interaction

**Coverage**: ~500K chemicals, ~9.6M proteins (all species)

**Data structure:**
```
Chemical (STITCH CID) ↔ Protein (STRING ID → UniProt)
```

**Scoring system:**
- `combined_score`: 0-999 (use ≥700 for high confidence)
- Evidence channels: experimental, database, textmining, prediction

**Downloads** (Human):
```
http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz
http://stitch.embl.de/download/chemical.aliases.v5.0.tsv.gz
```

**ID conversion:**
- STITCH CID format: `CIDm` or `CIDs` + PubChem CID (m=metabolite, s=stereo)
- STRING protein ID: `9606.ENSPXXXXXXXXX` (taxon + Ensembl protein)

### 2.5 DGIdb: Aggregated Drug-Gene Interactions

**Coverage**: ~100K interactions from 40+ sources

**Sources aggregated:**
- DrugBank, ChEMBL, PharmGKB, TTD, TEND
- Cancer-specific: CIViC, OncoKB, CGI

**REST API:**
```bash
# Search by drug name
curl "https://dgidb.org/api/v2/interactions.json?drugs=imatinib"

# Search by gene
curl "https://dgidb.org/api/v2/interactions.json?genes=ABL1"
```

**Response fields:**
- `gene_name`: HGNC symbol
- `drug_name`: Drug name
- `interaction_types`: inhibitor, activator, etc.
- `sources`: Contributing databases

---

## 3. Target Protein to Gene Links

### 3.1 UniProt to Ensembl/Entrez Mapping

**UniProt ID Mapping Service:**

```python
import requests

def uniprot_to_gene(uniprot_ids: list) -> dict:
    """Map UniProt accessions to Ensembl and Entrez Gene IDs."""

    # Submit mapping job
    url = "https://rest.uniprot.org/idmapping/run"
    data = {
        "from": "UniProtKB_AC-ID",
        "to": "GeneID",  # or "Ensembl"
        "ids": ",".join(uniprot_ids)
    }
    response = requests.post(url, data=data)
    job_id = response.json()["jobId"]

    # Poll for results
    result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    while True:
        result = requests.get(result_url)
        if result.status_code == 200:
            return result.json()
        time.sleep(1)
```

**Mapping file download:**
```
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
- idmapping_selected.tab.gz (full mapping table)
- by_organism/HUMAN_9606_idmapping.dat.gz
```

**Key ID types in idmapping:**
| ID Type | Column | Example |
|---------|--------|---------|
| UniProtKB-AC | 1 | P00533 |
| GeneID (Entrez) | 3 | 1956 |
| Ensembl | via dat file | ENSG00000146648 |
| HGNC | via dat file | HGNC:3236 |

### 3.2 HGNC as Gene Symbol Authority

**HGNC (HUGO Gene Nomenclature Committee)** provides official human gene symbols.

**Download:**
```
https://www.genenames.org/download/archive/
- hgnc_complete_set.txt (main file)
- withdrawn.txt (historical symbols)
```

**Key columns:**
- `hgnc_id`: HGNC:XXXXX (stable identifier)
- `symbol`: Official gene symbol
- `name`: Full gene name
- `alias_symbol`: Previous/alternative symbols
- `ensembl_gene_id`: ENSG...
- `entrez_id`: Entrez Gene ID
- `uniprot_ids`: UniProt accession(s)

**API access:**
```bash
# Search by symbol
curl "https://rest.genenames.org/search/symbol/EGFR" -H "Accept: application/json"

# Fetch by HGNC ID
curl "https://rest.genenames.org/fetch/hgnc_id/HGNC:3236" -H "Accept: application/json"
```

---

## 4. Gene/Protein to Pathway Links

### 4.1 Reactome: UniProt to Pathway

**Coverage**: ~2,700 human pathways, ~11K proteins

**REST API:**
```bash
# Get pathways for a protein
curl "https://reactome.org/ContentService/data/pathways/low/diagram/entity/UniProt:P00533"

# Get pathway details
curl "https://reactome.org/ContentService/data/query/R-HSA-177929"

# Get participating proteins in pathway
curl "https://reactome.org/ContentService/data/participants/R-HSA-177929/participatingPhysicalEntities"
```

**Analysis Service (batch):**
```python
import requests

def get_pathways_for_genes(gene_list: list) -> dict:
    """Submit gene list for pathway enrichment."""
    url = "https://reactome.org/AnalysisService/identifiers/projection"
    headers = {"Content-Type": "text/plain"}
    data = "\n".join(gene_list)

    response = requests.post(url, headers=headers, data=data)
    return response.json()
```

**Download files:**
```
https://reactome.org/download/current/
- UniProt2Reactome_All_Levels.txt
- UniProt2Reactome_PE_Pathway.txt
- ReactomePathways.txt
- ReactomePathwaysRelation.txt
```

### 4.2 KEGG: Gene to Pathway

**Coverage**: ~350 human pathways

**API (REST-like):**
```bash
# Get pathways for a gene (KEGG gene ID = Entrez)
curl "https://rest.kegg.jp/link/pathway/hsa:1956"

# Get genes in a pathway
curl "https://rest.kegg.jp/link/hsa/hsa00010"

# Get pathway info
curl "https://rest.kegg.jp/get/hsa00010"
```

**Mapping KEGG gene IDs:**
- KEGG uses Entrez Gene IDs with organism prefix: `hsa:1956` = EGFR
- Convert UniProt → Entrez → KEGG gene

**Download:**
```
https://rest.kegg.jp/link/pathway/hsa (all gene-pathway links)
https://rest.kegg.jp/list/pathway/hsa (all human pathways)
```

### 4.3 WikiPathways: Gene to Pathway

**Coverage**: ~3,000 human pathways (community-curated)

**SPARQL endpoint:**
```
https://sparql.wikipathways.org/sparql
```

**REST API:**
```bash
# Search pathways by gene
curl "https://webservice.wikipathways.org/findPathwaysByXref?ids=1956&codes=L&format=json"

# Get pathway info
curl "https://webservice.wikipathways.org/getPathwayInfo?pwId=WP254&format=json"
```

**ID codes for xref search:**
- `L` = Entrez Gene
- `En` = Ensembl
- `S` = UniProt

**Download:**
```
https://data.wikipathways.org/current/gpml/
- wikipathways-YYYYMMDD-gpml-Homo_sapiens.zip
```

### 4.4 GO: Gene to Biological Process

**Coverage**: ~20K GO terms, ~20K human genes with annotations

**QuickGO API:**
```bash
# Get GO terms for a gene (UniProt)
curl "https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=P00533&aspect=biological_process"

# Get genes for a GO term
curl "https://www.ebi.ac.uk/QuickGO/services/annotation/search?goId=GO:0006915&taxonId=9606"
```

**Download GAF files:**
```
http://current.geneontology.org/annotations/
- goa_human.gaf.gz
```

**GAF key columns:**
- Column 2: UniProt AC
- Column 5: GO ID
- Column 7: Evidence code
- Column 9: Aspect (P=Process, F=Function, C=Component)

**Evidence code filtering:**
- Experimental: EXP, IDA, IPI, IMP, IGI, IEP
- Computational: ISS, ISO, ISA, ISM, IGC, IBA, IBD, IKR, IRD, RCA
- Author statement: TAS, NAS
- Curator: IC, ND
- Electronic: IEA (lowest confidence)

---

## 5. Compound to Pathway Direct Links

### 5.1 SMPDB: Compound in Pathway

**Coverage**: ~48K metabolic/signaling pathways, ~5K metabolites

**Data structure:**
```
Pathway (SMPDB ID) contains Metabolite (HMDB ID/KEGG ID)
```

**Download:**
```
https://smpdb.ca/downloads
- smpdb_pathways.csv.zip
- smpdb_metabolites.csv.zip
- smpdb_proteins.csv.zip
```

**Metabolite-pathway mapping:**
```python
# smpdb_metabolites.csv columns
# SMPDB ID, Pathway Name, PathWhiz ID, KEGG ID, HMDB ID, PubChem CID, Name
```

### 5.2 Reactome: Small Molecule Participants

**Download:**
```
https://reactome.org/download/current/
- ChEBI2Reactome_All_Levels.txt
```

**Format:** `ChEBI_ID \t Reactome_Pathway_ID \t URL \t Pathway_Name \t Evidence \t Species`

**API query:**
```bash
# Get pathways containing a small molecule (ChEBI ID)
curl "https://reactome.org/ContentService/data/pathways/low/diagram/entity/ChEBI:15377"
```

**Note:** Reactome uses ChEBI IDs for small molecules. Map from PubChem CID → ChEBI first.

### 5.3 KEGG COMPOUND in KEGG PATHWAY

**API:**
```bash
# Get compounds in a pathway
curl "https://rest.kegg.jp/link/cpd/hsa00010"

# Get pathways containing a compound
curl "https://rest.kegg.jp/link/pathway/cpd:C00031"

# Get compound info
curl "https://rest.kegg.jp/get/cpd:C00031"
```

**Bulk download:**
```
https://rest.kegg.jp/link/pathway/compound (all compound-pathway links)
```

**ID mapping:**
- Map PubChem CID → KEGG COMPOUND: Use PubChem xrefs or KEGG `conv` API
```bash
curl "https://rest.kegg.jp/conv/compound/pubchem"
```

### 5.4 MetaCyc: Compound in Reaction

**Coverage**: ~18K pathways, ~17K reactions, ~16K compounds (all organisms)

**Data structure:**
```
Pathway → Reaction → Compound (substrate/product)
```

**Download** (requires license):
```
https://metacyc.org/downloads.shtml
- compounds.dat
- reactions.dat
- pathways.dat
```

**BioCyc API:**
```bash
# Get compound info
curl "https://websvc.biocyc.org/getxml?id=META:ATP"

# Get reactions involving compound
curl "https://websvc.biocyc.org/META/compound-rxns?id=ATP"
```

---

## 6. Integration APIs

### 6.1 Reactome Content Service

**Base URL:** `https://reactome.org/ContentService`

**Key endpoints:**

```bash
# Pathway analysis (gene/protein list)
POST /AnalysisService/identifiers/projection
Content-Type: text/plain
Body: UniProt IDs (one per line)

# Get participants in pathway
GET /data/participants/{pathway_id}/participatingPhysicalEntities

# Get pathways containing entity
GET /data/pathways/low/diagram/entity/{entity_id}

# Map IDs
GET /data/mapping/{resource}/{identifier}
# resource: UniProt, ChEBI, KEGG, etc.

# Cross-references
GET /data/xrefs/{identifier}
```

**Example: Complete pathway query**
```python
import requests

def get_full_pathway_info(pathway_id: str) -> dict:
    """Get comprehensive pathway information."""
    base = "https://reactome.org/ContentService"

    # Pathway details
    details = requests.get(f"{base}/data/query/{pathway_id}").json()

    # Participating proteins
    proteins = requests.get(
        f"{base}/data/participants/{pathway_id}/participatingPhysicalEntities"
    ).json()

    # Participating molecules
    molecules = requests.get(
        f"{base}/data/participants/{pathway_id}/referenceEntities"
    ).json()

    return {
        "details": details,
        "proteins": [p for p in proteins if "UniProt" in str(p)],
        "molecules": [m for m in molecules if "ChEBI" in str(m) or "COMPOUND" in str(m)]
    }
```

### 6.2 UniProt ID Mapping Service

**New REST API (2023+):**

```python
import requests
import time

class UniProtMapper:
    """UniProt ID mapping client."""

    BASE_URL = "https://rest.uniprot.org"

    VALID_FROM = [
        "UniProtKB_AC-ID", "UniParc", "UniRef50", "UniRef90", "UniRef100",
        "Gene_Name", "CRC64", "Ensembl", "Ensembl_Genomes", "GeneID",
        "KEGG", "PATRIC", "UCSC", "WBParaSite", "ChEMBL", "DrugBank",
        "GeneCards", "HGNC", "EMBL", "EMBL-CDS", "PIR", "RefSeq_Protein"
    ]

    VALID_TO = VALID_FROM + [
        "UniProtKB", "UniProtKB-Swiss-Prot"
    ]

    def map_ids(self, ids: list, from_db: str, to_db: str) -> dict:
        """Map IDs between databases."""

        # Submit job
        response = requests.post(
            f"{self.BASE_URL}/idmapping/run",
            data={
                "from": from_db,
                "to": to_db,
                "ids": ",".join(ids)
            }
        )
        job_id = response.json()["jobId"]

        # Poll for results
        while True:
            status = requests.get(
                f"{self.BASE_URL}/idmapping/status/{job_id}"
            ).json()

            if "results" in status or "failedIds" in status:
                break
            if status.get("jobStatus") == "FINISHED":
                break
            time.sleep(1)

        # Get results
        results = requests.get(
            f"{self.BASE_URL}/idmapping/results/{job_id}"
        ).json()

        return results

# Usage
mapper = UniProtMapper()
results = mapper.map_ids(
    ["P00533", "P04637", "P38398"],
    "UniProtKB_AC-ID",
    "GeneID"
)
```

### 6.3 PubChem PUG REST

**Base URL:** `https://pubchem.ncbi.nlm.nih.gov/rest/pug`

**Key endpoints:**

```bash
# Get compound by CID
GET /compound/cid/{cid}/JSON

# Get compound by name
GET /compound/name/{name}/cids/JSON

# Get compound by SMILES
GET /compound/smiles/{smiles}/cids/JSON

# Cross-references
GET /compound/cid/{cid}/xrefs/{xref_type}/JSON
# xref_type: PatentID, SourceName, RegistryID, etc.

# Get assay targets
GET /compound/cid/{cid}/assaysummary/JSON

# Get compound synonyms
GET /compound/cid/{cid}/synonyms/JSON
```

**Example: Get compound targets from bioassays**
```python
import requests

def get_pubchem_targets(cid: int) -> list:
    """Get protein targets from PubChem bioassays."""

    # Get assay summary
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    response = requests.get(url)

    if response.status_code != 200:
        return []

    data = response.json()
    targets = []

    for table in data.get("Table", {}).get("Row", []):
        cells = table.get("Cell", [])
        if len(cells) > 5:
            target_info = {
                "aid": cells[0],
                "target_gi": cells[4],
                "target_name": cells[5],
                "activity_outcome": cells[2]
            }
            if target_info["activity_outcome"] == "Active":
                targets.append(target_info)

    return targets
```

### 6.4 BridgeDb for ID Conversion

**Web service:** `https://webservice.bridgedb.org`

**Supported databases:**
- Ensembl (En), Entrez Gene (L), UniProt (S)
- ChEBI (Ce), KEGG Compound (Ck), PubChem (Cpc)
- HMDB (Ch), DrugBank (Dr)

**API:**
```bash
# Map single ID
GET /{organism}/xrefs/{systemCode}/{identifier}

# Map to specific database
GET /{organism}/xrefs/{systemCode}/{identifier}?dataSource={targetSystem}

# Batch mapping
POST /{organism}/xrefsBatch/{systemCode}/{targetSystemCode}
Content-Type: text/plain
Body: identifiers (one per line)
```

**Example:**
```bash
# Map UniProt to Ensembl (Human)
curl "https://webservice.bridgedb.org/Human/xrefs/S/P00533?dataSource=En"

# Map KEGG Compound to PubChem
curl "https://webservice.bridgedb.org/Human/xrefs/Ck/C00031?dataSource=Cpc"
```

**Python client:**
```python
import requests

class BridgeDbMapper:
    """BridgeDb ID mapping client."""

    BASE_URL = "https://webservice.bridgedb.org"

    SYSTEM_CODES = {
        "ensembl": "En",
        "entrez": "L",
        "uniprot": "S",
        "hgnc": "H",
        "chebi": "Ce",
        "kegg_compound": "Ck",
        "pubchem": "Cpc",
        "hmdb": "Ch",
        "drugbank": "Dr",
        "chembl": "Cl"
    }

    def map_id(self, identifier: str, from_system: str,
               to_system: str = None, organism: str = "Human") -> list:
        """Map single ID."""
        from_code = self.SYSTEM_CODES.get(from_system, from_system)

        url = f"{self.BASE_URL}/{organism}/xrefs/{from_code}/{identifier}"
        if to_system:
            to_code = self.SYSTEM_CODES.get(to_system, to_system)
            url += f"?dataSource={to_code}"

        response = requests.get(url)
        if response.status_code != 200:
            return []

        results = []
        for line in response.text.strip().split("\n"):
            if "\t" in line:
                mapped_id, system = line.split("\t")
                results.append({"id": mapped_id, "system": system})

        return results
```

---

## 7. Example Pipelines

### 7.1 Curcumin: Compound to Pathways

```
Curcumin (PubChem CID: 969516)
    ↓ PubChem → ChEMBL
ChEMBL ID: CHEMBL116438
    ↓ ChEMBL activities API
Targets: COX2, NF-kB, TNF-α, etc.
    ↓ ChEMBL target → UniProt
UniProt: P35354 (PTGS2), Q04206 (RELA), P01375 (TNF)
    ↓ UniProt → Reactome
Pathways:
  - R-HSA-6785807: Interleukin-4 and 13 signaling
  - R-HSA-9020702: Interleukin-1 signaling
  - R-HSA-448706: Interleukin-1 processing
```

**Code implementation:**
```python
import requests

def curcumin_to_pathways():
    """Trace curcumin from compound to pathways."""

    # Step 1: Get ChEMBL ID from PubChem
    pubchem_cid = 969516
    chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__canonical_smiles__flexmatch=COC1=C(C=CC(=C1)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC)O"

    # Or search by name
    chembl_search = requests.get(
        "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q=curcumin"
    ).json()
    chembl_id = "CHEMBL116438"

    # Step 2: Get targets from ChEMBL
    activities = requests.get(
        f"https://www.ebi.ac.uk/chembl/api/data/activity.json",
        params={
            "molecule_chembl_id": chembl_id,
            "pchembl_value__gte": 5,
            "limit": 100
        }
    ).json()

    target_chembl_ids = set()
    for act in activities.get("activities", []):
        if act.get("target_chembl_id"):
            target_chembl_ids.add(act["target_chembl_id"])

    # Step 3: Get UniProt IDs from targets
    uniprot_ids = []
    for target_id in list(target_chembl_ids)[:10]:
        target = requests.get(
            f"https://www.ebi.ac.uk/chembl/api/data/target/{target_id}.json"
        ).json()

        for component in target.get("target_components", []):
            if component.get("accession"):
                uniprot_ids.append(component["accession"])

    # Step 4: Get pathways from Reactome
    pathways = {}
    for uniprot_id in uniprot_ids[:5]:
        response = requests.get(
            f"https://reactome.org/ContentService/data/pathways/low/diagram/entity/UniProt:{uniprot_id}",
            headers={"Accept": "application/json"}
        )
        if response.status_code == 200:
            for pathway in response.json():
                pathways[pathway["stId"]] = pathway["displayName"]

    return {
        "compound": "Curcumin",
        "pubchem_cid": pubchem_cid,
        "chembl_id": chembl_id,
        "targets": uniprot_ids,
        "pathways": pathways
    }
```

### 7.2 CYP2D6: Gene to Drugs

```
CYP2D6 (HGNC:2625)
    ↓ HGNC → UniProt
UniProt: P10635
    ↓ UniProt → Reactome
Pathway: R-HSA-211981 (Xenobiotics)
    ↓ DrugBank/PharmGKB
Drugs metabolized: Codeine, Tamoxifen, Metoprolol, etc.
    ↓ Drug effects
Polymorphism impact on drug metabolism
```

**Code implementation:**
```python
import requests

def cyp2d6_drug_pathway():
    """Map CYP2D6 gene to drug metabolism pathway and affected drugs."""

    # Step 1: Get gene info from HGNC
    hgnc_response = requests.get(
        "https://rest.genenames.org/fetch/symbol/CYP2D6",
        headers={"Accept": "application/json"}
    ).json()

    gene_info = hgnc_response["response"]["docs"][0]
    uniprot_id = gene_info.get("uniprot_ids", ["P10635"])[0]
    entrez_id = gene_info.get("entrez_id", "1565")

    # Step 2: Get pathway from Reactome
    pathways = requests.get(
        f"https://reactome.org/ContentService/data/pathways/low/diagram/entity/UniProt:{uniprot_id}",
        headers={"Accept": "application/json"}
    ).json()

    # Step 3: Get drugs from DGIdb
    dgidb_response = requests.get(
        f"https://dgidb.org/api/v2/interactions.json?genes=CYP2D6"
    ).json()

    drugs = []
    for match in dgidb_response.get("matchedTerms", []):
        for interaction in match.get("interactions", []):
            drugs.append({
                "drug_name": interaction.get("drugName"),
                "interaction_types": interaction.get("interactionTypes"),
                "sources": interaction.get("sources")
            })

    return {
        "gene": "CYP2D6",
        "hgnc_id": "HGNC:2625",
        "uniprot_id": uniprot_id,
        "entrez_id": entrez_id,
        "pathways": [{"id": p["stId"], "name": p["displayName"]} for p in pathways],
        "drugs_affected": drugs[:20]
    }
```

---

## 8. SPARQL Federated Queries

### 8.1 Wikidata + WikiPathways

**Wikidata endpoint:** `https://query.wikidata.org/sparql`
**WikiPathways endpoint:** `https://sparql.wikipathways.org/sparql`

**Query: Find pathways for drugs targeting a specific protein**
```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX dcterms: <http://purl.org/dc/terms/>

# Find pathways containing genes targeted by imatinib
SELECT DISTINCT ?pathway ?pathwayLabel ?gene ?geneLabel
WHERE {
  # From Wikidata: imatinib targets
  SERVICE <https://query.wikidata.org/sparql> {
    wd:Q177094 wdt:P129 ?target .  # imatinib (Q177094) physically interacts with
    ?target wdt:P353 ?hgncSymbol .  # get HGNC gene symbol
    ?target rdfs:label ?geneLabel .
    FILTER(LANG(?geneLabel) = "en")
  }

  # From WikiPathways: find pathways with these genes
  ?gene a wp:GeneProduct ;
        rdfs:label ?geneLabel ;
        dcterms:isPartOf ?pathway .
  ?pathway a wp:Pathway ;
           dcterms:title ?pathwayLabel .

  FILTER(CONTAINS(LCASE(?geneLabel), LCASE(?hgncSymbol)))
}
LIMIT 100
```

**Query: Compounds in WikiPathways with Wikidata enrichment**
```sparql
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT ?metabolite ?metaboliteLabel ?pathway ?pathwayTitle ?wikidataItem ?chemblId
WHERE {
  # WikiPathways: metabolites in pathways
  ?metabolite a wp:Metabolite ;
              rdfs:label ?metaboliteLabel ;
              dcterms:isPartOf ?pathway ;
              wp:bdbWikidata ?wikidataItem .
  ?pathway a wp:Pathway ;
           dcterms:title ?pathwayTitle ;
           wp:organism <http://purl.obolibrary.org/obo/NCBITaxon_9606> .

  # Wikidata: get ChEMBL ID
  SERVICE <https://query.wikidata.org/sparql> {
    OPTIONAL { ?wikidataItem wdt:P592 ?chemblId . }
  }
}
LIMIT 100
```

### 8.2 Wikidata + UniProt (via Reactome)

**Query: Drug targets with pathway context**
```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>

SELECT ?drug ?drugLabel ?target ?targetLabel ?uniprotId ?reactomePathway
WHERE {
  # Wikidata: approved drugs with protein targets
  ?drug wdt:P31 wd:Q12140 ;  # instance of medication
        wdt:P129 ?target ;    # physically interacts with
        rdfs:label ?drugLabel .
  FILTER(LANG(?drugLabel) = "en")

  ?target wdt:P352 ?uniprotId ;  # UniProt ID
          rdfs:label ?targetLabel .
  FILTER(LANG(?targetLabel) = "en")

  # Get Reactome pathway
  OPTIONAL { ?target wdt:P3937 ?reactomePathway . }
}
LIMIT 100
```

**Query: Complete drug-gene-pathway chain**
```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT ?drug ?drugLabel ?pubchemCid ?target ?geneSymbol ?uniprotId ?ensemblGene ?reactomeId
WHERE {
  # Drug with identifiers
  ?drug wdt:P31 wd:Q12140 ;         # medication
        wdt:P662 ?pubchemCid ;       # PubChem CID
        wdt:P129 ?target ;           # target protein
        rdfs:label ?drugLabel .
  FILTER(LANG(?drugLabel) = "en")

  # Target protein with gene info
  ?target wdt:P353 ?geneSymbol ;     # HGNC symbol
          wdt:P352 ?uniprotId .      # UniProt ID

  # Gene identifiers
  OPTIONAL { ?target wdt:P594 ?ensemblGene . }  # Ensembl gene ID

  # Reactome pathway
  OPTIONAL { ?target wdt:P3937 ?reactomeId . }
}
ORDER BY ?drugLabel
LIMIT 500
```

---

## 9. Bulk Data Integration

### 9.1 Download All Mapping Files

**Download script:**
```bash
#!/bin/bash

# Create directory structure
mkdir -p data/{compounds,genes,proteins,pathways,mappings}

# === COMPOUNDS ===
# PubChem compound → external IDs
wget -P data/compounds/ \
  "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"

# ChEMBL compound mappings
wget -P data/compounds/ \
  "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_*_unichem.txt.gz"

# === GENES ===
# HGNC complete set
wget -P data/genes/ \
  "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"

# Ensembl ID mapping
wget -P data/genes/ \
  "https://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.*.entrez.tsv.gz"

# === PROTEINS ===
# UniProt ID mapping
wget -P data/proteins/ \
  "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"

# === PATHWAYS ===
# Reactome mappings
wget -P data/pathways/ \
  "https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"
wget -P data/pathways/ \
  "https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt"
wget -P data/pathways/ \
  "https://reactome.org/download/current/ReactomePathways.txt"

# GO annotations
wget -P data/pathways/ \
  "http://current.geneontology.org/annotations/goa_human.gaf.gz"

# === INTERACTIONS ===
# ChEMBL activities (requires processing)
wget -P data/mappings/ \
  "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_*_sqlite.tar.gz"

# STITCH (human only)
wget -P data/mappings/ \
  "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"

# DGIdb interactions
wget -P data/mappings/ \
  "https://dgidb.org/data/monthly_tsvs/interactions.tsv"

echo "Downloads complete. Extract and process files."
```

### 9.2 Build Local SQLite Mapping Database

**Schema:**
```sql
-- compounds table
CREATE TABLE compounds (
    pubchem_cid INTEGER PRIMARY KEY,
    chembl_id TEXT,
    drugbank_id TEXT,
    chebi_id TEXT,
    kegg_id TEXT,
    hmdb_id TEXT,
    inchikey TEXT,
    name TEXT
);
CREATE INDEX idx_compounds_chembl ON compounds(chembl_id);
CREATE INDEX idx_compounds_inchikey ON compounds(inchikey);

-- genes table
CREATE TABLE genes (
    hgnc_id TEXT PRIMARY KEY,
    symbol TEXT NOT NULL,
    name TEXT,
    entrez_id INTEGER,
    ensembl_id TEXT,
    uniprot_ids TEXT  -- comma-separated
);
CREATE INDEX idx_genes_symbol ON genes(symbol);
CREATE INDEX idx_genes_entrez ON genes(entrez_id);
CREATE INDEX idx_genes_ensembl ON genes(ensembl_id);

-- proteins table
CREATE TABLE proteins (
    uniprot_id TEXT PRIMARY KEY,
    gene_symbol TEXT,
    hgnc_id TEXT,
    entrez_id INTEGER,
    ensembl_gene TEXT,
    ensembl_protein TEXT,
    protein_name TEXT,
    FOREIGN KEY (hgnc_id) REFERENCES genes(hgnc_id)
);
CREATE INDEX idx_proteins_gene ON proteins(gene_symbol);

-- pathways table
CREATE TABLE pathways (
    pathway_id TEXT PRIMARY KEY,
    source TEXT,  -- reactome, kegg, go, wikipathways
    name TEXT,
    species TEXT
);
CREATE INDEX idx_pathways_source ON pathways(source);

-- compound_target interactions
CREATE TABLE compound_targets (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pubchem_cid INTEGER,
    chembl_id TEXT,
    uniprot_id TEXT,
    activity_type TEXT,
    activity_value REAL,
    source TEXT,
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
);
CREATE INDEX idx_ct_pubchem ON compound_targets(pubchem_cid);
CREATE INDEX idx_ct_uniprot ON compound_targets(uniprot_id);

-- protein_pathway associations
CREATE TABLE protein_pathways (
    uniprot_id TEXT,
    pathway_id TEXT,
    evidence TEXT,
    PRIMARY KEY (uniprot_id, pathway_id),
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id),
    FOREIGN KEY (pathway_id) REFERENCES pathways(pathway_id)
);
CREATE INDEX idx_pp_pathway ON protein_pathways(pathway_id);

-- compound_pathway (direct)
CREATE TABLE compound_pathways (
    compound_id TEXT,  -- chebi or kegg
    pubchem_cid INTEGER,
    pathway_id TEXT,
    role TEXT,  -- substrate, product, inhibitor, etc.
    source TEXT,
    PRIMARY KEY (compound_id, pathway_id),
    FOREIGN KEY (pathway_id) REFERENCES pathways(pathway_id)
);
CREATE INDEX idx_cp_pubchem ON compound_pathways(pubchem_cid);
```

**Load data script:**
```python
import sqlite3
import gzip
import csv
import os

def create_database(db_path: str):
    """Create and populate the mapping database."""

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create schema
    cursor.executescript("""
        -- [schema from above]
    """)

    # Load HGNC genes
    print("Loading HGNC...")
    with open("data/genes/hgnc_complete_set.txt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cursor.execute("""
                INSERT OR REPLACE INTO genes
                (hgnc_id, symbol, name, entrez_id, ensembl_id, uniprot_ids)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (
                row["hgnc_id"],
                row["symbol"],
                row["name"],
                row.get("entrez_id"),
                row.get("ensembl_gene_id"),
                row.get("uniprot_ids")
            ))

    # Load UniProt ID mapping
    print("Loading UniProt mappings...")
    with gzip.open("data/proteins/HUMAN_9606_idmapping.dat.gz", "rt") as f:
        proteins = {}
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 3:
                uniprot_id, id_type, value = parts
                if uniprot_id not in proteins:
                    proteins[uniprot_id] = {}
                proteins[uniprot_id][id_type] = value

        for uniprot_id, data in proteins.items():
            cursor.execute("""
                INSERT OR REPLACE INTO proteins
                (uniprot_id, gene_symbol, entrez_id, ensembl_gene, ensembl_protein)
                VALUES (?, ?, ?, ?, ?)
            """, (
                uniprot_id,
                data.get("Gene_Name"),
                data.get("GeneID"),
                data.get("Ensembl"),
                data.get("Ensembl_PRO")
            ))

    # Load Reactome pathways
    print("Loading Reactome...")
    with open("data/pathways/ReactomePathways.txt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                cursor.execute("""
                    INSERT OR REPLACE INTO pathways (pathway_id, name, species, source)
                    VALUES (?, ?, ?, 'reactome')
                """, (parts[0], parts[1], parts[2]))

    # Load UniProt-Reactome mappings
    print("Loading protein-pathway mappings...")
    with open("data/pathways/UniProt2Reactome_All_Levels.txt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2 and parts[5] == "Homo sapiens":
                cursor.execute("""
                    INSERT OR IGNORE INTO protein_pathways (uniprot_id, pathway_id, evidence)
                    VALUES (?, ?, ?)
                """, (parts[0], parts[1], parts[4] if len(parts) > 4 else None))

    # Load ChEBI-Reactome mappings
    print("Loading compound-pathway mappings...")
    with open("data/pathways/ChEBI2Reactome_All_Levels.txt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2 and "Homo sapiens" in line:
                cursor.execute("""
                    INSERT OR IGNORE INTO compound_pathways
                    (compound_id, pathway_id, source)
                    VALUES (?, ?, 'reactome')
                """, (parts[0], parts[1]))

    conn.commit()
    conn.close()
    print("Database created successfully!")

if __name__ == "__main__":
    create_database("compound_pathway_mappings.db")
```

### 9.3 Query Strategies

**Common queries:**

```sql
-- Get all pathways for a drug (by PubChem CID)
SELECT DISTINCT p.pathway_id, p.name, 'via_target' as link_type
FROM compound_targets ct
JOIN protein_pathways pp ON ct.uniprot_id = pp.uniprot_id
JOIN pathways p ON pp.pathway_id = p.pathway_id
WHERE ct.pubchem_cid = 2244  -- Aspirin

UNION

SELECT DISTINCT p.pathway_id, p.name, 'direct' as link_type
FROM compound_pathways cp
JOIN pathways p ON cp.pathway_id = p.pathway_id
WHERE cp.pubchem_cid = 2244;

-- Get all drugs affecting a pathway
SELECT DISTINCT c.name, c.pubchem_cid, ct.activity_type, ct.activity_value
FROM pathways p
JOIN protein_pathways pp ON p.pathway_id = pp.pathway_id
JOIN compound_targets ct ON pp.uniprot_id = ct.uniprot_id
JOIN compounds c ON ct.pubchem_cid = c.pubchem_cid
WHERE p.pathway_id = 'R-HSA-109582';  -- Hemostasis

-- Get genes/proteins in pathways affected by a compound
SELECT DISTINCT
    g.symbol,
    g.name,
    prot.uniprot_id,
    p.name as pathway_name
FROM compound_targets ct
JOIN proteins prot ON ct.uniprot_id = prot.uniprot_id
JOIN genes g ON prot.gene_symbol = g.symbol
JOIN protein_pathways pp ON prot.uniprot_id = pp.uniprot_id
JOIN pathways p ON pp.pathway_id = p.pathway_id
WHERE ct.chembl_id = 'CHEMBL25';  -- Aspirin

-- Cross-database ID lookup
SELECT
    c.pubchem_cid,
    c.chembl_id,
    c.drugbank_id,
    c.name
FROM compounds c
WHERE c.inchikey = 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N';  -- Aspirin InChIKey
```

**Index optimization:**
```sql
-- For frequent pathway queries
CREATE INDEX idx_pp_combined ON protein_pathways(pathway_id, uniprot_id);

-- For compound lookup
CREATE INDEX idx_ct_combined ON compound_targets(pubchem_cid, uniprot_id);

-- Analyze for query optimization
ANALYZE;
```

---

## 10. Code Examples

### 10.1 Python Functions for Each Link Type

```python
"""
Compound-Pathway Linking Library
================================
Functions for linking compounds to genes and pathways via multiple databases.
"""

import requests
import time
from typing import List, Dict, Optional, Any
from dataclasses import dataclass
from functools import lru_cache
import json


@dataclass
class CompoundInfo:
    """Compound information container."""
    pubchem_cid: Optional[int] = None
    chembl_id: Optional[str] = None
    name: Optional[str] = None
    smiles: Optional[str] = None
    inchikey: Optional[str] = None


@dataclass
class TargetInfo:
    """Target protein information container."""
    uniprot_id: str
    gene_symbol: Optional[str] = None
    protein_name: Optional[str] = None
    activity_type: Optional[str] = None
    activity_value: Optional[float] = None


@dataclass
class PathwayInfo:
    """Pathway information container."""
    pathway_id: str
    name: str
    source: str  # reactome, kegg, go, wikipathways


# ============================================================
# 1. Compound Identifier Resolution
# ============================================================

class CompoundResolver:
    """Resolve compound identifiers across databases."""

    @staticmethod
    def pubchem_to_chembl(pubchem_cid: int) -> Optional[str]:
        """Convert PubChem CID to ChEMBL ID."""
        url = f"https://www.ebi.ac.uk/unichem/rest/src_compound_id/{pubchem_cid}/22"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            for entry in data:
                if entry.get("src_id") == "1":  # ChEMBL
                    return f"CHEMBL{entry['src_compound_id']}"
        return None

    @staticmethod
    def name_to_pubchem(name: str) -> Optional[int]:
        """Resolve compound name to PubChem CID."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            return cids[0] if cids else None
        return None

    @staticmethod
    def smiles_to_pubchem(smiles: str) -> Optional[int]:
        """Resolve SMILES to PubChem CID."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{requests.utils.quote(smiles)}/cids/JSON"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            return cids[0] if cids else None
        return None

    @staticmethod
    def inchikey_to_pubchem(inchikey: str) -> Optional[int]:
        """Resolve InChIKey to PubChem CID."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            return cids[0] if cids else None
        return None

    @staticmethod
    def get_compound_info(identifier: str) -> CompoundInfo:
        """Get comprehensive compound info from any identifier."""
        info = CompoundInfo()

        # Try to identify the type of identifier
        if identifier.isdigit():
            info.pubchem_cid = int(identifier)
        elif identifier.startswith("CHEMBL"):
            info.chembl_id = identifier
        elif len(identifier) == 27 and "-" in identifier:
            info.inchikey = identifier
            info.pubchem_cid = CompoundResolver.inchikey_to_pubchem(identifier)
        else:
            # Assume it's a name
            info.name = identifier
            info.pubchem_cid = CompoundResolver.name_to_pubchem(identifier)

        # Enrich with PubChem data
        if info.pubchem_cid:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{info.pubchem_cid}/property/CanonicalSMILES,InChIKey,Title/JSON"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
                info.smiles = props.get("CanonicalSMILES")
                info.inchikey = props.get("InChIKey")
                info.name = props.get("Title")

            # Get ChEMBL ID
            if not info.chembl_id:
                info.chembl_id = CompoundResolver.pubchem_to_chembl(info.pubchem_cid)

        return info


# ============================================================
# 2. Compound → Target Links
# ============================================================

class CompoundTargetLinker:
    """Link compounds to protein targets."""

    @staticmethod
    def chembl_targets(chembl_id: str, pchembl_threshold: float = 6.0) -> List[TargetInfo]:
        """Get targets from ChEMBL with activity data."""
        targets = []

        url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
        params = {
            "molecule_chembl_id": chembl_id,
            "pchembl_value__gte": pchembl_threshold,
            "limit": 100
        }

        response = requests.get(url, params=params)
        if response.status_code != 200:
            return targets

        data = response.json()
        seen_targets = set()

        for activity in data.get("activities", []):
            target_id = activity.get("target_chembl_id")
            if not target_id or target_id in seen_targets:
                continue

            seen_targets.add(target_id)

            # Get target details
            target_url = f"https://www.ebi.ac.uk/chembl/api/data/target/{target_id}.json"
            target_response = requests.get(target_url)

            if target_response.status_code == 200:
                target_data = target_response.json()
                for component in target_data.get("target_components", []):
                    if component.get("accession"):
                        targets.append(TargetInfo(
                            uniprot_id=component["accession"],
                            gene_symbol=component.get("gene_name"),
                            protein_name=target_data.get("pref_name"),
                            activity_type=activity.get("standard_type"),
                            activity_value=activity.get("pchembl_value")
                        ))

        return targets

    @staticmethod
    def dgidb_targets(gene_or_drug: str) -> List[Dict[str, Any]]:
        """Get drug-gene interactions from DGIdb."""
        # Determine if input is gene or drug
        if gene_or_drug.isupper() and len(gene_or_drug) < 10:
            param = "genes"
        else:
            param = "drugs"

        url = f"https://dgidb.org/api/v2/interactions.json?{param}={gene_or_drug}"
        response = requests.get(url)

        if response.status_code != 200:
            return []

        data = response.json()
        interactions = []

        for match in data.get("matchedTerms", []):
            for interaction in match.get("interactions", []):
                interactions.append({
                    "gene": interaction.get("geneName"),
                    "drug": interaction.get("drugName"),
                    "types": interaction.get("interactionTypes", []),
                    "sources": interaction.get("sources", [])
                })

        return interactions

    @staticmethod
    def stitch_targets(pubchem_cid: int, score_threshold: int = 700) -> List[Dict[str, Any]]:
        """Get protein interactions from STITCH."""
        # STITCH uses its own CID format
        stitch_cid = f"CIDs{pubchem_cid:08d}"

        url = f"http://stitch.embl.de/api/json/interactors?identifiers={stitch_cid}&species=9606"
        response = requests.get(url)

        if response.status_code != 200:
            return []

        interactions = []
        for item in response.json():
            if item.get("score", 0) >= score_threshold / 1000:
                interactions.append({
                    "string_id": item.get("stringId"),
                    "name": item.get("preferredName"),
                    "score": item.get("score")
                })

        return interactions


# ============================================================
# 3. Protein → Gene Links
# ============================================================

class ProteinGeneLinker:
    """Link proteins to genes."""

    @staticmethod
    @lru_cache(maxsize=1000)
    def uniprot_to_gene(uniprot_id: str) -> Dict[str, Any]:
        """Get gene information for a UniProt ID."""
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        headers = {"Accept": "application/json"}

        response = requests.get(url, headers=headers)
        if response.status_code != 200:
            return {}

        data = response.json()

        result = {
            "uniprot_id": uniprot_id,
            "gene_symbol": None,
            "gene_name": None,
            "entrez_id": None,
            "ensembl_gene": None,
            "hgnc_id": None
        }

        # Extract gene info
        genes = data.get("genes", [])
        if genes:
            result["gene_symbol"] = genes[0].get("geneName", {}).get("value")

        # Extract cross-references
        for xref in data.get("uniProtKBCrossReferences", []):
            db = xref.get("database")
            if db == "GeneID":
                result["entrez_id"] = xref.get("id")
            elif db == "Ensembl":
                for prop in xref.get("properties", []):
                    if prop.get("key") == "GeneId":
                        result["ensembl_gene"] = prop.get("value")
            elif db == "HGNC":
                result["hgnc_id"] = xref.get("id")

        return result

    @staticmethod
    def batch_uniprot_to_gene(uniprot_ids: List[str]) -> Dict[str, Dict[str, Any]]:
        """Batch convert UniProt IDs to gene information."""
        results = {}

        # Use UniProt ID mapping service
        url = "https://rest.uniprot.org/idmapping/run"
        data = {
            "from": "UniProtKB_AC-ID",
            "to": "GeneID",
            "ids": ",".join(uniprot_ids)
        }

        response = requests.post(url, data=data)
        if response.status_code != 200:
            # Fallback to individual queries
            for uid in uniprot_ids:
                results[uid] = ProteinGeneLinker.uniprot_to_gene(uid)
            return results

        job_id = response.json()["jobId"]

        # Poll for results
        while True:
            status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
            status = requests.get(status_url).json()

            if "results" in status or status.get("jobStatus") == "FINISHED":
                break
            time.sleep(0.5)

        # Get results
        result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
        result_data = requests.get(result_url).json()

        for item in result_data.get("results", []):
            from_id = item.get("from")
            to_id = item.get("to")
            results[from_id] = {"uniprot_id": from_id, "entrez_id": to_id}

        return results

    @staticmethod
    def hgnc_lookup(symbol: str) -> Dict[str, Any]:
        """Look up gene information from HGNC."""
        url = f"https://rest.genenames.org/fetch/symbol/{symbol}"
        headers = {"Accept": "application/json"}

        response = requests.get(url, headers=headers)
        if response.status_code != 200:
            return {}

        data = response.json()
        docs = data.get("response", {}).get("docs", [])

        if not docs:
            return {}

        gene = docs[0]
        return {
            "hgnc_id": gene.get("hgnc_id"),
            "symbol": gene.get("symbol"),
            "name": gene.get("name"),
            "entrez_id": gene.get("entrez_id"),
            "ensembl_gene": gene.get("ensembl_gene_id"),
            "uniprot_ids": gene.get("uniprot_ids", [])
        }


# ============================================================
# 4. Gene/Protein → Pathway Links
# ============================================================

class GenePathwayLinker:
    """Link genes/proteins to pathways."""

    @staticmethod
    def reactome_pathways(uniprot_id: str) -> List[PathwayInfo]:
        """Get Reactome pathways for a UniProt ID."""
        url = f"https://reactome.org/ContentService/data/pathways/low/diagram/entity/UniProt:{uniprot_id}"
        headers = {"Accept": "application/json"}

        response = requests.get(url, headers=headers)
        if response.status_code != 200:
            return []

        pathways = []
        for p in response.json():
            pathways.append(PathwayInfo(
                pathway_id=p.get("stId"),
                name=p.get("displayName"),
                source="reactome"
            ))

        return pathways

    @staticmethod
    def kegg_pathways(entrez_id: str) -> List[PathwayInfo]:
        """Get KEGG pathways for an Entrez Gene ID."""
        url = f"https://rest.kegg.jp/link/pathway/hsa:{entrez_id}"

        response = requests.get(url)
        if response.status_code != 200 or not response.text.strip():
            return []

        pathways = []
        pathway_ids = []

        for line in response.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                pathway_id = parts[1].replace("path:", "")
                pathway_ids.append(pathway_id)

        # Get pathway names
        if pathway_ids:
            names_url = f"https://rest.kegg.jp/list/{'+'.join(pathway_ids)}"
            names_response = requests.get(names_url)

            if names_response.status_code == 200:
                for line in names_response.text.strip().split("\n"):
                    parts = line.split("\t")
                    if len(parts) == 2:
                        pid = parts[0]
                        name = parts[1].split(" - ")[0] if " - " in parts[1] else parts[1]
                        pathways.append(PathwayInfo(
                            pathway_id=pid,
                            name=name,
                            source="kegg"
                        ))

        return pathways

    @staticmethod
    def go_terms(uniprot_id: str, aspect: str = "P") -> List[PathwayInfo]:
        """Get GO terms for a UniProt ID (P=Process, F=Function, C=Component)."""
        url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
        params = {
            "geneProductId": uniprot_id,
            "aspect": "biological_process" if aspect == "P" else
                      "molecular_function" if aspect == "F" else "cellular_component",
            "limit": 100
        }

        response = requests.get(url, params=params)
        if response.status_code != 200:
            return []

        data = response.json()
        terms = {}

        for result in data.get("results", []):
            go_id = result.get("goId")
            if go_id and go_id not in terms:
                terms[go_id] = result.get("goName", "")

        return [
            PathwayInfo(pathway_id=go_id, name=name, source="go")
            for go_id, name in terms.items()
        ]

    @staticmethod
    def wikipathways(entrez_id: str) -> List[PathwayInfo]:
        """Get WikiPathways for an Entrez Gene ID."""
        url = f"https://webservice.wikipathways.org/findPathwaysByXref?ids={entrez_id}&codes=L&format=json"

        response = requests.get(url)
        if response.status_code != 200:
            return []

        data = response.json()
        pathways = []

        for result in data.get("result", []):
            pathways.append(PathwayInfo(
                pathway_id=result.get("id"),
                name=result.get("name"),
                source="wikipathways"
            ))

        return pathways


# ============================================================
# 5. Compound → Pathway Direct Links
# ============================================================

class CompoundPathwayLinker:
    """Link compounds directly to pathways."""

    @staticmethod
    def reactome_pathways(chebi_id: str) -> List[PathwayInfo]:
        """Get Reactome pathways containing a ChEBI compound."""
        url = f"https://reactome.org/ContentService/data/pathways/low/diagram/entity/ChEBI:{chebi_id}"
        headers = {"Accept": "application/json"}

        response = requests.get(url, headers=headers)
        if response.status_code != 200:
            return []

        return [
            PathwayInfo(
                pathway_id=p.get("stId"),
                name=p.get("displayName"),
                source="reactome"
            )
            for p in response.json()
        ]

    @staticmethod
    def kegg_pathways(kegg_compound_id: str) -> List[PathwayInfo]:
        """Get KEGG pathways containing a KEGG compound."""
        url = f"https://rest.kegg.jp/link/pathway/{kegg_compound_id}"

        response = requests.get(url)
        if response.status_code != 200 or not response.text.strip():
            return []

        pathways = []
        for line in response.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                pathway_id = parts[1].replace("path:", "")
                # Filter to human pathways (hsa prefix)
                if pathway_id.startswith("hsa") or pathway_id.startswith("map"):
                    pathways.append(PathwayInfo(
                        pathway_id=pathway_id,
                        name="",  # Would need separate lookup
                        source="kegg"
                    ))

        return pathways

    @staticmethod
    def pubchem_to_chebi(pubchem_cid: int) -> Optional[str]:
        """Convert PubChem CID to ChEBI ID."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/xrefs/RegistryID/JSON"

        response = requests.get(url)
        if response.status_code != 200:
            return None

        data = response.json()
        for registry in data.get("InformationList", {}).get("Information", [{}])[0].get("RegistryID", []):
            if registry.startswith("CHEBI:"):
                return registry.replace("CHEBI:", "")

        return None

    @staticmethod
    def pubchem_to_kegg(pubchem_cid: int) -> Optional[str]:
        """Convert PubChem CID to KEGG Compound ID."""
        # Use KEGG's conversion API
        url = f"https://rest.kegg.jp/conv/compound/pubchem:{pubchem_cid}"

        response = requests.get(url)
        if response.status_code == 200 and response.text.strip():
            parts = response.text.strip().split("\t")
            if len(parts) == 2:
                return parts[1]

        return None


# ============================================================
# 6. Complete Pipeline
# ============================================================

class CompoundPathwayPipeline:
    """Complete pipeline for compound → pathway linking."""

    def __init__(self):
        self.resolver = CompoundResolver()
        self.target_linker = CompoundTargetLinker()
        self.gene_linker = ProteinGeneLinker()
        self.pathway_linker = GenePathwayLinker()
        self.compound_pathway_linker = CompoundPathwayLinker()

    def compound_to_pathways(self, identifier: str,
                            include_direct: bool = True,
                            include_via_targets: bool = True) -> Dict[str, Any]:
        """
        Complete compound → pathway mapping.

        Args:
            identifier: Compound name, PubChem CID, ChEMBL ID, InChIKey, or SMILES
            include_direct: Include direct compound-pathway links
            include_via_targets: Include pathways via protein targets

        Returns:
            Dictionary with compound info, targets, and pathways
        """
        result = {
            "compound": None,
            "targets": [],
            "pathways_direct": [],
            "pathways_via_targets": []
        }

        # Step 1: Resolve compound
        print(f"Resolving compound: {identifier}")
        compound_info = self.resolver.get_compound_info(identifier)
        result["compound"] = {
            "pubchem_cid": compound_info.pubchem_cid,
            "chembl_id": compound_info.chembl_id,
            "name": compound_info.name,
            "smiles": compound_info.smiles,
            "inchikey": compound_info.inchikey
        }

        # Step 2: Get direct pathway links
        if include_direct and compound_info.pubchem_cid:
            print("Getting direct compound-pathway links...")

            # Try ChEBI → Reactome
            chebi_id = self.compound_pathway_linker.pubchem_to_chebi(compound_info.pubchem_cid)
            if chebi_id:
                reactome_pathways = self.compound_pathway_linker.reactome_pathways(chebi_id)
                result["pathways_direct"].extend([
                    {"id": p.pathway_id, "name": p.name, "source": "reactome"}
                    for p in reactome_pathways
                ])

            # Try KEGG Compound → KEGG Pathway
            kegg_id = self.compound_pathway_linker.pubchem_to_kegg(compound_info.pubchem_cid)
            if kegg_id:
                kegg_pathways = self.compound_pathway_linker.kegg_pathways(kegg_id)
                result["pathways_direct"].extend([
                    {"id": p.pathway_id, "name": p.name, "source": "kegg"}
                    for p in kegg_pathways
                ])

        # Step 3: Get targets and their pathways
        if include_via_targets and compound_info.chembl_id:
            print("Getting protein targets...")
            targets = self.target_linker.chembl_targets(compound_info.chembl_id)

            for target in targets[:10]:  # Limit to top 10 targets
                target_data = {
                    "uniprot_id": target.uniprot_id,
                    "gene_symbol": target.gene_symbol,
                    "activity_type": target.activity_type,
                    "activity_value": target.activity_value,
                    "pathways": []
                }

                # Get gene info
                gene_info = self.gene_linker.uniprot_to_gene(target.uniprot_id)
                target_data["entrez_id"] = gene_info.get("entrez_id")

                # Get pathways
                print(f"  Getting pathways for {target.gene_symbol or target.uniprot_id}...")

                # Reactome
                reactome = self.pathway_linker.reactome_pathways(target.uniprot_id)
                target_data["pathways"].extend([
                    {"id": p.pathway_id, "name": p.name, "source": "reactome"}
                    for p in reactome[:5]
                ])

                # KEGG (if we have Entrez ID)
                if gene_info.get("entrez_id"):
                    kegg = self.pathway_linker.kegg_pathways(gene_info["entrez_id"])
                    target_data["pathways"].extend([
                        {"id": p.pathway_id, "name": p.name, "source": "kegg"}
                        for p in kegg[:5]
                    ])

                result["targets"].append(target_data)

                # Aggregate pathways
                for pathway in target_data["pathways"]:
                    pathway_with_target = pathway.copy()
                    pathway_with_target["via_target"] = target.uniprot_id
                    result["pathways_via_targets"].append(pathway_with_target)

        return result

    def gene_to_compounds(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Map gene to affecting compounds.

        Args:
            gene_symbol: HGNC gene symbol

        Returns:
            Dictionary with gene info, protein, pathways, and affecting drugs
        """
        result = {
            "gene": None,
            "protein": None,
            "pathways": [],
            "compounds": []
        }

        # Step 1: Get gene info
        print(f"Looking up gene: {gene_symbol}")
        gene_info = self.gene_linker.hgnc_lookup(gene_symbol)
        result["gene"] = gene_info

        # Step 2: Get protein info
        if gene_info.get("uniprot_ids"):
            uniprot_id = gene_info["uniprot_ids"][0]
            result["protein"] = {"uniprot_id": uniprot_id}

            # Step 3: Get pathways
            print("Getting pathways...")
            reactome = self.pathway_linker.reactome_pathways(uniprot_id)
            result["pathways"].extend([
                {"id": p.pathway_id, "name": p.name, "source": "reactome"}
                for p in reactome
            ])

            if gene_info.get("entrez_id"):
                kegg = self.pathway_linker.kegg_pathways(str(gene_info["entrez_id"]))
                result["pathways"].extend([
                    {"id": p.pathway_id, "name": p.name, "source": "kegg"}
                    for p in kegg
                ])

        # Step 4: Get compounds from DGIdb
        print("Getting drug interactions...")
        interactions = self.target_linker.dgidb_targets(gene_symbol)
        result["compounds"] = [
            {
                "drug_name": i["drug"],
                "interaction_types": i["types"],
                "sources": i["sources"]
            }
            for i in interactions
        ]

        return result


# ============================================================
# Usage Examples
# ============================================================

if __name__ == "__main__":
    pipeline = CompoundPathwayPipeline()

    # Example 1: Compound → Pathways
    print("\n" + "="*60)
    print("Example 1: Curcumin → Pathways")
    print("="*60)

    result = pipeline.compound_to_pathways("curcumin")

    print(f"\nCompound: {result['compound']['name']}")
    print(f"PubChem CID: {result['compound']['pubchem_cid']}")
    print(f"ChEMBL ID: {result['compound']['chembl_id']}")

    print(f"\nDirect pathways: {len(result['pathways_direct'])}")
    for p in result['pathways_direct'][:5]:
        print(f"  - [{p['source']}] {p['id']}: {p['name']}")

    print(f"\nTargets: {len(result['targets'])}")
    for t in result['targets'][:3]:
        print(f"  - {t['gene_symbol']} ({t['uniprot_id']})")
        print(f"    Activity: {t['activity_type']} = {t['activity_value']}")
        print(f"    Pathways: {len(t['pathways'])}")

    # Example 2: Gene → Compounds
    print("\n" + "="*60)
    print("Example 2: CYP2D6 → Compounds")
    print("="*60)

    result = pipeline.gene_to_compounds("CYP2D6")

    print(f"\nGene: {result['gene']['symbol']}")
    print(f"HGNC ID: {result['gene']['hgnc_id']}")
    print(f"Entrez ID: {result['gene']['entrez_id']}")

    print(f"\nPathways: {len(result['pathways'])}")
    for p in result['pathways'][:5]:
        print(f"  - [{p['source']}] {p['name']}")

    print(f"\nAffecting compounds: {len(result['compounds'])}")
    for c in result['compounds'][:10]:
        print(f"  - {c['drug_name']}: {', '.join(c['interaction_types'])}")
```

### 10.2 Complete Pipeline Script

Save as `compound_pathway_link.py`:

```python
#!/usr/bin/env python3
"""
Compound-Pathway Linking CLI Tool
=================================
Command-line tool for linking compounds to pathways.

Usage:
    python compound_pathway_link.py compound "aspirin" --output json
    python compound_pathway_link.py gene "CYP2D6" --output table
    python compound_pathway_link.py pathway "R-HSA-211981" --output json
"""

import argparse
import json
import sys
from typing import Dict, Any

# Import the library functions (from previous section)
# from compound_pathway_lib import CompoundPathwayPipeline


def format_table(data: Dict[str, Any], indent: int = 0) -> str:
    """Format data as readable table."""
    lines = []
    prefix = "  " * indent

    for key, value in data.items():
        if isinstance(value, dict):
            lines.append(f"{prefix}{key}:")
            lines.append(format_table(value, indent + 1))
        elif isinstance(value, list):
            lines.append(f"{prefix}{key}: ({len(value)} items)")
            for i, item in enumerate(value[:10]):
                if isinstance(item, dict):
                    lines.append(f"{prefix}  [{i+1}]")
                    lines.append(format_table(item, indent + 2))
                else:
                    lines.append(f"{prefix}  - {item}")
            if len(value) > 10:
                lines.append(f"{prefix}  ... and {len(value) - 10} more")
        else:
            lines.append(f"{prefix}{key}: {value}")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Link compounds to genes and pathways"
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Compound command
    compound_parser = subparsers.add_parser(
        "compound", help="Map compound to pathways"
    )
    compound_parser.add_argument(
        "identifier", help="Compound name, PubChem CID, ChEMBL ID, or InChIKey"
    )
    compound_parser.add_argument(
        "--no-direct", action="store_true", help="Skip direct pathway links"
    )
    compound_parser.add_argument(
        "--no-targets", action="store_true", help="Skip target-mediated links"
    )
    compound_parser.add_argument(
        "--output", choices=["json", "table"], default="table"
    )

    # Gene command
    gene_parser = subparsers.add_parser(
        "gene", help="Map gene to compounds and pathways"
    )
    gene_parser.add_argument("symbol", help="HGNC gene symbol")
    gene_parser.add_argument(
        "--output", choices=["json", "table"], default="table"
    )

    # Pathway command
    pathway_parser = subparsers.add_parser(
        "pathway", help="Get compounds and genes in pathway"
    )
    pathway_parser.add_argument(
        "pathway_id", help="Pathway ID (Reactome, KEGG, or GO)"
    )
    pathway_parser.add_argument(
        "--output", choices=["json", "table"], default="table"
    )

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    pipeline = CompoundPathwayPipeline()

    try:
        if args.command == "compound":
            result = pipeline.compound_to_pathways(
                args.identifier,
                include_direct=not args.no_direct,
                include_via_targets=not args.no_targets
            )
        elif args.command == "gene":
            result = pipeline.gene_to_compounds(args.symbol)
        elif args.command == "pathway":
            # Implement pathway → compounds/genes lookup
            result = {"error": "Pathway lookup not yet implemented"}

        if args.output == "json":
            print(json.dumps(result, indent=2))
        else:
            print(format_table(result))

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
```

---

## Quick Reference

### ID Conversion Cheat Sheet

| From | To | Best Method |
|------|-----|-------------|
| PubChem CID | ChEMBL | UniChem API |
| PubChem CID | ChEBI | PubChem xrefs |
| PubChem CID | KEGG | KEGG conv API |
| ChEMBL | PubChem | UniChem API |
| UniProt | Entrez Gene | UniProt ID mapping |
| UniProt | Ensembl | UniProt ID mapping |
| Gene Symbol | UniProt | HGNC API |
| Gene Symbol | Entrez | HGNC API |
| Entrez Gene | UniProt | UniProt ID mapping |

### API Endpoints Summary

| Service | Base URL | Rate Limit |
|---------|----------|------------|
| PubChem PUG REST | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/` | 5 req/sec |
| ChEMBL API | `https://www.ebi.ac.uk/chembl/api/data/` | No limit |
| UniProt REST | `https://rest.uniprot.org/` | No limit |
| Reactome Content | `https://reactome.org/ContentService/` | No limit |
| KEGG REST | `https://rest.kegg.jp/` | 10 req/sec |
| BridgeDb | `https://webservice.bridgedb.org/` | No limit |
| DGIdb | `https://dgidb.org/api/v2/` | No limit |
| HGNC REST | `https://rest.genenames.org/` | No limit |

### Common Data Files

| Data | URL | Size |
|------|-----|------|
| HGNC complete | `genenames.org/download/archive/` | ~15 MB |
| UniProt human IDmapping | `ftp.uniprot.org/.../HUMAN_9606_idmapping.dat.gz` | ~200 MB |
| Reactome UniProt mapping | `reactome.org/download/.../UniProt2Reactome_All_Levels.txt` | ~50 MB |
| ChEMBL SQLite | `ftp.ebi.ac.uk/.../chembl_*_sqlite.tar.gz` | ~4 GB |
| GO human annotations | `current.geneontology.org/annotations/goa_human.gaf.gz` | ~50 MB |

---

## Troubleshooting

### Common Issues

1. **Rate limiting**: Add delays between requests (`time.sleep(0.2)`)
2. **ID not found**: Try alternative ID types or synonyms
3. **Empty results**: Check species filter (human = 9606)
4. **Timeout**: Increase timeout, use bulk downloads for large queries

### Validation Tips

1. **Cross-validate** IDs using multiple sources
2. **Check species** consistency (human UniProt starts with P, Q, or O)
3. **Verify activity values** are in expected ranges
4. **Use canonical** identifiers (reviewed UniProt, primary ChEMBL)

---

## See Also

- [PubChem Documentation](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
- [ChEMBL Web Services](https://chembl.gitbook.io/chembl-interface-documentation/)
- [UniProt Help](https://www.uniprot.org/help/)
- [Reactome Developer Guide](https://reactome.org/dev)
- [KEGG API Documentation](https://www.kegg.jp/kegg/rest/keggapi.html)

---

## Download

| Database | Method | URL |
|----------|--------|-----|
| PubChem | REST API | https://pubchem.ncbi.nlm.nih.gov/rest/pug/ |
| PubChem | ID Exchange | https://pubchem.ncbi.nlm.nih.gov/idexchange/ |
| ChEMBL | REST API | https://www.ebi.ac.uk/chembl/api/ |
| ChEMBL | SQL Dump | https://ftp.ebi.ac.uk/pub/databases/chembl/ |
| Reactome | REST API | https://reactome.org/ContentService/ |
| UniProt | REST API | https://rest.uniprot.org/ |
| UniProt | ID Mapping | https://www.uniprot.org/id-mapping/ |
| KEGG | REST API | https://www.kegg.jp/kegg/rest/keggapi.html |
| STRING | REST API | https://string-db.org/api/ |

---

## Data Format

| Format | Description | When Used |
|--------|-------------|-----------|
| JSON | Hierarchical structure, nested objects | API responses, linked data |
| TSV | Tab-separated records | Bulk downloads, line-oriented data |
| CSV | Comma-separated values | Spreadsheet-compatible format |
| SPARQL | Query language for RDF data | Federated queries across resources |
| XML | Markup structure with tags | Complex nested relationships |
| OWL/RDF | Semantic web ontologies | Ontology-based reasoning |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| compound_id | String | Primary compound identifier | CID:2662 |
| compound_name | String | IUPAC or common name | Curcumin |
| target_id | String | Protein target identifier | P35354 (UniProt) |
| target_name | String | Gene/protein symbol | COX2 |
| interaction_type | String | Nature of interaction | "inhibitor", "activator", "binder" |
| activity_value | Float | Quantitative measure of interaction | 5.2 |
| activity_unit | String | Unit of measurement | "μM", "nM", "IC50" |
| target_species | String | Organism of target | "Homo sapiens" |
| confidence_score | Float | Prediction confidence (0-1) | 0.85 |
| source_database | String | Data source | "ChEMBL", "PubChem", "KEGG" |

### Relationships

| From | To | Via | Field |
|------|-----|-----|-------|
| Compound | Target Protein | ChEMBL, PubChem | compound_id → target_id |
| Target Protein | Gene | UniProt | UniProt_ID → Ensembl_Gene |
| Gene | Pathway | Reactome, KEGG | Ensembl_ID → Pathway_ID |
| Compound | Pathway | Indirect | Via Target Protein → Pathway |

---

## Sample Data

### JSON Format

```json
{
  "compound_target_link": {
    "compound": {
      "id": "CHEMBL429661",
      "name": "Curcumin",
      "pubchem_cid": 969516,
      "inchi_key": "GHASVSINZRGABV-UHFFFAOYSA-N"
    },
    "target": {
      "uniprot_id": "P35354",
      "gene_symbol": "COX2",
      "gene_name": "Prostaglandin G/H synthase 2",
      "ensembl_id": "ENSG00000073756"
    },
    "interaction": {
      "type": "inhibitor",
      "activity": {
        "value": 5.2,
        "unit": "μM",
        "assay_type": "IC50"
      },
      "confidence": 0.92,
      "source": "ChEMBL"
    }
  }
}
```

### Query Result Table

| Compound | Target Gene | UniProt | Interaction | Value | Unit | Pathway | Confidence |
|----------|-------------|---------|-------------|-------|------|---------|------------|
| Curcumin | COX2 | P35354 | inhibitor | 5.2 | μM | Eicosanoid synthesis | 0.92 |
| Artemisinin | PfMPT1 | Q9Y5J8 | inhibitor | 2.1 | μM | Parasite survival | 0.88 |
| Resveratrol | SIRT1 | Q96EB6 | activator | 15.0 | μM | Stress response | 0.85 |

---

## License

| Resource | License | Commercial Use | Citation Required |
|----------|---------|-----------------|-------------------|
| PubChem | CC0 1.0 | Yes | Recommended |
| ChEMBL | CC BY 4.0 | Yes | Required |
| Reactome | CC BY 4.0 | Yes | Required |
| UniProt | CC BY 4.0 | Yes | Required |
| KEGG | Academic/Commercial | Yes (paid) | Required |
| STRING | CC BY 4.0 | Yes | Required |

---

## Data Set Size

| Resource | Records | Compressed Size | Uncompressed Size | Format |
|----------|---------|-----------------|-------------------|--------|
| PubChem Compounds | 120M+ | 15 GB | 150 GB | JSON |
| ChEMBL Activities | 15M+ | 2.5 GB | 25 GB | SQL |
| Reactome Pathways | 2,700+ | 500 MB | 3 GB | JSON |
| UniProt (Human) | 21K+ | 2 GB | 12 GB | FASTA/XML |
| KEGG Pathways | 600+ | 200 MB | 1.5 GB | KGML |
| STRING Interactions | 1B+ | 20 GB | 80 GB | TSV |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Compound-Target Linking | Mapping small molecules to their protein targets | Curcumin inhibits COX-2 |
| Pathway Enrichment | Statistical analysis identifying pathways overrepresented in a gene set | Fisher's exact test |
| Bioactivity | Measured effect of a compound on a biological target | IC50 = 100 nM |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| PubChem CID | Compound identifier in PubChem database | CID 969516 (curcumin) |
| ChEMBL ID | Compound identifier in ChEMBL database | CHEMBL116438 |
| InChIKey | 27-character hashed molecular structure identifier | VFLDPWHFBUODDF-FCXRPNKRSA-N |
| SMILES | Simplified Molecular Input Line Entry System notation | Chemical structure string |
| UniProt AC | Accession number identifying a protein in UniProt | P35354 (COX-2) |
| Reactome ID | Pathway identifier in Reactome database | R-HSA-109582 |
| IC50 | Half maximal inhibitory concentration | Potency measure |
| pChEMBL | Standardized -log10 activity value | Activity normalization |
| BridgeDb | ID mapping service for biological databases | Cross-reference tool |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST endpoints |
| CID | Compound Identifier | PubChem ID |
| DGIdb | Drug-Gene Interaction Database | Drug-target data |
| GO | Gene Ontology | Functional annotation |
| HGNC | HUGO Gene Nomenclature Committee | Gene symbols |
| IC50 | Half Maximal Inhibitory Concentration | Activity measure |
| InChI | International Chemical Identifier | Structure encoding |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |
| PPI | Protein-Protein Interaction | Network data |
| SMILES | Simplified Molecular Input Line Entry System | Chemical notation |
