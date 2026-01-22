# Wikipedia and Wikidata as Pharmaceutical Data Sources

**Last Updated:** January 2026
**Purpose:** Research documentation for leveraging Wikipedia/Wikidata ecosystem for drug, gene, and pathway data

---

## Overview

The Wikimedia ecosystem provides several complementary resources for pharmaceutical and biomedical data:

| Resource | Content | License | Access Method |
|----------|---------|---------|---------------|
| **Wikidata** | Structured knowledge graph | CC0 (Public Domain) | SPARQL, API, JSON dumps |
| **Wikipedia** | Article text + infoboxes | CC BY-SA 4.0 | API, dumps |
| **DBpedia** | Wikipedia extracted as RDF | CC BY-SA 3.0 | SPARQL endpoint |
| **WikiPathways** | Biological pathways | CC0 | SPARQL, REST API, downloads |

---

## 1. Wikidata Drug Items

### 1.1 How Drugs Are Represented

In Wikidata, drugs are represented as **items** with unique Q-identifiers (Q-IDs). The primary class for drugs is:

- **Q12140** - pharmaceutical drug (medication)
- **Q11173** - chemical compound (often used together with Q12140)
- **Q422248** - monoclonal antibody (subclass)

Example drug items:
- [Q188724](https://www.wikidata.org/wiki/Q188724) - Ibuprofen
- [Q407431](https://www.wikidata.org/wiki/Q407431) - Metformin
- [Q179452](https://www.wikidata.org/wiki/Q179452) - Amphetamine

### 1.2 Key Properties for Drugs (P-codes)

| Property | Name | Description | Example Value |
|----------|------|-------------|---------------|
| **P715** | DrugBank ID | Identifier in DrugBank database | DB00945 |
| **P231** | CAS Registry Number | Chemical Abstracts Service identifier | 15687-27-1 |
| **P267** | ATC code | WHO Anatomical Therapeutic Chemical code | M01AE01 |
| **P2868** | subject has role | Role of the drug (e.g., inhibitor, substrate) | enzyme inhibitor (Q427087) |
| **P2175** | medical condition treated | Disease/condition the drug treats | hypertension (Q41861) |
| **P129** | physically interacts with | Molecular targets (proteins, enzymes) | COX-2 (Q412415) |
| **P769** | significant drug interaction | Clinically significant drug-drug interactions | warfarin (Q407548) |
| **P176** | manufacturer | Company that produces the drug | Pfizer (Q206921) |
| **P274** | chemical formula | Molecular formula | C13H18O2 |
| **P233** | SMILES | Canonical SMILES notation | CC(C)Cc1ccc(cc1)C(C)C(=O)O |
| **P234** | InChI | IUPAC International Chemical Identifier | InChI=1S/C13H18O2/... |
| **P235** | InChIKey | Hashed InChI for lookup | HEFNNWSXXWATRW-UHFFFAOYSA-N |
| **P486** | MeSH descriptor ID | Medical Subject Headings identifier | D007052 |
| **P652** | FDA UNII | Unique Ingredient Identifier | WK2XYI10QM |
| **P662** | PubChem CID | PubChem Compound ID | 3672 |
| **P683** | ChEBI ID | Chemical Entities of Biological Interest | CHEBI:5855 |
| **P592** | ChEMBL ID | ChEMBL compound identifier | CHEMBL521 |

### 1.3 SPARQL Endpoint

**Endpoint URL:** `https://query.wikidata.org/sparql`

**Web Interface:** [https://query.wikidata.org](https://query.wikidata.org)

**HTTP Access:**
```bash
# GET request (URL-encoded query)
curl -G https://query.wikidata.org/sparql \
  --data-urlencode "query=SELECT * WHERE { ?s ?p ?o } LIMIT 10" \
  -H "Accept: application/sparql-results+json"

# POST request
curl -X POST https://query.wikidata.org/sparql \
  -H "Content-Type: application/sparql-query" \
  -H "Accept: application/sparql-results+json" \
  -d "SELECT * WHERE { ?s ?p ?o } LIMIT 10"
```

**Python Access:**
```python
from SPARQLWrapper import SPARQLWrapper, JSON

endpoint = "https://query.wikidata.org/sparql"
sparql = SPARQLWrapper(endpoint)
sparql.setQuery("""
    SELECT ?drug ?drugLabel WHERE {
        ?drug wdt:P31 wd:Q12140 .
        SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    LIMIT 100
""")
sparql.setReturnFormat(JSON)
results = sparql.query().convert()
```

### 1.4 Bulk Data Access Methods

#### JSON Dumps (Recommended for Large-Scale Analysis)

```bash
# Download latest full dump (~130 GB compressed)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.gz

# Incremental daily updates
wget https://dumps.wikimedia.org/other/incr/wikidatawiki/
```

**Dump Structure:**
- Each line is a separate JSON object (one entity per line)
- Can be processed line-by-line without loading entire file
- Available in JSON, RDF/Turtle, and N-Triples formats

#### Single Entity Access (Linked Data Interface)

```bash
# Get single item as JSON
curl https://www.wikidata.org/wiki/Special:EntityData/Q188724.json

# Available formats: .json, .rdf, .ttl, .nt, .jsonld
```

#### Python Libraries for Dumps

```python
# Using qwikidata
from qwikidata.json_dump import WikidataJsonDump

wjd = WikidataJsonDump("latest-all.json.gz")
for entity in wjd:
    if entity.get("type") == "item":
        claims = entity.get("claims", {})
        if "P267" in claims:  # Has ATC code
            print(entity["id"], entity["labels"].get("en", {}).get("value"))
```

---

## 2. Wikipedia Infobox Drug

### 2.1 Template Structure

The `{{Infobox drug}}` (also called `{{Drugbox}}`) template provides structured data in Wikipedia drug articles.

**Key Fields:**

| Field | Description | Example |
|-------|-------------|---------|
| `drugname` | Drug name | Ibuprofen |
| `IUPAC_name` | IUPAC chemical name | 2-(4-isobutylphenyl)propionic acid |
| `image` | Structure image file | Ibuprofen2.svg |
| `tradename` | Brand names | Advil, Motrin, Nurofen |
| `Drugs.com` | Drugs.com link | monograph/ibuprofen |
| `MedlinePlus` | MedlinePlus ID | a682159 |
| `pregnancy_AU` | Australian pregnancy category | C |
| `pregnancy_US` | US pregnancy category | C (D in 3rd trimester) |
| `routes_of_administration` | Administration routes | Oral, topical, IV |
| `ATC_prefix` | ATC classification prefix | M01 |
| `ATC_suffix` | ATC classification suffix | AE01 |
| `legal_AU` | Australian legal status | S2 |
| `legal_US` | US legal status | OTC / Rx-only |
| `bioavailability` | Oral bioavailability | 80-100% |
| `protein_bound` | Protein binding | 99% |
| `metabolism` | Metabolic pathway | Hepatic (CYP2C9) |
| `elimination_half-life` | Half-life | 2-4 hours |
| `excretion` | Excretion route | Renal |
| `CAS_number` | CAS Registry Number | 15687-27-1 |
| `PubChem` | PubChem CID | 3672 |
| `DrugBank` | DrugBank ID | DB00945 |
| `ChemSpider` | ChemSpider ID | 3544 |
| `UNII` | FDA UNII | WK2XYI10QM |
| `ChEBI` | ChEBI ID | 5855 |
| `ChEMBL` | ChEMBL ID | 521 |
| `chemical_formula` | Molecular formula | C<sub>13</sub>H<sub>18</sub>O<sub>2</sub> |
| `molecular_weight` | Molecular mass | 206.29 g/mol |
| `SMILES` | SMILES string | CC(C)Cc1ccc(cc1)C(C)C(=O)O |
| `InChI` | InChI identifier | 1S/C13H18O2/c1-9(2)... |
| `StdInChIKey` | InChIKey | HEFNNWSXXWATRW-UHFFFAOYSA-N |

### 2.2 Extracting Structured Data

#### Using MediaWiki API

```python
import requests

def get_infobox_data(title):
    """Extract infobox data via MediaWiki API parse action."""
    url = "https://en.wikipedia.org/w/api.php"
    params = {
        "action": "parse",
        "page": title,
        "prop": "wikitext",
        "format": "json"
    }
    response = requests.get(url, params=params)
    wikitext = response.json()["parse"]["wikitext"]["*"]
    # Parse infobox from wikitext
    return wikitext

# Example
wikitext = get_infobox_data("Ibuprofen")
```

#### Using wptools Library (Recommended)

```python
import wptools

def get_drug_infobox(drug_name):
    """Get parsed infobox data using wptools."""
    page = wptools.page(drug_name)
    page.get_parse()
    return page.data.get("infobox", {})

# Example
infobox = get_drug_infobox("Ibuprofen")
print(infobox.get("CAS_number"))  # 15687-27-1
print(infobox.get("DrugBank"))    # DB00945
```

#### Using BeautifulSoup for HTML Parsing

```python
import requests
from bs4 import BeautifulSoup

def get_infobox_html(title):
    """Extract infobox from rendered HTML."""
    url = "https://en.wikipedia.org/w/api.php"
    params = {
        "action": "parse",
        "page": title,
        "prop": "text",
        "format": "json"
    }
    response = requests.get(url, params=params)
    html = response.json()["parse"]["text"]["*"]
    soup = BeautifulSoup(html, "html.parser")
    infobox = soup.find("table", {"class": "infobox"})
    return infobox
```

### 2.3 DBpedia as Structured Wikipedia Access

DBpedia extracts structured data from Wikipedia infoboxes and makes it available as RDF.

**SPARQL Endpoint:** `http://dbpedia.org/sparql`

**Live Endpoint (real-time updates):** `http://live.dbpedia.org/sparql`

**Example Query - Get Drug Information:**

```sparql
PREFIX dbo: <http://dbpedia.org/ontology/>
PREFIX dbr: <http://dbpedia.org/resource/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT ?drug ?name ?casNumber ?drugbank ?atcCode WHERE {
    ?drug a dbo:Drug .
    ?drug rdfs:label ?name .
    OPTIONAL { ?drug dbo:casNumber ?casNumber }
    OPTIONAL { ?drug dbo:drugbank ?drugbank }
    OPTIONAL { ?drug dbo:atcCode ?atcCode }
    FILTER (lang(?name) = "en")
}
LIMIT 100
```

**DBpedia Drug Ontology Classes:**
- `dbo:Drug` - General drug class
- `dbo:ChemicalCompound` - Chemical compounds
- `dbo:ChemicalSubstance` - Chemical substances

---

## 3. WikiPathways

### 3.1 Overview

WikiPathways is a community-curated pathway database containing 3,000+ biological pathways with drug-pathway relationships.

**Website:** [https://www.wikipathways.org](https://www.wikipathways.org)

**License:** CC0 (Public Domain)

### 3.2 Drug-Pathway Relationships

WikiPathways includes:
- Metabolic pathways with drug targets
- Disease pathways showing drug intervention points
- Pharmacokinetic pathways (ADME)
- Drug metabolism pathways (CYP enzymes)

### 3.3 GPML Format

**GPML** (Graphical Pathway Markup Language) is WikiPathways' native XML format.

**Structure:**
```xml
<?xml version="1.0" encoding="UTF-8"?>
<Pathway xmlns="http://pathvisio.org/GPML/2013a"
         Name="Drug Metabolism"
         Organism="Homo sapiens">
  <DataNode TextLabel="CYP3A4" GraphId="abc123" Type="GeneProduct">
    <Xref Database="Ensembl" ID="ENSG00000160868"/>
  </DataNode>
  <DataNode TextLabel="Ibuprofen" GraphId="def456" Type="Metabolite">
    <Xref Database="ChEBI" ID="CHEBI:5855"/>
  </DataNode>
  <Interaction GraphId="ghi789">
    <Graphics ZOrder="12288" LineThickness="1.0">
      <Point X="100" Y="100" GraphRef="abc123" RelX="1.0" RelY="0.0"/>
      <Point X="200" Y="100" GraphRef="def456" RelX="-1.0" RelY="0.0" ArrowHead="mim-catalysis"/>
    </Graphics>
  </Interaction>
</Pathway>
```

### 3.4 API Access

#### REST API

```python
import requests

# Get pathway in GPML format
def get_pathway_gpml(pathway_id):
    """Retrieve pathway in GPML format."""
    url = f"https://webservice.wikipathways.org/getPathway?pwId={pathway_id}&format=json"
    response = requests.get(url)
    return response.json()

# Get cross-references (genes, metabolites)
def get_pathway_xrefs(pathway_id, system_code):
    """
    Get cross-references for a pathway.
    System codes:
      - 'L' = Entrez Gene
      - 'H' = HGNC
      - 'Ce' = ChEBI
      - 'Ik' = InChIKey
    """
    url = f"https://webservice.wikipathways.org/getXrefList?pwId={pathway_id}&code={system_code}&format=json"
    response = requests.get(url)
    return response.json()

# Example
gpml = get_pathway_gpml("WP554")  # ACE Inhibitor Pathway
genes = get_pathway_xrefs("WP554", "H")  # HGNC symbols
metabolites = get_pathway_xrefs("WP554", "Ce")  # ChEBI IDs
```

#### R Package (rWikiPathways)

```r
# Install
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rWikiPathways")

# Usage
library(rWikiPathways)

# Get GPML
gpml <- getPathway("WP554")

# Get metabolites (ChEBI)
metabolites <- getXrefList("WP554", "Ce")

# Get genes (HGNC)
genes <- getXrefList("WP554", "H")

# List all human pathways
human_pathways <- listPathways("Homo sapiens")
```

#### Python Package (pywikipathways)

```python
from pywikipathways import pywikipathways as pwpw

# Get pathway
pathway = pwpw.get_pathway("WP554")

# List pathways by organism
human_pathways = pwpw.list_pathways("Homo sapiens")
```

#### SPARQL Endpoint

**Endpoint:** `https://sparql.wikipathways.org/sparql`

**Example Query - Find Pathways Containing a Gene:**

```sparql
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT DISTINCT ?pathway ?pathwayLabel ?geneLabel
WHERE {
    ?gene a wp:GeneProduct ;
          rdfs:label ?geneLabel ;
          dcterms:isPartOf ?pathway .
    ?pathway a wp:Pathway ;
             rdfs:label ?pathwayLabel .
    FILTER regex(?geneLabel, "CYP3A4", "i")
}
LIMIT 50
```

### 3.5 Bulk Downloads

**Download Page:** [https://www.wikipathways.org/download.html](https://www.wikipathways.org/download.html)

**Available Formats:**
- GPML archives (by organism)
- GMT files (Gene Matrix Transposed) for enrichment analysis
- RDF dumps
- SVG/PNG images

**Zenodo Archive (Citeable):**
- Monthly releases with DOIs
- [https://zenodo.org/communities/wikipathways](https://zenodo.org/communities/wikipathways)

---

## 4. Wikidata Gene/Protein Items

### 4.1 Gene and Protein Representation

In Wikidata, genes and proteins are separate items linked by reciprocal properties:

| Entity Type | Class | Example |
|-------------|-------|---------|
| Gene | Q7187 (gene) | [Q17853226](https://www.wikidata.org/wiki/Q17853226) - BRCA1 (human) |
| Protein | Q8054 (protein) | [Q14864877](https://www.wikidata.org/wiki/Q14864877) - BRCA1 protein |

**Note:** Wikidata maintains separate items for the same gene/protein in different species.

### 4.2 Key Properties for Genes and Proteins

| Property | Name | Direction | Description |
|----------|------|-----------|-------------|
| **P702** | encoded by | Protein -> Gene | The gene that encodes this protein |
| **P688** | encodes | Gene -> Protein | The protein product of this gene |
| **P703** | found in taxon | Both | Species the gene/protein belongs to |
| **P351** | Entrez Gene ID | Gene | NCBI Gene identifier |
| **P353** | HGNC gene symbol | Gene | Official gene symbol |
| **P354** | HGNC ID | Gene | HGNC identifier |
| **P594** | Ensembl gene ID | Gene | Ensembl gene identifier |
| **P352** | UniProt protein ID | Protein | UniProt accession |
| **P637** | RefSeq protein ID | Protein | RefSeq protein accession |
| **P680** | molecular function | Protein | Gene Ontology molecular function |
| **P681** | cell component | Protein | Gene Ontology cellular component |
| **P682** | biological process | Protein | Gene Ontology biological process |

### 4.3 Linking Drugs to Genes

The key property for drug-gene/protein relationships is **P129 (physically interacts with)**.

**Chain: Drug -> Protein Target -> Gene**

```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>

SELECT ?drug ?drugLabel ?protein ?proteinLabel ?gene ?geneLabel
WHERE {
    # Drug physically interacts with protein
    ?drug wdt:P31 wd:Q12140 .           # instance of pharmaceutical drug
    ?drug wdt:P129 ?protein .            # physically interacts with

    # Protein is encoded by gene
    ?protein wdt:P31 wd:Q8054 .          # instance of protein
    ?protein wdt:P702 ?gene .            # encoded by

    # Gene is human
    ?gene wdt:P703 wd:Q15978631 .        # found in taxon: Homo sapiens

    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

---

## 5. Data Access Methods Summary

### 5.1 Wikidata Query Service (SPARQL)

| Feature | Details |
|---------|---------|
| **Endpoint** | https://query.wikidata.org/sparql |
| **Web UI** | https://query.wikidata.org |
| **Timeout** | 60 seconds |
| **Rate Limit** | Reasonable use policy (no hard limit) |
| **Formats** | JSON, XML, CSV, TSV |

### 5.2 Wikidata API

| Feature | Details |
|---------|---------|
| **Endpoint** | https://www.wikidata.org/w/api.php |
| **Actions** | wbgetentities, wbsearchentities, wbgetclaims |
| **Rate Limit** | 200 requests/second for bots |
| **Authentication** | Not required for reading |

**Example - Get Entity:**
```bash
curl "https://www.wikidata.org/w/api.php?action=wbgetentities&ids=Q188724&format=json"
```

### 5.3 Wikipedia API

| Feature | Details |
|---------|---------|
| **Endpoint** | https://en.wikipedia.org/w/api.php |
| **Actions** | parse, query, opensearch |
| **Infobox Access** | action=parse&prop=wikitext |

### 5.4 DBpedia SPARQL Endpoint

| Feature | Details |
|---------|---------|
| **Main Endpoint** | http://dbpedia.org/sparql |
| **Live Endpoint** | http://live.dbpedia.org/sparql |
| **Ontology** | http://dbpedia.org/ontology/ |
| **Resources** | http://dbpedia.org/resource/ |

### 5.5 Bulk JSON Dumps

| Resource | URL | Size | Update Frequency |
|----------|-----|------|------------------|
| Wikidata Full | dumps.wikimedia.org/wikidatawiki/entities/ | ~130 GB | Weekly |
| Wikidata Incremental | dumps.wikimedia.org/other/incr/wikidatawiki/ | Variable | Daily |
| Wikipedia | dumps.wikimedia.org/enwiki/ | ~20 GB | Monthly |
| DBpedia | downloads.dbpedia.org | ~10 GB | Quarterly |

---

## 6. Example SPARQL Queries

### 6.1 Get All Drugs with Gene Targets

```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>

SELECT ?drug ?drugLabel ?target ?targetLabel ?gene ?geneLabel ?geneSymbol
WHERE {
    # Drug is pharmaceutical drug
    ?drug wdt:P31 wd:Q12140 .

    # Drug interacts with target protein
    ?drug wdt:P129 ?target .

    # Target is a protein
    ?target wdt:P31 wd:Q8054 .

    # Protein is encoded by gene
    ?target wdt:P702 ?gene .

    # Get HGNC symbol
    OPTIONAL { ?gene wdt:P353 ?geneSymbol }

    # Human genes only
    ?gene wdt:P703 wd:Q15978631 .

    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
ORDER BY ?drugLabel
LIMIT 500
```

### 6.2 Get Drugs by ATC Code

```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>

SELECT ?drug ?drugLabel ?atc ?indication ?indicationLabel
WHERE {
    # Get drugs with ATC codes starting with "N05" (psycholeptics)
    ?drug wdt:P267 ?atc .
    FILTER(STRSTARTS(?atc, "N05"))

    # Get indications
    OPTIONAL { ?drug wdt:P2175 ?indication }

    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
ORDER BY ?atc
LIMIT 200
```

### 6.3 Get Drug-Disease Relationships

```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>

SELECT ?drug ?drugLabel ?disease ?diseaseLabel ?drugBankId ?atcCode
WHERE {
    # Drug treats disease
    ?drug wdt:P31 wd:Q12140 .
    ?drug wdt:P2175 ?disease .

    # Disease is subclass of disease
    ?disease wdt:P31/wdt:P279* wd:Q12136 .

    # Get identifiers
    OPTIONAL { ?drug wdt:P715 ?drugBankId }
    OPTIONAL { ?drug wdt:P267 ?atcCode }

    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
ORDER BY ?diseaseLabel
LIMIT 500
```

### 6.4 Link Drugs to Pathways (Federated Query with WikiPathways)

```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT ?drug ?drugLabel ?chebi ?pathwayId ?pathwayTitle
WHERE {
    # Get drug with ChEBI ID from Wikidata
    ?drug wdt:P31 wd:Q12140 .
    ?drug wdt:P683 ?chebiId .
    BIND(IRI(CONCAT("http://identifiers.org/chebi/CHEBI:", ?chebiId)) AS ?chebi)

    # Federated query to WikiPathways
    SERVICE <https://sparql.wikipathways.org/sparql> {
        ?metabolite a wp:Metabolite ;
                    wp:bdbChEBI ?chebi ;
                    dcterms:isPartOf ?pathway .
        ?pathway a wp:Pathway ;
                 dcterms:identifier ?pathwayId ;
                 dcterms:title ?pathwayTitle .
    }

    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

### 6.5 Get Drugs with Multiple Database Identifiers

```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>

SELECT ?drug ?drugLabel ?cas ?drugbank ?pubchem ?chebi ?chembl ?unii
WHERE {
    ?drug wdt:P31 wd:Q12140 .

    # Various identifiers
    OPTIONAL { ?drug wdt:P231 ?cas }        # CAS number
    OPTIONAL { ?drug wdt:P715 ?drugbank }   # DrugBank
    OPTIONAL { ?drug wdt:P662 ?pubchem }    # PubChem CID
    OPTIONAL { ?drug wdt:P683 ?chebi }      # ChEBI
    OPTIONAL { ?drug wdt:P592 ?chembl }     # ChEMBL
    OPTIONAL { ?drug wdt:P652 ?unii }       # FDA UNII

    # Must have at least DrugBank ID
    FILTER(BOUND(?drugbank))

    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
ORDER BY ?drugLabel
LIMIT 1000
```

### 6.6 Find Drug-Drug Interactions

```sparql
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>

SELECT ?drug1 ?drug1Label ?drug2 ?drug2Label
WHERE {
    ?drug1 wdt:P31 wd:Q12140 .
    ?drug1 wdt:P769 ?drug2 .  # significant drug interaction

    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 200
```

### 6.7 WikiPathways: Find All Pathways for a Gene

```sparql
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT DISTINCT ?pathway ?pathwayTitle ?organism ?geneLabel
WHERE {
    ?gene a wp:GeneProduct ;
          rdfs:label ?geneLabel ;
          dcterms:isPartOf ?pathway .
    ?pathway a wp:Pathway ;
             dcterms:title ?pathwayTitle ;
             wp:organismName ?organism .

    FILTER regex(?geneLabel, "^CYP3A4$", "i")
    FILTER (?organism = "Homo sapiens")
}
ORDER BY ?pathwayTitle
```

### 6.8 WikiPathways: Get All Drug Metabolites in Pathways

```sparql
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX dc: <http://purl.org/dc/elements/1.1/>

SELECT DISTINCT ?metabolite ?label ?chebi ?pathway ?pathwayTitle
WHERE {
    ?metabolite a wp:Metabolite ;
                rdfs:label ?label ;
                dcterms:isPartOf ?pathway .
    ?pathway a wp:Pathway ;
             dcterms:title ?pathwayTitle ;
             wp:organismName "Homo sapiens" .

    OPTIONAL { ?metabolite wp:bdbChEBI ?chebi }

    # Filter for drug-related pathways
    FILTER regex(?pathwayTitle, "drug|pharmacokinetics|metabolism", "i")
}
ORDER BY ?pathwayTitle ?label
LIMIT 500
```

---

## 7. Licensing Information

### 7.1 Wikidata - CC0 (Public Domain)

**License:** Creative Commons Zero (CC0 1.0)

**What this means:**
- No copyright restrictions
- No attribution required
- Can be used for any purpose including commercial
- Can be modified without restriction
- Maximum interoperability with other data sources

**Source:** [Wikidata:Licensing](https://www.wikidata.org/wiki/Wikidata:Licensing)

### 7.2 Wikipedia - CC BY-SA

**License:** Creative Commons Attribution-ShareAlike 4.0 (CC BY-SA 4.0)

**Requirements:**
- **Attribution:** Must credit Wikipedia and link to license
- **ShareAlike:** Derivatives must use same or compatible license
- Cannot add additional restrictions

**Source:** [Wikipedia:Reusing Wikipedia content](https://en.wikipedia.org/wiki/Wikipedia:Reusing_Wikipedia_content)

### 7.3 DBpedia - CC BY-SA

**License:** Creative Commons Attribution-ShareAlike 3.0 (CC BY-SA 3.0)

**Requirements:** Same as Wikipedia (derived from Wikipedia content)

### 7.4 WikiPathways - CC0 (Public Domain)

**License:** Creative Commons Zero (CC0 1.0)

**What this means:**
- Same freedoms as Wikidata
- No restrictions on reuse
- Ideal for integration with other databases

---

## 8. Integration Recommendations

### 8.1 For Drug Data Pipeline

1. **Primary Source:** Use Wikidata SPARQL for structured drug-target relationships
2. **Enrichment:** Supplement with Wikipedia infoboxes via wptools for additional details
3. **Pathways:** Query WikiPathways SPARQL for drug-pathway connections
4. **Bulk Processing:** Download Wikidata JSON dumps for large-scale analysis

### 8.2 Identifier Mapping Strategy

```
Wikidata Q-ID <-> DrugBank ID (P715)
                 |
                 +-> PubChem CID (P662)
                 |
                 +-> ChEBI ID (P683)
                 |
                 +-> ChEMBL ID (P592)
                 |
                 +-> CAS Number (P231)
                 |
                 +-> ATC Code (P267)
```

### 8.3 Data Freshness

| Source | Update Frequency | Lag |
|--------|------------------|-----|
| Wikidata SPARQL | Real-time | Minutes |
| Wikidata Dumps | Weekly | 1 week |
| Wikipedia API | Real-time | Minutes |
| DBpedia Main | Quarterly | 1-3 months |
| DBpedia Live | Real-time | Minutes |
| WikiPathways | Monthly | 1 month |

---

## 9. References

### Documentation
- [Wikidata:SPARQL query service](https://www.wikidata.org/wiki/Wikidata:SPARQL_query_service)
- [Wikidata:Data access](https://www.wikidata.org/wiki/Wikidata:Data_access)
- [Wikidata:Database download](https://www.wikidata.org/wiki/Wikidata:Database_download)
- [Template:Infobox drug](https://en.wikipedia.org/wiki/Template:Infobox_drug)
- [WikiPathways SPARQL queries](https://www.wikipathways.org/sparql.html)
- [WikiPathways RDF portal](https://www.wikipathways.org/rdf.html)
- [DBpedia SPARQL](https://www.dbpedia.org/resources/sparql/)

### Publications
- Waagmeester A, et al. (2016). Wikidata as a semantic framework for the Gene Wiki initiative. *Database*. doi:10.1093/database/baw015
- Slenter DN, et al. (2018). WikiPathways: a multifaceted pathway database bridging metabolomics to other omics research. *Nucleic Acids Research*. 46(D1):D661-D667

### Code Examples
- [sebotic/SPARQL - Drug query examples](https://github.com/sebotic/SPARQL/blob/master/drugs_sparql_examples.md)
- [wptools - Wikipedia parsing library](https://pypi.org/project/wptools/)
- [rWikiPathways - R package](https://bioconductor.org/packages/rWikiPathways/)
- [pywikipathways - Python package](https://pypi.org/project/pywikipathways/)
