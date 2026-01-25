---
id: schema-stitch
title: "STITCH Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, chemical-protein, drug-target, interactions, network, confidence-scores]
---

# STITCH Schema Documentation

**Document ID:** SCHEMA-STITCH
**Version:** 1.0
**Source Version:** STITCH 5.0 (2016)

---

## TL;DR

STITCH stores chemical-protein interactions with confidence scores from multiple evidence channels (experimental, database, text mining, prediction). Uses PubChem CID-based identifiers (CIDm/CIDs prefixes) for chemicals and STRING-compatible identifiers for proteins. Combined scores range 0-1000, with 700+ considered high confidence.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Chemicals | 500,000+ | STITCH 5.0 |
| Proteins | 9,643,763 | STRING v11 compatible |
| Organisms | 2,031 | Multi-species |
| Chemical-Protein Interactions | 1.6 billion | All species |
| Human Interactions | ~16 million | 9606 taxid |
| Drug-like Chemicals | 350,000+ | DrugBank overlap |
| Metabolites | 150,000+ | KEGG/HMDB overlap |

---

## Entity Relationship Overview

```
                    ┌─────────────────┐
                    │    Chemical     │
                    │  (CIDm/CIDs)    │
                    └────────┬────────┘
                             │
              ┌──────────────┼──────────────┐
              │              │              │
              ▼              ▼              ▼
    ┌─────────────┐  ┌─────────────┐  ┌─────────────┐
    │ Experimental│  │  Database   │  │ Text Mining │
    │   Score     │  │   Score     │  │   Score     │
    └──────┬──────┘  └──────┬──────┘  └──────┬──────┘
           │                │                │
           └────────────────┼────────────────┘
                            │
                            ▼
                   ┌────────────────┐
                   │ Combined Score │
                   │   (0-1000)     │
                   └───────┬────────┘
                           │
                           ▼
                   ┌────────────────┐
                   │    Protein     │
                   │ (taxid.ENSP)   │
                   └────────────────┘
                           │
                           ▼
                   ┌────────────────┐
                   │    Actions     │
                   │ (activation,   │
                   │  inhibition,   │
                   │  binding)      │
                   └────────────────┘
```

---

## Core Tables/Entities

### chemical_chemical.links

**Description:** Chemical-chemical similarity and interaction data.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| chemical1 | String | Yes | STITCH chemical ID (CIDm/CIDs format) |
| chemical2 | String | Yes | STITCH chemical ID (CIDm/CIDs format) |
| similarity | Integer | Yes | Tanimoto similarity (0-1000) |
| experimental | Integer | No | Experimental co-occurrence score |
| database | Integer | No | Pathway/database co-occurrence |
| textmining | Integer | No | Literature co-occurrence |
| combined_score | Integer | Yes | Probabilistic combination |

### protein_chemical.links

**Description:** Primary chemical-protein interaction table with evidence scores.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| chemical | String | Yes | STITCH chemical ID |
| protein | String | Yes | STRING protein ID (taxid.ENSP) |
| experimental | Integer | No | Binding assays, crystal structures (0-1000) |
| prediction | Integer | No | Binding site similarity, docking (0-1000) |
| database | Integer | No | DrugBank, KEGG, pathway DBs (0-1000) |
| textmining | Integer | No | Literature co-occurrence (0-1000) |
| combined_score | Integer | Yes | Probabilistic integration (0-1000) |

### protein_chemical.links.detailed

**Description:** Detailed breakdown including transfer scores.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| chemical | String | Yes | STITCH chemical ID |
| protein | String | Yes | STRING protein ID |
| experimental_direct | Integer | No | Direct experimental evidence |
| experimental_transferred | Integer | No | Evidence transferred from homologs |
| database_direct | Integer | No | Direct database annotation |
| database_transferred | Integer | No | Transferred annotation |
| textmining_direct | Integer | No | Direct text mining |
| textmining_transferred | Integer | No | Transferred text mining |
| prediction | Integer | No | Computational prediction |
| combined_score | Integer | Yes | Final integrated score |

### chemicals

**Description:** Chemical compound annotations and identifiers.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| chemical | String | Yes | STITCH chemical ID |
| name | String | Yes | Preferred chemical name |
| molecular_weight | Float | No | Molecular weight (Da) |
| SMILES_string | String | No | Canonical SMILES |
| InChIKey | String | No | InChI key identifier |

### chemical.sources

**Description:** Source attribution for chemical-protein interactions.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| chemical | String | Yes | STITCH chemical ID |
| protein | String | Yes | STRING protein ID |
| source | String | Yes | Evidence source database |
| evidence_type | String | Yes | Type of evidence |
| score | Integer | Yes | Source-specific score |

### actions

**Description:** Interaction mode/mechanism annotations.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| item_id_a | String | Yes | Chemical or protein ID |
| item_id_b | String | Yes | Protein ID |
| mode | String | Yes | Action type |
| action | String | No | Specific action |
| is_directional | Boolean | Yes | Has directionality |
| a_is_acting | Boolean | No | A acts on B |
| score | Integer | Yes | Confidence score |

---

## Identifier Formats

### Chemical Identifiers

| Prefix | Meaning | Example |
|--------|---------|---------|
| CIDm | Merged stereo | CIDm00002244 |
| CIDs | Stereospecific | CIDs00002244 |
| CID0 | PubChem exact | CID000002244 |

**Conversion to PubChem CID:**
- Remove prefix (CIDm, CIDs, CID0)
- Remove leading zeros
- Example: CIDm00002244 -> PubChem CID 2244 (Aspirin)

### Protein Identifiers

| Format | Pattern | Example |
|--------|---------|---------|
| STRING ID | taxid.ENSP | 9606.ENSP00000269305 |
| Human | 9606.{ENSP} | 9606.ENSP00000269305 |
| Mouse | 10090.{ENSP} | 10090.ENSMUSP00000001 |

---

## Score Channels

### Evidence Types

| Channel | Source | Description |
|---------|--------|-------------|
| experimental | ChEMBL, BindingDB, PDB | Direct binding evidence |
| prediction | SMAP, binding site similarity | Computational prediction |
| database | DrugBank, KEGG, MATADOR | Curated annotations |
| textmining | PubMed abstracts/full text | Literature co-occurrence |

### Score Calculation

Scores are combined using noisy-OR probability:
```
combined = 1 - (1-exp)(1-pred)(1-db)(1-tm)
```

Where each score is first converted to probability (score/1000).

### Confidence Thresholds

| Combined Score | Confidence | Interpretation |
|---------------|------------|----------------|
| 900-1000 | Highest | Multiple strong evidence types |
| 700-899 | High | Strong single or multiple moderate |
| 400-699 | Medium | Moderate evidence |
| 150-399 | Low | Weak evidence, mostly predicted |
| 0-149 | Lowest | Very weak, likely false positive |

---

## Action Types

| Action | Description | Example |
|--------|-------------|---------|
| activation | Chemical activates protein | Agonist |
| inhibition | Chemical inhibits protein | Antagonist, inhibitor |
| binding | Physical binding | Non-functional binding |
| catalysis | Enzyme-substrate relation | Metabolic substrate |
| reaction | Metabolic reaction | Product/substrate |
| expression | Affects protein expression | Inducer/repressor |
| ptmod | Post-translational modification | Phosphorylation |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| /api/json/interactions | GET | Get interactions for identifier |
| /api/json/resolve | GET | Map name to STITCH ID |
| /api/json/network | GET | Get mixed chem-protein network |
| /api/json/enrichment | GET | Functional enrichment |
| /api/image/network | GET | Network visualization image |

### API Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| identifier | String | Chemical name or STITCH ID |
| identifiers | String | Multiple IDs (newline-separated) |
| species | Integer | NCBI taxonomy ID (default: 9606) |
| limit | Integer | Max results (default: 10) |
| required_score | Integer | Min combined score (default: 400) |
| network_type | String | functional/physical |

---

## Data Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| TSV | .tsv | Tab-separated bulk data |
| JSON | .json | API responses |
| XML | .xml | Legacy API format |
| PNG/SVG | .png/.svg | Network images |

---

## Sample Record

### protein_chemical.links.v5.0.tsv

```tsv
chemical	protein	experimental	prediction	database	textmining	combined_score
CIDm00002244	9606.ENSP00000269305	900	0	700	650	976
CIDm00002244	9606.ENSP00000001008	0	450	600	800	923
CIDm00005090	9606.ENSP00000269305	800	0	900	700	987
```

### chemicals.v5.0.tsv

```tsv
chemical	name	molecular_weight	SMILES_string
CIDm00002244	aspirin	180.157	CC(=O)OC1=CC=CC=C1C(=O)O
CIDm00005090	ibuprofen	206.285	CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
CIDm00000000	caffeine	194.191	CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

### actions.v5.0.tsv

```tsv
item_id_a	item_id_b	mode	action	is_directional	a_is_acting	score
CIDm00002244	9606.ENSP00000269305	inhibition	inhibitor	t	t	900
CIDm00002244	9606.ENSP00000001008	binding		t	f	650
CIDm00005090	9606.ENSP00000269305	inhibition	inhibitor	t	t	800
```

### API Response (JSON)

```json
[
  {
    "chemicalId": "CIDm00002244",
    "stringId": "9606.ENSP00000269305",
    "preferredName_chemical": "aspirin",
    "preferredName_protein": "PTGS2",
    "ncbiTaxonId": 9606,
    "score": 0.976,
    "experimentalScore": 0.900,
    "predictionScore": 0.000,
    "databaseScore": 0.700,
    "textminingScore": 0.650,
    "actions": ["inhibition"]
  },
  {
    "chemicalId": "CIDm00002244",
    "stringId": "9606.ENSP00000001008",
    "preferredName_chemical": "aspirin",
    "preferredName_protein": "PTGS1",
    "ncbiTaxonId": 9606,
    "score": 0.923,
    "experimentalScore": 0.000,
    "predictionScore": 0.450,
    "databaseScore": 0.600,
    "textminingScore": 0.800,
    "actions": ["inhibition"]
  }
]
```

### Network API Response

```json
{
  "nodes": [
    {"id": "CIDm00002244", "name": "aspirin", "type": "chemical"},
    {"id": "9606.ENSP00000269305", "name": "PTGS2", "type": "protein"},
    {"id": "9606.ENSP00000001008", "name": "PTGS1", "type": "protein"}
  ],
  "edges": [
    {
      "source": "CIDm00002244",
      "target": "9606.ENSP00000269305",
      "score": 976,
      "label": "inhibition"
    },
    {
      "source": "CIDm00002244",
      "target": "9606.ENSP00000001008",
      "score": 923,
      "label": "inhibition"
    }
  ]
}
```

---

## Evidence Source Databases

| Source | Type | Content |
|--------|------|---------|
| ChEMBL | Experimental | Bioactivity assays |
| BindingDB | Experimental | Binding affinities |
| PDB | Experimental | Crystal structures |
| DrugBank | Database | Drug-target annotations |
| KEGG | Database | Metabolic pathways |
| MATADOR | Database | Drug-target relations |
| CTD | Database | Chemical-gene associations |
| SIDER | Database | Side effect targets |
| PubMed | Textmining | Literature co-occurrence |
| PMC | Textmining | Full-text mining |

---

## Glossary

| Term | Definition |
|------|------------|
| CIDm | Merged stereo chemical identifier (stereo-agnostic) |
| CIDs | Stereospecific chemical identifier |
| STRING ID | Protein identifier format: taxonomy.Ensembl_protein |
| Combined Score | Probabilistic integration of all evidence channels |
| Experimental Score | Evidence from binding assays and structures |
| Database Score | Evidence from curated pathway databases |
| Textmining Score | Evidence from literature co-occurrence |
| Prediction Score | Evidence from computational methods |
| Action | Mode of interaction (activation, inhibition, binding) |
| Transfer | Evidence from homologous proteins/chemicals |
| Noisy-OR | Probabilistic combination method |

---

## Cross-Reference Mappings

| External ID | STITCH Field | Mapping |
|-------------|--------------|---------|
| PubChem CID | chemical | CIDm + zero-padded |
| DrugBank | chemical | Via chemicals.aliases |
| ChEBI | chemical | Via chemicals.aliases |
| STRING | protein | Direct (same format) |
| UniProt | protein | Via STRING mapping |
| Ensembl Protein | protein | Embedded in STRING ID |
| Gene Symbol | protein | Via STRING aliases |

---

## Species Coverage

| Taxon ID | Species | Chemicals | Interactions |
|----------|---------|-----------|--------------|
| 9606 | Homo sapiens | 500,000+ | 16M+ |
| 10090 | Mus musculus | 400,000+ | 12M+ |
| 10116 | Rattus norvegicus | 350,000+ | 8M+ |
| 7955 | Danio rerio | 200,000+ | 5M+ |
| 6239 | C. elegans | 150,000+ | 3M+ |
| 7227 | D. melanogaster | 180,000+ | 4M+ |
| 4932 | S. cerevisiae | 120,000+ | 2M+ |

---

## References

1. http://stitch.embl.de/
2. Szklarczyk D et al. (2016) STITCH 5: augmenting protein-chemical interaction networks. Nucleic Acids Res. 44:D535-D539
3. Kuhn M et al. (2014) STITCH 4: integration of protein-chemical interactions with user data. Nucleic Acids Res. 42:D401-D407
4. Kuhn M et al. (2008) STITCH: interaction networks of chemicals and proteins. Nucleic Acids Res. 36:D684-D688
5. STRING v11: https://string-db.org/
