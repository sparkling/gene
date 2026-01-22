---
id: schema-batman-tcm
title: "BATMAN-TCM 2.0 Database Schema"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./_index.md)

# BATMAN-TCM 2.0 Database Schema

**Document ID:** BATMAN-TCM-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** BATMAN-TCM 2.0 (2023)

---

## TL;DR

BATMAN-TCM 2.0 (Bioinformatics Analysis Tool for Molecular mechANism of TCM) is the most comprehensive Traditional Chinese Medicine database with API access. Contains 54,832 formulas, 8,404 herbs, 39,171 ingredients, 17,068 known TTIs (Target-TCM Interactions), and 2,319,272 predicted TTIs. Provides REST API returning JSON data with bidirectional query support (ingredient->target and target->ingredient).

---

## Database Statistics (Version 2.0)

| Entity | Count | Change from v1.0 |
|--------|-------|-----------------|
| **TCM Formulas** | 54,832 | +16.9% |
| **Herbs** | 8,404 | +3% |
| **Ingredients (compounds)** | 39,171 | +215.9% (3.16x) |
| **Known TTIs** | 17,068 | +62.3x |
| **Predicted TTIs** | 2,319,272 | +3.23x |
| **Target Proteins** | ~15,000+ | (estimated) |

### Prediction Quality

| Metric | Value |
|--------|-------|
| ROC AUC | 0.9663 |
| Method | Weighted voting ensemble |
| Sources | ChEMBL, CTD, STITCH, DGIdb |

---

## Data Access

### Primary URL
```
http://bionet.ncpsb.org.cn/batman-tcm/
```

### API Access
```
http://bionet.ncpsb.org.cn/batman-tcm/api
```

**Note**: Server may have connectivity issues (timeout observed). Plan for retries.

### Download Availability
- Tab-delimited text files available
- JSON via API queries
- Bulk download documented

---

## API Documentation

### Overview

The BATMAN-TCM 2.0 API allows programmatic access to TCM ingredient-target interaction data.

**Key Features:**
- URL-based parameter building
- JSON format output (computer-readable)
- Hypertext format output (browser viewing)
- Bidirectional query support

### Query Modes

#### Mode 1: Traditional Query (TCM -> Targets)

Query TCM ingredients, herbs, or formulas to find target proteins.

| Input Type | Description |
|------------|-------------|
| Formula | TCM prescription/formula name |
| Herb | Single herb name |
| Compound | Chemical ingredient |

**Output includes:**
- Target protein list
- Interaction confidence scores
- Evidence source (known vs predicted)
- KEGG/GO/disease enrichment

#### Mode 2: Reverse Query (Targets -> TCM)

Query disease-specific signature genes to find candidate TCM ingredients.

**Input:**
- Differentially expressed genes
- Disease signature genes
- Gene list from omics studies

**Output:**
- Candidate TCM ingredients
- Candidate herbs
- Enrichment analysis results

### API URL Structure

```
http://bionet.ncpsb.org.cn/batman-tcm/api?
  request_type=[type]
  &output_format=[json|html]
  &input_item=[query]
```

### Response Format

**JSON Response Structure:**
```json
{
  "query": "herb_name or compound_name",
  "query_type": "herb|compound|formula",
  "results": [
    {
      "target_id": "UniProt accession or gene symbol",
      "target_name": "protein name",
      "interaction_type": "known|predicted",
      "confidence_score": 0.95,
      "evidence_source": "ChEMBL|CTD|STITCH|DGIdb",
      "pubmed_ids": ["12345678"],
      "binding_data": {
        "activity_type": "IC50|Ki|EC50",
        "activity_value": 1.5,
        "activity_unit": "nM|uM"
      }
    }
  ],
  "enrichment": {
    "kegg_pathways": [...],
    "go_terms": [...],
    "diseases": [...]
  }
}
```

---

## Data Model

### TTI (Target-TCM Interaction) Schema

The core data structure linking TCM components to molecular targets.

```
┌─────────────────────────────────────────────────────────────────┐
│                       TCM FORMULAS                               │
│  54,832 traditional prescriptions from ancient medical texts     │
└─────────────────────┬───────────────────────────────────────────┘
                      │ contains (many-to-many)
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                         HERBS                                    │
│  8,404 medicinal herbs from Chinese Pharmacopoeia + literature  │
└─────────────────────┬───────────────────────────────────────────┘
                      │ contains (many-to-many)
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                      INGREDIENTS                                 │
│  39,171 chemical compounds (3.16x increase from v1.0)           │
└─────────────────────┬───────────────────────────────────────────┘
                      │ interacts_with (many-to-many)
                      │ Known: 17,068 | Predicted: 2,319,272
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                     TARGET PROTEINS                              │
│  ~15,000+ human proteins with functional annotations             │
│  Linked to: KEGG pathways, GO terms, OMIM diseases, TTD         │
└─────────────────────────────────────────────────────────────────┘
```

### Entity Schemas

#### Formula Entity
```json
{
  "formula_id": "string (internal ID)",
  "formula_name_chinese": "string (Chinese characters)",
  "formula_name_pinyin": "string (romanized)",
  "formula_name_english": "string",
  "source_text": "string (classical text reference)",
  "herb_composition": [
    {
      "herb_id": "string",
      "herb_name": "string",
      "dosage": "string",
      "processing": "string (preparation method)"
    }
  ],
  "therapeutic_category": "string",
  "indications": ["string array"]
}
```

#### Herb Entity
```json
{
  "herb_id": "string (internal ID)",
  "herb_name_chinese": "string",
  "herb_name_pinyin": "string",
  "herb_name_english": "string",
  "herb_name_latin": "string (botanical name)",
  "pharmacopoeia_status": "CP2020|not_listed",
  "properties_tcm": {
    "nature": "cold|cool|neutral|warm|hot",
    "flavor": ["sweet", "bitter", "sour", "salty", "pungent"],
    "meridian_tropism": ["lung", "heart", "liver", "spleen", "kidney"]
  },
  "ingredient_count": "integer"
}
```

#### Ingredient Entity
```json
{
  "ingredient_id": "string (internal ID)",
  "compound_name": "string",
  "compound_name_chinese": "string (if applicable)",
  "cas_number": "string",
  "pubchem_cid": "integer",
  "chembl_id": "string",
  "smiles": "string",
  "inchi": "string",
  "inchi_key": "string",
  "molecular_formula": "string",
  "molecular_weight": "float",
  "source_herbs": ["herb_id array"],
  "drugbank_id": "string (if available)"
}
```

#### TTI (Target-TCM Interaction) Entity
```json
{
  "tti_id": "string (unique interaction ID)",
  "ingredient_id": "string",
  "target_id": "string (UniProt or gene symbol)",
  "interaction_type": "known|predicted",
  "confidence_score": "float (0-1)",
  "evidence_sources": {
    "chembl": {
      "assay_count": "integer",
      "activity_data": [
        {
          "activity_type": "IC50|Ki|EC50|Kd",
          "activity_value": "float",
          "activity_unit": "nM|uM|mM",
          "assay_id": "string"
        }
      ]
    },
    "ctd": {
      "interaction_types": ["increases|decreases|affects"],
      "pubmed_refs": ["string array"]
    },
    "stitch": {
      "combined_score": "integer (0-1000)"
    },
    "dgidb": {
      "interaction_types": ["string array"],
      "source_dbs": ["string array"]
    },
    "literature": {
      "pubmed_ids": ["string array"],
      "extraction_confidence": "float"
    }
  },
  "prediction_method": "weighted_voting_ensemble (if predicted)"
}
```

#### Target Entity
```json
{
  "target_id": "string (UniProt accession)",
  "gene_symbol": "string",
  "gene_name": "string",
  "uniprot_id": "string",
  "entrez_gene_id": "integer",
  "protein_class": "string (ChEMBL classification)",
  "functional_annotations": {
    "kegg_pathways": [
      {
        "pathway_id": "hsa:XXXXX",
        "pathway_name": "string"
      }
    ],
    "go_biological_process": ["GO:XXXXXXX"],
    "go_molecular_function": ["GO:XXXXXXX"],
    "go_cellular_component": ["GO:XXXXXXX"]
  },
  "disease_associations": {
    "omim_ids": ["integer array"],
    "ttd_diseases": ["string array"],
    "kegg_diseases": ["string array"]
  }
}
```

---

## Known TTI Data Sources

BATMAN-TCM 2.0 curates known interactions from:

| Source Database | Content Type |
|-----------------|--------------|
| **ChEMBL** | Bioactivity data with quantitative binding |
| **DrugBank** | Drug-target interactions |
| **KEGG** | Compound-enzyme relationships |
| **TTD** (Therapeutic Target Database) | Drug targets |
| **HIT** (Herbal Ingredients' Targets) | Herbal compound targets |
| **HERB** | TCM ingredient-target data |

---

## Predicted TTI Methodology

### Prediction Pipeline

1. **Feature extraction**: Compound structure descriptors + target sequence features
2. **Ensemble methods**: Multiple classifiers combined
3. **Weighted voting**: Confidence scores from each method
4. **Threshold filtering**: High-confidence predictions retained

### Performance Metrics

| Metric | Value |
|--------|-------|
| ROC AUC | 0.9663 |
| Validation | 5-fold cross-validation |
| Test set | Independent hold-out |

---

## Enrichment Analysis Features

### KEGG Pathway Enrichment

**Enrichment Ratio Formula:**
```
ER = (nquery/Nquery) / (nall/Nall)

Where:
- Nquery = total genes in query set
- nquery = genes interacting with TCM in query set
- Nall = total genes in database
- nall = genes interacting with TCM in database
```

### Available Analyses

| Analysis Type | Database Source |
|---------------|-----------------|
| Pathway enrichment | KEGG |
| GO enrichment | Gene Ontology |
| Disease association | OMIM, TTD |
| Network visualization | Cytoscape.js |

---

## Cross-References

### Identifier Mappings

| BATMAN-TCM Field | External Database | Notes |
|------------------|-------------------|-------|
| ingredient_id | PubChem CID | Via structure matching |
| ingredient_id | ChEMBL ID | Via InChIKey |
| target_id | UniProt | Primary identifier |
| target_id | Entrez Gene | Secondary identifier |
| herb_name | Chinese Pharmacopoeia | CP2020 standard |

### External Database Links

| Database | URL | Data Type |
|----------|-----|-----------|
| PubChem | https://pubchem.ncbi.nlm.nih.gov/ | Compound structures |
| ChEMBL | https://www.ebi.ac.uk/chembl/ | Bioactivity data |
| UniProt | https://www.uniprot.org/ | Protein sequences |
| KEGG | https://www.kegg.jp/ | Pathways |
| DrugBank | https://go.drugbank.com/ | Drug-target |
| TTD | https://db.idrblab.net/ttd/ | Therapeutic targets |

---

## Download Files

### Available Formats

| File Type | Content | Format |
|-----------|---------|--------|
| Formula data | All formulas with compositions | TSV |
| Herb data | Herb metadata and properties | TSV |
| Ingredient data | Compound information | TSV |
| Known TTI | Curated interactions | TSV |
| Predicted TTI | High-confidence predictions | TSV |
| Target data | Protein annotations | TSV |

### Bulk Download Structure

```
batman_tcm_2.0/
├── formulas.tsv
├── herbs.tsv
├── ingredients.tsv
├── known_tti.tsv
├── predicted_tti.tsv
├── targets.tsv
└── enrichment_annotations/
    ├── kegg_pathways.tsv
    ├── go_terms.tsv
    └── disease_associations.tsv
```

---

## Usage Examples

### API Query Example (cURL)

```bash
# Query targets for a specific herb
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=herb&query=Ginseng&output=json"

# Query targets for a compound
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=compound&query=ginsenoside_Rg1&output=json"

# Reverse query: find herbs targeting a gene set
curl "http://bionet.ncpsb.org.cn/batman-tcm/api?query_type=genes&genes=TP53,EGFR,MTOR&output=json"
```

### Python Integration Example

```python
import requests
import time

BASE_URL = "http://bionet.ncpsb.org.cn/batman-tcm/api"

def query_herb_targets(herb_name, output_format="json"):
    """Query targets for a TCM herb."""
    params = {
        "query_type": "herb",
        "query": herb_name,
        "output": output_format
    }
    response = requests.get(BASE_URL, params=params, timeout=30)
    return response.json()

def query_compound_targets(compound_name):
    """Query targets for a compound."""
    params = {
        "query_type": "compound",
        "query": compound_name,
        "output": "json"
    }
    response = requests.get(BASE_URL, params=params, timeout=30)
    return response.json()

def reverse_query_genes(gene_list):
    """Find TCM ingredients targeting a gene set."""
    params = {
        "query_type": "genes",
        "genes": ",".join(gene_list),
        "output": "json"
    }
    response = requests.get(BASE_URL, params=params, timeout=30)
    return response.json()

# Example usage
# Note: Add retry logic due to potential server timeouts
def query_with_retry(query_func, *args, retries=3, delay=5):
    for attempt in range(retries):
        try:
            return query_func(*args)
        except (requests.Timeout, requests.ConnectionError):
            if attempt < retries - 1:
                time.sleep(delay)
            else:
                raise

# Query ginseng targets
ginseng_targets = query_with_retry(query_herb_targets, "Ginseng")
print(f"Ginseng targets {len(ginseng_targets['results'])} proteins")

# Reverse query for cancer-related genes
cancer_herbs = query_with_retry(reverse_query_genes, ["TP53", "EGFR", "MTOR", "KRAS"])
print(f"Found {len(cancer_herbs['results'])} candidate herbs")
```

---

## Integration Recommendations

### Priority Integration Steps

1. **Download bulk files** for offline processing
2. **Use API for specific queries** (interactive use)
3. **Map compound IDs** to PubChem for structure standardization
4. **Map target IDs** to UniProt for protein integration

### Complementary Databases

| Database | Purpose | Priority |
|----------|---------|----------|
| KampoDB | Kampo (Japanese) medicine | High |
| HERB 2.0 | Gene expression data | High |
| SymMap 2.0 | Symptom mapping | Medium |
| TCMSID | ADME properties | Medium |
| HIT 2.0 | Validated interactions | High |

### Network Pharmacology Workflow

```
1. Query BATMAN-TCM for herb/formula targets
2. Perform enrichment analysis (KEGG, GO)
3. Validate key targets against HIT 2.0
4. Cross-reference with KampoDB for Japanese overlap
5. Add expression data from HERB 2.0
6. Map symptoms via SymMap 2.0
```

---

## License

**CC BY-NC 4.0** (Non-commercial use)

- **Permits**: Academic research, non-commercial use
- **Requires**: Attribution to BATMAN-TCM 2.0
- **Restricts**: Commercial use without agreement
- **Contact**: For commercial licensing, contact National Center for Protein Sciences Beijing

---

## Limitations

1. **Server reliability**: Timeouts observed; implement retry logic
2. **API rate limits**: Not documented; use reasonable delays (3+ sec between queries)
3. **Prediction confidence**: High-confidence only; validate critical targets
4. **Herb naming**: Chinese names may need standardization
5. **Commercial restrictions**: CC BY-NC 4.0 license

---

## Technical Notes

### Server Connectivity

The BATMAN-TCM server may experience intermittent connectivity issues:
- Implement exponential backoff for retries
- Consider bulk download for large-scale analysis
- Cache API responses locally

### Data Freshness

- Version 2.0 released: 2023
- Known TTI: Curated through 2023
- Predicted TTI: Generated 2023

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `formula_id` | Unique identifier for TCM prescription/formula | Internal string ID |
| `herb_id` | Unique identifier for medicinal herb | Internal string ID |
| `ingredient_id` | Unique identifier for chemical compound | Internal string ID |
| `tti_id` | Unique identifier for Target-TCM Interaction | Compound-target pair |
| `confidence_score` | Interaction prediction confidence (0-1) | 0.95 |
| `interaction_type` | Classification as known or predicted | known/predicted |
| `pubchem_cid` | PubChem Compound identifier | 5280343 |
| `chembl_id` | ChEMBL compound identifier | CHEMBL25 |
| `meridian_tropism` | TCM concept of organ system affinity | lung, heart, liver |
| `enrichment_ratio` | Ratio of observed to expected pathway hits | ER formula |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| TTI | Target-TCM Interaction linking ingredient to protein target | Core data model |
| TCM Formula | Traditional prescription combining multiple herbs | 54,832 formulas |
| Chinese Pharmacopoeia | Official standard for Chinese medicinal materials (CP2020) | Herb validation |
| Nature (TCM) | Thermal property of herb (cold/cool/neutral/warm/hot) | Properties_tcm |
| Flavor (TCM) | Taste classification (sweet/bitter/sour/salty/pungent) | Properties_tcm |
| Meridian | TCM organ system targeted by herb | Tropism annotation |
| Network Pharmacology | Systems-level analysis of multi-target drug effects | Integration workflow |
| Weighted Voting Ensemble | Machine learning method for TTI prediction | ROC AUC: 0.9663 |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| BATMAN-TCM | Bioinformatics Analysis Tool for Molecular mechANism of TCM | Database name |
| TCM | Traditional Chinese Medicine | Medical tradition |
| TTI | Target-TCM Interaction | Core interaction type |
| CP2020 | Chinese Pharmacopoeia 2020 | Standard reference |
| CTD | Comparative Toxicogenomics Database | Known TTI source |
| DGIdb | Drug-Gene Interaction Database | Known TTI source |
| STITCH | Search Tool for Interactions of Chemicals | TTI source |
| ChEMBL | Chemical database with bioactivity | TTI source |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway enrichment |
| GO | Gene Ontology | Functional enrichment |
| TTD | Therapeutic Target Database | Target annotation |
| OMIM | Online Mendelian Inheritance in Man | Disease association |
| CC BY-NC | Creative Commons Attribution-NonCommercial | License type |
| ROC AUC | Receiver Operating Characteristic Area Under Curve | Prediction quality |

---

## References

1. Kong X, et al. (2024). BATMAN-TCM 2.0: an enhanced integrative database for known and predicted interactions between traditional Chinese medicine ingredients and target proteins. Nucleic Acids Research, 52(D1):D1110-D1117.

2. Liu Z, et al. (2016). BATMAN-TCM: a Bioinformatics Analysis Tool for Molecular mechANism of Traditional Chinese Medicine. Scientific Reports, 6:21146.

3. Official Website: http://bionet.ncpsb.org.cn/batman-tcm/

4. PubMed: https://pubmed.ncbi.nlm.nih.gov/37904598/

5. NGDC Database Commons: https://ngdc.cncb.ac.cn/databasecommons/database/id/9483

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation from publication and database analysis |
