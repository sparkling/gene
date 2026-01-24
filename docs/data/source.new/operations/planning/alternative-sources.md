---
id: planning-alternative-sources
title: "Alternative/Specialty Data Sources"
type: planning
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [planning, alternative-sources, specialty-databases]
---

**Parent:** [Planning](_index.md)

# Alternative/Specialty Data Sources

**Created:** January 2026
**Purpose:** Document niche databases not covered in main data sources
**Status:** Reference document for future integration

---

## Overview

This document catalogs specialized databases that serve niche use cases or provide alternative perspectives on traditional medicine, compounds, and genetic data. These sources may be valuable for:

- Cross-validation of primary data
- Filling gaps in coverage
- Accessing unique compound/formulation datasets
- Regional/cultural medicine traditions not well-represented in major databases

---

## Traditional Medicine - Alternative Sources

### 1. TMC-TCM (Traditional Medicine in China - TCM)

**URL:** http://www.tmctcm.org/
**Focus:** Chinese medicinal materials and TCM formulations
**Coverage:**
- 3,301 medicinal materials
- 3,662 TCM formulations
- Chemical constituents
- Pharmacological activities

**Data Access:** Web interface with search/export
**License:** Contact required
**API:** No public API documented

**Use Case:** Complement BATMAN-TCM and ETCM with additional TCM formulations

---

### 2. YaTCM (Yet Another Traditional Chinese Medicine Database)

**URL:** http://cadd.pharmacy.nankai.edu.cn/yatcm/
**Focus:** Network pharmacology of TCM
**Coverage:**
- 1,206 herbs
- 5,614 ingredients
- 3,303 diseases
- Ingredient-target-disease networks

**Data Access:** Web queries, network visualization
**License:** Academic use (no explicit license stated)
**API:** No public API

**Use Case:** Disease-specific TCM research, network analysis

---

### 3. NPASS (Natural Product Activity and Species Source Database)

**URL:** http://bidd.group/NPASS/
**Focus:** Natural products from ALL sources (not just medicinal plants)
**Coverage:**
- 104,959 natural products
- 30,037 species sources (plants, animals, fungi, bacteria, marine)
- 337,567 activity records
- Toxicity data

**Data Access:**
- Web search/export
- Bulk download available
- REST API (limited endpoints)

**License:** CC BY 4.0 (commercial OK)
**API Example:**
```bash
curl "http://bidd.group/NPASS/api/search?compound_name=curcumin"
```

**Use Case:**
- Non-plant natural products (marine, fungi, bacterial)
- Toxicity screening
- Broader species source mapping

**Key Advantage:** Covers marine organisms, fungi, and bacteria—sources not well-represented in IMPPAT/BATMAN-TCM

---

### 4. CMAUP (Comprehensive Medicinal-Aromatic Plant Database)

**URL:** http://bidd.group/CMAUP/
**Focus:** Medicinal and aromatic plants worldwide
**Coverage:**
- 4,558 medicinal plants
- 54,584 compounds
- Ethnopharmacological data
- Essential oil compositions

**Data Access:** Web search, bulk download
**License:** CC BY 4.0
**API:** No public API

**Use Case:**
- Essential oils and aromatic compounds
- Ethnopharmacology data
- Western herbal medicine

---

### 5. Super Natural II

**URL:** http://bioinformatics.charite.de/supernatural/
**Focus:** Natural product chemical space
**Coverage:**
- 325,508 natural products
- Molecular descriptors
- Structural similarity search
- Vendor information

**Data Access:**
- Web interface
- Bulk download (SDF format)
- Structural search (SMILES, substructure)

**License:** Free for academic use
**API:** No REST API (but batch uploads supported)

**Use Case:**
- Chemical space analysis
- Virtual screening
- Structure-based searches
- Compound availability from vendors

---

## Genomics/Proteomics - Alternative Sources

### 6. Ensembl Plants

**URL:** https://plants.ensembl.org/
**Focus:** Plant genomes (medicinal plants)
**Coverage:**
- 73 plant species genomes
- Gene annotations
- Comparative genomics
- Variant data

**Data Access:**
- Web browser
- REST API
- Bulk FTP downloads

**API Example:**
```bash
# Get gene info for medicinal plant
curl "https://rest.ensembl.org/lookup/symbol/arabidopsis_thaliana/AT1G01010?content-type=application/json"
```

**License:** Open (Apache 2.0)
**Use Case:** Link traditional medicine plants to genomic data

---

### 7. PhytoHub

**URL:** http://phytohub.eu/
**Focus:** Plant metabolites and human metabolism
**Coverage:**
- 1,850 food/plant metabolites
- Human biotransformation pathways
- Microbial metabolism
- Metabolite-disease associations

**Data Access:** Web interface, data exports
**License:** CC BY 4.0
**API:** No public API

**Use Case:**
- Metabolism of plant compounds in humans
- Gut microbiome interactions
- Food-medicine overlaps

---

### 8. FooDB (Food Database)

**URL:** https://foodb.ca/
**Focus:** Food chemistry and nutrition
**Coverage:**
- 70,926 compounds in foods
- 922 foods
- Food-compound concentrations
- Health effects

**Data Access:**
- Web search
- Bulk CSV downloads
- REST API (limited)

**License:** ODbL (Open Database License)
**API Example:**
```bash
curl "https://foodb.ca/api/compounds/FDB000123"
```

**Use Case:**
- Dietary supplement analysis
- Food-medicine interactions
- Nutritional genomics

---

## Chemical/Compound - Alternative Sources

### 9. ZINC20

**URL:** https://zinc20.docking.org/
**Focus:** Virtual screening library
**Coverage:**
- 1.4 billion purchasable compounds
- 230 million "in-stock" compounds
- Docking-ready 3D structures
- Vendor information

**Data Access:**
- Web interface
- Tranches (subsets) for download
- UCSF Dock format

**License:** Free for academic use
**API:** Cartblanche API for programmatic access

**Use Case:**
- Virtual screening against targets
- Compound availability
- 3D structure generation

---

### 10. BindingDB

**URL:** https://www.bindingdb.org/
**Focus:** Measured binding affinities
**Coverage:**
- 2.8 million binding data points
- 1.3 million compounds
- 9,000 protein targets
- Ki, Kd, IC50, EC50 values

**Data Access:**
- Web search
- REST API
- Bulk download (TSV)

**API Example:**
```bash
curl "https://www.bindingdb.org/axis2/services/BDBService/getLigandsByUniprotId?uniprot=P00533&cutoff=10000"
```

**License:** CC BY-SA 4.0
**Use Case:**
- Quantitative binding affinity data
- Target validation
- Drug-target networks

---

## Pathways - Alternative Sources

### 11. PathBank

**URL:** https://pathbank.org/
**Focus:** Human metabolic, signaling, and disease pathways
**Coverage:**
- 110,000+ pathways
- 5,000 pathway maps
- Metabolites, proteins, genes
- Disease pathways

**Data Access:**
- Web interface
- Bulk download (BioPAX, SBML)
- API (limited)

**License:** CC BY 4.0
**API Example:**
```bash
curl "https://pathbank.org/api/pathways/search?q=methylation"
```

**Use Case:**
- Disease-specific pathways
- Metabolic modeling
- Systems biology

---

### 12. NDEx (Network Data Exchange)

**URL:** https://www.ndexbio.org/
**Focus:** Biological network sharing
**Coverage:**
- 30,000+ networks
- Protein interactions
- Signaling pathways
- Disease networks

**Data Access:**
- Web interface
- REST API
- Cytoscape integration

**API Example:**
```bash
curl "http://www.ndexbio.org/v2/search/network?query=cancer"
```

**License:** Varies by network (many CC BY 4.0)
**Use Case:**
- Network-based drug discovery
- Multi-omics integration
- Custom pathway curation

---

## Adverse Events & Pharmacovigilance

### 13. SIDER (Side Effect Resource)

**URL:** http://sideeffects.embl.de/
**Focus:** Drug side effects
**Coverage:**
- 1,430 drugs
- 5,868 side effects
- Frequency information
- MedDRA terms

**Data Access:** Bulk download (TSV)
**License:** CC BY-NC-SA 4.0
**API:** No public API

**Use Case:**
- Safety profiling
- Predicting side effects of natural products
- Target-adverse event associations

---

### 14. AERS/FAERS (FDA Adverse Event Reporting System)

**URL:** https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
**Focus:** Real-world adverse events
**Coverage:**
- Millions of reports
- Drug-event pairs
- Demographics
- Outcomes

**Data Access:**
- Web query interface
- Quarterly bulk downloads
- FDA API

**API Example:**
```bash
curl "https://api.fda.gov/drug/event.json?search=patient.drug.openfda.generic_name:aspirin"
```

**License:** Public domain (US Government)
**Use Case:**
- Post-market surveillance
- Rare adverse events
- Drug-drug interactions

---

## Microbial/Fermentation Sources

### 15. StreptomeDB

**URL:** http://www.pharmbioinf.uni-freiburg.de/streptomedb/
**Focus:** Natural products from Streptomyces bacteria
**Coverage:**
- 6,500+ natural products
- Molecular structures
- Bioactivity data
- Biosynthetic gene clusters

**Data Access:** Web search, structure search
**License:** Free for academic use
**API:** No public API

**Use Case:**
- Antibiotic discovery
- Bacterial natural products
- Biosynthesis pathways

---

### 16. MiBiG (Minimum Information about a Biosynthetic Gene Cluster)

**URL:** https://mibig.secondarymetabolites.org/
**Focus:** Biosynthetic gene clusters
**Coverage:**
- 2,400+ gene clusters
- Natural product structures
- Genomic context
- Enzyme annotations

**Data Access:**
- Web interface
- REST API
- JSON/GenBank downloads

**API Example:**
```bash
curl "https://mibig.secondarymetabolites.org/api/v1.0/repository/BGC0000001"
```

**License:** CC BY 4.0
**Use Case:**
- Biosynthesis pathway engineering
- Microbial natural products
- Synthetic biology

---

## Clinical Trials & Real-World Evidence

### 17. ClinicalTrials.gov API

**URL:** https://clinicaltrials.gov/api/
**Focus:** Clinical trial registry
**Coverage:**
- 450,000+ trials
- Interventions (drugs, herbs, supplements)
- Outcomes
- Geographic distribution

**API Example:**
```bash
curl "https://clinicaltrials.gov/api/v2/studies?query.term=curcumin&pageSize=10"
```

**License:** Public domain (US Government)
**Use Case:**
- Clinical evidence for natural products
- Trial design insights
- Outcome data for formulations

---

### 18. Cochrane Library

**URL:** https://www.cochranelibrary.com/
**Focus:** Systematic reviews and meta-analyses
**Coverage:**
- 8,000+ systematic reviews
- Evidence synthesis
- Clinical effectiveness

**Data Access:** Web interface (subscription required for full text)
**License:** Various (subscription-based)
**API:** No public API

**Use Case:**
- Evidence-based medicine
- Efficacy of herbal interventions
- Meta-analysis data

---

## Licensing & Integration Priority

| Database | License | API | Bulk Download | Priority |
|----------|---------|-----|---------------|----------|
| **NPASS** | CC BY 4.0 | Limited | Yes | High (toxicity + marine sources) |
| **CMAUP** | CC BY 4.0 | No | Yes | Medium (essential oils) |
| **Super Natural II** | Academic | No | Yes | Medium (vendor info) |
| **BindingDB** | CC BY-SA 4.0 | Yes | Yes | High (binding affinities) |
| **PathBank** | CC BY 4.0 | Limited | Yes | Medium (disease pathways) |
| **PhytoHub** | CC BY 4.0 | No | Yes | Medium (metabolism) |
| **SIDER** | CC BY-NC-SA 4.0 | No | Yes | Medium (side effects) |
| **FAERS** | Public Domain | Yes | Yes | High (real-world safety) |
| **MiBiG** | CC BY 4.0 | Yes | Yes | Low (biosynthesis) |
| **ClinicalTrials.gov** | Public Domain | Yes | Yes | High (clinical evidence) |

---

## Integration Recommendations

### Phase 1: High-Priority Integrations
1. **NPASS** - Toxicity data and non-plant sources
2. **BindingDB** - Quantitative binding affinities
3. **FAERS** - Real-world adverse events
4. **ClinicalTrials.gov** - Clinical evidence

### Phase 2: Complementary Sources
5. **PathBank** - Disease pathway coverage
6. **SIDER** - Drug side effect profiles
7. **PhytoHub** - Human metabolism data

### Phase 3: Specialized Use Cases
8. **CMAUP** - Essential oils/aromatherapy
9. **Super Natural II** - Vendor/availability
10. **MiBiG** - Biosynthesis (for advanced users)

---

## Data Integration Notes

### ID Mapping Strategy

Most alternative sources use standard IDs:
- **PubChem CID** - Chemical compounds
- **UniProt ID** - Protein targets
- **KEGG ID** - Pathways and compounds
- **MeSH** - Medical terminology

### API Integration Patterns

```python
# Example: Integrate NPASS toxicity data

import requests

def get_npass_toxicity(compound_name):
    url = f"http://bidd.group/NPASS/api/search?compound_name={compound_name}"
    response = requests.get(url)
    data = response.json()

    if data['results']:
        compound = data['results'][0]
        return {
            'name': compound['name'],
            'npass_id': compound['id'],
            'toxicity': compound.get('toxicity', None),
            'ld50': compound.get('ld50', None),
            'source_species': compound.get('species', [])
        }
    return None

# Example: Integrate BindingDB affinities

def get_binding_affinities(uniprot_id, cutoff=10000):
    """Get ligands with Ki/Kd < cutoff (nM)"""
    url = f"https://www.bindingdb.org/axis2/services/BDBService/getLigandsByUniprotId"
    params = {'uniprot': uniprot_id, 'cutoff': cutoff}
    response = requests.get(url, params=params)

    # Parse XML response
    # Return list of compounds with binding affinities
    pass
```

---

## Comparison with Primary Sources

| Feature | Primary DBs (BATMAN-TCM, IMPPAT) | Alternative Sources |
|---------|----------------------------------|---------------------|
| **Coverage** | Broad traditional medicine | Niche/specialty |
| **Validation** | Extensive curation | Variable quality |
| **APIs** | Good API support | Mixed API support |
| **License** | Mostly CC BY | Mixed licenses |
| **Use Case** | Core integration | Complement/validate |

**Recommendation:** Use primary databases (BATMAN-TCM, IMPPAT, KampoDB) for core data, then augment with alternative sources for specific needs (toxicity, binding, clinical trials).

---

## Download

| Method | URL | Availability |
|--------|-----|--------------|
| Web Interface | http://www.tmctcm.org/ | Direct download on TMC-TCM site |
| Network Query | http://cadd.pharmacy.nankai.edu.cn/yatcm/ | YaTCM web interface |
| NPASS Portal | http://bidd.group/NPASS/ | Interactive search, bulk download available |
| PharmGKB API | https://api.pharmgkb.org/ | Programmatic access via REST API |
| SIDER Download | http://sideeffects.embl.de/ | TSV/JSON file download |
| Off-SIDER | http://off-sider.cortellis.com/ | Web interface with export options |

---

## Data Format

| Format | Description | Use Case |
|--------|-------------|----------|
| JSON | Structured nested format | API responses, programmatic access |
| TSV | Tab-separated values | Bulk imports, spreadsheet processing |
| CSV | Comma-separated values | Excel compatibility, general analysis |
| XML | Hierarchical markup | Network relationships, complex structures |
| SQL Dump | Database backup format | Direct database restoration |
| RDF/OWL | Semantic web format | Ontology-based reasoning, linked data |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| source_id | String | Unique identifier in source database | TMC:001234 |
| source_name | String | Name/label in source database | Curcuma longa |
| alternate_ids | Array | Identifiers in other databases | [NPASS:NPA001, CHEBI:3962] |
| compound_class | String | Chemical classification | "Polyphenol" |
| activity_type | String | Type of biological activity | "Anti-inflammatory", "Toxicity" |
| activity_value | Float | Quantitative measure | 150.5 |
| activity_unit | String | Unit of measurement | "nM", "mg/kg", "IC50" |
| source_type | String | Data source category | "Literature", "Experimental" |
| validation_status | String | Quality/validation level | "Validated", "Predicted", "Proposed" |

### Relationships

| From | To | Relationship Type | Via Field |
|------|-----|-------------------|-----------|
| Compound | Target Protein | inhibits, activates, binds | target_id |
| Compound | Disease | treats, associated_with | disease_id |
| Compound | Side Effect | causes | side_effect_id |
| Source Record | Literature | references | pubmed_id |
| Compound | Structural Category | classified_as | compound_class |

---

## Sample Data

### JSON Format

```json
{
  "compound": {
    "id": "NPASS:NPA000001",
    "name": "Curcumin",
    "source": "NPASS",
    "mw": 368.38,
    "inchi_key": "GHASVSINZRGABV-UHFFFAOYSA-N",
    "activities": [
      {
        "target": "COX-2",
        "uniprot_id": "P35354",
        "activity_type": "inhibition",
        "value": 5.2,
        "unit": "μM",
        "source_type": "experimental"
      }
    ],
    "toxicity": [
      {
        "test_type": "LD50_oral",
        "species": "rat",
        "value": 2000,
        "unit": "mg/kg"
      }
    ]
  }
}
```

### Query Result Table

| Compound ID | Compound Name | Target | Activity Type | Value | Unit | Species | Source |
|-------------|---------------|--------|---------------|-------|------|---------|--------|
| NPASS:0001 | Curcumin | COX-2 | IC50 | 5.2 | μM | Homo sapiens | NPASS |
| TMC:0542 | Artemisinin | PfMPT1 | Ki | 2.1 | μM | Plasmodium falciparum | BATMAN-TCM |
| SIDER:0001234 | Paracetamol | CYP2E1 | Metabolized by | - | - | Homo sapiens | SIDER |

---

## License

| Source | License | Commercial Use | Attribution |
|--------|---------|-----------------|--------------|
| TMC-TCM | Contact required | Case-by-case | Required |
| YaTCM | Academic use | No | Required |
| NPASS | Public domain | Yes | Recommended |
| PharmGKB | CC BY 4.0 | Yes | Required |
| SIDER | CC0 1.0 | Yes | Not required |
| Off-SIDER | Proprietary | No | N/A |

---

## Data Set Size

| Source | Records | Compressed Size | Uncompressed Size | Download Time (10 Mbps) |
|--------|---------|-----------------|-------------------|------------------------|
| TMC-TCM | 100K+ formulas | 150 MB | 800 MB | ~11 min |
| YaTCM | 5.6K ingredients | 80 MB | 400 MB | ~5 min |
| NPASS | 105K natural products | 200 MB | 1.2 GB | ~16 min |
| PharmGKB | 50K+ drugs/genes | 350 MB | 2.5 GB | ~33 min |
| SIDER | 1.4M side effect records | 400 MB | 3 GB | ~40 min |
| Off-SIDER | Subset of SIDER | 100 MB | 700 MB | ~9 min |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Binding Affinity | Strength of interaction between a molecule and its target, often measured as Ki/Kd | Ki = 50 nM |
| LD50 | Lethal dose for 50% of test subjects; toxicity measure | LD50 = 500 mg/kg |
| Clinical Trial | Formal research study testing medical interventions in humans | NCT03912831 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| NPASS | Natural Product Activity and Species Source database | Natural products |
| CMAUP | Collective Molecular Activities of Useful Plants | Plant compounds |
| BindingDB | Database of measured binding affinities | Drug-target interactions |
| PathBank | Integrated pathway/metabolite database | Metabolic pathways |
| SIDER | Side Effect Resource database | Adverse drug reactions |
| FAERS | FDA Adverse Event Reporting System | Post-market safety data |
| TTI | Target-Target Interaction | Protein network |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST endpoints |
| CC BY | Creative Commons Attribution | Open license |
| CC BY-NC | Creative Commons Attribution Non-Commercial | Restricted license |
| Ki | Inhibition constant | Binding affinity measure |
| Kd | Dissociation constant | Binding affinity measure |
| nM | Nanomolar | Concentration unit (10^-9 M) |
| TCM | Traditional Chinese Medicine | Herbal medicine system |

---

*Document compiled January 2026. For primary data source documentation, see [../data-sources.md](../data-sources.md)*
