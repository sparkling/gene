# Natural Products Databases

**Document ID:** 43-52-NATURAL-PRODUCTS
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../index.md](./../index.md)

---

## TL;DR

Natural products databases provide structural, bioactivity, and source organism data for 750K+ unique compounds from plants, microbes, and marine organisms. Priority sources are COCONUT (695K structures, CC0), LOTUS (750K structure-organism pairs, CC0), and NPASS (204K compounds with quantitative activity data). Target prediction tools (SwissTargetPrediction, PharmMapper) and ADMET predictors (ProTox 3.0, pkCSM) complete the compound-to-mechanism pipeline.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary structure database | COCONUT 2.0 | Largest open (CC0), 695K structures, REST API | Jan 2026 |
| Structure-organism pairs | LOTUS/Wikidata | CC0 license, community-curated, SPARQL access | Jan 2026 |
| Quantitative bioactivity | NPASS | 1M+ activity records with IC50/Ki values | Jan 2026 |
| Microbial NPs | NPAtlas | Comprehensive microbial coverage, CC BY 4.0 | Jan 2026 |
| Target prediction | SwissTargetPrediction + PharmMapper | Complementary methods (2D/3D similarity + pharmacophore) | Jan 2026 |
| ADMET prediction | ProTox 3.0 + pkCSM | Free, comprehensive endpoints, no registration | Jan 2026 |
| Commercial databases | DNP excluded | Subscription-only, not suitable for open KB | Jan 2026 |

---

## Database Catalog

### Major Aggregated Natural Product Databases

#### COCONUT 2.0 (Collection of Open Natural Products)

| Field | Value |
|-------|-------|
| **URL** | https://coconut.naturalproducts.net/ |
| **Maintainer** | Steinbeck Lab, Friedrich Schiller University Jena, Germany |
| **Content** | Largest open-access aggregated NP database from 50+ sources |
| **Records** | 695,133 unique natural product structures |
| **License** | **CC0 (Creative Commons Zero)** - no attribution required |
| **API** | REST API with OpenAPI specification |
| **Documentation** | https://coconut.naturalproducts.net/api-documentation |
| **Update Frequency** | Regular updates with versioning (COCONUT 2.0: September 2024) |
| **Priority** | **Tier 1 (MVP)** |
| **Storage Estimate** | 32 GB (full SQL dump); 700 MB (SDF with metadata) |

**Download Options:**

| File Type | Size | Description |
|-----------|------|-------------|
| SDF 2D Lite | 287.6 MB | Basic 2D structures |
| SDF 2D Full | 691.7 MB | Full 2D with metadata |
| SDF 3D | 305.3 MB | RDKit-generated 3D coordinates |
| CSV Lite | 191 MB | Tabular format |
| CSV Full | 207.9 MB | Full metadata |
| SQL Dump | 31.91 GB | Complete PostgreSQL database |

**Key Fields:**
- coconut_id, SMILES, InChI, InChIKey
- molecular_formula, molecular_weight
- np_likeness_score, synthetic_accessibility_score, qed_score
- data_sources (array of contributing databases)
- stereochemistry: 82,220 without stereocenters, 539,350 with preserved stereochemistry

---

#### LOTUS (Natural Products Online)

| Field | Value |
|-------|-------|
| **URL** | https://lotus.naturalproducts.net/ and https://www.wikidata.org/ |
| **Maintainer** | Wikidata community / Originally Pierre-Marie Allard et al. |
| **Content** | Referenced structure-organism pairs from comprehensive literature curation |
| **Records** | 750,000+ structure-organism pairs |
| **License** | **CC0 (Creative Commons Zero)** - completely open |
| **API** | Wikidata SPARQL queries |
| **Documentation** | SPARQL endpoint: https://query.wikidata.org/ |
| **Update Frequency** | Community-driven continuous updates |
| **Priority** | **Tier 1 (MVP)** |
| **Storage Estimate** | ~500 MB (extracted subset) |

**SPARQL Access Example:**
```sparql
SELECT ?compound ?compoundLabel ?organism ?organismLabel
WHERE {
  ?compound wdt:P235 ?inchikey .
  ?compound p:P703 ?statement .
  ?statement ps:P703 ?organism .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 100
```

**Key Features:**
- First major NP database hosted on Wikidata
- Structure (chemical entity) linked to Organism (species) with literature reference
- Dual hosting: Wikidata (community curation) + lotus.naturalproducts.net (user-friendly search)
- Note: naturalproducts.net version being phased out in favor of Wikidata

---

#### NPASS (Natural Product Activity & Species Source)

| Field | Value |
|-------|-------|
| **URL** | https://bidd.group/NPASS/ |
| **Maintainer** | BIDD Group, National University of Singapore |
| **Content** | NPs with quantitative experimental activity data |
| **Records** | 204,023 NPs, 8,764 targets, 48,940 species, 1,048,756 activity records |
| **License** | Free for academic use |
| **API** | Web interface, download page |
| **Update Frequency** | Updated 2023 |
| **Priority** | **Tier 1 (MVP)** |
| **Storage Estimate** | ~2 GB |

**Quantitative Activity Types:**
- IC50, Ki, EC50, GI50, MIC values
- ~95,000 composition/concentration records
- ~66,600 NPs with estimated activity profiles

**Target Types:**
- Proteins and protein families
- Protein complexes
- Nucleic acids
- Subcellular targets
- Microbial/pathogenic organisms
- Cell lines, tissue targets

**ADMET Properties:** Computed for all NPs

---

#### NPAtlas (Natural Products Atlas)

| Field | Value |
|-------|-------|
| **URL** | https://www.npatlas.org/ |
| **Maintainer** | Linington Research Group, Simon Fraser University |
| **Content** | Microbially-derived natural products from peer-reviewed literature |
| **Records** | 36,545 compounds (database version 2024_09) |
| **License** | **CC BY 4.0** (Creative Commons Attribution 4.0) |
| **API** | RESTful API |
| **Documentation** | PostgreSQL with RDKit cartridge |
| **Update Frequency** | Regular releases with version tracking |
| **Priority** | **Tier 1 (MVP)** |
| **Storage Estimate** | ~500 MB |

**Organism Coverage:**
- Bacteria
- Fungi (including lichens, mushrooms)
- Cyanobacteria

**Exclusions:**
- Plant-derived compounds (unless also found in microbes)
- Marine invertebrates, diatoms
- Primary metabolites
- Patents and conference proceedings

**Features:**
- Molecular weight cutoff: 3000 Da
- SNAP-MS tool for MS/MS-based compound family annotation
- 590 structural corrections, 1,347 papers curated

---

#### SuperNatural 3.0

| Field | Value |
|-------|-------|
| **URL** | http://bioinf-applied.charite.de/supernatural_3/ |
| **Maintainer** | Charite - Universitatsmedizin Berlin |
| **Content** | Natural compounds with mechanism of action and pathway data |
| **Records** | 449,058 natural compounds |
| **License** | Free for academic use |
| **API** | Web interface (no login required) |
| **Update Frequency** | Static |
| **Priority** | **Tier 2** |
| **Storage Estimate** | ~1.5 GB |

**Seven Function Modules:**
1. Compound Search
2. Mechanism of Actions
3. Pathways
4. Target Library
5. Disease Indication
6. COVID-19 virtual screening
7. Organoleptic properties (taste prediction)

**Confidence Scoring:**
- Score 1.0: Has taxonomy info + linked to 3+ NP databases
- Score 0.5: No taxonomy but linked to 1+ NP database

**Disease-Specific Predictions:**
- Antiviral, antibacterial, antimalarial, anticancer
- CNS-targeted compounds

---

#### UNPD (Universal Natural Products Database)

| Field | Value |
|-------|-------|
| **URL** | Original site (pkuxxj.pku.edu.cn/UNPD) - **NO LONGER ACCESSIBLE** |
| **Maintainer** | Originally Peking University |
| **Content** | Natural product structures |
| **Records** | 229,358 structures |
| **License** | N/A (discontinued) |
| **Status** | **Discontinued** - data preserved in other databases |
| **Priority** | N/A (access via COCONUT) |

**Access UNPD Data Via:**
1. **COCONUT**: Integrated 156,865 UNPD entries
2. **ISDB (In-Silico MS/MS Database)**: https://oolonek.github.io/ISDB/
   - Contains 170,602 NPs from UNPD with predicted MS/MS spectra

---

### Traditional Medicine Integrated Databases

#### KNApSAcK Family

| Field | Value |
|-------|-------|
| **URL** | http://www.knapsackfamily.com/KNApSAcK/ |
| **Maintainer** | Nara Institute of Science and Technology (NAIST), Japan |
| **Content** | Metabolites, species, traditional medicines (KAMPO, JAMU) |
| **Records** | 50,048 metabolites, 20,741 species, 101,500 relationships |
| **License** | Academic use |
| **API** | Web interface, database download |
| **Update Frequency** | Periodic |
| **Priority** | **Tier 2** |
| **Storage Estimate** | ~800 MB |

**Database Family Components:**

| Database | Content |
|----------|---------|
| KNApSAcK Core | 101,500 species-metabolite relationships |
| KNApSAcK WorldMap | 41,548 geographic zone-plant pairs |
| KAMPO DB | 336 Japanese herbal formulae, 278 plants |
| JAMU DB | 5,310 Indonesian formulae, 550 plants |
| Biological Activity DB | 2,418 activities, 33,706 plant-activity relationships |

**Mass Spectrometry Features:**
- Search by accurate mass
- Multiple ionization modes
- Metabolomics research support

---

#### HERB 2.0 (High-throughput Experiment- and Reference-guided Database)

| Field | Value |
|-------|-------|
| **URL** | http://herb.ac.cn/ |
| **Maintainer** | Academic consortium (Chinese institutions) |
| **Content** | TCM herbs with transcriptomic validation |
| **Records** | 7,263 herbs, 49,258 ingredients, 12,933 targets, 28,212 diseases |
| **License** | Academic use |
| **API** | Web interface |
| **Update Frequency** | HERB 2.0 released 2024 |
| **Priority** | **Tier 2** |
| **Storage Estimate** | ~1 GB |

**Unique Features:**
- Re-analyzed 6,164 gene expression profiles from 1,037 high-throughput experiments
- Maps TCM herbs/ingredients to 2,837 modern drugs via pharmacotranscriptomics
- 8,558 clinical trials and 8,032 meta-analyses
- Knowledge graph representations

---

#### CMAUP (Collective Molecular Activities of Useful Plants)

| Field | Value |
|-------|-------|
| **URL** | https://bidd.group/CMAUP/ |
| **Maintainer** | BIDD Group, National University of Singapore |
| **Content** | Plants from 153 countries with targets and diseases |
| **Records** | 5,765 plants, 47,645 ingredients, 646 targets, 656 diseases |
| **License** | Free for academic use |
| **API** | Web interface, searchable by geography |
| **Update Frequency** | Updated 2024 |
| **Priority** | **Tier 2** |
| **Storage Estimate** | ~600 MB |

**Plant Categories:**
- 2,567 medicinal plants
- 170 food plants
- 1,567 edible plants
- 119 garden plants

**2024 Features:**
- Human transcriptomic changes overlapping with 1,152 targets
- 74 diseases from 20,027 patient samples
- Clinical information for 185 plants in 691 clinical trials

---

### Food and Dietary Compound Databases

#### FooDB

| Field | Value |
|-------|-------|
| **URL** | https://foodb.ca/ |
| **Maintainer** | The Metabolomics Innovation Centre (TMIC), University of Alberta |
| **Content** | World's largest food constituent database |
| **Records** | Comprehensive food-compound relationships |
| **License** | Free; citation required for publications |
| **API** | REST API (beta) |
| **Documentation** | https://foodb.ca/api_doc |
| **Download** | https://foodb.ca/downloads |
| **Update Frequency** | Periodic |
| **Priority** | **Tier 1 (MVP)** |
| **Storage Estimate** | ~2 GB |

**Content Categories:**
- Macronutrients, micronutrients
- Flavor, color, aroma compounds
- Texture compounds

---

#### Phenol-Explorer

| Field | Value |
|-------|-------|
| **URL** | http://phenol-explorer.eu/ |
| **Maintainer** | INRA (French National Institute for Agricultural Research) |
| **Content** | Polyphenol content in foods |
| **Records** | 35,000+ content values, 500 polyphenols, 400+ foods |
| **License** | Free academic access |
| **API** | Web interface, advanced search |
| **Update Frequency** | Version 3.6 (2015) latest |
| **Priority** | **Tier 2** |
| **Storage Estimate** | ~200 MB |

**Data Types:**
- Polyphenol content in raw foods (60,000+ original values from 1,300+ publications)
- 380 metabolites identified in biofluids
- Pharmacokinetic parameters
- Processing effects on 161 polyphenols

---

### Specialized Natural Product Databases

#### GNPS (Global Natural Products Social Molecular Networking)

| Field | Value |
|-------|-------|
| **URL** | https://gnps.ucsd.edu/ |
| **Maintainer** | Dorrestein Lab, UCSD |
| **Content** | MS/MS spectral data and molecular networks |
| **Records** | 1,800+ public datasets, 490,000+ MS files, 1.2 billion tandem mass spectra |
| **License** | Open access with sharing requirements |
| **API** | Web interface, data upload/download, API |
| **Data Format** | MS/MS spectra (MGF, mzML), molecular networks |
| **Update Frequency** | Continuous community contributions |
| **Priority** | **Tier 2** |
| **Storage Estimate** | Reference only (TB-scale full data) |

**Core Functions:**
1. Molecular networking
2. Spectral library search
3. Community-curated reference spectra
4. MASST search (search against all public datasets)
5. DEREPLICATOR for peptidic NPs

---

#### MIBiG (Minimum Information about a Biosynthetic Gene cluster)

| Field | Value |
|-------|-------|
| **URL** | https://mibig.secondarymetabolites.org/ |
| **Maintainer** | International consortium |
| **Content** | Experimentally validated biosynthetic gene clusters |
| **Records** | 2,000+ BGC entries |
| **License** | Open access |
| **API** | Web interface, bulk download, API |
| **Data Format** | JSON, standardized annotations |
| **Update Frequency** | Regular releases |
| **Priority** | **Tier 2** |
| **Storage Estimate** | ~500 MB |

**Data Fields:**
- Genomic locus coordinates
- Producing organism taxonomy
- Biosynthetic class
- Compound names and structures
- Gene functions
- Cross-links to NP Atlas, PubChem

---

#### NuBBEDB (Brazilian Biodiversity)

| Field | Value |
|-------|-------|
| **URL** | https://nubbe.iq.unesp.br/portal/nubbedb.html |
| **Maintainer** | NuBBE, UNESP, Brazil |
| **Content** | Brazilian biodiversity natural products |
| **Records** | 2,147 compounds (78% plant, 15% semi-synthetic, 5% microbial) |
| **License** | Free academic access |
| **API** | Web interface |
| **Update Frequency** | Static |
| **Priority** | **Tier 3** |
| **Storage Estimate** | ~100 MB |

**Unique Data:**
- NMR spectroscopic data
- Species geographic locations
- Toxicological information

---

#### NPCARE (Natural Products for Cancer Regulation)

| Field | Value |
|-------|-------|
| **URL** | http://silver.sejong.ac.kr/npcare |
| **Maintainer** | Sejong University, Korea |
| **Content** | Natural products with anticancer activity |
| **Records** | 6,578 compounds, 2,566 extracts, 1,952 species, 34 cancer types |
| **License** | Free academic access |
| **API** | Web interface |
| **Update Frequency** | Static |
| **Priority** | **Tier 3** |
| **Storage Estimate** | ~200 MB |

**Anticancer Criteria:**
- Growth inhibition of cancer cell lines
- Downregulation of oncogenes
- Upregulation of tumor suppressor genes

---

## Target Prediction Tools

### SwissTargetPrediction

| Field | Value |
|-------|-------|
| **URL** | http://www.swisstargetprediction.ch |
| **Maintainer** | SIB Swiss Institute of Bioinformatics |
| **Method** | 2D/3D molecular similarity |
| **Database** | 376,342 active compounds, 3,068 targets |
| **Species** | Human and other vertebrates |
| **Access** | Web interface (no API for batch processing) |
| **Output** | Ranked target predictions with probability scores |
| **Runtime** | 15-20 seconds per drug-like molecule |
| **Priority** | **Tier 1 (MVP)** |

**Performance:** >70% accuracy for at least one correct human target in top 15

**Input Methods:**
- SMILES string
- MarvinJS structure sketcher
- File upload

---

### SEA (Similarity Ensemble Approach)

| Field | Value |
|-------|-------|
| **URL** | https://sea.bkslab.org/ |
| **Maintainer** | Shoichet Lab, UCSF |
| **Method** | Set-wise chemical similarity among ligands |
| **Access** | Web interface |
| **Validation** | Only target prediction model with systematic experimental validation |
| **Priority** | **Tier 2** |

**Best Fingerprints (at significance 0.05):**

| Fingerprint | Precision |
|-------------|-----------|
| Atom pair | Highest |
| Topological | 83.7% |
| Morgan | Good |
| MACCS | Moderate |

---

### PharmMapper

| Field | Value |
|-------|-------|
| **URL** | http://lilab.ecust.edu.cn/pharmmapper/ |
| **Maintainer** | Li Lab, East China University of Science and Technology |
| **Method** | Reverse pharmacophore mapping |
| **Database** | 23,236 proteins (16,159 druggable + 51,431 ligandable models) |
| **Access** | Web interface (free, no login) |
| **Runtime** | ~1 hour for full database screen |
| **Input** | Mol2 file with 3D coordinates (required) |
| **Priority** | **Tier 2** |

**Applications:**
- Drug repurposing
- Side effect prediction
- Target identification for natural products

---

## ADMET and Toxicity Prediction Tools

### ProTox 3.0

| Field | Value |
|-------|-------|
| **URL** | https://tox.charite.de/ |
| **Maintainer** | Charite - Universitatsmedizin Berlin |
| **Models** | 61 toxicity prediction models |
| **Access** | Web interface (free, no login) |
| **Input** | SMILES or PubChem name |
| **Priority** | **Tier 1 (MVP)** |

**Toxicity Endpoints:**

| Category | Endpoints |
|----------|-----------|
| Oral Toxicity | 6 toxicity classes |
| Organ Toxicity | Hepato-, neuro-, respiratory, cardio-, nephrotoxicity |
| Toxicological | Mutagenicity, carcinogenicity, cytotoxicity, immunotoxicity |
| Pathways | 12 adverse outcome pathways (AOPs) |
| Targets | 15 toxicity targets |
| MIEs | 14 molecular initiating events |
| Metabolism | 6 molecular targets |

---

### pkCSM

| Field | Value |
|-------|-------|
| **URL** | https://biosig.lab.uq.edu.au/pkcsm/ |
| **Maintainer** | Biosig Lab, University of Queensland |
| **Method** | Graph-based molecular signatures |
| **Models** | 14 regression + 16 classification models |
| **Access** | Web interface (free, no data retention) |
| **Priority** | **Tier 1 (MVP)** |

**ADMET Categories:**
- **A**bsorption: Intestinal absorption, Caco-2 permeability, P-gp substrate
- **D**istribution: VDss, BBB penetration, CNS permeability
- **M**etabolism: CYP450 inhibition/substrate
- **E**xcretion: Total clearance, renal OCT2 substrate
- **T**oxicity: AMES test (83.8% accuracy), hERG inhibition, hepatotoxicity

---

### admetSAR 3.0

| Field | Value |
|-------|-------|
| **URL** | http://lmmd.ecust.edu.cn/admetsar3/ |
| **Maintainer** | East China University of Science and Technology |
| **Content** | 370,000+ experimental ADMET data, 104,652 compounds |
| **Models** | 119 prediction endpoints |
| **Access** | Web interface |
| **Priority** | **Tier 2** |

**Modules:**
1. **Search**: Experimental ADMET data + similarity search
2. **Prediction**: 119 endpoints including environmental/cosmetic risk
3. **ADMETopt**: Automatic ADMET optimization via transformation rules

---

### ADMETlab 3.0

| Field | Value |
|-------|-------|
| **URL** | https://admetmesh.scbdd.com/ |
| **Maintainer** | SCBDD (Computational Biology & Drug Design) |
| **Access** | Web interface, **API available** |
| **Features** | Decision support system for compound optimization |
| **Priority** | **Tier 2** |

---

## Bioactivity and Interaction Databases

### ChEMBL

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/chembl/ |
| **Maintainer** | EMBL-EBI |
| **Content** | Large-scale bioactivity database with NP flag |
| **Records** | ~64,000 molecules flagged as natural products (release 33) |
| **License** | **CC BY-SA 3.0** |
| **API** | REST API, full programmatic access |
| **Data Format** | SQL, SDF, TSV |
| **Priority** | **Tier 1 (MVP)** |
| **Storage Estimate** | ~50 GB (full database) |

**Natural Products Features:**
- NP flag for natural product molecules
- NP-likeness scores
- Mapping to COCONUT database
- Experimental bioactivity (Ki, Kd, IC50, EC50)

---

### BindingDB

| Field | Value |
|-------|-------|
| **URL** | https://www.bindingdb.org/ |
| **Maintainer** | BindingDB Team (UCSD) |
| **Content** | Experimental binding affinities |
| **Records** | 3.2M data points, 1.4M compounds, 11.4K targets |
| **License** | Free for academic use |
| **API** | RESTful API, bulk download, KNIME workflows |
| **Data Format** | TSV, SDF, FASTA |
| **Priority** | **Tier 1 (MVP)** |
| **Storage Estimate** | ~5 GB |

**Download Files:**
- BindingDBTargetSequences.fasta (7.23 MB)
- BindingDB_CID.txt - PubChem mapping (22.70 MB)
- BindingDB_UniProt.txt - UniProt mapping
- BindingDB_DrugBankID.txt - DrugBank mapping

---

### ChEBI (Chemical Entities of Biological Interest)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/chebi/ |
| **Maintainer** | EMBL-EBI |
| **Content** | Chemical entities with biological roles |
| **Records** | 195,000+ entries |
| **License** | **CC BY 4.0** (non-proprietary) |
| **API** | Web services + libChEBI API (Python, Java) |
| **Data Format** | Relational tables, flat files, SDF |
| **Priority** | **Tier 2** |
| **Storage Estimate** | ~3 GB |

---

### TTD (Therapeutic Target Database)

| Field | Value |
|-------|-------|
| **URL** | https://idrblab.net/ttd/ |
| **Maintainer** | IDRB, Zhejiang University |
| **Content** | Therapeutic targets and drugs including nature-derived |
| **Records** | 3,730 targets, 39,862 drugs |
| **License** | Free (no login required) |
| **API** | Web interface, download |
| **Priority** | **Tier 2** |
| **Storage Estimate** | ~1 GB |

**Natural Products Section:**
- Nature-derived drugs with species origins
- Cross-links to ICD disease codes

---

## Integration Recommendations

### Priority Tier Summary

**Tier 1 - Essential (MVP):**

| Database | Records | License | Use Case |
|----------|---------|---------|----------|
| COCONUT 2.0 | 695K structures | CC0 | Primary structure source |
| LOTUS | 750K structure-organism pairs | CC0 | Organism associations |
| NPASS | 204K NPs + 1M activities | Academic | Quantitative bioactivity |
| NPAtlas | 36K microbial NPs | CC BY 4.0 | Microbial coverage |
| ChEMBL | 64K NP-flagged | CC BY-SA 3.0 | Validated bioactivity |
| BindingDB | 1.4M compounds | Academic | Binding affinities |
| FooDB | Comprehensive | Citation required | Food compounds |
| SwissTargetPrediction | 3K targets | Free web | Target prediction |
| ProTox 3.0 | 61 models | Free web | Toxicity prediction |
| pkCSM | 30 models | Free web | ADMET prediction |

**Tier 2 - Important:**

| Database | Use Case |
|----------|----------|
| SuperNatural 3.0 | Mechanisms of action |
| KNApSAcK | Traditional medicine + MS |
| HERB 2.0 | TCM transcriptomics |
| Phenol-Explorer | Polyphenol content |
| GNPS | MS/MS spectral matching |
| MIBiG | Biosynthetic gene clusters |
| PharmMapper | 3D pharmacophore targets |
| admetSAR 3.0 | Extended ADMET |

**Tier 3 - Supplementary:**

| Database | Use Case |
|----------|----------|
| NuBBEDB | Brazilian biodiversity |
| NPCARE | Cancer-specific NPs |
| CMAUP | Geographic distribution |

---

### Target Prediction Pipeline

```
Input: Natural Product Structure (SMILES/InChI)
            |
            v
    +-------+-------+
    |               |
    v               v
SwissTarget    PharmMapper
(2D/3D sim)    (3D pharmacophore)
    |               |
    +-------+-------+
            |
            v
    Consensus Targets
            |
            v
    +-------+-------+-------+
    |       |       |       |
    v       v       v       v
  STRING  ChEMBL  BindingDB  TTD
  (PPI)   (activity) (Kd/Ki) (validated)
            |
            v
    Network Analysis (Cytoscape)
            |
            v
    Pathway Enrichment (Reactome, KEGG)
```

---

### Proposed Unified Schema

```json
{
  "compound": {
    "id": "string (internal)",
    "external_ids": {
      "coconut": "string",
      "lotus_wikidata": "string",
      "pubchem_cid": "integer",
      "chembl_id": "string",
      "inchikey": "string"
    },
    "structure": {
      "smiles": "string",
      "inchi": "string",
      "inchikey": "string",
      "mol_formula": "string",
      "mol_weight": "float"
    },
    "properties": {
      "np_likeness": "float",
      "drug_likeness": "float",
      "synthetic_accessibility": "float",
      "qed_score": "float"
    }
  },
  "source_organisms": [{
    "taxon_id": "integer (NCBI)",
    "species_name": "string",
    "common_names": ["string"],
    "part": "string (optional)",
    "geographic_region": "string (optional)"
  }],
  "bioactivity": [{
    "target_id": "string (UniProt)",
    "target_name": "string",
    "activity_type": "IC50|Ki|EC50|etc.",
    "activity_value": "float",
    "activity_unit": "nM|uM|etc.",
    "source": "string (database)",
    "reference": "string (PMID or DOI)"
  }],
  "predicted_targets": [{
    "target_id": "string",
    "target_name": "string",
    "prediction_method": "SwissTarget|SEA|PharmMapper",
    "probability": "float",
    "validation_status": "predicted|validated"
  }],
  "admet": {
    "absorption": {},
    "distribution": {},
    "metabolism": {},
    "excretion": {},
    "toxicity": {}
  }
}
```

---

### API Access Summary

| Database | API Type | Authentication | Batch Support |
|----------|----------|----------------|---------------|
| COCONUT | REST (OpenAPI) | None | Yes |
| LOTUS/Wikidata | SPARQL | None | Yes |
| ChEMBL | REST | None | Yes |
| NPAtlas | REST | None | Yes |
| BindingDB | REST | None | Yes |
| ChEBI | Web services + libChEBI | None | Yes |
| SwissTargetPrediction | Web only | None | No |
| PharmMapper | Web only | None | No |
| ADMETlab 3.0 | REST API | None | Yes |

---

### Licensing Summary

| License | Databases | Commercial Use |
|---------|-----------|----------------|
| **CC0** | COCONUT, LOTUS | Yes, unrestricted |
| **CC BY 4.0** | NPAtlas, ChEBI, STRING | Yes, with attribution |
| **CC BY-SA 3.0** | ChEMBL | Yes, share-alike |
| **Academic Only** | NPASS, SuperNatural, KNApSAcK | Check terms |
| **Commercial** | DNP (excluded) | Subscription required |

---

## Storage Estimates Summary

| Database | Size | Format |
|----------|------|--------|
| COCONUT 2.0 (full) | 32 GB | PostgreSQL |
| COCONUT 2.0 (SDF) | 700 MB | SDF |
| LOTUS subset | 500 MB | Extracted |
| NPASS | 2 GB | Download |
| NPAtlas | 500 MB | PostgreSQL |
| ChEMBL NP subset | 5 GB | Extracted |
| BindingDB | 5 GB | TSV/SDF |
| FooDB | 2 GB | Download |
| **Total Tier 1** | **~48 GB** | Mixed |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [pharmaceuticals.md](./pharmaceuticals.md) | Sister document for pharma interventions |
| [primary.md](./../pathways/primary.md) | Pathway databases for target enrichment |
| [tcm.md](./../traditional/tcm.md) | TCM-specific databases |
| [ayurveda.md](./../traditional/ayurveda.md) | Ayurveda-specific databases |

---

## References

### Primary Publications

1. **LOTUS**: Rutz A, et al. (2022) "The LOTUS initiative for open knowledge management in natural products research." eLife 11:e70780.

2. **COCONUT 2.0**: Venkata C, et al. (2024) "COCONUT 2.0: a comprehensive overhaul and curation of the collection of open natural products database." Nucleic Acids Res. gkae1063.

3. **NPASS**: Zeng X, et al. (2023) "NPASS database update 2023: quantitative natural product activity and species source database for biomedical research." Nucleic Acids Res. 51(D1):D621-D628.

4. **NPAtlas 3.0**: van Santen JA, et al. (2024) "Natural Products Atlas 3.0: extending the database of microbially derived natural products." Nucleic Acids Res. 53(D1):D691.

5. **SuperNatural 3.0**: Gallo K, et al. (2023) "SuperNatural 3.0-a database of natural products and natural product-based derivatives." Nucleic Acids Res. 51(D1):D654-D659.

6. **ChEMBL 2023**: Zdrazil B, et al. (2024) "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods." Nucleic Acids Res. 52(D1):D1180-D1192.

### Review Articles

- Sorokina M, Steinbeck C (2020) "Review on natural products databases: where to find data in 2020." J Cheminform. 12:20.

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial document with 20+ databases catalogued |
