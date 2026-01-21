# Database Schemas Index

**Document ID:** SCHEMAS-INDEX
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 3.0
**Parent:** [43-00-INDEX.md](../43-00-INDEX.md)

---

## Executive Summary

This directory contains **42 schema documentation files** covering databases across the THREE WORLDS framework (Modern Genetics, Traditional Medicine, Nutritional Science) plus cross-cutting resources. The research phase involved 8 parallel agents fetching actual API specifications, sample data, and statistics.

### Key Metrics

| Metric | Value |
|--------|-------|
| **Schema Documents** | 42 files |
| **Databases Documented** | 38 unique databases |
| **Total Records Accessible** | 1B+ variants, 3.6M compounds, 80K diseases |
| **CC0/Public Domain** | 12 databases (MVP-ready) |
| **CC BY/CC BY-SA** | 16 databases (commercial-friendly) |
| **Non-Commercial Only** | 5 databases |
| **Academic/Special License** | 5 databases |

---

## THREE WORLDS Coverage

### WORLD 1: Modern Genetics (11 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [DBSNP-SCHEMA.md](./DBSNP-SCHEMA.md) | dbSNP/ALFA R4 | 900M+ variants, 12 populations | Public Domain | Critical |
| [CLINVAR-SCHEMA.md](./CLINVAR-SCHEMA.md) | ClinVar | 2.5M+ clinical variants | Public Domain | Critical |
| [PHARMGKB-SCHEMA.md](./PHARMGKB-SCHEMA.md) | PharmGKB/ClinPGx | 33 CPIC guidelines | CC BY-SA 4.0 | Critical |
| [GWAS-CATALOG-SCHEMA.md](./GWAS-CATALOG-SCHEMA.md) | GWAS Catalog | 186K studies, 1M+ associations | CC0 | High |
| [GNOMAD-SCHEMA.md](./GNOMAD-SCHEMA.md) | gnomAD | 76K+ genomes, 730K exomes | Open Access | Critical |
| [ALPHAMISSENSE-SCHEMA.md](./ALPHAMISSENSE-SCHEMA.md) | AlphaMissense | 71M missense predictions | CC BY 4.0 | High |
| [DBNSFP-SCHEMA.md](./DBNSFP-SCHEMA.md) | dbNSFP | 84M+ variants, 46 prediction scores | Academic/Commercial | High |
| [DBVAR-SCHEMA.md](./DBVAR-SCHEMA.md) | dbVar | 6.5M+ structural variants | Public Domain | High |
| [SPLICEAI-SCHEMA.md](./SPLICEAI-SCHEMA.md) | SpliceAI | Pre-computed splice predictions | Illumina License | Medium |
| [ENCODE-SCHEMA.md](./ENCODE-SCHEMA.md) | ENCODE | 1M+ functional elements | Open Access | High |
| [GENE-ONTOLOGY-SCHEMA.md](./GENE-ONTOLOGY-SCHEMA.md) | Gene Ontology | 44K terms, 7M annotations | CC BY 4.0 | Critical |

### WORLD 2: Traditional Medicine (5 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [KAMPODB-SCHEMA.md](./KAMPODB-SCHEMA.md) | KampoDB (Kampo) | 481 formulas, 3K compounds, 63K targets | **CC BY-SA 4.0** | High |
| [DR-DUKES-SCHEMA.md](./DR-DUKES-SCHEMA.md) | Dr. Duke's (Western Herbal) | 50K entries, 1,900 activities | **CC0** | High |
| [BATMAN-TCM-SCHEMA.md](./BATMAN-TCM-SCHEMA.md) | BATMAN-TCM 2.0 (TCM) | 54K formulas, 39K compounds | CC BY-NC 4.0 | Medium |
| [IMPPAT-SCHEMA.md](./IMPPAT-SCHEMA.md) | IMPPAT 2.0 (Ayurveda) | 4K plants, 18K phytochemicals | CC BY-NC 4.0 | Medium |
| [LOTUS-SCHEMA.md](./LOTUS-SCHEMA.md) | LOTUS | 750K+ NP-organism pairs | **CC0** | High |

### WORLD 3: Nutritional Science (5 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [USDA-FOODDATA-CENTRAL.md](./USDA-FOODDATA-CENTRAL.md) | USDA FoodData Central | 20,900+ foods | **CC0** | Critical |
| [DSLD-NIH.md](./DSLD-NIH.md) | DSLD (NIH) | 200K+ supplement labels | Public Domain | High |
| [FOODB.md](./FOODB.md) | FooDB | 71K compounds, 778 foods | Free (citation) | High |
| [OPEN-FOOD-FACTS-SCHEMA.md](./OPEN-FOOD-FACTS-SCHEMA.md) | Open Food Facts | 3M+ food products | **ODbL** | High |
| [HMP-SCHEMA.md](./HMP-SCHEMA.md) | Human Microbiome Project | 5K+ metagenomes | CC BY 4.0 | Medium |

---

## Cross-Cutting Resources

### Natural Products & Interventions (5 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [CHEMBL-SCHEMA.md](./CHEMBL-SCHEMA.md) | ChEMBL 36 | 2.9M compounds, 24.3M activities | CC BY-SA 3.0 | Critical |
| [COCONUT-SCHEMA.md](./COCONUT-SCHEMA.md) | COCONUT | 716K natural products | **CC0** | Critical |
| [PUBCHEM-SCHEMA.md](./PUBCHEM-SCHEMA.md) | PubChem | 115M+ compounds | **Public Domain** | Critical |
| [CHEBI-SCHEMA.md](./CHEBI-SCHEMA.md) | ChEBI | 60K+ chemical entities | CC BY 4.0 | Critical |
| [DR-DUKES-SCHEMA.md](./DR-DUKES-SCHEMA.md) | Dr. Duke's | 50K ethnobotanical entries | **CC0** | High |

### Pathways & Networks (6 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [REACTOME-SCHEMA.md](./REACTOME-SCHEMA.md) | Reactome | 2,712 pathways, 11K proteins | CC BY 4.0 | Critical |
| [WIKIPATHWAYS-GPML-SCHEMA.md](./WIKIPATHWAYS-GPML-SCHEMA.md) | WikiPathways | 3,100+ pathways | CC BY 4.0 | High |
| [KGML-SCHEMA.md](./KGML-SCHEMA.md) | KEGG KGML | Reference format | Academic License | Reference |
| [DISGENET-SCHEMA.md](./DISGENET-SCHEMA.md) | DisGeNET | 628K gene-disease associations | CC BY-NC-SA 4.0 | High |
| [STRING-SCHEMA.md](./STRING-SCHEMA.md) | STRING | 68M+ protein interactions | CC BY 4.0 | Critical |
| [INTACT-SCHEMA.md](./INTACT-SCHEMA.md) | IntAct | 1.1M+ molecular interactions | CC BY 4.0 | High |

### Diseases & Phenotypes (3 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [MONDO-SCHEMA.md](./MONDO-SCHEMA.md) | MONDO | 26K diseases, 130K cross-refs | CC BY 4.0 | Critical |
| [HPO-SCHEMA.md](./HPO-SCHEMA.md) | HPO | 13K phenotypes, 156K annotations | Open Access | Critical |
| [ORPHANET-ORDO-SCHEMA.md](./ORPHANET-ORDO-SCHEMA.md) | Orphanet ORDO | 6.5K rare diseases | CC BY 4.0 | High |

### Health Domain Applications (2 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [CBIOPORTAL-SCHEMA.md](./CBIOPORTAL-SCHEMA.md) | cBioPortal | 423 cancer studies | Various | Medium |
| [GMREPO-SCHEMA.md](./GMREPO-SCHEMA.md) | GMrepo | 119K microbiome samples | Free academic | Medium |

### Integration & Cross-Reference (2 databases)

| Document | Database | Records | License | Priority |
|----------|----------|---------|---------|----------|
| [WIKIDATA-SCHEMA.md](./WIKIDATA-SCHEMA.md) | Wikidata | 59K genes, 200K diseases | **CC0** | Critical |
| [UNIPROT-IDMAPPING-SCHEMA.md](./UNIPROT-IDMAPPING-SCHEMA.md) | UniProt ID Mapping | **286 database cross-refs** | CC BY 4.0 | Critical |

### Analysis & Integration Documents

| Document | Description |
|----------|-------------|
| [UNIFIED-SCHEMA-ANALYSIS.md](./UNIFIED-SCHEMA-ANALYSIS.md) | Cross-database unified model with 5 core entities |
| [SAMPLE-DATA.md](./SAMPLE-DATA.md) | Actual API response samples from pathway databases |
| [INDEX.md](./INDEX.md) | Secondary index (pathway focus) |

---

## Prioritized Shortlist by License Tier

### Tier 1: CC0/Public Domain (Unrestricted - MVP Ready)

| Database | Category | Records | Use Case |
|----------|----------|---------|----------|
| **Wikidata** | Integration | 59K genes, 200K diseases | Knowledge graph backbone |
| **USDA FoodData** | Nutrition | 20,900 foods | Food composition |
| **dbSNP/ALFA** | Genetics | 900M variants | Population frequencies |
| **ClinVar** | Genetics | 2M+ variants | Clinical significance |
| **COCONUT** | Natural Products | 716K compounds | NP structures |
| **Dr. Duke's** | Traditional Medicine | 50K entries | Ethnobotanical context |
| **WikiPathways** | Pathways | 3,100 pathways | Community pathways |
| **PubChem** | Chemistry | 115M+ compounds | Universal compound reference |
| **dbVar** | Genetics | 6.5M+ SVs | Structural variants |
| **LOTUS** | Natural Products | 750K+ pairs | NP-organism links |
| **DSLD** | Nutrition | 200K labels | Supplement data |
| **GWAS Catalog** | Genetics | 186K studies | Trait associations |

### Tier 2: CC BY/CC BY-SA (Attribution Required - Commercial OK)

| Database | Category | Records | License |
|----------|----------|---------|---------|
| **UniProt ID Mapping** | Integration | 286 databases | CC BY 4.0 |
| **Reactome** | Pathways | 2,712 pathways | CC BY 4.0 |
| **MONDO** | Diseases | 26K diseases | CC BY 4.0 |
| **HPO** | Phenotypes | 13K terms | Open Access |
| **Orphanet** | Rare Diseases | 6.5K diseases | CC BY 4.0 |
| **ChEMBL** | Bioactivity | 24M activities | CC BY-SA 3.0 |
| **KampoDB** | Kampo Medicine | 481 formulas | **CC BY-SA 4.0** |
| **PharmGKB** | Pharmacogenomics | 33 guidelines | CC BY-SA 4.0 |
| **Gene Ontology** | Functional Annotation | 44K terms | CC BY 4.0 |
| **ChEBI** | Chemical Entities | 60K+ entities | CC BY 4.0 |
| **AlphaMissense** | Variant Pathogenicity | 71M predictions | CC BY 4.0 |
| **STRING** | Protein Networks | 68M+ interactions | CC BY 4.0 |
| **IntAct** | Molecular Interactions | 1.1M+ interactions | CC BY 4.0 |
| **gnomAD** | Population Genetics | 76K genomes | Open Access |
| **ENCODE** | Functional Elements | 1M+ elements | Open Access |
| **HMP** | Microbiome | 5K+ metagenomes | CC BY 4.0 |

### Tier 3: Non-Commercial Only (Academic/Research)

| Database | Category | Records | License |
|----------|----------|---------|---------|
| BATMAN-TCM 2.0 | TCM | 54K formulas | CC BY-NC 4.0 |
| IMPPAT 2.0 | Ayurveda | 4K plants | CC BY-NC 4.0 |
| DisGeNET | Gene-Disease | 628K associations | CC BY-NC-SA 4.0 |

### Tier 4: Special/Academic License

| Database | Category | Records | License |
|----------|----------|---------|---------|
| dbNSFP | Variant Annotation | 84M+ variants | Academic/Commercial |
| SpliceAI | Splice Prediction | Pre-computed | Illumina License |
| Open Food Facts | Food Products | 3M+ products | ODbL |

---

## Hub Identifier Strategy

Based on cross-database analysis, these identifiers should serve as canonical hubs:

| Entity Type | Primary Hub ID | Rationale | Coverage |
|-------------|----------------|-----------|----------|
| **Gene** | NCBI Gene ID (Entrez) | Universal, links to dbSNP/ClinVar/DisGeNET | 100% human genes |
| **Variant** | rsID + SPDI | rsID for SNVs, SPDI for normalization | 900M+ variants |
| **Compound** | **InChIKey** | Structure-based, 27-char standard | All chemistry DBs |
| **Disease** | **MONDO ID** | Best cross-ref (OMIM, Orphanet, DOID) | 130K mappings |
| **Pathway** | Reactome ID | Highest curation quality | 2.7K pathways |
| **Phenotype** | HPO ID | GA4GH Phenopackets standard | 13K terms |
| **Protein** | UniProt ID | Universal protein identifier | 286 DB mappings |
| **Interaction** | STRING/IntAct | Protein-protein interactions | 68M+ interactions |
| **Functional Term** | GO ID | Gene Ontology standard | 44K terms |

### Cross-Reference Chain

```
Gene (NCBI Gene ID)
  |
  +---> Variant (rsID) ---> Clinical (ClinVar VCV)
  |           |                    |
  |           +---> Pathogenicity (AlphaMissense, dbNSFP)
  |           |
  |           +---> Splicing (SpliceAI)
  |           |
  |           +---> Population (gnomAD)
  |
  +---> Protein (UniProt) -------+---> Disease (MONDO)
  |         |                              |
  |         +---> Pathway (Reactome)       +---> Phenotype (HPO)
  |         |
  |         +---> Function (GO)
  |         |
  |         +---> Interaction (STRING/IntAct)
  |
  +---> Compound (InChIKey)
            |
            +---> Bioactivity (ChEMBL)
            |
            +---> Structure (PubChem, ChEBI)
            |
            +---> Natural Source (COCONUT, LOTUS)
            |
            +---> Traditional Use (Dr. Duke's / KampoDB)
            |
            +---> Food Source (USDA, Open Food Facts, FooDB)
```

---

## Storage Estimates

### By Tier

| Tier | Content | Raw Size | Indexed |
|------|---------|----------|---------|
| **MVP Core** | 12 CC0 databases (CSV/JSON) | ~5 GB | ~12 GB |
| **Tier 1+2 Full** | All CC0/CC BY databases | ~150 GB | ~300 GB |
| **Complete KB** | All 38 databases | ~500 GB | ~1 TB |

### By Database (New Additions)

| Database | Format | Size | Notes |
|----------|--------|------|-------|
| gnomAD | VCF/TSV | 400+ GB | Full genomes/exomes |
| AlphaMissense | TSV | 5 GB | Pre-computed scores |
| dbNSFP | TSV | 90 GB | Comprehensive annotations |
| dbVar | VCF/TSV | 15 GB | Structural variants |
| STRING | TSV | 50 GB | Protein networks |
| IntAct | PSI-MI TAB | 5 GB | Interactions |
| PubChem | SDF/JSON | 100+ GB | Full compound set |
| ChEBI | OWL/SDF | 2 GB | Ontology + structures |
| Gene Ontology | OWL/OBO | 500 MB | Ontology |
| HMP | Various | 20 TB | Raw metagenomes |
| Open Food Facts | JSON/CSV | 10 GB | Product database |
| SpliceAI | TSV | 10 GB | Pre-computed scores |
| ENCODE | Various | 50+ GB | Functional data |
| LOTUS | CSV/JSON | 500 MB | NP-organism pairs |

### By Database (Existing)

| Database | Format | Size | Notes |
|----------|--------|------|-------|
| dbSNP/ALFA | VCF | 16 GB | Population frequencies |
| ClinVar | TSV/XML | 6 GB | Weekly updates |
| ChEMBL 36 | SQLite | 5.2 GB | Full bioactivity |
| COCONUT | PostgreSQL | 31.9 GB | Full dump |
| COCONUT | CSV | 207 MB | MVP format |
| Reactome | Neo4j | ~2 GB | Graph database |
| MONDO/HPO | OWL/OBO | ~300 MB | Ontologies |
| Wikidata | SPARQL | Query-based | No local storage needed |

---

## Implementation Roadmap (24 Weeks)

### Phase 1: Core Infrastructure (Weeks 1-4)

| Task | Databases | Deliverable |
|------|-----------|-------------|
| ID Resolution Service | UniProt ID Mapping, Wikidata | Cross-reference API |
| Gene/Variant Model | dbSNP, ClinVar | Variant entity store |
| Disease Ontology | MONDO, HPO | Disease/phenotype graph |
| Functional Annotation | Gene Ontology, ChEBI | GO/ChEBI integration |

### Phase 2: Clinical Layer (Weeks 5-8)

| Task | Databases | Deliverable |
|------|-----------|-------------|
| Clinical Annotations | ClinVar, PharmGKB | Clinical significance API |
| Gene-Disease Links | DisGeNET, MONDO | Association scoring |
| Population Frequencies | gnomAD, dbSNP | Frequency lookup |
| Pathogenicity Prediction | AlphaMissense, dbNSFP | Variant scoring |

### Phase 3: Variant Analysis (Weeks 9-12)

| Task | Databases | Deliverable |
|------|-----------|-------------|
| Structural Variants | dbVar | SV entity store |
| Splice Predictions | SpliceAI | Splice impact API |
| Functional Elements | ENCODE | Regulatory annotation |
| Comprehensive Scoring | dbNSFP integration | Multi-score API |

### Phase 4: Compounds & Targets (Weeks 13-16)

| Task | Databases | Deliverable |
|------|-----------|-------------|
| Compound Entity | ChEMBL, COCONUT, PubChem | Unified compound store |
| Chemical Ontology | ChEBI | Structure-function links |
| Bioactivity Data | ChEMBL | Target-activity API |
| Pathway Integration | Reactome, WikiPathways | Pathway visualization |
| Interaction Networks | STRING, IntAct | Network analysis |

### Phase 5: Traditional Medicine & Nutrition (Weeks 17-20)

| Task | Databases | Deliverable |
|------|-----------|-------------|
| Kampo Integration | KampoDB | Formula-compound links |
| Western Herbal | Dr. Duke's | Ethnobotanical context |
| NP-Organism Links | LOTUS | Natural source mapping |
| TCM (if NC OK) | BATMAN-TCM | TCM formulas |

### Phase 6: Nutrition & Completion (Weeks 21-24)

| Task | Databases | Deliverable |
|------|-----------|-------------|
| Food Composition | USDA FoodData, FooDB | Nutrient lookup |
| Food Products | Open Food Facts | Commercial food data |
| Supplements | DSLD | Supplement ingredients |
| Microbiome | HMP | Microbiome context |
| Full Integration | All databases | Unified API |

---

## Key Findings

### 1. Traditional Medicine Licensing

**KampoDB is the only traditional medicine database with commercial-friendly license (CC BY-SA 4.0)**. BATMAN-TCM and IMPPAT are NC-only, which may limit commercial use of TCM and Ayurveda data. **LOTUS (CC0)** provides excellent NP-organism mappings as a complement.

### 2. UniProt as Integration Hub

UniProt ID Mapping provides **286 database cross-references** - this is the critical resource for linking all biological entities. Essential for Phase 1.

### 3. InChIKey for Compounds

InChIKey serves as the universal compound identifier across ChEMBL, COCONUT, PubChem, ChEBI, and all chemistry databases. Must be the canonical compound ID.

### 4. MONDO for Disease Unification

MONDO provides **129,914 cross-references** to OMIM, Orphanet, DOID, MeSH - the best disease integration hub available.

### 5. Wikidata as CC0 Backbone

With 59K genes and 200K diseases under CC0, Wikidata can serve as the foundation for data that can't use NC-restricted sources.

### 6. Variant Pathogenicity Expansion

AlphaMissense (71M predictions), dbNSFP (46 scores), and SpliceAI provide comprehensive variant impact assessment, all with commercial-friendly or academic licenses.

### 7. Interaction Networks

STRING (68M+ protein interactions) and IntAct (1.1M+ curated interactions) provide complementary network data - STRING for coverage, IntAct for curation depth.

### 8. Gene Ontology Critical for Function

GO provides **44K terms** and **7M annotations** - essential for functional interpretation of genetic findings.

### 9. PubChem as Universal Chemistry Reference

With 115M+ compounds under Public Domain, PubChem is the most comprehensive chemistry reference, essential for compound normalization.

---

## Database Statistics Summary

### WORLD 1: Genetics (Expanded)

| Database | Variants/Records | Genes | Key Feature |
|----------|------------------|-------|-------------|
| dbSNP/ALFA | 900M+ | 30K | 12 population frequencies |
| ClinVar | 2.5M | - | Clinical significance |
| PharmGKB | 10K | 1,761 | CPIC guidelines |
| GWAS Catalog | 512K SNPs | - | 186K studies |
| gnomAD | 730K exomes | 19K | Population allele frequencies |
| AlphaMissense | 71M | 19K | Pathogenicity predictions |
| dbNSFP | 84M+ | - | 46 prediction scores |
| dbVar | 6.5M+ SVs | - | Structural variants |
| SpliceAI | Pre-computed | - | Splice site predictions |
| ENCODE | 1M+ elements | - | Functional annotations |
| Gene Ontology | 44K terms | - | 7M annotations |

### WORLD 2: Traditional Medicine (Expanded)

| Database | Formulas | Compounds | Targets | Tradition |
|----------|----------|-----------|---------|-----------|
| KampoDB | 481 | 3,002 | 62,906 | Kampo (Japan) |
| BATMAN-TCM | 54,832 | 39,171 | 2.3M TTIs | TCM (China) |
| IMPPAT | - | 17,967 | 5,042 | Ayurveda (India) |
| Dr. Duke's | - | 50K | 1,900 activities | Western Herbal |
| LOTUS | - | 750K+ pairs | - | NP-organism links |

### WORLD 3: Nutrition (Expanded)

| Database | Foods | Compounds | Key Feature |
|----------|-------|-----------|-------------|
| USDA FoodData | 20,900 | - | Authoritative nutrients |
| FooDB | 778 | 70,926 | Food metabolomics |
| DSLD | - | 200K labels | Supplement ingredients |
| Open Food Facts | 3M+ | - | Crowdsourced products |
| HMP | 5K+ metagenomes | - | Microbiome reference |

### Cross-Cutting (Expanded)

| Database | Primary Entity | Count | Cross-Refs |
|----------|----------------|-------|------------|
| ChEMBL | Bioactivities | 24.3M | UniProt, PubMed |
| COCONUT | Natural Products | 716K | NCBI Taxonomy |
| PubChem | Compounds | 115M+ | Universal |
| ChEBI | Chemical Entities | 60K+ | GO, Reactome |
| Reactome | Pathways | 2,712 | UniProt, ChEBI |
| STRING | Interactions | 68M+ | UniProt, GO |
| IntAct | Interactions | 1.1M+ | UniProt, ChEBI |
| MONDO | Diseases | 26K | 130K mappings |
| Wikidata | Entities | 100M+ | Universal |

---

## API Access Summary

| Database | API Type | Auth | Rate Limit | Batch |
|----------|----------|------|------------|-------|
| dbSNP | REST | None | 1/sec | 50K HGVS |
| ClinVar | E-utilities | API key | 10/sec | Yes |
| ChEMBL | REST (OpenAPI) | None | None stated | Yes |
| COCONUT | REST (Laravel) | Bearer (optional) | Unknown | Yes |
| Reactome | REST | None | None stated | Yes |
| WikiPathways | REST | None | None stated | Yes |
| Wikidata | SPARQL | None | Fair use | Yes |
| UniProt | REST | None | 1-3/sec | Yes |
| KampoDB | REST | None | None stated | No |
| GWAS Catalog | REST | None | None stated | Yes |
| gnomAD | REST/GraphQL | None | None stated | Yes |
| STRING | REST | None | None stated | Yes |
| IntAct | REST (PSICQUIC) | None | None stated | Yes |
| PubChem | REST (PUG) | None | 5/sec | Yes |
| ChEBI | REST/SOAP | None | None stated | Yes |
| Gene Ontology | REST (AmiGO) | None | None stated | Yes |
| Open Food Facts | REST | None | None stated | Yes |

---

## File Inventory (41 Documents)

```
schemas/
├── Analysis & Integration
│   ├── UNIFIED-SCHEMA-ANALYSIS.md    (898 lines)
│   ├── SCHEMAS-INDEX.md              (this file)
│   ├── INDEX.md                      (pathway index)
│   └── SAMPLE-DATA.md                (API samples)
│
├── WORLD 1: Genetics (11 files)
│   ├── DBSNP-SCHEMA.md               (502 lines)
│   ├── CLINVAR-SCHEMA.md             (620 lines)
│   ├── PHARMGKB-SCHEMA.md            (651 lines)
│   ├── GWAS-CATALOG-SCHEMA.md
│   ├── GNOMAD-SCHEMA.md              (Open Access)
│   ├── ALPHAMISSENSE-SCHEMA.md       (CC BY 4.0) [NEW]
│   ├── DBNSFP-SCHEMA.md              (Academic/Commercial) [NEW]
│   ├── DBVAR-SCHEMA.md               (Public Domain) [NEW]
│   ├── SPLICEAI-SCHEMA.md            (Illumina License) [NEW]
│   ├── ENCODE-SCHEMA.md              (Open Access) [NEW]
│   └── GENE-ONTOLOGY-SCHEMA.md       (CC BY 4.0) [NEW]
│
├── WORLD 2: Traditional Medicine (5 files)
│   ├── KAMPODB-SCHEMA.md             (CC BY-SA 4.0)
│   ├── DR-DUKES-SCHEMA.md            (CC0)
│   ├── BATMAN-TCM-SCHEMA.md          (CC BY-NC 4.0)
│   ├── IMPPAT-SCHEMA.md              (CC BY-NC 4.0)
│   └── LOTUS-SCHEMA.md               (CC0) [NEW]
│
├── WORLD 3: Nutrition (5 files)
│   ├── USDA-FOODDATA-CENTRAL.md      (CC0)
│   ├── FOODB.md
│   ├── DSLD-NIH.md
│   ├── OPEN-FOOD-FACTS-SCHEMA.md     (ODbL) [NEW]
│   └── HMP-SCHEMA.md                 (CC BY 4.0) [NEW]
│
├── Natural Products (3 files)
│   ├── CHEMBL-SCHEMA.md              (CC BY-SA 3.0)
│   ├── COCONUT-SCHEMA.md             (CC0)
│   └── PUBCHEM-SCHEMA.md             (Public Domain) [NEW]
│
├── Chemical Ontologies (1 file)
│   └── CHEBI-SCHEMA.md               (CC BY 4.0) [NEW]
│
├── Pathways & Networks (6 files)
│   ├── REACTOME-SCHEMA.md            (337 lines)
│   ├── WIKIPATHWAYS-GPML-SCHEMA.md   (420 lines)
│   ├── KGML-SCHEMA.md                (Academic License) [NEW]
│   ├── DISGENET-SCHEMA.md            (467 lines)
│   ├── STRING-SCHEMA.md              (CC BY 4.0) [NEW]
│   └── INTACT-SCHEMA.md              (CC BY 4.0) [NEW]
│
├── Diseases & Phenotypes (3 files)
│   ├── MONDO-SCHEMA.md               (459 lines)
│   ├── HPO-SCHEMA.md                 (582 lines)
│   └── ORPHANET-ORDO-SCHEMA.md
│
├── Health Domains (2 files)
│   ├── CBIOPORTAL-SCHEMA.md
│   └── GMREPO-SCHEMA.md
│
└── Integration (2 files)
    ├── WIKIDATA-SCHEMA.md            (412 lines)
    └── UNIPROT-IDMAPPING-SCHEMA.md   (574 lines)
```

---

## New Schema Files Summary (15 Additions)

| Schema File | Database | License | Category |
|-------------|----------|---------|----------|
| KGML-SCHEMA.md | KEGG KGML | Academic License | Pathways & Networks |
| GENE-ONTOLOGY-SCHEMA.md | Gene Ontology | CC BY 4.0 | WORLD 1: Genetics |
| CHEBI-SCHEMA.md | ChEBI | CC BY 4.0 | Chemical Ontologies |
| ALPHAMISSENSE-SCHEMA.md | AlphaMissense | CC BY 4.0 | WORLD 1: Genetics |
| DBNSFP-SCHEMA.md | dbNSFP | Academic/Commercial | WORLD 1: Genetics |
| GNOMAD-SCHEMA.md | gnomAD | Open Access | WORLD 1: Genetics |
| DBVAR-SCHEMA.md | dbVar | Public Domain | WORLD 1: Genetics |
| STRING-SCHEMA.md | STRING | CC BY 4.0 | Pathways & Networks |
| INTACT-SCHEMA.md | IntAct | CC BY 4.0 | Pathways & Networks |
| PUBCHEM-SCHEMA.md | PubChem | Public Domain | Natural Products |
| HMP-SCHEMA.md | Human Microbiome Project | CC BY 4.0 | WORLD 3: Nutrition |
| OPEN-FOOD-FACTS-SCHEMA.md | Open Food Facts | ODbL | WORLD 3: Nutrition |
| SPLICEAI-SCHEMA.md | SpliceAI | Illumina License | WORLD 1: Genetics |
| ENCODE-SCHEMA.md | ENCODE | Open Access | WORLD 1: Genetics |
| LOTUS-SCHEMA.md | LOTUS | CC0 | WORLD 2: Traditional Medicine |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-00-INDEX.md](../43-00-INDEX.md) | Parent navigation |
| [43-DATA-SOURCES.md](../../43-DATA-SOURCES.md) | Master data sources |
| [WORLD1-SCHEMA-RESEARCH.md](../research/WORLD1-SCHEMA-RESEARCH.md) | Genetics research notes |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial index with 3 NP schemas |
| 1.1 | January 2026 | Data Engineering | Added genetics, pathways, diseases |
| 1.2 | January 2026 | Data Engineering | Updated Tier 1 genetics with API data |
| 2.0 | January 2026 | Data Engineering | Complete overhaul: 27 schemas, THREE WORLDS coverage, prioritized shortlist, hub identifier strategy, implementation roadmap |
| **3.0** | January 2026 | Data Engineering | **Added 14 new schemas**: Gene Ontology, ChEBI, AlphaMissense, dbNSFP, gnomAD, dbVar, STRING, IntAct, PubChem, HMP, Open Food Facts, SpliceAI, ENCODE, LOTUS. Updated coverage to 41 files, 38 databases. Extended roadmap to 24 weeks. Added new license tiers and storage estimates. |
| **3.1** | January 2026 | Data Engineering | **Added KGML-SCHEMA.md** as reference documentation for KEGG pathway format. Updated coverage to 42 files. Added license notice and open alternatives recommendations. |
