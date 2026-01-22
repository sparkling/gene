# Environmental & Mitochondrial Data Sources

**Document ID:** 43-83-ENVIRONMENTAL-MITO
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [index.md](../index.md)

---

## TL;DR

Environmental and mitochondrial databases provide critical data for understanding gene-environment interactions and mitochondrial health. **CTD (50M+ relationships, free academic)** and **EPA CompTox (1M+ chemicals, public domain, free API)** are essential for toxicogenomics. **MITOMAP (62K+ sequences, MSeqDR API)** is the primary mtDNA variant resource. **Reactome (CC open, full API)** provides best-in-class detoxification pathway data. Prioritize Phase 1 integration with open-licensed sources (EPA CompTox, Reactome, T3DB, MITOMAP) before pursuing academic-licensed databases requiring commercial negotiations.

---

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **Toxicogenomics Primary** | CTD | 50M+ relationships, comprehensive chemical-gene-disease connections |
| **Environmental Exposure** | EPA CompTox | Public domain, 1M+ chemicals, free REST API |
| **Mitochondrial Variants** | MITOMAP + MSeqDR | Gold standard mtDNA resource, API via MSeqDR |
| **Detoxification Pathways** | Reactome | CC open license, comprehensive phase I/II coverage |
| **Heavy Metals** | CTD + T3DB | Integrated chemical-target-disease relationships |
| **Drug Metabolism** | DrugBank + SuperCYP | CYP450 focus, metabolism pathways |
| **Integration Order** | Open license first | Phase 1: public domain, Phase 2: academic, Phase 3: commercial |

---

## Database Catalog

### 1. Toxicogenomics Databases

#### 1.1 Comparative Toxicogenomics Database (CTD)

| Attribute | Value |
|-----------|-------|
| **Provider** | MDI Biological Laboratory and NC State University |
| **URL** | https://ctdbase.org |
| **Total Records** | 50+ million toxicogenomic relationships |
| **Annual Growth** | Monthly updates |
| **Date Range** | 1990s to present |
| **Content Type** | Chemical-gene-phenotype-disease relationships |
| **Relevance** | Critical - core toxicogenomics resource |
| **License** | Free for academic use; commercial licensing may apply |
| **Commercial Use** | REQUIRES LICENSE |

**Content Breakdown:**

| Entity Type | Count | Description |
|-------------|-------|-------------|
| Direct interactions | 3.8M | From 149K+ scientific articles |
| Chemicals | 17,700 | Environmental and pharmaceutical |
| Genes | 55,400 | Cross-species coverage |
| Phenotypes | 6,700 | GO-based phenotype classifications |
| Diseases | 7,200 | Disease associations |
| Exposures | 214,000 | Exposure study data |
| Species | 630+ | Multi-species coverage |

**Key Features:**
- CTD Tetramers (chemical-gene-phenotype-disease connections)
- Inferred relationships from curated data
- Batch query tools for bulk analysis
- Venn analysis for overlap identification

**API Access:**

| Method | Details |
|--------|---------|
| Batch Query | http://ctdbase.org/tools/batchQuery.go |
| Downloads | http://ctdbase.org/downloads/ |
| Web Services | Available (may require CAPTCHA) |
| Formats | CSV, TSV, XML |

**Use Cases:**
- Chemical exposure impact on gene expression
- Disease mechanisms from environmental exposures
- Drug-chemical interactions
- Air pollution/metal exposure linked to diseases

---

#### 1.2 T3DB - Toxin and Toxin Target Database

| Attribute | Value |
|-----------|-------|
| **Provider** | University of Alberta |
| **URL** | https://www.t3db.ca |
| **Total Records** | 3,678 toxins; 42,374 toxin-target associations |
| **Content Type** | Toxin data with target information |
| **Relevance** | High - essential toxin-target relationships |
| **License** | Non-proprietary; freely accessible; EULA for downloads |
| **Commercial Use** | YES (with EULA) |

**Content Breakdown:**

| Entity Type | Count | Description |
|-------------|-------|-------------|
| Toxins | 3,678 | Comprehensive toxin database |
| Synonyms | 41,602 | Alternative names |
| Target records | 2,073 | Protein target information |
| Associations | 42,374 | Toxin-target relationships |
| Data fields/toxin | 90+ | Comprehensive annotations |

**Toxin Categories:**
- Pollutants and pesticides
- Drugs and drug metabolites
- Food toxins and additives
- Household and workplace toxins
- Cigarette toxins
- Uremic toxins
- Heavy metals and metal salts

**ToxCard Content:**
- Chemical properties (structure, formula, weight)
- Toxicity values (LD50, LC50, etc.)
- Molecular interactions and targets
- Medical information and treatments
- Exposure routes and symptoms

---

### 2. Environmental Exposure Databases

#### 2.1 EPA CompTox Chemicals Dashboard

| Attribute | Value |
|-----------|-------|
| **Provider** | U.S. Environmental Protection Agency |
| **URL** | https://comptox.epa.gov/dashboard/ |
| **Total Records** | 1+ million chemicals |
| **Content Type** | Chemistry, toxicity, exposure data |
| **Relevance** | High - government-backed exposure data |
| **License** | PUBLIC DOMAIN |
| **Commercial Use** | YES - fully open |

**Content Breakdown:**

| Data Type | Coverage | Description |
|-----------|----------|-------------|
| Chemicals | 1M+ | Comprehensive chemical library |
| Chemical lists | 300+ | Structure/category groupings |
| Physicochemical properties | Full | Properties and predictions |
| Exposure data | Extensive | CPDat predictions |
| Bioactivity data | ToxCast | In vitro bioassay results |

**API Access (CTX APIs):**

| Endpoint | Function | Cost |
|----------|----------|------|
| Chemical API | Chemical search, properties | FREE |
| Exposure API | Exposure predictions (CPDat) | FREE |
| Hazard API | Hazard characterization | FREE |
| Bioactivity API | ToxCast bioassay data | FREE |

**API Details:**
- Base URL: Hosted on cloud.gov
- Authentication: Free API key (contact ccte_api@epa.gov)
- Format: JSON via REST
- Documentation: https://www.epa.gov/comptox-tools/computational-toxicology-and-exposure-apis

---

#### 2.2 HERCULES Exposome Research Center

| Attribute | Value |
|-----------|-------|
| **Provider** | Emory University (NIEHS-funded) |
| **URL** | https://emoryhercules.com |
| **Total Records** | 20,000+ chemical signals per sample |
| **Content Type** | Untargeted metabolomics data |
| **Relevance** | Medium - research collaboration focus |
| **License** | Academic research collaboration |
| **Commercial Use** | BY COLLABORATION |

**Coverage:**
- Environmental chemicals
- Food metabolites
- Microbiome-derived metabolites
- 20,000+ chemical signals per sample

**Access:** Collaborative access through facility cores; no public API.

---

#### 2.3 NEXUS (NIH Exposome Coordinating Center)

| Attribute | Value |
|-----------|-------|
| **Provider** | Columbia, Harvard, USC (NIH-funded) |
| **URL** | Via NIEHS |
| **Funding** | $7.7 million (September 2024) |
| **Content Type** | Environmental monitoring, wearables, biospecimens |
| **Relevance** | Medium - emerging resource |
| **License** | Expected publicly accessible |
| **Commercial Use** | TBD |

**Status:** Under development as part of NIH Real World Data Platform. Monitor for availability.

---

### 3. Mitochondrial Variant Databases

#### 3.1 MITOMAP

| Attribute | Value |
|-----------|-------|
| **Provider** | CHOP Center for Mitochondrial & Epigenomic Medicine |
| **URL** | https://www.mitomap.org |
| **Total Records** | 62,556 full-length GenBank sequences |
| **Content Type** | mtDNA polymorphisms and mutations |
| **Relevance** | Critical - primary mtDNA resource |
| **License** | Content property of contributing authors |
| **Commercial Use** | CHECK LICENSE |

**Content Breakdown:**

| Entity Type | Count | Description |
|-------------|-------|-------------|
| Full-length sequences | 62,556 | Complete mtDNA genomes |
| Control region sequences | 81,778 | D-loop sequences |
| SNVs | 19,892 | Single nucleotide variants |
| Disease mutations | Extensive | Published mutation reports |
| Deletions/rearrangements | Included | Complex mtDNA changes |

**Haplogroup Distribution (62,556 FL sequences):**

| Haplogroup | Percentage | Description |
|------------|------------|-------------|
| N | 68% | European/Asian lineages |
| L | 10% | African lineages |
| M | 20% | Asian lineages |

**Key Tools:**
- Allele Search: Variant lookup
- MITOMASTER: Sequence analysis, variant ID, haplogroup prediction

**API Access via MSeqDR mvTool:**

| Feature | Details |
|---------|---------|
| URL | https://mseqdr.org/mvtool.php |
| Input formats | VCF, HGVS, classical mtDNA nomenclature |
| Output formats | HTML, JSON, Excel, VCF |
| Integrations | MITOMAP, HmtDB, GeneDx, 1000 Genomes |
| Cost | FREE |

---

#### 3.2 HmtDB (Human Mitochondrial Database)

| Attribute | Value |
|-----------|-------|
| **Provider** | University of Bari, Italy |
| **URL** | https://www.hmtdb.uniba.it |
| **Total Records** | 3,700+ mitochondrial genomes (tripled since launch) |
| **Content Type** | mtDNA with phenotype annotations |
| **Relevance** | Medium - complementary to MITOMAP |
| **License** | Academic research use |
| **Commercial Use** | ACADEMIC ONLY |

**Key Features:**
- Healthy and disease phenotype sequences
- Site-specific nucleotide variability
- Continental subset analysis
- Phylotree-based haplogroup classification
- MToolBox pipeline integration

---

### 4. Detoxification Pathway Databases

#### 4.1 KEGG (Kyoto Encyclopedia of Genes and Genomes)

| Attribute | Value |
|-----------|-------|
| **Provider** | Kanehisa Laboratories, Kyoto University |
| **URL** | https://www.kegg.jp |
| **API URL** | https://rest.kegg.jp |
| **Content Type** | Pathway maps, compounds, enzymes, reactions |
| **Relevance** | High - essential pathway resource |
| **License** | Academic use only; commercial requires license |
| **Commercial Use** | REQUIRES LICENSE (Pathway Solutions) |

**Key Pathways:**

| Pathway ID | Name | Relevance |
|------------|------|-----------|
| hsa00980 | Drug metabolism - CYP450 | Critical |
| hsa00982 | Drug metabolism - other enzymes | Critical |
| hsa00480 | Glutathione metabolism | High |
| hsa00627 | Aminobenzoate degradation | Medium |

**API Examples:**
```
https://rest.kegg.jp/list/pathway/hsa    # List human pathways
https://rest.kegg.jp/get/hsa00980        # CYP450 pathway
https://rest.kegg.jp/find/compound/aspirin
```

**Libraries:**
- Python: Biopython Bio.KEGG.REST, kegg_pull
- R: KEGGREST (Bioconductor)

---

#### 4.2 Reactome

| Attribute | Value |
|-----------|-------|
| **Provider** | Ontario Institute for Cancer Research / EMBL-EBI |
| **URL** | https://reactome.org |
| **Downloads** | https://reactome.org/download-data |
| **Content Type** | Curated biological pathways |
| **Relevance** | High - open license, comprehensive |
| **License** | CREATIVE COMMONS (Open) |
| **Commercial Use** | YES |

**Key Pathway Coverage:**

| Pathway ID | Name | Description |
|------------|------|-------------|
| R-HSA-211859 | Biological oxidations | Phase I/II metabolism |
| R-HSA-211981 | Xenobiotics | Xenobiotic biotransformation |
| R-HSA-1430728 | Metabolism | Comprehensive metabolism |
| R-HSA-xxx | Detoxification of ROS | Reactive oxygen species |

**API Access:**

| Service | URL | Function |
|---------|-----|----------|
| Content Service | https://reactome.org/ContentService/ | Pathway data access |
| Analysis Service | https://reactome.org/AnalysisService/ | Pathway analysis |

**Data Formats:**
- BioPAX Level 3
- SBML Level 3 v1
- SBGN
- PSI-MITAB
- Neo4j GraphDB
- MySQL dumps

**Mapping Files:** UniProt, ChEBI, ENSEMBL, miRBase, NCBI, GtoP identifiers

---

#### 4.3 DrugBank

| Attribute | Value |
|-----------|-------|
| **Provider** | University of Alberta / OMx Personal Health Analytics |
| **URL** | https://go.drugbank.com |
| **Total Records** | 50,000+ drug entries |
| **Content Type** | Drug, drug-target, metabolism data |
| **Relevance** | High - gold standard drug reference |
| **License** | CC BY-NC 4.0 (datasets); CC0 (Open Data) |
| **Commercial Use** | ACADEMIC FREE / COMMERCIAL PAID |

**DrugBank 6.0 Content (2024):**

| Entity Type | Count | Description |
|-------------|-------|-------------|
| FDA-approved drugs | 4,563 | Approved medications |
| Investigational drugs | 6,231 | Pipeline compounds |
| Drug metabolites | 3,037 | With structures |
| Drug-drug interactions | 1.4M+ | Interaction pairs |
| Drug-food interactions | 2,475 | Food interactions |
| Drug-enzyme entries | 2,550 | Metabolism enzymes |
| Drug-transporter entries | 1,560 | Transport proteins |

**API Access:**
- Clinical API and Discovery API available
- Documentation: https://docs.drugbank.com/

---

#### 4.4 SuperCYP

| Attribute | Value |
|-----------|-------|
| **Provider** | Charite - Universitatsmedizin Berlin |
| **URL** | http://bioinformatics.charite.de/supercyp |
| **Prediction Tool** | http://insilico-cyp.charite.de/SuperCYPsPred/ |
| **Content Type** | CYP450 enzyme interactions |
| **Relevance** | Medium - specialized CYP resource |
| **License** | CC BY-NC-SA 3.0 |
| **Commercial Use** | NON-COMMERCIAL ONLY |

**Content:**

| Entity Type | Count | Description |
|-------------|-------|-------------|
| Drugs | 1,170 | Drug compounds |
| CYP interactions | 3,800+ | Drug-CYP relationships |
| SNPs/mutations | 2,000+ | Expression/activity effects |
| 3D homology models | 48 | All human CYPs |

**CYP Coverage:**
- CYP1A2, 2C9, 2C19, 2D6, 3A4 (90%+ of clinical drug metabolism)

**SuperCYPsPred Features:**
- Input: 2D chemical structure
- Output: CYP inhibition profile, confidence scores
- DDI tables, radar charts
- Free, no registration required

---

#### 4.5 PharmGKB / ClinPGx

| Attribute | Value |
|-----------|-------|
| **Provider** | Stanford University |
| **URL** | https://www.pharmgkb.org (redirects to https://www.clinpgx.org) |
| **Content Type** | Pharmacogenomics knowledge |
| **Relevance** | Medium - PGx focus |
| **License** | Academic use permitted |
| **Commercial Use** | CHECK CURRENT TERMS |

**Content:**
- Drug-gene associations
- Clinical annotations
- Dosing guidelines
- Pharmacokinetic pathways
- Variant annotations

---

### 5. Heavy Metal Interaction Databases

#### 5.1 CTD Heavy Metals Coverage

| Heavy Metal | Symbol | CTD Coverage |
|-------------|--------|--------------|
| Lead | Pb | Gene interactions, disease associations, exposure |
| Cadmium | Cd | Toxicogenomic relationships, phenotypes |
| Mercury | Hg | Methylmercury and inorganic forms |
| Arsenic | As | Pathway mechanisms to diseases |
| Aluminum | Al | Neurological associations |
| Chromium | Cr | Carcinogenicity data |

**CTD Tetramer Example:** Air pollution/metal exposure -> intermediate genes/phenotypes -> Alzheimer disease

---

#### 5.2 GWAS Data for Heavy Metals

| Attribute | Value |
|-----------|-------|
| **Source** | GWAS Catalog (https://www.ebi.ac.uk/gwas/) |
| **Content** | Genetic variants for metal blood levels |
| **Metals Studied** | Al, Cd, Co, Cu, Cr, Hg, Mn, Mo, Ni, Pb, Zn |
| **Notable Gene** | SLC39A8 (rs13107325) - Mn and Zn transport |

---

## Integration Priority

| Database | Priority | License | API | Key Use Case |
|----------|----------|---------|-----|--------------|
| EPA CompTox | HIGH | **Public Domain** | REST API (free) | Environmental exposure |
| Reactome | HIGH | **CC Open** | REST API | Detox pathways |
| T3DB | HIGH | Free (EULA) | Web/Download | Toxin targets |
| MITOMAP | HIGH | Author-owned | MSeqDR API | mtDNA variants |
| CTD | HIGH | Academic free | Batch/Download | Chemical-gene-disease |
| DrugBank | HIGH | CC BY-NC 4.0 | REST API | Drug metabolism |
| KEGG | HIGH | Academic only | REST API | Pathways (license needed) |
| HmtDB | MEDIUM | Academic | Web/Download | mtDNA population |
| SuperCYP | MEDIUM | CC BY-NC-SA | Web/Download | CYP450 focus |
| HERCULES | MEDIUM | Collaborative | None | Exposome research |

---

## Recommended Integration Order

### Phase 1 - Core (Open/Free License)

| Database | License | Method | Priority |
|----------|---------|--------|----------|
| EPA CompTox | Public domain | CTX REST APIs | Week 1-2 |
| Reactome | CC Open | Content/Analysis APIs | Week 2-3 |
| T3DB | Free (EULA) | Download + parse | Week 3-4 |
| MITOMAP | Author-owned | MSeqDR mvTool API | Week 4 |

### Phase 2 - Academic License

| Database | License | Method | Priority |
|----------|---------|--------|----------|
| CTD | Academic free | Batch downloads | Month 2 |
| KEGG | Academic free | REST API | Month 2 |
| HmtDB | Academic | Download + MToolBox | Month 3 |

### Phase 3 - Commercial Considerations

| Database | License | Method | Priority |
|----------|---------|--------|----------|
| DrugBank | Commercial paid | REST API | Month 4+ |
| SuperCYP | CC BY-NC-SA | Download | Month 4+ |
| PharmGKB | Academic | Download | Month 5+ |

---

## Data Volume Estimates

| Database | Records | Storage | Update Frequency |
|----------|---------|---------|------------------|
| CTD | 50M+ relationships | ~5 GB | Monthly |
| EPA CompTox | 1M+ chemicals | ~2 GB | Quarterly |
| MITOMAP | 62K sequences, 20K SNVs | ~500 MB | 4-6 months |
| Reactome | ~15,000 pathways | ~1 GB | Quarterly |
| T3DB | 42K associations | ~200 MB | Periodic |
| DrugBank | 50K+ entries | ~1 GB | Periodic |
| KEGG | Varies by DB | ~2 GB | Ongoing |
| HmtDB | 3,700+ genomes | ~100 MB | Periodic |

**Total Estimated Storage:** ~12 GB (compressed)

---

## API Integration Notes

### Free API Keys Required

| Database | Contact | Turnaround |
|----------|---------|------------|
| EPA CompTox | ccte_api@epa.gov | 1-2 weeks |

### No API Key Needed

| Database | Notes |
|----------|-------|
| Reactome | Content/Analysis Services open |
| KEGG REST | Academic use |
| MSeqDR mvTool | Open access |

### Download-Based Access

| Database | Requirement | Notes |
|----------|-------------|-------|
| CTD | Web download | May need CAPTCHA |
| T3DB | EULA agreement | Library files |
| DrugBank | Account required | Academic/commercial tiers |

---

## Dependencies

### Upstream Dependencies

| Dependency | Purpose | Risk if Unavailable |
|------------|---------|---------------------|
| EPA CompTox APIs | Exposure data | MEDIUM - alternatives exist |
| MSeqDR mvTool | mtDNA variant lookup | HIGH - primary MITOMAP API |
| Reactome API | Pathway data | LOW - downloads available |
| CTD Downloads | Toxicogenomics | HIGH - essential source |

### Downstream Dependents

| Dependent | Usage |
|-----------|-------|
| Environmental Risk Scoring | Chemical exposure assessment |
| Mitochondrial Health Module | mtDNA variant interpretation |
| Detox Pathway Analysis | CYP450 and phase I/II metabolism |
| Heavy Metal Impact | Metal-gene-disease relationships |
| Drug Interaction Checker | CYP metabolism predictions |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Initial document from research.old/data-sources-environmental-mito.md |
