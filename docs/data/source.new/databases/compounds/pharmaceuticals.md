---
id: compounds-pharmaceuticals
title: Pharmaceutical Databases
world: null
category: compounds
subcategory: pharmaceuticals
tier: 1
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
databases:
  - pharmgkb
  - cpic
  - drugbank
  - chembl
  - pubchem
  - drugcentral
  - open-targets
  - dgidb
  - ttd
  - bindingdb
  - gtopdb
  - dailymed
  - openfda
  - fda-pgx-table
  - rxnorm
  - civic
  - pgxdb
  - superdrug2
  - kegg-drug
  - uniprot
tags: [pharmaceuticals, drugs, interactions, pharmacogenomics]
---

# Pharmaceutical and Pharmacogenomics Data Sources

**Document ID:** 43-51-PHARMACEUTICALS
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

Pharmaceutical and pharmacogenomics databases provide the foundation for SNP-drug relationship mapping, clinical dosing recommendations, and drug-target interactions. PharmGKB/CPIC serve as the gold standard for pharmacogenomics with clinically validated dosing guidelines, while DrugBank and ChEMBL provide comprehensive drug information and bioactivity data. Integration of 22+ databases enables complete coverage from variant annotation through regulatory compliance.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary PGx source | PharmGKB + CPIC + DPWG | Gold standard with clinical validation | Jan 2026 |
| Drug information | DrugBank (academic) + ChEMBL | Comprehensive targets + open bioactivity | Jan 2026 |
| Drug-gene aggregator | DGIdb + Open Targets | Free aggregation from 40+ sources | Jan 2026 |
| Regulatory data | OpenFDA + DailyMed | Authoritative FDA labels and adverse events | Jan 2026 |
| Binding affinity | BindingDB + ChEMBL | Quantitative IC50/Ki/Kd data | Jan 2026 |
| Allele nomenclature | PharmVar | Official star allele definitions | Jan 2026 |

---

## Database Catalog

### Primary Pharmacogenomics Databases

#### 1. PharmGKB / ClinPGx

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org / https://www.clinpgx.org |
| **Maintainer** | Stanford University (NIH-funded) |
| **Content** | Clinical annotations linking variants to drug responses, dosing guidelines, pathways |
| **Records** | 700+ genes, 1,000+ drugs, 3,000+ clinical annotations, 20,000+ variant annotations, 150+ pathways |
| **License** | CC BY-SA 4.0 (free for academic and commercial with attribution) |
| **API** | RESTful API (JSON) - https://api.pharmgkb.org/ |
| **Documentation** | https://api.pharmgkb.org/swagger/ |
| **Bulk Download** | Yes (requires account and license agreement) |
| **Data Formats** | TSV/CSV, JSON (API), BioPax (pathways) |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB (full download) |

**Evidence Levels:**
- **1A:** CPIC/medical society guideline or implemented clinically
- **1B:** Preponderance of evidence shows association
- **2A:** Level 2B + Very Important Pharmacogene (VIP)
- **2B:** Moderate evidence, replicated but some non-significant studies
- **3:** Single significant study or lacking clear evidence
- **4:** Case report, non-significant study, or in vitro evidence only

**Schema:**
```
Gene <---> Variant <---> Clinical Annotation <---> Drug
                              |
                              v
                    Level of Evidence (1A-4)
                              |
                              v
                    Phenotype Categories
```

---

#### 2. CPIC (Clinical Pharmacogenetics Implementation Consortium)

| Field | Value |
|-------|-------|
| **URL** | https://cpicpgx.org |
| **Maintainer** | PharmGKB/PGRN consortium |
| **Content** | Evidence-based gene/drug clinical practice guidelines, allele tables, dosing recommendations |
| **Records** | 34 genes, 164 drugs, 27 published guidelines, 128+ healthcare institutions |
| **License** | CC0 (Public Domain) - completely free for any use |
| **API** | Via PharmGKB (80,000+ monthly queries) |
| **Bulk Download** | JSON format via PharmGKB |
| **Data Formats** | PDF (guidelines), Excel/TSV (tables), JSON |
| **Update Frequency** | As guidelines published |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |

**Key Tables:**
1. Allele Definition Table (variant positions, defining variants)
2. Allele Functionality Table (function assignments)
3. Diplotype-Phenotype Table
4. Dosing Recommendation Table

**Schema:**
```
Gene --> Allele Definitions --> Diplotype --> Phenotype --> Dosing Recommendation
                                                   |
                                                   v
                                        Strength (Strong/Moderate/Optional)
```

---

#### 3. DPWG (Dutch Pharmacogenetics Working Group)

| Field | Value |
|-------|-------|
| **URL** | https://www.knmp.nl |
| **Maintainer** | Royal Dutch Pharmacists Association (KNMP) |
| **Content** | Gene-drug interaction guidelines, therapeutic recommendations |
| **Records** | Guidelines for CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP1A2, DPYD, TPMT, NUDT15, UGT1A1, HLA-A, HLA-B |
| **License** | Open (published in peer-reviewed journals) |
| **API** | Via PharmGKB |
| **Bulk Download** | Via PharmGKB annotations |
| **Data Formats** | PDF, integrated into PharmGKB |
| **Update Frequency** | Regular updates (2024: antipsychotics, antidepressants, anti-epileptics) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 MB (via PharmGKB) |

**Recent Guidelines (2024):**
- CYP2D6, CYP3A4, CYP1A2 and antipsychotics (March 2024)
- CYP2D6, CYP2C19 and non-SSRI/non-TCA antidepressants (November 2024)
- CYP2C9, HLA-A, HLA-B with anti-epileptic drugs (August 2024)

---

#### 4. PharmVar (Pharmacogene Variation Consortium)

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmvar.org |
| **Maintainer** | PharmVar Consortium |
| **Content** | Official pharmacogene allele nomenclature, star allele definitions, haplotypes |
| **Records** | All major CYP450 enzymes (CYP2D6, CYP2C19, CYP2C9, CYP3A4), UGT1A1, DPYD, TPMT |
| **License** | Open access for research and clinical use |
| **API** | None (download-based) |
| **Bulk Download** | Gene-specific allele definition files |
| **Data Formats** | TSV/text files, Change logs |
| **Update Frequency** | Monthly releases (version 6.2+) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 MB |

**Recent Updates (2024-2025):**
- CYP1A2 completed transition from legacy nomenclature (December 2024)

**Schema:**
```
Gene --> Star Allele --> Defining Variants (rsIDs) --> Function Assignment
                              |
                              v
                    Reference Sequence Positions
```

---

### Comprehensive Drug Databases

#### 5. DrugBank

| Field | Value |
|-------|-------|
| **URL** | https://go.drugbank.com |
| **Maintainer** | DrugBank (OMx Personal Health Analytics) |
| **Content** | Comprehensive drug information, targets, enzymes, transporters, interactions, pharmacokinetics |
| **Records** | 2,700+ approved, 6,700+ experimental, 15,000+ total drugs, 5,000+ targets, 2.5M+ DDIs |
| **License** | Open Data (CC0) for vocabulary/structures; Academic (CC BY-NC 4.0) for full data; Commercial requires license |
| **API** | REST API (subscription for full access) |
| **Bulk Download** | XML, CSV, SDF, FASTA |
| **Data Formats** | XML (comprehensive), CSV, SDF (structures), FASTA (sequences) |
| **Update Frequency** | Regular (Version 6.0, 2024) |
| **Priority** | Tier 1 (MVP) - Academic license |
| **Storage Estimate** | ~5 GB (full XML) |

**Key Fields:**
- DrugBank ID, Name, Type, Description, Indication
- Mechanism of Action, Pharmacodynamics
- Targets (inhibitor, agonist, antagonist actions)
- Enzymes, Carriers, Transporters
- Drug-Drug Interactions, Food Interactions
- SNP Effects (pharmacogenomics)

**Schema:**
```
Drug --> Targets (proteins) --> Gene
  |          |
  |          v
  |     UniProt ID, Actions, Known Action
  |
  +--> Drug Interactions
  |
  +--> Metabolizing Enzymes
  |
  +--> Pharmacokinetics (ADMET)
```

---

#### 6. ChEMBL

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/chembl |
| **Maintainer** | EMBL-EBI |
| **Content** | Bioactivity data, drug-target binding, assays, ADMET properties, mechanisms |
| **Records** | 2.8M compounds, 21M+ bioactivities, 1.6M+ assays, 17,803 targets, 95,000+ publications |
| **License** | Open Access (Apache 2.0 for software) |
| **API** | REST API (full programmatic access) |
| **Bulk Download** | FTP/HTTPS |
| **Data Formats** | SQLite, MySQL, PostgreSQL, SDF, FASTA, RDF/Turtle, JSON/XML/YAML |
| **Update Frequency** | Regular (ChEMBL 36, July 2025) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~20 GB (full database) |

**Key Tables:**
- molecule_dictionary (compounds)
- assays (experimental setup)
- activities (IC50, Ki, EC50 measurements)
- target_dictionary (proteins, organisms)
- mechanism (drug mechanisms of action)

**Schema:**
```
Document --> Assay --> Activity --> Compound
               |                        |
               v                        v
            Target                  CHEMBL ID
               |                   Structure (SMILES, InChI)
               v
          UniProt ID
```

---

#### 7. PubChem

| Field | Value |
|-------|-------|
| **URL** | https://pubchem.ncbi.nlm.nih.gov |
| **Maintainer** | NIH/NCBI |
| **Content** | Chemical structures, biological activities, bioassays, patents, safety data |
| **Records** | 118M+ compounds, 326M+ substances, 1.5M+ bioassays, 19M+ patent compounds |
| **License** | Public Domain (no restrictions) |
| **API** | PUG REST API |
| **Bulk Download** | FTP site |
| **Data Formats** | SDF, JSON, XML, CSV, ASN.1 |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 GB (full) |

**Schema:**
```
Compound (CID) <--> Substance (SID) <--> Source
       |
       v
   BioAssay (AID) --> Activity Data
       |
       v
   Target (Gene/Protein)
```

---

#### 8. DrugCentral

| Field | Value |
|-------|-------|
| **URL** | https://drugcentral.org |
| **Maintainer** | University of New Mexico, Translational Informatics Division |
| **Content** | FDA/EMA/PMDA approved drugs, structures, mechanisms, indications, FAERS adverse events, veterinary drugs |
| **Records** | Active pharmaceutical ingredients with approval status |
| **License** | Open Access (free for all uses) |
| **API** | PostgreSQL database access |
| **Database Connection** | Host: unmtid-dbs.net, Port: 5433, DB: drugcentral, User: drugman |
| **Bulk Download** | PostgreSQL dump, Docker container |
| **Data Formats** | PostgreSQL, JSON (API) |
| **Update Frequency** | Regular (veterinary added 2023) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 GB |

**Schema:**
```
Drug --> Target --> Gene
  |          |
  |          v
  |     Action Type, Bioactivity
  |
  +--> Indications (diseases)
  |
  +--> Adverse Events (FAERS)
```

---

### Drug-Target Interaction Databases

#### 9. Open Targets Platform

| Field | Value |
|-------|-------|
| **URL** | https://platform.opentargets.org |
| **Maintainer** | Open Targets consortium (Wellcome Sanger, EMBL-EBI, GSK) |
| **Content** | Target-disease associations, drug-target relationships, genetic associations, pathways |
| **Records** | Genome-wide evidence aggregation from multiple sources |
| **License** | CC0 for data, open source tools |
| **API** | GraphQL API |
| **Bulk Download** | FTP, Google Cloud, AWS, Azure |
| **Data Formats** | Parquet (bulk), JSON (API), TSV (exports) |
| **Update Frequency** | Regular |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 GB (full Parquet) |

**Schema:**
```
Target (gene) <--> Disease <--> Drug
       |               |
       v               v
  Evidence Score   Evidence Sources
```

---

#### 10. DGIdb (Drug Gene Interaction Database)

| Field | Value |
|-------|-------|
| **URL** | https://dgidb.org |
| **Maintainer** | Washington University School of Medicine |
| **Content** | Drug-gene interactions aggregated from 40+ sources, druggable gene categories |
| **Records** | 10,000+ genes, 20,000+ drugs, 70,000+ interactions |
| **License** | Open Access (free for all uses) |
| **API** | GraphQL API |
| **Bulk Download** | TSV files |
| **Data Formats** | TSV, JSON, PostgreSQL backend |
| **Update Frequency** | Regular (v5.0, 2024) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~1 GB |

**Schema:**
```
Gene <--> Interaction <--> Drug
              |
              v
        Interaction Type
        (inhibitor, agonist, etc.)
              |
              v
         Source Database
```

---

#### 11. TTD (Therapeutic Target Database)

| Field | Value |
|-------|-------|
| **URL** | https://idrblab.net/ttd |
| **Maintainer** | Zhejiang University |
| **Content** | Therapeutic protein/nucleic acid targets, druggability, ligands, pathways, diseases |
| **Records** | 426 successful targets, 1,014 clinical trial targets, 212 preclinical, 1,479 literature |
| **License** | Free Access (no login required) |
| **API** | None |
| **Bulk Download** | https://idrblab.net/ttd/full-data-download |
| **Data Formats** | Downloadable files, MySQL backend |
| **Update Frequency** | Regular (2024) |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~500 MB |

**Schema:**
```
Target --> Drugs/Ligands
   |
   v
Disease --> Pathway
   |
   v
Druggability Metrics
```

---

#### 12. BindingDB

| Field | Value |
|-------|-------|
| **URL** | https://www.bindingdb.org |
| **Maintainer** | UC San Diego |
| **Content** | Experimentally measured binding affinities (IC50, Ki, Kd, EC50) |
| **Records** | 2.9M binding measurements, 1.3M compounds |
| **License** | Open Access (FAIR compliant) |
| **API** | REST API (JSON/XML) |
| **Bulk Download** | SDF, TSV files |
| **Data Formats** | SDF, TSV, FASTA, JSON/XML |
| **Update Frequency** | Regular (2024) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~10 GB |

**API Examples:**
```
# By UniProt IDs:
http://bindingdb.org/rest/getLigandsByUniprots?uniprot=P00176,P00183&cutoff=10000&response=application/json

# By SMILES structure:
http://bindingdb.org/rest/getLigandsBySmiles?smiles=...&similarity=0.8
```

**Schema:**
```
Compound (SMILES) <--> Binding Data <--> Target (UniProt)
                           |
                           v
                    Affinity Value (IC50, Ki, Kd)
                    Assay Type
                    Publication
```

---

#### 13. GtoPdb (IUPHAR/BPS Guide to PHARMACOLOGY)

| Field | Value |
|-------|-------|
| **URL** | https://www.guidetopharmacology.org |
| **Maintainer** | IUPHAR / BPS |
| **Content** | Expert-curated pharmacological targets, ligands, quantitative interactions |
| **Records** | 3,097 human targets, 1,758 with interactions, 13,131 ligands, 9,485 with interactions |
| **License** | Open Access (free for all uses) |
| **API** | REST API (JSON) |
| **Bulk Download** | CSV/TSV, PostgreSQL dump, RDF |
| **Data Formats** | CSV/TSV, PostgreSQL, RDF/N3, JSON |
| **Update Frequency** | Regular (Version 2025.4) |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 GB |

---

### Regulatory and Clinical Sources

#### 14. FDA DailyMed

| Field | Value |
|-------|-------|
| **URL** | https://dailymed.nlm.nih.gov |
| **Maintainer** | NIH/NLM |
| **Content** | FDA-approved drug labels (package inserts), SPL files, OTC labels, animal drug labels |
| **Records** | 155,000+ drug labels |
| **License** | Public Domain (no restrictions) |
| **API** | None (bulk download) |
| **Bulk Download** | Daily, weekly, monthly updates via FTP |
| **Data Formats** | XML (SPL format), ZIP archives |
| **Update Frequency** | Daily |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 GB |

---

#### 15. OpenFDA

| Field | Value |
|-------|-------|
| **URL** | https://open.fda.gov |
| **Maintainer** | FDA |
| **Content** | Drug adverse events (FAERS), drug labels, enforcement reports, NDC directory, shortages |
| **Records** | 200M+ API calls total, millions of adverse event reports |
| **License** | Public Domain (no restrictions) |
| **API** | REST API (free, rate-limited) |
| **Bulk Download** | Zipped JSON |
| **Data Formats** | JSON |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 GB |

**API Examples:**
```
# Adverse events for a drug:
https://api.fda.gov/drug/event.json?search=patient.drug.medicinalproduct:"aspirin"

# Drug labels:
https://api.fda.gov/drug/label.json?search=openfda.brand_name:"Prozac"
```

---

#### 16. FDA Table of Pharmacogenomic Biomarkers

| Field | Value |
|-------|-------|
| **URL** | https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling |
| **Maintainer** | FDA Division of Translational and Precision Medicine |
| **Content** | Drugs with pharmacogenomic biomarkers in labeling, gene-drug associations |
| **Records** | 541 drugs with PGx biomarkers (January 2025) |
| **License** | Public Domain |
| **API** | None |
| **Bulk Download** | PDF download |
| **Data Formats** | PDF, HTML table |
| **Update Frequency** | Periodic (last: June 2024) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~10 MB |

---

#### 17. RxNorm

| Field | Value |
|-------|-------|
| **URL** | https://www.nlm.nih.gov/research/umls/rxnorm |
| **Maintainer** | NIH/NLM |
| **Content** | Normalized drug names, vocabulary links, drug-drug interactions (via DrugBank) |
| **Records** | Links to First Databank, Micromedex, Gold Standard, Multum, DrugBank |
| **License** | Free for prescribable subset (SAB=RXNORM, SAB=MTHSPL) |
| **API** | RxNorm API, RxClass API via RxNav |
| **Web Browser** | https://lhncbc.nlm.nih.gov/RxNav/ |
| **Bulk Download** | RRF files (Rich Release Format) |
| **Data Formats** | RRF (pipe-delimited), UTF-8, MySQL/Oracle scripts |
| **Update Frequency** | Monthly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 GB |

---

### Specialized Databases

#### 18. CIViC (Clinical Interpretation of Variants in Cancer)

| Field | Value |
|-------|-------|
| **URL** | https://civicdb.org |
| **Maintainer** | Washington University, Griffith Lab |
| **Content** | Clinical interpretation of cancer variants, therapeutic/prognostic/diagnostic relevance |
| **Records** | 3,200+ variants, 470+ genes, 3,100+ publications, 300+ contributors |
| **License** | CC0 (Public Domain), MIT for software |
| **API** | REST API (full access) |
| **Bulk Download** | Nightly updates, monthly stable releases, AWS Open Data |
| **Data Formats** | JSON, TSV |
| **Update Frequency** | Nightly/Monthly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |

---

#### 19. PGxDB

| Field | Value |
|-------|-------|
| **URL** | https://pgx-db.org |
| **Maintainer** | Academic consortium (launched August 2024) |
| **Content** | Integrated pharmacogenomics data, drug-gene-variant associations, adverse reactions |
| **Records** | Integrated from DrugBank, PharmGKB, ChEMBL, UniProt, Ensembl, SIDER, DisGeNet, ClinicalTrials.gov |
| **License** | Open Access (FAIR compliant) |
| **API** | Programmatic access |
| **Bulk Download** | GitHub |
| **Data Formats** | Web tools, downloadable files |
| **Update Frequency** | Active development |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~1 GB |

---

#### 20. SuperDRUG2

| Field | Value |
|-------|-------|
| **URL** | http://cheminfo.charite.de/superdrug2 |
| **Maintainer** | Charite Berlin |
| **Content** | Approved/marketed drugs, 2D/3D structures, targets, PK data, DDIs, side effects |
| **Records** | 4,587 active pharmaceutical ingredients |
| **License** | Open Access |
| **API** | None |
| **Bulk Download** | Customized download links |
| **Data Formats** | MySQL, downloadable files |
| **Update Frequency** | Periodic |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~500 MB |

---

#### 21. KEGG Drug

| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/drug/ |
| **Maintainer** | Kanehisa Laboratories |
| **Content** | Approved drugs (Japan, USA, Europe), structures, targets, metabolizing enzymes |
| **Records** | Multi-regional drug coverage |
| **License** | API: Academic use only; Medicus: Free; Full FTP: Subscription |
| **API** | KEGG REST API (academic) |
| **Bulk Download** | Paid subscription for most data |
| **Data Formats** | KEGG flat file format |
| **Update Frequency** | Regular |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~1 GB |

---

#### 22. UniProt (Drug Target Reference)

| Field | Value |
|-------|-------|
| **URL** | https://www.uniprot.org |
| **Maintainer** | UniProt Consortium (EBI, SIB, PIR) |
| **Content** | Protein sequences, functions, PTMs, disease associations, drug binding sites |
| **Records** | 246M sequences in UniProtKB (2024) |
| **License** | CC BY 4.0 (Open Access) |
| **API** | REST API (full access) |
| **Bulk Download** | FTP site |
| **Data Formats** | FASTA, XML, JSON, RDF, TSV |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 (MVP) - for target mapping |
| **Storage Estimate** | ~100 GB (full) |

---

## Comparison Matrix

### Data Coverage

| Database | Drugs | Targets | Interactions | SNP-Drug | Dosing |
|----------|-------|---------|--------------|----------|--------|
| PharmGKB | 1,000+ | 700+ genes | 20,000+ annotations | YES | YES |
| CPIC | 164 | 34 genes | - | YES | YES |
| DrugBank | 15,000+ | 5,000+ | 2.5M DDI | YES | Partial |
| ChEMBL | 2.8M compounds | 17,803 | 21M bioactivities | No | No |
| DGIdb | 20,000+ | 10,000+ | 70,000+ | No | No |
| Open Targets | Integrated | Genome-wide | Evidence scores | Genetic | No |
| BindingDB | 1.3M | Thousands | 2.9M affinity | No | No |
| GtoPdb | 13,131 | 3,097 | Quantitative | No | No |
| TTD | Ligands | 3,131 | Target-drug | No | No |

### Access and Licensing

| Database | API | Bulk Download | Academic | Commercial | Format |
|----------|-----|---------------|----------|------------|--------|
| PharmGKB | REST/JSON | Yes (license) | Free | Free (CC BY-SA) | TSV, JSON |
| CPIC | Via PharmGKB | Yes | Free (CC0) | Free (CC0) | PDF, JSON |
| DrugBank | REST | Yes | Free | Paid | XML, CSV |
| ChEMBL | REST | Yes | Free | Free | SQL, SDF |
| PubChem | PUG REST | Yes | Free | Free | SDF, JSON |
| DGIdb | GraphQL | Yes | Free | Free | TSV, JSON |
| Open Targets | GraphQL | Yes | Free | Free | Parquet |
| OpenFDA | REST | Yes | Free | Free | JSON |
| BindingDB | REST | Yes | Free | Free | SDF, TSV |
| GtoPdb | REST | Yes | Free | Free | CSV, PostgreSQL |

### Pharmacogenomics Focus

| Database | SNP-Drug Pairs | Dosing Guidelines | Metabolizer Status | Clinical Implementation |
|----------|----------------|-------------------|--------------------|-----------------------|
| PharmGKB | PRIMARY | Yes (CPIC/DPWG) | Yes | Yes |
| CPIC | Yes | PRIMARY | Yes | Yes (EHR ready) |
| DPWG | Yes | Yes | Yes | Yes (Dutch systems) |
| PharmVar | Allele definitions | No | Function tables | Nomenclature |
| DrugBank | Partial | No | Enzyme info | Partial |
| FDA PGx Table | Biomarkers | In labels | No | FDA guidance |
| CIViC | Cancer variants | Therapeutic | No | Cancer treatment |

---

## Integration Recommendations

### Tier 1: Essential (MVP)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| PharmGKB/ClinPGx | PGx | Gold standard SNP-drug data | CC BY-SA 4.0 |
| CPIC | PGx | Clinical dosing guidelines | CC0 |
| DPWG | PGx | European guidelines | Open |
| PharmVar | PGx | Allele nomenclature | Open |
| DrugBank | Drug Info | Comprehensive drug profiles | Academic free |
| OpenFDA | Regulatory | Adverse events, labels | Public Domain |
| DailyMed | Regulatory | FDA drug labels | Public Domain |
| FDA PGx Table | Regulatory | Biomarker list | Public Domain |
| UniProt | Reference | Drug target proteins | CC BY 4.0 |

### Tier 2: Important (Post-MVP)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| ChEMBL | Bioactivity | 21M+ activity measurements | Open |
| DGIdb | Interactions | Aggregated from 40+ sources | Open |
| Open Targets | Target-Disease | Genetic evidence integration | CC0 |
| BindingDB | Binding | Quantitative affinity data | Open |
| GtoPdb | Pharmacology | Expert-curated targets | Open |
| CIViC | Oncology | Cancer variant interpretation | CC0 |
| RxNorm | Terminology | Drug name standardization | Free subset |
| DrugCentral | Drug Info | FAERS integration | Open |

### Tier 3: Supplementary (Future)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| PubChem | Chemistry | 118M+ compounds | Public Domain |
| TTD | Targets | Druggability assessment | Free |
| PGxDB | PGx | Integrated platform | Open |
| SuperDRUG2 | Drug Info | PK simulation | Open |
| KEGG Drug | Pathways | Drug-pathway links | Academic/Paid |

---

## Data Pipeline Architecture

```
User Query (Gene/Drug/Variant)
         |
         v
+------------------+
| Normalize IDs    | <-- RxNorm, PharmVar, dbSNP
+------------------+
         |
         v
+------------------+
| PharmGKB/CPIC    | --> SNP-Drug associations, Dosing
+------------------+
         |
         v
+------------------+
| DrugBank/ChEMBL  | --> Drug info, Targets, Mechanisms
+------------------+
         |
         v
+------------------+
| DGIdb/Open Targets| --> Drug-Gene interactions
+------------------+
         |
         v
+------------------+
| OpenFDA/DailyMed | --> Labels, Adverse events
+------------------+
         |
         v
    Integrated Response
```

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Tier 1 databases | ~260 GB (PharmGKB, CPIC, DPWG, PharmVar, DrugBank, OpenFDA, DailyMed, FDA PGx, UniProt) |
| Tier 2 databases | ~90 GB (ChEMBL, DGIdb, Open Targets, BindingDB, GtoPdb, CIViC, RxNorm, DrugCentral) |
| Tier 3 databases | ~510 GB (PubChem, TTD, PGxDB, SuperDRUG2, KEGG Drug) |
| Total storage estimate | ~860 GB |
| Last updated | January 2026 |

*Note: PubChem dominates Tier 3 at ~500 GB; selective compound import recommended.*

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [../_index.md](../_index.md) | Parent index |
| [../../genetics/primary.md](../../genetics/primary.md) | dbSNP/ClinVar variant data |
| [../../pathways/primary.md](../../pathways/primary.md) | Drug pathway integration |

---

## License

This document catalogs multiple databases with varying license terms:

| Database | License | Commercial Use | Attribution | Access |
|----------|---------|----------------|-------------|--------|
| PharmGKB/ClinPGx | CC BY-SA 4.0 | Yes (with attribution) | Required | Free (account required) |
| CPIC | CC0 (Public Domain) | Yes | None required | Open |
| DPWG | Open (published journals) | Yes | Citation | Via PharmGKB |
| PharmVar | Open Access | Yes | Citation | Open |
| DrugBank | CC0 (vocabulary), CC BY-NC 4.0 (academic), Commercial license | Requires license | Required | Academic free |
| ChEMBL | Open Access (Apache 2.0 for software) | Yes | Citation | Open |
| PubChem | Public Domain | Yes | None required | Open |
| DrugCentral | Open Access | Yes | Citation | Open |
| Open Targets | CC0 (data), Open Source (tools) | Yes | Citation | Open |
| DGIdb | Open Access | Yes | Citation | Open |
| TTD | Free Access | Yes | Citation | Open |
| BindingDB | Open Access (FAIR compliant) | Yes | Citation | Open |
| GtoPdb | Open Access | Yes | Citation | Open |
| DailyMed | Public Domain | Yes | None required | Open |
| OpenFDA | Public Domain | Yes | None required | Open |
| FDA PGx Table | Public Domain | Yes | None required | Open |
| RxNorm | Free (prescribable subset) | Yes (subset) | Citation | Open (subset) |
| CIViC | CC0 (data), MIT (software) | Yes | None required | Open |
| PGxDB | Open Access (FAIR compliant) | Yes | Citation | Open |
| SuperDRUG2 | Open Access | Yes | Citation | Open |
| KEGG Drug | Academic API only, Subscription for FTP | No (FTP), Yes (Medicus) | Required | Subscription (full) |
| UniProt | CC BY 4.0 | Yes | Required | Open |

**Key Considerations:**
- **Fully Open (Commercial OK):** CPIC, PubChem, OpenFDA, DailyMed, FDA PGx Table, CIViC, Open Targets
- **Commercial Friendly with Attribution:** PharmGKB (CC BY-SA), UniProt (CC BY), ChEMBL
- **Academic Only:** DrugBank (full data), KEGG Drug (FTP)
- **Public Domain:** CPIC, PubChem, DailyMed, OpenFDA, CIViC

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **ChEMBL** | FTP/API | `ftp://ftp.ebi.ac.uk/pub/databases/chembl/` |
| **PubChem** | FTP/API | `ftp://ftp.ncbi.nlm.nih.gov/pubchem/` |
| **DrugBank** | Download | `https://go.drugbank.com/releases/latest` (registration required) |
| **PharmGKB** | Bulk download | `https://www.pharmgkb.org/downloads` |
| **OpenFDA** | API | `https://open.fda.gov/apis/` |
| **DailyMed** | Download | `https://dailymed.nlm.nih.gov/dailymed/spl-resources.cfm` |
| **BindingDB** | Download | `https://www.bindingdb.org/bind/chemsearch/marvin/SDFdownload.jsp` |

**Access Requirements:** Most are freely accessible; DrugBank full data requires academic license; OpenFDA and DailyMed are public domain.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | SDF, TSV, JSON |
| Alternative | CSV, XML, RDF |
| Chemical structures | SMILES, InChI, MOL |
| Identifiers | DrugBank ID, ChEMBL ID, PubChem CID |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `drug_id` | string | Primary identifier | "CHEMBL25" |
| `name` | string | Generic drug name | "Aspirin" |
| `mechanism` | string | Mechanism of action | "COX inhibitor" |
| `indication` | array | Approved indications | ["pain", "inflammation"] |
| `targets` | array | Protein targets | ["PTGS1", "PTGS2"] |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `targets` | Protein | N:M |
| `has_indication` | Disease | N:M |
| `interacts_with` | Drug | N:M |

## Sample Data

### Example Drug Record
```json
{
  "chembl_id": "CHEMBL25",
  "drugbank_id": "DB00945",
  "name": "Aspirin",
  "synonyms": ["acetylsalicylic acid", "ASA"],
  "mechanism": "Irreversible COX-1/COX-2 inhibitor",
  "targets": [
    {"uniprot_id": "P23219", "gene": "PTGS1", "action": "inhibitor"},
    {"uniprot_id": "P35354", "gene": "PTGS2", "action": "inhibitor"}
  ],
  "max_phase": 4
}
```

### Sample Query Result
| drug | target | Ki_nM | indication |
|------|--------|-------|------------|
| Aspirin | PTGS1 | 1.67 | Pain/inflammation |
| Ibuprofen | PTGS2 | 6000 | Pain/inflammation |

## Data Set Size

| Metric | Value |
|--------|-------|
| ChEMBL compounds | 2.4M+ compounds |
| PubChem compounds | 115M+ compounds |
| DrugBank drugs | 15K+ approved/experimental |
| BindingDB interactions | 2.9M+ binding measurements |
| Total storage estimate | ~50-100 GB (combined sources) |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `pharmacogenomics (PGx)` | The study of how genetic variation affects drug response and metabolism | CYP2D6 poor metabolizer status affecting codeine efficacy |
| `star allele` | Standardized nomenclature for pharmacogene haplotypes using asterisk notation | CYP2D6*4 (non-functional allele) |
| `diplotype` | The combination of two haplotypes (star alleles) an individual carries | CYP2D6*1/*4 |
| `metabolizer phenotype` | Predicted drug metabolism capacity based on diplotype | Poor, Intermediate, Normal, Ultrarapid metabolizer |
| `IC50` | Half-maximal inhibitory concentration - potency measure for inhibitors | 10 nM IC50 indicates high potency |
| `Ki` | Inhibition constant - binding affinity of an inhibitor for its target | Lower Ki = tighter binding |
| `drug-drug interaction (DDI)` | When one drug affects the pharmacokinetics or pharmacodynamics of another | CYP3A4 inhibitor increasing statin levels |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `PharmGKB` | Pharmacogenomics Knowledge Base with curated SNP-drug associations | Clinical PGx |
| `CPIC` | Clinical Pharmacogenetics Implementation Consortium providing dosing guidelines | Dosing recommendations |
| `DPWG` | Dutch Pharmacogenetics Working Group with European guidelines | European PGx |
| `PharmVar` | Pharmacogene Variation Consortium defining official star allele nomenclature | Allele naming |
| `DrugBank` | Comprehensive drug database with targets, interactions, and PK data | Drug information |
| `ChEMBL` | EMBL-EBI bioactivity database with 21M+ activity measurements | Bioactivity data |
| `DGIdb` | Drug Gene Interaction Database aggregating 40+ sources | Drug-gene interactions |
| `Open Targets` | Platform integrating genetic evidence for drug target validation | Target identification |
| `FAERS` | FDA Adverse Event Reporting System for post-market drug safety | Adverse events |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| PGx | Pharmacogenomics | Genetic variation affecting drug response |
| CYP | Cytochrome P450 | Major drug-metabolizing enzyme family |
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity | Pharmacokinetic properties |
| DDI | Drug-Drug Interaction | Clinically significant drug combinations |
| SNP | Single Nucleotide Polymorphism | Common genetic variant type |
| VIP | Very Important Pharmacogene | PharmGKB-designated key PGx genes |
| EHR | Electronic Health Record | Clinical decision support integration target |
| FDA | Food and Drug Administration | US drug regulatory agency |
| NDC | National Drug Code | FDA drug identification system |
| SPL | Structured Product Labeling | FDA drug label XML format |
| CC0 | Creative Commons Zero | Public domain dedication license |
| CC BY | Creative Commons Attribution | Open license requiring attribution |

---

## Open Questions

- [ ] DrugBank commercial license - required for production deployment?
- [ ] KEGG licensing - academic use sufficient for platform?
- [ ] PubChem subset - which compounds to prioritize (size concerns)?
- [ ] Update automation - pipeline for monthly PharmGKB/CPIC updates?
- [ ] EHR integration - CPIC guideline format compatibility?

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial comprehensive catalog from interventions-pharma.md |
