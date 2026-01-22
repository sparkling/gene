# Drug Metabolism and Supplement Interaction Data Sources

**Document ID:** 43-53-DRUG-METABOLISM
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-00-INDEX](./../INDEX.md)

---

## TL;DR

Drug metabolism databases provide comprehensive coverage of CYP450 enzyme polymorphisms, drug-enzyme interactions, supplement-drug interactions, and transporter genetics. PharmVar serves as the gold standard for star allele nomenclature, while PharmGKB and CPIC provide clinical implementation guidelines. Integration of 25+ databases enables complete pharmacokinetic profiling from genotype through metabolizer status to personalized dosing recommendations.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary CYP allele nomenclature | PharmVar | Official star allele definitions used by PharmGKB/CPIC | Jan 2026 |
| Clinical PGx guidelines | PharmGKB + CPIC | Gold standard with EHR-ready implementation | Jan 2026 |
| Drug pathway integration | KEGG + DrugBank | Comprehensive pathway coverage with REST APIs | Jan 2026 |
| Supplement interactions | NatMed Pro + MSKCC | Industry-leading accuracy with free alternative | Jan 2026 |
| Transporter genetics | PharmGKB + CPIC | Curated SLCO1B1, ABCB1, ABCG2 annotations | Jan 2026 |
| PK modeling data | PK-DB | Full REST API with meta-analysis support | Jan 2026 |

---

## Database Catalog

### CYP450 Enzyme Databases

#### 1. PharmVar (Pharmacogene Variation Consortium)

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmvar.org/genes |
| **Maintainer** | PharmVar Consortium |
| **Content** | Star (*) allele nomenclature for 29 CYP enzymes plus NADPH cytochrome P450 oxidoreductase (POR). Functionally characterized variants mapped to RefSeqGene, LRG, GRCh37/hg19, GRCh38/hg38. SNPs linked to rs numbers/dbSNP identifiers. |
| **Records** | ~2,000+ SNPs and mutations across CYP genes |
| **License** | Free for academic/research use |
| **API** | None (download-based) |
| **Bulk Download** | FASTA sequence files per gene |
| **Data Formats** | FASTA, TSV |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- Official standard used by PharmGKB and CPIC
- Comprehensive coverage of CYP2D6, CYP2C9, CYP2C19, CYP2B6 (highly polymorphic genes)
- Transitioned from Human Cytochrome P450 Allele Nomenclature Database (2017)
- Cross-references with PharmGKB and dbSNP

**Schema:**
```
Gene --> Star Allele --> Defining Variants (rsIDs) --> Function Assignment
                              |
                              v
                    Reference Sequence Positions
```

---

#### 2. SuperCYP / SuperCYPsPred

| Field | Value |
|-------|-------|
| **URL** | https://insilico-cyp.charite.de/SuperCYPsPred/ |
| **Maintainer** | Charite - Universitatsmedizin Berlin |
| **Content** | 1,170 drugs with 3,800+ interactions. ~2,000 SNPs/mutations with effect on expression/activity. Covers CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4 (responsible for >90% clinical drug metabolism). Homology-modeled 3D structures. |
| **Records** | 1,170 drugs, 3,800+ interactions, ~2,000 SNPs |
| **License** | Free academic access |
| **API** | None (web interface) |
| **Bulk Download** | PDB files (CYP structures), MOL files (drug structures) |
| **Data Formats** | PDB, MOL, web interface |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~100 MB |

**Key Features:**
- Drug-cocktail mixing predictions
- SNP effect analysis
- Substrate/inhibitor/inducer classification

---

#### 3. Indiana University Flockhart CYP450 Table

| Field | Value |
|-------|-------|
| **URL** | https://drug-interactions.medicine.iu.edu/ |
| **Maintainer** | Indiana University |
| **Content** | Clinically relevant CYP-drug interactions organized by P450 isoform (8 columns). Lists substrates, inhibitors, and inducers. |
| **Records** | Hundreds of clinically significant interactions |
| **License** | Free academic/clinical use |
| **API** | None (mobile-friendly search interface) |
| **Bulk Download** | No |
| **Data Formats** | Web interface |
| **Update Frequency** | Twice yearly minimum |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~10 MB |

**Key Features:**
- Gold standard open-access CYP table
- Clinical focus
- Excellent teaching resource

---

#### 4. Curated CYP450 Interaction Dataset (2025)

| Field | Value |
|-------|-------|
| **URL** | https://www.nature.com/articles/s41597-025-05753-8 |
| **Maintainer** | Scientific Data (Nature) |
| **Content** | Substrates and non-substrates for 6 principal CYP450 isoforms: CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP3A4. Machine learning-ready dataset with GCN support. |
| **Records** | ~2,000 compounds per enzyme |
| **License** | Open access (Scientific Data journal) |
| **API** | Dataset download via Scientific Data repository |
| **Bulk Download** | Yes |
| **Data Formats** | CSV, ML-ready format |
| **Update Frequency** | Publication-based |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- Designed for ML/AI pharmacokinetic modeling
- GCN (Graph Convolutional Network) support
- Modern dataset for predictive model building

---

### Drug Metabolism Pathway Databases

#### 5. KEGG DRUG / KEGG PATHWAY

| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/drug/ |
| **Maintainer** | Kanehisa Laboratories |
| **Content** | Approved drugs from Japan, USA, Europe unified by chemical structure. KEGG DGROUP: Drug groups by metabolizing enzymes and transporters. 495 reference pathways across 4,700+ organisms. |
| **Records** | 11,000+ drugs, 495 reference pathways |
| **License** | Free for academic use; commercial requires license |
| **API** | KEGG REST API (https://rest.kegg.jp/) |
| **Bulk Download** | Subscription required for full FTP |
| **Data Formats** | KEGG flat file format |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~1 GB |

**Key Features:**
- Industry standard for pathway analysis
- KEGG Mapper for large-scale data integration
- Integrates with DrugBank for drug target downloads

---

#### 6. PharmGKB (Pharmacogenomics Knowledge Base)

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org/ |
| **Maintainer** | Stanford University (NIH-funded) |
| **Content** | 715+ drugs, 1,761 genes, 227 diseases, 165 clinical guidelines, 784 drug labels. 9,000+ manually curated literature annotations. 153 drug PK/PD pathways, 68 VIP gene summaries. 26,865 variant annotations. |
| **Records** | 775 drugs, 234 pathways, 26,865 variant annotations |
| **License** | Free for research (registration required) |
| **API** | REST web services |
| **Bulk Download** | Zipped spreadsheets (license agreement required) |
| **Data Formats** | TSV/CSV, JSON, BioPax |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB |

**Evidence Levels:**
- **1A:** CPIC/medical society guideline or implemented clinically
- **1B:** Preponderance of evidence shows association
- **2A:** Level 2B + Very Important Pharmacogene (VIP)
- **2B:** Moderate evidence, replicated
- **3:** Single significant study
- **4:** Case report, in vitro evidence only

---

#### 7. DrugBank

| Field | Value |
|-------|-------|
| **URL** | https://go.drugbank.com/ |
| **Maintainer** | DrugBank (OMx Personal Health Analytics) |
| **Content** | 2,832 drugs, 800 drug metabolites. Comprehensive drug-target, drug-enzyme, drug-transporter data. |
| **Records** | 2,832 drugs, 800 metabolites |
| **License** | CC BY-NC 4.0 (academic); commercial requires license |
| **API** | REST API (Clinical API for commercial) |
| **Bulk Download** | XML, CSV, JSON |
| **Data Formats** | XML, CSV, SDF, FASTA |
| **Update Frequency** | Regular |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~5 GB |

**Key Features:**
- Maps to MeSH, ATC, RxNorm
- Part of HMDB suite
- Custom exports available with license

---

#### 8. CPIC (Clinical Pharmacogenetics Implementation Consortium)

| Field | Value |
|-------|-------|
| **URL** | https://cpicpgx.org/ |
| **Maintainer** | PharmGKB/PGRN consortium |
| **Content** | Clinical practice guidelines for gene-drug pairs. Diplotype-to-phenotype mappings, dosing recommendations, CDS language. |
| **Records** | 34 genes, 164 drugs, 27 published guidelines |
| **License** | CC0 (Public Domain) |
| **API** | RESTful API (public access, 80,000+ monthly queries) |
| **Bulk Download** | JSON via PharmGKB, GitHub exports |
| **Data Formats** | PDF, Excel/TSV, JSON |
| **Update Frequency** | As guidelines published |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |

**Key Tables:**
- Allele Definition Table
- Allele Functionality Table
- Diplotype-Phenotype Table
- Dosing Recommendation Table

---

### Supplement-Drug Interaction Databases

#### 9. NatMed Pro (Natural Medicines Database)

| Field | Value |
|-------|-------|
| **URL** | https://naturalmedicines.therapeuticresearch.com/ |
| **Maintainer** | TRC Healthcare |
| **Content** | 1,400+ vitamins, minerals, botanicals, non-botanical supplements, foods. Drug-nutrient interactions, mechanism of action, safety data. |
| **Records** | 1,400+ ingredients |
| **License** | Subscription-based |
| **API** | RESTful JSON:API (licensed users) |
| **Bulk Download** | No (API access) |
| **Data Formats** | JSON |
| **Update Frequency** | Daily evidence updates |
| **Priority** | Tier 1 (if licensed) |
| **Storage Estimate** | ~500 MB |

**Key Features:**
- Industry-leading supplement database
- Validated for accuracy (2024 JCO Oncology Advances study)
- Multi-ingredient product checking

---

#### 10. MSKCC About Herbs

| Field | Value |
|-------|-------|
| **URL** | https://www.mskcc.org/cancer-care/diagnosis-treatment/symptom-management/integrative-medicine/herbs |
| **Maintainer** | Memorial Sloan Kettering Cancer Center |
| **Content** | ~284 dietary supplement monographs. Purported uses, adverse effects, herb-drug interactions. |
| **Records** | ~284 supplement monographs |
| **License** | Free public access |
| **API** | None (web access and mobile app) |
| **Bulk Download** | No |
| **Data Formats** | Web interface |
| **Update Frequency** | Continual |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- Oncology focus but broadly applicable
- Mobile app available
- Public and healthcare professional versions

---

#### 11. Stockley's Drug Interactions

| Field | Value |
|-------|-------|
| **URL** | https://www.medicinescomplete.com/mc/stockley/current/ |
| **Maintainer** | Pharmaceutical Press |
| **Content** | Comprehensive drug-drug and drug-herb interactions. 5 severity classifications from "avoid" to "no action needed." |
| **Records** | Thousands of interaction monographs |
| **License** | Commercial subscription |
| **API** | None (web access) |
| **Bulk Download** | No |
| **Data Formats** | Web interface |
| **Update Frequency** | Regular |
| **Priority** | Tier 3 |
| **Storage Estimate** | N/A (subscription) |

---

#### 12. FDB (First Databank) MedKnowledge

| Field | Value |
|-------|-------|
| **URL** | https://www.fdbhealth.com/solutions/medknowledge-drug-database/ |
| **Maintainer** | First Databank |
| **Content** | Drug-drug, drug-alternative therapy interactions. Three severity levels plus clinical effects subcategories. |
| **Records** | Comprehensive (enterprise-grade) |
| **License** | Commercial license required |
| **API** | Commercial API integration |
| **Bulk Download** | Enterprise agreement |
| **Data Formats** | Multiple (EHR integration) |
| **Update Frequency** | Regular |
| **Priority** | Tier 3 |
| **Storage Estimate** | N/A (commercial) |

---

### Natural Product-Drug Interaction Resources

#### 13. IMgateway (University of Sydney)

| Field | Value |
|-------|-------|
| **URL** | https://www.imgateway.net/ |
| **Maintainer** | University of Sydney |
| **Content** | 1,000+ drug interactions with foods, herbs, supplements. Evidence-based knowledge base. |
| **Records** | 1,000+ interactions |
| **License** | Partnership/license required |
| **API** | Available for healthcare application integration |
| **Bulk Download** | License required |
| **Data Formats** | API-based |
| **Update Frequency** | Regular |
| **Priority** | Tier 3 |
| **Storage Estimate** | ~100 MB |

---

#### 14. T3DB (Toxin and Toxin Target Database)

| Field | Value |
|-------|-------|
| **URL** | https://www.t3db.ca/ |
| **Maintainer** | HMDB suite |
| **Content** | 3,678 toxins (pollutants, pesticides, drugs, food toxins) with 42,374 toxin-target associations. 2,073 toxin target records. 90+ data fields per ToxCard. |
| **Records** | 3,678 toxins, 2,073 targets, 42,374 associations |
| **License** | Free access |
| **API** | Web interface with download options |
| **Bulk Download** | Yes |
| **Data Formats** | Web download |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~200 MB |

**Key Features:**
- Focus on toxicity mechanisms and target proteins
- Links to HMDB and DrugBank
- Useful for drug-toxin interaction prediction

---

#### 15. HMDB (Human Metabolome Database)

| Field | Value |
|-------|-------|
| **URL** | https://hmdb.ca/ |
| **Maintainer** | University of Alberta |
| **Content** | 220,945 metabolite entries (water-soluble and lipid-soluble). 8,610 protein sequences (enzymes and transporters). Spectral data (NMR, MS/MS, GC-MS). |
| **Records** | 220,945 metabolites, 8,610 proteins |
| **License** | CC BY-NC 4.0 |
| **API** | Contact required (academic: eponine@ualberta.ca) |
| **Bulk Download** | Yes (https://hmdb.ca/downloads) |
| **Data Formats** | XML, SDF, CSV |
| **Update Frequency** | Periodic |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~10 GB |

**Key Features:**
- Core metabolomics resource
- Links to DrugBank, T3DB, FooDB
- Essential for understanding drug metabolites

---

### Absorption/Bioavailability Databases

#### 16. PK-DB (Pharmacokinetics Database)

| Field | Value |
|-------|-------|
| **URL** | https://pk-db.com/ |
| **Maintainer** | Academic consortium |
| **Content** | Pharmacokinetic data from multiple studies. Experimental errors, normalized units, biological ontology annotations. |
| **Records** | Multi-study aggregated PK data |
| **License** | Free academic access |
| **API** | Full REST API (https://pk-db.com/api/v1/swagger/) |
| **Bulk Download** | Via API |
| **Data Formats** | JSON |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~500 MB |

**Key Features:**
- Designed for meta-analysis
- PBPK/PK-PD/population PK modeling support
- Collaborative data curation workflow

---

#### 17. SwissADME

| Field | Value |
|-------|-------|
| **URL** | http://www.swissadme.ch |
| **Maintainer** | Swiss Institute of Bioinformatics |
| **Content** | Predictive models for physicochemical properties, pharmacokinetics, drug-likeness, medicinal chemistry friendliness. |
| **Records** | Predictive tool (not static database) |
| **License** | Free academic/research use |
| **API** | None (web tool) |
| **Bulk Download** | No |
| **Data Formats** | Web interface |
| **Update Frequency** | Ongoing |
| **Priority** | Tier 3 |
| **Storage Estimate** | N/A (web tool) |

---

#### 18. FDA BCS (Biopharmaceutics Classification System)

| Field | Value |
|-------|-------|
| **URL** | https://www.fda.gov/media/148472/download |
| **Maintainer** | FDA |
| **Content** | Drug classification by solubility and permeability. Class I-IV classifications. Biowaiver guidance. |
| **Records** | Classification framework + drug examples |
| **License** | Public Domain |
| **API** | None (PDF guidance) |
| **Bulk Download** | PDF |
| **Data Formats** | PDF |
| **Update Frequency** | Periodic |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~5 MB |

---

#### 19. FDA Table of Pharmacogenomic Biomarkers

| Field | Value |
|-------|-------|
| **URL** | https://www.fda.gov/media/124784/download |
| **Maintainer** | FDA |
| **Content** | Approved products with PGx information in drug labeling. Biomarkers per HUGO nomenclature. |
| **Records** | 300+ drug labels with PGx biomarkers |
| **License** | Public Domain |
| **API** | PharmGKB parser: https://github.com/PharmGKB/fda-biomarker |
| **Bulk Download** | PDF, parseable |
| **Data Formats** | PDF |
| **Update Frequency** | Regular (June 2024) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~10 MB |

---

### Transporter Genetics Databases

#### 20. PharmGKB Transporter Annotations

| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org/ |
| **Maintainer** | Stanford University |
| **Content** | SLCO1B1 (OATP1B1), ABCB1 (MDR1/P-gp), ABCG2 (BCRP), OCT1, OAT1, OAT3, MATE1, MATE2-K genetic variants. Clinical annotations for transporter polymorphisms. |
| **Records** | Comprehensive transporter pharmacogenomics data |
| **License** | Free research use (registration required) |
| **API** | REST web services |
| **Bulk Download** | Via PharmGKB downloads |
| **Data Formats** | TSV, JSON |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | Included in PharmGKB |

**Key Features:**
- Best curated resource for transporter PGx
- SLCO1B1 521T>C and 388A>G extensively characterized
- Statin-SLCO1B1 associations documented

---

#### 21. CPIC Guidelines for Transporters

| Field | Value |
|-------|-------|
| **URL** | https://cpicpgx.org/guidelines/ |
| **Maintainer** | CPIC |
| **Content** | Clinical guidelines for SLCO1B1 and statin dosing. Transporter genotype-to-phenotype translations. |
| **Records** | Guideline-specific |
| **License** | CC0 (Public Domain) |
| **API** | CPIC API |
| **Bulk Download** | Via CPIC |
| **Data Formats** | PDF, JSON |
| **Update Frequency** | As guidelines published |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | Included in CPIC |

**Key Features:**
- SLCO1B1 genotyping for statin-induced myopathy prevention
- Model case for transporter pharmacogenetics

---

### Key Transporter Polymorphisms Reference

| Transporter | Gene | Key Variants | Clinical Relevance |
|-------------|------|--------------|-------------------|
| P-glycoprotein | ABCB1 (MDR1) | 1236C>T (rs1128503), 2677G>T/A (rs2032582), 3435C>T (rs1045642) | Digoxin, fexofenadine, cyclosporine PK |
| OATP1B1 | SLCO1B1 | 521T>C (V174A, *5), 388A>G (N130D, *1b), *15 (double) | Statin-induced myopathy, methotrexate clearance |
| OATP1B3 | SLCO1B3 | Multiple SNPs | Rifampin, telmisartan PK |
| BCRP | ABCG2 | 421C>A (Q141K) | Rosuvastatin, sulfasalazine PK |
| OCT1 | SLC22A1 | Multiple reduced-function variants | Metformin response |
| OAT1/OAT3 | SLC22A6/SLC22A8 | Various | Renal drug elimination |
| MATE1/MATE2-K | SLC47A1/SLC47A2 | Various | Metformin elimination |

---

## Comparison Matrix

### Data Coverage by Category

| Database | CYP Variants | Drug-Drug | Supplement | Transporter | PK Data |
|----------|--------------|-----------|------------|-------------|---------|
| PharmVar | PRIMARY | No | No | No | No |
| SuperCYP | 2,000 SNPs | 3,800+ | No | No | No |
| PharmGKB | Via alleles | Via pathways | No | YES | YES |
| CPIC | Via guidelines | No | No | YES | YES |
| DrugBank | Partial | 2.5M | No | YES | YES |
| KEGG | No | Pathway | No | YES | YES |
| NatMed Pro | No | No | PRIMARY | No | No |
| MSKCC | No | No | 284 | No | No |
| T3DB | No | No | Partial | Partial | No |
| HMDB | No | No | Metabolites | YES | No |
| PK-DB | No | No | No | No | PRIMARY |

### Access and Licensing

| Database | API | Bulk Download | Academic | Commercial | Priority |
|----------|-----|---------------|----------|------------|----------|
| PharmVar | No | Yes | Free | Free | Tier 1 |
| PharmGKB | REST | Yes | Free | CC BY-SA | Tier 1 |
| CPIC | REST | Yes | Free (CC0) | Free (CC0) | Tier 1 |
| DrugBank | REST | Yes | Free | Paid | Tier 1 |
| KEGG | REST | Partial | Free | Paid | Tier 1 |
| NatMed Pro | REST | No | Subscription | Subscription | Tier 1* |
| MSKCC | No | No | Free | Free | Tier 2 |
| HMDB | Contact | Yes | Free | Contact | Tier 1 |
| PK-DB | REST | Yes | Free | Free | Tier 2 |
| SuperCYP | No | Partial | Free | Free | Tier 2 |

---

## Integration Priority

### Tier 1: Essential (MVP)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| PharmVar | CYP Nomenclature | Official star allele definitions | Free |
| PharmGKB | PGx | Clinical annotations, pathways | CC BY-SA 4.0 |
| CPIC | Clinical | EHR-ready dosing guidelines | CC0 |
| DrugBank | Drug Metabolism | Drug-enzyme relationships | Academic free |
| KEGG | Pathways | Drug metabolism pathways | Academic free |
| HMDB | Metabolites | Drug metabolite data | CC BY-NC 4.0 |
| FDA PGx Table | Regulatory | Biomarker labeling | Public Domain |

### Tier 2: Important (Post-MVP)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| SuperCYP | CYP Interactions | Detailed CYP-drug data | Free |
| Flockhart Table | Clinical Reference | Quick clinical lookup | Free |
| MSKCC About Herbs | Supplements | Free herb-drug interactions | Free |
| T3DB | Toxicity | Toxin mechanisms | Free |
| PK-DB | Pharmacokinetics | PK modeling data | Free |
| FDA BCS | Bioavailability | Classification framework | Public Domain |

### Tier 3: Supplementary (Future)

| Source | Category | Key Value | License |
|--------|----------|-----------|---------|
| NatMed Pro | Supplements | Industry-leading accuracy | Subscription |
| Stockley's | Interactions | Comprehensive reference | Subscription |
| FDB MedKnowledge | EHR Integration | Enterprise DDI | Commercial |
| IMgateway | Natural Products | Academic research | License |
| SwissADME | ADME Prediction | Computational predictions | Free |

---

## Data Pipeline Architecture

```
User Query (Drug/Gene/Supplement)
         |
         v
+------------------+
| Normalize IDs    | <-- PharmVar (*alleles), RxNorm, dbSNP
+------------------+
         |
         v
+------------------+
| CYP Status       | --> PharmVar + CPIC (metabolizer phenotype)
+------------------+
         |
         v
+------------------+
| Drug Metabolism  | --> DrugBank + KEGG (enzymes, pathways)
+------------------+
         |
         v
+------------------+
| Interactions     | --> NatMed Pro, MSKCC (supplement-drug)
+------------------+
         |
         v
+------------------+
| Transporter PGx  | --> PharmGKB (SLCO1B1, ABCB1, ABCG2)
+------------------+
         |
         v
+------------------+
| PK Parameters    | --> PK-DB, HMDB (pharmacokinetics)
+------------------+
         |
         v
    Personalized Dosing Recommendation
```

---

## Storage Estimates Summary

| Tier | Databases | Estimated Storage |
|------|-----------|-------------------|
| Tier 1 | PharmVar, PharmGKB, CPIC, DrugBank, KEGG, HMDB, FDA PGx | ~19 GB |
| Tier 2 | SuperCYP, Flockhart, MSKCC, T3DB, PK-DB, FDA BCS | ~1 GB |
| Tier 3 | NatMed Pro, Stockley's, FDB, IMgateway, SwissADME | ~1 GB (subscriptions) |
| **Total** | All | **~21 GB** |

---

## Key Identifiers for Cross-Linking

| Identifier | Used By | Purpose |
|------------|---------|---------|
| rsID (dbSNP) | PharmVar, PharmGKB, CPIC | Variant identification |
| HGNC Symbol | All databases | Gene identification |
| Star Allele (*) | PharmVar, PharmGKB, CPIC | Haplotype naming |
| DrugBank ID | DrugBank, HMDB, Reactome | Drug cross-reference |
| KEGG D number | KEGG | Drug identification |
| ATC Code | DrugBank, KEGG | Drug classification |
| RxNorm | DrugBank, CPIC | Clinical drug vocabulary |
| CID (PubChem) | Multiple | Chemical structure |
| InChIKey | Multiple | Chemical identifier |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-00-INDEX](./../INDEX.md) | Parent index |
| [43-51-PHARMACEUTICALS](./../compounds/PHARMACEUTICALS.md) | Drug-gene core data |
| [43-11-GENETICS-PRIMARY](./../genetics/PRIMARY.md) | dbSNP/ClinVar variant data |
| [43-41-PATHWAYS-PRIMARY](./../pathways/PRIMARY.md) | Drug pathway integration |
| [43-81-BIOMARKERS-LABS](./../clinical/BIOMARKERS-LABS.md) | Lab test integration |
| [44-ARCHITECTURE](../44-ARCHITECTURE.md) | Database design |
| [45-DATA-MODEL](../45-DATA-MODEL.md) | Entity structure |

---

## Open Questions

- [ ] NatMed Pro licensing - subscription cost for production?
- [ ] PharmVar API - future development plans?
- [ ] KEGG commercial licensing - required for production?
- [ ] Transporter guidelines - CPIC expansion timeline?
- [ ] PK-DB coverage - sufficient for all major drugs?

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Initial catalog from data-sources-drug-metabolism.md |
