# Drug Metabolism and Supplement Interaction Data Sources

Research compiled for the Gene Platform to support pharmacogenomics, drug metabolism prediction, and supplement-drug interaction analysis.

---

## Table of Contents

1. [CYP450 Enzyme Databases](#1-cyp450-enzyme-databases)
2. [Drug Metabolism Pathway Databases](#2-drug-metabolism-pathway-databases)
3. [Supplement-Drug Interaction Databases](#3-supplement-drug-interaction-databases)
4. [Natural Product-Drug Interaction Resources](#4-natural-product-drug-interaction-resources)
5. [Absorption/Bioavailability Databases](#5-absorptionbioavailability-databases)
6. [Transporter Genetics Databases (OATP, P-gp)](#6-transporter-genetics-databases-oatp-p-gp)
7. [Summary Comparison Table](#7-summary-comparison-table)
8. [Integration Recommendations](#8-integration-recommendations)

---

## 1. CYP450 Enzyme Databases

### 1.1 PharmVar (Pharmacogene Variation Consortium)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pharmvar.org/genes |
| **Content** | Star (*) allele nomenclature for 29 CYP enzymes plus NADPH cytochrome P450 oxidoreductase (POR). Contains functionally characterized variants mapped to RefSeqGene, LRG, GRCh37/hg19, and GRCh38/hg38. SNPs linked to rs numbers/dbSNP identifiers. |
| **API** | No public REST API. Data available via web interface and downloadable sequence files (FASTA). Cross-references with PharmGKB and dbSNP. |
| **License** | Free for academic/research use. Star allele definitions are the official standard used by PharmGKB and CPIC. |
| **Size** | ~2000+ SNPs and mutations across CYP genes. Comprehensive coverage of CYP2D6, CYP2C9, CYP2C19, CYP2B6 (highly polymorphic genes). |
| **Notes** | Gold standard for pharmacogene nomenclature. Works closely with PharmGKB and CPIC for data curation. Transitioned from the Human Cytochrome P450 Allele Nomenclature Database in 2017. |

### 1.2 SuperCYP / SuperCYPsPred

| Attribute | Details |
|-----------|---------|
| **URL** | https://insilico-cyp.charite.de/SuperCYPsPred/ (Original: http://bioinformatics.charite.de/supercyp) |
| **Content** | 1,170 drugs with 3,800+ interactions. ~2,000 SNPs/mutations ordered by effect on expression/activity. Covers CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4 (responsible for >90% clinical drug metabolism). Includes homology-modeled 3D structures. |
| **API** | No REST API. Web interface for queries. PDB files downloadable for CYP structures, MOL files for drug structures. |
| **License** | Free academic access. Maintained by Charite - Universitatsmedizin Berlin. |
| **Size** | 1,170 drugs, 3,800+ interactions, ~2,000 SNPs |
| **Notes** | Useful for drug-cocktail mixing predictions and SNP effect analysis. Includes substrate/inhibitor/inducer classification. |

### 1.3 Indiana University Flockhart CYP450 Table

| Attribute | Details |
|-----------|---------|
| **URL** | https://drug-interactions.medicine.iu.edu/ |
| **Content** | Clinically relevant CYP-drug interactions organized by P450 isoform (8 columns). Lists substrates, inhibitors, and inducers. Considered the "gold standard" open-access CYP table. |
| **API** | No API. Mobile-friendly search interface. |
| **License** | Free academic/clinical use. Citation: Flockhart DA, Thacker D, McDonald C, Desta Z. The Flockhart Cytochrome P450 Drug-Drug Interaction Table. |
| **Size** | Curated selection of clinically significant interactions (hundreds of drugs) |
| **Notes** | Updated at least twice yearly. Focus on clinical relevance. Excellent teaching resource. |

### 1.4 P450Rdb

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.sciencedirect.com/science/article/pii/S2090123223003168 |
| **Content** | Manually curated database of reactions catalyzed by cytochrome P450 enzymes. Covers P450 reactions across all kingdoms of life. |
| **API** | Not specified in available documentation |
| **License** | Academic publication - data availability per journal requirements |
| **Size** | Not specified |
| **Notes** | Focus on catalytic reactions rather than drug interactions specifically. |

### 1.5 Curated CYP450 Interaction Dataset (2025)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.nature.com/articles/s41597-025-05753-8 (Scientific Data) |
| **Content** | Substrates and non-substrates for 6 principal CYP450 isoforms: CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP3A4. Machine learning-ready dataset with GCN (Graph Convolutional Network) support. |
| **API** | Dataset download via Scientific Data repository |
| **License** | Open access (Scientific Data journal) |
| **Size** | ~2,000 compounds per enzyme |
| **Notes** | Modern dataset designed for ML/AI pharmacokinetic modeling. Excellent for building predictive models. |

---

## 2. Drug Metabolism Pathway Databases

### 2.1 KEGG DRUG / KEGG PATHWAY

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.genome.jp/kegg/drug/ (DRUG) / https://www.genome.jp/kegg/pathway.html (PATHWAY) |
| **Content** | KEGG DRUG: Approved drugs from Japan, USA, Europe unified by chemical structure. KEGG DGROUP: Drug groups by metabolizing enzymes and transporters. 495 reference pathways across 4,700+ organisms. Drug metabolism annotations, therapeutic targets, molecular interaction networks. |
| **API** | KEGG REST API (https://rest.kegg.jp/) - free for academic use. Supports programmatic access to pathways, compounds, reactions. |
| **License** | Free for academic use. Commercial use requires license from Kanehisa Laboratories. |
| **Size** | 11,000+ drugs, 495 reference pathways, extensive cross-linking |
| **Notes** | Industry standard for pathway analysis. Integrates with DrugBank for drug target downloads. KEGG Mapper enables large-scale data integration. |

### 2.2 PharmGKB (Pharmacogenomics Knowledge Base)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | 715+ drugs, 1,761 genes, 227 diseases, 165 clinical guidelines, 784 drug labels. 9,000+ manually curated literature annotations. 153 drug PK/PD pathways, 68 VIP gene summaries. 26,865 variant annotations. |
| **API** | REST web services for bulk download. User registration required (free). Data packages as zipped spreadsheets. |
| **License** | Free for research. Registration required, agreement to non-redistribution terms. |
| **Size** | 775 drugs, 234 pathways, 26,865 variant annotations, 9,000+ curated papers |
| **Notes** | Premier pharmacogenomics resource. Integrates PharmVar star alleles, CPIC guidelines, FDA labels. Essential for clinical PGx. |

### 2.4 DrugBank

| Attribute | Details |
|-----------|---------|
| **URL** | https://go.drugbank.com/ (Academic: https://go.drugbank.com/academic_research) |
| **Content** | 2,832 drugs, 800 drug metabolites. Comprehensive drug-target, drug-enzyme, drug-transporter data. Cross-references to KEGG, PubChem, ChEBI, PDB. |
| **API** | REST API (Clinical API) for commercial license. Academic license provides XML/CSV/JSON downloads. CC BY-NC 4.0 for academic datasets. |
| **License** | Free academic license (non-commercial). Commercial requires paid license. CC0 1.0 for certain cross-linking datasets. |
| **Size** | 2,832 drugs, 800 metabolites, extensive annotations |
| **Notes** | Part of HMDB suite. Maps to MeSH, ATC, RxNorm. Custom exports available with license. |

### 2.5 CPIC (Clinical Pharmacogenetics Implementation Consortium)

| Attribute | Details |
|-----------|---------|
| **URL** | https://cpicpgx.org/ (API: https://cpicpgx.org/api-and-database/) |
| **Content** | Clinical practice guidelines for gene-drug pairs. Diplotype-to-phenotype mappings, dosing recommendations, CDS language. Allele frequencies, function assignments. Integrated with PharmGKB. |
| **API** | RESTful API - public access. Supports 80,000+ monthly queries. Whole database exports on GitHub. JSON downloads for EHR integration. |
| **License** | Free public access. Data documentation critical for proper use. |
| **Size** | Multiple gene-drug guidelines (growing), comprehensive allele definitions |
| **Notes** | Gold standard for clinical PGx implementation. Integrated into Epic's foundational genomics module. Guidelines peer-reviewed and updated. |

---

## 3. Supplement-Drug Interaction Databases

### 3.1 NatMed Pro (Natural Medicines Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://naturalmedicines.therapeuticresearch.com/ (TRC Healthcare) |
| **Content** | 1,400+ vitamins, minerals, botanicals, non-botanical supplements, foods. Drug-nutrient interactions, mechanism of action, safety data. Multi-ingredient product checking. |
| **API** | RESTful JSON:API for licensed users. Interaction checker endpoint, monograph endpoint. API access requires specific license type. |
| **License** | Subscription-based. Individual/institutional licenses. Available through many academic libraries. |
| **Size** | 1,400+ ingredients, comprehensive interaction database |
| **Notes** | Industry-leading supplement database. Validated for accuracy in 2024 JCO Oncology Advances study. Daily evidence updates. |

### 3.2 MSKCC About Herbs

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.mskcc.org/cancer-care/diagnosis-treatment/symptom-management/integrative-medicine/herbs |
| **Content** | ~284 dietary supplement monographs. Purported uses, adverse effects, herb-drug interactions. Public and healthcare professional versions. |
| **API** | No API. Free web access and mobile app. |
| **License** | Free public access. Maintained by Memorial Sloan Kettering Cancer Center Integrative Medicine Service. |
| **Size** | ~284 supplement monographs |
| **Notes** | Oncology focus but broadly applicable. Continually updated by pharmacist and botanicals experts. Mobile app available. |

### 3.3 Stockley's Drug Interactions

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.medicinescomplete.com/mc/stockley/current/ (subscription) |
| **Content** | Comprehensive drug-drug and drug-herb interactions. 5 severity classifications from "avoid" to "no action needed." Evidence-based recommendations. |
| **API** | No public API. Web access via subscription. |
| **License** | Commercial subscription (Pharmaceutical Press). |
| **Size** | Thousands of interaction monographs |
| **Notes** | Gold standard reference for clinical drug interactions. Used extensively in pharmacy practice. |

### 3.4 FDB (First Databank) MedKnowledge

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.fdbhealth.com/solutions/medknowledge-drug-database/ |
| **Content** | Drug-drug, drug-alternative therapy interactions. Three severity levels plus "clinical effects" subcategories. Includes herbals and dietary supplements. |
| **API** | Commercial API integration available. EHR/pharmacy system integration. |
| **License** | Commercial license required. Enterprise healthcare solutions. |
| **Size** | Comprehensive (enterprise-grade) |
| **Notes** | Major commercial drug database. Widely integrated in healthcare IT systems. |

---

## 4. Natural Product-Drug Interaction Resources

### 4.1 IMgateway (University of Sydney)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.imgateway.net/ |
| **Content** | 1,000+ drug interactions with foods, herbs, supplements. Evidence-based knowledge base. |
| **API** | API available for healthcare application integration. |
| **License** | Partnership/license required. |
| **Size** | 1,000+ interactions |
| **Notes** | Academic-developed, evidence-based approach. |

### 4.2 PHYDGI (Phytotherapy Drug-Interaction Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.synapse-medicine.com/blog/blogpost/herb-drug-interactions-database-phydgi |
| **Content** | Herbal entities and pharmacokinetic interaction strengths. Data from scientific literature and French pharmacovigilance database. |
| **API** | Integration via Synapse Medicine platform |
| **License** | Commercial/research partnerships |
| **Size** | Growing database of herb-drug interactions |
| **Notes** | Includes pharmacovigilance case reports. |

### 4.3 T3DB (Toxin and Toxin Target Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.t3db.ca/ |
| **Content** | 3,678 toxins (pollutants, pesticides, drugs, food toxins) with 42,374 toxin-target associations. 2,073 toxin target records. 90+ data fields per ToxCard including toxicity values, molecular interactions, medical information. |
| **API** | Web interface. Download options available (check website for current formats). |
| **License** | Free access. Part of HMDB suite. |
| **Size** | 3,678 toxins, 2,073 targets, 42,374 associations |
| **Notes** | Focus on toxicity mechanisms and target proteins. Useful for drug-toxin interaction prediction. Links to HMDB and DrugBank. |

### 4.4 HMDB (Human Metabolome Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://hmdb.ca/ (Downloads: https://hmdb.ca/downloads, API: https://hmdb.ca/simple/api) |
| **Content** | 220,945 metabolite entries (water-soluble and lipid-soluble). 8,610 protein sequences (enzymes and transporters). Spectral data (NMR, MS/MS, GC-MS). |
| **API** | API access via contact (academic: eponine@ualberta.ca, commercial: samackay@ualberta.ca). Bulk downloads available. |
| **License** | CC BY-NC 4.0. Free for non-commercial use. |
| **Size** | 220,945 metabolites, 8,610 proteins |
| **Notes** | Core metabolomics resource. Links to DrugBank, T3DB, FooDB. Essential for understanding drug metabolites. |

---

## 5. Absorption/Bioavailability Databases

### 5.1 PK-DB (Pharmacokinetics Database)

| Attribute | Details |
|-----------|---------|
| **URL** | https://pk-db.com/ (API: https://pk-db.com/api/v1/swagger/) |
| **Content** | Pharmacokinetic data from multiple studies. Experimental errors, normalized units, biological ontology annotations. PK parameter calculations from concentration-time profiles. |
| **API** | Full REST API with Swagger documentation. Supports programmatic queries independent of programming language. |
| **License** | Free academic access |
| **Size** | Multi-study aggregated PK data |
| **Notes** | Designed for meta-analysis and PBPK/PK-PD/population PK modeling. Collaborative data curation workflow. Strong validation rules. |

### 5.2 SwissADME

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.swissadme.ch |
| **Content** | Predictive models for physicochemical properties, pharmacokinetics, drug-likeness, medicinal chemistry friendliness. ADME property predictions. |
| **API** | No REST API. Free web tool with user-friendly interface. |
| **License** | Free academic/research use |
| **Size** | Predictive tool (not a static database) |
| **Notes** | Computational predictions, not experimental data. Good for early drug discovery screening. |

### 5.3 FDA BCS (Biopharmaceutics Classification System)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.fda.gov/media/148472/download (M9 Guidance) |
| **Content** | Drug classification by solubility and permeability: Class I (high/high), Class II (low/high), Class III (high/low), Class IV (low/low). Biowaiver guidance. |
| **API** | No API. Regulatory guidance documents (PDF). |
| **License** | Public domain (US Government) |
| **Size** | Classification framework + drug examples |
| **Notes** | Regulatory standard adopted by FDA, EMA, WHO. Used for bioequivalence study waivers. |

### 5.4 FDA Table of Pharmacogenomic Biomarkers in Drug Labeling

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.fda.gov/media/124784/download (Updated June 2024) |
| **Content** | Approved products with PGx information in drug labeling. Biomarkers per HUGO nomenclature. Germline/somatic variants, functional deficiencies, gene expression differences. |
| **API** | No official API. PharmGKB provides parser tool: https://github.com/PharmGKB/fda-biomarker |
| **License** | Public domain (US Government) |
| **Size** | 300+ drug labels with PGx biomarkers |
| **Notes** | Official FDA resource. Updated regularly. Indicates labeling sections containing biomarker info. |

---

## 6. Transporter Genetics Databases (OATP, P-gp)

### 6.1 PharmGKB Transporter Annotations

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | SLCO1B1 (OATP1B1), ABCB1 (MDR1/P-gp), ABCG2 (BCRP), OCT1, OAT1, OAT3, MATE1, MATE2-K genetic variants. Clinical annotations for transporter polymorphisms. Statin-SLCO1B1 associations extensively documented. |
| **API** | REST web services (registration required) |
| **License** | Free research use with registration |
| **Size** | Comprehensive transporter pharmacogenomics data |
| **Notes** | Best curated resource for transporter PGx. SLCO1B1 521T>C and 388A>G extensively characterized. |

### 6.2 CPIC Guidelines for Transporters

| Attribute | Details |
|-----------|---------|
| **URL** | https://cpicpgx.org/guidelines/ |
| **Content** | Clinical guidelines for SLCO1B1 and statin dosing. Transporter genotype-to-phenotype translations. |
| **API** | CPIC API (see Section 2.5) |
| **License** | Free public access |
| **Size** | Guideline-specific |
| **Notes** | SLCO1B1 genotyping for statin-induced myopathy prevention is a model case for transporter pharmacogenetics. |

### 6.3 International Transporter Consortium (ITC) Resources

| Attribute | Details |
|-----------|---------|
| **URL** | Publications in PMC (e.g., https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6348109/) |
| **Content** | Recommendations for transporter polymorphism investigation in drug development. Guidance on ABCG2, SLCO1B1, SLC22A1 (OCT1). FDA regulatory guidance alignment. |
| **API** | No API - literature/regulatory guidance |
| **License** | Open access publications |
| **Size** | Guidance documents |
| **Notes** | Provides scientific consensus on clinically relevant transporter variants. |

### 6.4 Key Transporter Polymorphisms Reference

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

## 7. Summary Comparison Table

| Database | Focus | API | License | Size | Best For |
|----------|-------|-----|---------|------|----------|
| **PharmVar** | CYP allele nomenclature | No | Free | 2000+ SNPs | Star allele definitions |
| **PharmGKB** | PGx knowledge | REST | Free (registration) | 775 drugs, 26K annotations | Clinical PGx |
| **CPIC** | Clinical guidelines | REST | Free | Guidelines | EHR integration |
| **DrugBank** | Drug data | REST (licensed) | Free academic | 2,832 drugs | Drug metabolism |
| **KEGG** | Pathways | REST | Free academic | 11K+ drugs | Pathway analysis |
| **SuperCYP** | CYP interactions | No | Free | 1,170 drugs | CYP-specific |
| **NatMed Pro** | Supplements | REST (licensed) | Subscription | 1,400+ ingredients | Clinical supplement checking |
| **MSKCC About Herbs** | Herbs/supplements | No | Free | 284 monographs | Oncology herb safety |
| **PK-DB** | Pharmacokinetics | REST | Free | Multi-study | PBPK modeling |
| **HMDB** | Metabolome | Contact | CC BY-NC 4.0 | 220K metabolites | Metabolomics |
| **T3DB** | Toxins | Web download | Free | 3,678 toxins | Toxicity data |
| **Flockhart Table** | CYP clinical | No | Free | Curated | Clinical teaching |
| **FDA PGx Labels** | Regulatory | Parser available | Public domain | 300+ labels | Regulatory compliance |

---

## 8. Integration Recommendations

### 8.1 Priority Data Sources for Gene Platform

**Tier 1 - Essential (Free, High-Quality, API Access)**
1. **PharmGKB** - Core PGx knowledge base, REST API
2. **CPIC** - Clinical guidelines, REST API, EHR-ready
3. **PharmVar** - Star allele standards
4. **KEGG** - Pathway integration, REST API
5. **FDA PGx Labels** - Regulatory data, parser available

**Tier 2 - Important (Free, Manual Integration)**
1. **DrugBank** (academic license) - Comprehensive drug data
2. **HMDB** - Metabolite data
3. **SuperCYP** - CYP interaction details
4. **Flockhart Table** - Clinical CYP reference

**Tier 3 - Supplementary (Specialized or Licensed)**
1. **NatMed Pro** - Supplement interactions (if licensed)
2. **MSKCC About Herbs** - Free herb data
3. **PK-DB** - PK modeling data
4. **T3DB** - Toxicity mechanisms

### 8.2 Data Integration Strategy

```
                    ┌─────────────────────────────────────┐
                    │       Gene Platform Core           │
                    └─────────────────────────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        │                           │                           │
        ▼                           ▼                           ▼
┌───────────────┐          ┌───────────────┐          ┌───────────────┐
│   Genotype    │          │  Drug/Pathway │          │  Interactions │
│    Layer      │          │    Layer      │          │    Layer      │
├───────────────┤          ├───────────────┤          ├───────────────┤
│ • PharmVar    │          │ • DrugBank    │          │ • NatMed Pro  │
│ • PharmGKB    │   ───►   │ • KEGG        │   ───►   │ • MSKCC       │
│ • CPIC        │          │ • Reactome    │          │ • SuperCYP    │
│ • FDA Labels  │          │ • HMDB        │          │ • Flockhart   │
└───────────────┘          └───────────────┘          └───────────────┘
        │                           │                           │
        └───────────────────────────┼───────────────────────────┘
                                    │
                                    ▼
                    ┌─────────────────────────────────────┐
                    │     Personalized Recommendations    │
                    │  • CYP metabolizer status          │
                    │  • Transporter function            │
                    │  • Drug-supplement interactions    │
                    │  • Dosing adjustments              │
                    └─────────────────────────────────────┘
```

### 8.3 Key Identifiers for Cross-Linking

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

### 8.4 Estimated Data Volumes

| Data Type | Estimated Records | Storage Estimate |
|-----------|-------------------|------------------|
| Drug entries | ~15,000 (deduplicated) | ~500 MB |
| CYP variants | ~5,000 | ~50 MB |
| Gene-drug annotations | ~30,000 | ~200 MB |
| Pathway diagrams | ~150,000 | ~2 GB (images) |
| Interaction records | ~50,000 | ~300 MB |
| PK parameters | ~100,000 | ~100 MB |
| **Total** | - | **~3-4 GB** |

### 8.5 API Integration Priority

1. **CPIC API** - Immediate clinical utility
2. **KEGG REST** - Pathway visualization
3. **PK-DB API** - Pharmacokinetic modeling
4. **PharmGKB Web Services** - Comprehensive PGx data
5. **DrugBank API** (if licensed) - Drug metabolism details

---

## References

1. PharmVar Consortium: https://www.pharmvar.org/
2. PharmGKB: https://www.pharmgkb.org/
3. CPIC: https://cpicpgx.org/
4. DrugBank: https://go.drugbank.com/
5. KEGG: https://www.genome.jp/kegg/
6. Reactome: https://reactome.org/
7. HMDB: https://hmdb.ca/
8. NatMed Pro: https://naturalmedicines.therapeuticresearch.com/
9. MSKCC About Herbs: https://www.mskcc.org/cancer-care/diagnosis-treatment/symptom-management/integrative-medicine/herbs
10. SuperCYPsPred: https://insilico-cyp.charite.de/SuperCYPsPred/
11. Flockhart Table: https://drug-interactions.medicine.iu.edu/
12. PK-DB: https://pk-db.com/
13. T3DB: https://www.t3db.ca/
14. FDA PGx Biomarkers: https://www.fda.gov/drugs/science-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling
15. SwissADME: http://www.swissadme.ch/

---

*Research compiled: January 2026*
*For Gene Platform pharmacogenomics integration*
