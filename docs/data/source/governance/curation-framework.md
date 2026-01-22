# Data Sources Pruning Analysis

**Created:** January 2026
**Updated:** January 2026 (added original research evaluation)
**Purpose:** Critical evaluation of ALL databases (original ~40 + new 200+) for removal

---

## Executive Summary

| Pruning Reason | New DBs | Original DBs | Total |
|----------------|---------|--------------|-------|
| **Commercial/Paid License** | 18 | 0 | 18 |
| **Controlled Access** | 15 | 0 | 15 |
| **Legal/Scraping Risk** | 12 | 0 | 12 |
| **Out of Scope** | 22 | 2 | 24 |
| **Redundant/Superseded** | 14 | 6 | 20 |
| **Defunct/Unavailable** | 6 | 1 | 7 |
| **Low Value** | 16 | 2 | 18 |
| **TOTAL PRUNABLE** | **103** | **11** | **114** |
| **KEEP (New)** | ~97 | - | - |
| **KEEP (Original)** | - | ~29 | - |
| **TOTAL KEEP** | - | - | **~126** |

---

## Original Research Databases Evaluation (~40 databases)

### PRUNE from Original Research (12 databases)

| Database | Source | Reason | Alternative |
|----------|--------|--------|-------------|
| **OpenSNP** | data-sources-public | DEFUNCT - Shut down April 2025 | Use gnomAD, archived Zenodo |
| ~~SMPDB~~ | ~~data-sources-public~~ | **MOVED TO KEEP** - Unique small molecule pathways | Complements Reactome |
| ~~STITCH~~ | ~~data-sources-public~~ | **MOVED TO KEEP** - Chemical-protein interactions | Complements STRING |
| ~~FooDB~~ | ~~western-herbal-models~~ | **MOVED TO KEEP** - Nutrigenomics value | Tier 2 |
| **Natural Products Atlas** | western-herbal-models | REDUNDANT - COCONUT 2.0 is larger | COCONUT 2.0 |
| **Health Canada LNHPD** | western-herbal-models | OUT OF SCOPE - Canadian regulatory only | DSLD (US, CC0) |
| **WAKANYAKU** | data-sources-public | LOW VALUE - niche Japanese academic | KampoDB covers Kampo |
| **TTD (standalone)** | drug-target-models | REDUNDANT - covered by DGIdb aggregator | DGIdb, DrugBank |
| **1000 Genomes (full)** | data-sources-public | LOW VALUE - 260TB, only need small subsets | gnomAD (derived frequencies) |
| **KEGG (FTP bulk)** | data-sources-public | LICENSE - Academic license for bulk | Reactome + WikiPathways |
| **TCMSP (standalone)** | data-sources-public | REDUNDANT - superseded by BATMAN-TCM 2.0 | BATMAN-TCM 2.0 |
| **TCMBank (duplicate)** | tcm-data-models | REDUNDANT - listed twice, keep one entry | Single entry in Tier 2 |

### KEEP from Original Research (28 databases)

#### Tier 1 - Already in MVP list
| Database | License | Notes |
|----------|---------|-------|
| ClinVar | Public domain | Core variant database |
| PharmGKB | CC BY-SA 4.0 | Drug-gene relationships |
| Reactome | CC0 | Pathway gold standard |
| HPO | Open | Phenotype ontology |

#### Tier 2 - High Value Original
| Database | License | Notes |
|----------|---------|-------|
| dbSNP | Public domain | 1.2B SNPs reference |
| gnomAD | ODC-ODbL | Allele frequencies |
| GWAS Catalog | Open | Trait associations |
| SNPedia | CC BY-NC-SA | Curated 109K SNPs |
| GeneCards | Free academic | Gene summaries |
| UniProt | CC BY 4.0 | Protein database |
| NCBI Gene | Public domain | Gene reference |
| Ensembl | Open | Genome annotations |
| WikiPathways | CC0 | Community pathways |
| ChEMBL | CC BY-SA 3.0 | Bioactivity data |
| Open Targets | Open | Gene-disease associations |
| DSLD | CC0 | Supplement labels |
| NIH ODS API | Public domain | Supplement facts |
| FDA DailyMed | Public domain | Drug labels |
| DisGeNET | CC BY-NC-SA | Gene-disease associations |
| OMIM | Academic API | Gene-disease (curated) |
| DGIdb | Open | Drug-gene aggregator |
| BindingDB | Open | Binding affinities |
| GtoPdb | Open | Pharmacology guide |
| Dr. Duke's | CC0 | Ethnobotanical |
| KampoDB | Academic | Japanese Kampo |
| PubMed | Public domain | Literature abstracts |
| PMC OA | Various open | Full-text articles |
| IMPPAT | CC BY 4.0 | Ayurveda (moved to Tier 1) |

---

## Category 1: Commercial/Paid License (REMOVE - 18 databases)

These require paid subscriptions or commercial licenses incompatible with an open platform.

| Database | Source File | License Issue | Alternative |
|----------|-------------|---------------|-------------|
| ~~COSMIC~~ | ~~cancer-oncology~~ | **MOVED TO KEEP** - Academic license free | Best cancer mutation DB |
| **GOBIOM Database** | biomarkers-labs | Commercial subscription required | Use MarkerDB (free academic) |
| **ConsumerLab** | scraping-targets | $99/2yr subscription, no scraping | Use Examine.com (partial free) |
| **Labdoor** | scraping-targets | Proprietary test data | Skip supplement testing data |
| **NatMed Pro** | drug-metabolism | $695/yr subscription | Use DrugBank free tier |
| **MSKCC About Herbs** | drug-metabolism | Limited free, full access paid | Use IMPPAT + BATMAN-TCM |
| **IPD-IMGT/HLA** | hormones-autoimmune | Complex license, commercial restrictions | Use minimal HLA data from ClinVar |
| **BRCA Exchange** | cancer-oncology | Complex data sharing terms | Use ClinVar BRCA variants |
| **Genographic Project** | genetics-expanded | Closed (National Geographic ended) | Use gnomAD + 1000 Genomes |
| **deCODE genetics** | genetics-expanded | Proprietary Amgen data | Use public GWAS Catalog |
| **InSiGHT (MMR)** | cancer-oncology | Database access fees | Use ClinVar for MMR variants |
| **NCBI RefSeq Select** | genetics-expanded | Use requires attribution complexity | Use Ensembl/MANE |
| **LabCorp/Quest ranges** | biomarkers-labs | Proprietary reference ranges | Use CALIPER + HMDB |
| **WebMD** | scraping-targets | Commercial site, heavy anti-scraping | Use MedlinePlus (public domain) |
| **Examine.com Premium** | scraping-targets | Free limited, premium gated | Use free tier only |
| **Mayo Clinic Lab** | biomarkers-labs | No bulk access, commercial | Use ABIM + CALIPER |
| **UpToDate** | scraping-targets | Commercial medical reference | Use PubMed abstracts |
| **Elsevier ClinicalKey** | papers-literature | Subscription required | Use OpenAlex (CC0) |

---

## Category 2: Controlled Access / Lengthy Approval (REMOVE - 15 databases)

These require formal applications, IRB approval, or data access agreements with uncertain timelines.

| Database | Source File | Access Issue | Alternative |
|----------|-------------|--------------|-------------|
| **UK Biobank** | rare-diseases | 6-12 month application, fees | Use gnomAD (derived data free) |
| **BioBank Japan** | rare-diseases | Japanese researcher priority | Use published GWAS summaries |
| **dbGaP** | genetics-expanded | NIH DAR process (3-6 months) | Use GWAS Catalog summaries |
| **All of Us (full)** | genetics-expanded | Requires researcher registration | Use public summary stats |
| **TKDL (Traditional Knowledge)** | trad-med-expanded | Government of India restricted | Use IMPPAT (open alternative) |
| **GDC/TCGA (controlled)** | cancer-oncology | DAA required for individual data | Use open-access tiers |
| **TOPMed BRAVO (full)** | genetics-expanded | dbGaP controlled access | Use gnomAD (overlapping data) |
| **NIMH Genetics** | mental-cognitive | NIH application process | Use PGC public summaries |
| **GenomeConnect** | womens-pediatric | Patient consent gated | Use ClinVar submissions |
| **PatientsLikeMe** | communities-biohackers | Partnership agreements only | Use Open Humans |
| **DECIPHER (bulk)** | rare-diseases | Data Access Agreement | Use open-access tier |
| **H3Africa** | genetics-expanded | Governance approval needed | Use gnomAD for African alleles |
| **DIAGRAM consortium** | cardio-metabolic | Membership required | Use GWAS Catalog T2D |
| **Estonian Biobank** | genetics-expanded | Estonian ethics approval | Use published summary stats |
| **Framingham Heart Study** | cardio-metabolic | dbGaP controlled | Use public GWAS results |

---

## Category 3: Legal/Scraping Risk (REMOVE - 12 databases)

These involve Terms of Service violations, lawsuit precedents, or significant legal exposure.

| Database | Source File | Legal Issue | Alternative |
|----------|-------------|-------------|-------------|
| **Reddit forums** | scraping-targets | ToS violation, Reddit sued Perplexity | Use formal research API if approved |
| **HealthUnlocked** | communities-biohackers | Patient privacy, no scraping allowed | Skip - partnership only |
| **MyFitnessPal** | scraping-targets | Under Armour ToS, login required | Skip - no alternative |
| **ConsumerLab** | scraping-targets | Copyright on test results | Skip (see paid above) |
| **Drugs.com** | drug-metabolism | Copyright, anti-scraping | Use DrugBank |
| **WebMD** | scraping-targets | Commercial, robot blocks | Use MedlinePlus |
| **RxList** | drug-metabolism | Copyright content | Use DrugBank |
| **PubMed Central (full text)** | papers-literature | PMC OA subset only legal | Use PMC OA + OpenAlex |
| **Google Scholar** | papers-literature | Explicit prohibition in ToS | Use Semantic Scholar API |
| **Phoenix Rising (bulk)** | communities-biohackers | Personal health data concerns | Use public posts only with care |
| **23andMe community** | communities-biohackers | Private platform, no data export | Skip |
| **Ancestry forums** | communities-biohackers | Private platform | Skip |

---

## Category 4: Out of Scope for Gene Platform (REMOVE - 22 databases)

These are tangential to the core mission (genetics + traditional medicine + personalized health).

| Database | Source File | Why Out of Scope |
|----------|-------------|------------------|
| **HOMD (oral microbiome)** | oral-skin-sensory | Dental-specific, not systemic health |
| **Skin Gene Database** | oral-skin-sensory | Dermatology-specific |
| **RetNet** | oral-skin-sensory | Ophthalmology-specific |
| **BitterDB** | oral-skin-sensory | Taste research, not therapeutics |
| **CircaDB** | sleep-longevity-nutri | Mouse circadian data, not human |
| **4D Nucleome** | genetics-expanded | Research tool, not clinical |
| **ENCODE (full)** | genetics-expanded | Keep only cCRE subset |
| **FANTOM5 (full)** | genetics-expanded | Keep only enhancer/promoter subset |
| **Roadmap Epigenomics** | genetics-expanded | Superseded by ENCODE 4 |
| **IHEC Data Portal** | genetics-expanded | Redundant with ENCODE |
| **MethBank 4.0** | genetics-expanded | Methylation atlas - niche |
| **PhosphoSitePlus** | drug-metabolism | Phosphorylation research focus |
| **BioGRID** | drug-metabolism | Interaction research, redundant with STRING |
| **IEDB (full)** | hormones-autoimmune | Keep epitope subset only |
| **ThyroidOmics** | hormones-autoimmune | Too specialized |
| **ImmunoBase (full)** | hormones-autoimmune | Use GWAS Catalog instead |
| **Probio-ichnos** | microbiome | Probiotic research niche |
| **PharmacoMicrobiomics** | microbiome | Research-stage interactions |
| **MetaHIT** | microbiome | Superseded by HMP/GMrepo |
| **SV4GD** | genetics-expanded | SV research database, niche |
| **HGSVC** | genetics-expanded | Assembly reference, not clinical |
| **DGV** | genetics-expanded | Use gnomAD-SV instead |

---

## Category 5: Redundant / Superseded (REMOVE - 14 databases)

Better open alternatives exist that cover the same data.

| Database | Source File | Superseded By | Reason |
|----------|-------------|---------------|--------|
| **DGV** | genetics-expanded | gnomAD-SV v4.1 | gnomAD has frequency + functional data |
| **SuperNatural 3.0** | trad-med-expanded | COCONUT 2.0 | COCONUT includes SuperNatural |
| **BIOFACQUIM** | trad-med-expanded | COCONUT 2.0 | Latin American NPs in COCONUT |
| **MetaHIT** | microbiome | GMrepo v3 | GMrepo is more current |
| **Roadmap Epigenomics** | genetics-expanded | ENCODE 4 | ENCODE 4 supersedes |
| **IHEC Data Portal** | genetics-expanded | ENCODE 4 | Redundant epigenome data |
| **CureTogether** | communities-biohackers | 23andMe (acquired 2012) | Data absorbed into 23andMe |
| **REVEL (standalone)** | genetics-expanded | dbNSFP v4.9 | dbNSFP includes REVEL scores |
| **EVE (standalone)** | genetics-expanded | dbNSFP v4.9 | dbNSFP includes EVE scores |
| **PrimateAI-3D** | genetics-expanded | AlphaMissense | AlphaMissense more comprehensive |
| **Europa PMC** | papers-literature | OpenAlex | OpenAlex has better coverage + API |
| **Semantic Scholar** | papers-literature | OpenAlex | OpenAlex is CC0, larger |
| **ChEBI** | mental-cognitive | ChEMBL | ChEMBL includes chemical ontology |
| ~~SMPDB~~ | ~~drug-metabolism~~ | **MOVED TO KEEP** | Unique metabolite pathways |

---

## Category 6: Defunct / Unavailable (REMOVE - 6 databases)

These are no longer maintained, have broken links, or have ceased operations.

| Database | Source File | Status |
|----------|-------------|--------|
| **Genographic Project** | genetics-expanded | Shut down 2020 |
| **CureTogether** | communities-biohackers | Absorbed into 23andMe 2012 |
| **p-ANAPL** | trad-med-expanded | Physical library only, not digital |
| **BSI Medicinal Plants** | trad-med-expanded | Limited digital access |
| **Several older TCM DBs** | trad-med-expanded | Not maintained, 404 errors |
| **HPGDB (histamine)** | allergy-pain | Outdated, sparse data |

---

## Category 7: Low Value / High Effort (REMOVE - 18 databases)

The effort to integrate exceeds the value for the Gene Platform.

| Database | Source File | Why Low Value |
|----------|-------------|---------------|
| **Wikipedia infoboxes** | wikipedia-wikidata | Unstructured, unreliable |
| **Pain genetics GWAS** | allergy-pain | Sparse, early-stage research |
| **HAGR (GenAge/GenDR)** | sleep-longevity-nutri | Model organism focus |
| **CircaDB** | sleep-longevity-nutri | Mouse data, not human |
| **Open Genes** | sleep-longevity-nutri | Overlaps HAGR |
| ~~FooDB~~ | ~~sleep-longevity-nutri~~ | **MOVED TO KEEP** - nutrigenomics |
| **LabTestsOnline** | biomarkers-labs | Consumer info, no data API |
| **BIONDA** | biomarkers-labs | Text-mined, lower quality |
| **TheMarker** | biomarkers-labs | No API, web interface only |
| **ABIM ranges** | biomarkers-labs | PDF only, use CALIPER |
| **UMLS full** | biomarkers-labs | 30GB, use LOINC directly |
| **Genspace projects** | communities-biohackers | Physical lab, minimal digital data |
| **r/StackAdvice** | communities-biohackers | Anecdotal, not evidence-based |
| **MEpedia** | communities-biohackers | Wiki format, unreliable |
| **Health Rising** | communities-biohackers | News site, not data source |
| **Quantified Self Archive** | communities-biohackers | Sparse project documentation |
| **Fitabase** | scraping-targets | Wearable research, out of scope |

---

## KEEP List: High-Value Databases (~95 databases)

### Tier 1 - Critical for MVP (16 databases)

| Database | Category | Why Keep |
|----------|----------|----------|
| **dbNSFP v4.9** | Genetics | 35+ pathogenicity scores in one file |
| **AlphaMissense** | Genetics | Best protein variant predictor |
| **gnomAD-SV v4.1** | Genetics | Only open SV frequency database |
| **ClinVar** | Genetics | Clinical variant gold standard |
| **COCONUT 2.0** | Traditional Med | 800K+ natural products, open API |
| **BATMAN-TCM 2.0** | Traditional Med | TCM compound-target mappings |
| **IMPPAT 2.0** | Traditional Med | Ayurveda CC BY 4.0 (commercial OK) |
| **GMrepo v3** | Microbiome | Gut microbiome-phenotype, REST API |
| **OncoKB** | Cancer | FDA-recognized, free academic API |
| **CIViC** | Cancer | Open cancer interpretation |
| **CPIC** | Pharmacogenomics | Clinical guidelines API |
| **PharmVar** | Pharmacogenomics | Star allele definitions |
| **PGC** | Mental Health | Psychiatric GWAS summaries |
| **Orphanet** | Rare Disease | Rare disease reference, CC BY 4.0 |
| **HPO** | Phenotypes | Phenotype ontology standard |
| **OpenAlex** | Literature | 271M papers, CC0, good API |

### Tier 2 - High Value (35 databases)

| Database | Category | Why Keep |
|----------|----------|----------|
| **SNPedia** | Genetics | Curated 109K SNPs |
| **PharmGKB** | Pharmacogenomics | Drug-gene annotations |
| **DrugBank** | Drugs | Comprehensive drug database |
| **Reactome** | Pathways | Pathway gold standard |
| **KEGG** | Pathways | Metabolic pathways (attribution OK) |
| **STRING** | Interactions | Protein interactions |
| **ENCODE 4 (cCRE)** | Epigenetics | Regulatory elements |
| **SpliceAI** | Genetics | Splice site predictions |
| **CADD v1.7** | Genetics | Variant deleteriousness |
| **SymMap 2.0** | TCM | Symptom-disease mapping |
| **ETCM v2.0** | TCM | Formula composition |
| **HERB 2.0** | TCM | Herb transcriptomics |
| **LOTUS** | Natural Products | Wikidata compound-taxon |
| **ANPDB** | African Med | African natural products |
| **SANCDB** | African Med | South African compounds |
| **NuBBEDB** | Latin Am Med | Brazilian biodiversity |
| **HMP Portal** | Microbiome | Human microbiome project |
| **gutMGene v2.0** | Microbiome | Gut microbiome-gene |
| **CTD** | Environmental | Chemical-gene-disease |
| **MITOMAP** | Mitochondrial | mtDNA variants |
| **MitoCarta** | Mitochondrial | Mitochondrial genes |
| **ClinGen** | Clinical | Gene curation authority |
| **DECIPHER (open)** | Rare Disease | Open-access CNV data |
| **Monarch Initiative** | Phenotypes | Cross-species phenotypes |
| **MarkerDB 2.0** | Biomarkers | Biomarker reference |
| **HMDB** | Metabolites | Metabolite concentrations |
| **FooDB** | Nutrigenomics | 100K+ food compounds, vitamins, minerals, bioactives - nutrient-gene interactions |
| **LOINC** | Labs | Lab test codes |
| **CALIPER** | Labs | Pediatric reference intervals |
| **ClinicalTrials.gov** | Trials | Trial registry, public API |
| **PubMed** | Literature | Biomedical abstracts |
| **Wikidata (full)** | Knowledge Graph | Universal hub: drugs, vitamins, genes, proteins, diseases, compounds, biological processes - downloadable dumps |
| **Allen Brain Atlas** | Brain | Brain expression |
| **GTEx** | Expression | Tissue eQTLs |
| **BrainSpan** | Brain Dev | Developmental expression |
| **CARDIoGRAMplusC4D** | Cardio | Heart disease GWAS |

### Tier 3 - Specialized / As Needed (~44 databases)

These are kept but lower priority for schema documentation.

| Category | Databases |
|----------|-----------|
| **Genetics Specialized** | RegulomeDB v2, MaveDB, gnomAD (main), ALFA |
| **TCM Specialized** | TCMID 2.0, TCMBank, CMAUP |
| **Microbiome** | mBodyMap, MASI, MDAD |
| **Drug/Supplement** | T3DB, PK-DB, SuperCYP |
| **Cancer** | Cancer Gene Census |
| **Rare Disease** | DDG2P, OMIM (partial) |
| **Biomarkers** | Metabolomics Workbench |
| **Cardio/Metabolic** | GLGC |
| **Community (careful use)** | Phoenix Rising (public posts), Open Humans |

---

## Pruning Summary

### Combined (Original + New Research)

| Status | Count | Source |
|--------|-------|--------|
| **PRUNE (New)** | 103 | New 200+ databases |
| **PRUNE (Original)** | 11 | Original ~40 databases |
| **TOTAL PRUNED** | **114** | |
| **KEEP (Tier 1)** | 16 | MVP critical |
| **KEEP (Tier 2)** | 35 | High value |
| **KEEP (Tier 3)** | 44 | Specialized |
| **KEEP (Original)** | 28 | Already documented |
| **TOTAL KEEP** | **~126** | |
| **TOTAL EVALUATED** | **~240** | |

### Pruning Rate
- New databases: 103/200 = **52% pruned**
- Original databases: 11/40 = **28% pruned**
- Overall: 114/240 = **48% pruned**

---

## Recommended Actions

1. **Remove 103 databases** from consideration
2. **Prioritize Tier 1 (16)** for immediate schema documentation
3. **Document Tier 2 (35)** after Tier 1 complete
4. **Keep Tier 3 (44)** in inventory for specialized needs
5. **Total effort reduced by 48%** by pruning

---

## Pruning Rationale by Category

| Reason | Count | Key Examples |
|--------|-------|--------------|
| Commercial/Paid | 17 | ConsumerLab, NatMed Pro, GOBIOM |
| Controlled Access | 15 | UK Biobank, dbGaP, TCGA controlled |
| Legal Risk | 12 | Reddit scraping, HealthUnlocked, Google Scholar |
| Out of Scope | 22 | HOMD, BitterDB, 4D Nucleome |
| Redundant | 11 | DGV→gnomAD-SV, REVEL→dbNSFP |
| Defunct | 6 | Genographic, CureTogether |
| Low Value | 16 | Wikipedia infoboxes, pain GWAS, CircaDB |
| **Moved to Keep** | 4 | STITCH, SMPDB, COSMIC, FooDB |

---

*Generated: January 2026*
