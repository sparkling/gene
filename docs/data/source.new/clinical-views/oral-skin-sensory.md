---
id: domain-oral-skin-sensory
title: "Oral Health, Skin Health, and Sensory Genetics Data Sources"
type: health-domain
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [disease, health-domain, databases]
---

**Parent:** [Health Domains](./_index.md)

# Oral Health, Skin Health, and Sensory Genetics Data Sources

Research compiled for the Gene Platform - Comprehensive database inventory for oral microbiome, dermatogenomics, and sensory genetics (vision, hearing, taste/smell).

---

## Table of Contents

1. [Oral Health / Dental Microbiome Databases](#1-oral-health--dental-microbiome-databases)
2. [Skin Condition Genetics / Dermatogenomics](#2-skin-condition-genetics--dermatogenomics)
3. [Vision / Eye Genetics Databases](#3-vision--eye-genetics-databases)
4. [Hearing Genetics Databases](#4-hearing-genetics-databases)
5. [Taste / Smell Genetics Databases](#5-taste--smell-genetics-databases)
6. [Summary Comparison Table](#6-summary-comparison-table)

---

## 1. Oral Health / Dental Microbiome Databases

### 1.1 Human Oral Microbiome Database (HOMD / eHOMD)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.homd.org/ |
| **Alternative URL** | http://www.ehomd.org/ |
| **Download Page** | https://v31.homd.org/download |
| **Description** | The first and most comprehensive human body-site specific microbiome database, part of the NIH Human Microbiome Project. Contains curated 16S rRNA gene reference sequences linked to genomes for species-level taxonomy assignment. |
| **Content** | 836 taxa total (525 primarily oral, 22 primarily nasal); 49% named species, 21% unnamed but cultivated, 29% uncultivated phylotypes |
| **Current Version** | Taxonomy V4.1, 16S rRNA RefSeq V16.03, Genomic RefSeq V11.02, Viruses V1.2 |
| **Data Formats** | FASTA (aligned/unaligned), Taxonomy files for MOTHUR and QIIME, SVG phylogenetic trees, PROKKA-annotated genomes |
| **API** | BLAST servers (genome and RefSeq), JBrowse genome viewer, FTP site access |
| **License** | Free academic use; cite when publishing |
| **Size** | ~836 taxa with full genome sequences for many species |
| **Citations** | >5,500 scientific publications |

**Key Features:**
- 16S rRNA training sequence sets (full-length and V1V3)
- Prophage predictions for all genomes (beta)
- Cross-links to PubMed, Entrez, and other databases
- Systematic naming scheme for unnamed/uncultivated taxa

---

### 1.2 Human Reference Oral Microbiome (HROM)

| Attribute | Details |
|-----------|---------|
| **URL** | Published in Cell Host & Microbe (2025) |
| **Source** | https://www.cell.com/cell-host-microbe/abstract/S1931-3128(25)00415-9 |
| **Description** | High-quality genomic catalog from 7,956 oral samples spanning 21 countries, covering diverse oral sites including saliva and dental samples. |
| **Content** | 72,641 high-quality genomes from 3,426 species (2,019 previously unidentified species) |
| **Key Finding** | 1,137 previously uncharacterized candidate phyla radiation (CPR) species; Patescibacteria identified as most prevalent oral phylum |
| **Data Formats** | Metagenome-assembled genomes (MAGs) |
| **API** | Supplementary data available with publication |
| **License** | Academic research use (journal terms) |
| **Size** | 72,641 genomes |

---

### 1.3 Cultivated Oral Bacteria Genome Reference (COGR)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.nature.com/articles/s41522-023-00414-3 |
| **Description** | Reference genomes based on large-scale cultivation of human oral bacteria isolated from dental plaques, tongue, and saliva. |
| **Content** | 1,089 high-quality genomes; 195 species-level clusters (95 with no taxonomic annotation); 315 genomes representing novel species |
| **Coverage** | 5 phyla; 111 person-specific clusters |
| **Data Formats** | Genome sequences (FASTA), annotations |
| **API** | Supplementary data with publication |
| **License** | Open access (Nature npj) |
| **Size** | 1,089 genomes |

---

### 1.4 Periodontal Disease GWAS Data

| Attribute | Details |
|-----------|---------|
| **URL** | GWAS Catalog (https://www.ebi.ac.uk/gwas/) |
| **Key Studies** | GLIDE Consortium, UK Biobank, All of Us Research Program |
| **Description** | Genome-wide association studies for dental caries and periodontitis |
| **Content** | 47+ novel risk loci for dental caries; Multiple ancestry-specific associations |
| **Heritability** | ~50% for periodontitis |
| **Key Loci** | CLEC19A, TRA, GGTA2P, TM9SF2, IFI16, RBMS3, SLC15A4, PKP2, SNRPN |
| **Data Formats** | Summary statistics (TSV/CSV) |
| **API** | GWAS Catalog REST API |
| **License** | CC0 or EMBL-EBI terms |

---

## 2. Skin Condition Genetics / Dermatogenomics

### 2.1 GWAS Catalog - Skin Diseases

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **REST API** | https://www.ebi.ac.uk/gwas/rest/docs/api |
| **Summary Statistics API** | https://www.ebi.ac.uk/gwas/summary-statistics/docs/ |
| **Description** | Curated collection of published GWAS for psoriasis, atopic dermatitis, eczema, acne, and other skin conditions |

#### Psoriasis GWAS Data
| Metric | Value |
|--------|-------|
| **Total Loci** | 109 distinct susceptibility loci (46 novel) |
| **Study Size** | 36,466 cases / 458,078 controls (European meta-analysis) |
| **Multi-Ancestry** | 28,869 cases / 443,950 controls; 74 genome-wide significant loci (32 novel) |
| **Heritability** | 66% (overall); 15% SNP-based |
| **Key Genes** | HLA-C*06:02, IL23R, IL23A, IL12B, TRAF3IP2, ERAP1, RNF114, IFIH1, LCE3B/3D |
| **Study Accessions** | Multiple; search "psoriasis" in GWAS Catalog |

#### Atopic Dermatitis GWAS Data
| Metric | Value |
|--------|-------|
| **Total Loci** | 91 genetic loci (multi-ancestry) |
| **Study Size** | >1 million individuals; 56,146 cases / 602,280 controls |
| **Novel Loci** | 29 novel (European); 16 novel (multi-ancestry) |
| **Key Drug Targets** | CSF1, CTSS, IL15, MMP12 |
| **Study Accession** | GCST003184 (example) |

**Data Formats:** Summary statistics (TSV), association data (JSON/TSV)
**API:** REST API with R package (gwasrapidd)
**License:** CC0 or EMBL-EBI standard terms

---

### 2.2 SKIOME Project

| Attribute | Details |
|-----------|---------|
| **URL** | https://github.com/giuliaago/SKIOMEMetadataRetrieval |
| **Publication** | https://academic.oup.com/database/article/doi/10.1093/database/baac033/6586378 |
| **Description** | Curated collection of skin microbiome datasets enriched with study-related metadata for meta-analysis |
| **Content** | Aggregated skin microbiome sequencing datasets from public databases |
| **Data Formats** | Data frames, supplementary files |
| **License** | Open access |

---

### 2.3 Expanded Skin Microbial Genome Collection (eSMGC)

| Attribute | Details |
|-----------|---------|
| **Source** | iMetaOmics (2025) |
| **Description** | Human skin microbiome reference catalog |
| **Content** | 675 prokaryotic genomes, 12 fungal genomes, 2,344 viral genomes, 4,935 plasmid genomes |
| **Source Data** | 2,264 publicly available metagenomic datasets |
| **Data Formats** | Genome sequences |
| **License** | Journal open access terms |

---

### 2.4 Early-Life Skin Genomes (ELSG) Catalog

| Attribute | Details |
|-----------|---------|
| **URL** | https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03090-w |
| **Description** | Genome catalog of early-life human skin microbiome |
| **Content** | 9,483 prokaryotic genomes (1,056 species), 206 fungal genomes (13 species), 39 eukaryotic viral sequences |
| **Improvement** | 21% improvement in sequencing data classification rate |
| **Data Formats** | Genome assemblies, annotations |
| **License** | Open access (BMC) |

---

### 2.5 Expression Atlas (Skin Gene Expression)

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ebi.ac.uk/gxa/home |
| **Description** | Gene expression data for skin and other tissues |
| **Content** | RNA-seq datasets, microarray experiments |
| **Data Formats** | Raw counts (RNA-seq), normalized intensities (microarray) |
| **API** | R package on Bioconductor |
| **License** | Open access |

---

## 3. Vision / Eye Genetics Databases

### 3.1 RetNet - Retinal Information Network

| Attribute | Details |
|-----------|---------|
| **URL** | https://retnet.org/ |
| **Alternative** | https://sph.uth.edu/retnet/ |
| **Maintainer** | Stephen P. Daiger, PhD (University of Texas Health Science Center Houston) |
| **Description** | Comprehensive catalog of mapped and cloned genes causing inherited retinal diseases |
| **Content** | >300 genes implicated in retinal degenerations |
| **Diseases Covered** | Retinitis pigmentosa, macular degeneration, Usher syndrome, Leber congenital amaurosis |
| **Update Frequency** | Weekly review of new research publications |
| **Data Formats** | Web tables with links to external databases |
| **API** | None (web-based reference); cross-references with BioMart, PubMed |
| **License** | Free academic use |
| **Funding** | Foundation Fighting Blindness, George Gund Foundation, Hermann Eye Fund |

---

### 3.2 eyeGENE - National Ophthalmic Disease Genotyping Network

| Attribute | Details |
|-----------|---------|
| **URL** | https://eyegene.nih.gov/ |
| **Data Commons** | https://neidatacommons.nei.nih.gov/eyegene |
| **Maintainer** | National Eye Institute (NEI), NIH |
| **Description** | Genomic medicine initiative for rare inherited eye diseases with DNA biorepository and phenotype/genotype database |
| **Content** | 6,403 participants from 5,385 families; >30 different inherited eye conditions |
| **Common Diagnoses** | Retinitis pigmentosa, Stargardt disease, choroideremia |
| **Top Genes** | ABCA4 (37%), USH2A (7%), RPGR (6%), CHM (5%), PRPH2 (3%) |
| **Diagnostic Rate** | 62.1% with pathogenic/likely pathogenic variants |
| **Data Formats** | Clinical and genetic data (controlled access) |
| **API** | Controlled access via Data Access Committee (DAC) |
| **License** | Data User Agreement (DUA) required |
| **RRID** | SCR_004523 |

---

### 3.3 GEDi Panel (Genetic Eye Disease)

| Attribute | Details |
|-----------|---------|
| **URL** | https://oculargenomics.meei.harvard.edu/services/genetic-diagnostic-testing-service/ |
| **Provider** | Ocular Genomics Institute, Mass. Eye and Ear / Harvard |
| **Description** | Comprehensive genetic diagnostic testing panel for inherited eye diseases |

#### Panel Types:
| Panel | Genes | Conditions |
|-------|-------|------------|
| **GEDi-R** | 267 genes | Inherited retinal degenerations + mitochondrial genome |
| **GEDi-O** | 22 genes | Optic nerve disease, early-onset glaucoma (ACO2, CYP1B1, MYOC, OPA1, etc.) |
| **GEDi-S** | 8 genes | Strabismus, cranial dysinnervation disorders |

| Metric | Value |
|--------|-------|
| **Sensitivity** | 97.9% |
| **Specificity** | 100% |
| **Clinical Sensitivity** | 51% diagnostic rate |
| **Technology** | SureSelect + Illumina MiSeq NGS |

---

### 3.4 RetinoGenetics Database

| Attribute | Details |
|-----------|---------|
| **Description** | Comprehensive mutation database for retinal degeneration genes |
| **Content** | 4,270 mutations in 186 genes; 164 phenotypes from 934 publications |
| **Data** | Functional annotations included |

---

### 3.5 Ocular Genomics Institute Resources

| Attribute | Details |
|-----------|---------|
| **URL** | https://oculargenomics.meei.harvard.edu/ |
| **Genomics Core Services** | Clinical genetic diagnostic testing, WES/WGS, SNP genotyping, CNV analysis |
| **Mission** | Translate genomic medicine into precision ophthalmic care |

---

## 4. Hearing Genetics Databases

### 4.1 Deafness Variation Database (DVD)

| Attribute | Details |
|-----------|---------|
| **URL** | https://deafnessvariationdatabase.org/ |
| **Maintainer** | Molecular Otolaryngology & Renal Research Laboratories (MORL), University of Iowa |
| **Description** | Expert-curated resource of genetic variation in hearing loss genes |
| **Version** | v9.2 (current) |
| **Content** | 876,139 variants across 224 deafness-associated genes |
| **Classifications** | Pathogenic/likely pathogenic: >8,100 variants; Benign/likely benign: >172,000 variants; VUS: >695,000 variants |
| **Platform** | OtoSCOPE genetic screening platform |
| **Data Formats** | Web-based search; downloadable variant data |
| **API** | Web interface queries |
| **License** | Freely available with citation requirement |
| **Copyright** | MORL, University of Iowa |

---

### 4.2 Hereditary Hearing Loss Homepage

| Attribute | Details |
|-----------|---------|
| **URL** | https://hereditaryhearingloss.org/ |
| **Description** | Up-to-date overview of genetics of hereditary hearing impairment |
| **Last Updated** | February 19, 2025 |
| **Content** | 156 nonsyndromic hearing loss genes total |

| Category | Gene Count |
|----------|------------|
| Autosomal Dominant | 64 |
| Autosomal Recessive | 88 |
| Sex-linked | 7 |
| Mitochondrial | 9 |
| Auditory Neuropathy | 5 |

| Additional Stats | Value |
|-----------------|-------|
| **Total NSHL Loci** | 170 |
| **Related Genes** | 124 |

**Data Formats:** Web tables with publication references
**License:** Academic reference use
**ClinGen:** Classifications included

---

### 4.3 Gene4HL

| Attribute | Details |
|-----------|---------|
| **URL** | http://www.genemed.tech/gene4hl/ |
| **Publication** | https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.773009/full |
| **Description** | Integrated genetic database for hearing loss with search, browse, and analysis capabilities |
| **Content** | 326 HL-related genes (170 NSHL + 156 SHL); 3,872 genetic variations; Data from 1,608 published studies |
| **Data Sources** | 62 integrated genetic databases |
| **Features** | Spatiotemporal expression patterns, functional networks, genotype-phenotype correlations |
| **Modules** | Search, Analysis, Gene, Upload, Download |
| **Data Formats** | Web interface; download module available |
| **License** | Open access (Frontiers) |

---

### 4.4 AudioGene

| Attribute | Details |
|-----------|---------|
| **Description** | Machine learning software for predicting genetic cause of hearing loss from audiogram phenotypes |
| **Target** | Autosomal dominant nonsyndromic hearing loss |
| **Accuracy** | 68% for top 3 predictions |
| **Method** | Phenotype-based gene prediction |

---

### 4.5 GeneReviews - Genetic Hearing Loss

| Attribute | Details |
|-----------|---------|
| **URL** | https://www.ncbi.nlm.nih.gov/books/NBK1434/ |
| **Description** | Clinical overview of genetic hearing loss |
| **Key Statistics** | 80% of prelingual hearing loss is genetic; 80% autosomal recessive, 19% autosomal dominant, <1% mitochondrial/X-linked |
| **Top Gene** | GJB2 (connexin 26) - up to 50% of genetic deafness in Europe/USA |

---

## 5. Taste / Smell Genetics Databases

### 5.1 BitterDB

| Attribute | Details |
|-----------|---------|
| **URL** | http://bitterdb.agri.huji.ac.il/dbbitter.php |
| **Maintainer** | Hebrew University of Jerusalem |
| **Description** | Central resource for bitter-tasting molecules and their receptors (TAS2Rs) |
| **Version** | 2024 update |
| **Content** | >2,200 bitter molecules; ~700 molecules with known TAS2R associations; ~1,800 ligand-TAS2R associations |
| **Species** | 66 species (human, chicken, mouse, cat, dog, birds, fishes, primates) |
| **Human Receptors** | 25 functional TAS2R isoforms |
| **Data Formats** | SDF, SMILES, images; compound identifiers (CAS, IUPAC) |
| **API** | Simple search, structure-based similarity search (2D), advanced search |
| **Cross-References** | PubChem, ZINC, DrugBank, IUPHAR/BPS Guide to PHARMACOLOGY |
| **License** | Creative Commons Attribution-NonCommercial (CC BY-NC 4.0) |

**Key Features:**
- Local BLAST for receptor sequence similarity
- Global alignment of 25 human bitter taste receptors
- Natural vs. synthetic compound classification

---

### 5.2 TAS2R38 Genetic Data

| Attribute | Details |
|-----------|---------|
| **Description** | Most studied bitter taste receptor gene |
| **Key Polymorphisms** | A49P, V262A, I296V |
| **Common Haplotypes** | PAV ("taster"), AVI ("non-taster") |
| **Phenotype** | Ability to taste PTC (phenylthiocarbamide) and PROP (6-n-propylthiouracil) |
| **Population Data** | 5,589 individuals from 105 populations (1000 Genomes Project) |
| **BitterDB Ligands** | 37 distinct ligands for T2R38 |

---

### 5.3 ORDB - Olfactory Receptor Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://ordb.biotech.ttu.edu/ORDB/ |
| **Alternative** | http://senselab.med.yale.edu/senselab/ordb |
| **Maintainer** | Texas Tech University / Yale University |
| **Description** | Central repository of olfactory and chemosensory receptor gene/protein sequences |
| **Content** | >18,000 entries; >70 organisms; 35 source tissues |
| **Receptor Types** | Olfactory receptors (ORs), taste papilla receptors (TPRs), vomeronasal receptors (VNRs), insect olfactory receptors (IORs), C. elegans chemosensory receptors, fungal pheromone receptors |
| **Data Formats** | Sequences automatically downloaded from GenBank, SWISS-PROT; EAV/CR database architecture |
| **API** | Search by name, chromosome, or any attribute |
| **License** | Academic use with citation |

**Citation:** Crasto C., Marenco L., Miller P.L., and Shepherd G.S. (2002) Nucleic Acids Research 1:354-360

---

### 5.4 OlfactionBase

| Attribute | Details |
|-----------|---------|
| **URL** | https://bioserver.iiita.ac.in/olfactionbase/ |
| **Description** | Repository for odors, odorants, olfactory receptors, and OR-odorant interactions |
| **Content** | 3,985 odorant compounds; 1,124 odorless compounds; 2,871 OBP/PBP proteins; 408 human odorant-OR pairs; 466 mouse odorant-OR pairs |
| **Human Data** | 197 odorants, 67 receptors, 419 associations |
| **Mouse Data** | 133 ligands, 81 receptors, 518 associations |
| **Features** | Physicochemical and ADMET properties; odor/subodor classification |
| **Data Formats** | Web interface with chemical information (SMILES, CAS, PubChem ID, ZINC ID) |
| **License** | Free, open access |

---

### 5.5 M2OR - Molecule to Olfactory Receptor Database

| Attribute | Details |
|-----------|---------|
| **URL** | https://m2or.fr/ |
| **Publication** | https://academic.oup.com/nar/article/52/D1/D1370/7327067 |
| **Description** | Curated database of OR-molecule interactions including mixtures |
| **Content** | OR responses to molecules and mixtures; receptor sequences; experimental details |
| **Data Formats** | Structured interaction data |
| **License** | Open access |

---

### 5.6 Human Olfactory Receptor Gene Family Statistics

| Metric | Value |
|--------|-------|
| **Total OR Genes** | 636 (339 intact + 297 pseudogenes) |
| **Chromosomal Distribution** | 51 loci on 21 chromosomes |
| **Functional Variation** | 63% of ORs have function-altering polymorphisms |
| **Individual Variation** | ~30% allelic differences between individuals |
| **Estimated Mammalian Range** | ~400 (human) to ~5,000 (elephant) |

---

## 6. Summary Comparison Table

### Oral Health Databases

| Database | URL | Content Size | API | License |
|----------|-----|--------------|-----|---------|
| **HOMD/eHOMD** | homd.org | 836 taxa, genomes | FTP, BLAST | Free academic |
| **HROM** | Cell Host & Microbe | 72,641 genomes | Supplementary | Journal terms |
| **COGR** | Nature npj | 1,089 genomes | Supplementary | Open access |

### Skin Genetics Databases

| Database | URL | Content Size | API | License |
|----------|-----|--------------|-----|---------|
| **GWAS Catalog** | ebi.ac.uk/gwas | 109+ psoriasis loci, 91+ AD loci | REST API | CC0/EMBL-EBI |
| **SKIOME** | GitHub | Aggregated datasets | GitHub | Open access |
| **eSMGC** | iMetaOmics | 675 prokaryotic genomes | - | Open access |
| **ELSG** | Genome Biology | 9,483 genomes | - | Open access |

### Eye Genetics Databases

| Database | URL | Content Size | API | License |
|----------|-----|--------------|-----|---------|
| **RetNet** | retnet.org | >300 genes | Cross-refs | Free academic |
| **eyeGENE** | eyegene.nih.gov | 6,403 participants | Controlled access | DUA required |
| **GEDi Panel** | MEEI/Harvard | 267-297 genes | Clinical service | Clinical |
| **RetinoGenetics** | - | 4,270 mutations, 186 genes | - | Academic |

### Hearing Genetics Databases

| Database | URL | Content Size | API | License |
|----------|-----|--------------|-----|---------|
| **DVD** | deafnessvariationdatabase.org | 876,139 variants, 224 genes | Web queries | Free with citation |
| **HHHL** | hereditaryhearingloss.org | 156 genes, 170 loci | - | Academic |
| **Gene4HL** | genemed.tech/gene4hl | 326 genes, 3,872 variants | Download module | Open access |

### Taste/Smell Databases

| Database | URL | Content Size | API | License |
|----------|-----|--------------|-----|---------|
| **BitterDB** | bitterdb.agri.huji.ac.il | 2,200+ compounds, 66 species | Search, BLAST | CC BY-NC |
| **ORDB** | ordb.biotech.ttu.edu | 18,000+ receptors, 70+ organisms | Search | Academic |
| **OlfactionBase** | bioserver.iiita.ac.in | 3,985 odorants, 2,871 proteins | Web interface | Open access |
| **M2OR** | m2or.fr | OR-molecule pairs | - | Open access |

---

## Data Integration Recommendations

### Priority 1 - Core Databases (Comprehensive, Well-maintained)
1. **HOMD/eHOMD** - Gold standard for oral microbiome
2. **GWAS Catalog** - Central resource for all skin disease genetics
3. **DVD** - Most comprehensive hearing variant database
4. **RetNet** - Authoritative retinal disease gene reference
5. **BitterDB** - Primary taste receptor resource

### Priority 2 - Supplementary Databases
1. **Gene4HL** - Complements DVD with integrated analysis
2. **eyeGENE** - Clinical validation data (requires DUA)
3. **OlfactionBase** - Comprehensive olfactory receptor data
4. **HROM** - Newest oral microbiome genomic catalog

### API Integration Notes
- **REST APIs available**: GWAS Catalog, HOMD (limited)
- **Programmatic access**: GWAS Catalog R package, FTP downloads
- **Controlled access**: eyeGENE (DAC approval required)
- **Web scraping may be needed**: RetNet, HHHL, some others

### License Considerations
- Most databases are free for academic/research use
- BitterDB requires non-commercial attribution (CC BY-NC)
- eyeGENE requires Data User Agreement
- Always cite original publications when using data

---

## References

1. Chen et al. (2010) The Human Oral Microbiome Database. Database (Oxford)
2. HOMD v4 announcement: https://forsyth.org/ada-forsyth-announces-version-4-update/
3. Tsoi LC et al. (2012) Identification of 15 new psoriasis susceptibility loci. Nat Genet
4. Paternoster L et al. (2015) Multi-ancestry GWAS of atopic dermatitis. Nat Genet
5. Shearer AE et al. (2014) Comprehensive genetic testing for hearing loss. Genet Med
6. Daiger SP et al. RetNet: https://retnet.org/
7. Wiener A et al. (2012) BitterDB: A database of bitter compounds. Nucleic Acids Res
8. Crasto C et al. (2002) Olfactory Receptor Database. Nucleic Acids Res

---

## Sample Data

### Example Record: HOMD Oral Microbiome Taxon

```json
{
  "homd_id": "HOT-096",
  "taxon_name": "Streptococcus mutans",
  "phylum": "Firmicutes",
  "genus": "Streptococcus",
  "cultivation_status": "Named",
  "site": "Oral cavity, dental plaque",
  "disease_association": "Dental caries",
  "genome_size_mb": 2.03,
  "gc_content": 36.8
}
```

### Sample Query Result: GWAS Catalog Psoriasis Loci

| rsid | gene | chromosome | p_value | odds_ratio | population | pmid |
|------|------|------------|---------|------------|------------|------|
| rs12191877 | HLA-C | 6 | 4.1e-298 | 4.66 | European | 27723756 |
| rs20541 | IL13 | 5 | 2.5e-15 | 1.18 | European | 27723756 |
| rs240993 | TRAF3IP2 | 6 | 3.8e-32 | 1.35 | European | 27723756 |

### Sample Query Result: Deafness Variation Database

| gene | variant | classification | evidence | condition |
|------|---------|----------------|----------|-----------|
| GJB2 | c.35delG | Pathogenic | Expert panel | Nonsyndromic hearing loss |
| MYO7A | c.6025delG | Pathogenic | Literature | Usher syndrome type 1B |
| SLC26A4 | c.919-2A>G | Pathogenic | Functional | Pendred syndrome |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| HOMD/eHOMD taxa | 836 taxa with genomes |
| HROM genomes | 72,641 high-quality genomes |
| COGR genomes | 1,089 cultivated oral bacteria |
| GWAS Catalog skin loci | 109+ psoriasis, 91+ atopic dermatitis |
| eSMGC skin genomes | 675 prokaryotic + 2,344 viral |
| ELSG early-life genomes | 9,483 prokaryotic genomes |
| DVD hearing variants | 876,139 variants across 224 genes |
| Gene4HL variations | 3,872 variations, 326 genes |
| BitterDB compounds | 2,200+ bitter molecules |
| ORDB entries | 18,000+ olfactory receptors |
| Total storage estimate | ~20-30 GB (combined sources) |
| Last updated | January 2026 |

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **HOMD** | FTP | `http://www.homd.org/ftp/` |
| **HROM** | Web | `https://www.ncbi.nlm.nih.gov/bioproject/PRJNA928779/` |
| **PsoGene** | Web | `https://www.psogene.cn/download.html` |
| **DVD Hearing** | Web | `https://deafnessvariationdatabase.org/` |
| **BitterDB** | Direct | `http://bitterdb.agri.huji.ac.il/dbbitter.php` |
| **ORDB** | Web | `https://senselab.med.yale.edu/ordb/` |

**Access Requirements:** Most databases are open access; some require academic affiliation.

---

## Data Format

| Format | Description | Used By |
|--------|-------------|---------|
| FASTA | Reference sequences | HOMD, HROM |
| TSV | Variant data, gene lists | DVD, Gene4HL |
| CSV | Compound data | BitterDB |
| XML | Detailed records | HOMD |
| JSON | API responses | Various APIs |
| BIOM | Microbiome tables | HOMD, HROM |

**Compression:** gzip (.gz) for sequence data
**Encoding:** UTF-8

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `domain` | string | Health domain | "Oral Health, Skin, Sensory" |
| `database` | string | Data source | "HOMD, GWAS Catalog, RetNet" |
| `category` | string | Data category | "Microbiome, Genetic variants, Gene-disease" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `covers` | Condition | 1:N |
| `sources` | Database | N:M |
| `associates` | Gene | N:M |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `16S rRNA` | Ribosomal RNA gene used for bacterial identification and phylogenetic classification | HOMD 16S RefSeq V16.03 |
| `MAG` | Metagenome-Assembled Genome - genome reconstructed computationally from metagenomic data | HROM 72,641 MAGs |
| `GWAS` | Genome-Wide Association Study - statistical approach identifying genetic variants associated with traits | Psoriasis 109 loci |
| `SNP` | Single Nucleotide Polymorphism - single base pair variation in DNA sequence | TAS2R38 A49P, V262A, I296V |
| `NGS` | Next-Generation Sequencing - high-throughput DNA sequencing technology | GEDi panel uses Illumina MiSeq |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `Phylotype` | Uncultivated microorganism identified by 16S rRNA sequence | HOMD 29% uncultivated |
| `Periodontitis` | Inflammatory disease affecting supporting structures of teeth | 50% heritability |
| `TAS2R` | Taste 2 Receptor - bitter taste receptor gene family | 25 functional human isoforms |
| `Olfactory receptor` | G-protein coupled receptor detecting odorant molecules | 636 human OR genes |
| `Retinitis pigmentosa` | Inherited retinal degeneration causing progressive vision loss | RetNet >300 genes |
| `Connexin 26` | Gap junction protein encoded by GJB2 gene | 50% of genetic deafness |
| `Patescibacteria` | Bacterial superphylum including CPR bacteria | Most prevalent oral phylum |
| `HLA-C*06:02` | HLA allele strongly associated with psoriasis risk | Key psoriasis susceptibility |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| HOMD | Human Oral Microbiome Database | 836 taxa, NIH HMP |
| eHOMD | Expanded HOMD | Current version |
| HROM | Human Reference Oral Microbiome | 72,641 genomes |
| COGR | Cultivated Oral Bacteria Genome Reference | 1,089 genomes |
| DVD | Deafness Variation Database | 876,139 variants |
| HHHL | Hereditary Hearing Loss Homepage | 156 genes |
| RetNet | Retinal Information Network | >300 retinal disease genes |
| eyeGENE | National Ophthalmic Disease Genotyping Network | NEI/NIH resource |
| GEDi | Genetic Eye Disease panel | Harvard/MEEI |
| ORDB | Olfactory Receptor Database | 18,000+ entries |
| TAS2R38 | Taste 2 Receptor Member 38 | PTC/PROP tasting gene |
| CPR | Candidate Phyla Radiation | Uncultivated bacteria |
| MORL | Molecular Otolaryngology & Renal Research Labs | University of Iowa |
| NEI | National Eye Institute | Part of NIH |
| CC0 | Creative Commons Zero (public domain) | GWAS Catalog |
| CC BY-NC 4.0 | Creative Commons Attribution-NonCommercial 4.0 | BitterDB license |
| DUA | Data User Agreement | eyeGENE access |

---

## References

1. Chen et al. (2010) The Human Oral Microbiome Database. Database (Oxford)
2. HOMD v4 announcement: https://forsyth.org/ada-forsyth-announces-version-4-update/
3. Tsoi LC et al. (2012) Identification of 15 new psoriasis susceptibility loci. Nat Genet
4. Paternoster L et al. (2015) Multi-ancestry GWAS of atopic dermatitis. Nat Genet
5. Shearer AE et al. (2014) Comprehensive genetic testing for hearing loss. Genet Med
6. Daiger SP et al. RetNet: https://retnet.org/
7. Wiener A et al. (2012) BitterDB: A database of bitter compounds. Nucleic Acids Res
8. Crasto C et al. (2002) Olfactory Receptor Database. Nucleic Acids Res

---

*Document compiled: January 2025*
*For Gene Platform data integration planning*
