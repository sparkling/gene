---
id: genetics-population
title: Population Genetics Databases
world: 1
category: genetics
tier: 2
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [genetics, population, diversity]
---

# Expanded Genetics & SNP Data Sources

**Last Updated:** January 2026
**Purpose:** New databases NOT already in data-sources-public.md
**Constraint:** Publicly available sources with open or academic licenses
**Parent:** [../_index.md](../_index.md)

---

## Summary of New Sources

| Category | New Sources | Key Highlights |
|----------|-------------|----------------|
| Population Genetics | 6 | TOPMed BRAVO, All of Us, ALFA, UK Biobank AFB, GenomeAsia 100K, H3Africa |
| Ancestry-Specific | 4 | African Genome Variation Project, AGenDA, H3ABioNet, GenomeAsia 100K |
| Functional Annotation | 8 | AlphaMissense, dbNSFP, EVE, MaveDB, PrimateAI-3D, RegulomeDB, SpliceAI, CADD |
| Epigenetics | 6 | ENCODE 4, IHEC, MethBank, 4D Nucleome, Roadmap Epigenomics, FANTOM5 |
| Structural Variants | 6 | gnomAD-SV v4.1, DGV, dbVar, DECIPHER, SV4GD, HGSVC |

---

## 1. Population Genetics Databases (Beyond 1000 Genomes, gnomAD)

### 1.1 TOPMed BRAVO
| Field | Value |
|-------|-------|
| **URL** | https://bravo.sph.umich.edu/ |
| **Content** | 868M+ variants from 150,000+ whole genomes |
| **Populations** | ~60% non-European ancestry |
| **Freeze** | Freeze 10 (current) |
| **API** | REST API via bravo_api (https://github.com/statgen/bravo_api) |
| **License** | Public (summary data), controlled access (individual) |
| **Access** | Web browser, BioData Catalyst for controlled data |
| **Size** | Summary data ~10 GB |
| **Note** | Best for rare variants in diverse populations |

### 1.2 All of Us Researcher Workbench
| Field | Value |
|-------|-------|
| **URL** | https://www.researchallofus.org/ |
| **Content** | 1B+ genetic variants from 245K+ individuals |
| **Novel Variants** | 275M+ previously unreported variants |
| **Diversity** | ~50% non-European ancestry |
| **API** | Researcher Workbench APIs (controlled tier) |
| **Formats** | VDS, Hail MT, VCF, BGEN, PLINK |
| **License** | Data Use Agreement required |
| **Access** | Researcher Workbench (free registration, $300 credits) |
| **Size** | ~50 TB (full), summary data available publicly |
| **Updates** | Fall 2025: Verily Workbench integration |

### 1.3 ALFA (Allele Frequency Aggregator) - NCBI
| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/ |
| **Content** | Allele frequencies from 200K+ dbGaP subjects |
| **Version** | Release 4 (June 2025) - doubled cohort size from R3 |
| **Populations** | 12 populations (European, African, Asian, Latin American, etc.) |
| **API** | E-utilities API, Variation Services API |
| **Formats** | VCF, JSON via API |
| **License** | Public domain |
| **FTP** | https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/ |
| **Size** | ~5 GB |
| **ClinVar Coverage** | 960K+ ClinVar RS IDs (74% increase from R3) |

### 1.4 UK Biobank Allele Frequency Browser (AFB)
| Field | Value |
|-------|-------|
| **URL** | https://afb.ukbiobank.ac.uk/ |
| **Content** | SNP/indel frequencies from 150,119 WGS individuals |
| **Full Dataset** | 490,640 WGS (Sept 2025 Nature publication) |
| **Browser** | Global Biobank Engine (https://biobankengine.stanford.edu/) |
| **API** | Limited public access, full via Research Analysis Platform |
| **License** | Public (AFB), controlled (full data) |
| **Access** | AFB publicly available, full data requires application |
| **Size** | ~30 TB (full), AFB summary ~10 GB |

### 1.5 GenomeAsia 100K
| Field | Value |
|-------|-------|
| **URL** | https://browser.genomeasia100k.org/ |
| **Content** | 1,739 WGS from 219 Asian population groups |
| **Countries** | 64 countries across Asia |
| **Variants** | 66M+ catalogued variants |
| **API** | Data access via application |
| **Formats** | VCF (individual-level) |
| **License** | Open for research (data access agreement) |
| **Contact** | dataaccess@genomeasia100k.org |
| **Size** | ~2 TB (full VCFs) |

### 1.6 H3Africa / African Genome Variation Project
| Field | Value |
|-------|-------|
| **URL** | https://h3africa.org/ |
| **Content** | WGS from 426+ individuals, 50 ethnolinguistic groups |
| **AGVP** | 1,481 dense genotypes + 320 WGS across sub-Saharan Africa |
| **Novel Variants** | 3M+ previously undescribed variants |
| **Repository** | EGA (EGAS00001002976) |
| **API** | H3ABioNet tools (https://github.com/h3abionet/h3agwas) |
| **License** | Controlled access (bona fide researchers) |
| **Contact** | dbac@h3africa.org |
| **Array** | H3Africa 2.3M SNP array available |

---

## 2. Ancestry-Specific Variant Databases

### 2.1 AGenDA (Assessing Genetic Diversity in Africa)
| Field | Value |
|-------|-------|
| **URL** | Part of H3Africa consortium |
| **Content** | WGS from under-represented African groups |
| **Countries** | 9 African countries |
| **Publication** | Nature 2025 (10.1038/s41586-025-09935-7) |
| **Purpose** | Enriching global datasets with African diversity |
| **License** | Controlled access |
| **Note** | Addresses gaps in existing reference panels |

### 2.2 gnomAD Ancestry-Stratified Data
| Field | Value |
|-------|-------|
| **URL** | https://gnomad.broadinstitute.org/ |
| **Populations** | African/African American, Admixed American, Ashkenazi Jewish, East Asian, European (Finnish), European (non-Finnish), Middle Eastern, South Asian |
| **Version** | v4.1 (807,162 individuals) |
| **API** | GraphQL API |
| **License** | Open access |
| **Note** | Use for ancestry-specific allele frequency filtering |

### 2.3 jMorp (Japanese Multi-Omics Reference Panel)
| Field | Value |
|-------|-------|
| **URL** | https://jmorp.megabank.tohoku.ac.jp/ |
| **Content** | Japanese population-specific variants |
| **Samples** | 38K+ Japanese individuals |
| **Data Types** | WGS, WES, proteomics, metabolomics |
| **API** | REST API available |
| **License** | Open access (summary), controlled (individual) |

---

## 3. Functional Annotation Databases

### 3.1 AlphaMissense (Google DeepMind)
| Field | Value |
|-------|-------|
| **URL** | https://alphamissense.hegelab.org/ |
| **Content** | Pathogenicity scores for 71M missense variants |
| **Classification** | 57% likely benign, 32% likely pathogenic, 11% uncertain |
| **Download** | Google Cloud Public Dataset |
| **R Package** | AlphaMissenseR (Bioconductor) |
| **Formats** | TSV, VEP plugin compatible |
| **License** | CC BY 4.0 (academic and commercial use) |
| **GitHub** | https://github.com/google-deepmind/alphamissense |
| **Size** | ~10 GB (all predictions) |
| **Integration** | Ensembl VEP plugin available |

### 3.2 dbNSFP v4.9
| Field | Value |
|-------|-------|
| **URL** | https://www.dbnsfp.org/home |
| **Content** | 83M nsSNVs + 2.4M ssSNVs with 35 prediction scores |
| **Scores Included** | SIFT, PolyPhen2, CADD, REVEL, AlphaMissense, SpliceAI, EVE, PrimateAI, etc. |
| **Conservation** | PhyloP, phastCons, GERP++, bStatistic |
| **Version** | v4.9a (academic), v4.9c (commercial) |
| **Download** | Amazon S3, Box |
| **Formats** | TSV, VCF-ready |
| **License** | Academic (v4.9a), Commercial (v4.9c - excludes some scores) |
| **Size** | ~35 GB (compressed) |
| **Updates** | October 2025 (GENCODE 49, Ensembl 115) |

### 3.3 CADD (Combined Annotation Dependent Depletion)
| Field | Value |
|-------|-------|
| **URL** | https://cadd.gs.washington.edu/ |
| **Content** | Deleteriousness scores for all possible SNVs |
| **Version** | v1.7 (current) |
| **Scores** | Raw and PHRED-scaled C-scores |
| **API** | Web API, offline scoring available |
| **Download** | Pre-scored files for coding/non-coding |
| **License** | Academic use free, commercial license available |
| **Size** | ~300 GB (genome-wide), ~15 GB (coding only) |
| **Integration** | VEP plugin, dbNSFP included |

### 3.4 SpliceAI
| Field | Value |
|-------|-------|
| **URL** | https://spliceailookup.broadinstitute.org/ |
| **Content** | Splice-altering predictions for all SNVs/indels |
| **Scores** | DS_AG, DS_AL, DS_DG, DS_DL (0-1 scale) |
| **Download** | Illumina BaseSpace (https://basespace.illumina.com/s/otSPW8hnhaZR) |
| **PyPI** | `pip install spliceai` |
| **Thresholds** | 0.2 (high recall), 0.5 (recommended), 0.8 (high precision) |
| **License** | Free for academic, commercial license for Illumina |
| **Size** | ~50 GB (pre-computed scores) |
| **Note** | 2025 update: SpliceAI-splint tool addresses precomputed score limitations |

### 3.5 REVEL
| Field | Value |
|-------|-------|
| **URL** | https://sites.google.com/site/revelgenomics/ |
| **Content** | Ensemble missense pathogenicity scores |
| **Components** | MutPred, FATHMM, VEST, PolyPhen, SIFT, PROVEAN, etc. |
| **Range** | 0-1 (higher = more pathogenic) |
| **Download** | Direct download available |
| **License** | Academic use |
| **Size** | ~5 GB |
| **Integration** | dbNSFP, VEP plugin |

### 3.6 EVE (Evolutionary model of Variant Effect)
| Field | Value |
|-------|-------|
| **URL** | https://evemodel.org/ |
| **Content** | Pathogenicity predictions for 36M+ variants across 3,219 genes |
| **VUS Classified** | 256K+ variants of unknown significance |
| **Method** | Bayesian VAE trained on evolutionary data |
| **Download** | evemodel.org, Supplementary Data |
| **GitHub** | https://github.com/OATML-Markslab/EVE |
| **License** | Academic use |
| **Publication** | Nature 2021 |

### 3.7 PrimateAI-3D
| Field | Value |
|-------|-------|
| **URL** | https://primateai3d.basespace.illumina.com/download |
| **Content** | Pathogenicity scores using primate variation + AlphaFold structures |
| **Method** | Semi-supervised 3D-CNN trained on 4.5M benign variants |
| **Download** | Illumina BaseSpace |
| **GitHub** | https://github.com/Illumina/PrimateAI |
| **License** | Free academic, commercial license from Illumina |
| **Thresholds** | >0.8 likely pathogenic, <0.6 likely benign |

### 3.8 MaveDB
| Field | Value |
|-------|-------|
| **URL** | https://www.mavedb.org/ |
| **Content** | 7M+ variant effect measurements from MAVE experiments |
| **Datasets** | 1,884 datasets (Nov 2024) |
| **Data Types** | DMS, MPRA, saturation genome editing |
| **API** | REST API (JSON format) |
| **Download** | CSV per score set |
| **GitHub** | https://github.com/VariantEffect/mavedb-api |
| **License** | Open access |
| **Publication** | Genome Biology Jan 2025 |

### 3.9 RegulomeDB v2
| Field | Value |
|-------|-------|
| **URL** | http://regulomedb.org |
| **Content** | Regulatory potential scores for non-coding variants |
| **Data Sources** | ENCODE, Roadmap Epigenomics, GTEx |
| **Scoring** | 1 (strongest) to 6 (weakest) regulatory evidence |
| **Intervals** | 650M (hg19), 1.5B (GRCh38) |
| **Builds** | hg19 and GRCh38 |
| **License** | Open access |
| **Integration** | SURF algorithm scores included |

---

## 4. Epigenetics Databases

### 4.1 ENCODE 4
| Field | Value |
|-------|-------|
| **URL** | https://www.encodeproject.org/ |
| **Content** | 926,535 human cCREs, 339,815 mouse cCREs |
| **Data Types** | ChIP-seq (H3K4me1/3, H3K27ac/me3, H3K36me3, H3K9me3, CTCF), DNase-seq, ATAC-seq, RNA-seq |
| **Biosamples** | 234 human biosamples, 1,794 experiments |
| **API** | REST API |
| **Download** | FTP, rsync (recommended) |
| **Browser** | SCREEN web interface |
| **License** | Open access |
| **Size** | ~5 TB (full), ~50 GB (processed) |

### 4.2 Roadmap Epigenomics
| Field | Value |
|-------|-------|
| **URL** | https://egg2.wustl.edu/roadmap/web_portal/ |
| **Content** | 111 reference human epigenomes |
| **Data Types** | Histone ChIP-seq (H3K4me1/3, H3K27me3, H3K36me3, H3K9me3, H3K27ac, H3K9ac), DNA methylation, RNA-seq |
| **Methylation** | 104 datasets (WGBS, RRBS, MeDIP-seq, MRE-seq) |
| **Download** | http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics |
| **Human Epigenome Atlas** | Release 9 |
| **License** | Open access |
| **Size** | ~2 TB |

### 4.3 IHEC Data Portal
| Field | Value |
|-------|-------|
| **URL** | https://epigenomesportal.ca/ihec/ |
| **Content** | 7,500+ reference epigenomic datasets |
| **Tissues** | 600+ tissue types |
| **Contributors** | ENCODE, Roadmap, CEEHRC, Blueprint, DEEP, AMED-CREST, KNIH |
| **Browser** | UCSC mirror with preloaded tracks |
| **Download** | Bulk download available |
| **License** | Open access |
| **Size** | ~10 TB (full) |

### 4.4 MethBank 4.0
| Field | Value |
|-------|-------|
| **URL** | https://ngdc.cncb.ac.cn/methbank/ |
| **Content** | 1,449 single-base resolution methylomes |
| **Species** | 23 species (human primary) |
| **Human Data** | 53,680 age-specific DMCs, 1,716 age-specific DMRs |
| **Tools** | DMR Toolkit for analysis |
| **Download** | Web interface, bulk download |
| **License** | Open access |
| **Updates** | Aug 2025: 30 new human projects |
| **Cancer Module** | Online Jan 2025 |

### 4.5 4D Nucleome Data Portal
| Field | Value |
|-------|-------|
| **URL** | https://data.4dnucleome.org/ |
| **Content** | 1,800+ experiment sets, 36,000+ files |
| **Data Types** | Hi-C, ChIA-PET, Capture Hi-C, PLAC-Seq |
| **High-Resolution** | 40+ datasets with >1B read pairs |
| **Visualization** | HiGlass, 3D Genome Browser, WashU Browser |
| **Download** | Bulk download, filtering interface |
| **License** | Open access |
| **Publication** | Nature 2025 (final consortium paper) |

### 4.6 FANTOM5
| Field | Value |
|-------|-------|
| **URL** | https://fantom.gsc.riken.jp/5/ |
| **Content** | 210K human promoters, 63K human enhancers |
| **Samples** | 1,800 human CAGE profiles, 1,000 mouse profiles |
| **Data Types** | CAGE (cap analysis gene expression) |
| **Tools** | BioMart, UCSC DataHub, Enhancer Selector |
| **Download** | https://dbarchive.biosciencedbc.jp/data/fantom5/ |
| **License** | CC BY 4.0 |
| **Update** | NAR 2025 update for lncRNA interfaces |

---

## 5. Structural Variant Databases

### 5.1 gnomAD-SV v4.1
| Field | Value |
|-------|-------|
| **URL** | https://gnomad.broadinstitute.org/ |
| **Content** | 1.2M SVs from 464,297 exomes + 63,046 genomes |
| **SV Types** | Deletions, duplications, insertions, inversions, mCNVs, translocations, complex |
| **Per-Person** | ~11,844 SVs average |
| **Detection** | GATK-SV (genomes), GATK-gCNV (exomes) |
| **API** | GraphQL API |
| **Download** | BED files, Hail Tables (AWS) |
| **License** | Open access |
| **Size** | ~5 GB (summary), ~500 GB (full) |

### 5.2 Database of Genomic Variants (DGV)
| Field | Value |
|-------|-------|
| **URL** | http://dgv.tcag.ca/ |
| **Content** | 2.5M+ SV entries from 55 studies |
| **SV Types** | CNVs >50 bp, inversions |
| **Population** | Control/healthy individuals |
| **Download** | FTP, web interface |
| **Formats** | GFF, BED, VCF |
| **License** | Open access |
| **Size** | ~500 MB |

### 5.3 dbVar (NCBI)
| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/dbvar/ |
| **Content** | 6M+ structural variants from 185+ studies |
| **Sources** | 1000 Genomes, gnomAD, ClinVar, ClinGen |
| **Common SVs** | nstd186 (curated common SVs) |
| **API** | E-utilities |
| **Download** | FTP |
| **License** | Public domain |
| **Size** | ~10 GB |

### 5.4 DECIPHER
| Field | Value |
|-------|-------|
| **URL** | https://www.deciphergenomics.org/ |
| **Content** | 51,894 patients, 51K+ variants, 172K+ phenotype terms |
| **SV Types** | CNVs, aneuploidy, UPD, inversions, insertions, STRs |
| **Phenotypes** | HPO terms with quantitative data |
| **Projects** | 250+ projects from 40 countries |
| **Download** | Bulk download (data access agreement) |
| **License** | Open (consented data), Data Display Agreement for third parties |
| **Integration** | Ensembl, UCSC display |

### 5.5 SV4GD (Structural Variation for Genetic Diseases)
| Field | Value |
|-------|-------|
| **URL** | https://sv4gd.com/ (NAR 2025) |
| **Content** | 10,305 germline SV records |
| **Diseases** | 58 neoplastic, 232 non-neoplastic genetic diseases |
| **Disease-Related** | 2,695 disease-related SVs |
| **Common Types** | Deletions (1,462), duplications (743) = 81.8% |
| **Curation** | Manual from ~300 publications |
| **License** | Open access |
| **Publication** | NAR Jan 2025 |

### 5.6 HGSVC (Human Genome Structural Variation Consortium)
| Field | Value |
|-------|-------|
| **URL** | https://www.hgsvc.org/ |
| **Content** | Complete haplotype sequences from 65 diverse individuals |
| **HGSVC3** | Near T2T completeness, fully resolved centromeres |
| **Data Types** | PacBio, Strand-Seq, Bionano maps |
| **Download** | https://www.internationalgenome.org/data-portal/data-collection/hgsvc2 |
| **GitHub** | https://github.com/hgsvc |
| **License** | Open access |
| **Publications** | Nature Aug 2025, Nature Jul 2025 |

### 5.7 ClinGen Dosage Sensitivity
| Field | Value |
|-------|-------|
| **URL** | https://dosage.clinicalgenome.org/ |
| **Content** | Haploinsufficiency and triplosensitivity scores |
| **Download** | https://search.clinicalgenome.org/kb/downloads |
| **Formats** | CSV, BED, TSV |
| **API** | REST API available |
| **License** | Open access |
| **Updates** | Daily |

---

## 6. Integration Tools & Annotation Pipelines

### 6.1 Ensembl VEP (Variant Effect Predictor)
| Field | Value |
|-------|-------|
| **URL** | https://www.ensembl.org/vep |
| **Version** | Release 115 (Sept 2025) |
| **Plugins** | AlphaMissense, SpliceAI, CADD, MaveDB, UTRAnnotator, NMD, OpenTargets |
| **API** | REST API |
| **Installation** | Conda, Docker, web interface |
| **License** | Open source (Apache 2.0) |
| **Updates 2025** | gnomAD-SV allele frequencies, ClinVar somatic classifications, All of Us frequencies |

### 6.2 ANNOVAR
| Field | Value |
|-------|-------|
| **URL** | https://annovar.openbioinformatics.org/ |
| **Content** | Gene-based, region-based, filter-based annotation |
| **Databases** | gnomAD, dbSNP, ClinVar, dbNSFP |
| **License** | Free for academic use |
| **Updates** | March 2024 version |

---

## 7. Recommended Ingestion Priority

### Tier 1 - Critical (Weeks 1-2)
| Database | Priority Reason |
|----------|-----------------|
| **dbNSFP v4.9** | 35 scores in one download, comprehensive |
| **AlphaMissense** | State-of-art missense predictions, CC BY 4.0 |
| **gnomAD-SV v4.1** | Population SV frequencies |
| **ALFA R4** | NCBI-integrated, public domain |

### Tier 2 - High Value (Weeks 3-4)
| Database | Priority Reason |
|----------|-----------------|
| **TOPMed BRAVO** | Diverse populations, rare variants |
| **ENCODE 4** | Regulatory element annotations |
| **MaveDB** | Experimental variant effects |
| **ClinGen Dosage** | Copy number interpretation |

### Tier 3 - Enrichment (Weeks 5-8)
| Database | Priority Reason |
|----------|-----------------|
| **GenomeAsia 100K** | Asian-specific variants |
| **H3Africa/AGVP** | African diversity |
| **4D Nucleome** | 3D genome context |
| **FANTOM5** | Enhancer/promoter annotations |

---

## 8. API Summary

| Database | API Type | Rate Limits | Authentication |
|----------|----------|-------------|----------------|
| gnomAD | GraphQL | None stated | None |
| TOPMed BRAVO | REST | Unknown | None (summary) |
| All of Us | Workbench APIs | Per-project | Registration required |
| ALFA | E-utilities | 3/sec (key), 10/sec (API key) | Optional API key |
| ENCODE | REST | None stated | None |
| MaveDB | REST | None stated | None |
| Ensembl VEP | REST | 15 req/sec | None |
| ClinGen | REST/FTP | None stated | None |

---

## 9. License Comparison

| Database | License | Commercial Use | Attribution |
|----------|---------|----------------|-------------|
| AlphaMissense | CC BY 4.0 | Yes | Required |
| dbNSFP | Academic/Commercial split | v4.9c only | Required |
| gnomAD | Open access | Yes | Citation |
| ENCODE | Open access | Yes | Citation |
| FANTOM5 | CC BY 4.0 | Yes | Required |
| MaveDB | Open access | Yes | Citation |
| All of Us | DUA | Research only | Required |
| TOPMed | Open/Controlled | Summary yes | Citation |

---

## 10. Size Estimates for Full Integration

| Category | Compressed Size | Uncompressed | Notes |
|----------|-----------------|--------------|-------|
| Population genetics | ~100 GB | ~500 GB | Summary-level only |
| Functional annotation | ~60 GB | ~300 GB | dbNSFP largest |
| Epigenetics | ~20 GB | ~100 GB | ENCODE cCREs focus |
| Structural variants | ~15 GB | ~75 GB | gnomAD-SV primary |
| **Total** | **~195 GB** | **~975 GB** | Summary/processed data |

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **gnomAD v4.1** | Download | `https://gnomad.broadinstitute.org/downloads` |
| **TOPMed BRAVO** | Web browser | `https://bravo.sph.umich.edu/` |
| **ALFA** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/snp/population_frequency/` |
| **1000 Genomes** | FTP | `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/` |
| **dbSNP** | FTP | `ftp://ftp.ncbi.nih.gov/snp/` |
| **ClinVar** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/` |

**Access Requirements:** Most are open access; All of Us requires Data Use Agreement; TOPMed individual data requires dbGaP application.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | VCF, TSV, JSON |
| Alternative | Parquet, BigQuery |
| Variant notation | rsID, HGVS, VCF format |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `rsid` | string | Reference SNP ID | "rs1801133" |
| `chrom` | string | Chromosome | "1" |
| `pos` | integer | Genomic position | 11856378 |
| `ref` | string | Reference allele | "G" |
| `alt` | string | Alternate allele | "A" |
| `af` | float | Global allele frequency | 0.35 |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `in_gene` | Gene | N:1 |
| `has_annotation` | Functional prediction | 1:N |
| `associated_with` | Phenotype | N:M |

## Sample Data

### Example Variant Record
```json
{
  "rsid": "rs1801133",
  "chrom": "1",
  "pos": 11856378,
  "ref": "G",
  "alt": "A",
  "gene": "MTHFR",
  "consequence": "missense_variant",
  "af_global": 0.35,
  "af_by_pop": {
    "EUR": 0.36,
    "AFR": 0.11,
    "EAS": 0.44,
    "SAS": 0.12
  }
}
```

### Sample Query Result
| rsid | gene | consequence | af_global |
|------|------|-------------|-----------|
| rs1801133 | MTHFR | missense | 0.35 |
| rs1801131 | MTHFR | missense | 0.29 |

## Data Set Size

| Metric | Value |
|--------|-------|
| gnomAD v4.1 variants | 786M+ variants from 807K individuals |
| TOPMed variants | 868M+ variants from 150K genomes |
| ALFA variants | 1B+ frequency records |
| dbSNP variants | 1.1B+ submitted variants |
| Total storage estimate | ~195 GB compressed |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `freeze` | A versioned snapshot of sequencing consortium data representing a stable release | TOPMed Freeze 10 |
| `ancestry stratification` | Division of genetic data by population or continental ancestry groups | European, African, East Asian |
| `allele frequency` | The proportion of a specific allele in a population | 0.15 (15% of alleles) |
| `imputation` | Statistical inference of unobserved genotypes using reference panel haplotypes | Imputing rare variants from array data |
| `reference panel` | A collection of haplotypes used for imputation or population frequency comparison | gnomAD v4.1 reference |
| `novel variant` | A genetic variant not previously catalogued in reference databases | Variant not in dbSNP |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `gnomAD` | Genome Aggregation Database with 807K+ individuals across 8 ancestry groups | Population frequencies |
| `TOPMed BRAVO` | Browser for TOPMed allele frequencies with 868M+ variants from 150K genomes | Rare variant frequencies |
| `ALFA` | NCBI Allele Frequency Aggregator from 200K+ dbGaP subjects | NCBI integration |
| `All of Us` | NIH precision medicine cohort with 50% non-European ancestry | Diverse populations |
| `H3Africa` | Human Heredity and Health in Africa consortium addressing African diversity gaps | African genetics |
| `GenomeAsia 100K` | Asian population reference with 219 population groups across 64 countries | Asian diversity |
| `jMorp` | Japanese Multi-Omics Reference Panel with 38K+ Japanese individuals | Japanese population |
| `AGVP` | African Genome Variation Project with 1,481 individuals from sub-Saharan Africa | African reference |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| WGS | Whole Genome Sequencing | Complete genome at 30x+ coverage |
| WES | Whole Exome Sequencing | Coding regions only (~2% of genome) |
| VDS | Variant Dataset | Hail's columnar variant storage format |
| AFB | Allele Frequency Browser | UK Biobank public frequency tool |
| EGA | European Genome-phenome Archive | Controlled access genetic data repository |
| dbGaP | Database of Genotypes and Phenotypes | NCBI controlled access repository |
| MAF | Minor Allele Frequency | Less common allele's frequency |
| SNP | Single Nucleotide Polymorphism | Point variant in the genome |
| SV | Structural Variant | Large genomic rearrangements (>50 bp) |
| CNV | Copy Number Variant | Deletions or duplications |
| DMS | Deep Mutational Scanning | Experimental saturation mutagenesis |
| MPRA | Massively Parallel Reporter Assay | High-throughput regulatory element testing |

---

## License

This document catalogs multiple databases with varying license terms:

| Database | License | Commercial Use | Attribution | Access |
|----------|---------|----------------|-------------|--------|
| TOPMed BRAVO | Public (summary), Controlled (individual) | Summary data: Yes | Citation | Open (summary) |
| All of Us | Data Use Agreement | Research only | Required | Registration required |
| ALFA R4 | Public Domain | Yes | Citation | Open |
| UK Biobank AFB | Public (AFB), Controlled (full) | Summary data: Yes | Required | AFB open |
| GenomeAsia 100K | Data Access Agreement | Research only | Required | Application required |
| H3Africa/AGVP | Controlled Access | Research only | Required | Application required |
| AGenDA | Controlled Access | Research only | Required | Application required |
| gnomAD v4.1 | Open Access | Yes | Citation | Open |
| jMorp | Open (summary), Controlled (individual) | Summary data: Yes | Citation | Open (summary) |
| AlphaMissense | CC BY 4.0 | Yes | Required | Open |
| dbNSFP v4.9 | Academic (v4.9a), Commercial (v4.9c) | v4.9c only | Required | Download |
| CADD v1.7 | Academic free, Commercial available | License required | Required | Download |
| SpliceAI | Academic free, Commercial from Illumina | License required | Required | Download |
| REVEL | Academic use | No | Required | Download |
| EVE | Academic use | No | Required | Download |
| PrimateAI-3D | Academic free, Commercial from Illumina | License required | Required | Download |
| MaveDB | Open Access | Yes | Citation | Open |
| RegulomeDB | Open Access | Yes | Citation | Open |
| ENCODE 4 | Open Access | Yes | Citation | Open |
| Roadmap Epigenomics | Open Access | Yes | Citation | Open |
| IHEC | Open Access | Yes | Citation | Open |
| MethBank 4.0 | Open Access | Yes | Citation | Open |
| 4D Nucleome | Open Access | Yes | Citation | Open |
| FANTOM5 | CC BY 4.0 | Yes | Required | Open |
| gnomAD-SV v4.1 | Open Access | Yes | Citation | Open |
| DGV | Open Access | Yes | Citation | Open |
| dbVar | Public Domain | Yes | Citation | Open |
| DECIPHER | Open (consented), Data Display Agreement | Agreement required | Required | Agreement required |
| SV4GD | Open Access | Yes | Citation | Open |
| HGSVC | Open Access | Yes | Citation | Open |
| ClinGen Dosage | Open Access | Yes | Citation | Open |
| VEP | Apache 2.0 | Yes | Required | Open Source |
| ANNOVAR | Academic free | No | Required | Download |

**Key Considerations:**
- **Fully Open (Commercial OK):** ALFA, AlphaMissense, FANTOM5, gnomAD, MaveDB, dbVar
- **Academic Only:** dbNSFP (v4.9a), CADD, REVEL, EVE, ANNOVAR
- **Controlled Access:** TOPMed (individual), All of Us, GenomeAsia 100K, H3Africa, AGenDA
- **Commercial License Required:** SpliceAI, PrimateAI-3D (Illumina)

---

## References

1. gnomAD v4.1: https://gnomad.broadinstitute.org/news/2025-03-cnv-sv-hts/
2. TOPMed BRAVO: https://bravo.sph.umich.edu/
3. All of Us: https://www.nature.com/articles/s41586-023-06957-x
4. ALFA R4: https://ncbiinsights.ncbi.nlm.nih.gov/2025/06/09/ncbi-alfa-release-4/
5. AlphaMissense: https://www.science.org/doi/10.1126/science.adg7492
6. dbNSFP v4: https://link.springer.com/article/10.1186/s13073-020-00803-9
7. ENCODE 4: https://www.nature.com/articles/s41586-020-2493-4
8. MaveDB 2024: https://link.springer.com/article/10.1186/s13059-025-03476-y
9. SV4GD: https://academic.oup.com/nar/article/53/D1/D1557/7889244
10. 4D Nucleome: https://www.nature.com/articles/s41586-025-09890-3
