# Microbiome and Gut Health Databases for Gene Platform

## Research Summary

This document catalogs key microbiome and gut health databases relevant to the Gene Platform, including data access methods, licensing, and integration considerations.

---

## 1. Human Microbiome Project (HMP)

### Overview
The NIH-funded Human Microbiome Project is a collaborative effort of over 300 scientists from 80+ organizations to comprehensively characterize microbial communities inhabiting the human body and their role in health and disease.

### URLs
- **Main Portal**: https://hmpdacc.org/
- **Data Portal**: https://portal.hmpdacc.org/
- **AWS Open Data**: https://registry.opendata.aws/human-microbiome-project/
- **GitHub (Portal API)**: https://github.com/IGS/portal-api
- **GitHub (Organization)**: https://github.com/ihmpdcc

### Content
- **Samples**: 11,000+ physical samples from 300 healthy adult subjects
- **Body Sites**: 18 specific sites across 5 body regions (oral cavity, airways, urogenital tract, skin, gut)
- **Sequencing**: 16S rRNA sequencing for all samples; ~2,300 samples with whole metagenomic shotgun sequencing
- **Data Types**: 16S, metagenomic WGS, transcriptome, proteome, metabolome (host and microbiome), host whole genomes, epigenomes

### API Access
| Method | Description |
|--------|-------------|
| **HMP Portal API** | Flask app using GraphQL with Neo4j backend; RESTful data transfer |
| **GA4GH DRS API** | Data Repository Service implementation for standardized access |
| **HMP16SData (Bioconductor)** | R package for programmatic access to 16S data with integrated phylogeny and taxonomy |
| **AWS S3** | Direct download via AWS Open Data Registry |

### License
- **Type**: Creative Commons Attribution (CC BY 4.0)
- **Terms**: Unrestricted reuse, distribution, and reproduction with proper citation
- **Controlled Data**: Some data types require dbGaP authorization

### Data Size
- **Total**: >48 TB (combined HMP and iHMP phases)
- **Public vs Controlled**: Both public and controlled access datasets

### Status
- **Active Support**: 2007-2016 (NIH Common Fund)
- **Current**: Maintained as archive; data still accessible

### Gene Platform Integration Notes
- Excellent for baseline healthy microbiome reference
- Multi-omic data enables correlation with host genetics
- AWS hosting facilitates cloud-based pipeline integration

---

## 2. MetaHIT (Metagenomics of the Human Intestinal Tract)

### Overview
European Commission-funded collaborative project (15 institutes, 8 countries) focused on establishing associations between human intestinal microbiota genes and health/disease states.

### URLs
- **ENA Browser**: https://www.ebi.ac.uk/ena/browser/view/PRJNA32811
- **MAHMI Database**: https://academic.oup.com/database/article/doi/10.1093/database/baw157/2884895
- **Project Info (CORDIS)**: https://cordis.europa.eu/project/id/201052/reporting

### Content
- **Samples**: 124 European individuals (healthy, overweight, obese, IBD patients) from Denmark and Spain
- **Gene Catalog**: 3.3 million unique open reading frames (ORFs)
- **Focus Diseases**: Inflammatory Bowel Disease (IBD) and obesity
- **Key Finding**: Discovered human enterotypes (3 distinct gut bacterial community clusters)

### API Access
| Method | Description |
|--------|-------------|
| **ENA Browser Tools** | Python utilities for accession-based downloads |
| **wget/FTP** | Direct file download from EBI FTP servers |
| **Globus** | Fast, reliable file transfer service |
| **MAHMI REST API** | Public API for immunomodulatory/antiproliferative peptide data |

### License
- **Type**: Open Access (specific terms vary by dataset)
- **Terms**: Data deposited in European Nucleotide Archive (ENA) follows EBI data sharing policies

### Data Size
- **Sequence Data**: 576.7 Gb of raw sequence data
- **Gene Catalog**: ~200x more data than all previous studies combined

### Gene Platform Integration Notes
- Foundation dataset for gut microbiome research
- Gene catalog useful for metagenomic annotation pipelines
- Enterotype classification relevant for personalization features

---

## 3. GMrepo (Gut Microbiome Repository)

### Overview
Curated database of human gut metagenomes designed to increase reusability and accessibility of gut metagenomic data, enabling cross-project and phenotype comparisons.

### URLs
- **Main Website**: https://gmrepo.humangut.info
- **Documentation**: https://evolgeniusteam.github.io/gmrepodocumentation/
- **GitHub (API Access)**: https://github.com/evolgeniusteam/GMrepoProgrammableAccess

### Content (v3 - Latest)
| Metric | Count |
|--------|-------|
| Projects | 890 |
| Samples/Runs | 118,965 (87,048 16S + 31,917 WGS) |
| Diseases | 302 annotated |
| Marker Taxa | 1,299 (726 species, 573 genera) |
| Phenotype Pairs | 167 |

### API Access
| Method | Description |
|--------|-------------|
| **REST API** | Full programmatic access to database contents |
| **R Examples** | Sample code for R-based data retrieval |
| **Perl Examples** | Sample code for Perl-based data retrieval |
| **Python Examples** | Sample code for Python-based data retrieval |

### API Endpoints (Examples)
```
GET /api/getRunsByProjectID
GET /api/getRunsByPhenotype
GET /api/getSpeciesAbundance
GET /api/getMarkerTaxa
```

### License
- **Type**: Creative Commons Attribution Non-Commercial 3.0 (CC BY-NC 3.0)
- **Terms**: Free for all academic users; commercial use requires permission

### Data Size
- **Note**: Raw sequence data not directly downloadable due to hardware limitations
- **Processed Data**: Relative abundances, taxonomic profiles, metadata available for download

### Gene Platform Integration Notes
- Excellent for disease-marker microbiome associations
- Pre-computed relative abundances reduce processing needs
- Cross-dataset comparison features valuable for meta-analyses
- REST API enables real-time integration

---

## 4. mBodyMap

### Overview
Curated database for microbes across the human body and their associations with health and diseases, focusing on reusability of human-associated metagenomic data.

### URLs
- **Main Website**: https://mbodymap.microbiome.cloud
- **GitHub**: https://github.com/whchenlab/mBodymap
- **API Documentation**: https://github.com/whchenlab/mBodymap/tree/main/programmable-access

### Content
| Metric | Count |
|--------|-------|
| Total Runs | 63,148 |
| Valid Runs | 61,913 |
| Metagenomes | 14,401 |
| Amplicons | 48,747 |
| Body Sites | 22 |
| Diseases | 56 |
| Projects | 136 |
| Species Identified | 6,247 (from 1,645 genera) |

### API Access
| Method | Description |
|--------|-------------|
| **REST API** | Programmatic access via documented endpoints |
| **R Access** | Sample code available on GitHub |
| **Python Access** | Sample code available on GitHub |

### API Endpoints (Examples)
```
GET /api/getAllBodySites
GET /api/getSpeciesOfOneBodySite
GET /api/getRelativeAbundanceOfOneBodySite
GET /api/getDiseaseAssociations
```

### License
- **Type**: Creative Commons Attribution Non-Commercial 3.0 (CC BY-NC 3.0)
- **Terms**: Free for academic users

### Data Size
- Processed abundance data and metadata available for download
- Bulk download available from Help page

### Gene Platform Integration Notes
- Multi-body-site coverage unique among databases
- Disease-microbe associations pre-curated
- Useful for systemic microbiome analysis beyond gut
- MeSH standardization enables integration with medical ontologies

---

## 5. gutMGene (Gut Microbiota-Gene Database)

### Overview
Manually curated database providing comprehensive resource of target genes of gut microbes and microbial metabolites in humans and mice.

### URLs
- **Main Website (v2.0)**: http://bio-computing.hrbmu.edu.cn/gutmgene
- **Legacy (v1.0)**: http://bio-annotation.cn/gutmgene
- **Database Commons Entry**: https://ngdc.cncb.ac.cn/databasecommons/database/id/8012

### Content
**Version 1.0:**
| Organism | Relationships | Microbes | Metabolites | Genes |
|----------|---------------|----------|-------------|-------|
| Human | 1,331 | 332 | 207 | 223 |
| Mouse | 2,349 | 209 | 149 | 544 |

**Version 2.0:**
| Source | Associations |
|--------|-------------|
| Literature-based | 4,860 |
| Metabolic Reconstructions | 2,476,040 |

### API Access
| Method | Description |
|--------|-------------|
| **Web Interface** | Search by microbes, disorders, intervention measures |
| **Autocomplete** | Query refinement for microbes, metabolites, genes |
| **Bulk Download** | Resource page for complete dataset download |

### License
- **Type**: Creative Commons Attribution Non-Commercial (CC BY-NC 4.0)
- **Terms**: Non-commercial use permitted with citation; commercial use requires permission from journals.permissions@oup.com

### Data Size
- ~2.5 million associations (v2.0 with metabolic reconstructions)
- Compact structured database (relationships, not raw sequences)

### Gene Platform Integration Notes
- Critical for microbiome-gene interaction mapping
- Metabolite data enables pathway analysis
- Human + mouse data supports translational research
- Causal vs correlational classification aids interpretation

---

## 6. Microbiome-Drug Interaction Databases

### 6.1 MASI (Microbiota-Active Substance Interactions)

#### URLs
- **Main Website**: https://www.aiddlab.com/MASI/
- **About Page**: https://www.aiddlab.com/MASI/about.html

#### Content
| Category | Count |
|----------|-------|
| Pharmaceutical Substances | 1,051 |
| Dietary Substances | 103 |
| Herbal Substances | 119 |
| Probiotic Substances | 46 |
| Environmental Substances | 142 |
| Microbiota Species | 806 |
| Linked Diseases | 56 |
| Microbiota-Disease Associations | 784 |

**Interactions:**
- Bacteria-Pharmaceutical: 11,215
- Bacteria-Herbal: 914
- Bacteria-Dietary: 309
- Bacteria-Environmental: 753

#### API Access
| Method | Description |
|--------|-------------|
| **Web Browse** | Browse all entries by category |
| **Download** | CSV, Excel, PDF export for table data |
| **Bulk Download** | Complete datasets in plain text or Excel |

#### License
- Not explicitly stated; academic use appears unrestricted

#### Data Size
- Structured database with ~13,000+ interaction entries
- Cross-linked to NCBI Taxonomy, ChEMBL, DrugBank, TTD, PubChem, NPASS

---

### 6.2 PharmacoMicrobiomics Portal

#### URLs
- **Main Website**: http://pharmacomicrobiomics.com/
- **Alternative**: http://www.pharmacomicrobiomics.org
- **Source Code**: http://sourceforge.net/projects/pharmacomicro

#### Content
- Drug-Microbiome Interactions: 60+ drugs
- Curated from: 100+ research and review articles
- Focus: Drug metabolism, toxicity, efficacy modification by microbiome

#### API Access
| Method | Description |
|--------|-------------|
| **Web Portal** | Search engine for drug-microbiome interactions |
| **Cross-Links** | Integration with PubMed, PubChem, Comparative Toxicogenomics |

#### License
- **Type**: Creative Commons Attribution 2.5

#### Technical Details
- Built with Django (Python), JQuery
- Currently hosted by Digital Ocean

---

### 6.3 MDAD (Microbe-Drug Association Database)

#### URLs
- **Main Publication**: https://www.frontiersin.org/articles/10.3389/fcimb.2018.00424/full

#### Content
| Metric | Count |
|--------|-------|
| Total Entries | 5,055 |
| Drugs | 1,388 (1,373 after deduplication) |
| Microbes | 180 (173 after deduplication) |
| Unique Associations | 2,470 |

#### Data Sources
- DrugBank, TTD, SuperTarget, MATADOR
- TDR targets, PDTD, ChEMBL, aBiofilm, Integrity

#### API Access
| Method | Description |
|--------|-------------|
| **Web Interface** | Browse and search functionality |
| **Data Submission** | User upload with manual verification |

#### License
- **Type**: Creative Commons Attribution (CC BY)
- **Terms**: Open access with citation requirement

---

## 7. Probiotic Databases

### 7.1 PROBIO Database

#### URLs
- **Main Website**: http://bidd.group/probio/homepage.htm
- **Reference**: Tao et al. 2017, Journal of Agricultural and Food Chemistry

#### Content
| Category | Strains |
|----------|---------|
| Marketed Probiotics | 448 |
| Clinical Trial/Field Trial | 167 |
| Research Probiotics | 382 |

- **Applications**: Humans, animals, plants

#### API Access
| Method | Description |
|--------|-------------|
| **Web Interface** | Browse by category and function |
| **Download** | Statistics and data downloads available |

#### License
- Academic use permitted (specific terms not explicitly stated)

---

### 7.2 Probio-ichnos

#### URLs
- **Main Publication**: https://www.mdpi.com/2076-2607/12/10/1955
- **PubMed**: https://pubmed.ncbi.nlm.nih.gov/39458265/

#### Content
| Metric | Count |
|--------|-------|
| Entries | 12,993 |
| Distinct Strains | 11,202 |
| Species | 470 |
| Genera | 143 |
| Source Studies | 2,207 |

**Probiotic Characteristics Covered:**
- Resistance to acid and bile
- Attachment to host epithelia
- Antimicrobial activity
- Immunomodulatory activity
- Antiproliferative activity
- Antioxidant activity

#### API Access
| Method | Description |
|--------|-------------|
| **Web Interface** | Browse by species, strain, or probiotic trait |
| **Search** | Strain-specific attribute queries |

#### License
- **Type**: Creative Commons Attribution (CC BY 4.0)
- **Terms**: Full reuse permitted with citation

---

### 7.3 ProbioMinServer

#### URLs
- **Main Website**: https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad153/7321112

#### Content
- **Species Covered**: 25 common probiotic species (11 genera)
- **Analysis Features**: ARGs, virulence factors, pathogenic genes, plasmid types, prophage regions, CAZy, metabolite gene clusters

#### API Access
| Method | Description |
|--------|-------------|
| **Web Interface** | Keyword search, browsing |
| **Download** | Tab-delimited file export |
| **Genome Upload** | WGS analysis for probiotic potential |

#### License
- Academic use (specific terms via publication)

---

### 7.4 IPDB (Integrated Probiotic DataBase)

#### URLs
- **Web Access**: http://probiogenomics.unipr.it/cmu/
- **Direct Download**: http://probiogenomics.unipr.it/files/Probiotic_Bifidobacteria_DataBase

#### Content
- 34 publicly available Bifidobacterium strains
- Focus: Commercialized health-promoting supplements
- Genomic compendium of bifidobacterial strains

#### License
- Academic research use

---

## Comparative Summary Table

| Database | Focus | Samples/Entries | API | License | Best For |
|----------|-------|-----------------|-----|---------|----------|
| **HMP** | Multi-omic human microbiome | 11,000+ samples, >48TB | REST, GraphQL, GA4GH | CC BY 4.0 | Reference baseline |
| **MetaHIT** | Gut metagenomics | 124 subjects, 3.3M genes | ENA tools, FTP | Open Access | Gene catalog |
| **GMrepo** | Gut microbiome curation | 118,965 samples | REST API | CC BY-NC 3.0 | Disease markers |
| **mBodyMap** | Multi-body-site microbiome | 63,148 runs | REST API | CC BY-NC 3.0 | Body-wide analysis |
| **gutMGene** | Microbe-gene interactions | 2.5M associations | Web/Download | CC BY-NC 4.0 | Pathway analysis |
| **MASI** | Drug-microbiome | 13,000+ interactions | Download | Not stated | Pharmacomicrobiomics |
| **MDAD** | Drug-microbe | 5,055 entries | Web | CC BY | Drug associations |
| **Probio-ichnos** | Probiotic strains | 12,993 entries | Web | CC BY 4.0 | Strain selection |

---

## Integration Recommendations for Gene Platform

### Priority 1: Core Databases
1. **GMrepo v3** - REST API ready, disease markers, academic-friendly license
2. **gutMGene v2.0** - Essential for gene-microbiome linking
3. **MASI** - Comprehensive drug interaction coverage

### Priority 2: Reference Data
4. **HMP** - Gold standard healthy reference via AWS
5. **mBodyMap** - Multi-site coverage with API

### Priority 3: Specialized
6. **Probio-ichnos** - For probiotic recommendation features
7. **MetaHIT** - Gene catalog for annotation
8. **MDAD/PharmacoMicrobiomics** - Supplementary drug data

### Technical Considerations
- **Commercial Use**: HMP (CC BY) and Probio-ichnos (CC BY) have fewest restrictions
- **Academic Only**: GMrepo, mBodyMap, gutMGene require NC license compliance
- **API Availability**: GMrepo, mBodyMap have best-documented REST APIs
- **Data Volume**: HMP requires significant storage (~50TB); others are more compact

---

## References

1. Human Microbiome Project Consortium. (2012). Nature, 486(7402), 207-214.
2. Qin, J., et al. (2010). A human gut microbial gene catalogue. Nature, 464(7285), 59-65.
3. Wu, S., et al. (2021). GMrepo v2. Nucleic Acids Research, 50(D1), D777-D784.
4. Jin, H., et al. (2022). mBodyMap. Nucleic Acids Research, 50(D1), D808-D816.
5. Cheng, L., et al. (2022). gutMGene. Nucleic Acids Research, 50(D1), D795-D800.
6. Dai, Z., et al. (2021). MASI. Nucleic Acids Research, 49(D1), D776-D782.
7. Sun, Y., et al. (2018). MDAD. Frontiers in Cellular and Infection Microbiology, 8, 424.
8. Tsifintaris, M., et al. (2024). Probio-ichnos. Microorganisms, 12(10), 1955.
9. Tao, W., et al. (2017). PROBIO. Journal of Agricultural and Food Chemistry.

---

*Document generated: January 2026*
*For Gene Platform data integration planning*
