# Public Data Sources Inventory

**Last Updated:** January 2026
**Constraint:** All sources must be publicly available - no commercial licenses, no partnerships required.

---

## Summary

| Category | Sources | Total Records | Est. Size |
|----------|---------|---------------|-----------|
| SNP/Variants | 6 | 1.2B+ SNPs | ~20 GB |
| Genes/Proteins | 4 | 250M+ proteins | ~120 GB |
| Pathways | 4 | 5K+ pathways | ~5 GB |
| Pharmacogenomics | 4 | 60K+ interactions | ~5 GB |
| Supplements | 3 | 140K+ products | ~1 GB |
| Traditional Medicine | 6 | 100K+ compounds | ~10 GB |
| Research | 2 | 36M+ citations | ~100 GB |
| Specialized | 4 | Various | ~5 GB |

---

## 1. SNP & Genetic Variant Databases

### 1.1 dbSNP (NCBI) - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/snp/ |
| **Content** | 1.2 billion Reference SNP (rs) records |
| **Version** | Build 157 (March 2025) |
| **API** | E-utilities API, Variation Services API |
| **Formats** | JSON, VCF, HGVS, SPDI |
| **License** | Public domain |
| **FTP** | ftp.ncbi.nlm.nih.gov/snp/ |
| **Size** | ~15 GB (VCF compressed) |

### 1.2 ClinVar (NCBI) - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | 3M+ variants with clinical significance |
| **API** | E-utilities API |
| **Formats** | XML, VCF, TSV |
| **License** | Public domain |
| **FTP** | ftp.ncbi.nlm.nih.gov/pub/clinvar/ |
| **Size** | ~500 MB compressed |
| **Updates** | Weekly (Mondays), monthly full release |

### 1.3 GWAS Catalog (EMBL-EBI)
| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Content** | 500K+ SNP-trait associations (p≤5×10⁻⁸) |
| **API** | REST API |
| **Formats** | TSV, JSON |
| **License** | Open access |
| **Size** | ~100 MB |

### 1.4 SNPedia
| Field | Value |
|-------|-------|
| **URL** | https://www.snpedia.com/ |
| **Content** | 109,729 SNPs with human-curated phenotype annotations |
| **API** | MediaWiki API, SNPediaR (R package) |
| **Bulk Download** | https://snpedia.com/index.php/Bulk |
| **License** | CC BY-NC-SA 3.0 |
| **Size** | ~1 GB |
| **Note** | Best for human-readable interpretations |

### 1.5 gnomAD
| Field | Value |
|-------|-------|
| **URL** | https://gnomad.broadinstitute.org/ |
| **Content** | Allele frequencies from 807,162 individuals (v4.1) |
| **API** | GraphQL API |
| **Download** | AWS Open Data Registry (s3://gnomad-public-us-east-1/) |
| **License** | Open access |
| **Size** | ~30 TB (full), ~50 GB (summary) |
| **Note** | Use AWS to avoid egress fees |

### 1.6 1000 Genomes Project
| Field | Value |
|-------|-------|
| **URL** | https://www.internationalgenome.org/ |
| **Content** | Population-level genetic variation |
| **Download** | FTP, AWS (s3://1000genomes/) |
| **License** | Open access |
| **Size** | ~260 TB |

---

## 2. Gene & Protein Databases

### 2.1 GeneCards
| Field | Value |
|-------|-------|
| **URL** | https://www.genecards.org/ |
| **Content** | Comprehensive gene info from ~200 sources |
| **API** | JSON dumps (v5.24, May 2025) |
| **License** | Free for academic use |
| **Size** | ~5 GB |
| **Note** | Best for gene summaries |

### 2.2 UniProt
| Field | Value |
|-------|-------|
| **URL** | https://www.uniprot.org/ |
| **Content** | 250M+ proteins |
| **API** | REST API, SPARQL endpoint |
| **Formats** | XML, FASTA, JSON |
| **License** | CC BY 4.0 |
| **FTP** | ftp.uniprot.org |
| **Size** | ~110 GB (full), ~5 GB (human only) |

### 2.3 NCBI Gene
| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/gene/ |
| **Content** | 200K+ genes with cross-references |
| **API** | E-utilities |
| **License** | Public domain |
| **FTP** | ftp.ncbi.nlm.nih.gov/gene/ |
| **Size** | ~5 GB |

### 2.4 Ensembl
| Field | Value |
|-------|-------|
| **URL** | https://www.ensembl.org/ |
| **Content** | Genome annotations, variants, regulation |
| **API** | REST API |
| **Formats** | GFF, VCF, FASTA |
| **License** | Open access |
| **FTP** | ftp.ensembl.org |
| **Size** | ~50 GB |

---

## 3. Biochemical Pathway Databases

### 3.1 Reactome - PRIMARY (CC0)
| Field | Value |
|-------|-------|
| **URL** | https://reactome.org/ |
| **Content** | 2,500+ expert-curated human pathways |
| **API** | REST API, SPARQL |
| **Formats** | BioPAX, SBML, SBGN |
| **License** | CC0 (public domain) |
| **Download** | reactome.org/download-data |
| **Size** | ~200 MB |
| **Note** | Preferred over KEGG (no license issues) |

### 3.2 WikiPathways - PRIMARY (CC0)
| Field | Value |
|-------|-------|
| **URL** | https://wikipathways.org/ |
| **Content** | 3,000+ community-curated pathways |
| **API** | REST API, SPARQL |
| **Formats** | GPML, BioPAX |
| **License** | CC0 (public domain) |
| **Size** | ~100 MB |
| **Note** | Largest metabolic pathway collection |

### 3.3 KEGG (Limited Access)
| Field | Value |
|-------|-------|
| **URL** | https://www.genome.jp/kegg/pathway.html |
| **Content** | 500+ metabolic, signaling, disease pathways |
| **API** | REST API (web access only) |
| **License** | Academic license required for FTP |
| **Note** | Use Reactome/WikiPathways for bulk data |

### 3.4 SMPDB
| Field | Value |
|-------|-------|
| **URL** | https://smpdb.ca/ |
| **Content** | Metabolic pathways with compound structures |
| **License** | Open access |

---

## 4. Pharmacogenomics & Drug Databases

### 4.1 PharmGKB - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://www.pharmgkb.org/ |
| **Content** | 11,000+ variant-drug interactions |
| **Download** | TSV files (registration required) |
| **License** | CC BY-SA 4.0 |
| **Size** | ~500 MB |
| **Note** | Best for drug-gene relationships |

### 4.2 Open Targets
| Field | Value |
|-------|-------|
| **URL** | https://platform.opentargets.org/ |
| **Content** | Gene-disease associations |
| **API** | GraphQL API |
| **Formats** | Parquet, JSON |
| **License** | Open access |
| **FTP** | ftp.ebi.ac.uk/pub/databases/opentargets/ |
| **Size** | ~10 GB |

### 4.3 ChEMBL
| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/chembl/ |
| **Content** | 2M+ compounds, bioactivity data |
| **API** | REST API |
| **Formats** | SQL, SDF |
| **License** | CC BY-SA 3.0 |
| **FTP** | ftp.ebi.ac.uk/pub/databases/chembl/ |
| **Size** | ~3 GB |
| **Note** | Use instead of DrugBank (no license) |

### 4.4 CPIC Guidelines
| Field | Value |
|-------|-------|
| **URL** | https://cpicpgx.org/ |
| **Content** | Clinical dosing guidelines by genotype |
| **Formats** | JSON, PDF |
| **License** | Open access |

---

## 5. Supplement & Nutrient Databases

### 5.1 NIH ODS API - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://ods.od.nih.gov/api/ |
| **Content** | Official fact sheets on vitamins, minerals, supplements |
| **API** | REST API (free, no registration) |
| **Formats** | JSON |
| **License** | Public domain |
| **Size** | ~50 MB |

### 5.2 DSLD (Dietary Supplement Label Database)
| Field | Value |
|-------|-------|
| **URL** | https://dsld.od.nih.gov/ |
| **Content** | 140,000+ supplement product labels |
| **API** | https://dsld.od.nih.gov/api-guide |
| **Formats** | JSON |
| **License** | Public domain |
| **Size** | ~500 MB |

### 5.3 FDA NDC / DailyMed
| Field | Value |
|-------|-------|
| **URL** | https://open.fda.gov / https://dailymed.nlm.nih.gov |
| **Content** | Drug products, labels |
| **API** | REST API |
| **Formats** | JSON, XML |
| **License** | Public domain |
| **Size** | ~5 GB |

---

## 6. Traditional Medicine Databases

### 6.1 TCMSP (TCM Systems Pharmacology) - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://tcmsp-e.com/ |
| **Content** | 499 herbs, 29,384 ingredients, 3,311 targets |
| **License** | Open Database License |
| **Priority** | HIGH - TCM primary source |

### 6.2 TCMBank
| Field | Value |
|-------|-------|
| **Content** | 9,192 herbs, 61,966 ingredients, 15,179 targets |
| **License** | Non-commercial |
| **Note** | Largest downloadable TCM database |

### 6.3 HERB Database
| Field | Value |
|-------|-------|
| **URL** | http://herb.ac.cn/ |
| **Content** | 7,263 herbs, 49,258 ingredients |
| **License** | Open |

### 6.4 IMPPAT 2.0 (Ayurveda) - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://cb.imsc.res.in/imppat/ |
| **Content** | 4,010 Indian plants, 17,967 phytochemicals |
| **Formats** | SDF, MOL, MOL2, PDB, PDBQT |
| **License** | Academic use |

### 6.5 KampoDB (Japanese Kampo)
| Field | Value |
|-------|-------|
| **URL** | Publication: nature.com/articles/s41598-018-29516-1 |
| **Content** | Kampo formulas, crude drugs, compounds, targets |
| **Method** | Docking + ML predictions |

### 6.6 WAKANYAKU Database Portal
| Field | Value |
|-------|-------|
| **URL** | https://www.inm.u-toyama.ac.jp/en/database/ |
| **Content** | Natural compound-protein interactions |
| **License** | Academic |

---

## 7. Research Literature

### 7.1 PubMed/PMC (NCBI) - PRIMARY
| Field | Value |
|-------|-------|
| **URL** | https://pubmed.ncbi.nlm.nih.gov/ |
| **Content** | 36M+ biomedical citations |
| **API** | E-utilities, Entrez Direct |
| **Bulk Download** | Annual baseline (1,274 XML files) + daily updates |
| **License** | Public domain (abstracts) |
| **FTP** | ftp.ncbi.nlm.nih.gov/pubmed/baseline/ |
| **Size** | ~100 GB (baseline XML) |
| **Hugging Face** | ncbi/pubmed dataset |

### 7.2 PMC Open Access Subset
| Field | Value |
|-------|-------|
| **Content** | Full-text articles under open licenses |
| **FTP** | ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/ |
| **Note** | Only source for legal bulk full-text |

---

## 8. Specialized Databases

### 8.1 OMIM
| Field | Value |
|-------|-------|
| **URL** | https://omim.org/ |
| **Content** | 16K+ gene-disease relationships |
| **API** | API with registration |
| **License** | Academic |

### 8.2 Human Phenotype Ontology (HPO)
| Field | Value |
|-------|-------|
| **URL** | https://hpo.jax.org/ |
| **Content** | 18K+ standardized phenotype terms |
| **License** | Open |
| **Use** | Symptom standardization |

### 8.3 DisGeNET
| Field | Value |
|-------|-------|
| **URL** | https://www.disgenet.org/ |
| **Content** | 1M+ gene-disease associations |
| **API** | REST API |
| **License** | CC BY-NC-SA |

### 8.4 STITCH
| Field | Value |
|-------|-------|
| **URL** | http://stitch.embl.de/ |
| **Content** | Chemical-protein interactions |
| **API** | REST API |

---

## Not Available (License Required)

| Source | Issue | Alternative |
|--------|-------|-------------|
| **DrugBank** (full) | Commercial license | Use **ChEMBL** |
| **KEGG** (FTP) | Academic license | Use **Reactome + WikiPathways** |
| **NatMed Pro** | Enterprise only | Use **NIH ODS + PubMed** |
| **Examine.com** | No API | Manual curation or PubMed |
| **OpenSNP** | Shut down April 2025 | Archived on Zenodo/GenomePrep |

---

## Ingestion Priority

### Tier 1 - Core (Weeks 1-4)
1. SNPedia (109K SNPs, human-curated)
2. dbSNP subset (1-10M key SNPs)
3. ClinVar (clinical variants)
4. PharmGKB (drug-gene interactions)
5. Reactome (pathways, CC0)

### Tier 2 - Enrichment (Weeks 5-8)
1. GeneCards (gene summaries)
2. GWAS Catalog (trait associations)
3. UniProt (protein function)
4. WikiPathways (community pathways)
5. TCMSP (TCM compounds)

### Tier 3 - Expansion (Weeks 9-12)
1. IMPPAT (Ayurveda)
2. PubMed abstracts (for RAG)
3. gnomAD (allele frequencies)
4. ChEMBL (compounds)
