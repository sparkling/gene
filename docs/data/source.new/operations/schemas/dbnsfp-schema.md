---
id: schema-dbnsfp
title: "dbNSFP Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./_index.md)

# dbNSFP Schema Documentation

**Document ID:** SCHEMA-DBNSFP
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** dbNSFP Official Documentation, Genome Medicine Publication

---

## TL;DR

dbNSFP is a comprehensive database integrating functional predictions from 36+ algorithms for all possible non-synonymous SNVs (~84M) in the human genome. It serves as a one-stop resource combining deleteriousness scores (SIFT, PolyPhen-2, CADD, REVEL, AlphaMissense), conservation metrics (PhyloP, GERP++, phastCons), and population frequencies (gnomAD, 1000 Genomes) with gene-level annotations.

**Key Facts:**
- **Coverage:** 84M+ nsSNVs and splice-site variants
- **Prediction Scores:** 36+ algorithms integrated
- **Conservation Scores:** 9 evolutionary metrics
- **Population Data:** gnomAD v4.1, 1000 Genomes, TOPMed, All of Us
- **License:** Academic free; commercial requires contact
- **Current Version:** 5.3.1 (January 2026)

---

## Overview

dbNSFP (database for Non-synonymous SNPs' Functional Predictions) is developed and maintained by Dr. Xiaoming Liu's lab. It provides pre-computed functional predictions and annotations for every potential protein-altering variant in the human genome, eliminating the need to run individual prediction tools.

**Primary Use Cases:**
- Variant pathogenicity assessment
- Clinical variant interpretation
- Research prioritization
- Annotation pipeline integration

**Official Website:** https://www.dbnsfp.org/

**Legacy Site:** https://sites.google.com/site/jpopgen/dbNSFP

---

## License Information

| Branch | License | Excluded Resources |
|--------|---------|-------------------|
| **dbNSFP Academic (a)** | Free for academic use | None - includes all resources |
| **dbNSFP Commercial (c)** | Contact required | PolyPhen-2, VEST, REVEL, ClinPred, CADD, LINSIGHT, GenoCanyon, MutScore |

**Commercial Contact:** collaboration@dbnsfp.org

**Citation Requirements:**
- v1.x: Paper 1 only
- v2.x: Papers 1 & 2
- v3.x: Papers 1 & 3
- v4.x+: Papers 1 & 4

---

## Database Statistics (v5.3.1)

| Metric | Value |
|--------|-------|
| **Non-synonymous SNVs** | 83,049,507 |
| **Splice-site SNVs** | 2,446,464 |
| **Total Variants** | 85,495,971 |
| **Prediction Algorithms** | 36+ |
| **Conservation Scores** | 9 |
| **Population Databases** | 6 |
| **GENCODE Version** | Release 49 (Ensembl 115) |
| **Reference Genome** | GRCh38/hg38 (primary) |

---

## Version History

| Version | Date | Key Additions |
|---------|------|---------------|
| **5.3.1** | Jan 2026 | GENCODE 49, gnomAD v4.1, All of Us |
| **4.9a** | Aug 2024 | MutScore, ClinVar 20240805 |
| **4.8a** | Jun 2024 | PHACTboost, MutFormer, GERP_91_mammals |
| **4.7a** | Mar 2024 | CADD v1.7, gnomAD v4.0.0 |
| **4.6a** | Feb 2024 | GTEx V8 sQTLs, eQTLGen phase I |
| **4.5a** | Sep 2023 | AlphaMissense, ESM1b, EVE |
| **4.4a** | Apr 2023 | VARITY, MPC |
| **4.1a** | Nov 2020 | Initial v4 release |

---

## Primary Variant Fields

### Position and Allele Columns

| Column | Type | Description |
|--------|------|-------------|
| `chr` | string | Chromosome number (1-22, X, Y, M) |
| `pos(1-based)` | integer | Physical position on chromosome (hg38, 1-based) |
| `ref` | string | Reference nucleotide allele (+ strand) |
| `alt` | string | Alternative nucleotide allele (+ strand) |
| `aaref` | string | Reference amino acid ("." for splice variants) |
| `aaalt` | string | Alternative amino acid ("." for splice variants) |
| `hg19_chr` | string | Chromosome as to hg19 ("." if missing) |
| `hg19_pos(1-based)` | integer | Position in hg19 coordinates |
| `hg18_pos(1-based)` | integer | Position in hg18 coordinates (deprecated) |

### Gene and Transcript Identifiers

| Column | Type | Description |
|--------|------|-------------|
| `genename` | string | HGNC gene symbol |
| `Ensembl_geneid` | string | Ensembl gene ID (ENSG) |
| `Ensembl_transcriptid` | string | Ensembl transcript IDs (semicolon-separated) |
| `Ensembl_proteinid` | string | Ensembl protein IDs (semicolon-separated) |
| `Uniprot_acc` | string | UniProt accession numbers (semicolon-separated) |
| `Uniprot_id` | string | UniProt entry names (semicolon-separated) |
| `Uniprot_aapos` | string | AA positions in UniProt (semicolon-separated) |

### Amino Acid Context

| Column | Type | Description |
|--------|------|-------------|
| `aapos` | integer | Amino acid position in protein (-1 for splice sites) |
| `refcodon` | string | Reference codon sequence |
| `codonpos` | integer | Position within codon (1, 2, or 3) |
| `cds_strand` | string | Coding sequence strand (+ or -) |
| `fold_degenerate` | integer | Codon degeneracy (0, 2, or 3) |
| `Ancestral_allele` | string | Inferred ancestral allele |

---

## Prediction Score Columns

### SIFT Family

| Column | Type | Range | Interpretation |
|--------|------|-------|----------------|
| `SIFT_score` | float | 0-1 | Lower = more damaging; <0.05 = Damaging |
| `SIFT_converted_rankscore` | float | 0-1 | Percentile rank (higher = more damaging) |
| `SIFT_pred` | string | D/T | D(amaging) or T(olerated) |
| `SIFT4G_score` | float | 0-1 | SIFT4G algorithm score |
| `SIFT4G_converted_rankscore` | float | 0-1 | SIFT4G percentile rank |
| `SIFT4G_pred` | string | D/T | D(amaging) or T(olerated) |

### PolyPhen-2

| Column | Type | Range | Interpretation |
|--------|------|-------|----------------|
| `Polyphen2_HDIV_score` | string | 0-1 | HumDiv-trained score (semicolon-separated) |
| `Polyphen2_HDIV_rankscore` | float | 0-1 | Percentile rank |
| `Polyphen2_HDIV_pred` | string | D/P/B | D(amaging), P(ossibly), B(enign) |
| `Polyphen2_HVAR_score` | string | 0-1 | HumVar-trained score (semicolon-separated) |
| `Polyphen2_HVAR_rankscore` | float | 0-1 | Percentile rank |
| `Polyphen2_HVAR_pred` | string | D/P/B | D(amaging), P(ossibly), B(enign) |

**Note:** PolyPhen-2 excluded from commercial branch (dbNSFP c).

### CADD

| Column | Type | Range | Interpretation |
|--------|------|-------|----------------|
| `CADD_raw` | float | -7 to 35+ | Raw C-score |
| `CADD_raw_rankscore` | float | 0-1 | Percentile rank |
| `CADD_phred` | float | 0-99 | PHRED-scaled; >20 = top 1% |

**Version:** CADD v1.7 (as of dbNSFP 4.7+)

**Note:** CADD excluded from commercial branch.

### REVEL

| Column | Type | Range | Interpretation |
|--------|------|-------|----------------|
| `REVEL_score` | float | 0-1 | Ensemble score; higher = more pathogenic |
| `REVEL_rankscore` | float | 0-1 | Percentile rank |

**Threshold:** >0.5 commonly used for pathogenic prediction.

**Note:** REVEL excluded from commercial branch.

### AlphaMissense

| Column | Type | Range | Interpretation |
|--------|------|-------|----------------|
| `AlphaMissense_score` | float | 0-1 | Pathogenicity probability |
| `AlphaMissense_rankscore` | float | 0-1 | Percentile rank |
| `AlphaMissense_pred` | string | LB/A/LP | likely_benign/ambiguous/likely_pathogenic |

**Thresholds:** <0.34 = LB, >0.564 = LP

### Meta-Predictors

| Column | Type | Range | Interpretation |
|--------|------|-------|----------------|
| `MetaSVM_score` | float | -2 to 3 | SVM ensemble; >0 = Damaging |
| `MetaSVM_rankscore` | float | 0-1 | Percentile rank |
| `MetaSVM_pred` | string | D/T | D(amaging) or T(olerated) |
| `MetaLR_score` | float | 0-1 | Logistic regression ensemble |
| `MetaLR_rankscore` | float | 0-1 | Percentile rank |
| `MetaLR_pred` | string | D/T | D(amaging) or T(olerated) |
| `MetaRNN_score` | float | 0-1 | Recurrent neural network ensemble |
| `MetaRNN_rankscore` | float | 0-1 | Percentile rank |
| `MetaRNN_pred` | string | D/T | D(amaging) or T(olerated) |

### Additional Predictors

| Column | Type | Description |
|--------|------|-------------|
| `LRT_score` | float | Likelihood ratio test p-value |
| `LRT_converted_rankscore` | float | Converted rank score |
| `LRT_pred` | string | D(eleterious), N(eutral), U(nknown) |
| `LRT_Omega` | float | Estimated dN/dS ratio |
| `MutationTaster_score` | float | MutationTaster probability |
| `MutationTaster_converted_rankscore` | float | Converted rank score |
| `MutationTaster_pred` | string | A(uto)/D(isease)/N(eutral)/P(olymorphism) |
| `MutationAssessor_score` | float | Functional impact score |
| `MutationAssessor_rankscore` | float | Percentile rank |
| `MutationAssessor_pred` | string | H(igh)/M(edium)/L(ow)/N(eutral) |
| `FATHMM_score` | float | FATHMM default score; <-1.5 = Damaging |
| `FATHMM_converted_rankscore` | float | Converted rank score |
| `FATHMM_pred` | string | D(AMAGING) or T(OLERATED) |
| `PROVEAN_score` | float | PROVEAN score; <-2.5 = Damaging |
| `PROVEAN_converted_rankscore` | float | Converted rank score |
| `PROVEAN_pred` | string | D(amaging) or N(eutral) |
| `VEST4_score` | float | VEST4 probability (0-1) |
| `VEST4_rankscore` | float | Percentile rank |
| `M-CAP_score` | float | M-CAP score; >0.025 = possibly pathogenic |
| `M-CAP_rankscore` | float | Percentile rank |
| `M-CAP_pred` | string | D(amaging) or T(olerated) |
| `MutPred_score` | float | MutPred probability |
| `MutPred_rankscore` | float | Percentile rank |
| `MVP_score` | float | MVP missense score |
| `MVP_rankscore` | float | Percentile rank |
| `PrimateAI_score` | float | PrimateAI score (0-1) |
| `PrimateAI_rankscore` | float | Percentile rank |
| `PrimateAI_pred` | string | D(amaging) or T(olerated) |
| `DEOGEN2_score` | float | DEOGEN2 probability |
| `DEOGEN2_rankscore` | float | Percentile rank |
| `DEOGEN2_pred` | string | D(amaging) or T(olerated) |
| `BayesDel_addAF_score` | float | BayesDel with AF features |
| `BayesDel_addAF_rankscore` | float | Percentile rank |
| `BayesDel_addAF_pred` | string | D(amaging) or T(olerated) |
| `BayesDel_noAF_score` | float | BayesDel without AF features |
| `BayesDel_noAF_rankscore` | float | Percentile rank |
| `BayesDel_noAF_pred` | string | D(amaging) or T(olerated) |
| `ClinPred_score` | float | ClinPred probability |
| `ClinPred_rankscore` | float | Percentile rank |
| `ClinPred_pred` | string | D(amaging) or T(olerated) |
| `LIST-S2_score` | float | LIST-S2 probability |
| `LIST-S2_rankscore` | float | Percentile rank |
| `LIST-S2_pred` | string | D(amaging) or T(olerated) |
| `VARITY_R_score` | float | VARITY score (regulatory) |
| `VARITY_ER_score` | float | VARITY score (enhanced) |
| `VARITY_R_LOO_score` | float | VARITY leave-one-out score |
| `VARITY_ER_LOO_score` | float | VARITY enhanced LOO score |
| `ESM1b_score` | float | ESM-1b language model score |
| `ESM1b_rankscore` | float | Percentile rank |
| `ESM1b_pred` | string | D(amaging) or T(olerated) |
| `EVE_score` | float | EVE variational autoencoder score |
| `EVE_rankscore` | float | Percentile rank |
| `EVE_Class10_pred` - `EVE_Class90_pred` | string | EVE predictions at various thresholds |
| `MutFormer_score` | float | MutFormer transformer score |
| `MutFormer_rankscore` | float | Percentile rank |
| `PHACTboost_score` | float | PHACTboost gradient boosting score |
| `PHACTboost_rankscore` | float | Percentile rank |
| `MutScore_score` | float | MutScore ensemble score |
| `MutScore_rankscore` | float | Percentile rank |

---

## Conservation Score Columns

| Column | Type | Description |
|--------|------|-------------|
| `phyloP100way_vertebrate` | float | PhyloP score (100 vertebrates) |
| `phyloP100way_vertebrate_rankscore` | float | Percentile rank |
| `phyloP30way_mammalian` | float | PhyloP score (30 mammals) |
| `phyloP30way_mammalian_rankscore` | float | Percentile rank |
| `phyloP17way_primate` | float | PhyloP score (17 primates) |
| `phyloP17way_primate_rankscore` | float | Percentile rank |
| `phastCons100way_vertebrate` | float | phastCons score (100 vertebrates) |
| `phastCons100way_vertebrate_rankscore` | float | Percentile rank |
| `phastCons30way_mammalian` | float | phastCons score (30 mammals) |
| `phastCons30way_mammalian_rankscore` | float | Percentile rank |
| `phastCons17way_primate` | float | phastCons score (17 primates) |
| `phastCons17way_primate_rankscore` | float | Percentile rank |
| `GERP++_NR` | float | GERP++ neutral rate |
| `GERP++_RS` | float | GERP++ rejected substitution score |
| `GERP++_RS_rankscore` | float | Percentile rank |
| `GERP_91_mammals` | float | GERP score (91 mammals) |
| `GERP_91_mammals_rankscore` | float | Percentile rank |
| `SiPhy_29way_logOdds` | float | SiPhy log-odds score |
| `SiPhy_29way_logOdds_rankscore` | float | Percentile rank |
| `SiPhy_29way_pi` | string | A:C:G:T stationary distribution |
| `bStatistic` | float | Background selection statistic |
| `bStatistic_rankscore` | float | Percentile rank |

**Interpretation:** Higher conservation scores indicate more conserved positions, suggesting functional importance.

---

## Population Frequency Columns

### gnomAD v4.1

| Column | Type | Description |
|--------|------|-------------|
| `gnomAD_exomes_AC` | integer | Exome alternative allele count |
| `gnomAD_exomes_AN` | integer | Exome total allele number |
| `gnomAD_exomes_AF` | float | Exome allele frequency |
| `gnomAD_exomes_nhomalt` | integer | Exome homozygous alt count |
| `gnomAD_exomes_AFR_AC` | integer | African/African American AC |
| `gnomAD_exomes_AFR_AF` | float | African/African American AF |
| `gnomAD_exomes_AMR_AC` | integer | Latino/Admixed American AC |
| `gnomAD_exomes_AMR_AF` | float | Latino/Admixed American AF |
| `gnomAD_exomes_ASJ_AC` | integer | Ashkenazi Jewish AC |
| `gnomAD_exomes_ASJ_AF` | float | Ashkenazi Jewish AF |
| `gnomAD_exomes_EAS_AC` | integer | East Asian AC |
| `gnomAD_exomes_EAS_AF` | float | East Asian AF |
| `gnomAD_exomes_FIN_AC` | integer | Finnish AC |
| `gnomAD_exomes_FIN_AF` | float | Finnish AF |
| `gnomAD_exomes_NFE_AC` | integer | Non-Finnish European AC |
| `gnomAD_exomes_NFE_AF` | float | Non-Finnish European AF |
| `gnomAD_exomes_SAS_AC` | integer | South Asian AC |
| `gnomAD_exomes_SAS_AF` | float | South Asian AF |
| `gnomAD_genomes_*` | various | Genome sequencing equivalents |

### Other Population Databases

| Column | Type | Description |
|--------|------|-------------|
| `1000Gp3_AC` | integer | 1000 Genomes Phase 3 alt count |
| `1000Gp3_AF` | float | 1000 Genomes Phase 3 frequency |
| `1000Gp3_*_AF` | float | Population-specific frequencies |
| `TOPMed_AC` | integer | TOPMed alt allele count |
| `TOPMed_AN` | integer | TOPMed total allele number |
| `TOPMed_AF` | float | TOPMed allele frequency |
| `ALFA_Total_AC` | integer | ALFA database alt count |
| `ALFA_Total_AF` | float | ALFA database frequency |
| `UK10K_AC` | integer | UK10K alt count |
| `UK10K_AF` | float | UK10K frequency |
| `ESP6500_AA_AC` | integer | ESP African American count |
| `ESP6500_AA_AF` | float | ESP African American frequency |
| `ESP6500_EA_AC` | integer | ESP European American count |
| `ESP6500_EA_AF` | float | ESP European American frequency |
| `ExAC_*` | various | ExAC legacy frequencies |

---

## Clinical and Functional Annotations

| Column | Type | Description |
|--------|------|-------------|
| `clinvar_id` | string | ClinVar variation ID |
| `clinvar_clnsig` | string | ClinVar clinical significance |
| `clinvar_trait` | string | ClinVar associated trait |
| `clinvar_review` | string | ClinVar review status |
| `clinvar_hgvs` | string | ClinVar HGVS nomenclature |
| `Interpro_domain` | string | InterPro domain annotation |
| `GTEx_V8_gene` | string | GTEx gene expression data |
| `GTEx_V8_tissue` | string | GTEx tissue expression |
| `Geuvadis_eQTL_target_gene` | string | eQTL target gene |
| `eQTLGen_target_gene` | string | eQTLGen eQTL targets |

---

## Gene-Level Annotations (dbNSFP_gene)

Separate file containing gene-level annotations:

| Column | Type | Description |
|--------|------|-------------|
| `Gene_name` | string | HGNC gene symbol |
| `Ensembl_gene` | string | Ensembl gene ID |
| `chr` | string | Chromosome |
| `Gene_old_names` | string | Previous gene symbols |
| `Gene_other_names` | string | Aliases |
| `Uniprot_acc` | string | UniProt accession |
| `Uniprot_id` | string | UniProt entry name |
| `Entrez_gene_id` | integer | NCBI Gene ID |
| `CCDS_id` | string | CCDS identifier |
| `Refseq_id` | string | RefSeq transcript ID |
| `ucsc_id` | string | UCSC gene ID |
| `MIM_id` | integer | OMIM gene ID |
| `MIM_phenotype_id` | string | OMIM phenotype IDs |
| `MIM_disease` | string | OMIM disease descriptions |
| `Orphanet_id` | string | Orphanet disorder IDs |
| `Orphanet_disorder` | string | Orphanet disorder names |
| `Orphanet_association_type` | string | Gene-disease association type |
| `GDI` | float | Gene damage index |
| `GDI_phred` | float | GDI PHRED-scaled |
| `LOEUF` | float | Loss-of-function observed/expected upper bound |
| `LOEUF_percentile` | float | LOEUF percentile |
| `MOEUF` | float | Missense observed/expected upper bound |
| `MOEUF_percentile` | float | MOEUF percentile |
| `pLI` | float | Probability of loss-of-function intolerance |
| `pRec` | float | Probability of recessive disease |
| `pNull` | float | Probability of null phenotype |
| `Essential_gene` | string | Essential gene designation |
| `Essential_gene_CRISPR` | string | CRISPR essentiality |
| `Essential_gene_CRISPR2` | string | CRISPR2 essentiality |
| `Essential_gene_gene-trap` | string | Gene-trap essentiality |
| `MGI_mouse_gene` | string | Mouse ortholog gene |
| `MGI_mouse_phenotype` | string | Mouse phenotype |
| `ZFIN_zebrafish_gene` | string | Zebrafish ortholog |
| `ZFIN_zebrafish_phenotype` | string | Zebrafish phenotype |
| `Expression(egenetics)` | string | Expression patterns |
| `Expression(GNF/Atlas)` | string | GNF Atlas expression |
| `GO_biological_process` | string | GO BP terms |
| `GO_cellular_component` | string | GO CC terms |
| `GO_molecular_function` | string | GO MF terms |
| `Pathway(BioCarta)_short` | string | BioCarta pathways |
| `Pathway(ConsensusPathDB)` | string | ConsensusPathDB pathways |
| `Pathway(KEGG)_full` | string | KEGG pathways |

---

## Sample Data

### Variant Record
```tsv
chr	pos(1-based)	ref	alt	aaref	aaalt	genename	Ensembl_geneid	SIFT_score	SIFT_pred	Polyphen2_HDIV_score	Polyphen2_HDIV_pred	CADD_phred	REVEL_score	AlphaMissense_score	AlphaMissense_pred
7	140753336	A	T	V	E	BRAF	ENSG00000157764	0	D	1	D	34	0.952	0.9972	LP
17	7673803	G	A	R	H	TP53	ENSG00000141510	0	D	1	D	29.5	0.918	0.9856	LP
11	5227002	T	A	E	V	HBB	ENSG00000244734	0	D	0.999	D	26.1	0.734	0.8234	LP
```

### Gene-Level Record
```tsv
Gene_name	Ensembl_gene	chr	pLI	LOEUF	MIM_disease
BRCA1	ENSG00000012048	17	0.00	0.55	Breast-ovarian cancer, familial 1
TP53	ENSG00000141510	17	0.98	0.24	Li-Fraumeni syndrome
CFTR	ENSG00000001626	7	0.00	1.02	Cystic fibrosis
```

---

## Download Instructions

### Official Sources

**Amazon S3:**
```bash
# Download academic version (all resources)
wget https://dbnsfp.s3.amazonaws.com/dbNSFP5.3.1a.zip

# Download commercial version (restricted resources removed)
wget https://dbnsfp.s3.amazonaws.com/dbNSFP5.3.1c.zip
```

**Box/Google Drive:** Links available at https://sites.google.com/site/jpopgen/dbNSFP

### File Sizes

| File | Compressed | Uncompressed |
|------|------------|--------------|
| dbNSFP5.3.1a.zip | ~40 GB | ~150 GB |
| dbNSFP_gene.gz | ~15 MB | ~50 MB |

### Verification

```bash
# Verify MD5 checksum
md5sum -c dbNSFP5.3.1a.zip.md5
```

---

## Integration Examples

### Python: Parse with Pandas

```python
import pandas as pd

def load_dbnsfp_region(filepath: str, chrom: str, start: int, end: int):
    """Load dbNSFP data for a genomic region."""
    # Use chunked reading for large files
    chunks = pd.read_csv(
        filepath,
        sep='\t',
        compression='gzip',
        chunksize=100000,
        usecols=['chr', 'pos(1-based)', 'ref', 'alt',
                 'SIFT_score', 'Polyphen2_HDIV_score',
                 'CADD_phred', 'REVEL_score', 'AlphaMissense_score']
    )

    results = []
    for chunk in chunks:
        mask = (
            (chunk['chr'] == chrom) &
            (chunk['pos(1-based)'] >= start) &
            (chunk['pos(1-based)'] <= end)
        )
        results.append(chunk[mask])

    return pd.concat(results, ignore_index=True)
```

### SnpSift Annotation

```bash
# Annotate VCF with dbNSFP
java -jar SnpSift.jar dbnsfp \
    -v \
    -db dbNSFP5.3.1a.txt.gz \
    -f SIFT_score,Polyphen2_HDIV_score,CADD_phred,REVEL_score,AlphaMissense_score \
    input.vcf > annotated.vcf
```

### VEP Plugin

```bash
# Run VEP with dbNSFP plugin
vep -i input.vcf \
    --plugin dbNSFP,dbNSFP5.3.1a.txt.gz,SIFT_score,Polyphen2_HDIV_score,CADD_phred \
    -o output.vcf
```

### ANNOVAR Integration

```bash
# Convert dbNSFP for ANNOVAR
perl prepare_annovar_user.pl \
    -dbtype dbnsfp \
    dbNSFP5.3.1a.txt.gz

# Annotate with ANNOVAR
perl table_annovar.pl \
    input.avinput \
    humandb/ \
    -buildver hg38 \
    -out output \
    -protocol dbnsfp \
    -operation f
```

---

## Score Interpretation Guide

### Ensemble Approach

For clinical interpretation, consider multiple scores:

| Scenario | Recommended Scores |
|----------|-------------------|
| **Pathogenicity screening** | REVEL, AlphaMissense, CADD |
| **Conservation assessment** | GERP++, PhyloP, phastCons |
| **Clinical validation** | ClinPred, BayesDel, MetaRNN |
| **Research prioritization** | MetaSVM, M-CAP, VEST4 |

### Concordance Analysis

High confidence when multiple predictors agree:

```python
def assess_pathogenicity(row):
    """Assess pathogenicity based on multiple scores."""
    damaging_count = 0

    if row.get('SIFT_pred') == 'D':
        damaging_count += 1
    if row.get('Polyphen2_HDIV_pred') == 'D':
        damaging_count += 1
    if row.get('CADD_phred', 0) >= 20:
        damaging_count += 1
    if row.get('REVEL_score', 0) >= 0.5:
        damaging_count += 1
    if row.get('AlphaMissense_pred') == 'LP':
        damaging_count += 1

    if damaging_count >= 4:
        return 'Likely Pathogenic'
    elif damaging_count >= 2:
        return 'Uncertain'
    else:
        return 'Likely Benign'
```

---

## Known Limitations

1. **Missense Focus:** Primarily covers non-synonymous SNVs; limited frameshift/indel support

2. **Transcript Variability:** Scores may differ by transcript; multiple values semicolon-separated

3. **Version Dependency:** Field names and meanings change across versions

4. **Missing Values:** Not all predictors available for all variants (use "." for missing)

5. **Commercial Restrictions:** Several popular scores excluded from commercial branch

6. **File Size:** Full database is very large (~40GB compressed, ~150GB uncompressed)

7. **Update Lag:** Prediction scores may be from older algorithm versions

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `chr` | Chromosome number (1-22, X, Y, M) | 17 |
| `pos(1-based)` | Physical position on chromosome in 1-based coordinates | 7673803 |
| `aaref` | Reference amino acid at variant position | R |
| `aaalt` | Alternative amino acid after variant | H |
| `SIFT_score` | Tolerance to amino acid substitution (0-1, lower=damaging) | 0.01 |
| `CADD_phred` | PHRED-scaled CADD score (>20 = top 1%) | 29.5 |
| `REVEL_score` | Ensemble pathogenicity score (0-1, higher=pathogenic) | 0.918 |
| `AlphaMissense_score` | Deep learning pathogenicity probability (0-1) | 0.9856 |
| `pLI` | Probability of loss-of-function intolerance (0-1) | 0.98 |
| `LOEUF` | Loss-of-function observed/expected upper bound | 0.24 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| nsSNV | Non-synonymous Single Nucleotide Variant changing amino acid | Core scope |
| Splice-site Variant | Variant affecting mRNA splicing | Core scope |
| Rankscore | Percentile ranking (0-1) for cross-algorithm comparison | Score normalization |
| Deleteriousness | Likelihood that variant impairs protein function | Prediction target |
| Conservation Score | Evolutionary constraint indicating functional importance | PhyloP, GERP++ |
| Meta-predictor | Ensemble combining multiple prediction algorithms | MetaSVM, MetaRNN |
| Gene Damage Index | Gene-level tolerance to damaging variants | GDI field |
| Constraint Metric | Measure of selective pressure against variants | pLI, LOEUF |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| dbNSFP | Database for Non-synonymous SNPs' Functional Predictions | Database name |
| SIFT | Sorting Intolerant From Tolerant | Prediction algorithm |
| CADD | Combined Annotation Dependent Depletion | Integrative score |
| REVEL | Rare Exome Variant Ensemble Learner | Ensemble predictor |
| GERP | Genomic Evolutionary Rate Profiling | Conservation score |
| PhyloP | Phylogenetic P-values | Conservation score |
| phastCons | Phylogenetic Analysis with Space/Time Models | Conservation score |
| gnomAD | Genome Aggregation Database | Population frequencies |
| TOPMed | Trans-Omics for Precision Medicine | Population source |
| ESP | Exome Sequencing Project | Population source |
| OMIM | Online Mendelian Inheritance in Man | Disease association |
| HGNC | HUGO Gene Nomenclature Committee | Gene naming |
| nsSNV | Non-synonymous SNV | Variant type |
| LoF | Loss of Function | Variant impact |
| VEP | Variant Effect Predictor | Annotation tool |

---

## References

1. Liu X, Li C, Mou C, et al. (2020). dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. Genome Med. 12:103.

2. Liu X, Wu C, Li C, Boerwinkle E. (2016). dbNSFP v3.0: A One-Stop Database of Functional Predictions and Annotations for Human Nonsynonymous and Splice-Site SNVs. Hum Mutat. 37:235-241.

3. Liu X, Jian X, Boerwinkle E. (2013). dbNSFP v2.0: a database of human non-synonymous SNVs and their functional predictions and annotations. Hum Mutat. 34:E2393-2402.

4. Liu X, Jian X, Boerwinkle E. (2011). dbNSFP: a lightweight database of human nonsynonymous SNPs and their functional predictions. Hum Mutat. 32:894-899.

---

## Related Resources

- **Official Website:** https://www.dbnsfp.org/
- **Legacy Site:** https://sites.google.com/site/jpopgen/dbNSFP
- **SnpSift Integration:** https://pcingola.github.io/SnpEff/snpsift/dbnsfp/
- **VEP Plugin:** https://github.com/Ensembl/VEP_plugins
- **ANNOVAR:** http://annovar.openbioinformatics.org/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
