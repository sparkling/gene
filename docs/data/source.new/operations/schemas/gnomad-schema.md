---
id: schema-gnomad
title: gnomAD Schema Documentation
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, genomics, variants]
---

**Parent:** [Schema Documentation](./_index.md)

# gnomAD Schema Documentation

**Document ID:** SCHEMA-GNOMAD
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** gnomAD v4.0 Documentation, GraphQL API, VCF Specifications

---

## TL;DR

gnomAD (Genome Aggregation Database) is the world's largest open-access collection of human genetic variation data, aggregating exome and genome sequencing from 807,162 individuals across 8+ genetic ancestry groups. The database provides allele frequencies, variant annotations, loss-of-function predictions, and constraint metrics for human disease research.

**Key Access Points:**
- GraphQL API: `https://gnomad.broadinstitute.org/api`
- Downloads: AWS, Google Cloud, Microsoft Azure
- Browser: `https://gnomad.broadinstitute.org`

---

## License & Access

| Attribute | Value |
|-----------|-------|
| **License** | Open Access (ODC-BY 1.0) |
| **Attribution** | Required - cite gnomAD consortium |
| **Commercial Use** | Permitted |
| **Data Sharing** | Permitted with attribution |
| **Terms** | https://gnomad.broadinstitute.org/terms |

**Citation:**
> Karczewski KJ, et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. Nature. 581(7809):434-443.

---

## Database Statistics (v4.0 - November 2023)

| Metric | Value | Source |
|--------|-------|--------|
| **Total Individuals** | 807,162 | v4.0 release |
| **Exome Samples** | 730,947 | Including UK Biobank |
| **Genome Samples** | 76,215 | GRCh38 aligned |
| **UK Biobank Subset** | 416,555 | Exomes |
| **Non-European Ancestry** | ~138,000 | Individuals |
| **Passing SNVs** | 62,901,592 | After QC filters |
| **Passing Indels** | 6,189,261 | After QC filters |
| **Reference Assembly** | GRCh38 | v4.0+ |
| **Storage Size (VDS)** | 18 TB | Variant Dataset format |

---

## Genetic Ancestry Groups

gnomAD stratifies allele frequencies by genetic ancestry:

| Code | Population | Description |
|------|------------|-------------|
| `afr` | African/African American | Sub-Saharan African ancestry |
| `amr` | Admixed American | Latino/Admixed American ancestry |
| `asj` | Ashkenazi Jewish | Ashkenazi Jewish ancestry |
| `eas` | East Asian | East Asian ancestry |
| `fin` | European (Finnish) | Finnish European ancestry |
| `mid` | Middle Eastern | Middle Eastern ancestry (v4+) |
| `nfe` | European (Non-Finnish) | Non-Finnish European ancestry |
| `sas` | South Asian | South Asian ancestry |
| `oth` | Other / Remaining | Remaining individuals |

**Note:** These categories are analytical groupings used for allele frequency calculations, not natural population boundaries.

---

## GraphQL API

### Endpoint

```
POST https://gnomad.broadinstitute.org/api
Content-Type: application/json
```

### Rate Limits

- Recommended: 1 request/second
- Batch queries supported for efficiency
- Large queries may be blocked

### Core Query Types

#### Variant Query

```graphql
query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    chrom
    pos
    ref
    alt
    exome {
      ac
      an
      af
      ac_hom
      populations {
        id
        ac
        an
        ac_hom
      }
    }
    genome {
      ac
      an
      af
      ac_hom
      populations {
        id
        ac
        an
        ac_hom
      }
    }
    joint {
      ac
      an
      af
    }
    flags
    lof
    lof_filter
    lof_flags
    transcript_consequences {
      gene_id
      gene_symbol
      transcript_id
      hgvsc
      hgvsp
      consequence
      lof
      lof_filter
      lof_flags
    }
    in_silico_predictors {
      id
      value
    }
  }
}
```

#### Gene Query

```graphql
query GeneQuery($geneSymbol: String!, $datasetId: DatasetId!) {
  gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
    gene_id
    symbol
    name
    chrom
    start
    stop
    strand
    canonical_transcript_id
    mane_select_transcript {
      ensembl_id
      refseq_id
    }
    gnomad_constraint {
      exp_lof
      obs_lof
      oe_lof
      oe_lof_lower
      oe_lof_upper
      lof_z
      pLI
    }
    variants(dataset: $datasetId) {
      variantId
      exome {
        ac
        an
        af
      }
    }
  }
}
```

### Available Dataset IDs

| Dataset | Description |
|---------|-------------|
| `gnomad_r4` | gnomAD v4.0 (current) |
| `gnomad_r3` | gnomAD v3.1 |
| `gnomad_r2_1` | gnomAD v2.1.1 |
| `gnomad_r2_1_controls` | Controls subset |
| `gnomad_r2_1_non_neuro` | Non-neurological subset |
| `gnomad_r2_1_non_cancer` | Non-cancer subset |
| `gnomad_r2_1_non_topmed` | Non-TOPMed subset |
| `exac` | ExAC (legacy) |

---

## VCF Schema

### INFO Fields - Allele Counts & Frequencies

| Field | Type | Description |
|-------|------|-------------|
| `AC` | Integer | Allele count in genotypes |
| `AN` | Integer | Total number of alleles |
| `AF` | Float | Allele frequency (AC/AN) |
| `AC_XX` | Integer | Allele count in XX samples |
| `AC_XY` | Integer | Allele count in XY samples |
| `AC_afr` | Integer | Allele count - African ancestry |
| `AC_amr` | Integer | Allele count - Admixed American |
| `AC_asj` | Integer | Allele count - Ashkenazi Jewish |
| `AC_eas` | Integer | Allele count - East Asian |
| `AC_fin` | Integer | Allele count - Finnish |
| `AC_mid` | Integer | Allele count - Middle Eastern |
| `AC_nfe` | Integer | Allele count - Non-Finnish European |
| `AC_sas` | Integer | Allele count - South Asian |
| `nhomalt` | Integer | Count of homozygous individuals |
| `nhomalt_XX` | Integer | Homozygous count in XX samples |

### INFO Fields - Filtering Allele Frequency

| Field | Type | Description |
|-------|------|-------------|
| `faf95` | Float | Filtering AF (95% CI upper bound) |
| `faf99` | Float | Filtering AF (99% CI upper bound) |
| `faf95_afr` | Float | Filtering AF - African (95% CI) |
| `fafmax` | Struct | Maximum filtering AF across populations |
| `fafmax_gnomad` | Float | Maximum FAF across genetic ancestry |

### INFO Fields - Quality Metrics

| Field | Type | Description |
|-------|------|-------------|
| `FS` | Float | Fisher strand bias (phred-scaled) |
| `MQ` | Float | Mapping quality |
| `MQRankSum` | Float | Mapping quality rank sum test |
| `QD` | Float | Quality by depth |
| `ReadPosRankSum` | Float | Read position rank sum test |
| `SOR` | Float | Symmetric odds ratio |
| `VarDP` | Integer | Variant depth |
| `InbreedingCoeff` | Float | Inbreeding coefficient |

### INFO Fields - Allele-Specific Quality

| Field | Type | Description |
|-------|------|-------------|
| `AS_FS` | Float | Allele-specific Fisher strand |
| `AS_MQ` | Float | Allele-specific mapping quality |
| `AS_MQRankSum` | Float | Allele-specific MQ rank sum |
| `AS_QD` | Float | Allele-specific QD |
| `AS_ReadPosRankSum` | Float | Allele-specific read position rank sum |
| `AS_SB_TABLE` | String | Allele-specific strand bias table |
| `AS_SOR` | Float | Allele-specific symmetric odds ratio |
| `AS_pab_max` | Float | Maximum p-value for allele balance |

### INFO Fields - Region Flags

| Field | Type | Description |
|-------|------|-------------|
| `lcr` | Flag | In low-complexity region |
| `decoy` | Flag | In decoy region |
| `segdup` | Flag | In segmental duplication |
| `nonpar` | Flag | In non-pseudoautosomal region |

### INFO Fields - VQSR

| Field | Type | Description |
|-------|------|-------------|
| `VQSLOD` | Float | VQSR log odds score |
| `culprit` | String | Most informative annotation |
| `POSITIVE_TRAIN_SITE` | Flag | Training site (positive) |
| `NEGATIVE_TRAIN_SITE` | Flag | Training site (negative) |

### INFO Fields - VRS (GA4GH)

| Field | Type | Description |
|-------|------|-------------|
| `VRS_Allele_IDs` | String | GA4GH VRS Allele identifiers |
| `VRS_Starts` | Integer | VRS interval start positions |
| `VRS_Ends` | Integer | VRS interval end positions |
| `VRS_States` | String | VRS sequence states |

### FORMAT Fields

| Field | Type | Description |
|-------|------|-------------|
| `GT` | String | Genotype |
| `GQ` | Integer | Genotype quality |
| `DP` | Integer | Read depth |
| `AD` | Integer,Integer | Allelic depths (ref,alt) |
| `PL` | Integer,Integer,Integer | Phred-scaled genotype likelihoods |
| `SB` | Integer,Integer,Integer,Integer | Strand bias components |
| `PGT` | String | Physical phasing haplotype |
| `PID` | String | Physical phasing ID |

---

## Constraint Metrics

### Gene Constraint Fields

| Field | Type | Description |
|-------|------|-------------|
| `exp_lof` | Float | Expected loss-of-function variants |
| `obs_lof` | Float | Observed loss-of-function variants |
| `oe_lof` | Float | Observed/expected ratio (LoF) |
| `oe_lof_lower` | Float | oe_lof 90% CI lower bound |
| `oe_lof_upper` | Float | oe_lof 90% CI upper bound (LOEUF) |
| `lof_z` | Float | LoF Z-score |
| `pLI` | Float | Probability of LoF intolerance |
| `exp_mis` | Float | Expected missense variants |
| `obs_mis` | Float | Observed missense variants |
| `oe_mis` | Float | Observed/expected ratio (missense) |
| `mis_z` | Float | Missense Z-score |
| `exp_syn` | Float | Expected synonymous variants |
| `obs_syn` | Float | Observed synonymous variants |
| `oe_syn` | Float | Observed/expected ratio (synonymous) |
| `syn_z` | Float | Synonymous Z-score |

### Constraint Interpretation

| Metric | Interpretation |
|--------|----------------|
| **LOEUF < 0.35** | LoF-intolerant gene (constrained) |
| **pLI > 0.9** | High probability of haploinsufficiency |
| **mis_z > 3.09** | Missense-constrained |
| **syn_z > 3.09** | Synonymous-constrained (unusual) |

---

## LOFTEE Annotations

Loss-Of-Function Transcript Effect Estimator predicts high-confidence LoF variants.

### LOFTEE Fields

| Field | Type | Description |
|-------|------|-------------|
| `lof` | String | LOFTEE prediction: HC (high confidence), LC (low confidence) |
| `lof_filter` | String | LOFTEE filter flags |
| `lof_flags` | String | Additional LOFTEE warnings |
| `lof_info` | String | LOFTEE information string |

### LOFTEE Filter Values

| Filter | Description |
|--------|-------------|
| `END_TRUNC` | End truncation |
| `INCOMPLETE_CDS` | Incomplete CDS |
| `EXON_INTRON_UNDEF` | Exon/intron boundary undefined |
| `SMALL_INTRON` | Intron too small |
| `ANC_ALLELE` | Ancestral allele |
| `NON_ACCEPTOR_DISRUPTING` | Not acceptor disrupting |
| `NON_DONOR_DISRUPTING` | Not donor disrupting |
| `RESCUE_DONOR` | Rescue donor site |
| `RESCUE_ACCEPTOR` | Rescue acceptor site |
| `GC_TO_GT_DONOR` | GC-to-GT donor variant |
| `5UTR_SPLICE` | 5' UTR splice variant |
| `3UTR_SPLICE` | 3' UTR splice variant |

---

## In Silico Prediction Scores

### Available Predictors

| Predictor | Range | Interpretation |
|-----------|-------|----------------|
| **REVEL** | 0-1 | >0.5 likely pathogenic |
| **CADD** | 0-99 (phred) | >20 top 1%, >30 top 0.1% |
| **phyloP** | -20 to 10 | >2 conserved, <-2 fast-evolving |
| **SIFT** | 0-1 | <0.05 damaging |
| **PolyPhen** | 0-1 | >0.85 probably damaging |
| **SpliceAI** | 0-1 | >0.5 likely splice-altering |
| **Pangolin** | 0-1 | Splice impact score |

### Score Fields in VCF

| Field | Type | Description |
|-------|------|-------------|
| `cadd_phred` | Float | CADD Phred-scaled score |
| `cadd_raw` | Float | CADD raw score |
| `revel` | Float | REVEL score |
| `sift` | Float | SIFT score |
| `polyphen` | Float | PolyPhen-2 score |
| `spliceai_ds_max` | Float | Maximum SpliceAI delta score |
| `spliceai_dp_max` | Integer | Position of max splice effect |
| `pangolin_max` | Float | Maximum Pangolin score |

---

## Quality Histograms

Stored for allele-specific quality assessment:

| Histogram | Description |
|-----------|-------------|
| `gq_hist_alt` | Genotype quality for alt allele carriers |
| `gq_hist_all` | Genotype quality for all genotypes |
| `dp_hist_alt` | Depth for alt allele carriers |
| `dp_hist_all` | Depth for all genotypes |
| `ab_hist_alt` | Allele balance for heterozygotes |
| `age_hist_het` | Age distribution of heterozygotes |
| `age_hist_hom` | Age distribution of homozygotes |

---

## Download Options

### Cloud Providers

| Provider | Location |
|----------|----------|
| **AWS** | `s3://gnomad-public-us-east-1/` |
| **Google Cloud** | `gs://gcp-public-data--gnomad/` |
| **Azure** | Available via Azure Open Datasets |

### File Types

| Format | Use Case | Size (v4 exomes) |
|--------|----------|------------------|
| `.vcf.bgz` | VCF with tabix index | ~300 GB |
| `.ht` | Hail Table format | Variable |
| `.vds` | Variant Dataset (efficient) | ~18 TB total |
| `.tsv.bgz` | Tab-separated summaries | Variable |

### Download URLs (v4.0)

```
# Exome VCF (GRCh38)
gs://gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.vcf.bgz

# Genome VCF (GRCh38)
gs://gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.vcf.bgz

# Constraint metrics
gs://gcp-public-data--gnomad/release/4.0/constraint/gnomad.v4.0.constraint_metrics.tsv
```

---

## Sample Data

### GraphQL Response

```json
{
  "data": {
    "variant": {
      "variantId": "11-5227002-T-A",
      "chrom": "11",
      "pos": 5227002,
      "ref": "T",
      "alt": "A",
      "exome": {
        "ac": 156,
        "an": 1420000,
        "af": 0.00011,
        "ac_hom": 12,
        "populations": [
          {"id": "afr", "ac": 89, "an": 45000, "ac_hom": 8},
          {"id": "nfe", "ac": 12, "an": 800000, "ac_hom": 0},
          {"id": "sas", "ac": 45, "an": 100000, "ac_hom": 4}
        ]
      },
      "genome": {
        "ac": 42,
        "an": 152000,
        "af": 0.00028
      },
      "joint": {
        "ac": 198,
        "an": 1572000,
        "af": 0.000126
      },
      "flags": ["lcr"],
      "lof": "HC",
      "lof_filter": null,
      "transcript_consequences": [
        {
          "gene_id": "ENSG00000244734",
          "gene_symbol": "HBB",
          "transcript_id": "ENST00000335295",
          "hgvsc": "c.20A>T",
          "hgvsp": "p.Glu7Val",
          "consequence": "missense_variant",
          "lof": null
        }
      ],
      "in_silico_predictors": [
        {"id": "revel", "value": "0.891"},
        {"id": "cadd_phred", "value": "25.3"},
        {"id": "sift", "value": "0.01"},
        {"id": "polyphen", "value": "0.998"}
      ]
    }
  }
}
```

### VCF Record

```
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AC_afr,Number=A,Type=Integer,Description="Allele count in African ancestry">
##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="Count of homozygous individuals">
##INFO=<ID=faf95,Number=A,Type=Float,Description="Filtering allele frequency (95% CI)">
##INFO=<ID=cadd_phred,Number=A,Type=Float,Description="CADD Phred score">
##INFO=<ID=revel,Number=A,Type=Float,Description="REVEL score">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
11	5227002	rs334	T	A	.	PASS	AC=198;AN=1572000;AF=0.000126;AC_afr=89;AC_nfe=12;nhomalt=12;faf95=0.00015;cadd_phred=25.3;revel=0.891
```

---

## Integration Examples

### Python: GraphQL Query

```python
import requests

def query_gnomad_variant(variant_id: str, dataset: str = "gnomad_r4") -> dict:
    """Query gnomAD for variant information via GraphQL."""
    query = """
    query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
      variant(variantId: $variantId, dataset: $datasetId) {
        variantId
        chrom
        pos
        ref
        alt
        exome { ac an af ac_hom }
        genome { ac an af ac_hom }
        joint { ac an af }
        lof
        lof_filter
        transcript_consequences {
          gene_symbol
          hgvsc
          hgvsp
          consequence
        }
        in_silico_predictors { id value }
      }
    }
    """

    response = requests.post(
        "https://gnomad.broadinstitute.org/api",
        json={
            "query": query,
            "variables": {
                "variantId": variant_id,
                "datasetId": dataset
            }
        }
    )
    response.raise_for_status()
    return response.json()

# Example: Query sickle cell variant
result = query_gnomad_variant("11-5227002-T-A")
variant = result["data"]["variant"]
print(f"Allele frequency: {variant['joint']['af']}")
```

### Python: Parse VCF

```python
import gzip

def parse_gnomad_vcf(filepath: str):
    """Parse gnomAD VCF file."""
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom, pos, id_, ref, alt, qual, filter_, info = fields[:8]

            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = True

            yield {
                'chrom': chrom,
                'pos': int(pos),
                'id': id_,
                'ref': ref,
                'alt': alt,
                'filter': filter_,
                'ac': int(info_dict.get('AC', 0)),
                'an': int(info_dict.get('AN', 0)),
                'af': float(info_dict.get('AF', 0)),
                'ac_afr': int(info_dict.get('AC_afr', 0)),
                'nhomalt': int(info_dict.get('nhomalt', 0)),
                'cadd_phred': float(info_dict.get('cadd_phred', 0)) if 'cadd_phred' in info_dict else None,
                'revel': float(info_dict.get('revel', 0)) if 'revel' in info_dict else None
            }
```

### Calculate Population-Specific Frequency

```python
def calculate_population_af(ac: int, an: int) -> float:
    """Calculate allele frequency from counts."""
    if an == 0:
        return 0.0
    return ac / an

# Note: gnomAD GraphQL returns AC/AN per population,
# but not AF directly for populations in some queries
# Calculate: AF = AC/AN
```

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `AC` | Allele Count - number of times the alternate allele is observed in the dataset | 198 |
| `AN` | Allele Number - total number of alleles genotyped at this position | 1572000 |
| `AF` | Allele Frequency - proportion of alleles that are the alternate (AC/AN) | 0.000126 |
| `nhomalt` | Number of individuals homozygous for the alternate allele | 12 |
| `faf95` | Filtering Allele Frequency at 95% confidence, upper bound for pathogenicity assessment | 0.00015 |
| `LOEUF` | Loss-of-function Observed/Expected Upper bound Fraction, key constraint metric | 0.35 (constrained) |
| `pLI` | Probability of Loss-of-function Intolerance, likelihood gene is haploinsufficient | 0.99 |
| `CADD` | Combined Annotation Dependent Depletion score for variant deleteriousness | 25.3 (phred-scaled) |
| `REVEL` | Rare Exome Variant Ensemble Learner score for missense pathogenicity | 0.891 |
| `lof` | Loss-of-function classification from LOFTEE (HC=high confidence, LC=low confidence) | HC |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Allele Frequency | Proportion of chromosomes carrying a specific allele in a population | AF, AC, AN fields |
| Constraint | Degree to which a gene tolerates variation, measured by observed vs expected variants | pLI, LOEUF, oe_lof |
| Loss-of-Function (LoF) | Variant that abolishes or severely impairs gene function (nonsense, frameshift, splice) | LOFTEE annotations |
| Haploinsufficiency | Condition where one functional gene copy is insufficient for normal function | pLI score |
| Filtering Allele Frequency | Conservative AF estimate used to filter variants in clinical settings | faf95, faf99 |
| Genetic Ancestry | Population grouping based on genetic similarity, not self-reported ethnicity | afr, eas, nfe, etc. |
| VQSR | Variant Quality Score Recalibration, machine learning filter for variant quality | VQSLOD field |
| Missense Constraint | Measure of selection against missense variants in a gene | mis_z, oe_mis |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| gnomAD | Genome Aggregation Database | Broad Institute resource |
| LoF | Loss-of-Function | Variant abolishing gene function |
| LOEUF | LoF Observed/Expected Upper bound Fraction | Constraint metric (lower = more constrained) |
| pLI | Probability of LoF Intolerance | Legacy constraint metric |
| LOFTEE | Loss-Of-Function Transcript Effect Estimator | LoF annotation tool |
| CADD | Combined Annotation Dependent Depletion | Deleteriousness score |
| REVEL | Rare Exome Variant Ensemble Learner | Missense pathogenicity |
| SpliceAI | Splice Artificial Intelligence predictor | Splice impact score |
| VCF | Variant Call Format | Standard variant file format |
| VQSR | Variant Quality Score Recalibration | Quality filtering method |
| VRS | Variant Representation Specification | GA4GH standard |
| GRCh38 | Genome Reference Consortium Human Build 38 | Current reference genome |
| QD | Quality by Depth | Quality metric |
| FS | Fisher Strand | Strand bias metric |
| MQ | Mapping Quality | Alignment quality |
| AF | Allele Frequency | Proportion of alleles |
| CI | Confidence Interval | Statistical range |

---

## References

1. gnomAD Browser: https://gnomad.broadinstitute.org
2. gnomAD v4 Announcement: https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/
3. gnomAD Methods Documentation: https://broadinstitute.github.io/gnomad_methods/
4. Karczewski KJ, et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. Nature. 581(7809):434-443.
5. LOFTEE: https://github.com/konradjk/loftee

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema from v4.0 documentation |
