# gnomAD - LLM Context Reference

> Population allele frequency database: 807K individuals, 8 ancestry groups, constraint metrics for variant interpretation.

## Quick Reference

| Field | Value |
|-------|-------|
| **URL** | https://gnomad.broadinstitute.org |
| **Maintainer** | Broad Institute |
| **License** | ODC-BY 1.0 (Open) |
| **Commercial OK** | Yes |
| **Update Freq** | Major releases |
| **Version** | v4.0 (Nov 2023) |

## Content Summary

gnomAD aggregates exome and genome sequencing from 807K individuals to provide population-level allele frequencies stratified by 8 genetic ancestry groups. Essential for filtering common variants and assessing gene constraint (pLI, LOEUF scores).

**Record Counts:**
- Individuals: 807,162
- SNVs: 62.9M (passing QC)
- Indels: 6.2M (passing QC)
- Genes with constraint: 19,704

## Key Identifiers

| ID Type | Format | Example |
|---------|--------|---------|
| Variant ID | chr-pos-ref-alt | 1-55516888-G-A |
| Gene ID | ENSG + digits | ENSG00000141510 |

**Cross-references:** rsID (dbSNP), Ensembl, HGNC, OMIM

## Core Schema

### Primary Entity: Variant

| Field | Type | Description |
|-------|------|-------------|
| variantId | string | chr-pos-ref-alt |
| ac | int | Allele count |
| an | int | Allele number (chromosomes) |
| af | float | Allele frequency (ac/an) |
| ac_hom | int | Homozygote count |
| flags | array | QC flags (lcr, segdup, etc.) |

### Key Relationships

```
Variant --[transcript_consequences]--> Transcript
Transcript --[gene_id]--> Gene
Gene --[gnomad_constraint]--> ConstraintMetrics
Variant --[populations]--> PopulationFrequencies
```

## Access Methods

### API (GraphQL)
```
POST https://gnomad.broadinstitute.org/api
Content-Type: application/json
```
Rate limit: 1 req/sec recommended

### Bulk Download
```
gs://gcp-public-data--gnomad/release/4.0/vcf/
```
Format: VCF, Hail VariantDataset (VDS)
Size: 18 TB (full VDS), ~50 GB (VCF per chromosome)

## Query Examples

### Get variant frequency
```graphql
query {
  variant(variantId: "1-55516888-G-A", dataset: gnomad_r4) {
    variantId
    exome { af ac an ac_hom }
    genome { af ac an ac_hom }
    populations { id ac an }
  }
}
```

### Get gene constraint
```graphql
query {
  gene(gene_symbol: "BRCA1", reference_genome: GRCh38) {
    gnomad_constraint {
      pLI
      oe_lof
      oe_lof_lower
      oe_lof_upper
    }
  }
}
```

## Sample Record

```json
{
  "variantId": "1-55516888-G-A",
  "exome": {
    "ac": 1234,
    "an": 1461870,
    "af": 0.000844,
    "ac_hom": 2
  },
  "populations": [
    {"id": "nfe", "ac": 800, "an": 900000},
    {"id": "afr", "ac": 150, "an": 120000},
    {"id": "eas", "ac": 50, "an": 80000}
  ],
  "transcript_consequences": [{
    "gene_symbol": "PCSK9",
    "consequence": "missense_variant",
    "hgvsp": "p.Arg46Leu"
  }]
}
```

## Integration Notes

- **Primary use:** Filter common variants, assess variant rarity
- **Best for:** Population frequency lookup, gene constraint scores
- **Limitations:** Primarily adult, healthy individuals; ancestry bias (60% European v3)
- **Pairs with:** ClinVar (pathogenicity), dbNSFP (predictions), OMIM (disease)

---
*Size: ~50GB (VCF) | Updated: January 2026*
