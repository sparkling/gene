# ClinVar - LLM Context Reference

> Clinical variant interpretations: pathogenicity classifications from 2000+ submitters with review status.

## Quick Reference

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar |
| **Maintainer** | NCBI/NIH |
| **License** | Public Domain |
| **Commercial OK** | Yes |
| **Update Freq** | Weekly (Mon) |
| **Version** | Current release |

## Content Summary

ClinVar aggregates clinical interpretations of genetic variants from diagnostic labs, research groups, and expert panels. Each variant has aggregated classification (VCV) from individual submissions (SCV) with review status indicating evidence strength.

**Record Counts:**
- Variant records (VCV): 2.5M+
- Submissions (SCV): 4M+
- Submitting orgs: 2,000+
- Genes: 43,000+

## Key Identifiers

| ID Type | Format | Example |
|---------|--------|---------|
| VCV | VCV + 9 digits + version | VCV000012345.3 |
| RCV | RCV + 9 digits + version | RCV000012345.6 |
| Variation ID | Integer | 12345 |

**Cross-references:** rsID (dbSNP), OMIM, MedGen CUI, Gene ID

## Core Schema

### Primary Entity: VariationArchive (VCV)

| Field | Type | Description |
|-------|------|-------------|
| VariationID | int | Internal variant ID |
| Name | string | HGVS expression |
| ClassifiedRecord.Classification | string | Aggregated interpretation |
| ReviewStatus | string | Evidence level (0-4 stars) |
| ConditionList | array | Associated diseases |

### Key Relationships

```
VCV (aggregate) --[contains]--> RCV (variant+condition)
RCV --[from]--> SCV (individual submission)
Variant --[in_gene]--> Gene
Variant --[associated_with]--> Condition (MedGen)
```

## Access Methods

### API (E-utilities)
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={query}
```
Rate limit: 3/sec (without API key), 10/sec (with key)

### Bulk Download
```
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
```
Format: VCF, XML, TSV (variant_summary.txt.gz)
Size: 177 MB (VCF), 2.5 GB (full XML)

## Query Examples

### Get variant by rsID (E-utilities)
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=rs121913279[variantid]&retmode=json
```

### VCF INFO fields
```
##INFO=<ID=CLNSIG,Description="Clinical significance">
##INFO=<ID=CLNREVSTAT,Description="Review status">
##INFO=<ID=CLNDN,Description="Disease name">
##INFO=<ID=GENEINFO,Description="Gene:GeneID">
```

## Sample Record

```json
{
  "VariationID": 17661,
  "Name": "NM_000546.6(TP53):c.743G>A (p.Arg248Gln)",
  "GeneSymbol": "TP53",
  "Classification": {
    "GermlineClassification": "Pathogenic",
    "ReviewStatus": "criteria provided, multiple submitters, no conflicts"
  },
  "Conditions": [
    {"Name": "Li-Fraumeni syndrome", "MedGenCUI": "C0085390"}
  ],
  "rsID": "rs28934576"
}
```

## Classification Values

| Classification | Meaning |
|---------------|---------|
| Pathogenic | Disease-causing |
| Likely pathogenic | >90% certainty pathogenic |
| Uncertain significance (VUS) | Insufficient evidence |
| Likely benign | >90% certainty benign |
| Benign | Not disease-causing |

## Review Status (Stars)

| Stars | Status | Meaning |
|-------|--------|---------|
| 4 | practice guideline | Expert panel guideline |
| 3 | reviewed by expert panel | Expert review |
| 2 | criteria provided, multiple submitters | Consensus |
| 1 | criteria provided, single submitter | Single lab |
| 0 | no assertion criteria | No criteria |

## Integration Notes

- **Primary use:** Clinical variant interpretation lookup
- **Best for:** Pathogenicity assessment, disease associations
- **Limitations:** Submitter disagreements exist; review status varies
- **Pairs with:** gnomAD (frequency), dbNSFP (predictions), OMIM (disease detail)

---
*Size: ~2.5GB | Updated: January 2026*
