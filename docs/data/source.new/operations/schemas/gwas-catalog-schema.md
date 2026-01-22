---
id: schema-gwas-catalog
title: GWAS Catalog REST API Schema
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, gwas, genomics]
---

**Parent:** [Schema Documentation](./_index.md)

# GWAS Catalog REST API Schema

**Document ID:** GWAS-CATALOG-SCHEMA
**Source:** EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/rest/api)
**Last Updated:** January 2026
**API Version:** REST API (current)

---

## Database Statistics (January 2026)

| Metric | Count |
|--------|-------|
| Studies | 186,237 |
| Associations | 1,058,471 |
| SNPs | 512,069 |
| Summary Statistics Datasets | 155,485 |
| EFO Traits | 21,004 |

**Technical Details:**
- Genome Build: GRCh38.p14
- EFO Version: v3.86.0
- Ensembl Build: 115
- dbSNP Build: 156

---

## REST API Endpoints

### Base URL
```
https://www.ebi.ac.uk/gwas/rest/api
```

### Available Resources

| Endpoint | URL | Description |
|----------|-----|-------------|
| Studies | `/studies` | GWAS study metadata and publications |
| Associations | `/associations` | Genetic association results |
| SNPs | `/singleNucleotidePolymorphisms` | Variant information |
| EFO Traits | `/efoTraits` | Experimental Factor Ontology traits |
| Unpublished Studies | `/unpublished-studies` | Pre-publication studies |
| Profile | `/profile` | API profile information |

### Common Query Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `page` | integer | Page number (0-indexed) |
| `size` | integer | Results per page (default: 20) |
| `sort` | string | Sort field and direction |
| `projection` | string | Field projection (not available for unpublished-studies) |

---

## Study Schema

### Core Fields

```json
{
  "accessionId": "GCST000854",
  "diseaseTrait": {
    "trait": "Suicide risk"
  },
  "initialSampleSize": "3,117 European ancestry Bipolar disorder cases...",
  "snpCount": 1922309,
  "imputed": true,
  "pooled": false,
  "gxe": false,
  "gxg": false,
  "fullPvalueSet": false
}
```

### Field Specifications

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `accessionId` | string | Yes | GCST accession number (e.g., "GCST000854") |
| `diseaseTrait.trait` | string | Yes | Reported trait name |
| `initialSampleSize` | string | Yes | Sample description text |
| `snpCount` | integer | No | Number of SNPs analyzed |
| `imputed` | boolean | Yes | Whether genotypes were imputed |
| `pooled` | boolean | Yes | Whether samples were pooled |
| `gxe` | boolean | Yes | Gene-environment interaction study |
| `gxg` | boolean | Yes | Gene-gene interaction study |
| `fullPvalueSet` | boolean | Yes | Full summary statistics available |

### Publication Information

```json
{
  "publicationInfo": {
    "pubmedId": "21041247",
    "publicationDate": "2010-11-01",
    "publication": "Am J Psychiatry",
    "title": "Genome-wide association study of...",
    "author": {
      "fullname": "Perlis RH",
      "orcid": "0000-0001-5608-2293"
    }
  }
}
```

| Field | Type | Description |
|-------|------|-------------|
| `pubmedId` | string | PubMed ID |
| `publicationDate` | string (ISO 8601) | Publication date |
| `publication` | string | Journal name |
| `title` | string | Article title |
| `author.fullname` | string | First author name |
| `author.orcid` | string/null | ORCID identifier |

### Ancestry Information

```json
{
  "ancestries": [
    {
      "type": "initial",
      "ancestralGroups": [
        {
          "ancestralGroup": "European"
        }
      ],
      "countryOfRecruitment": [
        {
          "countryName": "United Kingdom"
        }
      ],
      "numberOfIndividuals": 96598
    }
  ]
}
```

---

## Association Schema

### Core Fields

```json
{
  "pvalue": 7.0e-8,
  "pvalueDescription": "(progression)",
  "snpType": "novel",
  "riskFrequency": 0.41,
  "orPerCopyNum": 1.105,
  "betaNum": null,
  "betaUnit": null,
  "betaDirection": null,
  "range": "[1.07-1.14]"
}
```

### Field Specifications

| Field | Type | Description |
|-------|------|-------------|
| `pvalue` | float (scientific notation) | Association p-value |
| `pvalueDescription` | string | Context (e.g., "(SBP)", "(progression)") |
| `snpType` | string | "novel" or "known" classification |
| `riskFrequency` | float (0-1) | Risk allele frequency in population |
| `orPerCopyNum` | float | Odds ratio per allele copy |
| `betaNum` | float/null | Effect size (for quantitative traits) |
| `betaUnit` | string/null | Unit of effect size |
| `betaDirection` | string/null | Direction of effect |
| `range` | string | Confidence interval notation |

### Loci Structure

```json
{
  "loci": [
    {
      "description": "Single variant",
      "strongestRiskAlleles": [
        {
          "riskAlleleName": "rs9497975-?"
        }
      ],
      "authorReportedGenes": [
        {
          "geneName": "TAGAP",
          "entrezGeneIds": [
            {
              "entrezGeneId": "117289"
            }
          ],
          "ensemblGeneIds": [
            {
              "ensemblGeneId": "ENSG00000164691"
            }
          ]
        }
      ]
    }
  ]
}
```

---

## SNP (Single Nucleotide Polymorphism) Schema

### Core Fields

```json
{
  "rsId": "rs7329174",
  "merged": 0,
  "functionalClass": "intergenic_variant",
  "lastUpdateDate": "2024-03-15T00:00:00.000+00:00"
}
```

### Location Information

```json
{
  "locations": [
    {
      "chromosomeName": "13",
      "chromosomePosition": 40983974,
      "region": {
        "name": "13q14.11"
      }
    }
  ]
}
```

### Genomic Contexts

```json
{
  "genomicContexts": [
    {
      "isIntergenic": false,
      "isUpstream": true,
      "isDownstream": false,
      "distance": 4523,
      "source": "Ensembl",
      "mappingMethod": "Ensembl_pipeline",
      "isClosestGene": true,
      "gene": {
        "geneName": "ELF1",
        "entrezGeneIds": [
          {"entrezGeneId": "1997"}
        ],
        "ensemblGeneIds": [
          {"ensemblGeneId": "ENSG00000120690"}
        ]
      }
    }
  ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `isIntergenic` | boolean | Located between genes |
| `isUpstream` | boolean | Upstream of gene |
| `isDownstream` | boolean | Downstream of gene |
| `distance` | integer | Distance in base pairs |
| `source` | string | Annotation source (Ensembl/NCBI) |
| `mappingMethod` | string | Method used for mapping |
| `isClosestGene` | boolean | Nearest gene flag |

---

## EFO Trait Schema

### Sample Response

```json
{
  "trait": "Celiac disease",
  "uri": "http://www.ebi.ac.uk/efo/EFO_0001060",
  "shortForm": "EFO_0001060",
  "_links": {
    "studies": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/EFO_0001060/studies"
    },
    "associations": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/EFO_0001060/associations"
    }
  }
}
```

### Trait Categories (Sample)

| Category | Example Traits |
|----------|----------------|
| Disease Traits | Celiac disease (EFO_0001060), endometriosis (EFO_0001065), malaria (EFO_0001068), lung carcinoma (EFO_0001071) |
| Measurement Traits | MHPG measurement, thyroxine amount, antioxidant measurement, CCL5 measurement, IGFBP-1 measurement |
| Morphological Traits | Arm span, gestational age, head circumference, body composition measurement |
| Physiological Traits | Energy expenditure, fatty acid amount, folic acid amount, energy intake |

---

## HATEOAS Navigation Links

All responses include hypermedia links for navigation:

```json
{
  "_links": {
    "self": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000854"
    },
    "study": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000854"
    },
    "associations": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000854/associations"
    },
    "snps": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000854/snps"
    },
    "efoTraits": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000854/efoTraits"
    }
  }
}
```

---

## Pagination Response

```json
{
  "page": {
    "size": 20,
    "totalElements": 186237,
    "totalPages": 9312,
    "number": 0
  }
}
```

---

## Example API Calls

### Fetch Studies by PubMed ID
```bash
curl "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId?pubmedId=20686565"
```

### Fetch SNP Data
```bash
curl "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs7329174"
```

### Fetch Associations with Pagination
```bash
curl "https://www.ebi.ac.uk/gwas/rest/api/associations?page=0&size=100"
```

### Fetch EFO Traits
```bash
curl "https://www.ebi.ac.uk/gwas/rest/api/efoTraits"
```

---

## Integration Notes

### EFO Ontology Integration
- All traits use Experimental Factor Ontology (EFO) for standardization
- EFO URIs enable cross-database linking
- Supports hierarchical trait classification

### Related APIs
- Summary Statistics API: `https://www.ebi.ac.uk/gwas/summary-statistics/api`
- FTP Downloads: `http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/`

### Rate Limiting
- No explicit rate limits documented
- Reasonable use expected

### Data Formats
- Primary: JSON (HAL format with hypermedia links)
- Bulk: TSV, OWL/RDF via FTP

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `accessionId` | GWAS Catalog study accession number | GCST000854 |
| `pvalue` | Statistical significance of genetic association | 7.0e-8 |
| `orPerCopyNum` | Odds ratio per copy of risk allele | 1.105 |
| `betaNum` | Effect size for quantitative traits | 0.15 |
| `riskFrequency` | Frequency of risk allele in study population | 0.41 |
| `rsId` | RefSNP identifier from dbSNP | rs7329174 |
| `functionalClass` | Sequence Ontology functional annotation | intergenic_variant |
| `chromosomePosition` | Genomic coordinate on chromosome | 40983974 |
| `EFO_XXXXXXX` | Experimental Factor Ontology trait identifier | EFO_0001060 |
| `pubmedId` | PubMed literature reference identifier | 21041247 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| GWAS | Genome-Wide Association Study | Study type |
| Association | Statistical link between genetic variant and trait | Core data type |
| Risk Allele | Genetic variant associated with increased trait/disease risk | riskAlleleName |
| Odds Ratio | Relative odds of outcome given exposure to risk allele | orPerCopyNum |
| Effect Size | Magnitude of genetic variant's influence on trait | betaNum |
| Locus | Genomic region containing associated variant(s) | loci structure |
| Summary Statistics | Full statistical results from GWAS analysis | fullPvalueSet |
| EFO Trait | Standardized trait term from Experimental Factor Ontology | Trait classification |
| Ancestry | Population genetic background of study participants | ancestralGroups |
| Imputation | Statistical inference of ungenotyped variants | imputed field |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| GWAS | Genome-Wide Association Study | Study methodology |
| SNP | Single Nucleotide Polymorphism | Variant type |
| EFO | Experimental Factor Ontology | Trait vocabulary |
| OR | Odds Ratio | Effect measure |
| CI | Confidence Interval | Statistical range |
| MAF | Minor Allele Frequency | Population metric |
| LD | Linkage Disequilibrium | Variant correlation |
| GCST | GWAS Catalog Study | Study identifier prefix |
| GRCh38 | Genome Reference Consortium Human Build 38 | Reference genome |
| dbSNP | Database of Single Nucleotide Polymorphisms | Variant database |
| TSV | Tab-Separated Values | File format |
| OWL | Web Ontology Language | Ontology format |
| RDF | Resource Description Framework | Semantic web format |
| HAL | Hypertext Application Language | API response format |
| HATEOAS | Hypermedia as the Engine of Application State | REST principle |

---

## References

- GWAS Catalog: https://www.ebi.ac.uk/gwas/
- API Documentation: https://www.ebi.ac.uk/gwas/rest/api
- EFO Ontology: https://www.ebi.ac.uk/efo/
