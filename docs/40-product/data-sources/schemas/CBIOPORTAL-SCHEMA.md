# cBioPortal REST API Schema

**Document ID:** CBIOPORTAL-SCHEMA
**Source:** cBioPortal for Cancer Genomics (https://www.cbioportal.org/api)
**Last Updated:** January 2026
**API Version:** REST API v1

---

## Database Statistics (January 2026)

| Metric | Count |
|--------|-------|
| Cancer Studies | 423 |
| Cancer Types | 800+ |
| Genes | 4,000+ |
| Molecular Profiles | Multiple per study |

---

## REST API Endpoints

### Base URL
```
https://www.cbioportal.org/api
```

### Primary Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/studies` | GET | List all cancer studies |
| `/cancer-types` | GET | Cancer type taxonomy |
| `/genes` | GET | Gene catalog |
| `/molecular-profiles` | GET | Molecular data profiles |
| `/samples` | GET | Sample data |
| `/patients` | GET | Patient data |
| `/clinical-attributes` | GET | Clinical attribute definitions |
| `/mutations/fetch` | POST | Mutation data retrieval |

---

## Study Schema

### Core Fields

```json
{
  "name": "Glioblastoma (TCGA, Nature 2008)",
  "description": "Comprehensive genomic characterization...",
  "publicStudy": true,
  "pmid": "18772890",
  "citation": "TCGA, Nature 2008",
  "groups": "PUBLIC",
  "status": 0,
  "importDate": "2025-10-21 00:00:00",
  "allSampleCount": 206,
  "readPermission": true,
  "studyId": "gbm_tcga_pub",
  "cancerTypeId": "difg",
  "referenceGenome": "hg19"
}
```

### Field Specifications

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `studyId` | string | Yes | Unique study identifier |
| `name` | string | Yes | Study display name |
| `description` | string | Yes | Full study description |
| `publicStudy` | boolean | Yes | Public accessibility flag |
| `pmid` | string | No | PubMed ID(s), comma-separated |
| `citation` | string | No | Citation text |
| `groups` | string | No | Access groups, semicolon-separated |
| `status` | integer | Yes | Study status code |
| `importDate` | string | Yes | Data import timestamp (YYYY-MM-DD HH:MM:SS) |
| `allSampleCount` | integer | Yes | Total samples in study |
| `readPermission` | boolean | Yes | User read access |
| `cancerTypeId` | string | Yes | Cancer type classification |
| `referenceGenome` | string | Yes | Genomic reference (hg19/hg38) |

### Resource Counts (Nested)

```json
{
  "resourceCounts": [
    {
      "resourceId": "PATHOLOGY_SLIDE",
      "displayName": "Pathology Slide",
      "description": "Digital pathology images",
      "resourceType": "SAMPLE",
      "priority": "1",
      "openByDefault": true,
      "sampleCount": 143,
      "patientCount": 143,
      "studyId": "gbm_tcga_pub"
    }
  ]
}
```

---

## Cancer Type Schema

### Hierarchical Taxonomy

```json
{
  "name": "Acute Myeloid Leukemia",
  "shortName": "AML",
  "cancerTypeId": "aml",
  "dedicatedColor": "Orange",
  "parent": "myeloid"
}
```

### Field Specifications

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Full cancer type designation |
| `shortName` | string | Abbreviated identifier |
| `cancerTypeId` | string | Unique identifier (lowercase, hyphenated) |
| `dedicatedColor` | string | Visual representation color |
| `parent` | string | Parent category ID |

### Hierarchy Levels

| Level | Examples |
|-------|----------|
| Primary (Tissue) | Brain, Breast, Kidney, Liver, Lung, Lymphoid, Myeloid, Skin, Soft Tissue |
| Secondary (Family) | Acute Myeloid Leukemia, Hodgkin Lymphoma, Medulloblastoma |
| Tertiary (Subtype) | AML with Maturation, Classical Hodgkin Lymphoma, SHH/WNT subtypes |

### Subtype Specifications

The taxonomy supports:
- Genetic marker subtypes (e.g., "BCR-ABL1", "t(9;22)")
- Morphological variants (adenocarcinoma, squamous cell, carcinoid)
- Grade classifications (Grade 2, Grade 3, anaplastic)
- Site-specific lesions (in situ variants, borderline tumors)

---

## Gene Schema

```json
{
  "entrezGeneId": 1,
  "hugoGeneSymbol": "A1BG",
  "type": "protein-coding"
}
```

### Field Specifications

| Field | Type | Description |
|-------|------|-------------|
| `entrezGeneId` | integer | NCBI Entrez Gene ID |
| `hugoGeneSymbol` | string | Official HUGO gene symbol |
| `type` | string | Gene classification |

### Gene Type Categories

| Type | Description |
|------|-------------|
| `protein-coding` | Protein-coding genes |
| `ncRNA` | Non-coding RNA |
| `pseudogene` | Pseudogenes |
| `phosphoprotein` | Phosphoprotein markers |
| `other` | Other gene types |
| `unknown` | Unclassified |

---

## Molecular Profile Schema

### Profile Types

```json
{
  "molecularProfileId": "gbm_tcga_pub_mutations",
  "studyId": "gbm_tcga_pub",
  "molecularAlterationType": "MUTATION_EXTENDED",
  "datatype": "MAF",
  "name": "Mutations",
  "description": "Mutation data from whole exome sequencing",
  "showProfileInAnalysisTab": true,
  "patientLevel": false
}
```

### Molecular Alteration Types

| Type | Description | Data Format |
|------|-------------|-------------|
| `MUTATION_EXTENDED` | Somatic mutations | MAF |
| `COPY_NUMBER_ALTERATION` | CNV data | DISCRETE or LOG2-VALUE |
| `MRNA_EXPRESSION` | Gene expression | CONTINUOUS or Z-SCORE |
| `STRUCTURAL_VARIANT` | Structural variants | SV |
| `PROTEIN_LEVEL` | Protein expression | RPPA |
| `METHYLATION` | DNA methylation | HM27/HM450 |
| `GENERIC_ASSAY` | Other assays | Various |

### Copy Number Discrete Values

| Value | Interpretation |
|-------|----------------|
| -2 | Homozygous deletion |
| -1 | Hemizygous deletion |
| 0 | Neutral / no change |
| 1 | Gain |
| 2 | High-level amplification |

### Expression Data Types

| Type | Description |
|------|-------------|
| CONTINUOUS | Raw measurements (FPKM, TPM, RSEM, CPM) |
| Z-SCORE | Normalized relative to diploid or all samples |

---

## Clinical Attribute Schema

```json
{
  "displayName": "Age at Diagnosis",
  "description": "Patient age at initial diagnosis",
  "datatype": "NUMBER",
  "patientAttribute": true,
  "priority": "1",
  "clinicalAttributeId": "AGE",
  "studyId": "gbm_tcga_pub"
}
```

### Field Specifications

| Field | Type | Description |
|-------|------|-------------|
| `displayName` | string | User-friendly label |
| `description` | string | Attribute explanation |
| `datatype` | string | STRING or NUMBER |
| `patientAttribute` | boolean | Patient-level (true) or sample-level (false) |
| `priority` | string | Display priority (numeric) |
| `clinicalAttributeId` | string | Standardized ID |
| `studyId` | string | Associated study |

### Attribute Categories

| Category | Examples |
|----------|----------|
| Demographic | Age at diagnosis, age at procurement, age at metastasis |
| Genomic Markers | 1p/19q status, ploidy, purity, copy number alterations |
| Treatment History | Adjuvant chemotherapy, radiation, immunotherapy, targeted agents |
| Clinical Outcomes | RECIST response, lesion decline, ctDNA recurrence status |
| Tumor Characteristics | Histological findings, anatomical location, subtype |

---

## Mutation Data Schema (MAF Format)

### Standard MAF Fields

| Field | Type | Description |
|-------|------|-------------|
| `Hugo_Symbol` | string | HUGO gene symbol |
| `Entrez_Gene_Id` | integer | NCBI Gene ID |
| `Center` | string | Sequencing center |
| `NCBI_Build` | string | Reference genome build |
| `Chromosome` | string | Chromosome number |
| `Start_Position` | integer | Start genomic coordinate |
| `End_Position` | integer | End genomic coordinate |
| `Strand` | string | DNA strand (+/-) |
| `Variant_Classification` | string | Mutation type |
| `Variant_Type` | string | SNP/INS/DEL |
| `Reference_Allele` | string | Reference nucleotide(s) |
| `Tumor_Seq_Allele1` | string | Tumor allele 1 |
| `Tumor_Seq_Allele2` | string | Tumor allele 2 |
| `Tumor_Sample_Barcode` | string | Sample identifier |

### Variant Classifications

| Classification | Description |
|----------------|-------------|
| Frame_Shift_Del | Deletion causing frameshift |
| Frame_Shift_Ins | Insertion causing frameshift |
| In_Frame_Del | In-frame deletion |
| In_Frame_Ins | In-frame insertion |
| Missense_Mutation | Single amino acid change |
| Nonsense_Mutation | Premature stop codon |
| Splice_Site | Splice junction alteration |
| Silent | Synonymous change |
| 3'UTR | 3' untranslated region |
| 5'UTR | 5' untranslated region |
| Intron | Intronic variant |

---

## Sample Study: GBM TCGA

### Study Summary

| Attribute | Value |
|-----------|-------|
| Study ID | gbm_tcga_pub |
| Name | Glioblastoma (TCGA, Nature 2008) |
| Cancer Type | Diffuse Glioma (DIFG) |
| Reference Genome | hg19 |
| Total Samples | 206 |
| PMID | 18772890 |

### Available Molecular Profiles

| Profile Type | Sample Count |
|--------------|--------------|
| Targeted Sequencing | 91 |
| Copy Number Analysis | 206 |
| Microarray Expression | 206 |
| microRNA Analysis | 206 |
| DNA Methylation (HM27) | 206 |

---

## API Usage Examples

### Fetch All Studies
```bash
curl "https://www.cbioportal.org/api/studies"
```

### Fetch Cancer Types
```bash
curl "https://www.cbioportal.org/api/cancer-types"
```

### Fetch Study Details
```bash
curl "https://www.cbioportal.org/api/studies/gbm_tcga_pub"
```

### Fetch Molecular Profiles for Study
```bash
curl "https://www.cbioportal.org/api/molecular-profiles?studyId=gbm_tcga_pub"
```

### Fetch Clinical Attributes
```bash
curl "https://www.cbioportal.org/api/clinical-attributes"
```

### Fetch Genes
```bash
curl "https://www.cbioportal.org/api/genes"
```

---

## Client Libraries

| Language | Library | Description |
|----------|---------|-------------|
| R | cBioPortalData | Bioconductor integration with MultiAssayExperiment |
| R | cbioportalR | Tidyverse-compatible, clinical researcher-focused |
| Python | bravado | Direct Swagger/OpenAPI client generation |
| Python | cbio_py | Simplified wrapper API |

### Python Example (bravado)

```python
from bravado.client import SwaggerClient

cbioportal = SwaggerClient.from_url(
    'https://www.cbioportal.org/api/v2/api-docs',
    config={"validate_requests": False, "validate_responses": False}
)

# Fetch mutations
mutations = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
    molecularProfileId="gbm_tcga_pub_mutations",
    sampleListId="gbm_tcga_pub_all"
).result()
```

---

## Authentication

- **Public Data**: No authentication required
- **Private/Authenticated Portals**: Bearer token via Authorization header

```bash
curl -H "Authorization: Bearer YOUR_TOKEN" "https://www.cbioportal.org/api/studies"
```

---

## Integration Notes

### Data Harmonization
- All studies use standardized cancer type ontology
- Reference genomes: hg19 (legacy) and hg38 (current)
- Gene symbols follow HUGO nomenclature

### Swagger/OpenAPI
- Full specification available at `/api/v2/api-docs`
- Interactive documentation via Swagger UI

### Rate Limiting
- Generally unrestricted for reasonable use
- Contact for high-volume access

---

## References

- cBioPortal: https://www.cbioportal.org/
- API Documentation: https://docs.cbioportal.org/web-api-and-clients/
- Swagger UI: https://www.cbioportal.org/api/swagger-ui.html
- GitHub: https://github.com/cBioPortal/cbioportal
