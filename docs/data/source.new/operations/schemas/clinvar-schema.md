---
id: schema-clinvar
title: "ClinVar Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database]
---

**Parent:** [Schema Documentation](./_index.md)

# ClinVar Schema Documentation

**Document ID:** SCHEMA-CLINVAR
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** ClinVar XSD v2.5, README.txt, Tab-Delimited Files

---

## Overview

ClinVar is NCBI's public archive of interpretations of clinically relevant variants. This document details the actual XML schemas, file formats, and data structures from ClinVar's official documentation and XSD files.

**FTP Base:** `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/`
**Update Schedule:** Weekly (Mondays), Monthly (First Thursday)

---

## Database Statistics (January 2026)

| Metric | Value | Source |
|--------|-------|--------|
| **Top Submitter** | Ambry Genetics | 1.75M submissions |
| **Second Submitter** | LabCorp Genetics | 1.89M submissions |
| **Third Submitter** | GeneDx | 436K submissions |
| **VCF File Size** | 177 MB | GRCh38 compressed |
| **Submission History** | 2009-2026 | Active growth |

---

## Core Identifier Hierarchy

ClinVar uses a three-tier accession system:

### VCV (Variation Archive)

The aggregate interpretation across all submitters for a specific variant.

| Property | Description |
|----------|-------------|
| **Format** | `VCV000000123.4` (accession.version) |
| **Scope** | One variant, all conditions |
| **Aggregation** | Precedence to records with assertion criteria |

### RCV (Record)

Individual record representing one submitter's interpretation of a variant-condition pair.

| Property | Description |
|----------|-------------|
| **Format** | `RCV000012345.6` (accession.version) |
| **Scope** | One variant + one condition |
| **Content** | Clinical significance, review status, evidence |

### SCV (Submitted Record)

The accession assigned to each submitted interpretation.

| Property | Description |
|----------|-------------|
| **Format** | `SCV000123456.7` (accession.version) |
| **Scope** | Single submission |
| **Content** | Original submitter assertion |

---

## XML Schema Structure (XSD v2.5)

### Root Element

```xml
<ClinVarVariationRelease>
  <VariationArchive>
    <RecordType>classified | included</RecordType>
    <VariationID/>
    <VariationName/>
    <ClassifiedRecord/> <!-- or IncludedRecord -->
  </VariationArchive>
</ClinVarVariationRelease>
```

### ClassifiedRecord Structure

For variants with explicit submitted classifications:

```xml
<ClassifiedRecord>
  <SimpleAllele/>           <!-- or Haplotype or Genotype -->
  <RCVList/>                <!-- Multiple RCV accessions -->
  <Classifications>
    <GermlineClassification/>
    <SomaticClinicalImpact/>
    <OncogenicityClassification/>
  </Classifications>
  <ClinicalAssertionList/>  <!-- Individual assertions -->
  <TraitMappingList/>       <!-- MedGen CUI mappings -->
</ClassifiedRecord>
```

### IncludedRecord Structure

For alleles from complex submissions without independent classification:

```xml
<IncludedRecord>
  <SimpleAllele/>
  <SubmittedClassificationList/>  <!-- Component SCV accessions -->
  <ClassifiedVariationList/>      <!-- Referenced variants -->
</IncludedRecord>
```

---

## Classification Framework

### GermlineClassification

```xml
<GermlineClassification>
  <ReviewStatus>criteria provided, single submitter</ReviewStatus>
  <Description>Pathogenic</Description>
  <ConditionList>
    <TraitSet Type="Disease">
      <Trait Type="Disease">
        <Name>
          <ElementValue Type="Preferred">Sickle cell anemia</ElementValue>
        </Name>
        <XRef DB="MedGen" ID="C0002895"/>
        <XRef DB="OMIM" ID="603903"/>
      </Trait>
    </TraitSet>
  </ConditionList>
  <DateLastEvaluated>2025-06-15</DateLastEvaluated>
  <NumberOfSubmitters>5</NumberOfSubmitters>
  <NumberOfSubmissions>8</NumberOfSubmissions>
</GermlineClassification>
```

### Review Status Values

| Status | Evidence Level |
|--------|----------------|
| `practice guideline` | Highest - expert panel guideline |
| `reviewed by expert panel` | Expert panel review |
| `criteria provided, multiple submitters, no conflicts` | Strong - consensus |
| `criteria provided, single submitter` | Moderate - single source |
| `criteria provided, conflicting interpretations` | Conflicting data |
| `no assertion criteria provided` | Minimal - no criteria |
| `no classification provided` | None |

### Classification Values

| Category | Values |
|----------|--------|
| **Germline** | Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign |
| **Somatic** | Tier I - Strong, Tier II - Potential, Tier III - Unknown, Tier IV - Benign/Likely benign |
| **Oncogenicity** | Oncogenic, Likely oncogenic, Uncertain significance, Likely benign, Benign |

---

## Complex Type Definitions

### typeAllele

Represents single variants:

```xml
<SimpleAllele AlleleID="12345" VariationID="67890">
  <GeneList>
    <Gene Symbol="HBB" GeneID="3043">
      <Relationship>within single gene</Relationship>
    </Gene>
  </GeneList>
  <Name>NM_000518.5(HBB):c.20A>T (p.Glu7Val)</Name>
  <VariantType>single nucleotide variant</VariantType>
  <Location>
    <CytogeneticLocation>11p15.4</CytogeneticLocation>
    <SequenceLocation Assembly="GRCh38"
                      Chr="11"
                      start="5227002"
                      stop="5227002"
                      referenceAllele="T"
                      alternateAllele="A"/>
  </Location>
  <HGVSlist>
    <HGVS Type="genomic, top-level">NC_000011.10:g.5227002T>A</HGVS>
    <HGVS Type="coding">NM_000518.5:c.20A>T</HGVS>
    <HGVS Type="protein">NP_000509.1:p.Glu7Val</HGVS>
  </HGVSlist>
  <AlleleFrequencyList>
    <AlleleFrequency Value="0.0108" Source="gnomAD"/>
  </AlleleFrequencyList>
</SimpleAllele>
```

### typeLocation

Genomic positioning:

| Element | Description | Required |
|---------|-------------|----------|
| `CytogeneticLocation` | Chromosome band (e.g., 11p15.4) | No |
| `SequenceLocation@Assembly` | GRCh37 or GRCh38 | Yes |
| `SequenceLocation@Chr` | 1-22, X, Y, MT | Yes |
| `SequenceLocation@start` | 1-based start position | Yes |
| `SequenceLocation@stop` | 1-based stop position | Yes |
| `SequenceLocation@referenceAllele` | Reference sequence | Yes |
| `SequenceLocation@alternateAllele` | Variant sequence | Yes |
| `SequenceLocation@positionVCF` | VCF-style position (left-aligned) | No |

### typeHGVSExpression

HGVS nomenclature expressions:

```xml
<HGVSlist>
  <HGVS Type="genomic, top-level" Assembly="GRCh38">
    NC_000011.10:g.5227002T>A
  </HGVS>
  <HGVS Type="genomic">NG_000007.3:g.70644T>A</HGVS>
  <HGVS Type="coding">NM_000518.5:c.20A>T</HGVS>
  <HGVS Type="non-coding"/>
  <HGVS Type="protein">NP_000509.1:p.Glu7Val</HGVS>
</HGVSlist>
```

| Type | Description |
|------|-------------|
| `genomic, top-level` | Chromosome-level (NC_) |
| `genomic` | Gene-level (NG_) |
| `coding` | Transcript coding (NM_) |
| `non-coding` | Non-coding transcript (NR_) |
| `protein` | Protein change (NP_) |

### typeSample

Supporting evidence samples:

```xml
<Sample>
  <Origin>germline</Origin>
  <Species TaxonomyId="9606">human</Species>
  <AffectedStatus>yes</AffectedStatus>
  <NumberTested>1</NumberTested>
  <Gender>female</Gender>
  <FamilyData NumFamilies="1">
    <FamilyHistory>yes</FamilyHistory>
  </FamilyData>
  <Indication Type="Indication">
    <Trait Type="Finding">
      <Name>
        <ElementValue Type="Preferred">Sickle cell trait</ElementValue>
      </Name>
    </Trait>
  </Indication>
</Sample>
```

### Origin Values (typeOrigin)

| Value | Description |
|-------|-------------|
| `germline` | Inherited variant |
| `somatic` | Acquired variant |
| `de novo` | New mutation |
| `inherited` | Confirmed inherited |
| `maternal` | Maternal inheritance |
| `paternal` | Paternal inheritance |
| `uniparental` | Uniparental disomy |
| `biparental` | Both parents |
| `experimentally generated` | Lab-created |

---

## Tab-Delimited File Schemas

### variant_summary.txt

Core variant data by genomic location.

| Column | Type | Description |
|--------|------|-------------|
| `AlleleID` | integer | Unique simple variant identifier |
| `Type` | string | Variant type (SNV, deletion, etc.) |
| `Name` | string | Preferred name |
| `GeneID` | integer | NCBI Gene ID |
| `GeneSymbol` | string | HGNC symbol |
| `HGNC_ID` | string | HGNC identifier |
| `ClinicalSignificance` | string | Semicolon-separated classifications |
| `ClinSigSimple` | integer | 0=benign, 1=pathogenic, -1=conflicting |
| `LastEvaluated` | date | Most recent evaluation |
| `RS#(dbSNP)` | integer | dbSNP RS identifier |
| `nsv/esv(dbVar)` | string | dbVar SV identifier |
| `RCVaccession` | string | Semicolon-separated RCV IDs |
| `PhenotypeIDS` | string | MedGen CUI identifiers |
| `PhenotypeList` | string | Condition names |
| `Origin` | string | Variant origin |
| `OriginSimple` | string | Simplified origin |
| `Assembly` | string | GRCh37 or GRCh38 |
| `ChromosomeAccession` | string | RefSeq chromosome |
| `Chromosome` | string | Chromosome number |
| `Start` | integer | Start position (right-shifted) |
| `Stop` | integer | Stop position |
| `ReferenceAllele` | string | Reference sequence |
| `AlternateAllele` | string | Variant sequence |
| `Cytogenetic` | string | Cytogenetic location |
| `ReviewStatus` | string | Aggregate review status |
| `NumberSubmitters` | integer | Submitter count |
| `Guidelines` | string | Associated guidelines |
| `TestedInGTR` | string | GTR test availability |
| `OtherIDs` | string | External identifiers |
| `SubmitterCategories` | integer | Submitter category flags |
| `VariationID` | integer | ClinVar variation ID |
| `PositionVCF` | integer | VCF-style position (left-aligned) |
| `ReferenceAlleleVCF` | string | VCF-style reference |
| `AlternateAlleleVCF` | string | VCF-style alternate |
| `SomaticClinicalImpact` | string | Somatic classification |
| `SomaticClinicalImpactLastEvaluated` | date | Somatic evaluation date |
| `ReviewStatusSomatic` | string | Somatic review status |
| `Oncogenicity` | string | Oncogenicity classification |
| `OncogenicityLastEvaluated` | date | Oncogenicity evaluation date |
| `ReviewStatusOncogenicity` | string | Oncogenicity review status |

### submission_summary.txt

Individual submission details:

| Column | Type | Description |
|--------|------|-------------|
| `VariationID` | integer | ClinVar variation ID |
| `ClinicalSignificance` | string | Submitter's classification |
| `DateLastEvaluated` | date | Evaluation date |
| `Description` | string | Interpretation description |
| `SubmittedPhenotypeInfo` | string | Submitted condition |
| `ReportedPhenotypeInfo` | string | Mapped MedGen terms |
| `ReviewStatus` | string | Submission review level |
| `CollectionMethod` | string | How evidence was collected |
| `OriginCounts` | string | Origin breakdown |
| `Submitter` | string | Organization name |
| `SCV` | string | SCV accession |
| `SubmittedGeneSymbol` | string | Submitted gene |
| `ExplanationOfInterpretation` | string | Classification rationale |
| `SomaticClinicalImpact` | string | Somatic classification |
| `Oncogenicity` | string | Oncogenicity classification |

### gene_specific_summary.txt

Gene-level aggregations:

| Column | Type | Description |
|--------|------|-------------|
| `Symbol` | string | Gene symbol |
| `GeneID` | integer | NCBI Gene ID |
| `Total_submissions` | integer | All submissions |
| `Total_alleles` | integer | Unique alleles |
| `Submissions_reporting_this_gene` | integer | Gene-specific submissions |
| `Alleles_reported_Pathogenic_Likely_pathogenic` | integer | P/LP allele count |
| `Gene_MIM_number` | integer | OMIM gene ID |
| `Number_Uncertain` | integer | VUS count |
| `Number_with_conflicts` | integer | Conflicting interpretations |

### hgvs4variation.txt

HGVS nomenclature by variation:

| Column | Type | Description |
|--------|------|-------------|
| `Symbol` | string | Gene symbol |
| `GeneID` | integer | NCBI Gene ID |
| `VariationID` | integer | ClinVar variation ID |
| `AlleleID` | integer | Allele identifier |
| `Type` | string | coding, genomic, non-coding, protein |
| `Assembly` | string | Reference assembly |
| `NucleotideExpression` | string | Nucleotide HGVS |
| `NucleotideChange` | string | Change portion only |
| `ProteinExpression` | string | Protein HGVS |
| `ProteinChange` | string | Protein change only |
| `UsedForNaming` | string | yes/no |
| `Submitted` | string | Submitter-provided HGVS |
| `OnRefSeqGene` | string | yes/no |

### allele_gene.txt

Allele-to-gene relationships:

| Column | Type | Description |
|--------|------|-------------|
| `AlleleID` | integer | Allele identifier |
| `GeneID` | integer | NCBI Gene ID |
| `Symbol` | string | Gene symbol |
| `Name` | string | Gene name |
| `GeneRelationship` | string | Relationship type |
| `Source` | string | Submitted or calculated |

**Relationship Types:**
- `asserted, but not computed`
- `genes overlapped by variant`
- `near gene, upstream`
- `near gene, downstream`
- `within multiple genes by overlap`
- `within single gene`

### cross_references.txt

External database mappings:

| Column | Type | Description |
|--------|------|-------------|
| `AlleleID` | integer | ClinVar allele |
| `Database` | string | External database name |
| `ID` | string | External identifier |
| `last_updated` | date | Mapping update date |

**Supported Databases:** dbSNP, dbVar

---

## VCF File Format

### File Locations

| File | Size | Description |
|------|------|-------------|
| `clinvar.vcf.gz` | 177 MB | Current release (GRCh38) |
| `clinvar_papu.vcf.gz` | 68 KB | Pathogenic/Likely pathogenic subset |

### VCF INFO Fields

| Field | Type | Description |
|-------|------|-------------|
| `ALLELEID` | Integer | ClinVar allele ID |
| `CLNDN` | String | Disease name(s) |
| `CLNDISDB` | String | Disease database references |
| `CLNHGVS` | String | HGVS genomic expression |
| `CLNREVSTAT` | String | Review status |
| `CLNSIG` | String | Clinical significance |
| `CLNSIGCONF` | String | Conflicting interpretations |
| `CLNSIGINCL` | String | Included variant significance |
| `CLNVC` | String | Variant type |
| `CLNVCSO` | String | Sequence Ontology ID |
| `CLNVI` | String | Variation identifiers |
| `DBVARID` | String | dbVar ID |
| `GENEINFO` | String | Gene:ID pairs |
| `MC` | String | Molecular consequence |
| `ORIGIN` | Integer | Origin code |
| `RS` | Integer | dbSNP RS ID |

### CLNSIG Values

| Code | Meaning |
|------|---------|
| `Benign` | Benign |
| `Likely_benign` | Likely benign |
| `Uncertain_significance` | VUS |
| `Likely_pathogenic` | Likely pathogenic |
| `Pathogenic` | Pathogenic |
| `drug_response` | Drug response variant |
| `Conflicting_classifications_of_pathogenicity` | Conflicting |
| `not_provided` | Not provided |

### ORIGIN Codes

| Code | Meaning |
|------|---------|
| 0 | Unknown |
| 1 | Germline |
| 2 | Somatic |
| 4 | Inherited |
| 8 | Paternal |
| 16 | Maternal |
| 32 | De novo |
| 64 | Biparental |
| 128 | Uniparental |
| 256 | Not tested |
| 512 | Tested inconclusive |
| 1073741824 | Other |

---

## Phenotype Mapping Files

### disease_names

Condition name standardization:

| Column | Type | Description |
|--------|------|-------------|
| `DiseaseName` | string | Preferred name |
| `SourceName` | string | Source database |
| `ConceptID` | string | UMLS C-number or NCBI CN-number |
| `SourceID` | string | Source identifier |
| `DiseaseMIM` | integer | OMIM disease ID |
| `LastUpdated` | date | Update timestamp |
| `Category` | string | Classification category |

**Categories:** Blood group, Disease, Finding, Named protein variant, Pharmacological response, Phenotype instruction

### gene_condition_source_id

Gene-disease associations:

| Column | Type | Description |
|--------|------|-------------|
| `GeneID` | integer | NCBI Gene ID |
| `GeneSymbol` | string | HGNC symbol |
| `ConceptID` | string | MedGen concept |
| `DiseaseName` | string | Condition name |
| `SourceName` | string | Source database |
| `SourceID` | string | External identifier |
| `AssociatedGenes` | string | Related genes |
| `RelatedGenes` | string | Additional genes |

---

## Data File Sizes

| File | Compressed | Description |
|------|------------|-------------|
| `variant_summary.txt.gz` | 421 MB | Core variant data |
| `submission_summary.txt.gz` | 359 MB | Submission details |
| `hgvs4variation.txt.gz` | 484 MB | HGVS nomenclature |
| `summary_of_conflicting_interpretations.txt` | 1.6 GB | Conflicts (uncompressed) |
| `var_citations.txt` | 215 MB | Citations |
| `cross_references.txt` | 112 MB | External references |
| `allele_gene.txt.gz` | 87 MB | Gene mappings |
| `variation_allele.txt.gz` | 23 MB | Variation-allele links |
| `gene_specific_summary.txt` | 3.5 MB | Gene summaries |
| `organization_summary.txt` | 864 KB | Submitter info |

---

## Update Schedule

| Type | Frequency | Day | Archive |
|------|-----------|-----|---------|
| Tab-delimited | Weekly | Monday | Monthly |
| VCF | Weekly | Monday | Monthly |
| XML (weekly) | Weekly | Monday | No |
| XML (monthly) | Monthly | First Thursday | Yes |

---

## Integration Examples

### Python: Parse variant_summary.txt

```python
import gzip
import csv

def parse_variant_summary(filepath: str):
    """Parse ClinVar variant_summary.txt.gz file."""
    with gzip.open(filepath, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            yield {
                'allele_id': int(row['AlleleID']),
                'gene_symbol': row['GeneSymbol'],
                'clinical_significance': row['ClinicalSignificance'],
                'review_status': row['ReviewStatus'],
                'chromosome': row['Chromosome'],
                'start': int(row['Start']) if row['Start'] else None,
                'stop': int(row['Stop']) if row['Stop'] else None,
                'hgvs_c': None,  # From hgvs4variation
                'rsid': int(row['RS#(dbSNP)']) if row['RS#(dbSNP)'] != '-1' else None
            }
```

### Python: Parse VCF

```python
import gzip

def parse_clinvar_vcf(filepath: str):
    """Parse ClinVar VCF file."""
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

            yield {
                'chrom': chrom,
                'pos': int(pos),
                'id': id_,
                'ref': ref,
                'alt': alt,
                'allele_id': int(info_dict.get('ALLELEID', 0)),
                'clnsig': info_dict.get('CLNSIG', ''),
                'clnrevstat': info_dict.get('CLNREVSTAT', ''),
                'gene_info': info_dict.get('GENEINFO', '')
            }
```

---

## Download

### FTP Access (Primary)

ClinVar files are distributed via FTP and updated weekly:

**FTP Base:** https://ftp.ncbi.nlm.nih.gov/pub/clinvar/

```bash
# Download weekly VCF (current/GRCh38)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# Download GRCh37 (hg19) version
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

# Download all tab-delimited files
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/hgvs4variation.txt.gz
```

### File Organization

```
/pub/clinvar/
├── vcf_GRCh38/              # Current (hg38)
│   ├── clinvar.vcf.gz
│   ├── clinvar.vcf.gz.tbi
│   ├── clinvar_papu.vcf.gz  # P/LP subset
│   └── README.txt
├── vcf_GRCh37/              # Legacy (hg19)
│   ├── clinvar.vcf.gz
│   └── clinvar.vcf.gz.tbi
├── tab_delimited/           # TSV files
│   ├── variant_summary.txt.gz
│   ├── submission_summary.txt.gz
│   ├── hgvs4variation.txt.gz
│   ├── gene_specific_summary.txt.gz
│   ├── gene_condition_source_id
│   ├── disease_names
│   ├── allele_gene.txt.gz
│   ├── variation_allele.txt.gz
│   └── cross_references.txt.gz
├── xml/                     # XML files (monthly archive)
│   ├── clinvar_YYYYMM.xml.gz
│   └── ...
└── archive/                 # Historical releases by date
```

### Update Schedule

| Type | Frequency | Day | Archive |
|------|-----------|-----|---------|
| **VCF** | Weekly | Monday | Monthly |
| **Tab-delimited** | Weekly | Monday | Monthly |
| **XML (full)** | Monthly | First Thursday | Archived |

### File Sizes (as of January 2026)

| File | Compressed | Uncompressed | Description |
|------|-----------|--------------|-------------|
| clinvar.vcf.gz | 177 MB | ~750 MB | Current GRCh38 VCF |
| clinvar_papu.vcf.gz | 68 KB | ~250 KB | Pathogenic/Likely pathogenic subset |
| variant_summary.txt.gz | 421 MB | ~1.2 GB | All variants (core data) |
| submission_summary.txt.gz | 359 MB | ~980 MB | Individual submissions |
| hgvs4variation.txt.gz | 484 MB | ~1.4 GB | HGVS nomenclature |

**Total Dataset:** ~1.5 TB (compressed), ~5.6 TB (uncompressed)

### Alternative Access

**NCBI Entrez Direct:**
```bash
# Search and download ClinVar records
esearch -db clinvar -query "TP53" | efetch -format xml > tp53_variants.xml
```

**REST API:**
```bash
# Query by variation ID
curl "https://www.ncbi.nlm.nih.gov/clinvar/?term=VCV000000001"

# Query by gene
curl "https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA1[gene]"
```

---

## Data Format

| Format | Description |
|--------|-------------|
| Primary | VCF (.vcf.gz), XML |
| Alternative | TSV (tab-delimited text files) |
| Compression | gzip (.gz) |
| Encoding | UTF-8 |
| API Response | XML, JSON |

### Format Details

**VCF (Variant Call Format)**
- Text-based variant representation
- Standard #CHROM POS ID REF ALT columns
- INFO field with ClinVar-specific annotations
- Indexed with tabix (.tbi) for fast lookup
- Sorted by chromosome and position

**XML (eXtensible Markup Language)**
- Complete ClinVar submission records
- Hierarchical variant/interpretation structure
- Full classification history
- Evidence details included
- Complies with XSD v2.5

**TSV (Tab-Separated Values)**
- Flat file format
- One variant or submission per row
- Header row with column names
- Pipe-delimited for multi-value fields
- Text encoding UTF-8

---

## Data Set Size

| Metric | Value | Notes |
|--------|-------|-------|
| **Total Submissions** | 10,000,000+ | From hundreds of organizations |
| **Unique Variants** | 1,500,000+ | Distinct alleles |
| **Top Submitter** | LabCorp Genetics | 1.89M submissions |
| **Variant Types** | SNV, indel, SV, complex | All variation types |
| **Assembly Coverage** | GRCh38, GRCh37 | Both major builds |
| **Active Updates** | Weekly | Monday releases |
| **Archive Size** | ~6 TB (uncompressed) | Historical + current |
| **Submission Span** | 2009-2026 | 17 years of data |
| **Review Status Levels** | 7 categories | From practice guideline to unreviewed |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "VCV000000123.4" |
| `name` | string | Entity name | "Variant Archive" |
| `type` | string | Record type | "variation" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `associated_with` | Entity | N:M |

---

## License

| Resource | License | Commercial Use |
|----------|---------|----------------|
| ClinVar | CC BY (public domain) | Yes |

---

## Sample Data

### Example Record
```json
{
  "vcv": "VCV000000123.4",
  "rcv": "RCV000012345.6",
  "allele_id": 12345,
  "gene_symbol": "HBB",
  "clinical_significance": "Pathogenic",
  "review_status": "criteria provided, multiple submitters, no conflicts",
  "hgvs_c": "NM_000518.5:c.20A>T",
  "chromosome": "11"
}
```

### Sample Query Result
| vcv | rcv | gene_symbol | clinical_significance | review_status | hgvs_c |
|-----|-----|-------------|----------------------|----------------|--------|
| VCV000000123.4 | RCV000012345.6 | HBB | Pathogenic | criteria provided, multiple submitters, no conflicts | NM_000518.5:c.20A>T |
| VCV000000124.5 | RCV000012346.7 | BRCA1 | Likely pathogenic | reviewed by expert panel | NM_007294.4:c.68_69delAG |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `VCV` | Variation Archive accession representing aggregate interpretation across all submitters for a variant | VCV000000123.4 |
| `RCV` | Record accession for a single variant-condition pair with clinical interpretation | RCV000012345.6 |
| `SCV` | Submitted record accession for an individual submission from a single submitter | SCV000123456.7 |
| `AlleleID` | Unique identifier for a simple variant (single allele) in ClinVar | 12345 |
| `VariationID` | Identifier for any variation record including complex variants (haplotypes, genotypes) | 67890 |
| `CLNSIG` | Clinical significance classification assigned to a variant | Pathogenic |
| `CLNREVSTAT` | Review status indicating level of evidence and submitter consensus | criteria provided, multiple submitters |
| `HGVS` | Human Genome Variation Society nomenclature for describing sequence variants | NM_000518.5:c.20A>T |
| `GeneID` | NCBI Entrez Gene identifier linking variant to affected gene | 3043 (HBB gene) |
| `MedGen CUI` | Concept Unique Identifier from MedGen for disease/phenotype standardization | C0002895 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Pathogenic | Variant causes disease with strong evidence | CLNSIG classification |
| Likely Pathogenic | Variant probably causes disease (>90% certainty) | CLNSIG classification |
| VUS (Uncertain Significance) | Insufficient evidence to classify variant as benign or pathogenic | CLNSIG classification |
| Likely Benign | Variant probably does not cause disease (>90% certainty) | CLNSIG classification |
| Benign | Variant does not cause disease with strong evidence | CLNSIG classification |
| Germline Variant | Inherited genetic change present in egg/sperm cells | Origin field |
| Somatic Variant | Acquired genetic change occurring after conception | Origin field |
| De Novo Variant | New mutation not inherited from either parent | Origin field |
| Review Status | Level of evidence and agreement among submitters | CLNREVSTAT field |
| Practice Guideline | Highest review level with expert panel consensus | Review status value |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| ClinVar | Clinical Variation database | NCBI public archive |
| VCV | Variation Archive Accession | Aggregate variant record |
| RCV | Record Accession | Variant-condition pair |
| SCV | Submitted Record Accession | Individual submission |
| HGVS | Human Genome Variation Society | Variant nomenclature standard |
| VUS | Variant of Uncertain Significance | Classification category |
| P/LP | Pathogenic/Likely Pathogenic | Combined classification |
| B/LB | Benign/Likely Benign | Combined classification |
| GRCh37 | Genome Reference Consortium Human Build 37 | hg19 equivalent |
| GRCh38 | Genome Reference Consortium Human Build 38 | hg38 equivalent |
| VCF | Variant Call Format | Standard variant file format |
| OMIM | Online Mendelian Inheritance in Man | Disease database |
| MedGen | Medical Genetics database | NCBI phenotype resource |
| GTR | Genetic Testing Registry | NCBI test database |
| ACMG | American College of Medical Genetics | Classification guidelines |
| AMP | Association for Molecular Pathology | Classification guidelines |
| SNV | Single Nucleotide Variant | Point mutation |
| CUI | Concept Unique Identifier | UMLS/MedGen identifier |

---

## References

1. ClinVar FTP: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
2. ClinVar XSD: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/
3. ClinVar Documentation: https://www.ncbi.nlm.nih.gov/clinvar/docs/
4. Landrum MJ, et al. (2018). ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res. 46(D1):D1062-D1067.

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema from XSD v2.5 and README |
