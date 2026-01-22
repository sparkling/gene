---
id: schema-dbvar
title: "dbVar Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, database, structural-variation, ncbi]
---

**Parent:** [Schema Documentation](./_index.md)

# dbVar Schema Documentation

**Document ID:** SCHEMA-DBVAR
**Status:** Draft
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source:** NCBI dbVar Documentation, GitHub Repository, VCF Submission Guidelines

---

## TL;DR

dbVar is NCBI's public database of human genomic structural variation, archiving variants greater than 50 base pairs including deletions, duplications, insertions, inversions, and complex rearrangements. The database contains over 6 million structural variants from 185+ studies, available in TSV, BED, BEDPE, VCF, and GVF formats.

**Key Access Points:**
- FTP: `https://ftp.ncbi.nlm.nih.gov/pub/dbVar/`
- Browser: `https://www.ncbi.nlm.nih.gov/dbvar`
- GitHub: `https://github.com/ncbi/dbvar`

---

## License & Access

| Attribute | Value |
|-----------|-------|
| **License** | Public Domain (US Government Work) |
| **Attribution** | Recommended but not required |
| **Commercial Use** | Permitted |
| **Data Sharing** | Unrestricted |
| **Terms** | NCBI Data Usage Policy |

**Note:** dbVar data is in the public domain. Clinical structural variants should be submitted to ClinVar, which automatically forwards structural variants to dbVar.

---

## Database Statistics (January 2026)

| Metric | Value | Source |
|--------|-------|--------|
| **Total Structural Variants** | 6M+ submitted | dbVar statistics |
| **Deletions** | 1.9M | Non-redundant set |
| **Duplications** | 659K | Non-redundant set |
| **Insertions** | 1.7M | Non-redundant set |
| **Studies** | 185+ | Human studies |
| **Reference Assemblies** | GRCh38, GRCh37, NCBI36 | Multiple builds |
| **Minimum SV Size** | 50 bp | Definition threshold |

---

## Core Data Objects

### Studies (std)

Collections of variant regions and calls representing coherent methodological analyses.

| Prefix | Source | Description |
|--------|--------|-------------|
| `nstd` | NCBI | NCBI-assigned study ID |
| `estd` | EBI | EBI-assigned study ID |
| `dstd` | DDBJ | DDBJ-assigned study ID |

### Variant Regions (sv)

Genomic locations marking regions containing observed structural variation. These represent the submitter's assertions about variant locations.

| Prefix | Source | Description |
|--------|--------|-------------|
| `nsv` | NCBI | NCBI-assigned variant region ID |
| `esv` | EBI | EBI-assigned variant region ID |
| `dsv` | DDBJ | DDBJ-assigned variant region ID |

### Variant Calls (ssv)

Individual structural variant instances from experimental detection.

| Prefix | Source | Description |
|--------|--------|-------------|
| `nssv` | NCBI | NCBI-assigned variant call ID |
| `essv` | EBI | EBI-assigned variant call ID |
| `dssv` | DDBJ | DDBJ-assigned variant call ID |

---

## Structural Variant Types

### Major Categories

| Type | SO ID | Description |
|------|-------|-------------|
| **deletion** | SO:0000159 | Loss of sequence |
| **duplication** | SO:1000035 | Gain of sequence |
| **copy_number_gain** | SO:0001742 | Increased copy number |
| **copy_number_loss** | SO:0001743 | Decreased copy number |
| **copy_number_variation** | SO:0001019 | Variable copy number |
| **insertion** | SO:0000667 | Sequence addition |
| **inversion** | SO:1000036 | Sequence reversal |
| **translocation** | SO:0000199 | Chromosomal rearrangement |
| **complex** | SO:0001784 | Multiple changes |

### Insertion Subtypes

| Type | SO ID | Description |
|------|-------|-------------|
| **mobile_element_insertion** | SO:0001837 | Transposable element |
| **Alu_insertion** | SO:0002063 | Alu element insertion |
| **LINE1_insertion** | SO:0002064 | LINE-1 element insertion |
| **SVA_insertion** | SO:0002066 | SVA element insertion |
| **HERV_insertion** | SO:0002186 | HERV element insertion |
| **novel_sequence_insertion** | SO:0001838 | New sequence |
| **tandem_duplication** | SO:1000173 | Adjacent duplication |

### Complex Types

| Type | Description |
|------|-------------|
| **delins** | Deletion-insertion (indel) |
| **interchromosomal_translocation** | Between chromosomes |
| **intrachromosomal_translocation** | Within chromosome |

---

## Coordinate Systems

dbVar uses multiple coordinate representations:

### Standard Coordinates

| Field | Base | Description |
|-------|------|-------------|
| `start` | 1-based | Genomic start position |
| `stop` | 1-based | Genomic stop position |
| `length` | N/A | Variant size in bp |

### Uncertainty Coordinates

For imprecise breakpoints:

| Field | Description |
|-------|-------------|
| `inner_start` | Minimum affected start (breakpoint outside) |
| `inner_stop` | Minimum affected stop (breakpoint outside) |
| `outer_start` | Maximum boundary start (breakpoint inside) |
| `outer_stop` | Maximum boundary stop (breakpoint inside) |

### BED Format Coordinates

| Position | Base | Description |
|----------|------|-------------|
| `chrom` | N/A | Chromosome with "chr" prefix |
| `chromStart` | 0-based | Start position |
| `chromEnd` | 1-based | End position |

**Note:** BED files use 0-based starts and 1-based stops. Standard dbVar files use 1-based coordinates.

---

## File Formats

### TSV Format (Primary)

Tab-separated values with variant details.

#### Deletion/Duplication TSV Columns

| Col | Field | Type | Description |
|-----|-------|------|-------------|
| 1 | `chr` | String | Chromosome (e.g., chr1) |
| 2 | `outermost_start` | Integer | 1-based start position |
| 3 | `outermost_stop` | Integer | 1-based stop position |
| 4 | `variant_count` | Integer | Number of supporting variants |
| 5 | `variant_type` | String | SV classification |
| 6 | `method` | String | Detection method(s) |
| 7 | `analysis` | String | Analysis type(s) |
| 8 | `platform` | String | Sequencing/array platform(s) |
| 9 | `study` | String | Source study accession(s) |
| 10 | `variant` | String | dbVar call accession(s) |
| 11 | `clinical_assertion` | String | Clinical significance |
| 12 | `clinvar_accession` | String | ClinVar accession(s) |
| 13 | `bin_size` | String | Size category |

#### Insertion TSV Additional Columns

| Col | Field | Type | Description |
|-----|-------|------|-------------|
| 14 | `min_insertion_length` | Integer | Minimum insert size |
| 15 | `max_insertion_length` | Integer | Maximum insert size |

#### ACMG-Annotated TSV Additional Columns

| Col | Field | Type | Description |
|-----|-------|------|-------------|
| 14+ | `gene` | String | Overlapping gene symbol(s) |

### BED Format

Standard BED format for genome browser viewing.

| Col | Field | Type | Description |
|-----|-------|------|-------------|
| 1 | `chrom` | String | Chromosome (chr prefix) |
| 2 | `chromStart` | Integer | 0-based start |
| 3 | `chromEnd` | Integer | 1-based end |
| 4 | `name` | String | NR_SV_id (chr_start_stop_type) |

### BEDPE Format

Paired-end format for complex rearrangements.

| Col | Field | Type | Description |
|-----|-------|------|-------------|
| 1 | `chrom1` | String | First chromosome |
| 2 | `start1` | Integer | First start (0-based) |
| 3 | `end1` | Integer | First end (1-based) |
| 4 | `chrom2` | String | Second chromosome ("." if N/A) |
| 5 | `start2` | Integer | Second start ("." if N/A) |
| 6 | `end2` | Integer | Second end ("." if N/A) |
| 7 | `name` | String | Variant identifier |
| 8 | `score` | String | Score ("." placeholder) |
| 9 | `strand1` | String | First strand |
| 10 | `strand2` | String | Second strand |
| 11+ | Optional | String | Additional metadata |

### VCF Format

Standard Variant Call Format for structural variants.

#### VCF Fixed Columns

| Col | Field | Type | Description |
|-----|-------|------|-------------|
| 1 | `CHROM` | String | Chromosome identifier |
| 2 | `POS` | Integer | 1-based position |
| 3 | `ID` | String | Variant identifier |
| 4 | `REF` | String | Reference allele (first base) |
| 5 | `ALT` | String | Alternate allele or SV notation |
| 6 | `QUAL` | Float | Quality score |
| 7 | `FILTER` | String | Filter status |
| 8 | `INFO` | String | Semicolon-separated annotations |

#### VCF INFO Fields

| Field | Type | Description |
|-------|------|-------------|
| `SVTYPE` | String | Structural variant type |
| `END` | Integer | End position |
| `SVLEN` | Integer | Length of SV (negative for deletions) |
| `CIPOS` | Integer,Integer | Confidence interval around POS |
| `CIEND` | Integer,Integer | Confidence interval around END |
| `IMPRECISE` | Flag | Imprecise structural variant |
| `EXPERIMENT` | Integer | Experiment identifier |
| `SAMPLE` | String | Sample identifier |
| `REGION` | String | dbVar region accession |
| `CALLID` | String | dbVar call accession |

### GVF Format

Genome Variation Format (sequence ontology-based).

| Col | Field | Description |
|-----|-------|-------------|
| 1 | `seqid` | Sequence identifier |
| 2 | `source` | Data source |
| 3 | `type` | SO term for variant type |
| 4 | `start` | Start position |
| 5 | `end` | End position |
| 6 | `score` | Score |
| 7 | `strand` | Strand |
| 8 | `phase` | Phase |
| 9 | `attributes` | Key=value pairs |

---

## Clinical Significance Values

| Value | Description |
|-------|-------------|
| `Pathogenic` | Disease-causing |
| `Likely pathogenic` | Probable disease-causing |
| `Uncertain significance` | Unknown clinical impact |
| `Likely benign` | Probably not disease-causing |
| `Benign` | Not disease-causing |
| `not provided` | No assertion provided |

---

## Detection Methods

| Method | Description |
|--------|-------------|
| `Sequencing` | High-throughput sequencing |
| `SNP genotyping analysis` | SNP array analysis |
| `Oligo aCGH` | Oligonucleotide array CGH |
| `BAC aCGH` | BAC array CGH |
| `Optical mapping` | Optical genome mapping |
| `PCR` | Polymerase chain reaction |
| `FISH` | Fluorescence in situ hybridization |
| `Karyotyping` | Chromosome analysis |
| `Multiple` | Multiple methods |

---

## Bin Size Categories

| Category | Size Range |
|----------|------------|
| `small` | < 50 bp |
| `medium` | 50 bp - 1 Mb |
| `large` | >= 1 Mb |

---

## Non-Redundant (NR) Files

### File Naming Convention

```
GRCh{38|37}.nr_<type>.{tsv|bed|bedpe}.gz
```

### Available NR Files

| Type | Description |
|------|-------------|
| `deletions` | Non-redundant deletions |
| `duplications` | Non-redundant duplications |
| `insertions` | Non-redundant insertions |

### Non-Redundancy Definition

Variants are grouped by exact coordinate overlap:
- Same chromosome
- Same outermost_start
- Same outermost_stop

Partial overlaps are NOT merged.

---

## FTP Directory Structure

```
ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/
├── data/
│   └── Homo_sapiens/
│       ├── by_assembly/
│       │   ├── GRCh38/
│       │   │   ├── vcf/
│       │   │   ├── gvf/
│       │   │   └── tab/
│       │   └── GRCh37/
│       │       └── ...
│       └── by_study/
│           └── nstd*/
├── Structural_Variant_Sets/
│   └── Nonredundant_Structural_Variants/
│       ├── GRCh38/
│       │   ├── GRCh38.nr_deletions.tsv.gz
│       │   ├── GRCh38.nr_deletions.bed.gz
│       │   ├── GRCh38.nr_duplications.tsv.gz
│       │   ├── GRCh38.nr_duplications.bed.gz
│       │   ├── GRCh38.nr_insertions.tsv.gz
│       │   └── GRCh38.nr_insertions.bed.gz
│       └── GRCh37/
│           └── ...
└── sandbox/
```

---

## Sample Data

### TSV Record (Deletion)

```tsv
#chr	outermost_start	outermost_stop	variant_count	variant_type	method	analysis	platform	study	variant	clinical_assertion	clinvar_accession	bin_size
chr1	10001	20000	5	deletion	Sequencing	Sequence alignment	Illumina HiSeq	nstd166	nssv14580340;nssv14580341	Pathogenic	VCV000012345	medium
chr1	30000	35000	2	deletion	Oligo aCGH	Probe signal intensity	Agilent 244K	nstd102	nssv1234567	Benign	.	medium
```

### TSV Record (Insertion)

```tsv
#chr	outermost_start	outermost_stop	variant_count	variant_type	method	analysis	platform	study	variant	clinical_assertion	clinvar_accession	bin_size	min_insertion_length	max_insertion_length
chr2	50000	50001	3	insertion	Sequencing	Sequence alignment	PacBio	nstd186	nssv16789012	Uncertain significance	.	small	150	350
```

### BED Record

```bed
chr1	10000	20000	chr1_10001_20000_deletion
chr1	29999	35000	chr1_30000_35000_deletion
chr2	49999	50001	chr2_50000_50001_insertion
```

### VCF Record

```vcf
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=REGION,Number=1,Type=String,Description="dbVar region accession">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10001	nsv1234567	N	<DEL>	.	PASS	SVTYPE=DEL;END=20000;SVLEN=-9999;CIPOS=-50,50;CIEND=-50,50;REGION=nsv1234567
chr2	50000	nsv2345678	N	<INS>	.	PASS	SVTYPE=INS;END=50001;SVLEN=250;CIPOS=-10,10;REGION=nsv2345678
```

### GVF Record

```gvf
##gvf-version 1.10
##genome-build NCBI GRCh38
chr1	dbVar	copy_number_loss	10001	20000	.	.	.	ID=nsv1234567;Variant_seq=.;Reference_seq=.;Start_range=9951,10051;End_range=19950,20050
```

---

## Integration Examples

### Python: Parse NR TSV

```python
import gzip
import csv

def parse_dbvar_nr_tsv(filepath: str):
    """Parse dbVar non-redundant TSV file."""
    with gzip.open(filepath, 'rt') as f:
        # Skip comment lines
        for line in f:
            if not line.startswith('#'):
                break

        # Get header
        header = line.strip().split('\t')

        reader = csv.DictReader(f, fieldnames=header, delimiter='\t')
        for row in reader:
            yield {
                'chr': row['chr'],
                'start': int(row['outermost_start']),
                'stop': int(row['outermost_stop']),
                'variant_count': int(row['variant_count']),
                'variant_type': row['variant_type'],
                'methods': row['method'].split(';'),
                'studies': row['study'].split(';'),
                'clinical': row['clinical_assertion'],
                'clinvar': row['clinvar_accession'] if row['clinvar_accession'] != '.' else None,
                'size': int(row['outermost_stop']) - int(row['outermost_start']) + 1
            }

# Example usage
for sv in parse_dbvar_nr_tsv('GRCh38.nr_deletions.tsv.gz'):
    if sv['clinical'] == 'Pathogenic':
        print(f"{sv['chr']}:{sv['start']}-{sv['stop']} ({sv['size']} bp)")
```

### Python: Parse VCF

```python
import gzip

def parse_dbvar_vcf(filepath: str):
    """Parse dbVar VCF file."""
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
                'svtype': info_dict.get('SVTYPE', ''),
                'end': int(info_dict.get('END', pos)),
                'svlen': int(info_dict.get('SVLEN', 0)) if info_dict.get('SVLEN') else None,
                'region': info_dict.get('REGION', ''),
                'imprecise': 'IMPRECISE' in info_dict
            }
```

### Python: Convert BED to 1-Based

```python
def bed_to_1based(bed_start: int, bed_end: int) -> tuple:
    """Convert BED 0-based coordinates to 1-based."""
    return (bed_start + 1, bed_end)

# BED: chr1  10000  20000  (0-based start)
# 1-based: chr1:10001-20000
start_1based, end_1based = bed_to_1based(10000, 20000)
```

### Bedtools: Intersect with Genes

```bash
# Add track header and remove chrMT
echo 'track name="dbVar NR deletions" description="non-redundant deletions from dbVar"' > nr_deletions_track.bed
zcat GRCh38.nr_deletions.bed.gz | grep -v ^chrMT >> nr_deletions_track.bed

# Intersect with gene annotations
bedtools intersect -a nr_deletions_track.bed -b genes.bed -wa -wb > deletions_in_genes.bed
```

---

## Data Exchange

dbVar synchronizes with the European Variation Archive (EVA/DGVa):
- Monthly data exchange
- Shared ID prefixes (nsv/esv/dsv)
- Harmonized variant types

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `nsv` | NCBI-assigned variant region identifier | `nsv1234567` |
| `nssv` | NCBI-assigned variant call identifier for individual observations | `nssv14580340` |
| `nstd` | NCBI-assigned study identifier | `nstd166` |
| `outermost_start` | 1-based genomic start position of variant region | `10001` |
| `outermost_stop` | 1-based genomic stop position of variant region | `20000` |
| `SVTYPE` | Structural variant type in VCF format | `DEL` |
| `SVLEN` | Length of structural variant (negative for deletions) | `-9999` |
| `CIPOS` | Confidence interval around variant start position | `-50,50` |
| `CIEND` | Confidence interval around variant end position | `-50,50` |
| `clinical_assertion` | Clinical significance classification | `Pathogenic` |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Structural Variant | Genomic alteration greater than 50 base pairs | variant_type |
| Copy Number Variation | Gain or loss of genomic segments | copy_number_gain, copy_number_loss |
| Mobile Element Insertion | Insertion of transposable element sequence | Alu_insertion, LINE1_insertion |
| Non-Redundant Set | Deduplicated variants with exact coordinate matches | NR files |
| Breakpoint | Precise location where structural change occurs | CIPOS, CIEND |
| Imprecise Variant | Structural variant with uncertain exact boundaries | IMPRECISE flag |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| dbVar | Database of Genomic Structural Variation | NCBI SV archive |
| SV | Structural Variant | Variants >50bp |
| VCF | Variant Call Format | Standard variant file format |
| BED | Browser Extensible Data | Genomic interval format |
| BEDPE | BED Paired-End | Format for paired breakpoints |
| GVF | Genome Variation Format | SO-based variant format |
| SO | Sequence Ontology | Standardized variant terminology |
| NR | Non-Redundant | Deduplicated variant set |
| EVA | European Variation Archive | EBI structural variant database |
| aCGH | Array Comparative Genomic Hybridization | CNV detection method |
| FISH | Fluorescence In Situ Hybridization | Cytogenetic detection method |

---

## References

1. dbVar Browser: https://www.ncbi.nlm.nih.gov/dbvar
2. dbVar FTP: https://ftp.ncbi.nlm.nih.gov/pub/dbVar/
3. dbVar GitHub: https://github.com/ncbi/dbvar
4. dbVar Help: https://www.ncbi.nlm.nih.gov/dbvar/content/help/
5. Lappalainen I, et al. (2013). dbVar and DGVa: public archives for genomic structural variation. Nucleic Acids Res. 41(D1):D903-D908.

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema from dbVar documentation |
