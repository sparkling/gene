---
id: schema-refseq
title: "NCBI RefSeq Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: final
tags: [schema, database, sequences, ncbi, reference, proteins, transcripts]
---

# NCBI RefSeq Schema Documentation

**Document ID:** SCHEMA-REFSEQ
**Version:** 2026.01
**Source Version:** Release 224

---

## TL;DR

NCBI RefSeq provides non-redundant, well-annotated reference sequences for genomes, transcripts, and proteins across 120,000+ organisms. Records are distinguished by accession prefixes (NP_, NM_, NC_, etc.) indicating molecule type and curation status. Data links directly to Gene, taxonomy, and other NCBI resources.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Protein records | 400,000,000+ | RefSeq Statistics |
| Transcript records | 200,000,000+ | RefSeq Statistics |
| Organisms | 120,000+ | RefSeq Statistics |
| Human protein records | 110,000+ | RefSeq Homo sapiens |
| Complete genomes | 100,000+ | RefSeq Genomes |

---

## Entity Relationship Overview

```
RefSeq Genome (NC_/NW_/NT_)
  └── Gene (GeneID)
        ├── Transcript (NM_/NR_/XM_/XR_)
        │     └── CDS → Protein (NP_/XP_/YP_/WP_)
        └── Gene Features
              ├── Exons
              ├── UTRs
              └── Regulatory regions
```

---

## Accession Prefix System

### Nucleotide Sequences

| Prefix | Description | Status |
|--------|-------------|--------|
| NC_ | Complete genomic molecule | Curated |
| NG_ | Incomplete genomic region | Curated |
| NM_ | mRNA | Curated |
| NR_ | Non-coding RNA | Curated |
| NT_ | Contig (genomic) | Constructed |
| NW_ | WGS contig | Constructed |
| NZ_ | WGS supercontig | Constructed |
| XM_ | Predicted mRNA | Model |
| XR_ | Predicted ncRNA | Model |

### Protein Sequences

| Prefix | Description | Status |
|--------|-------------|--------|
| NP_ | Protein (from NM_) | Curated |
| XP_ | Predicted protein | Model |
| YP_ | Protein (bacterial/viral) | Curated |
| WP_ | Non-redundant protein | Prokaryotic |
| AP_ | Annotated on AC_ | Third-party |

---

## Core Tables/Entities

### Protein Record

**Description:** RefSeq protein sequence with annotations.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| accession | string | Yes | Accession with version (NP_000001.1) |
| gi | integer | Deprecated | GI number (legacy) |
| locus | string | No | Gene locus tag |
| definition | string | Yes | Protein description |
| organism | string | Yes | Source organism |
| taxonomy | array | Yes | Taxonomic lineage |
| gene | string | No | Gene symbol |
| gene_id | integer | No | NCBI Gene ID |
| coded_by | string | Yes | Source transcript |
| sequence | string | Yes | Amino acid sequence |
| length | integer | Yes | Sequence length |
| features | array | No | Annotated features |
| dbxrefs | array | No | Cross-references |

### Transcript Record

**Description:** RefSeq mRNA/RNA sequence.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| accession | string | Yes | Accession with version (NM_000001.3) |
| definition | string | Yes | Transcript description |
| organism | string | Yes | Source organism |
| gene | string | No | Gene symbol |
| gene_id | integer | No | NCBI Gene ID |
| chromosome | string | No | Chromosome location |
| map_location | string | No | Cytogenetic location |
| cds | object | No | CDS coordinates |
| protein_id | string | No | Encoded protein accession |
| exons | array | No | Exon coordinates |
| sequence | string | Yes | Nucleotide sequence |

### Feature Table Entry

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| key | string | Yes | Feature type |
| location | string | Yes | Feature location |
| qualifiers | object | No | Key-value annotations |

---

## API Endpoints (E-utilities)

| Endpoint | Method | Description |
|----------|--------|-------------|
| esearch.fcgi | GET | Search and get IDs |
| efetch.fcgi | GET | Retrieve records |
| einfo.fcgi | GET | Database information |
| elink.fcgi | GET | Find related records |
| esummary.fcgi | GET | Document summaries |

### Common Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| db | Database name | `protein`, `nuccore` |
| id | Record IDs | `NP_000001.1` |
| rettype | Return type | `fasta`, `gb`, `gp` |
| retmode | Return mode | `xml`, `json`, `text` |
| term | Search term | `TP53[gene] AND human[organism]` |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | GenBank/GenPept flat file |
| Alternative | FASTA, ASN.1, GFF3 |
| Encoding | UTF-8 |
| Compression | gzip on FTP |

---

## Sample Record

### GenPept Format (Abbreviated)

```
LOCUS       NP_000537                393 aa            linear   PRI 15-JAN-2026
DEFINITION  cellular tumor antigen p53 [Homo sapiens].
ACCESSION   NP_000537
VERSION     NP_000537.3
DBSOURCE    REFSEQ: accession NM_000546.6
KEYWORDS    RefSeq; MANE Select.
SOURCE      Homo sapiens (human)
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Mammalia; Primates; Hominidae; Homo.
REFERENCE   1  (residues 1 to 393)
  AUTHORS   Kastenhuber,E.R. and Lowe,S.W.
  TITLE     Putting p53 in Context
  JOURNAL   Cell 170 (6), 1062-1078 (2017)
COMMENT     VALIDATED REFSEQ: This record has undergone validation
FEATURES             Location/Qualifiers
     source          1..393
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
                     /chromosome="17"
     Protein         1..393
                     /product="cellular tumor antigen p53"
                     /note="tumor suppressor p53"
     Region          102..292
                     /region_name="P53"
                     /note="p53 DNA-binding domain"
                     /db_xref="CDD:cd08367"
     CDS             1..393
                     /gene="TP53"
                     /gene_synonym="P53"
                     /coded_by="NM_000546.6:203..1384"
                     /db_xref="CCDS:CCDS11118.1"
                     /db_xref="GeneID:7157"
ORIGIN
        1 meepqsdpsv epplsqetfs dlwkllpenn vlsplpsqam ddlmlspddi eqwftedpgp
       61 deaprmpeaa prvapapaap tpaapapaps wplsssvpsq ktyqgsygfr lgflhsgtak
      ...
//
```

### JSON Format (via E-utilities)

```json
{
  "id": "NP_000537",
  "version": 3,
  "accessionVersion": "NP_000537.3",
  "title": "cellular tumor antigen p53 [Homo sapiens]",
  "organism": "Homo sapiens",
  "taxid": 9606,
  "sourceDb": "refseq",
  "moltype": "aa",
  "slen": 393,
  "gene": "TP53",
  "geneId": 7157,
  "codedBy": "NM_000546.6",
  "createDate": "1999/02/10",
  "updateDate": "2026/01/15"
}
```

---

## Feature Types

| Feature | Description | Example |
|---------|-------------|---------|
| Protein | Full protein sequence | 1..393 |
| Region | Domain or region | P53 domain |
| Site | Binding/active site | Zn binding |
| mat_peptide | Mature peptide | After processing |
| sig_peptide | Signal peptide | Secretory signal |
| CDS | Coding sequence | Source transcript |

---

## Annotation Status

| Status | Description |
|--------|-------------|
| VALIDATED | Full manual review |
| REVIEWED | Curator reviewed |
| PROVISIONAL | Automated + limited review |
| PREDICTED | Computational only |
| MODEL | Gene model based |
| INFERRED | Inferred from related |

---

## Cross-References

| Database | ID Type | Example |
|----------|---------|---------|
| GeneID | NCBI Gene | 7157 |
| UniProtKB | Swiss-Prot/TrEMBL | P04637 |
| CCDS | Consensus CDS | CCDS11118.1 |
| Ensembl | Ensembl protein | ENSP00000269305 |
| HGNC | Gene nomenclature | HGNC:11998 |
| MIM | OMIM | 191170 |
| PDB | Structure | 1TUP |
| InterPro | Domain | IPR011615 |

---

## Glossary

| Term | Definition |
|------|------------|
| Accession | Stable identifier without version |
| Version | Sequence-specific version number |
| GI | GenInfo Identifier (deprecated) |
| MANE Select | Matched Annotation from NCBI and EBI |
| WP_ protein | Non-redundant prokaryotic protein |
| Coding sequence (CDS) | Region that encodes protein |

---

## References

1. O'Leary NA, et al. (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Research. https://doi.org/10.1093/nar/gkv1189
2. RefSeq FTP Site: https://ftp.ncbi.nlm.nih.gov/refseq/
3. RefSeq FAQ: https://www.ncbi.nlm.nih.gov/books/NBK50679/
