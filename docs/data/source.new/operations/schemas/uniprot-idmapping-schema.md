---
id: schemas-uniprot-idmapping-schema
title: "UniProt ID Mapping Schema"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, uniprot, id-mapping, cross-reference, protein, integration]
---

**Parent:** [Schema Documentation](./_index.md)

# UniProt ID Mapping Schema

**Document ID:** UNIPROT-IDMAPPING-SCHEMA
**Status:** Final
**Last Updated:** January 2026
**Version:** 1.0

---

## Overview

UniProt ID Mapping is the master cross-reference service linking 286 databases to UniProt protein entries. It provides the most comprehensive protein-centric ID mapping available, with pre-computed files updated every 8 weeks synchronized with UniProt releases.

### Service Endpoints

| Resource | URL |
|----------|-----|
| **Web Interface** | https://www.uniprot.org/id-mapping |
| **REST API** | https://rest.uniprot.org/idmapping/ |
| **FTP Downloads** | ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ |
| **Field Configuration** | https://rest.uniprot.org/configure/idmapping/fields |

---

## File Formats

### idmapping.dat (~15 GB compressed)

Three-column tab-delimited file containing all mappings.

| Column | Name | Description | Example |
|--------|------|-------------|---------|
| 1 | UniProtKB-AC | UniProt accession | P04637 |
| 2 | ID_type | External database name | GeneID |
| 3 | ID | External identifier | 7157 |

**Sample Records:**
```
P04637	UniProtKB-ID	P53_HUMAN
P04637	GeneID	7157
P04637	RefSeq	NP_000537.3
P04637	RefSeq_NT	NM_000546.6
P04637	PDB	1AIE
P04637	PDB	1C26
P04637	GO	GO:0000785
P04637	GO	GO:0001046
P04637	Ensembl	ENST00000269305
P04637	HGNC	HGNC:11998
```

### idmapping_selected.tab (~3 GB compressed)

Pre-computed 22-column table with most frequently requested mappings.

| Col | Field Name | Description | Example |
|-----|------------|-------------|---------|
| 1 | UniProtKB-AC | UniProt accession | P04637 |
| 2 | UniProtKB-ID | UniProt entry name | P53_HUMAN |
| 3 | GeneID | NCBI Entrez Gene | 7157 |
| 4 | RefSeq | RefSeq protein IDs | NP_000537.3 |
| 5 | GI | GenInfo Identifier (deprecated) | 120407068 |
| 6 | PDB | Protein Data Bank | 1AIE;1C26;1DT7 |
| 7 | GO | Gene Ontology terms | GO:0000785;GO:0001046 |
| 8 | UniRef100 | UniRef100 cluster | UniRef100_P04637 |
| 9 | UniRef90 | UniRef90 cluster | UniRef90_P04637 |
| 10 | UniRef50 | UniRef50 cluster | UniRef50_P04637 |
| 11 | UniParc | UniParc ID | UPI000002ED67 |
| 12 | PIR | PIR accession | A25224 |
| 13 | NCBI-taxon | Taxonomy ID | 9606 |
| 14 | MIM | OMIM ID | 191170 |
| 15 | UniGene | UniGene cluster | Hs.408312 |
| 16 | PubMed | PubMed IDs | 11433014;11488916 |
| 17 | EMBL | EMBL accession | AB082923;AF052180 |
| 18 | EMBL-CDS | EMBL CDS IDs | BAC16799.1;AAC12971.1 |
| 19 | Ensembl | Ensembl transcript | ENST00000269305 |
| 20 | Ensembl_TRS | Ensembl transcript | ENST00000269305 |
| 21 | Ensembl_PRO | Ensembl protein | ENSP00000269305 |
| 22 | Additional_PubMed | Extra PubMed refs | Additional references |

---

## Complete Database List (286 Databases)

### Sequence Databases (8)

| Database | ID Type | Description |
|----------|---------|-------------|
| EMBL | EMBL/GenBank/DDBJ | Nucleotide sequence |
| GenBank | GenBank | NCBI nucleotide |
| DDBJ | DDBJ | DNA Data Bank of Japan |
| RefSeq | RefSeq | NCBI Reference Sequences |
| RefSeq_NT | RefSeq Nucleotide | RefSeq nucleotide |
| PIR | PIR | Protein Information Resource |
| UniParc | UniParc | UniProt Archive |
| UniRef | UniRef clusters | UniRef100/90/50 |

### Gene Databases (12)

| Database | ID Type | Description |
|----------|---------|-------------|
| GeneID | NCBI Entrez Gene | Primary gene identifier |
| HGNC | HGNC | Human gene nomenclature |
| Ensembl | Ensembl | Ensembl gene/transcript |
| KEGG | KEGG | KEGG gene |
| MGI | MGI | Mouse Genome Informatics |
| RGD | RGD | Rat Genome Database |
| SGD | SGD | Saccharomyces Genome Database |
| FlyBase | FlyBase | Drosophila database |
| WormBase | WormBase | C. elegans database |
| ZFIN | ZFIN | Zebrafish database |
| Xenbase | Xenbase | Xenopus database |
| dictyBase | dictyBase | Dictyostelium database |

### Protein Databases (15)

| Database | ID Type | Description |
|----------|---------|-------------|
| PDB | PDB | Protein Data Bank |
| AlphaFoldDB | AlphaFold | AlphaFold structure predictions |
| SMR | SWISS-MODEL | Swiss-Model Repository |
| BioGRID | BioGRID | Biological interaction database |
| MINT | MINT | Molecular interactions |
| IntAct | IntAct | Molecular interactions |
| STRING | STRING | Protein interactions |
| ComplexPortal | Complex Portal | Protein complexes |
| CORUM | CORUM | Comprehensive Resource of Mammalian Complexes |
| DIP | DIP | Database of Interacting Proteins |
| PRIDE | PRIDE | Proteomics data |
| MassIVE | MassIVE | Mass spectrometry data |
| ProteomicsDB | ProteomicsDB | Proteomics database |
| MaxQB | MaxQB | MaxQuant database |
| PeptideAtlas | PeptideAtlas | Peptide identification |

### Disease & Variant Databases (18)

| Database | ID Type | Description |
|----------|---------|-------------|
| MIM | OMIM | Online Mendelian Inheritance in Man |
| Orphanet | Orphanet | Rare disease database |
| DisGeNET | DisGeNET | Gene-disease associations |
| MalaCards | MalaCards | Human disease database |
| GeneCards | GeneCards | Human gene compendium |
| PharmGKB | PharmGKB | Pharmacogenomics database |
| ClinVar | ClinVar | Clinical variants |
| COSMIC | COSMIC | Cancer mutations |
| dbSNP | dbSNP | Single nucleotide polymorphisms |
| OpenTargets | Open Targets | Drug target validation |
| Reactome | Reactome | Pathway database |
| WikiPathways | WikiPathways | Community pathway database |
| KEGG | KEGG Pathway | KEGG pathways |
| BioCyc | BioCyc | Pathway/genome databases |
| CTD | CTD | Comparative Toxicogenomics |
| DrugBank | DrugBank | Drug database |
| ChEMBL | ChEMBL | Bioactivity database |
| BindingDB | BindingDB | Binding affinity data |

### Ontology & Annotation (12)

| Database | ID Type | Description |
|----------|---------|-------------|
| GO | Gene Ontology | Gene Ontology terms |
| Pfam | Pfam | Protein families |
| InterPro | InterPro | Protein domains |
| PROSITE | PROSITE | Protein patterns |
| SMART | SMART | Domain architecture |
| SUPFAM | SUPERFAMILY | Structural domains |
| PANTHER | PANTHER | Protein classification |
| CDD | CDD | Conserved domains |
| eggNOG | eggNOG | Orthologous groups |
| OMA | OMA | Orthologous Matrix |
| OrthoDB | OrthoDB | Orthologs database |
| HOGENOM | HOGENOM | Homologous genes |

---

## REST API Usage

### Submit Job

```bash
curl --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'ids="P04637,P53_HUMAN,7157"' \
  --form 'from="UniProtKB_AC-ID"' \
  --form 'to="GeneID"'
```

### Check Job Status

```bash
curl 'https://rest.uniprot.org/idmapping/status/{jobId}'
```

### Get Results

```bash
curl 'https://rest.uniprot.org/idmapping/results/{jobId}'
```

---

## Integration Notes

- **Update Frequency**: Every 8 weeks with UniProt releases
- **Coverage**: 286 databases, 250M+ cross-references
- **Primary Hub**: UniProt accession (P#####)
- **Best for**: Protein-centric ID resolution
- **Limitations**: Gene-level mappings require NCBI Gene as intermediary

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `UniProtKB-AC` | UniProt Knowledge Base Accession (primary protein ID) | P04637 |
| `UniProtKB-ID` | UniProt entry name (mnemonic identifier) | P53_HUMAN |
| `GeneID` | NCBI Entrez Gene identifier | 7157 |
| `RefSeq` | NCBI Reference Sequence protein identifier | NP_000537.3 |
| `PDB` | Protein Data Bank structure identifier | 1AIE |
| `UniRef100` | UniRef cluster at 100% sequence identity | UniRef100_P04637 |
| `UniRef90` | UniRef cluster at 90% sequence identity | UniRef90_P04637 |
| `UniRef50` | UniRef cluster at 50% sequence identity | UniRef50_P04637 |
| `UniParc` | UniProt Archive unique protein sequence ID | UPI000002ED67 |
| `HGNC` | HUGO Gene Nomenclature Committee identifier | HGNC:11998 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| ID Mapping | Converting identifiers between different databases | Core service function |
| Cross-Reference | Link between databases using equivalent identifiers | idmapping.dat |
| Protein-Centric | Organized around protein as primary entity | UniProt design philosophy |
| Accession Number | Stable identifier that persists across releases | UniProtKB-AC |
| Entry Name | Human-readable mnemonic identifier | UniProtKB-ID |
| Sequence Identity | Percentage of matching amino acids between proteins | UniRef clustering |
| Selected Mappings | Pre-computed common ID mappings | idmapping_selected.tab |
| Master Cross-Reference | Comprehensive mapping hub | 286 databases |
| Batch Conversion | Converting multiple IDs in single request | REST API feature |
| Release Synchronization | Updates aligned with UniProt releases | Every 8 weeks |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| UniProt | Universal Protein Resource | Protein knowledge base |
| UniProtKB | UniProt Knowledge Base | Curated protein database |
| UniRef | UniProt Reference Clusters | Sequence clusters |
| UniParc | UniProt Archive | Comprehensive sequence archive |
| HGNC | HUGO Gene Nomenclature Committee | Human gene naming |
| MGI | Mouse Genome Informatics | Mouse gene database |
| RGD | Rat Genome Database | Rat gene database |
| SGD | Saccharomyces Genome Database | Yeast database |
| PDB | Protein Data Bank | Structure database |
| GO | Gene Ontology | Function annotation |
| OMIM | Online Mendelian Inheritance in Man | Disease database |
| MIM | Mendelian Inheritance in Man | OMIM identifier |
| Pfam | Protein Families | Domain database |
| InterPro | Integrated Protein Resource | Domain integration |
| CDD | Conserved Domain Database | NCBI domain database |
| OMA | Orthologous Matrix | Ortholog database |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
