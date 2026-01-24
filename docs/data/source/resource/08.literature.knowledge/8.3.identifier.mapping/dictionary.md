# 8.3 Identifier Mapping - Data Dictionary

**Subcategory ID:** 8.3
**Subcategory Name:** Identifier Mapping
**Data Sources:** NCBI E-Link, PMC ID Converter, UniProt ID Mapping

## Overview

This subcategory provides cross-reference mapping services that enable translation between identifiers from different biological databases. These services are essential for integrating data across the diverse landscape of genomic, proteomic, and literature databases.

---

## Unified Fields

### Core Mapping Fields

| Field | Data Type | Required | Description | Example |
|-------|-----------|----------|-------------|---------|
| `source_id` | string | Yes | Input identifier to be mapped | `12345678`, `P04637` |
| `source_database` | string | Yes | Database the source identifier belongs to | `pubmed`, `UniProtKB_AC-ID` |
| `target_ids` | array[object] | Yes | Mapped identifiers in target databases | See structure below |
| `mapping_score` | number | No | Confidence or relevance score for the mapping | `1000`, `95` |

**Target ID Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `id` | string | Target identifier value |
| `database` | string | Target database name |
| `score` | number | Relevance score (if available) |
| `link_name` | string | Named link type (NCBI E-Link) |

**Source Mappings - Core Fields:**

| Field | NCBI E-Link | PMC ID Converter | UniProt ID Mapping |
|-------|-------------|------------------|-------------------|
| `source_id` | Id | requested-id | from_id |
| `source_database` | DbFrom | idtype | from |
| `target_ids[].id` | Link/Id | pmid/pmcid/doi | to_id |
| `target_ids[].database` | DbTo | idtype | to |
| `mapping_score` | Score | - | - |

---

## NCBI E-Link Specific Fields

### Link Information

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `ncbi_elink_data.link_name` | string | Named link type | `pubmed_gene`, `gene_protein` |
| `ncbi_elink_data.has_linkout` | boolean | Whether record has external links | `true` |
| `ncbi_elink_data.has_neighbor` | boolean | Whether record has related records | `true` |
| `ncbi_elink_data.linkset_db` | array[string] | Array of link sets to target databases | `["gene", "protein", "snp"]` |
| `ncbi_elink_data.query_key` | string | History server query key | `1` |

### External LinkOut

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `ncbi_elink_data.linkouts` | array[object] | External LinkOut URLs | See structure below |

**LinkOut Object Structure:**

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| `url` | string | External resource URL | `https://www.ncbi.nlm.nih.gov/...` |
| `provider_name` | string | LinkOut provider name | `PubMed Central` |
| `category` | string | LinkOut category | `Full Text Sources`, `Medical` |
| `attribute` | string | Link attribute | `free resource`, `subscription` |

### Cross-Database Links

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `ncbi_elink_data.related_pubmed` | array[string] | Related PubMed article IDs | `["12345678", "87654321"]` |
| `ncbi_elink_data.linked_genes` | array[string] | Associated gene IDs | `["7157", "7124"]` |
| `ncbi_elink_data.linked_proteins` | array[string] | Associated protein IDs | `["P04637", "P53_HUMAN"]` |
| `ncbi_elink_data.linked_clinvar` | array[string] | Clinical variant IDs | `["VCV000012345"]` |
| `ncbi_elink_data.linked_snps` | array[string] | SNP variant IDs | `["rs1042522", "rs28934576"]` |

### Common NCBI Link Names

| Link Name | From Database | To Database | Description |
|-----------|---------------|-------------|-------------|
| pubmed_gene | PubMed | Gene | Articles to associated genes |
| pubmed_protein | PubMed | Protein | Articles to associated proteins |
| pubmed_pubmed | PubMed | PubMed | Related articles |
| gene_pubmed | Gene | PubMed | Gene to associated articles |
| gene_protein | Gene | Protein | Gene to protein products |
| gene_clinvar | Gene | ClinVar | Gene to clinical variants |
| gene_snp | Gene | dbSNP | Gene to SNP variants |
| protein_gene | Protein | Gene | Protein to encoding gene |
| protein_structure | Protein | Structure | Protein to 3D structures |

---

## PMC ID Converter Specific Fields

### Identifier Conversions

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `pmc_converter_data.pmid` | string | PubMed identifier | `12345678` |
| `pmc_converter_data.pmcid` | string | PubMed Central identifier (with PMC prefix) | `PMC1234567` |
| `pmc_converter_data.doi` | string | Digital Object Identifier | `10.1000/example` |
| `pmc_converter_data.mid` | string | NIH Manuscript ID (NIHMS format) | `NIHMS123456` |

### Article Status

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `pmc_converter_data.live` | boolean | Whether article is publicly accessible | `true` |
| `pmc_converter_data.status` | string | Article status | `current`, `retracted`, `removed`, `preprint` |
| `pmc_converter_data.release_date` | date | Public release date (YYYY-MM-DD) | `2024-01-15` |

### Version and Correction Information

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `pmc_converter_data.versions` | array[string] | Version history (PMCID.version format) | `["PMC1234567.1", "PMC1234567.2"]` |
| `pmc_converter_data.errata` | object | Correction/erratum information | `{pmid, doi}` |
| `pmc_converter_data.error_message` | string | Error message if mapping failed | `ID not found` |

**Errata Object Structure:**

| Property | Data Type | Description |
|----------|-----------|-------------|
| `pmid` | string | PubMed ID of correction |
| `doi` | string | DOI of correction |

### Valid ID Types for PMC Converter

| ID Type | Description | Format |
|---------|-------------|--------|
| pmid | PubMed ID | `[0-9]+` |
| pmcid | PubMed Central ID | `PMC[0-9]+` |
| doi | Digital Object Identifier | `10\.[0-9]+/.+` |
| mid | NIH Manuscript ID | `NIHMS[0-9]+` |

---

## UniProt ID Mapping Specific Fields

### Protein Identifiers

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `uniprot_mapping_data.uniprot_accession` | string | UniProtKB accession (primary protein ID) | `P04637` |
| `uniprot_mapping_data.uniprot_id` | string | UniProt entry name (mnemonic) | `P53_HUMAN` |
| `uniprot_mapping_data.gene_id` | string | NCBI Entrez Gene identifier | `7157` |
| `uniprot_mapping_data.refseq` | string | NCBI RefSeq protein accession | `NP_000537.3` |

### Sequence Archive Identifiers

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `uniprot_mapping_data.uniref100` | string | UniRef100 cluster ID | `UniRef100_P04637` |
| `uniprot_mapping_data.uniref90` | string | UniRef90 cluster ID | `UniRef90_P04637` |
| `uniprot_mapping_data.uniref50` | string | UniRef50 cluster ID | `UniRef50_P04637` |
| `uniprot_mapping_data.uniparc` | string | UniParc sequence archive ID | `UPI000002ED67` |
| `uniprot_mapping_data.embl` | string | EMBL/GenBank nucleotide accession | `M14695` |

### Structure and Domain Identifiers

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `uniprot_mapping_data.pdb_ids` | array[string] | Protein Data Bank structure IDs | `["1TUP", "2OCJ", "3KMD"]` |
| `uniprot_mapping_data.pfam` | array[string] | Pfam protein family IDs | `["PF00870", "PF07710"]` |
| `uniprot_mapping_data.interpro` | array[string] | InterPro domain IDs | `["IPR011615", "IPR012346"]` |

### Functional Annotations

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `uniprot_mapping_data.go_terms` | array[string] | Gene Ontology term accessions | `["GO:0003677", "GO:0005634"]` |
| `uniprot_mapping_data.reactome` | array[string] | Reactome pathway IDs | `["R-HSA-69563", "R-HSA-5357801"]` |

### Gene and Nomenclature

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `uniprot_mapping_data.hgnc` | string | HGNC identifier | `HGNC:11998` |
| `uniprot_mapping_data.ensembl_transcript` | string | Ensembl transcript ID | `ENST00000269305` |
| `uniprot_mapping_data.ensembl_protein` | string | Ensembl protein ID | `ENSP00000269305` |

### Taxonomy and Disease

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `uniprot_mapping_data.ncbi_taxon` | string | NCBI Taxonomy ID | `9606` |
| `uniprot_mapping_data.omim` | string | OMIM disease ID | `191170` |
| `uniprot_mapping_data.pubmed_refs` | array[string] | PubMed reference IDs | `["10617473", "8451856"]` |

### Drug/Compound Identifiers

| Field | Data Type | Description | Example |
|-------|-----------|-------------|---------|
| `uniprot_mapping_data.drugbank` | array[string] | DrugBank drug IDs | `["DB00945", "DB01050"]` |
| `uniprot_mapping_data.chembl` | array[string] | ChEMBL compound IDs | `["CHEMBL25", "CHEMBL941"]` |

### Supported Database Types for UniProt Mapping

| From/To | Database Name | Description |
|---------|---------------|-------------|
| UniProtKB_AC-ID | UniProtKB | UniProt Knowledgebase accession |
| UniProtKB-ID | UniProtKB | UniProt entry name |
| GeneID | NCBI Gene | Entrez Gene identifier |
| RefSeq_Protein | RefSeq | RefSeq protein accession |
| PDB | Protein Data Bank | 3D structure identifier |
| GO | Gene Ontology | GO term accession |
| UniRef100/90/50 | UniRef | UniRef cluster identifier |
| UniParc | UniParc | Sequence archive ID |
| EMBL | EMBL-ENA | Nucleotide accession |
| Ensembl_TRS | Ensembl | Transcript identifier |
| Ensembl_PRO | Ensembl | Protein identifier |
| HGNC | HGNC | Human gene nomenclature |
| Pfam | Pfam | Protein family |
| InterPro | InterPro | Domain annotation |
| DrugBank | DrugBank | Drug identifier |
| ChEMBL | ChEMBL | Compound identifier |
| Reactome | Reactome | Pathway identifier |
| OMIM | OMIM | Mendelian inheritance |

---

## Source Mappings Summary

### NCBI E-Link Field Mappings

| Unified Field | NCBI E-Link Source |
|---------------|-------------------|
| `source_id` | Id |
| `source_database` | DbFrom |
| `target_ids[].id` | Link/Id |
| `target_ids[].database` | DbTo |
| `mapping_score` | Score |
| `ncbi_elink_data.link_name` | link_name |
| `ncbi_elink_data.has_linkout` | HasLinkOut |
| `ncbi_elink_data.has_neighbor` | HasNeighbor |
| `ncbi_elink_data.linkset_db` | LinkSetDb |
| `ncbi_elink_data.query_key` | QueryKey |
| `ncbi_elink_data.linkouts` | ObjUrl |
| `ncbi_elink_data.related_pubmed` | pubmed_pubmed |
| `ncbi_elink_data.linked_genes` | pubmed_gene |
| `ncbi_elink_data.linked_proteins` | pubmed_protein |
| `ncbi_elink_data.linked_clinvar` | gene_clinvar |
| `ncbi_elink_data.linked_snps` | gene_snp |

### PMC ID Converter Field Mappings

| Unified Field | PMC ID Converter Source |
|---------------|------------------------|
| `source_id` | requested-id |
| `source_database` | idtype |
| `pmc_converter_data.pmid` | pmid |
| `pmc_converter_data.pmcid` | pmcid |
| `pmc_converter_data.doi` | doi |
| `pmc_converter_data.mid` | mid |
| `pmc_converter_data.live` | live |
| `pmc_converter_data.status` | status |
| `pmc_converter_data.errata` | errata |
| `pmc_converter_data.versions` | versions |
| `pmc_converter_data.release_date` | release-date |
| `pmc_converter_data.error_message` | errmsg |

### UniProt ID Mapping Field Mappings

| Unified Field | UniProt Source |
|---------------|----------------|
| `source_id` | from_id |
| `source_database` | from |
| `target_ids[].id` | to_id |
| `target_ids[].database` | to |
| `uniprot_mapping_data.uniprot_accession` | UniProtKB_AC-ID |
| `uniprot_mapping_data.uniprot_id` | UniProtKB-ID |
| `uniprot_mapping_data.gene_id` | GeneID |
| `uniprot_mapping_data.refseq` | RefSeq |
| `uniprot_mapping_data.pdb_ids` | PDB |
| `uniprot_mapping_data.go_terms` | GO |
| `uniprot_mapping_data.uniref100` | UniRef100 |
| `uniprot_mapping_data.uniref90` | UniRef90 |
| `uniprot_mapping_data.uniref50` | UniRef50 |
| `uniprot_mapping_data.uniparc` | UniParc |
| `uniprot_mapping_data.ncbi_taxon` | NCBI_TaxID |
| `uniprot_mapping_data.omim` | OMIM |
| `uniprot_mapping_data.pubmed_refs` | PubMed |
| `uniprot_mapping_data.embl` | EMBL |
| `uniprot_mapping_data.ensembl_transcript` | Ensembl_TRS |
| `uniprot_mapping_data.ensembl_protein` | Ensembl_PRO |
| `uniprot_mapping_data.hgnc` | HGNC |
| `uniprot_mapping_data.pfam` | Pfam |
| `uniprot_mapping_data.interpro` | InterPro |
| `uniprot_mapping_data.drugbank` | DrugBank |
| `uniprot_mapping_data.chembl` | ChEMBL |
| `uniprot_mapping_data.reactome` | Reactome |

---

## Data Source Metadata

| Field | Data Type | Description |
|-------|-----------|-------------|
| `_source.primary_source` | string | Name of the primary mapping service |
| `_source.source_id` | string | Original ID in the source system |
| `_source.extraction_date` | date | Date when data was extracted |
| `_source.source_version` | string | Version of the source data/API |

---

## API Endpoints

### NCBI E-Link

- **Base URL:** `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi`
- **Parameters:** `dbfrom`, `db`, `id`, `cmd`, `linkname`

### PMC ID Converter

- **Base URL:** `https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/`
- **Parameters:** `ids`, `idtype`, `format`, `versions`

### UniProt ID Mapping

- **Base URL:** `https://rest.uniprot.org/idmapping/`
- **Endpoints:** `/run`, `/status/{jobId}`, `/results/{jobId}`
- **Parameters:** `from`, `to`, `ids`
