# Biomedical Data Integration - Data Dictionary

**Document ID:** DATA-DICT-001
**Version:** 1.0.0
**Generated:** 2026-01-24
**Status:** Final

---

## Executive Summary

This data dictionary provides semantic definitions, datatypes, and cardinality specifications for all fields used across 9 biomedical data categories in the Gene Platform taxonomy. It enables consistent data integration across 180+ data sources.

**Key Statistics:**
- **47 Shared Fields** - Fields appearing across 2+ categories
- **5 Identifier Hubs** - Gene, Protein, Compound, Disease, Publication
- **89 Cross-Category Links** - Foreign key relationships
- **45 Unified Schemas** - One per 3rd-level subcategory

---

## Table of Contents

1. [Identifier Hubs](#identifier-hubs)
2. [Shared Fields Catalogue](#shared-fields-catalogue)
3. [Field Definitions by Domain](#field-definitions-by-domain)
4. [Cross-Reference Mapping Tables](#cross-reference-mapping-tables)
5. [Cardinality Patterns](#cardinality-patterns)
6. [Integration Guidelines](#integration-guidelines)

---

## Identifier Hubs

### Gene Hub

| Property | Value |
|----------|-------|
| **Primary Identifier** | `gene_symbol` |
| **Secondary Identifiers** | `entrez_gene_id`, `ensembl_gene_id`, `hgnc_id`, `refseq_id` |
| **Categories Connected** | 8 of 9 categories |
| **Primary Authority** | HGNC (https://www.genenames.org/) |
| **Wikidata Property** | P351 |

**Mapping Strategies:**
- Direct symbol match (case-insensitive)
- HGNC ID to symbol lookup
- Ensembl to Entrez cross-reference
- Alias/synonym resolution

### Protein Hub

| Property | Value |
|----------|-------|
| **Primary Identifier** | `uniprot_accession` |
| **Secondary Identifiers** | `refseq_id`, `pdb_id`, `ensembl_protein_id` |
| **Categories Connected** | 6 categories |
| **Primary Authority** | UniProt Consortium (https://www.uniprot.org/) |
| **Wikidata Property** | P352 |

**Mapping Strategies:**
- UniProt accession direct lookup
- Gene symbol to protein mapping
- PDB to UniProt cross-reference
- Isoform resolution via UniProt

### Compound Hub

| Property | Value |
|----------|-------|
| **Primary Identifier** | `inchi_key` |
| **Secondary Identifiers** | `pubchem_cid`, `chembl_id`, `drugbank_id`, `chebi_id`, `hmdb_id`, `kegg_compound_id`, `cas_number` |
| **Categories Connected** | 5 categories |
| **Primary Authority** | IUPAC InChI Trust |
| **Wikidata Property** | P662 |

**Mapping Strategies:**
- InChI Key exact match (structure-based)
- PubChem CID lookup
- SMILES canonicalization and match
- Cross-database ID mapping via UniChem

### Disease Hub

| Property | Value |
|----------|-------|
| **Primary Identifier** | `mondo_id` |
| **Secondary Identifiers** | `omim_id`, `hpo_id`, `mesh_id`, `orphanet_id`, `icd10_code`, `doid`, `umls_cui`, `efo_id` |
| **Categories Connected** | 6 categories |
| **Primary Authority** | Mondo Disease Ontology |
| **Wikidata Property** | P5270 |

**Mapping Strategies:**
- MONDO ID as unifying identifier
- OMIM to MONDO mapping
- HPO phenotype to disease association
- ICD-10 to MONDO/OMIM cross-reference
- UMLS CUI semantic linking

### Publication Hub

| Property | Value |
|----------|-------|
| **Primary Identifier** | `pmid` |
| **Secondary Identifiers** | `doi`, `pmcid` |
| **Categories Connected** | 6 categories |
| **Primary Authority** | NCBI PubMed |
| **Wikidata Property** | P698 |

**Mapping Strategies:**
- PMID direct lookup
- DOI resolution to PMID
- PMC to PubMed cross-reference
- Title/author fuzzy matching

---

## Shared Fields Catalogue

### Top 20 Most Shared Fields (by source count)

| Rank | Field | Type | Sources | Categories | Format/Pattern |
|------|-------|------|---------|------------|----------------|
| 1 | `gene_symbol` | string | 78 | 7 | `^[A-Z][A-Z0-9-]+$` |
| 2 | `pmid` | integer | 65 | 6 | `^[0-9]{1,8}$` |
| 3 | `uniprot_accession` | string | 58 | 6 | UniProt AC pattern |
| 4 | `pubchem_cid` | integer | 52 | 5 | `^[0-9]+$` |
| 5 | `ncbi_taxonomy_id` | integer | 48 | 4 | `^[0-9]+$` |
| 6 | `entrez_gene_id` | integer | 45 | 5 | `^[0-9]+$` |
| 7 | `inchi_key` | string | 42 | 4 | `^[A-Z]{14}-[A-Z]{10}-[A-Z]$` |
| 8 | `mondo_id` | string | 38 | 4 | `^MONDO:[0-9]{7}$` |
| 9 | `chembl_id` | string | 36 | 4 | `^CHEMBL[0-9]+$` |
| 10 | `doi` | string | 35 | 3 | `^10\.[0-9]{4,}/[^\s]+$` |
| 11 | `ensembl_gene_id` | string | 34 | 4 | `^ENSG[0-9]{11}$` |
| 12 | `chebi_id` | string | 32 | 4 | `^CHEBI:[0-9]+$` |
| 13 | `canonical_smiles` | string | 31 | 4 | SMILES notation |
| 14 | `mesh_id` | string | 30 | 4 | `^D[0-9]{6}$` |
| 15 | `hpo_id` | string | 28 | 3 | `^HP:[0-9]{7}$` |
| 16 | `drugbank_id` | string | 28 | 4 | `^DB[0-9]{5}$` |
| 17 | `omim_id` | integer | 26 | 3 | `^[0-9]{6}$` |
| 18 | `hmdb_id` | string | 25 | 3 | `^HMDB[0-9]{7}$` |
| 19 | `go_term_id` | string | 24 | 3 | `^GO:[0-9]{7}$` |
| 20 | `kegg_pathway_id` | string | 23 | 3 | `^(hsa\|map\|ko)[0-9]{5}$` |

---

## Field Definitions by Domain

### Gene Identifiers

| Field | Type | Definition | Example | Authority |
|-------|------|------------|---------|-----------|
| `gene_symbol` | string | Official HGNC gene symbol | BRCA1 | HGNC |
| `entrez_gene_id` | integer | NCBI Entrez Gene ID | 672 | NCBI Gene |
| `ensembl_gene_id` | string | Ensembl gene identifier | ENSG00000141510 | Ensembl |
| `hgnc_id` | string | HGNC identifier | HGNC:1100 | HGNC |
| `refseq_id` | string | RefSeq accession | NM_000546.6 | NCBI RefSeq |

**Aliases:** gene, symbol, hgnc_symbol, gene_name, geneSymbol, SYMBOL

### Protein Identifiers

| Field | Type | Definition | Example | Authority |
|-------|------|------------|---------|-----------|
| `uniprot_accession` | string | UniProt KB accession | P04637 | UniProt |
| `pdb_id` | string | Protein Data Bank ID | 1TUP | RCSB PDB |
| `interpro_id` | string | InterPro family/domain ID | IPR000980 | InterPro |

**Aliases:** uniprot_id, UniProt_ID, accession, uniprot

### Compound Identifiers

| Field | Type | Definition | Example | Authority |
|-------|------|------------|---------|-----------|
| `inchi_key` | string | InChI hash key (27 char) | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | IUPAC |
| `pubchem_cid` | integer | PubChem Compound ID | 2244 | NCBI PubChem |
| `chembl_id` | string | ChEMBL molecule ID | CHEMBL25 | EMBL-EBI |
| `drugbank_id` | string | DrugBank ID | DB00945 | DrugBank |
| `chebi_id` | string | ChEBI ontology ID | CHEBI:15377 | EMBL-EBI |
| `hmdb_id` | string | Human Metabolome DB ID | HMDB0000001 | HMDB |
| `kegg_compound_id` | string | KEGG Compound ID | C00031 | KEGG |
| `cas_number` | string | CAS Registry Number | 50-78-2 | CAS |
| `canonical_smiles` | string | Canonical SMILES | CC(=O)OC1=CC=CC=C1C(=O)O | Daylight |
| `molecular_weight` | number | Molecular weight (Da) | 180.16 | - |
| `molecular_formula` | string | Molecular formula | C6H12O6 | - |

### Disease Identifiers

| Field | Type | Definition | Example | Authority |
|-------|------|------------|---------|-----------|
| `mondo_id` | string | Mondo Disease Ontology ID | MONDO:0005015 | Monarch |
| `omim_id` | integer | OMIM identifier | 114480 | OMIM |
| `hpo_id` | string | Human Phenotype Ontology ID | HP:0001250 | HPO |
| `mesh_id` | string | MeSH descriptor ID | D001943 | NLM |
| `orphanet_id` | string | Orphanet rare disease ID | ORPHA:558 | Orphanet |
| `icd10_code` | string | ICD-10 diagnosis code | E11.9 | WHO |
| `doid` | string | Disease Ontology ID | DOID:9352 | DO |
| `umls_cui` | string | UMLS Concept ID | C0011849 | NLM UMLS |
| `efo_id` | string | Experimental Factor Ontology ID | EFO:0000400 | EBI |
| `clinical_significance` | string | Variant pathogenicity | pathogenic | ClinVar |

**Allowed Values for clinical_significance:** pathogenic, likely_pathogenic, uncertain_significance, likely_benign, benign

### Genomic Position

| Field | Type | Definition | Example | Notes |
|-------|------|------------|---------|-------|
| `chromosome` | string | Chromosome identifier | chr17 | Pattern: `^(chr)?([1-9]\|1[0-9]\|2[0-2]\|X\|Y\|MT?)$` |
| `position` | integer | Genomic coordinate (1-based) | 7577120 | - |
| `variant_id` | string | dbSNP rsID | rs7412 | Pattern: `^rs[0-9]+$` |
| `clinvar_id` | integer | ClinVar variation ID | 12345 | - |

### Pathway Identifiers

| Field | Type | Definition | Example | Authority |
|-------|------|------------|---------|-----------|
| `kegg_pathway_id` | string | KEGG pathway ID | hsa04151 | KEGG |
| `reactome_id` | string | Reactome pathway ID | R-HSA-109582 | Reactome |
| `go_term_id` | string | Gene Ontology term ID | GO:0008150 | GO Consortium |
| `kegg_orthology_id` | string | KEGG Orthology ID | K00001 | KEGG |

### Literature Identifiers

| Field | Type | Definition | Example | Authority |
|-------|------|------------|---------|-----------|
| `pmid` | integer | PubMed ID | 12345678 | NCBI PubMed |
| `doi` | string | Digital Object Identifier | 10.1038/nature12373 | DOI Foundation |
| `pmcid` | string | PubMed Central ID | PMC3531190 | NCBI PMC |

### Organism/Taxonomy

| Field | Type | Definition | Example | Authority |
|-------|------|------------|---------|-----------|
| `ncbi_taxonomy_id` | integer | NCBI Taxonomy ID | 9606 | NCBI Taxonomy |
| `organism` | string | Scientific name | Homo sapiens | - |

### Quality/Score Fields

| Field | Type | Definition | Example | Range |
|-------|------|------------|---------|-------|
| `confidence_score` | number | Interaction confidence | 0.95 | 0-1 |

---

## Cross-Reference Mapping Tables

### Gene → Other Domains

| From Field | To Domain | Link Type | Cardinality |
|------------|-----------|-----------|-------------|
| gene_symbol | Diseases | gene_disease_association | N:M |
| gene_symbol | Pathways | gene_pathway_membership | N:M |
| gene_symbol | Proteins | gene_protein_encoding | 1:N |
| gene_symbol | Traditional Medicine | target_compound_association | N:M |
| gene_symbol | Microbiome | microbiome_gene_interaction | N:M |

### Protein → Other Domains

| From Field | To Domain | Link Type | Cardinality |
|------------|-----------|-----------|-------------|
| uniprot_accession | Compounds | drug_target_binding | N:M |
| uniprot_accession | Pathways | protein_pathway_function | N:M |
| uniprot_accession | Traditional Medicine | herbal_target_interaction | N:M |

### Compound → Other Domains

| From Field | To Domain | Link Type | Cardinality |
|------------|-----------|-----------|-------------|
| inchi_key | Traditional Medicine | structure_based_identity | 1:1 |
| inchi_key | Nutrition | structure_based_identity | 1:1 |
| pubchem_cid | Traditional Medicine | compound_source_identification | 1:N |
| pubchem_cid | Nutrition | food_compound_identity | 1:N |
| pubchem_cid | Microbiome | microbiome_metabolite_identity | 1:N |

### Disease → Other Domains

| From Field | To Domain | Link Type | Cardinality |
|------------|-----------|-----------|-------------|
| mondo_id | Pathways | disease_pathway_mechanism | N:M |
| mondo_id | Literature | disease_literature_annotation | 1:N |
| mondo_id | Microbiome | microbiome_disease_association | N:M |
| hpo_id | Genetics | phenotype_genotype_association | N:M |

---

## Cardinality Patterns

### Cardinality Definitions

| Symbol | Meaning | Description |
|--------|---------|-------------|
| 1:1 | One-to-One | Exact identity match (e.g., InChI Key) |
| 1:N | One-to-Many | One source entity maps to multiple targets |
| N:M | Many-to-Many | Complex associations (e.g., gene-disease) |

### Common Patterns by Relationship Type

| Relationship | Cardinality | Example |
|--------------|-------------|---------|
| Structure identity | 1:1 | InChI Key matching |
| Gene to proteins | 1:N | One gene encodes multiple isoforms |
| Gene to diseases | N:M | Pleiotropy and polygenic diseases |
| Compound to targets | N:M | Drug polypharmacology |
| Disease to pathways | N:M | Multiple pathways per disease |
| Literature to entities | 1:N | One paper references multiple entities |

---

## Integration Guidelines

### Identifier Resolution Priority

1. **Compounds:** Use InChI Key for exact structural identity
2. **Genes:** Use gene_symbol with HGNC validation
3. **Proteins:** Use UniProt accession as primary
4. **Diseases:** Use MONDO ID for unification (maps to OMIM, Orphanet, etc.)
5. **Literature:** Use PMID (most comprehensive coverage)

### Cross-Reference Strategies

| Strategy | Method | Fallback |
|----------|--------|----------|
| Structure-based compound linking | Match InChI Keys for exact identity | Canonical SMILES, then PubChem CID |
| Gene-centric integration | HGNC symbol as primary | Alias resolution via HGNC |
| Disease ontology unification | MONDO as integration layer | MeSH or UMLS CUI |
| Literature evidence linking | PMID as primary key | DOI via CrossRef |

### Data Quality Recommendations

1. Always validate identifiers against primary authorities
2. Use confidence scores to filter low-quality associations
3. Track provenance through source database citations
4. Apply cardinality constraints during data integration
5. Implement version tracking for ontology terms

### API Endpoints for Identifier Resolution

| Domain | Endpoint |
|--------|----------|
| Gene | https://rest.genenames.org/search/ |
| Protein | https://rest.uniprot.org/uniprotkb/ |
| Compound | https://pubchem.ncbi.nlm.nih.gov/rest/pug/ |
| Disease | https://api.monarchinitiative.org/api/ |
| Literature | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ |

---

## Appendix: Field Alias Mappings

### Gene Symbol Aliases
`gene`, `symbol`, `hgnc_symbol`, `gene_name`, `geneSymbol`, `SYMBOL`

### UniProt Accession Aliases
`uniprot_id`, `UniProt_ID`, `uniprotId`, `accession`, `uniprot`, `UniProtKB`

### PubChem CID Aliases
`PubChem_CID`, `pubchem_compound_id`, `cid`, `compound_id`

### SMILES Aliases
`smiles`, `SMILES`, `isomeric_smiles`, `smiles_string`

### PMID Aliases
`pubmed_id`, `PubMed_ID`, `pubmedId`, `PMID`, `pubmed`

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0.0 | 2026-01-24 | Schema Unification Swarm | Initial data dictionary from 130 schema files |
