# NPASS - Data Dictionary

## Overview

This data dictionary documents the schema for NPASS (Natural Product Activity and Species Source) database.

| Property | Value |
|----------|-------|
| **Level** | source |
| **ID** | npass |
| **Name** | NPASS |
| **Parent** | 2.1.natural.products |
| **Total Fields** | 25+ |
| **Last Updated** | 2026-01-25 |

---

## Core Fields

### Compound Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| npc_id | string | 1:1 | Yes | NPASS compound identifier | NPC12345 |
| name | string | 1:1 | No | Compound name | Artemisinin |
| canonical_smiles | string | 1:1 | Yes | Canonical SMILES structure | CC1CCC2C(C)... |
| inchi | string | 1:1 | No | InChI identifier | InChI=1S/... |
| inchi_key | string | 1:1 | Yes | InChI Key for lookups | BLUAFEHZUWYNDE-... |
| molecular_formula | string | 1:1 | No | Molecular formula | C15H22O5 |
| molecular_weight | decimal | 1:1 | No | Molecular weight (Da) | 282.33 |
| compound_class | string | 1:1 | No | Natural product class | Terpenoid |

### Activity Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| activity_id | integer | 1:1 | Yes | Primary identifier | 1001 |
| activity_type | string | 1:1 | Yes | Measurement type | IC50 |
| activity_value | decimal | 1:1 | Yes | Numeric value | 4.8 |
| activity_unit | string | 1:1 | Yes | Unit of measurement | nM |
| activity_relation | string | 1:1 | No | Relation operator | =, <, >, ~ |
| assay_type | string | 1:1 | No | Assay category | Binding, Functional |
| reference | string | 1:1 | No | PubMed ID or DOI | PMID:12345678 |

### Target Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| target_id | integer | 1:1 | Yes | Primary identifier | 501 |
| target_name | string | 1:1 | Yes | Target name | Plasmodium falciparum |
| uniprot_id | string | 1:1 | No | UniProt accession | P12345 |
| gene_symbol | string | 1:1 | No | Gene symbol | PfATP6 |
| organism | string | 1:1 | No | Target organism | Plasmodium falciparum |
| target_type | string | 1:1 | No | Target classification | Protein, Enzyme, Receptor |

### Species Information

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| species_id | integer | 1:1 | Yes | Primary identifier | 301 |
| scientific_name | string | 1:1 | Yes | Binomial name | Artemisia annua |
| ncbi_taxon_id | integer | 1:1 | No | NCBI Taxonomy ID | 35608 |
| kingdom | string | 1:1 | No | Taxonomic kingdom | Plantae |
| family | string | 1:1 | No | Taxonomic family | Asteraceae |

### Compound-Species Link

| Field Name | Data Type | Cardinality | Required | Description | Examples |
|------------|-----------|-------------|----------|-------------|----------|
| npc_id | string | 1:1 | Yes | Foreign key to compounds | NPC12345 |
| species_id | integer | 1:1 | Yes | Foreign key to species | 301 |
| part | string | 1:1 | No | Plant part if applicable | Leaf |
| reference | string | 1:1 | No | Literature source | PMID:12345 |

---

## Identifiers

| ID Type | Pattern | Example | Description |
|---------|---------|---------|-------------|
| NPC ID | NPC + digits | NPC12345 | Primary compound identifier |
| InChI Key | 27 characters | BLUAFEHZUWYNDE-NNWCWBAJSA-N | Structure hash |
| UniProt ID | Accession | P12345 | Protein target |
| NCBI Taxon ID | Integer | 35608 | Organism identifier |
| PubMed ID | Integer | 12345678 | Literature reference |

---

## Enumerations

### Activity Types

| Type | Description | Typical Unit |
|------|-------------|--------------|
| IC50 | Half-maximal inhibitory concentration | nM, uM |
| EC50 | Half-maximal effective concentration | nM, uM |
| Ki | Inhibition constant | nM, uM |
| Kd | Dissociation constant | nM, uM |
| MIC | Minimum inhibitory concentration | ug/ml |
| LD50 | Lethal dose 50% | mg/kg |

### Activity Relations

| Relation | Meaning |
|----------|---------|
| = | Exactly equal |
| < | Less than |
| > | Greater than |
| ~ | Approximately |
| <= | Less than or equal |
| >= | Greater than or equal |

### Target Types

| Type | Description |
|------|-------------|
| Protein | Generic protein target |
| Enzyme | Catalytic protein |
| Receptor | Cell surface receptor |
| Ion Channel | Ion channel protein |
| Transporter | Membrane transporter |
| Organism | Whole organism (pathogen) |

### Taxonomic Kingdoms

| Kingdom | Description |
|---------|-------------|
| Plantae | Plant-derived |
| Fungi | Fungal metabolites |
| Bacteria | Bacterial products |
| Animalia | Animal-derived |
| Protista | Protist sources |

---

## Entity Relationships

### Compound to Activities
- **Cardinality:** 1:N
- **Description:** One compound can have multiple activity measurements
- **Key Fields:** npc_id, activity_id

### Activity to Target
- **Cardinality:** N:1
- **Description:** Multiple activity measurements per target
- **Key Fields:** activity_id, target_id

### Compound to Species
- **Cardinality:** N:M
- **Description:** Compounds from multiple species; species produce multiple compounds
- **Key Fields:** npc_id, species_id

### Compound-Target-Activity Triplet
- **Cardinality:** Core model
- **Description:** Activity measurement links compound to specific target
- **Key Fields:** npc_id, target_id, activity_type

---

## Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| NPASS | Natural Product Activity and Species Source | Database name |
| NPC | NPASS Compound | Identifier prefix |
| IC50 | Half-maximal Inhibitory Concentration | Activity type |
| EC50 | Half-maximal Effective Concentration | Activity type |
| Ki | Inhibition Constant | Binding affinity |
| Kd | Dissociation Constant | Binding affinity |
| MIC | Minimum Inhibitory Concentration | Antimicrobial |
| SAR | Structure-Activity Relationship | Analysis type |
| SMILES | Simplified Molecular Input Line Entry System | Structure |
| InChI | International Chemical Identifier | IUPAC standard |
| NCBI | National Center for Biotechnology Information | Taxonomy |

---

## Cross-References

| External Database | ID Type | Relationship |
|-------------------|---------|--------------|
| PubChem | InChI Key | Chemical data |
| ChEMBL | InChI Key | Bioactivity |
| UniProt | UniProt ID | Protein targets |
| NCBI Taxonomy | Taxon ID | Organism classification |
| PubMed | PMID | Literature |
| COCONUT | InChI Key | Natural products |

---

## Data Quality Notes

1. **Quantitative Focus:** Emphasizes numeric activity values over qualitative data
2. **Standardized Units:** All activity values converted to standard units
3. **Species Validation:** NCBI Taxonomy linkage for organism verification
4. **Target Annotation:** UniProt cross-references for protein targets
5. **Chemotaxonomy:** Enables analysis of chemical patterns across taxa
6. **Update Frequency:** Annual major releases with continuous additions
