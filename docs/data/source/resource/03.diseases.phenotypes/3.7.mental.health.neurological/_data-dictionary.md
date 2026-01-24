# 3.7 Mental Health/Neurological - Data Dictionary

## Overview

This subcategory contains data about psychiatric and neurological conditions, including brain expression data, GWAS associations for mental health disorders, and synaptic gene annotations.

**Data Sources:** Allen Brain Atlas, PGC (Psychiatric Genomics Consortium), SynGO

---

## Unified Fields

These fields are harmonized across multiple data sources.

| Field Name | Data Type | Cardinality | Description | Sources | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `gene_symbol` | string | Required (1:1) | HGNC gene symbol | Allen Brain Atlas, PGC, SynGO | `SYN1`, `DRD2`, `GRIN2A` |
| `brain_regions` | array[string] | Optional (N:M) | Brain anatomical regions | Allen Brain Atlas, SynGO | `["hippocampus", "prefrontal cortex"]` |
| `disorders` | array[string] | Optional (N:M) | Psychiatric or neurological disorders | PGC | `["schizophrenia", "bipolar disorder"]` |
| `snp_ids` | array[string] | Optional (1:N) | Associated genetic variants. Pattern: `rs[0-9]+` | PGC | `["rs1006737", "rs4648845"]` |

---

## Source-Specific Fields

### Allen Brain Atlas

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `structure_id` | integer | Optional | Allen Brain structure identifier | - | `4020`, `4219`, `4084` |
| `structure_acronym` | string | Optional | Structure abbreviation | - | `HIP`, `CTX`, `TH`, `AMY` |
| `donor_id` | string | Optional | Human brain donor identifier | `H[0-9]+\.[0-9]+` | `H0351.1001`, `H0351.2001` |
| `expression_level` | float | Optional | Gene expression value (log2 normalized) | - | `8.45`, `6.23`, `10.12` |
| `data_type` | enum | Optional | Expression profiling method | `microarray`, `RNA-seq`, `ISH` | `RNA-seq` |

---

### PGC (Psychiatric Genomics Consortium)

| Field Name | Data Type | Cardinality | Description | Example Values |
|------------|-----------|-------------|-------------|----------------|
| `gwas_study_id` | string | Optional | PGC study identifier | `PGC_SCZ3`, `PGC_BIP3`, `PGC_MDD3` |
| `disorder_code` | string | Optional | Psychiatric disorder abbreviation | `SCZ`, `BIP`, `MDD`, `ADHD`, `ASD` |
| `p_value` | float | Optional | GWAS association p-value | `5.2e-12`, `1.3e-8` |
| `odds_ratio` | float | Optional | Effect size (odds ratio) | `1.15`, `0.92` |
| `sample_size` | integer | Optional | Total GWAS sample size | `150000`, `500000` |
| `loci_count` | integer | Optional | Number of genome-wide significant loci | `108`, `287` |

---

### SynGO

| Field Name | Data Type | Cardinality | Description | Pattern | Example Values |
|------------|-----------|-------------|-------------|---------|----------------|
| `syngo_term` | string | Optional | SynGO ontology term | `SYNGO:[0-9]+` | `SYNGO:1234`, `SYNGO:5678` |
| `go_term` | string | Optional | Gene Ontology term | `GO:[0-9]{7}` | `GO:0099055`, `GO:0098794` |
| `uniprot_id` | string | Optional | UniProt protein identifier | `[A-Z0-9]{6,10}` | `P17600`, `Q9UQF2` |
| `synaptic_location` | string | Optional | Synaptic localization | - | `presynapse`, `postsynapse`, `synaptic cleft` |
| `synaptic_function` | string | Optional | Synaptic biological process | - | `neurotransmitter release`, `synaptic plasticity` |
| `evidence_code` | string | Optional | Evidence ontology code | `ECO:[0-9]+` | `ECO:0000314`, `ECO:0000315` |
| `supporting_pmids` | array[integer] | Optional | Supporting publication PMIDs | - | `[12345678, 23456789]` |

---

## Source Field Mappings

### Allen Brain Atlas Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `gene_symbol` | `gene_symbol` |
| `structure_id` | `structure_id` |
| `structure_acronym` | `structure_acronym` |
| `structure_name` | `brain_regions` |
| `donor_id` | `donor_id` |
| `expression_level` | `expression_level` |
| `data_type` | `data_type` |

### PGC Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `gene_symbol` | `gene_symbol` |
| `gwas_study_id` | `gwas_study_id` |
| `disorder` | `disorders` |
| `disorder_code` | `disorder_code` |
| `snp_id` | `snp_ids` |
| `p_value` | `p_value` |
| `odds_ratio` | `odds_ratio` |
| `sample_size` | `sample_size` |
| `loci_count` | `loci_count` |

### SynGO Mappings

| Source Field | Unified Field |
|--------------|---------------|
| `gene_symbol` | `gene_symbol` |
| `syngo_term` | `syngo_term` |
| `go_term` | `go_term` |
| `uniprot_id` | `uniprot_id` |
| `synaptic_location` | `synaptic_location` |
| `synaptic_function` | `synaptic_function` |
| `evidence_code` | `evidence_code` |
| `supporting_pmid` | `supporting_pmids` |
| `brain_region` | `brain_regions` |

---

## PGC Disorder Codes

| Code | Full Name | Description |
|------|-----------|-------------|
| `SCZ` | Schizophrenia | Psychotic disorder |
| `BIP` | Bipolar Disorder | Mood disorder with manic/depressive episodes |
| `MDD` | Major Depressive Disorder | Clinical depression |
| `ADHD` | Attention Deficit Hyperactivity Disorder | Neurodevelopmental disorder |
| `ASD` | Autism Spectrum Disorder | Developmental disorder |
| `AN` | Anorexia Nervosa | Eating disorder |
| `OCD` | Obsessive-Compulsive Disorder | Anxiety-related disorder |
| `PTSD` | Post-Traumatic Stress Disorder | Trauma-related disorder |
| `TS` | Tourette Syndrome | Tic disorder |
| `ED` | Eating Disorders | Combined eating disorder phenotypes |
| `ANX` | Anxiety Disorders | Various anxiety conditions |
| `SUD` | Substance Use Disorders | Addiction phenotypes |

---

## Allen Brain Atlas Structure Acronyms

| Acronym | Full Name | Description |
|---------|-----------|-------------|
| `CTX` | Cerebral Cortex | Outer layer of cerebrum |
| `HIP` | Hippocampus | Memory and spatial navigation |
| `AMY` | Amygdala | Emotion processing |
| `TH` | Thalamus | Sensory relay center |
| `HY` | Hypothalamus | Homeostatic regulation |
| `MB` | Midbrain | Motor control, sensory |
| `Cb` | Cerebellum | Motor coordination |
| `BS` | Brainstem | Vital functions |
| `STR` | Striatum | Movement, reward |
| `GP` | Globus Pallidus | Motor control |
| `SN` | Substantia Nigra | Dopamine production |
| `PFC` | Prefrontal Cortex | Executive function |

---

## SynGO Synaptic Locations

| Location | Description |
|----------|-------------|
| `presynapse` | The transmitting end of the synapse |
| `postsynapse` | The receiving end of the synapse |
| `synaptic cleft` | The gap between pre- and postsynapse |
| `synaptic vesicle` | Vesicles containing neurotransmitters |
| `active zone` | Presynaptic neurotransmitter release site |
| `postsynaptic density` | Protein-dense region of postsynapse |
| `synaptic membrane` | Membrane at the synapse |

---

## Evidence Codes (ECO)

| Code | Name | Description |
|------|------|-------------|
| `ECO:0000314` | Direct assay evidence | Evidence from direct experimental assay |
| `ECO:0000315` | Mutant phenotype evidence | Evidence from mutant phenotype analysis |
| `ECO:0000316` | Genetic interaction evidence | Evidence from genetic interaction studies |
| `ECO:0000353` | Physical interaction evidence | Evidence from protein-protein interaction |
| `ECO:0000269` | Experimental evidence | General experimental evidence |
| `ECO:0007001` | High-throughput evidence | Evidence from high-throughput experiments |

---

## Metadata Fields

| Field Name | Data Type | Description | Example Values |
|------------|-----------|-------------|----------------|
| `_source.database` | string | Name of the source database | `Allen_Brain_Atlas`, `PGC`, `SynGO` |
| `_source.version` | string | Version of the source data | `2024-01`, `SCZ3`, `1.2` |
| `_source.access_date` | date | Date the data was accessed | `2026-01-24` |
| `_source.original_id` | string | Original identifier in source | `SYNGO:1234` |
