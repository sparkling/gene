# 6.4 Metabolomics - Data Dictionary

## Overview

This data dictionary documents all fields for the Metabolomics subcategory, integrating data from Exposome-Explorer and HMDB (Human Metabolome Database).

| Attribute | Value |
|-----------|-------|
| Subcategory ID | 6.4 |
| Subcategory Name | Metabolomics |
| Data Sources | Exposome-Explorer, HMDB |
| Schema ID | https://gene.taxonomy/schemas/6.4-metabolomics |

---

## Unified Fields

These fields have been harmonized across all data sources in this subcategory.

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| metabolite_id | string | Required (1:1) | Unique identifier for the metabolite | `123`, `HMDB0000001` |
| metabolite_name | string | Required (1:1) | Name of the metabolite or biomarker | `4-O-methylgallic acid`, `1-Methylhistidine` |
| pubchem_cid | integer | Optional (1:1) | PubChem compound identifier | `13159`, `92105` |
| hmdb_id | string | Optional (1:1) | Human Metabolome Database identifier | `HMDB0002122`, `HMDB0000001` |
| chebi_id | string | Optional (1:1) | ChEBI chemical entities identifier | `CHEBI:50599`, `CHEBI:16243` |
| cas_number | string | Optional (1:1) | Chemical Abstracts Service registry number | `3934-84-7`, `332-80-9` |
| molecular_formula | string | Optional (1:1) | Chemical formula of the metabolite | `C8H8O5`, `C7H11N3O2` |
| smiles | string | Optional (1:1) | SMILES notation for molecular structure | `CN1C=NC(CC(N)C(O)=O)=C1` |
| inchi | string | Optional (1:1) | International Chemical Identifier | `InChI=1S/C7H11N3O2/c1-10-4-9-5...` |
| inchikey | string | Optional (1:1) | Hashed InChI key for structure lookup | `BRMWTNUJHUMWMS-LURJTMIESA-N` |
| specimen_types | array[string] | Optional (1:N) | Biological sample types where metabolite is measured | `["plasma", "urine", "serum"]` |
| concentrations | array[object] | Optional (1:N) | Measured concentrations in biological samples | See structure below |
| biomarker_type | string | Optional (1:1) | Type of biomarker | `metabolite`, `nutrient`, `Phenolic acid metabolite` |

### Concentrations Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| specimen_type | string | Type of biological sample | `Blood`, `Urine` |
| value | number | Measured concentration | `5.2` |
| unit | string | Unit of measurement | `uM`, `umol/L`, `ng/mL` |
| condition | string | Health condition | `Normal`, `Diabetes` |
| age | string | Age group of subjects | `Adult (>18 years)` |
| sex | string | Sex of study subjects | `Male`, `Female`, `Both` |
| subject_count | integer | Number of subjects | `120` |

---

## Source-Specific Fields

### Exposome-Explorer

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| biomarker_id | integer | Optional (1:1) | Unique identifier for biomarker | `123` |
| exposure_id | integer | Optional (1:1) | Unique identifier for exposure | `789` |
| exposure_name | string | Optional (1:1) | Name of dietary or environmental exposure | `Tea consumption`, `Coffee consumption` |
| exposure_category | string | Optional (1:1) | Category of exposure | `Dietary`, `Environmental` |
| exposure_subcategory | string | Optional (1:1) | Subcategory of exposure | `Beverages`, `Fruits`, `Pollutants` |
| association_type | string | Optional (1:1) | Direction of biomarker-exposure relationship | `positive`, `negative`, `correlation` |
| effect_size | decimal | Optional (1:1) | Magnitude of association | `0.42` |
| effect_unit | string | Optional (1:1) | Type of effect size metric | `Pearson correlation`, `Odds ratio` |
| p_value | decimal | Optional (1:1) | Statistical significance | `0.00001` |
| confidence_interval | string | Optional (1:1) | 95% confidence interval for effect size | `0.35-0.49` |
| sample_size | integer | Optional (1:1) | Number of subjects in study | `450` |
| population | string | Optional (1:1) | Study population description | `European adults`, `Free-living adults, US` |
| study_type | string | Optional (1:1) | Type of epidemiological study | `cross-sectional`, `cohort`, `intervention` |

### HMDB

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| accession | string | Optional (1:1) | HMDB accession number | `HMDB0000001` |
| version | integer | Optional (1:1) | Database version for the metabolite record | `5` |
| iupac_name | string | Optional (1:1) | IUPAC systematic name | `(2S)-2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid` |
| average_molecular_weight | decimal | Optional (1:1) | Average molecular weight in g/mol | `169.181` |
| monoisotopic_molecular_weight | decimal | Optional (1:1) | Monoisotopic mass for MS identification | `169.085126611` |
| state | string | Optional (1:1) | Physical state of the compound | `Solid`, `Liquid`, `Gas` |
| biofluid_locations | array[string] | Optional (1:N) | Biofluids where metabolite is found | `["Blood", "Urine", "CSF"]` |
| cellular_locations | array[string] | Optional (1:N) | Cellular compartments where metabolite is present | `["Cytoplasm", "Mitochondria"]` |
| tissue_locations | array[string] | Optional (1:N) | Tissues where metabolite is found | `["Muscle", "Liver"]` |
| disease_name | string | Optional (1:N) | Disease associated with metabolite levels | `Histidinemia`, `Diabetes` |
| omim_id | string | Optional (1:1) | Online Mendelian Inheritance in Man identifier | `235800` |
| concentration_change | string | Optional (1:1) | Direction of change in disease state | `Elevated`, `Reduced` |
| pathway_name | string | Optional (1:N) | Metabolic pathway involving the metabolite | `Histidine metabolism`, `Glycolysis` |
| kegg_id | string | Optional (1:1) | KEGG compound identifier | `C01152` |
| smpdb_id | string | Optional (1:1) | Small Molecule Pathway Database ID | `SMP00001` |
| drugbank_id | string | Optional (1:1) | DrugBank database identifier | `DB00118` |
| foodb_id | string | Optional (1:1) | FooDB database identifier for food connection | `FDB000001` |
| ms_spectrum | object | Optional (1:N) | Mass spectrometry spectral data | See structure below |
| nmr_spectrum | object | Optional (1:N) | NMR spectral data for identification | See structure below |
| taxonomy | object | Optional (1:1) | Chemical classification taxonomy | See structure below |
| condition | string | Optional (1:1) | Health condition | `Normal`, `Diabetes`, `Cancer` |
| age | string | Optional (1:1) | Age group of subjects | `Adult (>18 years)`, `Pediatric` |
| sex | string | Optional (1:1) | Sex of study subjects | `Male`, `Female`, `Both` |
| subject_count | integer | Optional (1:1) | Number of subjects in concentration study | `120` |
| pubmed_id | integer | Optional (1:1) | PubMed reference ID for source data | `12345678` |

#### MS Spectrum Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| spectrum_type | string | Type of MS analysis | `MS/MS`, `GC-MS` |
| ionization_mode | string | Ionization polarity | `Positive`, `Negative` |
| peaks | array | Mass/intensity pairs | `[[100.0, 1000], [150.5, 500]]` |

#### NMR Spectrum Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| nucleus | string | NMR active nucleus | `1H`, `13C` |
| frequency | string | Spectrometer frequency | `400 MHz`, `600 MHz` |

#### Taxonomy Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| kingdom | string | Chemical kingdom | `Organic compounds` |
| superclass | string | Chemical superclass | `Organic acids` |
| class | string | Chemical class | `Amino acids` |
| subclass | string | Chemical subclass | `Alpha amino acids` |

---

## Unified Complex Fields

### Exposures Array (Exposome-Explorer)

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| exposures | array[object] | Optional (1:N) | Dietary or environmental exposure associations | See structure below |

#### Exposure Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| exposure_id | integer | Unique exposure identifier | `789` |
| exposure_name | string | Name of exposure | `Tea consumption` |
| exposure_category | string | Category | `Dietary` |
| exposure_subcategory | string | Subcategory | `Beverages` |
| association_type | string | Relationship direction | `positive` |
| effect_size | number | Association magnitude | `0.42` |
| effect_unit | string | Effect measure type | `Pearson correlation` |
| p_value | number | Statistical significance | `0.00001` |
| confidence_interval | string | 95% CI | `0.35-0.49` |
| sample_size | integer | Number of subjects | `450` |
| population | string | Study population | `European adults` |
| study_type | string | Study design | `cross-sectional` |

### Disease Associations (HMDB)

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| disease_associations | array[object] | Optional (1:N) | Diseases associated with metabolite levels | See structure below |

#### Disease Association Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| disease_name | string | Name of disease | `Histidinemia` |
| omim_id | string | OMIM identifier | `235800` |
| concentration_change | string | Direction of change | `Elevated`, `Reduced` |

### Pathways (HMDB)

| Field Name | Data Type | Cardinality | Description | Example |
|------------|-----------|-------------|-------------|---------|
| pathways | array[object] | Optional (1:N) | Metabolic pathways involving the metabolite | See structure below |

#### Pathway Object Structure

| Property | Data Type | Description | Example |
|----------|-----------|-------------|---------|
| pathway_name | string | Name of pathway | `Histidine metabolism` |
| kegg_id | string | KEGG identifier | `C01152` |
| smpdb_id | string | SMPDB identifier | `SMP00001` |

---

## Field Mappings

### Exposome-Explorer to Unified Schema

| Source Field | Unified Field |
|--------------|---------------|
| biomarker_id | metabolite_id |
| biomarker_name | metabolite_name |
| biomarker_type | biomarker_type |
| hmdb_id | hmdb_id |
| pubchem_cid | pubchem_cid |
| chebi_id | chebi_id |
| cas_number | cas_number |
| molecular_formula | molecular_formula |
| specimen_type | specimen_types |
| mean_concentration | concentrations[].value |
| concentration_unit | concentrations[].unit |
| exposure_id | exposures[].exposure_id |
| exposure_name | exposures[].exposure_name |
| exposure_category | exposures[].exposure_category |
| exposure_subcategory | exposures[].exposure_subcategory |
| association_type | exposures[].association_type |
| effect_size | exposures[].effect_size |
| effect_unit | exposures[].effect_unit |
| p_value | exposures[].p_value |
| confidence_interval | exposures[].confidence_interval |
| sample_size | exposures[].sample_size |
| population | exposures[].population |
| study_type | exposures[].study_type |

### HMDB to Unified Schema

| Source Field | Unified Field |
|--------------|---------------|
| accession | metabolite_id |
| name | metabolite_name |
| version | version |
| iupac_name | iupac_name |
| chemical_formula | molecular_formula |
| average_molecular_weight | average_molecular_weight |
| monoisotopic_molecular_weight | monoisotopic_molecular_weight |
| smiles | smiles |
| inchi | inchi |
| inchikey | inchikey |
| cas_registry_number | cas_number |
| state | state |
| biofluid | specimen_types |
| biofluid_locations | biofluid_locations |
| cellular_locations | cellular_locations |
| tissue_locations | tissue_locations |
| concentration_value | concentrations[].value |
| concentration_units | concentrations[].unit |
| condition | concentrations[].condition |
| age | concentrations[].age |
| sex | concentrations[].sex |
| subject_count | concentrations[].subject_count |
| disease_name | disease_associations[].disease_name |
| omim_id | disease_associations[].omim_id |
| concentration_change | disease_associations[].concentration_change |
| pathway_name | pathways[].pathway_name |
| kegg_id | pathways[].kegg_id |
| smpdb_id | pathways[].smpdb_id |
| drugbank_id | drugbank_id |
| foodb_id | foodb_id |
| taxonomy | taxonomy |
| ms_spectrum | ms_spectra |
| nmr_spectrum | nmr_spectra |
| pubmed_id | pubmed_id |

---

## Cardinality Legend

| Symbol | Meaning |
|--------|---------|
| 1:1 | Exactly one value per record |
| 1:N | One or more values per record |
| Required | Field must be present |
| Optional | Field may be null or absent |

---

## Cross-Reference Databases

| Identifier | Database | URL Pattern |
|------------|----------|-------------|
| hmdb_id | Human Metabolome Database | https://hmdb.ca/metabolites/{id} |
| pubchem_cid | PubChem | https://pubchem.ncbi.nlm.nih.gov/compound/{id} |
| chebi_id | ChEBI | https://www.ebi.ac.uk/chebi/searchId.do?chebiId={id} |
| kegg_id | KEGG Compound | https://www.genome.jp/dbget-bin/www_bget?cpd:{id} |
| drugbank_id | DrugBank | https://go.drugbank.com/drugs/{id} |
| foodb_id | FooDB | https://foodb.ca/compounds/{id} |
| omim_id | OMIM | https://omim.org/entry/{id} |
| smpdb_id | SMPDB | https://smpdb.ca/view/{id} |

---

## Study Types Reference

| Study Type | Description |
|------------|-------------|
| cross-sectional | Single time-point measurement across population |
| cohort | Longitudinal follow-up of defined group |
| case-control | Comparison of cases vs controls |
| intervention | Controlled experimental manipulation |
| randomized controlled trial | Randomized experimental design |

---

## Data Quality Notes

1. **metabolite_id** - HMDB uses accession format (HMDB0000001); Exposome-Explorer uses numeric IDs
2. **hmdb_id** - Primary cross-reference identifier; links to comprehensive metabolite information
3. **foodb_id** - Links metabolites to food composition data (FooDB)
4. **exposures** - Exposome-Explorer specific; links biomarkers to dietary/environmental factors
5. **disease_associations** - HMDB specific; documents metabolite-disease relationships
6. **ms_spectra/nmr_spectra** - Experimental spectra for compound identification
7. **monoisotopic_molecular_weight** - More precise than average; used for MS identification
8. **concentration values** - Units vary; always check `unit` field for interpretation
9. **taxonomy** - Chemical classification from ClassyFire ontology
