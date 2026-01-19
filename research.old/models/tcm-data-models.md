# Traditional Chinese Medicine (TCM) Database Data Models

This document provides comprehensive data models, schemas, and field descriptions for major TCM databases used in systems pharmacology and drug discovery research.

---

## Table of Contents

1. [BATMAN-TCM 2.0](#1-batman-tcm-20)
2. [TCMBank](#2-tcmbank)
3. [HERB 2.0](#3-herb-20)
4. [SymMap](#4-symmap)
5. [TCMSID (Supplementary)](#5-tcmsid-supplementary)
6. [Cross-Database ID Mapping](#6-cross-database-id-mapping)

---

## 1. BATMAN-TCM 2.0

**URL:** http://bionet.ncpsb.org.cn/batman-tcm/
**Publication:** Nucleic Acids Research, 2024
**Primary Function:** Known and predicted TCM ingredient-target protein interactions (TTIs)

### 1.1 Database Statistics

| Entity | Count |
|--------|-------|
| Formulae | 54,832 |
| Herbs | 8,404 |
| Ingredients | 39,171 |
| Known TTIs | 17,068 |
| Predicted TTIs | 2,319,272 |

### 1.2 Download File Formats

Files are available as **tab-delimited text files** from the Download page.

#### 1.2.1 Formula-Herb-Ingredient Associations

Based on TCMID database integration.

| Field | Data Type | Description |
|-------|-----------|-------------|
| formula_id | string | Unique formula identifier |
| formula_name | string | Chinese/Pinyin formula name |
| herb_id | string | Herb identifier |
| herb_name | string | Herb name (Chinese/Latin) |
| ingredient_id | string | Compound identifier |
| ingredient_name | string | Chemical compound name |

#### 1.2.2 Ingredient-Target Interactions (TTI)

| Field | Data Type | Description |
|-------|-----------|-------------|
| ingredient_id | string | Compound identifier (cross-linked to PubChem/TCMID) |
| ingredient_name | string | Chemical compound name |
| target_gene_symbol | string | Gene symbol for target protein |
| target_id | string | Target identifier (cross-linked to GeneBank) |
| score | float | Prediction confidence score (0-1 range) |
| interaction_type | string | "known" or "predicted" |
| evidence_source | string | Source database (DrugBank, KEGG, TTD, or predicted) |
| likelihood_ratio | float | LR value for predicted interactions |

#### 1.2.3 Target Classification (TDL)

Target proteins are classified using the **Target Development Level (TDL)** scheme:

| TDL Category | Description |
|--------------|-------------|
| Tclin | Targets with approved drugs (clinical) |
| Tchem | Targets with small molecule ligands (chemistry) |
| Tbio | Targets with biological annotations |
| Tdark | Dark genome targets with minimal annotation |

### 1.3 API Structure

**Format:** RESTful API returning JSON or HTML

**URL Pattern:**
```
http://bionet.ncpsb.org.cn/batman-tcm/api?type={request_type}&format={output_format}&item={input_item}
```

**Parameters:**
| Parameter | Values | Description |
|-----------|--------|-------------|
| type | ingredient, target, herb, formula | Query type |
| format | json, html | Output format |
| item | identifier/name | Query input |

**JSON Response Structure (Ingredient Query):**
```json
{
  "ingredient_id": "string",
  "ingredient_name": "string",
  "pubchem_cid": "integer",
  "tcmid_id": "string",
  "known_targets": [
    {
      "gene_symbol": "string",
      "gene_id": "string",
      "source": "DrugBank|KEGG|TTD"
    }
  ],
  "predicted_targets": [
    {
      "gene_symbol": "string",
      "score": "float",
      "likelihood_ratio": "float",
      "similarity_features": {
        "compound_similarity": "float",
        "protein_similarity": "float"
      }
    }
  ]
}
```

### 1.4 Supported Gene Identifiers

| Identifier Type | Example Format |
|-----------------|----------------|
| Gene Symbol | SPP1, DNMT3B |
| Entrez Gene ID | 6696 |
| UniProtKB AC/ID | P10451 |
| Ensembl ID | ENSG00000118785 |
| RefSeq Accession | NM_001251830 |
| HGNC ID | HGNC:11255 |
| MIM ID | 166490 |
| IMGT/GENE-DB ID | varies |

### 1.5 Similarity Features for Scoring

| Feature | Description |
|---------|-------------|
| ATC-GO | ATC classification vs GO annotation similarity |
| FP2-closeness | Fingerprint-based structural closeness |
| STRING-sequence | STRING database sequence similarity |
| expression-closeness | Gene expression profile closeness |
| ATC-sequence | ATC vs protein sequence similarity |
| functional_group-sequence | Functional group vs sequence similarity |
| functional_group-GO | Functional group vs GO similarity |
| side_effect-sequence | Side effect profile vs sequence similarity |

---

## 2. TCMBank

**URL:** https://tcmbank.cn/
**Publication:** Chemical Science, 2023
**Primary Function:** Largest herb-ingredient-target-disease mapping

### 2.1 Database Statistics

| Entity | Count | Connected |
|--------|-------|-----------|
| Herbs | 9,192 | 9,010 |
| Ingredients | 61,966 | 54,676 |
| Targets | 15,179 | - |
| Diseases | 32,529 | - |

### 2.2 Download Files

Available formats:
- **XLSX files** - Spreadsheets with detailed entity information
- **Mol2 files** - 3D molecular structures (MM2 force field minimized)

### 2.3 Herb Schema (herbs.xlsx)

| Field | Data Type | Description |
|-------|-----------|-------------|
| TCMBank_ID | string | Unique TCMBank identifier (e.g., TCMB-H00001) |
| TCM_name | string | Chinese herb name |
| TCM_name_en | string | English herb name |
| Herb_pinyin_name | string | Pinyin transliteration |
| Herb_latin_name | string | Latin botanical name |
| Properties | string | TCM properties (cold/hot, etc.) |
| Meridian_tropism | string | Target meridians |
| Function | text | Therapeutic functions |
| Indication | text | Clinical indications |
| Therapeutic_class | string | Classification category |

### 2.4 Ingredient Schema (ingredients.xlsx)

| Field | Data Type | Description |
|-------|-----------|-------------|
| TCMBank_ID | string | Unique identifier (e.g., TCMB-I00001) |
| Ingredient_Name | string | Chemical compound name |
| Ingredient_Alias | string | Alternative names |
| SMILES | string | Canonical SMILES notation |
| Molecular_Formula | string | Chemical formula (e.g., C15H10O7) |
| Molecular_Weight | float | Molecular weight (Da) |
| Molecular_Volume | float | Molecular volume |
| mol2_path | string | Path to mol2 structure file |
| OB_score | float | Oral bioavailability score |
| ADMET | object | ADMET properties object |
| LogD | float | Distribution coefficient |
| AlogP | float | Lipophilicity |
| Solubility | float | Aqueous solubility |
| Polar_Surface_Area | float | TPSA value |

### 2.5 Target Schema (targets.xlsx)

| Field | Data Type | Description |
|-------|-----------|-------------|
| TCMBank_ID | string | Unique identifier (e.g., TCMB-T00001) |
| Gene_Name | string | Official gene name |
| Gene_Alias | string | Alternative gene names |
| Gene_Symbol | string | HGNC gene symbol |
| Description | text | Functional description |
| Chromosome | string | Chromosomal location |
| UniProt_ID | string | UniProt accession |
| HGNC_ID | string | HGNC identifier |
| OMIM_ID | string | OMIM entry number |

### 2.6 Disease Schema (diseases.xlsx)

| Field | Data Type | Description |
|-------|-----------|-------------|
| TCMBank_ID | string | Unique identifier (e.g., TCMB-D00001) |
| Disease_Name | string | Disease name |
| Disease_Type | string | Disease category |
| DisGeNET_Link | string | DisGeNET URL |
| MeSH_Link | string | MeSH URL |
| MeSH_Name | string | MeSH preferred term |
| DO_Name | string | Disease Ontology name |
| HPO_ID | string | Human Phenotype Ontology ID |

### 2.7 Relationship Tables

#### Herb-Ingredient Associations
| Field | Data Type | Description |
|-------|-----------|-------------|
| Herb_ID | string | TCMBank herb ID |
| Ingredient_ID | string | TCMBank ingredient ID |
| Evidence_Source | string | Data source |

#### Ingredient-Target Associations
| Field | Data Type | Description |
|-------|-----------|-------------|
| Ingredient_ID | string | TCMBank ingredient ID |
| Target_ID | string | TCMBank target ID |
| Interaction_Type | string | Known/Predicted |
| Confidence | float | Confidence score |

#### Target-Disease Associations
| Field | Data Type | Description |
|-------|-----------|-------------|
| Target_ID | string | TCMBank target ID |
| Disease_ID | string | TCMBank disease ID |
| Evidence_Level | string | Association evidence |
| Source | string | Data source (DisGeNET, OMIM, etc.) |

### 2.8 ID Format Patterns

| Entity | Pattern | Example |
|--------|---------|---------|
| Herb | TCMB-H##### | TCMB-H00001 |
| Ingredient | TCMB-I##### | TCMB-I00001 |
| Target | TCMB-T##### | TCMB-T00001 |
| Disease | TCMB-D##### | TCMB-D00001 |

---

## 3. HERB 2.0

**URL:** http://herb.ac.cn/
**Publication:** Nucleic Acids Research, 2025
**Primary Function:** High-throughput experiment and reference-guided TCM database with clinical evidence

### 3.1 Database Statistics

| Entity | Count |
|--------|-------|
| Herbs | 6,892 |
| Ingredients | 44,595 |
| Formulae | 6,743 |
| Gene Targets | 15,515 |
| Diseases | 30,170 |
| Clinical Trials | 8,558 |
| Meta-analyses | 8,032 |
| High-throughput Experiments | 2,231 |
| Curated References | 6,644 |

### 3.2 Knowledge Graph Structure

**Nine Entity Types:**
1. Herbs
2. Ingredients
3. Formulae
4. Targets
5. Diseases
6. Clinical Trials
7. Meta-analyses
8. High-throughput Experiments
9. Curated References

**28 Relationship Types** connecting entities with weighted edges (0-1)

### 3.3 Clinical Trial Data Schema

| Field | Data Type | Description |
|-------|-----------|-------------|
| Trial_ID | string | ClinicalTrials.gov identifier (NCT#) |
| Title | string | Trial title |
| Status | string | Current status |
| Phase | string | Phase 1, 2, 3, or 4 |
| Blinding | string | Open, Single, Double, Triple |
| Randomization | boolean | Randomized trial flag |
| TCM_Entity | string | Herb/Ingredient/Formula studied |
| Disease | string | Target disease |
| Conclusions | text | Clinical conclusions (for 1,941 trials) |
| Source | string | ClinicalTrials.gov |

### 3.4 Meta-analysis Data Schema

| Field | Data Type | Description |
|-------|-----------|-------------|
| MA_ID | string | PROSPERO identifier |
| Title | string | Review title |
| TCM_Entity | string | Herb/Ingredient/Formula reviewed |
| Disease | string | Target disease |
| Conclusions | text | Review conclusions (for 593 analyses) |
| Supporting_Papers | array | References supporting conclusions |
| Source | string | PROSPERO |

### 3.5 Gene Expression Data Format

#### Experiment Metadata
| Field | Data Type | Description |
|-------|-----------|-------------|
| EXP_ID | string | HERB experiment identifier |
| GSE_ID | string | GEO Series identifier |
| GSM_IDs | array | GEO Sample identifiers |
| Platform | string | Microarray/RNA-seq |
| Species | string | Homo sapiens / Mus musculus |
| TCM_Entity | string | Herb/Ingredient studied |
| Cell_Line | string | Cell line or tissue type |
| Treatment_Time | string | Duration of treatment |

#### Differential Expression Results
| Field | Data Type | Description |
|-------|-----------|-------------|
| Gene_Symbol | string | Official gene symbol |
| log2FC | float | Log2 fold change |
| P_value | float | Statistical p-value |
| FDR | float | False discovery rate |
| Direction | string | Up/Down regulated |

**Thresholds:** |log2(fold change)| >= 0.5, P <= 0.05

#### Gene Set Enrichment
| Field | Data Type | Description |
|-------|-----------|-------------|
| Term_ID | string | GO/KEGG term identifier |
| Term_Name | string | Term description |
| P_value | float | Enrichment p-value |
| FDR | float | Adjusted p-value |
| Gene_Count | integer | Genes in term |
| Gene_List | array | Enriched gene symbols |

### 3.6 Connectivity Mapping Interface

**Input Format (Gene Expression Signature):**
```
# Upregulated genes (ranked by significance)
GENE1
GENE2
GENE3
...

# Downregulated genes (ranked by significance)
GENEA
GENEB
GENEC
...
```

**Output:** Similarity scores to curated herb/ingredient profiles

### 3.7 Target Standardization

| Source Database | Cross-reference |
|-----------------|-----------------|
| GeneBank | Primary identifier |
| GeneCards | Alias resolution |
| TTD | Drug target validation |
| UniProt | Protein sequences |

### 3.8 Disease Standardization

| Source Database | Purpose |
|-----------------|---------|
| DisGeNET | Primary standardization |
| OMIM | Genetic disease mapping |
| HPO | Phenotype annotation |
| Disease Ontology | Hierarchical classification |

---

## 4. SymMap

**URL:** http://www.symmap.org/
**Publication:** Nucleic Acids Research, 2019
**Primary Function:** TCM-Modern medicine integration through symptom mapping

### 4.1 Database Statistics

| Entity | Count |
|--------|-------|
| TCM Symptoms | 1,717 |
| MM Symptoms | 961 |
| Herbs | 499 |
| Ingredients | 19,595 |
| Targets | 4,302 |
| Diseases | 5,235 |
| Total Nodes | 32,281 |
| Total Edges | 403,318 |

### 4.2 Network Architecture

**Direct Associations:** 106,721 edges
**Indirect Associations:** 296,597 edges (statistically inferred)

### 4.3 Six Direct Association Types

| Relationship | Count | Avg per Entity |
|--------------|-------|----------------|
| Herb - TCM Symptom | 6,638 | 13.30/herb, 3.87/symptom |
| TCM Symptom - MM Symptom | 2,978 | 1.74/TCM, 3.13/MM |
| Herb - Ingredient | 48,372 | - |
| MM Symptom - Disease | 12,107 | - |
| Ingredient - Target | 29,370 | - |
| Target - Disease | 7,256 | - |

### 4.4 Nine Indirect Association Types

Derived through statistical inference (Fisher's Exact Test):
- Herb - MM Symptom
- Herb - Disease
- Herb - Target
- TCM Symptom - Disease
- TCM Symptom - Ingredient
- TCM Symptom - Target
- Ingredient - Disease
- Ingredient - MM Symptom
- Target - MM Symptom

### 4.5 Entity Schemas

#### TCM Symptom
| Field | Data Type | Description |
|-------|-----------|-------------|
| Symptom_ID | string | SymMap symptom identifier |
| TCM_Term | string | Chinese symptom term |
| TCM_Term_Pinyin | string | Pinyin transliteration |
| Category | string | Symptom category |
| Body_Region | string | Associated body region |

#### MM Symptom (Modern Medicine)
| Field | Data Type | Description |
|-------|-----------|-------------|
| MM_Symptom_ID | string | SymMap MM symptom identifier |
| UMLS_ID | string | UMLS Concept Unique Identifier |
| HPO_ID | string | Human Phenotype Ontology ID |
| Symptom_Name | string | Standard symptom name |
| Definition | text | Symptom definition |

#### Herb
| Field | Data Type | Description |
|-------|-----------|-------------|
| Herb_ID | string | SymMap herb identifier |
| Chinese_Name | string | Chinese herb name |
| Pinyin_Name | string | Pinyin transliteration |
| Latin_Name | string | Botanical name |
| CP_Edition | string | Chinese Pharmacopoeia edition |

#### Ingredient
| Field | Data Type | Description |
|-------|-----------|-------------|
| Ingredient_ID | string | SymMap ingredient identifier |
| CAS | string | CAS Registry Number |
| PubChem_CID | integer | PubChem Compound ID |
| InChI_Key | string | InChI Key identifier |
| SMILES | string | SMILES notation |
| Source_DB | string | TCMID or TCM-ID |

#### Target
| Field | Data Type | Description |
|-------|-----------|-------------|
| Target_ID | string | SymMap target identifier |
| Gene_Symbol | string | HGNC gene symbol |
| UniProt_ID | string | UniProt accession |
| Source | string | HIT, HPO, DrugBank, NCBI |

#### Disease
| Field | Data Type | Description |
|-------|-----------|-------------|
| Disease_ID | string | SymMap disease identifier |
| OMIM_ID | string | OMIM entry number |
| Orphanet_ID | string | Orphanet identifier |
| Disease_Name | string | Disease name |
| Category | string | Disease category |

### 4.6 Download Data Formats

**Three dataset options per relationship:**
1. **Full Set** - All associations
2. **Loose Selection** - P-value < 0.05
3. **Stringent Selection** - FDR (Bonferroni & BH) < 0.05

### 4.7 Association Confidence Fields

| Field | Data Type | Description |
|-------|-----------|-------------|
| Entity1_ID | string | First entity identifier |
| Entity2_ID | string | Second entity identifier |
| Association_Type | string | Direct/Indirect |
| P_value | float | Fisher's exact test p-value |
| FDR_BH | float | Benjamini-Hochberg adjusted p-value |
| FDR_Bonferroni | float | Bonferroni adjusted p-value |
| Evidence_Count | integer | Supporting evidence instances |

### 4.8 Data Source Integration

| Component | Source Databases |
|-----------|------------------|
| TCM Knowledge | Chinese Pharmacopoeia (2015) |
| MM Symptoms | UMLS, MeSH, SIDER, HPO |
| Drugs | DrugBank |
| Diseases | OMIM, Orphanet |
| Ingredients | TCMID, TCM-ID |
| Targets | HIT, HPO, DrugBank, NCBI Gene |

### 4.9 Technical Infrastructure

- **Backend:** MySQL
- **Web Framework:** Python-Flask
- **Web Server:** Nginx
- **Network Visualization:** Interactive web interface

---

## 5. TCMSID (Supplementary)

**URL:** https://tcmsid.pku.edu.cn/
**Publication:** Journal of Cheminformatics, 2022
**Primary Function:** Simplified integrated database for TCM drug discovery

### 5.1 Database Statistics

| Entity | Count |
|--------|-------|
| TCM Herbs | 499 |
| Unique Ingredients | 20,015 |
| Target Proteins | 3,270 |
| Drugs | 10,487 (3,883 FDA-approved) |
| Herb-Ingredient Pairs | 50,053 |

### 5.2 Ingredient Quality Metrics

| Field | Data Type | Values | Description |
|-------|-----------|--------|-------------|
| Significance_Degree | integer | 0, 1, 2 | Herb-ingredient relationship confidence |
| Structural_Reliability | integer | 1-5 | Structure accuracy score |
| ClassyFire_Kingdom | string | varies | Chemical taxonomy kingdom |
| ClassyFire_Superclass | string | varies | Chemical taxonomy superclass |
| ClassyFire_Class | string | varies | Chemical taxonomy class |

### 5.3 Enhanced ADME/T Properties

| Property | Data Type | Description |
|----------|-----------|-------------|
| Caco2_Permeability | float | Intestinal absorption |
| F30_Bioavailability | float | 30% bioavailability threshold |
| PPB | float | Plasma protein binding |
| BBB_Penetration | float | Blood-brain barrier |
| T_half | float | Half-life |
| CL | float | Clearance rate |
| hERG_Inhibition | boolean | Cardiac toxicity risk |
| Hepatotoxicity | boolean | Liver toxicity risk |
| Drug_Likeness | float | DL score |
| LogP | float | Partition coefficient |
| LogS | float | Aqueous solubility |

### 5.4 Target Development Level

| TDL | Description |
|-----|-------------|
| Tclin | FDA-approved drug targets |
| Tchem | Chemical probe targets |
| Tbio | Biology-annotated targets |
| Tdark | Understudied targets |

### 5.5 Similarity Metrics

| Fingerprint | Algorithm | Description |
|-------------|-----------|-------------|
| FCFP6 | Morgan/Circular | Functional connectivity fingerprint (radius 3) |
| ECFP4 | Morgan/Circular | Extended connectivity fingerprint (radius 2) |
| Tanimoto | Coefficient | Similarity score (0-1) |

---

## 6. Cross-Database ID Mapping

### 6.1 Compound Identifiers

| Database | Primary ID | Cross-references |
|----------|------------|------------------|
| BATMAN-TCM | Internal ID | PubChem, TCMID |
| TCMBank | TCMB-I##### | SMILES, PubChem |
| HERB | Internal ID | PubChem, CAS |
| SymMap | Internal ID | CAS, PubChem, InChIKey |

### 6.2 Target Identifiers

| Database | Primary ID | Cross-references |
|----------|------------|------------------|
| BATMAN-TCM | GeneBank ID | UniProt, HGNC, Ensembl, Entrez |
| TCMBank | TCMB-T##### | UniProt, HGNC, OMIM |
| HERB | GeneBank | GeneCards, TTD, UniProt |
| SymMap | Internal ID | UniProt, Gene Symbol |

### 6.3 Disease Identifiers

| Database | Primary Sources |
|----------|-----------------|
| BATMAN-TCM | MeSH, OMIM |
| TCMBank | DisGeNET, MeSH, DO, HPO |
| HERB | DisGeNET, OMIM, HPO, DO |
| SymMap | OMIM, Orphanet |

### 6.4 Recommended ID Mapping Strategy

```python
# Example cross-database compound mapping
compound_mapping = {
    "pubchem_cid": "primary",  # Best for cross-database
    "inchikey": "secondary",    # Structure-based mapping
    "cas": "tertiary",          # Legacy support
    "smiles": "structure"       # Structural comparison
}

# Example target mapping
target_mapping = {
    "gene_symbol": "primary",   # Human-readable
    "uniprot_id": "secondary",  # Protein-level
    "entrez_id": "tertiary",    # NCBI standard
    "ensembl_id": "genomic"     # Genomic context
}
```

---

## References

1. **BATMAN-TCM 2.0:** Kong X, et al. (2024) Nucleic Acids Research, 52(D1):D1110-D1117. https://pubmed.ncbi.nlm.nih.gov/37904598/

2. **TCMBank:** Xu HY, et al. (2023) Chemical Science, 14(39):10684-10700. https://pubs.rsc.org/en/content/articlehtml/2023/sc/d3sc02139d

3. **HERB 2.0:** Fang S, et al. (2025) Nucleic Acids Research, 53(D1):D1404-D1416. https://pubmed.ncbi.nlm.nih.gov/39558177/

4. **SymMap:** Wu Y, et al. (2019) Nucleic Acids Research, 47(D1):D1110-D1117. https://pubmed.ncbi.nlm.nih.gov/30380087/

5. **TCMSID:** Xue R, et al. (2022) Journal of Cheminformatics, 14(1):89. https://pmc.ncbi.nlm.nih.gov/articles/PMC9805110/

---

## Appendix: Example Data Rows

### A1. TCMBank Herb Example

| TCMBank_ID | TCM_name | TCM_name_en | Herb_pinyin_name | Function |
|------------|----------|-------------|------------------|----------|
| TCMB-H00001 | 人参 | Ginseng | Ren Shen | Tonify Qi, generate fluids |

### A2. BATMAN-TCM Predicted TTI Example

| Ingredient | Target | Score | LR | Source |
|------------|--------|-------|-----|--------|
| Theophylline | SPP1 | 0.82 | 4.5 | predicted |
| Vitamin E | MMP9 | 0.70 | 3.2 | predicted |

### A3. SymMap Association Example

| TCM_Symptom | MM_Symptom | P_value | FDR_BH |
|-------------|------------|---------|--------|
| 心悸 (Palpitation) | Palpitations | 0.001 | 0.01 |

---

*Document generated: January 2026*
*Last updated: January 18, 2026*
