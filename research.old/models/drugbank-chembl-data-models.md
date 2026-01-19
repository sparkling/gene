# Drug Database Data Models: DrugBank, ChEMBL, PubChem, and DrugCentral

## Overview

This document provides comprehensive data model documentation for four major drug databases used in pharmaceutical research, drug discovery, and pharmacological analysis. These databases serve as foundational resources for understanding drug-target interactions, bioactivity measurements, and clinical pharmacology.

---

## 1. DrugBank

### 1.1 Overview

DrugBank is a comprehensive, freely accessible online database containing information on drugs and drug targets. As both a bioinformatics and cheminformatics resource, DrugBank combines detailed drug (chemical, pharmacological, and pharmaceutical) data with comprehensive drug target (sequence, structure, and pathway) information.

**Current Version:** DrugBank 6.0 (2024)
**License:** CC BY-NC 4.0
**Website:** https://go.drugbank.com
**Documentation:** https://docs.drugbank.com/xml/

### 1.2 XML Schema Structure

#### 1.2.1 Root Element

```xml
<drugbank xmlns="http://www.drugbank.ca"
          version="5.1"
          exported-on="2024-01-01"
          schemaLocation="http://www.drugbank.ca/docs/drugbank.xsd">
```

| Attribute | Type | Description |
|-----------|------|-------------|
| `version` | string | Schema version (current: 5.1) |
| `exported-on` | date | Export date of the database |

#### 1.2.2 Drug Element Structure

```xml
<drug type="small molecule" created="2005-06-13" updated="2024-01-01">
  <drugbank-id primary="true">DB00001</drugbank-id>
  <drugbank-id>BTD00024</drugbank-id>
  <name>Lepirudin</name>
  <description>...</description>
  <cas-number>120993-53-5</cas-number>
  <unii>Y43GF64R34</unii>
  <state>liquid</state>
  <groups>
    <group>approved</group>
  </groups>
  ...
</drug>
```

**Drug Element Attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `type` | enum | "small molecule", "biotech" |
| `created` | date | Record creation date |
| `updated` | date | Last update date |

**Drug Groups (Approval Status):**

| Value | Description |
|-------|-------------|
| `approved` | Approved for clinical use |
| `illicit` | Controlled/illegal substance |
| `experimental` | In experimental stage |
| `withdrawn` | Withdrawn from market |
| `nutraceutical` | Dietary supplement |
| `investigational` | Under investigation |
| `vet_approved` | Approved for veterinary use |

#### 1.2.3 Core Drug Data Elements

| Element | Type | Description |
|---------|------|-------------|
| `drugbank-id` | string | Primary identifier (e.g., DB00001) |
| `name` | string | Drug name |
| `description` | text | Detailed description |
| `cas-number` | string | CAS Registry Number |
| `unii` | string | FDA Unique Ingredient Identifier |
| `average-mass` | decimal | Average molecular mass |
| `monoisotopic-mass` | decimal | Monoisotopic mass |
| `state` | enum | solid, liquid, gas |
| `synthesis-reference` | text | Synthesis method reference |
| `indication` | text | Therapeutic indication |
| `pharmacodynamics` | text | Pharmacodynamic properties |
| `mechanism-of-action` | text | Mechanism of action description |
| `toxicity` | text | Toxicity information |
| `metabolism` | text | Metabolic pathways |
| `absorption` | text | Absorption characteristics |
| `half-life` | text | Half-life information |
| `protein-binding` | text | Protein binding percentage |
| `route-of-elimination` | text | Elimination routes |
| `volume-of-distribution` | text | Volume of distribution |
| `clearance` | text | Clearance rate |
| `fda-label` | URL | FDA label link |
| `msds` | URL | Material Safety Data Sheet link |

### 1.3 Targets, Enzymes, Carriers, Transporters

DrugBank distinguishes four types of protein interactions:

| Type | Element | Description |
|------|---------|-------------|
| **Targets** | `<targets>` | Disease-related proteins that the drug acts on |
| **Enzymes** | `<enzymes>` | Proteins that metabolize the drug |
| **Carriers** | `<carriers>` | Proteins that move drugs around the body |
| **Transporters** | `<transporters>` | Proteins that move drugs into/around cells |

#### 1.3.1 Common Structure for Interactants

```xml
<target position="1">
  <id>BE0000048</id>
  <name>Prothrombin</name>
  <organism>Humans</organism>
  <actions>
    <action>inhibitor</action>
  </actions>
  <references>
    <articles>
      <article>
        <pubmed-id>10505536</pubmed-id>
        <citation>...</citation>
      </article>
    </articles>
  </references>
  <known-action>yes</known-action>
  <polypeptide id="P00734" source="Swiss-Prot">
    <name>Prothrombin</name>
    <general-function>...</general-function>
    <specific-function>...</specific-function>
    <gene-name>F2</gene-name>
    <locus>11p11.2</locus>
    <cellular-location>Secreted</cellular-location>
    <transmembrane-regions></transmembrane-regions>
    <signal-regions>1-24</signal-regions>
    <theoretical-pi>5.64</theoretical-pi>
    <molecular-weight>70036.26</molecular-weight>
    <chromosome-location>11</chromosome-location>
    <organism ncbi-taxonomy-id="9606">Humans</organism>
    <external-identifiers>
      <external-identifier>
        <resource>UniProtKB</resource>
        <identifier>P00734</identifier>
      </external-identifier>
    </external-identifiers>
    <amino-acid-sequence format="FASTA">...</amino-acid-sequence>
    <gene-sequence format="FASTA">...</gene-sequence>
    <pfams>
      <pfam>
        <identifier>PF00051</identifier>
        <name>Kringle</name>
      </pfam>
    </pfams>
    <go-classifiers>
      <go-classifier>
        <category>function</category>
        <description>serine-type endopeptidase activity</description>
      </go-classifier>
    </go-classifiers>
  </polypeptide>
</target>
```

**Interactant Fields:**

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | DrugBank entity ID (BE#######) |
| `name` | string | Protein name |
| `organism` | string | Source organism |
| `actions` | list | Action types |
| `known-action` | enum | yes, no, unknown |
| `inhibition-strength` | enum | strong, moderate, weak, unknown (enzymes) |
| `induction-strength` | enum | strong, moderate, weak, unknown (enzymes) |

### 1.4 Drug-Target Action Types

| Action Type | Classification | Description |
|-------------|---------------|-------------|
| `inhibitor` | Negative | Reduces target activity |
| `antagonist` | Negative | Blocks receptor activation |
| `blocker` | Negative | Prevents ion/substrate passage |
| `suppressor` | Negative | Reduces expression/activity |
| `activator` | Positive | Increases target activity |
| `agonist` | Positive | Activates receptor |
| `potentiator` | Positive | Enhances activity |
| `inducer` | Positive | Increases expression |
| `partial agonist` | Mixed | Partial receptor activation |
| `inverse agonist` | Negative | Produces opposite effect |
| `modulator` | Neutral | Modifies activity (direction varies) |
| `binder` | Neutral | Binds without specified effect |
| `substrate` | Neutral | Metabolized by enzyme |
| `product of` | Neutral | Product of enzymatic reaction |
| `cofactor` | Neutral | Required for enzyme activity |
| `ligand` | Neutral | Binds to target |
| `cleavage` | Neutral | Cleaves substrate |

### 1.5 SNP Effects Schema

#### 1.5.1 SNP-Associated Drug Effects (`<snp-effects>`)

```xml
<snp-effects>
  <effect>
    <protein-name>Cytochrome P450 2D6</protein-name>
    <gene-symbol>CYP2D6</gene-symbol>
    <uniprot-id>P10635</uniprot-id>
    <rs-id>rs3892097</rs-id>
    <allele>CYP2D6*4</allele>
    <defining-change>1846G>A</defining-change>
    <description>Poor metabolizer phenotype</description>
    <pubmed-id>15625014</pubmed-id>
  </effect>
</snp-effects>
```

| Field | Type | Description |
|-------|------|-------------|
| `protein-name` | string | Affected protein name |
| `gene-symbol` | string | Gene symbol (e.g., CYP2D6) |
| `uniprot-id` | string | UniProt accession |
| `rs-id` | string | dbSNP reference ID |
| `allele` | string | Allele designation |
| `defining-change` | string | Nucleotide/amino acid change |
| `description` | text | Effect description |
| `pubmed-id` | integer | PubMed reference |

#### 1.5.2 SNP-Associated Adverse Drug Reactions (`<snp-adverse-drug-reactions>`)

```xml
<snp-adverse-drug-reactions>
  <reaction>
    <protein-name>ATP-binding cassette sub-family B member 1</protein-name>
    <gene-symbol>ABCB1</gene-symbol>
    <uniprot-id>P08183</uniprot-id>
    <rs-id>rs1045642</rs-id>
    <allele>T</allele>
    <adverse-reaction>Increased drug toxicity</adverse-reaction>
    <description>...</description>
    <pubmed-id>16220112</pubmed-id>
  </reaction>
</snp-adverse-drug-reactions>
```

### 1.6 Drug Interactions Format

#### 1.6.1 Simple Drug Interactions

```xml
<drug-interactions>
  <drug-interaction>
    <drugbank-id>DB00091</drugbank-id>
    <name>Cyclosporine</name>
    <description>The therapeutic efficacy of Lepirudin can be decreased
                 when used in combination with Cyclosporine.</description>
  </drug-interaction>
</drug-interactions>
```

#### 1.6.2 Structured Drug Interactions (DrugBank 6.0+)

```xml
<structured-drug-interactions>
  <structured-drug-interaction>
    <subject-drug>
      <drugbank-id>DB00001</drugbank-id>
      <name>Lepirudin</name>
      <role>object</role>
    </subject-drug>
    <affected-drug>
      <drugbank-id>DB00091</drugbank-id>
      <name>Cyclosporine</name>
    </affected-drug>
    <severity>moderate</severity>
    <description>...</description>
    <extended-description>...</extended-description>
    <management>Monitor closely</management>
    <mechanism-of-interaction>
      Cyclosporine inhibits P-glycoprotein efflux transporters
    </mechanism-of-interaction>
  </structured-drug-interaction>
</structured-drug-interactions>
```

**Severity Levels:**

| Level | Description |
|-------|-------------|
| `minor` | Minimal clinical significance |
| `moderate` | May require monitoring/dose adjustment |
| `major` | Avoid combination or use with extreme caution |
| `contraindicated` | Should not be used together |

### 1.7 External Identifiers

```xml
<external-identifiers>
  <external-identifier>
    <resource>UniProtKB</resource>
    <identifier>P01050</identifier>
  </external-identifier>
  <external-identifier>
    <resource>PubChem Compound</resource>
    <identifier>123456</identifier>
  </external-identifier>
</external-identifiers>
```

**Supported External Resources:**

| Resource | Description | Example |
|----------|-------------|---------|
| UniProtKB | Protein database | P00734 |
| PubChem Compound | PubChem CID | 123456 |
| PubChem Substance | PubChem SID | 46507042 |
| ChEMBL | Bioactivity database | CHEMBL221959 |
| ChEBI | Chemical entities | CHEBI:6443 |
| KEGG Drug | KEGG database | D00022 |
| KEGG Compound | KEGG compound | C07481 |
| PharmGKB | Pharmacogenomics | PA164748896 |
| PDB | Protein Data Bank | 1A28 |
| BindingDB | Binding data | 50001934 |
| Wikipedia | Wikipedia article | Lepirudin |
| ChemSpider | ChemSpider ID | 5145 |
| ZINC | ZINC database | ZINC000000005145 |
| RxCUI | RxNorm identifier | 120608 |
| GenBank | Nucleotide sequence | J00228 |
| Guide to Pharmacology | IUPHAR/BPS | 6321 |
| DrugProduct Database (DPD) | Health Canada | 02242077 |
| National Drug Code (NDC) | FDA NDC | 0002-3227 |

### 1.8 Products Element

```xml
<products>
  <product>
    <name>Refludan</name>
    <labeller>Bayer Healthcare Pharmaceuticals Inc</labeller>
    <ndc-id>50419-150-01</ndc-id>
    <ndc-product-code>50419-150</ndc-product-code>
    <dpd-id>02242963</dpd-id>
    <ema-product-code>EMEA/H/C/000122</ema-product-code>
    <ema-ma-number>EU/1/97/035/001</ema-ma-number>
    <started-marketing-on>1998-03-06</started-marketing-on>
    <ended-marketing-on></ended-marketing-on>
    <dosage-form>Powder for solution</dosage-form>
    <strength>50 mg</strength>
    <route>Intravenous</route>
    <fda-application-number>BLA103829</fda-application-number>
    <generic>false</generic>
    <over-the-counter>false</over-the-counter>
    <approved>true</approved>
    <country>US</country>
    <source>FDA NDC</source>
  </product>
</products>
```

### 1.9 Example DrugBank Record (JSON API Response)

```json
{
  "drugbank_id": "DB00001",
  "name": "Lepirudin",
  "type": "biotech",
  "description": "Lepirudin is a recombinant form of hirudin...",
  "cas_number": "120993-53-5",
  "unii": "Y43GF64R34",
  "state": "liquid",
  "groups": ["approved"],
  "indication": "For the treatment of heparin-induced thrombocytopenia...",
  "mechanism_of_action": "Lepirudin forms a stable non-covalent complex with alpha-thrombin...",
  "pharmacodynamics": "Lepirudin is used to break up clots...",
  "absorption": "Bioavailability is 100% following subcutaneous injection.",
  "half_life": "Approximately 1.3 hours",
  "protein_binding": "Not protein bound",
  "external_identifiers": [
    {"resource": "UniProtKB", "identifier": "P01050"},
    {"resource": "PubChem Compound", "identifier": "16132432"}
  ],
  "targets": [
    {
      "id": "BE0000048",
      "name": "Prothrombin",
      "organism": "Humans",
      "actions": ["inhibitor"],
      "known_action": "yes",
      "uniprot_id": "P00734",
      "gene_name": "F2"
    }
  ]
}
```

---

## 2. ChEMBL

### 2.1 Overview

ChEMBL is a manually curated database of bioactive molecules with drug-like properties, maintained by the European Bioinformatics Institute (EMBL-EBI). It brings together chemical, bioactivity, and genomic data to aid in the translation of genomic information into effective new drugs.

**Current Version:** ChEMBL 34 (2024)
**License:** CC BY-SA 3.0
**Website:** https://www.ebi.ac.uk/chembl/
**Schema Documentation:** https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html

### 2.2 Database Schema (Core Tables)

#### 2.2.1 MOLECULE_DICTIONARY

The primary repository for all compounds in ChEMBL.

```sql
CREATE TABLE molecule_dictionary (
    molregno            BIGINT PRIMARY KEY,
    pref_name           VARCHAR(255),
    chembl_id           VARCHAR(20) NOT NULL UNIQUE,
    max_phase           SMALLINT,
    therapeutic_flag    SMALLINT,
    dosed_ingredient    SMALLINT,
    structure_type      VARCHAR(10),
    chebi_par_id        INTEGER,
    molecule_type       VARCHAR(30),
    first_approval      SMALLINT,
    oral                SMALLINT,
    parenteral          SMALLINT,
    topical             SMALLINT,
    black_box_warning   SMALLINT,
    natural_product     SMALLINT,
    first_in_class      SMALLINT,
    chirality           SMALLINT,
    prodrug             SMALLINT,
    inorganic_flag      SMALLINT,
    usan_year           SMALLINT,
    availability_type   SMALLINT,
    usan_stem           VARCHAR(50),
    polymer_flag        SMALLINT,
    usan_substem        VARCHAR(50),
    usan_stem_definition VARCHAR(1000),
    indication_class    VARCHAR(1000),
    withdrawn_flag      SMALLINT,
    withdrawn_year      SMALLINT,
    withdrawn_country   VARCHAR(2000),
    withdrawn_reason    VARCHAR(2000)
);
```

| Column | Type | Description |
|--------|------|-------------|
| `molregno` | BIGINT | Internal primary key |
| `chembl_id` | VARCHAR(20) | Public ChEMBL identifier (e.g., CHEMBL25) |
| `pref_name` | VARCHAR(255) | Preferred compound name |
| `max_phase` | SMALLINT | Highest development phase (0-4) |
| `therapeutic_flag` | SMALLINT | Therapeutic use flag (1=yes) |
| `molecule_type` | VARCHAR(30) | Small molecule, Protein, Antibody, etc. |
| `structure_type` | VARCHAR(10) | MOL (small molecule), SEQ (protein), NONE |
| `first_approval` | SMALLINT | Year of first approval |
| `oral` | SMALLINT | Oral administration (1=yes) |
| `parenteral` | SMALLINT | Parenteral administration (1=yes) |
| `topical` | SMALLINT | Topical administration (1=yes) |
| `black_box_warning` | SMALLINT | Has black box warning (1=yes) |
| `natural_product` | SMALLINT | Natural product derived (1=yes) |
| `first_in_class` | SMALLINT | First in class drug (1=yes) |
| `prodrug` | SMALLINT | Is a prodrug (1=yes) |
| `chirality` | SMALLINT | Chirality: -1=achiral, 0=unknown, 1=single, 2=racemic |
| `withdrawn_flag` | SMALLINT | Withdrawn from any market (1=yes) |

**Max Phase Values:**

| Value | Description |
|-------|-------------|
| 0 | Preclinical / Unknown |
| 0.5 | Early Phase 1 Clinical Trials |
| 1 | Phase 1 Clinical Trials |
| 2 | Phase 2 Clinical Trials |
| 3 | Phase 3 Clinical Trials |
| 4 | Approved |

#### 2.2.2 COMPOUND_STRUCTURES

```sql
CREATE TABLE compound_structures (
    molregno            BIGINT PRIMARY KEY,
    molfile             TEXT,
    standard_inchi      TEXT,
    standard_inchi_key  VARCHAR(27),
    canonical_smiles    TEXT
);
```

| Column | Type | Description |
|--------|------|-------------|
| `molregno` | BIGINT | Foreign key to molecule_dictionary |
| `molfile` | TEXT | MDL MOL format structure |
| `standard_inchi` | TEXT | Standard InChI string |
| `standard_inchi_key` | VARCHAR(27) | InChI Key |
| `canonical_smiles` | TEXT | Canonical SMILES (RDKit) |

#### 2.2.3 COMPOUND_PROPERTIES

```sql
CREATE TABLE compound_properties (
    molregno            BIGINT PRIMARY KEY,
    mw_freebase         NUMERIC(9,2),
    alogp               NUMERIC(9,2),
    hba                 SMALLINT,
    hbd                 SMALLINT,
    psa                 NUMERIC(9,2),
    rtb                 SMALLINT,
    ro3_pass            VARCHAR(3),
    num_ro5_violations  SMALLINT,
    cx_most_apka        NUMERIC(9,2),
    cx_most_bpka        NUMERIC(9,2),
    cx_logp             NUMERIC(9,2),
    cx_logd             NUMERIC(9,2),
    molecular_species   VARCHAR(50),
    full_mwt            NUMERIC(9,2),
    aromatic_rings      SMALLINT,
    heavy_atoms         SMALLINT,
    qed_weighted        NUMERIC(3,2),
    mw_monoisotopic     NUMERIC(11,4),
    full_molformula     VARCHAR(100),
    hba_lipinski        SMALLINT,
    hbd_lipinski        SMALLINT,
    num_lipinski_ro5_violations SMALLINT
);
```

| Column | Type | Description |
|--------|------|-------------|
| `mw_freebase` | NUMERIC | Molecular weight (parent compound) |
| `alogp` | NUMERIC | ALogP partition coefficient |
| `hba` | SMALLINT | Hydrogen bond acceptors |
| `hbd` | SMALLINT | Hydrogen bond donors |
| `psa` | NUMERIC | Polar surface area (A^2) |
| `rtb` | SMALLINT | Rotatable bonds |
| `num_ro5_violations` | SMALLINT | Lipinski Rule of 5 violations |
| `aromatic_rings` | SMALLINT | Number of aromatic rings |
| `heavy_atoms` | SMALLINT | Number of heavy atoms |
| `qed_weighted` | NUMERIC | Quantitative Estimate of Drug-likeness |
| `cx_logp` | NUMERIC | ChemAxon LogP |
| `cx_logd` | NUMERIC | ChemAxon LogD at pH 7.4 |
| `cx_most_apka` | NUMERIC | Most acidic pKa |
| `cx_most_bpka` | NUMERIC | Most basic pKa |

#### 2.2.4 ASSAYS

```sql
CREATE TABLE assays (
    assay_id            BIGINT PRIMARY KEY,
    doc_id              BIGINT,
    description         TEXT,
    assay_type          VARCHAR(1),
    assay_test_type     VARCHAR(20),
    assay_category      VARCHAR(20),
    assay_organism      VARCHAR(250),
    assay_tax_id        BIGINT,
    assay_strain        VARCHAR(200),
    assay_tissue        VARCHAR(100),
    assay_cell_type     VARCHAR(100),
    assay_subcellular_fraction VARCHAR(100),
    tid                 BIGINT,
    relationship_type   VARCHAR(1),
    confidence_score    SMALLINT,
    curated_by          VARCHAR(32),
    src_id              SMALLINT,
    src_assay_id        VARCHAR(50),
    chembl_id           VARCHAR(20) UNIQUE,
    cell_id             BIGINT,
    bao_format          VARCHAR(11),
    tissue_id           BIGINT,
    variant_id          BIGINT,
    aidx                VARCHAR(200)
);
```

| Column | Type | Description |
|--------|------|-------------|
| `assay_id` | BIGINT | Primary key |
| `chembl_id` | VARCHAR(20) | Public identifier (CHEMBL#######) |
| `doc_id` | BIGINT | Foreign key to docs |
| `tid` | BIGINT | Foreign key to target_dictionary |
| `assay_type` | VARCHAR(1) | B=Binding, F=Functional, A=ADMET, T=Toxicity, P=Physicochemical, U=Unclassified |
| `assay_test_type` | VARCHAR(20) | In vivo, In vitro, Ex vivo |
| `assay_category` | VARCHAR(20) | Screening, Confirmatory, etc. |
| `confidence_score` | SMALLINT | Target mapping confidence (0-9) |
| `relationship_type` | VARCHAR(1) | D=Direct, H=Homologous, M=Molecular, N=Non-molecular, S=Subcellular, U=Uncurated |

**Assay Type Values:**

| Code | Description |
|------|-------------|
| B | Binding assay |
| F | Functional assay |
| A | ADMET assay |
| T | Toxicity assay |
| P | Physicochemical assay |
| U | Unclassified |

**Confidence Score:**

| Score | Meaning |
|-------|---------|
| 9 | Direct single protein target |
| 8 | Homologous single protein target |
| 7 | Direct protein complex |
| 6 | Homologous protein complex |
| 5 | Direct protein family |
| 4 | Homologous protein family |
| 3 | Direct protein selectivity group |
| 1 | Default (non-molecular target) |
| 0 | Uncurated/unassigned |

#### 2.2.5 ACTIVITIES

```sql
CREATE TABLE activities (
    activity_id             BIGINT PRIMARY KEY,
    assay_id                BIGINT NOT NULL,
    doc_id                  BIGINT,
    record_id               BIGINT NOT NULL,
    molregno                BIGINT,
    standard_relation       VARCHAR(50),
    standard_value          NUMERIC,
    standard_units          VARCHAR(100),
    standard_flag           SMALLINT,
    standard_type           VARCHAR(250),
    activity_comment        VARCHAR(4000),
    data_validity_comment   VARCHAR(30),
    potential_duplicate     SMALLINT,
    pchembl_value           NUMERIC(4,2),
    bao_endpoint            VARCHAR(11),
    uo_units                VARCHAR(10),
    qudt_units              VARCHAR(70),
    toid                    BIGINT,
    upper_value             NUMERIC,
    standard_upper_value    NUMERIC,
    src_id                  SMALLINT,
    type                    VARCHAR(250),
    relation                VARCHAR(50),
    value                   NUMERIC,
    units                   VARCHAR(100),
    text_value              VARCHAR(1000),
    standard_text_value     VARCHAR(1000),
    action_type             VARCHAR(50)
);
```

| Column | Type | Description |
|--------|------|-------------|
| `activity_id` | BIGINT | Primary key |
| `assay_id` | BIGINT | Foreign key to assays |
| `molregno` | BIGINT | Foreign key to molecule_dictionary |
| `standard_type` | VARCHAR(250) | Standardized activity type |
| `standard_value` | NUMERIC | Standardized numeric value |
| `standard_units` | VARCHAR(100) | Standardized units |
| `standard_relation` | VARCHAR(50) | Relation (=, <, >, <=, >=, ~) |
| `pchembl_value` | NUMERIC | -log10(molar) value for comparison |
| `data_validity_comment` | VARCHAR(30) | Data quality flag |
| `action_type` | VARCHAR(50) | Action type from ACTION_TYPE table |

### 2.3 Activity Types (IC50, Ki, EC50, Kd)

#### 2.3.1 Standard Activity Types

| Type | Description | Typical Units |
|------|-------------|---------------|
| IC50 | Half maximal inhibitory concentration | nM |
| Ki | Inhibition constant | nM |
| EC50 | Half maximal effective concentration | nM |
| Kd | Dissociation constant | nM |
| AC50 | Half maximal activity concentration | nM |
| XC50 | Half maximal concentration (unspecified) | nM |
| Potency | General potency measure | nM |
| ED50 | Effective dose 50% | mg/kg |
| MIC | Minimum inhibitory concentration | ug/mL |
| GI50 | Growth inhibition 50% | nM |
| LC50 | Lethal concentration 50% | uM |
| LD50 | Lethal dose 50% | mg/kg |

#### 2.3.2 pChEMBL Value

The pChEMBL value allows comparison across different activity types:

```
pChEMBL = -log10(molar value)

Example: IC50 = 1 nM = 1e-9 M
         pChEMBL = -log10(1e-9) = 9
```

| pChEMBL | Approximate Activity |
|---------|---------------------|
| 9 | 1 nM |
| 8 | 10 nM |
| 7 | 100 nM |
| 6 | 1 uM |
| 5 | 10 uM |

#### 2.2.6 TARGET_DICTIONARY

```sql
CREATE TABLE target_dictionary (
    tid                 BIGINT PRIMARY KEY,
    target_type         VARCHAR(30),
    pref_name           VARCHAR(200),
    tax_id              BIGINT,
    organism            VARCHAR(150),
    chembl_id           VARCHAR(20) UNIQUE,
    species_group_flag  SMALLINT
);
```

| Column | Type | Description |
|--------|------|-------------|
| `tid` | BIGINT | Primary key |
| `chembl_id` | VARCHAR(20) | Public identifier (CHEMBL#######) |
| `target_type` | VARCHAR(30) | Target classification |
| `pref_name` | VARCHAR(200) | Preferred name |
| `organism` | VARCHAR(150) | Target organism |
| `tax_id` | BIGINT | NCBI taxonomy ID |
| `species_group_flag` | SMALLINT | Representative species (1=yes) |

**Target Types:**

| Type | Description |
|------|-------------|
| SINGLE PROTEIN | Single protein target |
| PROTEIN COMPLEX | Multi-protein complex |
| PROTEIN FAMILY | Group of related proteins |
| PROTEIN-PROTEIN INTERACTION | PPI target |
| CHIMERIC PROTEIN | Engineered fusion protein |
| SELECTIVITY GROUP | Selectivity panel |
| NUCLEIC-ACID | DNA/RNA target |
| CELL-LINE | Cell line assay target |
| TISSUE | Tissue-based target |
| ORGANISM | Whole organism |
| SUBCELLULAR | Subcellular fraction |
| PHENOTYPE | Phenotypic endpoint |
| SMALL MOLECULE | Small molecule target |
| MACROMOLECULE | Other macromolecule |

#### 2.2.7 TARGET_COMPONENTS

```sql
CREATE TABLE target_components (
    tid                 BIGINT,
    component_id        BIGINT,
    targcomp_id         BIGINT PRIMARY KEY,
    homologue           SMALLINT
);
```

#### 2.2.8 COMPONENT_SEQUENCES

```sql
CREATE TABLE component_sequences (
    component_id        BIGINT PRIMARY KEY,
    component_type      VARCHAR(50),
    accession           VARCHAR(25),
    sequence            TEXT,
    sequence_md5sum     VARCHAR(32),
    description         VARCHAR(200),
    tax_id              BIGINT,
    organism            VARCHAR(150),
    db_source           VARCHAR(25),
    db_version          VARCHAR(10)
);
```

| Column | Type | Description |
|--------|------|-------------|
| `component_id` | BIGINT | Primary key |
| `accession` | VARCHAR(25) | UniProt accession |
| `sequence` | TEXT | Protein sequence |
| `component_type` | VARCHAR(50) | PROTEIN, DNA, RNA |
| `tax_id` | BIGINT | NCBI taxonomy ID |
| `organism` | VARCHAR(150) | Source organism |

### 2.4 Mechanism of Action Schema

#### 2.4.1 DRUG_MECHANISM

```sql
CREATE TABLE drug_mechanism (
    mec_id              BIGINT PRIMARY KEY,
    record_id           BIGINT NOT NULL,
    molregno            BIGINT,
    mechanism_of_action VARCHAR(250),
    tid                 BIGINT,
    site_id             BIGINT,
    action_type         VARCHAR(50),
    direct_interaction  SMALLINT,
    molecular_mechanism SMALLINT,
    disease_efficacy    SMALLINT,
    mechanism_comment   VARCHAR(2000),
    selectivity_comment VARCHAR(1000),
    binding_site_comment VARCHAR(1000)
);
```

| Column | Type | Description |
|--------|------|-------------|
| `mec_id` | BIGINT | Primary key |
| `molregno` | BIGINT | Foreign key to molecule_dictionary |
| `tid` | BIGINT | Foreign key to target_dictionary |
| `mechanism_of_action` | VARCHAR(250) | MoA description (e.g., "PDE5 inhibitor") |
| `action_type` | VARCHAR(50) | Action type (INHIBITOR, AGONIST, etc.) |
| `direct_interaction` | SMALLINT | Direct target engagement (1=yes) |
| `molecular_mechanism` | SMALLINT | Molecular (vs physiological) (1=yes) |
| `disease_efficacy` | SMALLINT | Target involved in efficacy (1=yes) |

#### 2.4.2 ACTION_TYPE Table

```sql
CREATE TABLE action_type (
    action_type         VARCHAR(50) PRIMARY KEY,
    description         VARCHAR(200),
    parent_type         VARCHAR(50)
);
```

**Action Types:**

| Action Type | Parent Type | Description |
|-------------|-------------|-------------|
| INHIBITOR | NEGATIVE MODULATOR | Enzyme/protein inhibitor |
| ANTAGONIST | NEGATIVE MODULATOR | Receptor antagonist |
| BLOCKER | NEGATIVE MODULATOR | Channel blocker |
| NEGATIVE ALLOSTERIC MODULATOR | NEGATIVE MODULATOR | Negative allosteric effect |
| AGONIST | POSITIVE MODULATOR | Receptor agonist |
| PARTIAL AGONIST | POSITIVE MODULATOR | Partial receptor agonist |
| ACTIVATOR | POSITIVE MODULATOR | Enzyme/protein activator |
| POSITIVE ALLOSTERIC MODULATOR | POSITIVE MODULATOR | Positive allosteric effect |
| OPENER | POSITIVE MODULATOR | Channel opener |
| INVERSE AGONIST | NEGATIVE MODULATOR | Inverse agonist |
| SUBSTRATE | OTHER | Enzyme substrate |
| RELEASING AGENT | OTHER | Neurotransmitter releaser |
| BINDING AGENT | OTHER | General binder |
| MODULATOR | OTHER | Non-specific modulation |

#### 2.4.3 MECHANISM_REFS

```sql
CREATE TABLE mechanism_refs (
    mecref_id           BIGINT PRIMARY KEY,
    mec_id              BIGINT NOT NULL,
    ref_type            VARCHAR(50),
    ref_id              VARCHAR(200),
    ref_url             VARCHAR(400)
);
```

| ref_type | Description | Example ref_id |
|----------|-------------|----------------|
| PubMed | PubMed article | 12345678 |
| DailyMed | Drug label | setid-xxxx |
| FDA | FDA document | NDA012345 |
| ISBN | Book reference | 978-0-xxxx |
| ClinicalTrials | Trial | NCT01234567 |
| DOI | DOI reference | 10.1000/xxxx |
| Wikipedia | Wikipedia | Article_title |

### 2.5 UniProt Relationship

ChEMBL targets are linked to UniProt through the component_sequences table:

```sql
-- Get UniProt accessions for a ChEMBL target
SELECT DISTINCT
    td.chembl_id AS target_chembl_id,
    td.pref_name AS target_name,
    cs.accession AS uniprot_accession,
    cs.description AS protein_description,
    cs.organism
FROM target_dictionary td
JOIN target_components tc ON td.tid = tc.tid
JOIN component_sequences cs ON tc.component_id = cs.component_id
WHERE td.chembl_id = 'CHEMBL203'
  AND cs.db_source = 'Swiss-Prot';
```

### 2.6 REST API Endpoints

**Base URL:** `https://www.ebi.ac.uk/chembl/api/data/`

#### 2.6.1 Available Resources

| Endpoint | Description |
|----------|-------------|
| `/molecule` | Compound information |
| `/molecule/{chembl_id}` | Single compound |
| `/target` | Target information |
| `/target/{chembl_id}` | Single target |
| `/assay` | Assay definitions |
| `/activity` | Bioactivity measurements |
| `/mechanism` | Mechanism of action |
| `/drug_indication` | Drug indications |
| `/drug_warning` | Drug warnings/withdrawals |
| `/document` | Publications |
| `/similarity/{smiles}/{similarity}` | Similarity search |
| `/substructure/{smiles}` | Substructure search |
| `/target_component` | Target protein components |
| `/image/{chembl_id}` | Structure image |

#### 2.6.2 Query Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `format` | Output format | json, xml, yaml |
| `limit` | Results per page (max 1000) | 20 |
| `offset` | Pagination offset | 0 |
| `molecule_chembl_id` | Filter by compound | CHEMBL25 |
| `target_chembl_id` | Filter by target | CHEMBL203 |
| `assay_chembl_id` | Filter by assay | CHEMBL1217643 |
| `max_phase` | Development phase | 4 |
| `standard_type` | Activity type | IC50 |

#### 2.6.3 Example API Responses

**Molecule Response:**
```json
{
  "molecule_chembl_id": "CHEMBL25",
  "pref_name": "ASPIRIN",
  "molecule_type": "Small molecule",
  "max_phase": 4,
  "structure_type": "MOL",
  "molecule_structures": {
    "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "standard_inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
    "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
  },
  "molecule_properties": {
    "mw_freebase": 180.16,
    "alogp": 1.31,
    "hba": 3,
    "hbd": 1,
    "psa": 63.60,
    "rtb": 3,
    "num_ro5_violations": 0,
    "qed_weighted": 0.55
  }
}
```

**Activity Response:**
```json
{
  "activity_id": 31863,
  "assay_chembl_id": "CHEMBL1217643",
  "assay_type": "B",
  "molecule_chembl_id": "CHEMBL25",
  "target_chembl_id": "CHEMBL2094253",
  "standard_type": "IC50",
  "standard_value": 26000.0,
  "standard_units": "nM",
  "standard_relation": "=",
  "pchembl_value": 4.59,
  "data_validity_comment": null
}
```

**Mechanism Response:**
```json
{
  "mechanism_of_action": "Cyclooxygenase inhibitor",
  "molecule_chembl_id": "CHEMBL25",
  "target_chembl_id": "CHEMBL2094253",
  "action_type": "INHIBITOR",
  "direct_interaction": true,
  "molecular_mechanism": true,
  "max_phase": 4,
  "mechanism_refs": [
    {
      "ref_type": "PubMed",
      "ref_id": "7678454",
      "ref_url": "http://europepmc.org/abstract/MED/7678454"
    }
  ]
}
```

---

## 3. PubChem

### 3.1 Overview

PubChem is a public repository for information on chemical substances and their biological activities, maintained by the National Center for Biotechnology Information (NCBI) at the National Institutes of Health (NIH).

**Website:** https://pubchem.ncbi.nlm.nih.gov
**API Documentation:** https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
**Launched:** 2004

### 3.2 Data Organization

PubChem organizes data into three interconnected databases:

| Database | Identifier | Description |
|----------|------------|-------------|
| **Substance** | SID | Depositor-provided sample information |
| **Compound** | CID | Unique standardized chemical structures |
| **BioAssay** | AID | Biological assay results |

### 3.3 Compound (CID) Schema

#### 3.3.1 Basic Compound Properties

| Property | Type | Description |
|----------|------|-------------|
| `CID` | integer | Unique compound identifier |
| `MolecularFormula` | string | Molecular formula (e.g., C9H8O4) |
| `MolecularWeight` | float | Molecular weight |
| `CanonicalSMILES` | string | Canonical SMILES notation |
| `IsomericSMILES` | string | Isomeric SMILES (with stereochemistry) |
| `InChI` | string | IUPAC International Chemical Identifier |
| `InChIKey` | string | Hashed InChI (27 characters) |
| `IUPACName` | string | IUPAC systematic name |
| `Title` | string | Preferred compound name |

#### 3.3.2 Physicochemical Properties

| Property | Type | Unit | Description |
|----------|------|------|-------------|
| `XLogP` | float | - | Computed octanol-water partition coefficient |
| `ExactMass` | float | Da | Exact mass |
| `MonoisotopicMass` | float | Da | Monoisotopic mass |
| `TPSA` | float | A^2 | Topological polar surface area |
| `Complexity` | float | - | Molecular complexity (Bertz formula) |
| `Charge` | integer | - | Formal charge |
| `HBondDonorCount` | integer | - | Hydrogen bond donors |
| `HBondAcceptorCount` | integer | - | Hydrogen bond acceptors |
| `RotatableBondCount` | integer | - | Rotatable bonds |
| `HeavyAtomCount` | integer | - | Non-hydrogen atoms |
| `IsotopeAtomCount` | integer | - | Isotope-labeled atoms |
| `AtomStereoCount` | integer | - | Stereocenters (atom) |
| `DefinedAtomStereoCount` | integer | - | Defined stereocenters |
| `UndefinedAtomStereoCount` | integer | - | Undefined stereocenters |
| `BondStereoCount` | integer | - | Stereo bonds |
| `DefinedBondStereoCount` | integer | - | Defined stereo bonds |
| `UndefinedBondStereoCount` | integer | - | Undefined stereo bonds |
| `CovalentUnitCount` | integer | - | Covalently bonded units |

#### 3.3.3 3D Properties (Conformer-Dependent)

| Property | Type | Description |
|----------|------|-------------|
| `Volume3D` | float | Molecular volume |
| `XStericQuadrupole3D` | float | X-axis steric quadrupole |
| `YStericQuadrupole3D` | float | Y-axis steric quadrupole |
| `ZStericQuadrupole3D` | float | Z-axis steric quadrupole |
| `FeatureCount3D` | integer | Total pharmacophore features |
| `FeatureAcceptorCount3D` | integer | Acceptor features |
| `FeatureDonorCount3D` | integer | Donor features |
| `FeatureAnionCount3D` | integer | Anionic features |
| `FeatureCationCount3D` | integer | Cationic features |
| `FeatureRingCount3D` | integer | Ring features |
| `FeatureHydrophobeCount3D` | integer | Hydrophobic features |
| `ConformerModelRMSD3D` | float | RMSD of conformer model |
| `EffectiveRotorCount3D` | float | Effective rotatable bonds |
| `ConformerCount3D` | integer | Number of conformers |

### 3.4 Substance (SID) Schema

```json
{
  "SID": 123456789,
  "SourceName": "ChEMBL",
  "SourceID": "CHEMBL25",
  "Synonyms": [
    "aspirin",
    "acetylsalicylic acid",
    "2-acetoxybenzoic acid"
  ],
  "RegistryID": "50-78-2",
  "Comment": "Depositor-provided description",
  "CreateDate": "2004-09-16",
  "ModifyDate": "2024-01-15",
  "CompoundCID": 2244,
  "Standardized": true
}
```

| Field | Type | Description |
|-------|------|-------------|
| `SID` | integer | Unique substance identifier |
| `SourceName` | string | Depositor name |
| `SourceID` | string | Depositor's internal ID |
| `Synonyms` | array | Chemical names |
| `RegistryID` | string | CAS number or registry ID |
| `CompoundCID` | integer | Linked compound ID |
| `CreateDate` | date | Submission date |
| `ModifyDate` | date | Last update date |

### 3.5 BioAssay (AID) Schema

#### 3.5.1 Assay Description

```json
{
  "AID": 1851,
  "SourceName": "ChEMBL",
  "Name": "Inhibition of human recombinant COX-2",
  "Description": [
    "Inhibition of cyclooxygenase-2 enzyme activity"
  ],
  "Protocol": [
    "Compounds were tested for inhibition of COX-2...",
    "IC50 values were determined using..."
  ],
  "AssayType": "Confirmatory",
  "ActivityOutcomeMethod": "Confirmatory",
  "ProjectCategory": 1,
  "Target": {
    "Name": "Prostaglandin G/H synthase 2",
    "Type": "protein",
    "GeneID": 5743,
    "GeneSymbol": "PTGS2",
    "TargetURL": "https://www.ncbi.nlm.nih.gov/gene/5743"
  },
  "XRef": [
    {
      "DBName": "GeneID",
      "ID": "5743"
    },
    {
      "DBName": "UniProt",
      "ID": "P35354"
    }
  ],
  "Results": [
    {
      "TID": 1,
      "Name": "IC50",
      "Description": "Half maximal inhibitory concentration",
      "Type": "float",
      "Unit": "um"
    },
    {
      "TID": 2,
      "Name": "Outcome",
      "Description": "Activity outcome",
      "Type": "string"
    }
  ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `AID` | integer | Unique assay identifier |
| `Name` | string | Assay name |
| `Description` | array | Detailed description |
| `Protocol` | array | Experimental protocol |
| `AssayType` | string | Screening, Confirmatory, Summary |
| `Target` | object | Target information |
| `Results` | array | Result type definitions (TID) |

#### 3.5.2 Assay Data

```json
{
  "AID": 1851,
  "Data": [
    {
      "SID": 12345678,
      "CID": 2244,
      "Outcome": 2,
      "OutcomeDescription": "Active",
      "Data": [
        {
          "TID": 1,
          "Value": 0.026
        }
      ]
    }
  ]
}
```

**Outcome Values:**

| Value | Description |
|-------|-------------|
| 1 | Inactive |
| 2 | Active |
| 3 | Inconclusive |
| 4 | Unspecified |
| 5 | Probe |

### 3.6 PUG REST API Response Formats

#### 3.6.1 URL Structure

```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/<input>/<operation>/<output>
```

**Components:**

| Part | Options | Example |
|------|---------|---------|
| Input | compound/cid, compound/name, compound/smiles, compound/inchi, compound/inchikey, substance/sid, substance/name, assay/aid | compound/cid/2244 |
| Operation | record, property/[properties], synonyms, sids, cids, aids, description, targets, classification | property/MolecularWeight |
| Output | JSON, XML, CSV, SDF, PNG, TXT | JSON |

#### 3.6.2 Example API Calls

**Get Properties by CID:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey/JSON
```

Response:
```json
{
  "PropertyTable": {
    "Properties": [
      {
        "CID": 2244,
        "MolecularFormula": "C9H8O4",
        "MolecularWeight": 180.16,
        "CanonicalSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
      }
    ]
  }
}
```

**Search by Name:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/cids/JSON
```

Response:
```json
{
  "IdentifierList": {
    "CID": [2244]
  }
}
```

**Get Synonyms:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/synonyms/JSON
```

Response:
```json
{
  "InformationList": {
    "Information": [
      {
        "CID": 2244,
        "Synonym": [
          "aspirin",
          "acetylsalicylic acid",
          "Ecotrin",
          "Acenterine",
          "2-Acetoxybenzoic acid",
          "..."
        ]
      }
    ]
  }
}
```

**Get BioAssay Results:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/1851/JSON?cid=2244
```

#### 3.6.3 Property Fields Reference

**Complete property list for PUG REST:**

```
MolecularFormula, MolecularWeight, CanonicalSMILES, IsomericSMILES,
InChI, InChIKey, IUPACName, XLogP, ExactMass, MonoisotopicMass,
TPSA, Complexity, Charge, HBondDonorCount, HBondAcceptorCount,
RotatableBondCount, HeavyAtomCount, IsotopeAtomCount, AtomStereoCount,
DefinedAtomStereoCount, UndefinedAtomStereoCount, BondStereoCount,
DefinedBondStereoCount, UndefinedBondStereoCount, CovalentUnitCount,
Volume3D, XStericQuadrupole3D, YStericQuadrupole3D, ZStericQuadrupole3D,
FeatureCount3D, FeatureAcceptorCount3D, FeatureDonorCount3D,
FeatureAnionCount3D, FeatureCationCount3D, FeatureRingCount3D,
FeatureHydrophobeCount3D, ConformerModelRMSD3D, EffectiveRotorCount3D,
ConformerCount3D
```

### 3.7 PUG-View API

PUG-View provides access to annotation data with a different data model:

```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/2244/JSON
```

Response includes sections like:
- Pharmacology and Biochemistry
- Drug and Medication Information
- Safety and Hazards
- Literature
- Patents
- Classification
- Related Compounds

---

## 4. DrugCentral

### 4.1 Overview

DrugCentral is an online drug compendium developed and maintained by the Translational Informatics Division at the University of New Mexico. It integrates structure, bioactivity, regulatory, pharmacologic action, and disease information for approved drugs.

**Website:** https://drugcentral.org
**Download:** https://drugcentral.org/download
**Database:** PostgreSQL
**License:** CC BY 4.0

### 4.2 PostgreSQL Schema

#### 4.2.1 Core Tables

##### STRUCTURES (Central Table)

```sql
CREATE TABLE structures (
    id              SERIAL PRIMARY KEY,
    cd_id           INTEGER,
    cd_formula      VARCHAR(100),
    cd_molweight    NUMERIC(12,4),
    name            VARCHAR(500),
    cas_reg_no      VARCHAR(50),
    inchi           TEXT,
    inchikey        VARCHAR(27),
    smiles          TEXT,
    stem            VARCHAR(50),
    mol2d           TEXT,
    mol3d           TEXT,
    mrdef           TEXT,
    no_formulations INTEGER,
    active_moiety_struct_id INTEGER
);
```

| Column | Type | Description |
|--------|------|-------------|
| `id` | SERIAL | Primary key (drug structure ID) |
| `cd_formula` | VARCHAR | Molecular formula |
| `cd_molweight` | NUMERIC | Molecular weight |
| `name` | VARCHAR | Drug name |
| `cas_reg_no` | VARCHAR | CAS Registry Number |
| `inchi` | TEXT | InChI string |
| `inchikey` | VARCHAR | InChI Key |
| `smiles` | TEXT | SMILES notation |
| `stem` | VARCHAR | USAN stem |
| `active_moiety_struct_id` | INTEGER | Parent compound ID |

##### IDENTIFIER

```sql
CREATE TABLE identifier (
    id              SERIAL PRIMARY KEY,
    identifier      VARCHAR(100) NOT NULL,
    id_type         VARCHAR(50) NOT NULL,
    struct_id       INTEGER REFERENCES structures(id),
    parent_match    BOOLEAN
);
```

| Column | Type | Description |
|--------|------|-------------|
| `identifier` | VARCHAR | External database ID |
| `id_type` | VARCHAR | Source database name |
| `struct_id` | INTEGER | Foreign key to structures |
| `parent_match` | BOOLEAN | Matches parent structure |

**Identifier Types:**

| id_type | Description | Example |
|---------|-------------|---------|
| DRUGBANK_ID | DrugBank identifier | DB00945 |
| CHEMBL_ID | ChEMBL identifier | CHEMBL25 |
| PUBCHEM_CID | PubChem Compound ID | 2244 |
| UNII | FDA UNII | R16CO5Y76E |
| RXCUI | RxNorm concept | 1191 |
| MESH_ID | MeSH identifier | D001241 |
| KEGG_ID | KEGG Drug ID | D00109 |
| CHEBI_ID | ChEBI identifier | CHEBI:15365 |
| SNOMEDCT_US | SNOMED CT code | 387458008 |
| NDF-RT | NDF-RT identifier | N0000145918 |
| IUPHAR_LIGAND_ID | IUPHAR ligand | 4139 |

##### SYNONYMS

```sql
CREATE TABLE synonyms (
    id              SERIAL PRIMARY KEY,
    syn_id          INTEGER NOT NULL,
    name            VARCHAR(500) NOT NULL,
    preferred_name  SMALLINT DEFAULT 0,
    lname           VARCHAR(500),
    struct_id       INTEGER REFERENCES structures(id)
);
```

#### 4.2.2 Drug-Target Relationships

##### ACT_TABLE_FULL (Bioactivity Data)

```sql
CREATE TABLE act_table_full (
    id              SERIAL PRIMARY KEY,
    struct_id       INTEGER REFERENCES structures(id),
    target_id       INTEGER REFERENCES targets(id),
    target_name     VARCHAR(500),
    target_class    VARCHAR(100),
    accession       VARCHAR(20),
    gene            VARCHAR(50),
    swissprot       VARCHAR(20),
    act_value       NUMERIC(12,4),
    act_unit        VARCHAR(20),
    act_type        VARCHAR(50),
    act_comment     TEXT,
    act_source      VARCHAR(100),
    relation        VARCHAR(10),
    moa             SMALLINT,
    moa_source      VARCHAR(200),
    moa_source_url  VARCHAR(500),
    action_type     VARCHAR(50),
    tdl             VARCHAR(10)
);
```

| Column | Type | Description |
|--------|------|-------------|
| `struct_id` | INTEGER | Drug structure ID |
| `target_id` | INTEGER | Target ID |
| `target_name` | VARCHAR | Target name |
| `accession` | VARCHAR | UniProt accession |
| `gene` | VARCHAR | Gene symbol |
| `act_value` | NUMERIC | Activity value |
| `act_unit` | VARCHAR | Activity units (nM, uM, etc.) |
| `act_type` | VARCHAR | Activity type (Ki, IC50, EC50, Kd) |
| `act_source` | VARCHAR | Data source |
| `moa` | SMALLINT | Mechanism of action (1=yes) |
| `action_type` | VARCHAR | Action type (INHIBITOR, AGONIST, etc.) |
| `tdl` | VARCHAR | Target Development Level |

**Activity Types (act_type):**

| Type | Description |
|------|-------------|
| Ki | Inhibition constant |
| Kd | Dissociation constant |
| IC50 | Half-maximal inhibitory concentration |
| EC50 | Half-maximal effective concentration |
| Kb | Binding constant |
| pKi | -log10(Ki) |
| pIC50 | -log10(IC50) |

**Activity Sources (act_source):**

| Source | Percentage | Description |
|--------|------------|-------------|
| ChEMBLdb | 58.73% | ChEMBL database |
| DrugMatrix | 15.95% | DrugMatrix/Iconix |
| WOMBAT-PK | 13.90% | WOMBAT-PK |
| Guide to Pharmacology | 5.29% | IUPHAR |
| PDSP | 3.41% | NIMH PDSP Ki Database |
| Literature | 2.35% | Primary literature |
| Drug labels | 0.36% | FDA labels |

##### TARGETS

```sql
CREATE TABLE targets (
    id              SERIAL PRIMARY KEY,
    name            VARCHAR(500) NOT NULL,
    target_class    VARCHAR(100),
    accession       VARCHAR(20),
    gene            VARCHAR(50),
    swissprot       VARCHAR(20),
    organism        VARCHAR(100),
    tdl             VARCHAR(10),
    idg_family      VARCHAR(100)
);
```

| Column | Type | Description |
|--------|------|-------------|
| `name` | VARCHAR | Target name |
| `target_class` | VARCHAR | Target classification |
| `accession` | VARCHAR | UniProt accession |
| `gene` | VARCHAR | Gene symbol |
| `tdl` | VARCHAR | Target Development Level (Tclin, Tchem, Tbio, Tdark) |
| `idg_family` | VARCHAR | IDG family classification |

**Target Development Level (TDL):**

| Level | Description |
|-------|-------------|
| Tclin | Clinical - approved drug target |
| Tchem | Chemical - small molecule activity |
| Tbio | Biological - UniProt function annotation |
| Tdark | Dark - little is known |

#### 4.2.3 OMOP Relationship Tables

##### OMOP_RELATIONSHIP (Indications & Contra-indications)

```sql
CREATE TABLE omop_relationship (
    id              SERIAL PRIMARY KEY,
    struct_id       INTEGER REFERENCES structures(id),
    concept_id      INTEGER NOT NULL,
    relationship_name VARCHAR(100) NOT NULL,
    concept_name    VARCHAR(500),
    umls_cui        VARCHAR(20),
    snomed_conceptid BIGINT,
    snomed_full_name TEXT
);
```

| Column | Type | Description |
|--------|------|-------------|
| `struct_id` | INTEGER | Drug structure ID |
| `concept_id` | INTEGER | OMOP concept ID |
| `relationship_name` | VARCHAR | Type of relationship |
| `concept_name` | VARCHAR | Condition name |
| `umls_cui` | VARCHAR | UMLS Concept Unique Identifier |
| `snomed_conceptid` | BIGINT | SNOMED CT concept ID |
| `snomed_full_name` | TEXT | SNOMED CT preferred term |

**Relationship Types:**

| relationship_name | Count (approx.) | Description |
|-------------------|-----------------|-------------|
| indication | 10,707 | Approved indication |
| contraindication | 27,851 | Contraindication |
| off-label use | 2,496 | Off-label indication |

#### 4.2.4 Adverse Event Data (FAERS)

##### FAERS (FDA Adverse Event Reporting System)

```sql
CREATE TABLE faers (
    id              SERIAL PRIMARY KEY,
    struct_id       INTEGER REFERENCES structures(id),
    meddra_code     VARCHAR(20) NOT NULL,
    meddra_name     VARCHAR(500),
    llr             NUMERIC(12,4),
    llr_threshold   NUMERIC(12,4),
    drug_ae_count   INTEGER,
    drug_count      INTEGER,
    ae_count        INTEGER,
    total_count     INTEGER
);
```

| Column | Type | Description |
|--------|------|-------------|
| `struct_id` | INTEGER | Drug structure ID |
| `meddra_code` | VARCHAR | MedDRA Preferred Term code |
| `meddra_name` | VARCHAR | MedDRA adverse event term |
| `llr` | NUMERIC | Log-likelihood ratio statistic |
| `llr_threshold` | NUMERIC | Significance threshold (p<0.05) |
| `drug_ae_count` | INTEGER | Reports with drug and event |
| `drug_count` | INTEGER | Total reports for drug |
| `ae_count` | INTEGER | Total reports for event |
| `total_count` | INTEGER | Total FAERS reports |

**LLR Statistical Method:**

The log-likelihood ratio (LLR) identifies drug-adverse event combinations with disproportionally high reporting rates:

```
LLR = 2 * [a*ln(a/E_a) + b*ln(b/E_b) + c*ln(c/E_c) + d*ln(d/E_d)]

Where:
a = drug_ae_count (reports with drug AND event)
b = ae_count - a (reports with event but NOT drug)
c = drug_count - a (reports with drug but NOT event)
d = total_count - a - b - c (reports with neither)
```

**Age-Specific FAERS Tables:**

| Table | Description |
|-------|-------------|
| `faers` | All ages |
| `faers_ped` | Pediatric (1 day - 17 years) |
| `faers_ger` | Geriatric (>65 years) |
| `faers_female` | Female patients |
| `faers_male` | Male patients |

##### FAERS_MED (FAERS processed with MedDRA hierarchy)

```sql
CREATE TABLE faers_med (
    id              SERIAL PRIMARY KEY,
    struct_id       INTEGER REFERENCES structures(id),
    meddra_code     VARCHAR(20),
    meddra_name     VARCHAR(500),
    meddra_level    VARCHAR(10),
    llr             NUMERIC(12,4),
    prr             NUMERIC(12,4),
    count           INTEGER
);
```

| meddra_level | Description |
|--------------|-------------|
| PT | Preferred Term |
| HLT | High Level Term |
| HLGT | High Level Group Term |
| SOC | System Organ Class |

#### 4.2.5 Additional Tables

##### ATC (Anatomical Therapeutic Chemical Classification)

```sql
CREATE TABLE atc (
    id              SERIAL PRIMARY KEY,
    struct_id       INTEGER REFERENCES structures(id),
    atc_code        VARCHAR(10) NOT NULL,
    atc_name        VARCHAR(500)
);
```

| ATC Level | Example | Description |
|-----------|---------|-------------|
| 1 | N | Nervous system |
| 2 | N02 | Analgesics |
| 3 | N02B | Other analgesics and antipyretics |
| 4 | N02BA | Salicylic acid and derivatives |
| 5 | N02BA01 | Aspirin |

##### APPROVAL (Regulatory Approval)

```sql
CREATE TABLE approval (
    id              SERIAL PRIMARY KEY,
    struct_id       INTEGER REFERENCES structures(id),
    approval_date   DATE,
    approval_type   VARCHAR(100),
    agency          VARCHAR(50)
);
```

| agency | Description |
|--------|-------------|
| FDA | US Food and Drug Administration |
| EMA | European Medicines Agency |
| PMDA | Japan Pharmaceuticals and Medical Devices Agency |
| TGA | Australian Therapeutic Goods Administration |

### 4.3 Database Connection

**Public Access Credentials:**

```yaml
host: unmtid-dbs.net
port: 5433
database: drugcentral
user: drugman
password: dosage
```

**Docker Deployment:**
```bash
docker pull drugcentral/drugcentral
docker run -p 5433:5432 drugcentral/drugcentral
```

### 4.4 Example Queries

#### Get Drug-Target Bioactivity:
```sql
SELECT
    s.name AS drug_name,
    s.smiles,
    a.target_name,
    a.gene,
    a.accession AS uniprot,
    a.act_type,
    a.act_value,
    a.act_unit,
    a.action_type,
    a.act_source
FROM structures s
JOIN act_table_full a ON s.id = a.struct_id
WHERE s.name = 'aspirin'
ORDER BY a.act_value;
```

#### Get Drug Indications:
```sql
SELECT
    s.name AS drug_name,
    o.relationship_name,
    o.concept_name AS indication,
    o.snomed_full_name,
    o.umls_cui
FROM structures s
JOIN omop_relationship o ON s.id = o.struct_id
WHERE s.name = 'aspirin'
  AND o.relationship_name = 'indication';
```

#### Get Significant Adverse Events:
```sql
SELECT
    s.name AS drug_name,
    f.meddra_name AS adverse_event,
    f.meddra_code,
    f.llr,
    f.llr_threshold,
    f.drug_ae_count
FROM structures s
JOIN faers f ON s.id = f.struct_id
WHERE s.name = 'aspirin'
  AND f.llr > f.llr_threshold
ORDER BY f.llr DESC
LIMIT 20;
```

#### Get Cross-References:
```sql
SELECT
    s.name AS drug_name,
    i.id_type,
    i.identifier
FROM structures s
JOIN identifier i ON s.id = i.struct_id
WHERE s.name = 'aspirin'
ORDER BY i.id_type;
```

---

## 5. Cross-Database ID Mapping

### 5.1 UniChem (EBI)

UniChem provides cross-references between chemistry databases using InChI Key matching.

**API:** https://www.ebi.ac.uk/unichem/api/v1/

| Source ID | Database |
|-----------|----------|
| 1 | ChEMBL |
| 2 | DrugBank |
| 3 | PDB |
| 4 | IUPHAR |
| 6 | KEGG Ligand |
| 7 | ChEBI |
| 10 | eMolecules |
| 14 | FDA/SPL |
| 17 | ZINC |
| 22 | PubChem |
| 31 | BindingDB |
| 36 | DrugCentral |

**Example Cross-Reference Query:**
```
GET https://www.ebi.ac.uk/unichem/rest/src_compound_id/CHEMBL25/1/2
```

### 5.2 Common Identifier Mappings

| DrugBank | ChEMBL | PubChem CID | DrugCentral | Name |
|----------|--------|-------------|-------------|------|
| DB00945 | CHEMBL25 | 2244 | 1967 | Aspirin |
| DB00316 | CHEMBL1771 | 5090 | 3696 | Acetaminophen |
| DB01050 | CHEMBL521 | 2519 | 2844 | Ibuprofen |
| DB00688 | CHEMBL34259 | 5281 | 3376 | Mycophenolic acid |
| DB00661 | CHEMBL103 | 54675776 | 5224 | Verapamil |

---

## 6. Summary Comparison

| Feature | DrugBank | ChEMBL | PubChem | DrugCentral |
|---------|----------|--------|---------|-------------|
| **Focus** | Approved drugs | Bioactivity data | Chemical substances | Approved drugs |
| **Primary ID** | DB###### | CHEMBL###### | CID (compound) | struct_id |
| **Format** | XML, CSV | SQL, JSON | JSON, XML, SDF | PostgreSQL |
| **Compounds** | ~15,000 drugs | ~2.4M molecules | ~120M compounds | ~5,000 drugs |
| **Targets** | ~5,000 | ~17,000 | Variable | ~2,200 |
| **Activities** | Qualitative | ~23M quantitative | ~300M | ~14,000 |
| **MoA Data** | Yes | Yes | Limited | Yes |
| **Indications** | Yes | Yes | No | Yes (OMOP) |
| **ADRs** | SNP-ADR | No | No | FAERS LLR |
| **API** | REST (paid) | REST (free) | REST (free) | PostgreSQL |
| **License** | CC BY-NC 4.0 | CC BY-SA 3.0 | Public Domain | CC BY 4.0 |

---

## References

### DrugBank
- Wishart DS, et al. DrugBank 6.0: the DrugBank Knowledgebase for 2024. Nucleic Acids Res. 2024;52(D1):D1265-D1275.
- XML Format Reference: https://docs.drugbank.com/xml/

### ChEMBL
- Zdrazil B, et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Res. 2024;52(D1):D1180-D1192.
- Schema Documentation: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html
- Web Services: https://www.ebi.ac.uk/chembl/api/data/docs

### PubChem
- Kim S, et al. PubChem Substance and Compound databases. Nucleic Acids Res. 2016;44(D1):D1202-D1213.
- PUG REST: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

### DrugCentral
- Avram S, et al. DrugCentral 2023 extends human clinical data and integrates veterinary drugs. Nucleic Acids Res. 2023;51(D1):D1276-D1287.
- Ursu O, et al. DrugCentral: online drug compendium. Nucleic Acids Res. 2017;45(D1):D932-D939.
- Download: https://drugcentral.org/download

### UniChem
- Chambers J, et al. UniChem: a unified chemical structure cross-referencing and identifier tracking system. J Cheminform. 2013;5(1):3.
- API: https://www.ebi.ac.uk/unichem/
