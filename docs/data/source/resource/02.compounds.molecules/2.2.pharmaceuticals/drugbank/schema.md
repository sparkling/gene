---
id: schema-drugbank
title: "DrugBank Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-23
status: migrated
tags: [schema, database, drugs, pharmacology, targets, interactions]
---

# DrugBank - Comprehensive Drug Database Schema

**Document ID:** SCHEMA-DRUGBANK
**Version:** 5.0
**Source Version:** 2024

---

## TL;DR

DrugBank provides comprehensive drug data including chemical properties, pharmacology, targets, interactions, and clinical information. The schema organizes drugs with detailed target proteins (sequences, pathways), drug-drug interactions, and mechanisms of action, supporting drug discovery, clinical practice, and pharmaceutical research.

---

## Database Statistics

| Metric | Value | Source |
|--------|-------|--------|
| Drug Entries | 16,000+ | Drug collection |
| Approved Small Molecules | 2,700+ | Approved drugs |
| Approved Biologics | 1,500+ | Biologic drugs |
| Experimental Drugs | 6,500+ | Investigational |
| Drug Targets | 5,200+ | Protein targets |
| Drug-Drug Interactions | 380,000+ | DDI database |

---

## Entity Relationship Overview

```
Drugs (1) ←→ (many) Drug_Targets (many) ←→ (1) Targets
  ↓                      ↓                      ↓
DrugBank ID         Actions, Ki          UniProt, Sequence

Drugs (1) ←→ (many) Drug_Interactions
                         ↓
                 Interacting drugs, severity

Drugs (1) ←→ (many) Drug_Pathways (many) ←→ (1) Pathways
                         ↓
                   KEGG, SMPDB links

Drugs (1) ←→ (1) Chemical_Properties
                    ↓
             SMILES, InChI, LogP
```

---

## Core Tables/Entities

### drugs

**Description:** Main drug entries with identification and classification.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| drugbank_id | string | Yes | DrugBank identifier (DB00001) |
| name | string | Yes | Drug name |
| type | string | Yes | small molecule, biotech |
| description | text | No | Drug description |
| cas_number | string | No | CAS registry number |
| unii | string | No | FDA UNII identifier |
| state | string | No | solid, liquid, gas |
| indication | text | No | Therapeutic indications |
| pharmacodynamics | text | No | PD description |
| mechanism_of_action | text | No | MOA description |
| toxicity | text | No | Toxicity information |
| metabolism | text | No | Metabolic pathways |
| absorption | text | No | Absorption profile |
| half_life | string | No | Elimination half-life |
| protein_binding | string | No | Protein binding % |
| route_of_elimination | text | No | Elimination routes |
| volume_of_distribution | string | No | Vd value |
| clearance | string | No | Clearance value |
| groups | array | Yes | approved, investigational, etc. |

### chemical_properties

**Description:** Chemical and physical properties of drugs.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| drugbank_id | string | Yes | Foreign key to drugs |
| molecular_formula | string | No | Chemical formula |
| molecular_weight | decimal | No | MW in Daltons |
| monoisotopic_mass | decimal | No | Exact mass |
| smiles | string | No | SMILES string |
| inchi | string | No | InChI identifier |
| inchikey | string | No | InChI Key |
| logp | decimal | No | Partition coefficient |
| pka | string | No | Ionization constant |

### targets

**Description:** Drug target proteins and other biomolecules.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| target_id | string | Yes | Target identifier |
| name | string | Yes | Target name |
| organism | string | Yes | Species |
| gene_name | string | No | Gene symbol |
| uniprot_id | string | No | UniProt accession |
| pdb_ids | array | No | PDB structure IDs |
| amino_acid_sequence | text | No | Protein sequence |

### drug_targets

**Description:** Links drugs to their targets with action types.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| drugbank_id | string | Yes | Foreign key to drugs |
| target_id | string | Yes | Foreign key to targets |
| actions | array | Yes | inhibitor, agonist, etc. |
| known_action | string | No | yes, no, unknown |
| pharmacological_action | string | No | yes, no |
| ki | decimal | No | Inhibition constant |
| kd | decimal | No | Dissociation constant |
| ic50 | decimal | No | IC50 value |
| ec50 | decimal | No | EC50 value |
| references | array | No | PubMed IDs |

### drug_interactions

**Description:** Drug-drug interactions.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| drugbank_id | string | Yes | First drug |
| interacting_drug_id | string | Yes | Interacting drug |
| description | text | Yes | Interaction description |
| severity | string | No | major, moderate, minor |

### pathways

**Description:** Metabolic and signaling pathways.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| smpdb_id | string | Yes | SMPDB identifier |
| name | string | Yes | Pathway name |
| category | string | No | Metabolic, signaling, etc. |
| drugs | array | Yes | Associated DrugBank IDs |

---

## API Endpoints (Commercial)

| Endpoint | Method | Description |
|----------|--------|-------------|
| /drugs/{id} | GET | Get drug by DrugBank ID |
| /drugs/search | GET | Search drugs |
| /targets/{id} | GET | Get target information |
| /interactions/{id} | GET | Get drug interactions |

---

## Data Formats

| Format | Description |
|--------|-------------|
| Primary | XML (full database) |
| Alternative | CSV (selected fields) |
| API Response | JSON |
| Encoding | UTF-8 |

---

## Sample Record

```xml
<drug type="small molecule" created="2005-06-13" updated="2024-01-02">
  <drugbank-id primary="true">DB00945</drugbank-id>
  <name>Aspirin</name>
  <description>
    Aspirin is an analgesic, antipyretic, and anti-inflammatory agent...
  </description>
  <cas-number>50-78-2</cas-number>
  <groups>
    <group>approved</group>
  </groups>
  <indication>
    For the treatment of pain, fever, and inflammation...
  </indication>
  <calculated-properties>
    <property>
      <kind>SMILES</kind>
      <value>CC(=O)Oc1ccccc1C(=O)O</value>
    </property>
    <property>
      <kind>InChIKey</kind>
      <value>BSYNRYMUTXBXSQ-UHFFFAOYSA-N</value>
    </property>
  </calculated-properties>
  <targets>
    <target>
      <id>BE0000017</id>
      <name>Prostaglandin G/H synthase 1</name>
      <organism>Humans</organism>
      <actions>
        <action>inhibitor</action>
      </actions>
      <polypeptide id="P23219" source="Swiss-Prot"/>
    </target>
  </targets>
</drug>
```

---

## Drug Groups

| Group | Description |
|-------|-------------|
| approved | FDA/EMA approved |
| investigational | Clinical trials |
| experimental | Preclinical |
| withdrawn | Removed from market |
| nutraceutical | Nutritional products |
| vet_approved | Veterinary use |

---

## Target Actions

| Action | Description |
|--------|-------------|
| inhibitor | Blocks target activity |
| agonist | Activates receptor |
| antagonist | Blocks receptor |
| substrate | Metabolized by enzyme |
| inducer | Increases expression |
| binder | Binds without action |

---

## Glossary

| Term | Definition |
|------|------------|
| DrugBank ID | Identifier format DB + 5 digits |
| MOA | Mechanism of Action |
| DDI | Drug-Drug Interaction |
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity |
| Ki | Inhibition constant |

---

## References

1. DrugBank: https://go.drugbank.com
2. Wishart DS, et al. (2018) Nucleic Acids Res. 46(D1):D1074-D1082
3. API Documentation: https://docs.drugbank.com
