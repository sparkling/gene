---
id: schema-hit2
title: "HIT 2.0 Database Schema"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [schema, herbal-medicine, validated-targets, experimental, binding-affinity]
---

# HIT 2.0 Database Schema

**Document ID:** HIT2-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** HIT 2.0

---

## TL;DR

HIT 2.0 (Herbal Ingredients' Targets) is the gold standard database for experimentally validated interactions between herbal ingredients and molecular targets. Contains 10,000+ herbal ingredients, 5,000+ validated targets, 100,000+ herb-target interactions, 50,000+ literature references, and 30,000+ binding affinity data points. Unlike prediction-based databases, HIT 2.0 exclusively includes interactions supported by experimental evidence.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **Herbal Ingredients** | 10,000+ |
| **Validated Targets** | 5,000+ |
| **Herb-Target Interactions** | 100,000+ |
| **Literature References** | 50,000+ |
| **Binding Affinity Data** | 30,000+ |
| **Traditional Systems** | Multiple |

---

## Data Model

### Entity Relationships

```
Herbal Sources (Multiple Systems)
+-- TCM Herbs
+-- Ayurvedic Plants
+-- Western Herbs
            |
            v
   Ingredients (10,000+)
            |
            v (experimental validation)
   Target Proteins (5,000+)
            |
            +-- Binding Affinity (30,000+)
            +-- Literature (50,000+)
            +-- Experimental Method
```

### Ingredient Entity (Primary)

```json
{
  "ingredient_id": "string (HIT internal, e.g., HIT001234)",
  "ingredient_name": "string",
  "identifiers": {
    "pubchem_cid": "integer",
    "chembl_id": "string",
    "cas_number": "string",
    "inchi_key": "string"
  },
  "structure": {
    "smiles": "string",
    "molecular_formula": "string",
    "molecular_weight": "float"
  },
  "herbal_sources": [
    {
      "herb_name": "string",
      "traditional_system": "string (TCM|Ayurveda|Kampo|Western)",
      "plant_part": "string"
    }
  ],
  "validated_targets": ["target_id array"],
  "interaction_count": "integer",
  "binding_affinity_count": "integer"
}
```

### Interaction Entity (Core Feature)

```json
{
  "interaction_id": "string (HIT internal, e.g., HIT_INT_001234)",
  "ingredient_id": "string",
  "target_id": "string",
  "interaction_type": "string (inhibitor|agonist|antagonist|modulator|binder)",
  "evidence": {
    "experimental_method": "string",
    "assay_type": "string",
    "binding_affinity": {
      "metric": "string (IC50|Ki|Kd|EC50|KD)",
      "value": "float",
      "unit": "string (nM|uM|mM)",
      "relation": "string (=|<|>|~)"
    },
    "confidence_level": "string (high|medium|low)"
  },
  "literature": {
    "pubmed_id": "integer",
    "doi": "string",
    "authors": "string",
    "journal": "string",
    "year": "integer"
  }
}
```

### Target Entity

```json
{
  "target_id": "string (UniProt accession or Gene Symbol)",
  "gene_symbol": "string",
  "gene_name": "string",
  "uniprot_id": "string",
  "protein_class": "string (Kinase|GPCR|Nuclear Receptor|Enzyme|Channel|etc.)",
  "organism": "string (Human|Mouse|Rat)",
  "pathway_associations": ["pathway_id array"],
  "disease_associations": ["disease_id array"],
  "ingredient_count": "integer",
  "interaction_count": "integer"
}
```

### Herbal Source Entity

```json
{
  "herb_id": "string",
  "herb_name": "string",
  "botanical_name": "string",
  "traditional_system": "string",
  "region": "string",
  "plant_parts": ["string array"],
  "ingredient_count": "integer"
}
```

---

## Evidence Types

### Experimental Methods

| Method | Description | Confidence |
|--------|-------------|------------|
| Binding Assay | Direct binding measurement (radioligand, SPR) | High |
| Enzyme Inhibition | IC50, Ki determination | High |
| Cell-Based Assay | Functional cellular readout | Medium |
| Receptor Binding | Ligand-receptor interaction | High |
| In Vivo | Animal model confirmation | Medium |
| Clinical | Human studies | Variable |

### Binding Affinity Metrics

| Metric | Full Name | Description |
|--------|-----------|-------------|
| IC50 | Half-maximal Inhibitory Concentration | Concentration for 50% inhibition |
| Ki | Inhibition Constant | Equilibrium dissociation constant for inhibitor |
| Kd | Dissociation Constant | Binding affinity constant |
| EC50 | Half-maximal Effective Concentration | Concentration for 50% effect |
| KD | Equilibrium Dissociation Constant | SPR-derived binding constant |

### Interaction Types

| Type | Description |
|------|-------------|
| Inhibitor | Reduces target activity |
| Agonist | Activates target |
| Antagonist | Blocks target activation |
| Modulator | Alters target function |
| Binder | Binds without known functional effect |
| Substrate | Acts as enzyme substrate |

---

## Multi-System Coverage

| Traditional System | Coverage Level | Examples |
|--------------------|----------------|----------|
| Traditional Chinese Medicine | Comprehensive | Ginseng, Astragalus |
| Ayurveda | Major herbs | Ashwagandha, Turmeric |
| Japanese Kampo | Selected | Glycyrrhiza, Ephedra |
| Western Herbal | Common herbs | St. John's Wort, Echinacea |
| African Traditional | Limited | Selected compounds |

---

## Validation Use Case

### Cross-Reference Workflow

```
Predicted Interaction (from BATMAN-TCM, KampoDB, IMPPAT)
                    |
                    v
            Query HIT 2.0
                    |
         +----------+----------+
         |                     |
         v                     v
    Found: VALIDATED      Not Found: UNVALIDATED
         |                     |
         v                     v
   Use with confidence   Requires experimental testing
```

### Integration Example

```python
import pandas as pd

# Load predicted interactions (e.g., from BATMAN-TCM)
predicted = pd.read_csv('batman_predicted_tti.tsv', sep='\t')

# Load HIT 2.0 validated interactions
validated = pd.read_csv('hit2_interactions.tsv', sep='\t')

# Create lookup set of validated pairs
validated_pairs = set(
    zip(validated['ingredient_pubchem_cid'],
        validated['target_uniprot'])
)

# Mark predictions as validated or not
predicted['experimentally_validated'] = predicted.apply(
    lambda row: (row['pubchem_cid'], row['uniprot_id']) in validated_pairs,
    axis=1
)

# Summary statistics
print(f"Predictions validated: {predicted['experimentally_validated'].sum()}")
print(f"Predictions unvalidated: {(~predicted['experimentally_validated']).sum()}")
print(f"Validation rate: {100*predicted['experimentally_validated'].mean():.1f}%")
```

---

## Cross-References

### Identifier Mappings

| HIT Field | External Database | Notes |
|-----------|-------------------|-------|
| pubchem_cid | PubChem | Compound structure |
| chembl_id | ChEMBL | Bioactivity data |
| uniprot_id | UniProt | Protein sequence |
| pubmed_id | PubMed | Literature reference |
| kegg_pathway | KEGG | Pathway annotation |

---

## Sample Data

### Sample Interaction Record

| Field | Value |
|-------|-------|
| interaction_id | HIT_INT_000001 |
| ingredient_name | Curcumin |
| target_gene | TP53 |
| interaction_type | Modulator |
| affinity_type | Ki |
| affinity_value | 250 |
| affinity_unit | nM |
| pubmed_id | 12345678 |

### Sample Binding Affinity Distribution

| Affinity Range | Interaction Count |
|----------------|-------------------|
| < 100 nM | 5,000+ |
| 100 nM - 1 uM | 12,000+ |
| 1 uM - 10 uM | 10,000+ |
| > 10 uM | 3,000+ |

---

## Integration Recommendations

### Primary Use: Validation

HIT 2.0 should be used to validate predictions from:

| Database | Validation Purpose |
|----------|---------------------|
| BATMAN-TCM | Confirm predicted TTIs |
| KampoDB | Validate docking predictions |
| IMPPAT | Cross-reference STITCH predictions |
| SymMap | Ground symptom-target links |
| TCMSID | Validate filtered compounds |

### Analysis Workflow

```
1. Generate predictions (BATMAN-TCM, KampoDB, etc.)
2. Query HIT 2.0 for experimental validation
3. Separate validated vs unvalidated predictions
4. Prioritize validated interactions for downstream analysis
5. Flag unvalidated high-confidence predictions for experimental testing
```

---

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

---

## Key Strengths

| Strength | Description |
|----------|-------------|
| Experimental Validation | No predictions, only confirmed interactions |
| Multi-System | Spans multiple traditional medicine systems |
| Quantitative | Binding affinity data included |
| Literature-Linked | Full citation support |
| Gold Standard | Reference for validating other databases |

---

## Limitations

1. Literature curation lag behind publications
2. Not all interactions have binding affinity data
3. Coverage biased toward well-studied compounds
4. English literature bias
5. Validation coverage varies by target class

---

## Glossary

| Term | Definition |
|------|------------|
| Validated Interaction | Compound-target relationship confirmed by experiment |
| Binding Affinity | Quantitative measure of interaction strength |
| IC50 | Concentration for 50% inhibition |
| Ki | Inhibition constant (binding affinity for inhibitors) |
| Kd | Dissociation constant (general binding affinity) |
| Gold Standard | Reference database for validation purposes |

---

## References

1. Ye H, et al. (2011). HIT: linking herbal active ingredients to targets. Nucleic Acids Research.
2. HIT 2.0 Official Website: http://hit2.badd-cao.net/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
