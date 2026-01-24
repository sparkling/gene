---
id: schema-tcmsid
title: "TCMSID Database Schema"
type: schema
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [schema, tcm, adme, pharmacokinetics, systems-pharmacology]
---

# TCMSID Database Schema

**Document ID:** TCMSID-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** TCMSID

---

## TL;DR

TCMSID (Traditional Chinese Medicine Systematic pharmacology database and analysis platform) provides comprehensive ADME/T properties for TCM compounds. Contains 499 herbs, 29,384 compounds with oral bioavailability and drug-likeness scores, 3,311 targets, 837 diseases, and 98,215 compound-target pairs. Emphasizes pharmacokinetic filtering for drug discovery.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **TCM Herbs** | 499 |
| **Compounds** | 29,384 |
| **Target Proteins** | 3,311 |
| **Diseases** | 837 |
| **Compound-Target Pairs** | 98,215 |
| **OB Data** | All compounds |
| **DL Data** | All compounds |

---

## Data Model

### Entity Relationships

```
Herbs (499)
    |
    v
Compounds (29,384) ---> ADME Properties
    |                   - Oral Bioavailability (OB)
    |                   - Drug-Likeness (DL)
    |                   - Caco-2 Permeability
    |                   - BBB Permeability
    v
Targets (3,311) -----> Diseases (837)
    |
    v
Pathway Networks
```

### Compound Entity (Primary Focus)

```json
{
  "compound_id": "string (MOL prefix, e.g., MOL000001)",
  "compound_name": "string",
  "identifiers": {
    "pubchem_cid": "integer",
    "cas_number": "string"
  },
  "structure": {
    "smiles": "string",
    "molecular_formula": "string",
    "molecular_weight": "float"
  },
  "adme_properties": {
    "ob": "float (Oral Bioavailability %)",
    "dl": "float (Drug-Likeness 0-1)",
    "caco2": "float (Caco-2 permeability)",
    "bbb": "float (Blood-Brain Barrier)",
    "hl": "string (Half-Life: Long/Short)"
  },
  "lipinski": {
    "mw_ok": "boolean",
    "logp_ok": "boolean",
    "hbd_ok": "boolean",
    "hba_ok": "boolean",
    "violations": "integer (0-4)"
  },
  "source_herbs": ["herb_id array"],
  "predicted_targets": ["target_id array"]
}
```

### ADME Property Definitions

```json
{
  "ob_definition": {
    "name": "Oral Bioavailability",
    "description": "Fraction of oral dose reaching systemic circulation",
    "unit": "percent",
    "threshold": 30,
    "interpretation": ">=30% considered favorable for oral drugs"
  },
  "dl_definition": {
    "name": "Drug-Likeness",
    "description": "Quantitative estimate of similarity to known drugs",
    "unit": "score 0-1",
    "threshold": 0.18,
    "interpretation": ">=0.18 considered drug-like"
  },
  "caco2_definition": {
    "name": "Caco-2 Permeability",
    "description": "Predicted intestinal epithelial permeability",
    "unit": "log value",
    "threshold": -0.4,
    "interpretation": ">=-0.4 indicates good intestinal absorption"
  },
  "bbb_definition": {
    "name": "Blood-Brain Barrier Permeability",
    "description": "Predicted ability to cross BBB",
    "unit": "log value",
    "threshold": -0.3,
    "interpretation": ">=-0.3 suggests BBB penetration"
  },
  "hl_definition": {
    "name": "Half-Life",
    "description": "Pharmacokinetic elimination half-life",
    "values": ["Long", "Short"],
    "interpretation": "Long HL suitable for once-daily dosing"
  }
}
```

### Herb Entity

```json
{
  "herb_id": "string (TCMSID internal)",
  "herb_name_chinese": "string",
  "herb_name_pinyin": "string",
  "herb_name_english": "string",
  "herb_name_latin": "string",
  "compound_count": "integer",
  "drug_like_compounds": "integer (passing OB>=30, DL>=0.18)"
}
```

### Target Entity

```json
{
  "target_id": "string (UniProt or Gene Symbol)",
  "gene_symbol": "string",
  "gene_name": "string",
  "uniprot_id": "string",
  "source_compounds": ["compound_id array"],
  "disease_associations": ["disease_id array"]
}
```

### Disease Entity

```json
{
  "disease_id": "string",
  "disease_name": "string",
  "associated_targets": ["target_id array"],
  "associated_herbs": ["herb_id array"]
}
```

---

## Filtering Workflow

### Standard TCMSID Filtering

| Step | Filter | Criteria | Purpose |
|------|--------|----------|---------|
| 1 | OB | >= 30% | Adequate oral absorption |
| 2 | DL | >= 0.18 | Drug-like properties |
| 3 | Caco-2 | >= -0.4 | Intestinal permeability |
| 4 | Target | Has prediction | Confirmed interactions |

### Example Filter Code

```python
import pandas as pd

# Load TCMSID compound data
compounds = pd.read_csv('tcmsid_compounds.csv')

# Standard TCMSID filtering
drug_like = compounds[
    (compounds['OB'] >= 30) &        # Oral bioavailability
    (compounds['DL'] >= 0.18) &      # Drug-likeness
    (compounds['Caco2'] >= -0.4)     # Intestinal permeability
]

print(f"Drug-like compounds: {len(drug_like)} of {len(compounds)}")
print(f"Percentage passing: {100*len(drug_like)/len(compounds):.1f}%")

# Further filter for CNS targets
cns_candidates = drug_like[drug_like['BBB'] >= -0.3]
print(f"CNS candidates: {len(cns_candidates)}")
```

---

## Network Pharmacology Integration

### Pathway Analysis

```
Drug-like Compounds (filtered)
         |
         v
    Target Proteins
         |
         v
    KEGG Pathways
         |
         v
    Disease Networks
```

### Enrichment Output

```json
{
  "herb_id": "TCMSID001",
  "drug_like_compounds": 45,
  "targets": 234,
  "enriched_pathways": [
    {
      "pathway_id": "hsa04151",
      "pathway_name": "PI3K-Akt signaling pathway",
      "p_value": 0.0001,
      "gene_count": 23
    }
  ],
  "enriched_diseases": [
    {
      "disease_id": "DOID:12345",
      "disease_name": "Example disease",
      "gene_overlap": 15
    }
  ]
}
```

---

## Cross-References

### Identifier Mappings

| TCMSID Field | External Database | Notes |
|--------------|-------------------|-------|
| compound_id | PubChem CID | Via structure |
| target_id | UniProt | Protein accession |
| target_id | Gene Symbol | HGNC standard |
| pathway_id | KEGG | hsa prefix |

---

## Sample Data

### Sample Compound Record

| Field | Value |
|-------|-------|
| compound_id | MOL000001 |
| compound_name | Quercetin |
| OB | 46.43% |
| DL | 0.28 |
| Caco-2 | 0.05 |
| BBB | -0.77 |
| HL | Long |

### Sample Filtering Result

| Herb | Total Compounds | OB>=30 | OB>=30 & DL>=0.18 |
|------|-----------------|--------|-------------------|
| Herb A | 156 | 67 | 23 |
| Herb B | 89 | 34 | 12 |
| Herb C | 234 | 98 | 41 |

---

## Integration Recommendations

### Complementary Databases

| Database | Purpose |
|----------|---------|
| BATMAN-TCM | Target prediction comparison |
| TCMBank | Additional compound structures |
| HERB | Expression validation |
| HIT 2.0 | Experimental target confirmation |

### Systems Pharmacology Workflow

```
1. Filter TCMSID compounds (OB, DL criteria)
2. Identify predicted targets
3. Map to KEGG pathways
4. Construct compound-target-pathway network
5. Identify hub targets/pathways
6. Validate with HIT 2.0 experimental data
```

---

## License

| Aspect | Value |
|--------|-------|
| License | Academic use free |
| Commercial Use | Contact maintainers |
| Attribution | Citation required |

---

## Limitations

1. ADME predictions are computational models
2. Not all herbs have comprehensive compound data
3. OB thresholds may exclude valid candidates
4. Target predictions require experimental validation

---

## Glossary

| Term | Definition |
|------|------------|
| OB | Oral Bioavailability - fraction reaching systemic circulation |
| DL | Drug-Likeness - similarity to known drugs |
| Caco-2 | Human intestinal epithelial cell permeability model |
| BBB | Blood-Brain Barrier permeability |
| HL | Half-Life - pharmacokinetic elimination parameter |
| Systems Pharmacology | Network-based analysis of multi-target drug effects |

---

## References

1. Ru J, et al. (2014). TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. Journal of Cheminformatics.
2. Official Website: http://lsp.nwu.edu.cn/tcmsid.php

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
