---
id: schema-herb
title: "HERB Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: active
tags: [schema, tcm, gene-expression, transcriptomics]
---

# HERB Database Schema

**Document ID:** HERB-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** HERB 2.0

---

## TL;DR

HERB (High-throughput Experiment and Reference Database for TCM) uniquely integrates gene expression data with TCM information. Contains 1,037 herbs, 12,933 ingredients, 2,064 targets, 866 diseases, and 2,000+ gene expression experiments from GEO. Provides transcriptomic signatures for herbs and connectivity analysis with disease signatures.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **TCM Herbs** | 1,037 |
| **Ingredients (Compounds)** | 12,933 |
| **Target Genes** | 2,064 |
| **Diseases** | 866 |
| **Gene Expression Experiments** | 2,000+ |
| **Herb-Ingredient Links** | 49,258 |
| **Ingredient-Target Links** | 28,212 |

---

## Data Model

### Entity Relationships

```
Gene Expression Experiments (GEO)
            |
            v
+-----> Herbs (1,037) <-----> Diseases (866)
|           |
|           v
|    Ingredients (12,933)
|           |
|           v
+------ Targets (2,064) <---- Experimental Validation
```

### Herb Entity

```json
{
  "herb_id": "string (HERB internal ID)",
  "herb_name_chinese": "string",
  "herb_name_pinyin": "string",
  "herb_name_english": "string",
  "herb_name_latin": "string (botanical name)",
  "tcm_properties": {
    "nature": "cold|cool|neutral|warm|hot",
    "flavor": ["sweet", "bitter", "sour", "salty", "pungent"],
    "meridian_tropism": ["lung", "heart", "liver", "spleen", "kidney"]
  },
  "ingredient_count": "integer",
  "expression_experiments": ["GEO accession array"]
}
```

### Ingredient Entity

```json
{
  "ingredient_id": "string (HERB internal ID)",
  "compound_name": "string",
  "pubchem_cid": "integer",
  "smiles": "string",
  "inchi_key": "string",
  "molecular_formula": "string",
  "molecular_weight": "float",
  "source_herbs": ["herb_id array"],
  "target_genes": ["gene_symbol array"]
}
```

### Target Entity

```json
{
  "target_id": "string (Gene Symbol or Entrez ID)",
  "gene_symbol": "string",
  "entrez_gene_id": "integer",
  "gene_name": "string",
  "uniprot_id": "string",
  "source_ingredients": ["ingredient_id array"],
  "expression_data": {
    "upregulated_by": ["herb_id array"],
    "downregulated_by": ["herb_id array"]
  }
}
```

### Expression Experiment Entity

```json
{
  "experiment_id": "string (GEO accession)",
  "geo_accession": "string (GSE format)",
  "platform": "string (GPL format)",
  "organism": "string",
  "herb_or_compound": "string",
  "treatment_description": "string",
  "differential_genes": {
    "upregulated": ["gene_symbol array"],
    "downregulated": ["gene_symbol array"]
  },
  "fold_change_threshold": "float",
  "significance_threshold": "float (p-value)"
}
```

### Disease Entity

```json
{
  "disease_id": "string (internal or MESH)",
  "disease_name": "string",
  "mesh_id": "string",
  "disease_signature": {
    "upregulated_genes": ["gene_symbol array"],
    "downregulated_genes": ["gene_symbol array"]
  },
  "connected_herbs": ["herb_id array"],
  "connectivity_scores": [
    {
      "herb_id": "string",
      "score": "float"
    }
  ]
}
```

---

## Expression Data Features

| Feature | Description |
|---------|-------------|
| Differential Expression | Up/down-regulated genes per treatment |
| Connectivity Scores | Similarity to disease signatures |
| Pathway Enrichment | KEGG/GO enrichment results |
| Dose-Response | Multiple concentration data (when available) |

### Connectivity Analysis

```
Herb Treatment Signature
    |
    v (comparison)
Disease Gene Signature
    |
    v
Connectivity Score (-1 to +1)
    - Positive: Similar expression pattern
    - Negative: Opposite expression pattern (potential therapeutic)
```

---

## Cross-References

### Identifier Mappings

| HERB Field | External Database | Notes |
|------------|-------------------|-------|
| herb_id | Internal | HERB-specific |
| ingredient_id | PubChem CID | Via structure |
| target_id | Entrez Gene ID | Primary identifier |
| experiment_id | GEO Accession | GSE prefix |
| disease_id | MESH | Disease concepts |

### External Links

| Database | Purpose |
|----------|---------|
| GEO | Gene expression data |
| PubChem | Compound structures |
| KEGG | Pathway annotations |
| Gene Ontology | Functional terms |
| MESH | Disease terminology |

---

## Data Quality

### Expression Data Standards

| Metric | Criteria |
|--------|----------|
| Platform | Microarray or RNA-seq |
| Fold Change | Typically >= 2 |
| P-value | Typically <= 0.05 |
| Replicates | Minimum 3 biological |

### Curation Standards

- Gene expression linked to GEO accessions
- Herb identities verified against pharmacopoeia
- Compound structures standardized via PubChem

---

## Sample Data

### Sample Herb Record

| Field | Value |
|-------|-------|
| herb_id | HERB001234 |
| herb_name_chinese | Example |
| herb_name_pinyin | li-zi |
| ingredient_count | 45 |
| experiments | 12 |

### Sample Expression Record

| Field | Value |
|-------|-------|
| geo_accession | GSE12345 |
| herb_treatment | Ginseng extract |
| upregulated_genes | 234 |
| downregulated_genes | 189 |

---

## Integration Recommendations

### Complementary Databases

| Database | Integration Purpose |
|----------|---------------------|
| BATMAN-TCM | Validate predicted targets with expression |
| SymMap | Add expression data to symptom mapping |
| DisGeNET | Connect herb effects to disease genes |
| CTD | Compare with chemical-gene interactions |

### Analysis Workflow

```
1. Query HERB for herb expression signature
2. Compare with disease signature
3. Identify therapeutic candidates (negative connectivity)
4. Validate targets via BATMAN-TCM predictions
5. Cross-reference with HIT 2.0 for experimental validation
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

1. Gene expression biased toward commonly studied herbs
2. Not all herbs have associated experiments
3. Western experimental conditions may not reflect traditional usage
4. Expression data quality varies by study

---

## Glossary

| Term | Definition |
|------|------------|
| GEO Accession | Gene Expression Omnibus identifier (GSE format) |
| Connectivity Score | Measure of expression signature similarity |
| Differential Expression | Genes significantly changed by treatment |
| Fold Change | Ratio of expression levels (treated/control) |
| Transcriptomic Signature | Set of differentially expressed genes |

---

## References

1. Fang S, et al. (2021). HERB: a high-throughput experiment- and reference-guided database of traditional Chinese medicine. Nucleic Acids Research.
2. Official Website: http://herb.ac.cn/
3. GEO: https://www.ncbi.nlm.nih.gov/geo/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
