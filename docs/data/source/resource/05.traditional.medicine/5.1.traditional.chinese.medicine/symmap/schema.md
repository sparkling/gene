---
id: schema-symmap
title: "SymMap 2.0 Database Schema"
type: schema
parent: README.md
last_updated: 2026-01-23
status: active
tags: [schema, tcm, symptom-mapping, phenotype, network-medicine]
---

# SymMap 2.0 Database Schema

**Document ID:** SYMMAP-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** SymMap 2.0

---

## TL;DR

SymMap 2.0 uniquely bridges TCM symptoms to modern phenotypes via Human Phenotype Ontology (HPO) mappings. Contains 1,717 TCM symptoms, 5,000+ HPO phenotypes, 961 herbs, 19,595 compounds, 4,302 targets, 5,235 diseases, 28,212 symptom-herb associations, and 322,073 symptom-gene associations. Enables phenotype-driven analysis of traditional medicine.

---

## Database Statistics

| Entity | Count |
|--------|-------|
| **TCM Symptoms** | 1,717 |
| **Modern Phenotypes (HPO)** | 5,000+ |
| **TCM Herbs** | 961 |
| **Compounds** | 19,595 |
| **Target Genes** | 4,302 |
| **Diseases** | 5,235 |
| **Symptom-Herb Associations** | 28,212 |
| **Symptom-Gene Associations** | 322,073 |

---

## Data Model

### Entity Relationships

```
TCM Symptoms (1,717) <----> Modern Phenotypes (HPO)
        |                           |
        v                           v
    Herbs (961)              Diseases (5,235)
        |                           |
        v                           v
  Compounds (19,595) --------> Targets (4,302)
```

### Symptom Entity (Core)

```json
{
  "symptom_id": "string (SMSY internal, e.g., SMSY00001)",
  "symptom_name_chinese": "string",
  "symptom_name_english": "string",
  "symptom_category": "string (TCM classification)",
  "tcm_pattern": "string (Zheng pattern)",
  "hpo_mappings": [
    {
      "hpo_id": "string (HP:XXXXXXX)",
      "hpo_term": "string",
      "mapping_confidence": "float"
    }
  ],
  "associated_herbs": ["herb_id array"],
  "associated_genes": ["gene_symbol array"],
  "associated_diseases": ["disease_id array"]
}
```

### HPO Phenotype Entity

```json
{
  "hpo_id": "string (HP:XXXXXXX)",
  "hpo_term": "string",
  "definition": "string",
  "synonyms": ["string array"],
  "parent_terms": ["hpo_id array"],
  "child_terms": ["hpo_id array"],
  "tcm_symptom_mappings": ["symptom_id array"],
  "disease_associations": ["disease_id array"]
}
```

### Herb Entity

```json
{
  "herb_id": "string (SMHB internal, e.g., SMHB00001)",
  "herb_name_chinese": "string",
  "herb_name_pinyin": "string",
  "herb_name_english": "string",
  "herb_name_latin": "string",
  "tcm_properties": {
    "nature": "string",
    "flavor": ["string array"],
    "meridian_tropism": ["string array"]
  },
  "symptom_associations": ["symptom_id array"],
  "compound_count": "integer"
}
```

### Compound Entity

```json
{
  "compound_id": "string (SMCP internal, e.g., SMCP00001)",
  "pubchem_cid": "integer",
  "compound_name": "string",
  "smiles": "string",
  "molecular_weight": "float",
  "source_herbs": ["herb_id array"],
  "target_genes": ["gene_symbol array"]
}
```

### Target Entity

```json
{
  "gene_symbol": "string",
  "gene_name": "string",
  "entrez_gene_id": "integer",
  "uniprot_id": "string",
  "symptom_associations": ["symptom_id array"],
  "compound_associations": ["compound_id array"],
  "pathway_annotations": ["pathway_id array"]
}
```

### Disease Entity

```json
{
  "disease_id": "string (MESH or OMIM)",
  "disease_name": "string",
  "mesh_id": "string",
  "omim_id": "integer",
  "symptom_associations": ["symptom_id array"],
  "gene_associations": ["gene_symbol array"],
  "hpo_phenotypes": ["hpo_id array"]
}
```

---

## Symptom-Phenotype Mapping

### TCM to HPO Bridging

```
TCM Symptom (Traditional)
    |
    v (manual curation + text mining)
HPO Phenotype (Modern Medical)
    |
    v
Genetic/Disease Associations
```

### Mapping Types

| Type | Description | Count |
|------|-------------|-------|
| Direct | Clear correspondence | ~60% |
| Partial | Related but not identical | ~30% |
| Composite | TCM maps to multiple HPO | ~10% |

### Example Mappings

| TCM Symptom | HPO Phenotype(s) |
|-------------|------------------|
| Headache (Tou Tong) | HP:0002315 Headache |
| Insomnia (Shi Mian) | HP:0100785 Insomnia |
| Fatigue (Pi Lao) | HP:0012378 Fatigue |
| Blood Stasis Pattern | HP:0001892 Abnormal bleeding + others |

---

## Network Layers

### Multi-Layer Network Structure

```
Layer 1: TCM Symptoms (1,717)
    |
    v (mapping)
Layer 2: HPO Phenotypes (5,000+)
    |
    v (association)
Layer 3: Disease/Gene Networks
    |
    v (therapeutic)
Layer 4: Herb-Compound Networks
```

### Traversal Paths

| From | To | Via | Use Case |
|------|-----|-----|----------|
| Symptom | Herb | Direct association | Find herbs for symptoms |
| Symptom | Gene | HPO mapping | Molecular basis of symptoms |
| HPO | TCM Herb | Symptom bridge | Modern-to-traditional discovery |
| Disease | Herb | Gene-compound | Drug repurposing |

---

## Cross-References

### Identifier Mappings

| SymMap Field | External Database | Notes |
|--------------|-------------------|-------|
| hpo_id | Human Phenotype Ontology | HP:XXXXXXX format |
| compound_id | PubChem CID | Via structure |
| gene_symbol | HGNC | Standard symbols |
| disease_id | MESH / OMIM | Disease concepts |

### External Links

| Database | Purpose |
|----------|---------|
| HPO | Phenotype ontology |
| PubChem | Compound structures |
| KEGG | Pathway annotations |
| OMIM | Genetic disease associations |
| MESH | Medical subject headings |

---

## Sample Data

### Sample Symptom Record

| Field | Value |
|-------|-------|
| symptom_id | SMSY00001 |
| symptom_name_english | Headache |
| hpo_mapping | HP:0002315 |
| associated_herbs | 45 |
| associated_genes | 234 |

### Sample Symptom-Gene Association

| Symptom | Gene | Evidence |
|---------|------|----------|
| SMSY00001 (Headache) | CALCA | HPO mapping |
| SMSY00001 (Headache) | CGRP | Literature |
| SMSY00002 (Insomnia) | MELATONIN | HPO mapping |

---

## Unique Features

### Phenotype-Driven Analysis

| Feature | Description |
|---------|-------------|
| TCM-to-Modern Bridge | First systematic mapping |
| Symptom-based Drug Discovery | Find herbs by phenotype |
| Molecular Mechanism | Link symptoms to genes |
| Cross-cultural Validation | Compare TCM with HPO |

### Network Analysis Capabilities

- Symptom-gene association networks
- Herb-compound-target networks
- Disease-symptom clustering
- Cytoscape visualization support

---

## Integration Recommendations

### Complementary Databases

| Database | Integration Purpose |
|----------|---------------------|
| BATMAN-TCM | Add compound-target predictions |
| HERB | Add gene expression validation |
| DisGeNET | Extend disease-gene associations |
| HPO | Full ontology hierarchy |

### Analysis Workflow

```
1. Query SymMap for TCM symptom pattern
2. Map to HPO phenotypes
3. Identify associated genes
4. Find therapeutic herbs via compounds
5. Validate with expression data (HERB)
6. Confirm targets (HIT 2.0)
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

1. Symptom-phenotype mappings may be approximate
2. Not all TCM concepts map cleanly to HPO
3. Manual curation limits coverage expansion
4. Western medical bias in HPO may miss TCM nuances

---

## Glossary

| Term | Definition |
|------|------------|
| TCM Symptom | Traditional diagnostic sign/symptom |
| Zheng | TCM pattern/syndrome |
| HPO | Human Phenotype Ontology |
| Phenotype | Observable characteristic |
| Meridian Tropism | TCM organ system affinity |
| Phenotype-driven | Analysis starting from symptoms |

---

## References

1. Wu Y, et al. (2019). SymMap: an integrative database of traditional Chinese medicine enhanced by symptom mapping. Nucleic Acids Research.
2. Official Website: http://www.symmap.org/
3. HPO: https://hpo.jax.org/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |
