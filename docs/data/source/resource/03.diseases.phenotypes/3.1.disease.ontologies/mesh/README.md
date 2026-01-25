---
id: mesh
title: "Medical Subject Headings (MeSH)"
type: source
parent: ../README.md
tier: 1
status: active
category: diseases.phenotypes
subcategory: disease.ontologies
tags:
  - vocabulary
  - medical-terms
  - nlm
  - pubmed
  - indexing
  - diseases
---

# Medical Subject Headings (MeSH)

## Overview

Medical Subject Headings (MeSH) is the National Library of Medicine's (NLM) controlled vocabulary thesaurus used for indexing articles in PubMed and cataloging books in the NLM catalog. MeSH provides a hierarchically-organized terminology for describing biomedical concepts and allows searching across varying levels of specificity.

MeSH contains over 30,000 descriptors arranged in a hierarchical tree structure with 16 major categories including Anatomy, Diseases, Chemicals and Drugs, and Phenomena and Processes. Each descriptor has a unique identifier and includes synonyms (entry terms), scope notes, and cross-references to related concepts.

For disease research, MeSH provides standardized disease terminology that enables consistent literature searches across millions of biomedical publications. The vocabulary is updated annually with new terms and structural revisions, making it an essential resource for systematic reviews, meta-analyses, and biomedical text mining applications.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Descriptors | ~30,000 |
| Disease Descriptors | ~5,000 |
| Supplementary Concepts | ~280,000 |
| Tree Categories | 16 |
| Annual Updates | Continuous |

## Primary Use Cases

1. PubMed literature indexing and search
2. Medical vocabulary standardization across databases
3. Disease classification cross-referencing
4. Biomedical text mining and NLP applications
5. Systematic review search strategy development

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| MeSH UI | `D[0-9]{6}` | D003920 (Diabetes Mellitus) |
| Tree Number | `C[0-9]{2}(\.[0-9]{3})+` | C18.452.394.750 |
| MeSH Term | Free text | Diabetes Mellitus, Type 2 |
| Entry Term | Free text | Type 2 Diabetes |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| MeSH Browser | https://meshb.nlm.nih.gov/ | Web interface |
| MeSH RDF | https://id.nlm.nih.gov/mesh/ | Linked data |
| FTP Download | https://www.nlm.nih.gov/databases/download/mesh.html | Bulk files |
| SPARQL | https://id.nlm.nih.gov/mesh/sparql | Query endpoint |

## Data Formats

| Format | File | Notes |
|--------|------|-------|
| XML | desc2026.xml | Full descriptors |
| ASCII | d2026.bin | Flat file format |
| RDF/N-Triples | mesh.nt | Linked data |
| JSON-LD | Available via API | Web-friendly |

## Limitations

- Designed for literature indexing; may not map to clinical terminologies
- Annual update cycle means new terms lag publications
- Non-English language support limited
- Supplementary concepts less structured than descriptors

## See Also

- [Schema Definition](./schema.json) - Data structure and field types
- [Field Dictionary](./dictionary.md) - Field semantics and definitions
- [License Terms](./license.md) - Usage rights and restrictions
- [ICD](../icd/README.md) - Clinical classification system
- [MONDO](../mondo/README.md) - Disease ontology with MeSH mappings
