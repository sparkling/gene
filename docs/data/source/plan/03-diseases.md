# 03. Diseases & Phenotypes - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 6 |
| Top Score | 26/27 |
| Key Identifiers | MONDO, OMIM, HPO, EFO, MeSH |

## Selected Sources

### 3.1 Disease Ontologies

#### MONDO | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [2/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: MONDO, OMIM, Orphanet, DOID, MeSH, ICD, UMLS
- **Why Include**: THE unified disease ontology, maps across all disease terminologies

#### MeSH | Score: 24/27
- **Tier 1**: Identifiers [2/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [3/3], Versioned [3/3], Size [2/3]
- **License**: Public Domain (US Gov)
- **Key IDs**: MeSH UI, Tree Numbers, UMLS
- **Why Include**: THE standard for PubMed literature indexing, SPARQL endpoint

### 3.2 Phenotype Databases

#### HPO | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [2/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [2/3], Size [2/3]
- **License**: Custom open access
- **Key IDs**: HPO, OMIM, ORPHA, Gene symbols, UMLS
- **Why Include**: THE phenotype ontology, 156K+ disease-phenotype annotations

### 3.3 Disease-Gene Associations

#### Open Targets | Score: 26/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [3/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: Ensembl Gene, EFO, ChEMBL, UniProt, Reactome, PubMed
- **Why Include**: Best-in-class integration of genetics, drugs, pathways (20+ sources)

### 3.4 Cancer Oncology

#### GDC/TCGA | Score: 21/27
- **Tier 1**: Identifiers [2/3], License [2/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [1/3], Versioned [3/3], Size [2/3]
- **License**: NIH GDS Policy (open-access tier)
- **Key IDs**: TCGA barcode, GDC UUID, gene symbols
- **Why Include**: THE cancer genomics reference dataset

### 3.5 Rare Diseases

#### Orphanet | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [2/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: ORPHA, OMIM, ICD-10, UMLS, gene symbols
- **Why Include**: THE rare disease portal, comprehensive nomenclature and gene associations

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| EFO | 21/27 | Use MONDO instead (more comprehensive) |
| OMIM | 19/27 | Restrictive license, use MONDO/HPO mappings |
| DisGeNET | 20/27 | CC BY-NC-SA (non-commercial) |
| ICD | 17/27 | CC BY-ND, access via MONDO mappings |
| DECIPHER | 17/27 | Limited open data |
| ImmunoBase | 18/27 | Outdated (2019), use GWAS Catalog |
| PGC | 20/27 | Specialized GWAS, use GWAS Catalog |
| Allen Brain Atlas | 18/27 | Expression resource, not disease ontology |

## Integration Points

```
MONDO (disease hub)
   ├── HPO (phenotypes)
   ├── Orphanet (rare diseases)
   └── MeSH (literature)
         ↓
   Open Targets (gene-disease)
         ↓
   GDC/TCGA (cancer specifics)
```
