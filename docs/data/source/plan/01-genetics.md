# 01. Genetics & Genomics - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 9 |
| Top Score | 26/27 |
| Key Identifiers | rsID, Gene ID, Ensembl, OMIM, UniProt |

## Selected Sources

### 1.1 Variant Repositories

#### ClinVar | Score: 26/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [1/3], Versioned [3/3], Size [3/3]
- **License**: CC0/Public Domain
- **Key IDs**: VCV, RCV, rsID, OMIM, MedGen CUI, Gene ID, UniProt, HPO
- **Why Include**: Clinical variant gold standard, excellent cross-references, weekly updates

#### dbSNP | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [1/3], Versioned [3/3], Size [1/3]
- **License**: CC0/Public Domain
- **Key IDs**: rsID (universal variant identifier)
- **Why Include**: Provides the universal rsID system essential for variant linking

### 1.2 Functional Prediction

#### dbNSFP | Score: 22/27
- **Tier 1**: Identifiers [3/3], License [2/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [1/3]
- **License**: Academic free, commercial contact required
- **Key IDs**: Position (chr:pos), Gene symbols, Ensembl IDs
- **Why Include**: Single file aggregates 36+ prediction scores (CADD, REVEL, SpliceAI, etc.)

#### AlphaMissense | Score: 24/27
- **Tier 1**: Identifiers [2/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [2/3], XRefs [2/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: Genomic position, UniProt, Ensembl transcript
- **Why Include**: State-of-the-art missense predictor (AUROC 0.94), small footprint (643MB)

### 1.3 Population Genetics

#### gnomAD | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [0/3]
- **License**: ODbL 1.0
- **Key IDs**: Variant ID, rsID, Ensembl gene, HGNC
- **Why Include**: Essential population frequencies (807K individuals), constraint metrics

### 1.4 Pharmacogenomics

#### PharmGKB | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [3/3]
- **License**: CC BY-SA 4.0
- **Key IDs**: PharmGKB accessions, rsID, Gene symbol, Drug name
- **Why Include**: THE pharmacogenomics knowledge base, includes CPIC guidelines

### 1.5 Expression & Regulation

#### GWAS Catalog | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [3/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: GCST study ID, rsID, EFO trait ID, PubMed ID
- **Why Include**: Canonical variant-trait associations, has RDF/knowledge graph

#### GTEx | Score: 22/27
- **Tier 1**: Identifiers [3/3], License [2/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [1/3]
- **License**: Open for aggregate data; dbGaP for individual
- **Key IDs**: Ensembl gene, rsID, tissue name
- **Why Include**: Essential tissue expression and eQTL data

### 1.6 Cancer Genomics

#### CIViC | Score: 26/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [2/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [3/3]
- **License**: CC0 (Public Domain)
- **Key IDs**: CIViC IDs, Entrez Gene, DOID, DrugBank
- **Why Include**: Community-curated cancer interpretation, public domain, small size (~35MB)

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| dbVar | 21/27 | SV focus - defer to Phase 2 |
| 1000 Genomes | 23/27 | Superseded by gnomAD |
| COSMIC | 18/27 | Commercial license required |
| CADD | 19/27 | Commercial license; covered by dbNSFP |
| SpliceAI | 17/27 | Commercial license; covered by dbNSFP |
| TOPMed | 16/27 | dbGaP required, no bulk download |
| UK Biobank | 14/27 | Application required |

## Integration Points

```
ClinVar ←→ dbSNP (rsID)
    ↓
gnomAD (frequencies)
    ↓
GWAS Catalog (trait associations)
    ↓
PharmGKB (drug response)
    ↓
GTEx (expression context)
```
