# 02. Compounds & Molecules - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 10 |
| Top Score | 26/27 |
| Key Identifiers | PubChem CID, InChIKey, ChEBI ID, SMILES |

## Selected Sources

### 2.1 Natural Products

#### LOTUS | Score: 26/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [2/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [3/3], Size [3/3]
- **License**: CC0 (Public Domain)
- **Key IDs**: InChIKey, Wikidata QID, NCBI Taxon ID
- **Why Include**: Structure-organism pairs (750K+), full RDF/Wikidata integration

#### COCONUT | Score: 23/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [1/3], Versioned [2/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: InChIKey, SMILES, CAS, PubChem CID
- **Why Include**: Largest open NP database (~450K compounds), aggregates 52+ sources

### 2.2 Pharmaceuticals

#### ChEMBL | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [3/3], Size [2/3]
- **License**: CC BY-SA 3.0
- **Key IDs**: ChEMBL ID, InChIKey, SMILES, UniProt
- **Why Include**: THE canonical bioactivity database (2.4M compounds, 20M activities)

#### RxNorm | Score: 21/27
- **Tier 1**: Identifiers [2/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [3/3], Size [2/3]
- **License**: UMLS (free, commercial OK)
- **Key IDs**: RXCUI, NDC, ATC
- **Why Include**: THE canonical drug nomenclature standard (300K+ concepts)

### 2.4 Food Compounds

#### USDA FoodData Central | Score: 21/27
- **Tier 1**: Identifiers [2/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [2/3]
- **License**: Public Domain (US Gov)
- **Key IDs**: FDC ID, NDB number, GTIN/UPC
- **Why Include**: Authoritative food composition (380K+ foods), continuous updates

### 2.6 Chemical Ontology

#### ChEBI | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [3/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: ChEBI ID, InChIKey, SMILES
- **Why Include**: THE chemical ontology (150K+ entities), OBO Foundry member, full RDF

#### PubChem | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [2/3], Size [1/3]
- **License**: Public Domain (US Gov)
- **Key IDs**: CID, SID, InChIKey, SMILES
- **Why Include**: THE largest chemical repository (115M+ compounds), central ID hub

### 2.7 Compound-Target Interactions

#### DGIdb | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: Gene symbol, Entrez ID, ChEMBL ID
- **Why Include**: Aggregates 40+ sources (90K+ interactions), small size (~100MB)

#### BindingDB | Score: 23/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [2/3], Size [2/3]
- **License**: CC BY 3.0
- **Key IDs**: MonomerID, InChIKey, UniProt
- **Why Include**: Best quantitative binding affinities (2.9M measurements)

#### GtoPdb | Score: 22/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [2/3], Size [3/3]
- **License**: CC BY-SA 4.0
- **Key IDs**: GtoPdb ID, HGNC, UniProt
- **Why Include**: IUPHAR official database, expert-curated pharmacology

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| DrugBank | 16/27 | CC BY-NC 4.0 (non-commercial) |
| NPASS | 15/27 | CC BY-NC 4.0 (non-commercial) |
| TTD | 17/27 | CC BY-NC 4.0 (non-commercial) |
| SuperCYP | 13/27 | Non-commercial, outdated (2013) |
| Phenol-Explorer | 15/27 | Non-commercial, outdated (2015) |
| PhytoHub | 14/27 | Restricted license, outdated |
| DailyMed | 20/27 | Too large (60GB), limited chemical IDs |
| NPAtlas | 21/27 | Covered by COCONUT |
| Dr. Duke's | 18/27 | Archived (2016) |

## Integration Points

```
PubChem (CID) ←→ ChEBI (ontology)
       ↓
   ChEMBL (bioactivity)
       ↓
  BindingDB (affinities) ←→ DGIdb (aggregated)
       ↓
   GtoPdb (pharmacology)
       ↓
   RxNorm (clinical naming)
```
