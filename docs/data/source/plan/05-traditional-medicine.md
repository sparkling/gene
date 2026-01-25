# 05. Traditional Medicine - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 4 |
| Top Score | 21/27 |
| Key Identifiers | PubChem CID, UniProt, NCBI Taxonomy |

## Selected Sources

### 5.1 Traditional Chinese Medicine

#### BATMAN-TCM 2.0 | Score: 21/27
- **Tier 1**: Identifiers [3/3], License [2/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [1/3]
- **License**: CC BY-NC 4.0 (contact for commercial)
- **Key IDs**: PubChem CID, DrugBank, TTD, UniProt, KEGG, ChEMBL
- **Why Include**: Only TCM DB with REST API, largest coverage (54K formulas, 39K compounds), high prediction AUC (0.97)

### 5.2 South & East Asian Systems

#### IMPPAT 2.0 | Score: 20/27
- **Tier 1**: Identifiers [3/3], License [2/3], Bulk [2/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [1/3]
- **License**: CC BY 4.0
- **Key IDs**: PubChem CID (via InChIKey/UniChem), ChEMBL, UniProt, HGNC
- **Why Include**: Largest Indian medicinal plant DB (4K plants, 18K compounds), CC BY 4.0, comprehensive ADMET

#### KampoDB | Score: 18/27
- **Tier 1**: Identifiers [3/3], License [1/3], Bulk [2/3]
- **Tier 2**: Canonical [3/3], Active [2/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [1/3], Size [3/3]
- **License**: CC BY-SA 4.0
- **Key IDs**: PubChem CID, NCBI Gene ID, KEGG
- **Why Include**: Only Kampo DB with REST API, 3M+ docking simulations, Japanese traditional medicine

### 5.4 Multi-System Integration

#### HIT 2.0 | Score: 17/27
- **Tier 1**: Identifiers [3/3], License [1/3], Bulk [2/3]
- **Tier 2**: Canonical [3/3], Active [2/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [2/3], Size [1/3]
- **License**: Academic use free
- **Key IDs**: PubChem CID, UniProt, Gene Symbol, PubMed ID
- **Why Include**: 100K+ experimentally validated interactions (not predictions), essential for validation

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| NAPRALERT | 10/27 | **Subscription required** - not open |
| HERB | 18/27 | Overlaps with BATMAN-TCM |
| ETCM | 12/27 | Limited bulk access, academic-only |
| TCMBank | 14/27 | Academic-only license |
| TCMSID | 14/27 | Overlaps with IMPPAT |
| SymMap | 14/27 | Academic-only license |
| NPACT | 11/27 | Narrow focus, limited access |
| TM-MC | 11/27 | No API, variable coverage |
| EMA Herbal | 15/27 | PDF format only |

## Integration Strategy

### Primary ID Bridge: PubChem CID
```
BATMAN-TCM → PubChem CID → other compound DBs
IMPPAT → InChIKey → UniChem → PubChem CID
KampoDB → PubChem CID (native)
HIT 2.0 → PubChem CID (native)
```

### Target ID Bridge: UniProt/Gene Symbol
All 4 sources provide UniProt or Gene Symbol mappings for linking to protein/gene databases.

### Validation Flow
```
BATMAN-TCM (predicted interactions)
        ↓
    HIT 2.0 (validation - which predictions have experimental support?)
        ↓
IMPPAT/KampoDB (extend to other traditional systems)
```

## Coverage by System

| System | Primary Source | Compounds | Targets |
|--------|---------------|-----------|---------|
| TCM | BATMAN-TCM 2.0 | 39K | 8K |
| Ayurveda | IMPPAT 2.0 | 18K | 4K |
| Kampo | KampoDB | 2K | 1K |
| Validation | HIT 2.0 | 100K+ interactions | - |
