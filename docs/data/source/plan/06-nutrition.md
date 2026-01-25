# 06. Nutrition & Food - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 4 |
| Top Score | 22/27 |
| Key Identifiers | HMDB ID, PubChem CID, ChEBI ID, FDC ID |

## Selected Sources

### 6.1 Food Composition

#### FooDB | Score: 20/27
- **Tier 1**: Identifiers [3/3], License [1/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [2/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [2/3], Size [3/3]
- **License**: CC BY-NC 4.0 (contact for commercial)
- **Key IDs**: PubChem CID, HMDB ID, ChEBI ID, KEGG ID, DrugBank ID
- **Why Include**: Best-in-class food compound database (70K compounds, 1K foods, 800K food-compound links)

### 6.2 Dietary Supplements

#### DSLD | Score: 22/27
- **Tier 1**: Identifiers [2/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [0/3], Versioned [2/3], Size [1/3]
- **License**: Public Domain (US Government)
- **Key IDs**: DSLD ID, UPC Barcode, NHP ID
- **Why Include**: THE canonical US supplement database (NIH), 150K+ products, public domain

### 6.4 Metabolomics

#### HMDB | Score: 21/27
- **Tier 1**: Identifiers [3/3], License [1/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [2/3], Size [3/3]
- **License**: CC BY-NC 4.0 (contact for commercial)
- **Key IDs**: HMDB ID, KEGG ID, PubChem CID, ChEBI ID, InChIKey
- **Why Include**: Gold standard human metabolome (220K metabolites, 650K MS/MS spectra)

#### Exposome-Explorer | Score: 19/27
- **Tier 1**: Identifiers [2/3], License [3/3], Bulk [2/3]
- **Tier 2**: Canonical [2/3], Active [2/3], XRefs [2/3]
- **Tier 3**: RDF [0/3], Versioned [2/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: PubChem CID, HMDB ID, CAS Number
- **Why Include**: Unique dietary biomarker-exposure data (IARC/WHO), complements HMDB

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| Open Food Facts | 18/27 | No chemical identifiers (product-focused) |
| ConsumerLab | 5/27 | Proprietary, subscription required |
| Natural Medicines | 4/27 | Proprietary, subscription required |
| eBASIS | 11/27 | EuroFIR agreement required |
| Phenol-Explorer | 15/27 | CC BY-NC, outdated (2015) |
| PhytoHub | 14/27 | Restricted license, outdated |

## Integration Points

```
FooDB (food compounds)
   ↓ HMDB ID
HMDB (metabolomics)
   ↓ PubChem CID
Exposome-Explorer (biomarkers)
   ↓
DSLD (supplement ingredients → need external ID mapping)
```

## Note on Licenses

FooDB and HMDB both use CC BY-NC 4.0 (non-commercial). For commercial MVP:
- Contact database maintainers for commercial license
- Or use PubChem/ChEBI mappings to get compound data from open sources
- DSLD and Exposome-Explorer are fully open (Public Domain and CC BY 4.0)

## Data Size Estimates

| Source | Size | Format |
|--------|------|--------|
| FooDB | ~1GB | CSV, XML, JSON, SDF |
| HMDB | ~5GB | XML, SDF, TSV |
| DSLD | ~500MB | JSON, CSV |
| Exposome-Explorer | ~50MB | TSV |
