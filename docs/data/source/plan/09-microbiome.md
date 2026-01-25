# 09. Microbiome - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 3 |
| Top Score | 24/27 |
| Key Identifiers | NCBI Taxonomy ID, SRA accession, UniProt |

## Selected Sources

### 9.1 Gut Microbiome / 9.2 Body Site Microbiomes

#### HMP (Human Microbiome Project) | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [1/3], Versioned [2/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: NCBI Taxonomy IDs, SRA accessions (SRR/ERR), BioProject (PRJNA), body site codes
- **Why Include**: THE canonical human microbiome reference, NIH-funded, multi-omics (16S, WGS, proteomics, metabolomics)

#### Key Features
- 3,000+ reference genomes across body sites
- 18 body site controlled vocabulary
- AWS S3 (no auth), HTTP, Aspera bulk downloads
- Active updates via iHMP studies

### 9.2 Body Site Microbiomes

#### HOMD (Human Oral Microbiome Database) | Score: 20/27
- **Tier 1**: Identifiers [2/3], License [2/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [3/3]
- **License**: Academic use open, commercial requires permission
- **Key IDs**: HOMT IDs (HMT-###), NCBI Taxonomy IDs, ATCC strain IDs
- **Why Include**: THE authoritative oral microbiome taxonomy (HOMT system), 770+ oral taxa

#### Key Features
- 1,800+ genomes
- Quarterly 16S database updates
- FTP and direct HTTP downloads

### 9.3 Microbe-Host Interactions

#### VMH (Virtual Metabolic Human) | Score: 23/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [1/3], Versioned [2/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: VMH Reaction/Metabolite IDs, NCBI Taxonomy IDs, KEGG Reaction IDs, UniProt IDs
- **Why Include**: THE canonical microbiome metabolic modeling resource, 7,200+ AGORA2 microbial models

#### Key Features
- SBML format (semantic interoperability)
- REST API available
- Links to KEGG, UniProt, PubChem

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| GMrepo | 19/27 | CC BY-NC 3.0 restricts commercial |
| MetaHIT | 18/27 | Legacy resource, no updates |
| gutMGene | 14/27 | Limited bulk access, academic-only |
| mBodyMap | 13/27 | Not canonical, limited bulk access |
| gutMDisorder | 14/27 | No bulk download |
| MASI | 8/27 | Very limited access |

## Integration Points

```
HMP (composition)
   ↓ NCBI Taxonomy ID
VMH (metabolic modeling)
   ↓ UniProt, KEGG
HOMD (oral specialization)
```

## Recommended MVP Strategy

1. **HMP** provides canonical microbiome composition data
2. **VMH** provides metabolic modeling framework for microbiome-host interactions
3. **HOMD** adds oral-specific detail if in scope

All three have permissive licenses (CC BY 4.0 for HMP/VMH) and excellent bulk download capabilities.

## Data Size Estimates

| Source | Size | Notes |
|--------|------|-------|
| HMP | ~500GB | Full multi-omics; subset for MVP |
| VMH | ~1GB | SBML models |
| HOMD | ~2GB | 16S + genomes |

## Linking to Other Categories

Microbiome data links to rest of knowledge graph via:

| Link Type | Source ID | Target Category |
|-----------|-----------|-----------------|
| Taxonomy | NCBI Taxonomy | Literature (PubMed) |
| Metabolites | KEGG, VMH | Compounds (ChEBI, PubChem) |
| Proteins | UniProt | Proteins |
| Genes | NCBI Gene | Genetics |
| Diseases | MeSH, MONDO | Diseases |
