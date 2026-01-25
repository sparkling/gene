# 04. Pathways & Networks - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 8 |
| Top Score | **27/27** (Gene Ontology - perfect score) |
| Key Identifiers | GO, UniProt, Reactome ID, KEGG (via Reactome) |

## Selected Sources

### 4.1 Metabolic Pathways

#### Reactome | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [3/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: Reactome Stable ID (R-HSA-), UniProt, Ensembl, ChEBI
- **Why Include**: Premier open pathway database (2,712 human pathways), BioPAX/SBML/Neo4j

#### WikiPathways | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [2/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [2/3], Size [2/3]
- **License**: CC0 (Public Domain)
- **Key IDs**: WikiPathways ID (WP####), Ensembl, Entrez Gene, UniProt
- **Why Include**: CC0 license, full SPARQL endpoint, complements Reactome

### 4.3 Protein-Protein Interactions

#### STRING | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [1/3], Versioned [3/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: STRING ID (taxid.ENSP), Ensembl Protein, UniProt
- **Why Include**: Gold standard for functional associations, scored evidence

#### IntAct | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [3/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: IntAct ID (EBI-), IMEx ID, UniProt, PSI-MI terms
- **Why Include**: Highest quality experimental PPI evidence, IMEx standards

### 4.4 Drug-Target Interactions

#### STITCH | Score: 22/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [2/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: STITCH Chemical ID, STRING Protein ID, PubChem CID
- **Why Include**: Essential chemical-protein networks (STRING sister database)

### 4.5 Gene Function Ontology

#### Gene Ontology | Score: 27/27 (PERFECT)
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [3/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: GO Term (GO:0000000), UniProt, Ensembl, Evidence Codes
- **Why Include**: THE universal standard for gene function annotation

#### MSigDB | Score: 21/27
- **Tier 1**: Identifiers [3/3], License [2/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [0/3], Versioned [3/3], Size [2/3]
- **License**: CC BY 4.0 (with KEGG exceptions)
- **Key IDs**: Gene Set Name, Systematic Name (M#####), HGNC Symbol
- **Why Include**: THE standard for gene set enrichment analysis, Hallmark sets

### 4.6 Regulatory Networks

#### JASPAR | Score: 25/27
- **Tier 1**: Identifiers [2/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [1/3], Versioned [3/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: JASPAR ID (MA#####.#), UniProt, TF Name
- **Why Include**: Gold standard for TF binding motifs, compact (~3MB)

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| KEGG | 15/27 | **Commercial license required** |
| BioGRID | 24/27 | Redundant with IntAct/STRING |
| Pathway Commons | 21/27 | Aggregator - use primary sources |
| Roadmap Epigenomics | 22/27 | Phase 2 - chromatin states only |

## Integration Points

```
Gene Ontology (function hub)
       ↓
Reactome ←→ WikiPathways (pathways)
       ↓
STRING ←→ IntAct (interactions)
       ↓
STITCH (chemical-protein)
       ↓
JASPAR (regulatory)
```

## Note on KEGG

KEGG was **excluded** despite being canonical because:
- Commercial license required for bulk data
- FTP requires subscription
- Reactome + WikiPathways provide open alternatives with good coverage
