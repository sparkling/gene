# 07. Proteins & Molecular Biology - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 3 |
| Top Score | 26/27 |
| Key Identifiers | UniProt AC (central hub), PDB ID, Ensembl |

## Selected Sources

### 7.1 Protein Sequences & Annotations

#### UniProt | Score: 26/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [3/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: UniProt AC, Gene ID, Ensembl, RefSeq, PDB, GO, InterPro, Pfam, KEGG, Reactome, OMIM, STRING, ChEMBL
- **Why Include**: THE gold standard for protein annotation, best-in-class ID mapping to 200+ databases

#### Key Features
- Swiss-Prot: 570K reviewed entries (ideal MVP starting point)
- TrEMBL: 250M unreviewed entries
- Full SPARQL endpoint at sparql.uniprot.org
- Monthly updates

### 7.2 Protein Structures

#### PDB | Score: 26/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [3/3], Size [3/3]
- **License**: CC0 (Public Domain)
- **Key IDs**: PDB ID, UniProt AC, Gene names, Ligand IDs
- **Why Include**: THE authoritative experimental structure source (220K structures)

#### Key Features
- Weekly updates (Wednesdays)
- Multiple formats: mmCIF, PDB, FASTA
- AWS Open Data, rsync mirroring
- Rich ligand cross-references

#### AlphaFold DB | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [2/3], XRefs [2/3]
- **Tier 3**: RDF [1/3], Versioned [3/3], Size [2/3]
- **License**: CC BY 4.0
- **Key IDs**: AlphaFold ID (AF-{UniProt}-F{fragment}), UniProt AC
- **Why Include**: Revolutionary AI predictions (200M+ structures), direct UniProt mapping

#### Key Features
- Human proteome: ~23K proteins (~15GB) - fits MVP scope
- pLDDT confidence scores for quality filtering
- Google Cloud + AWS Open Data

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| RefSeq | 23/27 | Redundant with UniProt for protein-centric MVP |
| SWISS-MODEL | 18/27 | Superseded by AlphaFold, CC BY-SA adds complexity |

## Integration Points

UniProt serves as the **central protein hub** for the entire knowledge graph:

```
                    ┌─── PDB (experimental structures)
                    │
UniProt ───────────┼─── AlphaFold DB (predicted structures)
(Central Hub)       │
                    ├─── Gene IDs (NCBI, Ensembl, HGNC)
                    │
                    ├─── Pathways (Reactome, KEGG, GO)
                    │
                    ├─── Diseases (OMIM, DisGeNET)
                    │
                    └─── Compounds (ChEMBL, DrugBank)
```

## UniProt ID Mapping Service

UniProt's ID mapping is critical for MVP integration:

| From | To | Coverage |
|------|-------|----------|
| UniProt | NCBI Gene | 99%+ |
| UniProt | Ensembl | 99%+ |
| UniProt | PDB | 50%+ (where structures exist) |
| UniProt | Reactome | 80%+ |
| UniProt | ChEMBL | 40%+ (druggable) |

## Data Size Estimates

| Source | MVP Subset | Full Size |
|--------|------------|-----------|
| UniProt (Swiss-Prot) | ~2GB | ~100GB (TrEMBL) |
| PDB | ~100GB (mmCIF) | ~100GB |
| AlphaFold (Human) | ~15GB | ~23TB (full) |

## Recommendation

For MVP, use:
1. **UniProt Swiss-Prot** (reviewed entries only) as central hub
2. **PDB** for experimental structures
3. **AlphaFold** human proteome for predicted structures

This provides comprehensive protein coverage while keeping data manageable.
