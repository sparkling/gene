# ADR-004: SPARQL Federation Strategy

## Status
Accepted

## Context

The Gene platform needs access to large biomedical datasets:
- **UniProt**: ~500M+ triples (all species), ~20K human proteins
- **dbSNP/ClinVar**: ~1.2M clinically relevant variants
- **Reactome/KEGG**: ~2K pathways
- **Gene Ontology**: Ontology with ~50K terms
- **Disease Ontology**: ~15K disease terms
- **PubChem/ChEMBL**: Millions of compounds (optional)

SPARQL federation allows querying remote endpoints, but external federation has reliability concerns.

## Decision

**Local import preferred**, with external SPARQL federation only for data sources too large to import.

### Strategy (In Order of Preference)

1. **Local Import** (Preferred)
   - Import all data we need into local QLever instance
   - Most reliable, fastest queries
   - Our dataset (~15-35M triples) easily fits

2. **Multiple QLever Instances** (Fallback 1)
   - If data exceeds single instance capacity
   - Internal federation between our own QLever instances
   - Still under our control

3. **External SPARQL Federation** (Fallback 2)
   - Only for data too large to import (e.g., full UniProt 500M+ triples)
   - Accept reliability tradeoffs
   - Cache aggressively

### Data Import Plan

| Source | Triples | Strategy | Priority |
|--------|---------|----------|----------|
| Human proteins (UniProt subset) | ~500K | Local import | P0 |
| dbSNP/ClinVar variants | ~5M | Local import | P0 |
| Reactome/KEGG pathways | ~500K | Local import | P0 |
| Gene Ontology | ~2M | Local import | P0 |
| Disease Ontology | ~200K | Local import | P0 |
| Herb databases (TCMSP, etc.) | ~1M | Local import | P1 |
| Full UniProt (all species) | ~500M | External federation | P2 |
| PubChem/ChEMBL | ~100M+ | External federation | P2 |

**Total local import: ~10M triples** - well within QLever's capacity.

## Consequences

### Positive
- Reliable, fast queries for core data
- No dependency on external endpoint availability
- Full control over data freshness (quarterly re-import)
- Simpler debugging and testing

### Negative
- Need to download and convert data to RDF
- Storage requirements (~5-10GB)
- Need periodic re-import for freshness

### Neutral
- External federation available when needed
- RDF conversion tooling required (rdflib, Apache Jena)

## Options Considered

### Option 1: Full External Federation
- **Pros**: Always fresh data, no storage
- **Cons**: Unreliable (external endpoints go down), slow, rate limits

### Option 2: Full Local Import
- **Pros**: Most reliable, fastest
- **Cons**: Storage for 500M+ triples if including all UniProt

### Option 3: Hybrid (Chosen)
- **Pros**: Reliable for core data, federation for edge cases
- **Cons**: More complex strategy

## RDF Conversion

Use **rdflib** (Python) or **Apache Jena** for converting biomedical data to N-Triples:

```python
from rdflib import Graph, URIRef, Literal, Namespace

# Example: Convert gene to RDF
GENE = Namespace("http://gene.example.org/")
g = Graph()
g.add((GENE.MTHFR, GENE.symbol, Literal("MTHFR")))
g.add((GENE.MTHFR, GENE.name, Literal("Methylenetetrahydrofolate reductase")))
g.serialize("genes.nt", format="nt")
```

## Related Decisions
- [ADR-001](./adr-001-three-database-architecture.md): Three-Database Architecture
- [ADR-002](./adr-002-knowledge-graph-engine.md): Knowledge Graph Engine Selection

## References
- [SPARQL 1.1 Federation](https://www.w3.org/TR/sparql11-federated-query/)
- [UniProt SPARQL Endpoint](https://sparql.uniprot.org/)
- [rdflib Documentation](https://rdflib.readthedocs.io/)
- [Apache Jena](https://jena.apache.org/)
