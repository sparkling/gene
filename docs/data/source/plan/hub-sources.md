# Hub Sources - Critical Identity Bridges

These sources serve as central integration anchors, connecting multiple domains via shared identifiers.

## Primary Hubs

### UniProt - Protein Hub
**Score**: 26/27 | **License**: CC BY 4.0

The central protein hub connecting genes to pathways to diseases to structures.

```
Genes (Ensembl, NCBI Gene, HGNC)
                ↓
            UniProt ←→ Structures (PDB, AlphaFold)
                ↓
Pathways (Reactome, KEGG, GO)
                ↓
Diseases (OMIM, DisGeNET)
                ↓
Compounds (ChEMBL, DrugBank)
```

**ID Mappings**: 200+ databases
**Key Service**: UniProt ID Mapping (idmapping.dat)

---

### PubChem - Chemical Hub
**Score**: 24/27 | **License**: Public Domain

The central chemical hub connecting compounds to targets to literature.

```
Chemical Ontology (ChEBI)
            ↓
        PubChem ←→ Bioactivity (ChEMBL, BindingDB)
            ↓
Natural Products (LOTUS, COCONUT)
            ↓
Literature (PubMed)
            ↓
Regulatory (FDA, RxNorm)
```

**ID Mappings**: InChIKey, SMILES, CAS, ChEMBL, DrugBank
**Key Feature**: Largest chemical repository (115M+ compounds)

---

### Wikidata - Universal Hub
**Score**: 26/27 | **License**: CC0

The universal knowledge graph connecting virtually all biological entities.

```
Genes (P351, P594, P354)
            ↓
        Wikidata ←→ Proteins (P352)
            ↓
Compounds (P662, P592, P715)
            ↓
Diseases (P492, P494)
            ↓
Literature (P698, DOI)
```

**Key Properties**:
| Property | ID Type | Domain |
|----------|---------|--------|
| P351 | NCBI Gene | Genes |
| P352 | UniProt | Proteins |
| P354 | HGNC | Gene symbols |
| P492 | OMIM | Diseases |
| P594 | Ensembl | Genes |
| P662 | PubChem | Compounds |
| P592 | ChEMBL | Bioactivity |
| P686 | GO | Functions |
| P698 | PubMed | Literature |
| P715 | DrugBank | Drugs |

**Key Feature**: Native RDF, SPARQL endpoint, continuous community updates

---

### MONDO - Disease Hub
**Score**: 25/27 | **License**: CC BY 4.0

The unified disease ontology connecting all disease terminologies.

```
OMIM (Mendelian)
            ↓
        MONDO ←→ Orphanet (Rare)
            ↓
DOID (Disease Ontology)
            ↓
ICD-10/11 (Clinical)
            ↓
MeSH (Literature)
```

**Cross-References**: 129,914+ mappings
**Key Feature**: Resolves inconsistencies between disease vocabularies

---

### Gene Ontology - Function Hub
**Score**: 27/27 (Perfect) | **License**: CC BY 4.0

The universal standard for functional annotation.

```
Molecular Function
            ↓
    Gene Ontology ←→ Biological Process
            ↓
Cellular Component
            ↓
All protein databases (UniProt, RefSeq, Ensembl)
```

**Annotations**: 45K terms, 700K human annotations
**Key Feature**: OWL/RDF format, SPARQL endpoint, evidence codes

---

## Hub Integration Strategy

### Step 1: Load Hub Sources First
```
1. UniProt (protein identities)
2. PubChem/ChEBI (chemical identities)
3. MONDO (disease identities)
4. Gene Ontology (function identities)
5. Wikidata (universal bridging)
```

### Step 2: Link Domain Sources to Hubs
```
Genetics → UniProt (via Gene ID, Ensembl)
Compounds → PubChem (via InChIKey, CID)
Diseases → MONDO (via OMIM, Orphanet)
Pathways → GO (via GO terms)
All → Wikidata (via QID properties)
```

### Step 3: Cross-Domain Queries
With hubs in place, queries like these become possible:

```sparql
# Find drugs targeting genes associated with diabetes
SELECT ?drug ?gene ?disease
WHERE {
  ?disease mondo:MONDO_0005015 .  # Diabetes
  ?gene rdfs:seeAlso ?disease .
  ?drug chembl:hasTarget ?gene .
}
```

## Hub Completeness by Domain

| Domain | Primary Hub | Secondary Hub | Coverage |
|--------|-------------|---------------|----------|
| Genetics | UniProt | Ensembl | 99%+ |
| Compounds | PubChem | ChEBI | 99%+ |
| Diseases | MONDO | MeSH | 95%+ |
| Pathways | GO | Reactome | 95%+ |
| Literature | PubMed | Wikidata | 99%+ |
| Proteins | UniProt | PDB | 99%+ |
