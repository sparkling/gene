# Integration Map - How Sources Connect

Visual map of how the 55 MVP data sources link together via shared identifiers.

## Master Integration Diagram

```
                                    LITERATURE
                                        │
                    ┌───────────────────┼───────────────────┐
                    │                   │                   │
                PubMed ◄────────► Wikidata ◄────────► OpenAlex
                (PMID)              (QID)              (OpenAlex ID)
                    │                   │                   │
                    └───────────────────┼───────────────────┘
                                        │
        ┌───────────────────────────────┼───────────────────────────────┐
        │                               │                               │
    GENETICS                        PROTEINS                       COMPOUNDS
        │                               │                               │
   ┌────┴────┐                    ┌─────┴─────┐                   ┌─────┴─────┐
   │         │                    │           │                   │           │
ClinVar   gnomAD              UniProt ◄──────────────────────► PubChem    ChEBI
(VCV)     (variant)           (AC)     ←──────────────────────  (CID)      (ID)
   │         │                    │                               │           │
   └────┬────┘                    │                               └─────┬─────┘
        │                         │                                     │
      rsID ◄──────────────────────┤                                     │
   (dbSNP)                        │                                     │
        │                         │                                     │
        ▼                         ▼                                     ▼
   ┌─────────┐              ┌─────────┐                           ┌─────────┐
   │  GWAS   │              │   PDB   │                           │ ChEMBL  │
   │ Catalog │              │ (PDB ID)│                           │(ChEMBL) │
   └─────────┘              └─────────┘                           └─────────┘
        │                         │                                     │
        ▼                         ▼                                     ▼
   ┌─────────┐              ┌─────────┐                           ┌─────────┐
   │  EFO    │              │AlphaFold│                           │BindingDB│
   │ (trait) │              │  (AF-)  │                           │  (Ki)   │
   └─────────┘              └─────────┘                           └─────────┘
        │                                                               │
        ▼                                                               ▼
   ┌─────────────────────────────────────────────────────────────────────────┐
   │                              DISEASES                                    │
   │                                                                          │
   │    MONDO ◄────────► HPO ◄────────► Orphanet ◄────────► Open Targets    │
   │   (MONDO)          (HP)           (ORPHA)              (target-disease) │
   │      │               │                │                       │          │
   │      └───────────────┴────────────────┴───────────────────────┘          │
   │                                  │                                       │
   │                              MeSH (UI)                                   │
   └─────────────────────────────────────────────────────────────────────────┘
                                      │
   ┌──────────────────────────────────┼──────────────────────────────────────┐
   │                                  │                                      │
   │                              PATHWAYS                                   │
   │                                                                         │
   │   Gene Ontology ◄────────► Reactome ◄────────► STRING ◄────────► IntAct│
   │      (GO)                 (R-HSA-)            (taxid.ENSP)      (EBI-)  │
   │         │                     │                    │                │   │
   │         └─────────────────────┴────────────────────┴────────────────┘   │
   │                                  │                                      │
   │                            WikiPathways                                 │
   │                              (WP####)                                   │
   └─────────────────────────────────────────────────────────────────────────┘
```

## Identifier Flow by Domain

### Genetics → Proteins
```
rsID (dbSNP/ClinVar)
    → Gene ID (NCBI Gene)
    → UniProt AC
    → Protein sequence/structure
```

### Compounds → Targets
```
PubChem CID / InChIKey
    → ChEMBL ID (bioactivity)
    → UniProt AC (target)
    → Gene ID
```

### Diseases → Genes
```
MONDO ID
    → OMIM (via xref)
    → Gene associations
    → HPO phenotypes
```

### Pathways → Everything
```
GO Term
    → UniProt (proteins)
    → Gene ID (genes)
    → Reactome (pathways)
    → MeSH (literature)
```

## Cross-Domain ID Bridges

### UniProt as Central Bridge
| From | Via UniProt To | Coverage |
|------|----------------|----------|
| Gene ID | PDB structures | 50%+ |
| Ensembl | Reactome pathways | 80%+ |
| HGNC | ChEMBL targets | 40%+ |
| RefSeq | GO annotations | 95%+ |

### Wikidata as Universal Bridge
| Property | Maps | Coverage |
|----------|------|----------|
| P352 | UniProt | 200K+ proteins |
| P351 | NCBI Gene | 50K+ genes |
| P662 | PubChem | 100M+ compounds |
| P492 | OMIM | 10K+ diseases |
| P698 | PubMed | 30M+ articles |

### PubChem as Chemical Bridge
| Source | ID Type | Linkable |
|--------|---------|----------|
| ChEBI | ChEBI ID | Yes |
| ChEMBL | ChEMBL ID | Yes |
| DrugBank | InChIKey | Yes |
| HMDB | HMDB ID | Yes |
| KEGG | KEGG ID | Yes |

## Integration Priority Order

### Phase 1: Load Hubs
```
1. UniProt (protein hub)
2. PubChem + ChEBI (chemical hub)
3. MONDO (disease hub)
4. Gene Ontology (function hub)
5. Wikidata (universal bridge)
```

### Phase 2: Load Domain Sources
```
6. ClinVar, dbSNP, gnomAD (genetics)
7. ChEMBL, BindingDB (bioactivity)
8. HPO, Orphanet, Open Targets (diseases)
9. Reactome, STRING, IntAct (pathways)
10. PubMed, OpenAlex (literature)
```

### Phase 3: Load Specialized
```
11. PDB, AlphaFold (structures)
12. Traditional medicine DBs
13. Nutrition/food DBs
14. Microbiome DBs
```

## Key Integration Patterns

### Pattern 1: Gene-centric
```sparql
Gene → Variants → Diseases → Drugs
(Ensembl) → (ClinVar) → (MONDO) → (ChEMBL)
```

### Pattern 2: Drug-centric
```sparql
Drug → Targets → Pathways → Side Effects
(PubChem) → (UniProt) → (Reactome) → (FDA)
```

### Pattern 3: Disease-centric
```sparql
Disease → Genes → Variants → Literature
(MONDO) → (Open Targets) → (gnomAD) → (PubMed)
```

### Pattern 4: Natural Product-centric
```sparql
Plant → Compound → Target → Disease
(NCBI Taxon) → (LOTUS) → (UniProt) → (MONDO)
```

## Missing Links (Gaps to Address)

| Gap | Workaround |
|-----|------------|
| Food compound → Disease | Via HMDB metabolites |
| Microbiome → Drug | Via VMH metabolic models |
| Traditional medicine → Clinical trials | Manual curation needed |
| Epigenetics → Variants | Roadmap in Phase 2 |
