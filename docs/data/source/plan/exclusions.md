# Excluded Sources - Rationale

Sources analyzed but not selected for MVP, organized by exclusion reason.

## License Restrictions

### Commercial License Required

| Source | Category | Score | License | Notes |
|--------|----------|-------|---------|-------|
| **KEGG** | Pathways | 15/27 | Commercial | FTP requires subscription |
| **COSMIC** | Genetics | 18/27 | Commercial | Academic free, commercial paid |
| **CADD** | Genetics | 19/27 | Commercial | Included in dbNSFP anyway |
| **SpliceAI** | Genetics | 17/27 | GPLv3/Commercial | Included in dbNSFP anyway |

### Non-Commercial Only (CC BY-NC)

| Source | Category | Score | License | Notes |
|--------|----------|-------|---------|-------|
| **DrugBank** | Compounds | 16/27 | CC BY-NC 4.0 | Contact for commercial |
| **DisGeNET** | Diseases | 20/27 | CC BY-NC-SA 4.0 | Large but restricted |
| **NPASS** | Compounds | 15/27 | CC BY-NC 4.0 | Natural products bioactivity |
| **TTD** | Compounds | 17/27 | CC BY-NC 4.0 | Drug targets |
| **SuperCYP** | Compounds | 13/27 | CC BY-NC-SA 3.0 | Also outdated |
| **Phenol-Explorer** | Nutrition | 15/27 | CC BY-NC 3.0 | Also outdated |
| **GMrepo** | Microbiome | 19/27 | CC BY-NC 3.0 | Good data though |

### Subscription Required

| Source | Category | Score | License | Notes |
|--------|----------|-------|---------|-------|
| **NAPRALERT** | Trad Med | 10/27 | Subscription | World's largest NP literature DB |
| **ConsumerLab** | Nutrition | 5/27 | Proprietary | $49.95/year |
| **Natural Medicines** | Nutrition | 4/27 | Proprietary | Enterprise only |

### Restrictive Custom License

| Source | Category | Score | License | Notes |
|--------|----------|-------|---------|-------|
| **OMIM** | Diseases | 19/27 | Custom | Registration + commercial license |
| **UK Biobank** | Genetics | 14/27 | Custom | Application required (months) |
| **TOPMed** | Genetics | 16/27 | dbGaP | Individual data controlled |

---

## Access Barriers

### No Bulk Download

| Source | Category | Score | Issue |
|--------|----------|-------|-------|
| **SwissADME** | Compounds | 12/27 | Web-only tool, not database |
| **MASI** | Microbiome | 8/27 | Contact maintainers only |
| **gutMDisorder** | Microbiome | 14/27 | Web interface only |

### Registration/Approval Required

| Source | Category | Score | Issue |
|--------|----------|-------|-------|
| **UK Biobank** | Genetics | 14/27 | Multi-month application |
| **TOPMed** | Genetics | 16/27 | dbGaP application |
| **DECIPHER** | Diseases | 17/27 | Patient data agreement |

---

## Redundancy / Overlap

### Covered by Other Sources

| Source | Category | Score | Covered By |
|--------|----------|-------|------------|
| **RefSeq** | Proteins | 23/27 | UniProt (for proteins) |
| **SWISS-MODEL** | Proteins | 18/27 | AlphaFold DB |
| **1000 Genomes** | Genetics | 23/27 | gnomAD (for frequencies) |
| **NPAtlas** | Compounds | 21/27 | COCONUT |
| **HERB** | Trad Med | 18/27 | BATMAN-TCM |
| **Europe PMC** | Literature | 22/27 | PubMed |
| **EFO** | Diseases | 21/27 | MONDO |
| **BioGRID** | Pathways | 24/27 | IntAct + STRING |
| **Pathway Commons** | Pathways | 21/27 | Primary sources |

---

## Outdated / Archived

| Source | Category | Score | Last Update |
|--------|----------|-------|-------------|
| **Dr. Duke's** | Compounds | 18/27 | 2016 (archived) |
| **Phenol-Explorer** | Nutrition | 15/27 | 2015 |
| **PhytoHub** | Nutrition | 14/27 | 2015 |
| **SuperCYP** | Compounds | 13/27 | 2013 |
| **ImmunoBase** | Diseases | 18/27 | 2019 |
| **MetaHIT** | Microbiome | 18/27 | 2015 (completed) |

---

## Too Large for MVP

| Source | Category | Score | Size | Notes |
|--------|----------|-------|------|-------|
| **DailyMed** | Compounds | 20/27 | ~60GB | Limited chemical IDs |
| **PubMed Central** | Literature | 20/27 | ~1TB | Mixed licensing |

---

## Specialized / Phase 2

### Defer to Later Phases

| Source | Category | Score | Reason |
|--------|----------|-------|--------|
| **dbVar** | Genetics | 21/27 | Structural variants (after SNV) |
| **Roadmap Epigenomics** | Pathways | 22/27 | Chromatin states only |
| **PGC** | Diseases | 20/27 | Psychiatric GWAS (use GWAS Catalog) |
| **Allen Brain Atlas** | Diseases | 18/27 | Brain expression (specialized) |
| **SynGO** | Diseases | 21/27 | Synaptic GO extension |
| **IPD-IMGT/HLA** | Diseases | 20/27 | HLA alleles (specialized) |

---

## Summary Statistics

| Exclusion Reason | Count |
|------------------|-------|
| Commercial license required | 4 |
| Non-commercial only (CC BY-NC) | 7 |
| Subscription required | 3 |
| Restrictive custom license | 3 |
| No bulk download | 3 |
| Registration/approval required | 3 |
| Covered by other sources | 9 |
| Outdated/archived | 6 |
| Too large for MVP | 2 |
| Specialized/Phase 2 | 6 |
| **Total Excluded** | **~46** |

---

## Reconsider for Phase 2

These high-scoring sources should be reconsidered after MVP:

| Source | Score | Barrier | Mitigation |
|--------|-------|---------|------------|
| DrugBank | 16/27 | CC BY-NC | Contact for commercial license |
| DisGeNET | 20/27 | CC BY-NC-SA | Academic use, or Open Targets covers similar |
| KEGG | 15/27 | Commercial | Budget for license |
| RefSeq | 23/27 | Redundant | Add for transcript-protein mapping |
| BioGRID | 24/27 | Redundant | MIT license - easy to add |
