---
id: research-interventions-priority
title: "Interventions Data Sources: Final Recommendations"
type: research-priority
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [research, priorities, swarm-synthesis, interventions, pharmaceuticals, tcm, kampo, ayurveda, herbal]
---

# Interventions Data Sources: Final Recommendations

**Document Type:** Synthesis & Recommendations
**Date:** January 2026
**Research Swarm:** 6 parallel research agents
**Total Research Output:** ~200KB across 6 detailed documents
**Parent:** [_index.md](./_index.md)

---

## Executive Summary

After comprehensive analysis by a 6-agent research swarm, here are the definitive recommendations for intervention data integration in the Gene platform's knowledge base.

### Key Decision: Multi-Source Integration Strategy

| Category | Primary Source | Coverage | Access | License |
|----------|---------------|----------|--------|---------|
| **Pharmaceuticals** | PharmGKB + CPIC | 1000+ drugs, 164 w/guidelines | REST API, Download | CC BY-SA 4.0 |
| **TCM** | BATMAN-TCM 2.0 | 54,832 formulas, 2.3M TTIs | REST API, Download | CC BY-NC |
| **Kampo** | KampoDB + TM-MC 2.0 | 298 + 192 JP herbs | Web + Download | CC BY-SA 4.0 |
| **Ayurveda** | IMPPAT 2.0 | 4,010 plants, 17,967 compounds | Export, Download | CC BY 4.0 |
| **Western Herbal** | DSLD + Dr. Duke's | 200K products + phytochemicals | REST API, CSV | CC0 |
| **Natural Products** | COCONUT + LOTUS | 695K + 750K compounds | REST API, SPARQL | CC0 |

---

## 1. Pharmaceutical/Pharmacogenomics Sources

### Recommended: PharmGKB + CPIC + DrugBank

| Source | Records | Focus | Access | Priority |
|--------|---------|-------|--------|----------|
| **PharmGKB** | 20K+ annotations | SNP-drug relationships | REST API | **CRITICAL** |
| **CPIC** | 164 drugs | Dosing guidelines | JSON download | **CRITICAL** |
| **DrugBank** | 15K+ drugs | Drug info, targets | Academic license | HIGH |
| **ChEMBL** | 2.8M compounds | Bioactivity data | REST API | HIGH |
| **DGIdb** | 70K+ interactions | Drug-gene interactions | GraphQL API | MEDIUM |
| **OpenFDA** | FAERS data | Adverse events | REST API | MEDIUM |

**Data Integration Flow:**
```
User SNP → PharmGKB/CPIC → Dosing Guidelines
         → DrugBank → Drug Mechanisms, Interactions
         → ChEMBL → Binding Affinity Data
         → OpenFDA → Adverse Event Warnings
```

### Storage Estimate (Pharmaceuticals)
- PharmGKB core data: ~500 MB
- CPIC guidelines: ~50 MB
- DrugBank (academic): ~1 GB
- ChEMBL (filtered): ~5 GB
- **Total: ~7 GB**

---

## 2. Traditional Chinese Medicine (TCM) Sources

### Recommended: BATMAN-TCM 2.0 + HERB 2.0

| Source | Formulas | Herbs | Compounds | Targets | Access |
|--------|----------|-------|-----------|---------|--------|
| **BATMAN-TCM 2.0** | 54,832 | 8,404 | 39,171 | 2.3M TTIs | REST API |
| **HERB 2.0** | 9 | 7,263 | 49,258 | 12,933 | Web |
| **TCMBank** | - | 9,192 | 61,966 | 15,179 | Download |
| **SymMap** | - | 499 | 19,595 | 4,302 | Download |

**Priority Integration Order:**
1. BATMAN-TCM 2.0 - Best API, largest target predictions
2. HERB 2.0 - Gene expression validation
3. SymMap - Symptom mapping (unique)

**Storage Estimate (TCM):**
- BATMAN-TCM target data: ~2 GB
- HERB expression data: ~1 GB
- **Total: ~3 GB**

---

## 3. Japanese Kampo Sources

### Recommended: KampoDB + TM-MC 2.0 + STORK

| Source | Formulas | Compounds | Targets | License |
|--------|----------|-----------|---------|---------|
| **KampoDB** | 298 | 3,002 | 62,906 | CC BY-SA 4.0 |
| **TM-MC 2.0** | 5,075 | 34,107 | 13,992 | Open |
| **STORK** | 148 (all approved) | - | - | Open reference |
| **KNApSAcK KAMPO** | 336 | - | - | Non-commercial |

**Integration Strategy:**
```
Layer 1: STORK (official 148 approved formulas)
Layer 2: KampoDB (target predictions)
Layer 3: TM-MC 2.0 (JP pharmacopoeia herbs)
Layer 4: KNApSAcK (metabolite-species relationships)
```

**Storage Estimate (Kampo): ~1 GB**

---

## 4. Ayurveda/Indian Traditional Medicine Sources

### Recommended: IMPPAT 2.0 + NPASS + CMAUP

| Source | Plants | Compounds | Targets | License |
|--------|--------|-----------|---------|---------|
| **IMPPAT 2.0** | 4,010 | 17,967 | 27,365 | CC BY 4.0 |
| **OSADHI** | 6,959 | 27,440 | - | Open |
| **GRAYU** | 12,743 | 129,542 | Indirect | Research only |
| **NPACT** | - | 1,574 | 1,980 | Open |
| **CMAUP** | 5,645 | 48,000 | 646 | Free |

**Note:** TKDL (454K+ formulations) is **restricted access** - not recommended.

**Storage Estimate (Ayurveda): ~2 GB**

---

## 5. Western Herbal/Supplements Sources

### Recommended: DSLD + Dr. Duke's + ODS API

| Source | Products/Plants | Compounds | Access | License |
|--------|-----------------|-----------|--------|---------|
| **DSLD** | 200,000+ labels | Ingredients | REST API | CC0 |
| **Dr. Duke's** | Thousands | Many | CSV download | CC0 |
| **ODS API** | 80+ fact sheets | - | REST API | Public domain |
| **Health Canada LNHPD** | All CA NHPs | - | REST API | Open Gov |
| **NatMed Pro** | Comprehensive | - | Enterprise API | Commercial |

**API Access Priority:**
1. DSLD - Best API, complete label data, CC0
2. Health Canada LNHPD - Free REST API, daily updates
3. Dr. Duke's - CC0 bulk download, phytochemicals
4. ODS - Authoritative fact sheets

**Note:** NatMed Pro ($182/year) and Examine.com (API in development) are valuable for efficacy ratings but require licensing.

**Storage Estimate (Western Herbal): ~3 GB**

---

## 6. Cross-System Natural Product Sources

### Recommended: COCONUT + LOTUS + NPAtlas

| Source | Compounds | Access | License | Focus |
|--------|-----------|--------|---------|-------|
| **COCONUT 2.0** | 695,133 | REST API | CC0 | Aggregated NPs |
| **LOTUS** | 750,000+ pairs | SPARQL | CC0 | Structure-organism |
| **NPAtlas** | 36,545 | REST API | CC BY 4.0 | Microbial NPs |
| **NPASS** | 204,023 | Download | Academic | Quantitative activity |
| **SuperNatural 3.0** | 449,058 | Web | Academic | MoA, targets |

**Integration Value:**
- COCONUT/LOTUS provide structure identification
- NPASS provides quantitative bioactivity
- NPAtlas fills microbial gap
- All can cross-reference via InChIKey

**Storage Estimate (Natural Products): ~10 GB**

---

## 7. Target Prediction Pipeline

### Recommended Tools (All Free)

| Tool | Method | Batch | Use Case |
|------|--------|-------|----------|
| **SwissTargetPrediction** | 2D/3D similarity | No | Initial prediction |
| **SEA** | Ligand set similarity | No | Validation |
| **PharmMapper** | 3D pharmacophore | No | Alternative targets |

**Workflow:**
```
Compound Structure → SwissTargetPrediction (2D/3D)
                   → PharmMapper (pharmacophore)
                   → Consensus targets
                   → Validate against STRING/ChEMBL
                   → Network analysis (Cytoscape/STRING)
```

---

## 8. Data Architecture

### Proposed Schema (RuVector-Compatible)

```typescript
// Compound collection
const compoundsCollection = {
  name: 'compounds',
  dimension: 384,
  distanceMetric: 'cosine',
  properties: {
    inchikey: { type: 'string', indexed: true },
    smiles: { type: 'string' },
    name: { type: 'string', indexed: true },
    mol_weight: { type: 'number' },
    sources: { type: 'string[]' },  // ['tcmsp', 'imppat', 'dsld']
    systems: { type: 'string[]' },  // ['tcm', 'ayurveda', 'western']
    drug_likeness: { type: 'number' },
    np_likeness: { type: 'number' },
  }
};

// Targets collection
const targetsCollection = {
  name: 'targets',
  dimension: 384,
  properties: {
    uniprot_id: { type: 'string', indexed: true },
    gene_symbol: { type: 'string', indexed: true },
    name: { type: 'string' },
    target_class: { type: 'string' },
  }
};

// Relationships (Graph)
// (:Compound)-[:TARGETS {activity_type, activity_value, source}]->(:Target)
// (:Compound)-[:FOUND_IN]->(:Plant)
// (:Plant)-[:USED_IN]->(:Formula)
// (:Formula)-[:TREATS]->(:Condition)
// (:SNP)-[:AFFECTS_RESPONSE]->(:Compound)
```

### ID Mapping Strategy

| System | Primary ID | Cross-Reference |
|--------|------------|-----------------|
| Compounds | InChIKey | PubChem CID, ChEMBL ID |
| Targets | UniProt ID | Gene Symbol, Ensembl |
| Plants | NCBI Taxon ID | Common names, TCM names |
| Diseases | MeSH ID | ICD-10, DOID |

---

## 9. Storage & Cost Estimates

### Total Database Size

| Category | Raw Data | With Embeddings | Notes |
|----------|----------|-----------------|-------|
| Pharmaceuticals | 7 GB | 10 GB | PharmGKB + DrugBank + ChEMBL |
| TCM | 3 GB | 5 GB | BATMAN-TCM + HERB |
| Kampo | 1 GB | 1.5 GB | KampoDB + TM-MC + STORK |
| Ayurveda | 2 GB | 3 GB | IMPPAT + OSADHI + GRAYU |
| Western Herbal | 3 GB | 5 GB | DSLD + Dr. Duke's |
| Natural Products | 10 GB | 15 GB | COCONUT + LOTUS + NPASS |
| **Total** | **27 GB** | **40-45 GB** | With RuVector compression |

### Infrastructure Cost (Self-Hosted)

| Configuration | Storage | RAM | Monthly Cost |
|---------------|---------|-----|--------------|
| MVP (subset) | 50 GB | 32 GB | ~$0 (existing server) |
| Full | 100 GB | 64 GB | ~$50-100 |
| Comprehensive | 200 GB | 128 GB | ~$100-200 |

---

## 10. Implementation Timeline

### Phase 1: Pharmaceuticals (Week 1-2)
- [ ] PharmGKB API integration
- [ ] CPIC guidelines download + parse
- [ ] DrugBank academic license application
- [ ] SNP-drug relationship mapping

### Phase 2: TCM + Kampo (Week 3-4)
- [ ] BATMAN-TCM 2.0 API integration
- [ ] HERB 2.0 data extraction
- [ ] KampoDB web scraping or contact maintainers
- [ ] Target predictions linkage

### Phase 3: Ayurveda + Western Herbal (Week 5-6)
- [ ] IMPPAT 2.0 export processing
- [ ] DSLD API integration
- [ ] Dr. Duke's CSV import
- [ ] Cross-system compound deduplication

### Phase 4: Natural Products Integration (Week 7-8)
- [ ] COCONUT REST API integration
- [ ] LOTUS SPARQL queries
- [ ] NPASS bioactivity data
- [ ] Target prediction pipeline

---

## 11. Licensing Summary

### Free for Commercial Use (CC0/Public Domain)

| Source | License | Notes |
|--------|---------|-------|
| DSLD | CC0 | Complete label database |
| Dr. Duke's | CC0 | Phytochemical data |
| COCONUT | CC0 | Largest NP database |
| LOTUS | CC0 | Wikidata-hosted |
| CPIC | CC0 | Dosing guidelines |
| OpenFDA | Public Domain | Adverse events |
| ODS API | Public Domain | Fact sheets |

### Free for Academic/Attribution Required

| Source | License | Commercial Notes |
|--------|---------|------------------|
| PharmGKB | CC BY-SA 4.0 | Free with attribution |
| KampoDB | CC BY-SA 4.0 | Free with attribution |
| NPAtlas | CC BY 4.0 | Free with attribution |
| ChEMBL | Open Access | Free for all |
| IMPPAT 2.0 | CC BY 4.0 | Free with attribution |

### Academic Only / Contact for Commercial

| Source | License | Notes |
|--------|---------|-------|
| DrugBank | CC BY-NC | Academic free, commercial license required |
| BATMAN-TCM 2.0 | CC BY-NC | Contact for commercial |
| TCMBank | Non-commercial | Contact authors |
| NPASS | Academic | Verify terms |

### Commercial Subscription Required

| Source | Cost | Value |
|--------|------|-------|
| NatMed Pro | $182/year | Efficacy ratings, interactions |
| Examine.com | TBD | Evidence grades |
| DNP | Subscription | 328K compounds |

---

## 12. Research Documents Summary

| Document | Size | Key Content |
|----------|------|-------------|
| `interventions-pharma.md` | 34 KB | 23+ databases, PharmGKB/CPIC deep dive |
| `interventions-tcm.md` | 33 KB | 30+ TCM databases, BATMAN-TCM API |
| `interventions-kampo.md` | 30 KB | 25+ databases, KampoDB schema |
| `interventions-ayurveda.md` | 30 KB | 16+ databases, IMPPAT focus |
| `interventions-western-herbal.md` | 30 KB | 25+ databases, DSLD API |
| `interventions-natural-products.md` | 48 KB | Cross-system databases, target prediction |
| **Total Research** | **~205 KB** | Comprehensive intervention coverage |

---

## 13. Final Recommendations

### For a Bootstrapped Genetics Platform:

1. **Start with CC0/public domain sources** - DSLD, CPIC, COCONUT, LOTUS
2. **Apply for academic licenses** - DrugBank, PharmGKB
3. **Prioritize pharmacogenomics** - SNP-drug relationships are core value
4. **Use BATMAN-TCM API** - Best programmatic access for TCM
5. **IMPPAT for Ayurveda** - Most comprehensive, CC BY 4.0
6. **Cross-reference via InChIKey** - Universal compound identifier

### Data Priority (for MVP):

| Priority | Data Type | Source | Why |
|----------|-----------|--------|-----|
| 1 | SNP-Drug relationships | PharmGKB/CPIC | Core platform value |
| 2 | Drug interactions | DrugBank + OpenFDA | Safety critical |
| 3 | TCM compounds + targets | BATMAN-TCM + HERB | Largest traditional medicine system |
| 4 | Supplement labels | DSLD | Western market focus |
| 5 | Natural product structures | COCONUT/LOTUS | Cross-reference all sources |

### Future Expansion Path:

1. **Phase 1 (MVP):** Pharmaceuticals + Western supplements
2. **Phase 2:** TCM + Kampo integration
3. **Phase 3:** Ayurveda + cross-system natural products
4. **Phase 4:** Target prediction pipeline
5. **Phase 5:** NatMed Pro/Examine.com licensing (if revenue supports)

---

## Data Set Size

| Metric | Value |
|--------|-------|
| Pharmaceuticals (raw) | ~7 GB (PharmGKB + DrugBank + ChEMBL) |
| TCM databases | ~3 GB (BATMAN-TCM + HERB) |
| Kampo databases | ~1 GB (KampoDB + TM-MC + STORK) |
| Ayurveda databases | ~2 GB (IMPPAT + OSADHI + GRAYU) |
| Western Herbal | ~3 GB (DSLD + Dr. Duke's) |
| Natural Products | ~10 GB (COCONUT + LOTUS + NPASS) |
| Total raw data | ~27 GB |
| With embeddings | ~40-45 GB (RuVector compression) |
| Last updated | January 2026 |

---

## Download

| Resource | Method | URL |
|----------|--------|-----|
| PharmGKB | REST API | https://api.clinpgx.org/ |
| CPIC Guidelines | Download | https://cpicpgx.org/guidelines/ |
| BATMAN-TCM 2.0 | REST API | https://batman.bi.a.u-tokyo.ac.jp/ |
| IMPPAT 2.0 | Export | https://cb.imppat.org/ |
| DSLD | REST API | https://dsld.nlm.nih.gov/ |
| COCONUT 2.0 | REST API | https://coconut.naturalproducts.net/ |
| LOTUS | SPARQL | https://lotus.naturalproducts.net/ |

**Access Requirements:** Academic or institutional credentials for some sources; most freely available.

## Data Format

| Format | Description |
|--------|-------------|
| JSON | Structured data from APIs |
| XML | Hierarchical data formats |
| CSV/TSV | Tabular data downloads |
| SDF | Chemical structure files |
| UTF-8 | Text encoding standard |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `compound_id` | string | Unique compound identifier | "CHEMBL123456" |
| `source_system` | string | Origin system (TCM, Ayurveda, etc.) | "tcm" |
| `targets` | array | Target proteins/genes | ["CYP2D6", "ABCB1"] |
| `clinical_significance` | string | Therapeutic relevance | "high" |
| `license` | string | Data usage rights | "CC BY-SA 4.0" |

## Sample Data

### Example Compound Record
```json
{
  "compound_id": "BATMAN-TCM-1234",
  "compound_name": "Curcumin",
  "inchikey": "GHASVSINZPUNUT-UHFFFAOYSA-N",
  "smiles": "O1C(=CC(OC)=CC1=O)C(=CC2=CC(OC)=C(OC)C=C2)O",
  "sources": ["BATMAN-TCM", "COCONUT", "LOTUS"],
  "systems": ["tcm", "ayurveda"],
  "targets": ["TNF", "IL6", "CYP3A4"],
  "mol_weight": 368.38,
  "drug_likeness": 0.85
}
```

## License

| Resource | License | Notes |
|----------|---------|-------|
| PharmGKB | CC BY-SA 4.0 | Attribution required |
| CPIC | CC0 | Public domain |
| BATMAN-TCM | CC BY-NC | Non-commercial use |
| IMPPAT | CC BY 4.0 | Attribution required |
| DSLD | CC0 | Public domain |
| COCONUT | CC0 | Public domain |

## Data Set Size

| Metric | Value |
|--------|-------|
| Pharmaceutical compounds | ~100K (PharmGKB + DrugBank) |
| TCM compounds | ~54K (BATMAN-TCM 2.0) |
| Kampo compounds | ~3K (KampoDB) |
| Ayurveda compounds | ~18K (IMPPAT 2.0) |
| Natural products | ~695K (COCONUT 2.0) |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Intervention | Therapeutic action including drugs, supplements, or herbal medicines | Metformin, curcumin |
| Pharmacogenomics | Study of how genes affect drug response | CYP2D6 metabolizer status |
| Target Prediction | Computational prediction of protein targets for compounds | BATMAN-TCM scoring |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| PharmGKB | Pharmacogenomics knowledge base | SNP-drug relationships |
| CPIC | Clinical Pharmacogenetics Implementation Consortium | Dosing guidelines |
| BATMAN-TCM | TCM target prediction database | Compound-target mapping |
| IMPPAT | Indian Medicinal Plants database | Ayurveda compounds |
| KampoDB | Japanese Kampo medicine database | Herbal formulas |
| DSLD | Dietary Supplement Label Database | Supplement products |
| COCONUT | COlleCtion of Open Natural prodUcTs | Natural products |
| LOTUS | Natural Products Online | Natural products |
| TTI | Target-Target Interaction | Protein networks |
| InChIKey | International Chemical Identifier hash | Compound matching |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | Data access |
| CC BY | Creative Commons Attribution | Open license |
| CC BY-NC | Creative Commons Attribution Non-Commercial | Restricted license |
| CC BY-SA | Creative Commons Attribution Share-Alike | Open license |
| CC0 | Creative Commons Zero (Public Domain) | Public domain |
| CPIC | Clinical Pharmacogenetics Implementation Consortium | Guidelines |
| FAERS | FDA Adverse Event Reporting System | Safety data |
| MVP | Minimum Viable Product | Development phase |
| SNP | Single Nucleotide Polymorphism | Genetic variant |
| TCM | Traditional Chinese Medicine | Herbal medicine |
| TTI | Target-Target Interaction | Network data |

---

*This recommendation synthesizes findings from 6 parallel research agents analyzing pharmaceutical, TCM, Kampo, Ayurveda, Western herbal, and natural product databases for intervention data integration.*
