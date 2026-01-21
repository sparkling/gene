# Microbiome Data Sources

**Document ID:** 43-77-MICROBIOME
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [43-00-INDEX.md](./43-00-INDEX.md)

---

## TL;DR

Microbiome databases span gut, oral, and skin communities with 130K+ samples across major repositories. Priority sources include GMrepo (REST API, 119K samples, disease markers), HMP (gold standard reference, 48TB), and gutMGene (2.5M microbe-gene associations). Commercial use requires attention to NC licenses on GMrepo, mBodyMap, and gutMGene.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary gut microbiome source | GMrepo v3 | REST API ready, disease markers, 119K samples | Jan 2026 |
| Reference baseline | HMP via AWS | Gold standard healthy reference, CC BY 4.0 | Jan 2026 |
| Microbe-gene linking | gutMGene v2.0 | 2.5M associations, metabolite data | Jan 2026 |
| Drug-microbiome interactions | MASI | 13K+ interactions, comprehensive coverage | Jan 2026 |
| Multi-body-site coverage | mBodyMap | 22 body sites, 63K runs | Jan 2026 |
| Probiotic recommendations | Probio-ichnos | 13K entries, CC BY 4.0 license | Jan 2026 |

---

## Database Catalog

### 1. Human Microbiome Project (HMP)

| Field | Value |
|-------|-------|
| **URL** | https://hmpdacc.org/ |
| **Data Portal** | https://portal.hmpdacc.org/ |
| **AWS Open Data** | https://registry.opendata.aws/human-microbiome-project/ |
| **Content** | Multi-omic human microbiome reference (16S, WGS, transcriptome, proteome, metabolome) |
| **Records** | 11,000+ samples from 300 healthy adults across 18 body sites |
| **License** | CC BY 4.0 (some controlled data requires dbGaP authorization) |
| **API** | HMP Portal API (GraphQL/REST), GA4GH DRS API, HMP16SData (Bioconductor R package), AWS S3 |
| **Update Frequency** | Archive maintained (active 2007-2016) |
| **Priority** | Tier 1 (Reference baseline) |
| **Storage Estimate** | ~50 TB (combined HMP and iHMP phases) |

**Gene Platform Integration Notes:**
- Excellent for baseline healthy microbiome reference
- Multi-omic data enables correlation with host genetics
- AWS hosting facilitates cloud-based pipeline integration

---

### 2. MetaHIT (Metagenomics of the Human Intestinal Tract)

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/ena/browser/view/PRJNA32811 |
| **MAHMI Database** | https://academic.oup.com/database/article/doi/10.1093/database/baw157/2884895 |
| **Content** | European gut metagenomics, gene catalog, enterotype discovery |
| **Records** | 124 individuals (Denmark/Spain), 3.3M unique ORFs |
| **License** | Open Access (EBI data sharing policies) |
| **API** | ENA Browser Tools (Python), wget/FTP, Globus, MAHMI REST API |
| **Update Frequency** | Static (foundational dataset) |
| **Priority** | Tier 2 (Gene catalog for annotation) |
| **Storage Estimate** | ~577 GB raw sequence data |

**Gene Platform Integration Notes:**
- Foundation dataset for gut microbiome research
- Gene catalog useful for metagenomic annotation pipelines
- Enterotype classification relevant for personalization features

---

### 3. GMrepo (Gut Microbiome Repository)

| Field | Value |
|-------|-------|
| **URL** | https://gmrepo.humangut.info |
| **Documentation** | https://evolgeniusteam.github.io/gmrepodocumentation/ |
| **GitHub** | https://github.com/evolgeniusteam/GMrepoProgrammableAccess |
| **Content** | Curated gut metagenomes with disease-marker associations |
| **Records** | 890 projects, 118,965 samples (87K 16S + 32K WGS), 302 diseases, 1,299 marker taxa |
| **License** | CC BY-NC 3.0 (commercial use requires permission) |
| **API** | REST API with R/Perl/Python examples |
| **Update Frequency** | Active curation |
| **Priority** | Tier 1 (Primary gut microbiome) |
| **Storage Estimate** | ~50 GB processed data (raw sequences not directly downloadable) |

**API Endpoints:**
```
GET /api/getRunsByProjectID
GET /api/getRunsByPhenotype
GET /api/getSpeciesAbundance
GET /api/getMarkerTaxa
```

**Gene Platform Integration Notes:**
- Excellent for disease-marker microbiome associations
- Pre-computed relative abundances reduce processing needs
- Cross-dataset comparison features valuable for meta-analyses
- REST API enables real-time integration

---

### 4. mBodyMap

| Field | Value |
|-------|-------|
| **URL** | https://mbodymap.microbiome.cloud |
| **GitHub** | https://github.com/whchenlab/mBodymap |
| **Content** | Multi-body-site microbiome associations with health/disease |
| **Records** | 63,148 runs (14K metagenomes + 49K amplicons), 22 body sites, 56 diseases, 6,247 species |
| **License** | CC BY-NC 3.0 (academic use) |
| **API** | REST API with R/Python examples |
| **Update Frequency** | Active curation |
| **Priority** | Tier 1 (Multi-site coverage) |
| **Storage Estimate** | ~30 GB processed abundance data |

**API Endpoints:**
```
GET /api/getAllBodySites
GET /api/getSpeciesOfOneBodySite
GET /api/getRelativeAbundanceOfOneBodySite
GET /api/getDiseaseAssociations
```

**Gene Platform Integration Notes:**
- Multi-body-site coverage unique among databases
- Disease-microbe associations pre-curated
- Useful for systemic microbiome analysis beyond gut
- MeSH standardization enables integration with medical ontologies

---

### 5. gutMGene (Gut Microbiota-Gene Database)

| Field | Value |
|-------|-------|
| **URL** | http://bio-computing.hrbmu.edu.cn/gutmgene |
| **Legacy** | http://bio-annotation.cn/gutmgene |
| **Content** | Microbe-metabolite-gene relationships in humans and mice |
| **Records** | v1.0: 3,680 relationships; v2.0: 2.48M associations (including metabolic reconstructions) |
| **License** | CC BY-NC 4.0 (commercial requires permission) |
| **API** | Web interface with autocomplete, bulk download |
| **Update Frequency** | Version updates (v2.0 current) |
| **Priority** | Tier 1 (Microbe-gene linking) |
| **Storage Estimate** | ~5 GB structured data |

**v1.0 Content:**
| Organism | Relationships | Microbes | Metabolites | Genes |
|----------|---------------|----------|-------------|-------|
| Human | 1,331 | 332 | 207 | 223 |
| Mouse | 2,349 | 209 | 149 | 544 |

**Gene Platform Integration Notes:**
- Critical for microbiome-gene interaction mapping
- Metabolite data enables pathway analysis
- Human + mouse data supports translational research
- Causal vs correlational classification aids interpretation

---

### 6. MASI (Microbiota-Active Substance Interactions)

| Field | Value |
|-------|-------|
| **URL** | https://www.aiddlab.com/MASI/ |
| **Content** | Comprehensive microbiota-substance interactions |
| **Records** | 1,051 pharmaceuticals, 103 dietary, 119 herbal, 46 probiotic, 142 environmental substances; 806 microbiota species; 13K+ interactions |
| **License** | Not explicitly stated (academic use appears unrestricted) |
| **API** | Web browse, CSV/Excel/PDF export, bulk download |
| **Update Frequency** | Active curation |
| **Priority** | Tier 1 (Drug-microbiome) |
| **Storage Estimate** | ~500 MB structured data |

**Interaction Categories:**
| Type | Count |
|------|-------|
| Bacteria-Pharmaceutical | 11,215 |
| Bacteria-Herbal | 914 |
| Bacteria-Dietary | 309 |
| Bacteria-Environmental | 753 |

**Cross-links:** NCBI Taxonomy, ChEMBL, DrugBank, TTD, PubChem, NPASS

**Gene Platform Integration Notes:**
- Comprehensive drug interaction coverage
- Multi-category substance coverage (pharma, dietary, herbal, environmental)
- Disease associations (56 diseases, 784 microbiota-disease links)

---

### 7. PharmacoMicrobiomics Portal

| Field | Value |
|-------|-------|
| **URL** | http://pharmacomicrobiomics.com/ |
| **Source Code** | http://sourceforge.net/projects/pharmacomicro |
| **Content** | Drug-microbiome interactions affecting metabolism, toxicity, efficacy |
| **Records** | 60+ drugs curated from 100+ research articles |
| **License** | CC BY 2.5 |
| **API** | Web search, cross-links to PubMed/PubChem/Comparative Toxicogenomics |
| **Update Frequency** | Community-maintained |
| **Priority** | Tier 2 (Supplementary drug data) |
| **Storage Estimate** | ~100 MB |

---

### 8. MDAD (Microbe-Drug Association Database)

| Field | Value |
|-------|-------|
| **URL** | https://www.frontiersin.org/articles/10.3389/fcimb.2018.00424/full |
| **Content** | Comprehensive microbe-drug associations |
| **Records** | 5,055 entries, 1,373 drugs, 173 microbes, 2,470 unique associations |
| **License** | CC BY |
| **API** | Web interface, user data submission |
| **Update Frequency** | Community-maintained |
| **Priority** | Tier 2 (Supplementary drug data) |
| **Storage Estimate** | ~200 MB |

**Data Sources:** DrugBank, TTD, SuperTarget, MATADOR, TDR targets, PDTD, ChEMBL, aBiofilm, Integrity

---

### 9. PROBIO Database

| Field | Value |
|-------|-------|
| **URL** | http://bidd.group/probio/homepage.htm |
| **Content** | Marketed and research probiotic strains |
| **Records** | 997 strains (448 marketed, 167 clinical trial, 382 research) |
| **License** | Academic use |
| **API** | Web interface, statistics download |
| **Update Frequency** | Periodic updates |
| **Priority** | Tier 2 (Probiotic data) |
| **Storage Estimate** | ~50 MB |

---

### 10. Probio-ichnos

| Field | Value |
|-------|-------|
| **URL** | https://www.mdpi.com/2076-2607/12/10/1955 |
| **Content** | Comprehensive probiotic strain characteristics |
| **Records** | 12,993 entries, 11,202 distinct strains, 470 species, 143 genera |
| **License** | CC BY 4.0 |
| **API** | Web interface, strain-specific queries |
| **Update Frequency** | Published 2024 |
| **Priority** | Tier 1 (Probiotic recommendations) |
| **Storage Estimate** | ~100 MB |

**Probiotic Characteristics Covered:**
- Resistance to acid and bile
- Attachment to host epithelia
- Antimicrobial activity
- Immunomodulatory activity
- Antiproliferative activity
- Antioxidant activity

---

### 11. ProbioMinServer

| Field | Value |
|-------|-------|
| **URL** | https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad153/7321112 |
| **Content** | Probiotic genome analysis (ARGs, virulence factors, CAZy) |
| **Records** | 25 common probiotic species (11 genera) |
| **License** | Academic use |
| **API** | Web interface, tab-delimited export, WGS upload for analysis |
| **Update Frequency** | Published 2023 |
| **Priority** | Tier 3 (Specialized analysis) |
| **Storage Estimate** | ~20 MB |

---

### 12. IPDB (Integrated Probiotic DataBase)

| Field | Value |
|-------|-------|
| **URL** | http://probiogenomics.unipr.it/cmu/ |
| **Download** | http://probiogenomics.unipr.it/files/Probiotic_Bifidobacteria_DataBase |
| **Content** | Bifidobacterium genomic compendium |
| **Records** | 34 publicly available Bifidobacterium strains |
| **License** | Academic research use |
| **API** | Direct download |
| **Update Frequency** | Static |
| **Priority** | Tier 3 (Specialized) |
| **Storage Estimate** | ~10 MB |

---

## Comparative Summary

| Database | Focus | Samples/Entries | API | License | Best For |
|----------|-------|-----------------|-----|---------|----------|
| **HMP** | Multi-omic human microbiome | 11,000+ samples, >48TB | REST, GraphQL, GA4GH | CC BY 4.0 | Reference baseline |
| **MetaHIT** | Gut metagenomics | 124 subjects, 3.3M genes | ENA tools, FTP | Open Access | Gene catalog |
| **GMrepo** | Gut microbiome curation | 118,965 samples | REST API | CC BY-NC 3.0 | Disease markers |
| **mBodyMap** | Multi-body-site microbiome | 63,148 runs | REST API | CC BY-NC 3.0 | Body-wide analysis |
| **gutMGene** | Microbe-gene interactions | 2.5M associations | Web/Download | CC BY-NC 4.0 | Pathway analysis |
| **MASI** | Drug-microbiome | 13,000+ interactions | Download | Not stated | Pharmacomicrobiomics |
| **MDAD** | Drug-microbe | 5,055 entries | Web | CC BY | Drug associations |
| **Probio-ichnos** | Probiotic strains | 12,993 entries | Web | CC BY 4.0 | Strain selection |

---

## Integration Recommendations

### Priority 1: Core Databases (MVP)

| Database | Rationale | Integration Method |
|----------|-----------|-------------------|
| **GMrepo v3** | REST API ready, disease markers, academic-friendly | REST API integration |
| **gutMGene v2.0** | Essential for gene-microbiome linking | Bulk download + periodic sync |
| **MASI** | Comprehensive drug interaction coverage | Bulk download |

### Priority 2: Reference Data (Post-MVP)

| Database | Rationale | Integration Method |
|----------|-----------|-------------------|
| **HMP** | Gold standard healthy reference | AWS S3 cloud integration |
| **mBodyMap** | Multi-site coverage with API | REST API integration |

### Priority 3: Specialized (Future)

| Database | Rationale | Integration Method |
|----------|-----------|-------------------|
| **Probio-ichnos** | Probiotic recommendation features | Bulk download |
| **MetaHIT** | Gene catalog for annotation | ENA/FTP download |
| **MDAD/PharmacoMicrobiomics** | Supplementary drug data | Bulk download |

---

## Technical Considerations

### Commercial Use Licensing

| License Type | Databases | Action Required |
|--------------|-----------|-----------------|
| **CC BY 4.0** (Unrestricted) | HMP, Probio-ichnos | None |
| **CC BY** (Unrestricted) | MDAD | None |
| **CC BY 2.5** | PharmacoMicrobiomics | None |
| **CC BY-NC 3.0** (Academic Only) | GMrepo, mBodyMap | Commercial license negotiation |
| **CC BY-NC 4.0** (Academic Only) | gutMGene | Commercial license negotiation |
| **Not Stated** | MASI | Legal review recommended |
| **Academic Use** | MetaHIT, PROBIO, ProbioMinServer, IPDB | Commercial license negotiation |

### API Availability

| Quality | Databases |
|---------|-----------|
| **Best-documented REST APIs** | GMrepo, mBodyMap |
| **Standard APIs** | HMP (GraphQL/REST), GA4GH DRS |
| **Download/FTP only** | MetaHIT, gutMGene, MASI, MDAD |
| **Web interface only** | Probiotic databases |

### Storage Requirements

| Category | Estimate |
|----------|----------|
| HMP (full) | ~50 TB |
| MetaHIT (raw) | ~577 GB |
| Processed databases (GMrepo, mBodyMap, gutMGene, MASI, MDAD, probiotics) | ~1 TB total |
| **Recommended MVP** | ~2 TB (processed data only, no raw sequences) |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [43-DATA-SOURCES](../43-DATA-SOURCES.md) | Parent overview |
| [43-00-INDEX](./43-00-INDEX.md) | Navigation index |
| [43-11-GENETICS-PRIMARY](./43-11-GENETICS-PRIMARY.md) | Host genetics integration |
| [43-51-PHARMACEUTICALS](./43-51-PHARMACEUTICALS.md) | Drug interaction cross-reference |
| [44-ARCHITECTURE](../44-ARCHITECTURE.md) | Informs database design |
| [45-DATA-MODEL](../45-DATA-MODEL.md) | Informs entity structure |

---

## Open Questions

- [ ] HMP controlled access data - is dbGaP authorization needed for platform use?
- [ ] GMrepo/mBodyMap commercial license - negotiation timeline and cost?
- [ ] Oral and skin microbiome - sufficient coverage in mBodyMap or need additional sources?
- [ ] Microbiome data update pipeline - how to sync with source databases?
- [ ] Raw sequence storage - needed or processed data sufficient?

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial document migrated from research.old/data-sources-microbiome.md |
