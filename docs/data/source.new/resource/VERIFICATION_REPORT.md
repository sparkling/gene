---
title: "Migration Verification Report"
date: 2026-01-23
status: complete
version: "1.0"
---

# Migration Verification Report

## Executive Summary

A comprehensive verification was performed comparing the source directories (`databases/`, `domains/`, `operations/`) against the destination (`resource/`). The core migration infrastructure is complete, but significant gaps exist in specialty database coverage.

| Source | Coverage | Status |
|--------|----------|--------|
| Symlinks | **100%** | ✅ All valid |
| operations/schemas/ | **97.6%** | ✅ Nearly complete |
| databases/ | **~50%** | ⚠️ Gaps in specialty areas |
| domains/ | **~40%** | ⚠️ View docs, different purpose |

---

## 1. Symlink Verification

**Result: 31/31 Valid (100%)**

All polyhierarchical symlinks are functioning correctly.

| # | Symlink Location | Target | Status |
|---|------------------|--------|--------|
| 1 | 01.genetics.genomics/1.5.expression.regulation/disgenet | 03.diseases.phenotypes/3.3.disease.gene.associations/disgenet | ✓ |
| 2 | 02.compounds.molecules/2.1.natural.products/imppat | 05.traditional.medicine/5.2.south.east.asian.systems/imppat | ✓ |
| 3 | 02.compounds.molecules/2.3.traditional.medicine.compounds/batman.tcm | 05.traditional.medicine/5.1.traditional.chinese.medicine/batman.tcm | ✓ |
| 4 | 02.compounds.molecules/2.3.traditional.medicine.compounds/herb | 05.traditional.medicine/5.1.traditional.chinese.medicine/herb | ✓ |
| 5 | 02.compounds.molecules/2.3.traditional.medicine.compounds/kampodb | 05.traditional.medicine/5.2.south.east.asian.systems/kampodb | ✓ |
| 6 | 02.compounds.molecules/2.4.food.compounds.nutrients/foodb | 06.nutrition.food/6.1.food.composition/foodb | ✓ |
| 7 | 02.compounds.molecules/2.5.drug.metabolism.pharmacokinetics/cpic | 01.genetics.genomics/1.4.pharmacogenomics/cpic | ✓ |
| 8 | 02.compounds.molecules/2.5.drug.metabolism.pharmacokinetics/pharmgkb | 01.genetics.genomics/1.4.pharmacogenomics/pharmgkb | ✓ |
| 9 | 02.compounds.molecules/2.7.compound.target.interactions/chembl | 02.compounds.molecules/2.2.pharmaceuticals/chembl | ✓ |
| 10 | 03.diseases.phenotypes/3.3.disease.gene.associations/clinvar | 01.genetics.genomics/1.1.variant.repositories/clinvar | ✓ |
| 11 | 03.diseases.phenotypes/3.3.disease.gene.associations/gwas.catalog | 01.genetics.genomics/1.5.expression.regulation/gwas.catalog | ✓ |
| 12 | 03.diseases.phenotypes/3.4.cancer.oncology/civic | 01.genetics.genomics/1.6.cancer.genomics/civic | ✓ |
| 13 | 03.diseases.phenotypes/3.4.cancer.oncology/cosmic | 01.genetics.genomics/1.6.cancer.genomics/cosmic | ✓ |
| 14 | 03.diseases.phenotypes/3.4.cancer.oncology/oncokb | 01.genetics.genomics/1.6.cancer.genomics/oncokb | ✓ |
| 15 | 03.diseases.phenotypes/3.5.rare.orphan.diseases/omim | 03.diseases.phenotypes/3.2.phenotype.databases/omim | ✓ |
| 16 | 03.diseases.phenotypes/3.7.mental.health.neurological/gtex | 01.genetics.genomics/1.5.expression.regulation/gtex | ✓ |
| 17 | 04.pathways.networks/4.2.signaling.pathways/kegg | 04.pathways.networks/4.1.metabolic.pathways/kegg | ✓ |
| 18 | 04.pathways.networks/4.2.signaling.pathways/reactome | 04.pathways.networks/4.1.metabolic.pathways/reactome | ✓ |
| 19 | 04.pathways.networks/4.4.drug.target.interactions/dgidb | 02.compounds.molecules/2.7.compound.target.interactions/dgidb | ✓ |
| 20 | 04.pathways.networks/4.4.drug.target.interactions/open.targets | 03.diseases.phenotypes/3.3.disease.gene.associations/open.targets | ✓ |
| 21 | 04.pathways.networks/4.6.regulatory.networks/encode | 01.genetics.genomics/1.5.expression.regulation/encode | ✓ |
| 22 | 05.traditional.medicine/5.3.western.global.herbal/dr.dukes | 02.compounds.molecules/2.1.natural.products/dr.dukes | ✓ |
| 23 | 05.traditional.medicine/5.4.multi.system.integration/wikidata | 08.literature.knowledge/8.2.knowledge.bases/wikidata | ✓ |
| 24 | 06.nutrition.food/6.1.food.composition/usda.fooddata | 02.compounds.molecules/2.4.food.compounds.nutrients/usda.fooddata | ✓ |
| 25 | 06.nutrition.food/6.3.bioactive.food.compounds/phenol.explorer | 02.compounds.molecules/2.4.food.compounds.nutrients/phenol.explorer | ✓ |
| 26 | 06.nutrition.food/6.3.bioactive.food.compounds/phytohub | 02.compounds.molecules/2.4.food.compounds.nutrients/phytohub | ✓ |
| 27 | 07.proteins.molecular.biology/7.3.molecular.interactions/intact | 04.pathways.networks/4.3.protein.protein.interactions/intact | ✓ |
| 28 | 07.proteins.molecular.biology/7.3.molecular.interactions/reactome | 04.pathways.networks/4.1.metabolic.pathways/reactome | ✓ |
| 29 | 07.proteins.molecular.biology/7.3.molecular.interactions/string | 04.pathways.networks/4.3.protein.protein.interactions/string | ✓ |
| 30 | 08.literature.knowledge/8.3.identifier.mapping/uniprot | 07.proteins.molecular.biology/7.1.protein.sequences.annotations/uniprot | ✓ |
| 31 | 08.literature.knowledge/8.4.regulatory.legal/dailymed | 02.compounds.molecules/2.2.pharmaceuticals/dailymed | ✓ |

---

## 2. Operations/Schemas Migration

**Result: 40/41 Schema Files Migrated (97.6%)**

### Successfully Migrated Schemas

| # | Source (operations/schemas/) | Destination (resource/) | Lines |
|---|------------------------------|-------------------------|-------|
| 1 | alphamissense-schema.md | 01.genetics.genomics/1.1.variant.databases/alphamissense/schema.md | 549 |
| 2 | batman-tcm-schema.md | 05.traditional.medicine/5.1.tcm.databases/batman-tcm/schema.md | 711 |
| 3 | cbioportal-schema.md | 03.diseases.phenotypes/3.2.clinical.databases/cbioportal/schema.md | 602 |
| 4 | chebi-schema.md | 02.compounds.molecules/2.6.chemical.ontology.classification/chebi/schema.md | 740 |
| 5 | chembl-schema.md | 02.compounds.molecules/2.2.pharmaceuticals/chembl/schema.md | 939 |
| 6 | clinvar-schema.md | 01.genetics.genomics/1.1.variant.databases/clinvar/schema.md | 877 |
| 7 | coconut-schema.md | 02.compounds.molecules/2.1.natural.products/coconut/schema.md | 651 |
| 8 | dbnsfp-schema.md | 01.genetics.genomics/1.1.variant.databases/dbnsfp/schema.md | 726 |
| 9 | dbsnp-schema.md | 01.genetics.genomics/1.1.variant.databases/dbsnp/schema.md | 648 |
| 10 | dbvar-schema.md | 01.genetics.genomics/1.1.variant.databases/dbvar/schema.md | 659 |
| 11 | dgidb-open-targets-schema.md | 03.diseases.phenotypes/3.4.drug.disease.targets/dgidb/schema.md | 736 |
| 12 | disgenet-schema.md | 03.diseases.phenotypes/3.2.clinical.databases/disgenet/schema.md | 577 |
| 13 | dr-dukes-schema.md | 05.traditional.medicine/5.4.ethnobotany/dr-dukes/schema.md | 606 |
| 14 | dsld-nih.md | 06.nutrition.food/6.2.supplements/dsld/schema.md | 665 |
| 15 | encode-schema.md | 01.genetics.genomics/1.3.functional.genomics/encode/schema.md | 775 |
| 16 | fda-regulatory-schema.md | 08.literature.knowledge/8.4.regulatory.legal/fda.openfda/schema.md | 2403 |
| 17 | foodb.md | 06.nutrition.food/6.1.food.composition/foodb/schema.md | 515 |
| 18 | gene-ontology-schema.md | 04.pathways.networks/4.5.gene.function.ontology/gene.ontology/schema.md | 707 |
| 19 | gmrepo-schema.md | 09.microbiome/9.1.gut.microbiome/gmrepo/schema.md | 556 |
| 20 | gnomad-schema.md | 01.genetics.genomics/1.1.variant.databases/gnomad/schema.md | 781 |
| 21 | gwas-catalog-schema.md | 01.genetics.genomics/1.2.gwas.qtl/gwas-catalog/schema.md | 542 |
| 22 | hmp-schema.md | 09.microbiome/9.1.gut.microbiome/hmp/schema.md | 706 |
| 23 | hpo-schema.md | 03.diseases.phenotypes/3.1.disease.ontologies/hpo/schema.md | 714 |
| 24 | imppat-schema.md | 05.traditional.medicine/5.2.ayurveda.databases/imppat/schema.md | 819 |
| 25 | intact-schema.md | 04.pathways.networks/4.3.protein.protein.interactions/intact/schema.md | 746 |
| 26 | kampodb-schema.md | 05.traditional.medicine/5.3.kampo.databases/kampodb/schema.md | 729 |
| 27 | kgml-schema.md | 04.pathways.networks/4.1.metabolic.pathways/kegg/schema.md | 342 |
| 28 | lotus-schema.md | 02.compounds.molecules/2.1.natural.products/lotus/schema.md | 713 |
| 29 | mondo-schema.md | 03.diseases.phenotypes/3.1.disease.ontologies/mondo/schema.md | 243 |
| 30 | open-food-facts-schema.md | 06.nutrition.food/6.1.food.composition/open-food-facts/schema.md | 267 |
| 31 | orphanet-ordo-schema.md | 03.diseases.phenotypes/3.1.disease.ontologies/orphanet/schema.md | 266 |
| 32 | pharmgkb-schema.md | 03.diseases.phenotypes/3.3.pharmacogenomics/pharmgkb/schema.md | 396 |
| 33 | pubchem-schema.md | 02.compounds.molecules/2.6.chemical.ontology.classification/pubchem/schema.md | 447 |
| 34 | reactome-schema.md | 04.pathways.networks/4.1.metabolic.pathways/reactome/schema.md | 695 |
| 35 | spliceai-schema.md | 01.genetics.genomics/1.1.variant.databases/spliceai/schema.md | 597 |
| 36 | string-schema.md | 04.pathways.networks/4.3.protein.protein.interactions/string/schema.md | 689 |
| 37 | uniprot-idmapping-schema.md | 08.literature.knowledge/8.3.identifier.mapping/uniprot.id.mapping/schema.md | 382 |
| 38 | usda-fooddata-central.md | 06.nutrition.food/6.1.food.composition/usda/schema.md | 355 |
| 39 | wikidata-schema.md | 08.literature.knowledge/8.2.knowledge.bases/wikidata/schema.md | 278 |
| 40 | wikipathways-gpml-schema.md | 04.pathways.networks/4.1.metabolic.pathways/wikipathways/schema.md | 416 |

### Requires Attention

| File | Issue | Resolution |
|------|-------|------------|
| binding-affinity-schema.md (1630 lines) | Combined schema for TTD, BindingDB, GtoPdb | Needs splitting into 3 individual schema.md files |

### Intentionally Not Migrated (Meta Documents)

| File | Purpose | Status |
|------|---------|--------|
| _index.md | Navigation index | Superseded by resource/ structure |
| schemas-index.md | Schema index | Superseded by resource/ structure |
| schemas-navigation.md | Navigation links | Superseded by resource/ structure |
| pathway-formats.md | Multi-format reference | Retained as reference document |
| sample-data.md | Sample data examples | Retained as reference document |
| ruvector-three-worlds-schema.md | Architecture design | Retained as architecture document |
| unified-schema-analysis.md | Schema analysis | Retained as analysis document |

---

## 3. Databases Folder Migration

**Result: ~50% Coverage**

### Source Files Analyzed

| Category | Files | Databases Documented |
|----------|-------|---------------------|
| Genetics | 3 | 30+ |
| Compounds | 4 | 25+ |
| Pathways | 4 | 20+ |
| Traditional Medicine | 6 | 60+ |
| Nutrition | 1 | 7+ |
| Literature | 7 | 15+ |
| **Total** | **26** | **150+** |

### Migration Status by Category

#### Genetics (~50% migrated)

| Migrated | Missing |
|----------|---------|
| dbSNP, ClinVar, gnomAD, dbNSFP, AlphaMissense, CADD, SpliceAI, ENCODE, UK Biobank, TOPMed, 1000 Genomes | All of Us, ALFA R4, GenomeAsia 100K, H3Africa, MaveDB, RegulomeDB, IHEC, MethBank, 4D Nucleome, FANTOM5, SV4GD, HGSVC |

#### Compounds (~48% migrated)

| Migrated | Missing |
|----------|---------|
| ChEMBL, DrugBank, PubChem, DailyMed, SuperCYP, SwissADME | PharmVar, Flockhart Table, BindingDB, GtoPdb, T3DB, PK-DB, NatMed Pro |

#### Traditional Medicine (~25% migrated)

| Migrated | Missing |
|----------|---------|
| BATMAN-TCM, HERB, SymMap, ETCM, TCMBank, TCMSID, IMPPAT, KampoDB, Dr. Duke's, NAPRALERT, EMA Herbal | TCMID, YaTCM, HIT 2.0, TCMSP, CMAUP, TCMM, TCM-Mesh, DCABM-TCM, LTM-TCM, CPMCP, TCMIO, TM-MC, TCMNP, SuperTCM, TCMPG, IGTCM, ANPDB, SANCDB, p-ANAPL, ETM-DB, AfroDb, NuBBEDB, BIOFACQUIM, MAMPDB, SuperNatural, DSLD, ODS API, LNHPD, HerbMed, USP, MSK Herbs, Commission E, ESCOP, WHO Monographs, PhytoHub, MPD3 (35+ missing) |

#### Pathways (~50% migrated)

| Migrated | Missing |
|----------|---------|
| Reactome, KEGG, WikiPathways, STRING, Gene Ontology, DisGeNET, HPO, MONDO, Orphanet, OMIM | MetaCyc, SMPDB, NDEx, PathBank, Rhea, BRENDA, IntEnz, SABIO-RK, BioModels, Human Protein Atlas, MalaCards, Monarch |

---

## 4. Domains Folder Analysis

**Result: ~40% Database Coverage (By Design)**

The domains/ folder contains **health-focused "view" documents** that aggregate databases by clinical application area. These are NOT meant to be 1:1 migrated but provide valuable cross-cutting perspectives.

### Domain Files Overview

| Domain File | Purpose | Databases Listed | In resource/ |
|-------------|---------|-----------------|--------------|
| mental-cognitive.md | Mental health & cognitive | 35+ | ~15 |
| cardio-metabolic.md | Cardiovascular & metabolic | 10 | ~5 |
| cancer-oncology.md | Cancer & oncology | 8 | 8 ✓ |
| autoimmune.md | Autoimmune & hormone | 23 | ~8 |
| rare.md | Rare diseases & biobanks | 9 | 9 ✓ |
| womens-pediatric.md | Women's health & pediatric | 30+ | ~10 |
| microbiome.md | Microbiome | 12 | 10 |
| allergy-pain.md | Allergy, histamine & pain | 26 | ~5 |
| oral-skin-sensory.md | Oral, skin & sensory | 20+ | ~3 |
| sleep-longevity-nutri.md | Sleep, longevity & nutrigenomics | 30+ | ~5 |
| clinical-environmental-mito.md | Environmental & mitochondrial | 10+ | ~2 |
| clinical-biomarkers-labs.md | Biomarkers & lab references | 21 | ~2 |
| community-patient-networks.md | Patient communities | 30+ | 0 |

### Major Gaps Identified

| Category | Missing Databases | Priority |
|----------|-------------------|----------|
| **Longevity/Aging** | HAGR (GenAge, LongevityMap, DrugAge, CellAge, GenDR), Open Genes | HIGH |
| **Circadian/Sleep** | CircaDB, CirGRDB, Sleep Disorders Portal | HIGH |
| **Environmental/Toxicology** | CTD, T3DB, EPA CompTox, HERCULES | HIGH |
| **Mitochondrial** | MITOMAP, HmtDB | HIGH |
| **Biomarkers/Labs** | MarkerDB, LOINC, CALIPER, TheMarker, BioGPS | MEDIUM |
| **Sensory** | RetNet, eyeGENE, DVD, Gene4HL, BitterDB | MEDIUM |
| **Cardiovascular GWAS** | CARDIoGRAMplusC4D, GLGC, DIAGRAM, GIANT, ICBP | MEDIUM |

---

## 5. Recommendations

### Immediate Actions

1. **Split binding-affinity-schema.md** into individual schema.md files for:
   - TTD (Therapeutic Target Database)
   - BindingDB
   - GtoPdb (Guide to Pharmacology)

2. **Create new resource categories** for missing domains:
   - `10.environmental.toxicology/` - CTD, T3DB, EPA CompTox, MITOMAP, HmtDB
   - `11.biomarkers.reference/` - MarkerDB, LOINC, CALIPER
   - `12.longevity.aging/` - HAGR suite, Open Genes, CircaDB

### High Priority Database Additions

| Database | Category | Reason |
|----------|----------|--------|
| CTD | Toxicology | 50M+ toxicogenomics relationships |
| MITOMAP | Mitochondrial | Primary mtDNA resource |
| HAGR suite | Aging | 2,500+ aging genes |
| MarkerDB | Biomarkers | Clinical biomarker reference |
| CircaDB | Circadian | Circadian gene expression |

### Documentation Retention

The following source files should be **retained** as complementary documents:

1. **domains/*.md** - Health-focused view documents
2. **operations/pathway-formats.md** - Technical format reference
3. **operations/unified-schema-analysis.md** - Architecture analysis
4. **databases/literature/*.md** - Pipeline design documents

---

## 6. Summary Statistics

### Current State

| Metric | Count |
|--------|-------|
| Category folders | 9 |
| Total .md files in resource/ | 300 |
| _index.md files | 190 |
| schema.md files | 40 |
| download.md files | 30 |
| xrefs.md files | 30 |
| Valid symlinks | 31 |

### Migration Coverage

| Source | Total Files | Content Coverage |
|--------|-------------|------------------|
| operations/schemas/ | 48 | **97.6%** |
| databases/ | 26 | **~50%** |
| domains/ | 14 | **~40%** (by design) |

### Gaps Summary

| Category | Missing Databases |
|----------|-------------------|
| Traditional Medicine | 35+ |
| Longevity/Aging | 30+ |
| Biomarkers/Labs | 19 |
| Environmental/Toxicology | 10+ |
| Sensory/Specialty | 10+ |
| Population Genetics | 4 |
| **Total** | **~100+** |

---

## Generated: 2026-01-23

This report was generated by a 5-agent verification swarm analyzing:
- Symlink integrity
- operations/schemas/ migration
- databases/ content coverage
- domains/ content coverage
- Cross-source reconciliation
