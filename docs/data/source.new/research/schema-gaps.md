---
id: research-schema-gaps
title: Data Sources Needing Schema Documentation
type: research
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [research, schema, documentation, gaps, integration]
---

**Parent:** [Research](./_index.md)

# Data Sources Needing Schema Documentation

**Created:** January 2026
**Updated:** January 2026 (pruned from 200+ to 95 databases)
**Purpose:** Track data sources that need detailed schemas, download code, and entity relationships

---

## Summary

After pruning analysis, **95 databases** remain from the original 200+. Pruned databases (105 total) are documented in `DATA-SOURCES-PRUNING-ANALYSIS.md` with rationale.

**Pruning removed:**
- 18 commercial/paid license databases
- 15 controlled access databases (lengthy approval)
- 12 legal/scraping risk sources
- 22 out-of-scope databases
- 14 redundant/superseded databases
- 6 defunct databases
- 18 low-value/high-effort sources

**Remaining databases need:**
1. **Detailed field schemas** (field names, data types, JSON examples)
2. **Download code examples** (curl/wget/Python scripts)
3. **Entity relationship documentation**

---

## Tier 1 - Critical for MVP (16 databases)

These must be documented first for minimum viable product.

| Database | Category | Schema Needed | License |
|----------|----------|---------------|---------|
| **dbNSFP v4.9** | Genetics | 35 score columns, field definitions | Academic free |
| **AlphaMissense** | Genetics | TSV fields, score interpretation | CC BY 4.0 |
| **gnomAD-SV v4.1** | Genetics | SV VCF fields, Hail Table schema | ODC-ODbL |
| **ClinVar** | Genetics | VCF fields, submission schema | Public domain |
| **COCONUT 2.0** | Natural Products | REST API schema, SDF fields | CC0 |
| **BATMAN-TCM 2.0** | TCM | Compound-target schema | CC BY-NC |
| **IMPPAT 2.0** | Ayurveda | Plant-compound schema | CC BY 4.0 |
| **GMrepo v3** | Microbiome | REST API endpoints/response | Free |
| **OncoKB** | Cancer | API response schema, evidence levels | Free academic |
| **CIViC** | Cancer | GraphQL schema, evidence types | CC0 |
| **CPIC** | Pharmacogenomics | API response, guideline format | CC BY-SA |
| **PharmVar** | Pharmacogenomics | Star allele definition schema | Free |
| **PGC** | Mental Health | GWAS summary statistics format | Public summaries |
| **Orphanet** | Rare Disease | Disease-gene schema (ORDO) | CC BY 4.0 |
| **HPO** | Phenotypes | Phenotype ontology schema | Open |
| **OpenAlex** | Literature | API response schema | CC0 |

---

## Tier 2 - High Value (35 databases)

Document after Tier 1 is complete.

### Genetics & Functional Annotation

| Database | Schema Needed | License |
|----------|---------------|---------|
| **SNPedia** | MediaWiki API schema | CC BY-NC-SA |
| **CADD v1.7** | Score fields, C-score vs PHRED | Free academic |
| **SpliceAI** | DS_AG/AL/DG/DL fields | Free |
| **MaveDB** | Variant effect measurement schema | CC0 |
| **RegulomeDB v2** | Scoring system fields | Free |
| **ENCODE 4 (cCRE)** | cCRE format, experiment metadata | CC BY 4.0 |
| **ALFA (NCBI)** | VCF format, E-utilities response | Public domain |

### Traditional Medicine

| Database | Schema Needed | License |
|----------|---------------|---------|
| **SymMap 2.0** | Symptom-disease mapping schema | Free |
| **ETCM v2.0** | Formula composition fields | Free |
| **HERB 2.0** | Transcriptomic data schema | Free |
| **LOTUS** | Wikidata SPARQL schema | CC0 |
| **ANPDB** | API response schema, compound fields | Free |
| **SANCDB** | REST API schema | Free |
| **NuBBEDB** | SDF fields, ADMET schema | Free |

### Microbiome

| Database | Schema Needed | License |
|----------|---------------|---------|
| **HMP Portal** | GA4GH DRS schema, GraphQL | Public |
| **gutMGene v2.0** | Gene-microbe relationship schema | Free |
| **mBodyMap** | REST API schema | Free |

### Drugs & Metabolism

| Database | Schema Needed | License |
|----------|---------------|---------|
| **PharmGKB** | Already documented - verify current | CC BY-SA 4.0 |
| **DrugBank** | Already documented - verify current | CC BY-NC 4.0 |
| **KEGG DRUG** | REST API response schema | Attribution |
| **PK-DB** | REST API schema | Free |
| **SuperCYP** | CYP interaction schema | Free |
| **HMDB** | Metabolite XML schema | Free academic |
| **FooDB** | Food compound schema - nutrients, vitamins, bioactives for nutrigenomics | Open |
| **T3DB** | Toxin record schema | Free |

### Pathways & Interactions

| Database | Schema Needed | License |
|----------|---------------|---------|
| **Reactome** | Already documented - verify current | CC0 |
| **STRING** | API response schema | CC BY 4.0 |

### Brain & Expression

| Database | Schema Needed | License |
|----------|---------------|---------|
| **Allen Brain Atlas** | Expression data schema | Free |
| **GTEx** | eQTL schema | dbGaP open |
| **BrainSpan** | Developmental expression schema | Free |

### Clinical & Rare Disease

| Database | Schema Needed | License |
|----------|---------------|---------|
| **ClinGen** | Gene curation schema | CC0 |
| **DECIPHER (open)** | Open-access CNV annotation schema | Free open tier |
| **Monarch Initiative** | Knowledge graph schema | CC BY 3.0 |
| **DDG2P** | Gene-disease schema | Open |

### Biomarkers & Labs

| Database | Schema Needed | License |
|----------|---------------|---------|
| **MarkerDB 2.0** | Biomarker annotation schema | Free academic |
| **LOINC** | Code structure, hierarchy | Free |
| **CALIPER** | Reference interval schema | Free |

### Environmental & Cardio

| Database | Schema Needed | License |
|----------|---------------|---------|
| **CTD** | Chemical-gene-disease schema | Free |
| **MITOMAP** | Variant annotation schema | Free |
| **MitoCarta** | Gene annotation schema | Free |
| **CARDIoGRAMplusC4D** | GWAS summary format | Public |
| **GLGC** | Lipid GWAS schema | Public |

### Trials & Literature

| Database | Schema Needed | License |
|----------|---------------|---------|
| **ClinicalTrials.gov** | API response schema | Public domain |
| **PubMed** | E-utilities API schema | Public domain |
| **Wikidata (full dump)** | SPARQL + JSON dump schema - drugs, vitamins, genes, proteins, diseases, compounds, bio processes | CC0 |

---

## Tier 3 - Specialized (44 databases)

Lower priority - document when specific needs arise.

### Genetics Specialized

| Database | Schema Needed | Notes |
|----------|---------------|-------|
| gnomAD (main) | Variant API schema | Large dataset |
| dbVar | VCF structure | SV reference |
| ClinGen Dosage | HI/TS score schema | CNV interpretation |

### TCM Specialized

| Database | Schema Needed | Notes |
|----------|---------------|-------|
| TCMID 2.0 | Herb-compound-target schema | Backup to BATMAN |
| TCMBank | Interaction prediction format | AI predictions |
| CMAUP | Plant-compound mapping | Plants focus |
| NPASS | Activity data schema | Activity assays |

### Microbiome Specialized

| Database | Schema Needed | Notes |
|----------|---------------|-------|
| MASI | Interaction table schema | Microbe-drug |
| MDAD | Drug-microbe association | Disease focus |

### Cancer Specialized

| Database | Schema Needed | Notes |
|----------|---------------|-------|
| Cancer Gene Census | Gene annotation schema | Gene list |
| GDC/TCGA (open) | MAF format, open tier only | Open access data |

### Rare Disease Specialized

| Database | Schema Needed | Notes |
|----------|---------------|-------|
| OMIM | Entry structure | Gene-disease |

### Mental Health Specialized

| Database | Schema Needed | Notes |
|----------|---------------|-------|
| PubChem | Assay result schema | Compound assays |
| KEGG (neuro) | KGML schema | Pathways |

### Biomarkers Specialized

| Database | Schema Needed | Notes |
|----------|---------------|-------|
| Metabolomics Workbench | Study data schema | Research data |

### Environmental Specialized

| Database | Schema Needed | Notes |
|----------|---------------|-------|
| EPA CompTox | Chemical data schema | Toxicology |

---

## Already Documented (Verify Current)

These databases have existing schemas in `models/` - verify they're up to date:

| Database | Schema File | Status |
|----------|-------------|--------|
| BATMAN-TCM | tcm-data-models.md | KEEP - verify current |
| IMPPAT | ayurveda-data-models.md | KEEP - verify current |
| KampoDB | kampo-data-models.md | KEEP - verify current |
| PharmGKB | pharmgkb-cpic-data-models.md | KEEP - verify current |
| DrugBank | drugbank-chembl-data-models.md | KEEP - verify current |
| ChEMBL | drugbank-chembl-data-models.md | KEEP - verify current |
| Reactome | pathway-formats-detailed.md | KEEP - verify current |
| ClinVar | (needs schema doc) | KEEP - needs documentation |
| DSLD | western-herbal-data-models.md | KEEP - verify current |
| DGIdb | drug-target-interaction-models.md | KEEP - verify current |
| BindingDB | drug-target-interaction-models.md | KEEP - verify current |

---

## Action Items

1. **Document Tier 1 schemas** (16 databases) in `models/` directory
2. **Add download code** to DATA-INTEGRATION-GUIDE.md for Tier 1
3. **Verify existing schemas** are current (10 databases)
4. **Document Tier 2** (35 databases) after Tier 1 complete
5. **Keep Tier 3** (44 databases) in inventory for future needs

---

## Database Count Summary

| Category | Count | Source |
|----------|-------|--------|
| Tier 1 (MVP Critical) | 16 | New research |
| Tier 2 (High Value) | 37 | New research |
| Tier 3 (Specialized) | 43 | New research |
| Already Documented | 10 | Original research |
| **TOTAL DATABASES** | **~106** | To document |

---

## Cross-References

- **Existing schemas:** `models/*.md` (460KB, ~33 databases)
- **Download code:** `DATA-INTEGRATION-GUIDE.md`
- **Pipeline design:** `bulk-downloads/data-processing-pipeline.md`

---

## Download

| Resource | Method | URL |
|----------|--------|-----|
| dbNSFP | Academic License | https://sites.google.com/site/jpopgen/dbNSFP |
| ClinVar | FTP | https://ftp.ncbi.nih.gov/pub/clinvar/ |
| gnomAD | Cloud (AWS/GCP) | https://gnomad.broadinstitute.org/downloads |
| COCONUT | REST API | https://coconut.naturalproducts.net/ |
| BATMAN-TCM | REST API | https://batman.bi.a.u-tokyo.ac.jp/ |

**Access Requirements:** Most freely available; some require academic institutional access.

## Data Format

| Format | Description |
|--------|-------------|
| VCF | Variant Call Format (genomics standard) |
| TSV/CSV | Tabular data |
| JSON | API responses |
| SQL | Database dumps |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `database_id` | string | Unique database identifier | "clinvar-001" |
| `tier_level` | string | Priority tier for documentation | "tier1" |
| `schema_status` | string | Documentation completeness | "complete" |
| `license_type` | string | License category | "open" |
| `record_count` | number | Total records in database | 1000000 |

## Sample Data

### Example Database Schema Entry
```json
{
  "database_name": "ClinVar",
  "tier_level": "tier1",
  "schema_status": "documented",
  "fields": [
    {
      "name": "VariationID",
      "type": "integer",
      "description": "Unique variant identifier",
      "example": 123456
    },
    {
      "name": "ClinicalSignificance",
      "type": "enum",
      "values": ["pathogenic", "likely_pathogenic", "uncertain", "likely_benign", "benign"],
      "example": "pathogenic"
    }
  ],
  "license": "public_domain",
  "total_records": 2500000
}
```

## License

| Resource | License | Notes |
|----------|---------|-------|
| ClinVar | Public Domain | Freely usable |
| dbNSFP | Academic | Requires registration |
| gnomAD | Open Access | CC0 |
| COCONUT | CC0 | Public domain |
| BATMAN-TCM | CC BY-NC | Non-commercial |

## Data Set Size

| Metric | Value |
|--------|-------|
| Tier 1 databases | 16 |
| Tier 2 databases | 35 |
| Tier 3 databases | 44 |
| Already documented | 10 |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Schema Documentation | Formal description of a database's data structure and fields | Field definitions, types |
| Tier Classification | Priority ranking for implementation order | Tier 1 = MVP critical |
| Data Model | Conceptual representation of data entities and relationships | Entity-relationship diagram |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Tier 1 | MVP-critical databases requiring immediate documentation | 16 databases |
| Tier 2 | High-value databases for post-MVP implementation | 35 databases |
| Tier 3 | Specialized databases for future consideration | 44 databases |
| Schema Gap | Missing formal documentation for a database's structure | Undocumented fields |
| Field Mapping | Correspondence between fields in different databases | ID mapping |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | Data access |
| DB | Database | Data source |
| ETL | Extract, Transform, Load | Data pipeline |
| JSON | JavaScript Object Notation | Data format |
| KB | Kilobyte | Storage unit |
| MVP | Minimum Viable Product | Development phase |

---

*Updated: January 2026*
