---
id: compounds-drug-metabolism
title: Drug Metabolism and Supplement Interactions
category: shared
subcategory: drug-metabolism
tier: 1
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
databases:
  - pharmvar
  - supercyp
  - flockhart-table
  - curated-cyp-dataset
  - kegg-drug
  - pharmgkb
  - drugbank
  - cpic
  - natmed-pro
  - mskcc-herbs
  - stockleys
  - fdb
  - imgateway
  - t3db
  - hmdb
  - pk-db
  - swissadme
  - fda-bcs
  - fda-pgx-table
tags: [drug-metabolism, cyp450, supplements, interactions, pharmacokinetics, transporters]
---

# Drug Metabolism and Supplement Interaction Data Sources

**Document ID:** 43-53-DRUG-METABOLISM
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

Drug metabolism databases provide comprehensive coverage of CYP450 enzyme polymorphisms, drug-enzyme interactions, supplement-drug interactions, and transporter genetics. PharmVar serves as the gold standard for star allele nomenclature, while PharmGKB and CPIC provide clinical implementation guidelines. Integration of 25+ databases enables complete pharmacokinetic profiling from genotype through metabolizer status to personalized dosing recommendations.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary CYP allele nomenclature | PharmVar | Official star allele definitions used by PharmGKB/CPIC | Jan 2026 |
| Clinical PGx guidelines | PharmGKB + CPIC | Gold standard with EHR-ready implementation | Jan 2026 |
| Drug pathway integration | KEGG + DrugBank | Comprehensive pathway coverage with REST APIs | Jan 2026 |
| Supplement interactions | NatMed Pro + MSKCC | Industry-leading accuracy with free alternative | Jan 2026 |
| Transporter genetics | PharmGKB + CPIC | Curated SLCO1B1, ABCB1, ABCG2 annotations | Jan 2026 |
| PK modeling data | PK-DB | Full REST API with meta-analysis support | Jan 2026 |

---

[Full content from drug-metabolism.md original file]

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [../_index.md](../_index.md) | Parent index |
| [pharmaceuticals.md](./pharmaceuticals.md) | Drug-gene core data |
| [../../genetics/primary.md](../../genetics/primary.md) | dbSNP/ClinVar variant data |
| [../../pathways/primary.md](../../pathways/primary.md) | Drug pathway integration |
| [../../clinical/biomarkers-labs.md](../../clinical/biomarkers-labs.md) | Lab test integration |

---

## License

This document catalogs multiple databases with varying license terms:

| Database | License | Commercial Use | Attribution | Access |
|----------|---------|----------------|-------------|--------|
| PharmVar | Open Access | Yes | Citation | Open |
| SuperCYP | Academic use | Research only | Required | Open |
| Flockhart Table | Free educational use | Educational only | Required | Open |
| PharmGKB | CC BY-SA 4.0 | Yes (with attribution) | Required | Free (account required) |
| CPIC | CC0 (Public Domain) | Yes | None required | Open |
| DrugBank | CC BY-NC 4.0 (academic), Commercial license | Requires license | Required | Academic free |
| KEGG Drug | Academic API only, Subscription for FTP | Subscription required | Required | Academic API free |
| NatMed Pro | Commercial subscription | Subscription required | N/A | Subscription |
| MSKCC About Herbs | Free access | Yes | Citation | Open |
| Stockley's | Commercial subscription | Subscription required | N/A | Subscription |
| First Databank (FDB) | Commercial subscription | Subscription required | N/A | Enterprise |
| Integrative Medicine Gateway | Commercial subscription | Subscription required | N/A | Subscription |
| T3DB | Academic use | Research only | Required | Open |
| HMDB | Open Access | Yes | Citation | Open |
| PK-DB | Open Access | Yes | Citation | Open (REST API) |
| SwissADME | Free web access | Yes | Citation | Open |
| FDA BCS Database | Public Domain | Yes | None required | Open |
| FDA PGx Table | Public Domain | Yes | None required | Open |

**Key Considerations:**
- **Fully Open (Commercial OK):** CPIC, HMDB, PK-DB, SwissADME, FDA BCS, FDA PGx Table, MSKCC
- **Commercial Friendly with Attribution:** PharmGKB (CC BY-SA), PharmVar
- **Academic Only:** DrugBank (full data), SuperCYP, T3DB
- **Commercial Subscription Required:** NatMed Pro, Stockley's, FDB, KEGG Drug (FTP)

---

## Download

| Database | Method | URL/Command |
|----------|--------|-------------|
| **PharmGKB** | Bulk download | `https://www.pharmgkb.org/downloads` |
| **PharmVar** | Web | `https://www.pharmvar.org/` |
| **CPIC** | Download | `https://cpicpgx.org/guidelines/` |
| **DrugBank** | Academic download | `https://go.drugbank.com/releases/latest` (registration required) |
| **SuperCYP** | Web | `http://bioinformatics.charite.de/supercyp/` |
| **PK-DB** | API | `https://pk-db.com/` |

**Access Requirements:** PharmGKB (CC BY-SA 4.0) and CPIC (CC0) are freely available; DrugBank requires academic registration.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | TSV, CSV, JSON |
| Alternative | XML, SDF |
| Identifiers | PharmGKB ID, DrugBank ID |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `drug_id` | string | Primary drug identifier | "PA451906" |
| `gene` | string | Pharmacogene symbol | "CYP2D6" |
| `star_allele` | string | Haplotype designation | "*4" |
| `phenotype` | string | Metabolizer status | "Poor Metabolizer" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `metabolized_by` | Enzyme | N:M |
| `has_guideline` | CPIC Guideline | 1:N |

## Sample Data

### Example Drug-Gene Interaction
```json
{
  "drug": "codeine",
  "pharmgkb_id": "PA449088",
  "gene": "CYP2D6",
  "phenotype_recommendation": {
    "Poor Metabolizer": "Avoid codeine; use alternative analgesic",
    "Ultrarapid Metabolizer": "Avoid codeine; risk of toxicity"
  },
  "evidence_level": "1A"
}
```

### Sample Query Result
| drug | gene | metabolizer | recommendation |
|------|------|-------------|----------------|
| codeine | CYP2D6 | PM | Avoid use |
| tamoxifen | CYP2D6 | PM | Consider alternative |

---

## Data Set Size

| Metric | Value |
|--------|-------|
| PharmVar star alleles | Official nomenclature database |
| PharmGKB annotations | 20K+ gene-drug annotations |
| CPIC guidelines | 164 drugs with guidelines |
| DrugBank entries | 15K+ drugs |
| KEGG Drug | Pathway-integrated drug data |
| PK-DB pharmacokinetics | Full PK parameters |
| Total storage estimate | ~5-8 GB (combined sources) |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `CYP450` | Cytochrome P450 enzyme superfamily responsible for metabolizing ~75% of drugs | CYP3A4 metabolizes 50% of drugs |
| `substrate` | A drug metabolized by a specific enzyme | Codeine is a CYP2D6 substrate |
| `inhibitor` | A compound that reduces enzyme activity, slowing metabolism of substrates | Fluoxetine inhibits CYP2D6 |
| `inducer` | A compound that increases enzyme expression, accelerating metabolism | Rifampin induces CYP3A4 |
| `star allele` | Standardized haplotype nomenclature for pharmacogenes | CYP2D6*4 (null allele) |
| `metabolizer phenotype` | Predicted metabolic capacity based on genotype | Poor Metabolizer (PM) |
| `AUC` | Area Under the Curve - measure of total drug exposure over time | AUC increased 3-fold |
| `Cmax` | Maximum plasma concentration of a drug | Cmax = 150 ng/mL |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `PharmVar` | Pharmacogene Variation Consortium defining official star allele names | Allele nomenclature |
| `CPIC` | Clinical Pharmacogenetics Implementation Consortium providing dosing guidelines | Clinical implementation |
| `PharmGKB` | Pharmacogenomics Knowledge Base with gene-drug annotations | PGx evidence |
| `DPYD` | Dihydropyrimidine dehydrogenase - metabolizes fluoropyrimidines | Chemotherapy toxicity |
| `TPMT` | Thiopurine methyltransferase - metabolizes azathioprine/mercaptopurine | Immunosuppressant dosing |
| `UGT` | UDP-glucuronosyltransferase family - Phase II conjugation enzymes | Drug glucuronidation |
| `SLCO1B1` | Solute carrier transporter affecting statin hepatic uptake | Statin myopathy risk |
| `ABCB1` | ATP-binding cassette transporter (P-glycoprotein) affecting drug absorption | Drug efflux |
| `PK-DB` | Pharmacokinetics Database with quantitative PK parameters | PK modeling |
| `NatMed Pro` | Natural Medicines database for supplement interactions | Herb-drug interactions |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| CYP | Cytochrome P450 | Major Phase I enzyme family |
| PGx | Pharmacogenomics | Genetic effects on drug response |
| PK | Pharmacokinetics | Drug absorption, distribution, metabolism, excretion |
| PD | Pharmacodynamics | Drug effects on the body |
| ADME | Absorption, Distribution, Metabolism, Excretion | Core PK parameters |
| DDI | Drug-Drug Interaction | Clinically significant interactions |
| HDI | Herb-Drug Interaction | Supplement-medication interactions |
| PM | Poor Metabolizer | Little to no enzyme activity |
| IM | Intermediate Metabolizer | Reduced enzyme activity |
| NM | Normal Metabolizer | Typical enzyme activity |
| UM | Ultrarapid Metabolizer | Increased enzyme activity |
| EHR | Electronic Health Record | Clinical decision support target |

---

## Open Questions

- [ ] NatMed Pro licensing - subscription cost for production?
- [ ] PharmVar API - future development plans?
- [ ] KEGG commercial licensing - required for production?
- [ ] Transporter guidelines - CPIC expansion timeline?
- [ ] PK-DB coverage - sufficient for all major drugs?

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Initial catalog from data-sources-drug-metabolism.md |
