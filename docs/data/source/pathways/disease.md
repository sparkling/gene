# Disease and Phenotype Databases

**Document ID:** 43-43-PATHWAYS-DISEASE
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../index.md](./../index.md)

---

## TL;DR

Disease and phenotype databases link genes, variants, and biological pathways to human diseases. Primary sources include DisGeNET (1.1M+ gene-disease associations), OMIM (25K+ Mendelian disorders), and Human Phenotype Ontology (HPO, 13K+ phenotype terms). These databases enable genotype-phenotype correlation, disease mechanism understanding, and clinical variant interpretation.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary gene-disease source | DisGeNET | Largest curated collection, integrates multiple sources | Jan 2026 |
| Mendelian disease reference | OMIM | Gold standard for genetic disorders, comprehensive | Jan 2026 |
| Phenotype ontology | HPO | Standard vocabulary, cross-species integration via Monarch | Jan 2026 |
| Disease ontology | MONDO | Unified ontology with precise equivalence axioms | Jan 2026 |
| Rare disease integration | Orphanet/ORDO | ELIXIR Core Resource, comprehensive rare disease coverage | Jan 2026 |
| Clinical variant interpretation | ClinGen + ClinVar | NIH-funded authoritative sources | Jan 2026 |

---

## Database Catalog

### DisGeNET

| Field | Value |
|-------|-------|
| **URL** | https://www.disgenet.org/ |
| **Content** | Gene-disease associations, variant-disease associations |
| **Records** | 1.1M+ gene-disease associations, 600K+ variant-disease associations |
| **License** | Free for academic research; commercial license required |
| **API** | REST API (registration required) |
| **Update Frequency** | Bi-annual |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~2 GB (full dataset) |

**Data Types:**
- Gene-disease associations with evidence scores
- Variant-disease associations
- Disease-disease associations
- Gene-gene associations (shared diseases)
- Integration from 20+ source databases

**Score Types:**
| Score | Description |
|-------|-------------|
| GDA Score | Gene-Disease Association strength (0-1) |
| VDA Score | Variant-Disease Association strength (0-1) |
| EI (Evidence Index) | Ratio of supporting vs total publications |
| DSI (Disease Specificity Index) | How specific gene is to disease |
| DPI (Disease Pleiotropy Index) | Number of disease classes associated |

**Download Formats:** TSV, SQLite database, RDF

---

### OMIM (Online Mendelian Inheritance in Man)

| Field | Value |
|-------|-------|
| **URL** | https://www.omim.org/ |
| **Content** | Mendelian disorders and genes |
| **Records** | 25,000+ entries (genes, phenotypes, gene-phenotype relationships) |
| **License** | Johns Hopkins University; requires registration for API/download |
| **API** | REST API (registration required, rate limited) |
| **Update Frequency** | Daily |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |

**Entry Types:**
| MIM Number Prefix | Type |
|-------------------|------|
| * (asterisk) | Gene description |
| + (plus) | Gene + phenotype description |
| # (number sign) | Phenotype, molecular basis known |
| % (percent) | Phenotype, molecular basis unknown |
| None | Phenotype or locus |
| ^ (caret) | Entry removed/moved |

**Key Fields:**
- Gene symbols (HGNC approved)
- Clinical synopsis
- Allelic variants
- Inheritance patterns
- Chromosomal location
- Literature references

---

### Human Phenotype Ontology (HPO)

| Field | Value |
|-------|-------|
| **URL** | https://hpo.jax.org/ |
| **Content** | Standardized phenotypic abnormality vocabulary |
| **Records** | 13,000+ phenotype terms, 156,000+ disease annotations |
| **License** | Open access with attribution |
| **API** | REST API (free) |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~100 MB |

**Ontology Structure:**
```
Phenotypic abnormality (HP:0000118)
├── Abnormality of the nervous system
│   ├── Abnormality of nervous system physiology
│   │   ├── Seizure
│   │   ├── Intellectual disability
│   │   └── Movement abnormality
│   └── Abnormality of nervous system morphology
├── Abnormality of the cardiovascular system
├── Abnormality of metabolism/homeostasis
└── [14 other top-level categories]
```

**Download Formats:**
| Format | URL |
|--------|-----|
| OWL | http://purl.obolibrary.org/obo/hp.owl |
| OBO | http://purl.obolibrary.org/obo/hp.obo |
| JSON | GitHub releases |

**Languages:** English, German, Spanish, French, Italian, Dutch, Portuguese, Turkish, Japanese, Chinese, and more

---

### MONDO Disease Ontology

| Field | Value |
|-------|-------|
| **URL** | https://mondo.monarchinitiative.org/ |
| **Content** | Unified disease ontology |
| **Records** | 25,938 disease concepts (22,977 human diseases) |
| **License** | CC-BY-4.0 |
| **API** | Via Monarch Initiative API |
| **Update Frequency** | Monthly |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~150 MB |

**Key Features:**
- Precise 1:1 equivalence axioms to OMIM, Orphanet, NCIT
- Validated by OWL reasoning
- Integrates disease definitions across sources
- Maintained by Monarch Initiative

**Download Formats:**
| Format | File | Description |
|--------|------|-------------|
| OWL | mondo-with-equivalents.owl | With equivalence axioms |
| OBO | mondo.obo | Using xrefs for linking |
| JSON | mondo-with-equivalents.json | Equivalent to OWL |

---

### Orphanet / ORDO (Orphanet Rare Disease Ontology)

| Field | Value |
|-------|-------|
| **URL** | https://www.orpha.net/ |
| **Content** | Rare diseases, orphan drugs, expert centres |
| **Records** | 6,528 rare diseases, 4,512 linked genes |
| **License** | CC-BY-4.0 |
| **API** | REST API (https://api.orphadata.com/) |
| **Update Frequency** | Bi-annual (July and December) |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~500 MB |

**Data Types:**
| Dataset | Description |
|---------|-------------|
| ORDO | Structured vocabulary for rare diseases (44-50 MB OWL) |
| Disease Classifications | Hierarchical organization |
| Gene-Disease Associations | Gene linkages to rare diseases |
| HPO-ORDO Module (HOOM) | Phenotype mappings |
| Epidemiological Data | Prevalence and incidence figures |
| Expert Centres | 8,722 specialized healthcare centers |
| Diagnostic Tests | 36,595 available tests |

**Cross-References:** OMIM, MONDO, UniProtKB, HGNC, Ensembl, Reactome, IUPHAR

---

### ClinVar

| Field | Value |
|-------|-------|
| **URL** | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **Content** | Clinical variant interpretations |
| **Records** | 2M+ submissions, 1.5M+ unique variants |
| **License** | Public Domain (NCBI) |
| **API** | E-utilities REST API |
| **Update Frequency** | Daily |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~5 GB |

**Clinical Significance Categories:**
| Category | Description |
|----------|-------------|
| Pathogenic | Causes disease |
| Likely pathogenic | Probably causes disease |
| Uncertain significance (VUS) | Insufficient evidence |
| Likely benign | Probably does not cause disease |
| Benign | Does not cause disease |

**Download Options:**
| File | Description |
|------|-------------|
| variant_summary.txt.gz | All variants with interpretations |
| clinvar.vcf.gz | VCF format (GRCh37/38) |
| clinvar_*.xml.gz | Full XML export |
| submission_summary.txt.gz | Individual lab submissions |

---

### ClinGen (Clinical Genome Resource)

| Field | Value |
|-------|-------|
| **URL** | https://clinicalgenome.org/ |
| **Content** | Gene-disease validity, variant pathogenicity, clinical actionability |
| **Records** | 3,289 gene-disease validity curations, 11,133 variant curations |
| **License** | Open access with citation |
| **API** | Multiple APIs (Evidence Repository, Allele Registry, etc.) |
| **Update Frequency** | Daily |
| **Priority** | Tier 1 (MVP) |
| **Storage Estimate** | ~1 GB |

**Data Types:**
| Curation Type | Description | Count |
|---------------|-------------|-------|
| Gene-Disease Validity | Strength of gene-disease relationships | 3,289 |
| Variant Pathogenicity | ACMG/AMP classifications | 11,133 |
| Dosage Sensitivity | Copy number effects | 1,557 genes |
| Clinical Actionability | Intervention evaluations | 447 pairs |

**Gene-Disease Validity Classifications:**
| Classification | Description |
|----------------|-------------|
| Definitive | Conclusive evidence |
| Strong | Strong evidence supporting association |
| Moderate | Moderate evidence |
| Limited | Limited evidence |
| Disputed | Conflicting evidence |
| Refuted | Evidence against association |
| No Known Disease Relationship | No evidence found |

**Download Formats:** CSV, TSV, BED, JSON (nested/flat)

---

### Disease Ontology (DO)

| Field | Value |
|-------|-------|
| **URL** | https://disease-ontology.org/ |
| **Content** | Standardized disease terminology |
| **Records** | 11,000+ disease terms |
| **License** | CC0 1.0 (Public Domain) |
| **API** | Via OBO Foundry / OLS |
| **Update Frequency** | Monthly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 MB |

**Key Features:**
- Integrates with MeSH, ICD, NCI Thesaurus, OMIM
- OBO Foundry compliant
- Used by major databases for disease annotation

**Cross-References:**
| Resource | Mapping Type |
|----------|--------------|
| MeSH | Disease terms |
| ICD-9/10 | Clinical codes |
| OMIM | Mendelian disorders |
| NCI Thesaurus | Cancer terminology |
| UMLS | Unified medical language |

---

### MalaCards

| Field | Value |
|-------|-------|
| **URL** | https://www.malacards.org/ |
| **Content** | Human disease database (from GeneCards suite) |
| **Records** | 22,000+ disease entries |
| **License** | Free for academic; commercial license required |
| **API** | Web services (registration required) |
| **Update Frequency** | Quarterly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~3 GB |

**Data Types:**
- Disease-gene associations
- Disease-drug relationships
- Disease-publication links
- Disease-pathway associations
- Anatomical classifications

**Integrated Sources:** OMIM, Orphanet, ClinVar, HPO, UniProt, Reactome, KEGG

---

### Monarch Initiative Knowledge Graph

| Field | Value |
|-------|-------|
| **URL** | https://monarchinitiative.org/ |
| **Content** | Cross-species phenotype-genotype integration |
| **Records** | 33 integrated data sources |
| **License** | Open access (individual sources may vary) |
| **API** | REST API v3 (https://api-v3.monarchinitiative.org/) |
| **Update Frequency** | Continuous |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~10 GB (full KG) |

**Integrated Sources Include:**
- Model organism databases: MGD, ZFIN, WormBase, FlyBase
- Human disease: OMIM, Orphanet
- Phenotypes: HPO, Gene Ontology
- Pathways: Reactome
- Interactions: STRING

**Available Datasets:**
| Dataset | Description |
|---------|-------------|
| Monarch KG | Main knowledge graph |
| Exomiser | Variant prioritization resources |
| Phenologs | Phenotype ortholog data |
| Semantic Similarity | Phenotype similarity metrics |
| UPHENO2 | Unified Phenotype Ontology |

---

### DECIPHER

| Field | Value |
|-------|-------|
| **URL** | https://www.deciphergenomics.org/ |
| **Content** | Clinical genomic variants and phenotypes |
| **Records** | 51,894 open-access patients, 250+ contributing centers |
| **License** | Data Access Agreement required for bulk access |
| **API** | Deposition API (with API keys) |
| **Update Frequency** | Frequent |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~2 GB |

**Data Types:**
- Copy Number Variants (CNVs)
- Sequence Variants (SNVs, indels)
- Phenotypes (HPO-linked)
- CNV Syndromes

**Access Tiers:**
| Tier | Requirements |
|------|--------------|
| Open Access | Explicit patient consent (public) |
| Bulk Research | Data Access Agreement |
| Display Agreement | Genome browser integration |

---

### GWAS Catalog

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/gwas/ |
| **Content** | Published genome-wide association studies |
| **Records** | 6,000+ publications, 500K+ associations |
| **License** | Open access |
| **API** | REST API + Summary Statistics API |
| **Update Frequency** | Weekly |
| **Priority** | Tier 2 |
| **Storage Estimate** | ~50 GB (with summary statistics) |

**Data Types:**
- SNP-trait associations
- Summary statistics (per study)
- Study metadata
- Ancestral breakdowns

**Download Options:**
| Format | Description |
|--------|-------------|
| TSV | Association data |
| OWL/RDF | Ontology format |
| VCF | Variant format |
| Summary Stats | Full GWAS results per study |

---

## Summary Comparison

### Database Focus Areas

| Database | Primary Focus | Records | License |
|----------|---------------|---------|---------|
| DisGeNET | Gene-disease associations | 1.1M+ associations | Academic free |
| OMIM | Mendelian disorders | 25K+ entries | Registration required |
| HPO | Phenotype ontology | 13K+ terms | Open |
| MONDO | Disease ontology | 26K+ concepts | CC-BY-4.0 |
| Orphanet | Rare diseases | 6,500+ diseases | CC-BY-4.0 |
| ClinVar | Variant interpretations | 2M+ submissions | Public Domain |
| ClinGen | Clinical curation | 14K+ curations | Open |
| Disease Ontology | Disease terminology | 11K+ terms | CC0 |
| MalaCards | Disease encyclopedia | 22K+ diseases | Academic free |
| Monarch | Cross-species integration | 33 sources | Open |
| DECIPHER | Clinical variants | 52K+ patients | Agreement required |
| GWAS Catalog | Association studies | 500K+ associations | Open |

### Access Methods

| Database | Web | REST API | SPARQL | Bulk Download |
|----------|-----|----------|--------|---------------|
| DisGeNET | Yes | Yes | Yes | Yes |
| OMIM | Yes | Yes | No | Yes (license) |
| HPO | Yes | Yes | No | Yes |
| MONDO | Yes | Yes | No | Yes |
| Orphanet | Yes | Yes | No | Yes |
| ClinVar | Yes | Yes | No | Yes |
| ClinGen | Yes | Yes | No | Yes |
| Disease Ontology | Yes | Yes | No | Yes |
| MalaCards | Yes | Limited | No | No |
| Monarch | Yes | Yes | No | Yes |
| DECIPHER | Yes | Yes | No | Agreement |
| GWAS Catalog | Yes | Yes | No | Yes |

### Update Frequency

| Database | Frequency |
|----------|-----------|
| ClinVar | Daily |
| ClinGen | Daily |
| GWAS Catalog | Weekly |
| HPO | Monthly |
| MONDO | Monthly |
| Disease Ontology | Monthly |
| DisGeNET | Bi-annual |
| Orphanet | Bi-annual |
| MalaCards | Quarterly |
| OMIM | Daily |

---

## Integration Recommendations

### Priority Integration Order

| Priority | Database | Rationale |
|----------|----------|-----------|
| 1 | HPO | Foundation for phenotype annotation |
| 2 | MONDO | Unified disease identifiers |
| 3 | ClinVar | Clinical variant gold standard |
| 4 | DisGeNET | Comprehensive gene-disease links |
| 5 | ClinGen | Expert-curated validity |
| 6 | OMIM | Mendelian disease reference |
| 7 | Orphanet | Rare disease coverage |
| 8 | GWAS Catalog | Population associations |

### Ontology Integration Strategy

1. **Load base ontologies:**
   - HPO for phenotypes
   - MONDO for diseases
   - ORDO for rare diseases

2. **Map equivalences:**
   - MONDO provides precise mappings to OMIM, Orphanet
   - Use OWL reasoning for inference

3. **Gene-disease layer:**
   - ClinGen for validated associations
   - DisGeNET for comprehensive coverage
   - Cross-reference with OMIM

4. **Variant layer:**
   - ClinVar for clinical interpretations
   - DECIPHER for structural variants
   - GWAS Catalog for common variants

### Identifier Mapping

| Database | Primary ID | Cross-References |
|----------|------------|------------------|
| HPO | HP:xxxxxxx | OMIM, Orphanet, MONDO |
| MONDO | MONDO:xxxxxxx | OMIM, Orphanet, DO, NCIT |
| OMIM | MIM:xxxxxx | HPO, Orphanet, ClinVar |
| Orphanet | ORPHA:xxxxx | OMIM, HPO, MONDO |
| DisGeNET | UMLS CUI | OMIM, DO, HPO |
| ClinVar | RCV/VCV | dbSNP, ClinGen, OMIM |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [primary.md](./primary.md) | Pathway-disease connections |
| [primary.md](./../genetics/primary.md) | Gene-disease variant data |
| [rare.md](./../diseases/rare.md) | Detailed rare disease sources |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial disease database catalog |
