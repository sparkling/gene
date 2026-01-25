# 08. Literature & Knowledge - MVP Sources

## Summary

| Metric | Value |
|--------|-------|
| Sources Selected | 8 |
| Top Score | 26/27 |
| Key Identifiers | PMID, DOI, Wikidata QID, ORCID |

## Selected Sources

### 8.1 Scientific Literature

#### PubMed | Score: 25/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [1/3], Versioned [3/3], Size [3/3]
- **License**: CC0/Public Domain
- **Key IDs**: PMID (universal), DOI, PMC ID, MeSH terms, ORCID
- **Why Include**: THE biomedical literature database (36M+ citations), daily updates

#### OpenAlex | Score: 26/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [2/3], Versioned [3/3], Size [2/3]
- **License**: CC0
- **Key IDs**: OpenAlex ID, DOI, PMID, ORCID, ROR
- **Why Include**: Best open bibliometric resource (250M+ works), successor to Microsoft Academic Graph

### 8.2 Knowledge Bases

#### Wikidata | Score: 26/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [3/3], Versioned [2/3], Size [2/3]
- **License**: CC0
- **Key IDs**: QID, plus 100+ external DB properties (UniProt P352, Ensembl P594, OMIM P492, ChEMBL P592, PubMed P698, NCBI Gene P351, DrugBank P715, PubChem P662, HGNC P354, GO P686)
- **Why Include**: THE universal ID bridge, native RDF/SPARQL

### 8.3 Identifier Mapping

#### UniProt ID Mapping | Score: 24/27
- **Tier 1**: Identifiers [3/3], License [2/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [1/3], Versioned [3/3], Size [3/3]
- **License**: CC BY 4.0
- **Key IDs**: UniProt AC → 100+ databases
- **Why Include**: Critical for protein/gene ID resolution

#### PMC ID Converter | Score: 22/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [1/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [1/3], Size [3/3]
- **License**: Public Domain
- **Key IDs**: PMID ↔ PMCID ↔ DOI
- **Why Include**: Essential for literature ID harmonization

#### NCBI E-Link | Score: 21/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [1/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [3/3]
- **Tier 3**: RDF [0/3], Versioned [1/3], Size [3/3]
- **License**: Public Domain
- **Key IDs**: All NCBI database UIDs
- **Why Include**: Links across 40+ NCBI databases

### 8.4 Regulatory & Legal

#### ClinicalTrials.gov | Score: 24/27
- **Tier 1**: Identifiers [2/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [1/3], Versioned [2/3], Size [3/3]
- **License**: Public Domain
- **Key IDs**: NCT Number, MeSH conditions, drug names
- **Why Include**: THE clinical trials registry (500K+ studies)

#### FDA OpenFDA | Score: 23/27
- **Tier 1**: Identifiers [3/3], License [3/3], Bulk [3/3]
- **Tier 2**: Canonical [3/3], Active [3/3], XRefs [2/3]
- **Tier 3**: RDF [0/3], Versioned [2/3], Size [2/3]
- **License**: Public Domain
- **Key IDs**: NDC Code, NDA/ANDA, RxNorm, UNII
- **Why Include**: Essential for drug safety and pharmacovigilance

## Excluded Sources

| Source | Score | Reason |
|--------|-------|--------|
| PubMed Central | 20/27 | Large size, mixed licensing |
| Europe PMC | 22/27 | Overlaps with PubMed |
| Semantic Scholar | 18/27 | ODC-BY adds compliance burden |
| Wikipedia | 17/27 | CC BY-SA ShareAlike; Wikidata sufficient |

## Integration Points

```
PubMed (PMID) ←→ PMC ID Converter ←→ DOI
       ↓
   OpenAlex (bibliometrics)
       ↓
   Wikidata (universal hub)
       ↓
UniProt ID Mapping ←→ NCBI E-Link
       ↓
ClinicalTrials.gov ←→ FDA OpenFDA
```

## Wikidata as Universal Bridge

Wikidata connects virtually all MVP sources via its properties:

| Property | ID Type | Links To |
|----------|---------|----------|
| P352 | UniProt | Proteins |
| P351 | NCBI Gene | Genes |
| P594 | Ensembl | Genes |
| P698 | PubMed | Literature |
| P492 | OMIM | Diseases |
| P662 | PubChem | Compounds |
| P592 | ChEMBL | Bioactivity |
| P715 | DrugBank | Drugs |
| P686 | GO | Functions |
| P354 | HGNC | Gene symbols |

## Data Size Estimates

| Source | Size | Notes |
|--------|------|-------|
| PubMed | ~350GB | FTP baseline |
| OpenAlex | ~500GB | S3 snapshots |
| Wikidata | ~100GB | Compressed JSON |
| UniProt ID Mapping | ~30GB | idmapping.dat.gz |
| ClinicalTrials.gov | ~5-20GB | Manageable |
| FDA OpenFDA | ~20GB | Can subset |
