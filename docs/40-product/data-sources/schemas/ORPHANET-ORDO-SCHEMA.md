# Orphanet/ORDO Rare Disease Schema

**Document ID:** ORPHANET-ORDO-SCHEMA
**Source:** Orphanet (https://www.orpha.net/) & ORDO Ontology (EBI OLS)
**Last Updated:** January 2026
**ORDO Version:** 4.7

---

## Database Statistics (January 2026)

### Orphanet Database

| Metric | Count |
|--------|-------|
| Rare Diseases | 6,528 |
| Linked Genes | 4,512 |
| Diagnostic Tests | 36,595 |
| Expert Centres | 8,722 |
| Languages | 8 (EN, CZ, NL, FR, DE, IT, ES, PL) |

### ORDO Ontology (via EBI OLS)

| Metric | Count |
|--------|-------|
| Total Terms | 15,788 |
| Properties | 35 |
| Individuals | 0 |
| Version | 4.7 |
| Last Updated | January 19, 2026 |

---

## Data Sources

| Resource | URL | Format | License |
|----------|-----|--------|---------|
| Orphanet Portal | https://www.orpha.net/ | Web | CC BY 4.0 |
| Orphadata Science | https://sciences.orphadata.com/ | XML, JSON | CC BY 4.0 |
| ORDO Ontology | https://www.ebi.ac.uk/ols4/ontologies/ordo | OWL | CC BY 4.0 |
| GitHub Archive | github.com/Orphanet | XML | CC BY 4.0 |

---

## Orphadata Products (December 2025)

### Available Datasets

| Product | File | Entries | Size | Description |
|---------|------|---------|------|-------------|
| Genes | en_product6.xml | 4,128 | 21.77 MB | Gene-disease associations |
| Phenotypes | en_product4.xml | 4,337 | 46+ MB | HPO-linked disease phenotypes |
| Disorders | en_product1.xml | 6,528 | Variable | Disease classifications |
| ORDO | ordo.owl | 15,788 | ~50 MB | Full ontology |
| HOOM | hoom.owl | Variable | Variable | HPO-ORDO Ontological Module |
| Orphapackets | JSON | Variable | Variable | Structured disease packets |

### Schema Documentation

Each dataset includes:
- XSD (XML Schema Definition)
- PDF description document
- JPEG visual diagram

---

## ORDO Ontology Schema

### Term Structure

```json
{
  "iri": "http://www.orpha.net/ORDO/Orphanet_90033",
  "label": "Warm autoimmune hemolytic anemia",
  "description": "Most prevalent form, defined by warm autoantibodies against red blood cells",
  "synonyms": ["WAIHA", "Warm antibody autoimmune hemolytic anemia"],
  "oboId": "Orphanet:90033",
  "shortForm": "ORDO:90033",
  "isDefiningOntology": true,
  "hasChildren": true,
  "isObsolete": false
}
```

### Key Properties

| Property | Description |
|----------|-------------|
| `iri` | Full URI identifier |
| `label` | Primary disease name |
| `description` | Clinical definition |
| `synonyms` | Alternative names |
| `oboId` | OBO format ID |
| `shortForm` | Short identifier |
| `hasChildren` | Has subtypes flag |
| `isObsolete` | Deprecation status |

---

### Disease Term Schema

```json
{
  "orpha_code": "90033",
  "name": "Warm autoimmune hemolytic anemia",
  "definition": "A form of autoimmune hemolytic anemia characterized by warm autoantibodies against red blood cells",
  "inheritance_pattern": "autosomal dominant",
  "age_of_onset": "All ages",
  "prevalence": {
    "geographic_area": "Europe",
    "prevalence_class": "1-9 / 100,000",
    "prevalence_type": "Point prevalence"
  },
  "cross_references": {
    "ICD-10": "D59.1",
    "ICD-11": "3A20.0",
    "UMLS": "C0272118",
    "MeSH": "D000744",
    "MedDRA": "10003553"
  }
}
```

### Prevalence Classes

| Class | Range |
|-------|-------|
| Unknown | No data available |
| <1 / 1,000,000 | Ultra-rare |
| 1-9 / 1,000,000 | Very rare |
| 1-9 / 100,000 | Rare |
| 1-5 / 10,000 | Less common |
| >1 / 1,000 | Common |

---

### Gene-Disease Association Schema

```json
{
  "orpha_code": "90033",
  "disease_name": "Warm autoimmune hemolytic anemia",
  "gene_associations": [
    {
      "gene_symbol": "TNFRSF13B",
      "gene_name": "TNF receptor superfamily member 13B",
      "hgnc_id": "HGNC:11924",
      "entrez_id": "23495",
      "ensembl_id": "ENSG00000240505",
      "association_type": "Disease-causing germline mutation(s) in",
      "association_status": "Assessed"
    }
  ]
}
```

### Association Types

| Type | Description |
|------|-------------|
| Disease-causing germline mutation(s) in | Causal germline variants |
| Disease-causing somatic mutation(s) in | Causal somatic variants |
| Major susceptibility factor in | Increased risk factor |
| Modifying germline mutation in | Disease modifier |
| Part of a fusion gene in | Fusion gene involvement |
| Role in the phenotype of | Phenotype contribution |
| Candidate gene tested in | Under investigation |

---

### Phenotype (HPO) Annotation Schema

```json
{
  "orpha_code": "90033",
  "disease_name": "Warm autoimmune hemolytic anemia",
  "hpo_annotations": [
    {
      "hpo_id": "HP:0001878",
      "hpo_term": "Hemolytic anemia",
      "frequency": "Very frequent (99-80%)",
      "diagnostic_criteria": true
    },
    {
      "hpo_id": "HP:0001744",
      "hpo_term": "Splenomegaly",
      "frequency": "Frequent (79-30%)",
      "diagnostic_criteria": false
    }
  ]
}
```

### Frequency Classes

| Class | Range |
|-------|-------|
| Obligate (100%) | Always present |
| Very frequent | 99-80% |
| Frequent | 79-30% |
| Occasional | 29-5% |
| Very rare | <4-1% |
| Excluded | 0% |

---

## Hierarchical Structure

### Disease Classification Levels

```
Rare Diseases
├── Rare genetic disease
│   ├── Rare developmental defect during embryogenesis
│   │   ├── Disorder of sex development
│   │   ├── Orofacial clefting syndrome
│   │   └── Multiple congenital anomalies/dysplastic syndrome
│   ├── Rare neurologic disease
│   │   ├── Rare movement disorder
│   │   └── Rare neurodegenerative disease
│   └── Rare inborn errors of metabolism
│       ├── Disorder of amino acid metabolism
│       └── Disorder of carbohydrate metabolism
├── Rare hematologic disease
│   ├── Rare acquired hemolytic anemia
│   │   ├── Autoimmune hemolytic anemia
│   │   │   ├── Warm autoimmune hemolytic anemia (ORDO:90033)
│   │   │   ├── Cold agglutinin disease
│   │   │   └── Mixed-type autoimmune hemolytic anemia
│   │   └── Drug-induced hemolytic anemia
│   └── Rare inherited hemolytic anemia
└── Rare oncologic disease
```

---

## Cross-Reference Systems

### External Database Links

| System | Prefix | Example |
|--------|--------|---------|
| OMIM | OMIM: | OMIM:614470 |
| ICD-10 | ICD-10: | ICD-10:D59.1 |
| ICD-11 | ICD-11: | ICD-11:3A20.0 |
| MeSH | MSH: | MSH:D000744 |
| UMLS | UMLS: | UMLS:C0272118 |
| MedDRA | MedDRA: | MedDRA:10003553 |
| MONDO | MONDO: | MONDO:0015893 |
| GARD | GARD: | GARD:5864 |
| SNOMED CT | SCTID: | SCTID:371096002 |

---

## XML Schema Example (Orphadata)

### Gene-Disease Product (en_product6.xml)

```xml
<?xml version="1.0" encoding="UTF-8"?>
<JDBOR>
  <DisorderList>
    <Disorder id="90033">
      <OrphaCode>90033</OrphaCode>
      <Name lang="en">Warm autoimmune hemolytic anemia</Name>
      <DisorderGeneAssociationList count="2">
        <DisorderGeneAssociation>
          <Gene id="14564">
            <Symbol>TNFRSF13B</Symbol>
            <Name lang="en">TNF receptor superfamily member 13B</Name>
            <SynonymList count="2">
              <Synonym lang="en">TACI</Synonym>
              <Synonym lang="en">CD267</Synonym>
            </SynonymList>
            <ExternalReferenceList>
              <ExternalReference>
                <Source>HGNC</Source>
                <Reference>11924</Reference>
              </ExternalReference>
              <ExternalReference>
                <Source>Ensembl</Source>
                <Reference>ENSG00000240505</Reference>
              </ExternalReference>
            </ExternalReferenceList>
          </Gene>
          <DisorderGeneAssociationType id="17949">
            <Name lang="en">Disease-causing germline mutation(s) in</Name>
          </DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus id="17991">
            <Name lang="en">Assessed</Name>
          </DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
    </Disorder>
  </DisorderList>
</JDBOR>
```

---

## API Access

### SPARQL Endpoint (via EBI OLS)
```sparql
PREFIX ordo: <http://www.orpha.net/ORDO/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT ?disease ?label ?parent
WHERE {
  ?disease rdfs:subClassOf ?parent .
  ?disease rdfs:label ?label .
  FILTER(STRSTARTS(STR(?disease), "http://www.orpha.net/ORDO/"))
}
LIMIT 100
```

### REST API (EBI OLS4)

```bash
# Get ORDO ontology info
curl "https://www.ebi.ac.uk/ols4/api/ontologies/ordo"

# Get disease term details
curl "https://www.ebi.ac.uk/ols4/api/ontologies/ordo/terms/http%3A%2F%2Fwww.orpha.net%2FORDO%2FOrphanet_90033"

# Search for disease
curl "https://www.ebi.ac.uk/ols4/api/search?q=hemolytic%20anemia&ontology=ordo"
```

### Orphadata Downloads

```bash
# Download gene-disease associations
wget https://www.orphadata.com/data/xml/en_product6.xml

# Download phenotype annotations
wget https://www.orphadata.com/data/xml/en_product4.xml

# Download ORDO ontology
wget http://www.orphadata.org/data/ORDO/ordo.owl
```

---

## Integration Notes

### OMIM Cross-References

Orphanet maintains mappings to OMIM for:
- Disease entries (phenotype MIM numbers)
- Gene entries (gene MIM numbers)
- Note: Not all rare diseases have OMIM entries

### HPO Integration (HOOM)

The HPO-ORDO Ontological Module provides:
- Disease-to-phenotype mappings
- Frequency annotations
- Diagnostic criteria flags
- Validated clinical features

### MONDO Alignment

MONDO Disease Ontology provides:
- Precise equivalence axioms to ORDO
- Additional cross-database mappings
- Unified disease classification

---

## Data Quality

### Expert Curation
- All disease entries expert-reviewed
- Gene-disease associations from literature
- Prevalence data from epidemiological studies
- Regular updates (bi-annual major releases)

### Recognition Status

| Designation | Description |
|-------------|-------------|
| ELIXIR Core Data Resource | Recognized as essential infrastructure |
| Global Core Biodata Resource | International recognition |
| IRDiRC Recognized Resource | Rare disease community standard |

---

## Use Cases

### Clinical Genomics
- Variant interpretation for rare diseases
- Differential diagnosis support
- Gene panel design

### Research
- Rare disease cohort identification
- Phenotype-genotype correlation studies
- Drug repurposing targets

### Public Health
- Rare disease epidemiology
- Healthcare resource planning
- Patient registry development

---

## License

**CC BY 4.0 (Creative Commons Attribution 4.0 International)**

Requirements:
- Attribution to Orphanet
- Indicate if changes were made
- Link to license

Permitted:
- Commercial use
- Derivative works
- Distribution

---

## References

1. Orphanet: https://www.orpha.net/

2. Orphadata Science: https://sciences.orphadata.com/

3. ORDO in OLS: https://www.ebi.ac.uk/ols4/ontologies/ordo

4. Rath A, et al. (2012). Representation of rare diseases in health information systems: the Orphanet approach to serve a wide range of end users. Hum Mutat. 33(5):803-8. https://doi.org/10.1002/humu.22078

5. Vasant D, et al. (2014). ORDO: An Ontology Connecting Rare Disease, Epidemiology and Genetic Data. In Proceedings of SWAT4LS.

6. HPO-ORDO: https://github.com/obophenotype/hoom
